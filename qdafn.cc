#include "qdafn.h"

// -----------------------------------------------------------------------------
//  Ziggurat Method standard normal pseudorandom number generator code from 
//  George Marsaglia and Wai Wan Tsang (2000).
// 
//  "The Ziggurat Method for Generating Random Variables". Journal of
//  Statistical Software 5 (8).
// -----------------------------------------------------------------------------
unsigned long jsr = 123456789;
unsigned long jz;
long hz;
unsigned long iz, kn[128], ke[256];
float wn[128],fn[128], we[256],fe[256];

// -----------------------------------------------------------------------------
//  nfix() generates variates from the residue when rejection in RNOR occurs
// -----------------------------------------------------------------------------
float nfix()
{
	const float r = 3.442620;		// The start of the right tail
	static float x, y;

	while (1) {
		x = hz * wn[iz];
		if (iz == 0) {				// iz == 0, handles the base strip
			do {
				x = -log(UNI)*0.2904764;
				y = -log(UNI);
			} while (y+y < x*x);

			return (hz > 0) ? r+x : -r-x;
		}
									// iz > 0, handle the wedges of other strips
		if (fn[iz] + UNI*(fn[iz-1]-fn[iz]) < exp(-0.5*x*x)) {
			return x;
		}

		hz = SHR3;					// initiate, try to exit loop
		iz = hz & 127;
		if (abs(hz) < kn[iz]) return hz*wn[iz];
	}
}

// -----------------------------------------------------------------------------
//  efix() generates variates from the residue when rejection in REXP occurs
// -----------------------------------------------------------------------------
float efix()
{
	float x;

	while (1) {
		if (iz == 0) return 7.69711F-log(UNI); // iz == 0

		x = jz * we[iz];
		if (fe[iz] + UNI*(fe[iz-1]-fe[iz]) < exp(-x)) return x;

		jz = SHR3;					// initiate, try to exit loop
		iz = (jz & 255);
		if (jz < ke[iz]) return jz*we[iz];
	}
}

// -----------------------------------------------------------------------------
void zigset(						// set the seed and create the tables
	unsigned long jsrseed)				// new seed
{
	const double m1 = 2147483648.0;
	const double m2 = 4294967296.0;

	double dn = 3.442619855899;
	double tn = dn;
	double vn = 9.91256303526217e-3;
	double q;

	double de = 7.697117470131487;
	double te = de;
	double ve = 3.949659822581572e-3;
	int i;

	jsr ^= jsrseed;

	// -------------------------------------------------------------------------
	//  set up tables for RNOR
	// -------------------------------------------------------------------------
	q = vn / exp(-0.5 * dn * dn);
	kn[0] = (dn / q) * m1;
	kn[1] = 0;

	wn[0] = q / m1;
	wn[127] = dn / m1;

	fn[0] = 1.0;
	fn[127] = exp(-0.5 * dn * dn);

	for (i = 126; i >= 1; i--) {
		dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
		kn[i+1] = (dn / tn) * m1;
		tn = dn;
		fn[i] = exp(-0.5 * dn * dn);
		wn[i] = dn / m1;
	}

	// -------------------------------------------------------------------------
	//  set up tables for REXP
	// -------------------------------------------------------------------------
	q = ve / exp(-de);
	ke[0] = (de / q) * m2;
	ke[1] = 0;

	we[0] = q / m2;
	we[255] = de / m2;

	fe[0] = 1.0;
	fe[255] = exp(-de);

	for (i = 254; i >= 1; i--) {
		de = -log(ve / de + exp(-de));
		ke[i+1] = (de / te) * m2;
		te = de;
		fe[i] = exp(-de);
		we[i] = de / m2;
	}
}

// -----------------------------------------------------------------------------
QDAFN::QDAFN()						// default constructor
{
	n_pts_   = -1;
	dim_     = -1;
	B_       = -1;
	l_       = -1;
	m_       = -1;
	page_io_ = -1;
	dist_io_ = -1;
	proj_    = NULL;
	table_   = NULL;
	trees_   = NULL;
}

// -----------------------------------------------------------------------------
QDAFN::~QDAFN()						// destructor
{
	for (int i = 0; i < l_; ++i) { delete[] proj_[i]; proj_[i] = NULL; }
	delete[] proj_; proj_ = NULL;
	g_memory -= SIZEFLOAT * l_ * dim_;

	if (trees_ != NULL) {
		for (int i = 0; i < l_; ++i) {
			if (trees_[i] != NULL) { delete trees_[i]; trees_[i] = NULL; }
		}
		delete[] trees_; trees_ = NULL;
	}
	if (table_ != NULL) {
		for (int i = 0; i < l_; ++i) {
			if (table_[i] != NULL) { delete table_[i]; table_[i] = NULL; }
		}
		delete[] table_; table_ = NULL;
		g_memory -= sizeof(Result) * l_ * m_;
	}
}

// -----------------------------------------------------------------------------
int QDAFN::build(					// build index
	int   n,							// number of data objects
	int   d,							// dimension of space
	int   B,							// page size
	int   l,							// number of projections
	int   m,							// number of candidates
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char  *path)					// index path
{
	// -------------------------------------------------------------------------
	//  set up parameters of QDAFN
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	B_     = B;

	strcpy(path_, path); create_dir(path_);

	if (l == 0 || m == 0) {
		l_ = 2 * (int) ceil(pow((float) n, 1.0F/(ratio*ratio)));
		if (l_ < 1) {
			printf("bad number of projection <l> %d\n", l_);
			return 1;
		}

		float x = pow(log((float) n), (ratio*ratio/2.0F - 1.0F/3.0F));
		m_ = 1 + (int) ceil(E * E * l_ * x);
		if (m_ < 1) {
			printf("bad number of candidates <m> %d\n", m_);
			return 1;
		}
	}
	else {
		l_ = l;
		m_ = m;
	}

	// -------------------------------------------------------------------------
	//  generate random projection directions
	// -------------------------------------------------------------------------
	zigset(MAGIC + 17); 			// use fix seed 
	// zigset(MAGIC + time(NULL));

	g_memory += SIZEFLOAT * l_ * dim_;
	proj_ = new float*[l_];
	for (int i = 0; i < l_; ++i) {
		proj_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			proj_[i][j] = RNOR / sqrt((float) dim_);
		}
	}

	// -------------------------------------------------------------------------
	//  write the "para" file
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "wb");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	fwrite(&n_pts_, SIZEINT, 1, fp);
	fwrite(&dim_,   SIZEINT, 1, fp);
	fwrite(&B_,     SIZEINT, 1, fp);
	fwrite(&l_,     SIZEINT, 1, fp);
	fwrite(&m_,     SIZEINT, 1, fp);
	for (int i = 0; i < l_; ++i) fwrite(proj_[i], SIZEFLOAT, dim_, fp);

	// -------------------------------------------------------------------------
	//  build index (QAB+tree or simply array)
	// -------------------------------------------------------------------------
	if (m_ > CANDIDATES) {
		Result *table = new Result[n_pts_];
		trees_ = new B_Tree*[l_];
		for (int i = 0; i < l_; ++i) {
			for (int j = 0; j < n_pts_; ++j) {
				table[j].id_  = j;
				table[j].key_ = calc_proj(i, data[j]);
			}
			qsort(table, n_pts_, sizeof(Result), ResultComp);

			// -----------------------------------------------------------------
			//  build index with QAB+trees
			// -----------------------------------------------------------------
			get_tree_filename(i, fname);

			trees_[i] = new B_Tree();
			trees_[i]->init(B_, fname);
			if (trees_[i]->bulkload(m_, table)) return 1;
		}
		delete[] table; table = NULL;
	}
	else {
		table_ = new Result*[l_];
		for (int i = 0; i < l_; ++i) {
			table_[i] = new Result[n_pts_];
			for (int j = 0; j < n_pts_; ++j) {
				table_[i][j].id_  = j;
				table_[i][j].key_ = calc_proj(i, data[j]);
			}
			qsort(table_[i], n_pts_, sizeof(Result), ResultComp);

			// -----------------------------------------------------------------
			//  store in the 'para' file
			// -----------------------------------------------------------------
			fwrite(table_[i], sizeof(Result), m_, fp);
		}
	}
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
float QDAFN::calc_proj(				// calc projection of input data object
	int   id,							// projection vector id
	const float *data)					// input data object
{
	return calc_inner_product(dim_, (const float*) proj_[id], data);
}

// -----------------------------------------------------------------------------
inline void QDAFN::get_tree_filename( // get file name of b-tree
	int  tid,							// tree id, from 0 to m-1
	char *fname)						// file name (return)
{
	sprintf(fname, "%s%d.qdafn", path_, tid);
}

// -----------------------------------------------------------------------------
void QDAFN::display()				// display parameters
{
	printf("Parameters of QDAFN (SISAP2015 paper):\n");
	printf("    n    = %d\n", n_pts_);
	printf("    d    = %d\n", dim_);
	printf("    B    = %d\n", B_);
	printf("    l    = %d\n", l_);
	printf("    m    = %d\n", m_);
	printf("    algo = %s\n", algoname[2]);
	printf("    path = %s\n", path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int QDAFN::load(					// load index
	const char *path)					// index path
{
	strcpy(path_, path);

	// -------------------------------------------------------------------------
	//  read the "para" file
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "rb");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	fread(&n_pts_, SIZEINT, 1, fp);
	fread(&dim_,   SIZEINT, 1, fp);
	fread(&B_,     SIZEINT, 1, fp);
	fread(&l_,     SIZEINT, 1, fp);
	fread(&m_,     SIZEINT, 1, fp);

	g_memory += SIZEFLOAT * l_ * dim_;
	proj_ = new float*[l_];
	for (int i = 0; i < l_; ++i) {
		proj_[i] = new float[dim_];
		fread(proj_[i], SIZEFLOAT, dim_, fp);
	}

	if (m_ > CANDIDATES) {
		// ---------------------------------------------------------------------
		//  load <l> B+ trees
		// ---------------------------------------------------------------------
		trees_ = new B_Tree*[l_];
		for (int i = 0; i < l_; ++i) {
			get_tree_filename(i, fname);

			trees_[i] = new B_Tree();
			trees_[i]->init_restore(fname);
		}
	}
	else {
		// ---------------------------------------------------------------------
		//  read from 'para' file
		// ---------------------------------------------------------------------
		g_memory += sizeof(Result) * l_ * m_;
		table_ = new Result*[l_];
		for (int i = 0; i < l_; ++i) {
			table_[i] = new Result[m_];
			fread(table_[i], sizeof(Result), m_, fp);
		}
	}
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
uint64_t QDAFN::search(				// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// query object
	const char *data_folder,			// new format data folder
	MaxK_List *list)					// top-k results (return)
{
	if (m_ > CANDIDATES) return ext_search(top_k, query, data_folder, list);
	else return int_search(top_k, query, data_folder, list);
}

// -----------------------------------------------------------------------------
uint64_t QDAFN::int_search(			// internal search
	int   top_k,						// top-k value
	const float *query,					// query object
	const char *data_folder,			// new format data folder
	MaxK_List *list)					// top-k results (return)
{
	// -------------------------------------------------------------------------
	//  allocation and initialize <proj_q>
	// -------------------------------------------------------------------------
	std::vector<int>  next(l_, 0);
	std::vector<bool> checked(n_pts_, false);

	float *proj_q = new float[l_];
	for (int i = 0; i < l_; ++i) proj_q[i] = calc_proj(i, query);;

	float *data = new float[dim_];
	for (int i = 0; i < dim_; ++i) data[i] = -1.0f;

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int cand = std::min(m_+top_k, n_pts_);
	
	std::priority_queue<Result, std::vector<Result>, Cmp> pri_queue;
	Result q_item;
	for (int i = 0; i < l_; ++i) {
		q_item.key_ = fabs(table_[i][next[i]].key_ - proj_q[i]);
		q_item.id_  = i;

		pri_queue.push(q_item);
	}

	dist_io_ = 0;
	for (int i = 0; i < cand; ++i) {
		// ---------------------------------------------------------------------
		//  get obj with largest proj dist and remove it from the queue
		// ---------------------------------------------------------------------
		if (pri_queue.empty()) break;
		q_item = pri_queue.top();
		pri_queue.pop();

		// ---------------------------------------------------------------------
		//  check candidate
		// ---------------------------------------------------------------------
		int pid = q_item.id_;
		int id  = table_[pid][next[pid]].id_;
		if (!checked[id]) {
			checked[id] = true;
			read_data_new_format(id, dim_, B_, data_folder, data);

			float dist = calc_l2_dist(dim_, (const float *) data, query);
			list->insert(dist, id + 1);
			++dist_io_;
		}
		// ---------------------------------------------------------------------
		//  update priority queue
		// ---------------------------------------------------------------------
		if (++next[pid] < m_) {
			q_item.key_ = fabs(table_[pid][next[pid]].key_ - proj_q[pid]);
			pri_queue.push(q_item);
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	while (!pri_queue.empty()) pri_queue.pop();

	delete[] proj_q; proj_q = NULL;
	delete[] data;   data   = NULL;

	return dist_io_;
}

// -----------------------------------------------------------------------------
uint64_t QDAFN::ext_search(			// external search
	int   top_k,						// top-k value
	const float *query,					// query object
	const char *data_folder,			// new format data folder
	MaxK_List *list)					// top-k results (return)
{
	// -------------------------------------------------------------------------
	//  allocation and initialization
	// -------------------------------------------------------------------------
	std::vector<bool> checked(n_pts_, false);

	float *data = new float[dim_];
	for (int i = 0; i < dim_; ++i) data[i] = 0.0f;

	float *proj_q = new float[l_];
	QDAFN_Page *page = new QDAFN_Page[l_];
	for (int i = 0; i < l_; ++i) {
		page[i].node_ = NULL;
		page[i].pos_  = -1;
	}

	// -------------------------------------------------------------------------
	//  compute hash value <proj_q> of query and init page buffers <page> 
	// -------------------------------------------------------------------------
	page_io_ = 0;					// page i/os for search
	dist_io_ = 0;					// i/os for distance computation
	init_buffer(query, page, proj_q);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int cand = std::min(m_ + top_k, n_pts_);

	std::priority_queue<Result, std::vector<Result>, Cmp> pri_queue;
	Result q_item;
	for (int i = 0; i < l_; ++i) {
		if (page[i].node_) {
			q_item.key_ = calc_dist(proj_q[i], &page[i]);
			q_item.id_  = i;

			pri_queue.push(q_item);
		}
	}

	for (int i = 0; i < cand; ++i) {
		// ---------------------------------------------------------------------
		//  get obj with largest proj dist and remove it from the queue
		// ---------------------------------------------------------------------
		if (pri_queue.empty()) break;
		q_item = pri_queue.top();	// get the object with largest proj dist
		pri_queue.pop();			// delete the object from the queue

		// ---------------------------------------------------------------------
		//  check candidate
		// ---------------------------------------------------------------------
		int j = q_item.id_;
		int id = page[j].node_->get_son(page[j].pos_);
		if (!checked[id]) {
			checked[id] = true;
			read_data_new_format(id, dim_, B_, data_folder, data);

			float dist = calc_l2_dist(dim_, (const float *) data, query);
			list->insert(dist, id + 1);
			++dist_io_;
		}

		// ---------------------------------------------------------------------
		//  update priority queue
		// ---------------------------------------------------------------------
		update_page(&page[j]);
		if (page[j].node_) {
			q_item.key_ = calc_dist(proj_q[j], &page[j]);
			pri_queue.push(q_item);
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	while (!pri_queue.empty()) pri_queue.pop();

	delete[] proj_q; proj_q = NULL;
	delete[] data;   data   = NULL;

	for (int i = 0; i < l_; ++i) {
		if (page[i].node_ != NULL) { 
			delete page[i].node_; page[i].node_ = NULL; 
		}
	}
	delete[] page; page = NULL;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
void QDAFN::init_buffer(			// init page buffer
	const float *query,					// query point
	QDAFN_Page *page,					// buffer page (return)
	float *proj_q)						// projection of query (return)
{
	int block = -1;
	B_Node *node = NULL;

	for (int i = 0; i < l_; ++i) {
		proj_q[i] = calc_proj(i, query);

		block = trees_[i]->root_;
		node = new B_Node();
		node->init_restore(trees_[i], block);
		++page_io_;

		if (node->get_level() == 0) {
			// -----------------------------------------------------------------
			//  leaf level
			// -----------------------------------------------------------------
			page[i].node_ = node;
			page[i].pos_ = 0;

			node = NULL;
		}
		else {
			// -----------------------------------------------------------------
			//  non-leaf level
			// -----------------------------------------------------------------
			while (node->get_level() > 1) {
				block = node->get_son(0);
				delete node; node = NULL;

				node = new B_Node();
				node->init_restore(trees_[i], block);
				++page_io_;
			}

			block = node->get_son(0);
			page[i].node_ = new B_Node();
			page[i].node_->init_restore(trees_[i], block);
			page[i].pos_ = 0;
			++page_io_;

			if (node != NULL) {
				delete node; node = NULL;
			}
		}
	}
}

// -----------------------------------------------------------------------------
void QDAFN::update_page(			// update right node info
	QDAFN_Page *page)					// page buffer
{
	B_Node *node = NULL;
	B_Node *old_node = NULL;

	++page->pos_;
	if (page->pos_ > page->node_->get_num_entries() - 1) {
		old_node = page->node_;
		node = page->node_->get_right_sibling();

		if (node != NULL) {
			page->node_ = node;
			page->pos_ = 0;
			++page_io_;
		}
		else {
			page->node_ = NULL;
			page->pos_ = -1;
		}
		delete old_node; old_node = NULL;
	}
}

// -----------------------------------------------------------------------------
inline float QDAFN::calc_dist(		// calc proj_dist
	float proj_q,						// projection of query
	const QDAFN_Page *page)				// page buffer
{
	float key = page->node_->get_key(page->pos_);
	return fabs(key - proj_q);
}
