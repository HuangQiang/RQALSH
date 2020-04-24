#include "rqalsh.h"

// -----------------------------------------------------------------------------
RQALSH::RQALSH()					// constructor
{
	n_pts_ = -1;
	dim_   = -1;
	B_     = -1;
	beta_  = -1.0f;
	delta_ = -1.0f;
	ratio_ = -1.0f;
	w_     = -1.0f;
	m_     = -1;
	l_     = -1;
	a_     = NULL;
	trees_ = NULL;

	dist_io_ = -1;
	page_io_ = -1;
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i]; a_[i] = NULL;
		delete trees_[i]; trees_[i] = NULL;
	}
	delete[] a_; a_ = NULL;
	delete[] trees_; trees_ = NULL;

	g_memory -= SIZEFLOAT * m_ * dim_;
}

// -----------------------------------------------------------------------------
int RQALSH::build(					// build index
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char *path)					// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	B_     = B;
	beta_  = (float) beta / (float) n;
	delta_ = delta;
	ratio_ = ratio;
	strcpy(path_, path); create_dir(path_);

	// -------------------------------------------------------------------------
	//  init <w_> <m_> and <l_> (auto tuning-w)
	// -------------------------------------------------------------------------
	w_ = sqrt((8.0f * log(ratio_)) / (ratio_ * ratio_ - 1.0f));
	
	float p1 = calc_l2_prob(w_ / 2.0f);
	float p2 = calc_l2_prob(w_ * ratio_ / 2.0f);

	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);
	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			a_[i][j] = gaussian(0.0f, 1.0f);
		}
	}

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1;

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	if (bulkload(data)) return 1;

	return 0;
}

// -----------------------------------------------------------------------------
inline float RQALSH::calc_l2_prob(	// calc prob <p1_> and <p2_> of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
int RQALSH::write_params()			// write parameters to disk
{
	char fname[200];
	strcpy(fname, path_); strcat(fname, "para");

	FILE *fp = fopen(fname, "rb");
	if (fp)	{ printf("Hash Tables Already Exist\n\n"); return 1; }
	fp = fopen(fname, "wb");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	fwrite(&n_pts_, SIZEINT,   1, fp);
	fwrite(&dim_,   SIZEINT,   1, fp);
	fwrite(&B_,     SIZEINT,   1, fp);
	fwrite(&beta_,  SIZEFLOAT, 1, fp);
	fwrite(&delta_, SIZEFLOAT, 1, fp);
	fwrite(&ratio_, SIZEFLOAT, 1, fp);
	fwrite(&w_,     SIZEFLOAT, 1, fp);
	fwrite(&m_,     SIZEINT,   1, fp);
	fwrite(&l_,     SIZEINT,   1, fp);

	for (int i = 0; i < m_; ++i) fwrite(a_[i], SIZEFLOAT, dim_, fp);
	fclose(fp);	

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::bulkload(				// build QAB+Trees by bulkloading
	const float** data)					// data set
{
	// -------------------------------------------------------------------------
	//  write hash tables (indexed by B+ Tree) to disk
	// -------------------------------------------------------------------------
	Result *table = new Result[n_pts_];

	trees_ = new QAB_Tree*[m_];
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < n_pts_; ++j) {
			table[j].id_  = j;
			table[j].key_ = calc_hash_value(i, data[j]);
		}
		qsort(table, n_pts_, sizeof(Result), ResultComp);

		char fname[200];
		get_tree_filename(i, fname);
		trees_[i] = new QAB_Tree();
		trees_[i]->init(B_, fname);
		if (trees_[i]->bulkload(n_pts_, (const Result *) table)) return 1;
	}
	delete[] table; table = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
float RQALSH::calc_hash_value( 		// calc hash value
	int   tid,							// hash table id
	const float *data)					// one data object
{
	return calc_inner_product(dim_, a_[tid], data);
}

// -----------------------------------------------------------------------------
inline void RQALSH::get_tree_filename( // get file name of QAB+Tree
	int  tid,							// tree id, from 0 to m-1
	char *fname)						// file name (return)
{
	sprintf(fname, "%s%d.rqalsh", path_, tid);
}

// -----------------------------------------------------------------------------
void RQALSH::display()				// display parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    B     = %d\n",   B_);
	printf("    beta  = %f\n",   beta_);
	printf("    delta = %f\n",   delta_);
	printf("    ratio = %.1f\n", ratio_);
	printf("    w     = %f\n",   w_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n",   l_);
	printf("    path  = %s\n",   path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH::load(					// load index
	const char *path)					// index path
{
	// -------------------------------------------------------------------------
	//  read parameters from disk
	// -------------------------------------------------------------------------
	strcpy(path_, path);
	if (read_params()) return 1;

	// -------------------------------------------------------------------------
	//  load qab-tree for k-FN search
	// -------------------------------------------------------------------------
	char fname[200];
	trees_ = new QAB_Tree*[m_];
	for (int i = 0; i < m_; ++i) {
		get_tree_filename(i, fname);

		trees_[i] = new QAB_Tree();
		trees_[i]->init_restore(fname);
	}
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::read_params()			// read parameters from disk
{
	char fname[200];
	strcpy(fname, path_); strcat(fname, "para");

	FILE *fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fread(&n_pts_, SIZEINT,   1, fp);
	fread(&dim_,   SIZEINT,   1, fp);
	fread(&B_,     SIZEINT,   1, fp);
	fread(&beta_,  SIZEFLOAT, 1, fp);
	fread(&delta_, SIZEFLOAT, 1, fp);
	fread(&ratio_, SIZEFLOAT, 1, fp);
	fread(&w_,     SIZEFLOAT, 1, fp);
	fread(&m_,     SIZEINT,   1, fp);
	fread(&l_,     SIZEINT,   1, fp);
	
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		fread(a_[i], SIZEFLOAT, dim_, fp);
	}
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
uint64_t RQALSH::kfn(				// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// query object
	const int *index,					// mapping index for data objects
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-FN results (return)
{
	int   *freq    = new int[n_pts_];
	bool  *checked = new bool[n_pts_];
	bool  *flag    = new bool[m_];
	float *q_val   = new float[m_];
	float *data    = new float[dim_];
	
	Page **lptrs = new Page*[m_];
	Page **rptrs = new Page*[m_];
	for (int i = 0; i < m_; ++i) {
		lptrs[i] = new Page();
		rptrs[i] = new Page();
	}

	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	memset(freq, 0, n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);

	init_search_params(query, q_val, lptrs, rptrs);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   candidates = CANDIDATES + top_k - 1; // threshold of candidates
	float kdist  = MINREAL;			// k-th furthest neighbor distance
	float radius = find_radius(q_val, (const Page**) lptrs, (const Page**) rptrs);
	float width  = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_flag = 0;
		memset(flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects (dynamic separation counting)
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!flag[i]) continue;

				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				Page *lptr = lptrs[i];
				Page *rptr = rptrs[i];

				float ldist = -1.0f;
				float rdist = -1.0f;
				if (lptr->size_ != -1) ldist = calc_dist(q_val[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val[i], rptr);

				// -------------------------------------------------------------
				//  step 2.2: determine the closer direction (left or right)
				//  and do separation counting to find frequent objects.
				//
				//  For the frequent object, we calc the Lp distance with
				//  query, and update the c-k-AFN results.
				// -------------------------------------------------------------
				if (ldist > width && ldist > rdist) {
					int count = lptr->size_;
					int start = lptr->leaf_pos_;
					int end   = start + count;
					
					for (int j = start; j < end; ++j) {
						int id = lptr->leaf_node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							if (index != NULL) id = index[id];
							read_data_new_format(id, dim_, B_, data_folder, data);

							float dist = calc_l2_dist(dim_, data, query);
							kdist = list->insert(dist, id + 1);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(rptr, lptr);
				}
				else if (rdist > width && ldist <= rdist) {
					int count = rptr->size_;
					int end   = rptr->leaf_pos_;
					int start = end - count;

					for (int j = end; j > start; --j) {
						int id = rptr->leaf_node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							if (index != NULL) id = index[id];
							read_data_new_format(id, dim_, B_, data_folder, data);

							float dist = calc_l2_dist(dim_, data, query);
							kdist = list->insert(dist, id + 1);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					flag[i] = false;
					++num_flag;
				}
				if (num_flag >= m_ || dist_io_ >= candidates) break;
			}
			if (num_flag >= m_ || dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist > radius / ratio_ && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = radius / ratio_;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete_tree_ptr(lptrs, rptrs);	

	delete[] freq;    freq    = NULL;
	delete[] checked; checked = NULL;
	delete[] flag;    flag    = NULL;
	delete[] q_val;   q_val   = NULL;
	delete[] data;    data    = NULL;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
void RQALSH::init_search_params(	// init parameters
	const float *query,					// query object
	float *q_val,						// hash values of query (return)
	Page  **lptrs,						// left buffer (return)
	Page  **rptrs)						// right buffer (return)
{
	page_io_ = 0;
	dist_io_ = 0;

	for (int i = 0; i < m_; ++i) {
		lptrs[i]->leaf_node_ = NULL;
		lptrs[i]->index_pos_ = -1;
		lptrs[i]->leaf_pos_  = -1;
		lptrs[i]->size_      = -1;

		rptrs[i]->leaf_node_ = NULL;
		rptrs[i]->index_pos_ = -1;
		rptrs[i]->leaf_pos_  = -1;
		rptrs[i]->size_      = -1;
	}

	QAB_IndexNode *index_node = NULL;
	QAB_LeafNode  *leaf_node  = NULL;

	int block       = -1;			// variables for index node
	int pos         = -1;			// variables for leaf node
	int increment   = -1;
	int num_entries = -1;
	int num_keys    = -1;

	for (int i = 0; i < m_; ++i) {
		float q_v = calc_hash_value(i, query);
		QAB_Tree   *tree = trees_[i];
		Page *lptr = lptrs[i];
		Page *rptr = rptrs[i];

		q_val[i] = q_v;
		block = tree->root_;
		if (block == 1) {
			// -----------------------------------------------------------------
			//  the qab+tree has only one leaf node
			// -----------------------------------------------------------------
			lptr->leaf_node_ = new QAB_LeafNode();
			lptr->leaf_node_->init_restore(trees_[i], block);
			++page_io_;

			leaf_node   = lptr->leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();
			if (num_keys > 1) {
				lptr->index_pos_ = 0;
				lptr->leaf_pos_  = 0;
				lptr->size_      = increment;

				rptr->leaf_node_ = lptr->leaf_node_;
				rptr->index_pos_ = num_keys - 1;
				rptr->leaf_pos_  = num_entries - 1;
				rptr->size_      = num_entries - (num_keys - 1) * increment;
			}
			else {
				lptr->index_pos_ = 0;
				lptr->leaf_pos_  = 0;
				lptr->size_      = num_entries;

				rptr->leaf_node_ = NULL;
				rptr->index_pos_ = -1;
				rptr->leaf_pos_  = -1;
				rptr->size_      = -1;
			}
		}
		else {
			// -----------------------------------------------------------------
			//  the qab+tree has index node
			// 
			//  (1) initialize left leaf node
			// -----------------------------------------------------------------
			block = tree->root_;
			index_node = new QAB_IndexNode();
			index_node->init_restore(tree, block);
			++page_io_;
									// find the left most leaf node
			while (index_node->get_level() > 1) {
				block = index_node->get_son(0);
				delete index_node; index_node = NULL;

				index_node = new QAB_IndexNode();
				index_node->init_restore(tree, block);
				++page_io_;			// access a new node (a new page)
			}

			block = index_node->get_son(0);
			lptr->leaf_node_ = new QAB_LeafNode();
			lptr->leaf_node_->init_restore(tree, block);
			++page_io_;

			lptr->index_pos_ = 0;
			lptr->leaf_pos_  = 0;

			increment   = lptr->leaf_node_->get_increment();
			num_entries = lptr->leaf_node_->get_num_entries();
			if (increment > num_entries) {
				lptr->size_ = num_entries;
			} else {
				lptr->size_ = increment;
			}

			if (index_node != NULL) {
				delete index_node; index_node = NULL;
			}

			// -----------------------------------------------------------------
			//  Initialize right leaf node
			// -----------------------------------------------------------------
			block = tree->root_;
			index_node = new QAB_IndexNode();
			index_node->init_restore(tree, block);
			++page_io_;
									// find the right most leaf node
			while (index_node->get_level() > 1) {
				num_entries = index_node->get_num_entries();
				block = index_node->get_son(num_entries - 1);
				delete index_node; index_node = NULL;

				index_node = new QAB_IndexNode();
				index_node->init_restore(tree, block);
				++page_io_;			// access a new node (a new page)
			}

			num_entries = index_node->get_num_entries();
			block = index_node->get_son(num_entries - 1);
			rptr->leaf_node_ = new QAB_LeafNode();
			rptr->leaf_node_->init_restore(tree, block);
			++page_io_;

			leaf_node   = rptr->leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();

			rptr->index_pos_ = num_keys - 1;
			rptr->leaf_pos_  = num_entries - 1;
			rptr->size_      = num_entries - (num_keys - 1) * increment;

			if (index_node != NULL) {
				delete index_node; index_node = NULL;
			}
		}
	}
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	const float *q_val,					// hash value of query
	const Page **lptrs,					// left buffer
	const Page **rptrs)					// right buffer
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	std::vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lptrs[i]->size_ != -1) {
			list.push_back(calc_dist(q_val[i], lptrs[i]));
		}
		if (rptrs[i]->size_ != -1) {
			list.push_back(calc_dist(q_val[i], rptrs[i]));
		}
	}

	// -------------------------------------------------------------------------
	//  sort the array in ascending order 
	// -------------------------------------------------------------------------
	std::sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int   num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(ratio_));
	return pow(ratio_, kappa);
}

// -----------------------------------------------------------------------------
void RQALSH::update_left_buffer(	// update left buffer
	const Page *rptr,					// right buffer
	Page *lptr)							// left buffer (return)
{
	QAB_LeafNode* leaf_node     = NULL;
	QAB_LeafNode* old_leaf_node = NULL;

	if (lptr->index_pos_ < lptr->leaf_node_->get_num_keys() - 1) {
		lptr->index_pos_++;

		int pos         = lptr->index_pos_;
		int increment   = lptr->leaf_node_->get_increment();
		lptr->leaf_pos_ = pos * increment;
		if (pos == lptr->leaf_node_->get_num_keys() - 1) {
			int num_entries = lptr->leaf_node_->get_num_entries();
			lptr->size_ = num_entries - pos * increment;
		} else {
			lptr->size_ = increment;
		}
	}
	else {
		old_leaf_node = lptr->leaf_node_;
		leaf_node     = lptr->leaf_node_->get_right_sibling();

		if (leaf_node) {
			lptr->leaf_node_ = leaf_node;
			lptr->index_pos_ = 0;
			lptr->leaf_pos_  = 0;

			int increment    = leaf_node->get_increment();
			int num_entries  = leaf_node->get_num_entries();
			if (increment > num_entries) {
				lptr->size_ = num_entries;
			} else {
				lptr->size_ = increment;
			}
			++page_io_;
		}
		else {
			lptr->leaf_node_ = NULL;
			lptr->index_pos_ = -1;
			lptr->leaf_pos_  = -1;
			lptr->size_      = -1;
		}

		if (rptr->leaf_node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
void RQALSH::update_right_buffer(	// update right buffer
	const Page* lptr,					// left buffer
	Page* rptr)							// right buffer (return)
{
	QAB_LeafNode* leaf_node     = NULL;
	QAB_LeafNode* old_leaf_node = NULL;

	if (rptr->index_pos_ > 0) {
		rptr->index_pos_--;

		int pos         = rptr->index_pos_;
		int increment   = rptr->leaf_node_->get_increment();
		rptr->leaf_pos_ = pos * increment + increment - 1;
		rptr->size_     = increment;
	}
	else {
		old_leaf_node = rptr->leaf_node_;
		leaf_node     = rptr->leaf_node_->get_left_sibling();

		if (leaf_node) {
			rptr->leaf_node_ = leaf_node;
			rptr->index_pos_ = leaf_node->get_num_keys() - 1;

			int pos          = rptr->index_pos_;
			int increment    = leaf_node->get_increment();
			int num_entries  = leaf_node->get_num_entries();
			rptr->leaf_pos_  = num_entries - 1;
			rptr->size_      = num_entries - pos * increment;
			++page_io_;
		}
		else {
			rptr->leaf_node_ = NULL;
			rptr->index_pos_ = -1;
			rptr->leaf_pos_  = -1;
			rptr->size_      = -1;
		}

		if (lptr->leaf_node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
inline float RQALSH::calc_dist(		// calc projected distance
	float q_val,						// hash value of query
	const Page *ptr)					// page buffer
{
	int   pos  = ptr->index_pos_;
	float key  = ptr->leaf_node_->get_key(pos);

	return fabs(key - q_val);
}

// -----------------------------------------------------------------------------
void RQALSH::delete_tree_ptr(		// delete the pointers of QAB+Trees
	Page **lptrs,						// left buffer (return)
	Page **rptrs)						// right buffer (return)
{
	for (int i = 0; i < m_; ++i) {
		// ---------------------------------------------------------------------
		//  CANNOT remove the condition
		//              <lptrs[i].leaf_node != rptrs[i].leaf_node>
		//  because <lptrs[i].leaf_node> and <rptrs[i].leaf_node> may point 
		//  to the same address, then we would delete it twice and receive 
		//  the runtime error or segmentation fault.
		// ---------------------------------------------------------------------
		if (lptrs[i]->leaf_node_ && lptrs[i]->leaf_node_!=rptrs[i]->leaf_node_) {
			delete lptrs[i]->leaf_node_; lptrs[i]->leaf_node_ = NULL;
		}
		if (rptrs[i]->leaf_node_) {
			delete rptrs[i]->leaf_node_; rptrs[i]->leaf_node_ = NULL;
		}

		if (lptrs[i] != NULL) { delete[] lptrs[i]; lptrs[i] = NULL; }
		if (rptrs[i] != NULL) { delete[] rptrs[i]; rptrs[i] = NULL; }
	}
	delete[] lptrs; lptrs = NULL;
	delete[] rptrs; rptrs = NULL;
}
