#include "headers.h"


// -----------------------------------------------------------------------------
//  RQALSH: structure of rqalsh indexed by b+ tree. RQALSH is used to solve
//  the problem of c-Approximate Furthest Neighbor (c-AFN) search.
// -----------------------------------------------------------------------------
RQALSH::RQALSH()					// constructor
{
	n_pts_ = dim_ = B_ = m_ = l_ = -1;

	appr_ratio_ = w_ = p1_ = p2_ = -1.0f;
	alpha_ = beta_ = delta_ = -1.0f;

	a_array_ = NULL;
	trees_ = NULL;
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	if (a_array_) {
		delete[] a_array_; a_array_ = NULL;
		g_memory -= SIZEFLOAT * m_ * dim_;
	}

	if (trees_) {
		for (int i = 0; i < m_; i++) {
			delete trees_[i]; trees_[i] = NULL;
		}
		delete[] trees_; trees_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void RQALSH::init(					// init params of rqalsh
	int   n,							// number of points
	int   d,							// dimension of space
	int   B,							// page size
	int   beta,							// false probability percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	char* output_folder)				// folder of info of rqalsh
{
	n_pts_ = n;						// init <n_pts_>
	dim_   = d;						// init <dim_>
	B_     = B;						// init <B_>
	appr_ratio_ = ratio;			// init <appr_ratio_>

	beta_ = (float) beta / n_pts_;
	delta_ = delta;
									// init <index_path_>
	strcpy(index_path_, output_folder);
	strcat(index_path_, "indices/");

									// init <w_> <delta_> <beta_> <p1_>
	calc_params();					// <p2_> <alpha_> <m_> and <l_>

	gen_hash_func();				// init <a_array_>
	display_params();				// display params
	trees_ = NULL;					// init <trees_>
}

// -----------------------------------------------------------------------------
void RQALSH::calc_params()			// calc params of rqalsh
{
	// -------------------------------------------------------------------------
	//  init <delta_> and <beta_>
	// -------------------------------------------------------------------------
	//beta_  = 100.0f / n_pts_;
	//delta_ = 0.49F;

	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	// -------------------------------------------------------------------------
	//w_  = sqrt((8.0f*log(appr_ratio_)) / (appr_ratio_*appr_ratio_-1.0f));
	w_ = 1.0f;
	
	p1_ = calc_l2_prob(w_ / 2.0f);
	p2_ = calc_l2_prob(w_ * appr_ratio_ / 2.0f);

	// -------------------------------------------------------------------------
	//  init <alpha_> <m_> and <l_>
	// -------------------------------------------------------------------------
	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);

	float eta = para1 / para2;
	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha_ * m_);
	//l_ = (int) ceil((p1_ * para1 + p2_ * para2) * (para1 + para2) / para3);
}

// -----------------------------------------------------------------------------
float RQALSH::calc_l2_prob(			// calc prob <p1_> and <p2_> of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
void RQALSH::display_params()		// display params of rqalsh
{
	printf("Parameters of RQALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    B          = %d\n", B_);
	printf("    ratio      = %0.2f\n", appr_ratio_);
	printf("    w          = %0.4f\n", w_);
	printf("    p1         = %0.4f\n", p1_);
	printf("    p2         = %0.4f\n", p2_);
	printf("    alpha      = %0.6f\n", alpha_);
	printf("    beta       = %0.6f\n", beta_);
	printf("    delta      = %0.6f\n", delta_);
	printf("    m          = %d\n", m_);
	printf("    l          = %d\n", l_);
	printf("    beta * n   = %d\n", 100);
	printf("    index path = %s\n\n", index_path_);
}

// -----------------------------------------------------------------------------
void RQALSH::gen_hash_func()		// generate hash function <a_array>
{
	g_memory += SIZEFLOAT * m_ * dim_;
	a_array_ = new float[m_ * dim_];

	for (int i = 0; i < m_ * dim_; i++) {
		a_array_[i] = gaussian(0.0f, 1.0f);
	}
}

// -----------------------------------------------------------------------------
int RQALSH::bulkload(				// build m b-trees by bulkloading
	float** data)						// data set
{
	// -------------------------------------------------------------------------
	//  Check whether the default maximum memory is enough
	// -------------------------------------------------------------------------
	g_memory += sizeof(HashValue) * n_pts_;
	if (check_mem()) {
		printf("*** memory = %.2f MB\n\n", g_memory / (1024.0f * 1024.0f));
		return 1;
	}

	// -------------------------------------------------------------------------
	//  Check whether the directory exists. If the directory does not exist, we
	//  create the directory for each folder.
	// -------------------------------------------------------------------------
	create_dir(index_path_);

	// -------------------------------------------------------------------------
	//  Write the file "para" where the parameters and hash functions are 
	//  stored in it.
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, index_path_);		// write the "para" file
	strcat(fname, "para");
	if (write_para_file(fname)) return 1;

	// -------------------------------------------------------------------------
	//  Write the hash tables (indexed by b+ tree) to the disk
	// -------------------------------------------------------------------------
	HashValue* hashtable = new HashValue[n_pts_];
	for (int i = 0; i < m_; i++) {
		printf("    Tree %3d (out of %d)\n", i + 1, m_);

		printf("        Computing Hash Values...\n");
		for (int j = 0; j < n_pts_; j++) {
			hashtable[j].id_ = j;
			hashtable[j].proj_ = calc_hash_value(i, data[j]);
		}

		printf("        Sorting...\n");
		qsort(hashtable, n_pts_, sizeof(HashValue), HashValueQsortComp);

		printf("        Bulkloading...\n");
		get_tree_filename(i, fname);

		BTree *bt = new BTree();
		bt->init(fname, B_);
		if (bt->bulkload(hashtable, n_pts_)) {
			return 1;
		}
		delete bt; bt = NULL;
	}

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (hashtable != NULL) {
		delete[] hashtable; hashtable = NULL;
		g_memory -= sizeof(HashValue) * n_pts_;
	}
	return 0;						// success to return
}

// -----------------------------------------------------------------------------
int RQALSH::write_para_file(		// write "para" file from disk
	char* fname)						// file name of "para" file
{
	FILE* fp = NULL;
	fp = fopen(fname, "r");
	if (fp) {						// ensure the file not exist
		error("Hash tables exist.\n", true);
	}

	fp = fopen(fname, "w");			// open "para" file to write
	if (!fp) {
		printf("Could not create %s.\n", fname);
		printf("Perhaps no such folder %s?\n", index_path_);
		return 1;					// fail to return
	}

	fprintf(fp, "n = %d\n", n_pts_);
	fprintf(fp, "d = %d\n", dim_);
	fprintf(fp, "B = %d\n", B_);
	fprintf(fp, "ratio = %f\n", appr_ratio_);

	fprintf(fp, "w = %f\n", w_);
	fprintf(fp, "p1 = %f\n", p1_);
	fprintf(fp, "p2 = %f\n", p2_);

	fprintf(fp, "alpha = %f\n", alpha_);
	fprintf(fp, "beta = %f\n", beta_);
	fprintf(fp, "delta = %f\n", delta_);

	fprintf(fp, "m = %d\n", m_);
	fprintf(fp, "l = %d\n", l_);

	int count = 0;
	for (int i = 0; i < m_; i++) {	// write <a_array_>
		fprintf(fp, "%f", a_array_[count++]);
		for (int j = 1; j < dim_; j++) {
			fprintf(fp, " %f", a_array_[count++]);
		}
		fprintf(fp, "\n");
	}
	if (fp) fclose(fp);				// close para file
	
	return 0;						// success to return
}

// -----------------------------------------------------------------------------
float RQALSH::calc_hash_value(		// calc hash value
	int table_id,						// hash table id
	float* point)						// a point
{
	float ret = 0.0f;
	for (int i = 0; i < dim_; i++) {
		ret += (a_array_[table_id * dim_ + i] * point[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
void RQALSH::get_tree_filename(		// get file name of b-tree
	int tree_id,						// tree id, from 0 to m-1
	char* fname)						// file name (return)
{
	char c[20];

	strcpy(fname, index_path_);
	sprintf(c, "%d", tree_id);
	strcat(fname, c);
	strcat(fname, ".rqalsh");
}

// -----------------------------------------------------------------------------
int RQALSH::restore(				// load existing b-trees.
	char* output_folder)				// folder of info of rqalsh
{
	strcpy(index_path_, output_folder);
	strcat(index_path_, "indices/");

	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");
	if (read_para_file(fname)) {
		return 1;
	}

	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; i++) {
		get_tree_filename(i, fname);

		trees_[i] = new BTree();
		trees_[i]->init_restore(fname);
	}
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::read_para_file(			// read "para" file
	char* fname)						// file name of "para" file
{
	FILE* fp = NULL;
	fp = fopen(fname, "r");
	if (!fp) {						// ensure we can open the file
		fprintf(stderr, "Could not open %s.\n", fname);
		exit(1);
	}

	fscanf(fp, "n = %d\n", &n_pts_);
	fscanf(fp, "d = %d\n", &dim_);
	fscanf(fp, "B = %d\n", &B_);
	fscanf(fp, "ratio = %f\n", &appr_ratio_);

	fscanf(fp, "w = %f\n", &w_);
	fscanf(fp, "p1 = %f\n", &p1_);
	fscanf(fp, "p2 = %f\n", &p2_);

	fscanf(fp, "alpha = %f\n", &alpha_);
	fscanf(fp, "beta = %f\n", &beta_);
	fscanf(fp, "delta = %f\n", &delta_);

	fscanf(fp, "m = %d\n", &m_);
	fscanf(fp, "l = %d\n", &l_);

	a_array_ = new float[m_ * dim_];// read <a_array_>
	g_memory += SIZEFLOAT * m_ * dim_;
	int count = 0;
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < dim_; j++) {
			fscanf(fp, "%f", &a_array_[count++]);
		}
		fscanf(fp, "\n");
	}
	if (fp) fclose(fp);				// close para file
	
	display_params();				// display params
	return 0;						// success to return
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// c-k-AFN search
	int top_k,							// top-k value
	float* query,						// one query point
	char*  data_folder,					// folder to store new format of data
	MaxK_List* list)					// c-k-AFN results (returned)
{
	int num_tables = m_;
	int col_threshold = l_;

	// -------------------------------------------------------------------------
	//  Space allocation and initialization
	// -------------------------------------------------------------------------
	g_memory += ((SIZEBOOL + SIZEINT) * n_pts_ + SIZEFLOAT * dim_);
	int* freq = new int[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		freq[i] = 0;				// objects frequency
	}

	bool* checked = new bool[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		checked[i] = false;			// whether an object is checked
	}

	float* data = new float[dim_];
	for (int i = 0; i < dim_; i++) {
		data[i] = 0.0f;				// one data object
	}

	g_memory += (SIZEFLOAT + SIZEBOOL) * num_tables;
	bool* flag = new bool[num_tables];
	for (int i = 0; i < num_tables; i++) {
		flag[i] = true;				// flag to stop search the hash table
	}

	float* q_val = new float[num_tables];
	for (int i = 0; i < num_tables; i++) {
		q_val[i] = -1.0f;			// hash values of the query
	}

	g_memory += (SIZECHAR * B_ * num_tables * 2 + SIZEINT * num_tables * 6);
	PageBuffer* lptr = new PageBuffer[num_tables];
	PageBuffer* rptr = new PageBuffer[num_tables];
	for (int i = 0; i < num_tables; i++) {
		lptr[i].leaf_node_ = NULL;	// left page buffer
		lptr[i].index_pos_ = -1;
		lptr[i].leaf_pos_  = -1;
		lptr[i].size_      = -1;

		rptr[i].leaf_node_ = NULL;	// right page buffer
		rptr[i].index_pos_ = -1;
		rptr[i].leaf_pos_  = -1;
		rptr[i].size_      = -1;
	}

	// -------------------------------------------------------------------------
	//  Compute hash values <q_val> of query and init the page buffers 
	//  <lptr> and <rptr>.
	// -------------------------------------------------------------------------
	page_io_ = 0;					// num of page i/os
	dist_io_ = 0;					// num of dist cmpt
	for (int i = 0; i < num_tables; i++) {
		q_val[i] = calc_hash_value(i, query);
	}
	init_buffer(num_tables, lptr, rptr);

	// -------------------------------------------------------------------------
	//  Determine the basic <radius> and <bucket_width> 
	// -------------------------------------------------------------------------
	float radius = find_radius(num_tables, lptr, rptr, q_val);
	float bucket_width = (w_ * radius / 2.0f);

	// -------------------------------------------------------------------------
	//  c-k-Approximate Furthest Neighbor (c-k-AFN) search
	// -------------------------------------------------------------------------
	bool again      = true;			// stop flag
	int  candidates = 99 + top_k;	// number of candidates
	int  flag_num   = 0;			// used for bucket bound

	int id    = -1;					// current object id
	int count = -1;					// count size in one page
	int start = -1;					// start position
	int end   = -1;					// end position

	float left_dist  = -1.0f;		// left dist with query
	float right_dist = -1.0f;		// right dist with query
	float kfn_dist   = MIN_FLT;		// k-th furthest neighbor distance

	while (again) {
		// ---------------------------------------------------------------------
		//  Step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		flag_num = 0;
		for (int i = 0; i < num_tables; i++) flag[i] = true;

		// ---------------------------------------------------------------------
		//  Step 2: find frequent objects (dynamic separation counting)
		// ---------------------------------------------------------------------
		while (true) {
			for (int i = 0; i < num_tables; i++) {
				if (!flag[i]) continue;

				// -------------------------------------------------------------
				//  Step 2.1: compute <left_dist> and <right_dist>
				// -------------------------------------------------------------
				left_dist = -1.0f;
				if (lptr[i].size_ != -1) {
					left_dist = calc_proj_dist(&lptr[i], q_val[i]);
				}

				right_dist = -1.0f;
				if (rptr[i].size_ != -1) {
					right_dist = calc_proj_dist(&rptr[i], q_val[i]);
				}

				// -------------------------------------------------------------
				//  Step 2.2: determine the closer direction (left or right)
				//  and do separation counting to find frequent objects.
				//
				//  For the frequent object, we calc the Lp distance with
				//  query, and update the c-k-AFN results.
				// -------------------------------------------------------------
				if (left_dist > bucket_width && left_dist > right_dist) {
					count = lptr[i].size_;
					start = lptr[i].leaf_pos_;
					end   = start + count;
					for (int j = start; j < end; j++) {
						id = lptr[i].leaf_node_->get_entry_id(j);
						freq[id]++;
						if (freq[id] > col_threshold && !checked[id]) {
							checked[id] = true;
							read_data(id, dim_, B_, data, data_folder);

							float dist = calc_l2_dist(data, query, dim_);
							int obj_id = id + 1;
							list->insert(dist, obj_id);

							// -------------------------------------------------
							//  Terminating condition 2
							// -------------------------------------------------
							dist_io_++;
							if (dist_io_ >= candidates) {
								again = false;
								flag_num += num_tables;
								break;
							}
						}
					}
					update_left_buffer(&lptr[i], &rptr[i]);
				}
				else if (right_dist > bucket_width && right_dist > left_dist) {
					count = rptr[i].size_;
					end   = rptr[i].leaf_pos_;
					start = end - count;
					for (int j = end; j > start; j--) {
						id = rptr[i].leaf_node_->get_entry_id(j);
						freq[id]++;
						if (freq[id] > col_threshold && !checked[id]) {
							checked[id] = true;
							read_data(id, dim_, B_, data, data_folder);

							float dist = calc_l2_dist(data, query, dim_);
							int obj_id = id + 1;
							list->insert(dist, obj_id);

							// -------------------------------------------------
							//  Terminating condition 2
							// -------------------------------------------------
							dist_io_++;
							if (dist_io_ >= candidates) {
								again = false;
								flag_num += num_tables;
								break;
							}
						}
					}
					update_right_buffer(&lptr[i], &rptr[i]);
				}
				else {
					flag[i] = false;
					flag_num++;
				}
				if (flag_num >= num_tables) break;
			}
			if (flag_num >= num_tables) break;
		}
		// ---------------------------------------------------------------------
		//  Terminating condition 1 & 2
		// ---------------------------------------------------------------------
		kfn_dist = list->min_key();
		if (kfn_dist > radius/appr_ratio_ || dist_io_ >= candidates) {
			again = false;
			break;
		}

		// ---------------------------------------------------------------------
		//  Step 3: auto-update <radius>
		// ---------------------------------------------------------------------
		//radius = radius / appr_ratio_;
		radius = update_radius(num_tables, lptr, rptr, q_val);
		bucket_width = radius * w_ / 2.0f;
	}

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete[] data; data = NULL;
	delete[] freq; freq = NULL;
	delete[] checked; checked = NULL;
	g_memory -= ((SIZEBOOL + SIZEINT) * n_pts_ + SIZEFLOAT * dim_);

	delete[] q_val; q_val = NULL;
	delete[] flag; flag = NULL;
	g_memory -= (SIZEFLOAT + SIZEBOOL) * num_tables;

	for (int i = 0; i < num_tables; i++) {
		// ---------------------------------------------------------------------
		//  CANNOT remove the condition
		//              <lptrs[i].leaf_node_ != rptrs[i].leaf_node_>
		//  Because <lptrs[i].leaf_node> and <rptrs[i].leaf_node> may point 
		//  to the same address, then we would delete it twice and receive 
		//  the runtime error or segmentation fault.
		// ---------------------------------------------------------------------
		if (lptr[i].leaf_node_ && lptr[i].leaf_node_ != rptr[i].leaf_node_) {
			delete lptr[i].leaf_node_; lptr[i].leaf_node_ = NULL;
		}
		if (rptr[i].leaf_node_) {
			delete rptr[i].leaf_node_; rptr[i].leaf_node_ = NULL;
		}
	}
	delete[] lptr; lptr = NULL;
	delete[] rptr; rptr = NULL;
	g_memory -= (SIZECHAR * B_ * num_tables * 2 + SIZEINT * num_tables * 6);

	return (page_io_ + dist_io_);
}

// -----------------------------------------------------------------------------
void RQALSH::init_buffer(			// init page buffer (loc pos of b-treee)
	int num_tables,						// numb of hash tables used for search
	PageBuffer* lptr,					// left buffer page (return)
	PageBuffer* rptr)					// right buffer page (return)
{
	int block   = -1;				// variables for index node
	BIndexNode* index_node = NULL;
	BLeafNode* leaf_node = NULL;

	int pos = -1;					// variables for leaf node
	int increment = -1;
	int num_entries = -1;
	int num_keys = -1;

	for (int i = 0; i < num_tables; i++) {
		block = trees_[i]->root_;
		if (block == 1) {
			// -----------------------------------------------------------------
			// the b-tree has only one leaf node
			// -----------------------------------------------------------------
			lptr[i].leaf_node_ = new BLeafNode();
			lptr[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			leaf_node   = lptr[i].leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();
			if (num_keys > 1) {
				lptr[i].index_pos_ = 0;
				lptr[i].leaf_pos_  = 0;
				lptr[i].size_      = increment;

				rptr[i].leaf_node_ = lptr[i].leaf_node_;
				rptr[i].index_pos_ = num_keys - 1;
				rptr[i].leaf_pos_  = num_entries - 1;
				rptr[i].size_      = num_entries - (num_keys - 1) * increment;
			}
			else {
				lptr[i].index_pos_ = 0;
				lptr[i].leaf_pos_  = 0;
				lptr[i].size_      = num_entries;

				rptr[i].leaf_node_ = NULL;
				rptr[i].index_pos_ = -1;
				rptr[i].leaf_pos_  = -1;
				rptr[i].size_      = -1;
			}
		}
		else {
			// -----------------------------------------------------------------
			// the b-tree has index node
			// -----------------------------------------------------------------

			// -----------------------------------------------------------------
			//  Initialize left leaf node
			// -----------------------------------------------------------------
			block = trees_[i]->root_;
			index_node = new BIndexNode();
			index_node->init_restore(trees_[i], block);
			page_io_++;
									// find the left most leaf node
			while (index_node->get_level() > 1) {
				block = index_node->get_son(0);
				delete index_node; index_node = NULL;

				index_node = new BIndexNode();
				index_node->init_restore(trees_[i], block);
				page_io_++;			// access a new node (a new page)
			}

			block = index_node->get_son(0);
			lptr[i].leaf_node_ = new BLeafNode();
			lptr[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			lptr[i].index_pos_ = 0;
			lptr[i].leaf_pos_  = 0;

			increment   = lptr[i].leaf_node_->get_increment();
			num_entries = lptr[i].leaf_node_->get_num_entries();
			if (increment > num_entries) {
				lptr[i].size_ = num_entries;
			} else {
				lptr[i].size_ = increment;
			}

			if (index_node != NULL) {
				delete index_node; index_node = NULL;
			}

			// -----------------------------------------------------------------
			//  Initialize right leaf node
			// -----------------------------------------------------------------
			block = trees_[i]->root_;
			index_node = new BIndexNode();
			index_node->init_restore(trees_[i], block);
			page_io_++;
									// find the right most leaf node
			while (index_node->get_level() > 1) {
				num_entries = index_node->get_num_entries();
				block = index_node->get_son(num_entries - 1);
				delete index_node; index_node = NULL;

				index_node = new BIndexNode();
				index_node->init_restore(trees_[i], block);
				page_io_++;			// access a new node (a new page)
			}

			num_entries = index_node->get_num_entries();
			block = index_node->get_son(num_entries - 1);
			rptr[i].leaf_node_ = new BLeafNode();
			rptr[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			leaf_node   = rptr[i].leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();

			rptr[i].index_pos_ = num_keys - 1;
			rptr[i].leaf_pos_  = num_entries - 1;
			rptr[i].size_      = num_entries - (num_keys - 1) * increment;

			if (index_node != NULL) {
				delete index_node; index_node = NULL;
			}
		}
	}
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	int num_tables,						// numb of hash tables used for search
	PageBuffer* lptr,					// left page buffer
	PageBuffer* rptr,					// right page buffer
	float* q_val)						// hash value of query
{
	return update_radius(num_tables, lptr, rptr, q_val);
}

// -----------------------------------------------------------------------------
float RQALSH::update_radius(		// update radius
	int num_tables,						// numb of hash tables used for search
	PageBuffer* lptr,					// left page buffer
	PageBuffer* rptr,					// right page buffer
	float* q_val)						// hash value of query
{
	float dist = 0.0f;				// tmp vars
	vector<float> list;

									// find an array of proj dist on each table
	for (int i = 0; i < num_tables; i++) {
		if (lptr[i].size_ != -1) {
			dist = calc_proj_dist(&lptr[i], q_val[i]);
			list.push_back(dist);
		}
		if (rptr[i].size_ != -1) {
			dist = calc_proj_dist(&rptr[i], q_val[i]);
			list.push_back(dist);
		}
	}
	sort(list.begin(), list.end());	// sort the array

	int num = (int) list.size();
	if (num % 2 == 0) {				// find median dist
		dist = (list[num/2 - 1] + list[num/2]) / 2.0f;
	} else {
		dist = list[num/2];
	}
	list.clear();

	int kappa = (int) floor(log(2.0f * dist / w_) / log(appr_ratio_));
	dist = pow(appr_ratio_, kappa);
	return dist;
}

// -----------------------------------------------------------------------------
void RQALSH::update_left_buffer(	// update left buffer
	PageBuffer* lptr,					// left buffer
	const PageBuffer* rptr)				// right buffer
{
	BLeafNode* leaf_node     = NULL;
	BLeafNode* old_leaf_node = NULL;

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
			page_io_++;
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
	const PageBuffer* lptr,				// left buffer
	PageBuffer* rptr)					// right buffer
{
	BLeafNode* leaf_node     = NULL;
	BLeafNode* old_leaf_node = NULL;

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
			page_io_++;
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
float RQALSH::calc_proj_dist(		// calc proj dist
	const PageBuffer* ptr,				// page buffer
	float q_val)						// hash value of query
{
	int pos = ptr->index_pos_;
	float key = ptr->leaf_node_->get_key(pos);

	return fabs(key - q_val);
}

// -----------------------------------------------------------------------------
//  Comparison function for qsort called by RQALSH::bulkload()
// -----------------------------------------------------------------------------
int HashValueQsortComp(				// compare function for qsort
	const void* e1,						// 1st element
	const void* e2)						// 2nd element
{
	int ret = 0;
	HashValue* value1 = (HashValue*) e1;
	HashValue* value2 = (HashValue*) e2;

	if (value1->proj_ < value2->proj_) {
		ret = -1;
	} else if (value1->proj_ > value2->proj_) {
		ret = 1;
	} else {
		if (value1->id_ < value2->id_) ret = -1;
		else if (value1->id_ > value2->id_) ret = 1;
	}
	return ret;
}
