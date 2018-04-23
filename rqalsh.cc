#include "headers.h"


// -----------------------------------------------------------------------------
//  RQALSH: structure of rqalsh indexed by b+ tree. RQALSH is used to solve
//  the problem of c-Approximate Furthest Neighbor (c-AFN) search.
// -----------------------------------------------------------------------------
RQALSH::RQALSH()					// constructor
{
	n_pts_      = -1;
	dim_        = -1;
	B_          = -1;
	beta_       = -1.0f;
	delta_      = -1.0f;
	appr_ratio_ = -1.0f;
	
	w_          = -1.0f;
	p1_         = -1.0f;
	p2_         = -1.0f;
	alpha_      = -1.0f;
	beta_       = -1.0f;
	delta_      = -1.0f;
	m_          = -1;
	l_          = -1;
	a_array_    = NULL;
	trees_      = NULL;

	dist_io_    = -1;
	page_io_    = -1;
	freq_       = NULL;
	checked_    = NULL;
	flag_       = NULL;
	data_       = NULL;
	q_val_      = NULL;
	lptr_       = NULL;
	rptr_       = NULL;
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	if (a_array_ != NULL) {
		delete[] a_array_; a_array_ = NULL;
	}
	if (trees_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete trees_[i]; trees_[i] = NULL;
		}
		delete[] trees_; trees_ = NULL;
	}

	if (freq_ != NULL) {
		delete[] freq_; freq_ = NULL;
	}
	if (checked_ != NULL) {
		delete[] checked_; checked_ = NULL;
	}
	if (flag_ != NULL) {
		delete[] flag_; flag_ = NULL;
	}
	if (data_ != NULL) {
		delete[] data_; data_ = NULL;
	}
	if (q_val_ != NULL) {
		delete[] q_val_; q_val_ = NULL;
	}
	if (lptr_ != NULL) {
		delete[] lptr_; lptr_ = NULL;
	}
	if (rptr_ != NULL) {
		delete[] rptr_; rptr_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int RQALSH::build(					// build index
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *index_path,				// index path
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	B_          = B;
	beta_       = (float) beta / (float) n;
	delta_      = delta;
	appr_ratio_ = ratio;

	strcpy(index_path_, index_path);
	create_dir(index_path_);

	// -------------------------------------------------------------------------
	//  calc parameters and generate hash functions
	// -------------------------------------------------------------------------
	calc_params();
	gen_hash_func();

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	bulkload(data);
	// display();

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH::calc_params()			// calc params of rqalsh
{
	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	// -------------------------------------------------------------------------
	w_  = sqrt((8.0f*log(appr_ratio_)) / (appr_ratio_*appr_ratio_-1.0f));
	
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
	
	freq_    = new int[n_pts_];
	checked_ = new bool[n_pts_];
	flag_    = new bool[m_];
	data_    = new float[dim_];
	q_val_   = new float[m_];
	lptr_    = new PageBuffer[m_];
	rptr_    = new PageBuffer[m_];
}

// -----------------------------------------------------------------------------
float RQALSH::calc_l2_prob(			// calc prob <p1_> and <p2_> of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
void RQALSH::display()				// display parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    B          = %d\n", B_);
	printf("    ratio      = %f\n", appr_ratio_);
	printf("    w          = %f\n", w_);
	printf("    p1         = %f\n", p1_);
	printf("    p2         = %f\n", p2_);
	printf("    alpha      = %f\n", alpha_);
	printf("    beta       = %f\n", beta_);
	printf("    delta      = %f\n", delta_);
	printf("    m          = %d\n", m_);
	printf("    l          = %d\n", l_);
	printf("\n");
}

// -----------------------------------------------------------------------------
void RQALSH::gen_hash_func()		// generate hash functions
{
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; ++i) {
		a_array_[i] = gaussian(0.0f, 1.0f);
	}
}

// -----------------------------------------------------------------------------
int RQALSH::bulkload(				// build m b-trees by bulkloading
	const float** data)					// data set
{
	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1;

	// -------------------------------------------------------------------------
	//  write hash tables (indexed by B+ Tree) to disk
	// -------------------------------------------------------------------------
	char fname[200];
	Result *table = new Result[n_pts_];
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < n_pts_; ++j) {
			table[j].id_  = j;
			table[j].key_ = calc_hash_value(i, data[j]);
		}

		qsort(table, n_pts_, sizeof(Result), ResultComp);

		BTree *bt = new BTree();
		get_tree_filename(i, fname);
		bt->init(B_, fname);
		if (bt->bulkload(n_pts_, (const Result *) table)) return 1;

		delete bt; bt = NULL;
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] table; table = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::write_params()			// write parameters to disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (fp)	{
		printf("Hash tables exist.\n");
		exit(1);
	}

	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		printf("Perhaps no such folder %s?\n", index_path_);
		return 1;
	}

	fprintf(fp, "n     = %d\n", n_pts_);
	fprintf(fp, "d     = %d\n", dim_);
	fprintf(fp, "B     = %d\n", B_);
	fprintf(fp, "ratio = %f\n", appr_ratio_);
	fprintf(fp, "w     = %f\n", w_);
	fprintf(fp, "p1    = %f\n", p1_);
	fprintf(fp, "p2    = %f\n", p2_);
	fprintf(fp, "alpha = %f\n", alpha_);
	fprintf(fp, "beta  = %f\n", beta_);
	fprintf(fp, "delta = %f\n", delta_);
	fprintf(fp, "m     = %d\n", m_);
	fprintf(fp, "l     = %d\n", l_);

	int count = 0;
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			fprintf(fp, "%f ", a_array_[count++]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);	

	return 0;
}

// -----------------------------------------------------------------------------
float RQALSH::calc_hash_value(		// calc hash value
	int table_id,						// hash table id
	const float *point)					// one data object
{
	float ret = 0.0f;
	for (int i = 0; i < dim_; ++i) {
		ret += (a_array_[table_id * dim_ + i] * point[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
void RQALSH::get_tree_filename(		// get file name of b-tree
	int tree_id,						// tree id, from 0 to m-1
	char *fname)						// file name (return)
{
	sprintf(fname, "%s%d.rqalsh", index_path_, tree_id);
}

// -----------------------------------------------------------------------------
int RQALSH::load(					// load index
	const char *index_path)				// index path
{
	strcpy(index_path_, index_path);
	if (read_params()) return 1;

	char fname[200];
	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; ++i) {
		get_tree_filename(i, fname);

		trees_[i] = new BTree();
		trees_[i]->init_restore(fname);
	}

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::read_params()			// read parameters from disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s.\n", fname);
		return 1;
	}

	fscanf(fp, "n     = %d\n", &n_pts_);
	fscanf(fp, "d     = %d\n", &dim_);
	fscanf(fp, "B     = %d\n", &B_);
	fscanf(fp, "ratio = %f\n", &appr_ratio_);
	fscanf(fp, "w     = %f\n", &w_);
	fscanf(fp, "p1    = %f\n", &p1_);
	fscanf(fp, "p2    = %f\n", &p2_);
	fscanf(fp, "alpha = %f\n", &alpha_);
	fscanf(fp, "beta  = %f\n", &beta_);
	fscanf(fp, "delta = %f\n", &delta_);
	fscanf(fp, "m     = %d\n", &m_);
	fscanf(fp, "l     = %d\n", &l_);
	
	a_array_ = new float[m_ * dim_];
	int count = 0;
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			fscanf(fp, "%f ", &a_array_[count++]);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	freq_    = new int[n_pts_];
	checked_ = new bool[n_pts_];
	flag_    = new bool[m_];
	data_    = new float[dim_];
	q_val_   = new float[m_];
	lptr_    = new PageBuffer[m_];
	rptr_    = new PageBuffer[m_];

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-NN results (return)
{
	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	init_search_params(query);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   candidates = CANDIDATES + top_k - 1; // threshold of candidates
	int   num_flag   = -1;			// used for bucket width
	int   id         = -1;			// data object id
	int   count      = -1;			// count size in one page
	int   start      = -1;			// start position
	int   end        = -1;			// end position

	float dist       = -1.0f;		// real distance between data and query
	float ldist      = -1.0f;		// left distance with query
	float rdist      = -1.0f;		// right distance with query
	float kfn_dist   = MINREAL;		// k-th furthest neighbor distance

	float radius     = find_radius(q_val_, lptr_, rptr_);
	float bucket     = w_ * radius / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		num_flag = 0;
		for (int i = 0; i < m_; ++i) flag_[i] = true;

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects (dynamic separation counting)
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!flag_[i]) continue;

				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				ldist = -1.0f;
				if (lptr_[i].size_ != -1) ldist = calc_dist(q_val_[i], &lptr_[i]);

				rdist = -1.0f;
				if (rptr_[i].size_ != -1) rdist = calc_dist(q_val_[i], &rptr_[i]);

				// -------------------------------------------------------------
				//  step 2.2: determine the closer direction (left or right)
				//  and do separation counting to find frequent objects.
				//
				//  For the frequent object, we calc the Lp distance with
				//  query, and update the c-k-AFN results.
				// -------------------------------------------------------------
				if (ldist > bucket && ldist > rdist) {
					count = lptr_[i].size_;
					start = lptr_[i].leaf_pos_;
					end   = start + count;
					for (int j = start; j < end; ++j) {
						id = lptr_[i].leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data_);

							dist = calc_l2_dist(dim_, (const float*) data_, query);
							kfn_dist = list->insert(dist, id + 1);

							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(&rptr_[i], &lptr_[i]);
				}
				else if (rdist > bucket && ldist <= rdist) {
					count = rptr_[i].size_;
					end   = rptr_[i].leaf_pos_;
					start = end - count;
					for (int j = end; j > start; --j) {
						id = rptr_[i].leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data_);

							dist = calc_l2_dist(dim_, (const float*) data_, query);
							kfn_dist = list->insert(dist, id + 1);

							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(&lptr_[i], &rptr_[i]);
				}
				else {
					flag_[i] = false;
					num_flag++;
				}
				if (num_flag >= m_ || dist_io_ >= candidates) break;
			}
			if (num_flag >= m_ || dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kfn_dist > radius / appr_ratio_ && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	delete_tree_ptr();	

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	const int *object_id,				// object id mapping
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	init_search_params(query);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   candidates = CANDIDATES + top_k - 1; // threshold of candidates
	int   num_flag   = -1;			// used for bucket width
	int   id         = -1;			// current object id
	int   oid        = -1;			// read data object id
	int   count      = -1;			// count size in one page
	int   start      = -1;			// start position
	int   end        = -1;			// end position

	float dist       = -1.0f;		// real distance between data and query
	float ldist      = -1.0f;		// left distance with query
	float rdist      = -1.0f;		// right distance with query
	float kfn_dist   = MINREAL;		// k-th furthest neighbor distance

	float radius     = find_radius(q_val_, lptr_, rptr_);
	float bucket     = w_ * radius / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		num_flag = 0;
		for (int i = 0; i < m_; ++i) flag_[i] = true;

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects (dynamic separation counting)
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!flag_[i]) continue;

				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				ldist = -1.0f;
				if (lptr_[i].size_ != -1) ldist = calc_dist(q_val_[i], &lptr_[i]);

				rdist = -1.0f;
				if (rptr_[i].size_ != -1) rdist = calc_dist(q_val_[i], &rptr_[i]);

				// -------------------------------------------------------------
				//  step 2.2: determine the closer direction (left or right)
				//  and do separation counting to find frequent objects.
				//
				//  For the frequent object, we calc the Lp distance with
				//  query, and update the c-k-AFN results.
				// -------------------------------------------------------------
				if (ldist > bucket && ldist > rdist) {
					count = lptr_[i].size_;
					start = lptr_[i].leaf_pos_;
					end   = start + count;
					for (int j = start; j < end; ++j) {
						id = lptr_[i].leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							oid = object_id[id];
							assert(oid >= 0 && oid < 59000);
							read_data_new_format(oid, dim_, B_, data_folder, data_);

							dist = calc_l2_dist(dim_, (const float*) data_, query);
							kfn_dist = list->insert(dist, oid + 1);

							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(&rptr_[i], &lptr_[i]);
				}
				else if (rdist > bucket && ldist <= rdist) {
					count = rptr_[i].size_;
					end   = rptr_[i].leaf_pos_;
					start = end - count;
					for (int j = end; j > start; --j) {
						id = rptr_[i].leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							oid = object_id[id];
							assert(oid >= 0 && oid < 59000);
							read_data_new_format(oid, dim_, B_, data_folder, data_);

							dist = calc_l2_dist(dim_, (const float*) data_, query);
							kfn_dist = list->insert(dist, oid + 1);

							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(&lptr_[i], &rptr_[i]);
				}
				else {
					flag_[i] = false;
					num_flag++;
				}
				if (num_flag >= m_ || dist_io_ >= candidates) break;
			}
			if (num_flag >= m_ || dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kfn_dist > radius / appr_ratio_ && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	delete_tree_ptr();	
	
	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
void RQALSH::init_search_params(	// init parameters for k-FN search
	const float *query)					// query object
{
	page_io_ = 0;
	dist_io_ = 0;

	for (int i = 0; i < n_pts_; ++i) {
		freq_[i]    = 0;
		checked_[i] = false;
	}

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);

		lptr_[i].leaf_node_ = NULL;
		lptr_[i].index_pos_ = -1;
		lptr_[i].leaf_pos_  = -1;
		lptr_[i].size_      = -1;

		rptr_[i].leaf_node_ = NULL;
		rptr_[i].index_pos_ = -1;
		rptr_[i].leaf_pos_  = -1;
		rptr_[i].size_      = -1;
	}

	BIndexNode *index_node = NULL;
	BLeafNode  *leaf_node = NULL;

	int block       = -1;			// variables for index node
	int pos         = -1;			// variables for leaf node
	int increment   = -1;
	int num_entries = -1;
	int num_keys    = -1;

	for (int i = 0; i < m_; ++i) {
		block = trees_[i]->root_;
		if (block == 1) {
			// -----------------------------------------------------------------
			//  the b-tree has only one leaf node
			// -----------------------------------------------------------------
			lptr_[i].leaf_node_ = new BLeafNode();
			lptr_[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			leaf_node   = lptr_[i].leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();
			if (num_keys > 1) {
				lptr_[i].index_pos_ = 0;
				lptr_[i].leaf_pos_  = 0;
				lptr_[i].size_      = increment;

				rptr_[i].leaf_node_ = lptr_[i].leaf_node_;
				rptr_[i].index_pos_ = num_keys - 1;
				rptr_[i].leaf_pos_  = num_entries - 1;
				rptr_[i].size_      = num_entries - (num_keys - 1) * increment;
			}
			else {
				lptr_[i].index_pos_ = 0;
				lptr_[i].leaf_pos_  = 0;
				lptr_[i].size_      = num_entries;

				rptr_[i].leaf_node_ = NULL;
				rptr_[i].index_pos_ = -1;
				rptr_[i].leaf_pos_  = -1;
				rptr_[i].size_      = -1;
			}
		}
		else {
			// -----------------------------------------------------------------
			//  the b-tree has index node
			// -----------------------------------------------------------------

			// -----------------------------------------------------------------
			//  (1) initialize left leaf node
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
			lptr_[i].leaf_node_ = new BLeafNode();
			lptr_[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			lptr_[i].index_pos_ = 0;
			lptr_[i].leaf_pos_  = 0;

			increment   = lptr_[i].leaf_node_->get_increment();
			num_entries = lptr_[i].leaf_node_->get_num_entries();
			if (increment > num_entries) {
				lptr_[i].size_ = num_entries;
			} else {
				lptr_[i].size_ = increment;
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
			rptr_[i].leaf_node_ = new BLeafNode();
			rptr_[i].leaf_node_->init_restore(trees_[i], block);
			page_io_++;

			leaf_node   = rptr_[i].leaf_node_;
			num_keys    = leaf_node->get_num_keys();
			increment   = leaf_node->get_increment();
			num_entries = leaf_node->get_num_entries();

			rptr_[i].index_pos_ = num_keys - 1;
			rptr_[i].leaf_pos_  = num_entries - 1;
			rptr_[i].size_      = num_entries - (num_keys - 1) * increment;

			if (index_node != NULL) {
				delete index_node; index_node = NULL;
			}
		}
	}
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	const float *q_val,					// hash value of query
	const PageBuffer *lptr,				// left page buffer
	const PageBuffer *rptr)				// right page buffer
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lptr[i].size_ != -1) {
			list.push_back(calc_dist(q_val[i], &lptr[i]));
		}
		if (rptr[i].size_ != -1) {
			list.push_back(calc_dist(q_val[i], &rptr[i]));
		}
	}

	// -------------------------------------------------------------------------
	//  sort the array in ascending order 
	// -------------------------------------------------------------------------
	sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int   num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(appr_ratio_));
	return pow(appr_ratio_, kappa);
}

// -----------------------------------------------------------------------------
void RQALSH::update_left_buffer(	// update left buffer
	const PageBuffer *rptr,				// right buffer
	PageBuffer *lptr)					// left buffer (return)
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
	PageBuffer* rptr)					// right buffer (return)
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
float RQALSH::calc_dist(			// calc projected distance
	float q_val,						// hash value of query
	const PageBuffer *ptr)				// page buffer
{
	int   pos  = ptr->index_pos_;
	float key  = ptr->leaf_node_->get_key(pos);
	float dist = fabs(key - q_val);

	return dist;
}

// -----------------------------------------------------------------------------
void RQALSH::delete_tree_ptr()		// delete the pointers of B+ Trees
{
	for (int i = 0; i < m_; ++i) {
		// ---------------------------------------------------------------------
		//  CANNOT remove the condition
		//              <lptrs[i].leaf_node != rptrs[i].leaf_node>
		//  because <lptrs[i].leaf_node> and <rptrs[i].leaf_node> may point 
		//  to the same address, then we would delete it twice and receive 
		//  the runtime error or segmentation fault.
		// ---------------------------------------------------------------------
		if (lptr_[i].leaf_node_ && lptr_[i].leaf_node_ != rptr_[i].leaf_node_) {
			delete lptr_[i].leaf_node_; lptr_[i].leaf_node_ = NULL;
		}
		if (rptr_[i].leaf_node_) {
			delete rptr_[i].leaf_node_; rptr_[i].leaf_node_ = NULL;
		}
	}
}