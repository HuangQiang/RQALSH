#include "rqalsh_star.h"

// -----------------------------------------------------------------------------
RQALSH_STAR::RQALSH_STAR()			// default constructor
{
	n_pts_ = -1;
	dim_   = -1;
	B_     = -1;
	L_     = -1;
	M_     = -1;
	cand_  = NULL;
	lsh_   = NULL;
}

// -----------------------------------------------------------------------------
RQALSH_STAR::~RQALSH_STAR()			// destructor
{
	delete[] cand_; cand_ = NULL; 
	g_memory -= SIZEINT * L_ * M_;
	
	if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::build(				// build index	
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char  *path)					// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	B_          = B;
	L_          = L;
	M_          = M;
	strcpy(path_, path); strcat(path_, "indices/");
	create_dir(path_);

	// -------------------------------------------------------------------------
	//  get representative data from data dependent selection
	// -------------------------------------------------------------------------
	int n_cand = L_ * M_;
	g_memory += SIZEINT * n_cand;
	cand_ = new int[n_cand];

	data_dependent_select(data, cand_);

	// -------------------------------------------------------------------------
	//  build rqalsh for representative data if necessary
	// -------------------------------------------------------------------------
	if (n_cand > CANDIDATES) {
		float **cand_data = new float*[n_cand];
		for (int i = 0; i < n_cand; ++i) {
			cand_data[i] = new float[dim_];

			int id = cand_[i];
			for (int j = 0; j < dim_; ++j) cand_data[i][j] = data[id][j];
		}

		lsh_ = new RQALSH();
		lsh_->build(n_cand, dim_, B_, beta, delta, ratio, 
			(const float **) cand_data, path_);

		for (int i = 0; i < n_cand; ++i) {
			delete[] cand_data[i]; cand_data[i] = NULL;
		}
		delete[] cand_data; cand_data = NULL;
	}

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, path_); strcat(fname, "rqalsh_star_para");
	
	FILE *fp = fopen(fname, "wb");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }

	fwrite(&n_pts_, SIZEINT, 1,      fp);
	fwrite(&dim_,   SIZEINT, 1,      fp);
	fwrite(&B_,     SIZEINT, 1,      fp);
	fwrite(&L_,     SIZEINT, 1,      fp);
	fwrite(&M_,     SIZEINT, 1,      fp);
	fwrite(cand_,   SIZEINT, n_cand, fp);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::data_dependent_select( // data dependent selection
	const float **data,					// data objects
	int   *cand)						// candidate id (return)
{
	// -------------------------------------------------------------------------
	//  calc the shift data
	// -------------------------------------------------------------------------
	int   max_id = -1;
	float max_norm = -1.0f;
	float *norm = new float[n_pts_];
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) shift_data[i] = new float[dim_];

	calc_shift_data(data, max_id, max_norm, norm, shift_data);

	// -------------------------------------------------------------------------
	//  data dependent selection
	// -------------------------------------------------------------------------
	float  *proj  = new float[dim_];
	Result *score = new Result[n_pts_];

	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			proj[j] = shift_data[max_id][j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] >= 0.0f) {
				float offset = calc_inner_product(dim_, shift_data[j], proj);

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j][k] - offset * proj[k]);
				}
				score[j].id_  = j;
				score[j].key_ = offset * offset - distortion;
			}
			else {
				score[j].id_  = j;
				score[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < M_; ++j) {
			int id = score[j].id_;

			cand[i * M_ + j] = id;
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > max_norm) { max_norm = norm[j]; max_id = j; }
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;  norm  = NULL;
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::calc_shift_data( 	// calculate shift data objects
	const float **data,					// data objects
	int   &max_id,						// data id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float **shift_data) 				// shift data (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	std::vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid[i] /= n_pts_;

	// -------------------------------------------------------------------------
	//  calc shift data and their l2-norm and find max l2-norm and its id
	// -------------------------------------------------------------------------
	max_id   = -1;
	max_norm = MINREAL;

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = 0.0f;
		for (int j = 0; j < dim_; ++j) {
			float tmp = data[i][j] - centroid[j];
			shift_data[i][j] = tmp;
			norm[i] += SQR(tmp);
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::display()			// display parameters
{
	printf("Parameters of RQALSH*:\n");
	printf("    n    = %d\n", n_pts_);
	printf("    d    = %d\n", dim_);
	printf("    B    = %d\n", B_);
	printf("    L    = %d\n", L_);
	printf("    M    = %d\n", M_);
	printf("    path = %s\n", path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::load(				// restore parameters
	const char *path)					// index path
{
	strcpy(path_, path); strcat(path_, "indices/");

	// -------------------------------------------------------------------------
	//  read parameters from disk
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, path_); strcat(fname, "rqalsh_star_para");

	FILE* fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fread(&n_pts_, SIZEINT, 1, fp);
	fread(&dim_,   SIZEINT, 1, fp);
	fread(&B_,     SIZEINT, 1, fp);
	fread(&L_,     SIZEINT, 1, fp);
	fread(&M_,     SIZEINT, 1, fp);
	
	int n_cand = L_ * M_;
	g_memory += SIZEINT * n_cand;
	cand_ = new int[n_cand];
	fread(cand_, SIZEINT, n_cand, fp);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  load rqalsh if necessary
	// -------------------------------------------------------------------------
	if (n_cand > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->load(path_);
	}
	return 0;
}

// -----------------------------------------------------------------------------
uint64_t RQALSH_STAR::kfn(			// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  use index to speed up c-k-AFN search
	// -------------------------------------------------------------------------
	int n_cand = L_ * M_;
	int candidates = CANDIDATES + top_k - 1;

	if (n_cand > candidates) {
		return lsh_->kfn(top_k, query, (const int*) cand_, data_folder, list);
	}

	// -------------------------------------------------------------------------
	//  otherwise, linear scan directly
	// -------------------------------------------------------------------------
	float *data = new float[dim_];		
	for (int i = 0; i < n_cand; ++i) {
		int id  = cand_[i];
		read_data_new_format(id, dim_, B_, data_folder, data);

		float dist = calc_l2_dist(dim_, (const float*) data, query);
		list->insert(dist, id + 1);
	}
	delete[] data; data = NULL;
	
	return (uint64_t) n_cand;
}
