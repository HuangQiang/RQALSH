#include "headers.h"

// -----------------------------------------------------------------------------
RQALSH_Star::RQALSH_Star()			// constructor
{
	n_pts_       = -1;
	dim_         = -1;
	B_           = -1;
	L_           = -1;
	M_           = -1;
	beta_        = -1;
	delta_       = -1.0f;
	appr_ratio_  = -1.0f;
	sample_size_ = -1;
	sample_id_   = NULL;
	lsh_         = NULL;
}

// -----------------------------------------------------------------------------
RQALSH_Star::~RQALSH_Star()			// destructor
{
	if (sample_id_ != NULL) {
		delete[] sample_id_; sample_id_ = NULL;
	}
	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int RQALSH_Star::build(				// build index	
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
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
	L_          = L;
	M_          = M;
	beta_       = beta;
	delta_      = delta;
	appr_ratio_ = ratio;

	strcpy(index_path_, index_path);
	create_dir(index_path_);

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload(data);
	// display();
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::bulkload(			// bulkloading for each block
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  calculate shift data
	// -------------------------------------------------------------------------
	float *shift_data = new float[n_pts_ * dim_];
	calc_shift_data(data, shift_data);

	// -------------------------------------------------------------------------
	//  get sample data from data dependent selection
	// -------------------------------------------------------------------------
	sample_size_ = L_ * M_;
	sample_id_   = new int[sample_size_];
	data_dependent_select(shift_data);

	float **sample_data = new float*[sample_size_];
	for (int i = 0; i < sample_size_; ++i) {
		int id = sample_id_[i];
		
		sample_data[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			sample_data[i][j] = data[id][j];
		}
	}

	// -------------------------------------------------------------------------
	//  build hash tables for objects from drusilla select using RQALSH
	// -------------------------------------------------------------------------
	if (sample_size_ > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->build(sample_size_, dim_, B_, beta_, delta_, appr_ratio_, 
			index_path_, (const float **) sample_data);
	}

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	write_params();
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] shift_data; shift_data = NULL;
	for (int i = 0; i < sample_size_; ++i) {
		delete[] sample_data[i]; sample_data[i] = NULL;
	}
	delete[] sample_data; sample_data = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::calc_shift_data( 	// calc shift data
	const float **data,					// data objects
	float *shift_data)					// shift data objects (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) {
		centroid[i] /= (float) n_pts_;
	}

	// -------------------------------------------------------------------------
	//  make a copy of data objects which move to the centroid of data objects
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		int base = i * dim_;
		for (int j = 0; j < dim_; ++j) {
			shift_data[base + j] = data[i][j] - centroid[j];
		}
	}
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::data_dependent_select( // drusilla select
	const float *shift_data)			// shift data
{
	// -------------------------------------------------------------------------
	//  calc the norm of data objects and find the data object with max norm
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = -1.0f;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		int base = i * dim_;
		for (int j = 0; j < dim_; ++j) {
			float x = shift_data[base + j];
			norm[i] += x * x;
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id   = i;
		}
	}

	vector<bool>  close_angle(n_pts_);
	vector<float> projection(dim_);
	Result *score_pair = new Result[n_pts_];

	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			projection[j] = shift_data[max_id * dim_ + j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			int base = j * dim_;

			if (norm[j] >= 0.0f) {
				float offset = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					offset += (shift_data[base + k] * projection[k]);
				}

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					float x = shift_data[base + k] - offset * projection[k];
					distortion += x * x;
				}

				score_pair[j].id_ = j;
				score_pair[j].key_ = offset * offset - distortion;
			}
			else {
				score_pair[j].id_ = j;
				score_pair[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score_pair, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < M_; ++j) {
			int id = score_pair[j].id_;

			sample_id_[i * M_ + j] = id;
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > max_norm) {
				max_norm = norm[j];
				max_id = j;
			}
		}
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] score_pair; score_pair = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH_Star::display()			// display parameters
{
	printf("Parameters of RQALSH_Star:\n");
	printf("    n           = %d\n", n_pts_);
	printf("    d           = %d\n", dim_);
	printf("    B           = %d\n", B_);
	printf("    L           = %d\n", L_);
	printf("    M           = %d\n", M_);
	printf("    beta        = %d\n", beta_);
	printf("    delta       = %f\n", delta_);
	printf("    c           = %f\n", appr_ratio_);
	printf("    sample_size = %d\n", sample_size_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH_Star::write_params()		// write parameters to disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "rqalsh_star_para");
	
	FILE *fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		exit(1);
	}

	fprintf(fp, "n     = %d\n", n_pts_);
	fprintf(fp, "d     = %d\n", dim_);
	fprintf(fp, "B     = %d\n", B_);
	fprintf(fp, "L     = %d\n", L_);
	fprintf(fp, "M     = %d\n", M_);
	fprintf(fp, "beta  = %d\n", beta_);
	fprintf(fp, "delta = %f\n", delta_);
	fprintf(fp, "c     = %f\n", appr_ratio_);

	for (int i = 0; i < sample_size_; ++i) {
		fprintf(fp, "%d ", sample_id_[i]);
	}
	fprintf(fp, "\n");
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::load(				// restore parameters
	const char *index_path)				// output folder
{
	strcpy(index_path_, index_path);	
	if (read_params() == 1) return 1;

	if (sample_size_ > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->load(index_path_);
	}

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::read_params()		// read parameters from disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "rqalsh_star_para");

	FILE* fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s.\n", fname);
		return 1;
	}

	fscanf(fp, "n     = %d\n", &n_pts_);
	fscanf(fp, "d     = %d\n", &dim_);
	fscanf(fp, "B     = %d\n", &B_);
	fscanf(fp, "L     = %d\n", &L_);
	fscanf(fp, "M     = %d\n", &M_);
	fscanf(fp, "beta  = %d\n", &beta_);
	fscanf(fp, "delta = %f\n", &delta_);
	fscanf(fp, "c     = %f\n", &appr_ratio_);

	sample_size_ = L_ * M_;
	sample_id_   = new int[sample_size_];
	for (int i = 0; i < sample_size_; ++i) {
		fscanf(fp, "%d ", &sample_id_[i]);
	}
	fscanf(fp, "\n");
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::kfn(				// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  use index to speed up c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1;
	if (sample_size_ > candidates) {
		return lsh_->kfn(top_k, query, (const int*) sample_id_, data_folder, list);
	}

	// -------------------------------------------------------------------------
	//  if the number of samples is small enough, linear scan directly
	// -------------------------------------------------------------------------
	float *data = new float[dim_];		
	for (int i = 0; i < sample_size_; ++i) {
		int id  = sample_id_[i];
		read_data_new_format(id, dim_, B_, data_folder, data);

		float dist = calc_l2_dist(dim_, (const float*) data, query);
		list->insert(dist, id + 1);
	}
	delete [] data; data = NULL;
	
	return sample_size_;
}
