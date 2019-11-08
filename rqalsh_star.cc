#include <algorithm>
#include <cassert>
#include <cstring>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"
#include "rqalsh.h"
#include "rqalsh_star.h"

// -----------------------------------------------------------------------------
//  RQALSH* is used to solve the problem of c-k-Approximate Furthest Neighbor 
//  (c-k-AFN) search
// -----------------------------------------------------------------------------
RQALSH_Star::RQALSH_Star()			// default constructor
{
	n_pts_       = -1;
	dim_         = -1;
	B_           = -1;
	L_           = -1;
	M_           = -1;
	beta_        = -1;
	delta_       = -1.0f;
	appr_ratio_  = -1.0f;
	n_cand_      = -1;
	cand_        = NULL;
	lsh_         = NULL;
}

// -----------------------------------------------------------------------------
RQALSH_Star::~RQALSH_Star()			// destructor
{
	delete[] cand_; cand_ = NULL;
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
	beta_       = beta;
	delta_      = delta;
	appr_ratio_ = ratio;

	strcpy(path_, path);
	strcat(path_, "indices/");
	create_dir(path_);

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload(data);
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::bulkload(			// bulkloading for each block
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  calculate shift data
	// -------------------------------------------------------------------------
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) shift_data[i] = new float[dim_];
	calc_shift_data(data, shift_data);

	// -------------------------------------------------------------------------
	//  get sample data from data dependent selection
	// -------------------------------------------------------------------------
	n_cand_ = L_ * M_;
	cand_   = new int[n_cand_];
	data_dependent_select((const float **) shift_data);

	float **cand_data = new float*[n_cand_];
	for (int i = 0; i < n_cand_; ++i) {
		int id = cand_[i];
		
		cand_data[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			cand_data[i][j] = data[id][j];
		}
	}

	// -------------------------------------------------------------------------
	//  build hash tables for objects from drusilla select using RQALSH
	// -------------------------------------------------------------------------
	if (n_cand_ > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->build(n_cand_, dim_, B_, beta_, delta_, appr_ratio_, 
			(const float **) cand_data, path_);
	}

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	write_params();
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;

	for (int i = 0; i < n_cand_; ++i) {
		delete[] cand_data[i]; cand_data[i] = NULL;
	}
	delete[] cand_data; cand_data = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::calc_shift_data( 	// calc shift data
	const float **data,					// data objects
	float **shift_data)					// shift data objects (return)
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
	for (int i = 0; i < dim_; ++i) {
		centroid[i] /= (float) n_pts_;
	}

	// -------------------------------------------------------------------------
	//  make a copy of data objects which move to the centroid of data objects
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			shift_data[i][j] = data[i][j] - centroid[j];
		}
	}
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::data_dependent_select( // drusilla select
	const float **shift_data)			// shift data
{
	// -------------------------------------------------------------------------
	//  calc the norm of data objects and find the data object with max norm
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = -1.0f;
	std::vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, shift_data[i], shift_data[i]));
		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id   = i;
		}
	}

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

			cand_[i * M_ + j] = id;
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
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH_Star::display()			// display parameters
{
	printf("Parameters of RQALSH*:\n");
	printf("    n      = %d\n",   n_pts_);
	printf("    d      = %d\n",   dim_);
	printf("    B      = %d\n",   B_);
	printf("    L      = %d\n",   L_);
	printf("    M      = %d\n",   M_);
	printf("    beta   = %d\n",   beta_);
	printf("    delta  = %f\n",   delta_);
	printf("    c      = %.1f\n", appr_ratio_);
	printf("    n_cand = %d\n\n", n_cand_);
}

// -----------------------------------------------------------------------------
int RQALSH_Star::write_params()		// write parameters to disk
{
	char fname[200];
	strcpy(fname, path_);
	strcat(fname, "rqalsh_star_para");
	
	FILE *fp = fopen(fname, "wb");
	if (!fp) {
		printf("Could not create %s\n", fname);
		return 1;
	}

	fwrite(&n_pts_,      SIZEINT,   1,       fp);
	fwrite(&dim_,        SIZEINT,   1,       fp);
	fwrite(&B_,          SIZEINT,   1,       fp);
	fwrite(&L_,          SIZEINT,   1,       fp);
	fwrite(&M_,          SIZEINT,   1,       fp);
	fwrite(&beta_,       SIZEINT,   1,       fp);
	fwrite(&delta_,      SIZEFLOAT, 1,       fp);
	fwrite(&appr_ratio_, SIZEFLOAT, 1,       fp);
	fwrite(cand_,        SIZEINT,   n_cand_, fp);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::load(				// restore parameters
	const char *path)					// index path
{
	strcpy(path_, path);
	strcat(path_, "indices/");
	if (read_params() == 1) return 1;

	if (n_cand_ > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->load(path_);
	}

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::read_params()		// read parameters from disk
{
	char fname[200];
	strcpy(fname, path_);
	strcat(fname, "rqalsh_star_para");

	FILE* fp = fopen(fname, "rb");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	fread(&n_pts_,      SIZEINT,   1, fp);
	fread(&dim_,        SIZEINT,   1, fp);
	fread(&B_,          SIZEINT,   1, fp);
	fread(&L_,          SIZEINT,   1, fp);
	fread(&M_,          SIZEINT,   1, fp);
	fread(&beta_,       SIZEINT,   1, fp);
	fread(&delta_,      SIZEFLOAT, 1, fp);
	fread(&appr_ratio_, SIZEFLOAT, 1, fp);
	
	n_cand_ = L_ * M_;
	cand_   = new int[n_cand_];
	fread(cand_, SIZEINT, n_cand_, fp);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
long long RQALSH_Star::kfn(			// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  use index to speed up c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1;
	if (n_cand_ > candidates) {
		return lsh_->kfn(top_k, query, (const int*) cand_, data_folder, list);
	}

	// -------------------------------------------------------------------------
	//  if the number of samples is small enough, linear scan directly
	// -------------------------------------------------------------------------
	float *data = new float[dim_];		
	for (int i = 0; i < n_cand_; ++i) {
		int id  = cand_[i];
		read_data_new_format(id, dim_, B_, data_folder, data);

		float dist = calc_l2_dist(dim_, (const float*) data, query);
		list->insert(dist, id + 1);
	}
	delete[] data; data = NULL;
	
	return (long long) n_cand_;
}
