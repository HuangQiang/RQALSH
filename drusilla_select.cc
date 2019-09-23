#include "headers.h"

// -----------------------------------------------------------------------------
Drusilla_Select::Drusilla_Select()	// default constructor
{
	n_pts_ = -1;
	dim_   = -1;
	l_     = -1;
	m_     = -1;
	B_     = -1;
	cand_  = NULL;
}

// -----------------------------------------------------------------------------
Drusilla_Select::~Drusilla_Select()	// destructor
{
	delete[] cand_; cand_ = NULL;
}

// -----------------------------------------------------------------------------
int Drusilla_Select::build(			// build index
	int   n,							// number of data points
	int   d,							// number of dimensions
	int   l,							// number of projections
	int   m,							// number of candidates on each proj
	int   B,							// page size
	const float **data,					// data objects
	const char  *path)					// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	l_     = l;
	m_     = m;
	B_     = B;

	strcpy(path_, path);
	strcat(path_, "drusilla.index");

	// -------------------------------------------------------------------------
	//  build hash tables
	// -------------------------------------------------------------------------
	bulkload(data);

	// -------------------------------------------------------------------------
	//  write parameter to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1; 

	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Select::bulkload(		// build hash tables
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  calculate centroid
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
	//  calc the shift data
	// -------------------------------------------------------------------------
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		shift_data[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			shift_data[i][j] = data[i][j] - centroid[j];
		}
	}

	// -------------------------------------------------------------------------
	//  find the idect with maximum Euclidean norm
	// -------------------------------------------------------------------------
	int   max_id = -1;
	float max_norm = -1.0f;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, shift_data[i], shift_data[i]));
		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id = i;
		}
	}

	vector<bool>  close_angle(n_pts_, false);
	Result *score = new Result[n_pts_];
	float  *proj  = new float[dim_];

	cand_ = new int[l_ * m_];
	for (int i = 0; i < l_; ++i) {
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
			score[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = calc_inner_product(dim_, shift_data[j], proj);

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j][k] - offset * proj[k]);
				}
				distortion = sqrt(distortion);

				score[j].key_ = fabs(offset) - fabs(distortion);
				if (atan(distortion / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < FLOATZERO) {
				score[j].key_ = MINREAL + 1.0f;
			}
			else {
				score[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the idects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < m_; ++j) {
			int id = score[j].id_;
			cand_[i * m_ + j] = id;
			
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding idect
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) {
				norm[j] = 0.0f;
			}
			if (norm[j] > max_norm) {
				max_norm = norm[j];
				max_id = j;
			}
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void Drusilla_Select::display()		// display parameters
{
	printf("Parameters of Drusilla_Select (SISAP2016 paper):\n");
	printf("    n    = %d\n",   n_pts_);
	printf("    d    = %d\n",   dim_);
	printf("    l    = %d\n",   l_);
	printf("    m    = %d\n",   m_);
	printf("    B    = %d\n",   B_);
	printf("    path = %s\n\n", path_);
}

// -----------------------------------------------------------------------------
int Drusilla_Select::write_params()	// write parameters to disk
{
	FILE *fp = fopen(path_, "wb");
	if (!fp) {
		printf("Culd not create %s\n", path_);
		return 1;
	}

	fwrite(&n_pts_, SIZEINT, 1,     fp);
	fwrite(&dim_,   SIZEINT, 1,     fp);
	fwrite(&B_,     SIZEINT, 1,     fp);
	fwrite(&l_,     SIZEINT, 1,     fp);
	fwrite(&m_,     SIZEINT, 1,     fp);
	fwrite(cand_,   SIZEINT, l_*m_, fp);
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Select::load(			// load index
	const char *path)					// index path
{
	strcpy(path_, path);
	strcat(path_, "drusilla.index");

	// -------------------------------------------------------------------------
	//  read index file from disk
	// -------------------------------------------------------------------------
	if (read_params() == 1) return 1;

	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Select::read_params()	// read parameters from disk
{
	FILE *fp = fopen(path_, "rb");
	if (!fp) {
		printf("Could not open %s\n", path_);
		return 1;
	}

	fread(&n_pts_, SIZEINT, 1, fp);
	fread(&dim_,   SIZEINT, 1, fp);
	fread(&B_,     SIZEINT, 1, fp);
	fread(&l_,     SIZEINT, 1, fp);
	fread(&m_,     SIZEINT, 1, fp);

	int size = l_ * m_;
	cand_ = new int[size];
	fread(cand_, SIZEINT, size, fp);
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
long long Drusilla_Select::search(	// c-k-AFN search
	const float *query,					// query point
	const char  *data_folder,			// new format data folder
	MaxK_List   *list)					// top-k results (return)
{
	float *data = new float[dim_];	
	long long size = l_ * m_;
	for (int i = 0; i < size; ++i) {
		int id = cand_[i];
		read_data_new_format(id, dim_, B_, data_folder, data);

		float dist = calc_l2_dist(dim_, (const float *) data, query);
		list->insert(dist, id + 1);
	}
	delete[] data; data = NULL;

	return size;
}
