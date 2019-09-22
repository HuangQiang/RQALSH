#include "headers.h"

// -----------------------------------------------------------------------------
Drusilla_Index::Drusilla_Index()	// default constructor
{
	n_pts_   = -1;
	dim_     = -1;
	l_       = -1;
	m_       = -1;
	B_       = -1;
	fn_cand_ = NULL;
}

// -----------------------------------------------------------------------------
Drusilla_Index::~Drusilla_Index()	// destructor
{
	delete[] fn_cand_; fn_cand_ = NULL;
}

// -----------------------------------------------------------------------------
int Drusilla_Index::build(			// build index
	int   n,							// number of data points
	int   d,							// number of dimensions
	int   l,							// number of projections
	int   m,							// number of candidates on each proj
	int   B,							// page size
	const float **data,					// data objects
	const char *index_path)				// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	l_     = l;
	m_     = m;
	B_     = B;

	strcpy(fname_, index_path);
	strcat(fname_, "drusilla.para");

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
int Drusilla_Index::bulkload(		// build hash tables
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
	vector<vector<float> > shift_data(n_pts_, vector<float>(dim_, 0.0f));
	for (int i = 0; i < n_pts_; ++i) {
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
		for (int j = 0; j < dim_; ++j) {
			float x = shift_data[i][j];
			norm[i] += x * x;
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id = i;
		}
	}

	vector<float> projection(dim_, 0.0f);
	vector<bool>  close_angle(n_pts_, false);
	Result *score_pair = new Result[n_pts_];

	fn_cand_ = new int[l_ * m_];
	for (int i = 0; i < l_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			projection[j] = shift_data[max_id][j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			score_pair[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = 0.0f;
				for (int k = 0; k < dim_; ++k) {
					offset += shift_data[j][k] * projection[k];
				}

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					float x = shift_data[j][k] - offset * projection[k];
					distortion += x * x;
				}
				distortion = sqrt(distortion);

				score_pair[j].key_ = fabs(offset) - fabs(distortion);
				if (atan(distortion / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < FLOATZERO) {
				score_pair[j].key_ = MINREAL + 1.0f;
			}
			else {
				score_pair[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the idects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score_pair, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < m_; ++j) {
			int id = score_pair[j].id_;
			fn_cand_[i * m_ + j] = id;
			
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
	delete[] score_pair; score_pair = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void Drusilla_Index::display()		// display parameters
{
	printf("Parameters of Drusilla_Select (SISAP2016 paper):\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    l     = %d\n",   l_);
	printf("    m     = %d\n",   m_);
	printf("    B     = %d\n",   B_);
	printf("    fname = %s\n\n", fname_);
}

// -----------------------------------------------------------------------------
int Drusilla_Index::write_params()	// write parameters to disk
{
	FILE* fp = fopen(fname_, "w");
	if (!fp) {
		printf("Culd not create %s.\n", fname_);
		return 1;
	}

	fprintf(fp, "n = %d\n", n_pts_);
	fprintf(fp, "d = %d\n", dim_);
	fprintf(fp, "l = %d\n", l_);
	fprintf(fp, "m = %d\n", m_);
	fprintf(fp, "B = %d\n", B_);

	int size = l_ * m_;
	for (int i = 0; i < size; ++i) {
		fprintf(fp, "%d\n", fn_cand_[i]);
	}
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Index::load(			// load index
	const char *index_path)				// index path
{
	strcpy(fname_, index_path);
	strcat(fname_, "drusilla.para");

	// -------------------------------------------------------------------------
	//  read index file from disk
	// -------------------------------------------------------------------------
	if (read_params() == 1) return 1;

	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Index::read_params()	// read parameters from disk
{
	FILE* fp = fopen(fname_, "r");
	if (!fp) {
		fprintf(stderr, "Could not open %s\n", fname_);
		return 1;
	}

	fscanf(fp, "n = %d\n", &n_pts_);
	fscanf(fp, "d = %d\n", &dim_);
	fscanf(fp, "l = %d\n", &l_);
	fscanf(fp, "m = %d\n", &m_);
	fscanf(fp, "B = %d\n", &B_);

	int size = l_ * m_;
	fn_cand_ = new int[size];
	for (int i = 0; i < size; ++i) {
		fscanf(fp, "%d\n", &fn_cand_[i]);
	}
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int Drusilla_Index::search(			// c-k-AFN search
	const float *query,					// query point
	const char *data_folder,			// new format data folder
	MaxK_List *list)					// top-k results (return)
{
	float *data = new float[dim_];	
	int size = l_ * m_;
	for (int i = 0; i < size; ++i) {
		int id = fn_cand_[i];
		read_data_new_format(id, dim_, B_, data_folder, data);

		float dist = calc_l2_dist(dim_, (const float *) data, query);
		list->insert(dist, id + 1);
	}
	delete[] data; data = NULL;

	return size;
}
