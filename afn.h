#ifndef __AFN_H
#define __AFN_H

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	const float  **query,				// query set
	const Result **R,					// truth set
	const char   *data_folder,			// data folder
	const char   *output_folder);		// output folder

// -----------------------------------------------------------------------------
int indexing_of_rqalsh_star(		// indexing of RQALSH*
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int kfn_of_rqalsh_star(				// c-k-AFN search of RQALSH*
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float  **query,				// query set
	const Result **R,					// truth set
	const char   *data_folder,			// data folder
	const char   *output_folder);		// output folder

// -----------------------------------------------------------------------------
int indexing_of_rqalsh(				// indexing of RQALSH
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int kfn_of_rqalsh(					// c-k-AFN search of RQALSH
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float  **query,				// query set
	const Result **R,					// truth set
	const char   *data_folder,			// data folder
	const char   *output_folder);		// output folder

// -----------------------------------------------------------------------------
int indexing_of_drusilla_select(	// indexing of Drusilla_Select
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	const float **data,					// data set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int kfn_of_drusilla_select(			// c-k-AFN via Drusilla_Select
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float  **query,				// query set
	const Result **R,					// truth set
	const char   *data_folder,			// data folder
	const char   *output_folder);		// output folder

// -----------------------------------------------------------------------------
int indexing_of_qdafn(				// indexing of QDAFN
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   L,							// number of projections
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int kfn_of_qdafn(					// c-k-AFN via QDAFN
	int   qn,							// number of query points
	int   d,							// dimensionality
	const float  **query,				// query set
	const Result **R,					// truth set
	const char   *data_folder,			// data folder
	const char   *output_folder);		// output folder

#endif // __AFN_H
