#ifndef __AFN_H
#define __AFN_H

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set);				// address of truth set

// -----------------------------------------------------------------------------
int indexing_of_rqalsh_star(		// indexing of RQALSH_Star
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int kfn_of_rqalsh_star(				// c-k-AFN search of RQALSH_Star
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection
	int   M,							// number of candidates
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int indexing_of_rqalsh(				// indexing of RQALSH
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int kfn_of_rqalsh(					// c-k-AFN search of RQALSH
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

#endif // __AFN_H
