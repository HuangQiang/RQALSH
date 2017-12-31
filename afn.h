#ifndef __APPROXIMATE_FURTHEST_NEIGHBOR_H
#define __APPROXIMATE_FURTHEST_NEIGHBOR_H

// -----------------------------------------------------------------------------
int ground_truth(					// output ground truth
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	char* data_set,						// address of data set
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* output_folder);				// folder to store info of rqalsh
	
// -----------------------------------------------------------------------------
int indexing(						// build hash tables for the dataset
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   beta,
	float delta,
	float ratio,						// approximation ratio
	char* data_set,						// address of data set
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// folder to store info of rqalsh

// -----------------------------------------------------------------------------
int rqalsh_afn(						// c-k-AFN search
	int   qn,							// number of query points
	int   d,							// dimensionality
	char* query_set,					// path of query set
	char* truth_set,					// groundtrue file
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// output folder

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   B,							// page size
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// output folder

#endif
