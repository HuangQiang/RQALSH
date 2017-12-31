#ifndef __UTIL_H
#define __UTIL_H

class MaxK_List;

// -----------------------------------------------------------------------------
//  Global variables
// -----------------------------------------------------------------------------
extern long g_memory;

// -----------------------------------------------------------------------------
//  Uitlity functions
// -----------------------------------------------------------------------------
int compfloats(						// compare two float values
	float v1,							// 1st float value
	float v2);							// 2nd float value

// -----------------------------------------------------------------------------
void error(							// an error message
	char* msg,							// an message
	bool is_exit);						// whether exit the program

// -----------------------------------------------------------------------------
int check_mem();					// check memory is enough

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char* data_path);					// input path

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	float* vec1,						// 1st point
	float* vec2,						// 2nd point
	int    dim);						// dimension

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall between two ID list
	int n,								// length of ID list
	MaxK_List* list,					// ID list of k-nn method
	int* ID);							// ground truth ID list

// -----------------------------------------------------------------------------
void calc_relative_constrast(		// calc relative constrast
	int n,								// number of data objects
	int qn,								// number of queries
	int d,								// dimensionality
	FILE* fp,							// file pointer
	float** data,						// data objects
	float** query);						// queries


// -----------------------------------------------------------------------------
//  Functions used for the input/output of datasets and query sets.
// -----------------------------------------------------------------------------
int read_set(						// read (data or query) set from disk
	int n,								// number of data points
	int d,								// dimensionality
	char* set,							// address of dataset
	float** points);					// data or queries (return)

// -----------------------------------------------------------------------------
int write_data_new_form(			// write dataset with new format
	int n,								// cardinality
	int d,								// dimensionality
	int B,								// page size
	float** data,						// data set
	char* output_path);					// output path

// -----------------------------------------------------------------------------
void get_data_filename(				// get file name of data
	int data_id,						// data file id
	char* data_path,					// path to store data in new format
	char* fname);						// file name of data (return)

// -----------------------------------------------------------------------------
void write_data_to_buffer(			// write data to buffer
	int d,								// dimensionality
	int left,							// left data id
	int right,							// right data id
	float** data,						// data set
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write data to one page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
int read_data(						// read data from page
	int id,								// index of data
	int d,								// dimensionality
	int B,								// page size
	float* data,						// real data (return)
	char* output_path);					// output path

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read data from page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
void read_data_from_buffer(			// read data from buffer
	int index,							// index of data in buffer
	int d,								// dimensionality
	float* data,						// data set
	char* buffer);						// buffer to store data


#endif
