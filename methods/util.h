#ifndef __UTIL_H
#define __UTIL_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "def.h"
#include "pri_queue.h"

struct Result;
class  MaxK_List;

extern timeval  g_start_time;		// global parameter: start time
extern timeval  g_end_time;			// global parameter: end time

extern float    g_runtime;			// global parameter: running time
extern float    g_ratio;			// global parameter: overall ratio
extern float    g_recall;			// global parameter: recall
extern uint64_t g_io;				// global parameter: i/o cost
extern uint64_t g_memory;			// global parameter: memory usage

// -------------------------------------------------------------------------
void create_dir(					// create directory
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_txt_data(					// read data (text) from disk
	int   n,							// number of data/query objects
	int   d,							// dimensionality
	const char *fname,					// address of data/query set
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data points
	int   d,							// dimensionality
	const char *fname,					// address of data
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int write_data_new_form(			// write dataset with new format
	int   n,							// cardinality
	int   d,							// dimensionality
	int   B,							// page size
	const float **data,					// data set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
void get_data_filename(				// get file name of data
	int   data_id,						// data file id
	const char *data_path,				// path to store data in new format
	char  *fname);						// file name of data (return)

// -----------------------------------------------------------------------------
void write_data_to_buffer(			// write data to buffer
	int   d,							// dimensionality
	int   left,							// left  data id
	int   right,						// right data id
	const float **data,					// data set
	char  *buffer);						// buffer to store data (return)

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write data to one page
	int   B,							// page size
	const char *fname,					// file name of data
	const char *buffer);				// buffer to store data

// -----------------------------------------------------------------------------
int read_data_new_format(			// read data with new format from disk
	int   id,							// index of data
	int   d,							// dimensionality
	int   B,							// page size
	const char *output_path,			// output path
	float *data);						// real data (return)

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read data from page
	int   B,							// page size
	const char *fname,					// file name of data
	char  *buffer);						// buffer to store data (return)

// -----------------------------------------------------------------------------
void read_data_from_buffer(			// read data from buffer
	int   index,						// index of data in buffer
	int   d,							// dimensionality
	const char *buffer,					// buffer to store data
	float *data);						// data object (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	const Result *result);				// results returned by algorithms

// -----------------------------------------------------------------------------
uint64_t linear(					// linear scan search
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   top_k,						// top-k value
	const float *query,					// query object
	const char  *data_folder,			// data folder
	MaxK_List   *list);					// k-FN results (return)

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);			// address of truth set

#endif
