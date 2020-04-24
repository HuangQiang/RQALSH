#ifndef __RQALSH_STAR_H
#define __RQALSH_STAR_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"
#include "rqalsh.h"

class RQALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  RQALSH* is used for c-k-Approximate Furthest Neighbor (c-k-AFN) search
// -----------------------------------------------------------------------------
class RQALSH_STAR {
public:
	RQALSH_STAR();					// default constructor
	~RQALSH_STAR();					// destructor

	// -------------------------------------------------------------------------
	int build(						// build index		
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   B,						// page size
		int   L,						// number of projection
		int   M,						// number of candidates
		int   beta,						// false positive percentage
		float delta,					// error probability
		float ratio,					// approximation ratio
		const float **data, 			// data objects
		const char  *path);				// index path

	// -------------------------------------------------------------------------
	int load(   					// load index
		const char *path);				// index path

	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	uint64_t kfn(					// c-k-AFN search
		int   top_k,					// top-k value
		const float *query,				// query objects
		const char  *data_folder,		// data folder
		MaxK_List   *list);				// k-FN results (return)

protected:
	int    n_pts_;					// number of data objects
	int    dim_;					// dimensionality
	int    B_;						// page size
	int    L_;						// number of projection
	int    M_;						// number of candidates
	char   path_[200];				// index path
	int    *cand_;				    // candidate id
	RQALSH *lsh_;					// index of sample data objects

	// -------------------------------------------------------------------------
	void data_dependent_select(		// data dependent selection
		const float **data,				// data objects
		int  *cand);					// candidate id (return)

	// -------------------------------------------------------------------------
	void calc_shift_data(			// calculate shift data objects
		const float **data,				// data objects
		int   &max_id,					// data id with max l2-norm (return)
		float &max_norm,				// max l2-norm (return)
		float *norm,					// l2-norm of shift data (return)
		float **shift_data); 			// shift data (return)
}; 

#endif // __RQALSH_STAR_H
