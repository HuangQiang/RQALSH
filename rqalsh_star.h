#ifndef __RQALSH_STAR_H
#define __RQALSH_STAR_H

class RQALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  RQALSH_Star is used to solve the problem of c-k-Approximate Furthest 
//  Neighbor (c-k-AFN) search
// -----------------------------------------------------------------------------
class RQALSH_Star {
public:
	RQALSH_Star();					// default constructor
	~RQALSH_Star();					// destructor

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
		const char *index_path,			// index path
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int load(   					// load index
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int top_k,						// top-k value
		const float *query,				// query objects
		const char *data_folder,		// data folder
		MaxK_List *list);				// k-FN results (return)

protected:
	int    n_pts_;					// number of data objects
	int    dim_;					// dimensionality
	int    B_;						// page size
	int    L_;						// number of projection
	int    M_;						// number of candidates
	int    beta_;                   // false positive percentage
    float  delta_;                  // error probability
    float  appr_ratio_;				// approximation ratio
	char   index_path_[200];		// index path

	int    sample_size_;			// number of sample data objects
	int    *sample_id_;			    // sample data objects id
	RQALSH *lsh_;					// index of sample data objects

    // -------------------------------------------------------------------------
	int bulkload(					// bulkloading
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int calc_shift_data(			// calc shift data
		const float **data,  			// data objects
		float *shift_data);  			// shift data objects (return)

	// -------------------------------------------------------------------------
	int data_dependent_select(		// data dependent selection
		const float *shift_data);		// shift data objects

	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk
}; 

#endif // __RQALSH_STAR_H