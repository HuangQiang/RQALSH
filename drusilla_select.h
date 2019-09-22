#ifndef __DRUSILLA_SELECT_H
#define __DRUSILLA_SELECT_H

// -----------------------------------------------------------------------------
//  Drusilla_Index: data structure of Drusilla_Select for c-k-AFN search
// -----------------------------------------------------------------------------
class Drusilla_Index {
public:
	Drusilla_Index();				// default constructor
	~Drusilla_Index();				// destructor

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// number of data objects
		int   d,						// number of dimensions
		int   l,						// number of projections
		int   m,						// number of candidates on each proj
		int   B,						// page size
		const float **data,				// data objects
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int search(						// c-k-AFN search
		const float *query,				// query object
		const char *data_folder,		// new format data folder
		MaxK_List* list);				// top-k results (return)

protected:
	int  n_pts_;					// number of data objects
	int  dim_;						// dimensionality
	int  l_;						// number of random projections
	int  m_;						// number of candidates
	int  B_;						// page size

	char fname_[200];				// address of index
	int  *fn_cand_;					// candidates on each projection

	// -------------------------------------------------------------------------
	int bulkload(					// build hash tables
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk
};

#endif // __DRUSILLA_SELECT_H
