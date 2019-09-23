#ifndef __DRUSILLA_SELECT_H
#define __DRUSILLA_SELECT_H

// -----------------------------------------------------------------------------
//  Drusilla_Select: data structure of Drusilla_Select for c-k-AFN search
// -----------------------------------------------------------------------------
class Drusilla_Select {
public:
	Drusilla_Select();				// default constructor
	~Drusilla_Select();				// destructor

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// number of data objects
		int   d,						// number of dimensions
		int   l,						// number of projections
		int   m,						// number of candidates on each proj
		int   B,						// page size
		const float **data,				// data objects
		const char  *path);				// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *path);				// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	long long search(				// c-k-AFN search
		const float *query,				// query object
		const char  *data_folder,		// new format data folder
		MaxK_List   *list);				// top-k results (return)

protected:
	int  n_pts_;					// number of data objects
	int  dim_;						// dimensionality
	int  l_;						// number of random projections
	int  m_;						// number of candidates
	int  B_;						// page size

	char path_[200];				// address of index
	int  *cand_;					// candidates on each projection

	// -------------------------------------------------------------------------
	int bulkload(					// build hash tables
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk
};

#endif // __DRUSILLA_SELECT_H
