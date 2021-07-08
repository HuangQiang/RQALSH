#ifndef __RQALSH_H
#define __RQALSH_H

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
#include "qab_node.h"
#include "qab_tree.h"

class QAB_Node;
class QAB_LeafNode;
class QAB_Tree;
class MaxK_List;

// -----------------------------------------------------------------------------
//  Page: a buffer of one page for c-k-AFN search
// -----------------------------------------------------------------------------
struct Page {
	QAB_LeafNode *leaf_node_;		// leaf node (level = 0)
	int index_pos_;					// cur pos of key in leaf node
	int leaf_pos_;					// cur pos of object id in leaf node
	int size_;						// size for one scan
};

// -----------------------------------------------------------------------------
//  RQALSH: structure of RQALSH indexed by Query-Aware B+tree (QAB+Tree), which
//  is used for c-Approximate Furthest Neighbor (c-AFN) search.
// -----------------------------------------------------------------------------
class RQALSH {
public:
	RQALSH();						// default constructor
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// number of data objects
		int   d,						// dimension of space
		int   B,						// page size
		int   beta,						// false positive percentage
		float delta,					// error probability
		float ratio,					// approximation ratio
		const float **data,				// data objects
		const char  *path);				// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *path);				// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters
	
	// -------------------------------------------------------------------------
	uint64_t kfn(					// c-k-AFN search
		int   top_k,					// top-k value
		const float *query,				// query object
		const int   *index,				// mapping index for data objects
		const char  *data_folder,		// data folder
		MaxK_List   *list);				// k-FN results (return)

protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	int   B_;						// page size
	float beta_;					// false positive percentage
	float delta_;					// error probability
	float ratio_;					// approximation ratio
	float w_;						// bucket width
	int   m_;						// number of hashtables
	int   l_;						// collision threshold
	char  path_[200];				// index path

	float **a_;						// hash functions
	QAB_Tree **trees_;				// query-aware b+ trees
	uint64_t dist_io_;				// io for computing distance
	uint64_t page_io_;				// io for scanning pages

	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L2 distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	int bulkload(					// build QAB+Trees by bulkloading
		const float** data);			// data set

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// hash table id
		const float *data);				// one data object

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk

	// -------------------------------------------------------------------------
	void get_tree_filename(			// get file name of QAB+Tree
		int  tid,						// tree id
		char *fname);					// file name (return)

	// -------------------------------------------------------------------------
	void init_search_params(		// init parameters
		const float *query,				// query object
		float *q_val,					// hash values of query (return)
		Page  **lptrs,					// left buffer (return)
		Page  **rptrs);					// right buffer (return)

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius
		const float *q_val,				// hash value of query
		const Page **lptrs,				// left buffer
		const Page **rptrs);			// right buffer

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		const Page *rptr,				// right buffer
		Page *lptr);					// left buffer (return)

	// -------------------------------------------------------------------------
	void update_right_buffer(		// update right buffer
		const Page *lptr,				// left buffer
		Page *rptr);					// right buffer (return)

	// -------------------------------------------------------------------------
	float calc_dist(				// calc projected distance
		float q_val,					// hash value of query
		const Page *ptr);				// page buffer
	
	// -------------------------------------------------------------------------
	void delete_tree_ptr(			// delete the pointers of QAB+Trees
		Page **lptrs,					// left buffer (return)
		Page **rptrs);					// right buffer (return)
};

#endif // __RQALSH_H
