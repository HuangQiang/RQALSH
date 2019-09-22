#ifndef __RQALSH_H
#define __RQALSH_H

class QAB_Node;
class QAB_Tree;
class MaxK_List;

// -----------------------------------------------------------------------------
//  PageBuffer: a buffer of one page for c-k-AFN search
// -----------------------------------------------------------------------------
struct PageBuffer {
	QAB_LeafNode *leaf_node_;		// leaf node (level = 0)
	int index_pos_;					// cur pos of key in leaf node
	int leaf_pos_;					// cur pos of object id in leaf node
	int size_;						// size for one scan
};

// -----------------------------------------------------------------------------
//  RQALSH is used to solve the problem of c-k-Approximate Furthest Neighbor 
//  (c-k-AFN) search
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
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	long long kfn(					// c-k-AFN search
		int top_k,						// top-k value
		const float *query,				// query object
		const char *data_folder,		// data folder
		MaxK_List *list);				// k-FN results (return)
	
	// -------------------------------------------------------------------------
	long long kfn(					// c-k-AFN search
		int top_k,						// top-k value
		const float *query,				// query object
		const int *object_id,			// object id mapping
		const char *data_folder,		// data folder
		MaxK_List *list);				// k-FN results (return)

protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	int   B_;						// page size
	float beta_;					// false positive percentage
	float delta_;					// error probability
	float appr_ratio_;				// approximation ratio
	char  index_path_[300];			// folder path of index

	float w_;						// bucket width
	float p1_;						// positive probability
	float p2_;						// negative probability
	float alpha_;					// collision threshold percentage
	int   m_;						// number of hashtables
	int   l_;						// collision threshold
	float *a_array_;				// hash functions
	QAB_Tree **trees_;					// b-trees

	int   dist_io_;					// io for computing distance
	int   page_io_;					// io for scanning pages
	int   *freq_;					// frequency of data objects
	bool  *checked_;				// whether the data objects are checked
	bool  *flag_;					// flag of bucket width
	float *data_;					// one data object
	float *q_val_;					// hash value of query
	
	PageBuffer **lptr_;				// left  pointer of B+ Tree
	PageBuffer **rptr_;				// right pointer of B+ Tree

	// -------------------------------------------------------------------------
	void calc_params();				// calc parameters

	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L2 distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	void gen_hash_func();			// generate hash functions

	// -------------------------------------------------------------------------
	int bulkload(					// build B+ Trees by bulkloading
		const float** data);			// data set

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int tid,						// hash table id
		const float *point);			// one data object

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk

	// -------------------------------------------------------------------------
	void get_tree_filename(			// get file name of tree
		int  tree_id,					// tree id
		char *fname);					// file name of tree (return)

	// -------------------------------------------------------------------------
	void init_search_params(		// init parameters for k-NN search
		const float *query);			// query object

	// -------------------------------------------------------------------------
	float find_radius();			// find proper radius

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		const PageBuffer *rptr,			// right buffer
		PageBuffer *lptr);				// left buffer (return)

	// -------------------------------------------------------------------------
	void update_right_buffer(		// update right buffer
		const PageBuffer *lptr,			// left buffer
		PageBuffer* rptr);				// right buffer (return)

	// -------------------------------------------------------------------------
	float calc_dist(				// calc projected distance
		float q_val,					// hash value of query
		const PageBuffer *ptr);			// page buffer
	
	// -------------------------------------------------------------------------
	void delete_tree_ptr();			// delete the pointers of B+ Trees
};

#endif // __RQALSH_H
