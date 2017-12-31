#ifndef __RQALSH_H
#define __RQALSH_H

class BNode;
class BTree;
class MaxK_List;

// -----------------------------------------------------------------------------
//  PageBuffer: buffer of a page for the AFN search of rqalsh
// -----------------------------------------------------------------------------
struct PageBuffer {
	BLeafNode* leaf_node_;			// leaf node (level = 0)
	int index_pos_;					// cur pos of key in leaf node
	int leaf_pos_;					// cur pos of object id in leaf node
	int size_;						// size for one scan
};

// -----------------------------------------------------------------------------
//  HashValue: structure of hash table to store hash values
// -----------------------------------------------------------------------------
struct HashValue {
	int id_;						// object id
	float proj_;					// projection of the object
};

// -----------------------------------------------------------------------------
//  RQALSH: structure of rqalsh indexed by b+ tree. RQALSH is used to solve
//  the problem of c-Approximate Furthest Neighbor (c-AFN) search.
// -----------------------------------------------------------------------------
class RQALSH {
public:
	RQALSH();						// constructor
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init params of rqalsh
		int   n,						// number of data points
		int   d,						// dimension of space
		int   B,						// page size
		int   beta,						// false positive percentage
		float delta,					// error probability
		float ratio,					// approximation ratio
		char* output_folder);			// folder of info of rqalsh

	// -------------------------------------------------------------------------
	int restore(					// restore params of rqalsh
		char* output_folder);			// folder of info of rqalsh

	// -------------------------------------------------------------------------
	int bulkload(					// build b+ trees by bulkloading
		float** data);					// data set

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int top_k,						// top-k value
		float* query,					// one query point
		char*  data_folder,				// folder to store new format of data
		MaxK_List* list);				// c-k-AFN results (returned)
private:
	// -------------------------------------------------------------------------
	int   n_pts_;					// number of points
	int   dim_;						// dimensionality of space
	int   B_;						// page size in words
	float appr_ratio_;				// approximation ratio

	// -------------------------------------------------------------------------
	float w_;						// bucket width
	float p1_;						// positive probability
	float p2_;						// negative probability

	float alpha_;					// separation threshold percentage
	float beta_;					// false positive percentage
	float delta_;					// error probability

	int m_;							// number of hashtables
	int l_;							// separation threshold

	float* a_array_;				// hash functions
	char index_path_[200];			// folder path of index

	int dist_io_;					// io for calculating Euclidean distance
	int page_io_;					// io for scanning pages by rqalsh
	BTree** trees_;					// b-trees

	// -------------------------------------------------------------------------
	void calc_params();				// calc parama of rqalsh

	float calc_l2_prob(				// calc <p1> and <p2> for L2 distance
		float x);						// para

	void display_params();			// display params

	void gen_hash_func();			// generate hash functions

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int table_id,					// hash table id
		float* point);					// one point

	int write_para_file(			// write file of para
		char* fname);					// file name of para

	int read_para_file(				// read file of para
		char* fname);					// file name of para

	void get_tree_filename(			// get file name of tree
		int tree_id,					// tree id
		char* fname);					// file name of tree (return)

	// -------------------------------------------------------------------------
	void init_buffer(				// init page buffer (loc pos of b-treee)
		int num_tables,					// num of hash tables used for search
		PageBuffer* lptr,				// left buffer page (return)
		PageBuffer* rptr);				// right buffer page (return)

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius
		int num_tables,					// num of hash tables used for search
		PageBuffer* lptr,				// left page buffer
		PageBuffer* rptr,				// right page buffer
		float* q_val);					// hash value of query

	// -------------------------------------------------------------------------
	float update_radius(			// update radius
		int num_tables,					// num of hash tables used for search
		PageBuffer* lptr,				// left page buffer
		PageBuffer* rptr,				// right page buffer
		float* q_val);					// hash value of query

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		PageBuffer* lptr,				// left buffer
		const PageBuffer* rptr);		// right buffer

	void update_right_buffer(		// update right buffer
		const PageBuffer* lptr,			// left buffer
		PageBuffer* rptr);				// right buffer

	// -------------------------------------------------------------------------
	float calc_proj_dist(			// calc proj dist
		const PageBuffer* ptr,			// page buffer
		float q_val);					// hash value of query
};


// -----------------------------------------------------------------------------
//  Comparison function for qsort called by RQALSH::bulkload()
// -----------------------------------------------------------------------------
int HashValueQsortComp(				// compare function for qsort
	const void* e1,						// 1st element
	const void* e2);					// 2nd element

#endif
