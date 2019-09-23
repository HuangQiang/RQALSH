#ifndef __QAB_TREE_H
#define __QAB_TREE_H

class  BlockFile;
class  QAB_Node;
struct Result;

// -----------------------------------------------------------------------------
//  QAB_Tree: query-aware b-tree to index hash tables produced by rqalsh
// -----------------------------------------------------------------------------
class QAB_Tree {
public:
	int root_;						// address of disk for root
	QAB_Node *root_ptr_;			// pointer of root
	BlockFile *file_;				// file in disk to store
	
	// -------------------------------------------------------------------------
	QAB_Tree();						// constructor
	~QAB_Tree();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init a new b-tree
		int   b_length,					// block length
		const char *fname);				// file name	

	// -------------------------------------------------------------------------
	void init_restore(				// load an exist b-tree
		const char *fname);				// file name

	// -------------------------------------------------------------------------
	int bulkload(					// bulkload b-tree from hash table in mem
		int   n,						// number of entries
		const Result *table);		// hash table

protected:
	// -------------------------------------------------------------------------
	inline int read_header(const char *buf) { // read <root> from buffer
		memcpy(&root_, buf, SIZEINT);
		return SIZEINT;
	}

	// -------------------------------------------------------------------------
	inline int write_header(char *buf) { // write <root> into buffer
		memcpy(buf, &root_, SIZEINT);
		return SIZEINT;
	}

	// -------------------------------------------------------------------------
	inline void load_root() {		// load root of b-tree
		if (root_ptr_ == NULL) {
			root_ptr_ = new QAB_IndexNode();
			root_ptr_->init_restore(this, root_);
		}
	}

	// -------------------------------------------------------------------------
	inline void delete_root() {		// delete root of b-tree
		if (root_ptr_ != NULL) {
			delete root_ptr_; root_ptr_ = NULL;
		}
	}
};

#endif // __QAB_TREE_H
