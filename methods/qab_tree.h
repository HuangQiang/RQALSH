#ifndef __QAB_TREE_H
#define __QAB_TREE_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "def.h"
#include "util.h"
#include "block_file.h"
#include "qab_node.h"

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
	void load_root(); 				// load root of b-tree

	// -------------------------------------------------------------------------
	void delete_root();				// delete root of b-tree
};

#endif // __QAB_TREE_H
