#ifndef __QAB_NODE_H
#define __QAB_NODE_H

class BlockFile;
class QAB_Tree;

// -----------------------------------------------------------------------------
//  QAB_Node: query-aware node in query-aware b-tree
// -----------------------------------------------------------------------------
class QAB_Node {
public:
	QAB_Node();						// constructor
	virtual ~QAB_Node();			// destructor

	// -------------------------------------------------------------------------
	virtual void init(				// init a new node, which not exist
		int   level,					// level (depth) in b-tree
		QAB_Tree *btree);					// b-tree of this node

	// -------------------------------------------------------------------------
	virtual void init_restore(		// load an exist node from disk to init
		QAB_Tree *btree,					// b-tree of this node
		int   block);					// address of file of this node

	// -------------------------------------------------------------------------
	virtual void read_from_buffer(	// read a b-node from buffer
		const char *buf);				// store info of a b-node

	virtual void write_to_buffer(	// write a b-node into buffer
		char *buf);						// store info of a b-node (return)

	// -------------------------------------------------------------------------
	virtual int get_entry_size();	// get entry size in b-node

	virtual int find_position_by_key(// find pos just less than input key
		float key);						// input key

	virtual float get_key(			// get <key> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	virtual QAB_Node* get_left_sibling(); // get left sibling node

	virtual QAB_Node* get_right_sibling(); // get right sibling node

	// -------------------------------------------------------------------------
	int get_block();				// get <block_>

	int get_num_entries();			// get <num_entries_>

	int get_level();				// get <level_>

	// -------------------------------------------------------------------------
	int get_header_size();			// get header size in b-node

	float get_key_of_node();		// get key of this node

	bool isFull();					// whether is full?

	// -------------------------------------------------------------------------
	void set_left_sibling(			// set <left_sibling>
		int left_sibling);				// addr of left sibling node

	void set_right_sibling(			// set <right sibling>
		int right_sibling);				// addr of right sibling node

protected:
	char  level_;					// level of b-tree (level > 0)
	int   num_entries_;				// number of entries in this node
	int   left_sibling_;			// addr in disk for left  sibling
	int   right_sibling_;			// addr in disk for right sibling
	float *key_;					// keys

	bool  dirty_;					// if dirty, write back to file
	int   block_;					// addr of disk for this node
	int   capacity_;				// max num of entries can be stored
	QAB_Tree *btree_;					// b-tree of this node
};

// -----------------------------------------------------------------------------
//  QAB_IndexNode: query-aware index node in query-aware b-tree
// -----------------------------------------------------------------------------
class QAB_IndexNode : public QAB_Node {
public:
	QAB_IndexNode();					// constructor
	virtual ~QAB_IndexNode();			// destructor

	// -------------------------------------------------------------------------
	virtual void init(				// init a new node, which not exist
		int   level,					// level (depth) in b-tree
		QAB_Tree *btree);					// b-tree of this node

	virtual void init_restore(		// load an exist node from disk to init
		QAB_Tree *btree,					// b-tree of this node
		int   block);					// address of file of this node

	// -------------------------------------------------------------------------
	virtual void read_from_buffer(	// read a b-node from buffer
		const char *buf);				// store info of a b-node

	virtual void write_to_buffer(	// write a b-node into buffer
		char *buf);						// store info of a b-node (return)

	// -------------------------------------------------------------------------
	virtual int get_entry_size();	// get entry size in b-node

	virtual int find_position_by_key(// find pos just less than input key
		float key);						// input key

	virtual float get_key(			// get <key_> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	virtual QAB_IndexNode* get_left_sibling(); // get left sibling node

	virtual QAB_IndexNode* get_right_sibling(); // get right sibling node

	// -------------------------------------------------------------------------
	int get_son(					// get <son_> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	void add_new_child(				// add new child by its child node
		float key,						// input key
		int son);						// input son

protected:
	int *son_;						// addr of son node
};


// -----------------------------------------------------------------------------
//  QAB_LeafNode: query-aware leaf node in query-aware b-tree
// -----------------------------------------------------------------------------
class QAB_LeafNode : public QAB_Node {
public:
	QAB_LeafNode();					// constructor
	virtual ~QAB_LeafNode();		// destructor

	// -------------------------------------------------------------------------
	virtual void init(				// init a new node, which not exist
		int   level,					// level (depth) in b-tree
		QAB_Tree *btree);				// b-tree of this node

	virtual void init_restore(		// load an exist node from disk to init
		QAB_Tree *btree,				// b-tree of this node
		int   block);					// address of file of this node

	// -------------------------------------------------------------------------
	virtual void read_from_buffer(	// read a b-node from buffer
		const char *buf);				// store info of a b-node

	virtual void write_to_buffer(	// write a b-node into buffer
		char *buf);						// store info of a b-node (return)

	// -------------------------------------------------------------------------
	virtual int get_entry_size();	// get entry size in b-node

	virtual int find_position_by_key( // find pos just less than input key
		float key);						// input key

	virtual float get_key(			// get <key_> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	virtual QAB_LeafNode* get_left_sibling(); // get left sibling node

	virtual QAB_LeafNode* get_right_sibling(); // get right sibling node

	// -------------------------------------------------------------------------
	int get_key_size(				// get key size of this node
		int block_length);				// block length

	int get_increment();			// get <increment>

	int get_num_keys();				// get <num_keys_>

	int get_entry_id(				// get entry id indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	void add_new_child(				// add new child by input id and key
		int id,							// input object id
		float key);						// input key

protected:
	int num_keys_;					// number of keys
	int *id_;						// object id

	int capacity_keys_;				// max num of keys can be stored
};

#endif // __QAB_NODE_H
