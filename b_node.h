#ifndef __B_NODE_H
#define __B_NODE_H

class BlockFile;
class B_Tree;

// -----------------------------------------------------------------------------
//  B_Node: basic structure of node in b-tree
// -----------------------------------------------------------------------------
class B_Node {
public:
	B_Node();						// default constructor
	~B_Node();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init a new node, which not exist
		int   level,					// level (depth) in b-tree
		B_Tree *btree);					// b-tree of this node

	// -------------------------------------------------------------------------
	void init_restore(				// load an exist node from disk to init
		B_Tree *btree,					// b-tree of this node
		int   block);					// address of file of this node

	// -------------------------------------------------------------------------
	void read_from_buffer(			// read a b-node from buffer
		const char *buf);				// store info of a b-node

	// -------------------------------------------------------------------------
	void write_to_buffer(			// write a b-node into buffer
		char *buf);						// store info of a b-node (return)

	// -------------------------------------------------------------------------
	B_Node* get_left_sibling();		// get left sibling node

	// -------------------------------------------------------------------------
	B_Node* get_right_sibling();		// get right sibling node

	// -------------------------------------------------------------------------
	int get_son(					// get <son_> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	void add_new_child(				// add new child by its child node
		float key,						// input key
		int   son);						// input son

	// -------------------------------------------------------------------------
	int find_position_by_key(		// find pos just less than input key
		float key);						// input key

	// -------------------------------------------------------------------------
	float get_key(					// get <key> indexed by <index>
		int index);						// index

	// -------------------------------------------------------------------------
	int get_entry_size();			// get entry  size in b-node
	int get_header_size();			// get header size in b-node

	// -------------------------------------------------------------------------
	int get_block();				// get <block_>
	int get_num_entries();			// get <num_entries_>
	int get_level();				// get <level_>

	// -------------------------------------------------------------------------
	float get_key_of_node();		// get key of this node

	// -------------------------------------------------------------------------
	bool isFull();					// whether is full?

	// -------------------------------------------------------------------------
	void set_left_sibling(			// set <left_sibling>
		int left_sibling);				// left sibling node address

	// -------------------------------------------------------------------------
	void set_right_sibling(			// set <right sibling>
		int right_sibling);				// right sibling node address

protected:
	char  level_;					// level of b-tree (level > 0)
	int   num_entries_;				// number of entries in this node
	int   left_sibling_;			// addr in disk for left  sibling
	int   right_sibling_;			// addr in disk for right sibling
	float *key_;					// keys
	int   *son_;					// son node address or object id

	bool  dirty_;					// if dirty, write back to file
	int   block_;					// disk address for this node
	int   capacity_;				// max num of entries can be stored
	B_Tree *btree_;					// b-tree of this node
};

#endif // __B_NODE_H
