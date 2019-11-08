#include <algorithm>
#include <cassert>
#include <cstring>

#include "def.h"
#include "util.h"
#include "block_file.h"
#include "b_node.h"
#include "b_tree.h"

// -----------------------------------------------------------------------------
B_Node::B_Node()						// constructor
{
	level_         = -1;
	num_entries_   = -1;
	left_sibling_  = -1;
	right_sibling_ = -1;
	key_           = NULL;
	son_           = NULL;

	dirty_         = false;
	block_         = -1;
	capacity_      = -1;
	btree_         = NULL;
}

// -----------------------------------------------------------------------------
B_Node::~B_Node()						// destructor
{
	// -------------------------------------------------------------------------
	//  if dirty, write back to disk
	// -------------------------------------------------------------------------
	if (dirty_) {
		int block_length = btree_->file_->get_blocklength();
		char *buf = new char[block_length];
		write_to_buffer(buf);
		btree_->file_->write_block(buf, block_);
		delete[] buf;buf = NULL;
	}

	if (key_ != NULL) {
		delete[] key_; key_ = NULL;
	}
	if (son_ != NULL) {
		delete[] son_; son_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void B_Node::init(					// init a new node, which not exist
	int   level,						// level (depth) in b-tree
	B_Tree *btree)						// b-tree of this node
{
	btree_         = btree;
	level_         = (char) level;
	num_entries_   = 0;
	left_sibling_  = -1;
	right_sibling_ = -1;
	dirty_         = true;

	int b_len = btree_->file_->get_blocklength();
	capacity_ = (b_len - get_header_size()) / get_entry_size();
	if (capacity_ < 50) {			// at least 50 entries
		printf("capacity = %d, which is too small\n", capacity_);
		exit(1);
	}

	key_ = new float[capacity_];
	son_ = new int[capacity_];
	for (int i = 0; i < capacity_; i++) {
		key_[i] = MINREAL;
		son_[i] = -1;
	}

	char* blk = new char[b_len];	// init <block_>, get new addr
	block_ = btree_->file_->append_block(blk);
	delete[] blk; blk = NULL;
}

// -----------------------------------------------------------------------------
void B_Node::init_restore(			// load an exist node from disk to init
	B_Tree *btree,						// b-tree of this node
	int   block)						// addr of disk for this node
{
	btree_ = btree;
	block_ = block;
	dirty_ = false;

	int b_len = btree_->file_->get_blocklength();
	capacity_ = (b_len - get_header_size()) / get_entry_size();
	if (capacity_ < 50) {			// at least 50 entries
		printf("capacity = %d, which is too small\n", capacity_);
		exit(1);
	}

	key_ = new float[capacity_];
	son_ = new int[capacity_];
	for (int i = 0; i < capacity_; i++) {
		key_[i] = MINREAL;
		son_[i] = -1;
	}

	// -------------------------------------------------------------------------
	//  read the buffer <blk> to init <level_>, <num_entries_>, 
	//  <left_sibling_>, <right_sibling_>, <key_> and <son_>.
	// -------------------------------------------------------------------------
	char* blk = new char[b_len];
	btree_->file_->read_block(blk, block);
	read_from_buffer(blk);

	delete[] blk; blk = NULL;
}

// -----------------------------------------------------------------------------
void B_Node::read_from_buffer(		// read a b-node from buffer
	const char *buf)					// store info of a b-index node
{
	int i = 0;
	// -------------------------------------------------------------------------
	//  read header: <level_> <num_entries_> <left_sibling_> <right_sibling_>
	// -------------------------------------------------------------------------
	memcpy(&level_, &buf[i], SIZECHAR);
	i += SIZECHAR;

	memcpy(&num_entries_, &buf[i], SIZEINT);
	i += SIZEINT;

	memcpy(&left_sibling_, &buf[i], SIZEINT);
	i += SIZEINT;

	memcpy(&right_sibling_, &buf[i], SIZEINT);
	i += SIZEINT;

	// -------------------------------------------------------------------------
	//  read entries: <key_> and <son_>
	// -------------------------------------------------------------------------
	for (int j = 0; j < num_entries_; j++) {
		memcpy(&key_[j], &buf[i], SIZEFLOAT);
		i += SIZEFLOAT;

		memcpy(&son_[j], &buf[i], SIZEINT);
		i += SIZEINT;
	}
}

// -----------------------------------------------------------------------------
void B_Node::write_to_buffer(		// write info of node into buffer
	char *buf)							// store info of this node (return)
{
	int i = 0;
	// -------------------------------------------------------------------------
	//  write header: <level_> <num_entries_> <left_sibling_> <right_sibling_>
	// -------------------------------------------------------------------------
	memcpy(&buf[i], &level_, SIZECHAR);
	i += SIZECHAR;

	memcpy(&buf[i], &num_entries_, SIZEINT);
	i += SIZEINT;

	memcpy(&buf[i], &left_sibling_, SIZEINT);
	i += SIZEINT;

	memcpy(&buf[i], &right_sibling_, SIZEINT);
	i += SIZEINT;

	// -------------------------------------------------------------------------
	//  write entries: <key_> and <son_>
	// -------------------------------------------------------------------------
	for (int j = 0; j < num_entries_; j++) {
		memcpy(&buf[i], &key_[j], SIZEFLOAT);
		i += SIZEFLOAT;

		memcpy(&buf[i], &son_[j], SIZEINT);
		i += SIZEINT;
	}
}

// -----------------------------------------------------------------------------
B_Node* B_Node::get_left_sibling()	// get the left-sibling node
{
	B_Node *node = NULL;
	if (left_sibling_ != -1) {		// left sibling node exist
		node = new B_Node();			// read left-sibling from disk
		node->init_restore(btree_, left_sibling_);
	}
	return node;
}

// -----------------------------------------------------------------------------
B_Node* B_Node::get_right_sibling()	// get the right-sibling node
{
	B_Node *node = NULL;
	if (right_sibling_ != -1) {		// right sibling node exist
		node = new B_Node();			// read right-sibling from disk
		node->init_restore(btree_, right_sibling_);
	}
	return node;
}

// -----------------------------------------------------------------------------
int B_Node::get_son(					// get son indexed by <index>
	int index)							// input index
{
	assert(index >= 0 && index < num_entries_);
	return son_[index];
}

// -----------------------------------------------------------------------------
void B_Node::add_new_child(			// add a new entry from its child node
	float key,							// input key
	int son)							// input son
{
	assert(num_entries_ < capacity_);

	key_[num_entries_] = key;		// add new entry into its pos
	son_[num_entries_] = son;

	num_entries_++;					// update <num_entries_>
	dirty_ = true;					// node modified, <dirty_> is true
}

// -----------------------------------------------------------------------------
//  Finds position of entry that is just less than or equal to input entry.
//  If input entry is smaller than all entry in this node, we'll return -1. 
//  The scan order is from right to left.
// -----------------------------------------------------------------------------
int B_Node::find_position_by_key(	// find pos just less than input entry
	float key)							// input key
{
	int pos = -1;
	for (int i = num_entries_ - 1; i >= 0; --i) {
		if (key_[i] <= key) {
			pos = i;
			break;
		}
	}
	return pos;
}

// -----------------------------------------------------------------------------
float B_Node::get_key(				// get <key> indexed by <index>
	int index)							// input index
{
	assert(index >= 0 && index < num_entries_);
	return key_[index];
}

// -----------------------------------------------------------------------------
//  entry: <key_> (SIZEFLOAT) and <son_> (SIZEINT)
// -----------------------------------------------------------------------------
int B_Node::get_entry_size()			// get entry size of b-node
{
	int entry_size = SIZEFLOAT + SIZEINT;
	return entry_size;
}

// -----------------------------------------------------------------------------
//  <level>: SIZECHAR
//  <num_entries> <left_sibling> and <right_sibling>: SIZEINT
// -----------------------------------------------------------------------------
int B_Node::get_header_size()		// get header size of b-node
{
	int header_size = SIZECHAR + SIZEINT * 3;
	return header_size;
}

// -----------------------------------------------------------------------------
int B_Node::get_block()				// get <block_> (address of this node)
{
	return block_;
}

// -----------------------------------------------------------------------------
int B_Node::get_num_entries()		// get <num_entries_>
{
	return num_entries_;
}

// -----------------------------------------------------------------------------
int B_Node::get_level()				// get <level_>
{
	return level_;
}

// -----------------------------------------------------------------------------
float B_Node::get_key_of_node()		// get key of this node
{
	if (key_ != NULL) return key_[0];
	else return -1.0f;
}

// -----------------------------------------------------------------------------
bool B_Node::isFull()				// whether is full?
{
	if (num_entries_ >= capacity_) return true;
	else return false;
}

// -----------------------------------------------------------------------------
void B_Node::set_left_sibling(		// set addr of left sibling node
	int left_sibling)					// addr of left sibling node
{
	left_sibling_ = left_sibling;
}

// -----------------------------------------------------------------------------
void B_Node::set_right_sibling(		// set addr of right sibling node
	int right_sibling)					// addr of right sibling node
{
	right_sibling_ = right_sibling;
}
