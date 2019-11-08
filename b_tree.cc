#include <algorithm>
#include <cassert>
#include <cstring>

#include "def.h"
#include "util.h"
#include "block_file.h"
#include "b_node.h"
#include "b_tree.h"

// -----------------------------------------------------------------------------
//  B_Tree: structure to index the projection produced by QDAFN
// -----------------------------------------------------------------------------
B_Tree::B_Tree()						// default constructor
{
	root_     = -1;
	file_     = NULL;
	root_ptr_ = NULL;
}

// -----------------------------------------------------------------------------
B_Tree::~B_Tree()						// destructor
{
	char *header = new char[file_->get_blocklength()];
	write_header(header);			// write <root_> to <header>
	file_->set_header(header);		// write back to disk
	delete[] header; header = NULL;

	if (root_ptr_ != NULL) {
		delete root_ptr_; root_ptr_ = NULL;
	}
	if (file_ != NULL) {
		delete file_; file_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void B_Tree::init(					// init a new tree
	int   b_length,						// block length
	const char *fname)					// file name
{
	FILE *fp = fopen(fname, "r");
	if (fp) {
		fclose(fp);					// ask whether replace?
		printf("The file \"%s\" exists. Replace? (y/n)", fname);

		char c = getchar();			// input 'Y' or 'y' or others
		getchar();					// input <ENTER>
		assert(c == 'y' || c == 'Y');
		remove(fname);				// otherwise, remove existing file
	}
	file_ = new BlockFile(b_length, fname); // b-tree store here

	// -------------------------------------------------------------------------
	//  Init the first node: to store <blocklength> (page size of a node),
	//  <number> (number of nodes including both index node and leaf node), and
	//  <root> (address of root node)
	// -------------------------------------------------------------------------
	//root_ptr_ = new B_Node();
	//root_ptr_->init(0, this);
	//root_ = root_ptr_->get_block();
	//delete_root();
}

// -----------------------------------------------------------------------------
void B_Tree::init_restore(			// load the tree from a tree file
	const char *fname)					// file name
{
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("index file %s does not exist\n", fname);
		exit(1);
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  It doesn't matter to initialize blocklength to 0.
	//  after reading file, <blocklength> will be reinitialized by file.
	// -------------------------------------------------------------------------
	file_ = new BlockFile(0, fname);
	root_ptr_ = NULL;
	
	// -------------------------------------------------------------------------
	//  Read the content after first 8 bytes of first block into <header>
	// -------------------------------------------------------------------------
	char *header = new char[file_->get_blocklength()];
	file_->read_header(header);		// read remain bytes from header
	read_header(header);			// init <root> from <header>

	delete[] header; header = NULL;
}

// -----------------------------------------------------------------------------
int B_Tree::read_header(				// read <root> from buffer
	const char *buf)					// buffer
{
	int i = 0;
	memcpy(&root_, buf, SIZEINT);
	i += SIZEINT;

	return i;
}

// -----------------------------------------------------------------------------
int B_Tree::write_header(			// write <root> to buffer
	char *buf)							// buffer (return)
{
	int i = 0;
	memcpy(buf, &root_, SIZEINT);
	i += SIZEINT;

	return i;
}

// -----------------------------------------------------------------------------
void B_Tree::load_root()				// load <root_ptr> of b-tree
{
	if (root_ptr_ == NULL)  {
		root_ptr_ = new B_Node();
		root_ptr_->init_restore(this, root_);
	}
}

// -----------------------------------------------------------------------------
void B_Tree::delete_root()			// delete <root_ptr>
{
	delete root_ptr_; root_ptr_ = NULL;
}

// -----------------------------------------------------------------------------
int B_Tree::bulkload(				// bulkload b-tree
	int   n,							// number of entries
	const Result *table)				// table of projected distance
{
	B_Node *child   = NULL;
	B_Node *prev_nd = NULL;
	B_Node *act_nd  = NULL;

	int   son = -1;					// address of node or object id
	float key = MINREAL;			// key of node

	bool first_node  = false;		// determine relationship of sibling
	int  start_block = -1;			// position of first node
	int  end_block   = -1;			// position of last node

	// -------------------------------------------------------------------------
	//  Stop consition: lastEndBlock == lastStartBlock (only one node, as root)
	// -------------------------------------------------------------------------
	std::vector<float> vec;			// store the key of node
	vec.push_back(0.0);				// address of b-node starts from 1
									// node 0 store the header info of b-tree

	int current_level    = 0;		// current level (leaf level is 0)
	int last_start_block = 0;		// to build b-tree level by level
	int last_end_block   = n - 1;	// to build b-tree level by level

	while (last_end_block > last_start_block) {
		first_node = true;
		for (int i = last_start_block; i <= last_end_block; i++) {
			if (current_level == 0) {
				son = table[i].id_;
				key = table[i].key_;
			}
			else {
				son = i;			// get address of son node
				key = vec[i];
			}

			if (!act_nd) {
				act_nd = new B_Node();
				act_nd->init(current_level, this);

				if (first_node) {
					first_node = false;
					start_block = act_nd->get_block();
				}
				else {
					act_nd->set_left_sibling(prev_nd->get_block());
					prev_nd->set_right_sibling(act_nd->get_block());
									// write the node to disk
					delete prev_nd; prev_nd = NULL;
				}
				end_block = act_nd->get_block();
			}

			act_nd->add_new_child(key, son);
			if (act_nd->isFull()) {	// get key of this node
				vec.push_back(act_nd->get_key_of_node());

				prev_nd = act_nd;	// full, change next node
				act_nd = NULL;
			}
		}
		if (prev_nd != NULL) {
			delete prev_nd; prev_nd = NULL;
		}
		if (act_nd != NULL) {		// get key of this node
			vec.push_back(act_nd->get_key_of_node());
			delete act_nd; act_nd = NULL;
		}
									// update info
		last_start_block = start_block;
		last_end_block = end_block;
		current_level++;
	}
	root_ = last_start_block;		// update the <root>

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (prev_nd != NULL) {
		delete prev_nd; prev_nd = NULL;
	}
	if (act_nd != NULL) {
		delete act_nd;  act_nd = NULL;
	}
	if (child != NULL) {
		delete child;   child = NULL;
	}

	return 0;
}
