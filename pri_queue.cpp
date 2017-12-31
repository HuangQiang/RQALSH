#include "headers.h"


// -----------------------------------------------------------------------------
//  MaxK_List: the structure is one which maintains the largest k values (of 
//  type float) and associated object id (of type int).
//
//  It is currently implemented using an array with k items. Items are stored 
//  in descending order by key, and its insertion is made through standard 
//  insertion mode. (It is quite inefficient if k is large, but current 
//  applications are able to handle since the value of k is small.)
//
//  Note that the priority queue contains k + 1 entries, while the last entry 
//  is used as a simple place holder and is otherwise ignored.
// -----------------------------------------------------------------------------
MaxK_List::MaxK_List(				// constructor (given max size)
	int max)							// max size
{
	num_ = 0;
	k_ = max;
	mk_ = new MaxK_Node[max + 1];
}

// -----------------------------------------------------------------------------
MaxK_List::~MaxK_List() 			// destructor
{
	if (mk_ != NULL) {
		delete[] mk_; mk_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void MaxK_List::reset()				// make exist queue empty
{
	num_ = 0;						// setting active num = 0
}

// -----------------------------------------------------------------------------
float MaxK_List::max_key()			// get maximum key
{
	return (num_ > 0 ? mk_[0].key : MIN_FLT);
}

// -----------------------------------------------------------------------------
float MaxK_List::min_key()			// get minimum key
{
	return (num_ == k_ ? mk_[k_-1].key : MIN_FLT);
}

// -----------------------------------------------------------------------------
float MaxK_List::ith_largest_key(	// return ith key
	int i)								// ith position
{
	return (i < num_ ? mk_[i].key : MIN_FLT);
}

// -----------------------------------------------------------------------------
int MaxK_List::ith_largest_id(		// return ith info
	int i)								// ith position
{
	return (i < num_ ? mk_[i].id : MIN_INT);
}

// -----------------------------------------------------------------------------
bool MaxK_List::isFull()			// is full?
{
	if (num_ >= k_) return true;
	else return false;
}

// -----------------------------------------------------------------------------
void MaxK_List::insert(				// insert item (inline for speed)
	float key,							// key of item
	int id)								// id of item
{
	int i = 0;
	for (i = num_; i > 0; i--) {
		if (mk_[i-1].key < key) mk_[i] = mk_[i - 1];
		else break;
	}
	mk_[i].key = key;				// store new item here
	mk_[i].id = id;
	if (num_ < k_) num_++;			// increase the number of items
}


// -----------------------------------------------------------------------------
//  Pri_Queue: this structure is a list of items, where the items are indexed
//  in a descending order by key.
//
//  The priority queue is maintained using a standard binary heap.
//  (Implementation Note: indexing is performed from [1..max] rather than the 
//  standard of [0..max-1], which simplifies parent/child computation.)
// -----------------------------------------------------------------------------
Pri_Queue::Pri_Queue(				// constructor (given max size)
	int max)							// max size
{
	num_ = 0;
	max_size_ = max;
	pq_ = new PQ_Node[max + 1];
}

// -----------------------------------------------------------------------------
Pri_Queue::~Pri_Queue()				// destructor
{
	if (pq_ != NULL) {
		delete[] pq_; pq_ = NULL;
	}
}

// -----------------------------------------------------------------------------
bool Pri_Queue::empty()				// is empty?
{
	if (num_ == 0) return true;
	else return false;
}

// -----------------------------------------------------------------------------
bool Pri_Queue::non_empty()			// is not empty?
{
	if (num_ == 0) return false;
	else return true;
}

// -----------------------------------------------------------------------------
void Pri_Queue::reset()				// make exist queue empty
{
	num_ = 0;						// setting current active num = 0
}

// -----------------------------------------------------------------------------
void Pri_Queue::insert(				// insert item (inline for speed)
	float key,							// key of item
	int id)								// id of item
{
	if (++num_ > max_size_) {
		error("Priority Queue Overflow!", false);
	}

	int r = num_;					// new position
	while (r > 1) {					// find proper place
		int p = r / 2;
		if (pq_[p].key >= key) break;
		else {						// swap with parent
			pq_[r] = pq_[p];
			r = p;
		}
	}
	pq_[r].key = key;				// insert new item
	pq_[r].id = id;
}

// -----------------------------------------------------------------------------
void Pri_Queue::extr_max(			// extract max item (then delete it)
	float& key,						// key of max item (returned)
	int& id)							// info of max item (returned)
{
	key = pq_[1].key;				// extract max item
	id  = pq_[1].id;
									// delete max item from pri queue
	float kn = pq_[num_--].key;		// last item in queue
	int p = 1;						// p point to item out of position
	int r = p << 1;					// left child of p

	while (r <= num_) {				// set r to larger child of p
		if (r < num_ && pq_[r].key < pq_[r + 1].key) r++;
		if (kn >= pq_[r].key) {		// in proper order
			break;
		} else {					// swap with child
			pq_[p] = pq_[r];
			p = r;
			r = p << 1;
		}
	}
	pq_[p].key = pq_[num_ + 1].key;	// insert last item
	pq_[p].id  = pq_[num_ + 1].id;
}
