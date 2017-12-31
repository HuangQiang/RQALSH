#ifndef __PRI_QUEUE_H
#define __PRI_QUEUE_H


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
class MaxK_List {
public:
	MaxK_List(int max);				// constructor (given max size)
	~MaxK_List();					// destructor

	void reset();					// make exist list empty

	float min_key();				// return minimum key
	float max_key();				// return maximum key

	float ith_largest_key(int i);	// return ith key
	int ith_largest_id(int i);		// return ith id

	bool isFull();					// is full?

	void insert(					// insert item (inline for speed)
		float key,						// key of item
		int id);						// id of item

private:
	struct MaxK_Node {				// node in MaxK_List
		float key;						// key value
		int id;							// object id
	};

	int k_;							// max numner of keys
	int num_;						// number of key current active
	MaxK_Node* mk_;					// the list itself
};


// -----------------------------------------------------------------------------
//  Pri_Queue: this structure is a list of items, where the items are indexed
//  in a descending order by key.
//
//  The priority queue is maintained using a standard binary heap.
//  (Implementation Note: indexing is performed from [1..max] rather than the 
//  standard of [0..max-1], which simplifies parent/child computation.)
// -----------------------------------------------------------------------------
class Pri_Queue {
public:
	Pri_Queue(int max);				// constructor (given max size)
	~Pri_Queue();					// destructor

	bool empty();					// is empty?
	bool non_empty();				// is not empty?

	void reset();					// make exist queue empty

	void insert(					// insert item
		float key,						// key of item
		int id);						// id of item

	void extr_max(					// extract max item
		float& key,						// key of max item (returned)
		int& id);						// id of max item (returned)

private:
	struct PQ_Node {				// node of priority queue
		float key;						// key of node
		int id;							// id of node
	};

	int num_;						// number of items in queue
	int max_size_;					// maximum queue size
	PQ_Node* pq_;					// priority queue (array)
};

#endif
