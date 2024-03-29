#ifndef __PRI_QUEUE_H
#define __PRI_QUEUE_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "def.h"

// -----------------------------------------------------------------------------
//  basic data structures
// -----------------------------------------------------------------------------
struct Result {						// basic data structure 
	float key_;
	int   id_;
};

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
//  MaxK_List: the structure is one which maintains the largest k values (of 
//  type float) and associated object id (of type int).
// -----------------------------------------------------------------------------
class MaxK_List {
public:
	MaxK_List(int max);				// constructor (given max size)
	~MaxK_List();					// destructor

	// -------------------------------------------------------------------------
	inline void reset() { num_ = 0; }

	// -------------------------------------------------------------------------
	inline float max_key()	{ return num_ > 0 ? list_[0].key_ : MINREAL; }

	// -------------------------------------------------------------------------
	inline float min_key() { return num_ == k_ ? list_[k_-1].key_ : MINREAL; }

	// -------------------------------------------------------------------------
	inline float ith_key(int i) { return i < num_ ? list_[i].key_ : MINREAL; }

	// -------------------------------------------------------------------------
	inline int ith_id(int i) { return i < num_ ? list_[i].id_ : MININT; }

	// -------------------------------------------------------------------------
	inline int size() { return num_; }

	// -------------------------------------------------------------------------
	inline bool isFull() { if (num_ >= k_) return true; else return false; }

	// -------------------------------------------------------------------------
	float insert(					// insert item
		float key,						// key of item
		int id);						// id of item

private:
	int k_;							// max numner of keys
	int num_;						// number of key current active
	Result *list_;					// the list itself
};

#endif // __PRI_QUEUE_H
