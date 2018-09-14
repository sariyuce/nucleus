#ifndef __BUCKET__H
#define __BUCKET__H

#include <vector>
#include <unordered_map>
typedef long long lol;

typedef int item;
//typedef lol item; // for big graph runs

struct Naive_Bucket_element
{
	Naive_Bucket_element * prev;
	Naive_Bucket_element * next;
	Naive_Bucket_element();
};

struct Naive_Bucket
{
	Naive_Bucket_element **buckets; /* actual pointers to bucket heads */
	Naive_Bucket_element *elements; /* for direct access to bucket elements elements[id] is the id-th element */
	item nb_elements;
	item max_value;
	item *values; /* needed for update, in case bucket head changed. */
	item current_min_value;

	Naive_Bucket();
	~Naive_Bucket();
	void Initialize(item max_value, item nb_element);
	void Free (); /* value == INT_MAX means not present in bucket */
	void Insert(item id, item value);
	void Update(item id, item new_value);
	void DecVal(item id);
	item PopMin(item* ret_id, item* ret_value); /* returns -1 if empty */
	item CurrentValue(item id);
};

#endif
