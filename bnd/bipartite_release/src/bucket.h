#ifndef __BUCKET__H
#define __BUCKET__H

#include "larray.h"
#include <vector>
#include <unordered_map>


using namespace std;


struct Bucket_element
{
    int id;
    long long value;
    Bucket_element * next;
    Bucket_element * prev;
    Bucket_element(int id, long long value);
};

struct Bucket
{
    long long max_value;
    long long current_min_value;
//    std::vector<Bucket_element*> buckets;
    unordered_map<long long, Bucket_element*> buckets;

    HashMap<Bucket_element*> elements;
    Bucket();
    ~Bucket();
    void Initialize(long long max_value);
    void Insert(int id, long long value);
    void Update(int id, long long new_value);
    void DecVal(int id);
    void DecTo (int id, long long val);
    void IncVal(int id);
    // returns -1 if empty
    int PopMin(int* ret_id, long long* ret_value);
    // grab the current value of id
    long long CurrentValue(int id);
    void Free();
};
















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
	int nb_elements;
	int max_value;
	int *values; /* needed for update, in case bucket head changed. */
	int current_min_value;

	Naive_Bucket();
	~Naive_Bucket();
	void Initialize(int max_value, int nb_element);
	void Free ();
	/* value == INT_MAX means not present in bucket */
	void Insert(int id, int value);
	void Update(int id, int new_value);
	void DecTo(int id, int val);
	void DecVal(int id);
	void Adjust(int id);
	void Sulk (int id);
	/*returns -1 if empty*/
	int PopMin(int* ret_id, int* ret_value);
	int CurrentValue(int id);
};






struct Max_Bucket_element
{
	Max_Bucket_element * prev;
	Max_Bucket_element * next;
	Max_Bucket_element();
};

struct Max_Bucket
{
	Max_Bucket_element **buckets; /* actual pointers to bucket heads */
	Max_Bucket_element *elements; /* for direct access to bucket elements elements[id] is the id-th element */
	int nb_elements;
	int max_value;
	int *values; /* needed for update, in case bucket head changed. */
	int current_max_value;

	Max_Bucket();
	~Max_Bucket();
	void Initialize(int max_value, int nb_element);
	void Free ();
	/* value == INT_MAX means not present in bucket */
	void Insert(int id, int value);
	void Update(int id, int new_value);
//	void DecVal(int id);
	/*returns -1 if empty*/
	int PopMax(int* ret_id, int* ret_value);
	int CurrentValue(int id);
};



#endif
