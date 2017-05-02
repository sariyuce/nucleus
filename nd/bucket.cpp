#include "bucket.h"

/******* ARRAY BUCKET ********/

Naive_Bucket_element::Naive_Bucket_element()
: next(NULL), prev(NULL)
{}

Naive_Bucket::Naive_Bucket()
: max_value(0), current_min_value(1), elements(NULL), buckets(NULL), values(NULL), nb_elements(0)
{}

Naive_Bucket::~Naive_Bucket() {
	if (buckets != NULL)
		free (buckets);
	if (elements != NULL)
		free (elements);
	if (values != NULL)
		free (values);
}
void Naive_Bucket::Free () {
	free (buckets);
	buckets = NULL;
	free (elements);
	elements = NULL;
	free (values);
	values = NULL;
}
void Naive_Bucket::Initialize(int max_v, int nb_element) {
	int i;
	max_value = max_v;
	buckets = (Naive_Bucket_element **) malloc(sizeof(Naive_Bucket_element *) * (max_value+1));
	elements = (Naive_Bucket_element *) malloc(sizeof(Naive_Bucket_element) * nb_element);
	//    printf("nb_element: %d\n", nb_element);
	values = (int *) malloc(sizeof(int) * nb_element);
	nb_elements = nb_element;
	if (buckets == NULL || elements == NULL || values == NULL) {
		free(values);
		free(buckets);
		free(elements);
	}
	else {
		for (i = 0; i <= max_value; i++)
			buckets[i] = NULL;
		for (i = 0; i < nb_element; i++) {
			elements[i].prev = NULL;
			elements[i].next = NULL;
		}
	}
	current_min_value = max_value + 1;
}

void Naive_Bucket::Insert (int id, int value) {
	values[id] = value;
	elements[id].prev = NULL;
	elements[id].next = buckets[value];
	if (buckets[value] != NULL)
		buckets[value]->prev = &(elements[id]);
	else if (current_min_value > value)
		current_min_value = value;
	buckets[value] = &(elements[id]);
}

int Naive_Bucket::PopMin(int* id, int* ret_value) {
	for (; current_min_value <= max_value; current_min_value++) {
		if (buckets[current_min_value] != NULL) {
			*id = buckets[current_min_value] - elements; // pointer arithmetic. finds the index of element that buckets[current_min_value] points to
			buckets[current_min_value] = buckets[current_min_value]->next; // adjust the pointer to the new head of the list that has same degree elements
			if (buckets[current_min_value] != NULL)
				buckets[current_min_value]->prev = NULL;
			*ret_value = values[*id];
			values[*id] = -1;
			return 0;
		}
	}
	return -1; // if the bucket is empty
}

int Naive_Bucket::CurrentValue(int id) {
	return values[id];
}

void Naive_Bucket::DecVal(int id) {
//	printf ("%d is decremented to %d\n", id, values[id]-1);
	int old_value = values[id];
	// adjust the prev and next pointers
	if (elements[id].prev == NULL)
		buckets[old_value] = elements[id].next;
	else
		elements[id].prev->next = elements[id].next;
	if (elements[id].next != NULL)
		elements[id].next->prev = elements[id].prev;
	Naive_Bucket::Insert(id, values[id]-1);
	return;
}

// bring the id to head of list
void Naive_Bucket::Adjust(int id) {
//	printf ("%d is adjusted\n", id);
	int value = values[id];
	// adjust the prev and next pointers
	if (elements[id].prev == NULL) // at the head of list
		return;
	else {
		if (elements[id].next != NULL)
			elements[id].next->prev = elements[id].prev;
		elements[id].prev->next = elements[id].next;
		elements[id].next = buckets[value];
		buckets[value]->prev = &(elements[id]);
		elements[id].prev = NULL;
		buckets[value] = &(elements[id]);
	}
	return;
}

// bring the id to tail of list, not constant time
void Naive_Bucket::Sulk (int id) {
//	printf ("%d is sulked\n", id);
	int value = values[id];
	// adjust the prev and next pointers
	if (elements[id].next == NULL) // at the tail of list
		return;
	// find the tail
	Naive_Bucket_element* tail = &(elements[id]);
	while (tail->next != NULL)
		tail = tail->next;
	if (elements[id].prev == NULL) {// at the head of list
		buckets[value] = elements[id].next;
		buckets[value]->prev = NULL;
	}
	else {// in the middle of list
		elements[id].prev->next = elements[id].next;
		elements[id].next->prev = elements[id].prev;
	}

	elements[id].prev = tail;
	elements[id].prev->next = &(elements[id]);
	elements[id].next = NULL;

	return;
}



// MAX BUCKET

Max_Bucket_element::Max_Bucket_element()
: next(NULL), prev(NULL)
{}

Max_Bucket::Max_Bucket()
: max_value(0), current_max_value(1), elements(NULL), buckets(NULL), values(NULL), nb_elements(0)
{}

Max_Bucket::~Max_Bucket() {
	if (buckets != NULL)
		free (buckets);
	if (elements != NULL)
		free (elements);
	if (values != NULL)
		free (values);
}
void Max_Bucket::Free () {
	free (buckets);
	buckets = NULL;
	free (elements);
	elements = NULL;
	free (values);
	values = NULL;
}
void Max_Bucket::Initialize(int max_v, int nb_element) {
	int i;
	max_value = max_v;
	buckets = (Max_Bucket_element **) malloc(sizeof(Max_Bucket_element *) * (max_value+1));
	elements = (Max_Bucket_element *) malloc(sizeof(Max_Bucket_element) * nb_element);
	//    printf("nb_element: %d\n", nb_element);
	values = (int *) malloc(sizeof(int) * nb_element);
	for (int i = 0; i < nb_element; i++)
		values[i] = -1;
//	nb_elements = nb_element;
	if (buckets == NULL || elements == NULL || values == NULL) {
		free(values);
		free(buckets);
		free(elements);
	}
	else {
		for (i = 0; i <= max_value; i++)
			buckets[i] = NULL;
		for (i = 0; i < nb_element; i++) {
			elements[i].prev = NULL;
			elements[i].next = NULL;
		}
	}
	current_max_value = 0;
}

void Max_Bucket::Insert (int id, int value) {
	values[id] = value;
	elements[id].prev = NULL;
	elements[id].next = buckets[value];
	if (buckets[value] != NULL)
		buckets[value]->prev = &(elements[id]);
	else if (current_max_value < value)
		current_max_value = value;
	buckets[value] = &(elements[id]);
	nb_elements++;
}

int Max_Bucket::PopMax(int* id, int* ret_value) {
	for (; current_max_value >= 0; current_max_value--) {
		if (buckets[current_max_value] != NULL) {
			*id = buckets[current_max_value] - elements;
			buckets[current_max_value] = buckets[current_max_value]->next;
			if (buckets[current_max_value] != NULL)
				buckets[current_max_value]->prev = NULL;
			*ret_value = values[*id];
			values[*id] = INT_MAX;
			nb_elements--;
			return 0;
		}
	}
	*id = -19;
	return -1;
}

int Max_Bucket::CurrentValue(int id) {
	return values[id];
}


void Max_Bucket::Update(int id, int new_value) {
	int old_value = values[id];

	if (old_value == INT_MAX)
			return;

	values[id] = new_value;

	if (elements[id].prev == NULL)
		buckets[old_value] = elements[id].next;
	else
		elements[id].prev->next = elements[id].next;

	if (elements[id].next != NULL)
		elements[id].next->prev = elements[id].prev;

	Max_Bucket::Insert(id, new_value);
}


//void Max_Bucket::DecVal(int id) {
////	printf ("%d is decremented to %d\n", id, values[id]-1);
//	int old_value = values[id];
//	// adjust the prev and next pointers
//	if (elements[id].prev == NULL)
//		buckets[old_value] = elements[id].next;
//	else
//		elements[id].prev->next = elements[id].next;
//	if (elements[id].next != NULL)
//		elements[id].next->prev = elements[id].prev;
//	Max_Bucket::Insert(id, values[id]-1);
//	return;
//}

