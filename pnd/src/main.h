#include <algorithm>
#include <iostream>
#include <queue>
#include <sstream>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <utility>
#include <random>
#include <chrono>
#include <tuple>
#include <fstream>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <omp.h>

#include "bucket.h"
#include "timestamp.cpp"
#include "larray.h"

#define PRIME 251231
using namespace std;
using namespace util;

typedef long long lol;
typedef int vertex; //vertices are 32 bytes
typedef int edge; //edges are 32 bytes

typedef tuple<vertex, vertex> couple;
typedef tuple<vertex, vertex, vertex> triple;

struct triangle_id {
	tuple<vertex, vertex, vertex> triple;
	vertex id;
	bool operator==(const triangle_id &other) const {
		if (triple == other.triple)
			return 1;
		else {
			return 0;
		}
	}
	triangle_id () {
		triple = make_tuple (-1, -1, -1);
		id = -1;
	}

};

namespace std {

template <>
struct hash<triangle_id> {
	std::size_t operator()(const triangle_id& t) const {
		return (get<0>(t.triple) * PRIME * PRIME + get<1>(t.triple)) * PRIME + get<2>(t.triple);
	}
};

template <>
struct hash<couple> {
	std::size_t operator()(const couple& c) const {
		return (get<0>(c) * PRIME + get<1>(c));
	}
};

template <>
struct hash<triple> {
	std::size_t operator()(const triple& c) const
	{
		return (get<0>(c) * PRIME * PRIME + get<1>(c) * PRIME + get<2>(c));
	}
};

}

typedef vector<vector<vertex> > Graph;

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
    template<typename S, typename T> struct hash<pair<S, T>>
    {
        inline size_t operator()(const pair<S, T> & v) const
        {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };
}

inline vertex findInd (vertex a, vertex i, vertex* ordered_adj, edge* ordered_xadj) {
	for (vertex j = ordered_xadj[a]; j < ordered_xadj[a+1]; j++)
		if (ordered_adj[j] == i)
			return j;
	return -1;
}

inline bool isSmaller (vertex* xadj, vertex u, vertex v) {
	vertex deg_u = xadj[u+1] - xadj[u];
	vertex deg_v = xadj[v+1] - xadj[v];
	return (deg_u < deg_v || (deg_u == deg_v && u < v));
}

typedef tuple<int, int> eda;

inline bool kksort (eda i, eda j) { return (get<1>(i) > get<1>(j)); }

inline void print_Ks (edge nVtx, vertex* T, const char* vfile, int H = -1) {
	string st (vfile);
	if (H == -1)
		st += "_FINAL_K";
	else
		st += "_H_" + to_string(H);
	FILE* pp = fopen (st.c_str(), "w");
	for (edge i = 0; i < nVtx; i++)
		fprintf (pp, "%d\n", T[i]);
	fclose (pp);
}

inline void read_Ks (size_t sz, const char* fl, vertex** P) {
	string st (fl);
	*P = (vertex *) malloc (sizeof(vertex) * sz);
	FILE* fp = fopen (st.c_str(), "r");
	vertex num;
	for (size_t i = 0; i < sz; i++) {
		fscanf (fp, "%d", &num);
		if (num == -1)
			(*P)[i] = 0;
		else
			(*P)[i] = num;
	}
	fclose (fp);
}

inline void intersection2 (vertex* adj, edge* xadj, vertex u, vertex v, vector<vertex>& intersection) {
	vertex i = xadj[u];
	vertex j = xadj[v];
	vertex gu = xadj[u+1];
	vertex gv = xadj[v+1];

	while (i < gu && j < gv) {
		if (adj[i] < adj[j])
			i++;
		else if (adj[j] < adj[i])
			j++;
		else {
			intersection.push_back(adj[i]);
			i++;
			j++;
		}
	}
}

inline void intersection3 (vertex* adj, edge* xadj, vertex u, vertex v, vertex w, vector<vertex>& intersection) {
	vertex i = xadj[u];
	vertex j = xadj[v];
	vertex k = xadj[w];
	vertex gu = xadj[u+1];
	vertex gv = xadj[v+1];
	vertex gw = xadj[w+1];
	while (i < gu && j < gv && k < gw) {
		vertex a = adj[i];
		vertex b = adj[j];
		vertex c = adj[k];
		if (a == b && a == c) {
			intersection.push_back(a);
			i++; j++; k++;
		}
		else {
			vertex m = max (a, max (b, c));
			if (a != m)
				i++;
			if (b != m)
				j++;
			if (c != m)
				k++;
		}
	}
}

inline void createOrdered (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, couple* el, edge* xel, vertex* ordered_adj, edge* ordered_xadj) {
	edge xi = 0;
	vertex i = 0;
	xel[xi++] = 0;
	edge oxi = 0;
	vertex oi = 0;
	ordered_xadj[oxi++] = 0;
	for (vertex u = 0; u < nVtx; u++) {
		for (vertex j = xadj[u]; j < xadj[u+1]; j++) {
			vertex v = adj[j];
			if (isSmaller (xadj, u, v)) {
				ordered_adj[oi++] = v;
				couple c = make_tuple(u, v);
				el[i++] = c;
			}
		}
		ordered_xadj[oxi++] = oi;
		xel[xi++] = i;
	}
}

void baseLocal12 (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile);
void nmLocal12 (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile);
void kcore (vertex nVtx, vertex* adj, edge* xadj, vertex* K, const char* vfile);

void baseLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void ktruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);

void baseLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void k34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);

template <typename VtxType, typename EdgeType>
void readGraph (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj);


