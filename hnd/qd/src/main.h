#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <chrono>
#include <sys/stat.h>
#include "bucket.h"

using namespace std;

#define DEG_DIST 0
#define COUNT_ONLY 0 // make it 1 to count motifs only (and terminate)
#define LOWERBOUND 0
#define UPPERBOUND 5000 // compute densities of subgraphs with at most this size, set to INT_MAX to compute all but it could take a lot of time
#define THRESHOLD 0.0
#define PRIME 251231 // for hash function

typedef long long lol;
typedef int vertex; // vertices are 32 bytes
typedef int edge; // edges are 32 bytes
typedef unordered_multimap<int, int> mmap;
typedef pair<vertex, vertex> vp;
typedef tuple<vertex, vertex, vertex> vt;
typedef vector<vector<vertex> > Graph;

typedef chrono::duration<double> tms;
typedef tuple<vertex, vertex> couple;
typedef tuple<vertex, vertex, vertex> triple;

#ifdef SIGNS
	extern Graph signs;
	extern int MODE;
#endif

struct subcore {
	bool visible;
	vertex rank;
	vertex K;
	vertex parent;
	vertex root;
	vector<vertex> children;
	vertex size;
	edge nEdge;
	double ed;

	subcore (vertex k) {
		K = k;
		rank = 0;
		parent = -1;
		root = -1;
		visible = true;
		size = 0;
		nEdge = 0;
		ed = -1;
	}
};

struct helpers {
	helpers (vector<vp>* ael) {
		el = ael;
	}
	helpers (vector<vt>* atris) {
		tris = atris;
	}
	helpers () {}

	vector<vp>* el;
	vector<vt>* tris;
};


inline vertex P2M_ (vertex a) {return (-1 * a) - 1;}

inline vertex M2P (vertex a) {	return (-1 * (a + 1));}

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std{
namespace
{

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     https://stackoverflow.com/questions/4948780

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
	seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// Recursive template code derived from Matthieu M.
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct HashValueImpl
{
	static void apply(size_t& seed, Tuple const& tuple)
	{
		HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
		hash_combine(seed, get<Index>(tuple));
	}
};

template <class Tuple>
struct HashValueImpl<Tuple,0>
{
	static void apply(size_t& seed, Tuple const& tuple)
	{
		hash_combine(seed, get<0>(tuple));
	}
};
}

template <typename ... TT>
struct hash<std::tuple<TT...>>
{
	size_t
	operator()(std::tuple<TT...> const& tt) const
	{
		size_t seed = 0;
		HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
		return seed;
	}

};
}

inline bool uniquify (vector<vertex>& vertices) {
	unordered_map<vertex, bool> hermap;
	for (size_t i = 0; i < vertices.size(); i++) {
		int t = vertices[i];
		if (hermap.find (t) == hermap.end())
			hermap[t] = true;
		else {
			vertices.erase (vertices.begin() + i);
			i--;
		}
	}
	sort (vertices.begin(), vertices.end());
	return true;
}


inline bool hashUniquify (vector<vertex>& vertices) {
	unordered_map<vertex, bool> hermap;
	for (size_t i = 0; i < vertices.size(); i++) {
		int t = vertices[i];
		if (hermap.find (t) == hermap.end())
			hermap[t] = true;
		else {
			vertices.erase (vertices.begin() + i);
			i--;
		}
		if (i > UPPERBOUND)
			return false;
	}
	sort (vertices.begin(), vertices.end());
	return true;
}

inline void print_time (FILE* fp, const string& str, tms t) {
	fprintf (fp, "%s %.6lf\n", str.c_str(), t.count());
	fflush(fp);
}

inline bool less_than (vertex u, vertex v, Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
}

inline bool orderedConnected (Graph& graph, Graph& orderedGraph, vertex u, vertex v) {
	vertex a = u, b = v;
	if (less_than (v, u, graph))
		swap (a, b);
	for (auto c : orderedGraph[a])
		if (c == b)
			return true;
	return false;
}

inline void createOrderedIndexEdges (Graph& graph, vector<vp>& el, vector<vertex>& xel, Graph& orderedGraph) {
	xel.push_back(0);
	orderedGraph.resize(graph.size());
	for (size_t u = 0; u < graph.size(); u++) {
		for (size_t j = 0; j < graph[u].size(); j++) {
			vertex v = graph[u][j];
			if (less_than (u, v, graph)) {
				orderedGraph[u].push_back(v);
				vp c (u, v);
				el.push_back(c);
			}
		}
		xel.push_back(el.size());
	}
}

inline void createOrdered (Graph& orderedGraph, Graph& graph) {
	orderedGraph.resize (graph.size());
	for (vertex i = 0; i < graph.size(); ++i)
		for (auto u : graph[i])
			if (less_than (i, u, graph))
				orderedGraph[i].push_back(u);
}

inline vertex getEdgeId (vertex u, vertex v, vector<vertex>& xel, vector<vp>& el, Graph& graph) {

	vertex a = u, b = v;

	for (vertex i = xel[a]; i < xel[a+1]; i++)
		if (el[i].second == b)
			return i;

	printf ("getEdgeId returns -1\n");
	exit(1);
}

inline void intersection (vector<vertex>& a, vector<vertex>& b, vector<vertex>& c) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			c.push_back(a[i]);
			i++;
			j++;
		}
	}
}

inline int order_compare (vertex* Reals, vertex* T, vertex e, vertex f) {
	if ((Reals[e] < Reals[f] && T[e] < T[f]) || (Reals[e] > Reals[f] && T[e] > T[f]) || (Reals[e] == Reals[f] && T[e] == T[f]))
		return 1;
	else
		return -1;
}

inline void EdgeKendallTau (vertex* Reals, vertex* T, vertex a, vertex b, int* count, int* score) {
	int s = 0;
	s += order_compare (Reals, T, a, b);
	if (s == 1) // concordant -- at least one 1, no -1
		(*score)++;
	else
		(*score)--;
	(*count)++;
}

inline void TriangleKendallTau (vertex* Reals, vertex* T, vertex a, vertex b, vertex c, int* count, int* score) {
	int s = 0;
	s += order_compare (Reals, T, a, b);
	s += order_compare (Reals, T, a, c);
	s += order_compare (Reals, T, b, c);
	if (s == 3) // concordant -- at least one 1, no -1
		(*score)++;
	else
		(*score)--;
	(*count)++;
}

// a, b, c, d are ids
inline void FourCliqueKendallTau (vertex* Reals, vertex* T, vertex a, vertex b, vertex c, vertex d, int* count, int* score) {
	int s = 0;
	s += order_compare (Reals, T, a, b);
	s += order_compare (Reals, T, a, c);
	s += order_compare (Reals, T, a, d);
	s += order_compare (Reals, T, b, c);
	s += order_compare (Reals, T, b, d);
	s += order_compare (Reals, T, c, d);
	if (s == 6) // concordant -- at least one 1, no -1
		(*score)++;
	else
		(*score)--;
	(*count)++;
}


inline double kt (int count, int score) {
	return (double) score / count;
}

inline bool isSmaller (vertex* xadj, vertex u, vertex v) {
	vertex deg_u = xadj[u+1] - xadj[u];
	vertex deg_v = xadj[v+1] - xadj[v];
	return (deg_u < deg_v || (deg_u == deg_v && u < v));
}

inline bool exists (int val, vector<int>& v) {
	for (size_t i = 1; i < v.size(); i++) {
		if (v[i] == val)
			return true;
		else if (v[i] < 0)
			return false;
	}
	return false;
}

inline int ind (int val, vector<int>& v) {
	for (size_t i = 1; i < v.size(); i++)
		if (v[i] == val)
			return i;
	return -1;
}


inline void checkAndDec (vertex kw, vertex kv, vertex w, vertex v, Naive_Bucket* nBucket, vertex tc) {
	if (kv == -1 && kw == -1) {
		if ((*nBucket).CurrentValue(v) > tc)
			(*nBucket).DecVal(v);
		if ((*nBucket).CurrentValue(w) > tc)
			(*nBucket).DecVal(w);
	}
}

inline void checkAndDecAndHier (vertex w, vertex v, Naive_Bucket* nBucket, vertex tc, vertex u,
		bool hierarchy, vertex* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,
		vector<vertex>& component, vector<vertex>& unassigned, vector<vp>& relations) {
	if (K[v] == -1 && K[w] == -1) {
		if ((*nBucket).CurrentValue(v) > tc)
			(*nBucket).DecVal(v);
		if ((*nBucket).CurrentValue(w) > tc)
			(*nBucket).DecVal(w);
	}
	else if (hierarchy)
		createSkeleton (u, {v, w}, nSubcores, K, skeleton, component, unassigned, relations);
}



void inte (vector<vertex>& a, vector<vertex>& b, vector<vertex>& c, vector<vertex>& d, vector<vertex>& ret);
void inter (int F, int G, Graph& dg, vertex a, vertex b, vector<vertex>& ret);

template <typename VtxType, typename EdgeType>
void readGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge);

template <typename VtxType, typename EdgeType>
void readDirectedGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge);

void createSkeleton (vertex u, initializer_list<vertex> neighbors, vertex* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,	vector<vertex>& component, vector<vertex>& unassigned, vector<vp>& relations);
void updateUnassigned (vertex t, vector<vertex>& component, vertex* cid, vector<vp>& relations, vector<vertex>& unassigned);
void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex nVtx);
void presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, edge nEdge, helpers& ax, string vfile, FILE* gp);

void outgoings (vector<vertex>& b, vector<vertex>& ret);
void outgoings_and_asymmetric_undirecteds (vertex u, vector<vertex>& b, vector<vertex>& ret);
void asymmetric_undirecteds (vertex u, vector<vertex>& b, vector<vertex>& ret);
void undirecteds (vector<vertex>& b, vector<vertex>& ret);
void incomings (vector<vertex>& a, vector<vertex>& ret);

void cycle_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);
void acyclic_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);
void outp_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);
void cyclep_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);
void inp_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);
void cyclepp_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp);

double simple_count_acyclics (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
double simple_count_cycles (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
double simple_count_cycleps (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
double simple_count_cyclepps (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
double simple_count_outps (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
double simple_count_inps (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing);
