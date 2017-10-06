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
#include <sys/stat.h>

#include "timestamp.cpp"
#include "bucket.h"
#include "larray.h"


using namespace std;
using namespace util;


#define UPPERBOUND 200 // compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time
#define LOWERBOUND 0 // show subgraphs with at least this size
#define PRIME 251231 // for hash function


typedef int vertex; // vertices are 32 bytes
typedef int edge; // edges are 32 bytes
typedef unordered_multimap<int, int> mmap;
//typedef pair<vertex, vertex> vp;
typedef pair<vertex, vertex> vp;
typedef tuple<vertex, vertex, vertex> vt;
typedef vector<vector<vertex> > Graph;


struct subcore {
	bool visible;
	vertex rank;
	vertex K;
	vertex parent;
	vertex root;
	vector<vertex> children;
	vertex size;
	edge nedge;
	double ed;

	subcore (vertex k) {
		K = k;
		rank = 0;
		parent = -1;
		root = -1;
		visible = true;
		size = 0;
		nedge = 0;
		ed = -1;
	}
};
struct triangleId {
	vt triple;
	vertex id;

	triangleId () {
		triple = make_tuple (-1, -1, -1);
		id = -1;
	}

	triangleId (vertex a, vertex b, vertex c) {
		triple = make_tuple (a, b, c);
		id = -1;
	}

	bool operator==(const triangleId &other) const {
		if (triple == other.triple)
			return 1;
		else {
			return 0;
		}
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

namespace std {
template <>
struct hash<triangleId> {
	std::size_t operator()(const triangleId& t) const {
		return (get<0>(t.triple) * PRIME * PRIME + get<1>(t.triple)) * PRIME + get<2>(t.triple);
	}
};
}


inline void print_time (FILE* fp, const string& str, const timestamp& t) {
	fprintf (fp, "%s %d.%06d\n", str.c_str(), t.seconds, t.microseconds);
	fflush(fp);
}
inline void findRepresentative (vertex* child, vector<subcore>& skeleton) {
	vertex u = *child;
	if (skeleton[u].parent != -1) {
		vertex pr = skeleton[u].parent;
		while (skeleton[u].K == skeleton[pr].K) {

			u = pr;
			if (skeleton[u].parent != -1)
				pr = skeleton[u].parent;
			else
				break;
		}
	}
	*child = u;
}
inline void assignToRepresentative (vertex* ch, vector<subcore>& skeleton) { // 2-pass path compression
	vertex u = *ch;
	vector<int> vs;
	while (skeleton[u].parent != -1) {
		vertex n = skeleton[u].parent;
		if (skeleton[n].K == skeleton[u].K) {
			vs.push_back (u);
			u = n;
		}
		else
			break;
	}
	*ch = u;
	for (vertex i : vs) {
		if (i != u)
			skeleton[i].parent = u;
	}
}
inline void store (vertex uComp, vertex vComp, vector<vertex>& unassigned, vector<vp>& relations) {
	vp c (vComp, uComp);
	if (uComp == -1) // it is possible that u didn't get an id yet
		unassigned.push_back (relations.size()); // keep those indices to process after the loop, below
	relations.push_back (c);
}
inline void merge (vertex u, vertex v, vector<vertex>& component, vector<subcore>& skeleton, int* nSubcores) {

	if (component[u] == -1) {
		component[u] = component[v];
		skeleton.erase (skeleton.end() - 1);
	}
	else { // merge component[u] and component[v] nodes
		vertex child = component[u];
		vertex parent = component[v];
		assignToRepresentative (&child, skeleton);
		assignToRepresentative (&parent, skeleton);
		if (child != parent) {
			if (skeleton[child].rank > skeleton[parent].rank)
				swap (child, parent);
			skeleton[child].parent = parent;
			skeleton[child].visible = false;
			if (skeleton[parent].rank == skeleton[child].rank)
				skeleton[parent].rank++;
			*nSubcores--;
		}
	}
}
inline void createSkeleton (vertex u, initializer_list<vertex> neighbors, vertex* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,
		vector<vertex>& component, vector<vertex>& unassigned, vector<vp>& relations) {
	vertex smallest = -1, minK = INT_MAX;
	for (auto i : neighbors)
		if (K[i] != -1 && K[i] < minK) {
			smallest = i;
			minK = K[i];
		}
	if (smallest == -1)
		return;

	if (K[smallest] == K[u])
		merge (u, smallest, component, skeleton, nSubcores);
	else
		store (component[u], component[smallest], unassigned, relations);
}
inline void updateUnassigned (vertex t, vector<vertex>& component, vertex* cid, vector<vp>& relations, vector<vertex>& unassigned) {
	if (component[t] == -1) { // if e didn't get a component, give her a new one
		component[t] = *cid;
		++(*cid);
	}

	// update the unassigned components that are in the relations
	for (vertex i : unassigned)
		relations[i] = make_pair (relations[i].first, component[t]);
}
inline bool less_than (vertex u, vertex v, Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
}
inline bool orientedConnected (Graph& graph, Graph& orientedGraph, vertex u, vertex v) {
	vertex a = u, b = v;
	if (less_than (v, u, graph))
		swap (a, b);
	for (auto c : orientedGraph[a])
		if (c == b)
			return true;
	return false;
}
inline void createOriented (Graph& orientedGraph, Graph& graph) {
	orientedGraph.resize (graph.size());
	for (vertex i = 0; i < graph.size(); ++i)
		for (auto u : graph[i])
			if (less_than (i, u, graph))
				orientedGraph[i].push_back(u);
}
inline vertex getEdgeId (vertex u, vertex v, vector<vertex>& xel, vector<vp>& el, Graph& graph) {

	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);

	for (vertex i = xel[a]; i < xel[a+1]; i++)
		if (el[i].second == b)
			return i;

	printf ("getEdgeId returns -1\n");
	exit(1);
}
inline void indexEdges (Graph& graph, vector<vp>& el, vector<vertex>& xel, Graph& orientedGraph) {
	xel.push_back(0);
	orientedGraph.resize(graph.size());
	for (size_t u = 0; u < graph.size(); u++) {
		for (size_t j = 0; j < graph[u].size(); j++) {
			vertex v = graph[u][j];
			if (less_than (u, v, graph)) {
				orientedGraph[u].push_back(v);
				vp c (u, v);
				el.push_back(c);
			}
		}
		xel.push_back(el.size());
	}
}
inline bool incrementTCIfConnected (Graph& graph, Graph& orientedGraph, Graph& TC, vertex u, vertex v) {
	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);
	for (size_t k = 0; k < orientedGraph[a].size(); k++)
		if (orientedGraph[a][k] == b) {
			TC[a][k]++;
			return true;
		}
	return false;
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

template <typename VtxType, typename EdgeType>
void ReadGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge);
void base_kcore (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp);
void base_k13 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp);
void base_k14 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp);
void base_ktruss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp);
void base_k24 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* max24, string vfile, FILE* fp);
void base_k34 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* max34, string vfile, FILE* fp);
void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex nVtx);
void presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, edge nEdge, helpers& ax, string vfile, FILE* gp);
