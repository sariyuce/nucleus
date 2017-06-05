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

#define VERTEXUPPERBOUND 1000 // for TIP variants, compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time
#define EDGEUPPERBOUND 500 // for WING, compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time
#define LOWERBOUND 0 // show subgraphs with at least this size
#define PRIME 251231 // for hash function

typedef int vertex; //vertices are 32 bytes
typedef int edge; //edges are 32 bytes
typedef unordered_multimap<int, int> mmap;
typedef pair<vertex, vertex> vp;
typedef vector<vector<vertex> > Graph;

struct subcore {
	bool visible;
	vertex rank;
	vertex K;
	vertex parent;
	vertex root;
	vector<vertex> children;
	vertex primarySize;
	vertex secondarySize;
	edge nEdge;
	double ed;

	subcore(vertex k) {
		K = k;
		rank = 0;
		parent = -1;
		root = -1;
		visible = true;
	}
};

struct helpers {
	helpers (vector<vp>* ael) {
		el = ael;
	}
	helpers () {}

	vector<vp>* el;
};

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

inline void print_time (FILE* fp, const string& str, const timestamp& t) {
	fprintf (fp, "%s %d.%06d\n", str.c_str(), t.seconds, t.microseconds);
	fflush(fp);
}

inline vertex find_ind (vector<vertex>& rg, vertex g) {
	for (vertex i = 0; i < rg.size(); i++)
		if (rg[i] == g)
			return i;
	printf ("couldn't find\n");
	exit(1);
}


inline void prefixSum (vector<vertex>& xRight, Graph& rightGraph, vector<vp>& el) {
	vertex sum = 0;
	xRight.push_back (sum);
	for (vertex u = 0; u < rightGraph.size(); u++) {
		for (vertex v = 0; v < rightGraph[u].size(); v++) {
			vp p (u, rightGraph[u][v]);
			el.push_back (p);
		}
		sum += rightGraph[u].size();
		xRight.push_back (sum);
	}
}

inline bool hashUniquify(vector<vertex>& ch, bool fl = true, string variant = "NONE") {
	HashMap<bool> hermap(false);
	for (size_t i = 0; i < ch.size(); i++) {
		int t = ch[i];
		if (hermap.hasDefaultValue(t)) {
			hermap[t] = true;
		}
		else {
			ch.erase(ch.begin() + i);
			i--;
		}
		if (fl && (variant == "WING" && i > EDGEUPPERBOUND)
				|| ((variant == "RIGHT_TIP" || variant == "LEFT_TIP") && i > VERTEXUPPERBOUND))
			return false;

//		if (i % 100000 == 0)
//			printf ("i : %d\n", i);
	}

//	printf ("ch.size: %d\n", ch.size());
	if (fl)
		sort(ch.begin(), ch.end());
	return true;
}

inline void findRepresentative(vertex* child, vector<subcore>& skeleton) {
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

inline void assignToRepresentative(vertex* ch, vector<subcore>& skeleton) { // 2-pass path compression
	vertex u = *ch;
	vector<int> vs;
	while (skeleton[u].parent != -1) {
		int n = skeleton[u].parent;
		if (skeleton[n].K == skeleton[u].K) {
			vs.push_back(u);
			u = n;
		} else
			break;
	}
	*ch = u;
	for (vertex i : vs) {
		if (i != u)
			skeleton[i].parent = u;
	}
}

inline void neighborsOfNeighbors(vector<vertex>& vset, Graph& graph,
		vector<vertex>& allNeighbors, edge* nEdge) {
	for (vertex i = 0; i < vset.size(); i++) {
		vertex u = vset[i];
		allNeighbors.insert(allNeighbors.end(), graph[u].begin(),
				graph[u].end());
		*nEdge += graph[u].size();
	}
	hashUniquify(allNeighbors);
}

inline void store(vertex uComp, vertex vComp, vector<vertex>& unassigned,
		vector<vp>& relations) {
	vp c(vComp, uComp);
	if (uComp == -1) // it is possible that u didn't get an id yet
		unassigned.push_back(relations.size()); // keep those indices to process after the loop, below
	relations.push_back(c);
}

inline void merge(vertex u, vertex v, vector<vertex>& component,
		vector<subcore>& skeleton, int* nSubcores) {

	if (component[u] == -1) {
		component[u] = component[v];
		skeleton.erase(skeleton.end() - 1);
	} else { // merge component[u] and component[v] nodes
		vertex child = component[u];
		vertex parent = component[v];
		assignToRepresentative(&child, skeleton);
		assignToRepresentative(&parent, skeleton);
		if (child != parent) {
			if (skeleton[child].rank > skeleton[parent].rank)
				swap(child, parent);
			skeleton[child].parent = parent;
			skeleton[child].visible = false;
			if (skeleton[parent].rank == skeleton[child].rank)
				skeleton[parent].rank++;
			*nSubcores--;
		}
	}
}

inline void createSkeleton(vertex u, initializer_list<vertex> neighbors,
		vertex* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,
		vector<vertex>& component, vector<vertex>& unassigned,
		vector<vp>& relations) {
	vertex smallest = -1, minK = INT_MAX;
	for (auto i : neighbors)
		if (K[i] != -1 && K[i] < minK) {
			smallest = i;
			minK = K[i];
		}
	if (smallest == -1)
		return;

	if (K[smallest] == K[u])
		merge(u, smallest, component, skeleton, nSubcores);
	else
		store(component[u], component[smallest], unassigned, relations);
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

inline void indicesIntersection(vector<vertex>& a, vector<vertex>& b,
		vector<vertex>& res, vertex g) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
//			if (a[i] != g)
			{
				res.push_back(i);
				res.push_back(j);
			}
			i++;
			j++;
		}
	}
}

inline int nChoosek(int n, int k) {
	if (k > n)
		return 0;
	if (k * 2 > n)
		k = n - k;
	if (k == 0)
		return 1;

	int result = n;
	for (int i = 2; i <= k; ++i) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}


inline void intersectionM1 (vector<vertex>& a, vector<vertex>& b, vector<vertex>& c) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {

		bool fl = false;
		if (b[j] == -1) {
			j++;
			fl = true;
		}
		if (a[i] == -1) {
			i++;
			fl = true;
		}
		if (fl)
			continue;

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
void ReadBipartiteGraph(char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph);
void tipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);
void hopefullyFastertipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);
void insallahFasterTipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);

void wingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);
void insallahFasterwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);
void hopefullyFasterwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount);
void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex rightnVtx, vertex leftnVtx);
void presentNuclei (string variant, vector<subcore>& skeleton, vector<vertex>& component, edge nEdge, helpers& ax, string vfile, FILE* gp, Graph& leftGraph, Graph& rightGraph, vector<vertex>* xRight);
