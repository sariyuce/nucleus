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

#define VERTEXUPPERBOUND 200 // for TIP variants, compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time
#define EDGEUPPERBOUND 3000 // for WING, compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time
#define LOWERBOUND 0 // show subgraphs with at least this size
#define PRIME 251231 // for hash function

typedef int vertex; //vertices are 32 bytes
typedef int edge; //edges are 32 bytes
typedef unordered_multimap<int, int> mmap;
typedef pair<vertex, vertex> vp;
typedef vector<vector<vertex> > Graph;
struct wv {
    int n;
    double w;
};
typedef vector<vector<wv>> Wraph;

typedef long long ll;
typedef pair<ll, ll> llp;







// User defined class, Point
class vd
{
public:
   int v;
   long long deg;
   vd (int va, long long dega)
   {
	   v = va;
	   deg = dega;
   }

};

// To compare two points
class myComparator
{
public:
    int operator() (const vd& p1, const vd& p2)
    {
        return p1.deg > p2.deg;
    }
};









struct subcore {
	bool visible;
	vertex rank;
	vertex K;
	ll parent;
	ll root;
	vector<ll> children;
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
//		primarySize = secondarySize = -1;
//		ed = -1;
//		children.clear();
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

inline bool less_than (vertex u, vertex v, Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
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
		if (fl &&
				((variant == "WING" && i > EDGEUPPERBOUND) ||
						(variant == "TIP" && i > VERTEXUPPERBOUND) ||
						(variant == "CORE" && i > VERTEXUPPERBOUND) ||
						(variant == "TRUSS" && i > VERTEXUPPERBOUND)))
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

inline void assignToRepresentative(ll* ch, vector<subcore>& skeleton) { // 2-pass path compression
	ll u = *ch;
	vector<ll> vs;
	while (skeleton[u].parent != -1) {
		ll n = skeleton[u].parent;
		if (skeleton[n].K == skeleton[u].K) {
			vs.push_back(u);
			u = n;
		} else
			break;
	}
	*ch = u;
	for (ll i : vs) {
		if (i != u)
			skeleton[i].parent = u;
	}
}

inline void neighborsOfNeighbors(vector<vertex>& vset, Graph& graph, vector<vertex>& allNeighbors, edge* nEdge) {
	for (vertex i = 0; i < vset.size(); i++) {
		vertex u = vset[i];
		allNeighbors.insert(allNeighbors.end(), graph[u].begin(), graph[u].end());
		*nEdge += graph[u].size();
	}
	hashUniquify(allNeighbors);
}

inline void store(ll uComp, ll vComp, vector<ll>& unassigned,
		vector<llp>& relations) {
	llp c(vComp, uComp);
	if (uComp == -1) // it is possible that u didn't get an id yet
		unassigned.push_back(relations.size()); // keep those indices to process after the loop, below
	relations.push_back(c);
}

inline void merge(vertex u, vertex v, vector<ll>& component,
		vector<subcore>& skeleton, ll* nSubcores) {

	if (component[u] == -1) {
		component[u] = component[v];
		skeleton.erase(skeleton.end() - 1);
	}
	else { // merge component[u] and component[v] nodes
		ll child = component[u];
		ll parent = component[v];
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
		ll* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,
		vector<ll>& component, vector<ll>& unassigned,
		vector<llp>& relations) {
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

inline void updateUnassigned (vertex t, vector<ll>& component, ll* cid, vector<llp>& relations, vector<ll>& unassigned) {
	if (component[t] == -1) { // if e didn't get a component, give her a new one
		component[t] = *cid;
		++(*cid);
	}

	// update the unassigned components that are in the relations
	for (ll i : unassigned)
		relations[i] = make_pair (relations[i].first, component[t]);
}

inline void indicesIntersection (vector<vertex>& a, vector<vertex>& b,
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

inline void indicesIntersectionOld (vector<vertex>& a, vector<vertex>& b,
		vector<vertex>& res, vertex g) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			if (a[i] != g) {
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

inline void ff_vertex_ind_density (string file, Graph& graph) {
	FILE* fp = fopen (file.c_str(), "r");
	string temp = file + "_temp";
	FILE* tp = fopen (temp.c_str(), "w");
    int t;
    int ln = 0;
    double db;

    while (1) {
        fscanf (fp, "%d %d %d %d %d %lf %d %d", &t, &t, &t, &t, &t, &db, &t, &t);
        vector<int> vs;
        bool exit_flag = false;
        while (1) {
            if (fscanf (fp, "%d", &t) == EOF)
                exit_flag = true;
            if (t == -1)
                break;
            vs.push_back(t);
        }
        if (exit_flag)
        	break;
        vector<vertex> neigsOfvs;
        int nedge = 0;
        neighborsOfNeighbors (vs, graph, neigsOfvs, &nedge);
        fprintf (tp, " %d %d %d %lf\n", vs.size(), neigsOfvs.size(), nedge, vs.size() == 0 ? -1.0 : ((double) nedge / (vs.size() * neigsOfvs.size())));
    }

    fclose (fp);
    fclose (tp);
    string paste = "paste " + file + " " + temp + " > " + file + temp;
//    cout << paste << endl;
    system (paste.c_str());
    string mv = "mv " + file + temp + " " + file;
//    cout << mv << endl;
    system (mv.c_str());
    string rm = "rm " + temp;
//    cout << rm << endl;
    system (rm.c_str());
}

inline void ff_edge_ind_density (string file, Graph& graph) {
	FILE* fp = fopen (file.c_str(), "r");
	string temp = file + "_temp";
	FILE* tp = fopen (temp.c_str(), "w");

    int u, v;
    int ln = 0;
    double db;


    while (1) {
    	fscanf (fp, "%d %d %d %d %d %lf %d %d", &u, &u, &u, &u, &u, &db, &u, &u);
        vector<int> vs;
        bool exit_flag = false;
        while (1) {
            if (fscanf (fp, "%d", &u) == EOF)
                exit_flag = true;
            if (u == -1)
                break;
            fscanf (fp, "%d", &v);
            vs.push_back(u);
            vs.push_back(v);
        }
        if (exit_flag)
        	break;
        hashUniquify (vs);
        vector<vertex> neigsOfvs;
        int nedge = 0;
        neighborsOfNeighbors (vs, graph, neigsOfvs, &nedge);
        fprintf (tp, " %d %d %d %lf\n", vs.size(), neigsOfvs.size(), nedge, vs.size() == 0 ? -1.0 : ((double) nedge / (vs.size() * neigsOfvs.size())));
    }
    fclose (fp);
    fclose (tp);
    string paste = "paste " + file + " " + temp + " > " + file + temp;
//    cout << paste << endl;
    system (paste.c_str());
    string mv = "mv " + file + temp + " " + file;
//    cout << mv << endl;
    system (mv.c_str());
    string rm = "rm " + temp;
//    cout << rm << endl;
    system (rm.c_str());
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

template <typename VtxType, typename EdgeType>
void ReadBipartiteGraph(char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph);
void ReadRegularGraph (char *filename, Graph& graph, int* nEdge);
void ReadWeightedGraph (char *filename, Wraph& wraph, int* nEdges);

void writeToRegularBinary (string filename, vertex nVtx, edge nEdge, Graph& graph);

void writeToWeightedRegularBinary (string filename, vertex nVtx, edge nEdge, Wraph& wraph);

void tipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, long long* bCount);
void oldtipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, long long* bCount);

void wingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, long long* bCount);
void oldwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, long long* bCount);
void buildHierarchy (vertex cn, vector<llp>& relations, vector<subcore>& skeleton, ll* nSubcores, edge nEdge, vertex rightnVtx, vertex leftnVtx = -1);
//void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex rightnVtx, vertex leftnVtx = -1);
void presentNuclei (string variant, vector<subcore>& skeleton, vector<vertex>& component, edge nEdge, helpers& ax, string vfile, FILE* gp, Graph& leftGraph, Graph& rightGraph, vector<vertex>* xRight);

void unweighted_projection (Graph& left, Graph& right, string filename);
void weighted_projection (Graph& left, Graph& right, string filename);
void weighted_base_kcore (Wraph& wraph, int nEdge, vector<vertex>& K, bool hierarchy, int* maxCore, const char* vfile, FILE* fp);
void base_kcore (Graph& graph, int nEdge, vector<vertex>& K, bool hierarchy, int* maxCore, string vfile, FILE* fp);
void base_ktruss (Graph& graph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxtruss, string vfile, FILE* fp);
