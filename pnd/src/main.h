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

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <sys/stat.h>
#include <chrono>
#include <omp.h>
#include <tuple>
#include <fstream>
#include <string.h>
#include "bucket.h"
#include "timestamp.cpp"
#include "larray.h"


#define PRIME 251231
#define THRESHOLD 1

#define SUBGRAPHSIZELIMIT 500

using namespace std;
using namespace util;

typedef long long lol;
typedef int vertex; //vertices are 32 bytes
typedef int edge; //edges are 32 bytes

typedef tuple<vertex, vertex> couple;
typedef tuple<vertex, vertex, vertex> triple;
typedef unordered_multimap<int, int> mmap;


struct stats {
	int find_op;
	int union_op;
	int adjust_op;
	int total_size;
};

struct subcore {
	bool visible;
	int rank;
	int K;
	int parent;
	int root;
	vector<int> children;
	int size;
	int nedge;
	double ed;
};

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

struct report {
	vertex id;
	vertex max;
	vertex nvtx;
	vertex nedge;
	double edge_density;
	double triangle_density;
	bool isMaximal;
	int mapid;
	vector<vertex> vs;
};

struct recai {
	bool isMaximal;
	vertex id;
	vertex max;
	set<vertex> vs;
	set<couple> edges;
};

typedef vector<couple> EdgeList2;
typedef vector<vector<vertex> > Graph;

struct auxies {
	vector<vertex> xel;
	EdgeList2 el;
	vector<vertex> xtl;
	vector<triangle_id> tl;
	unordered_map<triangle_id, int> tlist;
};

struct p_auxies {
	vector<vertex>* xel;
	EdgeList2* el;
	vector<vertex>* xtl;
	vector<triangle_id>* tl;
	unordered_map<triangle_id, int>* tlist;
};

struct co {
	int id;
	int c;
};

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <class Tuple, std::size_t Index = std::tuple_size<Tuple>::value - 1>
struct tuple_hash_impl
{
    static inline void apply(std::size_t & seed, Tuple const & tuple)
    {
        tuple_hash_impl<Tuple, Index - 1>::apply(seed, tuple);
        hash_combine(seed, std::get<Index>(tuple));
    }
};

template <class Tuple>
struct tuple_hash_impl<Tuple, 0>
{
    static inline void apply(std::size_t & seed, Tuple const & tuple)
    {
        hash_combine(seed, std::get<0>(tuple));
    }
};

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

    template<typename ...Args> struct hash<tuple<Args...>>
    {
        inline size_t operator()(const tuple<Args...> & v) const
        {
            size_t seed = 0;
            tuple_hash_impl<tuple<Args...>>::apply(seed, v);
            return seed;
        }
    };
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


inline bool hashUniquify (vector<vertex>& vertices) {
	HashMap<bool> hermap (false);
	for (size_t i = 0; i < vertices.size(); i++) {
		int t = vertices[i];
		if (hermap.hasDefaultValue (t))
			hermap[t] = true;
		else {
			vertices.erase (vertices.begin() + i);
			i--;
		}
	}
	sort (vertices.begin(), vertices.end());
	return true;
}


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

void kcore_levels (vertex nVtx, vertex* adj, edge* xadj, vertex* L, const char* vfile);
void kcore_Sesh_levels (vertex nVtx, vertex* adj, edge* xadj, vertex* K, vertex* L, const char* vfile);

void fast12DegeneracyNumber (vertex nVtx, vertex* adj, edge* xadj, vertex* P, vertex topK);





void baseLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void ktruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);

void ktruss_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* L, const char* vfile);
void ktruss_Sesh_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, vertex* L, const char* vfile);

void fast23DegeneracyNumber (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, vertex topK);


void baseLocal23_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal23_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void ktruss_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);






void baseLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void k34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);

void k34_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* L, const char* vfile);

void fast34DegeneracyNumber (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, vertex topK);


void baseLocal34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void nmLocal34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
void k34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);


template <typename VtxType, typename EdgeType>
void readGraph (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj);






//void baseLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
//void nmLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile);
////template <typename VtxType, typename EdgeType>
////void readGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge);
//
//
////void ReadGraph(char *filename, vertex *numofvertex, edge **pxadj, vertex **padjncy, vertex **padjncyw, vertex **ppvw);
//
//

//
//
//
//void justpi34 (Graph& graph, int nEdge, volatile vertex* P, vector<vertex>& Reals, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp, string fname);
//
//void TryFasterPiTruss (Graph& graph, int nEdge, volatile vertex* P, vertex* Reals, util::timestamp& totaltime, int *pi_number, const char* vfile, FILE* fp);
//void TryFasterPi34 (Graph& graph, int nEdge, volatile vertex* P, vector<vertex>& Reals, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp, string fname);
//
//void k34_Seshlevels (Graph& graph, int nEdge, vector<vertex>& F, vector<vertex>& cono, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp);
//void ktruss_Seshlevels (Graph& graph, int nEdge, vector<vertex>& T, vertex* cono, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//void kcore_Seshlevels (Graph& graph, int nEdge, vertex *K, vertex *cono, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp);
//void k34_levels (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp);
//void ktruss_levels (Graph& graph, int nEdge, vector<vertex>& T, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//void kcore_levels (Graph& graph, int nEdge, vertex *K, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp);
//void base_pitruss_LessSpace (Graph& graph, int nEdge, vertex* T, vertex* Reals, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//void base_pi34_LessSpace (Graph& graph, int nEdge, vertex* F, vector<vertex>& Reals, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp, string s);
//long intersect3_count_classic_efficient (Graph& graph, Graph& ordered_graph, EdgeList2& el, vector<vertex>& xel, vector<vertex>& xtl, vector<triangle_id>& tl);
//
//void superTriangleList (Graph& graph, Graph& ordered_graph, vector<int>& xel, vector<int>& xtl, EdgeList2& el, vector<triangle_id>& tl, Graph& TF);
//void base_k34_Store4c (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp);
//void create_classic_triangleList (Graph& ordered_graph, EdgeList2& el, vector<triangle_id>& tl, vector<vertex>& xtl);
//
//
//void superTriangleList (Graph& graph, Graph& ordered_graph, EdgeList2& el, Fcounts& tlist, Graph& TF);
//int count_triangles_adv (Graph& graph, Graph& ordered_graph, vertex* T, EdgeList2& el, vector<vertex>& xel, Graph& tris);
//void base_ktruss_StoreTri (Graph& graph, int nEdge, vector<vertex>& T, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//
//void base_k34 (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp);
//void base_pi34 (Graph& graph, int nEdge, vertex* F, vector<vertex>& Reals, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp);
//vertex hm_getEdgeId3 (vertex u, vertex v, vertex w, unordered_map<triangle_id, int>& tlist, Graph& graph);
//
//
//void intersect3_new (Graph& graph, vertex u, vertex v, vertex w, vector<vertex>& intersection);
//
//void create_triangleList (Graph& ordered_graph, EdgeList2& el, Fcounts& tlist, vector<vertex>& xtl);
//
//void smart_intersection (Graph& orgraph, Graph& gr, int u, int v, vector<vertex>& new_vertex);
//void base_ktruss (Graph& graph, int nEdge, vector<vertex>& T, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//void base_kcore (Graph& graph, int nEdge, vertex *K, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp);
//void base_pitruss (Graph& graph, int nEdge, vertex* T, vertex* R, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp);
//bool connected (Graph& us_graph, vertex u, vertex v);
//bool connected_deg (Graph& graph, Graph& ordered_graph, vertex u, vertex v);
//void intersection (vector<vertex>& vertices, vector<vertex>& edgelist, vector<vertex>& new_vertex);
//void report_stat (set<vertex>& vs, Graph& graph, report& rp);
//void report_stat2 (recai& rc, Graph& graph, report& rp);
//void dump_vertices (set<vertex>& vertices, char* vfile);
//void dump_edges (set<couple>& edges, char* vfile);
//void intersect2 (Graph& graph, vertex u, vertex v, vector<vertex>& intersection);
//vertex create_ordered_graph (Graph& graph, EdgeList2& el, vector<vertex>& xel, Graph& ordered_graph);
//void report_dump (vector<recai>& all_vertices, Graph& graph, char* vfile, bool v_based);
//bool less_than (vertex u, vertex v, Graph& graph);
//void uniquify_vs (vector<int>& ch);
//void uniquify_es (vector<couple>& ch);
//void uniquify_ts (vector<triple>& ch);
//void print_tree (vector<subcore>& hrc, vector<subcore>& cpc, int ind, bool print, set<int>& tolabel, int len, FILE* fp);
//void report_tree (vector<subcore>& hrc, vector<subcore>& cpc, char* fl);
//void report_nested_circle (vector<subcore>& hrc, vector<subcore>& cpc, FILE* fp, string cfl);
//void find_subcore (vertex root, Graph& graph, vertex* K, HashMap<bool>& visited, vector<vertex>& subcore_ids, vector<subcore>& hrc, util::timestamp& mergetime, util::timestamp& looptime, util::timestamp& dftime, stats& op);
//void find_subtruss (vertex root, Graph& graph, vertex* T, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& visited,
//		vector<vertex>& subtruss_ids, vector<subcore>& hrc, stats& op);
//bool intersects13 (vertex v, vertex w, Graph& graph);
//bool intersects14 (vertex v, vertex w, Graph& graph);
//bool sort_rep (report s, report t);
//void report_all_vertex (int variant, Graph& graph, int cn, HashMap<bool>& exists, vector<vertex>& K, char* vfile);
//void report_all_edge (int variant, Graph& graph, int cn, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& exists, vector<vertex>& T, char* vfile);
//void report_all_triangle (int variant, Graph& graph, int cn, unordered_map<triangle_id, int>& tlist,
//		vector<vertex>& xtl, vector<triangle_id>& tl, HashMap<bool>& exists, vector<vertex>& F, char* vfile);
//void report_all_stuff (int variant, Graph& graph, int cn, auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile);
//void report_all_stuff2 (int variant, Graph& graph, int nEdge, int cn, p_auxies& ax, HashMap<bool>& exists, vertex* K, const char* vfile, FILE* ffp);
//void find_sub34 (vertex root, Graph& graph, vector<vertex>& F, unordered_map<triangle_id, int>& tlist,
//		vector<triangle_id>& tl, HashMap<bool>& visited,
//		vector<vertex>& sub34_ids, vector<subcore>& hrc, stats& op);
//void print_time (FILE* fp, const string& str, const timestamp& t);
//void print_time_new (FILE* fp, const string& str, const timestamp& t);
//void print_time_new2 (FILE* fp, const string& str, const timestamp& t, const string& str2);
//void edge_sampling (Graph& graph, double srate);
//void cleanminusones (Graph& graph);
//void cleanminusones2 (Graph& graph, Graph& ordered_graph);
//void wedge_sampling (Graph& ordered_graph, Graph& graph, double srate);
//void count_triangles (Graph& graph, Graph& ordered_graph, vector<vertex>& xel, vertex* TD);
//int count_triangles2 (Graph& graph, Graph& ordered_graph, Graph& TC);
//void compute_densities (int variant, HashMap<int>& K_comps, vector<int>& compid, HashMap<int>& parent, mmap& new_kids,
//		Graph& graph, p_auxies& ax, HashMap<double>& densities,
//		HashMap<int>& sizes, int id, const char* vfile, FILE* ffp, util::timestamp* cdtimes);
//void traverse (int ind, vector<vector<int> >& all_matches, HashMap<bool>& visited, vector<int>& result);
//bool sortBy(const vector<int>& vec1, const vector<int>& vec2);
//void find_my_family (int K, int sc, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct, util::timestamp& p1, util::timestamp& p2, util::timestamp& p3);
////void find_my_family (int K, int sc, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct);
//void compid_set (vector<int>& compid, HashMap<int>& direct);
//void build_hierarchy (int id, int tn, mmap& cmps, HashMap<int>& count_cmps, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct, FILE* fp);
////void build_hierarchy (int j, int tn, mmap& cmps, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct, FILE* fp);
//void feed_forward (HashMap<int>& direct);
//void inverse_anc2kid (int id, mmap& ancestors, mmap& kids, HashMap<bool>& visited, HashMap<int>& direct);
//void handle_matches (vector<vector<vector<int> > >& all_matches, timestamp& t02, timestamp& t03, HashMap<bool>& visited, HashMap<int>& direct);
//void compute_TCP (Graph& graph, vector<vertex>& T, EdgeList2& el, vector<vertex>& xel, Graph& ordered_graph, vector<vector<vector<couple> > >& TCP, FILE* fp);
//void report_by_TCP (int variant, Graph& graph, int tn, tcp& TCP, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& T, const char* vfile, FILE* fp);
//
//void report_results (int variant, int ind, HashMap<int>& ordermap, vector<vertex>& subx_ids, p_auxies& ax, vector<subcore>& hrc, Graph& graph, FILE* fp,
//		timestamp& ph1, timestamp& ph2, timestamp& ph3, timestamp& ph4);
//bool maxcore (vertex root, vertex core_of_root, Graph& graph, vector<vertex>& K, vector<vertex>& vertices, HashMap<bool>& visited);
//bool maxtruss (vertex root, vertex truss_of_root, Graph& graph, vector<vertex>& T, vector<vertex>& xel, EdgeList2& el, vector<couple>& eset, HashMap<bool>& visited);
//bool max34 (vertex root, vertex threefour_of_root, Graph& graph, vector<vertex>& F, unordered_map<triangle_id, int>& tlist,
//		vector<triangle_id>& tl, vector<triple>& tset, HashMap<bool>& visited);
//void report_ZERO (int variant, Graph& graph, int cn, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile, FILE* ffp);
//void report_subgraphs (int variant, vector<subcore>& hrc, vector<int>& subx_ids, Graph& graph, int nEdge, p_auxies& ax, const char* vfile, FILE* ffp);
//void build (int cn, vector<couple>& relations, vector<subcore>& ecr, int* num_subcores, int nEdge, int nVtx);
//void report_LCPS (int variant, Graph& graph, int cn, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile, FILE* ffp);
//
//void old_find_subcore (vertex root, Graph& graph, vector<vertex>& K, HashMap<bool>& visited, FILE* ffp);
//void old_find_subtruss (vertex root, Graph& graph, vector<vertex>& T, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& visited, FILE* ffp);
//void old_find_sub34 (vertex root, Graph& graph, vector<vertex>& F, unordered_map<triangle_id, int>& tlist,
//		vector<triangle_id>& tl, HashMap<bool>& visited, FILE* ffp);
