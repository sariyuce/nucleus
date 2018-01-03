#include "main.h"
////
////void create_classic_triangleList (Graph& ordered_graph, EdgeList2& el, vector<triangle_id>& tl, vector<vertex>& xtl) {
////	xtl.push_back(0);
////	vertex id = 0;
////	for (size_t i = 0; i < el.size(); i++) {
////		vector<vertex> is;
////		intersect2 (ordered_graph, get<0>(el[i]), get<1>(el[i]), is);
////		for (size_t j = 0; j < is.size(); j++) {
////			triangle_id t;
////			t.triple = make_tuple(get<0>(el[i]), get<1>(el[i]), is[j]);
////			t.id = 0;
////			tl.push_back(t);
////		}
////		xtl.push_back(tl.size());
////	}
////	return;
////}
////
//
////
////// parallel
////long intersect3_count_classic_efficient (Graph& graph, Graph& ordered_graph, EdgeList2& el, vector<vertex>& xel, vector<vertex>& xtl, vector<triangle_id>& tl) {
////	long count = 0;
////#pragma omp parallel for
////	for (size_t i = 0; i < el.size(); i++) { // u < v
////		int u = get<0>(el[i]);
////		int v = get<1>(el[i]);
////		vector<vertex> is;
////		intersect_2 (ordered_graph, u, v, is);
////		for (size_t j = 0; j < is.size(); j++) { // u < v < w
////			int wcand = is[j];
////			for (size_t k = j+1; k < is.size(); k++) { // u < v < x
////				int xcand = is[k];
////				int w = wcand, x = xcand;
////				if (x == smaller (graph, w, x)) {
////					int temp = w;
////					w = x;
////					x = temp;
////				}
////				if (find_ind (graph, ordered_graph, w, x) != -1) {
////					vertex in1 = fast_getTriId (u, v, w, xel, xtl, tl, ordered_graph, graph);
////					vertex in2 = fast_getTriId (u, v, x, xel, xtl, tl, ordered_graph, graph);
////					vertex in3 = fast_getTriId (u, w, x, xel, xtl, tl, ordered_graph, graph);
////					vertex in4 = fast_getTriId (v, w, x, xel, xtl, tl, ordered_graph, graph);
////#pragma omp atomic
////					(tl[in1].id)++;
////
////#pragma omp atomic
////					(tl[in2].id)++;
////
////#pragma omp atomic
////					(tl[in3].id)++;
////
////#pragma omp atomic
////					(tl[in4].id)++;
////
////#pragma omp atomic
////					count++;
////				}
////			}
////		}
////	}
////	return count;
////}
////
////bool comptl (triangle_id a, triangle_id b) {
////	if (get<0>(a.triple) == get<0>(b.triple)) {
////		if (get<1>(a.triple) == get<1>(b.triple))
////			return (get<2>(a.triple) < get<2>(b.triple));
////		else
////			return (get<1>(a.triple) < get<1>(b.triple));
////	}
////	else
////		return (get<0>(a.triple) < get<0>(b.triple));
////}
////
////
////inline void print_Ks (vector<triangle_id>& tl, vector<vertex>& F, const char* vfile) {
////	string st (vfile);
////	st += "_K";
////	FILE* pp = fopen (st.c_str(), "w");
////	sort (tl.begin(), tl.end(), comptl);
////	for (vertex i = 0; i < tl.size(); i++) {
////		int u = get<0>(tl[i].triple);
////		int v = get<1>(tl[i].triple);
////		int w = get<2>(tl[i].triple);
////		fprintf (pp, "%d %d %d   %d\n", u, v, w, (F[i] == -1) ? 0 : F[i]);
////	}
////	fclose (pp);
////}
////
//void base_k34 (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp) {
//
//	util::timestamp t_begin;
//	char *p = (char*) &totaltime;
//	util::timestamp *pout = (util::timestamp*) p;
//
//	// Ordered (directed) graph creation
//	EdgeList2 el;
//	vector<vertex> xel;
//	Graph ordered_graph;
//	create_ordered_graph (graph, el, xel, ordered_graph);
//
//	util::timestamp ttl;
//	cout << "create triangle list: " << ttl - t_begin << endl;
//
//	// Enumerate triangles in sorted order
//	vector<triangle_id> tl; // stores triangles
//	vector<vertex> xtl; // end of list indices
//	create_classic_triangleList (ordered_graph, el, tl, xtl);
//
//	// 4-clique counting for each triangle
//	intersect3_count_classic_efficient (graph, ordered_graph, el, xel, xtl, tl);
//
//	util::timestamp t_4c;
//	cout << "4c counting: " << t_4c - t_begin << endl;
//
//	// Peeling
//	F.resize(tl.size(), -1);
//	Naive_Bucket na_bs;
//	int cid = 0;
//	na_bs.Initialize(tl.size(), tl.size());
//	for (size_t i = 0; i < tl.size(); i++) {
//		if (tl[i].id > 0) // this is 4c count
//			na_bs.Insert (cid, tl[i].id);
//		cid++;
//	}
//	vertex fc_of_uvw = 0;
//	int order = 0;
//
//	while (1) {
//		int uvw, value_of_uvw ;
//		int ret = na_bs.PopMin(&uvw, &value_of_uvw);
//		if (ret == -1)
//			break;
//
////		printf ("order %d %d\n", order++, uvw);
//
//		if (value_of_uvw == 0)
//			continue;
//
//		fc_of_uvw = F[uvw] = value_of_uvw;
//		vertex u = get<0> (tl[uvw].triple);
//		vertex v = get<1> (tl[uvw].triple);
//		vertex w = get<2> (tl[uvw].triple);
//
//		vector<vertex> intersection;
//		intersect3_new (graph, u, v, w, intersection);
//		for (size_t j = 0; j < intersection.size(); j++) { // decrease the fc of the neighbor triangles with greater fc
//			vertex x = intersection[j];
//
//			vertex id1 = getTriId (u, v, x, xel, xtl, tl, ordered_graph, graph);
//			vertex id2 = getTriId (u, w, x, xel, xtl, tl, ordered_graph, graph);
//			vertex id3 = getTriId (v, w, x, xel, xtl, tl, ordered_graph, graph);
//			if (F[id1] == -1 && F[id2] == -1 && F[id3] == -1) {
//				int cur_id1 = na_bs.CurrentValue(id1);
//				int cur_id2 = na_bs.CurrentValue(id2);
//				int cur_id3 = na_bs.CurrentValue(id3);
//
//				if (cur_id1 > fc_of_uvw)
//					na_bs.DecVal(id1);
//				if (cur_id2 > fc_of_uvw)
//					na_bs.DecVal(id2);
//				if (cur_id3 > fc_of_uvw)
//					na_bs.DecVal(id3);
//			}
//		}
//	}
//
//	na_bs.Free();
//
//	util::timestamp td1;
//
//	print_Ks (tl, F, vfile);
//	cout << "Peeling K time: " << td1 - t_begin << endl;
//	printf ("max K: %d\n", fc_of_uvw);
//
//	return;
//}
//
//void base_k34_Store4c (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp) {
//	util::timestamp t_begin;
//	char *p = (char*) &totaltime;
//	util::timestamp *pout = (util::timestamp*) p;
//
//	// Ordered (directed) graph creation
//	EdgeList2 el;
//	vector<vertex> xel;
//	Graph ordered_graph;
//	create_ordered_graph (graph, el, xel, ordered_graph);
//
//	util::timestamp ttl;
//	cout << "create triangle list: " << ttl - t_begin << endl;
//
//	// Enumerate triangles in sorted order
//	vector<triangle_id> tl; // stores triangles
//	vector<vertex> xtl; // end of list indices
//	create_classic_triangleList (ordered_graph, el, tl, xtl);
//
//	// 4-clique counting & enumerating for each triangle
//	Graph TF (tl.size()); // triangle's neigs in each 4-clique, final size is 4 * #4-cliques
//	superTriangleList (graph, ordered_graph, xel, xtl, el, tl, TF);
//
//#pragma omp parallel for
//	for (int i = 0; i < TF.size(); i++)
//		tl[i].id = TF[i].size() / 3;
//
//	util::timestamp t_4c;
//	cout << "4c counting & enumeration: " << t_4c - t_begin << endl;
//
//	// Peeling
//	F.resize(tl.size(), -1);
//	Naive_Bucket na_bs;
//	int cid = 0;
//	na_bs.Initialize(tl.size(), tl.size());
//	for (size_t i = 0; i < tl.size(); i++) {
//		if (tl[i].id > 0) // this is 4c count
//			na_bs.Insert (cid, tl[i].id);
//		cid++;
//	}
//	vertex fc_of_uvw = 0;
//
//	while (1) {
//		int uvw, value_of_uvw ;
//		int ret = na_bs.PopMin(&uvw, &value_of_uvw);
//		if (ret == -1)
//			break;
//
//		if (value_of_uvw == 0)
//			continue;
//
//		fc_of_uvw = F[uvw] = value_of_uvw;
//		vertex u = get<0> (tl[uvw].triple);
//		vertex v = get<1> (tl[uvw].triple);
//		vertex w = get<2> (tl[uvw].triple);
//
//		for (int i = 0; i < TF[uvw].size(); i+=3) {
//			int id1 = TF[uvw][i];
//			int id2 = TF[uvw][i+1];
//			int id3 = TF[uvw][i+2];
//
//			if (F[id1] == -1 && F[id2] == -1 && F[id3] == -1) {
//				int cur_id1 = na_bs.CurrentValue(id1);
//				int cur_id2 = na_bs.CurrentValue(id2);
//				int cur_id3 = na_bs.CurrentValue(id3);
//
//				if (cur_id1 > fc_of_uvw)
//					na_bs.DecVal(id1);
//				if (cur_id2 > fc_of_uvw)
//					na_bs.DecVal(id2);
//				if (cur_id3 > fc_of_uvw)
//					na_bs.DecVal(id3);
//			}
//		}
//	}
//
//	na_bs.Free();
//
//	util::timestamp td1;
//#ifdef DEBUG_00
//	print_Ks (tl, F, vfile);
//#endif
//	cout << "Peeling K time storing 4cs: " << td1 - t_begin << endl;
//
//	return;
//}
//
//
//
//void k34_levels (Graph& graph, int nEdge, vector<vertex>& F, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp) {
//	util::timestamp t_begin;
//	char *p = (char*) &totaltime;
//	util::timestamp *pout = (util::timestamp*) p;
//
//	// Ordered (directed) graph creation
//	EdgeList2 el;
//	vector<vertex> xel;
//	Graph ordered_graph;
//	create_ordered_graph (graph, el, xel, ordered_graph);
//
//	util::timestamp ttl;
//	cout << "create triangle list: " << ttl - t_begin << endl;
//
//	// Enumerate triangles in sorted order
//	vector<triangle_id> tl; // stores triangles
//	vector<vertex> xtl; // end of list indices
//	create_classic_triangleList (ordered_graph, el, tl, xtl);
//
//	// 4-clique counting for each triangle
//	intersect3_count_classic_efficient (graph, ordered_graph, el, xel, xtl, tl);
//
//	util::timestamp t_4c;
//	cout << "4c counting: " << t_4c - t_begin << endl;
//
//	// Peeling
//	F.resize(tl.size(), -1);
//	Naive_Bucket na_bs;
//	int cid = 0;
//	na_bs.Initialize(tl.size(), tl.size());
//	for (size_t i = 0; i < tl.size(); i++) {
//		if (tl[i].id > 0) // this is 4c count
//			na_bs.Insert (cid, tl[i].id);
//		cid++;
//	}
//	vertex fc_of_uvw = 0;
//
//	int level = 1;
//	while (1) {
//		int uvw, value_of_uvw ;
//		int ret;
//
//		int val;
//		vector<int> mins;
//		ret = na_bs.PopMin(&uvw, &value_of_uvw);
//		if (ret == -1)
//			break;
//
//		if (value_of_uvw == 0)
//			continue;
//
//		mins.push_back (uvw);
//		while (1) {
//			ret = na_bs.PopMin(&uvw, &val);
//			if (ret == -1) // if the bucket is empty
//				break;
//			if (val > value_of_uvw) {
//				na_bs.Insert (uvw, val);
//				break;
//			}
//			mins.push_back (uvw);
//		}
//
//
//
//		for (int i = 0; i < mins.size(); i++) {
//			int uvw = mins[i];
//
//			F[uvw] = value_of_uvw;
//			vertex u = get<0> (tl[uvw].triple);
//			vertex v = get<1> (tl[uvw].triple);
//			vertex w = get<2> (tl[uvw].triple);
//
//			vector<vertex> intersection;
//			intersect3_new (graph, u, v, w, intersection);
//			for (size_t j = 0; j < intersection.size(); j++) { // decrease the fc of the neighbor triangles with greater fc
//				vertex x = intersection[j];
//
//				vertex id1 = getTriId (u, v, x, xel, xtl, tl, ordered_graph, graph);
//				vertex id2 = getTriId (u, w, x, xel, xtl, tl, ordered_graph, graph);
//				vertex id3 = getTriId (v, w, x, xel, xtl, tl, ordered_graph, graph);
//				if (F[id1] < 0 && F[id2] < 0 && F[id3] < 0) {
//					int cur_id1 = na_bs.CurrentValue(id1);
//					int cur_id2 = na_bs.CurrentValue(id2);
//					int cur_id3 = na_bs.CurrentValue(id3);
//
//					if (cur_id1 > value_of_uvw)
//						na_bs.DecVal(id1);
//					if (cur_id2 > value_of_uvw)
//						na_bs.DecVal(id2);
//					if (cur_id3 > value_of_uvw)
//						na_bs.DecVal(id3);
//				}
//			}
//		}
//		printf ("Level %d K: %d size: %d\n", level, value_of_uvw, mins.size());
//		level++;
//
//	}
//
//
//	na_bs.Free();
//
//	util::timestamp td1;
//#ifdef DEBUG_00
//	print_Ks (tl, F, vfile);
//#endif
//	cout << "Peeling K time: " << td1 - t_begin << endl;
//
//	return;
//}
//
//
//
//
//void k34_Seshlevels (Graph& graph, int nEdge, vector<vertex>& F, vector<vertex>& cono, util::timestamp& totaltime, vertex* threefour_number, const char* vfile, FILE* fp) {
//	util::timestamp t_begin;
//	char *p = (char*) &totaltime;
//	util::timestamp *pout = (util::timestamp*) p;
//
//	// Ordered (directed) graph creation
//	EdgeList2 el;
//	vector<vertex> xel;
//	Graph ordered_graph;
//	create_ordered_graph (graph, el, xel, ordered_graph);
//
//	util::timestamp ttl;
//	cout << "create triangle list: " << ttl - t_begin << endl;
//
//	// Enumerate triangles in sorted order
//	vector<triangle_id> tl; // stores triangles
//	vector<vertex> xtl; // end of list indices
//	create_classic_triangleList (ordered_graph, el, tl, xtl);
//
//	// 4-clique counting for each triangle
//	intersect3_count_classic_efficient (graph, ordered_graph, el, xel, xtl, tl);
//
//	util::timestamp t_4c;
//	cout << "4c counting: " << t_4c - t_begin << endl;
//
//	// Peeling
//	F.resize(tl.size(), -1);
//
//	vector<vertex> FC (tl.size(), 0);
//	for (size_t i = 0; i < tl.size(); i++) {
//		FC[i] = tl[i].id;
//	}
//
//	int level = 1;
//	while (1) {
//		int uvw, value_of_uvw ;
//		int ret;
//
//		vector<int> torem;
//		for (int i = 0; i < tl.size(); i++) {
//			if (FC[i] != -1 && FC[i] <= cono[i]) {
//				torem.push_back(i);
//			}
//		}
//
//		for (int i = 0; i < torem.size(); i++) {
//			int uvw = torem[i];
//
//			vertex u = get<0> (tl[uvw].triple);
//			vertex v = get<1> (tl[uvw].triple);
//			vertex w = get<2> (tl[uvw].triple);
//
//			vector<vertex> intersection;
//			intersect3_new (graph, u, v, w, intersection);
//			for (size_t j = 0; j < intersection.size(); j++) { // decrease the fc of the neighbor triangles with greater fc
//				vertex x = intersection[j];
//
//				vertex id1 = getTriId (u, v, x, xel, xtl, tl, ordered_graph, graph);
//				vertex id2 = getTriId (u, w, x, xel, xtl, tl, ordered_graph, graph);
//				vertex id3 = getTriId (v, w, x, xel, xtl, tl, ordered_graph, graph);
//				if (FC[id1] != -1 && FC[id2] != -1 && FC[id3] != -1) {
//					FC[id1]--;
//					FC[id2]--;
//					FC[id3]--;
//				}
//			}
//			FC[uvw] = -1;
//		}
//
//		if (torem.empty())
//			break;
//		printf ("Level %d   size: %d\n", level, torem.size());
//		level++;
//	}
//
//	util::timestamp td1;
//#ifdef DEBUG_00
//	print_Ks (tl, F, vfile);
//#endif
//	cout << "Peeling K time: " << td1 - t_begin << endl;
//
//	return;
//}
//
