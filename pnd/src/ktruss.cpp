#include "main.h"
//
//void count_triangles (Graph& graph, Graph& ordered_graph, vector<vertex>& xel, vertex* TD) {
//
//#pragma omp parallel for
//	for (size_t x = 0; x < ordered_graph.size(); x++) {
//		vertex i = x;
//		for (size_t j = 0; j < ordered_graph[i].size(); j++) {
//			for (size_t k = j + 1; k < ordered_graph[i].size(); k++) {
//				vertex a = ordered_graph[i][j];
//				vertex b = ordered_graph[i][k];
//				if (a == -1 || b == -1)
//					continue;
//
//				vertex aa = a, bb = b;
//				if (smaller (graph, bb, aa) == bb) {
//					aa = b;
//					bb = a;
//				}
//				vertex l = -1;
//				for (size_t ii = 0; ii < ordered_graph[aa].size(); ii++) {
//					if (ordered_graph[aa][ii] == bb) {
//						l = ii;
//						break;
//					}
//				}
//				if (l != -1) {
//					int id1 = xel[i] + j;
//					int id2 = xel[i] + k;
//					int id3 = xel[aa] + l;
//
//#pragma omp atomic
//					TD[id1]++;
//
//#pragma omp atomic
//					TD[id2]++;
//
//#pragma omp atomic
//					TD[id3]++;
//				}
//			}
//		}
//	}
//}

inline vertex getEdgeId (vertex u, vertex v, vector<vertex>& xel, EdgeList2& el, Graph& graph) {
	vertex a, b;
	if (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v)) {
		a = u;
		b = v;
	}
	else {
		a = v;
		b = u;
	}

	for (vertex i = xel[a]; i < xel[a+1]; i++)
		if (get<1>(el[i]) == b)
			return i;

	return -1;
}

//void base_ktruss (Graph& graph, int nEdge, vector<vertex>& T, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp) {


void base_ktruss_StoreTri (Graph& graph, int nEdge, vector<vertex>& T, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp) {
	util::timestamp t_begin;
	char *p = (char*) &totaltime;
	util::timestamp *pout = (util::timestamp*) p;

	// Ordered (directed) graph creation
	EdgeList2 el;
	vector<vertex> xel;
	Graph ordered_graph;
	create_ordered_graph (graph, el, xel, ordered_graph);

	util::timestamp ttl;
	cout << "create ordered graph: " << ttl - t_begin << endl;

	// Triangle counting for each edge
	vertex* TD = (vertex *) malloc (el.size() * sizeof(vertex));
	for (int i = 0; i < el.size(); i++)
		TD[i] = 0;
	vector<vector<vertex>> tris (el.size());
//	count_triangles_adv (graph, ordered_graph, TD, el, xel, tris); todo: fix later

	util::timestamp t_tc;
	cout << "Triangle counting & enumeration: " << t_tc - t_begin << endl;

	// Peeling
	T.resize (el.size(), -1);
	Naive_Bucket na_bs;
	na_bs.Initialize (graph.size(), nEdge);
	vertex id = 0;
	for (size_t i = 0; i < ordered_graph.size(); i++)
		for (size_t j = 0; j < ordered_graph[i].size(); j++)
			na_bs.Insert (id++, TD[xel[i] + j]);
	vertex tc_of_uv = 0;

	while (1) {
		int uv, value_of_uv;
		int ret = na_bs.PopMin(&uv, &value_of_uv);
		if (ret == -1) // if the bucket is empty
			break;

		if (value_of_uv == 0)
			continue;

		tc_of_uv = T[uv] = value_of_uv;
		for (size_t j = 0; j < tris[uv].size(); j+=2) {
			vertex id1 = tris[uv][j];
			vertex id2 = tris[uv][j+1];
			if (T[id1] == -1 && T[id2] == -1) {
				int cur_id1 = na_bs.CurrentValue(id1);
				int cur_id2 = na_bs.CurrentValue(id2);

				if (cur_id1 > tc_of_uv)
					na_bs.DecVal(id1);
				if (cur_id2 > tc_of_uv)
					na_bs.DecVal(id2);
			}
		}
	}

	na_bs.Free();

	util::timestamp td1;
#ifdef DEBUG_00
	string zt (vfile);
	zt += "_K";
	FILE* zp = fopen (zt.c_str(), "w");
	for (vertex i = 0; i < el.size(); i++)
		fprintf (zp, "%d\n", T[i]);
	fclose (zp);
#endif
	cout << "Peeling K time storing tris: " << td1 - t_begin << endl;

	return;
}







//
//void ktruss_Seshlevels (Graph& graph, int nEdge, vector<vertex>& T, vertex* cono, util::timestamp& totaltime, int* truss_number, const char* vfile, FILE* fp) {
//
//	util::timestamp t_begin;
//	char *p = (char*) &totaltime;
//	util::timestamp *pout = (util::timestamp*) p;
//	vertex nVtx = graph.size();
//
//	EdgeList2 el;
//	vector<vertex> xel;
//	Graph ordered_graph;
//	create_ordered_graph (graph, el, xel, ordered_graph);
//
//	util::timestamp tcog;
//	cout << "create ordered graph: " << tcog - t_begin << endl;
//	printf ("el.size: %d\n", el.size());
//	// Triangle counting for each edge
//	vertex* TD = (vertex *) malloc (sizeof (vertex) * el.size());
//	vertex* mark = (vertex *) malloc (sizeof (vertex) * el.size());
//	for (size_t i = 0; i < el.size(); i++) {
//		TD[i] = 0;
//		mark[i] = 0;
//	}
//
//	count_triangles (graph, ordered_graph, xel, TD);
//
//	util::timestamp ttl;
//	cout << "Triangle counting: " << ttl - t_begin << endl;
//
//	int level = 1;
//	while (1) {
//
//		vector<int> torem;
//		for (int i = 0; i < el.size(); i++)
//			if (mark[i] == 0 && TD[i] <= cono[i]) {
//				torem.push_back(i);
////				mark[i] = 1;
//			}
//
//
//		for (int i = 0; i < torem.size(); i++) {
//			int uv = torem[i];
//			vertex u = get<0>(el[uv]);
//			vertex v = get<1>(el[uv]);
//
//			vector<vertex> intersection;
//			intersect2(graph, u, v, intersection);
//			for (size_t j = 0; j < intersection.size(); j++) { /* decrease the tc of the neighbor edges with greater tc */
//				vertex w = intersection[j];
//				vertex id1 = getEdgeId (u, w, xel, el, graph);
//				vertex id2 = getEdgeId (v, w, xel, el, graph);
//
////				if (mark[id1] == 0 && mark[id2] == 0) {
//					if (cono[id1] <= cono[id2])
//						TD[id1]--;
//					if (cono[id2] <= cono[id1])
//						TD[id2]--;
////				}
//			}
//			mark[uv] = 1;
//		}
//
//		if (torem.empty())
//			break;
//		printf ("Level %d   size: %d\n", level, torem.size());
//
//
//
//		level++;
//	}
//
//	for (int i = 0; i < el.size(); i++)
//		if (mark[i] == 0)
//			printf ("K: %d TD: %d\n", cono[i], TD[i]);
//
//
////
////
////	for (int i = 0; i < el.size(); i++) {
////		int uv = i;
////		vertex u = get<0>(el[uv]);
////		vertex v = get<1>(el[uv]);
////
////		vector<vertex> intersection;
////		intersect2(graph, u, v, intersection);
////		for (size_t j = 0; j < intersection.size(); j++) {
////			vertex w = intersection[j];
////			vertex id1 = getEdgeId (u, w, xel, el, graph);
////			vertex id2 = getEdgeId (v, w, xel, el, graph);
////
//////			if (cono[uv] > cono[id1] && LN[uv] < LN[id1])
//////				printf ("%d has K: %d and L: %d, but its neig %d has K: %d and L: %d\n", uv, cono[uv], LN[uv], id1, cono[id1], LN[id1]);
////
//////			if (cono[uv] > cono[id2] && LN[uv] < LN[id2])
//////				printf ("%d has K: %d and L: %d, but its neig %d has K: %d and L: %d\n", uv, cono[uv], LN[uv], id2, cono[id2], LN[id2]);
////
////		}
////	}
//
//
//	return;
//}


//
//
//util::timestamp tcog;
//cout << "create ordered graph: " << tcog - t_begin << endl;
//printf ("el.size: %d\n", el.size());
//// Triangle counting for each edge
//vertex* TD = (vertex *) malloc (sizeof (vertex) * el.size());
//
//int level = 1;
//
//while (1) {
//	for (size_t i = 0; i < el.size(); i++)
//		if (TD[i] != INT_MIN)
//			TD[i] = 0;
//
//	count_triangles (graph, ordered_graph, xel, TD);
//
//	int count = 0;
//	for (int i = 0; i < el.size(); i++)
//		if (TD[i] != INT_MIN && TD[i] <= cono[i]) {
//			count++;
//			vertex u = get<0>(el[i]);
//			vertex v = get<1>(el[i]);
//			int a = u;
//			int b = v;
//			if (v == smaller (graph, u, v)) {
//				a = v;
//				b = u;
//			}
//			for (int j = 0; j < ordered_graph[a].size(); j++)
//				if (ordered_graph[a][j] == b) {
//					ordered_graph[a][j] = -1;
//					break;
//				}
//			TD[i] = INT_MIN;
//		}
//
//	if (count == 0)
//		break;
//	printf ("Level %d   size: %d\n", level, count);
//	level++;
//}
//
//
//
//
//
//
//
//
//
