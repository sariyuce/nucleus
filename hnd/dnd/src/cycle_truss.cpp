#include "main.h"

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

inline vertex smaller (Graph& graph, vertex u, vertex v) {
	if (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v))
		return u;
	return v;
}

vertex count_cycles (Graph& dgraph, Graph& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		for (vertex j = 1; j < dgraph[u].size(); j++) {
			vertex v = dgraph[u][j];
			if (v < 0)
				break;

			// another option is that you can check incoming edges of u
			if (exists (u, dgraph[v])) // u-v must be directed because no edges in cycle is reciprocal
				continue;

			vector<vertex> ints;
			inter (1, 2, dgraph, u, v, ints); // todo: two items written to ints although ints[k] is not used. because inter is generic, can be fixed later

			for (vertex k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[v][ints[k+1]]; // equal to M2P (dgraph[u][k])
				TC[u][j]++; // u-v
				TC[v][ints[k+1]]++; // v-w
				TC[w][ind (u, dgraph[w])]++; // w-u
				count++;
			}
		}
	}

	for (vertex i = 0; i < TC.size(); i++)
		for (vertex j = 1; j < TC[i].size(); j++)
			TC[i][j] /= 3; // 3 same type edges in a cycle
	count /= 3;

	printf ("total cycle count: %d\n", count);
	return count;
}


void cycle_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	// Cycle counting for each edge
	vertex nVtx = graph.size();
	Graph TC;
	TC.resize(nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (graph[i].size(), 0);
	lol tric = count_cycles (graph, TC);

	fprintf (fp, "# cycles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Cycle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (graph.size(), nEdge);
	vertex id = 0;
	vector<vp> el;
	vector<vertex> xel;
	xel.push_back(0);
	// each non-reciprocal edge and its cycle-count is inserted to bucket
	for (size_t i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings (graph[i], ret);
		for (size_t r = 0; r < ret.size(); r++) {
			vp c (i, graph[i][ret[r]]);
			el.push_back(c); // first is source, second is target
			printf ("cycle count of %d is %d\n", id, TC[i][ret[r]]);
			if (TC[i][ret[r]] > 0)
				nBucket.Insert (id++, TC[i][ret[r]]);
			else
				K[id++]= 0;
		}
		xel.push_back(el.size());
	}

//	for (auto i = 0; i < xel.size(); i++) {
//		printf ("xel[%d]: %d\n", i, xel[i]);
//		for (auto j = xel[i]; j < xel[i+1]; j++)
//			printf ("%d -> %d\n", el[j].first, el[j].second);
//	}

	vertex tc_e = 0;

	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;

		tc_e = K[e] = val;

		vertex u = el[e].first; // source
		vertex v = el[e].second; // target

		vector<vertex> ints;
		inter (1, 2, graph, u, v, ints); // todo: we don't need ints[k+1]s

		for (auto k = 0; k < ints.size(); k+=2) {
			vertex w = graph[v][ints[k+1]]; // equal to M2P (graph[u][k])
			vertex id1 = getEdgeId (w, u, xel, el, graph);
			vertex id2 = getEdgeId (v, w, xel, el, graph);
			if (K[id1] == -1 && K[id2] == -1) {
				if (nBucket.CurrentValue(id1) > tc_e)
					nBucket.DecVal(id1);
				if (nBucket.CurrentValue(id2) > tc_e)
					nBucket.DecVal(id2);
			}
		}
	}

	nBucket.Free();
	*maxK = tc_e;

	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}

	for (auto i = 0; i < K.size(); i++)
		printf ("truss of %d is %d\n", i, K[i]);
	return;
}

/*
int count_triangles2 (Graph& graph, Graph& ordered_graph, Graph& TC) {

	int count = 0;
	for (size_t x = 0; x < ordered_graph.size(); x++) {
		vertex i = x;
		for (size_t j = 0; j < ordered_graph[i].size(); j++) {
			for (size_t k = j + 1; k < ordered_graph[i].size(); k++) {
				vertex a = ordered_graph[i][j];
				vertex b = ordered_graph[i][k];
				if (a != -1 && b != -1) {
					vertex aa = a, bb = b;
					if (smaller (graph, bb, aa) == bb) {
						aa = b;
						bb = a;
					}
					vertex l = -1;
					for (size_t ii = 0; ii < ordered_graph[aa].size(); ii++) {
						if (ordered_graph[aa][ii] == bb) {
							l = ii;
							break;
						}
					}
					if (l != -1) {
						TC[i][j]++;
						TC[i][k]++;
						TC[aa][l]++;
						count++;
					}
				}
			}
		}
	}
	return count;
}


inline tuple<int, int> my_make_tuple (vertex a, vertex b) {
	return make_tuple (a<b?a:b, a<b?b:a);
}

inline void see_largerK (int s, int tn, HashMap<bool>& marked, vector<subcore>& hrc, vector<int>& merge_sc, stats& op) {
	if (marked.hasDefaultValue(s)) {
		vector<int> acc; // accumulate during find
		int next_id = hrc.size() - 1;
		int old_s = s;

		while (hrc[s].root != -1) {
			acc.push_back (s);
			s = hrc[s].root;
		}
		for (int i : acc)
			hrc[i].root = s;

//		while (hrc[s].root != -1) {
//			int r = hrc[s].root;
//			if (hrc[r].root != -1)
//				hrc[s].root = hrc[r].root;
//			s = r;
//		}
//		op.find_op++;

		if (marked.hasDefaultValue(s) && s != next_id) { // now you have the root as s
			if (hrc[s].K > tn) { // s becomes kid of next_id-th subcore
				hrc[s].parent = next_id;
//				hrc[next_id].children.push_back (s);
//				hrc[next_id].rank += hrc[s].rank;
				hrc[s].root = next_id;
//				op.adjust_op += hrc[s].children.size();
//				for (int ch : hrc[s].children)
//					hrc[ch].root = next_id;
			}
			else // hrc[s].K == tn
				merge_sc.push_back (s); // defer the merge processing to the end: need to know the ranks

			marked[s] = true;
		}
		marked[old_s] = true;
	}
}

inline void see_equalK (vertex u, vertex v, vertex w, vector<int>& es, HashMap<bool>& visited, queue<couple1>& bfsorder) {
	es.push_back (u);
	if (visited.hasDefaultValue(u)) {
		couple1 cp = my_make_tuple (v, w);
		bfsorder.push(cp);
		visited[u] = true;
	}
}


void find_subtruss (vertex root, Graph& graph, vector<vertex>& T, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& visited,
		vector<vertex>& subtruss_ids, vector<subcore>& hrc, stats& op) {
	vertex tn = T[root];
	int next_id = hrc.size();
	subcore sc;
	sc.rank = 0;
	sc.K = tn;
	sc.parent = -1;
	sc.root = -1;
	sc.visible = true;
	hrc.push_back (sc);
	vector<int> es; // put the edge ids in subtruss
	vector<int> merge_sc;
	HashMap<bool> marked (false); // to keep track of the neighbor subcores (and their ancestors)
	queue<couple1> bfsorder;
	bfsorder.push(el[root]);
	es.push_back (root);
	visited[root] = true;
	while (!bfsorder.empty()) {
		couple1 c = bfsorder.front();
		bfsorder.pop();
		vertex u = get<0>(c);
		vertex v = get<1>(c);
		vector<vertex> inter;
		intersection (graph[u], graph[v], inter);
		for (size_t j = 0; j < inter.size(); j++) {
			vertex w = getEdgeId(u, inter[j], xel, el, graph);
			vertex x = getEdgeId(v, inter[j], xel, el, graph);
			if (T[w] >= tn && T[x] >= tn) {
				if (T[w] == tn) {
//					if (visited.hasDefaultValue(w) && subtruss_ids[w] != -1) {
//						printf ("VIOLATION-2\n");
//						exit(1);
//					}

					see_equalK (w, u, inter[j], es, visited, bfsorder);
				}
				else // > tn
					see_largerK (subtruss_ids[w], tn, marked, hrc, merge_sc, op);
				if (T[x] == tn) {
//					if (visited.hasDefaultValue(x) && subtruss_ids[x] != -1) {
//						printf ("VIOLATION-2\n");
//						exit(1);
//					}

					see_equalK (x, v, inter[j], es, visited, bfsorder);
				}
				else // > tn
					see_largerK (subtruss_ids[x], tn, marked, hrc, merge_sc, op);
			}
		}
	}

	for (size_t i = 0; i < es.size(); i++)
		subtruss_ids[es[i]] = next_id;

#ifndef SHORT
	hrc[next_id].size = es.size();
#endif
	op.total_size += es.size();

	if (!merge_sc.empty()) {
		merge_sc.push_back (next_id);
		vector<bool> seen (merge_sc.size(), false);
		for (int i = 0; i < merge_sc.size(); i++) {
			if (!seen[i]) {
				for (int j = i + 1; j < merge_sc.size(); j++) {
					if (!seen[i] && !seen[j]) {
						op.union_op++;
						int ch = merge_sc[i];
						int pr = merge_sc[j];
						bool fl = false;
						if (hrc[ch].rank > hrc[pr].rank) {
							int t = pr;
							pr = ch;
							ch = t;
							seen[j] = true;
						}
						else {
							seen[i] = true;
							fl = true;
						}
						hrc[ch].parent = pr;
						hrc[ch].root = pr;
//						op.adjust_op += hrc[ch].children.size();
//						for (int kid : hrc[ch].children)
//							hrc[kid].root = pr;
//						hrc[pr].children.push_back (ch);
//						hrc[pr].rank += hrc[ch].rank;
						if (hrc[pr].rank == hrc[ch].rank)
							hrc[pr].rank++;

						if (fl)
							break;
					}
				}
			}
		}
	}
}

bool maxtruss (vertex root, vertex truss_of_root, Graph& graph, vector<vertex>& T,
		vector<vertex>& xel, EdgeList2& el, vector<couple1>& eset, HashMap<bool>& visited) {

	queue<couple1> bfsorder;
	visited[root] = true;
	bfsorder.push(el[root]);
	bool maximal = true;

	while (!bfsorder.empty()) {

		couple1 c = bfsorder.front();
		bfsorder.pop();
		vertex u = get<0>(c);
		vertex v = get<1>(c);
//		rc.vs.insert(u);
//		rc.vs.insert(v);
//		rc.edges.insert(c);
		eset.push_back(c);

		vector<vertex> inter;
		intersection (graph[u], graph[v], inter);
		vector<vertex> vs;
		vs.push_back(u);
		vs.push_back(v);

		for (size_t j = 0; j < inter.size(); j++) {
			vertex w = getEdgeId(vs[0], inter[j], xel, el, graph);
			vertex x = getEdgeId(vs[1], inter[j], xel, el, graph);
			if (T[w] >= truss_of_root && T[x] >= truss_of_root) {
				if (visited.hasDefaultValue(w)) {
					couple1 cp = my_make_tuple (vs[0], inter[j]);
					bfsorder.push(cp);
					visited[w] = true;
				}
				if (visited.hasDefaultValue(x)) {
					couple1 cp = my_make_tuple (vs[1], inter[j]);
					bfsorder.push(cp);
					visited[x] = true;
				}
				if (T[w] > truss_of_root || T[x] > truss_of_root)
					maximal = false;
			}
		}
	}

	return maximal;
}


inline void anc_up (int co, HashMap<bool>& seen, vector<int>& anc) {
	if (seen.hasDefaultValue (co)) {
		seen[co] = true;
		anc.push_back (co);
	}
}


inline void match_up (int co, int* min_compid, HashMap<bool>& seen, vector<int>& matches) {
	if (seen.hasDefaultValue (co)) {
		seen[co] = true;
		matches.push_back (co);
		if (co < *min_compid)
			*min_compid = co;
	}
}



//just traverses edge-triangle fashion. experimental purposes.

void traversal_truss (Graph& graph, vector<vertex>& xel, EdgeList2& el) {
	// printf ("traversal truss, el.size(): %d\n", el.size());
	HashMap<bool> visited (false);

	int cn = 0;
//	for (int i = 0; i < el.size(); i++) {
//		couple1 c = el[i];
//		vertex u = get<0>(c);
//		vertex v = get<1>(c);
//		printf ("u: %d, v: %d\n", u, v);
//		if (u > graph.size() || v > graph.size()) {
//			printf ("NOLUYO!\n");
//			exit(1);
//		}
//	}
//	exit(1);
	for (int i = 0; i < el.size(); i++) {
//		int sm = smaller (graph, get<0>(el[i]), get<1>(el[i]));
//		int bg = (sm == get<0>(el[i])) ? get<1>(el[i]) : get<0>(el[i]);
//		vertex ii = getEdgeId(sm, bg, xel, el, graph);

		if (visited.hasDefaultValue(i)) {
			visited[i] = true;
			queue<couple1> bfsorder;
//			couple1 cp = my_make_tuple (sm, bg);
//			bfsorder.push(cp);
			bfsorder.push(el[i]);

//			printf ("init u: %d, v: %d\n", get<0>(el[i]), get<1>(el[i]));
			while (!bfsorder.empty()) {
				couple1 c = bfsorder.front();
				bfsorder.pop();
				cn++;
				vertex u = get<0>(c);
				vertex v = get<1>(c);
//				printf ("u: %d, v: %d\n", u, v);
				vector<vertex> inter;
				intersection (graph[u], graph[v], inter);
				for (size_t j = 0; j < inter.size(); j++) {
					vertex w = getEdgeId(u, inter[j], xel, el, graph);
					vertex x = getEdgeId(v, inter[j], xel, el, graph);
					if (visited.hasDefaultValue(w)) {
						visited[w] = true;
						couple1 cp = my_make_tuple (u, inter[j]);
						bfsorder.push(cp);
					}
					if (visited.hasDefaultValue(x)) {
						visited[x] = true;
						couple1 cp = my_make_tuple (v, inter[j]);
						bfsorder.push(cp);
					}
				}
			}
		}
	}
	// printf ("cn: %d\n", cn);
}





void old_find_subtruss (vertex root, Graph& graph, vector<vertex>& T, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& visited, FILE* ffp) {
	int cnt = 0;
	vertex truss_of_root = T[root];
	queue<couple1> bfsorder;
	bfsorder.push(el[root]);
	visited[root] = true;

	while (!bfsorder.empty()) {
		couple1 c = bfsorder.front();
		bfsorder.pop();
		cnt++;
		vertex u = get<0>(c);
		vertex v = get<1>(c);
		vector<vertex> inter;
		intersection (graph[u], graph[v], inter);
		for (size_t j = 0; j < inter.size(); j++) {
			vertex w = getEdgeId(u, inter[j], xel, el, graph);
			vertex x = getEdgeId(v, inter[j], xel, el, graph);
			if (T[w] >= truss_of_root && T[x] >= truss_of_root) {
				if (T[w] == truss_of_root && visited.hasDefaultValue(w)) {
					couple1 cp = my_make_tuple (u, inter[j]);
					bfsorder.push(cp);
					visited[w] = true;
				}
				if (T[x] == truss_of_root && visited.hasDefaultValue(x)) {
					couple1 cp = my_make_tuple (v, inter[j]);
					bfsorder.push(cp);
					visited[x] = true;
				}
			}
		}
	}
	fprintf (ffp, "T: %d size: %d\n", truss_of_root, cnt);
}

*/
