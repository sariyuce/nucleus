#include "main.h"

// UTILITY FUNCTIONS, MIGHT BE NEEDED

inline vertex old_getEdgeId (vertex u, vertex v, vector<vertex>& xel, EdgeList2& el, Graph& graph) {

	vertex a, b;
	if (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v)) {
		a = u; b = v;
	}
	else {
		a = v; b = u;
	}

	for (vertex i = xel[a]; i < xel[a+1]; i++) {
		if (get<1>(el[i]) == b) {
			return i;
		}
	}
	return -1;
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
		if (marked.hasDefaultValue(s) && s != next_id) { // now you have the root as s
			if (hrc[s].K > tn) { // s becomes kid of next_id-th subcore
				hrc[s].parent = next_id;
				hrc[s].root = next_id;
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
		couple1 cp = make_tuple (v, w);
		bfsorder.push(cp);
		visited[u] = true;
	}
}

void find_subtruss (vertex root, Graph& graph, vertex* T, vector<vertex>& xel, EdgeList2& el, HashMap<bool>& visited,
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
			vertex w = old_getEdgeId(u, inter[j], xel, el, graph);
			vertex x = old_getEdgeId(v, inter[j], xel, el, graph);
			if (T[w] >= tn && T[x] >= tn) {
				if (T[w] == tn) {
					if (visited.hasDefaultValue(w) && subtruss_ids[w] != -1) {
						printf ("VIOLATION-2\n");
						exit(1);
					}

					see_equalK (w, u, inter[j], es, visited, bfsorder);
				}
				else // > tn
					see_largerK (subtruss_ids[w], tn, marked, hrc, merge_sc, op);
				if (T[x] == tn) {
					if (visited.hasDefaultValue(x) && subtruss_ids[x] != -1) {
						printf ("VIOLATION-2\n");
						exit(1);
					}

					see_equalK (x, v, inter[j], es, visited, bfsorder);
				}
				else // > tn
					see_largerK (subtruss_ids[x], tn, marked, hrc, merge_sc, op);
			}
		}
	}

	for (size_t i = 0; i < es.size(); i++)
		subtruss_ids[es[i]] = next_id;

	hrc[next_id].size = es.size();
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

inline void print_Ks (edge nEdge, couple1* el, volatile vertex* T, const char* vfile) {
	string st (vfile);
	st += "_FINAL_H_values";
	FILE* pp = fopen (st.c_str(), "w");
	for (int i = 0; i < nEdge; i++)
		fprintf (pp, "%d %d   %d\n", get<0>(el[i]), get<1>(el[i]), T[i]);
	fclose (pp);
}

inline void print_Hs (Graph& graph, int nEdge, vector<vertex>& xel, EdgeList2& el, volatile vertex* T, const char* vfile, int oc) {
#ifdef CONSTRUCT
	int mT = 0;
	HashMap<bool> exists (false);
	if (oc <= OC_LIMIT) {
#endif
		string st (vfile);
		st += "_H_values_" + to_string (oc);
		FILE* pp = fopen (st.c_str(), "w");
		for (int i = 0; i < el.size(); i++) {
#ifdef CONSTRUCT
			exists[T[i]] = true;
			if (T[i] > mT)
				mT = T[i];
#endif
			fprintf (pp, "%d\n", T[i]);
		}
		fclose (pp);

#ifdef CONSTRUCT

		p_auxies ax;
		ax.xel = &xel;
		ax.el = &el;

		string ss (vfile);
		string out_file = ss + "_HIER_out";
		string aa = ss + "_OC_" + to_string (oc);
		FILE* ffp = fopen (out_file.c_str(), "w");
		report_all_stuff2 (23, graph, nEdge, mT, ax, exists, T, aa.c_str(), ffp);
	}
#endif
}

inline void simple_distance (int nEdge, vertex* Reals, vertex* T, int oc) {
	double score = 0;
	int count = 0;
	for (vertex i = 0; i < nEdge; i++)
		if (T[i] > 0 && Reals[i] > 0) {
//			score += (double) Reals[i] / T[i];
			if (Reals[i] == T[i])
				score++;
			count++;
		}

	score /= count;

	printf ("H %d , similarity: %lf\n", oc, score);
}

// Kendall-Tau
inline void compute_KT (vertex* Reals, EdgeList2& el, Graph& tris, vertex* T, int oc) {
	int count = 0, score = 0;
	for (vertex i = 0; i < el.size(); i++) {
		int u = get<0>(el[i]);
		int v = get<1>(el[i]);
		for (vertex j = 0; j < tris[i].size(); j+=2) {
			int x = get<0>(el[tris[i][j]]);
			int y = get<0>(el[tris[i][j+1]]);
			if ((u == x || u == y) && (v == y || v == x)) {
				TriangleKendallTau (Reals, T, i, tris[i][j], tris[i][j+1], &count, &score);
			}
		}
	}
	printf ("H %d , KendallTau: %lf\n", oc, kt (count, score));
}








// ESSENTIALS

/* store edge-triangle graphs */

void store_count_triangles (vertex* T, vector<vertex>* tris, vertex nVtx, edge nEdge, vertex* adj, edge* xadj, couple1* el, edge* xel, vertex* ordered_adj, edge* ordered_xadj) {

#pragma omp parallel for
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = ordered_xadj[i]; j < ordered_xadj[i+1]; j++) {
			vertex a = ordered_adj[j];
			edge x = j;
			for (edge k = j + 1; k < ordered_xadj[i+1]; k++) {
				edge y = k;
				vertex b = ordered_adj[k];
				vertex v = a, w = b;
				if (lessThan (w, v, xadj))
					swap (v, w);
				edge l = -1;
				for (edge m = ordered_xadj[v]; m < ordered_xadj[v+1]; m++) {
					if (ordered_adj[m] == w) {
						l = m;
						break;
					}
				}
				if (l != -1) {
					edge z = l;
					T[x]++;
					tris[x].push_back (y);
					tris[x].push_back (z);
					T[y]++;
					tris[y].push_back (x);
					tris[y].push_back (z);
				}
			}
		}
	}

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = xadj[i]; j < xadj[i+1]; j++) {
			vertex a = adj[j];
			vertex index = findInd (a, i, ordered_adj, ordered_xadj);
			if (index != -1) {
				edge x = index;
				for (edge m = ordered_xadj[a]; m < ordered_xadj[a+1]; m++) {
					if (lessThan (i, ordered_adj[m], xadj)) {
						edge l = -1;
						for (edge k = ordered_xadj[i]; k < ordered_xadj[i+1]; k++) {
							if (ordered_adj[k] == ordered_adj[m]) {
								l = k;
								break;
							}
						}
						if (l != -1) {
							edge y = l;
							edge z = m;
							T[y]++;
							tris[y].push_back (x);
							tris[y].push_back (z);
						}
					}
				}
			}
		}
	}
}

// this is faster than sort-based computation
inline int mapInitialHI (edge ind, vector<vertex>* tris, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	HashMap<vertex> gmap (0);
	vertex greaters = 0;
	vertex equals = 0;
	vertex H = 0;

	for (vertex ii = 0; ii < tris[ind].size(); ii+=2) {
		edge id1 = tris[ind][ii];
		edge id2 = tris[ind][ii+1];
		vertex sm = (P[id1] <= P[id2] ? P[id1] : P[id2]);

		if (sm == H + 1) {
			if (equals > 0) {
				equals--;
				greaters++;
				gmap[sm]++;
			}
			else { // equals = 0
				H++;
				vertex gH = 0;
				if (!gmap.hasDefaultValue (H))
					gH = gmap[H];
				equals = gH + 1;
				greaters -= gH;
				gmap.erase (H);
			}
		}
		else if (sm > H + 1) {
			if  (equals > 0) {
				equals--;
				greaters++;
				gmap[sm]++;
			}
			else { // equals = 0
				greaters++;
				H++;
				vertex gH = 0;
				if (!gmap.hasDefaultValue (H))
					gH = gmap[H];
				equals = gH;
				greaters -= gH;
				gmap[sm]++;
				gmap.erase (H);
			}
		}
	}

	vertex oP = P[ind];
#ifdef SYNC
	Q[ind] = H;
#else
	P[ind] = H;
#endif

	if (H < oP) {
		return 1;
	}
	else
		return 0;
}

inline int regularUpdateHI (edge ind, vector<vertex>* tris, vertex* T
#ifdef SYNC
		, vertex* U
#endif
) {
	vertex previous_T = T[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_T, 0);
	bool yep_set = false;

	for (vertex j = 0; j < tris[ind].size(); j+=2) { // this is only for sequential
		vertex id1 = tris[ind][j];
		vertex id2 = tris[ind][j+1];
		vertex cur_T = (T[id1] <= T[id2]) ? T[id1] : T[id2];

		if (cur_T >= previous_T)
			greaterequals++;
		else
			smallers[cur_T]++;

		if (greaterequals == previous_T) {
			yep_set = true;
			break;
		}
	}
	if (!yep_set && tris[ind].size() > 0) {
		vertex j;
		for (j = previous_T - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}
#ifdef SYNC
		if (U[ind] != j) {
			U[ind] = j;
			return 1;
		}
#else
		if (T[ind] != j) {
			T[ind] = j;
			return 1;
		}
#endif
	}
	return 0;
}

inline void updateAndNotify (edge ind, vertex* P, int newP, vector<vertex>& neigs, bool* changed) {
	P[ind] = newP;
	changed[ind] = true; // *THIS*
	for (int k = 0; k < neigs.size(); k++) {
		if (P[neigs[k]] >= P[ind]) {
			changed[neigs[k]] = true;
		}
	}
}

inline int efficientUpdateHI (edge ind, vector<vertex>* tris, vertex* T, bool* changed) {

	vertex previous_T = T[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_T, 0);
	bool yep_set = false;

	for (vertex j = 0; j < tris[ind].size(); j+=2) { // this is only for sequential
		vertex id1 = tris[ind][j];
		vertex id2 = tris[ind][j+1];
		vertex cur_T = (T[id1] <= T[id2]) ? T[id1] : T[id2];

		if (cur_T >= previous_T)
			greaterequals++;
		else
			smallers[cur_T]++;

		if (greaterequals == previous_T) {
			yep_set = true;
			break;
		}
	}
	if (!yep_set && tris[ind].size() > 0) {
		vertex j;
		for (j = previous_T - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}

		updateAndNotify (ind, T, j, tris[ind], changed);
		return 1;
	}
	return 0;
}



// stores the edge-triangle graph; basic AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void pitruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	// Ordered (directed) graph creation
	timestamp t_begin;
	couple1* el = (couple1 *) malloc (sizeof(couple1) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));

	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

#ifdef SYNC
	vertex* U = (vertex *) calloc (nEdge, sizeof(vertex));
#endif

	vector<vertex> tris[nEdge];
	store_count_triangles (T, tris, nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);
	free (xel);
	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHI (ind, tris, T
#ifdef SYNC
				, U
#endif
		);
	}

#ifdef SYNC
	memcpy (T, U, sizeof(vertex) * nEdge);
#endif

	timestamp ts5;
	cout << "H 0 time: " << ts5 - t_tc << endl;
	oc++;

	while (flag) {
		flag = false;

		timestamp td1;

		timestamp td2;
		td += td2 - td1;


#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			int fl = regularUpdateHI (ind, tris, T
#ifdef SYNC
					, U
#endif
			);

			if (fl == 1)
				flag = true;
		}

#ifdef SYNC
		memcpy (T, U, sizeof(vertex) * nEdge);
#endif

		timestamp t_end;
		cout << "H " << oc << " time: " << t_end - t_begin - td << endl;
		oc++;
	}


	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

//	print_Ks (nEdge, el, T, vfile);

	free (T);
	free (el);
#ifdef SYNC
	free (U);
#endif

	return;
}

// stores the edge-triangle graph; AND algorithm with the notification mechanism, each vertex is examined to see if active -- no collection in a batch
void NMpitruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {
#ifndef SYNC
	// Ordered (directed) graph creation
	timestamp t_begin;
	couple1* el = (couple1 *) malloc (sizeof(couple1) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));

	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

	vector<vertex> tris[nEdge];
	store_count_triangles (T, tris, nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);
	free (xel);
	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHI (ind, tris, T);
	}
	oc++;

	timestamp ts4;

	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	timestamp ts5;
	cout << "H 0 time: " << ts5 - t_tc << endl;

	while (flag) {
		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (edge i = 0; i < nEdge; i++) {
			edge ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, tris, T, changed);
			if (a == 1)
				flag = true;
		}

		timestamp it3;
		cout << "H " << oc << " time: " << it3 - it1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

//	print_Ks (nEdge, el, T, vfile);

	free (T);
	free (el);

	return;
#else
#endif
}












/* does NOT store edge-triangle graphs, space efficient but slower of course */

inline vertex getEdgeId (vertex u, vertex v, edge* xel, couple1* el, vertex deg_u, vertex deg_v) {
	vertex a = v, b = u;
	if (deg_u < deg_v || (deg_u == deg_v && u < v))
		swap (a, b);
	for (edge i = xel[a]; i < xel[a+1]; i++)
		if (get<1>(el[i]) == b)
			return i;
	return -1;
}

// two loops, single lock
void inefficient_count_triangles (vertex* T, vertex nVtx, vertex* adj, edge* xadj, vertex* ordered_adj, edge* ordered_xadj) { // HERE

#pragma omp parallel for
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = ordered_xadj[i]; j < ordered_xadj[i+1]; j++) {
			vertex a = ordered_adj[j];
			edge x = j;
			for (edge k = j + 1; k < ordered_xadj[i+1]; k++) {
				edge y = k;
				vertex b = ordered_adj[k];
				vertex v = a, w = b;
				if (lessThan (w, v, xadj))
					swap (v, w);
				edge l = -1;
				for (edge m = ordered_xadj[v]; m < ordered_xadj[v+1]; m++) {
					if (ordered_adj[m] == w) {
						l = m;
						break;
					}
				}
				if (l != -1) {
					edge z = l;
					T[x]++;
					T[y]++;
				}
			}
		}
	}

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = xadj[i]; j < xadj[i+1]; j++) {
			vertex a = adj[j];
			vertex index = findInd (a, i, ordered_adj, ordered_xadj);
			if (index != -1) {
				edge x = index;
				for (edge m = ordered_xadj[a]; m < ordered_xadj[a+1]; m++) {
					if (lessThan (i, ordered_adj[m], xadj)) {
						edge l = -1;
						for (edge k = ordered_xadj[i]; k < ordered_xadj[i+1]; k++) {
							if (ordered_adj[k] == ordered_adj[m]) {
								l = k;
								break;
							}
						}
						if (l != -1) {
							edge y = l;
							edge z = m;
							T[y]++;
						}
					}
				}
			}
		}
	}
}

// single loop, three locks
void count_triangles (vertex* T, vertex nVtx, edge* xadj, vertex* ordered_adj, edge* ordered_xadj) {

#pragma omp parallel for
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = ordered_xadj[i]; j < ordered_xadj[i+1]; j++) {
			vertex a = ordered_adj[j];
			edge x = j;
			for (edge k = j + 1; k < ordered_xadj[i+1]; k++) {
				edge y = k;
				vertex b = ordered_adj[k];
				vertex v = a, w = b;
				if (lessThan (w, v, xadj))
					swap (v, w);
				edge l = -1;
				for (edge m = ordered_xadj[v]; m < ordered_xadj[v+1]; m++) {
					if (ordered_adj[m] == w) {
						l = m;
						break;
					}
				}
				if (l != -1) {
					edge z = l;
#pragma omp atomic
					T[x]++;
#pragma omp atomic
					T[y]++;
#pragma omp atomic
					T[z]++;
				}
			}
		}
	}
}

inline int mapInitialHISpaceEfficient (edge ind, edge* xadj, vertex* adj, edge* xel, couple1* el, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	HashMap<vertex> gmap (0);
	vertex greaters = 0;
	vertex equals = 0;
	vertex H = 0;
	vertex u = get<0>(el[ind]);
	vertex v = get<1>(el[ind]);
	vertex deg_u = xadj[u+1] - xadj[u];
	vertex deg_v = xadj[v+1] - xadj[v];

	vertex i = xadj[u], j = xadj[v];
	while (i < xadj[u+1] && j < xadj[v+1]) {
		if (adj[i] < adj[j])
			i++;
		else if (adj[j] < adj[i])
			j++;
		else {
			vertex w = adj[i];
			vertex deg_w = xadj[w+1] - xadj[w];
			vertex id1 = getEdgeId (u, w, xel, el, deg_u, deg_w);
			vertex id2 = getEdgeId (v, w, xel, el, deg_v, deg_w);
			vertex sm = (P[id1] <= P[id2] ? P[id1] : P[id2]);

			if (sm == H + 1) {
				if (equals > 0) {
					equals--;
					greaters++;
					gmap[sm]++;
				}
				else { // equals = 0
					H++;
					vertex gH = 0;
					if (!gmap.hasDefaultValue (H))
						gH = gmap[H];
					equals = gH + 1;
					greaters -= gH;
					gmap.erase (H);
				}
			}
			else if (sm > H + 1) {
				if  (equals > 0) {
					equals--;
					greaters++;
					gmap[sm]++;
				}
				else { // equals = 0
					greaters++;
					H++;
					vertex gH = 0;
					if (!gmap.hasDefaultValue (H))
						gH = gmap[H];
					equals = gH;
					greaters -= gH;
					gmap[sm]++;
					gmap.erase (H);
				}
			}
			i++;
			j++;
		}
	}

	vertex oP = P[ind];
#ifdef SYNC
	Q[ind] = H;
#else
	P[ind] = H;
#endif

	if (H < oP) {
		return 1;
	}
	else
		return 0;
}

inline int regularUpdateHISpaceEfficient (edge ind, edge* xadj, vertex* adj, edge* xel, couple1* el, vertex* T
#ifdef SYNC
		, vertex* U
#endif
) {
	vertex previous_T = T[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_T, 0);
	bool yep_set = false;

	vertex u = get<0>(el[ind]);
	vertex v = get<1>(el[ind]);
	vertex deg_u = xadj[u+1] - xadj[u];
	vertex deg_v = xadj[v+1] - xadj[v];
	vertex w = -1;
	vertex i = xadj[u], j = xadj[v];
	while (i < xadj[u+1] && j < xadj[v+1]) {
		if (adj[i] < adj[j])
			i++;
		else if (adj[j] < adj[i])
			j++;
		else {
			w = adj[i];
			vertex deg_w = xadj[w+1] - xadj[w];
			vertex id1 = getEdgeId (u, w, xel, el, deg_u, deg_w);
			vertex id2 = getEdgeId (v, w, xel, el, deg_v, deg_w);
			vertex cur_T = (T[id1] <= T[id2] ? T[id1] : T[id2]);

			if (cur_T >= previous_T)
				greaterequals++;
			else
				smallers[cur_T]++;

			if (greaterequals == previous_T) {
				yep_set = true;
				break;
			}
			i++;
			j++;
		}
	}

	if (!yep_set && w != -1) {
		vertex j;
		for (j = previous_T - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}
#ifdef SYNC
		if (U[ind] != j) {
			U[ind] = j;
			return 1;
		}
#else
		if (T[ind] != j) {
			T[ind] = j;
			return 1;
		}
#endif
	}
	return 0;
}

inline int efficientUpdateHISpaceEfficient (edge ind, edge* xadj, vertex* adj, edge* xel, couple1* el, vertex* T, bool* changed) {

	vector<edge> neigs;
	vertex previous_T = T[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_T, 0);
	bool yep_set = false;

	vertex u = get<0>(el[ind]);
	vertex v = get<1>(el[ind]);
	vertex deg_u = xadj[u+1] - xadj[u];
	vertex deg_v = xadj[v+1] - xadj[v];
	vertex w = -1;
	vertex i = xadj[u], j = xadj[v];
	while (i < xadj[u+1] && j < xadj[v+1]) {
		if (adj[i] < adj[j])
			i++;
		else if (adj[j] < adj[i])
			j++;
		else {
			w = adj[i];
			vertex deg_w = xadj[w+1] - xadj[w];
			vertex id1 = getEdgeId (u, w, xel, el, deg_u, deg_w);
			vertex id2 = getEdgeId (v, w, xel, el, deg_v, deg_w);
			neigs.push_back (id1);
			neigs.push_back (id2);
			vertex cur_T = (T[id1] <= T[id2]) ? T[id1] : T[id2];

			if (cur_T >= previous_T)
				greaterequals++;
			else
				smallers[cur_T]++;

			if (greaterequals == previous_T) {
				yep_set = true;
				break;
			}
			i++;
			j++;
		}
	}

	if (!yep_set && w != -1) {
		vertex j;
		for (j = previous_T - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}

		updateAndNotify (ind, T, j, neigs, changed);
		return 1;
	}
	return 0;
}



// does NOT store the edge-triangle graph; basic AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void pitrussSpaceEfficient (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	// Ordered (directed) graph creation
	timestamp t_begin;
	couple1* el = (couple1 *) malloc (sizeof(couple1) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));

	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

#ifdef SYNC
	vertex* U = (vertex *) calloc (nEdge, sizeof(vertex));
#endif

//	inefficient_count_triangles (T, nVtx, adj, xadj, ordered_adj, ordered_xadj); // two loops with no lock
	count_triangles (T, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_begin << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHISpaceEfficient (ind, xadj, adj, xel, el, T
#ifdef SYNC
				, U
#endif
		);
	}

#ifdef SYNC
	memcpy (T, U, sizeof(vertex) * nEdge);
#endif

	timestamp ts5;
	cout << "H 0 time: " << ts5 - t_tc << endl;
	oc++;

	while (flag) {
		flag = false;

		timestamp td1;

		timestamp td2;
		td += td2 - td1;


#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			int fl = regularUpdateHISpaceEfficient (ind, xadj, adj, xel, el, T
#ifdef SYNC
					, U
#endif
			);

			if (fl == 1)
				flag = true;
		}

#ifdef SYNC
		memcpy (T, U, sizeof(vertex) * nEdge);
#endif

		timestamp t_end;
		cout << "H " << oc << " time: " << t_end - t_begin - td << endl;
		oc++;
	}


	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

//	print_Ks (nEdge, el, T, vfile);

	free (T);
	free (xel);
	free (el);
#ifdef SYNC
	free (U);
#endif

	return;
}


// does NOT store the edge-triangle graph; AND algorithm with the notification mechanism, each vertex is examined to see if active -- no collection in a batch
void NMpitrussSpaceEfficient (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {
#ifndef SYNC
	// Ordered (directed) graph creation
	timestamp t_begin;
	couple1* el = (couple1 *) malloc (sizeof(couple1) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));

	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "creating ordered graph: " << t_cog - t_begin << endl;
//
//	// Triangle counting for each edge
//	T = (vertex *) calloc (nEdge, sizeof(vertex));

//	inefficient_count_triangles (T, nVtx, adj, xadj, ordered_adj, ordered_xadj); // two loops with no lock
	count_triangles (T, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHISpaceEfficient (ind, xadj, adj, xel, el, T);
	}
	oc++;

	timestamp ts4;

	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	timestamp ts5;
	cout << "H 0 time: " << ts5 - t_tc << endl;

	while (flag) {
		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (edge i = 0; i < nEdge; i++) {
			edge ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHISpaceEfficient (ind, xadj, adj, xel, el, T, changed);
			if (a == 1)
				flag = true;
		}

		timestamp it3;
		cout << "H " << oc << " time: " << it3 - it1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

//	print_Ks (nEdge, el, T, vfile);

//	free (T);
	free (el);

	return;
#else
#endif
	printf ("truss done\n");
}


