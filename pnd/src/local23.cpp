#include "main.h"

/* STORES EDGE-TRIANGLE GRAPHS */

lol storeCountTriangles (vertex* T, vector<vertex>* tris, vertex nVtx, edge nEdge, vertex* adj, edge* xadj, couple* el, vertex* ordered_adj, edge* ordered_xadj) {
	lol count = 0;
#pragma omp parallel for
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = ordered_xadj[i]; j < ordered_xadj[i+1]; j++) {
			vertex a = ordered_adj[j];
			edge x = j;
			for (edge k = j + 1; k < ordered_xadj[i+1]; k++) {
				edge y = k;
				vertex b = ordered_adj[k];
				vertex v = a, w = b;
				if (isSmaller (xadj, w, v))
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
#ifndef FAST
#pragma omp atomic
					count++;
#endif
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
					if (isSmaller (xadj, i, ordered_adj[m])) {
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
	return count;
}

// this is faster than sort-based computation
inline int mapInitialHI_ST (edge ind, vector<vertex>* tris, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	HashMap<vertex> gmap (-1);
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
				if (!gmap.hasDefaultValue (H)) {
					gH = gmap[H];
					gmap.erase (H);
				}
				equals = gH + 1;
				greaters -= gH;
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
				if (!gmap.hasDefaultValue (H)) {
					gH = gmap[H];
					gmap.erase (H);
				}
				equals = gH;
				greaters -= gH;
				gmap[sm]++;
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

inline int regularUpdateHI_ST (edge ind, vector<vertex>* tris, vertex* T
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

inline int efficientUpdateHI_ST (edge ind, vector<vertex>* tris, vertex* T, bool* changed) {

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
void baseLocal23_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	// Ordered (directed) graph creation
	timestamp t_begin;
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));

	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting and storing for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

#ifdef SYNC
	printf ("It is SYNC\n");
	vertex* U = (vertex *) calloc (nEdge, sizeof(vertex));
#else
	printf ("It is ASYNC\n");
#endif

	vector<vertex> tris[nEdge];
	lol tricount = storeCountTriangles (T, tris, nVtx, nEdge, adj, xadj, el, ordered_adj, ordered_xadj);
	free (xel);
	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif
	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHI_ST (ind, tris, T
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

#ifdef DUMP_Hs
	print_Ks (nEdge, P, vfile, oc);
#endif

	oc++;

	while (flag) {
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			int fl = regularUpdateHI_ST (ind, tris, T
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

#ifdef DUMP_Hs
		print_Ks (nEdge, P, vfile, oc);
#endif

		timestamp t_end;
		cout << "H " << oc << " time: " << t_end - t_begin - td << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif

	free (T);
	free (el);
#ifdef SYNC
	free (U);
#endif

	return;
}

// stores the edge-triangle graph; AND algorithm with the notification mechanism
void nmLocal23_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting and storing for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

	vector<vertex> tris[nEdge];
	lol tricount = storeCountTriangles (T, tris, nVtx, nEdge, adj, xadj, el, ordered_adj, ordered_xadj);
	free (xel);
	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHI_ST (ind, tris, T);
	}
	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_tc - td << endl;
#ifdef DUMP_Hs
	print_Ks (nEdge, T, vfile, oc);
#endif
	oc++;

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI_ST (ind, tris, T, changed);
			if (a == 1)
				flag = true;
		}

		timestamp td2;

#ifdef DUMP_Hs
		print_Ks (nEdge, T, vfile, oc);
#endif
		timestamp td3;
		td += td3 - td2;

		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif

	free (T);
	free (el);
	return;
#endif
}


// Partly parallel k-truss computation (only triangle counting is parallel)
void ktruss_ST (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting and storing for each edge
	edge* tc = (edge *) calloc (nEdge, sizeof(vertex));

	vector<vertex> tris[nEdge];
	lol tricount = storeCountTriangles (tc, tris, nVtx, nEdge, adj, xadj, el, ordered_adj, ordered_xadj);

	timestamp t_tc;
	cout << "Triangle counting & edge-triangle graph construction time: " << t_tc - t_begin << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif

	T = (vertex *) malloc (nEdge * sizeof(vertex));
	memset (T, -1, sizeof(vertex) * nEdge);

	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize (nVtx, nEdge);
	vertex id = 0;
	for (size_t i = 0; i < nVtx; i++) {
		vertex ord_deg = ordered_xadj[i+1] - ordered_xadj[i];
		for (size_t j = 0; j < ord_deg; j++)
			na_bs.Insert (id++, tc[xel[i] + j]);
	}
	free (ordered_adj);
	free (ordered_xadj);

	vertex tc_of_uv = 0;
	while (1) {
		edge uv, val;
		if (na_bs.PopMin(&uv, &val) == -1) // if the bucket is empty
			break;

		if (val == 0)
			continue;

		tc_of_uv = T[uv] = val;

		for (vertex j = 0; j < tris[uv].size(); j+=2) { // this is only for sequential
			vertex id1 = tris[uv][j];
			vertex id2 = tris[uv][j+1];
			if (T[id1] == -1 && T[id2] == -1) {
				if (na_bs.CurrentValue(id1) > tc_of_uv)
					na_bs.DecVal(id1);
				if (na_bs.CurrentValue(id2) > tc_of_uv)
					na_bs.DecVal(id2);
			}
		}
	}

	na_bs.Free();
	cout << "Max truss number: " << tc_of_uv << endl;
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif
	return;
}




/* DOES NOT STORE EDGE_TRIANGLE GRAPHS; SPACE-EFFICIENT (AND SLOWER) */

inline vertex getEdgeId (vertex u, vertex v, edge* xel, couple* el, vertex deg_u, vertex deg_v) {
	vertex a = v, b = u;
	if (deg_u < deg_v || (deg_u == deg_v && u < v))
		swap (a, b);
	for (edge i = xel[a]; i < xel[a+1]; i++)
		if (get<1>(el[i]) == b)
			return i;
	return -1;
}

// single loop, three locks
lol count_triangles (vertex* T, vertex nVtx, edge* xadj, vertex* ordered_adj, edge* ordered_xadj) {
	lol count = 0;
#pragma omp parallel for
	for (vertex i = 0; i < nVtx; i++) {
		for (edge j = ordered_xadj[i]; j < ordered_xadj[i+1]; j++) {
			vertex a = ordered_adj[j];
			edge x = j;
			for (edge k = j + 1; k < ordered_xadj[i+1]; k++) {
				edge y = k;
				vertex b = ordered_adj[k];
				vertex v = a, w = b;
				if (isSmaller (xadj, w, v))
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
#ifndef FAST
#pragma omp atomic
					count++;
#endif
				}
			}
		}
	}
	return count;
}

inline int mapInitialHI (edge ind, edge* xadj, vertex* adj, edge* xel, couple* el, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	HashMap<vertex> gmap (-1);
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
			vertex sm = min (P[id1], P[id2]);

			if (sm == H + 1) {
				if (equals > 0) {
					equals--;
					greaters++;
					gmap[sm]++;
				}
				else { // equals = 0
					H++;
					vertex gH = 0;
					if (!gmap.hasDefaultValue (H)) {
						gH = gmap[H];
						gmap.erase (H);
					}
					equals = gH + 1;
					greaters -= gH;
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
					if (!gmap.hasDefaultValue (H)) {
						gH = gmap[H];
						gmap.erase (H);
					}
					equals = gH;
					greaters -= gH;
					gmap[sm]++;
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

inline int regularUpdateHI (edge ind, edge* xadj, vertex* adj, edge* xel, couple* el, vertex* T
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
			vertex cur_T = min (T[id1], T[id2]);

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
			if (greaterequals >= j)
				break;
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

inline int efficientUpdateHI (edge ind, edge* xadj, vertex* adj, edge* xel, couple* el, vertex* T, bool* changed) {

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



// base AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void baseLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

#ifdef SYNC
	printf ("It is SYNC\n");
	vertex* U = (vertex *) calloc (nEdge, sizeof(vertex));
#else
	printf ("It is ASYNC\n");
#endif

	lol tricount = count_triangles (T, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_cog << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif
	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

#ifdef DEBUG_000
	timestamp ts1;
	int nt = 1, tn = 0;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
		tn = omp_get_thread_num();
	}
	int* ccs = (int *) calloc (nt, sizeof(int));
	timestamp ts2;
	td += ts2 - ts1;
#endif

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
#ifdef DEBUG_000
		ccs[tn] += 	mapInitialHI (ind, xadj, adj, xel, el, T
#ifdef SYNC
				, U
#endif
		);
#else
		mapInitialHI (ind, xadj, adj, xel, el, T
#ifdef SYNC
				, U
#endif
		);
#endif
	}

#ifdef SYNC
	memcpy (T, U, sizeof(vertex) * nEdge);
#endif

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_tc - td << endl;
#ifdef DUMP_Hs
	print_Ks (nEdge, T, vfile, oc);
#endif
	oc++;

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			int fl = regularUpdateHI (ind, xadj, adj, xel, el, T
#ifdef SYNC
					, U
#endif
			);

			if (fl == 1)
				flag = true;
		}

		timestamp td2;
#ifdef SYNC
		memcpy (T, U, sizeof(vertex) * nEdge);
#endif

#ifdef DUMP_Hs
	print_Ks (nEdge, T, vfile, oc);
#endif

		timestamp td3;
		td += td3 - td2;

		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif

	free (T);
	free (xel);
	free (el);
#ifdef SYNC
	free (U);
#endif

	return;
}

// AND algorithm with the notification mechanism
void nmLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (T, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_begin << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

#ifdef DEBUG_000
	timestamp ts1;
	int nt = 1, tn = 0;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
		tn = omp_get_thread_num();
	}
	int* ccs = (int *) calloc (nt, sizeof(int));
	timestamp ts2;
	td += ts2 - ts1;
#endif


#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
#ifdef DEBUG_000
		ccs[tn] += 	mapInitialHI (ind, xadj, adj, xel, el, T);
#else
		mapInitialHI (ind, xadj, adj, xel, el, T);
#endif
	}
	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_tc - td << endl;
#ifdef DUMP_Hs
	print_Ks (nEdge, T, vfile, oc);
#endif
	oc++;

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, xadj, adj, xel, el, T, changed);
#ifdef DEBUG_000
			ccs[tn] += a;
#endif
			if (a == 1)
				flag = true;
		}

		timestamp td2;
#ifdef DEBUG_000
		timestamp a1;
		int degisenler = 0;
		for(int i = 0; i < nt; i++)
			degisenler += ccs[i];
		memset (ccs, 0, nt * sizeof(int));
		cout << "CHANGEDS: " << degisenler << endl;
		timestamp a2;
		td += a2 - a1;
#endif

#ifdef DUMP_Hs
		print_Ks (nEdge, T, vfile, oc);
#endif
		timestamp td3;
		td += td3 - td2;

		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif

	free (T);
	free (el);
	return;
#endif
}


// Partly parallel k-truss computation (only triangle counting is parallel)
void ktruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	edge* tc = (edge *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (tc, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_cog << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif
	T = (vertex *) malloc (nEdge * sizeof(vertex));
	memset (T, -1, sizeof(vertex) * nEdge);

	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize (nVtx, nEdge);
	vertex id = 0;
	for (size_t i = 0; i < nVtx; i++) {
		vertex ord_deg = ordered_xadj[i+1] - ordered_xadj[i];
		for (size_t j = 0; j < ord_deg; j++)
			na_bs.Insert (id++, tc[xel[i] + j]);
	}
	free (ordered_adj);
	free (ordered_xadj);

	vertex tc_of_uv = 0;
	while (1) {
		edge uv, val;
		if (na_bs.PopMin(&uv, &val) == -1) // if the bucket is empty
			break;

		if (val == 0)
			continue;

		tc_of_uv = T[uv] = val;

		vertex u = get<0>(el[uv]);
		vertex v = get<1>(el[uv]);
		vector<vertex> intersection;

		intersection2 (adj, xadj, u, v, intersection);
		for (auto w : intersection) { /* decrease the tc of the neighbor edges with greater tc */
			vertex id1 = getEdgeId (u, w, xel, el, xadj[u+1] - xadj[u], xadj[w+1] - xadj[w]);
			vertex id2 = getEdgeId (v, w, xel, el, xadj[v+1] - xadj[v], xadj[w+1] - xadj[w]);
			if (T[id1] == -1 && T[id2] == -1) {
				if (na_bs.CurrentValue(id1) > tc_of_uv)
					na_bs.DecVal(id1);
				if (na_bs.CurrentValue(id2) > tc_of_uv)
					na_bs.DecVal(id2);
			}
		}
	}

	na_bs.Free();
	cout << "Max truss number: " << tc_of_uv << endl;
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (nEdge, T, vfile);
#endif
	return;
}







// ONGOING WORK

// Finds the max K value by iterating on top-number high degree vertices
void fast23DegeneracyNumber (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, vertex topK) {

	int number = topK; // only top-number vertices in degree are executed
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	P = (edge *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (P, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_cog << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	vector<eda> KK;
	for (vertex i = 0; i < nEdge; i++)
		KK.push_back (make_tuple(i, P[i]));

	sort (KK.begin(), KK.end(), kksort);

	edge start = 0;
	for (edge i = start; i < number; i++)
		changed[get<0>(KK[i])] = true;


	while (flag) {

		memset (changed, 0, sizeof(bool) * nEdge); // set all false
		for (edge i = start; i < number; i++)
			changed[get<0>(KK[i])] = true; // only top-k true

		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, xadj, adj, xel, el, P, changed);
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


	vertex maxT = 0;
	for (edge i = start; i < number; i++) {
		int curP = P[get<0>(KK[i])];
		if (curP > maxT)
			maxT = curP;
//		printf ("%d goes from %d to %d\n", get<0>(KK[i]), get<1>(KK[i]), curP);
	}

	printf ("max T: %d\n", maxT);
#endif
}


/*
- L_1: All those edges with the minimum triangle count.
- To find L_i: delete all edges in L_j (j < i). Then find all edges with the minimum triangle count in the remaining graph.
*/
void ktruss_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* L, const char* vfile) {

	timestamp t_begin;
	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	edge* tc = (edge *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (tc, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_cog << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif
	L = (vertex *) malloc (nEdge * sizeof(vertex));
	memset (L, -1, sizeof(vertex) * nEdge);

	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize (nVtx, nEdge);
	vertex id = 0;
	for (size_t i = 0; i < nVtx; i++) {
		vertex ord_deg = ordered_xadj[i+1] - ordered_xadj[i];
		for (size_t j = 0; j < ord_deg; j++)
			na_bs.Insert (id++, tc[xel[i] + j]);
	}
	free (ordered_adj);
	free (ordered_xadj);

	vertex tc_of_uv = 0;

	int level = 1;
	while (1) {
		edge uv, min_val, val;
		vector<vertex> mins;

		if (na_bs.PopMin(&uv, &min_val) == -1) // if the bucket is empty
			break;

		if (min_val == 0)
			continue;

		mins.push_back (uv);

		while (1) {
			if (na_bs.PopMin(&uv, &val) == -1) // if the bucket is empty
				break;
			if (val > min_val) {
				na_bs.Insert (uv, val);
				break;
			}
			mins.push_back (uv);
		}

		for (auto uv : mins) {
			L[uv] = min_val;
			vertex u = get<0>(el[uv]);
			vertex v = get<1>(el[uv]);

			vector<vertex> intersection;
			intersection2 (adj, xadj, u, v, intersection);
			for (auto w : intersection) { /* decrease the tc of the neighbor edges with greater tc */
				vertex id1 = getEdgeId (u, w, xel, el, xadj[u+1] - xadj[u], xadj[w+1] - xadj[w]);
				vertex id2 = getEdgeId (v, w, xel, el, xadj[v+1] - xadj[v], xadj[w+1] - xadj[w]);
				if (L[id1] == -1 && L[id2] == -1) {
					if (na_bs.CurrentValue(id1) > tc_of_uv)
						na_bs.DecVal(id1);
					if (na_bs.CurrentValue(id2) > tc_of_uv)
						na_bs.DecVal(id2);
				}
			}
		}

		printf ("Level %d K: %d size: %d\n", level, min_val, mins.size());
		level++;
	}
	na_bs.Free();
#ifdef DUMP_K
	print_Ks (nEdge, L, vfile);
#endif
	return;
}


/*
- L_1: All those edges whose triangle count is at most T-value.
- To find L_i: delete all edges in L_j (j < i). Then find all edges whose triangle counts in the remaining graph is at most (original) T-value.
*/
void ktruss_Sesh_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, vertex* L, const char* vfile) {

	timestamp t_begin;
	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Triangle counting for each edge
	edge* tc = (edge *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (tc, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	timestamp t_tc;
	cout << "Triangle counting time: " << t_tc - t_cog << endl;
#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif
	L = (vertex *) malloc (nEdge * sizeof(vertex));

	int level = 1;
	while (1) {
		vector<edge> atMostT;
		for (vertex i = 0; i < nEdge; i++) {
			if (tc[i] != -1 && tc[i] <= T[i]) {
//				printf ("for %d   tc: %d    T: %d\n", i, tc[i], T[i]);
				atMostT.push_back(i);
			}
		}

		if (atMostT.empty())
			break;

		for (edge uv : atMostT) {
			vertex u = get<0>(el[uv]);
			vertex v = get<1>(el[uv]);

			vector<vertex> intersection;
			intersection2 (adj, xadj, u, v, intersection);
			for (auto w : intersection) { /* decrease the tc of the neighbor edges with greater tc */
				vertex id1 = getEdgeId (u, w, xel, el, xadj[u+1] - xadj[u], xadj[w+1] - xadj[w]);
				vertex id2 = getEdgeId (v, w, xel, el, xadj[v+1] - xadj[v], xadj[w+1] - xadj[w]);
				if (tc[id1] > 0 && tc[id2] > 0) {
						tc[id2]--;
						tc[id1]--;
				}
			}
			tc[uv] = -1;
			L[uv] = level;
		}

		printf ("Level %d size: %d\n", level, atMostT.size());
		level++;
	}
#ifdef DUMP_K
	print_Ks (nEdge, L, vfile);
#endif
	return;
}




