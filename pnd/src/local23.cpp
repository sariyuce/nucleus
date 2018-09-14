#include "main.h"

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

// this is faster than sort-based computation
inline int mapInitialHI (edge ind, edge* xadj, vertex* adj, edge* xel, couple* el, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	unordered_map<vertex, vertex> gmap;
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
					auto it = gmap.find(H);
					if (it != gmap.end()) {
						gH = it->second;
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
					auto it = gmap.find(H);
					if (it != gmap.end()) {
						gH = it->second;
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

inline void updateAndNotify (edge ind, vertex* P, int newP, vector<vertex>& neigs, bool* changed) {
	P[ind] = newP;
	changed[ind] = true;
	for (int k = 0; k < neigs.size(); k++) {
		if (P[neigs[k]] >= P[ind]) {
			changed[neigs[k]] = true;
		}
	}
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

	const auto t_begin = chrono::steady_clock::now();
	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	const auto t_cog = chrono::steady_clock::now();
	tms t1 = t_cog - t_begin;
	printf ("Creating ordered graph: %.6lf secs\n", t1.count());

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
	const auto t_tc = chrono::steady_clock::now();
	tms t2 = t_tc - t_cog;
	printf ("Triangle counting time: %.6lf secs\n", t2.count());

#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif

	int oc = 0;
	bool flag = true;

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		mapInitialHI (ind, xadj, adj, xel, el, T
#ifdef SYNC
				, U
#endif
		);
	}

#ifdef SYNC
	memcpy (T, U, sizeof(vertex) * nEdge);
#endif

	const auto t_init = chrono::steady_clock::now();
	tms t3 = t_init - t_tc;
	printf ("H %d time: %.6lf secs\n", oc, t3.count());
	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
	const auto ts1 = chrono::steady_clock::now();
	print_Ks (nEdge, T, vfile, oc);
	const auto ts2 = chrono::steady_clock::now();
	td += ts2 - ts1;
#endif
	oc++;

	while (flag) {
		const auto td1 = chrono::steady_clock::now();
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

#ifdef SYNC
		memcpy (T, U, sizeof(vertex) * nEdge);
#endif

		const auto td2 = chrono::steady_clock::now();
#ifdef DUMP_Hs
		const auto ts1 = chrono::steady_clock::now();
		print_Ks (nEdge, T, vfile, oc);
		const auto ts2 = chrono::steady_clock::now();
		td += ts2 - ts1;
#endif

		tms step = td2 - td1;
		printf ("H %d time: %.6lf secs\n", oc, step.count());
		oc++;
	}

#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nEdge, T, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

	free (T);
	free (xel);
	free (el);
#ifdef SYNC
	free (U);
#endif
	printf ("Converges at %d\n", oc);
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin - td;
	printf ("Total time: %.6lf secs\n", total.count());

	return;
}

// AND algorithm with the notification mechanism
void nmLocal23 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else

	const auto t_begin = chrono::steady_clock::now();
	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	const auto t_cog = chrono::steady_clock::now();
	tms t1 = t_cog - t_begin;
	printf ("Creating ordered graph: %.6lf secs\n", t1.count());

	// Triangle counting for each edge
	T = (vertex *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (T, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	free (ordered_adj);
	free (ordered_xadj);
	const auto t_tc = chrono::steady_clock::now();
	tms t2 = t_tc - t_cog;
	printf ("Triangle counting time: %.6lf secs\n", t2.count());

#ifndef FAST
	cout << "# triangles: " << tricount << endl;
#endif

	int oc = 0;
	bool flag = true;

	int nt;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
	}

	int counters[nt];
	for (int i = 0; i < nt; i++)
		counters[i] = 0;

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++) {
		counters[omp_get_thread_num()]++;
		mapInitialHI (ind, xadj, adj, xel, el, T);
	}
	bool changed[nEdge];
	memset (changed, 255, sizeof(bool) * nEdge); // set all true

	const auto t_init = chrono::steady_clock::now();;
	tms t3 = t_init - t_tc	;
	printf ("H %d time: %.6lf secs\n", oc, t3.count());
	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
	const auto ts1 = chrono::steady_clock::now();
	print_Ks (nEdge, T, vfile, oc);
	const auto ts2 = chrono::steady_clock::now();
	td += ts2 - ts1;
#endif
	oc++;

	while (flag) {
		const auto td1 = chrono::steady_clock::now();
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			if (!changed[ind])
				continue;
			counters[omp_get_thread_num()]++;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, xadj, adj, xel, el, T, changed);
			if (a == 1)
				flag = true;
		}

		const auto td2 = chrono::steady_clock::now();
#ifdef DUMP_Hs
		const auto ts1 = chrono::steady_clock::now();
		print_Ks (nEdge, T, vfile, oc);
		const auto ts2 = chrono::steady_clock::now();
		td += ts2 - ts1;
#endif

		tms step = td2 - td1;
		printf ("H %d time: %.6lf secs\n", oc, step.count());
		oc++;
	}

#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nEdge, T, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

	free (T);
	free (el);
	printf ("Converges at %d\n", oc);
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin - td;
	printf ("Total time: %.6lf secs\n", total.count());
	return;
#endif
}


// Partly parallel k-truss computation (only triangle counting is parallel)
void ktruss (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* T, const char* vfile) {

	const auto t_begin = chrono::steady_clock::now();

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	const auto t_cog = chrono::steady_clock::now();
	tms t1 = t_cog - t_begin;
	printf ("Creating ordered graph: %.6lf secs\n", t1.count());

	// Triangle counting for each edge
	edge* tc = (edge *) calloc (nEdge, sizeof(vertex));

	lol tricount = count_triangles (tc, nVtx, xadj, ordered_adj, ordered_xadj); // single loop with three locks

	const auto t_tc = chrono::steady_clock::now();
	tms t2 = t_tc - t_cog;
	printf ("Triangle counting time: %.6lf secs\n", t2.count());

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
	for (size_t i = 0; i < nEdge; i++)
		if (T[i] == -1)
			T[i] = 0;

	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nEdge, T, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

	free (T);
	cout << "Max truss number: " << tc_of_uv << endl;
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin;
	printf ("Total time: %.6lf secs\n", total.count());
	return;
}


// For partialnmLocal23 function
inline void AnotherupdateAndNotify (edge ind, vertex* P, int newP, vector<vertex>& neigs, int* changed) {
	P[ind] = newP;
	if (changed[ind] == 0)
		changed[ind] = 1;
	for (int k = 0; k < neigs.size(); k++) {
		if (P[neigs[k]] >= P[ind]) {
			if (changed[ind] == 0)
				changed[ind] = 1;
		}
	}
}

inline int AnotherefficientUpdateHI (edge ind, edge* xadj, vertex* adj, edge* xel, couple* el, vertex* T, int* changed) {

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

		AnotherupdateAndNotify (ind, T, j, neigs, changed);
		return 1;
	}
	return 0;
}


// Function in Sec. 5.3
void partialnmLocal23 (int* changed, couple* el, edge* xel, vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P) {

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < nEdge; ind++)
		if (changed[ind] == 1)
			mapInitialHI (ind, xadj, adj, xel, el, P);

	bool flag = true;
	while (flag) {
		const auto td1 = chrono::steady_clock::now();
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < nEdge; ind++) {
			if (changed[ind] != 1)
				continue;
			changed[ind] = 0;
			int a = AnotherefficientUpdateHI (ind, xadj, adj, xel, el, P, changed);
			if (a == 1)
				flag = true;
		}
	}
	return;
}


// For Fig. 14b. Applies partialAND on the ego network of a given sample edge. 4 * TEST_SIZE is the number of samples.
// Reads triangle-counts and truss-numbers from files
void converge23onEgo (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* K, string kfile) {

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;
	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);


	K = (vertex *) malloc (nEdge * sizeof(vertex));
	vertex* T = (vertex *) malloc (nEdge * sizeof(vertex));

	string fk = kfile + "_2300_TRUSS_H_0"; // Triangle counts are read
	FILE* fp = fopen (fk.c_str(), "r");
	vector<tuple<int, int>> trv;
	int i = 0;
	while (fscanf (fp, "%d", &(T[i])) != EOF) {
		trv.push_back (make_tuple(i, T[i]));
		i++;
	}
	sort (trv.begin(), trv.end(), kksort);
	fclose (fp);

	fk = kfile + "_2300_TRUSS_FINAL_K"; // Truss numbers are read
	fp = fopen (fk.c_str(), "r");

	vector<tuple<int, int>> cv;
	i = 0;
	while (fscanf (fp, "%d", &(K[i])) != EOF) {
		cv.push_back (make_tuple(i, K[i]));
		i++;
	}
	sort (cv.begin(), cv.end(), kksort);
	fclose (fp);

	vector<int> sample_edges;

	srand(time(NULL));

	// samples vertices from triangle-count-sorted array and truss-number-sorted-array
	for (int i = 1; i <= TEST_SIZE; i++) {
		int u = rand() % (i * nEdge / TEST_SIZE);
		sample_edges.push_back (get<0>(trv[u]));
		sample_edges.push_back (get<0>(cv[u]));
		sample_edges.push_back (get<0>(trv[i-1]));
		sample_edges.push_back (get<0>(cv[i-1]));
	}

	vertex* P = (vertex *) calloc (nEdge, sizeof(vertex));
	int changed[nEdge];
	for (int i = 0; i < nEdge; i++)
		changed[i] = -1;

	for (auto e : sample_edges) {
		memcpy (P, T, nEdge * sizeof(vertex));
		const auto t_begin = chrono::steady_clock::now();
		changed[e] = 1;

		vector<vertex> neigs;
		{
			int ind = e;
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
					changed[id1] = 1;
					changed[id2] = 1;
					neigs.push_back (id1);
					neigs.push_back (id2);
					i++;
					j++;
				}
			}
		}

		partialnmLocal23 (changed, el, xel, nVtx, nEdge, adj, xadj, P);
		const auto t_end = chrono::steady_clock::now();
		tms t1 = t_end - t_begin;
		printf ("Local convergence time: %.6lf secs\n", t1.count());

		int pe = P[e], ke = K[e], te = T[e];

		vector<vertex> gt, rs;
		for (int i = 0; i < neigs.size(); i+=2) {
			int id1 = neigs[i];
			int id2 = neigs[i+1];

			// easier
			int lk = max (K[id1], K[id2]);
			int lp = max (P[id1], P[id2]);

			if (lk >= ke) {
				gt.push_back (id1);
				gt.push_back (id2);
			}

			if (lp >= pe) {
				rs.push_back (id1);
				rs.push_back (id2);
			}
		}

		printf ("23_EGO -- P: %d  vs  K: %d\t", P[e], K[e]);
		sort (rs.begin(), rs.end());
		sort (gt.begin(), gt.end());

		int is = commons (rs, gt);
		double prec = (double) is / rs.size(), recall = (double) is / gt.size();
		double f1 = 2 * prec * recall / (prec + recall), jaccard = (double) is / (gt.size() + rs.size() - is);

		// Comparison of partialAND results (P) and the ground-truth truss numbers (T)
		printf ("base: %d\tjac: %lf\tprec: %lf\trec: %lf\tf1: %lf\t\t", gt.size(), jaccard, prec, recall, f1);

		{
			for (int i = 0; i < nEdge; i++)
				P[i] = T[i];

			int pe = P[e], ke = K[e];
			vector<vertex> gt, rs;
			for (int i = 0; i < neigs.size(); i+=2) {
				int id1 = neigs[i];
				int id2 = neigs[i+1];

				// easier
				int lk = max (K[id1], K[id2]);
				int lp = max (P[id1], P[id2]);

				if (lk >= ke) {
					gt.push_back (id1);
					gt.push_back (id2);
				}

				if (lp >= pe) {
					rs.push_back (id1);
					rs.push_back (id2);
				}
			}

			printf ("triCount_EGO -- P: %d  vs  K: %d\t", P[e], K[e]);
			sort (rs.begin(), rs.end());
			sort (gt.begin(), gt.end());

			int is = commons (rs, gt);
			double prec = (double) is / rs.size(), recall = (double) is / gt.size();
			double f1 = 2 * prec * recall / (prec + recall), jaccard = (double) is / (gt.size() + rs.size() - is);

			// Comparison of triangle-counts (T) and the ground-truth truss numbers (T)
			printf ("base: %d\tjac: %lf\tprec: %lf\trec: %lf\tf1: %lf\n", gt.size(), jaccard, prec, recall, f1);

		}
	}
}
