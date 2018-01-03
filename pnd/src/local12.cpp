#include "main.h"

inline int sortInitialHI (vertex ind, vertex* adj, edge* xadj, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	vector<vertex> counts;
	for (edge j = xadj[ind]; j < xadj[ind+1]; j++) {
		vertex w = adj[j];
		counts.push_back (P[w]);
	}
	sort (counts.begin(), counts.end(), greater<vertex>());
	vertex j;
	for (j = 0; j < counts.size(); j++)
		if (counts[j] < j+1)
			break;

#ifdef SYNC
	if (Q[ind] != j) {
		Q[ind] = j;
		return 1;
	}
#else
	if (P[ind] != j) {
		P[ind] = j;
		return 1;
	}
#endif
	return 0;
}

inline int regularUpdateHI (vertex ind, vertex* adj, edge* xadj, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {
	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;
	vertex ne = xadj[ind+1] - xadj[ind];
	for (edge j = xadj[ind]; j < xadj[ind+1]; j++) {
		vertex w = adj[j];
		vertex pw = P[w];
		if (pw >= previous_P)
			greaterequals++;
		else
			smallers[pw]++;
		if (greaterequals == previous_P) {
			yep_set = true;
			break;
		}
	}

	if (!yep_set && ne > 0) { // watch for degree zeros
		vertex j = 0;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j)
				break;
		}
#ifdef SYNC
		if (Q[ind] != j) {
			Q[ind] = j;
			return 1;
		}
#else
		if (P[ind] != j) {
			P[ind] = j;
			return 1;
		}
#endif
	}
	return 0;
}

inline void updateAndNotify (vertex ind, vertex* P, vertex newP, vector<vertex>& neigs, bool* changed) {
	P[ind] = newP;
	changed[ind] = true; // *THIS*
	for (vertex k = 0; k < neigs.size(); k++)
		if (P[neigs[k]] >= P[ind])
			changed[neigs[k]] = true;
}

// this is faster than sort-based computation
inline int mapInitialHI (vertex ind, vertex* adj, edge* xadj, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {

	HashMap<vertex> gmap (-1);
	vertex greaters = 0;
	vertex equals = 0;
	vertex H = 0;

	for (vertex j = xadj[ind]; j < xadj[ind+1]; j++) {
		vertex w = adj[j];
		vertex sm = P[w];
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

	if (H < oP)
		return 1;
	else
		return 0;
}

inline int efficientUpdateHI (vertex ind, vertex* adj, edge* xadj, vertex* P, bool* changed) {
	vector<vertex> neigs;
	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;
	vertex ne = xadj[ind+1] - xadj[ind];
	for (edge j = xadj[ind]; j < xadj[ind+1]; j++) {
		vertex w = adj[j];
		vertex pw = P[w];
		if (pw <= previous_P)
			neigs.push_back (w);
		if (pw >= previous_P)
			greaterequals++;
		else
			smallers[pw]++;
		if (greaterequals == previous_P) {
			yep_set = true;
			break;
		}
	}

	if (!yep_set && ne > 0) { // watch for degree zeros
		vertex j;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j)
				break;
		}
		updateAndNotify (ind, P, j, neigs, changed);
		return 1;
	}
	return 0;
}


// base AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void baseLocal12 (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile) {

	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

#ifdef SYNC
	printf ("It is SYNC\n");
	vertex* Q = (vertex *) malloc (sizeof(vertex) * nVtx);
#else
	printf ("It is ASYNC\n");
#endif

	timestamp t_deg;
	cout << "Degree finding time: " << t_deg - t_begin << endl;

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
	for (vertex ind = 0; ind < nVtx; ind++) {
#ifdef DEBUG_000
		ccs[tn] += mapInitialHI (ind, adj, xadj, P
#ifdef SYNC
				, Q
#endif
		);
#else
		mapInitialHI (ind, adj, xadj, P
#ifdef SYNC
				, Q
#endif
		);
#endif
	}

#ifdef SYNC
	memcpy (P, Q, sizeof(vertex) * nVtx);
#endif

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_deg - td << endl;
#ifdef DUMP_Hs
	print_Ks (nVtx, P, vfile, oc);
#endif
	oc++;


#ifdef DEBUG_000
	vertex degisenler = 0;
	for(int i = 0; i < nt; i++)
		degisenler += ccs[i];
	memset (ccs, 0, nt * sizeof(int));
#endif

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex ind = 0; ind < nVtx; ind++) {
			int fl = regularUpdateHI (ind, adj, xadj, P
#ifdef SYNC
					, Q
#endif
			);

			if (fl == 1)
				flag = true;
		}

		timestamp td2;
#ifdef SYNC
		memcpy (P, Q, sizeof(vertex) * nVtx);
#endif

#ifdef DUMP_Hs
		print_Ks (nVtx, P, vfile, oc);
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
	print_Ks (nVtx, P, vfile);
#endif

	free (P);
#ifdef SYNC
	free (Q);
#endif

	return;
}


// AND algorithm with the notification mechanism
void nmLocal12 (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	timestamp t_deg;
	cout << "degreeFinding time: " << t_deg - t_begin << endl;

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


	bool changed[nVtx];
	memset (changed, 255, sizeof(bool) * nVtx); // set all true

#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++) {
#ifdef DEBUG_000
		ccs[tn] += mapInitialHI (ind, adj, xadj, P);
#else
		mapInitialHI (ind, adj, xadj, P);
#endif
	}

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_deg - td << endl;
#ifdef DUMP_Hs
	print_Ks (nVtx, P, vfile, oc);
#endif
	oc++;


#ifdef DEBUG_000
	vertex degisenler = 0;
	for(int i = 0; i < nt; i++)
		degisenler += ccs[i];
	memset (ccs, 0, nt * sizeof(int));
#endif

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, adj, xadj, P, changed);
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
		print_Ks (nVtx, P, vfile, oc);
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
	print_Ks (nVtx, P, vfile);
#endif

	free (P);
	return;
#endif
}

// Sequential k-core computation
void kcore (vertex nVtx, vertex* adj, edge* xadj, vertex* K, const char* vfile) {

	timestamp t_begin;
	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		auto deg = xadj[i+1] - xadj[i];
		if (deg > max_degree)
			max_degree = deg;
	}
	K = (vertex *) malloc (nVtx * sizeof(vertex));


	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize(max_degree, nVtx);
	for (size_t i = 0; i < nVtx; i++) {
		auto deg = xadj[i+1] - xadj[i];
		if (deg > 0)
			na_bs.Insert (i, deg);
	}

	vertex degree_of_u = 0;
	while (1) {
		vertex u, val;
		if (na_bs.PopMin(&u, &val) == -1) // if the bucket is empty
			break;

		if (val == 0)
			continue;

		degree_of_u = K[u] = val;

		for (int i = xadj[u]; i < xadj[u+1]; i++) { /* decrease the degree of the neighbors with greater degree */
			vertex v = adj[i];
			if (na_bs.CurrentValue(v) > val)
				na_bs.DecVal(v);
		}
	}

	na_bs.Free();
	cout << "Max core number: " << degree_of_u << endl;
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (nVtx, K, vfile);
#endif
	return;
}







// ONGOING WORK

// Finds the max K value by iterating on top-number high degree vertices
void fast12DegeneracyNumber (vertex nVtx, vertex* adj, edge* xadj, vertex* P, vertex topK) {

	int number = topK; // only top-number vertices in degree are executed
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	timestamp t_tc;
	cout << "Degree finding time: " << t_tc - t_begin << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

	bool changed[nVtx];
	memset (changed, 255, sizeof(bool) * nVtx); // set all true

	vector<eda> KK;
	for (vertex i = 0; i < nVtx; i++)
		KK.push_back (make_tuple(i, P[i]));

	sort (KK.begin(), KK.end(), kksort);

	int start = 0;
	for (int i = start; i < number; i++)
		changed[get<0>(KK[i])] = true;



	while (flag) {

		memset (changed, 0, sizeof(bool) * nVtx); // set all false
		for (int i = start; i < number; i++)
			changed[get<0>(KK[i])] = true; // only top-k true

		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, adj, xadj, P, changed);
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


	int maxK = 0;
	for (int i = start; i < number; i++) {
		int curP = P[get<0>(KK[i])];
		if (curP > maxK)
			maxK = curP;
//		printf ("%d goes from %d to %d\n", get<0>(KK[i]), get<1>(KK[i]), curP);
	}

	printf ("max K: %d\n", maxK);
#endif
}


/*
- L_1: All those vertices with the minimum degree.
- To find L_i: delete all vertices in L_j (j < i). Then find all vertices with the minimum degree in the remaining graph.
*/
void kcore_levels (vertex nVtx, vertex* adj, edge* xadj, vertex* L, const char* vfile) {

	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		auto deg = xadj[i+1] - xadj[i];
		if (deg > max_degree)
			max_degree = deg;
	}
	L = (vertex *) malloc (nVtx * sizeof(vertex));

	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize(max_degree, nVtx);
	for (size_t i = 0; i < nVtx; i++) {
		auto deg = xadj[i+1] - xadj[i];
		if (deg > 0)
			na_bs.Insert (i, deg);
	}
	vertex degree_of_u = 0;

	int level = 1;
	while (1) {
		vertex u, min_val, val;
		vector<vertex> mins;

		if (na_bs.PopMin(&u, &min_val) == -1) // if the bucket is empty
			break;

		if (min_val == 0)
			continue;

		mins.push_back (u);

		while (1) {
			if (na_bs.PopMin(&u, &val) == -1) // if the bucket is empty
				break;
			if (val > min_val) {
				na_bs.Insert (u, val); // put it back
				break;
			}
			mins.push_back (u);
		}

		for (auto v : mins) {
			L[v] = min_val;
			for (int i = xadj[v]; i < xadj[v+1]; i++) { /* decrease the degree of the neighbors with greater degree */
				vertex w = adj[i];
				if (na_bs.CurrentValue(w) > min_val)
					na_bs.DecVal(w);
			}
		}
		printf ("Level %d K: %d size: %d\n", level, min_val, mins.size());
		level++;
	}
	na_bs.Free();
#ifdef DUMP_K
	print_Ks (nVtx, L, vfile);
#endif
	return;
}


/*
- L_1: All those vertices whose degree is at most K-value.
- To find L_i: delete all vertices in L_j (j < i). Then find all vertices whose degree in the remaining graph is at most (original) K-value.
*/
void kcore_Sesh_levels (vertex nVtx, vertex* adj, edge* xadj, vertex* K, vertex* L, const char* vfile) {

	L = (vertex *) malloc (nVtx * sizeof(vertex));
	vertex* degs = (vertex *) malloc (nVtx * sizeof(vertex));
	for (size_t i = 0; i < nVtx; i++)
		degs[i] = xadj[i+1] - xadj[i];

	int level = 1;
	while (1) {
		vector<vertex> atMostK;
		for (vertex i = 0; i < nVtx; i++) {
			if (degs[i] != -1 && degs[i] <= K[i])
				atMostK.push_back(i);
		}

		if (atMostK.empty())
			break;

		for (vertex v : atMostK) {
			for (vertex j = xadj[v]; j < xadj[v+1]; j++) {
				vertex w = adj[j];
				if (degs[w] != -1)
					degs[w]--;
			}
			degs[v] = -1;
			L[v] = level;
		}

		printf ("Level %d  size: %d\n", level, atMostK.size());
		level++;
	}
#ifdef DUMP_K
	print_Ks (nVtx, L, vfile);
#endif
	return;
}


