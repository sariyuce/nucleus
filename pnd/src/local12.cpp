#include "main.h"

// UTILITY FUNCTIONS, MIGHT BE NEEDED

inline void print_Ks (int nVtx, volatile vertex* T, const char* vfile, int H = -1) {
	string st (vfile);
	if (H == -1)
		st += "_FINAL_K";
	else
		st += "_H_" + to_string(H);
	FILE* pp = fopen (st.c_str(), "w");
	for (int i = 0; i < nVtx; i++)
		fprintf (pp, "%d\n", T[i]);
	fclose (pp);
}

inline void compute_KT (vertex* Reals, Graph& graph, vertex* T, int oc) {
	int count = 0, score = 0;
	for (vertex i = 0; i < graph.size(); i++) {
		int u = i;
		for (vertex j = 0; j < graph[i].size(); j++) {
			int v = graph[i][j];
			if (u == smaller (graph, u, v))
				EdgeKendallTau (Reals, T, u, v, &count, &score);
		}
	}
	printf ("H %d , KendallTau: %lf\n", oc, kt (count, score));
}

inline void simple_distance (vertex* Reals, Graph& graph, vertex* T, int oc) {
	double score = 0;
	int count = 0;
	for (vertex i = 0; i < graph.size(); i++)
		if (T[i] > 0 && Reals[i] > 0) {
			//			printf ("sc: %lf\n", (double) Reals[i] / T[i]);
			//			score += (double) Reals[i] / T[i];
			if (Reals[i] == T[i])
				score++;
			//			else
			//				printf ("%d: Reals: %d   K: %d\n", i, Reals[i], T[i]);

			count++;
		}

	//	score /= count;
	printf ("H %d , similarity: %lf = %d / %d\n", oc, score/count, (int) score, count);
}








// ESSENTIALS

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

	HashMap<vertex> gmap (0);
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
void baseLocal12 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {

	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

#ifdef SYNC
	printf ("it is SYNC\n");
	vertex* Q = (vertex *) malloc (sizeof(vertex) * nVtx);
#else
	printf ("it is ASYNC\n");
#endif

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


#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++) {
#ifdef DEBUG_000
		ccs[tn] += mapInitialHI (ind, adj, xadj, P);
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

#pragma omp parallel for default (shared)
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
void nmLocal12 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	printf ("it is ASYNC\n");

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

#pragma omp parallel for schedule (dynamic, 1)
	for (vertex ind = 0; ind < nVtx; ind++) {
#ifdef DEBUG_000
		ccs[tn] += mapInitialHI (ind, adj, xadj, P);
#else
		mapInitialHI (ind, adj, xadj, P);
//		TRY_mapInitialHI (ind, adj, xadj, P, changed);
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

#pragma omp parallel for schedule (dynamic, 1)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
#ifdef DEBUG_000
			ccs[tn] += efficientUpdateHI (ind, adj, xadj, P, changed);
#else
			int a = efficientUpdateHI (ind, adj, xadj, P, changed);
			if (a == 1)
				flag = true;
#endif
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


// AND algorithm with the notification mechanism and no wait in the outer while loop
void NoWaitnmLocal12 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for No Wait notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));
	int nThreads = 1;
	volatile int* counter;
#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
		counter = (volatile int *) calloc (nThreads, sizeof(int));
	}

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];
	printf ("it is ASYNC\n");
	timestamp t_deg;
	cout << "degreeFinding time: " << t_deg - t_begin << endl;

	bool changed[nVtx];
	memset (changed, 255, sizeof(bool) * nVtx); // set all true
#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++)
		mapInitialHI (ind, adj, xadj, P);
	timestamp t_init;
	cout << "H 0 time: " << t_init - t_deg << endl;
#ifdef DUMP_Hs
	print_Ks (nVtx, P, vfile, oc);
#endif

	bool flags[nThreads];
#pragma omp parallel
	{
		int nt = omp_get_thread_num();
		flags[nt] = true;
		while (flags[nt]) {
			flags[nt] = false;
			counter[omp_get_thread_num()]++;

#pragma omp for nowait schedule (dynamic, 1000) // (static, nVtx/nThreads)
			for (vertex ind = 0; ind < nVtx; ind++) {
				if (!changed[ind])
					continue;
				changed[ind] = false;
				int a = efficientUpdateHI (ind, adj, xadj, P, changed);
				if (a == 1)
					flags[nt] = true;
			}
		}
	}

	for (int i = 0; i < nThreads; i++)
		printf ("counter[%d]: %d\n", i, counter[i]);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (nVtx, P, vfile);
#endif
	free (P);
	return;
#endif
}




typedef tuple<int, int> eda;
bool kksort (eda i, eda j) { return (get<1>(i) > get<1>(j)); }

// Finds the max K value by iterating on top-number high degree vertices
void fast12DegeneracyNumber (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {

	int number = 200; // only top-number vertices in degree are executed
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	printf ("it is ASYNC\n");

	timestamp t_tc;
	cout << "degreeFinding time: " << t_tc - t_begin << endl;

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

		memset (changed, 255, sizeof(bool) * nVtx); // set all true
		for (int i = start; i < number; i++)
			changed[get<0>(KK[i])] = true;

		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (!changed[ind])
				continue;
			changed[ind] = false;
#ifdef DEBUG_000
			ccs[tn] += efficientUpdateHI (ind, adj, xadj, P, changed);
#else
			int a = efficientUpdateHI (ind, adj, xadj, P, changed);
			if (a == 1)
				flag = true;
#endif
		}

#ifdef DEBUG_000
		int degisenler = 0;
		for(int i = 0; i < nt; i++)
			degisenler += ccs[i];
		memset (ccs, 0, nt * sizeof(int));
		cout << "CHANGEDS: " << degisenler << endl;
#endif

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

