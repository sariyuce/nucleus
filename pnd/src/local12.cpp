#include "main.h"

#define TRIAL

// this is faster than sort-based computation
inline int mapInitialHI (vertex ind, vertex* adj, edge* xadj, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {

	unordered_map<vertex, vertex> gmap;
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

	const auto t_begin = chrono::steady_clock::now();
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

	const auto t_deg = chrono::steady_clock::now();
	tms t1 = t_deg - t_begin;
	printf ("Degree finding time: %.6lf secs\n", t1.count());

	int oc = 0;
	bool flag = true;

#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++) {
		mapInitialHI (ind, adj, xadj, P
#ifdef SYNC
				, Q
#endif
		);
	}

#ifdef SYNC
	memcpy (P, Q, sizeof(vertex) * nVtx);
#endif

	const auto t_init = chrono::steady_clock::now();;
	tms t2 = t_init - t_deg	;
	printf ("H %d time: %.6lf secs\n", oc, t2.count());
	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
	const auto ts1 = chrono::steady_clock::now();
	print_Ks (nVtx, P, vfile, oc);
	const auto ts2 = chrono::steady_clock::now();
	td += ts2 - ts1;
#endif
	oc++;

	while (flag) {
		const auto td1 = chrono::steady_clock::now();
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


#ifdef SYNC
		memcpy (P, Q, sizeof(vertex) * nVtx);
#endif

		const auto td2 = chrono::steady_clock::now();
#ifdef DUMP_Hs
		const auto ts1 = chrono::steady_clock::now();
		print_Ks (nVtx, P, vfile, oc);
		const auto ts2 = chrono::steady_clock::now();
		td += ts2 - ts1;
#endif

		tms step = td2 - td1;
		printf ("H %d time: %.6lf secs\n", oc, step.count());
		oc++;
	}

#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nVtx, P, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

	free (P);
#ifdef SYNC
	free (Q);
#endif
	printf ("Converges at %d\n", oc);
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin - td;
	printf ("Total time: %.6lf secs\n", total.count());
	return;
}

// AND algorithm with the notification mechanism
void nmLocal12 (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	const auto t_begin = chrono::steady_clock::now();
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	const auto t_deg = chrono::steady_clock::now();
	tms t1 = t_deg - t_begin;
	printf ("Degree finding time: %.6lf secs\n", t1.count());

	int oc = 0;
	bool flag = true;

	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
	const auto ts1 = chrono::steady_clock::now();
	print_Ks (nVtx, P, vfile, oc);
	oc++;
	const auto ts2 = chrono::steady_clock::now();
	td += ts2 - ts1;
#endif

	int nt;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
	}

	int counters[nt];

	for (int i = 0; i < nt; i++)
		counters[i] = 0;

	bool changed[nVtx];
	memset (changed, 255, sizeof(bool) * nVtx); // set all true

#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++) {
		counters[omp_get_thread_num()]++;
		mapInitialHI (ind, adj, xadj, P);
	}

	const auto t_init = chrono::steady_clock::now();;
	tms t2 = t_init - t_deg	;
	printf ("H %d time: %.6lf secs\n", oc, t2.count());
//	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
//	const auto ts1 = chrono::steady_clock::now();
	print_Ks (nVtx, P, vfile, oc);
//	const auto ts2 = chrono::steady_clock::now();
	td += ts2 - ts1;
#endif
	oc++;

	while (flag) {
		const auto td1 = chrono::steady_clock::now();
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (!changed[ind])
				continue;
			counters[omp_get_thread_num()]++;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, adj, xadj, P, changed);
			if (a == 1)
				flag = true;
		}
		const auto td2 = chrono::steady_clock::now();

#ifdef DUMP_Hs
		const auto ts1 = chrono::steady_clock::now();
		print_Ks (nVtx, P, vfile, oc);
		const auto ts2 = chrono::steady_clock::now();
		td += ts2 - ts1;
#endif

		tms step = td2 - td1;
		printf ("H %d time: %.6lf secs\n", oc, step.count());
		oc++;
	}

#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nVtx, P, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

#ifdef TRIAL
	int maxK = 0;
	for (int i = 0; i < nVtx; i++)
		if (P[i] > maxK)
			maxK = P[i];
	cout << "maxK: " << maxK << endl;

	int total_count = 0;
	for (int i = 0; i < nt; i++)
		total_count += counters[i];
	cout << " # operations: " << total_count << " = " << (double) total_count/nVtx << " * |V|" << endl;

#endif

	free (P);
	printf ("Converges at %d\n", oc);
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin - td;
	printf ("Total time: %.6lf secs\n", total.count());





	return;
#endif
}

// Sequential k-core computation
void kcore (vertex nVtx, vertex* adj, edge* xadj, vertex* K, const char* vfile) {

	const auto t_begin = chrono::steady_clock::now();
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
		else
			K[i] = 0;
	}

	vertex degree_of_u = 0;
	while (1) {
		vertex u, val;
		if (na_bs.PopMin(&u, &val) == -1) // if the bucket is empty
			break;

		if (val == 0)
			continue;

		degree_of_u = K[u] = val;

		for (auto i = xadj[u]; i < xadj[u+1]; i++) { /* decrease the degree of the neighbors with greater degree */
			vertex v = adj[i];
			if (na_bs.CurrentValue(v) > val)
				na_bs.DecVal(v);
		}
	}

	na_bs.Free();
	tms td  = chrono::duration<double>::zero();
#ifdef DUMP_K
	const auto ts3 = chrono::steady_clock::now();
	print_Ks (nVtx, K, vfile);
	const auto ts4 = chrono::steady_clock::now();
	td += ts4 - ts3;
#endif

	free(K);
	cout << "Max core number: " << degree_of_u << endl;
	const auto t_end = chrono::steady_clock::now();
	tms total = t_end - t_begin - td;
	printf ("Total time: %.6lf secs\n", total.count());
	return;
}















// tries to find degeneracy number
void topKs (vertex nVtx, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	const auto t_begin = chrono::steady_clock::now();
	P = (vertex *) calloc (nVtx, sizeof(vertex));

#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	const auto t_deg = chrono::steady_clock::now();
	tms t1 = t_deg - t_begin;
	printf ("Degree finding time: %.6lf secs\n", t1.count());

	vector<tuple<int, int>> b_degs;
	for (vertex i = 0; i < nVtx; i++)
		b_degs.push_back (make_tuple(i, P[i]));

	sort (b_degs.begin(), b_degs.end(), kksort);



	int nt;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
	}

	int counters[nt];



	int topKs[7];// = {round (nVtx*0.0001), round (nVtx*0.0005), round (nVtx*0.001), round (nVtx*0.005), round (nVtx*0.01), round (nVtx*0.05), round (nVtx*0.1)};
	for (int tk : topKs) {

#pragma omp parallel for default (shared)
		for (vertex i = 0; i < nVtx; i++)
			P[i] = xadj[i+1] - xadj[i];

		for (int i = 0; i < nt; i++)
			counters[i] = 0;

		vector<tuple<int, int>> degs (b_degs);

		int oc = 0;
		bool flag = true;

		bool visited[nVtx];
		bool changed[nVtx];
		for (int i = 0; i < nVtx; i++) {
			changed[i] = false;
			visited[i] = false;
		}

		for (int i = 0; i < tk; i++) {
			changed[get<0>(degs[i])] = true;
			visited[get<0>(degs[i])] = true;
		}


#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex ind = 0; ind < nVtx; ind++) {
			if (changed[ind]) {
				counters[omp_get_thread_num()]++;
//				visited[ind] = true;
				mapInitialHI (ind, adj, xadj, P);
			}
		}



		const auto t_init = chrono::steady_clock::now();;
		tms t2 = t_init - t_deg	;
		printf ("H %d time: %.6lf secs\n", oc, t2.count());
		tms td  = chrono::duration<double>::zero();
#ifdef DUMP_Hs
		const auto ts1 = chrono::steady_clock::now();
		print_Ks (nVtx, P, vfile, oc);
		const auto ts2 = chrono::steady_clock::now();
		td += ts2 - ts1;
#endif
		oc++;




		{
			int total_count = 0;
			for (int i = 0; i < nt; i++)
				total_count += counters[i];


			int maxK = 0;
			for (int i = 0; i < nVtx; i++) {
				if (visited[i]) {
					if (P[i] > maxK)
						maxK = P[i];
				}
			}
			cout << "top " << tk << " vertices; maxK: " << maxK << " @ step-" << oc;
			cout << " , # operations: " << total_count << " = " << (double) total_count/nVtx << " * |V|" << endl;


		}



		while (flag) {
			const auto td1 = chrono::steady_clock::now();
			flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
			for (vertex i = 0; i < nVtx; i++) {
				vertex ind = i;
				if (!visited[ind] || !changed[ind])
					continue;
				counters[omp_get_thread_num()]++;
//				visited[ind] = true;
				changed[ind] = false;
				int a = efficientUpdateHI (ind, adj, xadj, P, changed);
				if (a == 1)
					flag = true;
			}
			const auto td2 = chrono::steady_clock::now();

	#ifdef DUMP_Hs
			const auto ts1 = chrono::steady_clock::now();
			print_Ks (nVtx, P, vfile, oc);
			const auto ts2 = chrono::steady_clock::now();
			td += ts2 - ts1;
	#endif

			tms step = td2 - td1;
			printf ("H %d time: %.6lf secs\n", oc, step.count());
			oc++;


			int maxK = 0;
			for (int i = 0; i < nVtx; i++) {
				if (visited[i]) {
					if (P[i] > maxK)
						maxK = P[i];
				}
			}

			int total_count = 0;
			for (int i = 0; i < nt; i++)
				total_count += counters[i];

			cout << "top " << tk << " vertices; maxK: " << maxK << " @ step-" << oc;
			cout << " , # operations: " << total_count << " = " << (double) total_count/nVtx << " * |V|" << endl;
		}

	#ifdef DUMP_K
		const auto ts3 = chrono::steady_clock::now();
		print_Ks (nVtx, P, vfile);
		const auto ts4 = chrono::steady_clock::now();
		td += ts4 - ts3;
	#endif

		int processed = 0;
		int maxK = 0;
		for (int i = 0; i < nVtx; i++) {
			if (visited[i]) {
				processed++;
				if (P[i] > maxK)
					maxK = P[i];
			}
		}





		vector<tuple<int, int>> topK_nodes;

		for (int i = 0; i < nVtx; i++) {
			if (visited[i]) {
				topK_nodes.push_back (make_tuple(i, P[i]));
			}
		}
		sort (topK_nodes.begin(), topK_nodes.end(), kksort);

		for (int i = 0; i < topK_nodes.size(); i++)
			printf ("Top nodes: P[ %d ]: %d\n", get<0>(topK_nodes[i]), get<1>(topK_nodes[i]));





		int total_count = 0;
		for (int i = 0; i < nt; i++)
			total_count += counters[i];


		cout << "top " << tk << " vertices, maxK: " << maxK;
		cout << " # visiteds: " << processed << " = " << (double) processed / nVtx << " * |V|";
		cout << " # operations: " << total_count << " = " << (double) total_count/nVtx << " * |V|" << endl;


		printf ("Converges at %d\n", oc);
		const auto t_end = chrono::steady_clock::now();
		tms total = t_end - t_begin - td;
		printf ("Total time: %.6lf secs\n", total.count());
	}
	free (P);
	return;
#endif
}












void maxcore (vertex v, vertex nVtx, vertex* adj, edge* xadj, vertex* K, vector<vertex>& result) {
	queue<vertex> bfs;
	bfs.push (v);
	unordered_map<vertex, bool> visited;
	visited[v] = true;
	vertex core_number = K[v];

	while (!bfs.empty()) {
		vertex u = bfs.front();
		bfs.pop();
		result.push_back (u);
		for (int i = xadj[u]; i < xadj[u+1]; i++) {
			int w = adj[i];
			if (visited.find(w) == visited.end() && K[w] >= core_number) {
				visited[w] = true;
				bfs.push (w);
			}
		}
	}
}


//#define TEST_SIZE 5
#define TEST_SIZE 100


void find_mcore (vertex nVtx, vertex* adj, edge* xadj, vertex* K, string kfile) {

//	vector<int> sample_vertices; // select random vertices and put here

	string fk = kfile + "_12_NUCLEI";
	ifstream f (fk);
//	FILE* fp = fopen (fk.c_str(), "r");
	int i = 0, a;
	vector<vector<vertex>> maxes;
	vector<vertex> mc;

	cout << fk << endl;
	string line;
	while (getline(f, line)) {
		istringstream iss(line);
		int a;
		for (int i = 0; i < 4; i++) {
			iss >> a;
//			printf ("a: %d\n", a);
		}
		double dd;
		iss >> dd;
		iss >> a;
		iss >> a;

		iss >> a;
		while (a != -1) {
			mc.push_back (a);
			iss >> a;
		}
		maxes.push_back (mc);
		mc.clear();
//		for (auto a : mc)
//			printf ("%d ", a);
//		printf ("\n");
	}

	printf ("done reading nuclei\n");


//	string fk = kfile + "_1200_CORE_FINAL_K";
//	FILE* fp = fopen (fk.c_str(), "r");
	K = (vertex *) malloc (nVtx * sizeof(vertex));

	/*
	vector<tuple<int, int>> cv;
	int i = 0;
	while (fscanf (fp, "%d", &(K[i])) != EOF) {
		cv.push_back (make_tuple(i, K[i]));
		i++;
	}
	sort (cv.begin(), cv.end(), kksort);

*/

//	for (int i = 0; i < cv.size(); i++) {
//		printf ("In cv, core of %d is %d\n", get<0>(cv[i]), get<1>(cv[i]));
//	}

/*
	fk = kfile + "_1200_CORE_H_0";
	fp = fopen (fk.c_str(), "r");
	vector<tuple<int, int>> dv;
	i = 0;
	while (fscanf (fp, "%d", &(K[i])) != EOF) {
		dv.push_back (make_tuple(i, K[i]));
		i++;
	}
	sort (dv.begin(), dv.end(), kksort);

*/


//	for (int i = 1; i <= TEST_SIZE; i++) {
//		sample_vertices.push_back (get<0>(cv[i-1]));
//	}
//
//	for (int i = 1; i <= TEST_SIZE; i++) {
//		sample_vertices.push_back (get<0>(dv[i-1]));
//	}

//	srand(time(NULL));
//	for (int i = 1; i <= TEST_SIZE; i++) {
//		int u = rand() % (i * nVtx / TEST_SIZE);
//		sample_vertices.push_back (get<0>(cv[u]));
//	}
//
//	for (int i = 1; i <= TEST_SIZE; i++) {
//		int u = rand() % (i * nVtx / TEST_SIZE);
//		sample_vertices.push_back (get<0>(dv[u]));
//	}
//
//	printf ("random vertices are selected\n");


//	sample_vertices.push_back (10);
//	sample_vertices.push_back (100);
//	sample_vertices.push_back (1000);
//	sample_vertices.push_back (4000);
//
//	sample_vertices.push_back (1912);
//	sample_vertices.push_back (1917);
//	sample_vertices.push_back (1918);
//	sample_vertices.push_back (1929);


	vector<vector<vector<vector<vertex>>>> all_results (maxes.size()); // each item is a list of subgraphs in each iteration
	// all_results :  list of all leaves
	// all_results [i] : for a leaf, list of maxcores that its each vertex participates
	// all_results [i][j] : a vector, maxcore of a vertex
	string fl = kfile + "_1200_CORE_H_";
	int iteration = 0;
	vector<vector<vertex>> degs (maxes.size()), core_numbers (maxes.size());
	for (int i = 0; i < maxes.size(); i++) {
		all_results[i].resize (maxes[i].size());
//			degs[i].resize (maxes[i].size());
//			core_numbers[i].resize (maxes[i].size());
	}

	while (true) {
		// each iteration is for an H file
		string kf = fl + to_string (iteration);
		FILE* fp = fopen (kf.c_str(), "r");
		if (fp == NULL) {
			for (int i = 0; i < maxes.size(); i++) {
				for (int j = 0; j < maxes[i].size(); j++) {
					vertex v = maxes[i][j];
					core_numbers[i].push_back (K[v]);
				}
			}
			break;
		}
		int i = 0;
		while (fscanf (fp, "%d", &(K[i])) != EOF)
			i++;
		for (int i = 0; i < maxes.size(); i++) {
			for (int j = 0; j < maxes[i].size(); j++) {
				int v = maxes[i][j];
				vector<vertex> result;
				maxcore (v, nVtx, adj, xadj, K, result);
				sort (result.begin(), result.end());
				all_results[i][j].push_back (result);
				if (iteration == 0)
					degs[i].push_back(K[v]);
			}
		}
		iteration++;
		printf ("iter %d done\n", iteration);
	}

	string of = kfile + "_12_LEAVES_EVOL";
	FILE* op = fopen (of.c_str(), "w");
	for (int i = 0; i < all_results.size(); i++) {
		printf ("subgraph %d / %d\n", i+1, all_results.size());
		for (int j = 0; j < all_results[i].size(); j++) {
			fprintf (op, "LEAF %d\tVERTEX %d: Degree: %d  Core_Number: %d\n", i, maxes[i][j], degs[i][j], core_numbers[i][j]);
			auto vv = all_results[i][j];
			int iter = 0;
			for (auto v : vv) {
				int is = commons (v, vv[vv.size()-1]);
				fprintf (op, "Leaf %d Iteration %d, size: %d ( wrt %d ) jaccard: %lf\t\t", i, iter++, v.size(), vv[vv.size()-1].size(), (double) (is) / (v.size() + vv[vv.size()-1].size() - is));
				double prec = (double) is / v.size(), recall = (double) is / vv[vv.size()-1].size();
				fprintf (op, "Precision: %lf Recall: %lf F1: %lf\n", prec, recall, 2 * prec * recall / (prec + recall));
			}
			fprintf (op, "\n");
		}
		fprintf (op, "\n\n");
	}
	fclose (op);

//
//	for (auto vv : all_results) {
//		int iter = 0;
//		for (auto v : vv) {
//			printf ("Iteration %d:  ", iter++);
//			for (auto i : v) {
//				printf ("%d ", i);
//			}
//			printf ("\n");
//		}
//		printf ("\n\n");
//	}


}

inline void AnotherupdateAndNotify (vertex ind, vertex* P, vertex newP, vector<vertex>& neigs, int* changed) {
	P[ind] = newP;
	if (changed[ind] == 0)
		changed[ind] = 1;
	for (vertex k = 0; k < neigs.size(); k++)
		if (P[neigs[k]] >= P[ind]) {
			if (changed[neigs[k]] == 0) // THIS is important!!!
				changed[neigs[k]] = 1;
		}
}

inline int AnotherefficientUpdateHI (vertex ind, vertex* adj, edge* xadj, vertex* P, int* changed) {
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
		AnotherupdateAndNotify (ind, P, j, neigs, changed);
		return 1;
	}
	return 0;
}


void partialnmLocal12 (int* changed, vertex nVtx, vertex* adj, edge* xadj, vertex* P) {



#pragma omp parallel for default (shared)
	for (vertex i = 0; i < nVtx; i++)
		P[i] = xadj[i+1] - xadj[i];

	int oc = 0;
	bool flag = true;

	// NO changed, it's coming from argument

#pragma omp parallel for schedule (dynamic, 1000)
	for (vertex ind = 0; ind < nVtx; ind++) {
		if (changed[ind] == 1)
			mapInitialHI (ind, adj, xadj, P);
	}

	oc++;

	while (flag) {
//		const auto td1 = chrono::steady_clock::now();
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (vertex i = 0; i < nVtx; i++) {
			vertex ind = i;
			if (changed[ind] != 1)
				continue;
			changed[ind] = 0;
			int a = AnotherefficientUpdateHI (ind, adj, xadj, P, changed);
			if (a == 1)
				flag = true;
		}
//		const auto td2 = chrono::steady_clock::now();


//		tms step = td2 - td1;
//		printf ("H %d time: %.6lf secs\n", oc, step.count());
		oc++;
	}

//	printf ("Converges at %d\n", oc);
	return;
}




void converge12onEgo (vertex nVtx, vertex* adj, edge* xadj, vertex* K, string kfile) {

	string 	fk = kfile + "_1200_CORE_H_0";
	FILE* fp = fopen (fk.c_str(), "r");
	K = (vertex *) malloc (nVtx * sizeof(vertex));
	vertex* D = (vertex *) malloc (nVtx * sizeof(vertex));

	vector<tuple<int, int>> dv;
	int i = 0;
	while (fscanf (fp, "%d", &(D[i])) != EOF) {
		dv.push_back (make_tuple(i, D[i]));
		i++;
	}
	sort (dv.begin(), dv.end(), kksort);

//	for (int i = 0; i < cv.size(); i++) {
//		printf ("In cv, core of %d is %d\n", get<0>(cv[i]), get<1>(cv[i]));
//	}


	fk = kfile + "_1200_CORE_FINAL_K";
	fp = fopen (fk.c_str(), "r");
	vector<tuple<int, int>> cv;
	i = 0;
	while (fscanf (fp, "%d", &(K[i])) != EOF) {
		cv.push_back (make_tuple(i, K[i]));
		i++;
	}
	sort (cv.begin(), cv.end(), kksort);


	vector<vertex> sample_vertices;

//	for (int i = 0; i <= TEST_SIZE; i++) {
//		printf ("deg of %d is %d\n", get<0>(dv[i]), get<1>(dv[i]));
//		sample_vertices.push_back (get<0>(dv[i]));
//	}


	srand(time(NULL));

	for (int i = 1; i <= TEST_SIZE; i++) {
		int u = rand() % (i * nVtx / TEST_SIZE);
		sample_vertices.push_back (get<0>(dv[u]));
		sample_vertices.push_back (get<0>(cv[u]));
		sample_vertices.push_back (get<0>(dv[i-1]));
		sample_vertices.push_back (get<0>(cv[i-1]));
	}

	vertex* P = (vertex *) calloc (nVtx, sizeof(vertex));
	int changed[nVtx];
	for (int i = 0; i < nVtx; i++)
		changed[i] = -1;

	for (auto v : sample_vertices) {

 		memset (P, 0, nVtx);
 		const auto t_begin = chrono::steady_clock::now();

		changed[v] = 1;
		for (int i = xadj[v]; i < xadj[v+1]; i++)
			changed[adj[i]] = 1;
		partialnmLocal12 (changed, nVtx, adj, xadj, P);

		const auto t_end = chrono::steady_clock::now();
		tms t1 = t_end - t_begin;
		printf ("Local convergence time: %.6lf secs\n", t1.count());

		int pv = P[v], kv = K[v];
//		int point = 0, base = 0;
		vector<vertex> gt, rs;
		for (int i = xadj[v]; i < xadj[v+1]; i++) {
			vertex w = adj[i];
			if (K[w] >= kv)
				gt.push_back (w);

			if (P[w] >= pv)
				rs.push_back (w);

//			if ((K[w] >= kv && P[w] >= pv) || (K[w] < kv && P[w] < pv))
//				point++;
//			if (K[w] >= kv) {
//					printf ("%d is neig with K: %d P : %d\t\t", w, K[w], P[w]);
//				base++;
//				if (P[w] >= pv) {
//						printf ("GOT THE POINT\n");
//					point++;
//				}
//				else {
//						printf ("\n");
//				}
//			}
		}

		printf ("12_EGO -- P: %d  vs  K: %d\t", P[v], K[v]);
//		printf ("score: %d / %d = %lf\t", point, base, (double) point / base);
		sort (rs.begin(), rs.end());
		sort (gt.begin(), gt.end());

		int is = commons (rs, gt);
		double prec = (double) is / rs.size(), recall = (double) is / gt.size();
		double f1 = 2 * prec * recall / (prec + recall), jaccard = (double) is / (gt.size() + rs.size() - is);
		printf ("base: %d\tjac: %lf\tprec: %lf\trec: %lf\tf1: %lf\t\t", gt.size(), jaccard, prec, recall, f1);

		{
			for (vertex i = 0; i < nVtx; i++)
				P[i] = D[i];

			int pv = P[v], kv = K[v];
//			int point = 0, base = 0;
			vector<vertex> gt, rs;
//			printf ("for vertex %d with K: %d and P(D): %d  \n", v, kv, pv);
			for (int i = xadj[v]; i < xadj[v+1]; i++) {
				vertex w = adj[i];

				if (K[w] >= kv)
					gt.push_back (w);

				if (P[w] >= pv)
					rs.push_back (w);
//				if ((K[w] >= kv && P[w] >= pv) || (K[w] < kv && P[w] < pv)) {
//					printf ("GOT THE POINT\n");
//					point++;
//				}

//				if (K[w] >= kv) {
////					printf ("%d is neig with K: %d P : %d\t\t", w, K[w], P[w]);
//					base++;
//					if (P[w] >= pv) {
////						printf ("GOT THE POINT\n");
//						point++;
//					}
//					else {
////						printf ("\n");
//					}
//				}
			}

			printf ("deg_EGO -- P: %d  vs  K: %d\t", P[v], K[v]);
//			printf ("score: %d / %d = %lf\t", point, base, (double) point / base);
			sort (rs.begin(), rs.end());
			sort (gt.begin(), gt.end());

			int is = commons (rs, gt);
			double prec = (double) is / rs.size(), recall = (double) is / gt.size();
			double f1 = 2 * prec * recall / (prec + recall), jaccard = (double) is / (gt.size() + rs.size() - is);
			printf ("base: %d\tjac: %lf\tprec: %lf\trec: %lf\tf1: %lf\n", gt.size(), jaccard, prec, recall, f1);

		}


	}


}
























