#include "main.h"

inline void print_Hs (int nVtx, vertex* T, const char* vfile) {
	string st (vfile);
	st += "_K_values";
	FILE* pp = fopen (st.c_str(), "w");
	for (int i = 0; i < nVtx; i++)
		fprintf (pp, "%d\n", T[i]);
	fclose (pp);
}

void base_kcore (Graph& graph, int nEdge, vertex *K, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp) {

	util::timestamp t_begin;
	// Initial declarations
	char *p = (char*) &totaltime;
	util::timestamp *pout = (util::timestamp*) p;
	size_t nVtx = graph.size();

	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		if (graph[i].size() > max_degree)
			max_degree = graph[i].size();
	}

	K = (vertex *) calloc (nVtx, sizeof(vertex));
//	vector<int> K (graph.size(), 0);
	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize(max_degree, (graph.size()));
	for (size_t i = 0; i < nVtx; i++)
		if (graph[i].size() > 0) {
			na_bs.Insert (i, graph[i].size());
			printf ("deg[%d]: %d\n", i, graph[i].size());
		}
	vertex degree_of_u = 0;

	HashMap<bool> exists (false);

	int order = 0;
	while (1) {
		int u, value_of_u;
		int ret = na_bs.PopMin(&u, &value_of_u);
		if (ret == -1) // if the bucket is empty
			break;

//		printf ("order %d %d\n", order++, u);
		if (value_of_u == 0)
			continue;


		degree_of_u = K[u] = value_of_u;
		exists[value_of_u] = true;

		for (size_t j = 0; j < graph[u].size(); j++) { /* decrease the degree of the neighbors with greater degree */
			vertex v = graph[u][j];
			int curval = na_bs.CurrentValue(v);
			if (curval > degree_of_u)
				na_bs.DecVal(v);
		}
	}

	na_bs.Free();

	*core_number = degree_of_u; // degree_of_u is degree of last popped vertex

	util::timestamp td1;
	cout << "K time: " << td1 - t_begin << endl;

	cout << "max K: " << *core_number << endl;

	print_Hs (nVtx, K, vfile);

	return;
}







void kcore_levels (Graph& graph, int nEdge, vertex *K, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp) {

	util::timestamp t_begin;
	// Initial declarations
	char *p = (char*) &totaltime;
	util::timestamp *pout = (util::timestamp*) p;
	size_t nVtx = graph.size();

	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		if (graph[i].size() > max_degree)
			max_degree = graph[i].size();
	}

	K = (vertex *) malloc (nVtx * sizeof(vertex));
	// Peeling
	Naive_Bucket na_bs;
	na_bs.Initialize(max_degree, (graph.size()));
	for (size_t i = 0; i < nVtx; i++)
		if (graph[i].size() > 0)
			na_bs.Insert (i, graph[i].size());
	vertex degree_of_u = 0;

	int level = 1;
	while (1) {
		int u, value_of_u;
		int ret;

		int val;
		vector<int> mins;
		ret = na_bs.PopMin(&u, &value_of_u);
		if (ret == -1) // if the bucket is empty
			break;

		if (value_of_u == 0)
			continue;

		mins.push_back (u);
		while (1) {
			ret = na_bs.PopMin(&u, &val);
			if (ret == -1) // if the bucket is empty
				break;
			if (val > value_of_u) {
				na_bs.Insert (u, val);
				break;
			}
			mins.push_back (u);
		}

//		if (ret == -1) // if the bucket is empty
//			break;


		for (int i = 0; i < mins.size(); i++) {
			int v = mins[i];
			K[v] = value_of_u;
			for (size_t j = 0; j < graph[v].size(); j++) { /* decrease the degree of the neighbors with greater degree */
				vertex w = graph[v][j];
				int curval = na_bs.CurrentValue(w);
				if (curval > value_of_u)
					na_bs.DecVal(w);
			}
		}
		printf ("Level %d K: %d size: %d\n", level, value_of_u, mins.size());
		level++;
	}

	na_bs.Free();

	*core_number = degree_of_u; // degree_of_u is degree of last popped vertex

	util::timestamp td1;

#ifdef DEBUG_00
	print_Hs (nVtx, K, vfile);
#endif
	cout << "K time: " << td1 - t_begin << endl;

	return;
}











void kcore_Seshlevels (Graph& graph, int nEdge, vertex *K, vertex *cono, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp) {

	util::timestamp t_begin;
	// Initial declarations
	char *p = (char*) &totaltime;
	util::timestamp *pout = (util::timestamp*) p;
	size_t nVtx = graph.size();

	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		if (graph[i].size() > max_degree)
			max_degree = graph[i].size();
	}
	K = (vertex *) malloc (nVtx * sizeof(vertex));

	vertex* LN = (vertex *) malloc (sizeof (vertex) * nVtx);
	int* D = (vertex *) malloc (nVtx * sizeof(vertex));
	for (size_t i = 0; i < nVtx; i++)
		D[i] = graph[i].size();

	int level = 1;
	while (1) {

		vector<int> torem;
		for (int i = 0; i < nVtx; i++)
			if (D[i] != INT_MIN && D[i] <= cono[i])
				torem.push_back(i);

		for (int i = 0; i < torem.size(); i++) {
			int v = torem[i];
			int sma = 0, geq = 0, prevs = 0;
			for (size_t j = 0; j < graph[v].size(); j++) {
				vertex w = graph[v][j];
				if (D[w] != INT_MIN) {
					if (D[w] < D[v])
						sma++;
					else
						geq++;
				}
				else {
					if (cono[w] >= cono[v])
						prevs++;
				}
			}
			printf ("K: %d  D: %d geq: %d   sma: %d   prevs: %d\n", cono[v], D[v], geq, sma, prevs);
		}


		for (int i = 0; i < torem.size(); i++) {
			int v = torem[i];
			for (size_t j = 0; j < graph[v].size(); j++) {
				vertex w = graph[v][j];
				if (D[w] != INT_MIN) {
//					if (D[w] < D[v])
//						printf ("WTF\n");
					D[w]--;
				}
			}
			D[v] = INT_MIN;
			LN[v] = level;
		}

		if (torem.empty())
			break;
		printf ("Level %d   size: %d\n", level, torem.size());
		level++;
	}



	for (size_t u = 0; u < nVtx; u++) {
//		if (LN[u] == 1) {
			int Pgt = 0, Psm = 0, Peq = 0, Cgt = 0, Csm = 0, Ceq = 0, Fgt = 0, Fsm = 0, Feq = 0;
			int smallers = 0, equals = 0, greaters = 0;
			int kk = 0, tt = 0;
			for (size_t j = 0; j < graph[u].size(); j++) {
				int v = graph[u][j];

				if (cono[v] < cono[u]) {
					int df = LN[v] - LN[u];
					if (df > 0) {
						kk++;
						tt += df;
					}
				}

//				if (LN[v] < LN[u]) {
//					if (cono[v] < cono[u])
//						Psm++;
//					else if (cono[v] > cono[u])
//						Pgt++;
//					else if (cono[u] == cono[v])
//						Peq++;
//				}
//				else if (LN[v] == LN[u]) {
//					if (cono[v] < cono[u])
//						Csm++;
//					else if (cono[v] > cono[u])
//						Cgt++;
//					else if (cono[u] == cono[v])
//						Ceq++;
//				}
//				else if (LN[v] > LN[u]) {
//					if (cono[v] < cono[u])
//						Fsm++;
//					else if (cono[v] > cono[u])
//						Fgt++;
//					else if (cono[u] == cono[v])
//						Feq++;
//				}
//
//				if (graph[v].size() < graph[u].size())
//					smallers++;
//				else if (graph[v].size() == graph[u].size())
//					equals++;
//				else if (graph[v].size() > graph[u].size())
//					greaters++;
			}
			printf ("for %d   D: %d  K: %d  L: %d, tt/kk : %lf\n", u, graph[u].size(), cono[u], LN[u], (double) tt / kk);
//			printf ("for %d   D: %d  K: %d  L: %d      Psm: %d  Peq: %d  Pgt: %d     Csm: %d  Ceq: %d  Cgt: %d     Fsm: %d  Feq: %d  Fgt: %d\n",
//					u, graph[u].size(), cono[u], LN[u], Psm, Peq, Pgt, Csm, Ceq, Cgt, Fsm, Feq, Fgt);
//			printf ("eda %d D: %d  K: %d   smD: %d   eqD: %d   gtD: %d     smK: %d   eqK: %d   gtK: %d\n",
//					u, graph[u].size(), cono[u], smallers, equals, greaters, Psm+Csm+Fsm, Peq+Ceq+Feq, Pgt+Cgt+Fgt);
//		}
	}


//	for (size_t u = 0; u < nVtx; u++) {
//		int gt = 0, sm = 0, eq = 0;
//		for (size_t j = 0; j < graph[u].size(); j++) {
//			int v = graph[u][j];
//			if (LN[u] > LN[v])
//				gt++;
//			else if (LN[u] < LN[v])
//				sm++;
//			else if (LN[u] == LN[v])
//				eq++;
//		}
//		printf ("for %d   D: %d  K: %d  gLNs: %d  eqLNs: %d sLNs: %d\n", u, graph[u].size(), cono[u], gt, eq, sm);
//	}

	return;
}
