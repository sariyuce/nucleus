#include "main.h"

int main(int argc, char *argv[]) {

	vertex nVtx, *adj;
	edge nEdge, *xadj;

	if (argc < 2) {
		fprintf(stderr, "usage: %s "
				"\n <filename>"
				"\n <mode>\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	int idx = tmp.find_last_of("/");
	string asdf = tmp.substr(idx+1);

	int depth = atoi (argv[2]);
	string vfile = asdf + "_" + argv[2];

	// dummy, just for backward compatibility
	int order_choice;
	if (depth % 10 == 0) {
		order_choice = atoi (argv[3]);
		vfile = vfile + "_" + argv[3];
	}

	// read the graph
	Graph graph;
	readGraph<vertex, edge>(filename, &nVtx, &nEdge, &adj, &xadj);

	nEdge = xadj[nVtx];

	printf ("|V|: %d   |E|: %d\n", nVtx, nEdge);


	string out_file;
	FILE* fp;


	vertex* Reals;
	vector<int> vReals;
#ifdef DEBUG_0
	if (depth == 340 || depth == 3400) {

		FILE* aa = fopen (argv[4], "r");
		int nn, ii = 0;
		while (fscanf (aa, "%d", &nn) != EOF)
			vReals.push_back ((nn == -1) ? 0 : nn);
		fclose (aa);
	}
	else if (depth == 120 || depth == 1200 || depth == 230 || depth == 2300) {
		int sz;
		if (depth == 120 || depth == 1200)
			sz = nVtx;
		else if (depth == 230 || depth == 2300)
			sz = nEdge / 2;


		Reals = (vertex *) malloc (sizeof(vertex) * sz);
		FILE* aa = fopen (argv[4], "r");
		int nn, ii = 0;
		while (fscanf (aa, "%d", &nn) != EOF)
			Reals[ii++] = (nn == -1) ? 0 : nn;
		fclose (aa);
	}
#endif


	string step_no = "";
	if (depth < -900 || depth == 2300 || depth == 3400) {
		out_file= vfile + "_out";
		fopen (out_file.c_str(), "w");
	}
	else {
		string str (argv[3]);
		step_no = str.substr(str.find_last_of("_") + 1);
		if (step_no == "values")
			step_no = "K";
		else {
			string ss = str.substr(str.find_first_of("/") + 1, 5);
			if (ss == "ASYNC")
				step_no = "ASYNC_H_" + step_no;
			else
				step_no = "H_" + step_no;
		}

	}

	timestamp totaltime (0, 0);
	timestamp totaltime_1;
	vertex maxP;

	// ./picore graph depth KvalueFile
	if (depth == 912 ) { // read K values for kcore and construct hierarchy with subgraphs
		vertex* P = (vertex *) malloc (sizeof(vertex) * nVtx);
		HashMap<bool> exists (false);
		FILE* fzp = fopen (argv[3], "r");
		int cn = 0, num;
		for (int i = 0; i < nVtx; i++) {
			fscanf (fzp, "%d", &num);
			if (num < 0 || num > nVtx)
				P[i] = 0;
			else
				P[i] = num;

			if (P[i] > cn)
				cn = P[i];

			if (P[i] > 0)
				exists[P[i]] = true;

		}
		fclose (fzp);

		p_auxies ax;
		vfile = asdf + "_" + argv[2] + "_HIER_" + step_no;
		cout << "vfile: " << vfile << endl;
		string out_file = vfile + "_out";
		FILE* ffp = fopen (out_file.c_str(), "w");
		report_all_stuff2 (12, graph, nEdge, cn, ax, exists, P, vfile.c_str(), ffp);
	}
	else if (depth == 923) { // read K values for ktruss and construct hierarchy with subgraphs
		vertex* P = (vertex *) malloc (sizeof(vertex) * nEdge);
		HashMap<bool> exists (false);
		FILE* fzp = fopen (argv[3], "r");
		int cn = 0, num;
		for (int i = 0; i < nEdge; i++) {
			fscanf (fzp, "%d", &num);
			if (num < 0 || num > nEdge)
				P[i] = 0;
			else
				P[i] = num;
//			printf ("P[%d]: %d\n", i, P[i]);

			if (P[i] > cn)
				cn = P[i];

			if (P[i] > 0)
				exists[P[i]] = true;

		}
		fclose (fzp);


		EdgeList2 el;
		vector<vertex> xel;
		Graph og;
		create_ordered_graph (graph, el, xel, og);
		p_auxies ax;
		ax.xel = &xel;
		ax.el = &el;


		Graph ordered_graph;

		vfile = asdf + "_" + argv[2] + "_HIER_" + step_no;
		cout << "vfile: " << vfile << endl;
		string out_file = vfile + "_out";
		cout << "outfile: " << out_file << endl;
		FILE* ffp = fopen (out_file.c_str(), "w");
		report_all_stuff2 (23, graph, nEdge, cn, ax, exists, P, vfile.c_str(), ffp);

	}


	vertex* P; // = (vertex *) malloc (sizeof(vertex) * nVtx);
	if (depth == 712) {
		vfile += "_asdf";
		kcore_levels (graph, nEdge/2, P, totaltime, &maxP, vfile.c_str(), fp);
	}
	else if (depth == 612) {
		vfile += "_asdf";
		vertex* cono = (vertex *) malloc (sizeof(vertex) * nVtx);
		FILE* fzp = fopen (argv[3], "r");
		int num;
		for (int i = 0; i < nVtx; i++) {
			fscanf (fzp, "%d", &num);
			if (num < 0 || num > nVtx)
				cono[i] = 0;
			else
				cono[i] = num;
		}
		fclose (fzp);
		kcore_Seshlevels (graph, nEdge/2, P, cono, totaltime, &maxP, vfile.c_str(), fp);
	}



	else if (depth == 19) {
		vfile += "_PI_CORE";
		fast12DegeneracyNumber (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}


	else if (depth == 120) {
		vfile += "_PI_CORE";

		baseLocal12 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 1200) {
		vfile += "_PI_CORE";
		nmLocal12 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 12000) {
		vfile += "_PI_CORE";
		NoWaitnmLocal12 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
//		TRY_NMpicore (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 121) {
		vfile += "_K_CORE";
		base_kcore (graph, nEdge/2, P, totaltime, &maxP, vfile.c_str(), fp);
	}


	else if (depth == 723) {
		vfile += "_K_TRUSS";
		vector<vertex> T;
		ktruss_levels (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}
	else if (depth == 623) {
		vector<vertex> T;
		vertex* cono = (vertex *) malloc (sizeof(vertex) * nEdge);
		FILE* fzp = fopen (argv[3], "r");
		int num;
		int i = 0;
		while (	fscanf (fzp, "%d", &num) != EOF) {
			if (num < 0)
				cono[i] = 0;
			else
				cono[i] = num;
			i++;
		}
		printf ("i : %d  nEdge: %d\n", i, nEdge);
		fclose (fzp);
		ktruss_Seshlevels (graph, nEdge/2, T, cono, totaltime, &maxP, vfile.c_str(), fp);

	}



	else if (depth == 230) {
		vfile += "_PI_TRUSS";
		baseLocal23 (nVtx, nEdge/2, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 2300) {
		vfile += "_PI_TRUSS";
		nmLocal23 (nVtx, nEdge/2, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 23000) {
		vfile += "_PI_TRUSS";
		baseLocal23_ST (nVtx, nEdge/2, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 230000) {
		vfile += "_PI_TRUSS";
		nmLocal23_ST (nVtx, nEdge/2, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 2300000) {
		vfile += "_PI_TRUSS";
		NoWaitnmLocal23 (nVtx, nEdge/2, adj, xadj, P, vfile.c_str());
	}

	else if (depth == 231) {
		vfile += "_K_TRUSS";
		vector<vertex> T;
		base_ktruss (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}
	else if (depth == 232) {
		vfile += "_STORE_TRI_K_TRUSS";
		vector<vertex> T;
		base_ktruss_StoreTri (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}

	else if (depth == 634) {
		vector<vertex> T;
		vector<vertex> cono;
		FILE* fzp = fopen (argv[3], "r");
		int num;
		while (	fscanf (fzp, "%d", &num) != EOF) {
			if (num < 0 || num > nVtx)
				cono.push_back (0);
			else
				cono.push_back (num);
		}
		fclose (fzp);

		k34_Seshlevels (graph, nEdge/2, T, cono, totaltime, &maxP, vfile.c_str(), fp);
	}

	else if (depth == 734) {
		vfile += "_K_34";
		vector<vertex> T;
		k34_levels (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}
	else if (depth == 340) {
		vfile += "_PI_34";
//		base_pi34 (graph, nEdge/2, P, vReals, totaltime, &maxP, vfile.c_str(), fp);
		justpi34 (graph, nEdge/2, P, vReals, totaltime, &maxP, vfile.c_str(), fp, asdf);
	}
	else if (depth == 3400) {
		vertex* P;
		vfile += "_PI_34_LS";
		base_pi34_LessSpace (graph, nEdge/2, P, vReals, totaltime, &maxP, vfile.c_str(), fp, asdf);
//		TryFasterPi34 (graph, nEdge/2, P, vReals, totaltime, &maxP, vfile.c_str(), fp, asdf);
	}
	else if (depth == 341) {
		vfile += "_K_34";
		vector<vertex> T;
		base_k34 (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}
	else if (depth == 342) {
		vfile += "_STORE_4C_K_34";
		vector<vertex> T;
		base_k34_Store4c (graph, nEdge/2, T, totaltime, &maxP, vfile.c_str(), fp);
	}

//	if (depth / 10 == 12) {
//		// for each edge, check direction and difference
//		for (vertex i = 0; i < nVtx; i++)
//			for (vertex j = 0; j < graph[i].size(); j++) {
//				int u = i, v = graph[i][j];
//				if (u < v)
//					printf ("u,v: %d\t-\t%d\tPs: %d\t-\t%d\n", u, v, P[u], P[v]);
//			}
//	}
//	else {
		// for each triangle, check direction and difference
//		for (vertex i = 0; i < nVtx; i++)
//			for (vertex j = 0; j < graph[i].size(); j++) {
//				int u = i, v = graph[i][j];
//				if (u < v) {
//					vector<vertex> g;
//					smart_intersection (graph, u, v, g);
//					for (int w : g)
//						printf ("u,v,w: %d\t-\t%d\t-\t%d\tPs: %d\t-\t%d\t-\t%d\n", u, v, w, P[u], P[v]);
//				}
//	}
//	if (depth < 900)
//		fclose (fp);

//	for (vertex i = 0; i < nVtx; ++i)
//		graph[i].clear();
//	graph.clear();

#ifdef DEBUG_0
	if (depth == 340)
		vReals.clear();
	else if (depth == 120 || depth == 230)
		free (Reals);
#endif

	return 0;
}
