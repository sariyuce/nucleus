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

	nEdge = xadj[nVtx] / 2;

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

	timestamp totaltime (0, 0);
	timestamp totaltime_1;
	vertex maxP;
	vertex* P; // = (vertex *) malloc (sizeof(vertex) * nVtx);

	if (depth == 121) {
		vfile += "_K_CORE";
		kcore (nVtx, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 120) {
		vfile += "_CORE";
		baseLocal12 (nVtx, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 1200) {
		vfile += "_CORE";
		nmLocal12 (nVtx, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 712) {
		vfile += "_LEVELS";
		vertex* L;
		kcore_levels (nVtx, adj, xadj, L, vfile.c_str());
	}
	else if (depth == 612) {
		if (argc < 5) { // todo: update to 4
			cout << "K file is needed\n";
			exit(1);
		}
		read_Ks (nVtx, argv[4], &P); // todo: update to 3
		vfile += "_SESH_LEVELS";
		vertex* L;
		kcore_Sesh_levels (nVtx, adj, xadj, P, L, vfile.c_str());
	}
	else if (depth == -12) {
		if (argc < 5) { // todo: update to 4
			cout << "Top k number is needed\n";
			exit(1);
		}
		vertex top_k = atoi(argv[4]);
		fast12DegeneracyNumber (nVtx, adj, xadj, P, top_k);
	}











	if (depth == 231) {
		vfile += "_K_TRUSS";
		ktruss (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 230) {
		vfile += "_TRUSS";
		baseLocal23 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 2300) {
		vfile += "_TRUSS";
		nmLocal23 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 723) {
		vfile += "_LEVELS";
		vertex* L;
		ktruss_levels (nVtx, nEdge, adj, xadj, L, vfile.c_str());
	}
	else if (depth == 623) {
		if (argc < 5) { // todo: update to 4
			cout << "T file is needed\n";
			exit(1);
		}
		read_Ks (nEdge, argv[4], &P); // todo: update to 3
		vfile += "_SESH_LEVELS";
		vertex* L;
		ktruss_Sesh_levels (nVtx, nEdge, adj, xadj, P, L, vfile.c_str());
	}
	else if (depth == -23) {
		if (argc < 5) { // todo: update to 4
			cout << "Top k number is needed\n";
			exit(1);
		}
		vertex top_k = atoi(argv[4]);
		fast23DegeneracyNumber (nVtx, nEdge, adj, xadj, P, top_k);
	}


	else if (depth == 232) {
		vfile += "_STORE_TRI_K_TRUSS";
		ktruss_ST (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 23000) {
		vfile += "_ST_TRUSS";
		baseLocal23_ST (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 230000) {
		vfile += "_ST_TRUSS";
		nmLocal23_ST (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}

















	if (depth == 341) {
		vfile += "_34";
		k34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 340) {
		vfile += "_34";
		baseLocal34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 3400) {
		vfile += "_34";
		nmLocal34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 734) {
		vfile += "_LEVELS";
		vertex* L;
		k34_levels (nVtx, nEdge, adj, xadj, L, vfile.c_str());
	}
	else if (depth == -34) {
		if (argc < 5) { // todo: update to 4
			cout << "Top k number is needed\n";
			exit(1);
		}
		vertex top_k = atoi(argv[4]);
		fast34DegeneracyNumber (nVtx, nEdge, adj, xadj, P, top_k);
	}


	else if (depth == 342) {
		vfile += "_STORE_4C_34";
		k34_SF (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 34000) {
		vfile += "_SF_34";
		baseLocal34_SF (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}
	else if (depth == 340000) {
		vfile += "_SF_34";
		nmLocal34_SF (nVtx, nEdge, adj, xadj, P, vfile.c_str());
	}















#ifdef DEBUG_0
	if (depth == 340)
		vReals.clear();
	else if (depth == 120 || depth == 230)
		free (Reals);
#endif

	return 0;
}
