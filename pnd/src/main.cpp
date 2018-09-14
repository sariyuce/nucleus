#include "main.h"
bool connected (int u, int v, int* xadj, int* adj) {
	for (int i = xadj[v]; i < xadj[v+1]; i++)
		if (adj[i] == u)
			return true;
	return false;
}



int main(int argc, char *argv[]) {

	vertex nVtx, *adj;
	edge nEdge, *xadj;
	if (argc < 3) {
		fprintf(stderr, "usage: %s \n <filename> \n <mode (12, 120, 1200, 23, 230, 2300, 34, 340, 3400)>\n", argv[0]);
		exit(1);
	}
	string gr (argv[1]);
	int depth = atoi (argv[2]);
	string vfile = gr.substr (gr.find_last_of("/")+1) + "_" + argv[2] ;

	// read the graph
	Graph graph;
	readGraph<vertex, edge>(argv[1], &nVtx, &nEdge, &adj, &xadj);
	nEdge = xadj[nVtx] / 2;
	printf ("|V|: %d   |E|: %d\n", nVtx, nEdge);
	vertex* P;

	switch (depth) {
	case 12:
		vfile += "_K_CORE";
		kcore (nVtx, adj, xadj, P, vfile.c_str());
		break;
	case 120:
		vfile += "_CORE";
		baseLocal12 (nVtx, adj, xadj, P, vfile.c_str());
		break;
	case 1200:
		vfile += "_CORE";
		nmLocal12 (nVtx, adj, xadj, P, vfile.c_str());
		break;
	case 129:
		vfile += "_CORE";
		topKs (nVtx, adj, xadj, P, vfile.c_str());
		break;
	case 612:
		find_mcore (nVtx, adj, xadj, P, gr.substr (gr.find_last_of("/")+1));
		break;
	case 912:
		converge12onEgo (nVtx, adj, xadj, P, gr.substr (gr.find_last_of("/")+1));
		break;
	case 23:
		vfile += "_K_TRUSS";
		ktruss (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 230:
		vfile += "_TRUSS";
		baseLocal23 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 2300:
		vfile += "_TRUSS";
		nmLocal23 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 239:
		vfile += "_TRUSS";
		topKs23 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 623:
		find_mtruss (nVtx, nEdge, adj, xadj, P, gr.substr (gr.find_last_of("/")+1));
		break;
	case 923:
		converge23onEgo (nVtx, nEdge, adj, xadj, P, gr.substr (gr.find_last_of("/")+1));
		break;
	case 34:
		vfile += "_34";
		k34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 340:
		vfile += "_34";
		baseLocal34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	case 3400:
		vfile += "_34";
		nmLocal34 (nVtx, nEdge, adj, xadj, P, vfile.c_str());
		break;
	}
	return 0;
}
