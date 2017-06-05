#include "main.h"

int main(int argc, char *argv[]) {

	if (argc < 4) {
		fprintf(stderr, "usage: %s\n filename\n algorithm: RIGHT_TIP, LEFT_TIP, WING\n hierarchy?: 1 yes, 0 no\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);  // RIGHT_TIP, LEFT_TIP, WING 910 for rightVertex based, 911 for leftVertex based, 92 for edge
	if (!(nd == "RIGHT_TIP" || nd == "LEFT_TIP" || nd == "WING")) {
		printf ("Invalid algorithm, options are RIGHT_TIP, LEFT_TIP and WING\n");
		exit(1);
	}

	// read the graph, gives sorted edges in left and rightGraph
	edge nEdge;
	Graph leftGraph, rightGraph;

	ReadBipartiteGraph<vertex, vertex> (filename, &nEdge, leftGraph, rightGraph);

	string vfile = gname + "_" + nd;
	string out_file;

	bool hierarchy = (atoi (argv[3]) == 1);
	if (hierarchy)
		out_file = vfile + "-Hierarchy";
	else
		out_file = vfile + "-K";

	FILE* fp = fopen (out_file.c_str(), "w");

	timestamp peelingTime (0, 0);

	timestamp t1;

	vertex bCount = 0;
	vertex maxK; // maximum K value in the graph
	vector<vertex> K;
	vector<vertex> xRight;
	vector<vp> el;

	if (nd == "RIGHT_TIP") {
//		insallahFasterTipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount);
		hopefullyFastertipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount);
//		 tipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
	}
	else if (nd == "LEFT_TIP") {
//		insallahFasterTipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
		hopefullyFastertipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
//		 tipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
	}
	else if (nd == "WING") { // todo: should be same when you swap leftGraph and rightGraph
		prefixSum (xRight, rightGraph, el);
		printf ("summed\n");
		insallahFasterwingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
//		wingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
	}
	timestamp t2;

	printf ("%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d nButterflies: %d\tTotal time: ", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	cout << t2 - t1 << endl;
	fprintf (fp, "%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d\tnButterflies: %d\t", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	print_time (fp, "Total Time: ", t2 - t1);
	fclose (fp);
//
	FILE* dd = fopen (vfile.c_str(), "w");
	for (int i = 0; i < K.size(); i++)
		fprintf (dd, "K[%d]: %d\n", i, K[i]);
	fclose (dd);
	return 0;
}

