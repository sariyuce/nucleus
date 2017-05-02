#include "main.h"

int main (int argc, char *argv[]) {

	if (argc < 3) {
		fprintf(stderr, "usage: %s "
				"\n <filename>"
				"\n <nucleus type: 12, 13, 14, 23, 24, 34>"
				"\n <hierarchy?: 1 yes, 0 no>\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	int nd = atoi (argv[2]);
	string nds (argv[2]);
	bool hierarchy = atoi (argv[3]) == 1 ? true : false;

	// read the graph
	edge nEdge;
	Graph graph;
	ReadGraph<vertex, edge> (filename, graph, &nEdge);

	string vfile = gname + "_" + nds;
	string out_file;

	if (hierarchy)
		out_file = vfile + "-Hierarchy";
	else
		out_file = vfile + "-K";

	FILE* fp = fopen (out_file.c_str(), "w");

	timestamp peelingTime (0, 0);

	timestamp t1;

	vertex maxK; // maximum K value in the graph
	vector<vertex> K;

	switch (nd) {
	case 12: base_kcore (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	case 13: base_k13 (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	case 14: base_k14 (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	case 23: base_ktruss (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	case 24: base_k24 (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	case 34: base_k34 (graph, hierarchy, nEdge/2, K, &maxK, vfile, fp); break;
	default: cout << "Invalid nd: " << nd << endl; exit (1);
	}

	timestamp t2;

	printf ("%s\t|V|: %d\t|E|: %d\tmaxK: %d\n", gname.c_str(), graph.size(), nEdge, maxK);
	cout << t2 - t1 << endl;
	fprintf (fp, "%s\t|V|: %d\t|E|: %d\tmaxK: %d\n", gname.c_str(), graph.size(), nEdge, maxK);
	print_time (fp, "Total Time: ", t2 - t1);
	fclose (fp);

	return 0;
}
