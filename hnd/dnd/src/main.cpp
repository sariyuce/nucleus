#include "main.h"

int main (int argc, char *argv[]) {



//	vector<couple> coords;
//	coords.push_back (make_tuple(2, 5));
//	coords.push_back (make_tuple(3, 4));
//	coords.push_back (make_tuple(2, 4));
//	coords.push_back (make_tuple(1, 4));
//	coords.push_back (make_tuple(1, 2));
//	coords.push_back (make_tuple(2, 3));
//
//
//	sort (coords.begin(), coords.end());
//	for (auto i : coords)
//		printf ("%d %d\n", get<0>(i), get<1>(i));
//
//	exit (1);

	const auto t1 = chrono::steady_clock::now();
	if (argc < 2) {
		fprintf(stderr, "usage: %s "
				"\n filename"
				"\n alg: \n"
				"cycle\n"
				"acyclic\n"
				"out-p\n"
				"cycle-p\n"
				"in-p\n"
				"cycle-pp\n"
				, argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);

	if (!(
			nd == "cycle-truss" ||
			nd == "cycle-core" ||
			nd == "acyclic" ||
			nd == "out-p" ||
			nd == "cycle-p" ||
			nd == "in-p" ||
			nd == "cycle-pp")) {
		printf ("Invalid algorithm, options are cycle\n"
				"acyclic\n"
				"out-p\n"
				"cycle-p\n"
				"in-p\n"
				"cycle-pp\n");
		exit(1);
	}

	// read the graph, give sorted edges in graph
	edge nEdge = 0;
	Graph graph;
	readDirectedGraph<vertex, edge> (filename, graph, &nEdge);

//	for (auto u = 0; u < graph.size(); u++) {
//		printf ("%d: ", u);
//		for (auto v = 0; v < graph[u].size(); v++)
//			printf ("%d ", graph[u][v]);
//		printf ("\n");
//	}

	string hrc (argv[3]);
	string vfile = gname + "_" + nd;
	string out_file;

	bool hierarchy = (hrc == "YES" ? true : false);
	if (hierarchy)
		out_file = vfile + "_Hierarchy";
	else
		out_file = vfile + "_K";

	FILE* fp = fopen (out_file.c_str(), "w");

	vertex maxK; // maximum K value in the graph
	vector<vertex> K;

	if (nd == "cycle-truss")
		cycle_truss (graph, hierarchy, nEdge, K, &maxK, fp);
	if (nd == "cycle-core")
		cycle_core (graph, hierarchy, nEdge, K, &maxK, fp);

//	else if (nd == "acyclic")
//		acyclic_truss (graph, nEdge, K, &maxK, vfile.c_str());
//	else if (nd == "out-p")
//		out_p_truss (graph, nEdge, K, &maxK, vfile.c_str());
//	else if (nd == "cycle-p")
//		cycle_p_truss (graph, nEdge, K, &maxK, vfile.c_str());
//	else if (nd == "in-p")
//		in_p_truss (graph, nEdge, K, &maxK, vfile.c_str());
//	else if (nd == "cycle-pp")
//		cycle_p_p_truss (graph, nEdge, K, &maxK, vfile.c_str());


#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%lld\n", K[i]);
	fclose (kf);
#endif

	const auto t2 = chrono::steady_clock::now();
	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s: %d\n", gname.c_str(), graph.size(), nEdge, nd.c_str(), maxK);
	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);

	return 0;
}
