#include "main.h"

int main (int argc, char *argv[]) {

	const auto t1 = chrono::steady_clock::now();
	if (argc < 3) {
		fprintf(stderr, "usage: %s "
				"\n <filename>"
				"\n 23"
				"\n <hierarchy?: YES or NO>\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);

	// read the graph, give sorted edges in graph
	edge nEdge = 0;
	Graph rawgraph;
	readGraph<vertex, edge> (filename, rawgraph, &nEdge);
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
	Graph graph (rawgraph);

	// read quark/subgraph/cluster lists, construct induced subgraphs, compute avg. motif degree and conductance metrics
	if (COUNT_ONLY) {
		string line;
		ifstream myfile(argv[4]);
		int num;
		while (getline( myfile, line )) {
			unordered_map<int, bool> numbers;
			istringstream iss(line);
			while (iss >> num)
				numbers.emplace(num, true);

			Graph boundary (rawgraph);
			unordered_map<int, bool> crossing;
			// find the ones near the crossing line
			for (int i = 0; i < boundary.size(); i++) {
				if (!boundary[i].empty() && numbers.find(i) != numbers.end()) {
					// go over edges
					for (int j = 0; j < boundary[i].size(); j++) {
						vertex v = boundary[i][j];
						if (numbers.find(v) == numbers.end()) {
							crossing.emplace (v, true);
						}
					}
				}
			}

			unordered_map<int, bool> total (numbers);
			for (auto& x: crossing)
				total.emplace (x.first, true);


			int totalE = 0;
			// first ditch all the adj lists of non-existing nodes
			for (int i = 0; i < boundary.size(); i++)
				if (total.find(i) == total.end())
					boundary[i].clear();


			// then ditch the non-existing nodes from adj lists of existing nodes
			for (int i = 0; i < boundary.size(); i++) {
				if (!boundary[i].empty()) {
					vector<int> newneigs;
					int newEndOfOut = 1;
					// go over edges
					for (int j = 0; j < boundary[i].size(); j++) {
						vertex v = boundary[i][j];
						if (total.find(v) != total.end()) {
							newneigs.push_back (boundary[i][j]);
							newEndOfOut++;
						}
					}

					totalE += newneigs.size();
					newneigs.insert (newneigs.begin(), newEndOfOut);
					boundary[i].assign (newneigs.begin(), newneigs.end());
				}
			}


			if (total.size() != numbers.size() + crossing.size())
				printf ("sth is wrong: %d + %d != %d\n", numbers.size(), crossing.size(), total.size());
			printf ("total \n");
			printf ("|V|: %d\t |E|: %d\n", total.size(), totalE);
			double cond;


			if (nd == "23") {
				vertex nVtx = boundary.size();

				// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
				vector<vp> el;
				vector<vertex> xel;
				Graph orderedGraph;
				createOrderedIndexEdges (boundary, el, xel, orderedGraph);

				// Triangle counting for each edge
				vector<vector<vertex> > TC (nVtx);
				for (vertex i = 0; i < nVtx; i++)
					TC[i].resize (orderedGraph[i].size(), 0);

				simple_count_triangles (boundary, orderedGraph, TC, numbers, crossing);
			}
		}
	}
	else {
		if (nd == "23") {
			base_ktruss (graph, hierarchy, nEdge, K, &maxK, vfile, fp);
	}

#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%lld\n", K[i]);
	fclose (kf);
#endif

	const auto t2 = chrono::steady_clock::now();
	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s-nucleus: %d\n", gname.c_str(), rawgraph.size(), nEdge, nd.c_str(), maxK);
	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);
	return 0;
}
