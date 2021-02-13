#include "main.h"
#define MAXLINE 1000000
#ifdef SIGNS
	Graph signs;
	int MODE;
#endif


int main (int argc, char *argv[]) {

	const auto t1 = chrono::steady_clock::now();
	if (argc < 4) {
		fprintf(stderr, "usage: %s "
				"\n filename"
				"\n motif: \n"
				"cycle\n"
				"acyclic\n"
				"acyclic-RA\n"
				"outp\n"
				"cyclep\n"
				"inp\n"
				"cyclepp\n"
				"hierarchy: YES or NO\n"
				, argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);

	// read the graph, give sorted edges in graph
	edge nEdge = 0;

	Graph rawgraph;
	readDirectedGraph<vertex, edge> (filename, rawgraph, &nEdge);

	string hrc (argv[3]);
#ifdef SIGNS
	string sm (argv[4]);
	string vfile = gname + "_" + nd + "_" + sm;
	MODE = atoi (argv[4]);
#else
	string vfile = gname + "_" + nd;
#endif

	string out_file;

	bool hierarchy = (hrc == "YES" ? true : false);
	if (hierarchy)
		out_file = vfile + "_Hierarchy";
	else
		out_file = vfile + "_K";


	FILE* fp = fopen (out_file.c_str(), "w");

	vertex maxK; // maximum K value in the graph
	vector<vertex> K;


	// to read subgraph/quark/cluster lists and compute avg. motif degree and conductance metrics
	if (COUNT_ONLY) {
		printf ("asdfa\n");
		string line;
		ifstream myfile(argv[4]);
		int num;
		while (getline( myfile, line )) {
			Graph graph (rawgraph);
			unordered_map<int, bool> numbers;
			istringstream iss(line);
			while (iss >> num) {
				numbers.emplace(num, true);
			}

			Graph boundary (rawgraph);
			unordered_map<int, bool> crossing;
			// find the ones near the crossing line
			for (int i = 0; i < boundary.size(); i++) {
				if (!boundary[i].empty() && numbers.find(i) != numbers.end()) {
					vector<int> newneigs;
					int endOfOut = boundary[i][0];
					// go over outgoing edges
					for (int j = 1; j < endOfOut; j++) {
						vertex v = boundary[i][j];
						if (numbers.find(v) == numbers.end())
							crossing.emplace (v, true);
					}

					//go over incoming edges
					for (int j = endOfOut; j < boundary[i].size(); j++) {
						vertex v = M2P(boundary[i][j]);
						if (numbers.find(v) == numbers.end())
							crossing.emplace (v, true);
					}
				}
			}


			unordered_map<int, bool> total (numbers);
			for (auto& x: crossing)
				total.emplace (x.first, true);


			int totalE = 0;
			// first ditch all the adj lists of non-existing nodes
			for (int i = 0; i < boundary.size(); i++) {
				if (total.find(i) == total.end())
					boundary[i].clear();
			}

			// then ditch the non-existing nodes from adj lists of existing nodes
			for (int i = 0; i < boundary.size(); i++) {
				if (!boundary[i].empty()) {
					vector<int> newneigs;
					int endOfOut = boundary[i][0];
					// go over outgoing edges
					int newEndOfOut = 1;
					for (int j = 1; j < endOfOut; j++) {
						vertex v = boundary[i][j];
						if (total.find(v) != total.end()) {
							newneigs.push_back (boundary[i][j]);
							newEndOfOut++;
						}
					}

					//go over incoming edges
					for (int j = endOfOut; j < boundary[i].size(); j++) {
						vertex v = M2P(boundary[i][j]);
						if (total.find(v) != total.end())
							newneigs.push_back (boundary[i][j]);
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

			if (nd == "cycle")
				cond = simple_count_cycles (boundary, numbers, crossing);
			else if (nd == "acyclic")
				cond = simple_count_acyclics (boundary, numbers, crossing);
			else if (nd == "outp")
				cond = simple_count_outps (boundary, numbers, crossing);
			else if (nd == "cyclep")
				cond = simple_count_cycleps (boundary, numbers, crossing);
			else if (nd == "inp")
				cond = simple_count_inps (boundary, numbers, crossing);
			else if (nd == "cyclepp")
				cond = simple_count_cyclepps (boundary, numbers, crossing);

			printf ("cond: %lf\n", cond);
		}
	}
	else {
		Graph graph (rawgraph);

		if (nd == "cycle")
			cycle_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "acyclic")
			acyclic_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "acyclic-RA")
			acyclic_truss_roleAware (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "outp")
			outp_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "cyclep")
			cyclep_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "inp")
			inp_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		else if (nd == "cyclepp")
			cyclepp_truss (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
	}

	const auto t2 = chrono::steady_clock::now();

	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s: %d\n", gname.c_str(), rawgraph.size(), nEdge, nd.c_str(), maxK);

	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);

#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%d\n", K[i]);
	fclose (kf);
#endif

	return 0;
}



