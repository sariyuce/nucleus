#include "main.h"
#define MAXLINE 1000000
#ifdef SIGNS
	Graph signs;
	int MODE;
#endif

#ifdef ORBITS
	unordered_map<couple, int> orb1, orb2, orb3;
#endif

int main (int argc, char *argv[]) {

//	char* line = (char*) malloc (sizeof (char) * MAXLINE);
//	FILE* matfp = fopen(argv[1], "r");
//
//	unordered_map<string, int> mp;
//	int id = 0;
//
//	string u, v;
//	while (fgets(line, MAXLINE, matfp)) {
//		stringstream ss (line);
//		ss >> u >> v;
//		if (mp.find (u) == mp.end())
//			mp[u] = id++;
//		if (mp.find (v) == mp.end())
//			mp[v] = id++;
//		printf ("ids: %d %d\n", mp[u], mp[v]);
//	}
//	fclose(matfp);
//
//	for (auto& x: mp) {
//	    std::cout << "map " <<x.first << " " << x.second << std::endl;
//	  }
//
//	return 0;

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
//
//	if (!(
//			nd == "cycle-truss" ||
//			nd == "cycle-core" ||
//			nd == "acyclic-truss" ||
//			nd == "acyclic-core" ||
//			nd == "outp-truss" ||
//			nd == "outp-core" ||
//			nd == "cyclep-truss" ||
//			nd == "cyclep-core" ||
//			nd == "inp-truss" ||
//			nd == "inp-core" ||
//			nd == "cyclepp-truss" ||
//			nd == "cyclepp-core")) {
//		printf ("Invalid algorithm, options are cycle\n"
//				"acyclic\n"
//				"out-p\n"
//				"cycle-p\n"
//				"in-p\n"
//				"cycle-pp\n");
//		exit(1);
//	}

	// read the graph, give sorted edges in graph
	edge nEdge = 0;
//	Graph graph;
//	readDirectedGraph<vertex, edge> (filename, graph, &nEdge);





					Graph rawgraph;
					readDirectedGraph<vertex, edge> (filename, rawgraph, &nEdge);
#ifdef SIGNS
	if (0){
		for (int i = 0; i < rawgraph.size(); i++) {
			vector<vertex> ret;
                	outgoings (rawgraph[i], ret);
                	for (vertex r = 0; r < ret.size(); r++) {
                        	vertex v = rawgraph[i][ret[r]];
				vertex s = signs[i][ind (v, rawgraph[i])];
				printf ("edge %d %d %d\n", i, v, s);
			}
		}

		for (int i = 0; i < rawgraph.size(); i++) {
                        vector<vertex> ret;
                        incomings (rawgraph[i], ret);
                        for (vertex r = 0; r < ret.size(); r++) {
                                vertex v = M2P (rawgraph[i][ret[r]]);
                                vertex s = signs[v][ind (i, rawgraph[v])];
                                printf ("rev %d %d %d\n", v, i, s);
                        }
                }
		exit (1);

	}
#endif





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


	// READ CLUSTER LISTS, CONSTRUCT INDUCED SUBGRAPHS

	if (COUNT_ONLY) {
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

			if (nd == "cycle-truss")
				cond = simple_count_cycles (boundary, numbers, crossing);
			else if (nd == "acyclic-truss")
				cond = simple_count_acyclics (boundary, numbers, crossing);
			else if (nd == "outp-truss")
				cond = simple_count_outps (boundary, numbers, crossing);
			else if (nd == "cyclep-truss")
				cond = simple_count_cycleps (boundary, numbers, crossing);
			else if (nd == "inp-truss")
				cond = simple_count_inps (boundary, numbers, crossing);
			else if (nd == "cyclepp-truss")
				cond = simple_count_cyclepps (boundary, numbers, crossing);

			printf ("cond: %lf\n", cond);
		}
	}
	else {
		Graph graph (rawgraph);
//		if (nd == "cycle-core")
//				cycle_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		cycle_core (graph, hierarchy, nEdge, K, &maxK, fp);
//			else if (nd == "acyclic-core")
//				acyclic_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		acyclic_core (graph, hierarchy, nEdge, K, &maxK, fp);
//			else if (nd == "outp-core")
//				outp_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		outp_core (graph, hierarchy, nEdge, K, &maxK, fp);
//			else if (nd == "cyclep-core")
//				cyclep_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		cyclep_core (graph, hierarchy, nEdge, K, &maxK, fp);
//			else if (nd == "inp-core")
//				inp_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		inp_core (graph, hierarchy, nEdge, K, &maxK, fp);
//			else if (nd == "cyclepp-core")
//				cyclepp_core_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
//		//		cyclepp_core (graph, hierarchy, nEdge, K, &maxK, fp);

			if (nd == "cycle-truss")
				cycle_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		cycle_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "acyclic-truss")
				acyclic_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		acyclic_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "acyclic-truss-RA")
				acyclic_truss_SUBS_roleAware (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		acyclic_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "outp-truss")
				outp_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		outp_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "cyclep-truss")
				cyclep_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		cyclep_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "inp-truss")
				inp_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		inp_truss (graph, hierarchy, nEdge, K, &maxK, fp);
			else if (nd == "cyclepp-truss")
				cyclepp_truss_SUBS (graph, hierarchy, nEdge, K, &maxK, fp, vfile);
		//		cyclepp_truss (graph, hierarchy, nEdge, K, &maxK, fp);
	}

	const auto t2 = chrono::steady_clock::now();
//	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s: %d\n", gname.c_str(), graph.size(), nEdge, nd.c_str(), maxK);


			printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s: %d\n", gname.c_str(), rawgraph.size(), nEdge, nd.c_str(), maxK);










	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);


	if (DEG_DIST) {
		vector<int> dist (maxK+1, 0);
		for (vertex i = 0; i < K.size(); i++)
			if (K[i] >= 0)
				dist[K[i]]++;

		for (vertex i = 0; i < dist.size(); i++)
			if (dist[i] > 0)
				printf ("K %d %d\n", i, dist[i]);
	}


#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%d\n", K[i]);
	fclose (kf);
#endif

	return 0;
}

/*
 *
 * if (0) {
		vector<double> pr (graph.size(), 1.0/graph.size());
		for (int count = 0; count < 50; count++) {
			for (int i = 0; i < graph.size(); i++) {
				double sum = 0;
				if (graph[i].empty())
					continue;
				for (int j = graph[i][0]; j < graph[i].size(); j++) {
					int v = M2P (graph[i][j]);
					if (!graph[v].empty() && graph[v][0] >= 1) {
						int sz = graph[v][0] - 1;
						if (sz > 1)
							sum += pr[v]/(sz*sz-sz);
					}
				}
				pr[i] = 0.75 * sum + 0.25 * (1.0/graph.size());
			}
		}

		for (int i = 0; i < graph.size(); i++)
			printf ("pr of %d is %lf\n", i, pr[i]);
		exit (1);
	}

	if (0) {
		FILE* fp = fopen ("undirected", "w");
		for (auto u = 0; u < graph.size(); u++) {
			vector<vertex> ret;
			undirecteds (graph[u], ret);
			for (vertex r = 0; r < ret.size(); r++) {
				vertex v = graph[u][ret[r]];
				if (u < v)
					fprintf (fp, "%d %d\n", u, v);
			}
		}
		fclose (fp);
		exit (1);
	}

	if (0) {
		vertex outs = 0;
		vertex ins = 0;
		vertex recs = 0;

		vertex NE = 0;
		vertex naiveNE = 0;
		for (auto u = 0; u < graph.size(); u++) {
			if (graph[u].size() > 0)
				naiveNE += graph[u][0] - 1;
			vector<vertex> ret;
			outgoings (graph[u], ret);
			outs += ret.size();


			ret.clear();
			undirecteds (graph[u], ret);
			recs += ret.size();




		}
//		recs /= 2;
		NE = outs + recs;
		recs /= 2;
		printf ("#reciprocals: %d\t #edges: %d\t ratio: %lf\n", recs, recs+outs, (double) recs/(recs+outs));
		printf ("NE: %d\t nEdge: %d\n", NE, nEdge);
		printf ("naiveNE: %d\n", naiveNE);


//		for (auto v = 0; v < graph[u].size(); v++)
//			printf ("%d ", graph[u][v]);
//		printf ("\n");
//	}

	}
 *
 *
 */

