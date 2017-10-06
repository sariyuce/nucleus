#include "main.h"

int main(int argc, char *argv[]) {

//	vector<int> aa = {3,5,7,1,7,187,1,8,8,10};
//
//	for (auto i = 0; i < aa.size(); i++) {
//		printf ("i: %d\n", i);
//		printf ("aa[%d]: %d\n", i, aa[i]);
//	}
//
////	for (auto i = aa.size() - 1; i >= 0; i--) {
////		printf ("i: %d\n", i);
////		printf ("aa[%d]: %d\n", i, aa[i]);
////	}
//	exit(1);

//
//	{
//		FILE* fp = fopen (argv[1], "r");
//		FILE* gp = fopen (argv[2], "w");
//
//		vector<vector<int>> vv;
//		vv.resize (10);
//		int a, b, pb = -1;
//		int id = -1;
//		while (fscanf (fp, "%d %d", &a, &b) != EOF) {
//
//			if (b!= pb)
//				id++;
//
//			vv[id].push_back (a);
//
//			pb = b;
//		}
//
//		fclose (fp);
//		unordered_map<long long, bool> mm;
//		long long count = 0;
//		for (auto v : vv) {
//			for (int i = 0; i < v.size(); i++) {
//				for (int j = i+1; j < v.size(); j++) {
//					long long s = v[i];
//					long long b = v[j];
//					if (s > b)
//						swap (s, b);
//					long long aaa = (s * 2000000) + b;
//					if (mm.find (aaa) == mm.end()) {
//						count++;
//
//						if (count % 100000000 == 0)
//							printf ("count: %lld\n", count);
////						fprintf (gp, "%d %d\n", v[i], v[j]);
//						mm[aaa] = true;
//					}
//				}
//			}
//		}
//		fclose (gp);
//		printf ("count: %lld\n", count);
//		return 0;
//	}
//
//
//
//
//	long long aa = 4000000000000;
//	aa+=1;
//	printf ("aa: %lld\n", aa);

	if (argc < 4) {
		fprintf(stderr, "usage: %s\n filename\n algorithm: RIGHT_TIP, LEFT_TIP, WING\n hierarchy?: 1 yes, 0 no\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);  // RIGHT_TIP, LEFT_TIP, WING 910 for rightVertex based, 911 for leftVertex based, 92 for edge
	if (!(nd == "RIGHT_TIP" || nd == "LEFT_TIP" || nd == "WING" || nd.find ("BUILD") != string::npos || nd.find ("RUN") != string::npos || nd.find ("MEASURE") != string::npos) ) {
		printf ("Invalid algorithm, options are RIGHT_TIP, LEFT_TIP and WING for bipartite graphs\n");
		printf ("Or, RUN_WEIGHTED_CORE and RUN_(CORE|TRUSS) for weighted and unweighted projections\n");
		printf ("Creating those projections are is possible via BUILD_(UNW|W)_(LEFT|RIGHT) in weighted and unweighted ways from either side\n");
		printf ("Lastly, MEASURE_(CORE|TRUSS)_(LEFT|RIGHT) computes the densities of the bipartite nuclei which are induced on the nuclei of LEFT|RIGHT projections, given as the third argument\n");
		exit(1);
	}

	// read the graph, gives sorted edges in left and rightGraph
	edge nEdge = 0;
	Graph leftGraph, rightGraph;
	Wraph wg;
	Graph gr;
	if (nd == "RUN_CORE" || nd == "RUN_TRUSS" || nd == "RUN_WEIGHTED_CORE") {
//		ReadRegularGraph (filename, gr, &nEdge);
		ReadWeightedGraph (filename, wg, &nEdge);
		printf ("nVtx: %d   nEdge:%d\n", wg.size(), nEdge);
		if (nd == "RUN_CORE" || nd == "RUN_TRUSS") {
			gr.resize (wg.size());
			for (int i = 0; i < wg.size(); i++) {
				gr[i].resize (wg[i].size());
				for (int j = 0; j < wg[i].size(); j++) {
					gr[i][j] = wg[i][j].n;
				}
			}

			for (int i = 0; i < gr.size(); i++) {
//				printf ("%d: ", i);
				for (int j = 0; j < gr[i].size(); j++) {

					int v = gr[i][j];
//					printf ("%d ", v);
					bool df = false;
					for (int k = 0; k < gr[v].size(); k++) {
						if (gr[v][k] == i) {
							df = true;
							break;
						}
					}
					if (!df) {
						printf ("%d and %d are not reciprocal\n", i, v);
						exit(1);
					}
				}
//				printf ("\n");
			}



		}
		printf ("nVtx: %d   nEdge:%d, looks good\n", gr.size(), nEdge);
	}
	else
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

	long long bCount = 0;
	vertex maxK; // maximum K value in the graph
	vector<vertex> K;
	vector<vertex> xRight;
	vector<vp> el;


//	for (int i = 0; i < leftGraph.size(); i++)
//		for (int j = 0; j < leftGraph[i].size(); j++)
//			printf ("gr %d %d\n", i, leftGraph[i][j]);
//
//	exit(1);
	// for projections
//	Graph gr;

//	if (nd == "RIGHT_TIP") {
////		vector<long long> K;
////		long long maxK;
//		tipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
////		oldtipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
//		FILE* dd = fopen (vfile.c_str(), "w");
//		for (int i = 0; i < K.size(); i++)
//			fprintf (dd, "K[%d]: %lld\n", i, K[i]);
//		fclose (dd);
//	}
//	else if (nd == "LEFT_TIP") {
////		vector<long long> K;
////		long long maxK;
//		tipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
////		oldtipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
//		FILE* dd = fopen (vfile.c_str(), "w");
//		for (int i = 0; i < K.size(); i++)
//			fprintf (dd, "K[%d]: %lld\n", i, K[i]);
//		fclose (dd);
//	}
//	else if (nd == "WING") { // todo: should be same when you swap leftGraph and rightGraph
//		if (leftGraph.size() < rightGraph.size()) {
//			fprintf (fp, "LEFT is the PRIMARY SET\n");
//			prefixSum (xRight, leftGraph, el);
//			if (hierarchy)
//				oldwingDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
//			else
//				wingDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
//		}
//		else {
//			fprintf (fp, "RIGHT is the PRIMARY SET\n");
//			prefixSum (xRight, rightGraph, el);
//			if (hierarchy)
//				oldwingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
//			else
//				wingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
//		}
//	}
////	else if (nd == "BUILD_UNW_LEFT") {
////		unweighted_projection (leftGraph, rightGraph, vfile); // construct unweighted projected graph from LEFT
////		return 0;
////	}
////	else if (nd == "BUILD_UNW_RIGHT") {
////		unweighted_projection (rightGraph, leftGraph, vfile); // construct unweighted projected graph from RIGHT
////		return 0;
////	}
//	else if (nd == "BUILD_LEFT") {
//		weighted_projection (leftGraph, rightGraph, vfile); // construct weighted projected graph from LEFT
//		return 0;
//	}
//	else if (nd == "BUILD_RIGHT") {
//		weighted_projection (rightGraph, leftGraph, vfile); // construct weighted projected graph from RIGHT
//		return 0;
//	}
//	else if (nd == "RUN_CORE") {
//		base_kcore (gr, nEdge, K, hierarchy, &maxK, vfile, fp); // run k-core on the given graph
////		string ff = vfile + "_WER";
////		FILE* dd = fopen (ff.c_str(), "w");
////		for (int i = 0; i < gr.size(); i++) {
////			fprintf (dd, "%d: ", i);
////			for (int j = 0; j < gr[i].size(); j++)
////				fprintf (dd, "%d ", gr[i][j]);
////			fprintf (dd, "\n");
////		}
////		fclose (dd);
//	}
//	else
		if (nd == "RUN_TRUSS") {
		base_ktruss (gr, nEdge, K, hierarchy, &maxK, vfile, fp); // run k-truss on the given graph
	}
//	else if (nd == "RUN_WEIGHTED_CORE") {
//		weighted_base_kcore (wg, nEdge, K, hierarchy, &maxK, vfile.c_str(), fp); // run fractional core on the given weighted graph
//	}
//    else if (nd == "MEASURE_CORE_LEFT" || nd == "MEASURE_CORE_RIGHT") {
//    	string fl (argv[3]);
//		if (nd == "MEASURE_CORE_RIGHT")
//			ff_vertex_ind_density (fl, rightGraph); // database_conf
//		else if (nd == "MEASURE_CORE_LEFT")
//			ff_vertex_ind_density (fl, leftGraph); // condmat
//		printf ("densities updated\n");
//		return 0;
//	}
    else if (nd == "MEASURE_TRUSS_LEFT" || nd == "MEASURE_TRUSS_RIGHT") {
		string fl (argv[3]);
		if (nd == "MEASURE_TRUSS_RIGHT")
			ff_edge_ind_density (fl, rightGraph); // database_conf
		else if (nd == "MEASURE_TRUSS_LEFT")
			ff_edge_ind_density (fl, leftGraph); // condmat
		printf ("densities updated\n");
		return 0;
	}



	timestamp t2;

	printf ("%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d nButterflies: %d\tTotal time: ", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	cout << t2 - t1 << endl;
	fprintf (fp, "%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d\tnButterflies: %d\t", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	print_time (fp, "Total Time: ", t2 - t1);
	fclose (fp);


	return 0;
}

