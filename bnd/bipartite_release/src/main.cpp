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
		if (nd == "RUN_CORE" || nd == "RUN_TRUSS") {
			gr.resize (wg.size());
			for (int i = 0; i < wg.size(); i++) {
				gr[i].resize (wg[i].size());
				for (int j = 0; j < wg[i].size(); j++) {
					gr[i][j] = wg[i][j].n;
				}
			}
		}
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

	// for projections
//	Graph gr;

	if (nd == "RIGHT_TIP") {
		 tipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
//		 oldtipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
	}
	else if (nd == "LEFT_TIP") {
		 tipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
//		 oldtipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
	}
	else if (nd == "WING") { // todo: should be same when you swap leftGraph and rightGraph
		if (leftGraph.size() < rightGraph.size()) {
			fprintf (fp, "LEFT is the PRIMARY SET\n");
			prefixSum (xRight, leftGraph, el);
			if (hierarchy)
				oldwingDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
			else
				wingDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
		}
		else {
			fprintf (fp, "RIGHT is the PRIMARY SET\n");
			prefixSum (xRight, rightGraph, el);
			if (hierarchy)
				oldwingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
			else
				wingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, el, xRight, &maxK, vfile, fp, &bCount);
		}
	}
//	else if (nd == "BUILD_UNW_LEFT") {
//		unweighted_projection (leftGraph, rightGraph, vfile); // construct unweighted projected graph from LEFT
//		return 0;
//	}
//	else if (nd == "BUILD_UNW_RIGHT") {
//		unweighted_projection (rightGraph, leftGraph, vfile); // construct unweighted projected graph from RIGHT
//		return 0;
//	}
	else if (nd == "BUILD_LEFT") {
		weighted_projection (leftGraph, rightGraph, vfile); // construct weighted projected graph from LEFT
		return 0;
	}
	else if (nd == "BUILD_RIGHT") {
		weighted_projection (rightGraph, leftGraph, vfile); // construct weighted projected graph from RIGHT
		return 0;
	}
	else if (nd == "RUN_CORE") {
		base_kcore (gr, nEdge, K, hierarchy, &maxK, vfile, fp); // run k-core on the given graph
//		string ff = vfile + "_WER";
//		FILE* dd = fopen (ff.c_str(), "w");
//		for (int i = 0; i < gr.size(); i++) {
//			fprintf (dd, "%d: ", i);
//			for (int j = 0; j < gr[i].size(); j++)
//				fprintf (dd, "%d ", gr[i][j]);
//			fprintf (dd, "\n");
//		}
//		fclose (dd);
	}
	else if (nd == "RUN_TRUSS") {
		base_ktruss (gr, nEdge, K, hierarchy, &maxK, vfile, fp); // run k-truss on the given graph
	}
	else if (nd == "RUN_WEIGHTED_CORE") {
		weighted_base_kcore (wg, nEdge, K, hierarchy, &maxK, vfile.c_str(), fp); // run fractional core on the given weighted graph
	}
    else if (nd == "MEASURE_CORE_LEFT" || nd == "MEASURE_CORE_RIGHT") {
    	string fl (argv[3]);
		if (nd == "MEASURE_CORE_RIGHT")
			ff_vertex_ind_density (fl, rightGraph); // database_conf
		else if (nd == "MEASURE_CORE_LEFT")
			ff_vertex_ind_density (fl, leftGraph); // condmat
		printf ("densities updated\n");
		return 0;
	}
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

	printf ("%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d nButterflies: %ld\tTotal time: ", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	cout << t2 - t1 << endl;
	fprintf (fp, "%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d\tnButterflies: %ld\t", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);
	print_time (fp, "Total Time: ", t2 - t1);
	fclose (fp);

	FILE* dd = fopen (vfile.c_str(), "w");
	for (int i = 0; i < K.size(); i++)
		fprintf (dd, "K[%d]: %d\n", i, K[i]);
	fclose (dd);
	return 0;
}

