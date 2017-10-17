#include "main.h"

int main(int argc, char *argv[]) {

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);  // RIGHT_TIP, LEFT_TIP, WING 910 for rightVertex based, 911 for leftVertex based, 92 for edge
	if (!(nd == "RIGHT_TIP" || nd == "LEFT_TIP" || nd == "WING" || nd.find ("BUILD") != string::npos || nd.find ("RUN") != string::npos || nd.find ("MEASURE") != string::npos) ) {
		printf ("Invalid algorithm, options are RIGHT_TIP, LEFT_TIP and WING for bipartite graphs\n");
		printf ("Or, RUN_WEIGHTED_CORE and RUN_(CORE|TRUSS) for weighted and unweighted projections\n");
		printf ("Creating those projections is possible via BUILD_(UNW|W)_(LEFT|RIGHT) in weighted and unweighted ways from either side\n");
		printf ("Lastly, MEASURE_(CORE|TRUSS)_(LEFT|RIGHT) computes the densities of the bipartite subgraphs which are induced on the nuclei of LEFT|RIGHT projections. Third argument should be resulting nuclei file for the projections.\n");
		exit(1);
	}

	// read the graph, gives sorted edges in left and rightGraph
	edge nEdge = 0;
	Graph leftGraph, rightGraph;

	string hrc;
#ifdef EXPS
	Graph gr;
	Wraph wg;
	if (nd == "RUN_CORE" || nd == "RUN_TRUSS" || nd == "RUN_WEIGHTED_CORE") {
		hrc = argv[3];
		readWeightedBinary (filename, wg, &nEdge);
		printf ("nVtx: %d   nEdge:%d\n", wg.size(), nEdge);
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
#endif
	if (nd == "RIGHT_TIP" || nd == "LEFT_TIP" || nd == "WING" || nd.find("BUILD") != string::npos ||  nd.find("MEASURE") != string::npos) {
		if (nd == "RIGHT_TIP" || nd == "LEFT_TIP" || nd == "WING")
			hrc = argv[3];
		readBipartite<vertex, vertex> (filename, &nEdge, leftGraph, rightGraph);
	}


	string vfile = gname + "_" + nd;
	string out_file;


	bool hierarchy = (hrc == "YES" ? true : false);
	if (hierarchy)
		out_file = vfile + "_Hierarchy";
	else
		out_file = vfile + "_K";

	FILE* fp;
	if (nd.find ("BUILD") == string::npos)
		fp = fopen (out_file.c_str(), "w");
	timestamp peelingTime (0, 0);

	timestamp t1;

	lol bCount = 0; // number of butterflies
	lol maxK; // maximum K value in the graph
	vector<vertex> K; // K values of vertices (or edges)

	if (nd == "RIGHT_TIP")
		tipDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Right vertices are primary
	else if (nd == "LEFT_TIP")
		tipDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount); // Left vertices are primary
	else if (nd == "WING") {
		if (leftGraph.size() < rightGraph.size()) { // select the smaller set for faster computation
			fprintf (fp, "LEFT is chosen as the PRIMARY set\n");
			if (hierarchy)
				wingDecompositionHrc (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount);
			else
				wingDecomposition (rightGraph, leftGraph, nEdge, K, hierarchy, &maxK, fp, &bCount);
		}
		else {
			fprintf (fp, "RIGHT is chosen as the PRIMARY set\n");
			if (hierarchy)
				wingDecompositionHrc (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, vfile, fp, &bCount);
			else
				wingDecomposition (leftGraph, rightGraph, nEdge, K, hierarchy, &maxK, fp, &bCount);
		}
	}

#ifdef EXPS
	else if (nd == "BUILD_LEFT") {
		weighted_projection (leftGraph, rightGraph, vfile); // construct weighted projected graph from LEFT
		return 0;
	}
	else if (nd == "BUILD_RIGHT") {
		weighted_projection (rightGraph, leftGraph, vfile); // construct weighted projected graph from RIGHT
		return 0;
	}
	else if (nd == "RUN_CORE")
		base_kcore (gr, hierarchy, nEdge, K, &maxK, vfile, fp); // run k-core on the given projection graph
	else if (nd == "RUN_TRUSS")
		base_ktruss (gr, hierarchy, nEdge, K, &maxK, vfile, fp); // run k-truss on the given projection graph
	else if (nd == "RUN_WEIGHTED_CORE")
		weighted_base_kcore (wg, hierarchy, nEdge, K, &maxK, vfile, fp); // run fractional core on the given weighted projection graph
    else if (nd == "MEASURE_CORE_LEFT" || nd == "MEASURE_CORE_RIGHT") {
    		string density_file (argv[3]);
		if (nd == "MEASURE_CORE_RIGHT")
			ff_vertex_ind_density (density_file, rightGraph); // database_conf
		else if (nd == "MEASURE_CORE_LEFT")
			ff_vertex_ind_density (density_file, leftGraph); // condmat
		printf ("densities updated\n");
	}
    else if (nd == "MEASURE_TRUSS_LEFT" || nd == "MEASURE_TRUSS_RIGHT") {
    		string density_file (argv[3]);
		if (nd == "MEASURE_TRUSS_RIGHT")
			ff_edge_ind_density (density_file, rightGraph); // database_conf
		else if (nd == "MEASURE_TRUSS_LEFT")
			ff_edge_ind_density (density_file, leftGraph); // condmat
		printf ("densities updated\n");
	}
#endif

	timestamp t2;

#ifdef K_VALUES
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "K[%d]: %lld\n", i, K[i]);
	fclose (kf);
#endif

	fclose (fp);
	printf ("%s\t|L|: %d\t|R|: %d\t|E|: %d\tmaxK: %d nButterflies: %d\n", gname.c_str(), leftGraph.size(), rightGraph.size(), nEdge, maxK, bCount);

	return 0;
}

