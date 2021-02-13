#include "main.h"

int main (int argc, char *argv[]) {

	const auto t1 = chrono::steady_clock::now();
	if (argc < 3) {
		fprintf(stderr, "usage: %s "
				"\n <filename>"
				"\n <23 for triangles, 34 for four-cliques>"
				"\n <FFM FMM for 23; FFFM FFMM FMMM for 34"
				"\n <hierarchy?: YES or NO>\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);
	string rest = tmp.substr (0, tmp.find_last_of("/"));
	string nd (argv[2]);
	if (!(nd == "23" || nd == "34")) {
		printf ("Invalid algorithm, options are 23 and 34\n");
		exit(1);
	}

	// read the graph, give sorted edges in graph
	edge nEdge = 0;
	Graph graph;
	readGraph<vertex, edge> (filename, graph, &nEdge);

	string vfile = gname + "_" + nd;
	string out_file;

	// read the metadata (genders)
	vector<int> gender (graph.size());
	string hname = gname.substr (0, gname.find_last_of("."));
	string metafile = rest + "/fb100META/" + hname + ".mat.txt";
	printf ("metafile: %s\n", metafile.c_str());
	printf ("graph.size: %d\n", graph.size());
	FILE* p = fopen (metafile.c_str(), "r");

	int a, b, c, d, e, f, g;
	char dummy;
	int i = 1;

	while (1) {
		if (fscanf (p, "%d %c %d %c %d %c %d %c %d %c %d %c %d",
			&a, &dummy, &b, &dummy, &c, &dummy, &d, &dummy, &e, &dummy, &f, &dummy, &g) == EOF)
			break;
		gender[i++] = b;
	}

	string type (argv[3]);
	vfile += "_" + type;
	int labelType;
	if (type == "FFM")		 labelType = FFM;
	else if (type == "FMM")	 labelType = FMM;
	else if (type == "FFF")	 labelType = FFF;
	else if (type == "MMM")	 labelType = MMM;
	else if (type == "FFFF") labelType = FFFF;
	else if (type == "FFFM") labelType = FFFM;
	else if (type == "FFMM") labelType = FFMM;
	else if (type == "FMMM") labelType = FMMM;
	else if (type == "MMMM") labelType = MMMM;
	else if (type == "XXX")	 labelType = XXX;
	else {
		printf ("invalid type\n");
		exit(1);
	}

	string hrc (argv[4]);
	bool hierarchy = (hrc == "YES" ? true : false);
	if (hierarchy)
		out_file = vfile + "_Hierarchy";
	else
		out_file = vfile + "_K";

	FILE* fp = fopen (out_file.c_str(), "w");

	vertex maxK; // maximum K value in the graph
	vector<vertex> K;

	if (nd == "23")
		base_ktruss (graph, hierarchy, nEdge, K, &maxK, vfile, fp, gender, labelType);
	else if (nd == "34")
		base_k34 (graph, hierarchy, nEdge, K, &maxK, vfile, fp, gender, labelType);

#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%lld\n", K[i]);
	fclose (kf);
#endif

	const auto t2 = chrono::steady_clock::now();
	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s-nucleus: %d\n", gname.c_str(), graph.size(), nEdge, nd.c_str(), maxK);
	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);

	return 0;
}
