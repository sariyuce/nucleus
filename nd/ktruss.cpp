#include "main.h"

inline int ifConnected (Graph& graph, Graph& orientedGraph, Graph& TC, vertex u, vertex v) {
	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);
	for (size_t k = 0; k < orientedGraph[a].size(); k++)
		if (orientedGraph[a][k] == b) {
			TC[a][k]++;
			return b;
		}
	return -1;
}

// per edge
long long countTriangles (Graph& graph, Graph& orientedGraph, Graph& TC) {
	FILE* fp = fopen ("tris", "w");
	for (size_t i = 0; i < orientedGraph.size(); i++) {
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orientedGraph[i].size(); k++) {
				vertex a = orientedGraph[i][j];
				vertex b = orientedGraph[i][k];
				vertex c = ifConnected (graph, orientedGraph, TC, a, b);
				if (c != -1) {
					fprintf (fp, "%d %d %d\n", a, b, c);
					TC[i][j]++;
					TC[i][k]++;
				}
			}
		}
	}
	fclose (fp);



	long long tc = 0;
	for (size_t i = 0; i < orientedGraph.size(); i++)
		for (size_t j = 0; j < orientedGraph[i].size(); j++)
			tc += TC[i][j];

	return tc;
}

void base_ktruss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();

	// Prepare a CSR-like structure to index edges and create directed graph from low degree vertices to higher degree vertices
	vector<vp> el;
	vector<vertex> xel;
	Graph orientedGraph;
	indexEdges (graph, el, xel, orientedGraph);
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orientedGraph[i].size(), 0);


//	// Triangle counting for each edge

 	long long tric = countTriangles (graph, orientedGraph, TC);

//	FILE* aa = fopen ("tcs", "r");
//	int a, b, t;
//	long long tric = 0;
//	char c1[2];
//	while (fscanf (aa, "%d %d %d", &a, &b, &t) != EOF) {
//		TC[a][b] = t;
//		tric += t;
////		if (tric > 100)
////			break;
//	}
//	fclose (fp);
//	cout << "tric: " << tric << endl;
//	exit(1);
	tric /= 3;
	fprintf (fp, "TC: %lld\n", tric);
	timestamp tcEnd;
	print_time (fp, "Triangle counting: ", tcEnd - peelingStart);


	// exit (1);
	// Peeling
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orientedGraph.size(); i++)
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
//			if (id % 100000 == 0) {
//				fprintf (fp, "id: %d  nEdge: %d\n", id, nEdge);
//				fflush (fp);
//			}
//			fprintf (fp, "tc  %d %d    %d\n", i, j, TC[i][j]);
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	fprintf (fp, "All inserted\n");
	vertex tc_e = 0;

//	exit(1);


	// required for hierarchy
//	vertex cid; // subcore id number
	vertex cid;
	vector<subcore> skeleton; // equal K valued cores
//	vector<vertex> component; // subcore ids for each vertex
	vector<vertex> component; // subcore ids for each vertex
	vector<vp> relations;
	vector<vertex> unassigned;
	vertex nSubcores;

	if (hierarchy) {
		cid = 0;
		nSubcores = 0;
		component.resize (nEdge, -1);
	}

	fprintf (fp, "Starts\n");
	fflush (fp);
	int counter = 0;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;


		if (counter % 10000 == 0) {
			fprintf (fp, "e: %d    val: %d    counter: %d  nEdge: %d\n", e, val, counter, nEdge);
			fflush (fp);
		}


		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			sc.children.clear();
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		vector<vertex> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // decrease the TC of the neighbor edges with greater TC
			edge f = getEdgeId (u, w, xel, el, graph);
			edge g = getEdgeId (v, w, xel, el, graph);
			if (K[f] == -1 && K[g] == -1) {
				if (nBucket.CurrentValue(f) > tc_e)
					nBucket.DecVal(f);
				if (nBucket.CurrentValue(g) > tc_e)
					nBucket.DecVal(g);
			}
			else if (hierarchy)
				createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
		counter++;
	}

	nBucket.Free();
	*maxtruss = tc_e; // tc_e is tc of the last popped edge

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "Nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}
