#include "main.h"

// per vertex
lol countTriangles (Graph& graph, Graph& orientedGraph, vector<vertex>& TC) {
	lol tc = 0;
	for (size_t i = 0; i < orientedGraph.size(); i++) {
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orientedGraph[i].size(); k++) {
				vertex a = orientedGraph[i][j];
				vertex b = orientedGraph[i][k];
				if (orientedConnected (graph, orientedGraph, a, b)) {
					TC[a]++;
					TC[b]++;
					TC[i]++;
					tc++;
				}
			}
		}
	}
	return tc;
}

void base_k13 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* max13, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices
	Graph orientedGraph;
	createOriented (orientedGraph, graph);

	// Triangle counting for each vertex
	vector<vertex> TC (nVtx, 0);
	lol tric = countTriangles (graph, orientedGraph, TC);

	fprintf (fp, "# triangles: %lld\n", tric);
	timestamp tcEnd;
	print_time (fp, "Triangle storing and counting: ", tcEnd - peelingStart);

	// Peeling
	K.resize (nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, nVtx); // maximum triangle count of a vertex is nEdge
	for (vertex i = 0; i < graph.size(); i++) {
		if (TC[i] > 0)
			nBucket.Insert (i, TC[i]);
		else
			K[i] = 0;
	}

	vertex tc_u = 0;

	// required for hierarchy
	vertex cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<vertex> component; // subcore ids for each vertex
	vector<vp> relations;
	vector<vertex> unassigned;
	vertex nSubcores;

	if (hierarchy) {
		cid = 0;
		nSubcores = 0;
		component.resize (graph.size(), -1);
	}

	while (true) {
		vertex u, val;
		if (nBucket.PopMin(&u, &val) == -1) // if the bucket is empty
			break;

		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_u = K[u] = val;

		for (vertex j = 0; j < graph[u].size(); j++) {
			vertex a = graph[u][j];
			for (vertex k = j + 1; k < graph[u].size(); k++) {
				vertex b = graph[u][k];
				if (orientedConnected (graph, orientedGraph, a, b)) {
					if (K[a] == -1 && K[b] == -1) {
						if (nBucket.CurrentValue(a) > tc_u)
							nBucket.DecVal(a);
						if (nBucket.CurrentValue(b) > tc_u)
							nBucket.DecVal(b);
					}
					else if (hierarchy)
						createSkeleton (u, {a, b}, &nSubcores, K, skeleton, component, unassigned, relations);
				}
			}
		}

		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*max13 = tc_u; // tc_u is tc of last popped vertex

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*max13, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "1-3 nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp;
		presentNuclei (13, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;
		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}
