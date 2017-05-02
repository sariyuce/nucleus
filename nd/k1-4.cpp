#include "main.h"

// per vertex
void count4cliques (Graph& graph, Graph& orientedGraph, vector<vertex>& FC) {
	for (vertex i = 0; i < orientedGraph.size(); i++) {
		for (vertex j = 0; j < orientedGraph[i].size(); j++) {
			for (vertex k = j + 1; k < orientedGraph[i].size(); k++) {
				for (vertex l = k + 1; l < orientedGraph[i].size(); l++) {
					vertex a = orientedGraph[i][j];
					vertex b = orientedGraph[i][k];
					vertex c = orientedGraph[i][l];
					if (orientedConnected(graph, orientedGraph, a, b) &&
							orientedConnected(graph, orientedGraph, b, c) &&
							orientedConnected(graph, orientedGraph, a, c)) {
						FC[a]++;
						FC[b]++;
						FC[c]++;
						FC[i]++;
					}
				}
			}
		}
	}
}

void base_k14 (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* max14, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices
	Graph orientedGraph;
	createOriented (orientedGraph, graph);

	// 4-clique counting for each vertex
	vector<vertex> FC (nVtx, 0);
	count4cliques (graph, orientedGraph, FC);

	vertex maxFC = 0;
	for (vertex i = 0; i < graph.size(); i++)
		if (FC[i] > maxFC)
			maxFC = FC[i];

	// Peeling
	K.resize(nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize(maxFC, nVtx);
	for (vertex i = 0; i < graph.size(); i++) {
		if (FC[i] > 0)
			nBucket.Insert (i, FC[i]);
		else
			K[i] = 0;
	}

	vertex fc_u = 0;

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

		fc_u = K[u] = val;

		for (vertex j = 0; j < graph[u].size(); j++) {
			vertex a = graph[u][j];
			for (vertex k = j + 1; k < graph[u].size(); k++) {
				vertex b = graph[u][k];
				if (orientedConnected(graph, orientedGraph, a, b)) {
					for (vertex l = k + 1; l < graph[u].size(); l++) {
						vertex c = graph[u][l];
						if (orientedConnected (graph, orientedGraph, a, c) && orientedConnected (graph, orientedGraph, b, c)) {
							if (K[a] == -1 && K[b] == -1 && K[c] == -1) {
								if (nBucket.CurrentValue(a) > fc_u)
									nBucket.DecVal(a);
								if (nBucket.CurrentValue(b) > fc_u)
									nBucket.DecVal(b);
								if (nBucket.CurrentValue(c) > fc_u)
									nBucket.DecVal(c);
							}
							else if (hierarchy)
								createSkeleton (u, {a, b, c}, &nSubcores, K, skeleton, component, unassigned, relations);
						}
					}
				}
			}
		}

		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*max14 = fc_u; // fc_u is fc of last popped vertex

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*max14, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "Nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp;
		presentNuclei (14, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}

}
