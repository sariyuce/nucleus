#include "main.h"


vertex count_acyclics (Graph& dgraph, vector<vertex>& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		outgoings (dgraph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			vector<vertex> ints;
			inter (2, 2, dgraph, u, v, ints); // blue orbit

			for (size_t k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[u][ints[k]]; // equal to dgraph[v][ints[k+1]]
				TC[u]++;
				TC[v]++;
				TC[w]++;
				// todo: node orbit counting is as follows
//				circle_orbit[u]++;
//				square_orbit[v]++;
//				triangle_orbit[w]++;
				count++;
			}
		}
	}

	// do not divide each TC (and count) by 3 since there is only one blue orbit
//	count /= 3;

	printf ("total acyclic count: %d\n", count);
	return count;
}

void acyclic_core (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	// Acyclic counting for each edge
	vertex nVtx = graph.size();
	vector<vertex> TC (nVtx, 0);
	lol tric =count_acyclics (graph, TC);

	fprintf (fp, "# acyclis: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Acyclic counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, nVtx);
	// each node and its acyclic-count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		printf ("acyclic count of %d is %d\n", i, TC[i]);
		if (TC[i] > 0)
			nBucket.Insert (i, TC[i]);
		else
			K[i]= 0;
	}

	vertex tc_u = 0;

	while (true) {
		vertex u, val;
		if (nBucket.PopMin(&u, &val) == -1) // if the bucket is empty
			break;

		tc_u = K[u] = val;

		// three node orbits
		// circle node orbit
		for (vertex j = 1; j < graph[u][0]; j++) { // outgoing edges
			vertex v = graph[u][j]; // square node orbit
			if (K[v] == -1) { // v might not be a part of acyclic, but must have -1 kappa if so anyways
				vector<vertex> ints;
				inter (2, 2, graph, u, v, ints); // todo: we don't need ints[k+1]
				for (auto k = 0; k < ints.size(); k+=2) {
					vertex w = graph[v][ints[k+1]]; // equal to graph[u][k]
					if (K[w] == -1) {
						if (nBucket.CurrentValue(v) > tc_u)
							nBucket.DecVal(v);
						if (nBucket.CurrentValue(w) > tc_u)
							nBucket.DecVal(w);
					}
				}
			}
		}

		//  square node orbit
		for (vertex j = 1; j < graph[u][0]; j++) { // outgoing edges
			vertex v = graph[u][j]; // triangle node orbit
			if (K[v] == -1) { // v might not be a part of cycle, but must have -1 kappa if so anyways
				vector<vertex> ints;
				inter (1, 1, graph, u, v, ints); // todo: we don't need ints[k+1]
				for (auto k = 0; k < ints.size(); k+=2) {
					vertex w = M2P (graph[v][ints[k+1]]); // equal to M2P (graph[u][k])
					if (K[w] == -1) {
						if (nBucket.CurrentValue(v) > tc_u)
							nBucket.DecVal(v);
						if (nBucket.CurrentValue(w) > tc_u)
							nBucket.DecVal(w);
					}
				}
			}
		}

		// triangle node orbit
		for (vertex j = graph[u][0]; j < graph[u].size(); j++) { // incoming edges
			vertex v = M2P (graph[u][j]); // circle node orbit
			if (K[v] == -1) { // v might not be a part of cycle, but must have -1 kappa if so anyways
				vector<vertex> ints;
				inter (1, 2, graph, u, v, ints); // todo: we don't need ints[k+1]
				for (auto k = 0; k < ints.size(); k+=2) {
					vertex w = graph[v][ints[k+1]]; // equal to M2P (graph[u][k])
					if (K[w] == -1) {
						if (nBucket.CurrentValue(v) > tc_u)
							nBucket.DecVal(v);
						if (nBucket.CurrentValue(w) > tc_u)
							nBucket.DecVal(w);
					}
				}
			}
		}
	}

	nBucket.Free();
	*maxK = tc_u;

	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}

	for (auto i = 0; i < K.size(); i++)
		printf ("core of %d is %d\n", i, K[i]); // the ones with -1 kappa do not participate in any acyclic
	return;
}


