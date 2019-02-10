#include "main.h"

vertex count_cycles (Graph& dgraph, vector<vertex>& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		for (vertex j = 1; j < dgraph[u].size(); j++) {
			vertex v = dgraph[u][j];
			if (v < 0)
				break;

			// another option is that you can check incoming edges of u
			if (exists (u, dgraph[v])) // u-v must be directed because no edge in cycle is reciprocal
				continue;

			vector<vertex> ints;
			inter (1, 2, dgraph, u, v, ints); // todo: two items written to ints although ints[k] is not used. because inter is generic, can be fixed later

			for (vertex k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[v][ints[k+1]]; // equal to M2P (dgraph[u][k])
				TC[u]++;
				TC[v]++;
				TC[w]++;
				count++;
			}
		}
	}

	for (vertex i = 0; i < TC.size(); i++)
			TC[i] /= 3; // 3 circle node orbits in a cycle
	count /= 3;

	printf ("total cycle count: %d\n", count);
	return count;
}


void cycle_core (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	// Cycle counting for each node
	vertex nVtx = graph.size();
	vector<vertex> TC (nVtx, 0);
	lol tric = count_cycles (graph, TC);

	fprintf (fp, "# cycles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Cycle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, nVtx);
	// each node and its cycle-count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		printf ("cycle count of %d is %d\n", i, TC[i]);
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

		// only one node orbit; has one incoming one outgoing
		for (vertex j = 1; j < graph[u][0]; j++) { // outgoing edges
			vertex v = graph[u][j];

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
		printf ("core of %d is %d\n", i, K[i]);
	return;
}

