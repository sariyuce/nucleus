#include "main.h"

vertex count_cycles (Graph& dgraph, Graph& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		outgoings (dgraph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			vector<vertex> ints;
			inter (1, 2, dgraph, u, v, ints); // green orbit
			// todo: two items written to ints although ints[k] is not used.
			// because inter is generic, can be fixed later
			for (vertex k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[v][ints[k+1]]; // equal to M2P (dgraph[u][ints[k]])
				TC[u][ret[r]]++; // u-v
				TC[v][ints[k+1]]++; // v-w
				TC[w][ind (u, dgraph[w])]++; // w-u
				count++;
			}
		}
	}

	for (vertex i = 0; i < TC.size(); i++)
		for (vertex j = 1; j < TC[i].size(); j++)
			TC[i][j] /= 3; // 3 green edge orbits in a cycle
	count /= 3;
	printf ("total cycle count: %d\n", count);
	return count;
}

void cycle_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {
	const auto t1 = chrono::steady_clock::now();
	// Cycle counting for each edge
	vertex nVtx = graph.size();
	Graph TC;
	TC.resize(nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (graph[i].size(), 0);
	lol tric = count_cycles (graph, TC);
	fprintf (fp, "# cycles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Cycle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge);
	vertex id = 0;
	vector<vp> el;
	vector<vertex> xel;
	xel.push_back(0);
	// each non-reciprocal edge and its cycle-count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings (graph[i], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vp c (i, graph[i][ret[r]]);
			el.push_back(c); // first is source, second is target
			printf ("cycle count of %d is %d\n", id, TC[i][ret[r]]);
			if (TC[i][ret[r]] > 0)
				nBucket.Insert (id++, TC[i][ret[r]]);
			else
				K[id++]= 0;
		}
		xel.push_back(el.size());
	}
	vertex tc_e = 0;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
		tc_e = K[e] = val;
		vertex u = el[e].first; // source
		vertex v = el[e].second; // target
		// only one edge orbit
		vector<vertex> ints;
		inter (1, 2, graph, u, v, ints);
		for (auto k = 0; k < ints.size(); k+=2) {
			vertex w = graph[v][ints[k+1]]; // equal to M2P (graph[u][k])
			vertex id1 = getEdgeId (w, u, xel, el, graph);
			vertex id2 = getEdgeId (v, w, xel, el, graph);
			if (K[id1] == -1 && K[id2] == -1) {
				if (nBucket.CurrentValue(id1) > tc_e)
					nBucket.DecVal(id1);
				if (nBucket.CurrentValue(id2) > tc_e)
					nBucket.DecVal(id2);
			}
		}
	}
	nBucket.Free();
	*maxK = tc_e;
	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	for (auto i = 0; i < K.size(); i++)
		printf ("truss of %d is %d\n", i, K[i]); // the ones with -1 kappa do not participate in any cycle
	return;
}

