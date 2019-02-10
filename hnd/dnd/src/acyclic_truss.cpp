#include "main.h"


vertex count_acyclics (Graph& dgraph, Graph& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		for (vertex j = 1; j < dgraph[u].size(); j++) {
			vertex v = dgraph[u][j];
			if (v < 0)
				break;

			if (exists (u, dgraph[v])) // u-v must be directed because no edges in acyclic is reciprocal
				continue;

			vector<vertex> ints;
			inter (2, 2, dgraph, u, v, ints); // blue orbit

			for (size_t k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[u][ints[k]]; // equal to dgraph[v][ints[k+1]]
				TC[u][j]++; // u-v
				TC[u][ints[k]]++; // u-w
				TC[v][ints[k+1]]++; // v-w
				count++;
			}
		}
	}

	// do not divide each TC (and count) by 3 since there is only one blue orbit
//	count /= 3;

	printf ("total acyclic count: %d\n", count);
	return count;
}

void acyclic_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	// Acyclic counting for each edge
	vertex nVtx = graph.size();
	Graph TC;
	TC.resize(nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (graph[i].size(), 0);
	lol tric =count_acyclics (graph, TC);

	fprintf (fp, "# acyclis: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Acyclic counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (graph.size(), nEdge);
	vertex id = 0;
	vector<vp> el;
	vector<vertex> xel;
	xel.push_back(0);

	for (vertex i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings (graph[i], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vp c (i, graph[i][ret[r]]);
			el.push_back(c); // first is source, second is target
			printf ("acyclic count of %d is %d\n", id, TC[i][ret[r]]);
			if (TC[i][ret[r]] > 0)
				nBucket.Insert (id++, TC[i][ret[r]]);
			else
				K[id++]= 0;
		}
		xel.push_back(el.size());
	}

//	for (auto i = 0; i < xel.size(); i++) {
//		printf ("xel[%d]: %d\n", i, xel[i]);
//		for (auto j = xel[i]; j < xel[i+1]; j++)
//			printf ("%d -> %d\n", el[j].first, el[j].second);
//	}

	vertex tc_e = 0;

	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;

		tc_e = K[e] = val;

		vertex u = el[e].first; // source
		vertex v = el[e].second; // target

//		printf ("u-v:   %d-%d\n", u, v);
		vector<vertex> ints;

		inter (2, 2, graph, u, v, ints); // blue orbit
		inter (1, 1, graph, u, v, ints); // green orbit
		inter (2, 1, graph, u, v, ints); // purple orbit

		vertex id1, id2;
		for (auto k = 0; k < ints.size(); k+=2) {
			id1 = id2 = -1;
			vertex w = graph[u][ints[k]];
			vertex x = graph[v][ints[k+1]];
//			printf ("w: %d   x: %d\t\t", w, x);
			if (w >= 0 && x >= 0) { // blue orbit
				id1 = getEdgeId (u, w, xel, el, graph);
				id2 = getEdgeId (v, w, xel, el, graph);
			}
			else if (w < 0 && x < 0) { // green orbit
				w = M2P (w);
				id1 = getEdgeId (w, u, xel, el, graph);
				id2 = getEdgeId (w, v, xel, el, graph);
			}
			else if (w >= 0 && x < 0) { // purple orbit
				id1 = getEdgeId (u, w, xel, el, graph);
				id2 = getEdgeId (w, v, xel, el, graph);
			}

//			printf ("id1: %d    id2: %d\n", id1, id2);
			if (id1 == -1 || id2 == -1) {
				printf ("I shouldn't see this\n");
				exit (1);
				continue;
			}

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
		printf ("truss of %d is %d\n", i, K[i]); // the ones with -1 kappa do not participate in any acyclic
	return;
}


