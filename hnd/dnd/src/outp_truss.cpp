#include "main.h"

vertex count_outps (Graph& dgraph, Graph& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		undirecteds (dgraph[u], ret); // u-v edge is undirected
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			// u-v edge is maintained by u if u < v. Otherwise, v should maintain it.
			if (u > v)
				continue;
			vector<vertex> ints;
			inter (1, 1, dgraph, u, v, ints);

			for (size_t k = 0; k < ints.size(); k+=2) {
				vertex w = M2P (dgraph[u][ints[k]]); // equal to M2P (dgraph[v][ints[k+1]])
				TC[u][ret[r]]++; // u-v: undirected and u < v
				TC[w][ind (u, dgraph[w])]++; // w-u
				TC[w][ind (v, dgraph[w])]++; // w-v
				count++;
			}
		}
	}

	printf ("total out_p count: %d\n", count);
	return count;
}


void outp_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	// Acyclic counting for each edge
	vertex nVtx = graph.size();
	Graph TC;
	TC.resize(nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (graph[i].size(), 0);
	lol tric = count_outps (graph, TC);

	fprintf (fp, "# outps: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Outp counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (graph.size(), nEdge);
	vertex id = 0;
	vector<vp> el;
	vector<vertex> xel;
	xel.push_back(0);
	// each non-reciprocal edge and undirected edges are inserted with their outp counts
	// in el, it's undirected edge if the second item is negative. Note that always u < v, so node 0 cannot exist as the second item.
	for (size_t i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings_and_undirecteds (i, graph[i], ret); // undirecteds if only i < v
		for (vertex r = 0; r < ret.size(); r++) {
			vertex ind = abs (ret[r]);
			vertex v = graph[i][ind];
			vp c;
			if (ret[r] < 0)
				c = make_pair (i, -v); // undirected
			else
				c = make_pair (i, v); // directed
			el.push_back(c);
			printf ("outp count of %d (%d-%d) is %d\n", id, i, v, TC[i][ind]);
			if (TC[i][ind] > 0)
				nBucket.Insert (id++, TC[i][ind]);
			else
				K[id++] = 0;
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

		vector<vertex> ints;
		vertex id1, id2;
		if (v < 0) { // undirected. green neighborhood
			v *= -1;
			inter (1, 1, graph, u, v, ints); // green neighbors.
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = M2P (graph[u][ints[k]]);
				id1 = getEdgeId (w, u, xel, el, graph); // directed
				id2 = getEdgeId (w, v, xel, el, graph); // directed
			}
		}
		else { // directed. blue neighborhood
			vector<vertex> ints;
			inter (2, 0, graph, u, v, ints); // blue neighbors.
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = graph[u][ints[k]];
				id1 = getEdgeId (u, w, xel, el, graph); // directed
				id2 = getEdgeId (min (w, v), -1 * max (w, v), xel, el, graph); // undirected
			}
		}

		if (K[id1] == -1 && K[id2] == -1) {
			if (nBucket.CurrentValue(id1) > tc_e)
				nBucket.DecVal(id1);
			if (nBucket.CurrentValue(id2) > tc_e)
				nBucket.DecVal(id2);
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
		printf ("truss of %d (%d-%d) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any outp or u > v for the corresponding u-v edge
	return;
}
