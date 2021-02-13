#include "main.h"

vertex count_inps (Graph& dgraph, vector<vertex>& TC) {
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		asymmetric_undirecteds (u, dgraph[u], ret); // u-v edge is undirected s.t. u < v
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			vector<vertex> ints;
			inter (2, 2, dgraph, u, v, ints);
			for (size_t k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[u][ints[k]];
				TC[u]++;
				TC[v]++;
				TC[w]++;
				count++;
			}
		}
	}
	printf ("total inp count: %d\n", count);
	return count;
}

void inp_core (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {
	const auto t1 = chrono::steady_clock::now();
	// Inp counting for each node
	vertex nVtx = graph.size();
	vector<vertex> TC (nVtx, 0);
	lol tric = count_inps (graph, TC);
	fprintf (fp, "# outps: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Outp counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, nVtx);
	// each node and its inp count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		printf ("inp count of %d is %d\n", i, TC[i]);
		if (TC[i] > 0)
			nBucket.Insert (i, TC[i]);
		else
			K[i] = 0;
	}
	vertex tc_u = 0;
	while (true) {
		vertex u, val;
		if (nBucket.PopMin(&u, &val) == -1) // if the bucket is empty
			break;
		tc_u = K[u] = val;
		// u as the circle node
		vector<vertex> ret;
		undirecteds (graph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = graph[u][ret[r]]; // circle node
			if (K[v] == -1) { // v might not be a part of inp, but must have -1 kappa if so anyways
				vector<vertex> ints;
				inter (2, 2, graph, u, v, ints);
				for (auto k = 0; k < ints.size(); k+=2) {
					vertex w = graph[u][ints[k]]; // triangle node
					checkAndDec (K[w], -1, w, v, &nBucket, tc_u);
				}
			}
		}
		// u as the triangle node
		ret.clear();
		incomings (graph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = M2P (graph[u][ret[r]]); // circle node orbit
			if (K[v] == -1) { // v might not be a part of outp, but must have -1 kappa if so anyways
				vector<vertex> ints;
				inter (1, 0, graph, u, v, ints);
				for (auto k = 0; k < ints.size(); k+=2) {
					vertex w = M2P (graph[u][ints[k]]); // circle node
					if (v < w) // otherwise v-w pair is processed twice
						checkAndDec (K[w], -1, w, v, &nBucket, tc_u);
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

void inp_core_SUBS (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp, string vfile) {
	const auto t1 = chrono::steady_clock::now();
	// Inp counting for each node
	vertex nVtx = graph.size();
	vector<vertex> TC (nVtx, 0);
	lol tric = count_inps (graph, TC);
	fprintf (fp, "# outps: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Outp counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nVtx, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, nVtx);
	// each node and its inp count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		printf ("inp count of %d is %d\n", i, TC[i]);
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
		// u as the circle node
		vector<vertex> ret;
		undirecteds (graph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = graph[u][ret[r]]; // circle node
			vector<vertex> ints;
			inter (2, 2, graph, u, v, ints);
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = graph[u][ints[k]]; // triangle node
//				checkAndDec (K[w], -1, w, v, &nBucket, tc_u);
				checkAndDecAndHier (w, v, &nBucket, tc_u, u, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}
		// u as the triangle node
		ret.clear();
		incomings (graph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = M2P (graph[u][ret[r]]); // circle node orbit
			vector<vertex> ints;
			inter (1, 0, graph, u, v, ints);
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = M2P (graph[u][ints[k]]); // circle node
				if (v < w) // otherwise v-w pair is processed twice
					checkAndDecAndHier (w, v, &nBucket, tc_u, u, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
//					checkAndDec (K[w], -1, w, v, &nBucket, tc_u);
			}
		}

		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);

	}
	nBucket.Free();
	*maxK = tc_u;
	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxK, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total inp-core nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp;
		presentNuclei (12, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total inp-core nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}

	for (auto i = 0; i < K.size(); i++)
		printf ("core of %d is %d\n", i, K[i]);
	return;
}

