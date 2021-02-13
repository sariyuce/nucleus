#include "main.h"

double simple_count_inps (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing) {
	vertex cut = 0, vol = 0, count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		asymmetric_undirecteds (u, dgraph[u], ret); // u-v edge is undirected s.t. u < v
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			vector<vertex> ints;
			inter (2, 2, dgraph, u, v, ints);
			for (size_t k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[u][ints[k]];

				int asdf = 0;
				if (crossing.find(u) != crossing.end())
					asdf++;
				if (crossing.find(v) != crossing.end())
					asdf++;
				if (crossing.find(w) != crossing.end())
					asdf++;

//				printf ("u: %d v: %d w: %d -- asdf: %d\n", u, v, w, asdf);

				if (asdf < 3) {
					if (numbers.find(u) != numbers.end())
						vol++;
					if (numbers.find(v) != numbers.end())
						vol++;
					if (numbers.find(w) != numbers.end())
						vol++;
					if (asdf > 0)
						cut++;
					if (asdf == 0)
						count++;
				}

			}
		}
	}

	printf ("cut: %d vol: %d -- cond: %lf\n", cut, vol, ((double) cut) / vol);

	int nv = numbers.size();
	printf ("count: %d numOfNodes: %d -- avg. motif degree: %lf\n", count, nv, ((double) count) / nv);

	return ((double) cut) / vol;
}

vertex count_inps (Graph& dgraph, Graph& TC) {
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
				TC[u][ret[r]]++; // u-v is undirected and u < v
				TC[u][ints[k]]++; // u-w
				TC[v][ints[k+1]]++; // v-w
				count++;
#ifdef ORBITS
				vector<vertex> a = {u, v};
				sort (a.begin(), a.end());
				auto t = make_tuple (a[0], a[1]);
				orb1[t]++;

				a = {v, w};
				sort (a.begin(), a.end());
				t = make_tuple (a[0], a[1]);
				orb2[t]++;

				a = {u, w};
				sort (a.begin(), a.end());
				t = make_tuple (a[0], a[1]);
				orb2[t]++;
#endif

			}
		}
	}
	printf ("total inp count: %d\n", count);

#ifdef ORBITS
	for (auto it = orb1.begin(); it != orb1.end(); it++)
		printf ("orb1 of %d %d is %d\n", get<0>(it->first), get<1>(it->first), it->second);

	for (auto it = orb2.begin(); it != orb2.end(); it++)
		printf ("orb2 of %d %d is %d\n", get<0>(it->first), get<1>(it->first), it->second);

	for (auto it = orb3.begin(); it != orb3.end(); it++)
		printf ("orb3 of %d %d is %d\n", get<0>(it->first), get<1>(it->first), it->second);

	exit(1);
#endif

	return count;
}
//
//void inp_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {
//	const auto t1 = chrono::steady_clock::now();
//	// Inp counting for each edge
//	vertex nVtx = graph.size();
//	Graph TC;
//	TC.resize(nVtx);
//	for (vertex i = 0; i < nVtx; i++)
//		TC[i].resize (graph[i].size(), 0);
//	lol tric = count_inps (graph, TC);
//	fprintf (fp, "# inps: %lld\n", tric);
//	const auto t2 = chrono::steady_clock::now();
//	print_time (fp, "Inp counting: ", t2 - t1);
//
//	// Peeling
//	const auto p1 = chrono::steady_clock::now();
//	K.resize (nEdge, -1);
//	Naive_Bucket nBucket;
//	nBucket.Initialize (nVtx, nEdge);
//	vertex id = 0;
//	vector<vp> el;
//	vector<vertex> xel;
//	xel.push_back(0);
//	// both non-reciprocal and undirected edges are inserted
//	// in el, if it's undirected edge, the second item is negative.
//	// Note that always u < v, so node 0 cannot exist as the second item.
//	for (vertex i = 0; i < graph.size(); i++) {
//		vector<vertex> ret;
//		outgoings_and_asymmetric_undirecteds (i, graph[i], ret); // undirecteds if only i < v
//		for (vertex r = 0; r < ret.size(); r++) {
//			vertex ind = abs (ret[r]);
//			vertex v = graph[i][ind];
//			vp c;
//			if (ret[r] < 0)
//				c = make_pair (i, -v); // undirected
//			else
//				c = make_pair (i, v); // directed
//			el.push_back(c);
//			printf ("inp count of %d (%d-%d) is %d\n", id, i, v, TC[i][ind]);
//			if (TC[i][ind] > 0)
//				nBucket.Insert (id++, TC[i][ind]);
//			else
//				K[id++] = 0;
//		}
//		xel.push_back(el.size());
//	}
//	vertex tc_e = 0;
//	while (true) {
//		edge e;
//		vertex val;
//		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
//			break;
//		tc_e = K[e] = val;
//		vertex u = el[e].first; // source
//		vertex v = el[e].second; // target
//		vector<vertex> ints;
//		vertex id1, id2;
//		if (v < 0) { // u-v is undirected
//			v *= -1;
//			inter (2, 2, graph, u, v, ints);
//			for (auto k = 0; k < ints.size(); k+=2) {
//				vertex w = graph[u][ints[k]];
//				id1 = getEdgeId (u, w, xel, el, graph); // directed
//				id2 = getEdgeId (v, w, xel, el, graph); // directed
//				checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
//			}
//		}
//		else { // directed
//			ints.clear();
//			inter (0, 1, graph, u, v, ints);
//			for (auto k = 0; k < ints.size(); k+=2) {
//				vertex w = M2P (graph[v][ints[k+1]]);
//				id1 = getEdgeId (w, v, xel, el, graph); // directed
//				id2 = getEdgeId (min (u, w), -1 * max (u, w), xel, el, graph); // undirected
//				checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
//			}
//		}
//	}
//	nBucket.Free();
//	*maxK = tc_e;
//	const auto p2 = chrono::steady_clock::now();
//
//	if (!hierarchy) {
//		print_time (fp, "Only peeling time: ", p2 - p1);
//		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
//	}
//	for (auto i = 0; i < el.size(); i++)
//		printf ("truss of %d (%d-%d) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any inp or u > v for the corresponding u-v edge
//	return;
//}

void inp_truss_SUBS (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp, string vfile) {
	const auto t1 = chrono::steady_clock::now();
	// Inp counting for each edge
	vertex nVtx = graph.size();
	Graph TC;
	TC.resize(nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (graph[i].size(), 0);
	lol tric = count_inps (graph, TC);
	fprintf (fp, "# inps: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Inp counting: ", t2 - t1);

	if (COUNT_ONLY)
		return;

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge);
	vertex id = 0;
	vector<vp> el;
	vector<vertex> xel;
	xel.push_back(0);

	vector<int> dist (nVtx+1, 0);

	// both non-reciprocal and undirected edges are inserted
	// in el, if it's undirected edge, the second item is negative.
	// Note that always u < v, so node 0 cannot exist as the second item.
	for (vertex i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings_and_asymmetric_undirecteds (i, graph[i], ret); // undirecteds if only i < v
		for (vertex r = 0; r < ret.size(); r++) {
			vertex ind = abs (ret[r]);
			vertex v = graph[i][ind];
			vp c;
			if (ret[r] < 0)
				c = make_pair (i, -v); // undirected
			else
				c = make_pair (i, v); // directed
			el.push_back(c);
//			printf ("inp count of %d (%d-%d) is %d\n", id, i, v, TC[i][ind]);
			if (TC[i][ind] > 0) {
				nBucket.Insert (id++, TC[i][ind]);
				dist[TC[i][ind]]++;
			}
			else
				K[id++] = 0;
		}
		xel.push_back(el.size());
	}
	vertex tc_e = 0;

	const auto c1 = chrono::steady_clock::now();
	if (!hierarchy && DEG_DIST) {
		for (vertex i = 0; i < dist.size(); i++)
			if (dist[i] > 0)
				printf ("deg %d %d\n", i, dist[i]);
	}
	const auto c2 = chrono::steady_clock::now();

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
		component.resize (nEdge, -1);
	}

	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;

		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;
		vertex u = el[e].first; // source
		vertex v = el[e].second; // target
		vector<vertex> ints;
		vertex id1, id2;
		if (v < 0) { // u-v is undirected
			v *= -1;
			inter (2, 2, graph, u, v, ints);
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = graph[u][ints[k]];
				id1 = getEdgeId (u, w, xel, el, graph); // directed
				id2 = getEdgeId (v, w, xel, el, graph); // directed
//				checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
				checkAndDecAndHier (id1, id2, &nBucket, tc_e, e, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}
		else { // directed
			ints.clear();
			inter (0, 1, graph, u, v, ints);
			for (auto k = 0; k < ints.size(); k+=2) {
				vertex w = M2P (graph[v][ints[k+1]]);
				id1 = getEdgeId (w, v, xel, el, graph); // directed
				id2 = getEdgeId (min (u, w), -1 * max (u, w), xel, el, graph); // undirected
//				checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
				checkAndDecAndHier (id1, id2, &nBucket, tc_e, e, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);

	}
	nBucket.Free();
	*maxK = tc_e;
	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		auto tm = p2 - p1 - (c2 - c1);
		print_time (fp, "Only peeling time: ", tm);
		print_time (fp, "Total time: ", tm + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxK, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total inp-truss nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total inp-truss nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}

	for (auto i = 0; i < el.size(); i++)
		if (K[i] == *maxK)
			printf ("truss of %d ( %d %d ) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any inp or u > v for the corresponding u-v edge
	return;
}
