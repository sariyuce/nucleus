#include "main.h"

double simple_count_cycles (Graph& dgraph, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing) {
	vertex cut = 0, vol = 0, count = 0;
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

	vol /= 3;
	cut /= 3;
	count /= 3;

	printf ("cut: %d vol: %d -- cond: %lf\n", cut, vol, ((double) cut) / vol);

	int nv = numbers.size();
	printf ("count: %d numOfNodes: %d -- avg. motif degree: %lf\n", count, nv, ((double) count) / nv);

	return ((double) cut) / vol;
}

vertex count_cycles (Graph& dgraph, Graph& TC) {
#ifdef SIGNS
	unordered_map<triple, bool> mp;
	vertex s1, s2, s3, sum;
#endif
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
#ifndef SIGNS
				TC[u][ret[r]]++; // u-v
				TC[v][ints[k+1]]++; // v-w
				TC[w][ind (u, dgraph[w])]++; // w-u
				count++;
#else
				s1 = signs[u][ret[r]];
				s2 = signs[v][ints[k+1]];
				s3 = signs[w][ind (u, dgraph[w])];
				sum = s1 + s2 + s3;
				if ((MODE == 3 && sum == 3) ||
						(MODE == 1 && sum == 1) ||
						(MODE == -1 && sum == -1) ||
						(MODE == -3 && sum == -3)) {
					TC[u][ret[r]]++; // u-v
					TC[v][ints[k+1]]++; // v-w
					TC[w][ind (u, dgraph[w])]++; // w-u
					count++;
//					printf ("cycle: %d %d\t%d %d\t%d %d\t %d %d %d\n", u, v, v, w, w, u, signs[u][ret[r]], signs[v][ints[k+1]], signs[w][ind (u, dgraph[w])]);
				}
//				vector<vertex> a = {u, v, w};
//				sort (a.begin(), a.end());
//				auto t = make_tuple (a[0], a[1], a[2]);
//
//				if (mp.find (t) == mp.end()) {
//					mp[t] = true;
//				}
#endif

			}
		}
	}

	int cc = 0;
	for (vertex i = 0; i < TC.size(); i++)
		for (vertex j = 1; j < TC[i].size(); j++) {
//			cc += TC[i][j];
			TC[i][j] /= 3; // 3 green edge orbits in a cycle
		}
	count /= 3;
	printf ("total cycle count: %d\n", count);
//	printf ("total cc count: %d\n", cc);
	return count;
}
//
//void cycle_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {
//	const auto t1 = chrono::steady_clock::now();
//	// Cycle counting for each edge
//	vertex nVtx = graph.size();
//	Graph TC;
//	TC.resize(nVtx);
//	for (vertex i = 0; i < nVtx; i++)
//		TC[i].resize (graph[i].size(), 0);
//	lol tric = count_cycles (graph, TC);
//	fprintf (fp, "# cycles: %lld\n", tric);
//	const auto t2 = chrono::steady_clock::now();
//	print_time (fp, "Cycle counting: ", t2 - t1);
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
//	// each non-reciprocal edge and its cycle-count is inserted to bucket
//	for (vertex i = 0; i < graph.size(); i++) {
//		vector<vertex> ret;
//		outgoings (graph[i], ret);
//		for (vertex r = 0; r < ret.size(); r++) {
//			vp c (i, graph[i][ret[r]]);
//			el.push_back(c); // first is source, second is target
//			printf ("cycle count of %d is %d\n", id, TC[i][ret[r]]);
//			if (TC[i][ret[r]] > 0)
//				nBucket.Insert (id++, TC[i][ret[r]]);
//			else
//				K[id++]= 0;
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
//		// only one edge orbit
//		vector<vertex> ints;
//		inter (1, 2, graph, u, v, ints);
//		for (auto k = 0; k < ints.size(); k+=2) {
//			vertex w = graph[v][ints[k+1]]; // equal to M2P (graph[u][k])
//			vertex id1 = getEdgeId (w, u, xel, el, graph);
//			vertex id2 = getEdgeId (v, w, xel, el, graph);
//			checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
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
//		printf ("truss of %d (%d-%d) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any outp or u > v for the corresponding u-v edge
//	return;
//}

void cycle_truss_SUBS (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp, string vfile) {
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

	// each non-reciprocal edge and its cycle-count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings (graph[i], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vp c (i, graph[i][ret[r]]);
			el.push_back(c); // first is source, second is target
//			printf ("cycle count of %d is %d\n", id, TC[i][ret[r]]);
			if (TC[i][ret[r]] > 0) {
				nBucket.Insert (id++, TC[i][ret[r]]);
				dist[TC[i][ret[r]]]++;
			}
			else
				K[id++]= 0;
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
		// only one edge orbit
		vector<vertex> ints;
		inter (1, 2, graph, u, v, ints);
		for (auto k = 0; k < ints.size(); k+=2) {
			vertex w = graph[v][ints[k+1]]; // equal to M2P (graph[u][k])

#ifdef SIGNS
			vertex s1 = signs[u][ind (v, graph[u])];
			vertex s2 = signs[v][ind (w, graph[v])];
			vertex s3 = signs[w][ind (u, graph[w])];
			vertex sum = s1 + s2 + s3;
#endif

			vertex id1 = getEdgeId (w, u, xel, el, graph);
			vertex id2 = getEdgeId (v, w, xel, el, graph);
//			checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
#ifndef SIGNS
			checkAndDecAndHier (id2, id1, &nBucket, tc_e, e, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
#else
			if ((MODE == 3 && sum == 3) ||
					(MODE == 1 && sum == 1) ||
					(MODE == -1 && sum == -1) ||
					(MODE == -3 && sum == -3))
				checkAndDecAndHier (id2, id1, &nBucket, tc_e, e, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);

#endif

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
		print_time (fp, "Total cycle-truss nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total cycle-truss nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}

//	for (auto i = 0; i < el.size(); i++)
//		printf ("truss of %d (%d-%d) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any outp or u > v for the corresponding u-v edge
	return;
}

