#include "main.h"

vertex count_acyclics (Graph& dgraph, Graph& TC) {
#ifdef SIGNS
	unordered_map<triple, bool> mp;
	vertex s1, s2, s3, sum, _11, _12, _13, _m11, _m12, _m13;
#endif
	vertex count = 0;
	for (vertex u = 0; u < dgraph.size(); u++) {
		vector<vertex> ret;
		outgoings (dgraph[u], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vertex v = dgraph[u][ret[r]];
			vector<vertex> ints;
			inter (2, 2, dgraph, u, v, ints); // blue orbit
			for (vertex k = 0; k < ints.size(); k+=2) {
				vertex w = dgraph[u][ints[k]]; // equal to dgraph[v][ints[k+1]]
#ifndef SIGNS
				TC[u][ret[r]]++; // u-v
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
				orb3[t]++;
#endif
#else
				_m11 = _11 = s1 = signs[u][ret[r]]; // u-v
				_m12 = _12 = s2 = signs[u][ints[k]]; // u-w
				_m13 = _13 = s3 = signs[v][ints[k+1]]; // v-w
				sum = s1 + s2 + s3;
				if ((MODE == 3 && sum == 3) ||
						(MODE == 11 && _11 == -1 && sum == 1) ||
						(MODE == 12 && _12 == -1 && sum == 1) ||
						(MODE == 13 && _13 == -1 && sum == 1) ||
						(MODE == -11 && _m11 == 1 && sum == -1) ||
						(MODE == -12 && _m12 == 1 && sum == -1) ||
						(MODE == -13 && _m13 == 1 && sum == -1) ||
						(MODE == -3 && sum == -3)) {
					TC[u][ret[r]]++; // u-v
					TC[u][ints[k]]++; // u-w
					TC[v][ints[k+1]]++; // v-w
					count++;
//					vector<vertex> a = {u, v, w};
//					sort (a.begin(), a.end());
//					auto t = make_tuple (a[0], a[1], a[2]);
//
//					if (mp.find (t) == mp.end()) {
//						mp[t] = true;
//						printf ("acyclic: %d %d\t%d %d\t%d %d\t %d %d %d\n", u, v, v, w, w, u, signs[u][ret[r]], signs[u][ints[k]], signs[v][ints[k+1]]);
//					}
				}

#endif

			}
		}
	}
	printf ("total acyclic count: %d\n", count);

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
//void acyclic_truss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp) {
//	const auto t1 = chrono::steady_clock::now();
//	// Acyclic counting for each edge
//	vertex nVtx = graph.size();
//	Graph TC;
//	TC.resize(nVtx);
//	for (vertex i = 0; i < nVtx; i++)
//		TC[i].resize (graph[i].size(), 0);
//	lol tric =count_acyclics (graph, TC);
//	fprintf (fp, "# acyclis: %lld\n", tric);
//	const auto t2 = chrono::steady_clock::now();
//	print_time (fp, "Acyclic counting: ", t2 - t1);
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
//	// each non-reciprocal edge and its acyclic-count is inserted to bucket
//	for (vertex i = 0; i < graph.size(); i++) {
//		vector<vertex> ret;
//		outgoings (graph[i], ret);
//		for (vertex r = 0; r < ret.size(); r++) {
//			vp c (i, graph[i][ret[r]]);
//			el.push_back(c); // first is source, second is target
//			printf ("acyclic count of %d is %d\n", id, TC[i][ret[r]]);
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
//		vector<vertex> ints;
//		inter (2, 2, graph, u, v, ints); // blue orbit
//		inter (1, 1, graph, u, v, ints); // green orbit
//		inter (2, 1, graph, u, v, ints); // purple orbit
//		vertex id1, id2;
//		for (auto k = 0; k < ints.size(); k+=2) {
//			vertex w = graph[u][ints[k]];
//			vertex x = graph[v][ints[k+1]];
//			if (w >= 0 && x >= 0) { // blue orbit
//				id1 = getEdgeId (u, w, xel, el, graph);
//				id2 = getEdgeId (v, w, xel, el, graph);
//			}
//			else if (w < 0 && x < 0) { // green orbit
//				w = M2P (w);
//				id1 = getEdgeId (w, u, xel, el, graph);
//				id2 = getEdgeId (w, v, xel, el, graph);
//			}
//			else if (w >= 0 && x < 0) { // purple orbit
//				id1 = getEdgeId (u, w, xel, el, graph);
//				id2 = getEdgeId (w, v, xel, el, graph);
//			}
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

void acyclic_truss_SUBS (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxK, FILE* fp, string vfile) {
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

	// each non-reciprocal edge and its acyclic-count is inserted to bucket
	for (vertex i = 0; i < graph.size(); i++) {
		vector<vertex> ret;
		outgoings (graph[i], ret);
		for (vertex r = 0; r < ret.size(); r++) {
			vp c (i, graph[i][ret[r]]);
			el.push_back(c); // first is source, second is target
//			printf ("acyclic count of %d is %d\n", id, TC[i][ret[r]]);
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
		vector<vertex> ints;
		inter (2, 2, graph, u, v, ints); // blue orbit
		inter (1, 1, graph, u, v, ints); // green orbit
		inter (2, 1, graph, u, v, ints); // purple orbit
		vertex id1, id2;
#ifdef SIGNS
		vertex s1, s2, s3, sum, _11, _12, _13, _m11, _m12, _m13;
#endif
		for (auto k = 0; k < ints.size(); k+=2) {
			vertex w = graph[u][ints[k]];
			vertex x = graph[v][ints[k+1]];
			if (w >= 0 && x >= 0) { // blue orbit
				id1 = getEdgeId (u, w, xel, el, graph);
				id2 = getEdgeId (v, w, xel, el, graph);
#ifdef SIGNS
				_m11 = _11 = s1 = signs[u][ind (v, graph[u])]; // MODE 11: s1 is -1; MODE -11: s1 is 1
				_m12 = _12 = s2 = signs[u][ind (w, graph[u])]; // MODE 12: s2 is -1; MODE -12: s2 is 1
				_m13 = _13 = s3 = signs[v][ind (w, graph[v])]; // MODE 13: s3 is -1; MODE -13: s3 is 1
				sum = s1 + s2 + s3;
#endif
			}
			else if (w < 0 && x < 0) { // green orbit
				w = M2P (w);
				id1 = getEdgeId (w, u, xel, el, graph);
				id2 = getEdgeId (w, v, xel, el, graph);
#ifdef SIGNS
				_m13 = _13 = s1 = signs[u][ind (v, graph[u])]; // MODE 13: s1 is -1; MODE -13: s1 is 1
				_m11 = _11 = s2 = signs[w][ind (u, graph[w])]; // MODE 11: s2 is -1; MODE -11: s2 is 1
				_m12 = _12 = s3 = signs[w][ind (v, graph[w])]; // MODE 12: s3 is -1; MODE -12: s3 is 1
				sum = s1 + s2 + s3;
#endif
			}
			else if (w >= 0 && x < 0) { // purple orbit
				id1 = getEdgeId (u, w, xel, el, graph);
				id2 = getEdgeId (w, v, xel, el, graph);
#ifdef SIGNS
				_m12 = _12 = s1 = signs[u][ind (v, graph[u])];// MODE 12: s1 is -1; MODE -12: s1 is 1;
				_m11 = _11 = s2 = signs[u][ind (w, graph[u])];// MODE 11: s2 is -1; MODE -11: s2 is 1
				_m13 = _13 = s3 = signs[w][ind (v, graph[w])];// MODE 13: s3 is -1; MODE -13: s3 is 1
				sum = s1 + s2 + s3;
#endif

			}
//			checkAndDec (K[id1], K[id2], id1, id2, &nBucket, tc_e);
#ifndef SIGNS
			checkAndDecAndHier (id2, id1, &nBucket, tc_e, e, hierarchy, &nSubcores, K, skeleton, component, unassigned, relations);
#else
			if ((MODE == 3 && sum == 3) ||
					(MODE == 11 && _11 == -1 && sum == 1) ||
					(MODE == 12 && _12 == -1 && sum == 1) ||
					(MODE == 13 && _13 == -1 && sum == 1) ||
					(MODE == -11 && _m11 == 1 && sum == -1) ||
					(MODE == -12 && _m12 == 1 && sum == -1) ||
					(MODE == -13 && _m13 == 1 && sum == -1) ||
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
		print_time (fp, "Total acyclic-truss nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total acyclic-truss nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}

//	for (auto i = 0; i < el.size(); i++)
//		printf ("truss of %d (%d-%d) is %d\n", i, el[i].first, abs(el[i].second), K[i]); // the ones with -1 kappa either do not participate in any outp or u > v for the corresponding u-v edge
	return;
}

