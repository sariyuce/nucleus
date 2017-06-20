#include "main.h"

#define EL(A, B) (xRight[A]+B)

inline void limitedIntersection (Graph& rightGraph, vertex b, vertex c, Graph& butterflyCounts, vertex a, long long* bCount) {

	vertex i = 0, j = 0;
	while (rightGraph[b][i] < a) // b - a index
		i++;
	vertex ba = i;
	i++;
	while (rightGraph[c][j] < a) // c - a index
		j++;
	vertex ca = j;
	j++;

	while (i < rightGraph[b].size() && j < rightGraph[c].size()) {
		if (rightGraph[b][i] < rightGraph[c][j])
			i++;
		else if (rightGraph[c][j] < rightGraph[b][i])
			j++;
		else {
			// d is rightGraph[c][j] = rightGraph[b][i] b - d and c - d indices
			vertex bd = i;
			vertex cd = j;
			butterflyCounts[b][ba]++;
			butterflyCounts[b][bd]++;
			butterflyCounts[c][ca]++;
			butterflyCounts[c][cd]++;
			(*bCount)++;
			i++;
			j++;
		}
	}
}

void countButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, long long* bCount) {

	timestamp t1;
	int nedge = 0;
	for (vertex i = 0; i < leftGraph.size(); i++) {
		nedge += leftGraph[i].size();
		vertex a = i;
		for (vertex j = 0; j < leftGraph[i].size(); j++) {
			vertex b = leftGraph[i][j];
			for (vertex k = j + 1; k < leftGraph[i].size(); k++) {
				// intersection set with greater than a's
				int c = leftGraph[i][k];
				limitedIntersection (rightGraph, b, c, butterflyCounts, a, bCount);
			}
		}
	}
}

inline vertex getEdgeIndex (vertex a, vertex b, vector<vp>& el, vector<vertex>& xRight) {
	for (vertex i = xRight[a]; i < xRight[a+1]; i++)
		if (el[i].second == b)
			return i - xRight[a];
}

// no hierarchy construction
void wingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, long long* bCount) {

	timestamp countingStart;
	vector<vp> left_el;
	vector<vertex> xLeft;

	prefixSum (xLeft, leftGraph, left_el);

	Graph butterflyCounts (rightGraph.size());
	for (vertex i = 0; i < rightGraph.size(); i++)
		butterflyCounts[i].resize (rightGraph[i].size(), 0);
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge

//	FasterCountButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge

	vertex maxBc = 0;
#ifdef BUTTERFLIES
	string cc = vfile + "-butterflies";
	FILE* ffp = fopen (cc.c_str(), "w");
#endif
	for (vertex i = 0; i < rightGraph.size(); i++) {
		for (vertex j = 0; j < rightGraph[i].size(); j++) {
			if (butterflyCounts[i][j] > maxBc)
				maxBc = butterflyCounts[i][j];
#ifdef BUTTERFLIES
			fprintf (ffp, "BF %d %d    %d\n", i, rightGraph[i][j], butterflyCounts[i][j]);
#endif
		}
	}
#ifdef BUTTERFLIES
	fprintf (ffp, "bCount: %ld\n", *bCount);
	fclose (ffp);
#endif

	timestamp peelingStart;
	printf ("# bflys: %ld\n", *bCount);
	cout << "Counting butterflies per edge time (includes prefix sum): " << peelingStart - countingStart << endl;
	print_time (fp, "Counting butterflies per edge time (includes prefix sum): ", peelingStart - countingStart);
	printf ("maxBc: %d\n", maxBc);
	printf ("nEdge: %d\n", nEdge);

	// peeling
	K.resize (el.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (maxBc, nEdge);
	vertex bid = 0;

	for (vertex i = 0; i < butterflyCounts.size(); i++) {
		for (vertex j = 0; j < butterflyCounts[i].size(); j++) {
			vertex c = butterflyCounts[i][j];
			if (c > 0) {
				nBucket.Insert (bid++, c);
			}
			else
				K[bid++] = 0;
		}
	}

	vertex bf_e = 0;

	timestamp tint (0, 0);
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1)
			break;

		bf_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		vertex uvInd = getEdgeIndex (u, v, el, xRight);
		vertex vuInd = getEdgeIndex (v, u, left_el, xLeft);

		for (vertex vtInd = 0; vtInd < leftGraph[v].size(); vtInd++) {
			vertex t = leftGraph[v][vtInd];
			if (t == u || t == -1)
				continue;
			vertex tvInd = getEdgeIndex (t, v, el, xRight);
			vertex g = EL(t, tvInd);
			size_t i = 0, j = 0;
			while (i < rightGraph[u].size() && j < rightGraph[t].size()) {
				bool fl = false;
				if (rightGraph[t][j] == -1) {
					j++;
					fl = true;
				}
				if (rightGraph[u][i] == -1) {
					i++;
					fl = true;
				}
				if (fl)
					continue;

				if (rightGraph[u][i] < rightGraph[t][j])
					i++;
				else if (rightGraph[t][j] < rightGraph[u][i])
					j++;
				else if (rightGraph[u][i] != v) {
					vertex uwInd = i;
					vertex twInd = j;

					vertex f = EL(u, uwInd);
					vertex h = EL(t, twInd);

					if (nBucket.CurrentValue (f) > val)
						nBucket.DecVal (f);
					if (nBucket.CurrentValue (g) > val)
						nBucket.DecVal (g);
					if (nBucket.CurrentValue (h) > val)
						nBucket.DecVal (h);

					i++;
					j++;
				}
				else {
					i++;
					j++;
				}
			}
		}

		rightGraph[u][uvInd] = leftGraph[v][vuInd] = -1;
	}

	nBucket.Free();
	*maxbicore = bf_e;

	timestamp peelingEnd;
	cout << "Only peeling time: " << peelingEnd - peelingStart << endl;
	print_time (fp, "Only peeling time: ", peelingEnd - peelingStart);
	cout << "Total time: " << peelingEnd - countingStart << endl;
	print_time (fp, "Total time: ", peelingEnd - countingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif
}

void oldwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, long long* bCount) {

	timestamp peelingStart;

	Graph butterflyCounts (rightGraph.size());
	for (vertex i = 0; i < rightGraph.size(); i++)
		butterflyCounts[i].resize (rightGraph[i].size(), 0);
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right

	vertex maxBc = 0;
	for (auto g: butterflyCounts)
		for (auto c: g)
			if (c > maxBc)
				maxBc = c;

	// peeling
	K.resize (el.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (maxBc, nEdge);
	vertex bid = 0;

	for (vertex i = 0; i < butterflyCounts.size(); i++) {
		for (vertex j = 0; j < butterflyCounts[i].size(); j++) {
			vertex c = butterflyCounts[i][j];
			if (c > 0)
				nBucket.Insert (bid++, c);
			else
				K[bid++] = 0;
		}
	}

	vertex bf_e = 0;

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
		component.resize (el.size(), -1);
	}

	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1)
			break;

		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		bf_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		vertex uvInd = getEdgeIndex (u, v, el, xRight);

		// v is left vertex, u is right vertex, t is a right vertex and v's neighbor, w is a left vertex and neighbor to u and t
		// xyInd is the number i s.t. blaGraph[x][i] = y
		for (vertex vtInd = 0; vtInd < leftGraph[v].size(); vtInd++) {
			vertex t = leftGraph[v][vtInd];
			if (t == u)
				continue;
			vertex tvInd = getEdgeIndex (t, v, el, xRight);

			vector<vertex> ds;
			indicesIntersectionOld (rightGraph[u], rightGraph[t], ds, v); // each intersection is vertex w, on the left
			for (vertex j = 0; j < ds.size(); j += 2) {
				vertex uwInd = ds[j];
				vertex twInd = ds[j+1];
				// we have u, t on right; v, w on left; indices from u, t to v, w
				vertex f = EL(u, uwInd);
				vertex g = EL(t, tvInd);
				vertex h = EL(t, twInd);
				if (K[f] == -1 && K[g] == -1 && K[h] == -1) {
					if (nBucket.CurrentValue (f) > val)
						nBucket.DecVal (f);
					if (nBucket.CurrentValue (g) > val)
						nBucket.DecVal (g);
					if (nBucket.CurrentValue (h) > val)
						nBucket.DecVal (h);
				}
				else if (hierarchy)
					createSkeleton (e, {f, g, h}, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxbicore = bf_e;

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxbicore, relations, skeleton, &nSubcores, nEdge, rightGraph.size(), leftGraph.size());
		timestamp nucleusEnd;

		print_time (fp, "Wing decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d (in edges) \t\t |E|: %d\n", nSubcores, skeleton.size(), nEdge);

		helpers hp (&el);
		presentNuclei ("WING", skeleton, component, nEdge, hp, vfile, fp, leftGraph, rightGraph, &xRight);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}

//inline void addIn (int v, int k, HashMap<int>& hm, int w, int vi, Graph& wList, int* id, Graph& rightGraph) {
//	if (hm.hasDefaultValue (w)) {
//		hm[w] = *id;
//		(*id)++;
//		vector<int> ll;
//		wList.push_back(ll);
//	}
//	wList[hm[w]].push_back (v);
//	wList[hm[w]].push_back (vi);
//	wList[hm[w]].push_back (k);
//}
//
//inline void flyIn (int i, int x, int ivInd, int wvInd, Graph& butterflyCounts, int c) {
//	butterflyCounts[i][ivInd] += c;
//	butterflyCounts[x][wvInd] += c;
//}
//
//
//// best, use it
//vertex FasterCountButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, vertex* bCount) {
//
//	timestamp t1;
//
//	int nedge = 0, cost = 0;
//	HashMap<int> dup (0);
//	for (vertex i = 0; i < rightGraph.size(); i++) {
//		nedge += rightGraph[i].size();
//		dup.reset (0);
//		Graph wList;
//		HashMap<int> hm (-1);
//		int id = 0;
//		for (int j = 0; j < rightGraph[i].size(); j++) {
//			int v = rightGraph[i][j];
//			int vi = find_ind (leftGraph[v], i);
//			for (int k = 0; k < leftGraph[v].size(); k++) {
//				int w = leftGraph[v][k];
//				if (i < w) { // todo: other orders?
//					dup[w]++;
//					addIn (v, k, hm, w, vi, wList, &id, rightGraph);
//					cost += 2;
//				}
//			}
//		}
//
//		for (auto it = dup.begin(); it != dup.end(); it++) {
//			int x = it->first;
//			int count = it->second;
//			if (x != i && count > 1) {
//				int c = count - 1;
//				for (int j = 0; j < wList[hm[x]].size(); j+=3) {
//					int v = wList[hm[x]][j];
//					int vi = wList[hm[x]][j+1];
//					int vk = wList[hm[x]][j+2];
//					butterflyCounts[v][vi] += c;
//					butterflyCounts[v][vk] += c;
//					cost += 2;
// 				}
//				*bCount += nChoosek (count, 2);
//			}
//		}
//
////		if (i % 100 == 0) {
////			timestamp t2;
////			cout << "i: " << i << " / " << rightGraph.size() << "  time: " << t2 - t1 << endl;
////		}
//	}
//
////	FILE* ffp = fopen ("butterflies", "w");
////	for (vertex i = 0; i < Graph.size(); i++) {
////		for (vertex j = 0; j < rightGraph[i].size(); j++) {
////			vp p = make_pair (i, rightGraph[i][j]);
////			fprintf (ffp, "BF %d %d    %d\n", i, rightGraph[i][j], butterflyCounts[i][j]);
////		}
////	}
////	fprintf (ffp, "bCount: %d\n", *bCount);
////	fclose (ffp);
//	printf ("bFly: %d, Left: %d, Right: %d, nEdge: %d\n", *bCount, leftGraph.size(), rightGraph.size(), nedge);
//	printf ("count cost: %d\n", cost);
//	return 0;
//}
