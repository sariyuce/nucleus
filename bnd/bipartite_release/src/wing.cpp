#include "main.h"

#define EL(A, B) (xRight[A]+B)

inline void limitedIntersection (Graph& rightGraph, vertex b, vertex c, Graph& butterflyCounts, vertex a, lol* bCount) {

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

void countButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, lol* bCount) {
	for (vertex i = 0; i < leftGraph.size(); i++) {
		vertex a = i;
		for (vertex j = 0; j < leftGraph[i].size(); j++) {
			vertex b = leftGraph[i][j];
			for (vertex k = j + 1; k < leftGraph[i].size(); k++) {
				// intersection set with greater than a's
				vertex c = leftGraph[i][k];
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
void wingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, lol* maxbicore, string vfile, FILE* fp, lol* bCount) {

	timestamp pr1;
	vector<vp> left_el;
	vector<vertex> xLeft;
	prefixSum (xLeft, leftGraph, left_el);
	vector<vp> right_el;
	vector<vertex> xRight;
	prefixSum (xRight, rightGraph, right_el);
	timestamp pr2;
	cout << "Prefix computation time (both sides): " << pr2 - pr1 << endl;

	timestamp c1;
	Graph butterflyCounts (rightGraph.size());
	for (vertex i = 0; i < rightGraph.size(); i++)
		butterflyCounts[i].resize (rightGraph[i].size(), 0);
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge
	vertex maxBc = 0;
	for (auto g : butterflyCounts)
		for (auto c : g)
			if (c > maxBc)
				maxBc = c;
	timestamp c2;
	printf ("# bflys: %lld\t\t maxBc: %lld\n", *bCount, maxBc);
	cout << "Counting butterflies per edge time: " << c2 - c1 << endl;


	// peeling
	timestamp p1;
	K.resize (right_el.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (maxBc, nEdge);
	vertex bid = 0;

	for (vertex i = 0; i < butterflyCounts.size(); i++)
		for (vertex j = 0; j < butterflyCounts[i].size(); j++) {
			vertex c = butterflyCounts[i][j];
			if (c > 0)
				nBucket.Insert (bid++, c);
			else
				K[bid++] = 0;
		}

	vertex bf_e = 0;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1)
			break;

		bf_e = K[e] = val;

		vertex u = right_el[e].first;
		vertex v = right_el[e].second;
		vertex uvInd = getEdgeIndex (u, v, right_el, xRight);
		vertex vuInd = getEdgeIndex (v, u, left_el, xLeft);

		// v is left vertex, u is right vertex, t is a right vertex and v's neighbor, w is a left vertex and neighbor to u and t
		// xyInd is the number i s.t. blaGraph[x][i] = y
		for (vertex vtInd = 0; vtInd < leftGraph[v].size(); vtInd++) {
			vertex t = leftGraph[v][vtInd];
			if (t == u || t == -1)
				continue;
			vertex tvInd = getEdgeIndex (t, v, right_el, xRight);
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

	timestamp p2;
	cout << "Only peeling time: " << p2 - p1 << endl;
	cout << "Total time: " << (p2 - p1) + (c2 - c1) + (pr2 - pr1)<< endl;
}

inline void indicesIntersectionHrc (vector<vertex>& a, vector<vertex>& b, vector<vertex>& res, vertex g) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			if (a[i] != g) {
				res.push_back(i);
				res.push_back(j);
			}
			i++;
			j++;
		}
	}
}

void wingDecompositionHrc (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, lol* maxbicore, string vfile, FILE* fp, long long* bCount) {

	timestamp pr1;
	vector<vp> right_el;
	vector<vertex> xRight;
	prefixSum (xRight, rightGraph, right_el);
	timestamp pr2;
	cout << "Prefix computation time (single side): " << pr2 - pr1 << endl;

	timestamp c1;
	Graph butterflyCounts (rightGraph.size());
	for (vertex i = 0; i < rightGraph.size(); i++)
		butterflyCounts[i].resize (rightGraph[i].size(), 0);
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right
	vertex maxBc = 0;
	for (auto g : butterflyCounts)
		for (auto c : g)
			if (c > maxBc)
				maxBc = c;
	timestamp c2;
	printf ("# bflys: %lld\t\t maxBc: %lld\n", *bCount, maxBc);
	cout << "Counting butterflies per edge time: " << c2 - c1 << endl;


	// peeling and hierarchy construction
	timestamp p1;
	K.resize (right_el.size(), -1);
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
		component.resize (right_el.size(), -1);
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

		vertex u = right_el[e].first;
		vertex v = right_el[e].second;
		vertex uvInd = getEdgeIndex (u, v, right_el, xRight);

		// v is left vertex, u is right vertex, t is a right vertex and v's neighbor, w is a left vertex and neighbor to u and t
		// xyInd is the number i s.t. blaGraph[x][i] = y
		for (vertex vtInd = 0; vtInd < leftGraph[v].size(); vtInd++) {
			vertex t = leftGraph[v][vtInd];
			if (t == u)
				continue;
			vertex tvInd = getEdgeIndex (t, v, right_el, xRight);
			vector<vertex> ds;
			indicesIntersectionHrc (rightGraph[u], rightGraph[t], ds, v); // each intersection is vertex w, on the left
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

	timestamp p2;
	cout << "Peeling + on-the-fly hierarchy construction time: " << p2 - p1 << endl;
	timestamp b1;
	buildHierarchy (*maxbicore, relations, skeleton, &nSubcores, nEdge, rightGraph.size(), leftGraph.size());
	timestamp b2;

	cout << "Building hierarchy time: " << b2 - b1 << endl;
	cout << "Total time (excluding density computation): " << (pr2 - pr1) + (c2 - c1) + (p2 - p1) + (b2 - b1) << endl;

	timestamp d1;
	helpers hp (&right_el);
	presentNuclei ("WING", skeleton, component, nEdge, hp, vfile, fp, leftGraph, rightGraph, &xRight);
	timestamp d2;

	cout << "Total time: " << (pr2 - pr1) + (c2 - c1) + (p2 - p1) + (b2 - b1) + (d2 - d1) << endl;
}
