#include "main.h"

#define EL(A, B) (xRight[A]+B)

inline void limitedIntersection (Graph& rightGraph, vertex b, vertex c, Graph& butterflyCounts, vertex a, vertex* bCount) {

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

void countButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, vertex* bCount) {

	for (vertex i = 0; i < leftGraph.size(); i++) {
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

void wingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {

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
			indicesIntersection (rightGraph[u], rightGraph[t], ds, v); // each intersection is vertex w, on the left
			for (vertex j = 0; j < ds.size(); j += 2) {
				vertex uwInd = ds[j];
				vertex twInd = ds[j+1];
				// we have u, t on right; v, w on left; indices from u, t to v, w
				vertex f = EL(u, uwInd);
				vertex g = EL(t, tvInd);
				vertex h = EL(t, twInd);
				if (K[f] == -1 && K[g] == -1 && K[h]) {
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

