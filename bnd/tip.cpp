#include "main.h"

// per vertex (on the given side)
vertex countButterflies (Graph& rightGraph, Graph& leftGraph, vertex* butterflyCounts, vertex* bCount) {

	vertex maxBc = 0;
	for (vertex i = 0; i < rightGraph.size(); i++) {
		// combine neigs of neigs of i in d2Neighbors
		vector<vertex> d2Neighbors;
		for (auto v : rightGraph[i])
			d2Neighbors.insert (d2Neighbors.begin(), leftGraph[v].begin(), leftGraph[v].end());

		// count duplicates, each count indicates count choose 2 butterflies
		sort (d2Neighbors.begin(), d2Neighbors.end());
		vertex dupCount = 1, bc = 0;
		for (vertex j = 1; j < d2Neighbors.size(); j++) {
			if (d2Neighbors[j] == d2Neighbors[j-1])
				dupCount++;
			else {
				bc += nChoosek (dupCount, 2); // (dupCount choose 2)
				dupCount = 1;
			}
		}
		bc += nChoosek (dupCount, 2); // for last d2Neighbors[i]

		bc -= nChoosek (rightGraph[i].size(), 2); // remove the duplicate 'i's for vertex i
		butterflyCounts[i] = bc;
		if (bc > maxBc)
			maxBc = bc;
		*bCount += bc;
	}

	*bCount /= 2;
	return maxBc;
}

// rightGraph is primary, leftGraph is secondary
void tipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {

	timestamp peelingStart;

	vertex* butterflyCounts = (vertex *) malloc (sizeof(vertex) * rightGraph.size());
	vertex maxBc = countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right

	// peeling
	K.resize (rightGraph.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (maxBc+1, rightGraph.size());
	for (size_t i = 0; i < rightGraph.size(); i++)
		if (butterflyCounts[i] > 0)
			nBucket.Insert (i, butterflyCounts[i]);
		else
			K[i] = 0;

	vertex bf_u = 0;

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
		component.resize (rightGraph.size(), -1);
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

		bf_u = K[u] = val;

		for (size_t i = 0; i < rightGraph[u].size(); i++) {
			vertex v = rightGraph[u][i];
			for (size_t j = i + 1; j < rightGraph[u].size(); j++) {
				vertex w = rightGraph[u][j];
				vector<vertex> commons;
				intersection (leftGraph[v], leftGraph[w], commons);

				// each x in commons (except u) is the 4th vertex of the butterfly u, x (rights), v, w (lefts)
				for (vertex k = 0; k < commons.size(); k++) {
					vertex x = commons[k];
					if (x == u)
						continue;
					if (K[x] == -1) {
						if (nBucket.CurrentValue(x) > val)
							nBucket.DecVal(x);
					}
					else if (hierarchy)
						createSkeleton (u, {x}, &nSubcores, K, skeleton, component, unassigned, relations);
				}
			}
		}

		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxbicore = bf_u;

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxbicore, relations, skeleton, &nSubcores, nEdge, rightGraph.size(), leftGraph.size());
		timestamp nucleusEnd;

		print_time (fp, "Tip decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), rightGraph.size());

		helpers dummy;
		presentNuclei ("TIP", skeleton, component, nEdge, dummy, vfile, fp, leftGraph, rightGraph, NULL);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}
