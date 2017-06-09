#include "main.h"

vertex countButterflies (Graph& rightGraph, Graph& leftGraph, vertex* butterflyCounts, long* bCount) {

	timestamp t1;
	vertex maxBc = 0;
	HashMap<int> dup (0);
	for (vertex i = 0; i < rightGraph.size(); i++) {
		dup.reset (0);
		for (auto v : rightGraph[i])
			for (auto w : leftGraph[v]) {
				if (i < w) { // other order definitions
					dup[w]++;
				}
			}

		for (auto it = dup.begin(); it != dup.end(); it++) {
			int x = it->first;
			int count = it->second;
			if (x != i && count > 1) {
				int c = nChoosek (count, 2);
				butterflyCounts[i] += c;
				butterflyCounts[x] += c;
			}
		}
		if (butterflyCounts[i] > maxBc)
			maxBc = butterflyCounts[i];
		*bCount += butterflyCounts[i];
		timestamp t3;
	}

	*bCount /= 2;

	printf ("bFly: %ld, Left: %d, Right: %d\n", *bCount, leftGraph.size(), rightGraph.size());
	return maxBc;
}

void tipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, long* bCount) {

	timestamp countingStart;

	vertex* butterflyCounts = (vertex *) calloc (sizeof(vertex), rightGraph.size());
	vertex maxBc = countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right

	timestamp peelingStart;
	print_time (fp, "Counting bflies (per vertex) time: ", peelingStart - countingStart);
	cout << "Counting bflies (per vertex) time: " << peelingStart - countingStart << endl;
	printf ("# bflys: %d\n", *bCount);

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

		// combine neigs of neigs of i in distance-2 neighbors
		HashMap<int> dup (0);
		for (auto v : rightGraph[u])
			for (auto w : leftGraph[v])
				dup[w]++;

		for (auto it = dup.begin(); it != dup.end(); it++) {
			int x = it->first;
			int count = it->second;
			if (x != u && count > 1) {
				if (K[x] == -1) {
					int decreaseAmount = nChoosek (count, 2);
					if (nBucket.CurrentValue(x) - decreaseAmount <= val)
						nBucket.DecTo (x, val);
					else
						nBucket.DecTo (x, nBucket.CurrentValue(x) - decreaseAmount);
				}
				else if (hierarchy)
					createSkeleton (u, {x}, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}

		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);

	}

	nBucket.Free();
	*maxbicore = bf_u;

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);
	cout << "Peeling time: " << peelingEnd - peelingStart << endl;

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




// rightGraph is primary, leftGraph is secondary
void oldtipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, long* bCount) {

	timestamp peelingStart;

	vertex* butterflyCounts = (vertex *) calloc (sizeof(vertex), rightGraph.size());
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


