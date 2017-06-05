#include "main.h"



inline int intersection (vector<vertex>& a, vector<vertex>& b) {
	size_t i = 0, j = 0, c = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			c++;
			i++;
			j++;
		}
	}
	return c;
}

// bad, don't use
vertex hopefullyFastercountButterflies (Graph& rightGraph, Graph& leftGraph, vertex* butterflyCounts, vertex* bCount) {

	timestamp t1;
	vertex maxBc = 0;
	long int cost = 0;
	for (vertex i = 0; i < rightGraph.size(); i++) {
		for (vertex j = i + 1; j < rightGraph.size(); j++) {

			cost += min (rightGraph[i].size(), rightGraph[j].size());
			int inter = intersection (rightGraph[i], rightGraph[j]);
			if (inter > 1) {
				int c = nChoosek (inter, 2);
				butterflyCounts[i] += c;
				butterflyCounts[j] += c;

				if (butterflyCounts[i] > maxBc)
					maxBc = butterflyCounts[i];
				if (butterflyCounts[j] > maxBc)
					maxBc = butterflyCounts[j];
				*bCount += c;
			}
		}
	}

	printf ("cost: %ld\n", cost);
	return maxBc;
}


// best, use it
vertex insallahFastercountButterflies (Graph& rightGraph, Graph& leftGraph, vertex* butterflyCounts, vertex* bCount) {

	timestamp t1;
//	long int cost = 0;
	vertex maxBc = 0;
	HashMap<int> dup (0);
	for (vertex i = 0; i < rightGraph.size(); i++) {
		dup.reset (0);
		for (auto v : rightGraph[i])
			for (auto w : leftGraph[v]) {
				if (i < w) { // other order definitions
					dup[w]++;
//					cost++;
				}
			}

		for (auto it = dup.begin(); it != dup.end(); it++) {
//			cost++;
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
//		if (i % 100 == 0)
//			cout << i << " / " << rightGraph.size() << "   after hasmap: " << t3 - t1 << endl;
	}

	*bCount /= 2;

	printf ("bfly: %d\n", *bCount);
	printf ("Left: %d, Right: %d\n", leftGraph.size(), rightGraph.size());
	return maxBc;
}

// old, don't use
vertex countButterflies (Graph& rightGraph, Graph& leftGraph, vertex* butterflyCounts, vertex* bCount) {

	timestamp t1;
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

		timestamp t3;
//		if (i % 100 == 0)
//			cout << i << "after sort: " << t3 - t1 << endl;
	}

	*bCount /= 2;
	return maxBc;
}

// rightGraph is primary, leftGraph is secondary
void tipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {

	timestamp countingStart;

	vertex* butterflyCounts = (vertex *) malloc (sizeof(vertex) * rightGraph.size());
	vertex maxBc = countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right

	timestamp peelingStart;
	print_time (fp, "Counting bflies (per vertex) time: ", peelingStart - countingStart);
	cout << "Counting bflies (per vertex) time: " << peelingStart - countingStart << endl;

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
	cout << "Peeling time: " << peelingEnd - peelingStart << endl;

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxbicore, relations, skeleton, &nSubcores, nEdge, rightGraph.size(), leftGraph.size());
		timestamp nucleusEnd;

		print_time (fp, "Tip decomposition time with hierarchy construction: ", nucleusEnd - countingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), rightGraph.size());

		helpers dummy;
		presentNuclei ("TIP", skeleton, component, nEdge, dummy, vfile, fp, leftGraph, rightGraph, NULL);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - countingStart);
	}
}

// bad don't use
void insallahFasterTipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {

	timestamp countingStart;


	vertex* butterflyCounts = (vertex *) malloc (sizeof(vertex) * rightGraph.size());
	vertex maxBc = countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right
	timestamp peelingStart;

	print_time (fp, "Counting bflies (per vertex) time: ", peelingStart - countingStart);
	cout << "Counting bflies (per vertex) time: " << peelingStart - countingStart << endl;


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
				intersectionM1 (leftGraph[v], leftGraph[w], commons);

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


		for (size_t i = 0; i < rightGraph[u].size(); i++) {
			int v = rightGraph[u][i];
			rightGraph[u][i] = -1;
			for (int j = 0; j < leftGraph[v].size(); j++)
				if (leftGraph[v][j] == u)
					leftGraph[v][j] = -1;
		}

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

		print_time (fp, "Tip decomposition time with hierarchy construction: ", nucleusEnd - countingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), rightGraph.size());

		helpers dummy;
		presentNuclei ("TIP", skeleton, component, nEdge, dummy, vfile, fp, leftGraph, rightGraph, NULL);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - countingStart);
	}
}

// best, use it
void hopefullyFastertipDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {

	timestamp countingStart;

	vertex* butterflyCounts = (vertex *) calloc (sizeof(vertex), rightGraph.size());

//	vertex maxBc = countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right
	vertex maxBc = insallahFastercountButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each vertex on the right

	timestamp peelingStart;

//	FILE* ff = fopen ("bflys", "w");
//	for (int i = 0; i < rightGraph.size(); i++)
//		fprintf (ff, "%d\n", butterflyCounts[i]);
//	fclose (ff);
//
	print_time (fp, "Counting bflies (per vertex) time: ", peelingStart - countingStart);
	cout << "Counting bflies (per vertex) time: " << peelingStart - countingStart << endl;
	printf ("# bflys: %d\n", *bCount);
//	exit(1);

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

		// combine neigs of neigs of i in d2Neighbors
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

}
