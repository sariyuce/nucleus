#include "main.h"

#define EL(A, B) (xRight[A]+B)

#define REL(A, B) (xRight[A]+B)

#define LEL(A, B) (xLeft[A]+B)

inline int Lel (int v, int ind, Graph& leftGraph, Graph& rightGraph, vector<vertex>& xRight) {
	int u = leftGraph[v][ind];
	for (int i = 0; i < rightGraph[u].size(); i++)
		if (rightGraph[u][i] == v)
			return REL(u,i);
	printf ("OOOOO\n");
	exit(1);
}

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

	timestamp t1;
	int co = 0;
	for (vertex i = 0; i < leftGraph.size(); i++) {
		vertex a = i;
		for (vertex j = 0; j < leftGraph[i].size(); j++) {
			vertex b = leftGraph[i][j];
			for (vertex k = j + 1; k < leftGraph[i].size(); k++) {
				// intersection set with greater than a's
				int c = leftGraph[i][k];
				limitedIntersection (rightGraph, b, c, butterflyCounts, a, bCount);
				co++;
//				if (co % 100 == 0) {
//					timestamp t2;
//					cout << "co: " << co << "  time: " << t2 - t1 << " i: " << i << endl;
//				}
			}
		}
	}
}


void insallahFasterCountButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, vertex* bCount) {

	timestamp t1;
	int co = 0;
	for (vertex i = 0; i < rightGraph.size(); i++) {
		if (rightGraph[i].size() > 1) {
			for (vertex j = i + 1; j < rightGraph.size(); j++) {
				if (rightGraph[j].size() > 1) {
					vector<vertex> commons;
					indicesIntersection (rightGraph[i], rightGraph[j], commons, -1);
					if (commons.empty())
						continue;
					int sz = commons.size() / 2;
					for (int k = 0; k < commons.size(); k+=2) {
						butterflyCounts[i][commons[k]] += sz - 1;
						butterflyCounts[j][commons[k+1]] += sz - 1;
						co++;
//						if (co % 100 == 0) {
//							timestamp t2;
//							cout << "co: " << co << "  time: " << t2 - t1 << " i: " << i << endl;
//						}

					}
				}
			}
		}
	}
}

void hopefullyFasterCountButterflies (Graph& rightGraph, Graph& leftGraph, Graph& butterflyCounts, vertex* bCount) {

	timestamp t1;
	HashMap<int> marked (-1);
	for (vertex i = 0; i < rightGraph.size(); i++) {
		marked.reset (-1);
		// mark all nu's
		int u = i;
		for (vertex j = 0; j < rightGraph[i].size(); j++) {
			int nu = rightGraph[i][j];
			marked[nu] = 0;
		}

		for (vertex j = 0; j < rightGraph[i].size(); j++) {
			int v = rightGraph[i][j];
			int count = 0;
			// get distance-2 neighbors of v with duplicate counts
			for (auto nv : leftGraph[v]) {
				if (nv != u) {
					for (auto n2v : rightGraph[nv]) {
						if (!(marked.hasDefaultValue(n2v)) && n2v != v) {
							count++;
//							if (count % 100000 == 0) {
//								timestamp t2;
//								cout << "co: " << count << "  time: " << t2 - t1 << " i: " << i << endl;
//							}
						}
					}
				}
			}

			butterflyCounts[u][j] = count;
			*bCount += count;
		}
	}

	*bCount /= 4;

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
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge

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



void insallahFasterwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {


	vector<vp> left_el;
	vector<vertex> xLeft;

	prefixSum (xLeft, leftGraph, left_el);

	timestamp countingStart;


	Graph butterflyCounts (rightGraph.size());
	for (vertex i = 0; i < rightGraph.size(); i++)
		butterflyCounts[i].resize (rightGraph[i].size(), 0);
	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge
//	insallahFasterCountButterflies (leftGraph, rightGraph, butterflyCounts, bCount); // counts butterflies for each edge
//	hopefullyFasterCountButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge

	timestamp countingEnd;
//	cout << "total time: " << countingEnd - countingStart << endl;
//	FILE* ffp = fopen ("butterflies", "w");
//	for (vertex i = 0; i < rightGraph.size(); i++) {
//		for (vertex j = 0; j < rightGraph[i].size(); j++) {
//			fprintf (ffp, "BF %d %d    %d\n", i, rightGraph[i][j], butterflyCounts[i][j]);
//		}
//	}
//	fprintf (ffp, "bCount: %d\n", *bCount);
//	fclose (fp);
//	exit(1);




	vertex maxBc = 0;
	for (auto g: butterflyCounts)
		for (auto c: g)
			if (c > maxBc)
				maxBc = c;

	timestamp peelingStart;
	cout << "counting time: " << peelingStart - countingStart << endl;
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

	timestamp tint (0, 0);
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

//		printf ("e: %d, val: %d\n", e, val);
		vertex u = el[e].first;
		vertex v = el[e].second;
		vertex uvInd = getEdgeIndex (u, v, el, xRight);
		vertex vuInd = getEdgeIndex (v, u, left_el, xLeft);


		// v is left vertex, u is right vertex, t is a right vertex and v's neighbor, w is a left vertex and neighbor to u and t
		// xyInd is the number i s.t. blaGraph[x][i] = y
		for (vertex vtInd = 0; vtInd < leftGraph[v].size(); vtInd++) {
			vertex t = leftGraph[v][vtInd];
			if (t == u || t == -1)
				continue;
//			printf ("v: %d, t: %d\n", v, t);
			vertex tvInd = getEdgeIndex (t, v, el, xRight);
			vertex g = EL(t, tvInd);
//			if (K[g] == -1) {
				size_t i = 0, j = 0;
				while (i < rightGraph[u].size() && j < rightGraph[t].size()) {

//					printf ("first: %d   i:%d, second: %d   j: %d\n", rightGraph[u][i], i, rightGraph[t][j], j);

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
						// we have u, t on right; v, w on left; indices from u, t to v, w
						vertex f = EL(u, uwInd);
						vertex h = EL(t, twInd);

//						if (K[f] == -1 && K[h] == -1) {
							if (nBucket.CurrentValue (f) > val)
								nBucket.DecVal (f);
							if (nBucket.CurrentValue (g) > val)
								nBucket.DecVal (g);
							if (nBucket.CurrentValue (h) > val)
								nBucket.DecVal (h);
//						}

						i++;
						j++;
					}
					else {
						i++;
						j++;
					}
				}
			}
//			}

		rightGraph[u][uvInd] = leftGraph[v][vuInd] = -1;
	}

	nBucket.Free();
	*maxbicore = bf_e;

	timestamp peelingEnd;
	cout << "peeling time: " << peelingEnd - peelingStart << endl;
	cout << "intersection time: " << tint << endl;
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

//
//void hopefullyFasterwingDecomposition (Graph& leftGraph, Graph& rightGraph, edge nEdge, vector<vertex>& K, bool hierarchy, vector<vp>& el, vector<vertex>& xRight, vertex* maxbicore, string vfile, FILE* fp, vertex* bCount) {
//
//
////	vector<vp> left_el;
////	vector<vertex> xLeft;
////
////	prefixSum (xLeft, leftGraph, left_el);
//
//
//	timestamp peelingStart;
//
//	Graph butterflyCounts (rightGraph.size());
//	for (vertex i = 0; i < rightGraph.size(); i++)
//		butterflyCounts[i].resize (rightGraph[i].size(), 0);
//	countButterflies (rightGraph, leftGraph, butterflyCounts, bCount); // counts butterflies for each edge
//
//	vertex maxBc = 0;
//	for (auto g: butterflyCounts)
//		for (auto c: g)
//			if (c > maxBc)
//				maxBc = c;
//
//	printf ("el.size(): %d\n", el.size());
//	// peeling
//	K.resize (el.size(), -1);
//	Naive_Bucket nBucket;
//	nBucket.Initialize (maxBc, nEdge);
//	vertex bid = 0;
//
//	for (vertex i = 0; i < butterflyCounts.size(); i++) {
//		for (vertex j = 0; j < butterflyCounts[i].size(); j++) {
//			vertex c = butterflyCounts[i][j];
//			if (c > 0)
//				nBucket.Insert (bid++, c);
//			else
//				K[bid++] = 0;
//		}
//	}
//
//	vertex bf_e = 0;
//
//	// required for hierarchy
//	vertex cid; // subcore id number
//	vector<subcore> skeleton; // equal K valued cores
//	vector<vertex> component; // subcore ids for each vertex
//	vector<vp> relations;
//	vector<vertex> unassigned;
//	vertex nSubcores;
//
//	if (hierarchy) {
//		cid = 0;
//		nSubcores = 0;
//		component.resize (el.size(), -1);
//	}
//
//	timestamp lt (0, 0);
//
//	while (true) {
//		edge e;
//		vertex val;
//		if (nBucket.PopMin(&e, &val) == -1)
//			break;
//
//		if (hierarchy) {
//			unassigned.clear();
//			subcore sc (val);
//			skeleton.push_back (sc);
//		}
//
//		bf_e = K[e] = val;
//
//		vertex u = el[e].first;
//		vertex v = el[e].second;
//		vertex uvInd = getEdgeIndex (u, v, el, xRight);
//
//		// mark & collect u's neighbors
//		HashMap<int> nvMap (-1);
//		vector<vertex> nvs;
//		for (int i = 0; i < leftGraph[v].size(); i++) {
////			if (K[REL(u, i)] == -1) {
//				nvMap[leftGraph[v][i]] = i;
////			}
//		}
//
//
//		unordered_multimap<int, int> eda;
//		for (int i = 0; i < rightGraph[u].size(); i++) {
//			if (K[REL(u, i)] == -1) {
//				int nu = rightGraph[u][i];
//				for (int j = 0; j < leftGraph[nu].size(); j++) {
//					int nv = leftGraph[nu][j]
//					if (!(nvMap.hasDefaultValue (nv))) {
//						eda.emplace (make_pair (nv, i));
//					}
//				}
//			}
//		}
//
//
//
//
//		// First phase
//		HashMap<int> myMap (0);
//		for (int i = 0; i < nus.size(); i++) {
//			int nu = rightGraph[u][nus[i]];
//			for (int j = 0; j < leftGraph[nu].size(); j++) {
//				int nnu = leftGraph[nu][j];
//				if (nnu != u) {
////					timestamp t1;
////					int x = Lel(nu, j, leftGraph, rightGraph, xRight);
////					timestamp t2;
////					lt += t2 - t1;
////					if (K[x] == -1) {
////						printf ("index 3 is %d\n", Lel(nu, j, leftGraph, rightGraph, xRight));
//						myMap[nnu]++;
////					}
//				}
//			}
//		}
//
//		int willDec = 0;
//		for (int i = 0; i < nvs.size(); i++) {
//			int nv = leftGraph[v][nvs[i]];
//			if (!(myMap.hasDefaultValue (nv))) {
//
//				timestamp t1;
//				int x = Lel(v, nvs[i], leftGraph, rightGraph, xRight);
//				timestamp t2;
//				lt += t2 - t1;
//
////				printf ("index 4 is %d with %d\n", x, myMap[nv]);
//				 // V - NV EDGE , WILL BE NV - V EDGE
////				if (nBucket.CurrentValue (x) == -1) {
////					printf ("HOOOOOO, x: %d and -1", x);
////				}
//				if (nBucket.CurrentValue (x) - myMap[nv] <= val)
//					nBucket.DecTo (x, val);
//				else
//					nBucket.DecTo (x, nBucket.CurrentValue(x) - myMap[nv]);
//			}
//		}
//
//		// Second phase
//		myMap.reset (0);
//		for (int i = 0; i < nvs.size(); i++) {
//			int nv = leftGraph[v][nvs[i]];
//			for (int j = 0; j < rightGraph[nv].size(); j++) {
//				int nnv = rightGraph[nv][j];
//				if (nnv != v) {
//					if (K[REL(nv, j)] == -1) {
////						printf ("index 6 is %d\n", REL(nv, j));
//						myMap[nnv]++;
//					}
//				}
//			}
//		}
//
//		for (int i = 0; i < nus.size(); i++) {
//			int nu = rightGraph[u][nus[i]];
//			if (!(myMap.hasDefaultValue (nu))) {
//
//				int x = REL(u, nus[i]);
////				printf ("index 7 is %d with %d\n", x, myMap[nu]);
//				 // U - NU EDGE
////				if (nBucket.CurrentValue (x) == -1) {
////					printf ("HOOOOOO, x: %d and -1", x);
////				}
//
//				if (nBucket.CurrentValue (x) - myMap[nu] <= val)
//					nBucket.DecTo (x, val);
//				else
//					nBucket.DecTo (x, nBucket.CurrentValue(x) - myMap[nu]);
//			}
//		}
//
//
//		// Third phase
//		for (int i = 0; i < nus.size(); i++) {
//			int nu = rightGraph[u][nus[i]];
//			for (int j = 0; j < nvs.size(); j++) {
//				int nv = leftGraph[v][nvs[j]];
//
//				bool itsthere = false;
//				int k = 0;
//				for (k = 0; k < rightGraph[nv].size(); k++)
//					if (rightGraph[nv][k] == nu) {
//						itsthere = true;
//						break;
//					}
//				if (itsthere) {
//					int x = REL(nv, k);
//					if (K[x] == -1) {
////						if (nBucket.CurrentValue (x) == -1) {
////							printf ("HOOOOOO, x: %d and -1", x);
////						}
//
//						if (nBucket.CurrentValue (x) > val) {
//							nBucket.DecVal (x); // NU - NV EDGE
////							printf ("index 8 is %d\n", x);
//						}
//					}
//				}
//			}
//		}
//
////		printf ("%d is assigned %d\n", e, val);
//	}
//
//	nBucket.Free();
//	*maxbicore = bf_e;
//
//	timestamp peelingEnd;
//	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);
//
//	cout << "lt time: " << lt << endl;
//	cout << "all time: " << peelingEnd - peelingStart << endl;
// #ifdef K_VALUES
//	for (int i = 0; i < K.size(); i++)
//		printf ("K[%d]: %d\n", i, K[i]);
//#endif
//
//	if (hierarchy) {
//		buildHierarchy (*maxbicore, relations, skeleton, &nSubcores, nEdge, rightGraph.size(), leftGraph.size());
//		timestamp nucleusEnd;
//
//		print_time (fp, "Wing decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
//		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d (in edges) \t\t |E|: %d\n", nSubcores, skeleton.size(), nEdge);
//
//		helpers hp (&el);
//		presentNuclei ("WING", skeleton, component, nEdge, hp, vfile, fp, leftGraph, rightGraph, &xRight);
//		timestamp totalEnd;
//
//		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
//	}
//}
