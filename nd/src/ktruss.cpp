#include "main.h"

inline int checkConnectedness (Graph& graph, Graph& orientedGraph, Graph& TC, vertex u, vertex v, vector<vertex>* xel = NULL) {
	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);
	for (size_t k = 0; k < orientedGraph[a].size(); k++)
		if (orientedGraph[a][k] == b) {
			TC[a][k]++;
			if (xel == NULL)
				return b;
			else
				return (*xel)[a] + k;
		}
	return -1;
}

// per edge
lol countTriangles (Graph& graph, Graph& orientedGraph, Graph& TC) {
	lol tc = 0;
#ifdef SAVE_TRIS
	FILE* fp = fopen ("tris", "w");
#endif
	for (size_t i = 0; i < orientedGraph.size(); i++) {
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orientedGraph[i].size(); k++) {
				vertex a = orientedGraph[i][j];
				vertex b = orientedGraph[i][k];
				vertex c = checkConnectedness (graph, orientedGraph, TC, a, b);
				if (c != -1) {
#ifdef SAVE_TRIS
					fprintf (fp, "%d %d %d\n", a, b, c);
#endif
					TC[i][j]++;
					TC[i][k]++;
					tc++;
				}
			}
		}
	}
#ifdef SAVE_TRIS
	fclose (fp);
#endif
	return tc;
}

void base_ktruss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orientedGraph;
	createOrientedIndexEdges (graph, el, xel, orientedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orientedGraph[i].size(), 0);
	lol tric;
#ifndef LOAD_TRIS
	// Compute
 	tric = countTriangles (graph, orientedGraph, TC);
#else
 	// Load
	FILE* aa = fopen ("tris", "r");
	int a, b, t;
	tric = 0;
	while (fscanf (aa, "%d %d %d", &a, &b, &t) != EOF) {
		TC[a][b] = t;
		tric += t;
	}
	tric /= 3;
	fclose (aa);
#endif

	fprintf (fp, "# triangles: %lld\n", tric);
	timestamp tcEnd;
	print_time (fp, "Triangle counting: ", tcEnd - peelingStart);

	// Peeling
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orientedGraph.size(); i++)
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	vertex tc_e = 0;

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

	vertex monitor = 0;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
#ifdef MONITOR
		if (monitor % 10000 == 0)
			printf ("e: %d    val: %d    counter: %d  nEdge: %d\n", e, val, monitor, nEdge);
		monitor++;
#endif
		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		vector<vertex> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // decrease the TC of the neighbor edges with greater TC
			edge f = getEdgeId (u, w, xel, el, graph);
			edge g = getEdgeId (v, w, xel, el, graph);
			if (K[f] == -1 && K[g] == -1) {
				if (nBucket.CurrentValue(f) > tc_e)
					nBucket.DecVal(f);
				if (nBucket.CurrentValue(g) > tc_e)
					nBucket.DecVal(g);
			}
			else if (hierarchy)
				createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxtruss = tc_e; // tc_e is tc of the last popped edge

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "2-3 nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;
		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}

// per edge
lol storeCountTriangles (Graph& graph, Graph& orientedGraph, Graph& TC, vector<vertex>& xel, Graph& tris) {
	lol tc = 0;
	for (vertex i = 0; i < orientedGraph.size(); i++) {
		for (vertex j = 0; j < orientedGraph[i].size(); j++) {
			vertex a = orientedGraph[i][j];
			edge x = xel[i] + j;
			for (vertex k = j + 1; k < orientedGraph[i].size(); k++) {
				edge y = xel[i] + k;
				vertex b = orientedGraph[i][k];
				vertex c = checkConnectedness (graph, orientedGraph, TC, a, b, &xel);

				if (c != -1) {
					tc++;
					TC[i][j]++;
					TC[i][k]++;
					vertex z = c;
					tris[x].push_back (y);
					tris[x].push_back (z);
					tris[y].push_back (x);
					tris[y].push_back (z);
					tris[z].push_back (x);
					tris[z].push_back (y);
				}
			}
		}
	}
	return tc;
}

void base_ktruss_storeTriangles (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orientedGraph;
	createOrientedIndexEdges (graph, el, xel, orientedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orientedGraph[i].size(), 0);

	vector<vector<vertex>> tris (el.size());
	lol tric = storeCountTriangles (graph, orientedGraph, TC, xel, tris);

	fprintf (fp, "# triangles: %lld\n", tric);
	timestamp tcEnd;
	print_time (fp, "Triangle storing and counting: ", tcEnd - peelingStart);


	// Peeling
	K.resize (el.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge);
	vertex id = 0;
	for (size_t i = 0; i < orientedGraph.size(); i++)
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	vertex tc_e = 0;

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

	vertex monitor = 0;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
#ifdef MONITOR
		if (monitor % 10000 == 0)
			printf ("e: %d    val: %d    counter: %d  nEdge: %d\n", e, val, monitor, nEdge);
		monitor++;
#endif

		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;

		for (size_t j = 0; j < tris[e].size(); j+=2) {
			vertex f = tris[e][j];
			vertex g = tris[e][j+1];
			if (K[f] == -1 && K[g] == -1) {
				if (nBucket.CurrentValue(f) > tc_e)
					nBucket.DecVal(f);
				if (nBucket.CurrentValue(g) > tc_e)
					nBucket.DecVal(g);
			}
			else if (hierarchy)
				createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();

	*maxtruss = tc_e; // tc_e is tc of the last popped edge

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "2-3 nucleus decomposition (storing triangles) time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;
		print_time (fp, "Total time (storing triangles), including the density computations: ", totalEnd - peelingStart);
	}
}

