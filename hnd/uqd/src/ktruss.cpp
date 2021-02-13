#include "main.h"

inline int checkConnectedness (Graph& graph, Graph& orderedGraph, Graph& TC, vertex u, vertex v, vector<vertex>* xel = NULL) {

	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);
	for (size_t k = 0; k < orderedGraph[a].size(); k++)
		if (orderedGraph[a][k] == b) {
			TC[a][k]++;
			if (xel == NULL)
				return b;
			else
				return (*xel)[a] + k;
		}
	return -1;
}

double simple_count_triangles (Graph& graph, Graph& orderedGraph, Graph& TC, unordered_map<int, bool>& numbers, unordered_map<int, bool>& crossing) {
	vertex cut = 0, vol = 0, count = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				vertex a = orderedGraph[i][j];
				vertex b = orderedGraph[i][k];
				vertex c = checkConnectedness (graph, orderedGraph, TC, a, b);
				if (c != -1) {
					int asdf = 0;
					if (crossing.find(a) != crossing.end())
						asdf++;
					if (crossing.find(b) != crossing.end())
						asdf++;
					if (crossing.find(c) != crossing.end())
						asdf++;

	//				printf ("u: %d v: %d w: %d -- asdf: %d\n", u, v, w, asdf);

					if (asdf < 3) {
						if (numbers.find(a) != numbers.end())
							vol++;
						if (numbers.find(b) != numbers.end())
							vol++;
						if (numbers.find(c) != numbers.end())
							vol++;
						if (asdf > 0)
							cut++;
						if (asdf == 0)
							count++;
					}
				}
			}
		}
	}
	printf ("cut: %d vol: %d -- cond: %lf\n", cut, vol, ((double) cut) / vol);

	int nv = numbers.size();
	printf ("count: %d numOfNodes: %d -- avg. motif degree: %lf\n", count, nv, ((double) count) / nv);

	return ((double) cut) / vol;
}


//#define SAVE_TRIS_FOR_NODE_ATTRIBUTES
// per edge
lol countTriangles (Graph& graph, Graph& orderedGraph, Graph& TC) {

	lol tc = 0;
#ifdef SAVE_TRIS_FOR_NODE_ATTRIBUTES
	FILE* gens = fopen("brown_meta", "r");
	int* gender = (int *) malloc (sizeof(int) * graph.size());
	int i=1, a;
	while (fscanf (gens, "%d %d %d %d %d %d %d", &a, &(gender[i++]), &a, &a, &a, &a, &a) != EOF)
		;
	fclose (gens);
	FILE* fp = fopen ("tris", "w");
#endif
	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				vertex a = orderedGraph[i][j];
				vertex b = orderedGraph[i][k];
				vertex c = checkConnectedness (graph, orderedGraph, TC, a, b);
				if (c != -1) {
#ifdef SAVE_TRIS_FOR_NODE_ATTRIBUTES
					fprintf (fp, "%d %d %d\n", gender[a], gender[b], gender[c]);
#endif
					TC[i][j]++;
					TC[i][k]++;
					tc++;
				}
			}
		}
	}
#ifdef SAVE_TRIS_FOR_NODE_ATTRIBUTES
	fclose (fp);
	exit(1);
#endif
	return tc;
}

void base_ktruss (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);


	//for (int i = 0; i < el.size(); i++)
	//	printf ("el: %d %d\n", el[i].first, el[i].second);


	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	lol tric;
#ifndef LOAD_TRIS
	// Compute
	tric = countTriangles (graph, orderedGraph, TC);
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
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Triangle counting: ", t2 - t1);
	
	printf ("triangle count: %lld\n", tric);

	if (COUNT_ONLY)
		return;
#ifdef DUMP_K
	string kfile = vfile + "_freq_values";
        FILE* dp = fopen (kfile.c_str(), "w");
#endif

    	vector<int> dist (nVtx+1, 0);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	printf ("nEdge: %d\n", nEdge);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
#ifdef DUMP_K
	        fprintf (dp, "%d\n", TC[i][j]);
#endif
			if (TC[i][j] > 0) {
				nBucket.Insert (id++, TC[i][j]);
				dist[TC[i][j]]++;
			}
			else
				K[id++] = 0;
		}

#ifdef DUMP_K
	fclose (dp);
#endif
	vertex tc_e = 0;

	const auto c1 = chrono::steady_clock::now();
	if (!hierarchy) {
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
	*maxtruss = tc_e;

	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		auto tm = p2 - p1 - (c2 - c1);
		print_time (fp, "Only peeling time: ", tm);
		print_time (fp, "Total time: ", tm + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 2,3 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 2,3 nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}
}


void tcBasedHierarchy (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	const auto t1 = chrono::steady_clock::now();
	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);


	//for (int i = 0; i < el.size(); i++)
	//	printf ("el: %d %d\n", el[i].first, el[i].second);


	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	lol tric;
#ifndef LOAD_TRIS
	// Compute
	tric = countTriangles (graph, orderedGraph, TC);
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
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Triangle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	printf ("nEdge: %d\n", nEdge);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
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
					; // nBucket.DecVal(f);
				if (nBucket.CurrentValue(g) > tc_e)
					; // nBucket.DecVal(g);
			}
			else if (hierarchy)
				createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxtruss = tc_e;

	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 2,3 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 2,3 nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}
}


// per edge
lol storeCountTriangles (Graph& graph, Graph& orderedGraph, Graph& TC, vector<vertex>& xel, Graph& tris) {

	lol tc = 0;
	for (vertex i = 0; i < orderedGraph.size(); i++) {
		for (vertex j = 0; j < orderedGraph[i].size(); j++) {
			vertex a = orderedGraph[i][j];
			edge x = xel[i] + j;
			for (vertex k = j + 1; k < orderedGraph[i].size(); k++) {
				edge y = xel[i] + k;
				vertex b = orderedGraph[i][k];
				vertex c = checkConnectedness (graph, orderedGraph, TC, a, b, &xel);

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

	const auto t1 = chrono::steady_clock::now();
	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	vector<vector<vertex>> tris (el.size());
	lol tric = storeCountTriangles (graph, orderedGraph, TC, xel, tris);

	fprintf (fp, "# triangles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Triangle storing and counting: ", t2 - t1);


	// Peeling
	const auto p1 = chrono::steady_clock::now();
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge);
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
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

	const auto p2 = chrono::steady_clock::now();

	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 2,3 nucleus decomposition (with stored triangles) time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		presentNuclei (23, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 2,3 nucleus decomposition (with stored triangles) time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}
}


