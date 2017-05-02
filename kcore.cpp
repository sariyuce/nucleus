#include "main.h"

void base_kcore (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp) {

	timestamp peelingStart;

	vertex nVtx = graph.size();
	vertex maxDeg = 0;
	for (auto g : graph)
		if (g.size() > maxDeg)
			maxDeg = g.size();

	// Peeling
	K.resize (graph.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize(maxDeg, nVtx);
	for (vertex i = 0; i < graph.size(); i++) {
		if (graph[i].size() > 0)
			nBucket.Insert (i, graph[i].size());
		else
			K[i] = 0;
	}

	vertex deg_u = 0;

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
		component.resize (graph.size(), -1);
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

		deg_u = K[u] = val;

		for (auto v : graph[u]) { // decrease the degree of the neighbors with greater degree
			vertex deg_v = nBucket.CurrentValue(v);
			if (deg_v > deg_u)
				nBucket.DecVal(v);
			else if (hierarchy) // hierarchy related
				createSkeleton (u, {v}, &nSubcores, K, skeleton, component, unassigned, relations);
		}
		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxCore = deg_u; // deg_u is degree of the last popped vertex

	timestamp peelingEnd;
	print_time (fp, "Peeling time: ", peelingEnd - peelingStart);

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif

	if (hierarchy) {
		buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "Nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp;
		presentNuclei (12, skeleton, component, graph, nEdge, hp, vfile, fp);
		timestamp totalEnd;
		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}
