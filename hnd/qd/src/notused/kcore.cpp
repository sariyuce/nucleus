#include "main.h"

void base_kcore (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp) {

	const auto p1 = chrono::steady_clock::now();
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

	const auto p2 = chrono::steady_clock::now();
	if (!hierarchy) {
		print_time (fp, "Peeling time: ", p2 - p1);
	}
	else  {
		print_time (fp, "Peeling + on-the-fly hiearchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 1,2 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp;
		presentNuclei (12, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 1,2 nucleus decomposition time: ", (p2 - p1) + (b2 - b1) + (d2 - d1));
	}
}



void degreeBasedHierarchy (Graph& graph, bool hierarchy, edge nEdge, vector<vertex>& K, vertex* maxCore, string vfile, FILE* fp) {

	const auto p1 = chrono::steady_clock::now();
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
				; // nBucket.DecVal(v);
			else if (hierarchy) // hierarchy related
				createSkeleton (u, {v}, &nSubcores, K, skeleton, component, unassigned, relations);
		}
		if (hierarchy)
			updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxCore = deg_u; // deg_u is degree of the last popped vertex

	const auto p2 = chrono::steady_clock::now();
	if (!hierarchy) {
		print_time (fp, "Peeling time: ", p2 - p1);
	}
	else  {
		print_time (fp, "Peeling + on-the-fly hiearchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 1,2 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp;
		presentNuclei (12, skeleton, component, graph, nEdge, hp, vfile, fp);
		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 1,2 nucleus decomposition time: ", (p2 - p1) + (b2 - b1) + (d2 - d1));
	}



}

//	report_all_stuff2 (12, graph, nEdge, *core_number, ax, exists, K, vfile, fp);
//}

/*


void find_subcore (vertex root, Graph& graph, vector<vertex>& K, bool* visited, vector<vertex>& subcore_ids, vector<subcore>& hrc) {

	vertex cn = K[root];
	int next_id = hrc.size();
	subcore sc (cn);
	sc.rank = 0;
	hrc.push_back (sc);

	vector<int> acc; // accumulate during find
	vector<vertex> vs; // put the vertices in subcore
	vector<int> merge_sc; // put the subcores to merge
	HashMap<bool> marked (false); // to keep track of the neighbor subcores (and their ancestors)


	queue<vertex> bfsorder; // we are doing bfs
	bfsorder.push(root);
	visited[root] = true;

	while (!bfsorder.empty()) {
		vertex v = bfsorder.front();
		bfsorder.pop();
		vs.push_back(v);
		for (size_t i = 0; i < graph[v].size(); i++) {
			vertex w = graph[v][i];
			// discovering an unvisited vertex with same K value: put it to subcore
			if (K[w] == cn && !visited[w]) {
				bfsorder.push(w);
				visited[w] = true;
			}
			// discovering a vertex with larger K value. It's visited for sure since we visit them in decreasing order of K values
			else if (K[w] > cn) {
				int s = subcore_ids[w];
				if (marked.find(s) == marked.end()) { // check if you see this subcore before in this traversal
					int old_s = s;
					// find the root at that point, it  and compress same K subcores in the mean time
					acc.clear();
					while (hrc[s].root != -1) {
						acc.push_back (s);
						s = hrc[s].root;
					}
					for (int i : acc)
						hrc[i].root = s;

					if (marked.hasDefaultValue(s) && s != next_id) { // now you have the root as s
						if (hrc[s].K > cn) { // s becomes kid of next_id-th subcore
							hrc[s].parent = next_id;
							hrc[s].root = next_id;
						}
						else // hrc[s].K == cn
							merge_sc.push_back (s); // defer the merge processing to the end: need to know the ranks

						marked[s] = true;
					}
					marked[old_s] = true;
				}
			}
		}
	}

//	printf ("subcore K: %d |V|: %d\n", cn, vs.size());
	for (size_t i = 0; i < vs.size(); i++)
		subcore_ids[vs[i]] = next_id;
#ifndef SHORT
	hrc[next_id].size = vs.size();
#endif
	op.total_size += vs.size();

	if (!merge_sc.empty()) {
		merge_sc.push_back (next_id);
		vector<bool> seen (merge_sc.size(), false);
		for (int i = 0; i < merge_sc.size(); i++) {
			if (!seen[i]) {
				for (int j = i + 1; j < merge_sc.size(); j++) {
					if (!seen[i] && !seen[j]) {
						op.union_op++;
						int ch = merge_sc[i];
						int pr = merge_sc[j];
						bool fl = false;
						// note that we already found the topmost ancestors above in **
						if (hrc[ch].rank > hrc[pr].rank) {
							int t = pr;
							pr = ch;
							ch = t;
							seen[j] = true;
						}
						else {
							seen[i] = true;
							fl = true;
						}
						hrc[ch].parent = pr;
						hrc[ch].visible = false;
						hrc[ch].root = pr;
						if (hrc[pr].rank == hrc[ch].rank)
							hrc[pr].rank++;
						if (fl)
							break;
					}
				}
			}
		}
	}
}

void report_all_stuff2 (int variant, Graph& graph, int nEdge, int cn, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile, FILE* ffp) {

	vector<subcore> hrc;
	size_t sz;
//	if (variant == 34)
//		sz = (*ax.tl).size();
//	else if (variant == 23)
//		sz = (*ax.el).size();
//	else
	sz = graph.size();

	vector<vertex> subx_ids (sz, -1);

	bool* visited = (bool *) malloc (sizeof(bool) * sz);
	double prev_rt = 0;
	int old_size = 0;
	int iter = 0;

	for (vertex i = cn; i > 0; i--) { // iterate over all node values

		if (exists.hasDefaultValue(i))
			continue;

		int cnt = 0;

		iter++;

		for (size_t j = 0; j < sz; j++) {
			if (!visited[j] && K[j] == i) {
				cnt++;

				double rt = visited.size() * 100 / (double) sz;
				if (rt - prev_rt >= 0.1) {
					prev_rt = rt;
					printf ("%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					fprintf (ffp, "%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					time_t now;
					time(&now);
					printf("CENTRAL TIME: %s", ctime(&now));
					fflush (ffp);
				}

//				if (variant == 34)
//					find_sub34 (j, graph, K, (*ax.tlist), (*ax.tl), visited, subx_ids, hrc, op);
//				else if (variant == 23)
//					find_subtruss (j, graph, K, (*ax.xel), (*ax.el), visited, subx_ids, hrc, op);
//				else
				find_subcore (j, graph, K, visited, subx_ids, hrc, tmerge, tloop, tdirect, op);

			}
		}
		timestamp n2;
		cout << "K: " << i << " count: " << cnt << " size: " << op.total_size - old_size << " total size: " << op.total_size << " time: " << n2 - n1 << " total time: " << n2 - t1 << endl;
	}

	fprintf (ffp, "iter: %d\n", iter);
	timestamp t2;
	print_time (ffp, "all subthings finding time: ", t2 - t1);
	cout << "all subthings finding time: " << t2 - t1 << endl;

	// insert the root as a parent of all parentless subcores
	int nid = hrc.size();
	subcore sc;
	sc.rank = 1;
	sc.K = 0;
	sc.parent = -1;
	sc.root = -1;
	sc.visible = true;
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].parent == -1) {
			hrc[i].parent = nid;
			hrc[i].root = nid;
			sc.children.push_back(i);
		}
	}


#ifndef SHORT
	sc.nedge = nEdge;
	sc.ed = nEdge / double (graph.size() * (graph.size() - 1) / 2);
#endif
	hrc.push_back (sc);

	timestamp t3;
	print_time (ffp, "root core insertion time: ", t3 - t2);
	cout << "root core insertion time: " << t3 - t2 << endl;

#if 1
	long totalmem = 0;
	for (int i = 0; i < hrc.size(); i++) {
		totalmem += sizeof(bool) + 4 * sizeof(int) + hrc[i].children.size() * sizeof(int);
	}
	printf ("hrc.size: %d , and total mem: %ld\n", hrc.size(), totalmem);
#endif


	report_subgraphs (variant, hrc, subx_ids, graph, nEdge, ax, vfile, ffp);
}


*/
