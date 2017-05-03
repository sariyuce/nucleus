#include "main.h"
//#define IMPR

inline bool intersect (int variant, vertex u, vertex v, Graph& graph) {
	if (variant == 12)
		return true;
	else if (variant == 13)
		return intersects13 (u, v, graph);
	else if (variant == 14)
		return intersects14 (u, v, graph);
	else {
		printf("invalid variant\n");
		exit(1);
	}
}

bool sort_rep (report s, report t) {
	return s.max < t.max;
}

inline bool exists (int val, vector<int>& v) {
	for (size_t i = 0; i < v.size(); i++)
		if (v[i] == val)
			return true;
	return false;
}

void find_subcore (vertex root, Graph& graph, vector<vertex>& K,
		HashMap<bool>& visited, vector<vertex>& subcore_ids, vector<subcore>& hrc,
		util::timestamp& asstime, util::timestamp& findtime, util::timestamp& adjusttime, stats& op) {
	vertex cn = K[root];
	int next_id = hrc.size();
	subcore sc;
	sc.rank = 0;
	sc.K = cn;
	sc.parent = -1;
	sc.root = -1;
	sc.visible = true;
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
			if (K[w] == cn && visited.hasDefaultValue(w)) {
				bfsorder.push(w);
				visited[w] = true;
			}
			// discovering a vertex with larger K value
			else if (K[w] > cn) {
				int s = subcore_ids[w];
				if (marked.hasDefaultValue(s)) { // check if you see this subcore before in this traversal
					int old_s = s;
					// find the root and compress same K subcores in the mean time
#ifdef NOROOT
					while (hrc[s].parent != -1) {
						int lcn = hrc[s].K; // the K value of s
						int n = hrc[s].parent;
						if (hrc[n].K == lcn) {
							if (hrc[n].parent != -1) {
								int g = hrc[n].parent;
								if (hrc[g].K == lcn) // if you have 3 level of subcores with last_cn
									hrc[s].parent = g; // *the extra line of code*
								s = g;
							}
							else {
								s = n; // n is the root
								break;
							}
						}
						else // hrc[n].K < lcn
							s = n;
					}
#else
//					op.find_op++;

					acc.clear();
					while (hrc[s].root != -1) {
						acc.push_back (s);
						s = hrc[s].root;
					}
					for (int i : acc)
						hrc[i].root = s;

//					while (hrc[s].root != -1) { // **
//						int r = hrc[s].root;
//						if (hrc[r].root != -1)
//							hrc[s].root = hrc[r].root;
//						s = r;
//					}
#endif
					if (marked.hasDefaultValue(s) && s != next_id) { // now you have the root as s
						if (hrc[s].K > cn) { // s becomes kid of next_id-th subcore
							hrc[s].parent = next_id;
//							hrc[next_id].children.push_back (s);
//							hrc[next_id].rank += hrc[s].rank;
#ifndef NOROOT
							hrc[s].root = next_id;
//							op.adjust_op += hrc[s].children.size();
//							for (int ch : hrc[s].children)
//								hrc[ch].root = next_id;
#endif
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
#ifndef NOROOT
						hrc[ch].root = pr;
//						op.adjust_op += hrc[ch].children.size();
//						for (int kid : hrc[ch].children)
//							hrc[kid].root = pr;
#endif
//						hrc[pr].children.push_back (ch);
						if (hrc[pr].rank == hrc[ch].rank)
							hrc[pr].rank++;
//						hrc[pr].rank += hrc[ch].rank;
						if (fl)
							break;
					}
				}
			}
		}
	}
}

bool maxcore (vertex root, vertex core_of_root, Graph& graph, vector<vertex>& K, vector<vertex>& vertices, HashMap<bool>& visited) {
	queue<vertex> bfsorder;
	visited[root] = true;
	bfsorder.push(root);
	bool maximal = true;
	while (!bfsorder.empty()) {
		vertex v = bfsorder.front();
		bfsorder.pop();
		vertices.push_back(v);
		for (size_t i = 0; i < graph[v].size(); i++) {
			vertex w = graph[v][i];
			if (K[w] >= core_of_root) {
				if (visited.hasDefaultValue(w)) {
					bfsorder.push(w);
					visited[w] = true;
					if (K[w] > core_of_root)
						maximal = false;
				}
			}
		}
	}
	return maximal;
}

bool myfunction (co i, co j) {
	if (i.c != j.c)
		return (i.c < j.c);
	else
		return (i.id < j.id);
}

inline void find_root (int* ch, vector<subcore>& ecr) {
	int s = *ch;
	while (ecr[s].root != -1) { // **
		int r = ecr[s].root;
		if (ecr[r].root != -1)
			ecr[s].root = ecr[r].root;
		s = r;
	}
	*ch = s;
}

//inline void find_same_K_parent (int* ch, vector<subcore>& ecr) {
//	int u = *ch;
//	while (ecr[u].parent != -1) {
//		int lcn = ecr[u].K; // the K value of ch
//		int n = ecr[u].parent;
//		if (ecr[n].K == lcn) {
//			if (ecr[n].parent != -1) {
//				int g = ecr[n].parent;
//				if (ecr[g].K == lcn) { // if you have 3 level of subcores with lcn
//					ecr[u].parent = g; // *the extra line of code*
//					u = g;
//				}
//				else {
//					u = n;
//					break;
//				}
//			}
//			else {
//				u = n; // n is the root
//				break;
//			}
//		}
//		else // hrc[n].K < lcn
//			break;
//	}
//	*ch = u;
//}

// 2-pass path compression
inline void find_same_K_parent1 (int* ch, vector<subcore>& ecr) {
	int u = *ch;
	vector<int> vs;
	while (ecr[u].parent != -1) {
		int n = ecr[u].parent;
		if (ecr[n].K == ecr[u].K) {
			vs.push_back (u);
			u = n;
		}
		else
			break;
	}
	*ch = u;
	for (int i : vs) {
		if (i != u)
			ecr[i].parent = u;
	}
}


inline void find_same_K_root (int* ch, vector<subcore>& ecr) {
	int u = *ch;
	while (ecr[u].root != -1) {
		int lcn = ecr[u].K; // the K value of ch
		int n = ecr[u].root;
		if (ecr[n].K == lcn) {
			if (ecr[n].root != -1) {
				int g = ecr[n].root;
				if (ecr[g].K == lcn) { // if you have 3 level of subcores with lcn
					ecr[u].root = g; // *the extra line of code*
					u = g;
				}
				else {
					u = n;
					break;
				}
			}
			else {
				u = n; // n is the root
				break;
			}
		}
		else // hrc[n].K < lcn
			break;
	}
	*ch = u;
}

void traversal_core (Graph& graph) {
	timestamp q1;
	printf ("traversal truss, graph.size(): %d\n", graph.size());
	vector<bool> visited (graph.size(), false);

	int cn = 0;
	for (int i = 0; i < graph.size(); i++) {
//		if (graph[i].size() > 0 && visited.hasDefaultValue(i)) {
		if (visited[i] == false) {
			visited[i] = true;
			queue<int> bfsorder;
			bfsorder.push(i);
			while (!bfsorder.empty()) {
				int c = bfsorder.front();
				bfsorder.pop();
				cn++;
				for (size_t j = 0; j < graph[c].size(); j++) {
					vertex w = graph[c][j];
					if (visited[w] == false) {
						visited[w] = true;
						bfsorder.push(w);
					}
				}
			}
		}
	}
	printf ("cn: %d\n", cn);
	timestamp q2;
	cout << "plain traverse: " << q2 - q1 << endl;

}



void base_kcore (Graph& graph, int nEdge, vector<vertex>& K, util::timestamp& totaltime, int *core_number, const char* vfile, FILE* fp) {

	util::timestamp t1;
	time_t now;
	time(&now);
	fprintf (fp, "%s\n", vfile);
	fflush (fp);
	// Initial declarations
	char *p = (char*) &totaltime;
	util::timestamp *pout = (util::timestamp*) p;

	p_auxies ax;
	size_t nVtx = graph.size();
	size_t max_degree = 0;
	for (size_t i = 0; i < nVtx; i++) {
		if (graph[i].size() > max_degree)
			max_degree = graph[i].size();
	}

	// Peeling
	K.resize (graph.size());
	Naive_Bucket na_bs;
	na_bs.Initialize(max_degree, (graph.size()));
	for (size_t i = 0; i < nVtx; i++)
		na_bs.Insert (i, graph[i].size());

	vertex degree_of_u = 0;
	HashMap<bool> exists (false);

	int id = 0; // subcore id number
	vector<subcore> ecr; // equal K valued cores
	vector<int> compid (graph.size(), -1); // subcore ids for each vertex
	vector<couple1> relations;
	HashMap<bool> seen (false);
	vector<int> minusones;
	int num_subcores = 0;
	while (1) {
		int u, value_of_u;
		int ret = na_bs.PopMin(&u, &value_of_u);
		if (ret == -1) // if the bucket is empty
			break;

		if (value_of_u == 0)
			continue;

		seen.reset (false);
		minusones.clear();
		subcore sc;
		sc.rank = 0;
		sc.K = value_of_u;
		sc.parent = -1;
		sc.root = -1;
		sc.visible = true;
		ecr.push_back (sc);

		degree_of_u = K[u] = value_of_u;
		exists[value_of_u] = true;

		for (size_t j = 0; j < graph[u].size(); j++) { /* decrease the degree of the neighbors with greater degree */
			vertex v = graph[u][j];
			int curval = na_bs.CurrentValue(v);
			if (curval > degree_of_u)
				na_bs.DecVal(v);
			else if (curval != degree_of_u) { // **
				if (K[v] == value_of_u) { // @      // if v had got same K value with u
					if (compid[u] == -1) {
						compid[u] = compid[v];
						ecr.erase (ecr.begin() + ecr.size() - 1);
					}
					else {
						// merge compid[u] and compid[v] nodes
						int ch = compid[u];
						int pr = compid[v];
						find_same_K_parent1 (&ch, ecr);
						find_same_K_parent1 (&pr, ecr);
						if (ch != pr) {
							if (ecr[ch].rank > ecr[pr].rank) {
								int t = pr;
								pr = ch;
								ch = t;
							}
							ecr[ch].parent = pr;
							ecr[ch].visible = false;
							if (ecr[pr].rank == ecr[ch].rank)
								ecr[pr].rank++;
							num_subcores--;
						}
					}
				}
				else if (K[v] < value_of_u) { // @@
					couple1 c = make_tuple(compid[v], compid[u]);
					if (compid[u] == -1) { // it is possible that u didn't get an id yet
						minusones.push_back (relations.size()); // keep those indices to process after the loop, below
						relations.push_back (c); // ecr[compid[v]].K < ecr[compid[u]].K
					}
					else
						relations.push_back (c);
				}
			}
		}

		if (compid[u] == -1) // if u didn't get a compid
			compid[u] = id++; // give him a new one

		// update the minusones and put into relations, unless already there.
		// leave as it is if they are there, remember to skip those in post-processing
		for (int i : minusones) {
			couple1 c = make_tuple(get<0>(relations[i]), compid[u]);
			relations[i] = c;
		}
	}

	na_bs.Free();
	*core_number = degree_of_u; // degree_of_u is degree of last popped vertex
	util::timestamp t33;
	print_time_new (fp, "peeling time: ", t33 - t1);
	build (*core_number, relations, ecr, &num_subcores);
	printf ("number of subcores: %d   number of subsubcores: %d   graph.size: %d\n", num_subcores, ecr.size(), graph.size());
	report_subgraphs (12, ecr, compid, graph, nEdge, ax, vfile, fp);
	return;
}

void old_find_subcore (vertex root, Graph& graph, vector<vertex>& K, HashMap<bool>& visited, FILE* ffp) {
	int cnt = 0;
	vertex cr = K[root];
	queue<int> bfsorder;
	bfsorder.push(root);
	visited[root] = true;

	while (!bfsorder.empty()) {
		int v = bfsorder.front();
		bfsorder.pop();
		cnt++;
		for (size_t i = 0; i < graph[v].size(); i++) {
			vertex w = graph[v][i];
			// discovering an unvisited vertex with same K value: put it to subcore
			if (K[w] == cr && visited.hasDefaultValue(w)) {
				bfsorder.push(w);
				visited[w] = true;
			}
		}
	}

	fprintf (ffp, "K: %d size: %d\n", cr, cnt);
}
