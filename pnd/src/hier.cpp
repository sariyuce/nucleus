#include "main.h"

inline void find_root (int* ch, vector<subcore>& ecr) {
	vector<int> acc;
	int s = *ch;
	while (ecr[s].root != -1) {
		acc.push_back (s);
		s = ecr[s].root;
	}
	for (int i : acc)
		ecr[i].root = s;
	*ch = s;
}

inline void find_same_K_parent4 (int* ch, vector<subcore>& ecr) {
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

void build (int cn, vector<couple>& relations, vector<subcore>& ecr, int* num_subcores, int nEdge, int nVtx) {
	// regularization
	vector<vector<couple>> mm;
	for (int i = 0; i < cn+1; i++) {
		vector<couple> t;
		mm.push_back(t);
	}

	printf ("relations.size: %d\n", relations.size());
	for (int i = 0; i < relations.size(); i++) {

		int b = get<1>(relations[i]);
		if (b == -1) // not necessary todo
			continue;  // not necessary todo
		int a = get<0>(relations[i]);

		find_same_K_parent4 (&a, ecr);
		find_same_K_parent4 (&b, ecr);

		if (a == b)
			continue;

		couple c = make_tuple (a, b);
		ecr[a].visible = ecr[b].visible = true; // not necessary, todo
		mm[ecr[a].K].push_back (c);
	}

	// process mm in reverse order
	for (int i = mm.size() - 1; i >= 0; i--) {
		vector<couple> merge_sc;
		for (int j = 0; j < mm[i].size(); j++) { // all have K of ecr[b].K
			int a = get<0>(mm[i][j]);
			int b = get<1>(mm[i][j]);
			// find b's top parent. no need to find a's top parent.
			int root = b;
			find_root (&root, ecr);
			if (a != root) {
				if (ecr[a].K < ecr[root].K) {
					ecr[root].parent = a;
					ecr[root].root = a;
				}
				else { // ecr[root].K == ecr[a].K
					int e = root;
					int f = a;
					couple c = (e<f) ? make_tuple (e, f) : make_tuple (f, e);
						merge_sc.push_back (c);
				}
			}
		}

		// handle merges
		for (int ii = 0; ii < merge_sc.size(); ii++) {
			int ch = get<0>(merge_sc[ii]);
			int pr = get<1>(merge_sc[ii]);
			find_same_K_parent4 (&ch, ecr);
			find_same_K_parent4 (&pr, ecr);
			if (ch == pr)
				continue;
			if (ecr[ch].rank > ecr[pr].rank) {
				int t = pr;
				pr = ch;
				ch = t;
			}
			ecr[ch].parent = pr;
			ecr[ch].root = pr;
			ecr[ch].visible = false;
			if (ecr[pr].rank == ecr[ch].rank)
				ecr[pr].rank++;
		}
	}

	*num_subcores += ecr.size();

	// root core
	int nid = ecr.size();
	subcore sc;
	sc.rank = 0;
	sc.K = 0;
	sc.parent = -1;
	sc.root = -1;
	sc.visible = true;
	for (size_t i = 0; i < ecr.size(); i++) {
		if (ecr[i].visible && ecr[i].parent == -1) { // first condition is not necessary todo
			ecr[i].parent = nid;
		}
	}
#ifndef SHORT
	sc.nedge = nEdge;
	sc.ed = nEdge / double (nVtx * (nVtx - 1) / 2);
#endif
	ecr.push_back (sc);
}

inline bool exists (int val, vector<int>& v) {
	for (size_t i = 0; i < v.size(); i++)
		if (v[i] == val)
			return true;
	return false;
}

bool sortBy(const vector<int>& vec1, const vector<int>& vec2) {
	if (vec1[0] == vec2[0]) {
		if (vec1.size() != vec2.size())
			return (vec1.size() < vec2.size());
		else {
			for(size_t i = 0; i < vec1.size() && i < vec2.size(); i++) {
				 if (vec1[i] > vec2[i])
					 return false;
				 else if (vec1[i] < vec2[i])
					 return true;
			}
			return false;
		}
	}
	else
		return (vec1[0] < vec2[0]);
}

void traverse (int ind, vector<vector<int> >& matches, HashMap<bool>& visited, vector<int>& result) {
	queue<int> bfsorder;
	visited[ind] = true;
	bfsorder.push(ind);
	while (!bfsorder.empty()) {
		vertex v = bfsorder.front();
		bfsorder.pop();
		result.push_back(v);
		for (size_t i = 1; i < matches[v].size(); i++) {
			int val = matches[v][i];
			for (size_t k = 0; k < matches.size(); k++) {
				if (visited.hasDefaultValue (k) && (exists (val, matches[k]))) {
					bfsorder.push (k);
					visited[k] = true;
				}
			}
		}
	}
}

void smart_insert (mmap& ancestors, int key, int value) {
	auto range = ancestors.equal_range(key);
	for (auto it = range.first; it != range.second; it++) {
		if (it->second == value)
			return;
	}
	ancestors.emplace (key, value);
}

void find_my_family (int K, int sc, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct, util::timestamp& p1, util::timestamp& p2, util::timestamp& p3) {
	timestamp t1;
	HashMap<bool> seen (false);
	vector<int> merge_sc;
	auto range = kids.equal_range (sc);
	for (auto i_it = range.first; i_it != range.second; i_it++) {
		int u = i_it->second;
		if (seen.hasDefaultValue(u)) {
			timestamp e1;
			int init_u = u;
			while (!direct.hasDefaultValue(u)) {
				u = direct[u];
			}
			timestamp e2;
			int mid_u = u;
			while (!parent.hasDefaultValue(u))
				u = parent[u];
			timestamp e3;
			p2 += e2 - e1;
			p3 += e3 - e2;
			if (seen.hasDefaultValue(u)) {
				if (K_comps[u] > K)
					smart_insert (new_kids, i_it->first, u);
				else if (K_comps[u] == K)
					merge_sc.push_back (u);
				seen[u] = true;
			}
			seen[init_u] = true;
			seen[mid_u] = true;
		}
	}
	timestamp t2;
	p1 += t2 - t1;

	range = new_kids.equal_range (sc);
	for (auto it = range.first; it != range.second; it++)
		parent[it->second] = sc;

	int new_parent = sc;
	vector<int> toemplace;
	for (size_t i = 0; i < merge_sc.size(); i++) {
		int ns = merge_sc[i];
		direct[ns] = sc;
		auto iter = new_kids.find(ns);
		if (iter != new_kids.end()) {
			auto range = new_kids.equal_range (ns);
			if (range.first != range.second) {
				for (auto it = range.first; it != range.second; ++it) {
					parent[it->second] = new_parent; // assign new_parent to new kid
					toemplace.push_back (it->second);
				}
				new_kids.erase (ns);
			}
		}
	}

	for (size_t i = 0; i < toemplace.size(); i++)
		smart_insert (new_kids, new_parent, toemplace[i]);
}

void handle_matches (vector<vector<vector<int> > >& all_matches, timestamp& t02, timestamp& t03, HashMap<bool>& visited, HashMap<int>& direct) {

	for (size_t i = 0; i < all_matches.size(); i++)
		for (size_t j = 0; j < all_matches[i].size(); j++)
			sort (all_matches[i][j].begin() + 1, all_matches[i][j].end()); // exclude K value while sorting
	timestamp t2;
	t02 = t2;

	for (size_t i = 0; i < all_matches.size(); i++) {
		sort (all_matches[i].begin(), all_matches[i].end(), sortBy);

		for (size_t j = 1; j < all_matches[i].size(); j++) {
			if (all_matches[i][j] == all_matches[i][j-1]) {
				all_matches[i].erase (all_matches[i].begin() + j);
				j--;
			}
		}
	}
	timestamp t3;
	t03 = t3;
	vector<int> result;

	for (size_t n = 0; n < all_matches.size(); n++) {
		visited.reset (false);
		for (size_t i = 0; i < all_matches[n].size(); i++) {
			if (!visited.hasDefaultValue (i))
				continue;
			traverse (i, all_matches[n], visited, result);
			int minid = INT_MAX;
			for (size_t j = 0; j < result.size(); j++) {
				for (size_t k = 1; k < all_matches[n][result[j]].size(); k++) {
					int sc = all_matches[n][result[j]][k]; // sc is subcore!
					if (sc < minid)
						minid = sc;
				}
			}

			for (size_t j = 0; j < result.size(); j++) {
				for (size_t k = 1; k < all_matches[n][result[j]].size(); k++) {
					int sc = all_matches[n][result[j]][k]; // sc is subcore!
					if (direct.hasDefaultValue (sc)) { // checks only once assignment of new comp id
						if (sc != minid)
							direct[sc] = minid;
					}
				}
			}
			result.clear();
		}
	}
}

void inverse_anc2kid (int id, mmap& ancestors, mmap& kids, HashMap<bool>& visited, HashMap<int>& direct) {
	for (vertex sc = 0; sc < id; ++sc) {
		auto rn = ancestors.equal_range (sc);
		visited.reset (false);
		for (auto local_it = rn.first; local_it != rn.second; ++local_it) {
			int u = local_it->second;
			if (!direct.hasDefaultValue(u))
				u = direct[u];
			if (visited.hasDefaultValue (u)) {
				visited[u] = true;
				kids.emplace (u, sc);
			}
		}
	}
	ancestors.clear();
}

void feed_forward (HashMap<int>& direct) {
	for (auto it = direct.begin(); it != direct.end(); it++) {
		int u = it->second;
		while (!direct.hasDefaultValue (u))
			u = direct[u];
		it->second = u;
	}
}



void build_hierarchy (int id, int tn, mmap& cmps, HashMap<int>& count_cmps, mmap& new_kids, mmap& kids, HashMap<int>& K_comps, HashMap<int>& parent, HashMap<int>& direct, FILE* fp) {
	timestamp p1 (0, 0);
	timestamp p2 (0, 0);
	timestamp p3 (0, 0);
	timestamp p4 (0, 0);
	int in = 0;
	double prev_rt = 0;
	for (int i = tn; i >= 0; i--) {

		timestamp n1;

		auto range = cmps.equal_range (i);
		int y = 0;
		int ch = 2;
		int pr = 10000;
		int chunk = count_cmps[i] / ch;
		chunk = (chunk < pr)? pr : chunk;
		for (auto it = range.first; it != range.second; it++, y++) {
			find_my_family (it->first, it->second, new_kids, kids, K_comps, parent, direct, p1, p2, p3);
			in++;
			double rt = (double) in * 100 / cmps.size();
			if (rt - prev_rt >= 1) {
				prev_rt = rt;
				fprintf (fp, "%lf of build hierarchy (%d / %d)\n", rt, in, cmps.size());
				fflush (fp);
			}
		}
		feed_forward (direct);
		timestamp n2;
		cout << "K: " << i << " count: " << count_cmps[i] << " time " << n2 - n1 << " p1 so far " << p1 << endl;
	}
	kids.clear();
	cout << "p1 in build hier " << p1 << endl;
	cout << "p2 in build hier " << p2 << endl;
	cout << "p3 in build hier " << p3 << endl;
	cout << "p4 in build hier " << p4 << endl;
}



void find_subcore (vertex root, Graph& graph, vertex* K,
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

	for (size_t i = 0; i < vs.size(); i++)
		subcore_ids[vs[i]] = next_id;
	hrc[next_id].size = vs.size();
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

