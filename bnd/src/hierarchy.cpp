#include "main.h"

inline void assignToRoot (int* ch, vector<subcore>& skeleton) {
	vector<int> acc;
	int s = *ch;
	while (skeleton[s].root != -1) {
		acc.push_back (s);
		s = skeleton[s].root;
	}
	for (auto i : acc)
		skeleton[i].root = s;
	*ch = s;
}

void buildHierarchy (vertex cn, vector<llp>& relations, vector<subcore>& skeleton, int* nSubcores, edge nEdge,
		vertex rightnVtx, vertex leftnVtx) {

	// bin the relations w.r.t. first's K
	vector<vector<llp>> binnedRelations (cn + 1);

	for (int i = 0; i < relations.size(); i++) {
		int a = relations[i].first;
		int b = relations[i].second;
		assignToRepresentative (&a, skeleton);
		assignToRepresentative (&b, skeleton);
		if (a == b)
			continue;
		llp c (a, b);
		binnedRelations[skeleton[a].K].push_back (c);
	}

	// process binnedRelations in reverse order
	for (int i = binnedRelations.size() - 1; i >= 0; i--) {
		vector<llp> mergeList;
		for (llp br : binnedRelations[i]) {
			int a = br.first;
			int root = br.second;
			assignToRoot (&root, skeleton);
			if (a != root) {
				if (skeleton[a].K < skeleton[root].K) {
					skeleton[root].parent = a;
					skeleton[root].root = a;
				}
				else { // skeleton[root].K == skeleton[a].K
					llp c = (root < a) ? make_pair (root, a) : make_pair (a, root);
					mergeList.push_back (c);
				}
			}
		}

		// handle merges
		for (auto sc : mergeList) {
			int child = sc.first;
			int parent = sc.second;
			assignToRepresentative (&child, skeleton);
			assignToRepresentative (&parent, skeleton);
			if (child != parent) {
				if (skeleton[child].rank > skeleton[parent].rank)
					swap (child, parent);
				skeleton[child].parent = parent;
				skeleton[child].root = parent;
				skeleton[child].visible = false;
				if (skeleton[parent].rank == skeleton[child].rank)
					skeleton[parent].rank++;
			}
		}
	}

	*nSubcores += skeleton.size();

	// root core
	int nid = skeleton.size();
	subcore sc (0);
	for (auto i = 0; i < skeleton.size(); i++)
		if (skeleton[i].parent == -1)
			skeleton[i].parent = nid;

	sc.primarySize = rightnVtx;
	sc.secondarySize = leftnVtx;
	sc.nEdge = nEdge;
	sc.ed = nEdge / double (rightnVtx * leftnVtx);
	skeleton.push_back (sc);
}
