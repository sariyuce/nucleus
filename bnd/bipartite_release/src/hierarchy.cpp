#include "main.h"

inline void assignToRoot (vertex* ch, vector<subcore>& skeleton) {
	vector<vertex> acc;
	vertex s = *ch;
	while (skeleton[s].root != -1) {
		acc.push_back (s);
		s = skeleton[s].root;
	}
	for (vertex i : acc)
		skeleton[i].root = s;
	*ch = s;
}

void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex rightnVtx, vertex leftnVtx) {

	// bin the relations w.r.t. first's K
	vector<vector<vp>> binnedRelations (cn + 1);

	for (int i = 0; i < relations.size(); i++) {
		vertex a = relations[i].first;
		vertex b = relations[i].second;

		assignToRepresentative (&a, skeleton);
		assignToRepresentative (&b, skeleton);

		if (a == b)
			continue;

		vp c (a, b);
		binnedRelations[skeleton[a].K].push_back (c);
	}

	// process binnedRelations in reverse order
	for (int i = binnedRelations.size() - 1; i >= 0; i--) {
		vector<vp> mergeList;
		for (vertex j = 0; j < binnedRelations[i].size(); j++) { // each binnedRelations[i] has K of skeleton[b].K
			vertex a = binnedRelations[i][j].first;
			vertex root = binnedRelations[i][j].second;
			assignToRoot (&root, skeleton);
			if (a != root) {
				if (skeleton[a].K < skeleton[root].K) {
					skeleton[root].parent = a;
					skeleton[root].root = a;
				}
				else { // skeleton[root].K == skeleton[a].K
					vp c = (root < a) ? make_pair (root, a) : make_pair (a, root);
					mergeList.push_back (c);
				}
			}
		}

		// handle merges
		for (auto sc : mergeList) {
			vertex child = sc.first;
			vertex parent = sc.second;
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
	vertex nid = skeleton.size();
	subcore sc (0);
	for (size_t i = 0; i < skeleton.size(); i++)
		if (skeleton[i].parent == -1)
			skeleton[i].parent = nid;

	sc.primarySize = rightnVtx;
	sc.secondarySize = leftnVtx;
	sc.nEdge = nEdge;
	sc.ed = nEdge / double (rightnVtx * leftnVtx);
	skeleton.push_back (sc);
}
