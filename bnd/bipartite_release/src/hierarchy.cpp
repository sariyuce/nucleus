#include "main.h"

inline void assignToRoot (ll* ch, vector<subcore>& skeleton) {
	vector<ll> acc;
	ll s = *ch;
	while (skeleton[s].root != -1) {
		acc.push_back (s);
		s = skeleton[s].root;
	}
	for (auto i : acc)
		skeleton[i].root = s;
	*ch = s;
}

void buildHierarchy (vertex cn, vector<llp>& relations, vector<subcore>& skeleton, ll* nSubcores, edge nEdge, vertex rightnVtx, vertex leftnVtx) {

	// bin the relations w.r.t. first's K
	vector<vector<llp>> binnedRelations (cn + 1);

	for (int i = 0; i < relations.size(); i++) {
		ll a = relations[i].first;
		ll b = relations[i].second;

//		printf ("a: %d sksz: %d    %d/%d\n", a, skeleton.size(), i , relations.size());
		assignToRepresentative (&a, skeleton);
//		printf ("b: %d\n", b);
		assignToRepresentative (&b, skeleton);

		if (a == b)
			continue;

		llp c (a, b);
		binnedRelations[skeleton[a].K].push_back (c);
//		printf ("push %lld %lld in to bR's %d\n", a, b, skeleton[a].K);
	}

	// process binnedRelations in reverse order
	for (int i = binnedRelations.size() - 1; i >= 0; i--) {
		vector<llp> mergeList;
//		for (auto j = 0; j < binnedRelations[i].size(); j++) { // each binnedRelations[i] has K of skeleton[b].K
//		printf ("process i: %d\n", i);
		for (llp br : binnedRelations[i]) {
			ll a = br.first; // binnedRelations[i][j].first;
			ll root = br.second; // binnedRelations[i][j].second;
			assignToRoot (&root, skeleton);
			if (a != root) {
//				printf ("a: %d, root: %d\n", a, root);
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
			ll child = sc.first;
			ll parent = sc.second;
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
	ll nid = skeleton.size();
//	printf ("sz before root: %d\n", skeleton.size());
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
