#include "main.h"

// finding the nuclei subgraphs and their densities

inline void tipDensity(vector<vertex>& vset, Graph& rightGraph, subcore* sc) {
	vector<vertex> allNeighbors;
	edge nEdgesSubgraph = 0;
	neighborsOfNeighbors(vset, rightGraph, allNeighbors, &nEdgesSubgraph);
	sc->nEdge = nEdgesSubgraph;
	sc->ed = (double) nEdgesSubgraph / (vset.size() * allNeighbors.size()); // edge density
	sc->primarySize = vset.size();
	sc->secondarySize = allNeighbors.size();
}

inline void wingDensity(vector<vertex>& eset, vector<vp>& el, subcore* sc) {

	vector<vertex> lefties, righties;
	for (vertex i = 0; i < eset.size(); i++) {
		vertex r = el[eset[i]].first;
		vertex l = el[eset[i]].second;
		righties.push_back(r);
		lefties.push_back(l);
	}

	hashUniquify (lefties, false);
	hashUniquify (righties, false);

	sc->nEdge = eset.size();
	sc->ed = (double) eset.size() / (lefties.size() * righties.size()); // edge density
	sc->primarySize = righties.size();
	sc->secondarySize = lefties.size();
}

inline bool pullChildrenSets (string variant, FILE* fp, vector<ll>& children, HashMap<ll>& orderInFile, vector<vertex>& vset, vector<subcore>& skeleton, Graph& rightGraph, bool edgeList, vector<vertex>* xRight) {

	int limit = INT_MAX;
	if (variant == "TIP" || variant == "CORE" || variant == "TRUSS")
		limit = VERTEXUPPERBOUND;
	else if (variant == "WING")
		limit = EDGEUPPERBOUND;
	char c;
	for (auto eda : children) {

		if (skeleton[eda].nEdge == -1)
			return false;

		if (orderInFile.hasDefaultValue (eda)) {
			printf ("PROBLEM: %d has -1 as order\n", eda);
			exit(1);
		}

		auto sc = orderInFile[eda];
		fseek (fp, 0, SEEK_SET);
		ll ln = 0;
		if (sc != 0) {
			do {
				c = fgetc (fp);
				if (c == '\n') {
					ln++;
					if (ln == sc)
						break;
				}
			} while (c != EOF);
		}

		// now you are at the correct line of the file and can get the vertices of eda
		int u, v, d;
		ll uu;
		double f;
		// check reportSubgraph() below to see what's what
		fscanf (fp, "%d %d %d %d %d %lf %d %d", &uu, &u, &u, &u, &u, &f, &u, &uu); // id, K, |LV|, |RV|, |E|, ed, LEAF?, parent-id
		while (fscanf (fp, "%d", &u) != EOF) {
			if (u != -1) {
				d = u;
				if (edgeList) {
					fscanf (fp, "%d", &v);
					d = (*xRight)[u] + find_ind (rightGraph[u], v);
				}
				vset.push_back (d);
				if (vset.size() > limit)
					return false;
			}
			else
				break;
		}
		fseek (fp, 0, SEEK_END);
	}
	return true;
}

inline void dummyLine (subcore* sc, FILE* fp, vertex ind) {
	sc->primarySize = -1;
	sc->secondarySize = -1;
	sc->nEdge = -1;
	sc->ed = -1;
	sc->parent = -19;
	fprintf(fp, "%lld %d %d %d %d %lf %d %lld ", ind, sc->K, sc->primarySize, sc->secondarySize, sc->nEdge, sc->ed,
			sc->children.empty()?1:0, sc->parent);

	fprintf(fp, "-1\n");
}

inline void removeChild (ll i, vector<subcore>& backup) {
	if (backup[i].parent != -1) {
		ll pr = backup[i].parent;
		for (auto j = 0; j < backup[pr].children.size(); j++)
			if (backup[pr].children[j] == i) {
				backup[pr].children.erase (backup[pr].children.begin() + j);
				break;
			}
	}
}

double color (double ed) {
	double a = (1 - ed) / 0.25;
	int X = floor(a);
	int Y = floor(255 * (a-X));
	if (X < 4)
		return (3 - X) + (255 - (double) Y) / 255;
	else
		return 0;
}

inline vertex subgraphSize (string variant, vertex i, vector<subcore>& skeleton) {
	if (variant == "WING")
		return skeleton[i].nEdge;
	else
		return skeleton[i].primarySize;
}

inline void rearrange (vector<subcore>& skeleton) { // rearrange children and parents based on visibility

	for (size_t i = 0; i < skeleton.size(); i++)
		skeleton[i].children.clear();

	for (size_t i = 0; i < skeleton.size() - 1; i++) {
		if (skeleton[i].visible) {
			int pr = skeleton[i].parent;
			while (!skeleton[pr].visible)
				pr = skeleton[pr].parent;
			skeleton[i].parent = pr;
			skeleton[pr].children.push_back (i);
		}
	}
}

// find the r-cliques whose component is ind, append the vertices in those cliques to the vertices of its all children, sort, compute its density
void reportSubgraph (string variant, ll index, HashMap<ll>& orderInFile, vector<ll>& component,
		helpers& ax, vector<subcore>& skeleton, Graph& leftGraph, Graph& rightGraph, edge nEdge, FILE* fp, vector<vertex>* xRight = NULL) {

	if (skeleton[index].parent == -1) {
		skeleton[index].primarySize = rightGraph.size();
		skeleton[index].secondarySize = leftGraph.size();
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		fprintf(fp, "%lld %d %d %d %d %lf %d %lld\t", index, skeleton[index].K, skeleton[index].primarySize, skeleton[index].secondarySize, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0, skeleton[index].parent);
		fprintf(fp, "-1\n");
		return;
	}

	vector<vertex> vset;
	for (auto i = 0; i < component.size(); i++)
		if (component[i] == index)
			vset.push_back (i);
	bool pass = true;
	bool edgeList = (variant == "WING");
	pass = pullChildrenSets (variant, fp, skeleton[index].children, orderInFile, vset, skeleton, rightGraph, edgeList, xRight);

	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	pass = hashUniquify (vset, true, variant);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	if (variant == "TIP")
		tipDensity (vset, rightGraph, &(skeleton[index]));
	else if (variant == "WING")
		wingDensity (vset, *(ax.el), &(skeleton[index]));
	else {
		skeleton[index].nEdge = -1;
		skeleton[index].ed = -1;
		skeleton[index].primarySize = vset.size();
		skeleton[index].secondarySize = -1;
	}

	fprintf(fp, "%lld %d %d %d %d %.2lf %d %lld\t", index, skeleton[index].K, skeleton[index].primarySize, skeleton[index].secondarySize, skeleton[index].nEdge,
					skeleton[index].ed,
					skeleton[index].children.empty()?1:0,
							skeleton[index].parent);

	if (variant == "WING" || variant == "TRUSS")
		for (size_t i = 0; i < vset.size(); i++)
			fprintf(fp, "%d %d  ", get<0>((*ax.el)[vset[i]]), get<1>((*ax.el)[vset[i]]));
	else
		for (size_t i = 0; i < vset.size(); i++)
			fprintf(fp, "%d ", vset[i]);

	fprintf(fp, "-1\n");

}

inline void bfsHierarchy (vector<subcore>& skeleton, stack<ll>& scs) {
	rearrange (skeleton);
	queue<ll> bfsorder; // we are doing bfs on the hierarchy tree and push the dequeued nodes to the stack
	bfsorder.push(skeleton.size() - 1);
	while (!bfsorder.empty()) {
		ll s = bfsorder.front();
		bfsorder.pop();
		scs.push (s);
		for (auto r : skeleton[s].children)
			bfsorder.push (r);
	}
}

void presentNuclei (string variant, vector<subcore>& skeleton, vector<ll>& component, edge nEdge, helpers& ax, string vfile, FILE* gp, Graph& leftGraph, Graph& rightGraph, vector<vertex>* xRight) {

	// assign unassigned items to top subcore
	for (auto i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	HashMap<ll> replace (-1);
	for (auto i = 0; i < skeleton.size(); i++) {
		auto sc = i;
		auto original = sc;
		findRepresentative (&sc, skeleton);
		if (original != sc)
			skeleton[original].visible = false;
		replace[original] = sc;
	}

	// each component takes its representative's component number
	for (auto i = 0; i < component.size(); i++)
		if (!replace.hasDefaultValue (component[i]))
			component[i] = replace[component[i]];

	stack<ll> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	string nFile = vfile + "_NUCLEI";
	FILE* fp = fopen (nFile.c_str(), "w+");
	vector<subcore> backup (skeleton);
	HashMap<ll> orderInFile (-1); // key is the skeleton index, value is the order
	ll o = 0; // order of subcores in file

	while (!subcoreStack.empty()) {
		ll i = subcoreStack.top();
		subcoreStack.pop();
		if (backup[i].visible) {
			orderInFile[i] = o++;
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, leftGraph, rightGraph, nEdge, fp, xRight);
			removeChild (i, backup);
		}
	}
	fclose (fp);
}





// building the projections

void unweighted_projection (Graph& left, Graph& right, string filename) {

	Graph graph;
	graph.resize (left.size());
	int nedge = 0;
	for (size_t i = 0; i < left.size(); i++) {
		HashMap<bool> mp (false);
		for (int j = 0; j < left[i].size(); j++) {
			int v = left[i][j];
			for (auto w : right[v])
				mp[w] = true;
		}
		for (auto it = mp.begin(); it != mp.end(); it++)
			if (it->first != i)
				graph[i].push_back (it->first);

		sort (graph[i].begin(), graph[i].end());
		nedge += graph[i].size();
	}

	for (int i = 0; i < graph.size(); i++)
		for (int j = 0; j < graph[i].size(); j++)
			if (graph[i][j] == i) {
				printf ("self loop   i: %d\n", i);
				exit(1);
			}
	filename += "_UW_projection";

	writeRegularBinary (filename, graph.size(), nedge, graph);
}

inline void add (vector<wv>& node, double weight, int u) {
	for (int i = 0; i < node.size(); i++) {
		if (node[i].n == u) {
			node[i].w += weight;
			return;
		}
	}
	wv t;
	t.n = u;
	t.w = weight;
	node.push_back (t);
	return;
}

bool wgc (wv i, wv j) { return (i.n < j.n); }

void weighted_projection (Graph& left, Graph& right, string filename) {
	Wraph wg;
	wg.resize (left.size());
	for (int i = 0; i < right.size(); i++) {
		if (right[i].size() <= 1)
			continue;

		double weight = 1.0 / right[i].size();
		for (int j = 0; j < right[i].size(); j++) {
			for (int k = j + 1; k < right[i].size(); k++) {
				int u = right[i][j];
				int v = right[i][k];
				add (wg[u], weight, v);
				add (wg[v], weight, u);
			}
		}
	}
	edge nedge = 0;
	for (int i = 0; i < wg.size(); i++) {
		sort (wg[i].begin(), wg[i].end(), wgc);
		nedge += wg[i].size();
	}
	string fn = filename + "_W_projection";
    writeWeightedRegularBinary (fn, wg.size(), nedge, wg);
}

