#include "main.h"

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

//#define INDSIZE subgraphSize(variant,ind,skeleton)
//#define PARENTSIZE subgraphSize(variant,skeleton[ch].parent,skeleton)
//#define CHSIZE subgraphSize(variant,ch,skeleton)
//
//void nestedCircles (vector<subcore>& skeleton, vertex ind, FILE* fp, string cfl, string variant) {
//	double parent_color = color(skeleton[ind].ed);
//	fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"index\": \"%d\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld",
//			parent_color, cfl.c_str(), ind, INDSIZE, skeleton[ind].ed, skeleton[ind].K, INDSIZE);
//	if (skeleton[ind].children.size() == 1) {
//		fprintf(fp, ", \"children\": [\n");
//		int ch = skeleton[ind].children[0];
//		// ghost child
//		fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ",
//				parent_color, PARENTSIZE - CHSIZE);
//		// real child
//		nestedCircles (skeleton, ch, fp, cfl, variant);
//		fprintf(fp, "\n]\n");
//	}
//	else if (skeleton[ind].children.size() > 1) {
//		fprintf(fp, ", \n\"children\": [\n");
//		size_t i;
//		for (i = 0; i < skeleton[ind].children.size() - 1; i++) {
//			nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl, variant);
//			fprintf(fp, ",\n");
//		}
//		nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl, variant);
//		fprintf(fp, "\n]");
//	}
//	fprintf(fp, "}\n");
//}

void rearrange (vector<subcore>& skeleton) { // rearrange children and parents based on visibility

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

void bfsHierarchy (vector<subcore>& skeleton, stack<ll>& scs) {

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

//	{
//		string cFile = vfile + "_COMPS";
//		FILE* ap = fopen (cFile.c_str(), "w");
//		for (vertex i = 0; i < component.size(); i++) {
//			fprintf (af, "%d\n", component[i]);
//		}
//		fclose (af);
//		printf("done");
//		exit(1);
//	}

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
//		if (backup[i].visible) {
//			if (backup[i].children.empty() || backup[i].children.size() > 1)
//		}
		if (backup[i].visible) {
//			printf ("i: %d   parent: %d\n", i, skeleton[i].parent);
			orderInFile[i] = o++;
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, leftGraph, rightGraph, nEdge, fp, xRight);
			removeChild (i, backup);
		}
	}
	fclose (fp);

//	for (size_t i = 0; i < skeleton.size(); i++)
//		if (skeleton[i].visible) {
//			if (skeleton[i].primarySize >=5 && skeleton[i].secondarySize >= 5 && skeleton[i].ed >= 0.8 && skeleton[i].children.empty())
//				printf ("%d: K: %d |PV|: %d |SV|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].primarySize, skeleton[i].secondarySize, skeleton[i].nEdge, skeleton[i].ed, skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
//			fprintf (gp, "%d: K: %d |PV|: %d |SV|: %d |E|: %d ed: %.5lf  #children: %d   pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].primarySize, skeleton[i].secondarySize, skeleton[i].nEdge, skeleton[i].ed, skeleton[i].children.size(), skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
//			fflush (gp);
//		}







//	// build the json file for visualization -- based on the interactive nested circles
//	string cfl = vfile + "_circle.json";
//	FILE* ncp = fopen (cfl.c_str(), "w");
//	nestedCircles (skeleton, skeleton.size() - 1, ncp, cfl, variant);
//	fclose (ncp);
}










// Projection building codes

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

//	writeToRegularBinary (filename, graph.size(), nedge, graph);
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

	string fan = "/scratch1/scratchdirs/aerdem/bipartite/newProjs/NEW_" + filename + "_W_projection";
	FILE* fff = fopen (fan.c_str(), "w");
	for (int i = 0; i < right.size(); i++) {
		if (right[i].size() <= 1)
			continue;

		double weight = 1.0 / right[i].size();
		for (int j = 0; j < right[i].size(); j++) {
			for (int k = j + 1; k < right[i].size(); k++) {
				int u = right[i][j];
				int v = right[i][k];
				fprintf (fff, "%d %d %lf\n", u<v?u:v, u<v?v:u, weight);
				//					add (wg[u], weight, v);
				//					add (wg[v], weight, u);
			}
		}
	}
	fclose (fff);
	return;






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
//		sort (gr[i].begin(), gr[i].end());
		nedge += wg[i].size();
	}


//	Graph gr;
//	gr.resize (left.size());

    string fn = "/scratch1/scratchdirs/aerdem/bipartite/newProjs/" + filename + "_W_projection";
    writeToWeightedRegularBinary (fn, wg.size(), nedge, wg);

//    string fasdf = "/scratch1/scratchdirs/aerdem/bipartite/newProjs/" + filename + "_UNW_projection";
//    writeToRegularBinary (fasdf, gr.size(), nedge, gr);
//    FILE* fp = fopen (fn.c_str(), "w");
//    fprintf (fp, "%d %d\n", wg.size(), nedge);
//    for (int i = 0; i < wg.size(); i++) {
//        for (int j = 0; j < wg[i].size(); j++)
//            fprintf (fp, "%d %.3lf ", wg[i][j].n, wg[i][j].w);
//        fprintf (fp, "\n");
//    }
//    fclose (fp);
}








//void weighted_base_kcore (Wraph& wraph, int nEdge, vector<vertex>& K, bool hierarchy, int* maxCore, const char* vfile, FILE* fp) {
//
//	util::timestamp peelingStart;
//	size_t nVtx = wraph.size();
//	double max_degree = 0;
//
//
////	for (size_t i = 0; i < nVtx; i++) {
////		printf ("i: %d   ", i);
////		for (int j = 0; j < wraph[i].size(); j++) {
////			printf ("%d %lf   ", wraph[i][j].n, wraph[i][j].w);
////		}
////		printf ("\n");
////	}
//
//
//	for (size_t i = 0; i < nVtx; i++) {
//		double wsum = 0;
//		for (int j = 0; j < wraph[i].size(); j++)
//			wsum += wraph[i][j].w;
//		if (wsum > max_degree)
//			max_degree = wsum;
//	}
//	max_degree *= 10;
//	int maxdeg = (int) max_degree + 1;
//
//	// Peeling
//	K.resize (wraph.size(), -1);
//	Naive_Bucket na_bs;
//	na_bs.Initialize(maxdeg, wraph.size());
//
//	for (size_t i = 0; i < wraph.size(); i++) {
//		int wsum = 0;
//		for (int j = 0; j < wraph[i].size(); j++)
//			wsum += (int) (10 * wraph[i][j].w);
////		if (wsum > 0)
//			na_bs.Insert (i, wsum);
////		else
////			K[i] = 0;
//	}
//
//	vertex degree_of_u = 0;
//	vertex cid;
//	vector<subcore> skeleton;
//	vector<vertex> component;
//	vector<vp> relations;
//	vector<vertex> unassigned;
//	vertex nSubcores;
//
//	if (hierarchy) {
//		cid = 0;
//		nSubcores = 0;
//		component.resize (wraph.size(), -1);
//	}
//
//	while (1) {
//		int u, value_of_u;
//		int ret = na_bs.PopMin(&u, &value_of_u);
//		if (ret == -1) // if the bucket is empty
//			break;
//
////		printf ("u: %d  val: %d\n", u, value_of_u);
//		if (value_of_u == 0) {
//			K[u] = 0;
//			continue;
//		}
//
//		if (hierarchy) {
//			unassigned.clear();
//			subcore sc (value_of_u);
////			printf ("skeleton.size: %d and extended\n", skeleton.size());
//			skeleton.push_back (sc);
//		}
//
////		printf ("in cid: %d\n", cid);
//		degree_of_u = K[u] = value_of_u;
//
//		for (size_t j = 0; j < wraph[u].size(); j++) { /* decrease the degree of the neighbors with greater degree */
//			vertex v = wraph[u][j].n;
//			int weig = (int) (wraph[u][j].w * 10);
//			int curval = na_bs.CurrentValue(v);
//			if ((K[v] == -1) && (curval - weig >= degree_of_u)) {
//				for (int k = 0; k < weig; k++)
//					na_bs.DecVal(v);
//			}
//			else if ((K[v] == -1) && (curval - weig < degree_of_u))
//				;
//			else if (hierarchy)
//				createSkeleton (u, {v}, &nSubcores, K, skeleton, component, unassigned, relations);
//		}
//
//
//		if (hierarchy)
//			updateUnassigned (u, component, &cid, relations, unassigned);
//
//	}
//
//	na_bs.Free();
//	*maxCore = degree_of_u; // degree_of_u is degree of last popped vertex
//	util::timestamp t33;
//	Graph graph (wraph.size());
//	for (int i = 0; i < wraph.size(); i++) {
//		for (int j = 0; j < wraph[i].size(); j++)
//			graph[i].push_back (wraph[i][j].n);
//	}
//
////	printf ("cid: %d\n", cid);
//	timestamp peelingEnd;
//	cout << "Total time: " << peelingEnd - peelingStart << endl;
//	print_time (fp, "Total time: ", peelingEnd - peelingStart);
//
//	if (hierarchy) {
//		buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, graph.size());
//		timestamp nucleusEnd;
//
//		print_time (fp, "Fractional k-core decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
//		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());
//
//		helpers dummy;
//		Graph dum;
//		presentNuclei ("CORE", skeleton, component, nEdge, dummy, vfile, fp, graph, dum, NULL);
//		timestamp totalEnd;
//
//		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
//	}
//
//	return;
//}
//
//


//
//void base_kcore (Graph& graph, int nEdge, vector<vertex>& K, bool hierarchy, int* maxCore, string vfile, FILE* fp) {
//
//	timestamp peelingStart;
//
//	vertex nVtx = graph.size();
//	vertex maxDeg = 0;
//	for (auto g : graph) {
//		if (g.size() > maxDeg)
//			maxDeg = g.size();
//	}
//
//	// Peeling
//	K.resize (graph.size(), -1);
//	Naive_Bucket nBucket;
//	nBucket.Initialize(maxDeg, nVtx);
//	for (vertex i = 0; i < graph.size(); i++) {
//		if (graph[i].size() > 0)
//			nBucket.Insert (i, graph[i].size());
//		else
//			K[i] = 0;
//	}
//
//	vertex deg_u = 0;
//
//	vertex cid;
//	vector<subcore> skeleton;
//	vector<vertex> component;
//	vector<vp> relations;
//	vector<vertex> unassigned;
//	vertex nSubcores;
//
//	if (hierarchy) {
//		cid = 0;
//		nSubcores = 0;
//		component.resize (graph.size(), -1);
//	}
//
//	while (true) {
//		vertex u, val;
//		if (nBucket.PopMin(&u, &val) == -1) // if the bucket is empty
//			break;
//
//		if (hierarchy) {
//			unassigned.clear();
//			subcore sc (val);
//			skeleton.push_back (sc);
//		}
//
//		deg_u = K[u] = val;
//
//		for (auto v : graph[u]) {
//			vertex deg_v = nBucket.CurrentValue(v);
//			if (deg_v > deg_u)
//				nBucket.DecVal(v);
//			else if (hierarchy)
//				createSkeleton (u, {v}, &nSubcores, K, skeleton, component, unassigned, relations);
//		}
//		if (hierarchy)
//			updateUnassigned (u, component, &cid, relations, unassigned);
//	}
//
//	nBucket.Free();
//	*maxCore = deg_u;
//
//	printf ("maxCore: %d\n", *maxCore);
//
//	timestamp peelingEnd;
//	cout << "Total time: " << peelingEnd - peelingStart << endl;
//	print_time (fp, "Total time: ", peelingEnd - peelingStart);
//
//
//	if (hierarchy) {
//		buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, nVtx);
//		timestamp nucleusEnd;
//
//		print_time (fp, "Nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
//		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());
//
//		helpers dummy;
//		Graph dum;
//		presentNuclei ("CORE", skeleton, component, nEdge, dummy, vfile, fp, graph, dum, NULL);
//		timestamp totalEnd;
//
//		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
//	}
//}
//
//



// per edge
void countTriangles (Graph& graph, Graph& orientedGraph, Graph& TC) {
	for (size_t i = 0; i < orientedGraph.size(); i++) {
		for (size_t j = 0; j < orientedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orientedGraph[i].size(); k++) {
				vertex a = orientedGraph[i][j];
				vertex b = orientedGraph[i][k];
				if (incrementTCIfConnected (graph, orientedGraph, TC, a, b)) {
					TC[i][j]++;
					TC[i][k]++;
				}
			}
		}
	}
}

void base_ktruss (Graph& graph, edge nEdge, vector<vertex>& K, bool hierarchy, vertex* maxtruss, string vfile, FILE* fp) {

	timestamp countingStart;

	vertex nVtx = graph.size();

	// Prepare a CSR-like structure to index edges and create directed graph from low degree vertices to higher degree vertices
	vector<vp> el;
	vector<vertex> xel;
	Graph orientedGraph;
	indexEdges (graph, el, xel, orientedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orientedGraph[i].size(), 0);
	countTriangles (graph, orientedGraph, TC);

	timestamp peelingStart;
	cout << "Triangle counting per edge time: " << peelingStart - countingStart << endl;
	print_time (fp, "Triangle counting per edge time: ", peelingStart - countingStart);

	for (int i = 0; i < TC.size(); i++)
		for (int j = 0; j < TC[i].size(); j++)
			printf ("TC: %d\n", TC[i][j]);

	exit(1);

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


	cout << "Bucket insertions done" << endl;
	vertex tc_e = 0;

	// required for hierarchy
	ll cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<ll> component; // subcore ids for each vertex
	vector<llp> relations;
	vector<ll> unassigned;
	ll nSubcores;

	if (hierarchy) {
		cid = 0;
		nSubcores = 0;
		component.resize (nEdge, -1);
	}

	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;

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
	cout << "Only peeling time: " << peelingEnd - peelingStart << endl;
	print_time (fp, "Only peeling time: ", peelingEnd - peelingStart);
	cout << "Total time: " << peelingEnd - countingStart << endl;
	print_time (fp, "Total time: ", peelingEnd - countingStart);

	cout << "maxtruss: " << *maxtruss << endl;

#ifdef K_VALUES
	for (int i = 0; i < K.size(); i++)
		printf ("K[%d]: %d\n", i, K[i]);
#endif


	if (hierarchy) {
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
		timestamp nucleusEnd;

		print_time (fp, "Nucleus decomposition time with hierarchy construction: ", nucleusEnd - peelingStart);
		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		helpers hp (&el);
		Graph dum;
		presentNuclei ("TRUSS", skeleton, component, nEdge, hp, vfile, fp, graph, dum, NULL);
		timestamp totalEnd;

		print_time (fp, "Total time, including the density computations: ", totalEnd - peelingStart);
	}
}






