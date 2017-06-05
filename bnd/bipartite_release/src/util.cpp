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

inline bool pullChildrenSets (FILE* fp, vector<vertex>& children, HashMap<vertex>& orderInFile, vector<vertex>& vset, vector<subcore>& skeleton, Graph& rightGraph, bool edgeList, vector<vertex>* xRight) {

	char c;
	for (vertex eda : children) {

		if (skeleton[eda].nEdge == -1)
			return false;

		if (orderInFile.hasDefaultValue (eda)) {
			printf ("PROBLEM: %d has -1 as order\n", eda);
			exit(1);
		}

		vertex sc = orderInFile[eda];
		fseek (fp, 0, SEEK_SET);
		vertex ln = 0;
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
		double f;
		// check reportSubgraph() below to see what's what
		fscanf (fp, "%d %d %d %d %d %lf %d %d", &u, &u, &u, &u, &u, &f, &u, &u); // id, K, |LV|, |RV|, |E|, ed, LEAF?, parent-id
		while (fscanf (fp, "%d", &u) != EOF) {
			if (u != -1) {
				d = u;
				if (edgeList) {
					fscanf (fp, "%d", &v);
					d = (*xRight)[u] + find_ind (rightGraph[u], v);
				}
				vset.push_back (d);
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
	fprintf(fp, "%d %d %d %d %lf %d ", ind, sc->K, -1, -1, -1, sc->children.empty()?1:0);
	fprintf(fp, "-1\n");
}

inline void removeChild (vertex i, vector<subcore>& backup) {
	if (backup[i].parent != -1) {
		vertex pr = backup[i].parent;
		for (vertex j = 0; j < backup[pr].children.size(); j++)
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

#define INDSIZE subgraphSize(variant,ind,skeleton)
#define PARENTSIZE subgraphSize(variant,skeleton[ch].parent,skeleton)
#define CHSIZE subgraphSize(variant,ch,skeleton)

void nestedCircles (vector<subcore>& skeleton, vertex ind, FILE* fp, string cfl, string variant) {
	double parent_color = color(skeleton[ind].ed);
	fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"index\": \"%d\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld",
			parent_color, cfl.c_str(), ind, INDSIZE, skeleton[ind].ed, skeleton[ind].K, INDSIZE);
	if (skeleton[ind].children.size() == 1) {
		fprintf(fp, ", \"children\": [\n");
		int ch = skeleton[ind].children[0];
		// ghost child
		fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ",
				parent_color, PARENTSIZE - CHSIZE);
		// real child
		nestedCircles (skeleton, ch, fp, cfl, variant);
		fprintf(fp, "\n]\n");
	}
	else if (skeleton[ind].children.size() > 1) {
		fprintf(fp, ", \n\"children\": [\n");
		size_t i;
		for (i = 0; i < skeleton[ind].children.size() - 1; i++) {
			nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl, variant);
			fprintf(fp, ",\n");
		}
		nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl, variant);
		fprintf(fp, "\n]");
	}
	fprintf(fp, "}\n");
}

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
void reportSubgraph (string variant, vertex index, HashMap<vertex>& orderInFile, vector<vertex>& component,
		helpers& ax, vector<subcore>& skeleton, Graph& leftGraph, Graph& rightGraph, edge nEdge, FILE* fp, vector<vertex>* xRight = NULL) {

	if (skeleton[index].parent == -1) {
		skeleton[index].primarySize = rightGraph.size();
		skeleton[index].secondarySize = leftGraph.size();
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		fprintf(fp, "%d %d %d %d %d %lf %d %d\t", index, skeleton[index].K, skeleton[index].primarySize, skeleton[index].secondarySize, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0, skeleton[index].parent);
		fprintf(fp, "-1\n");
		return;
	}

	vector<vertex> vset;
	for (vertex i = 0; i < component.size(); i++)
		if (component[i] == index)
			vset.push_back (i);
	bool pass = true;

	pass = pullChildrenSets (fp, skeleton[index].children, orderInFile, vset, skeleton, rightGraph, variant == "WING", xRight);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	pass = hashUniquify (vset);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	if (variant == "TIP")
		tipDensity (vset, rightGraph, &(skeleton[index]));
	else if (variant == "WING")
		wingDensity (vset, *(ax.el), &(skeleton[index]));

	fprintf(fp, "%d %d %d %d %d %.2lf %d %d\t", index, skeleton[index].K, skeleton[index].primarySize, skeleton[index].secondarySize, skeleton[index].nEdge,
					skeleton[index].ed, skeleton[index].children.empty()?1:0, skeleton[index].parent);

	if (variant == "WING")
		for (size_t i = 0; i < vset.size(); i++)
			fprintf(fp, "%d %d  ", get<0>((*ax.el)[vset[i]]), get<1>((*ax.el)[vset[i]]));
	else
		for (size_t i = 0; i < vset.size(); i++)
			fprintf(fp, "%d ", vset[i]);

	fprintf(fp, "-1\n");
}

void bfsHierarchy (vector<subcore>& skeleton, stack<vertex>& scs) {

	rearrange (skeleton);
	queue<vertex> bfsorder; // we are doing bfs on the hierarchy tree and push the dequeued nodes to the stack
	bfsorder.push(skeleton.size() - 1);
	while (!bfsorder.empty()) {
		vertex s = bfsorder.front();
		bfsorder.pop();
		scs.push (s);
		for (vertex r : skeleton[s].children)
			bfsorder.push (r);
	}
}

void presentNuclei (string variant, vector<subcore>& skeleton, vector<vertex>& component, edge nEdge, helpers& ax, string vfile, FILE* gp, Graph& leftGraph, Graph& rightGraph, vector<vertex>* xRight) {

	// assign unassigned items to top subcore
	for (vertex i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	HashMap<vertex> replace (-1);
	for (vertex i = 0; i < skeleton.size(); i++) {
		vertex sc = i;
		vertex original = sc;
		findRepresentative (&sc, skeleton);
		if (original != sc)
			skeleton[original].visible = false;
		replace[original] = sc;
	}

	// each component takes its representative's component number
	for (vertex i = 0; i < component.size(); i++)
		if (!replace.hasDefaultValue (component[i]))
			component[i] = replace[component[i]];

	stack<vertex> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	string nFile = vfile + "_NUCLEI";
	FILE* fp = fopen (nFile.c_str(), "w+");
	vector<subcore> backup (skeleton);
	HashMap<vertex> orderInFile (-1); // key is the skeleton index, value is the order
	vertex o = 0; // order of subcores in file

	while (!subcoreStack.empty()) {
		vertex i = subcoreStack.top();
		subcoreStack.pop();
		if (backup[i].visible) {
			orderInFile[i] = o++;
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, leftGraph, rightGraph, nEdge, fp, xRight);
			removeChild (i, backup);
		}
	}
	fclose (fp);

	for (size_t i = 0; i < skeleton.size(); i++)
		if (skeleton[i].visible) {
			if (skeleton[i].primarySize >=5 && skeleton[i].secondarySize >= 5 && skeleton[i].ed >= 0.8 && skeleton[i].children.empty())
				printf ("%d: K: %d |PV|: %d |SV|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].primarySize, skeleton[i].secondarySize, skeleton[i].nEdge, skeleton[i].ed, skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
			fprintf (gp, "%d: K: %d |PV|: %d |SV|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].primarySize, skeleton[i].secondarySize, skeleton[i].nEdge, skeleton[i].ed, skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
			fflush (gp);
		}

	// build the json file for visualization -- based on the interactive nested circles
	string cfl = vfile + "_circle.json";
	FILE* ncp = fopen (cfl.c_str(), "w");
	nestedCircles (skeleton, skeleton.size() - 1, ncp, cfl, variant);
	fclose (ncp);
}
