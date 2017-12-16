#include "main.h"

inline vertex commons (vector<vertex>& a, vector<vertex>& b) {
	vertex i = 0, j = 0;
	vertex count = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			count++;
			i++;
			j++;
		}
	}
	return count;
}

bool pullChildrenSets (FILE* fp, vector<vertex>& children, HashMap<vertex>& orderInFile, vector<vertex>& vset, vector<subcore>& skeleton) {

	int limit = UPPERBOUND;
	char c;
	for (vertex eda : children) {
		if (skeleton[eda].size == -1)
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
		int d;
		double f;
		// check line 163 to see what's what
		fscanf (fp, "%d %d %d %d %lf %d %d", &d, &d, &d, &d, &f, &d, &d); // id, K, |V|, |E|, ed, LEAF?
		while (fscanf (fp, "%d", &d) != EOF) {
			if (d != -1) {
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

inline void dummyLine (subcore* sc, FILE* fp, vertex index) {
	sc->size = -1;
	sc->ed = -1;
	sc->nEdge = -1;
	sc->parent = -1;
	fprintf(fp, "%d %d %d %d %lf %d %d ", index,	sc->K, sc->size, sc->nEdge, sc->ed, sc->children.empty()?1:0, sc->parent);
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

// find the r-cliques whose component is index, append the vertices in those cliques to the vertices of its all children, sort, compute its density
void reportSubgraph (int variant, vertex index, HashMap<vertex>& orderInFile, vector<vertex>& component, helpers& ax, vector<subcore>& skeleton, Graph& graph, edge nEdge, FILE* fp, FILE* gp) {

	if (skeleton[index].parent == -1) {
		skeleton[index].size = graph.size();
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		fprintf(fp, "%d %d %d %d %lf %d ", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0);
		fprintf(fp, "-1\n");
		return;
	}

	vector<vertex> vset;
	if (variant == 12 || variant == 13 || variant == 14) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == index)
				vset.push_back (i);
	}
	else if (variant == 23 || variant == 24) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == index) {
				vset.push_back (get<0>((*ax.el)[i]));
				vset.push_back (get<1>((*ax.el)[i]));
			}
	}
	else if (variant == 34) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == index) {
				vset.push_back (get<0>((*ax.tris)[i]));
				vset.push_back (get<1>((*ax.tris)[i]));
				vset.push_back (get<2>((*ax.tris)[i]));
			}
	}

	bool pass = true;
	pass = pullChildrenSets (fp, skeleton[index].children, orderInFile, vset, skeleton);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	pass = hashUniquify (vset);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	// edge density
	edge edge_count = 0;
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++)
			edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;
	skeleton[index].nEdge = edge_count;
	skeleton[index].size = vset.size();
	if (vset.size() > 1)
		skeleton[index].ed = (double) edge_count / (skeleton[index].size * (skeleton[index].size - 1) / 2);

	bool highlight = (skeleton[index].children.empty() && skeleton[index].ed >= THRESHOLD && skeleton[index].size >= LOWERBOUND) ? true : false;
	if (highlight)
		fprintf(gp, "id: %lld  K: %d  |V|: %d  |E|: %d  ed: %.2lf  LEAF?: %d  parent id: %lld\t", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge,
				skeleton[index].ed,	skeleton[index].children.empty()?1:0, skeleton[index].parent);

	fprintf(fp, "%d %d %d %d %lf %d %d\t", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0, skeleton[index].parent);
	for (size_t i = 0; i < vset.size(); i++) {
		fprintf(fp, "%d ", vset[i]);
		if (highlight)
			fprintf(gp, "%d ", vset[i]);
	}
	fprintf(fp, "-1\n");
	if (highlight)
		fprintf(gp, "-1\n");
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

void presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, edge nEdge, helpers& ax, string vfile, FILE* gp) {

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
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, graph, nEdge, fp, gp);
			removeChild (i, backup);
		}
	}
	fclose (fp);
}
