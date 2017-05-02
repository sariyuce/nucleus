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

	char c;
	for (vertex eda : children) {
		if (skeleton[eda].size == -19)
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
		fscanf (fp, "%d %d %d %d %lf %d", &d, &d, &d, &d, &f, &d); // id, K, |V|, |E|, ed, LEAF?
		while (fscanf (fp, "%d", &d) != EOF)
			if (d != -1)
				vset.push_back (d);
			else
				break;
		fseek (fp, 0, SEEK_END);
	}
	return true;
}

bool hashUniquify (vector<vertex>& vertices) {
	HashMap<bool> hermap (false);
	for (size_t i = 0; i < vertices.size(); i++) {
		int t = vertices[i];
		if (hermap.hasDefaultValue (t))
			hermap[t] = true;
		else {
			vertices.erase (vertices.begin() + i);
			i--;
		}
		if (i > UPPERBOUND)
			return false;
	}
	sort (vertices.begin(), vertices.end());
	return true;
}

inline void dummyLine (subcore* sc, FILE* fp, vertex ind) {
	sc->size = -1;
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

inline double color (double ed) {
	double a = (1 - ed) / 0.25;
	int X = floor(a);
	int Y = floor(255 * (a-X));
	if (X < 4)
		return (3 - X) + (255 - (double) Y) / 255;
	else
		return 0;
}

void nestedCircles (vector<subcore>& skeleton, vertex ind, FILE* fp, string cfl) {
	double parent_color = color(skeleton[ind].ed);
	fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"index\": \"%d\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld",
			parent_color, cfl.c_str(), ind, skeleton[ind].size, skeleton[ind].ed, skeleton[ind].K, skeleton[ind].size);
	if (skeleton[ind].children.size() == 1) {
		fprintf(fp, ", \"children\": [\n");
		int ch = skeleton[ind].children[0];
		// ghost child
		fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ",
				parent_color, skeleton[skeleton[ch].parent].size - skeleton[ch].size);
		// real child
		nestedCircles (skeleton, ch, fp, cfl);
		fprintf(fp, "\n]\n");
	}
	else if (skeleton[ind].children.size() > 1) {
		fprintf(fp, ", \n\"children\": [\n");
		size_t i;
		for (i = 0; i < skeleton[ind].children.size() - 1; i++) {
			nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl);
			fprintf(fp, ",\n");
		}
		nestedCircles (skeleton, skeleton[ind].children[i], fp, cfl);
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
void reportSubgraph (int variant, vertex ind, HashMap<vertex>& orderInFile, vector<vertex>& component, helpers& ax, vector<subcore>& skeleton, Graph& graph, edge nEdge, FILE* fp) {

	if (skeleton[ind].parent == -1) {
		skeleton[ind].size = graph.size();
		skeleton[ind].nedge = nEdge;
		skeleton[ind].ed = 0;
		fprintf(fp, "%d %d %d %d %lf %d ", ind, skeleton[ind].K, skeleton[ind].size, skeleton[ind].nedge, skeleton[ind].ed, skeleton[ind].children.empty()?1:0);
		fprintf(fp, "-1\n");
		return;
	}

	vector<vertex> vset;
	if (variant == 12 || variant == 13 || variant == 14) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == ind)
				vset.push_back (i);
	}
	else if (variant == 23 || variant == 24) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == ind) {
				vset.push_back (get<0>((*ax.el)[i]));
				vset.push_back (get<1>((*ax.el)[i]));
			}
	}
	else if (variant == 34) {
		for (vertex i = 0; i < component.size(); i++)
			if (component[i] == ind) {
				vset.push_back (get<0>((*ax.tris)[i]));
				vset.push_back (get<1>((*ax.tris)[i]));
				vset.push_back (get<2>((*ax.tris)[i]));
			}
	}


	bool pass = true;
	pass = pullChildrenSets (fp, skeleton[ind].children, orderInFile, vset, skeleton);
	if (!pass) {
		dummyLine (&skeleton[ind], fp, ind);
		return;
	}

	pass = hashUniquify (vset);
	if (!pass) {
		dummyLine (&skeleton[ind], fp, ind);
		return;
	}

	// edge density
	edge edge_count = 0;
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++)
				edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;
	skeleton[ind].nedge = edge_count;
	skeleton[ind].size = vset.size();
	if (vset.size() > 1)
		skeleton[ind].ed = (double) edge_count / (skeleton[ind].size * (skeleton[ind].size - 1) / 2);

	fprintf(fp, "%d %d %d %d %lf %d ", ind, skeleton[ind].K, skeleton[ind].size, skeleton[ind].nedge, skeleton[ind].ed, skeleton[ind].children.empty()?1:0);
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
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, graph, nEdge, fp);
			removeChild (i, backup);
		}
	}
	fclose (fp);

	for (size_t i = 0; i < skeleton.size(); i++)
		if (skeleton[i].visible) {
			printf ("%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].size, skeleton[i].nedge, skeleton[i].ed, skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
			fprintf (gp, "%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, skeleton[i].K, skeleton[i].size, skeleton[i].nedge, skeleton[i].ed, skeleton[i].parent, skeleton[i].children.empty()?"LEAF":"INTERMEDIATE");
			fflush (gp);
		}

	// build the json file for visualization -- based on the interactive nested circles
	string cfl = vfile + "_circle.json";
	FILE* ncp = fopen (cfl.c_str(), "w");
	nestedCircles (skeleton, skeleton.size() - 1, ncp, cfl);
	fclose (ncp);
}
