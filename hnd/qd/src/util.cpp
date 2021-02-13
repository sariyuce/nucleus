#include "main.h"


double color (double ed) {
        double a = (1 - ed) / 0.25;
        int X = floor(a);
        int Y = floor(255 * (a-X));
        if (X < 4)
                return (3 - X) + (255 - (double) Y) / 255;
        else
                return 0;
}

void print_nested_circle (vector<subcore>& hrc, int ind, FILE* fp, string cfl) {
	if (hrc[ind].size < LOWERBOUND)
            return;
        double parent_color = color(hrc[ind].ed);
        fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"index\": \"%d\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld",
                        parent_color, cfl.c_str(), ind, hrc[ind].size, hrc[ind].ed, hrc[ind].K, hrc[ind].size);
        if (hrc[ind].children.size() == 1) {
                fprintf(fp, ", \"children\": [\n");
                int ch = hrc[ind].children[0];
                // ghost child
                fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ",
                                parent_color, hrc[hrc[ch].parent].size - hrc[ch].size);
                // real child
                print_nested_circle (hrc, ch, fp, cfl);
                fprintf(fp, "\n]\n");
        }
        else if (hrc[ind].children.size() > 1) {
                fprintf(fp, ", \n\"children\": [\n");
                size_t i;
                for (i = 0; i < hrc[ind].children.size() - 1; i++) {
                        print_nested_circle (hrc, hrc[ind].children[i], fp, cfl);
                        fprintf(fp, ",\n");
                }
                print_nested_circle (hrc, hrc[ind].children[i], fp, cfl);
                fprintf(fp, "\n]");
        }
        fprintf(fp, "}\n");
}

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

bool pullChildrenSets (FILE* fp, vector<vertex>& children, unordered_map<vertex, vertex>& orderInFile, vector<vertex>& vset, vector<subcore>& skeleton) {

	int limit = UPPERBOUND;
	char c;
	for (vertex eda : children) {
		if (skeleton[eda].size == -1)
			return false;

		if (orderInFile.find (eda) == orderInFile.end()) {
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
				if (vset.size() > limit) {
					fseek (fp, 0, SEEK_END);
					return false;
				}
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
	fprintf(fp, "%d %d %d %d %lf %d %d -1 \n", index,	sc->K, sc->size, sc->nEdge, sc->ed, sc->children.empty()?1:0, sc->parent);
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
void reportSubgraph (int variant, vertex index, unordered_map<vertex, vertex>& orderInFile, vector<vertex>& component, helpers& ax, vector<subcore>& skeleton, Graph& graph, edge nEdge, FILE* fp, FILE* gp) {

	if (skeleton[index].parent == -1) {
		skeleton[index].size = graph.size();
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		fprintf(fp, "%d %d %d %d %lf %d ", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0);
		fprintf(fp, "-1\n");
		return;
	}

//	printf ("LEAF: ");
	vector<vertex> vset;
	if (variant == 12 || variant == 13 || variant == 14) {
		for (vertex i = 0; i < component.size(); i++) {
			if (component[i] == index)
				vset.push_back (i);
			if (skeleton[index].children.empty()) {
//				printf ("%d ", i);
			}
		}
	}
	else if (variant == 23 || variant == 24) {
		for (vertex i = 0; i < component.size(); i++) {
			if (component[i] == index) {
				vset.push_back (get<0>((*ax.el)[i]));
				if (get<1>((*ax.el)[i]) > 0)
					vset.push_back (get<1>((*ax.el)[i]));
				else
					vset.push_back (-1 * get<1>((*ax.el)[i]));

				if (skeleton[index].children.empty()) {
//					printf ("%d %d ", get<0>((*ax.el)[i]), get<1>((*ax.el)[i]));
				}
			}
		}
	}
	else if (variant == 34) {
		for (vertex i = 0; i < component.size(); i++) {
			if (component[i] == index) {
				vset.push_back (get<0>((*ax.tris)[i]));
				vset.push_back (get<1>((*ax.tris)[i]));
				vset.push_back (get<2>((*ax.tris)[i]));

				if (skeleton[index].children.empty()) {
//					printf ("%d %d %d ", get<0>((*ax.tris)[i]), get<1>((*ax.tris)[i]), get<2>((*ax.tris)[i]));
				}
			}
		}
	}

//	if (skeleton[index].children.empty()) {
//		printf ("\n");
//	}
	bool pass = true;
	pass = pullChildrenSets (fp, skeleton[index].children, orderInFile, vset, skeleton);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	vector<vertex> backup_vset (vset);
	pass = hashUniquify (vset);
	if (!pass) {
		dummyLine (&skeleton[index], fp, index);
		return;
	}

	// edge density
	edge edge_count = 0;
//	printf ("vset.size: %d\n", vset.size());
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++) {
			vertex u = vset[i];
			vector<vertex> ret;
			outgoings (graph[u], ret);
			vector<vertex> abc;
			for (auto v: ret)
				abc.push_back(graph[u][v]);
			sort (abc.begin(), abc.end());
			edge_count += commons (vset, ret);

			ret.clear();
			abc.clear();
			incomings (graph[u], ret);
			for (auto v: ret)
				abc.push_back (M2P (graph[u][v]));
			sort (abc.begin(), abc.end());
			edge_count += commons (vset, abc);

			ret.clear();
			abc.clear();
			undirecteds (graph[u], ret);
			for (auto v: ret)
				abc.push_back (graph[u][v]);
			sort (abc.begin(), abc.end());
			edge_count += commons (vset, abc);

		}

//	edge_count /= 2;
	skeleton[index].nEdge = edge_count;
	skeleton[index].size = vset.size();
	if (vset.size() > 1)
		skeleton[index].ed = (double) edge_count / (skeleton[index].size * (skeleton[index].size - 1) / 2);

	//bool highlight = (skeleton[index].children.empty() && skeleton[index].ed >= THRESHOLD && skeleton[index].size >= LOWERBOUND) ? true : false;
	//if (highlight)
	//	fprintf(gp, "id: %lld  K: %d  |V|: %d  |E|: %d  ed: %.2lf  LEAF?: %d  parent id: %lld\t", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge,
	//			skeleton[index].ed,	skeleton[index].children.empty()?1:0, skeleton[index].parent);

	fprintf(fp, "%d %d %d %d %lf %d %d\t", index, skeleton[index].K, skeleton[index].size, skeleton[index].nEdge, skeleton[index].ed, skeleton[index].children.empty()?1:0, skeleton[index].parent);

	for (size_t i = 0; i < vset.size(); i++) {
		fprintf(fp, "%d ", vset[i]);
		//if (highlight)
		//	fprintf(gp, "%d ", vset[i]);
	}

	fprintf(fp, "-1\n");
	//if (highlight)
	//	fprintf(gp, "-1\n");
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

inline void findRepresentative (vertex* child, vector<subcore>& skeleton) {
	vertex u = *child;
	if (skeleton[u].parent != -1) {
		vertex pr = skeleton[u].parent;
		while (skeleton[u].K == skeleton[pr].K) {
			u = pr;
			if (skeleton[u].parent != -1)
				pr = skeleton[u].parent;
			else
				break;
		}
	}
	*child = u;
}

void presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, edge nEdge, helpers& ax, string vfile, FILE* gp) {

	// assign unassigned items to top subcore
	for (vertex i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	unordered_map<vertex, vertex> replace;
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
		if (replace.find (component[i]) != replace.end())
			component[i] = replace[component[i]];

	stack<vertex> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	string nFile = vfile + "_NUCLEI";
	FILE* fp = fopen (nFile.c_str(), "w+");
	vector<subcore> backup (skeleton);

	unordered_map<vertex, vertex> orderInFile; // key is the skeleton index, value is the order
	vertex o = 0; // order of subcores in file



	while (!subcoreStack.empty()) {
		vertex i = subcoreStack.top();
		subcoreStack.pop();
		if (backup[i].visible) { // && backup[i].children.empty()) {
			orderInFile[i] = o++;
			reportSubgraph (variant, i, orderInFile, component, ax, skeleton, graph, nEdge, fp, gp);
			removeChild (i, backup);
		}
	}
	fclose (fp);


    string temp (vfile);
    string cfl = temp + "_circle.json";
    FILE* gip = fopen (cfl.c_str(), "w");
    print_nested_circle (skeleton, skeleton.size() - 1, gip, cfl);
    fclose(gip);
}






// ids of only incoming edges
void incomings (vector<vertex>& a, vector<vertex>& ret) {
	if (a.empty())
		return;
	vertex i = 1, ni = a[0];
	while (1) {
		if (i == a[0]) {
			while (ni < a.size()) {
				ret.push_back (ni);
				ni++;
			}
			break;
		}
		else if (ni == a.size())
			break;
		if (M2P (a[ni]) < a[i]) {
			ret.push_back (ni);
			ni++;
		}
		else if (a[i] < M2P (a[ni]))
			i++;
		else {
			i++;
			ni++;
		}
	}
}

// undirecteds for which u is the smaller node
void outgoings_and_asymmetric_undirecteds (vertex u, vector<vertex>& b, vector<vertex>& ret) {
	if (b.empty())
		return;
	vertex j = 1, nj = b[0];
	while (1) {
		if (j == b[0])
			break;
		else if (nj == b.size()) {
			while (j < b[0]) {
				ret.push_back (j);
				j++;
			}
			break;
		}
		if (M2P (b[nj]) < b[j])
			nj++;
		else if (b[j] < M2P (b[nj])) {
			ret.push_back (j);
			j++;
		}
		else {
			if (u < b[j])
				ret.push_back(-j);
			j++;
			nj++;
		}
	}
}

void outgoings_and_undirecteds (vector<vertex>& b, vector<vertex>& ret) {
	if (b.empty())
		return;
	vertex j = 1, nj = b[0];
	while (1) {
		if (j == b[0])
			break;
		else if (nj == b.size()) {
			while (j < b[0]) {
				ret.push_back (j);
				j++;
			}
			break;
		}
		if (M2P (b[nj]) < b[j])
			nj++;
		else if (b[j] < M2P (b[nj])) {
			ret.push_back (j);
			j++;
		}
		else {
			ret.push_back(j);
			j++;
			nj++;
		}
	}
}

// ids of only outgoing edges
void outgoings (vector<vertex>& b, vector<vertex>& ret) {
	if (b.empty())
		return;
	vertex j = 1, nj = b[0];
	while (1) {
		if (j == b[0])
			break;
		else if (nj == b.size()) { // means no more incoming edges exist
			while (j < b[0]) {
				ret.push_back (j);
				j++;
			}
			break;
		}
		// takes the diff of all outgoings from all incomings
		if (M2P (b[nj]) < b[j])
			nj++;
		else if (b[j] < M2P (b[nj])) {
			ret.push_back (j);
			j++;
		}
		else {
			j++;
			nj++;
		}
	}
}


// ids of bidirectional edges
void undirecteds (vector<vertex>& b, vector<vertex>& ret) {
	if (b.empty())
		return;
	vertex j = 1, nj = b[0];
	while (j < nj && nj < b.size()) {
		if (M2P (b[nj]) < b[j])
			nj++;
		else if (b[j] < M2P (b[nj]))
			j++;
		else {
			ret.push_back (j);
			j++;
			nj++;
		}
	}
}


// ids of bidirectional edges
void asymmetric_undirecteds (vertex u, vector<vertex>& b, vector<vertex>& ret) {
	if (b.empty())
		return;
	vertex j = 1, nj = b[0];
	while (j < nj && nj < b.size()) {
		if (M2P (b[nj]) < b[j])
			nj++;
		else if (b[j] < M2P (b[nj]))
			j++;
		else {
			if (u < b[j])
				ret.push_back (j);
			j++;
			nj++;
		}
	}
}


//// ids of bidirectional edges
//void undirecteds (vertex b, Graph& dg, vector<vertex>& ret) {
//	vertex j = 1, nj = dg[b][0];
//	while (j < nj && nj < dg[b].size()) {
//		if (M2P (dg[b][nj]) < dg[b][j])
//			nj++;
//		else if (dg[b][j] < M2P (dg[b][nj]))
//			j++;
//		else {
//			ret.push_back (j);
//			j++;
//			nj++;
//		}
//	}
//}


void inte (vector<vertex>& a, vector<vertex>& b, vector<vertex>& c, vector<vertex>& d, vector<vertex>& ret) {
	vertex i, j;
	i = j = 0;
	vertex u, v;
	while (i < c.size() && j < d.size()) {
		u = (a[c[i]] < 0) ? M2P (a[c[i]]) : a[c[i]];
		v = (b[d[j]] < 0) ? M2P (b[d[j]]) : b[d[j]];;
		if (u < v)
			i++;
		else if (v < u)
			j++;
		else {
			ret.push_back(c[i]);
			ret.push_back(d[j]);
			i++;
			j++;
		}
	}
}


void inter (int F, int G, Graph& dg, vertex a, vertex b, vector<vertex>& ret) {
	if (dg[a].size() == 1 || dg[b].size() == 1)
		return;
	vector<vertex> c, d;
	if (F == 0)	undirecteds (dg[a], c);
	else if (F == 1) incomings (dg[a], c);
	else if (F == 2) outgoings (dg[a], c);

	if (G == 0)	undirecteds (dg[b], d);
	else if (G == 1) incomings (dg[b], d);
	else if (G == 2) outgoings (dg[b], d);

	inte (dg[a], dg[b], c, d, ret); // intersection of c and d
}

