#include "main.h"

void cleanminusones (Graph& graph) {

	//		printf ("*** BEFORE CLEAN ***\n");
	//		for (vertex i = 0; i < nVtx; ++i) {
	//			printf ("%d: ", i);
	//			for (vertex j = 0; j < graph[i].size(); j++)
	//				printf ("%d ", graph[i][j]);
	//			printf ("\n");
	//		}
	//		printf ("**********\n");

	// clean -1s
	for (size_t i = 0; i < graph.size(); i++) {
		sort (graph[i].begin(), graph[i].end());
		bool fl = false;
		for (size_t j = 0; j < graph[i].size(); j++) {
			if (graph[i][j] > -1) {
				graph[i].erase (graph[i].begin(), graph[i].begin() + j);
				fl = true;
				break;
			}
		}
		if (!fl)
			graph[i].erase (graph[i].begin(), graph[i].end());
	}

	//		printf ("*** AFTER CLEAN ***\n");
	//		for (vertex i = 0; i < nVtx; ++i) {
	//			printf ("%d: ", i);
	//			for (vertex j = 0; j < graph[i].size(); j++)
	//				printf ("%d ", graph[i][j]);
	//			printf ("\n");
	//		}
	//		printf ("**********\n");
}

void print_time (FILE* fp, const string& str, const timestamp& t) {
	fprintf (fp, "%s %d.%06d\n", str.c_str(), t.seconds, t.microseconds);
	fflush(fp);
}

void print_time_new (FILE* fp, const string& str, const timestamp& t) {
	cout << str << t << endl;
	fprintf (fp, "%s %d.%06d\n", str.c_str(), t.seconds, t.microseconds);
	fflush(fp);
}

void print_time_new2 (FILE* fp, const string& str, const timestamp& t, const string& str2) {
	cout << str << t << str2 << endl;
	fprintf (fp, "%s %d.%06d %s\n", str.c_str(), t.seconds, t.microseconds, str2.c_str());
	fflush(fp);
}

bool sr (couple1 a, couple1 b) {
	int u = get<0>(a);
	int v = get<1>(a);

	int w = get<0>(b);
	int x = get<1>(b);

	if (u != w)
		return (u < w);
	else
		return (v < x);
}

void uniquify_es (vector<couple1>& ch) {
	sort (ch.begin(), ch.end(), sr);
	for (size_t i = 1; i < ch.size(); i++) {
		if (ch[i] == ch[i-1]) {
			ch.erase (ch.begin() + i);
			i--;
		}
	}
}

bool tr (triple2 a, triple2 b) {
	int u = get<0>(a);
	int v = get<1>(a);
	int w = get<2>(a);

	int x = get<0>(b);
	int y = get<1>(b);
	int z = get<2>(b);

	if (u != x)
		return (u < x);
	else {
		if (v != y)
			return (v < y);
		else
			return (w < z);
	}
}

void uniquify_ts (vector<triple2>& ch) {
	sort (ch.begin(), ch.end(), tr);
	for (size_t i = 1; i < ch.size(); i++) {
		if (ch[i] == ch[i-1]) {
			ch.erase (ch.begin() + i);
			i--;
		}
	}
}

void uniquify_vs (vector<int>& ch) {
	sort (ch.begin(), ch.end());
	for (size_t i = 1; i < ch.size(); i++) {
		if (ch[i] == ch[i-1]) {
			ch.erase (ch.begin() + i);
			i--;
		}
	}
}

bool hash_uniquify_vs (vector<int>& ch) {
	HashMap<bool> hermap (false);
	for (size_t i = 0; i < ch.size(); i++) {
		int t = ch[i];
		if (hermap.hasDefaultValue (t)) {
			hermap[t] = true;
		}
		else {
			ch.erase (ch.begin() + i);
			i--;
		}
		if (i > SUBGRAPHSIZELIMIT)
			return false;
	}

	sort (ch.begin(), ch.end());
	return true;
}


bool connected (Graph& us_graph, vertex u, vertex v) {
	vertex a = u<v?u:v;
	vertex b = u>=v?u:v;
	for (size_t i = 0; i < us_graph[a].size(); i++) {
		if (us_graph[a][i] == b)
			return true;
	}
	return false;
}

bool connected_deg (Graph& graph, Graph& ordered_graph, vertex u, vertex v) {
	vertex a, b;
	if (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v)) {
		a = u; b = v;
	}
	else {
		a = v; b = u;
	}

	for (size_t i = 0; i < ordered_graph[a].size(); i++) {
		if (ordered_graph[a][i] == b)
			return true;
	}
	return false;
}

void intersection (vector<vertex>& vertices, vector<vertex>& edgelist, vector<vertex>& new_vertex) {
	size_t i = 0, j = 0;
	while (i < vertices.size() && j < edgelist.size()) {
		if (vertices[i] < edgelist[j])
			i++;
		else if (edgelist[j] < vertices[i])
			j++;
		else {
			new_vertex.push_back(vertices[i]);
			i++;
			j++;
		}
	}
}

inline int intersect_count (vector<vertex>& a, vector<vertex>& b) {
	size_t i = 0, j = 0;
	int count = 0;
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

bool less_than (vertex u, vertex v, Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
}


inline bool exists (int val, vector<int>& v) {
	for (size_t i = 0; i < v.size(); i++)
		if (v[i] == val)
			return true;
	return false;
}

bool take_vsets (FILE* fp, vector<int>& toget, HashMap<int>& ordermap, vector<vertex>& vset, vector<subcore>& hrc) {
	char c;
	for (int er : toget) {
		if (hrc[er].size == -19)
			return false;

		if (ordermap.hasDefaultValue (er)) {
			printf ("PROBLEM: %d has -1 as order\n", er);
			exit(1);
		}

		int sc = ordermap[er];
		fseek (fp, 0, SEEK_SET);
		int ln = 0;
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
		// now you are on the correct place of file and can get the vlist
		int d;
		double f;
		fscanf (fp, "%d", &d);
		fscanf (fp, "%d", &d);
		fscanf (fp, "%d", &d);
		fscanf (fp, "%d", &d);
		fscanf (fp, "%lf", &f);
		fscanf (fp, "%d", &d);
		int v;
		while (fscanf (fp, "%d", &v) != EOF) {
			if (v != -1)
				vset.push_back (v);
			else
				break;
		}
		fseek (fp, 0, SEEK_END);
	}

	return true;
}

#ifndef SHORT
void report_results (int variant, int ind, HashMap<int>& ordermap, vector<vertex>& subx_ids, p_auxies& ax, vector<subcore>& hrc, Graph& graph, FILE* fp,
		timestamp& ph1, timestamp& ph2, timestamp& ph3, timestamp& ph4) {
	timestamp t1;
	if (hrc[ind].parent == -1) {
		hrc[ind].size = graph.size();
		int all_possible_edges = hrc[ind].size * (hrc[ind].size - 1) / 2;
		hrc[ind].nedge = 0;
		for (vector<int> g : graph)
			hrc[ind].nedge += g.size();
		hrc[ind].nedge /= 2;
		hrc[ind].ed = (double) hrc[ind].nedge / all_possible_edges;
		fprintf(fp, "%d %d %d %d %lf %d ", ind, hrc[ind].K, hrc[ind].size, hrc[ind].nedge, hrc[ind].ed, hrc[ind].children.empty()?1:0);
		fprintf(fp, "-1 %d\n", hrc[ind].children.size());
		return;
	}


	vector<vertex> vset;
	if (variant == 12) {
		for (int i = 0; i < subx_ids.size(); i++) {
			if (subx_ids[i] == ind)
				vset.push_back (i);
		}
	}
	else if (variant == 23) {
		for (int i = 0; i < subx_ids.size(); i++) {
			if (subx_ids[i] == ind) {
				vset.push_back (get<0>((*ax.el)[i]));
				vset.push_back (get<1>((*ax.el)[i]));
			}
		}
	}
	else if (variant == 34) {
		for (int i = 0; i < subx_ids.size(); i++) {
			if (subx_ids[i] == ind) {
				vset.push_back (get<0>((*ax.tl)[i].triple));
				vset.push_back (get<1>((*ax.tl)[i].triple));
				vset.push_back (get<2>((*ax.tl)[i].triple));
			}
		}
	}

	bool pass = true;
	pass = take_vsets (fp, hrc[ind].children, ordermap, vset, hrc);
	if (!pass) {
		hrc[ind].size = -19;
		fprintf(fp, "%d %d %d %d %lf %d ", ind, hrc[ind].K, hrc[ind].size, hrc[ind].nedge, hrc[ind].ed, hrc[ind].children.empty()?1:0);
		fprintf(fp, "-1 %d\n", hrc[ind].children.size());
		return;
	}
	pass = hash_uniquify_vs (vset);
	if (!pass) {
		hrc[ind].size = -19;
		fprintf(fp, "%d %d %d %d %lf %d ", ind, hrc[ind].K, hrc[ind].size, hrc[ind].nedge, hrc[ind].ed, hrc[ind].children.empty()?1:0);
		fprintf(fp, "-1 %d\n", hrc[ind].children.size());
		return;
	}
	timestamp t2;
	ph1 += t2 - t1;

	// edge density
	int edge_count = 0;
#ifndef NODENSITY
	if (vset.size() < SUBGRAPHSIZELIMIT) {
		for (size_t i = 0; i < vset.size(); i++)
			edge_count += intersect_count (vset, graph[vset[i]]);
	}
#endif
	edge_count /= 2;
	timestamp t3;
	ph2 += t3 - t2;

	hrc[ind].size = vset.size();
	int all_possible_edges = hrc[ind].size * (hrc[ind].size - 1) / 2;
	hrc[ind].nedge = edge_count;
	hrc[ind].ed = -1;
	if (all_possible_edges > 0)
		hrc[ind].ed = (double) edge_count / all_possible_edges;

	timestamp t4;
	ph3 += t4 - t3;

	fprintf(fp, "%d %d %d %d %lf %d ", ind, hrc[ind].K, hrc[ind].size, hrc[ind].nedge, hrc[ind].ed, hrc[ind].children.empty()?1:0);
	for (size_t i = 0; i < vset.size(); i++)
		fprintf(fp, "%d ", vset[i]);

	fprintf(fp, "-1 %d\n", hrc[ind].children.size());

	timestamp t5;
	ph4 += t5 - t4;
}
#endif

void intersect2 (Graph& graph, vertex u, vertex v, vector<vertex>& intersection) {
	size_t i = 0, j = 0;
	while (i < graph[u].size() && j < graph[v].size()) {
		if (graph[u][i] < graph[v][j])
			i++;
		else if (graph[v][j] < graph[u][i])
			j++;
		else {
			intersection.push_back(graph[u][i]);
			i++;
			j++;
		}
	}
}

void createOrdered (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, couple1* el, edge* xel, vertex* ordered_adj, edge* ordered_xadj) {
	edge xi = 0;
	vertex i = 0;
	xel[xi++] = 0;

	edge oxi = 0;
	vertex oi = 0;
	ordered_xadj[oxi++] = 0;

	for (vertex u = 0; u < nVtx; u++) {
		for (vertex j = xadj[u]; j < xadj[u+1]; j++) {
			vertex v = adj[j];
			if (lessThan (u, v, xadj)) {
				ordered_adj[oi++] = v;
				couple1 c = make_tuple(u, v);
				el[i++] = c;
			}
		}
		ordered_xadj[oxi++] = oi;
		xel[xi++] = i;
	}
}


vertex create_ordered_graph (Graph& graph, EdgeList2& el, vector<vertex>& xel, Graph& ordered_graph) {
	vertex nEdge = 0;
	xel.push_back(0);
	ordered_graph.resize(graph.size());
	for (size_t u = 0; u < graph.size(); u++) {
		nEdge += graph[u].size();
		for (size_t j = 0; j < graph[u].size(); j++) {
			vertex v = graph[u][j];
			if (less_than (u, v, graph)) {
				ordered_graph[u].push_back(v);
				couple1 c = make_tuple(u, v);
				el.push_back(c);
			}
		}
		ordered_graph[u].shrink_to_fit();
		xel.push_back(el.size());
	}
	ordered_graph.shrink_to_fit();
	el.shrink_to_fit();
	xel.shrink_to_fit();
	return nEdge;
}

#if 0
// vsets from cpc, childrens from hrc
void print_tree (vector<subcore>& hrc, vector<subcore>& cpc, int ind, bool print, set<int>& tolabel, int len, FILE* fp) {

	if (cpc[ind].vset.size() >= THRESHOLD) {
		if (hrc[ind].children.size() == 1) {
			if (print) {
				fprintf(fp, "\t\t %d -- ", ind);
				tolabel.insert(ind);
			}
			if (hrc[hrc[ind].children[0]].children.size() != 1) {
				if (len != 1)
					fprintf(fp, "%d [penwidth=%d];\n", hrc[ind].children[0], (int) ceil (len/ (double) 5)); // (int) ceil(log(len)));
				else
					fprintf(fp, "%d;\n", hrc[ind].children[0]);

				tolabel.insert(hrc[ind].children[0]);
				print_tree (hrc, cpc, hrc[ind].children[0], true, tolabel, 1, fp);
			}
			else
				print_tree (hrc, cpc, hrc[ind].children[0], false, tolabel, ++len, fp);
		}
		else {
			for (size_t i = 0; i < hrc[ind].children.size(); i++) {
				if (cpc[hrc[ind].children[i]].vset.size() >= THRESHOLD) {
					fprintf(fp, "\t\t %d -- %d;\n", ind, hrc[ind].children[i]);
					tolabel.insert(ind);
					tolabel.insert(hrc[ind].children[i]);
					print_tree (hrc, cpc, hrc[ind].children[i], true, tolabel, 1, fp);
				}
			}
		}
	}
}

// vsets from cpc, childrens from hrc
void report_tree (vector<subcore>& hrc, vector<subcore>& cpc, FILE* fp) {

	set<int> tolabel;
	fprintf(fp, "graph graphname {\n");
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].parent == -1) {
			print_tree (hrc, cpc, i, true, tolabel, 1, fp);
		}
	}

	for (set<int>::iterator it = tolabel.begin(); it != tolabel.end(); it++) {
		int size = cpc[*it].vset.size();
		string shape;
		if (size > 1 && size <= 100)
			shape = "circle";
		else if (size > 100 && size <= 1000)
			shape = "hexagon";
		else if (size > 1000 && size <= 10000)
			shape = "square";
		else if (size > 10000)
			shape = "triangle";

		int sz = floor(log (cpc[*it].vset.size()));
		double a = (1 - hrc[*it].ed) / 0.25;
		int X = floor(a);
		int Y = floor(255 * (a-X));
		int r, g, b;
		r = g = b = -1;

		switch(X) {
		case 0: r=255;g=Y;b=0;break;
		case 1: r=255-Y;g=255;b=0;break;
		case 2: r=0;g=255;b=Y;break;
		case 3: r=0;g=255-Y;b=255;break;
		case 4: r=0;g=0;b=255;break;
		}

		fprintf(fp, "%d [style=filled,label=\"\",fillcolor=\"#%.2x%.2x%.2x\",width=%d,height=%d,fontsize=32,shape=\"%s\"];\n",
				*it, r, g, b, sz, sz, shape.c_str());
	}

	fprintf(fp, "\t}\n");
}
#endif

double color (double ed) {
	double a = (1 - ed) / 0.25;
	int X = floor(a);
	int Y = floor(255 * (a-X));
	if (X < 4)
		return (3 - X) + (255 - (double) Y) / 255;
	else
		return 0;
}

#ifndef SHORT
// vsets from cpc, childrens from hrc
void print_nested_circle (vector<subcore>& hrc, int ind, FILE* fp, string cfl) {
	double parent_color = color(hrc[ind].ed);
	fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld", parent_color, cfl.c_str(), hrc[ind].size, hrc[ind].ed, hrc[ind].K, hrc[ind].size);
	if (hrc[ind].children.size() == 1) {
		fprintf(fp, ", \"children\": [\n");
		int ch = hrc[ind].children[0];
//		while (hrc[ch].children.size() == 1)
//			ch = hrc[ch].children[0];
		// ghost child
		fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ", parent_color, hrc[hrc[ch].parent].size - hrc[ch].size);
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
#endif

#if 0
// vsets from cpc, childrens and eds from hrc
void report_nested_circle (vector<subcore>& hrc, vector<subcore>& cpc, FILE* fp, string cfl) {
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].parent == -1 && (cpc[i].eset.size() > 0 || cpc[i].vset.size() > 0))
			print_nested_circle (hrc, i, fp, cfl);
	}
}
#endif

inline vertex getEdgeId (vertex u, vertex v, vector<vertex>& xel, EdgeList2& el, Graph& graph) {

	vertex a, b;
	if (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v)) {
		a = u;
		b = v;
	}
	else {
		a = v;
		b = u;
	}

	for (int i = xel[a]; i < xel[a+1]; i++) {
		if (get<1>(el[i]) == b)
			return i;
	}
}

void report_all_stuff (int variant, Graph& graph, int cn, auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile) {

}

string now () {
	srand (time(NULL));
	int r = rand();
	string str = to_string (r);
	return str;
}

inline bool edge_equal (couple1 i, couple1 j) {
	if ((get<0>(i) == get<0>(j) && get<1>(i) == get<1>(j)) || (get<0>(i) == get<1>(j) && get<1>(i) == get<0>(j)))
		return true;
	else
		return false;
}

inline bool tri_equal (triple2 i, triple2 j) {
	vector<int> a;
	a.push_back(get<0>(i));
	a.push_back(get<1>(i));
	a.push_back(get<2>(i));
	vector<int> b;
	b.push_back(get<0>(j));
	b.push_back(get<1>(j));
	b.push_back(get<2>(j));

	sort (a.begin(), a.end());
	sort (b.begin(), b.end());

	if (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
		return true;
	else
		return false;
}

#ifdef ZERO
bool nonemptyintersection (int variant, subcore& a, subcore& b) {
	if (variant == 12) {
		for (int i : a.vset) {
			for (int j : b.vset) {
				if (i == j)
					return true;
			}
		}
		return false;
	}
	else if (variant == 23) {
		for (couple1 i : a.eset) {
			for (couple1 j : b.eset) {
				if (edge_equal (i, j))
					return true;
			}
		}
		return false;
	}
	else if (variant == 34) {
		for (triple2 i : a.tset) {
			for (triple2 j : b.tset) {
				if (tri_equal (i, j))
					return true;
			}
		}
		return false;
	}
}

#ifndef SHORT
void dens_comp_ZERO (int variant, vector<subcore>& all_cores, int i, Graph& graph) {
	all_cores[i].nedge = 0;
	if (variant == 12) {
		sort (all_cores[i].vset.begin(), all_cores[i].vset.end());
		for (size_t k = 0; k < all_cores[i].vset.size(); k++)
			all_cores[i].nedge += intersect_count (all_cores[i].vset, graph[all_cores[i].vset[k]]);
		all_cores[i].nedge /= 2;
		all_cores[i].size = all_cores[i].vset.size();
		int all_possible_edges = all_cores[i].size * (all_cores[i].size - 1) / 2;
		if (all_possible_edges > 0)
			all_cores[i].ed = (double) all_cores[i].nedge / all_possible_edges;
	}
	else if (variant == 23) {
		vector<int> vset;
		for (couple1 c : all_cores[i].eset) {
			vset.push_back (get<0>(c));
			vset.push_back (get<1>(c));
		}
		uniquify_vs (vset);
		for (size_t k = 0; k < vset.size(); k++)
			all_cores[i].nedge += intersect_count (vset, graph[vset[k]]);
		all_cores[i].nedge /= 2;
		all_cores[i].size = vset.size();
		int all_possible_edges = all_cores[i].size * (all_cores[i].size - 1) / 2;
		if (all_possible_edges > 0)
			all_cores[i].ed = (double) all_cores[i].nedge / all_possible_edges;
	}
	else if (variant == 34) {
		vector<int> vset;
		for (triple2 c : all_cores[i].tset) {
			vset.push_back (get<0>(c));
			vset.push_back (get<1>(c));
			vset.push_back (get<2>(c));
		}
		uniquify_vs (vset);
		for (size_t k = 0; k < vset.size(); k++)
			all_cores[i].nedge += intersect_count (vset, graph[vset[k]]);
		all_cores[i].nedge /= 2;
		all_cores[i].size = vset.size();
		int all_possible_edges = all_cores[i].size * (all_cores[i].size - 1) / 2;
		if (all_possible_edges > 0)
			all_cores[i].ed = (double) all_cores[i].nedge / all_possible_edges;

	}
}
#endif

void report_ZERO (int variant, Graph& graph, int cn, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile, FILE* ffp) {
	timestamp ee1;
	size_t sz;
	if (variant == 34)
		sz = (*ax.tl).size();
	else if (variant >= 23)
		sz = (*ax.el).size();
	else
		sz = graph.size();

	double prev_rt = 0;
	vector<subcore> all_cores;
	for (int i = cn; i > 0; i--) {
		HashMap<bool> visited (false);
		if (exists.hasDefaultValue(i))
			continue;
		for (size_t j = 0; j < sz; j++) {
			if (visited.hasDefaultValue(j) && K[j] == i) {
				double rt = visited.size() * 100 / (double) sz;
				if (rt - prev_rt >= 1) {
					prev_rt = rt;
					printf ("%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					fprintf (ffp, "%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					fflush (ffp);
				}
				subcore c;
				c.K = i;
				c.parent = -1;
				if (variant == 34)
					max34 (j, i, graph, K, (*ax.tlist), (*ax.tl), c.tset, visited);
				else if (variant == 23)
					maxtruss (j, i, graph, K, (*ax.xel), (*ax.el), c.eset, visited);
				else if (variant == 12)
					maxcore (j, i, graph, K, c.vset, visited);
				all_cores.push_back (c);
			}
		}
	}
	timestamp ee2;
	cout << "Maximum k-core finding: " << ee2 - ee1 << endl;

	for (int i = 0; i < all_cores.size(); i++) {
		bool flag = false;
		for (int j = i + 1; j < all_cores.size(); j++) {
			if (all_cores[i].K > all_cores[j].K && nonemptyintersection (variant, all_cores[i], all_cores[j])) {
				all_cores[i].parent = j;
				all_cores[j].children.push_back (i);
#ifndef SHORT
				dens_comp_ZERO (variant, all_cores, i, graph);
#endif
				flag = true;
				break;
			}
		}
		if (!flag) {
#ifndef SHORT
			dens_comp_ZERO (variant, all_cores, i, graph);
#endif
		}
	}
	timestamp ee3;
	cout << "Hierarchy finding: " << ee3 - ee2 << endl;

#ifndef SHORT
	for (int i = 0; i < all_cores.size(); i++)
		printf("CORE-%d:  %d, MAX: %3d, size: (%4d, %4d), edge_density: %5lf\n", i, all_cores[i].children.empty()?1:0, all_cores[i].K,
				all_cores[i].size, all_cores[i].nedge, all_cores[i].ed);

	for (int i = 0; i < all_cores.size(); i++) {
		printf ("CORE %d contains ", i);
		for (int j : all_cores[i].children)
			printf ("%d ", j);
		printf ("\n");
	}
#endif

}
#endif

inline void remove_invisibles (vector<subcore>& cpc, int i) {
	for (int j = 0; j < cpc[i].children.size(); j++) {
		int ch = cpc[i].children[j];
		if (!cpc[ch].visible) {
			cpc[i].children.erase (cpc[i].children.begin() + j);
			j--;
		}
	}
}

void rearrange_children (vector<subcore>& hrc) {
	// rearrange children and parents based on visibility
	for (size_t i = 0; i < hrc.size(); i++)
		hrc[i].children.clear();

	for (size_t i = 0; i < hrc.size() - 1; i++) {
		if (hrc[i].visible) {
			int pr = hrc[i].parent;
			while (!hrc[pr].visible)
				pr = hrc[pr].parent;
			hrc[i].parent = pr;
			hrc[pr].children.push_back (i);
		}
	}
}

void report_all_stuff2 (int variant, Graph& graph, int nEdge, int cn, p_auxies& ax, HashMap<bool>& exists, vertex* K, const char* vfile, FILE* ffp) {

	timestamp t1;
	fprintf (ffp, "before subthings are being found\n");
	fflush (ffp);
	vector<subcore> hrc;
	size_t sz;

	if (variant == 34)
		sz = (*ax.tl).size();
	else if (variant == 23)
		sz = (*ax.el).size();
	else
		sz = graph.size();

	vector<vertex> subx_ids (sz, -1);

	timestamp tmerge (0, 0), tloop (0, 0), tdirect (0, 0);
	HashMap<bool> visited (false);
	double prev_rt = 0;
	stats op;
	op.find_op = op.union_op = op.adjust_op = 0; op.total_size = 0;
	int old_size = 0;
	int iter = 0;
	for (int i = cn; i > 0; i--) {
		int cnt = 0;
		if (exists.hasDefaultValue(i))
			continue;
		iter++;
		timestamp ff (0, 0);
		for (size_t j = 0; j < sz; j++) {
			if (visited.hasDefaultValue(j) && K[j] == i) {
				cnt++;
				double rt = visited.size() * 100 / (double) sz;
				if (rt - prev_rt >= 1) {
					prev_rt = rt;
//					printf ("%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					fprintf (ffp, "%.1lf subthing traversal (%d / %d)\n", rt, visited.size(), sz);
					fflush (ffp);
				}
//				if (variant == 34)
//					find_sub34 (j, graph, K, (*ax.tlist), (*ax.tl), visited, subx_ids, hrc, op);
//				else
				if (variant == 23)
					find_subtruss (j, graph, K, (*ax.xel), (*ax.el), visited, subx_ids, hrc, op);
				else {
//					printf ("find %d subcore from %d\n", i, j);
					find_subcore (j, graph, K, visited, subx_ids, hrc, tmerge, tloop, tdirect, op);
				}
			}
		}
		old_size = op.total_size;
	}

	timestamp t2;
	print_time (ffp, "all subthings finding time: ", t2 - t1);
	cout << "all subthings finding time: " << t2 - t1 << endl;

	// insert the root as a parent of all parentless subcores
	int nid = hrc.size();
	subcore sc;
	sc.rank = 1;
	sc.K = 0;
	sc.parent = -1;
	sc.root = -1;
	sc.visible = true;
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].parent == -1) {
			hrc[i].parent = nid;
			hrc[i].root = nid;
			sc.children.push_back(i);
		}
	}
	sc.nedge = nEdge;
	sc.ed = nEdge / double (graph.size() * (graph.size() - 1) / 2);
	hrc.push_back (sc);

	timestamp t3;
	print_time (ffp, "root core insertion time: ", t3 - t2);
	cout << "root core insertion time: " << t3 - t2 << endl;

	report_subgraphs (variant, hrc, subx_ids, graph, nEdge, ax, vfile, ffp);
}

inline bool doesitexistin (int val, vector<int>& v) {
	for (size_t i = 0; i < v.size(); i++)
		if (v[i] == val)
			return true;
	return false;
}

double report_density (vector<int>& vs, Graph& graph, int ismaximal, int curK, FILE* fp, FILE* ffp) {
	int edge_count = 0;
#ifndef SHORT
	for (size_t i = 0; i < vs.size(); i++)
		edge_count += intersect_count (vs, graph[vs[i]]);
	edge_count /= 2;
#endif
	int all_possible_edges = vs.size() * (vs.size() - 1) / 2;
	double ed = 0;
	if (all_possible_edges > 0)
		ed = (double) edge_count / all_possible_edges;
	if (vs.size() >= THRESHOLD) {
		// LeafFlag, K, |V|, |E|, density
		fprintf(fp, "%d %d %ld %d %lf    ", ismaximal, curK, vs.size(), edge_count, ed);
		fprintf(ffp, "%d, MAX: %3d, size: (%4ld, %4d), ratio: %5lf, edge_density: %5lf\n",
				ismaximal, curK, vs.size(), edge_count, (double) vs.size() / graph.size(), ed);
		for (size_t i = 0; i < vs.size(); i++)
			fprintf(fp, "%d ", vs[i]);
		fprintf(fp, "-1\n");
	}
	return ed;
}

void nested_circle (int co, mmap& new_kids, HashMap<int>& K_comps, HashMap<int>& parent,
		HashMap<double>& densities,	HashMap<int>& sizes, FILE* gip, string cfl) {
	double ed = densities[co];
	double parent_color = color(ed);
	int sz = sizes[co];
	fprintf(gip, "{\"color\": %lf, \"fl\": \"%s\", \"name\": \"%d %.2lf (%d)\", \"size\": %d", parent_color, cfl.c_str(), sz, ed, K_comps[co], sz);
	auto rn = new_kids.equal_range (co);
	int numch = distance (rn.first, rn.second);
	if (numch == 1) {
		fprintf(gip, ", \"children\": [\n");
		int ch = rn.first->second;
		// you may need no chains in the future
//			while (1) {
//				auto rn = new_kids.find (ch);
//				if (distance (rn.start, rn.end) == 1)
//					ch = rn.first->second;
//				else
//					break;
//			}
		// ghost child
		fprintf(gip, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %d}, ", parent_color, sz - sizes[ch]);
		// real child
		nested_circle (ch, new_kids, K_comps, parent, densities, sizes, gip, cfl);
		fprintf(gip, "\n]\n");
	}
	else if (numch > 1) {
		fprintf(gip, ", \n\"children\": [\n");
		auto it = rn.first;
		while (it != rn.second) {
			nested_circle (it->second, new_kids, K_comps, parent, densities, sizes, gip, cfl);
			it++;
			if (it == rn.second)
				fprintf(gip, "]");
			else
				fprintf(gip, ",\n");
		}
	}
	fprintf(gip, "}\n");
}

void compute_densities (int variant, HashMap<int>& K_comps, vector<int>& compid, HashMap<int>& parent, mmap& new_kids,
		Graph& graph, p_auxies& ax, HashMap<double>& densities,
		HashMap<int>& sizes, int id, const char* vfile, FILE* ffp, timestamp* cdtimes) {

	timestamp t1;
	for (size_t j = 0; j < compid.size(); j++)
		sizes[compid[j]]++;

	// do a bfs, collect comp ids in a 'bfs' list
	queue<int> order;
	int nedge = 0;
	for (vector<int> i : graph)
		nedge += i.size();
	densities[id] = (double) nedge / (graph.size() * (graph.size() - 1));
	K_comps[id] = -19;
	order.push (id);
	vector<int> vv;
	for (auto it = sizes.begin(); it != sizes.end(); it++) {
		int sc = it->first;
		if (parent.hasDefaultValue (sc)) {
			parent[sc] = id;
			vv.push_back (sc);
		}
	}
	for (int sc : vv)
		new_kids.emplace (id, sc);
	timestamp t2;
	cdtimes[0] += t2 - t1;

	vector<int> bfs;
	while (!order.empty()) {
		int co = order.front();
		order.pop();
		auto rn = new_kids.equal_range (co);
		for (auto iter = rn.first; iter != rn.second; iter++)
			order.push (iter->second);
		bfs.push_back (co);
	}
	// process 'bfs' in reverse order, compute density of each nucleus
	mmap all_kids;
	string tmp (vfile);
	string tfl = tmp + "_nuclei";
	char *a = new char[tfl.length() + 1];
	strcpy (a, tfl.c_str());
	FILE* fp = fopen (a, "w");
	vector<int> vs;
	vector<int> current_kids;
	timestamp t3;
	cdtimes[1] += t3 - t2;

	for (size_t i = bfs.size() - 1; i > 0; i--) {
		timestamp t31;
		current_kids.clear();
		int co = bfs[i];
		// accumulate current_kids
		current_kids.push_back (co);
		int ismaximal = 1;
		auto rn = all_kids.equal_range (co);
		for (auto iter = rn.first; iter != rn.second; iter++) {
			ismaximal = 0;
			current_kids.push_back (iter->second);
		}
		HashMap<bool> markcurkid (false);
		for (int u: current_kids)
			markcurkid[u] = true;
		timestamp t32;
		cdtimes[2] += t32 - t31;

		// find vertices/edges belonging to current subcore/subtruss
		int curK = K_comps[co];
		vs.clear();
		vs.shrink_to_fit();
		for (size_t j = 0; j < compid.size(); j++) {
			if (!markcurkid.hasDefaultValue (compid[j]))
				vs.push_back (j);
		}
		timestamp t33;
		cdtimes[3] += t33 - t32;

		// do density computation
		double ed = -1;
		if (variant == 12) {
			timestamp t1;
			sort (vs.begin(), vs.end());
			ed = report_density (vs, graph, ismaximal, curK, fp, ffp);
			timestamp t2;
			cdtimes[4] += t2 - t1;
		}
		else if (variant == 23) {
			HashMap<bool> seen (false);
			vector<vertex> vset;
			for (size_t j = 0; j < vs.size(); j++) {
				int u = get<0>((*ax.el)[vs[j]]);
				if (seen.hasDefaultValue(u)) {
					seen[u] = true;
					vset.push_back (u);
				}
				u = get<1>((*ax.el)[vs[j]]);
				if (seen.hasDefaultValue(u)) {
					seen[u] = true;
					vset.push_back (u);
				}
			}
			sizes[co] = vset.size();
			timestamp t1;
			sort (vset.begin(), vset.end());
			ed = report_density (vset, graph, ismaximal, curK, fp, ffp);
			timestamp t2;
			cdtimes[4] += t2 - t1;
		}
		else if (variant == 34) {
			HashMap<bool> seen (false);
			vector<vertex> vset;
			for (size_t j = 0; j < vs.size(); j++) {
				int u = get<0>((*ax.tl)[vs[j]].triple);
				if (seen.hasDefaultValue(u)) {
					seen[u] = true;
					vset.push_back (u);
				}
				u = get<1>((*ax.tl)[vs[j]].triple);
				if (seen.hasDefaultValue(u)) {
					seen[u] = true;
					vset.push_back (u);
				}
				u = get<2>((*ax.tl)[vs[j]].triple);
				if (seen.hasDefaultValue(u)) {
					seen[u] = true;
					vset.push_back (u);
				}
			}
			sizes[co] = vset.size();
			timestamp t1;
			sort (vset.begin(), vset.end());
			ed = report_density (vset, graph, ismaximal, curK, fp, ffp);
			timestamp t2;
			cdtimes[4] += t2 - t1;
		}

		densities[co] = ed;
		timestamp t34;
		cdtimes[5] += t34 - t33;

		// contribute to all_kids of parents
		if (!parent.hasDefaultValue (co)) {
			int pr = parent[co];
			if (variant == 12)
				sizes[pr] += sizes[co]; // todo: probably not needed if you do 'sizes[co] = vs.size()' along the sort above
			all_kids.emplace (pr, co); // put the current nucleus to all_kids of its parent
			// put all_kids of current nucleus to all_kids of its parent
			for (int cc : current_kids)
				all_kids.emplace (pr, cc);
		}
		timestamp t35;
		cdtimes[6] += t35 - t34;
	}
	timestamp t4;
	sizes[id] = graph.size();
//	chmod (a, S_IRUSR);
	fclose(fp);
	// build the nested circles for visualization
	string temp (vfile);
	string cfl = temp + "_circle.json";
	char *bb = new char[cfl.length() + 1];
	strcpy (bb, cfl.c_str());
	FILE* gip = fopen (bb, "w");
	nested_circle (id, new_kids, K_comps, parent, densities, sizes, gip, cfl);
//	chmod (bb, S_IRUSR);
	fclose(gip);
	timestamp t5;
	cdtimes[7] += t5 - t4;
}



void report_subgraphs (int variant, vector<subcore>& hrc, vector<int>& subx_ids, Graph& graph, int nEdge, p_auxies& ax, const char* vfile, FILE* ffp) {

#ifdef SHORT // just the hierarchy and K values
	timestamp p1;
	rearrange_children (hrc);
//	queue<vertex> bfsorder; // we are doing bfs
//	bfsorder.push(hrc.size()-1);
//	int ind = 0;
//	while (!bfsorder.empty()) {
//		int s = bfsorder.front();
//		bfsorder.pop();
//		fprintf (ffp, "%d: K: %d pr: %d\n", s, hrc[s].K, hrc[s].parent);
//		fflush (ffp);
//		printf ("%d: K: %d pr: %d\n", s, hrc[s].K, hrc[s].parent);
//		for (int r : hrc[s].children)
//			bfsorder.push (r);
//	}
//
//	timestamp p2;
//	cout << "printing time: " << p2 - p1 << endl;
//	print_time (ffp, "print time: ", p2 - p1);

#ifdef STAT
	// stats for reference
	timestamp s1;
	int sum = 0, max = 0, sel;
	for (int i = 0; i < hrc.size(); i++) {
		int init = i;
		int c = 0;

		while (hrc[i].root != -1) {
			i = hrc[i].root;
			c++;
		}
		sum += c;
		if (c > max) {
			max = c;
			sel = init;
		}
	}
	fprintf (ffp, "max-depth: %d, avg-depth: %lf\n", max, sum / (double) hrc.size());
#endif

	vector<int> tebaaNums (hrc.size(), -1);
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].visible) {

			queue<int> order;
			order.push (i);
			int sum = 0;
			while (!order.empty()) {
				int co = order.front();
				order.pop();

				int k = hrc[co].K;
				for (int j = 0; j < hrc[co].children.size(); j++) {
					int ch = hrc[co].children[j];
					if (hrc[ch].K == k)
						order.push (ch);
					else if (hrc[ch].K > k)
						sum++;
				}
			}
			tebaaNums[i] = sum;
		}
	}

	HashMap<int> dist (0);
	int avg = 0, max = 0, min = INT_MAX, n = 0;
	for (size_t i = 0; i < tebaaNums.size(); i++) {
		if (tebaaNums[i] == -1)
			continue;
		dist[tebaaNums[i]]++;
		if (tebaaNums[i] > max)
			max = tebaaNums[i];
		if (tebaaNums[i] < min)
			min = tebaaNums[i];
		avg += tebaaNums[i];
		n++;
	}
	printf ("Avg-branching: %lf  max: %d  min: %d\n", (double (avg) / n), max, min);

	for (auto iter = dist.begin(); iter != dist.end(); iter++)
		printf ("%d with %d children\n", iter->second, iter->first);


#else
	// set subx_ids of each vertex/edge to the id of its representative subcore (the top subcore with same consecutive K)
	timestamp q1;

	// assign unassigned items to top subcore
	int nd = hrc.size() - 1;
//	int nd = 0;
	for (int i = 0; i < subx_ids.size(); i++) {
		if (subx_ids[i] == -1)
			subx_ids[i] = nd;
	}

	// find the representative core, e.g., nearest visible ancestor
	HashMap<int> replace (-1);
	for (int i = 0; i < hrc.size(); i++) {
		int sc = i;
		int init = sc;
		if (hrc[sc].parent != -1) {
			int pr = hrc[sc].parent;
			while (hrc[sc].K == hrc[pr].K) {
				sc = pr;
				if (hrc[sc].parent != -1)
					pr = hrc[sc].parent;
				else
					break;
			}
		}
		if (init != sc)
			hrc[init].visible = false;
		replace[init] = sc;
	}

	// and get its id
	for (int i = 0; i < subx_ids.size(); i++) {
		if (!replace.hasDefaultValue (subx_ids[i]))
			subx_ids[i] = replace[subx_ids[i]];
	}

	rearrange_children (hrc);

	timestamp q2;
	print_time (ffp, "visibility setting and subxid assign time: ", q2 - q1);
	cout << "visibility setting and subxid assign time: " << q2 - q1 << endl;

	timestamp ph1 (0, 0), ph2 (0, 0), ph3 (0, 0), ph4 (0, 0);
	string tmp (vfile);
	string tfl = tmp + "_fakedb";
	char *a = new char[tfl.length() + 1];
	strcpy (a, tfl.c_str());
	FILE* fp = fopen (a, "w+");
	HashMap<bool> marked (false);
	vector<subcore> cpc (hrc);
	HashMap<int> ordermap (-1); // key is the hrc index, value is the order
	int order = 0; // order of subcores in file

	while (1) {
		bool flag = false;
		for (size_t ind = 0; ind < cpc.size(); ind++) {
			int i = ind;
			if (marked.hasDefaultValue(i) && cpc[i].visible && cpc[i].children.empty()) {
				flag = true;
				marked[i] = true;
				ordermap[i] = (order)++;

				// finds the items whose subx_id is i, appends the vertices of those items to the vlists combination of its children, sort, find the number of edges
				report_results (variant, i, ordermap, subx_ids, ax, hrc, graph, fp, ph1, ph2, ph3, ph4); // LEAVES

				if (cpc[i].parent != -1) {
					int pr = cpc[i].parent;
					for (int j = 0; j < cpc[pr].children.size(); j++)
						if (cpc[pr].children[j] == i) {
							cpc[pr].children.erase (cpc[pr].children.begin() + j);
							break;
						}
				}
			}
		}
		if (!flag)
			break;
	}



#ifndef SHORT
	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].visible) {
			printf ("%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, hrc[i].K, hrc[i].size, hrc[i].nedge, hrc[i].ed, hrc[i].parent, hrc[i].children.empty()?"LEAF":"INTERMEDIATE");
			fprintf (ffp, "%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, hrc[i].K, hrc[i].size, hrc[i].nedge, hrc[i].ed, hrc[i].parent, hrc[i].children.empty()?"LEAF":"INTERMEDIATE");
			fflush (ffp);
		}
	}
#endif

	timestamp q3;
	print_time (ffp, "printing (with densities) time: ", q3 - q2);
	cout << "printing (with densities) time: " << q3 - q2 << endl;

	// build the nested circles for visualization
	string temp (vfile);
	string cfl = temp + "_circle.json";
	char *bb = new char[cfl.length() + 1];
	strcpy (bb, cfl.c_str());
	FILE* gip = fopen (bb, "w");
	print_nested_circle (hrc, hrc.size() - 1, gip, cfl);
//	chmod (aa, S_IRUSR);
//	chmod (bb, S_IRUSR);
	fclose(gip);
	timestamp q4;
	print_time (ffp, "circle generation time: ", q4 - q3);
	cout << "circle generation time: " << q4 - q3 << endl;
#endif
}

#ifdef LCPS
void report_LCPS (int variant, Graph& graph, int cn, p_auxies& ax, HashMap<bool>& exists, vector<vertex>& K, const char* vfile, FILE* ffp) {
	timestamp q1;
	size_t sz;
	sz = graph.size();

	double prev_rt = 0;

	int max_K = 0;
	for (size_t i = 0; i < graph.size(); i++) {
		if (K[i] > max_K)
			max_K = K[i];
	}

	vector<subcore> hrc;
	vector<vertex> subx_ids (graph.size(), -1);
	int u = -19; // current vertex id
	int cr = -1; // current hrc id

	// LCPS-1
	int i = 0;
	int k = -1;
	int p;
	Max_Bucket bs;
	bs.Initialize(max_K, (graph.size()));

	vector<bool> visited (graph.size(), false);
	int cnt = 0;
	while (1) {
		cnt++;
		// LCPS-2
		if (i == graph.size())
			break;
		i++;

		if (u == -19) { // LCPS-3
			i--;
			k = -1;
			int j, myj;
			int maxk = 0;
			for (j = 0; j < subx_ids.size(); j++) {
				if (subx_ids[j] == -1 && K[j] > maxk) {
					maxk = K[j];
					myj = j;
				}
			}
			u = myj;
			if (visited[u] == true)
				break;
			visited[u] = true;
//			printf ("1- u: %d\n", u);
			for (int l = 0; l < k+1; l++)
				cr = hrc[cr].parent;

			for (int l = 0; l < K[u]; l++) {
				subcore sc;
				sc.K = l + 1;
				sc.parent = cr;
				hrc.push_back (sc);
				cr = hrc.size() - 1;
			}

			subx_ids[u] = cr;
			k = K[u];

			for (int jj = 0; jj < graph[u].size(); jj++) {
				if (subx_ids[graph[u][jj]] == -1) {
					int v = graph[u][jj];
					if (bs.CurrentValue(v) == -1) // not in the bucket
						bs.Insert (v, K[v]);
				}
			}

			int ret = bs.PopMax(&u, &p); // LCPS-4
			if (u != -19) {
				visited[u] = true;
//				printf ("2- u: %d\n", u);
			}
//			if (ret == -1 || p == 0)
//				break;
		}
		else {
			// LCPS-5
			for (int l = 0; l < k - p; l++)
				cr = hrc[cr].parent;

			for (int l = 0; l < p - k; l++) {
				subcore sc;
				sc.K = k + l + 1;
				sc.parent = cr;
				hrc.push_back (sc);
				cr = hrc.size() - 1;
			}

			subx_ids[u] = cr;
			k = K[u];

			// LCPS-6
			for (int jj = 0; jj < graph[u].size(); jj++) {
				if (subx_ids[graph[u][jj]] == -1) {
					int v = graph[u][jj];
					if (bs.CurrentValue(v) == -1) // not in the bucket
						bs.Insert (v, K[v]);
				}
			}

			int ret = bs.PopMax(&u, &p); // LCPS-4
			if (u != -19) {
				visited[u] = true;
//				printf ("3- u: %d\n", u);
			}
//			if (ret == -1 || p == 0)
//				break;
		}

	}

	int nEdge = 0;
	for (size_t i = 0; i < graph.size(); i++)
		nEdge += graph[i].size();

	// rearrange children and parents based on visibility
	for (size_t i = 0; i < hrc.size(); i++)
		hrc[i].children.clear();

	for (size_t i = 0; i < hrc.size(); i++) {
		int pr = hrc[i].parent;
		if (pr != -1)
			hrc[pr].children.push_back (i);
	}

	timestamp q2;
	int rr = 0;
	for (int i = 0; i < graph.size(); i++)
		if (graph[i].size() == 0)
			rr++;
	printf ("cnt: %d, graph.size: %d, zerodegrees: %d\n", cnt, graph.size(), rr);
	cout << "lcps traverse: " << q2 - q1 << endl;
//	for (size_t i = 0; i < graph.size(); i++) {
//		if (visited[i] == false)
//			printf ("%d is not visited, K: %d, deg: %d\n", i, K[i], graph[i].size());
//	}

#ifndef SHORT
	timestamp ph1 (0, 0), ph2 (0, 0), ph3 (0, 0), ph4 (0, 0);
	string tmp (vfile);
	string tfl = tmp + "_fakedb";
	char *a = new char[tfl.length() + 1];
	strcpy (a, tfl.c_str());
	FILE* fp = fopen (a, "w+");
	HashMap<bool> marked (false);
	vector<subcore> cpc (hrc);
	HashMap<int> ordermap (-1); // key is the hrc index, value is the order
	int order = 0; // order of subcores in file

	while (1) {
		bool flag = false;
		for (size_t ind = 0; ind < cpc.size(); ind++) {
			int i = ind;
			if (marked.hasDefaultValue(i) && cpc[i].children.empty()) {
				flag = true;
				marked[i] = true;
				ordermap[i] = order++;
				// finds the items whose subx_id is i, appends the vertices of those items to the vlists combination of its children, sort, find the number of edges
				report_results (variant, i, ordermap, subx_ids, ax, hrc, graph, fp, ph1, ph2, ph3, ph4);
				if (cpc[i].parent != -1) {
					int pr = cpc[i].parent;
					for (int j = 0; j < cpc[pr].children.size(); j++)
						if (cpc[pr].children[j] == i) {
							cpc[pr].children.erase (cpc[pr].children.begin() + j);
							break;
						}
				}
			}
		}
		if (!flag)
			break;
	}

	for (size_t i = 0; i < hrc.size(); i++) {
		if (hrc[i].visible) {
			printf ("%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, hrc[i].K, hrc[i].size, hrc[i].nedge, hrc[i].ed, hrc[i].parent, hrc[i].children.empty()?"LEAF":"INTERMEDIATE");
			fprintf (ffp, "%d: K: %d |V|: %d |E|: %d ed: %.5lf pr: %d (%s)\n", i, hrc[i].K, hrc[i].size, hrc[i].nedge, hrc[i].ed, hrc[i].parent, hrc[i].children.empty()?"LEAF":"INTERMEDIATE");
			fflush (ffp);
		}
	}

	timestamp q3;
#ifdef NODENSITY
	print_time (ffp, "printing (without densities) time: ", q3 - q2);
	cout << "printing (without densities) time: " << q3 - q2 << endl;
#else
	print_time (ffp, "printing (with densities) time: ", q3 - q2);
	cout << "printing (with densities) time: " << q3 - q2 << endl;
#endif

#ifndef NODENSITY
	// build the nested circles for visualization
	string temp (vfile);
	string cfl = temp + "_circle.json";
	char *bb = new char[cfl.length() + 1];
	strcpy (bb, cfl.c_str());
	FILE* gip = fopen (bb, "w");
	print_nested_circle (hrc, hrc.size() - 1, gip, cfl);
//	chmod (aa, S_IRUSR);
//	chmod (bb, S_IRUSR);
	fclose(gip);
	timestamp q4;
	print_time (ffp, "circle generation time: ", q4 - q3);
	cout << "circle generation time: " << q4 - q3 << endl;
#endif

#endif

}
#endif

