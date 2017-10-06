#include "main.h"

#define MAXLINE 1000000

typedef struct asdf {
	int f;
	int s;
} pprr;

int pcmp(const void *v1, const void *v2) {
        long long diff = (((pprr *)v1)->f - ((pprr *)v2)->f);
        if (diff != 0)
                return diff;
        else
                return (((pprr *)v1)->s - ((pprr *)v2)->s);
}

static int really_read(std::istream& is, char* buf, size_t global_size) {
	char* temp2 = buf;
	while (global_size != 0) {
		is.read(temp2, global_size);
		size_t s = is.gcount();
		if (!is)
			return -1;

		global_size -= s;
		temp2 += s;
	}
	return 0;
}

template <typename VtxType, typename EdgeType>
void ReadBinary(char *filename, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph, EdgeType* nEdge) {

	timestamp t1;
	ifstream in (filename);
	int vtxsize; //in bytes
	int edgesize; //in bytes

	//reading header
	in.read((char *)&vtxsize, sizeof(int));
	in.read((char *)&edgesize, sizeof(int));

	if (!in) {
			cerr<<"IOError"<<std::endl;
			return;
	}

	if (vtxsize != sizeof(VtxType)) {
			cerr<<"Incompatible VertexSize."<<endl;
			return;
	}

	if (edgesize != sizeof(EdgeType)) {
			cerr<<"Incompatible EdgeSize."<<endl;
			return;
	}

	//reading should be fine from now on.
	vertex leftVtx, rightVtx;
	edge nn;
	really_read (in, (char *)&leftVtx, sizeof(VtxType));
	really_read(in, (char*)&rightVtx, sizeof(VtxType));
	really_read(in, (char*)&nn, sizeof(EdgeType));
	printf ("nn: %d\n", nn);
	*nEdge = nn;

	leftGraph.resize (leftVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * leftVtx);
	really_read (in, (char*)pxadj, sizeof(EdgeType) * leftVtx);
	for (vertex i = 0; i < leftVtx; i++) {
		leftGraph[i].resize (pxadj[i]);
		really_read (in, (char*)&(leftGraph[i][0]), sizeof(VtxType) * pxadj[i]);
		if (i % 1000000 == 0)
			printf ("left i: %d\n", i);
	}

	rightGraph.resize (rightVtx);
	EdgeType *pxadj2 = (EdgeType*) malloc (sizeof(EdgeType) * rightVtx);
	really_read(in, (char*)pxadj2, sizeof(EdgeType) * rightVtx);
	for (vertex i = 0; i < rightVtx; i++) {
		rightGraph[i].resize (pxadj2[i]);
		really_read (in, (char*)&(rightGraph[i][0]), sizeof(VtxType) * pxadj2[i]);
		if (i % 1000000 == 0)
			printf ("right i: %d\n", i);
	}
	timestamp t2;
	cout << "time: " << t2 - t1 << endl;

	int ss = 0;
	for (vertex i = 0; i < rightVtx; i++) {
		ss += rightGraph[i].size();
	}
	*nEdge = ss;
	printf ("Left.sz: %d   Right.sz: %d\n", leftGraph.size(), rightGraph.size());
	int ls = 0;
	for (vertex i = 0; i < leftVtx; i++) {
		ls += leftGraph[i].size();
	}

	printf ("rE: %d   lE: %d\n", ss, ls);



//	for (vertex i = 0; i < rightGraph.size(); i++) {
//		for (vertex j = 0; j < rightGraph[i].size(); j++) {
//			int u = i;
//			int v = rightGraph[i][j];
//			bool f = false;
//			for (vertex k = 0; k < leftGraph[v].size(); k++) {
//				if (leftGraph[v][k] == u) {
//					f = true;
//					break;
//				}
//			}
//			if (!f) {
//				printf ("@1 HOOOOOOO   %d-%d is there but %d-%d is NOT\n", u, v, v, u);
//				exit(1);
//			}
//		}
//		if (i % 1000 == 0)
//			printf ("%d/%d\n", i, rightGraph.size());
//	}
//
//	printf ("passed @1\n");
//	for (vertex i = 0; i < leftGraph.size(); i++) {
//		for (vertex j = 0; j < leftGraph[i].size(); j++) {
//			int u = i;
//			int v = leftGraph[i][j];
//			bool f = false;
//			for (vertex k = 0; k < rightGraph[v].size(); k++) {
//				if (rightGraph[v][k] == u) {
//					f = true;
//					break;
//				}
//			}
//			if (!f) {
//				printf ("@2 HOOOOOOO   %d-%d is there but %d-%d is NOT\n", u, v, v, u);
//				exit(1);
//			}
//		}
//		if (i % 1000 == 0)
//			printf ("%d/%d\n", i, leftGraph.size());
//	}
//
//	printf ("passed @2\n");
//
//


	return;
}


template <typename VtxType, typename EdgeType>
void writeToBinary (char* filename, EdgeType nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	string str(filename);
	string fl = str + ".bin";
	FILE* filep = fopen (fl.c_str(), "w");
	int vtxt = sizeof (VtxType);
	int edget = sizeof (EdgeType);
	fwrite (&vtxt, sizeof(int), 1, filep);
	fwrite (&edget, sizeof(int), 1, filep);

	int leftVtx = leftGraph.size();
	int rightVtx = rightGraph.size();
	fwrite (&leftVtx, sizeof(VtxType), 1, filep);
	fwrite (&rightVtx, sizeof(VtxType), 1, filep);
	fwrite (&nEdge, sizeof(EdgeType), 1, filep);

	leftGraph.resize (leftVtx);
	for (VtxType i = 0; i < leftGraph.size(); i++) {
		VtxType sz = leftGraph[i].size();
		fwrite (&sz, sizeof(VtxType), 1, filep);
		if (i % 1000000 == 0)
			printf ("left sz i: %d\n", i);
	}

	for (VtxType i = 0; i < leftGraph.size(); i++) {
		size_t sz = leftGraph[i].size();
		fwrite (&(leftGraph[i][0]), sizeof(VtxType), sz, filep);
		if (i % 1000000 == 0)
			printf ("left i: %d\n", i);
	}

	rightGraph.resize (rightVtx);
	for (VtxType i = 0; i < rightGraph.size(); i++) {
		VtxType sz = rightGraph[i].size();
		fwrite (&sz, sizeof(VtxType), 1, filep);
		if (i % 1000000 == 0)
			printf ("right sz i: %d\n", i);
	}

	for (VtxType i = 0; i < rightGraph.size(); i++) {
		size_t sz = rightGraph[i].size();
		fwrite (&(rightGraph[i][0]), sizeof(VtxType), sz, filep);
		if (i % 1000000 == 0)
			printf ("right i: %d\n", i);
	}

	fclose (filep);
}

void ReadRegularBinary(char *filename, vector<vector<vertex>>& graph, edge* nEdge) {

	ifstream in (filename);
	int vtxsize; //in bytes
	int edgesize; //in bytes

	//reading header
	in.read((char *)&vtxsize, sizeof(int));
	in.read((char *)&edgesize, sizeof(int));

	if (!in) {
		cerr<<"IOError"<<std::endl;
		return;
	}

	if (vtxsize != sizeof(vertex)) {
		cerr<<"Incompatible VertexSize."<<endl;
		return;
	}

	if (edgesize != sizeof(edge)) {
		cerr<<"Incompatible EdgeSize."<<endl;
		return;
	}

	//reading should be fine from now on.
	vertex nVtx;
	in.read((char*)&nVtx, sizeof(vertex));
	in.read((char*)nEdge, sizeof(edge));

	printf ("nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
	graph.resize (nVtx);
	edge *pxadj = (edge*) malloc (sizeof(edge) * nVtx);
	really_read(in, (char*)pxadj, sizeof(edge) * nVtx);
	for (vertex i = 0; i < nVtx; i++) {
		graph[i].resize (pxadj[i]);
		really_read (in, (char*)&(graph[i][0]), sizeof(vertex) * pxadj[i]);
	}

//	debug
	cout << nVtx << " " << *nEdge << endl;
	int ne = 0;
	for (int i = 0; i < nVtx; i++) {
		ne += graph[i].size();
	}

	cout << "ne : " << ne << endl;
//		for (int j = 0; j < graph[i].size(); j++)
//			cout << graph[i][j] << " ";
//		cout << endl;
//	}

	return;
}

void ReadWeightedBinary(char *filename, Wraph& wraph, edge* nEdge) {

	ifstream in (filename);
	int vtxsize; //in bytes
	int edgesize; //in bytes

	//reading header
	in.read((char *)&vtxsize, sizeof(int));
	in.read((char *)&edgesize, sizeof(int));

	if (!in) {
		cerr<<"IOError"<<std::endl;
		return;
	}

	if (vtxsize != sizeof(vertex)) {
		cerr<<"Incompatible VertexSize."<<endl;
		return;
	}

	if (edgesize != sizeof(edge)) {
		cerr<<"Incompatible EdgeSize."<<endl;
		return;
	}

	//reading should be fine from now on.
	vertex nVtx;
	in.read((char*)&nVtx, sizeof(vertex));
	in.read((char*)nEdge, sizeof(edge));

	printf ("nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
	wraph.resize (nVtx);
	edge *pxadj = (edge*) malloc (sizeof(edge) * nVtx);
	really_read(in, (char*)pxadj, sizeof(edge) * nVtx);
	for (vertex i = 0; i < nVtx; i++) {
		wraph[i].resize (pxadj[i]);
		really_read (in, (char*)&(wraph[i][0]), sizeof(wv) * pxadj[i]);
	}

//	debug
	cout << nVtx << " " << *nEdge << endl;
//	for (int i = 0; i < nVtx; i++) {
//		for (int j = 0; j < graph[i].size(); j++)
//			cout << graph[i][j] << " ";
//		cout << endl;
//	}

	return;
}

template <typename VtxType, typename EdgeType>
void ReadGraphFromChacoFile(char *filename, Graph& graph, EdgeType* nEdge) {

//	printf ("I'll only create a binary version, ok? [Just Enter]\n");
//	char ch;
//	scanf ("%c", &ch);
//	if (ch != '\n')
//		exit (1);

	char line[MAXLINE];
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, MAXLINE, matfp);
	} while (line[0] == '%');

	VtxType nVtx, neig;
	string s = line;
	stringstream ss (s);
	ss >> nVtx >> *nEdge;

	graph.resize (nVtx);
	// read each edge list
	for (VtxType i = 0; i < nVtx; i++) {
		fgets(line, MAXLINE, matfp);
		string s = line;
		stringstream ss (s);
		while (ss >> neig)
			graph[i].push_back (neig);
	}

	// sort each neighbor list
	for(VtxType i = 0; i < nVtx; i++)
		sort (graph[i].begin(), graph[i].end());

	cout << nVtx << " " << *nEdge << endl;
	int ne = 0;
	for (int i = 0; i < nVtx; i++) {
		ne += graph[i].size();
	}

	cout << "ne : " << ne << endl;
	// write binary
//	writeToBinary (filename, nVtx, *nEdge, graph);

	return;
}

/* reads the Matrix Market format graph */
template <typename VtxType, typename EdgeType>
void ReadGraphFromMMFile (char *filename, Graph& graph, EdgeType* nEdge) {

//	printf ("I'll only create a binary version, ok? [Just Enter]\n");
//	char ch;
//	scanf ("%c", &ch);
//	if (ch != '\n')
//		exit (1);

	VtxType nVtx;
	char line[1000000];
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%');

	string s (line);
	stringstream ss (s);
	ss >> nVtx >> *nEdge;

	// remove duplicate edges, take one direction
	pprr* coords = (pprr*) malloc (sizeof(pprr) * 2 * *nEdge);
	VtxType itemp, jtemp, index = 0;


	for (EdgeType i = 0; i < *nEdge; i++) {
		fgets(line, 1000000, matfp);
		string s (line);
		stringstream ss (s);
		ss >> itemp >> jtemp;
		if(itemp != jtemp) {
			coords[index].f = coords[index + 1].s = itemp;
			coords[index + 1].f = coords[index].s = jtemp;
			index += 2;
		}
	}

	// onnz is # of edges
	qsort(coords, index, sizeof(pprr), pcmp);

	VtxType onnz = 1;
	for(EdgeType i = 1; i < index; i++) {
		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
			coords[onnz].f = coords[i].f;
			coords[onnz++].s = coords[i].s;
		}
	}

	// begin constructing graph
	graph.resize (nVtx);
	for(EdgeType i = 0; i < onnz; i++)
		graph[coords[i].f].push_back(coords[i].s);

	// sort each neighbor list
	edge numedge = 0;
	for(VtxType i = 0; i < nVtx; i++) {
		sort (graph[i].begin(), graph[i].end());
		numedge += graph[i].size();
	}
	*nEdge = numedge;

	// write binary
//	writeToBinary (filename, nVtx, *nEdge, graph);

	return;
}

// reads the Chcao format bipartite graph, vertices are the ones on the left
template <typename VtxType, typename EdgeType>
void ReadBipartiteGraphFromChacoFile (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	char line[MAXLINE];
	FILE* fp = fopen(filename, "r");

	// skip comments
	do {
		fgets (line, MAXLINE, fp);
	} while (line[0] == '%');

	VtxType leftnVtx, rightnVtx;
	stringstream ss (line);
	ss >> *nEdge >> leftnVtx >> rightnVtx;

	leftGraph.resize (leftnVtx);
	rightGraph.resize (rightnVtx);
	VtxType u, v;

	// read each edge list
	for (VtxType u = 0; u < leftnVtx; u++) {
		fgets (line, MAXLINE, fp);
		stringstream ss (line);
		while (ss >> v) {
			leftGraph[u].push_back (v);
			rightGraph[v].push_back (u);
		}
	}

	for (size_t i = 0; i < leftGraph.size(); i++)
		hashUniquify (leftGraph[i]);

	for (size_t i = 0; i < rightGraph.size(); i++)
		hashUniquify (rightGraph[i]);

	fclose (fp);
}

// reads the Matrix Market format bipartite graph, file can have multiple edges, they'll be filtered
template <typename VtxType, typename EdgeType>
void ReadBipartiteGraphFromMMFile (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* fp = fopen(filename, "r");

	// skip comments
	do {
		fgets (line, MAXLINE, fp);
	} while (line[0] == '%');

	VtxType leftnVtx, rightnVtx;
	stringstream ss (line);
	ss >> *nEdge >> leftnVtx >> rightnVtx;

	leftGraph.resize (leftnVtx);
	rightGraph.resize (rightnVtx);
	VtxType u, v;

	for (VtxType i = 0; i < *nEdge; i++) {
		fgets (line, MAXLINE, fp);
		stringstream ss (line);
		ss >> u >> v;
		leftGraph[u].push_back (v);
		rightGraph[v].push_back (u);
	}

	for (size_t i = 0; i < leftGraph.size(); i++)
		hashUniquify (leftGraph[i]);

	for (size_t i = 0; i < rightGraph.size(); i++)
		hashUniquify (rightGraph[i]);

	fclose (fp);
}

template <typename VtxType, typename EdgeType>
void ReadReuters (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	timestamp t1;
	char line[MAXLINE];
	FILE* fp = fopen(filename, "r");

	// skip comments
	do {
		fgets (line, MAXLINE, fp);
	} while (line[0] == '%');

	int co = 0;
	int u, v;
	string dum, sum;
	unordered_map<pair<int, int>, int> mp;
	int maxU = 0, maxV = 0;
//	leftGraph.resize (137695);
//	rightGraph.resize (2255877);
	while (fgets (line, MAXLINE, fp) != NULL) {
		stringstream ss (line);
		ss >> u >> v >> dum >> sum;

//		printf ("u: %d, v: %d\n", u, v);
		auto at = mp.find (make_pair(u, v));
		if (at == mp.end()) {
			mp.emplace (make_pair(make_pair(u, v), 1));

			if (u > maxU) {
				leftGraph.resize (u+1);
				maxU = u;
			}

			if (v > maxV) {
				rightGraph.resize (v+1);
				maxV = v;
			}
			leftGraph[u].push_back (v);
			rightGraph[v].push_back (u);
			(*nEdge)++;
		}

		co++;
//		if (co == 10) {
//			printf ("maxU: %d  maxV: %d asdf\n", maxU, maxV);
//			exit(1);
//		}
		if (co % 1000000 == 0) {
			timestamp t2;
			cout << "co: " << co << "  time: " << t2 - t1 << endl;

		}
	}

	printf ("maxU: %d  maxV: %d\n", maxU, maxV);
//	exit(1);
	for (size_t i = 0; i < leftGraph.size(); i++) {

//		if (i % 1000 == 0) {
//			timestamp t2;
//			cout << "left sort: " << i << "  size: "<< leftGraph[i].size() << " time: " << t2 - t1 << endl;
//		}
		sort (leftGraph[i].begin(), leftGraph[i].end());
	}

	for (size_t i = 0; i < rightGraph.size(); i++) {

//		if (i % 1000 == 0) {
//			timestamp t2;
//			cout << "right sort: " << i << "  size: "<< rightGraph[i].size() << " time: " << t2 - t1 << endl;
//		}
		sort (rightGraph[i].begin(), rightGraph[i].end());
	}

}




template <typename VtxType, typename EdgeType>
void ReadDelicious (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	timestamp t1;
	char line[MAXLINE];
	FILE* fp = fopen(filename, "r");

	leftGraph.resize (4511973);
	rightGraph.resize (33777769);
	// skip comments
	do {
		fgets (line, MAXLINE, fp);
	} while (line[0] == '%');

	int co = 0;
	int u, v, dum, sum;
	unordered_map<pair<int, int>, int> mp;
	while (fgets (line, MAXLINE, fp) != NULL) {
		stringstream ss (line);
		ss >> u >> v >> dum >> sum;

		auto at = mp.find (make_pair(u, v));
		if (at == mp.end()) {
			mp.emplace (make_pair(make_pair(u, v), 1));

			leftGraph[u].push_back (v);
			rightGraph[v].push_back (u);
		}

		co++;
		if (co % 1000000 == 0) {
			timestamp t2;
			cout << "co: " << co << "  time: " << t2 - t1 << endl;

		}
	}

	for (size_t i = 0; i < leftGraph.size(); i++) {

		if (i % 1000 == 0) {
			timestamp t2;
			cout << "left sort: " << i << "  size: "<< leftGraph[i].size() << " time: " << t2 - t1 << endl;
		}
			hashUniquify (leftGraph[i]);
	}

	for (size_t i = 0; i < rightGraph.size(); i++) {

		if (i % 1000 == 0) {
			timestamp t2;
			cout << "right sort: " << i << "  size: "<< rightGraph[i].size() << " time: " << t2 - t1 << endl;
		}
			hashUniquify (rightGraph[i]);
	}

}




template <typename VtxType, typename EdgeType>
void ReadBipartiteGraph(char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	string st (filename);
	int idx = st.find_last_of(".");
	string ext = st.substr(idx);

	if (ext == ".bin")
		ReadBinary<VtxType, VtxType> (filename, leftGraph, rightGraph, nEdge);
	else if (ext == ".deli")
		ReadDelicious<VtxType> (filename, nEdge, leftGraph, rightGraph);
	else if (ext == ".reuters")
		ReadReuters<VtxType> (filename, nEdge, leftGraph, rightGraph);
	else if (ext == ".graph") // Chaco
		ReadBipartiteGraphFromChacoFile<VtxType> (filename, nEdge, leftGraph, rightGraph);
	else if (ext == ".amazon") {
		FILE* fp = fopen (st.c_str(), "r");
		FILE* gp = fopen ("user_product.out", "w");
		char user[1000];
		char product[1000];
		unordered_map<string, int> userMap;
		unordered_map<string, int> productMap;
		int userId = 0, productId = 0;

		while (fscanf (fp, "%s %s", user, product) != EOF) {
			string us (user);
			string ps (product);

			if (userMap.find (us) == userMap.end())
				userMap[us] = userId++;

			if (productMap.find (ps) == productMap.end())
				productMap[ps] = productId++;

			fprintf (gp, "%d %d\n", userMap[us], productMap[ps]);
		}
		fclose (fp);
		fclose (gp);

		FILE* up = fopen ("userMap", "w");
		for (auto it = userMap.begin(); it != userMap.end(); it++)
			fprintf (up, "%s %d\n", (it->first).c_str(), it->second);
		fclose (up);

		FILE* pp = fopen ("productMap", "w");
		for (auto it = productMap.begin(); it != productMap.end(); it++)
			fprintf (pp, "%s %d\n", (it->first).c_str(), it->second);
		fclose (pp);
	}
	else if (st.find("out") != string::npos) {
		ReadReuters<VtxType> (filename, nEdge, leftGraph, rightGraph);
	}
	else // MatrixMarket
		ReadBipartiteGraphFromMMFile<VtxType> (filename, nEdge, leftGraph, rightGraph);

//
//
//
//
//
//
//
//	for (int i = 0; i < leftGraph.size(); i++) {
//		printf ("L %d:  ", i);
//		for (int j = 0; j < leftGraph[i].size(); j++)
//			printf ("%d ", leftGraph[i][j]);
//		printf ("\n");
//	}
//
//
//	for (int i = 0; i < rightGraph.size(); i++) {
//		printf ("R %d:  ", i);
//		for (int j = 0; j < rightGraph[i].size(); j++)
//			printf ("%d ", rightGraph[i][j]);
//		printf ("\n");
//	}


//	writeToBinary<vertex, vertex> (filename, *nEdge, leftGraph, rightGraph);
//	printf ("binary is written, leftGraph.size: %d, rightGraph.size(): %d, nEdge: %d\n",
//			leftGraph.size(), rightGraph.size(), *nEdge);
//	exit(1);

}

void ReadWeightedGraph (char *filename, Wraph& wraph, int* nEdges) {

	string fn (filename);
	if (fn.find(".bin") != string::npos) {
		ReadWeightedBinary (filename, wraph, nEdges);
	}
	else {
		ifstream infile (filename);
		string line;
		int a;
		double we;

		int nVtx;

		getline(infile, line);
		stringstream ss (line);
		ss >> nVtx >> *nEdges;
		wraph.resize(nVtx);
		int u = 0;
		while (getline(infile, line))
		{
			stringstream iss(line);
			while (iss >> a) {
				iss >> we;
				wv aa;
				aa.n = a;
				aa.w = we;
				wraph[u].push_back(aa);
			}
			u++;
		}

		infile.close();


	//    for (int i = 0; i < wraph.size(); i++) {
	//    	printf ("%d:   ", i);
	//    	for (int j = 0; j < wraph[i].size(); j++) {
	//    		printf ("%d %lf   ", wraph[i][j].n, wraph[i][j].w);
	//    	}
	//    	printf ("\n");
	//    }
	}
}

void ReadRegularGraph (char *filename, Graph& graph, int* nEdge) {

	string fn (filename);
	if (fn.find(".bin") != string::npos) {
		ReadRegularBinary (filename, graph, nEdge);
	}
	else {
	//	printf ("I'll only create a binary version, ok? [Just Enter]\n");
	//	char ch;
	//	scanf ("%c", &ch);
	//	if (ch != '\n')
	//		exit (1);

		char line[MAXLINE];
		FILE* matfp = fopen(filename, "r");

		// skip comments
		do {
			fgets(line, MAXLINE, matfp);
		} while (line[0] == '%');

		int nVtx, neig;
		string s = line;
		stringstream ss (s);
		ss >> nVtx >> *nEdge;

		graph.resize (nVtx);
		// read each edge list
		for (int i = 0; i < nVtx; i++) {
			fgets(line, MAXLINE, matfp);
			string s = line;
			stringstream ss (s);
			while (ss >> neig)
				graph[i].push_back (neig);
		}

		// sort each neighbor list
		for(int i = 0; i < nVtx; i++)
			sort (graph[i].begin(), graph[i].end());

		cout << nVtx << " " << *nEdge << endl;
		int ne = 0;
		for (int i = 0; i < nVtx; i++) {
			ne += graph[i].size();
		}

		cout << "ne : " << ne << endl;
		// write binary
	//	writeToBinary (filename, nVtx, *nEdge, graph);
	}

	return;
}

void writeToWeightedRegularBinary (string filename, vertex nVtx, edge nEdge, Wraph& wraph) {

	string str = filename;
	string fl = str + ".bin";
	FILE* filep = fopen (fl.c_str(), "w");
	int vtxt = sizeof (vertex);
	int edget = sizeof (edge);
	fwrite (&vtxt, sizeof(int), 1, filep);
	fwrite (&edget, sizeof(int), 1, filep);

	fwrite (&nVtx, sizeof(vertex), 1, filep);
	fwrite (&nEdge, sizeof(edge), 1, filep);

	for (vertex i = 0; i < nVtx; i++) {
		vertex sz = wraph[i].size();
		fwrite (&sz, sizeof(vertex), 1, filep);
	}

	for (vertex i = 0; i < nVtx; i++) {
		size_t sz = wraph[i].size();
		fwrite (&(wraph[i][0]), sizeof(wv), sz, filep);
	}

	fclose (filep);
}


void writeToRegularBinary (string filename, vertex nVtx, edge nEdge, vector<vector<vertex>>& graph) {

	string str = filename;
	string fl = str + ".bin";
	FILE* filep = fopen (fl.c_str(), "w");
	int vtxt = sizeof (vertex);
	int edget = sizeof (edge);
	fwrite (&vtxt, sizeof(int), 1, filep);
	fwrite (&edget, sizeof(int), 1, filep);

	fwrite (&nVtx, sizeof(vertex), 1, filep);
	fwrite (&nEdge, sizeof(edge), 1, filep);

	for (vertex i = 0; i < nVtx; i++) {
		vertex sz = graph[i].size();
		fwrite (&sz, sizeof(vertex), 1, filep);
	}

	for (vertex i = 0; i < nVtx; i++) {
		size_t sz = graph[i].size();
		fwrite (&(graph[i][0]), sizeof(vertex), sz, filep);
	}

	fclose (filep);
}

template void ReadBipartiteGraph (char *filename, edge* nEdge, Graph& leftGraph, Graph& rightGraph);
//void writeToRegularBinary (string filename, vertex nVtx, edge nEdge, Graph& graph);
//template void writeToRegularBinary (string filename, vertex nVtx, edge nEdge, Graph& graph);
//template void writeToWeightedRegularBinary (string filename, vertex nVtx, edge nEdge, Wraph& wraph);

//template void writeToRegularBinary (char* filename, vertex nVtx, edge nEdge, vector<vector<vertex>>& graph);
//template void writeToWeightedRegularBinary (char* filename, vertex nVtx, edge nEdge, vector<wv>& wraph);












