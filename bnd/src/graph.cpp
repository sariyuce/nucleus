#include "main.h"

#define MAXLINE 1000000
//#define WRITE_BINARY

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


// Regular graphs

void writeWeightedRegularBinary (string filename, vertex nVtx, edge nEdge, Wraph& wraph) {

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

void writeRegularBinary (string filename, vertex nVtx, edge nEdge, vector<vector<vertex>>& graph) {

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

void readWeightedBinary(char *filename, Wraph& wraph, edge* nEdge) {

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
}




// Bipartite graphs

template <typename VtxType, typename EdgeType>
void readBipartiteBinary(char *filename, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph, EdgeType* nEdge) {

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
	*nEdge = nn;

	leftGraph.resize (leftVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * leftVtx);
	really_read (in, (char*)pxadj, sizeof(EdgeType) * leftVtx);
	for (vertex i = 0; i < leftVtx; i++) {
		leftGraph[i].resize (pxadj[i]);
		really_read (in, (char*)&(leftGraph[i][0]), sizeof(VtxType) * pxadj[i]);
		if (i % 1000000 == 0)
			printf ("left i: %d / %d\n", i, leftVtx);
	}

	rightGraph.resize (rightVtx);
	EdgeType *pxadj2 = (EdgeType*) malloc (sizeof(EdgeType) * rightVtx);
	really_read(in, (char*)pxadj2, sizeof(EdgeType) * rightVtx);
	for (vertex i = 0; i < rightVtx; i++) {
		rightGraph[i].resize (pxadj2[i]);
		really_read (in, (char*)&(rightGraph[i][0]), sizeof(VtxType) * pxadj2[i]);
		if (i % 1000000 == 0)
			printf ("right i: %d / %d\n", i, rightVtx);
	}

	int ss = 0;
	for (vertex i = 0; i < rightVtx; i++) {
		ss += rightGraph[i].size();
	}
	*nEdge = ss;
}

template <typename VtxType, typename EdgeType>
void writeBipartiteBinary (char* filename, EdgeType nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

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
	}

	for (VtxType i = 0; i < rightGraph.size(); i++) {
		size_t sz = rightGraph[i].size();
		fwrite (&(rightGraph[i][0]), sizeof(VtxType), sz, filep);
		if (i % 1000000 == 0)
			printf ("right i: %d\n", i);
	}

	fclose (filep);
}

template <typename VtxType, typename EdgeType>
void readBipartiteChaco (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

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

//	writeBipartiteBinary (filename, nEdge, leftGraph, rightGraph);
}

template <typename VtxType, typename EdgeType>
void readBipartiteMM (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

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

//	writeBipartiteBinary (filename, nEdge, leftGraph, rightGraph);
}

template <typename VtxType, typename EdgeType>
void readBipartiteOut (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

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
	while (fgets (line, MAXLINE, fp) != NULL) {
		stringstream ss (line);
		ss >> u >> v;
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
		if (co % 1000000 == 0) {
			timestamp t2;
			cout << "read " << co << "  edges, time: " << t2 - t1 << endl;

		}
	}
	for (size_t i = 0; i < leftGraph.size(); i++)
		sort (leftGraph[i].begin(), leftGraph[i].end());

	for (size_t i = 0; i < rightGraph.size(); i++)
		sort (rightGraph[i].begin(), rightGraph[i].end());

	fclose (fp);

}

template <typename VtxType, typename EdgeType>
void readBipartite (char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	string st (filename);
	int idx = st.find_last_of(".");
	string ext = st.substr(idx);

	if (ext == ".bin") // Binary format
		readBipartiteBinary<VtxType, VtxType> (filename, leftGraph, rightGraph, nEdge);
	else if (ext == ".graph") // Chaco format
		readBipartiteChaco<VtxType> (filename, nEdge, leftGraph, rightGraph);
	else if (ext == ".amazon") { // For the user product data in Amazon ratings, each line has user id and product id (http://jmcauley.ucsd.edu/data/amazon/)
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
	else if (st.find("out") == 0) { // SNAP format
		readBipartiteOut<VtxType> (filename, nEdge, leftGraph, rightGraph);
	}
	else // MatrixMarket format
		readBipartiteMM<VtxType> (filename, nEdge, leftGraph, rightGraph);


#ifdef WRITE_BINARY
	if (ext != ".bin")
		writeBipartiteBinary (filename, nEdge, leftGraph, rightGraph);
#endif

	printf ("|Left|: %d		|Right|: %d		|Edge|: %d\n", leftGraph.size(), rightGraph.size(), *nEdge);
}

template void readBipartite (char *filename, edge* nEdge, Graph& leftGraph, Graph& rightGraph);

