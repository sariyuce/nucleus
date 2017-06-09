#include "main.h"

#define MAXLINE 1000000

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

	printf ("asdsdg\n");
	int co = 0;
	int u, v, dum, sum;
	unordered_map<pair<int, int>, int> mp;
	int maxU = 0, maxV = 0;
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

	printf ("asdsdg\n");
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
	else if (st.find("out") != string::npos) {
//		printf ("wiki read\n");
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



template void ReadBipartiteGraph (char *filename, edge* nEdge, Graph& leftGraph, Graph& rightGraph);















