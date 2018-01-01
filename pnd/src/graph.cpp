#include "main.h"

#define MAXLINE 1000000
#define WRITE_BINARY

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
void writeBinary (char* filename, VtxType nVtx, EdgeType nEdge, VtxType* adj, EdgeType* xadj) {

	string str(filename);
	string fl = str + ".bin";
	FILE* filep = fopen (fl.c_str(), "w");
	int vtxt = sizeof (VtxType);
	int edget = sizeof (EdgeType);
	fwrite (&vtxt, sizeof(int), 1, filep);
	fwrite (&edget, sizeof(int), 1, filep);

	fwrite (&nVtx, sizeof(VtxType), 1, filep);
	fwrite (&nEdge, sizeof(EdgeType), 1, filep);

	for (VtxType i = 0; i < nVtx; i++) {
		EdgeType sz = xadj[i+1] - xadj[i];
		fwrite (&sz, sizeof(EdgeType), 1, filep);
	}

	fwrite (adj, sizeof(VtxType), xadj[nVtx], filep);

	fclose (filep);
}

template <typename VtxType, typename EdgeType>
void readBinary(char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {

	ifstream in (filename);
	int vtxsize; //in bytes
	int edgesize; //in bytes

	//reading header
	in.read((char *)&vtxsize, sizeof(int));
	in.read((char *)&edgesize, sizeof(int));

	if (!in) {
		cerr << "IOError" << endl;
		return;
	}

	if (vtxsize != sizeof(VtxType)) {
		cerr << "Incompatible VertexSize" << endl;
		return;
	}

	if (edgesize != sizeof(EdgeType)) {
		cerr << "Incompatible EdgeSize" << endl;
		return;
	}

	in.read((char*)nVtx, sizeof(VtxType)); // we already write this as +1
	in.read((char*)nEdge, sizeof(EdgeType));

	(*xadj) = (EdgeType *) malloc (sizeof(EdgeType) * (*nVtx + 1));
	(*adj) = (VtxType *) malloc (sizeof(VtxType) * (*nEdge * 2));

	(*xadj)[0] = 0;
	for (VtxType i = 0; i < *nVtx; i++) {
		EdgeType nt = 0;
		really_read (in, (char*) &nt, sizeof(EdgeType));
		(*xadj)[i+1] = (*xadj)[i] + nt;
	}
	really_read (in, (char*)(*adj), sizeof(VtxType) * (*xadj)[*nVtx]);
}

template <typename VtxType, typename EdgeType>
void VV2CRS (vector<vector<VtxType>>& graph, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {
	EdgeType ei = 0;
	*xadj = (EdgeType *) malloc (sizeof(EdgeType) * (*nVtx + 1));
	*adj = (VtxType *) malloc (sizeof(VtxType) * (*nEdge * 2));
	(*xadj)[0] = 0;
	for (VtxType i = 0; i < *nVtx; i++) {
		(*xadj)[i+1] = graph[i].size() + (*xadj)[i];
		for (VtxType j = 0; j < graph[i].size(); j++)
			(*adj)[ei++] = graph[i][j];
		graph[i].clear();
	}
	graph.clear();
}

template <typename VtxType, typename EdgeType>
void readChaco (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, MAXLINE, matfp);
	} while (line[0] == '%' || line[0] == '#');

	VtxType neig;

	string s = line;
	stringstream ss (s);
	ss >> *nVtx >> *nEdge;

	*nVtx += 1; // Since our graphs are zero-based

	vector<vector<VtxType>> graph (*nVtx);

	// read each edge list
	for (VtxType i = 0; i < *nVtx; i++) {
		fgets(line, MAXLINE, matfp);
		stringstream ss (line);
		while (ss >> neig)
			graph[i].push_back (neig);
	}
	fclose (matfp);

	// sort each neighbor list
	for(VtxType i = 0; i < *nVtx; i++)
		hashUniquify (graph[i]);

	VV2CRS (graph, nVtx, nEdge, adj, xadj);
}

template <typename VtxType, typename EdgeType>
void readMM (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%' || line[0] == '#');

	stringstream ss (line);
//	int a;
//	ss >> a;
	ss >> *nVtx >> *nEdge;


	printf ("|V|: %d   |E|: %d\n", *nVtx, *nEdge);

	*nVtx += 1; // Since our graphs are zero-based

	printf ("|V|: %d   |E|: %d\n", *nVtx, *nEdge);
	// remove duplicate edges, take one direction
	pprr* coords = (pprr*) malloc (sizeof(pprr) * 2 * *nEdge);
	VtxType itemp, jtemp, index = 0;

	for (EdgeType i = 0; i < *nEdge; i++) {
		fgets(line, MAXLINE, matfp);
		stringstream ss (line);
		ss >> itemp >> jtemp;
		if(itemp != jtemp) {
			coords[index].f = coords[index + 1].s = itemp;
			coords[index + 1].f = coords[index].s = jtemp;
			index += 2;
		}
	}
	fclose (matfp);

	qsort(coords, index, sizeof(pprr), pcmp);

	VtxType onnz = 1; // onnz is # of edges
	for(EdgeType i = 1; i < index; i++) {
		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
			coords[onnz].f = coords[i].f;
			coords[onnz++].s = coords[i].s;
		}
	}

	// begin constructing graph
	vector<vector<VtxType>> graph (*nVtx);
	for(EdgeType i = 0; i < onnz; i++)
		graph[coords[i].f].push_back(coords[i].s);

	// sort each neighbor list
	EdgeType numedge = 0;
	for(VtxType i = 0; i < *nVtx; i++) {
		sort (graph[i].begin(), graph[i].end());
		numedge += graph[i].size();
	}

	if (numedge != *nEdge * 2)
		printf ("nEdge in header is wrong: %d (must be %d)\n", *nEdge, numedge/2);

	VV2CRS (graph, nVtx, nEdge, adj, xadj);
}

template <typename VtxType, typename EdgeType>
void readOut (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%' || line[0] == '#');

	VtxType u, v;
	string dum, sum;
	unordered_map<pair<VtxType, VtxType>, bool> mp;
	VtxType maxV = 0;
	vector<vector<VtxType>> graph;
	while (fgets (line, MAXLINE, matfp) != NULL) {
		stringstream ss (line);
		ss >> u >> v;
		if (u > v)
			swap (u, v);
		auto at = mp.find (make_pair (u, v));
		if (at == mp.end()) {
			mp.emplace (make_pair(make_pair(u, v), true));
			if (v > maxV) {
				graph.resize (v+1);
				maxV = v;
			}
			graph[u].push_back (v);
			graph[v].push_back (u);
			(*nEdge)++;
		}
	}
	fclose (matfp);

	for (auto w : graph)
		sort (w.begin(), w.end());

	*nVtx = graph.size();
	VV2CRS (graph, nVtx, nEdge, adj, xadj);
}

template <typename VtxType, typename EdgeType>
void readGraph (char *filename, VtxType* nVtx, EdgeType* nEdge, VtxType** adj, EdgeType** xadj) {

	string st (filename);
	string gname = st.substr (st.find_last_of("/") + 1);
	int idx = gname.find_last_of(".");
	string ext = gname.substr(idx);

	if (ext == ".bin")
		readBinary<VtxType, EdgeType> (filename, nVtx, nEdge, adj, xadj);
	else if (ext == ".graph")
		readChaco<VtxType, EdgeType> (filename, nVtx, nEdge, adj, xadj);
	else if (gname.find("out") == 0) {
		printf ("it's out\n");
		readOut<VtxType, EdgeType> (filename, nVtx, nEdge, adj, xadj);
	}
	else // .mtx or .txt
		readMM<VtxType, EdgeType> (filename, nVtx, nEdge, adj, xadj);

#ifdef WRITE_BINARY
	if (ext != ".bin") {
		writeBinary (filename, *nVtx, *nEdge, *adj, *xadj);
		printf ("Binary graph is written\n");
	}
#endif

	return;
}

template void readGraph (char *filename, vertex* nVtx, edge* nEdge, vertex** adj, edge** xadj);

