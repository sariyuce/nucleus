#include "main.h"

#define MAXLINE 1000000
#define WRITE_BINARY

typedef struct asdf {
	int f;
	int s;
} pprr;

int pcmp(const void *v1, const void *v2) {
	vertex diff = (((pprr *)v1)->f - ((pprr *)v2)->f);
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
void writeBinary (char* filename, VtxType nVtx, EdgeType nEdge, vector<vector<VtxType>>& graph) {

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
		VtxType sz = graph[i].size();
		fwrite (&sz, sizeof(VtxType), 1, filep);
	}

	for (VtxType i = 0; i < nVtx; i++) {
		size_t sz = graph[i].size();
		fwrite (&(graph[i][0]), sizeof(VtxType), sz, filep);
	}

	fclose (filep);
}

template <typename VtxType, typename EdgeType>
void readBinary(char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {

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

	vertex nVtx;
	in.read((char*)&nVtx, sizeof(VtxType));
	in.read((char*)nEdge, sizeof(EdgeType));

	printf ("nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
	graph.resize (nVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * nVtx);
	really_read(in, (char*)pxadj, sizeof(EdgeType) * nVtx);
	for (vertex i = 0; i < nVtx; i++) {
		graph[i].resize (pxadj[i]);
		really_read (in, (char*)&(graph[i][0]), sizeof(VtxType) * pxadj[i]);
	}
}

template <typename VtxType, typename EdgeType>
void readChaco (char *filename, Graph& graph, EdgeType* nEdge) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
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
		stringstream ss (line);
		while (ss >> neig)
			graph[i].push_back (neig);
	}

	// sort each neighbor list
	for(VtxType i = 0; i < nVtx; i++)
		hashUniquify (graph[i]);

	fclose (matfp);
}

template <typename VtxType, typename EdgeType>
void readMM (char *filename, Graph& graph, EdgeType* nEdge) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%');

	VtxType nVtx;
	stringstream ss (line);
	ss >> nVtx >> *nEdge;

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

	qsort(coords, index, sizeof(pprr), pcmp);

	VtxType onnz = 1; // onnz is # of edges
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

	fclose (matfp);

}

template <typename VtxType, typename EdgeType>
void readOut (char *filename, Graph& graph, EdgeType* nEdge) {

	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%');

	VtxType u, v;
	string dum, sum;
	unordered_map<pair<int, int>, bool> mp;
	VtxType maxV = 0;
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
	for (auto w : graph)
		sort (w.begin(), w.end());
	fclose (matfp);

}

template <typename VtxType, typename EdgeType>
void readGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {

	string tmp (filename);
	string filetype = tmp.substr (tmp.find_last_of("."));

	if (filetype == ".bin")
		readBinary<VtxType, EdgeType> (filename, graph, nEdge);
	else if (filetype == ".graph")
		readChaco<VtxType, EdgeType> (filename, graph, nEdge);
	else if (tmp.find("out") == 0)
		readOut<VtxType, EdgeType> (filename, graph, nEdge);
	else // .mtx or .txt
		readMM<VtxType, EdgeType> (filename, graph, nEdge);

#ifdef WRITE_BINARY
	if (filetype != ".bin") {
		vertex nVtx = graph.size();
		writeBinary (filename, nVtx, *nEdge, graph);
		printf ("Binary graph is written\n");
	}
#endif

	return;
}

template void readGraph (char *filename, vector<vector<vertex>>& graph, edge* nEdge);

