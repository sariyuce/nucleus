#include "main.h"

typedef struct pair {
	int f;
	int s;
} Pair;

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

int pcmp(const void *v1, const void *v2) {
	vertex diff = (((Pair *)v1)->f - ((Pair *)v2)->f);
	if (diff != 0)
		return diff;
	else
		return (((Pair *)v1)->s - ((Pair *)v2)->s);
}

template <typename VtxType, typename EdgeType>
void writeToBinary (char* filename, VtxType nVtx, EdgeType nEdge, vector<vector<VtxType>>& graph) {

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
void ReadBinary(char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {

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
	vertex nVtx;
	in.read((char*)&nVtx, sizeof(VtxType));
	in.read((char*)nEdge, sizeof(EdgeType));

	graph.resize (nVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * nVtx);
	really_read(in, (char*)pxadj, sizeof(EdgeType) * nVtx);
	for (vertex i = 0; i < nVtx; i++) {
		graph[i].resize (pxadj[i]);
		really_read (in, (char*)&(graph[i][0]), sizeof(VtxType) * pxadj[i]);
	}

//	debug
//	cout << nVtx << " " << *nEdge << endl;
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
	Pair* coords = (Pair*) malloc (sizeof(Pair) * 2 * *nEdge);
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
	qsort(coords, index, sizeof(Pair), pcmp);

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

template <typename VtxType, typename EdgeType>
void ReadGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {

	string tmp (filename);
	string filetype = tmp.substr (tmp.find_last_of("."));

	if (filetype == ".graph")
		ReadGraphFromChacoFile<VtxType, EdgeType> (filename, graph, nEdge);
	else
	if (filetype == ".bin")
		ReadBinary<VtxType, EdgeType> (filename, graph, nEdge);
	else // .mtx or .txt
		ReadGraphFromMMFile<VtxType, EdgeType> (filename, graph, nEdge);

	return;
}

template void ReadGraph (char *filename, vector<vector<vertex>>& graph, edge* nEdge);

