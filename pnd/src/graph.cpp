//#ifndef _GRAPH_H_
//#define _GRAPH_H_
//
//#include "ulib.h"
//#include <string>
//#include <sstream>
//#include <iostream>
//#include <fstream>
//#include <map>
//#include <vector>

#include "main.h"
typedef struct pair {
	int f;
	int s;
} Pair;

#define MAXLINE 128*(1024*1024)

using namespace std;

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
	long long diff = (((Pair *)v1)->f - ((Pair *)v2)->f);
	if (diff != 0)
		return diff;
	else
		return (((Pair *)v1)->s - ((Pair *)v2)->s);
}


template <typename VtxType, typename EdgeType>
void SubReadBinary(char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {

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
	VtxType nVtx;
	in.read((char*)&nVtx, sizeof(VtxType));
	in.read((char*)nEdge, sizeof(EdgeType));

	graph.resize (nVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * nVtx);
	really_read(in, (char*)pxadj, sizeof(EdgeType) * nVtx);
	for (VtxType i = 0; i < nVtx; i++) {
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

template <typename VtxType, typename EdgeType, typename WeightType>
void ReadBinary(char *filename, VtxType *numofvertex_r, VtxType *numofvertex_c, EdgeType **pxadj, VtxType **padj, WeightType **padjw, WeightType **ppvw) {

	vector<vector<VtxType>> graph;
	EdgeType nEdge;
	SubReadBinary<VtxType, EdgeType> (filename, graph, &nEdge);

	*pxadj = (EdgeType*) malloc (sizeof(EdgeType) * (graph.size() + 1));
	*padj =  (VtxType*) malloc (sizeof(VtxType) * (nEdge));

	(*pxadj)[0] = 0;
	int c  = 0;
	for (int i = 0; i < graph.size(); i++) {
		(*pxadj)[i+1] = graph[i].size();
		for (const auto w: graph[i])
			(*padj)[c++] = w;
	}
	*numofvertex_c = *numofvertex_r = graph.size();


//	if (ppvw != NULL) {
//		cerr<<"vertex weight is unsupported"<<std::endl;
//		return;
//	}
//
//	std::ifstream in (filename);
//
//	if (!in.is_open()) {
//		cerr<<"can not open file:"<<filename<<std::endl;
//		return;
//	}
//
//	int vtxsize; //in bytes
//	int edgesize; //in bytes
//	int weightsize; //in bytes
//
//	//reading header
//	in.read((char *)&vtxsize, sizeof(int));
//	in.read((char *)&edgesize, sizeof(int));
//	in.read((char *)&weightsize, sizeof(int));
//
//	if (!in) {
//		cerr<<"IOError"<<std::endl;
//		return;
//	}
//
//	if (vtxsize != sizeof(VtxType)) {
//		cerr<<"Incompatible VertexSize."<<endl;
//		return;
//	}
//
//	if (edgesize != sizeof(EdgeType)) {
//		cerr<<"Incompatible EdgeSize."<<endl;
//		return;
//	}
////
////	if (weightsize != sizeof(WeightType)) {
////		cerr<<"Incompatible WeightType."<<endl;
////		return;
////	}
//
//	//reading should be fine from now on.
//	in.read((char*)numofvertex_r, sizeof(VtxType));
////	in.read((char*)numofvertex_c, sizeof(VtxType));
//	EdgeType nnz;
//	in.read((char*)&nnz, sizeof(EdgeType));
////	if ((int) numofvertex_c <=0 || (int) numofvertex_r <=0 || nnz <= 0) {
////		cerr<<"graph makes no sense"<<endl;
////		return;
////	}
//
//	*pxadj = (EdgeType*) malloc (sizeof(EdgeType) * (*numofvertex_r+1));
//	*padj =  (VtxType*) malloc (sizeof(VtxType) * (nnz));
//
//
//	if (padjw) {
//		*padjw = new WeightType[nnz];
//	}
//
//	*pxadj[0] = 0;
//	int err = really_read(in, (char*)&(*pxadj[1]), sizeof(EdgeType)*(*numofvertex_r));
//
//	for (VtxType i = 0; i < *numofvertex_r; i++) {
//		int x = *pxadj[i];
//		really_read(in, (char*) &(padj[x]), sizeof(VtxType)*(*pxadj[i+1] - *pxadj[i]));
////		fwrite (&(entire_graph[i][0]), sizeof(VtxType), sz, filep);
//	}
////	err += really_read(in, (char*)*padj, sizeof(VtxType)*(nnz));
////	if (padjw)
////		err += really_read(in, (char*)*padjw, sizeof(WeightType)*(nnz));
//	if (!in || err != 0) {
//		cerr<<"IOError"<<endl;
//	}
//
//	return;
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

	printf ("read:   nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
	graph.resize (nVtx);
	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * nVtx);
	really_read(in, (char*)pxadj, sizeof(EdgeType) * nVtx);
	printf ("read:   nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
	for (vertex i = 0; i < nVtx; i++) {
		graph[i].resize (pxadj[i]);

		really_read (in, (char*)&(graph[i][0]), sizeof(VtxType) * pxadj[i]);
//		printf ("i %d: %d\n", i, pxadj[i]);
	}
}









template <typename VtxType, typename EdgeType, typename WeightType>
void ReadGraphFromChacoFile(char *filename, VtxType *numofvertex, EdgeType **pxadj, VtxType **padjncy, WeightType **padjncyw, WeightType **ppvw) {

	EdgeType *xadj, nedges;
	VtxType *adjncy, nvtxs, edge;
	int fmt, readew, readvw;
	EdgeType i, k, ncon;
	WeightType *adjncyw=NULL, *pvw=NULL;
	char *line;

	line = (char *) malloc(sizeof(char) * (MAXLINE+1));

	FILE* fpin = fopen(filename, "r");
	// skip comments
	do {
		fgets(line, MAXLINE, fpin);
	} while (line[0] == '%' && !feof(fpin));

	if (feof(fpin))
		errexit("empty graph!!!");

	// read header
	fmt = 0;
	{
		std::string s = line;
		std::stringstream ss (s);
		ss>>nvtxs>>nedges>>fmt>>ncon;
	}
	*numofvertex = nvtxs;

	// check if it is weighted
	readew = (fmt % 10 > 0);
	readvw = ((fmt / 10) % 10 > 0);
	if (fmt >= 100)
		errexit("invalid format");

	// allocations
	nedges *=2;

	xadj = *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * (nvtxs+2));
	adjncy = *padjncy = (VtxType*) malloc (sizeof(VtxType) * (nedges));
	if (padjncyw)
		adjncyw = *padjncyw = (VtxType*) malloc (sizeof(VtxType) * (nedges));
	if (ppvw)
		pvw = *ppvw = (WeightType*) malloc (sizeof(WeightType) * (nedges));

	// read each edge list
	for (xadj[0] = 0, k = 0, i = 0; i < nvtxs; i++) {
		char *oldstr=line, *newstr=NULL;
		int  ewgt=1, vw=1;

		do {
			fgets(line, MAXLINE, fpin);
		} while (line[0] == '%' && !feof(fpin));

		if (strlen(line) >= MAXLINE-5)
			errexit("\nBuffer for fgets not big enough!\n");

		if (readvw) {
			vw = (WeightType)strtol(oldstr, &newstr, 10);
			oldstr = newstr;
		}

		if (ppvw)
			pvw[i] = vw;

		for (;;) {
			edge = (VtxType)strtol(oldstr, &newstr, 10) -1;
			oldstr = newstr;

			if (readew) {
				ewgt = (WeightType)strtol(oldstr, &newstr, 10);
				oldstr = newstr;
			}

			if (edge < 0) {
				break;
			}

			if (edge == i)
				errexit("Self loop in the graph for vertex %d\n", i);

			bool flag = false;
			for (EdgeType j = xadj[i]; j < k; j++) {
				if (adjncy[j] == edge) {
					flag = true;
					break;
				}
			}

			if (!flag) {
				adjncy[k] = edge;
				if (padjncyw)
					adjncyw[k] = ewgt;
				k++;
			}
		}
		xadj[i+1] = k;
	}

	free(line);
	return;
}

/* reads the Matrix Market format graph, does mapping from VtxType numbers to VtxType numbers from 0 to nVtx */
template <typename VtxType, typename EdgeType, typename WeightType>
void ReadGraphFromMMFile(char *filename, VtxType *numofvertex, EdgeType **pxadj, VtxType **padjncy, WeightType **padjncyw, WeightType **ppvw) {
	Pair *coords, *new_coords;
	EdgeType nnz, tnnz, onnz;
	VtxType m, n;
	VtxType *adj, *adjncyw;
	VtxType itemp, jtemp;
	EdgeType *xadj;
	WeightType *pvw;

	int maxLineLength = 1000000;
	char line[maxLineLength];

	/* set return null parameter values, in case we exit with errors */
	m = nnz = 0;

	FILE* matfp = fopen(filename, "r");

	/* now continue scanning until you reach the end-of-comments */
	do {
		fgets(line, 1000000, matfp);
	} while (line[0] == '%');


	/* line[] is either blank or has M,N, nz */
	{
		std::string s = line;
		std::stringstream ss (s);
		ss>>n>>nnz>>nnz;
	}

	coords = (Pair*) malloc(sizeof(Pair) * 2 * nnz);

	// begin symmetrization

	// nnz is given # of edges
	tnnz = 0;
	VtxType max = 0;
	for(EdgeType i = 0; i < nnz; i++) {
		fgets(line, 1000000, matfp);
		std::string s = line;
		std::stringstream ss (s);
		ss>>itemp>>jtemp;
		if(itemp != jtemp) {
			coords[tnnz].f = itemp;
			coords[tnnz++].s = jtemp;
			coords[tnnz].f = jtemp;
			coords[tnnz++].s = itemp;
		}

		if (itemp > max)
			max = itemp;
		if (jtemp > max)
			max = jtemp;
//		if (i % 1000000 == 0)
//			printf("nnz: %d, i: %d\n", nnz, i);
	}

	m = n = max + 1;
	// tnnz is symmetrized # of edges
	qsort(coords, tnnz, sizeof(Pair), pcmp);

	onnz = 1;
	for(EdgeType i = 1; i < tnnz; i++) {
		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
			coords[onnz].f = coords[i].f;
			coords[onnz++].s = coords[i].s;
		}
	}

	// end symmetrization

	// allocations (onnz is duplicate-eliminated # of edges)
	*numofvertex = n;
//	printf("numofvertex: %d\n", *numofvertex);
	xadj = *pxadj = (EdgeType *) malloc((n+1) * sizeof(EdgeType));
	adj = *padjncy = (VtxType *) malloc(onnz * sizeof(VtxType));
	if (padjncyw)
		adjncyw = *padjncyw = imalloc (nnz, "ReadGraph: adjncyw");
	if (ppvw)
		pvw = *ppvw = imalloc (n+1, "ReadGraph: adjncyw");



	// begin mapping
	map<VtxType, VtxType> reed;
	typename map<VtxType, VtxType>::iterator reed_it; //Todo: CHANGE THE TYPES TO VtxType EXPLICITLY, C++ doesn't allow to say VtxType
	new_coords = (Pair*) malloc(sizeof(Pair) * 2 * nnz);
	VtxType vno = 0;


	//    // map the ids
	//    for(EdgeType i = 0; i < onnz; i++) {
	//    	VtxType temp = coords[i].f;
	//    	reed_it = reed.find(temp);
	//    	if (reed_it == reed.end()) {
	//    	    reed.insert (make_pair (temp, vno));
	//    	    new_coords[i].f = vno++;
	//    	}
	//    	else
	//    		new_coords[i].f = reed_it->second;
	//
	//    	temp = coords[i].s;
	//    	reed_it = reed.find(temp);
	//    	if (reed_it == reed.end()) {
	//    	    reed.insert (make_pair (temp, vno));
	//    	    new_coords[i].s = vno++;
	//    	}
	//    	else
	//    		new_coords[i].s = reed_it->second;
	//
	//		if (i % 1000000 == 0)
	//			printf("in map: onnz: %ld, i: %ld\n", onnz, i);
	//    }
	//    // end mapping




	// begin constructing graph
	vector<vector<VtxType> > entire_graph;
	entire_graph.resize(n);
	for(EdgeType i = 0; i < onnz; i++) {
		//    	entire_graph[new_coords[i].f].push_back(new_coords[i].s);
		entire_graph[coords[i].f].push_back(coords[i].s);
	}
//
//
//	xadj[0] = 0;
//	EdgeType j = 0;
//	for(VtxType i = 1; i < n+1; i++) {
//		xadj[i] = xadj[i-1] + entire_graph[i-1].size();
//		for (VtxType k = 0; k < entire_graph[i-1].size(); k++) {
//			adj[j++] = entire_graph[i-1][k];
//		}
//	}
//
//	free(coords);
//	for(VtxType i = 0; i < m; i++)
//		entire_graph[i].clear();
//
//	entire_graph.clear();
	free (new_coords);
	fclose (matfp);



	writeBinary (filename, *numofvertex, onnz, entire_graph);

	printf ("binary written\n");
	exit(1);

	// write binary
//	char* str = (char*) malloc (sizeof(char) * 50);
//	str = filename;
//	strcat (str, ".bin");
//	FILE* filep = fopen (str, "w");
//	int vtxt = sizeof (VtxType);
//	int edget = sizeof (EdgeType);
//	int wt = sizeof (WeightType);
//	fwrite (&vtxt, sizeof(int), 1, filep);
//	fwrite (&edget, sizeof(int), 1, filep);
////	fwrite (&wt, sizeof(int), 1, filep);
//
//	VtxType nvtx = *numofvertex;
//	EdgeType nedge = j;
////	fwrite (&nvtx, sizeof(VtxType), 1, filep);
//	fwrite (&nvtx, sizeof(VtxType), 1, filep);
//	fwrite (&nedge, sizeof(EdgeType), 1, filep);
//
//	for (VtxType i = 0; i < nvtx; i++) {
//		VtxType sz = entire_graph[i].size();
//		fwrite (&sz, sizeof(VtxType), 1, filep);
//	}
//
//	for (VtxType i = 0; i < nvtx; i++) {
//		size_t sz = entire_graph[i].size();
//		fwrite (&(entire_graph[i][0]), sizeof(VtxType), sz, filep);
//	}
//
////	fwrite (xadj, sizeof(EdgeType), n+1, filep);
////	fwrite (adj, sizeof(VtxType), j, filep);
//	fclose (filep);
////
//	exit(1);
	return;
}

template <typename VtxType, typename EdgeType, typename WeightType>
void ReadGraph (Graph& graph, EdgeType* nEdge, char *filename, VtxType *numofvertex, EdgeType **pxadj, VtxType **padjncy, WeightType **padjncyw, WeightType **ppvw) {

	const char * pch;
	//    pch = strstr (filename,".");
	string tmp (filename);
	int idx = tmp.find_last_of(".");
	string asdf = tmp.substr(idx);
	pch = asdf.c_str();
	//	cout << "pch: " << pch << endl;
	char t1[10];
	char t2[10];
	char t3[10];
	strcpy (t1, ".graph");
	strcpy (t2, ".mtx");
	strcpy (t3, ".bin");

	if (strcmp (pch, t1) == 0)
		ReadGraphFromChacoFile<VtxType, EdgeType, VtxType> (filename, numofvertex, pxadj, padjncy, padjncyw, ppvw);
	else if (strcmp (pch, t3) == 0) {
		readBinary<VtxType, EdgeType> (filename, graph, nEdge);
//		ReadBinary<VtxType, EdgeType, VtxType> (filename, numofvertex, numofvertex, pxadj, padjncy, NULL, NULL);

		printf ("|V|: %d   |E|: %d\n", graph.size(), *nEdge);
		(*pxadj) = (VtxType *) malloc (sizeof(vertex) * (graph.size() + 1));
		(*padjncy) = (VtxType *) malloc (sizeof(vertex) * (*nEdge));

		vertex c = 0;
		(*pxadj)[0] = 0;
		for (vertex i = 0; i < graph.size(); ++i) {
			(*pxadj)[i+1] = (*pxadj)[i] + graph[i].size();
			sort (graph[i].begin(), graph[i].end());
			for (vertex w : graph[i])
				(*padjncy)[c++] = w;
		}
		*numofvertex = graph.size();
	}
	else// if (strcmp (pch, t2) == 0)
		ReadGraphFromMMFile<VtxType, EdgeType, VtxType> (filename, numofvertex, pxadj, padjncy, padjncyw, ppvw);



//	VtxType* adj = *padjncy;
//	EdgeType* xadj = *pxadj;
//	VtxType n = *numofvertex;
	//    int* tadj = *ptadj;

	return;
}

template void ReadGraph (Graph& graph, edge* nEdge, char *filename, vertex *numofvertex, edge **pxadj, vertex **padjncy, vertex **padjncyw, vertex **ppvw);


//#endif
//void ReadGraph(char *filename, VtxType *numofvertex, EdgeType **pxadj, VtxType **padjncy, WeightType **padjncyw, WeightType **ppvw);



//#include "main.h"
//
//#define MAXLINE 1000000
//
//typedef struct asdf {
//	int f;
//	int s;
//} pprr;
//
//inline bool hashUniquify (vector<vertex>& vertices) {
//	HashMap<bool> hermap (false);
//	for (size_t i = 0; i < vertices.size(); i++) {
//		int t = vertices[i];
//		if (hermap.hasDefaultValue (t))
//			hermap[t] = true;
//		else {
//			vertices.erase (vertices.begin() + i);
//			i--;
//		}
////		if (i > UPPERBOUND)
////			return false;
//	}
//	sort (vertices.begin(), vertices.end());
//	return true;
//}
//
//int pcmp(const void *v1, const void *v2) {
//	vertex diff = (((pprr *)v1)->f - ((pprr *)v2)->f);
//	if (diff != 0)
//		return diff;
//	else
//		return (((pprr *)v1)->s - ((pprr *)v2)->s);
//}
//
//static int really_read(std::istream& is, char* buf, size_t global_size) {
//	char* temp2 = buf;
//	while (global_size != 0) {
//		is.read(temp2, global_size);
//		size_t s = is.gcount();
//		if (!is)
//			return -1;
//
//		global_size -= s;
//		temp2 += s;
//	}
//	return 0;
//}
//
//template <typename VtxType, typename EdgeType>
//void writeBinary (char* filename, VtxType nVtx, EdgeType nEdge, vector<vector<VtxType>>& graph) {
//
//	string str(filename);
//	string fl = str + ".bin";
//	FILE* filep = fopen (fl.c_str(), "w");
//	int vtxt = sizeof (VtxType);
//	int edget = sizeof (EdgeType);
//	fwrite (&vtxt, sizeof(int), 1, filep);
//	fwrite (&edget, sizeof(int), 1, filep);
//
//	fwrite (&nVtx, sizeof(VtxType), 1, filep);
//	fwrite (&nEdge, sizeof(EdgeType), 1, filep);
//
//	for (VtxType i = 0; i < nVtx; i++) {
//		VtxType sz = graph[i].size();
//		fwrite (&sz, sizeof(VtxType), 1, filep);
//	}
//
//	for (VtxType i = 0; i < nVtx; i++) {
//		size_t sz = graph[i].size();
//		fwrite (&(graph[i][0]), sizeof(VtxType), sz, filep);
//	}
//
//	fclose (filep);
//}
//
//template <typename VtxType, typename EdgeType>
//void readBinary(char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {
//
//	ifstream in (filename);
//	int vtxsize; //in bytes
//	int edgesize; //in bytes
//
//	//reading header
//	in.read((char *)&vtxsize, sizeof(int));
//	in.read((char *)&edgesize, sizeof(int));
//
//	if (!in) {
//		cerr<<"IOError"<<std::endl;
//		return;
//	}
//
//	if (vtxsize != sizeof(VtxType)) {
//		cerr<<"Incompatible VertexSize."<<endl;
//		return;
//	}
//
//	if (edgesize != sizeof(EdgeType)) {
//		cerr<<"Incompatible EdgeSize."<<endl;
//		return;
//	}
//
//	vertex nVtx;
//	in.read((char*)&nVtx, sizeof(VtxType));
//	in.read((char*)nEdge, sizeof(EdgeType));
//
//	printf ("nVtx: %d   nEdge:%d\n", nVtx, *nEdge);
//	graph.resize (nVtx);
//	EdgeType *pxadj = (EdgeType*) malloc (sizeof(EdgeType) * nVtx);
//	really_read(in, (char*)pxadj, sizeof(EdgeType) * nVtx);
//	for (vertex i = 0; i < nVtx; i++) {
//		graph[i].resize (pxadj[i]);
//		really_read (in, (char*)&(graph[i][0]), sizeof(VtxType) * pxadj[i]);
//	}
//}
//
//template <typename VtxType, typename EdgeType>
//void readChaco (char *filename, Graph& graph, EdgeType* nEdge) {
//
//	char* line = (char*) malloc (sizeof (char) * MAXLINE);
//	FILE* matfp = fopen(filename, "r");
//
//	// skip comments
//	do {
//		fgets(line, MAXLINE, matfp);
//	} while (line[0] == '%');
//
//	VtxType nVtx, neig;
//	string s = line;
//	stringstream ss (s);
//	ss >> nVtx >> *nEdge;
//
//	graph.resize (nVtx);
//	// read each edge list
//	for (VtxType i = 0; i < nVtx; i++) {
//		fgets(line, MAXLINE, matfp);
//		stringstream ss (line);
//		while (ss >> neig)
//			graph[i].push_back (neig);
//	}
//
//	// sort each neighbor list
//	for(VtxType i = 0; i < nVtx; i++)
//		hashUniquify (graph[i]);
//
//	fclose (matfp);
////	writeBinary (filename, nVtx, *nEdge, graph);
//}
//
//template <typename VtxType, typename EdgeType>
//void readMM (char *filename, Graph& graph, EdgeType* nEdge) {
//
//	char* line = (char*) malloc (sizeof (char) * MAXLINE);
//	FILE* matfp = fopen(filename, "r");
//
//	// skip comments
//	do {
//		fgets(line, 1000000, matfp);
//	} while (line[0] == '%');
//
//
//	VtxType nVtx;
//	stringstream ss (line);
//	ss >> nVtx >> *nEdge;
//
//	// remove duplicate edges, take one direction
//	pprr* coords = (pprr*) malloc (sizeof(pprr) * 2 * *nEdge);
//	VtxType itemp, jtemp, index = 0;
//
//
//	for (EdgeType i = 0; i < *nEdge; i++) {
//		fgets(line, MAXLINE, matfp);
//		stringstream ss (line);
//		ss >> itemp >> jtemp;
//		if(itemp != jtemp) {
//			coords[index].f = coords[index + 1].s = itemp;
//			coords[index + 1].f = coords[index].s = jtemp;
//			index += 2;
//		}
//	}
//
//	// onnz is # of edges
//	qsort(coords, index, sizeof(pprr), pcmp);
//
//	VtxType onnz = 1;
//	for(EdgeType i = 1; i < index; i++) {
//		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
//			coords[onnz].f = coords[i].f;
//			coords[onnz++].s = coords[i].s;
//		}
//	}
//
//	// begin constructing graph
//	graph.resize (nVtx);
//	for(EdgeType i = 0; i < onnz; i++)
//		graph[coords[i].f].push_back(coords[i].s);
//
//	// sort each neighbor list
//	edge numedge = 0;
//	for(VtxType i = 0; i < nVtx; i++) {
//		sort (graph[i].begin(), graph[i].end());
//		numedge += graph[i].size();
//	}
//	*nEdge = numedge;
//
//	fclose (matfp);
//
////	writeBinary (filename, nVtx, *nEdge, graph);
//}
//
//template <typename VtxType, typename EdgeType>
//void readOut (char *filename, Graph& graph, EdgeType* nEdge) {
//
//	char* line = (char*) malloc (sizeof (char) * MAXLINE);
//	FILE* matfp = fopen(filename, "r");
//
//	// skip comments
//	do {
//		fgets(line, 1000000, matfp);
//	} while (line[0] == '%');
//
//	VtxType u, v;
//	string dum, sum;
//	unordered_map<pair<int, int>, bool> mp;
//	VtxType maxV = 0;
//	while (fgets (line, MAXLINE, matfp) != NULL) {
//		stringstream ss (line);
//		ss >> u >> v;
//		if (u > v)
//			swap (u, v);
//		auto at = mp.find (make_pair (u, v));
//		if (at == mp.end()) {
//			mp.emplace (make_pair(make_pair(u, v), true));
//
//			if (v > maxV) {
//				graph.resize (v+1);
//				maxV = v;
//			}
//			graph[u].push_back (v);
//			graph[v].push_back (u);
//		}
//	}
//
//	*nEdge = 0;
//	for (auto w : graph) {
//		*nEdge += w.size();
//		sort (w.begin(), w.end());
//	}
//	fclose (matfp);
//
////	writeBinary (filename, nVtx, *nEdge, graph);
//}
//
//template <typename VtxType, typename EdgeType>
//void readGraph (char *filename, vector<vector<VtxType>>& graph, EdgeType* nEdge) {
//
//	string tmp (filename);
//	string filetype = tmp.substr (tmp.find_last_of("."));
//
//	if (filetype == ".bin")
//		readBinary<VtxType, EdgeType> (filename, graph, nEdge);
//	else if (filetype == ".graph")
//		readChaco<VtxType, EdgeType> (filename, graph, nEdge);
//	else if (tmp.find("out") != string::npos)
//		readOut<VtxType, EdgeType> (filename, graph, nEdge);
//	else // .mtx or .txt
//		readMM<VtxType, EdgeType> (filename, graph, nEdge);
//
//	return;
//}
//
//template void readGraph (char *filename, vector<vector<vertex>>& graph, edge* nEdge);
//


