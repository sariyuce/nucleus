#include "main.h"

#define MAXLINE 1000000

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
void ReadBipartiteGraph(char *filename, EdgeType* nEdge, vector<vector<VtxType>>& leftGraph, vector<vector<VtxType>>& rightGraph) {

	string st (filename);
	int idx = st.find_last_of(".");
	string ext = st.substr(idx);

	if (ext == ".graph") // Chaco
		ReadBipartiteGraphFromChacoFile<VtxType> (filename, nEdge, leftGraph, rightGraph);
	else // MatrixMarket
		ReadBipartiteGraphFromMMFile<VtxType> (filename, nEdge, leftGraph, rightGraph);
}

template void ReadBipartiteGraph (char *filename, edge* nEdge, Graph& leftGraph, Graph& rightGraph);
