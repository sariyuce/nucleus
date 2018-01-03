#include "main.h"

inline int getTriangleId (vertex u, vertex v, vertex w, edge* xtl, triangle_id* tlist, vertex* ordered_adj, edge* ordered_xadj) {
	edge ind = findInd (u, v, ordered_adj, ordered_xadj);
	for (vertex i = xtl[ind]; i < xtl[ind + 1]; i++)
		if (get<2>(tlist[i].triple) == w)
			return i;
	return -1;
}

void enumTriangles (edge nEdge, vertex* ordered_adj, edge* ordered_xadj, couple* el, vector<triangle_id>& tl, edge* xtl) {
	vertex ind = 0;
	xtl[ind++] = 0;
	for (size_t i = 0; i < nEdge; i++) {
		vector<vertex> is;
		intersection2 (ordered_adj, ordered_xadj, get<0>(el[i]), get<1>(el[i]), is);
		for (size_t j = 0; j < is.size(); j++) {
			triangle_id t;
			t.triple = make_tuple(get<0>(el[i]), get<1>(el[i]), is[j]);
			t.id = 0;
			tl.push_back(t);
		}
		xtl[ind++] = tl.size();
	}
	return;
}

/* STORES TRIANGLE_4CLIQUE GRAPHS */

inline void conditionalIntersection2 (vertex* adj, edge* xadj, vertex u, vertex v, vector<vertex>& intersection) {
	vertex i = xadj[u];
	vertex j = xadj[v];
	vertex gu = xadj[u+1];
	vertex gv = xadj[v+1];

	while (i < gu && j < gv) {
		if (adj[i] < adj[j])
			i++;
		else if (adj[j] < adj[i])
			j++;
		else {
			vertex w = adj[i];
			if (isSmaller (xadj, w, u) && isSmaller (xadj, w, v))
				intersection.push_back (w);
			i++;
			j++;
		}
	}
}

void storeCount4cliques (edge nEdge, edge* xadj, vertex* ordered_adj, edge* ordered_xadj, couple* el, triangle_id* tl, edge* xtl, vector<edge>* fours) {

#pragma omp parallel for
	for (size_t i = 0; i < nEdge; i++) { // u < v
		vertex u = get<0>(el[i]);
		vertex v = get<1>(el[i]);
		vector<vertex> is;
		intersection2 (ordered_adj, ordered_xadj, u, v, is);
		for (size_t j = 0; j < is.size(); j++) { // u < v < w
			vertex wcand = is[j];
			for (size_t k = j+1; k < is.size(); k++) { // u < v < x
				vertex xcand = is[k];
				vertex w = wcand, x = xcand;
				if (isSmaller (xadj, x, w))
					swap (x, w);
				if (findInd (w, x, ordered_adj, ordered_xadj) != -1) {
					vertex in1 = getTriangleId (u, v, w, xtl, tl, ordered_adj, ordered_xadj);
					vertex in2 = getTriangleId (u, v, x, xtl, tl, ordered_adj, ordered_xadj);
					vertex in3 = getTriangleId (u, w, x, xtl, tl, ordered_adj, ordered_xadj);
					vertex in4 = getTriangleId (v, w, x, xtl, tl, ordered_adj, ordered_xadj);

					fours[in1].push_back (in2);
					fours[in1].push_back (in3);
					fours[in1].push_back (in4);
					fours[in2].push_back (in1);
					fours[in2].push_back (in3);
					fours[in2].push_back (in4);
				}
			}
		}
	}

#pragma omp parallel for
	for (size_t i = 0; i < nEdge; i++) { // x < w
		vertex x = get<0>(el[i]);
		vertex w = get<1>(el[i]);
		vector<vertex> is;
		conditionalIntersection2 (ordered_adj, ordered_xadj, x, w, is);
		for (size_t j = 0; j < is.size(); j++) { // u < x < w
			vertex ucand = is[j];
			for (size_t k = j+1; k < is.size(); k++) { // v < x < w
				vertex vcand = is[k];
				vertex u = ucand, v = vcand;
				if (isSmaller (xadj, v, u))
					swap (v, u);
				if (findInd (u, v, ordered_adj, ordered_xadj) != -1) {
					vertex in1 = getTriangleId (u, x, w, xtl, tl, ordered_adj, ordered_xadj);
					vertex in2 = getTriangleId (v, x, w, xtl, tl, ordered_adj, ordered_xadj);
					vertex in3 = getTriangleId (u, v, w, xtl, tl, ordered_adj, ordered_xadj);
					vertex in4 = getTriangleId (u, v, x, xtl, tl, ordered_adj, ordered_xadj);

					fours[in1].push_back (in2);
					fours[in1].push_back (in3);
					fours[in1].push_back (in4);
					fours[in2].push_back (in1);
					fours[in2].push_back (in3);
					fours[in2].push_back (in4);
				}
			}
		}
	}
	return;
}

inline void updateAndNotify (int ind, vertex* P, int newP, vector<vertex>& neigs, bool* changed) {
	P[ind] = newP;
	changed[ind] = true; // *THIS* todo: Are we sure on this?
	for (auto w : neigs)
		if (P[w] >= P[ind])
			changed[w] = true;
}

inline int mapInitialHI_SF (edge ind, vector<edge>* fours, triangle_id* tlist, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {

	HashMap<vertex> gmap (-1);
	vertex greaters = 0;
	vertex equals = 0;
	vertex H = 0;

	for (vertex ii = 0; ii < fours[ind].size(); ii+=3) {
		vertex id1 = fours[ind][ii];
		vertex id2 = fours[ind][ii+1];
		vertex id3 = fours[ind][ii+2];
		vertex sm = min (min (P[id1], P[id2]), P[id3]);

		if (sm == H + 1) {
			if (equals > 0) {
				equals--;
				greaters++;
				gmap[sm]++;
			}
			else { // equals = 0
				H++;
				int gH = 0;
				if (!gmap.hasDefaultValue (H)) {
					gH = gmap[H];
					gmap.erase (H);
				}
				equals = gH + 1;
				greaters -= gH;
			}
		}
		else if (sm > H + 1) {
			if  (equals > 0) {
				equals--;
				greaters++;
				gmap[sm]++;
			}
			else { // equals = 0
				greaters++;
				H++;
				int gH = 0;
				if (!gmap.hasDefaultValue (H))  {
					gH = gmap[H];
					gmap.erase (H);
				}
				equals = gH;
				greaters -= gH;
				gmap[sm]++;
			}
		}
	}

	vertex oP = P[ind];
#ifdef SYNC
	Q[ind] = H;
#else
	P[ind] = H;
#endif

	if (H < oP) {
		return 1;
	}
	else
		return 0;
}

inline int regularUpdateHI_SF (edge ind, vector<edge>* fours, triangle_id* tlist, vertex* P
#ifdef SYNC
			, vertex* Q
#endif
	) {

	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;

	for (vertex ii = 0; ii < fours[ind].size(); ii+=3) {
		vertex id1 = fours[ind][ii];
		vertex id2 = fours[ind][ii+1];
		vertex id3 = fours[ind][ii+2];
		vertex cur_F = min (min (P[id1], P[id2]), P[id3]);

		if (cur_F >= previous_P)
			greaterequals++;
		else
			smallers[cur_F]++;

		if (greaterequals == previous_P) {
			yep_set = true;
			break;
		}
	}
	if (!yep_set && fours[ind].size() > 0) {
		vertex j;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j)
				break;
		}
#ifdef SYNC
		if (Q[ind] != j) {
			Q[ind] = j;
			return 1;
		}
#else
		if (P[ind] != j) {
			P[ind] = j;
			return 1;
		}
#endif
	}
	return 0;
}

inline int efficientUpdateHI_SF (edge ind, vector<edge>* fours, triangle_id* tlist, vertex* P, bool* changed) {

	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;

	for (vertex ii = 0; ii < fours[ind].size(); ii+=3) {
		vertex id1 = fours[ind][ii];
		vertex id2 = fours[ind][ii+1];
		vertex id3 = fours[ind][ii+2];
		vertex cur_F = min (min (P[id1], P[id2]), P[id3]);

		if (cur_F >= previous_P)
			greaterequals++;
		else
			smallers[cur_F]++;
		if (greaterequals == previous_P) {
			yep_set = true;
			break;
		}
	}

	if (!yep_set && fours[ind].size() > 0) {
		vertex j;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}
		updateAndNotify (ind, P, j, fours[ind], changed);
		return 1;
	}
	return 0;
}



// stores the triangle-4clique graph; base AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void baseLocal34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;
	cout << "# triangles: " << tlist.size() << endl;

	// 4-clique counting
	P = (vertex *) malloc (sizeof(vertex) * tlist.size());

#ifdef SYNC
	printf ("It is SYNC\n");
	vertex* Q = (vertex *) malloc (sizeof(vertex) * tlist.size());
#else
	printf ("It is ASYNC\n");
#endif

	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}
	vector<edge> fours[tlist.size()];
	storeCount4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, Ntlist, xtl, fours);

	lol fccount = 0;
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++) {
		P[i] = fours[i].size() / 3;
#ifndef FAST
#pragma omp atomic
		fccount += P[i];
#endif
	}


	timestamp t_fc;
#ifndef FAST
	fccount /= 4;
	cout << "# 4-cliques: " << fccount << endl;
#endif
	cout << "4-clique counting: " << t_fc - t_tri << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

	edge sz = tlist.size();
#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < sz; ind++) {
		mapInitialHI_SF (ind, fours, Ntlist, P
#ifdef SYNC
				, Q
#endif
		);
	}

#ifdef SYNC
	memcpy (P, Q, sizeof(vertex) * sz);
#endif

	timestamp ts5;
	cout << "H 0 time: " << ts5 - t_fc << endl;

#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif

	oc++;

	while (flag) {
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < sz; ind++) {
			int fl = regularUpdateHI_SF (ind, fours, Ntlist, P
#ifdef SYNC
					, Q
#endif
			);
			if (fl == 1)
				flag = true;
		}

#ifdef SYNC
		memcpy (P, Q, sizeof(vertex) * sz);
#endif

#ifdef DUMP_Hs
		print_Ks (sz, P, vfile, oc);
#endif

		timestamp t_end;
		cout << "H " << oc << " time: " << t_end - t_begin - td << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (sz, P, vfile);
#endif

	free (P);
	free (xel);
	free (el);
	free (Ntlist);
	free (xtl);
	free (ordered_adj);
	free (ordered_xadj);
#ifdef SYNC
	free (Q);
#endif
	return;
}

// stores the triangle-4clique graph; AND algorithm with the notification mechanism
void nmLocal34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;

	// 4-clique counting
	P = (vertex *) malloc (sizeof(vertex) * tlist.size());

	cout << "# triangles: " << tlist.size() << endl;
	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	vector<edge> fours[tlist.size()];
	storeCount4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, Ntlist, xtl, fours);

	lol fccount = 0;
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++) {
		P[i] = fours[i].size() / 3;
#ifndef FAST
#pragma omp atomic
		fccount += P[i];
#endif
	}


	timestamp t_fc;
#ifndef FAST
	fccount /= 4;
	cout << "# 4-cliques: " << fccount << endl;
#endif
	cout << "4-clique counting: " << t_fc - t_tri << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

	edge sz = tlist.size();

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < sz; ind++) {
		mapInitialHI_SF (ind, fours, Ntlist, P);
	}

	bool changed[sz];
	memset (changed, 255, sizeof(bool) * sz); // set all true

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_fc - td << endl;
#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif
	oc++;
	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < sz; ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI_SF (ind, fours, Ntlist, P, changed);
			if (a == 1)
				flag = true;
		}
		timestamp td2;

#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif
		timestamp td3;
		td += td3 - td2;

		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}
	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (sz, P, vfile);
#endif


	free (P);
	free (xel);
	free (el);
	free (Ntlist);
	free (xtl);
	free (ordered_adj);
	free (ordered_xadj);
	return;
#endif
}



// Partly parallel 3,4 computation (only 4clique counting is parallel)
void k34_SF (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* F, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);



	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;
	cout << "# triangles: " << tlist.size() << endl;

	// 4-clique counting
	vertex* fc = (vertex *) malloc (sizeof(vertex) * tlist.size());
	F = (vertex *) malloc (sizeof(vertex) * tlist.size());
	memset (F, -1, sizeof(vertex) * tlist.size());

	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}
	vector<edge> fours[tlist.size()];
	storeCount4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, Ntlist, xtl, fours);

	lol fccount = 0;
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++) {
		fc[i] = fours[i].size() / 3;
#ifndef FAST
#pragma omp atomic
		fccount += fc[i];
#endif
	}

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif

	// Peeling
	Naive_Bucket na_bs;
	edge cid = 0;
	na_bs.Initialize (tlist.size(), tlist.size());
	for (size_t i = 0; i < tlist.size(); i++) {
		if (fc[i] > 0) // this is 4c count
			na_bs.Insert (cid, fc[i]);
		cid++;
	}
	vertex fc_of_uvw = 0;
	int order = 0;

	while (1) {
		edge uvw, value_of_uvw ;
		edge ret = na_bs.PopMin(&uvw, &value_of_uvw);
		if (ret == -1)
			break;

		if (value_of_uvw == 0)
			continue;

		fc_of_uvw = F[uvw] = value_of_uvw;

		for (vertex ii = 0; ii < fours[uvw].size(); ii+=3) {
			vertex id1 = fours[uvw][ii];
			vertex id2 = fours[uvw][ii+1];
			vertex id3 = fours[uvw][ii+2];

			if (F[id1] == -1 && F[id2] == -1 && F[id3] == -1) {
				int cur_id1 = na_bs.CurrentValue(id1);
				int cur_id2 = na_bs.CurrentValue(id2);
				int cur_id3 = na_bs.CurrentValue(id3);

				if (na_bs.CurrentValue(id1) > fc_of_uvw)
					na_bs.DecVal(id1);
				if (na_bs.CurrentValue(id2) > fc_of_uvw)
					na_bs.DecVal(id2);
				if (na_bs.CurrentValue(id3) > fc_of_uvw)
					na_bs.DecVal(id3);
			}
		}
	}

	na_bs.Free();
	cout << "Max 34 number: " << fc_of_uvw << endl;
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (tlist.size(), F, vfile);
#endif
	return;
}



/* DOES NOT STORE TRIANGLE_4CLIQUE GRAPHS; SPACE-EFFICIENT (AND SLOWER) */

inline vertex getComplexTriangleId (vertex a, vertex b, vertex c, edge* xtl, triangle_id* tlist, edge* xadj, vertex* adj, edge* ordered_xadj, vertex* ordered_adj) {
	vertex u, v, w;
	if (isSmaller (xadj, c, a)) {
		u = c;
		v = a;
		w = b;
	}
	else if (isSmaller (xadj, c, b)) {
		u = a;
		v = c;
		w = b;
	}
	else {
		u = a;
		v = b;
		w = c;
	}
	return getTriangleId (u, v, w, xtl, tlist, ordered_adj, ordered_xadj);
}

lol intersection3for4cliques (edge nEdge, edge* xadj, vertex* ordered_adj, edge* ordered_xadj, couple* el, edge* xtl, triangle_id* tlist) {
	lol count = 0;
#pragma omp parallel for
	for (size_t i = 0; i < nEdge; i++) { // u < v
		vertex u = get<0>(el[i]);
		vertex v = get<1>(el[i]);
		vector<vertex> is;
		intersection2 (ordered_adj, ordered_xadj, u, v, is);
		for (size_t j = 0; j < is.size(); j++) { // u < v < w
			vertex wcand = is[j];
			for (size_t k = j+1; k < is.size(); k++) { // u < v < x
				vertex xcand = is[k];
				vertex w = wcand, x = xcand;
				if (isSmaller (xadj, x, w))
					swap (w, x);
				if (findInd (w, x, ordered_adj, ordered_xadj) != -1) {
					vertex in1 = getTriangleId (u, v, w, xtl, tlist, ordered_adj, ordered_xadj);
					vertex in2 = getTriangleId (u, v, x, xtl, tlist, ordered_adj, ordered_xadj);
					vertex in3 = getTriangleId (u, w, x, xtl, tlist, ordered_adj, ordered_xadj);
					vertex in4 = getTriangleId (v, w, x, xtl, tlist, ordered_adj, ordered_xadj);
#pragma omp atomic
					(tlist[in1].id)++;

#pragma omp atomic
					(tlist[in2].id)++;

#pragma omp atomic
					(tlist[in3].id)++;

#pragma omp atomic
					(tlist[in4].id)++;
#ifndef FAST
#pragma omp atomic
					count++;
#endif
				}
			}
		}
	}
	return count;
}

inline int mapInitialHI (edge ind, vertex* adj, edge* xadj, triangle_id* tlist, edge* xtl, vertex* ordered_adj, edge* ordered_xadj, vertex* P
#ifdef SYNC
		, vertex* Q
#endif
) {

	HashMap<vertex> gmap (-1);
	vertex greaters = 0;
	vertex equals = 0;
	vertex H = 0;
	vertex u = get<0> (tlist[ind].triple);
	vertex v = get<1> (tlist[ind].triple);
	vertex w = get<2> (tlist[ind].triple);

	vertex ii = xadj[u];
	vertex jj = xadj[v];
	vertex kk = xadj[w];
	vertex u_size = xadj[u+1];
	vertex v_size = xadj[v+1];
	vertex w_size = xadj[w+1];

	while (ii < u_size && jj < v_size && kk < w_size) {
		vertex a = adj[ii];
		vertex b = adj[jj];
		vertex c = adj[kk];
		if (a == b && a == c) {
			vertex x = a;
			vertex id1 = getComplexTriangleId (u, v, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id2 = getComplexTriangleId (u, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id3 = getComplexTriangleId (v, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);

			vertex sm = min (min (P[id1], P[id2]), P[id3]);

			if (sm == H + 1) {
				if (equals > 0) {
					equals--;
					greaters++;
					gmap[sm]++;
				}
				else { // equals = 0
					H++;
					int gH = 0;
					if (!gmap.hasDefaultValue (H)) {
						gH = gmap[H];
						gmap.erase (H);
					}
					equals = gH + 1;
					greaters -= gH;
				}
			}
			else if (sm > H + 1) {
				if  (equals > 0) {
					equals--;
					greaters++;
					gmap[sm]++;
				}
				else { // equals = 0
					greaters++;
					H++;
					int gH = 0;
					if (!gmap.hasDefaultValue (H))  {
						gH = gmap[H];
						gmap.erase (H);
					}
					equals = gH;
					greaters -= gH;
					gmap[sm]++;
				}
			}
			ii++;
			jj++;
			kk++;
		}
		else {
			vertex m = max (max (a, b), c);
			if (a != m)
				ii++;
			if (b != m)
				jj++;
			if (c != m)
				kk++;
		}
	}

	vertex oP = P[ind];
#ifdef SYNC
	Q[ind] = H;
#else
	P[ind] = H;
#endif

	if (H < oP) {
		return 1;
	}
	else
		return 0;
}

inline int regularUpdateHI (edge ind, vertex* adj, edge* xadj, triangle_id* tlist, edge* xtl, vertex* ordered_adj, edge* ordered_xadj, vertex* P
#ifdef SYNC
			, vertex* Q
#endif
	) {

	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;

	vertex u = get<0> (tlist[ind].triple);
	vertex v = get<1> (tlist[ind].triple);
	vertex w = get<2> (tlist[ind].triple);

	vertex ii = xadj[u];
	vertex jj = xadj[v];
	vertex kk = xadj[w];

	vertex u_size = xadj[u+1];
	vertex v_size = xadj[v+1];
	vertex w_size = xadj[w+1];

	bool has4Clique = false;
	while (ii < u_size && jj < v_size && kk < w_size) {
		vertex a = adj[ii];
		vertex b = adj[jj];
		vertex c = adj[kk];
		if (a == b && a == c) {
			has4Clique = true;
			vertex x = a;
			vertex id1 = getComplexTriangleId (u, v, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id2 = getComplexTriangleId (u, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id3 = getComplexTriangleId (v, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex cur_F = min (min (P[id1], P[id2]), P[id3]);

			if (cur_F >= previous_P)
				greaterequals++;
			else
				smallers[cur_F]++;
			if (greaterequals == previous_P) {
				yep_set = true;
				break;
			}
			ii++; jj++; kk++;
		}
		else {
			vertex m = max (max (a, b), c);
			if (a != m)
				ii++;
			if (b != m)
				jj++;
			if (c != m)
				kk++;
		}
	}

	if (!yep_set && has4Clique) {
		vertex j;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j)
				break;
		}
#ifdef SYNC
		if (Q[ind] != j) {
			Q[ind] = j;
			return 1;
		}
#else
		if (P[ind] != j) {
			P[ind] = j;
			return 1;
		}
#endif
	}
	return 0;
}

inline int efficientUpdateHI (edge ind, vertex* adj, edge* xadj, triangle_id* tlist, edge* xtl, vertex* ordered_adj, edge* ordered_xadj, vertex* P, bool* changed) {

	vector<vertex> neigs;
	vertex previous_P = P[ind];
	vertex greaterequals = 0;
	vector<vertex> smallers (previous_P, 0);
	bool yep_set = false;

	vertex u = get<0> (tlist[ind].triple);
	vertex v = get<1> (tlist[ind].triple);
	vertex w = get<2> (tlist[ind].triple);

	vertex ii = xadj[u];
	vertex jj = xadj[v];
	vertex kk = xadj[w];

	vertex u_size = xadj[u+1];
	vertex v_size = xadj[v+1];
	vertex w_size = xadj[w+1];

	bool has4Clique = false;
	while (ii < u_size && jj < v_size && kk < w_size) {
		vertex a = adj[ii];
		vertex b = adj[jj];
		vertex c = adj[kk];
		if (a == b && a == c) {
			has4Clique = true;
			vertex x = a;
			vertex id1 = getComplexTriangleId (u, v, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id2 = getComplexTriangleId (u, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id3 = getComplexTriangleId (v, w, x, xtl, tlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex cur_F = min (min (P[id1], P[id2]), P[id3]);

			if (P[id1] <= previous_P)
				neigs.push_back (id1);
			if (P[id2] <= previous_P)
				neigs.push_back (id2);
			if (P[id3] <= previous_P)
				neigs.push_back (id3);

			if (cur_F >= previous_P)
				greaterequals++;
			else
				smallers[cur_F]++;
			if (greaterequals == previous_P) {
				yep_set = true;
				break;
			}
			ii++; jj++; kk++;
		}
		else {
			vertex m = max (max (a, b), c);
			if (a != m)
				ii++;
			if (b != m)
				jj++;
			if (c != m)
				kk++;
		}
	}


	if (!yep_set && has4Clique) {
		vertex j;
		for (j = previous_P - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j) {
				break;
			}
		}
		updateAndNotify (ind, P, j, neigs, changed);
		return 1;
	}
	return 0;
}



// base AND and SND algorithms, no notification mechanism. compile with SYNC=yes to get the synchronous mode (SND)
void baseLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;
	cout << "# triangles: " << tlist.size() << endl;

	// 4-clique counting
	P = (vertex *) malloc (sizeof(vertex) * tlist.size());

#ifdef SYNC
	printf ("It is SYNC\n");
	vertex* Q = (vertex *) malloc (sizeof(vertex) * tlist.size());
#else
	printf ("It is ASYNC\n");
#endif

	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	lol fccount = intersection3for4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, xtl, Ntlist);
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++)
		P[i] = Ntlist[i].id;

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

#ifdef DEBUG_000
	timestamp ts1;
	int nt = 1, tn = 0;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
		tn = omp_get_thread_num();
	}
	int* ccs = (int *) calloc (nt, sizeof(int));
	timestamp ts2;
	td += ts2 - ts1;
#endif

	edge sz = tlist.size();
#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < sz; ind++) {
#ifdef DEBUG_000
		ccs[tn] += 	mapInitialHI (ind, adj, xadj, Ntlist, xel, xtl, ordered_adj, ordered_xadj, P
#ifdef SYNC
				, Q
#endif
		);
#else
		mapInitialHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P
#ifdef SYNC
				, Q
#endif
		);
#endif
	}

#ifdef SYNC
	memcpy (P, Q, sizeof(vertex) * sz);
#endif

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_fc - td << endl;
#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif
	oc++;

	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < sz; ind++) {
			int fl = regularUpdateHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P
#ifdef SYNC
				, Q
#endif
			);
			if (fl == 1)
				flag = true;
		}
		timestamp td2;

#ifdef SYNC
	memcpy (P, Q, sizeof(vertex) * sz);
#endif

#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif

		timestamp td3;
		td += td3 - td2;
		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (sz, P, vfile);
#endif

	free (P);
	free (xel);
	free (el);
	free (Ntlist);
	free (xtl);
	free (ordered_adj);
	free (ordered_xadj);
#ifdef SYNC
	free (Q);
#endif
	return;
}

// AND algorithm with the notification mechanism
void nmLocal34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, const char* vfile) {
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;

	// 4-clique counting
	P = (vertex *) malloc (sizeof(vertex) * tlist.size());

	cout << "# triangles: " << tlist.size() << endl;
	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	lol fccount = intersection3for4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, xtl, Ntlist);
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++)
		P[i] = Ntlist[i].id;

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif

	int oc = 0;
	bool flag = true;
	timestamp td (0, 0);
#ifdef DEBUG_000
	timestamp ts1;
	int nt = 1, tn = 0;
#pragma omp parallel
	{
		nt = omp_get_num_threads();
		tn = omp_get_thread_num();
	}
	int* ccs = (int *) calloc (nt, sizeof(int));
	timestamp ts2;
	td += ts2 - ts1;
#endif

	edge sz = tlist.size();

#pragma omp parallel for schedule (dynamic, 1000)
	for (edge ind = 0; ind < sz; ind++) {
#ifdef DEBUG_000
		ccs[tn] += 	mapInitialHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P);
#else
		mapInitialHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P);
#endif
	}

	bool changed[sz];
	memset (changed, 255, sizeof(bool) * sz); // set all true

	timestamp t_init;
	cout << "H 0 time: " << t_init - t_fc - td << endl;
#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif
	oc++;
	while (flag) {
		timestamp td1;
		flag = false;

#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < sz; ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P, changed);
#ifdef DEBUG_000
			ccs[tn] += a;
#endif
			if (a == 1)
				flag = true;
		}
		timestamp td2;
#ifdef DEBUG_000
		int degisenler = 0;
		for(int i = 0; i < nt; i++)
			degisenler += ccs[i];
		memset (ccs, 0, nt * sizeof(int));
		cout << "CHANGEDS: " << degisenler << endl;
#endif


#ifdef DUMP_Hs
	print_Ks (sz, P, vfile, oc);
#endif

		timestamp td3;
		td += td3 - td2;
		cout << "H " << oc << " time: " << td2 - td1 << endl;
		oc++;
	}
	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin - td << endl;

#ifdef DUMP_K
	print_Ks (sz, P, vfile);
#endif


	free (P);
	free (xel);
	free (el);
	free (Ntlist);
	free (xtl);
	free (ordered_adj);
	free (ordered_xadj);
	return;
#endif
}


// Partly parallel 3,4 computation (only 4clique counting is parallel)
void k34 (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* F, const char* vfile) {

	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;

	// 4-clique counting
	F = (vertex *) malloc (sizeof(vertex) * tlist.size());
	memset (F, -1, sizeof(vertex) * tlist.size());

	cout << "# triangles: " << tlist.size() << endl;
	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	lol fccount = intersection3for4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, xtl, Ntlist);

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif

	// Peeling
	Naive_Bucket na_bs;
	edge cid = 0;
	na_bs.Initialize (tlist.size(), tlist.size());
	for (size_t i = 0; i < tlist.size(); i++) {
		if (Ntlist[i].id > 0) // this is 4c count
			na_bs.Insert (cid, Ntlist[i].id);
		cid++;
	}
	vertex fc_of_uvw = 0;
	int order = 0;

	while (1) {
		edge uvw, value_of_uvw ;
		edge ret = na_bs.PopMin(&uvw, &value_of_uvw);
		if (ret == -1)
			break;

		if (value_of_uvw == 0)
			continue;

		fc_of_uvw = F[uvw] = value_of_uvw;
		vertex u = get<0> (tlist[uvw].triple);
		vertex v = get<1> (tlist[uvw].triple);
		vertex w = get<2> (tlist[uvw].triple);

		vector<vertex> intersection;
		intersection3 (adj, xadj, u, v, w, intersection);
		for (auto x: intersection) {
			vertex id1 = getComplexTriangleId (u, v, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id2 = getComplexTriangleId (u, w, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);
			vertex id3 = getComplexTriangleId (v, w, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);

			if (F[id1] == -1 && F[id2] == -1 && F[id3] == -1) {
				int cur_id1 = na_bs.CurrentValue(id1);
				int cur_id2 = na_bs.CurrentValue(id2);
				int cur_id3 = na_bs.CurrentValue(id3);

				if (na_bs.CurrentValue(id1) > fc_of_uvw)
					na_bs.DecVal(id1);
				if (na_bs.CurrentValue(id2) > fc_of_uvw)
					na_bs.DecVal(id2);
				if (na_bs.CurrentValue(id3) > fc_of_uvw)
					na_bs.DecVal(id3);
			}
		}
	}

	na_bs.Free();
	cout << "Max 34 number: " << fc_of_uvw << endl;
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;

#ifdef DUMP_K
	print_Ks (tlist.size(), F, vfile);
#endif
	return;
}



// ONGOING WORK

// Finds the max K value by iterating on top-number high degree vertices
void fast34DegeneracyNumber (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* P, vertex topK) {

	int number = topK; // only top-number vertices in degree are executed
#ifdef SYNC
	printf ("No SYNC for notification-mechanism\n");
	exit(1);
#else
	timestamp t_begin;

	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;

	// 4-clique counting
	P = (vertex *) malloc (sizeof(vertex) * tlist.size());

	cout << "# triangles: " << tlist.size() << endl;
	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	lol fccount = intersection3for4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, xtl, Ntlist);
#pragma omp parallel for
	for (edge i = 0; i < tlist.size(); i++)
		P[i] = Ntlist[i].id;

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif

	timestamp td (0, 0);
	int oc = 0;
	bool flag = true;

	bool changed[tlist.size()];
	memset (changed, 255, sizeof(bool) * tlist.size()); // set all true

	vector<eda> KK;
	for (vertex i = 0; i < tlist.size(); i++)
		KK.push_back (make_tuple(i, P[i]));

	sort (KK.begin(), KK.end(), kksort);

	edge start = 0;
	for (edge i = start; i < number; i++)
		changed[get<0>(KK[i])] = true;


	while (flag) {

		memset (changed, 0, sizeof(bool) * tlist.size()); // set all false
		for (edge i = start; i < number; i++)
			changed[get<0>(KK[i])] = true; // only top-k true

		timestamp it1;
		flag = false;
#pragma omp parallel for schedule (dynamic, 1000)
		for (edge ind = 0; ind < tlist.size(); ind++) {
			if (!changed[ind])
				continue;
			changed[ind] = false;
			int a = efficientUpdateHI (ind, adj, xadj, Ntlist, xtl, ordered_adj, ordered_xadj, P, changed);
			if (a == 1)
				flag = true;
		}

		timestamp it3;
		cout << "H " << oc << " time: " << it3 - it1 << endl;
		oc++;
	}

	printf ("Converges at %d\n", oc);
	timestamp t_end;
	cout << "Total time: " << t_end - t_begin << endl;


	vertex maxK = 0;
	for (edge i = start; i < number; i++) {
		int curP = P[get<0>(KK[i])];
		if (curP > maxK)
			maxK = curP;
//		printf ("%d goes from %d to %d\n", get<0>(KK[i]), get<1>(KK[i]), curP);
	}

	printf ("max K: %d\n", maxK);
#endif
}


/*
- L_1: All those triangles with the minimum 4clique count.
- To find L_i: mark all triangles in L_j deleted (j < i). Then find all triangles with the minimum 4clique count in the remaining graph.
*/
void k34_levels (vertex nVtx, edge nEdge, vertex* adj, edge* xadj, vertex* L, const char* vfile) {

	timestamp t_begin;
	// Ordered (directed) graph creation
	couple* el = (couple *) malloc (sizeof(couple) * nEdge);
	edge* xel = (edge *) malloc (sizeof(edge) * (nVtx+1));
	vertex* ordered_adj = (vertex *) malloc (sizeof(vertex) * nEdge );
	edge* ordered_xadj = (edge *) malloc (sizeof(edge) * (nVtx+1));;

	createOrdered (nVtx, nEdge, adj, xadj, el, xel, ordered_adj, ordered_xadj);

	timestamp t_cog;
	cout << "Creating ordered graph: " << t_cog - t_begin << endl;

	// Enumerating triangles
	edge* xtl = (edge *) malloc (sizeof(edge) * (nEdge + 1));
	vector<triangle_id> tlist;
	enumTriangles (nEdge, ordered_adj, ordered_xadj, el, tlist, xtl);

	timestamp t_tri;
	cout << "Enumerating triangles: " << t_tri - t_cog << endl;
	cout << "# triangles: " << tlist.size() << endl;
	triangle_id* Ntlist = (triangle_id *) malloc (sizeof (triangle_id) * tlist.size());
	for (edge i = 0; i < tlist.size(); i++) {
		Ntlist[i].id = tlist[i].id;
		Ntlist[i].triple = tlist[i].triple;
	}

	// 4-clique counting
	lol fccount = intersection3for4cliques (nEdge, xadj, ordered_adj, ordered_xadj, el, xtl, Ntlist);

	timestamp t_fc;
	cout << "4-clique counting: " << t_fc - t_tri << endl;
#ifndef FAST
	cout << "# 4-cliques: " << fccount << endl;
#endif


	L = (vertex *) malloc (sizeof(vertex) * tlist.size());
	memset (L, -1, sizeof(vertex) * tlist.size());

	// Peeling
	Naive_Bucket na_bs;
	edge cid = 0;
	na_bs.Initialize (tlist.size(), tlist.size());
	for (size_t i = 0; i < tlist.size(); i++) {
		if (Ntlist[i].id > 0) // this is 4c count
			na_bs.Insert (cid, Ntlist[i].id);
		cid++;
	}
	vertex fc_of_uvw = 0;

	int level = 1;
	while (1) {
		edge uvw, min_val, val;
		vector<vertex> mins;

		if (na_bs.PopMin(&uvw, &min_val) == -1) // if the bucket is empty
			break;

		if (min_val == 0)
			continue;

		mins.push_back (uvw);

		while (1) {
			if (na_bs.PopMin(&uvw, &val) == -1) // if the bucket is empty
				break;
			if (val > min_val) {
				na_bs.Insert (uvw, val);
				break;
			}
			mins.push_back (uvw);
		}

		for (auto uvw : mins) {
			L[uvw] = min_val;
			vertex u = get<0> (tlist[uvw].triple);
			vertex v = get<1> (tlist[uvw].triple);
			vertex w = get<2> (tlist[uvw].triple);

			vector<vertex> intersection;
			intersection3 (adj, xadj, u, v, w, intersection);
			for (auto x: intersection) {
				vertex id1 = getComplexTriangleId (u, v, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);
				vertex id2 = getComplexTriangleId (u, w, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);
				vertex id3 = getComplexTriangleId (v, w, x, xtl, Ntlist, xadj, adj, ordered_xadj, ordered_adj);

				if (L[id1] == -1 && L[id2] == -1 && L[id3] == -1) {
					int cur_id1 = na_bs.CurrentValue(id1);
					int cur_id2 = na_bs.CurrentValue(id2);
					int cur_id3 = na_bs.CurrentValue(id3);

					if (na_bs.CurrentValue(id1) > fc_of_uvw)
						na_bs.DecVal(id1);
					if (na_bs.CurrentValue(id2) > fc_of_uvw)
						na_bs.DecVal(id2);
					if (na_bs.CurrentValue(id3) > fc_of_uvw)
						na_bs.DecVal(id3);
				}
			}
		}

		printf ("Level %d K: %d size: %d\n", level, min_val, mins.size());
		level++;
	}
	na_bs.Free();
#ifdef DUMP_K
	print_Ks (nEdge, L, vfile);
#endif
	return;
}
