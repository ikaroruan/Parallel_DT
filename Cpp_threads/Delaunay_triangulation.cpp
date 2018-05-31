#include "Delaunay_triangulation.h"

#ifndef NDEBUG
	std::chrono::high_resolution_clock::time_point start_critical;
	std::chrono::high_resolution_clock::time_point start_region1;
	std::chrono::high_resolution_clock::time_point start_region2;
	std::chrono::high_resolution_clock::time_point start_region3;
	std::chrono::high_resolution_clock::time_point start_region4;
	std::chrono::high_resolution_clock::time_point end_critical;
	std::chrono::high_resolution_clock::time_point end_region1;
	std::chrono::high_resolution_clock::time_point end_region2;
	std::chrono::high_resolution_clock::time_point end_region3;
	std::chrono::high_resolution_clock::time_point end_region4;
	std::chrono::duration<double> elapsed_critical;
	std::chrono::duration<double> elapsed_region1;
	std::chrono::duration<double> elapsed_region2;
	std::chrono::duration<double> elapsed_region3;
	std::chrono::duration<double> elapsed_region4;

	std::chrono::high_resolution_clock::time_point start_locate;
	std::chrono::high_resolution_clock::time_point start_conflict;
	std::chrono::high_resolution_clock::time_point start_update;
	std::chrono::high_resolution_clock::time_point end_locate;
	std::chrono::high_resolution_clock::time_point end_conflict;
	std::chrono::high_resolution_clock::time_point end_update;
#endif


//					  ===================== 					 //
//                                              INSERTION                                                //
//					  ===================== 					 //

void Delaunay_triangulation::insert(std::vector<Point>::iterator start, std::vector<Point>::iterator end)
{
	typedef std::vector<Point>::iterator	Point_iterator;

	Spatial_sort_traits sst;
	CGAL::spatial_sort(start, end, sst);
	for(Point_iterator it = start; it != end; ++it)
		insert(*it);
}

Vertex_iterator Delaunay_triangulation::insert(Point& p)
{
	switch(dimension()){
		case -1:{
			return insert_first(p);
			break;
		}
		case 0:{
			return insert_second(p);
			break;
		}
		case 1:{
			return insert_dimension_1(p);
			break;
		}
		case 2:{
			return insert_dimension_2(p);
			break;
		}
		default:
			std::cerr << "ERROR: wrong dimension when trying to insert.\n";
	}
	// Should not reach this line.
	return nullptr;
}

Vertex_iterator Delaunay_triangulation::insert_first(Point& p)
{
	Vertex_iterator v = vertices().insert(Vertex(p));
	Face fc(v, infinite_vertex(), infinite_vertex());
	Face_iterator fcit = faces().insert(fc);
	v->incident_face(fcit);
	infinite_vertex()->incident_face(fcit);
	dimension( dimension() + 1 );

	// Setting all face neighbors pointing to itself
	for(int i = 0; i <= 2; ++i)
		fcit->neighbor(fcit, i);

	return v;
}

Vertex_iterator Delaunay_triangulation::insert_second(Point& p)
{
	Vertex_iterator v = vertices().insert(Vertex(p));
	Face_iterator fc = faces().begin();
	dimension(dimension() + 1);
	
	if(fc->vertex(dimension()) != infinite_vertex())
		std::cerr << "Infinite vertex not set on correct position. (insert_second) \n";

	if(fc->vertex(dimension()) == infinite_vertex())
		fc->vertex(v, dimension());
	else if(fc->vertex(dimension() - 1) == infinite_vertex())
		fc->vertex(v, dimension() - 1);
	else
		std::cerr << "ERROR: wrong vertex position for insertion.\n";
	v->incident_face(fc);

	return v;
}

Vertex_iterator Delaunay_triangulation::insert_dimension_1(Point& p)
{
	Vertex_location lt;
	int li = -1;
	Face_iterator fc = locate(p, lt, li);

	switch(lt){
		case ON_EDGE:
			return insert_on_edge_1(fc, p, lt, li);
			break;
		case OUTSIDE_CONVEX_HULL:
			return insert_outside_convex_hull_1(p, fc, li);
			break;
		case ON_VERTEX:
			return fc->vertex(li);
			break;
		case ON_INTERIOR:
			std::cerr << "ERROR: Locate function returned wrong result for dimension 1.\n";
			break;
		default:
			std::cerr << "ERROR: Wrong location when inserting dimension 1.\n";
			break;
	}
	
	// May not reach this line
	return nullptr;
}

Vertex_iterator Delaunay_triangulation::insert_on_edge_1(Face_iterator fc, Point& p, Vertex_location lt, int li)
{
	Vertex_iterator v = vertices().insert(Vertex(p));
	Face_iterator ff = nullptr;
	// li == -1 means that p is collinear between two points of fc
	if(li == -1){
		int infinite_index = fc->index(infinite_vertex());
		//Vertex_iterator v0 = fc->vertex(ccw(infinite_index));
		Vertex_iterator v1 = fc->vertex(cw(infinite_index));
		Face_iterator n0 = fc->neighbor(ccw(infinite_index));
		Face_iterator n1 = fc->neighbor(cw(infinite_index)); 
		Face ffobj(v, v1, infinite_vertex());
		ff = faces().insert(ffobj);
		fc->vertex(v, fc->index(v1));
		v->incident_face(fc); v1->incident_face(ff);

		fc->neighbor(ff, ccw(infinite_index));
		ff->neighbor(fc, cw(infinite_index));
		int i0 = cw(n0->index(infinite_vertex())); 
		int i1 = ccw(n1->index(infinite_vertex()));
		if(n0 != nullptr){
			ff->neighbor(n0, ccw(infinite_index)); 
			ff->neighbor(ff, ff->index(infinite_vertex()));
			n0->neighbor(ff, i0);
		}
		if(vertices().size() > 4){
			fc->neighbor(n1, cw(infinite_index)); 
			n1->neighbor(fc, i1);
		}
	}
	else{ // locate returns the closest vertex to the new point and one incident face
		// If the new vertex is to be inserted to the left: li == 0 && ui == 1,
		// but if to the right: li == 1 && ui == 0
		int ui = -1;
		if(li == 0)	 ui = 1;
		else if(li == 1) ui = 0;
		else	std::cerr << "ERROR: indices not match (insert_on_edge_1).\n";

		Vertex_iterator u = fc->vertex(li);
		Face_iterator n = fc->neighbor(ui);

		Face ffobj(infinite_vertex(), infinite_vertex(), infinite_vertex());
		ff = faces().insert(ffobj);
		ff->vertex(u, ui); ff->vertex(v, li); ff->vertex(infinite_vertex(), 2);
		v->incident_face(ff);
		
		ff->neighbor(fc, li); fc->neighbor(ff, ui);  
		int ni = li; //n->index(fc);
		if(n != nullptr){
			ff->neighbor(n, ui);  
			n->neighbor(ff, ni);
		}
	}
	int infinite_index = ff->index(infinite_vertex());
	ff->neighbor(ff, infinite_index);
	return v;
}

Vertex_iterator Delaunay_triangulation::insert_outside_convex_hull_1(Point& p, Face_iterator fc, int li)
{
	Vertex_iterator vmax = fc->vertex(li);
	std::list<Face_iterator> infinite_created;
	std::list<Face_iterator> to_copy;
	Vertex_iterator v = vertices().insert(Vertex(p));

	// FIXME: Adding new faces to the compact container while iterating along its elements
	// arises in infinite loop.
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		to_copy.push_back(it);
	}

	Face_iterator ff = Faces_begin();
	for(std::list<Face_iterator>::iterator it = to_copy.begin(); it != to_copy.end(); ++it){
		Face_iterator current_face = *it;
		Face fi(*current_face);
		Face_iterator fiit = faces().insert(fi);
		current_face->vertex(v, current_face->index(infinite_vertex()));
		ff = current_face;

		fiit->neighbor(current_face, fiit->index(infinite_vertex()));
		current_face->neighbor(fiit, current_face->index(v));
		infinite_created.push_back(fiit);
	}
	
	Face_iterator fi = *(infinite_created.begin());
	infinite_vertex()->incident_face(fi);
	v->incident_face(ff);

	// Checking if reorientation of vertices is needed
	Point pa = ff->vertex( ccw(ff->index(v)) )->point();
	Point pb = ff->vertex( cw(ff->index(v)) )->point();

	if(orientation_test(pa, pb, p) < 0){
		for(std::list<Face_iterator>::iterator it = infinite_created.begin(); it != infinite_created.end(); ++it)
			((*it)->neighbor((*it)->index(infinite_vertex())))->reorient();
	}
	else{
		for(std::list<Face_iterator>::iterator it = infinite_created.begin(); it != infinite_created.end(); ++it){
			(*it)->reorient();
		}
	}

	Face_circulator* ft = new Face_circulator(v);
	Face_iterator done = ft->current_face();
	do{
		Face_iterator fl = ft->current_face();
		ft->next();
		Face_iterator fr = ft->current_face();
		Face_iterator nl = fl->neighbor(fl->index(v));
		Face_iterator nr = fr->neighbor(fr->index(v));
		
		nl->neighbor(nr, cw(nl->index(infinite_vertex())));
		nr->neighbor(nl, ccw(nr->index(infinite_vertex())));
	}while(ft->current_face() != done);
	delete ft;
	
	// Creating new infinitex faces
	Face_iterator fh;
	int i = fc->index(v);
	int fci, fhi, ii;
	li = fc->index(vmax); // The face reorientation may change its indice.
	if(ccw(i) == li){  
		fci = cw(i);  fh = fc->neighbor(fci);
		ii = fh->index(v); fhi = ccw(ii);
	}
	else{  
		fci = ccw(i); fh = fc->neighbor(fci);
		ii = fh->index(v); fhi = cw(ii);
	}
	
	//fh = fc->neighbor(fci);

	if(vertices().size() == 4 && fci == fhi){
		if(ccw(i) == i)  fhi = cw(i);
		else		 fhi = ccw(i);
	}

	Face f1_obj(fc->vertex(cw(fci)), fc->vertex(ccw(fci)), infinite_vertex());
	Face f2_obj(fh->vertex(cw(fhi)), fh->vertex(ccw(fhi)), infinite_vertex());
	Face_iterator f1 = faces().insert(f1_obj);
	Face_iterator f2 = faces().insert(f2_obj);

	// Updating adjacencies for the new created faces
	Face_iterator nc = fc->neighbor(i); Face_iterator nh = fh->neighbor(ii);
	fc->neighbor(f1, fci); fh->neighbor(f2, fhi);
	f1->neighbor(fc, f1->index(infinite_vertex())); f2->neighbor(fh, f2->index(infinite_vertex()));
	f1->neighbor(nc, f1->index(v)); f2->neighbor(nh, f2->index(v));
	nc->neighbor(f1, nc->index(fc->vertex(fci)));
	nh->neighbor(f2, nh->index(fh->vertex(fhi)));
	f1->neighbor(f2, f1->index(static_cast<Face_iterator>(nullptr)));
	f2->neighbor(f1, f2->index(static_cast<Face_iterator>(nullptr)));

	dimension(dimension() + 1);

	return v;
}

Vertex_iterator Delaunay_triangulation::insert_dimension_2(Point& p)
{
	Vertex_location lc;
	int li = -1; // used for locating the triangle edge, when point lies on edge
	Face_iterator fc = locate(p, lc, li);
	
	switch(lc){
		case ON_EDGE:
			return insert_on_edge_2(fc, p, lc, li);
			break;
		case ON_INTERIOR:
			return insert_in_face(fc, p, lc, li);
			break;
		case ON_VERTEX:
			return is_vertex(fc, p);
			break;
		case OUTSIDE_CONVEX_HULL:
			return insert_outside_convex_hull_2(fc, p, lc, li);
			break;
		default:
			std::cerr << "ERROR: wrong location when inserting dimension 2.\n";
	}

	// FOR NOW, should never reach this line
	return nullptr;
}

Vertex_iterator Delaunay_triangulation::insert_in_hole(Face_iterator fc, Point& p, int li, Vertex_location lc)
{
	std::vector<Face_iterator> region;
	std::vector<std::pair<Face_iterator, int>> cavity;
	Vertex_iterator new_vertex = vertices().insert(Vertex(p));
	std::vector<Face_iterator> faces_created;

	conflict_region(fc, p, region, cavity, lc, li);

	for(unsigned int i = 0; i < cavity.size(); ++i){
		Face_iterator old_face = std::get<0>(cavity[i]);
		int old_int = std::get<1>(cavity[i]);

		// create new face and set vertices according to old face and new vertex to be inserted
		Face_iterator new_face = faces().insert(Face());
		new_face->vertex(old_face->vertex(cw(old_int)), cw(old_int));
		new_face->vertex(old_face->vertex(ccw(old_int)), ccw(old_int));
		new_face->vertex(new_vertex, old_int);

		// setting adjacencies for neighbors outside cavity
		Face_iterator old_neighbor = old_face->neighbor(old_int);
		old_neighbor->neighbor(new_face, old_neighbor->index(old_face));
		new_face->neighbor(old_neighbor, old_int);

		Vertex_iterator v = new_face->vertex(ccw(old_int));
		v->incident_face(new_face);
		if(i == 0)
			new_vertex->incident_face(new_face);
		faces_created.push_back(new_face);
	}

	// adjacencies among new faces
	for(unsigned int i = 0; i < faces_created.size(); ++i){
		Face_iterator fc = faces_created[i];
		int j = fc->index(new_vertex);
		Vertex_iterator vi = fc->vertex(cw(j));
		Face_iterator ff = vi->incident_face();

		if(fc == ff)
			std::cerr << "ERROR: Identical faces when setting neighbor pointers.\n";

		// indices for setting neighbors
		int m, n;
		if(fc->vertex(cw(j)) == vi)	m = ccw(j);
		else				m = cw(j);

		int k = ff->index(new_vertex);
		if(ff->vertex(cw(k)) == vi)	n = ccw(k);
		else				n = cw(k);

		fc->neighbor(ff, m);
		ff->neighbor(fc, n);
	}

	// removing old faces from the conflict region
	for(unsigned int i = 0; i < region.size(); ++i){
		Face_iterator to_remove = region[i];
		faces().erase(to_remove);
		if(key_face() == to_remove)
			key_face(nullptr);
	}
	
	return new_vertex;

}

Vertex_iterator Delaunay_triangulation::insert_in_face(Face_iterator fc, Point& p, Vertex_location lc, int li)
{
	return insert_in_hole(fc, p, -1, lc);
}

Vertex_iterator Delaunay_triangulation::insert_on_edge_2(Face_iterator fc, Point& p, Vertex_location lc, int li)
{	
	return insert_in_hole(fc, p, li, lc);
}

Vertex_iterator Delaunay_triangulation::insert_outside_convex_hull_2(Face_iterator fc, Point& p, Vertex_location lc, int li)
{
	return insert_in_hole(fc, p, -1, lc);
}

//					   ===================						 //
//                                               REMOVAL                                                 //
//					   ===================	 					 //

void Delaunay_triangulation::remove(Vertex_iterator v)
{
	switch(dimension()){
		case -1:
			std::cerr << "ERROR: There is no vertex to be removed.\n";
			break;
		case 0:
			remove_dim_0(v);
			break;
		case 1:
			remove_dim_1(v);
			break;
		case 2:
			remove_dim_2(v);
			break;
		default:
			std::cerr << "ERROR: wrong dimension when trying to remove. (remove point)\n";
			break;
	}
	
}

void Delaunay_triangulation::remove_dim_0(Vertex_iterator v)
{
	Face_iterator fc = faces().begin();
	faces().erase(fc);
	dimension(dimension() - 1);
	if(fc == key_face())
		key_face(nullptr);
	vertices().erase(v);
}

void Delaunay_triangulation::remove_dim_1(Vertex_iterator v)
{
	// If there is only one face to be removed.
	if(vertices().size() == 3){
		Face_iterator fc = faces().begin();
		int ii = fc->index(infinite_vertex());
		
		if(fc->vertex(ccw(ii)) == v)
			fc->vertex(infinite_vertex(), ccw(ii));
		else
			fc->vertex(infinite_vertex(), cw(ii));
		dimension(dimension() - 1);
	}
	// Here, triangulation has more than one face.
	else{
		bool need_reorientation = false;
		Face_iterator fc = v->incident_face();
		int ii = fc->index(infinite_vertex());
		
		int i = -1;
		if(fc->vertex(ccw(ii)) == v){
			i = cw(ii);
			need_reorientation = true;
		}
		else
			i = ccw(ii);
		Face_iterator fn = fc->neighbor(i);
		int j = -1;

		// If fn has v as vertex, then v has two incident faces.
		if(fn->is_vertex(v)){
			j = mirror_index(fc, i);
			Face ff_obj(fc->vertex(i), fn->vertex(j), infinite_vertex());
			Face_iterator ff = faces().insert(ff_obj);
			fc->vertex(i)->incident_face(ff); fn->vertex(j)->incident_face(ff);
			infinite_vertex()->incident_face(ff);

			Face_iterator nc = fc->neighbor(fc->index(v));
			Face_iterator nn = fn->neighbor(fn->index(v));

			ff->neighbor(nn, ff->index(fc->vertex(i)));
			ff->neighbor(nc, ff->index(fn->vertex(j)));
			ff->neighbor(ff, ff->index(infinite_vertex()));
			nc->neighbor(ff, nc->index(fc));
			nn->neighbor(ff, nn->index(fn));
			if(need_reorientation)
				ff->reorient();
			
			faces().erase(fc); faces().erase(fn);
			if(fc == key_face() || fn == key_face())
				key_face(nullptr);
		}
		// Otherwise, v has only one incindent face.
		else{
			if(faces().size() == 2){
				for(int i = 0; i <= 2; ++i)
					fn->neighbor(fn, i);
				if(fc->vertex(i)->incident_face() == fc)
					fc->vertex(i)->incident_face(fn);
			}
			else{
				j = fn->index(fc);
				Face_iterator nc = fc->neighbor(fc->index(v));
				int nci = nc->index(fc);
				if(fc->vertex(i)->incident_face() == fc)
					fc->vertex(i)->incident_face(nc);
			
				nc->neighbor(fn, nci); fn->neighbor(nc, j);
			}
			faces().erase(fc); 
			if(fc == key_face())
				key_face(nullptr);
			if(infinite_vertex()->incident_face() == fc)
				infinite_vertex()->incident_face(fn);
		}
	}
	
	vertices().erase(v);
}

void Delaunay_triangulation::remove_dim_2(Vertex_iterator v)
{
	if(v == nullptr)
		std::cerr << "ERROR: vertex to be removed should be different than nullptr.\n";
	
	if(dim_goes_down(v))
		remove_dim_decrease(v);
	else{
		std::vector<Face_iterator> old_faces;
		std::vector<Face_iterator> new_faces;
		std::queue<std::pair<Face_iterator, int>> queue;
	
		Face_circulator* fh = new Face_circulator(v);
		while(is_infinite(fh->current_face())) 
			fh->next();
	
		Face_iterator done = fh->current_face();
		do{
			old_faces.push_back(fh->current_face());
			fh->next();
		}while(done != fh->current_face());
		delete fh;
	
		Face_iterator fc = old_faces[0]; Face_iterator fn = nullptr;
		int ii = fc->index(v); 
		Vertex_iterator vi = fc->vertex(ccw(ii));

		// Creating new faces based on the cavity edges. Skips first and last edges.
		for(unsigned int i = 1; i < old_faces.size() - 1; ++i){
			if(vertices().size() == 4) break;
			fc = old_faces[i]; ii = fc->index(v);
			Face ff_obj(vi, fc->vertex(ccw(ii)), fc->vertex(cw(ii)));
			Face_iterator ff = faces().insert(ff_obj);
			new_faces.push_back(ff);
			ff->vertex(0)->incident_face(ff); ff->vertex(1)->incident_face(ff); ff->vertex(2)->incident_face(ff); 
	
			fn = fc->neighbor(ii);
			ff->neighbor(fn, ff->index(vi)); fn->neighbor(ff, fn->index(fc));
	
			// Updating adjacencies for first and last face created
			if(i == 1){
				fc = old_faces[i-1]; ii = fc->index(v);
				fn = fc->neighbor(ii);
				ff->neighbor(fn, cw(ff->index(vi))); fn->neighbor(ff, fn->index(fc));
			}
			if(i == old_faces.size() - 2){
				fc = old_faces[i+1]; ii = fc->index(v);
				fn = fc->neighbor(ii);
				ff->neighbor(fn, ccw(ff->index(vi))); fn->neighbor(ff, fn->index(fc));
			}

		}

		// Deleting old vertex and faces
		vertices().erase(v);
		for(std::vector<Face_iterator>::iterator it = old_faces.begin(); it != old_faces.end(); ++it){
			if(infinite_vertex()->incident_face() == *it)
				infinite_vertex()->incident_face((*it)->neighbor(ccw((*it)->index(infinite_vertex()))));
			if((*it) == key_face())
				key_face(nullptr);
			faces().erase(*it);
		}

		// If only one face remains, then it should have all neighbors pointers set to itself.
		if(faces().size() == 1){
			fc = faces().begin();
			fc->neighbor(fc, 0); fc->neighbor(fc, 1); fc->neighbor(fc, 2);
		}

		// Setting adjacencies among new faces.
		// In this approach, last face is skipped in order to avoid redundance.
		if(new_faces.size() > 0){
			for(unsigned int i = 0; i < new_faces.size() - 1; ++i){
				fc = new_faces[i]; fn = new_faces[i+1];
				fc->neighbor(fn, ccw(fc->index(vi))); fn->neighbor(fc, cw(fn->index(vi)));
			
				// Storing edges on queue for flip checking
				std::pair<Face_iterator, int> pair;
				if(i == 0){
					pair = std::make_pair(fc, fc->index(fn));
					queue.push(pair);
				}
				pair = std::make_pair(fn, fn->index(fc));
				queue.push(pair);
			}
		}
	
		flip_check(queue);
	}
}

void Delaunay_triangulation::remove_dim_decrease(Vertex_iterator v)
{
	std::vector<Face_iterator> to_delete;
	std::vector<Face_iterator> vector;
	Face_iterator fc; Face_iterator ff;
	Face_circulator* fh = new Face_circulator(v);
	Face_iterator done = fh->current_face();
	do{
		to_delete.push_back(fh->current_face());
		fh->next();
		if(is_infinite(fh->current_face())) vector.push_back(fh->current_face());
	}while(fh->current_face() != done);
	delete fh;

	// Updating adjacencies
	if(vector.size() > 2)
		std::cerr << "ERROR: finite vertex has more than two infinite faces.\n";
	
	fc = vector[0]; ff = vector[1];
	fc = fc->neighbor(fc->index(v)); ff = ff->neighbor(ff->index(v));
	fc->neighbor(ff, fc->index(vector[0]));
	ff->neighbor(fc, ff->index(vector[1]));

	// Deleting v's incident faces
	for(unsigned int i = 0; i < to_delete.size(); ++i){
		fc = to_delete[i];
		if(infinite_vertex()->incident_face() == fc)
			infinite_vertex()->incident_face(fc->neighbor(ccw(fc->index(infinite_vertex()))));
		if(fc == key_face())
			key_face(nullptr);
		faces().erase(fc);

	}
	vertices().erase(v);

	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		int ii = it->index(infinite_vertex());
		it->neighbor(it, ii);
		fc = it;
	}
	infinite_vertex()->incident_face(fc);
	
	// Updating vertices incident faces
	fh = new Face_circulator(infinite_vertex());
	done = fh->current_face();
	do{
		int ii = fh->current_face()->index(infinite_vertex());
		Vertex_iterator vi = fh->current_face()->vertex(ccw(ii));
		vi->incident_face(fh->current_face());
		fh->next();
		if(fh->current_face() == done){
			fh->previous();
			vi = fh->current_face()->vertex(cw(ii));
			vi->incident_face(fh->current_face());
			fh->next();
		}
	}while(fh->current_face() != done);
	delete fh;

	dimension(dimension() - 1);

}

//				     ================================ 					 //
//                                         GEOMETRIC PREDICATES                                          //
//				     ================================ 					 //


// Return value > 0 if vd lies inside the circumcircle, value < 0 lies outside the circumcircle 
// and value == 0 if points are cocircular
double Delaunay_triangulation::incircle_test(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc, Vertex_iterator vd)
{
	return incircle_test(va->point(), vb->point(), vc->point(), vd->point());
}

double Delaunay_triangulation::incircle_test(Point& a, Point& b, Point& c, Point& d)
{
	double pa[2]; double pb[2]; double pc[2]; double pd[2];

	pa[0] = a.get_x(); pa[1] = a.get_y();
	pb[0] = b.get_x(); pb[1] = b.get_y();
	pc[0] = c.get_x(); pc[1] = c.get_y();
	pd[0] = d.get_x(); pd[1] = d.get_y();

	return incircle(pa, pb, pc, pd);

}

//					  =====================						 //
//                                              UTILITIES                                                //
//					  =====================						 //

void Delaunay_triangulation::conflict_region(Face_iterator fc, Point& pp, std::vector<Face_iterator>& region, std::vector<std::pair<Face_iterator, int>>& cavity, Vertex_location lc, int li)
{
	std::stack<Face_iterator> stack;
	Face_iterator ff = nullptr;
	Face_iterator fn = nullptr;
	Face_iterator key = nullptr;
	Vertex_iterator v0; Vertex_iterator v1; Vertex_iterator v2;
	Point p0; Point p1; Point p2; 
	// if a point lies on the convex hull boundary and infinite face need to be added to cavity
	bool boundary_edge = false;

	if(lc == ON_EDGE){
		if(is_infinite(fc->neighbor(li))){
			boundary_edge = true;
			key = fc->neighbor(li);
		}
	}
	
	std::map<Face_iterator, bool> face_visited;
	stack.push(fc);
	face_visited[fc] = true;
	region.push_back(fc);
	while(!stack.empty()){
		ff = stack.top();
		stack.pop();
		face_visited[ff] = true;

		for (int i = 0; i < dimension() + 1; ++i){
			fn = ff->neighbor(i);
			if(!is_infinite(fn)){
				p0 = fn->vertex(0)->point();
       		     	   	p1 = fn->vertex(1)->point();
       		       		p2 = fn->vertex(2)->point();

				if(incircle_test(p0, p1, p2, pp) > 0){
					if(!face_visited[fn]){
						region.push_back(fn);
						stack.push(fn);
					}
				}
				else{
					std::pair<Face_iterator, int> pair = std::make_pair(ff, i);
					cavity.push_back(pair);
				}
			}
			else{ // infinite face, then boundary edge found --> ??
				if(fn != key){
	 				int pos = fn->index(infinite_vertex());
					p0 = fn->vertex(ccw(pos))->point();
					p1 = fn->vertex(cw(pos))->point();
	
					if(orientation_test(p0, p1, pp) > 0){
						if(!face_visited[fn]){
							region.push_back(fn);
							stack.push(fn);
						}
					}
					else{
						std::pair<Face_iterator, int> pair = std::make_pair(ff, i);
						cavity.push_back(pair);
					}
				}
			}
		}
	}

	if(boundary_edge){
		Face_iterator ff = fc->neighbor(li);
		region.push_back(ff);
		int i = ff->index(infinite_vertex());
		std::pair<Face_iterator, int> pair1 = std::make_pair(ff, cw(i));
		std::pair<Face_iterator, int> pair2 = std::make_pair(ff, ccw(i));
		cavity.push_back(pair1); cavity.push_back(pair2);
	}
}

void Delaunay_triangulation::flip_check(std::queue<std::pair<Face_iterator, int>>& queue)
{



	typedef std::pair<Face_iterator, int> Pair;
	std::map<Face_iterator, bool> face_visited;
	std::map<int, bool> index_visited;

	while(!queue.empty()){
		Pair pair = queue.front();
		queue.pop();

		Face_iterator fc = std::get<0>(pair);
		int i = std::get<1>(pair);
		Face_iterator fn = fc->neighbor(i); int j = mirror_index(fc, i);//fn->index(fc);
					
		if(is_infinite(fc) || is_infinite(fn)){
			if(is_infinite(fc)){
				int ii = fc->index(infinite_vertex());
	
				if(fn->vertex(j) == infinite_vertex()) continue;
				double o = orientation_test(fc->vertex(ccw(ii)), fc->vertex(cw(ii)), fn->vertex(j));
				if(o >= 0){
					if(o == 0 && is_infinite(fc) && is_infinite(fn)) continue;
					Pair pair1 = std::make_pair(fc->neighbor(ccw(i)), fc->neighbor(ccw(i))->index(fc));
					Pair pair2 = std::make_pair(fc->neighbor(cw(i)), fc->neighbor(cw(i))->index(fc));
					Pair pair3 = std::make_pair(fn->neighbor(ccw(j)), fn->neighbor(ccw(j))->index(fn));
					Pair pair4 = std::make_pair(fn->neighbor(cw(j)), fn->neighbor(cw(j))->index(fn));
	
					queue.push(pair1); queue.push(pair2); queue.push(pair3); queue.push(pair4);

					flip(fc, i);
				}
			}
			else{
				int ii = fn->index(infinite_vertex());
				j = mirror_index(fn, ii);

				if(fc->vertex(j) == infinite_vertex()) continue;
				double o = orientation_test(fn->vertex(ccw(ii)), fn->vertex(cw(ii)), fc->vertex(j));
				if(o >= 0){
					if(o == 0 && is_infinite(fc) && is_infinite(fn)) continue;	
					Pair pair1 = std::make_pair(fc->neighbor(ccw(i)), fc->neighbor(ccw(i))->index(fc));
					Pair pair2 = std::make_pair(fc->neighbor(cw(i)), fc->neighbor(cw(i))->index(fc));
					Pair pair3 = std::make_pair(fn->neighbor(ccw(j)), fn->neighbor(ccw(j))->index(fn));
					Pair pair4 = std::make_pair(fn->neighbor(cw(j)), fn->neighbor(cw(j))->index(fn));
	
					queue.push(pair1); queue.push(pair2); queue.push(pair3); queue.push(pair4);

					flip(fc, i);
				}
			}
		}
		else if(!is_infinite(fc) && !is_infinite(fn)){
			if(incircle_test(fc->vertex(0), fc->vertex(1), fc->vertex(2), fn->vertex(j)) > 0){
				Pair pair1 = std::make_pair(fc->neighbor(ccw(i)), fc->neighbor(ccw(i))->index(fc));
				Pair pair2 = std::make_pair(fc->neighbor(cw(i)), fc->neighbor(cw(i))->index(fc));
				Pair pair3 = std::make_pair(fn->neighbor(ccw(j)), fn->neighbor(ccw(j))->index(fn));
				Pair pair4 = std::make_pair(fn->neighbor(cw(j)), fn->neighbor(cw(j))->index(fn));

				queue.push(pair1); queue.push(pair2); queue.push(pair3); queue.push(pair4);

				flip(fc, i);
			}
		}
	}
}




/////////////////////////////////// PARALLEL //////////////////////////////////

void Delaunay_triangulation::parallel_insert(std::vector<Point>& p)
{
	Spatial_sort_traits sst;
	CGAL::spatial_sort(p.begin(), p.end(), sst);
	
	int sequential_size = std::floor(p.size() * 0.3);
	for(int i = 0; i < sequential_size; ++i){
		insert(p[i]);
	}

	if(dimension() != 2){
		std::cerr << "ERROR: Dimension should be 2 on parallel insertion.\n";
		return;
	}

	int tsize[NUMBER_OF_THREADS];
	int tindex[NUMBER_OF_THREADS + 1];
	int work_size = p.size() - sequential_size;
	int remainder = work_size % NUMBER_OF_THREADS;
	std::thread t[NUMBER_OF_THREADS];
	std::queue<Point> queue[NUMBER_OF_THREADS];
	Vertex_iterator start_vertex[NUMBER_OF_THREADS];
	work_size = work_size/NUMBER_OF_THREADS;
	for(int i = 0; i < NUMBER_OF_THREADS; ++i){
		start_vertex[i] = nullptr;
		tsize[i] = work_size + ((i < remainder) ? 1 : 0);
		tindex[i] = (i == 0 ? sequential_size : (tindex[i-1] + tsize[i-1]) );
	}
	tindex[NUMBER_OF_THREADS] = p.size();	

	// Interleaving is the number of chunks the work per thread will be subdivided
	int interleaving = 1;
	if(work_size >= 24)
		interleaving = 4; // USING INTERLEAVING 1
	else if(work_size >= 12)
		interleaving = 2;
	else
		interleaving = 1;
	
	Work_chunk workload[interleaving * NUMBER_OF_THREADS];
	std::queue<Work_chunk> work_queue[NUMBER_OF_THREADS];
		int count = 0;
	int index = 0;
	for(int tid = 0; tid < NUMBER_OF_THREADS; ++tid){
		int wsize = tsize[tid]/interleaving;
		int wremainder = tsize[tid] % interleaving;
		for(int i = 0; i < interleaving; ++i){
			workload[index].size = wsize + (i < wremainder ? 1 : 0);
			workload[index].beginning = ((index == 0) ? sequential_size : workload[index-1].end);
			workload[index].current = workload[index].beginning;
			workload[index].end = workload[index].beginning + workload[index].size;
			work_queue[tid].push(workload[index]);
			count += workload[index].size;
			index++;
		}

	}

	// Begin of Parallel region.
	for(int i = 0; i < NUMBER_OF_THREADS; ++i)
		t[i] = std::thread([&, i](){parallel_try_insert(p, queue, start_vertex, tindex, i, work_queue[i]);});

	for(int i = 0; i < NUMBER_OF_THREADS; ++i)
		t[i].join();
	// End of Parallel region.
	
#ifndef NDEBUG
#ifdef DETAILS
	std::cout << "\nPrinting Statistics\n\n";
	std::cout << "Locate: " << elapsed_locate.count() << "\n";
	std::cout << "Conflict: " << elapsed_conflict.count() << "\n";
	std::cout << "Update: " << elapsed_update.count() << "\n";
	std::cout << "\n=================================\n\n";
	std::cout << "Locate retreats: " << locate_retreats << "\n";
	std::cout << "Conflict retreats: " << conflict_retreats << "\n";
	std::cout << "\n=================================\n\n";
#endif
#endif

}

void Delaunay_triangulation::parallel_try_insert(std::vector<Point>& p, std::queue<Point>* queue, Vertex_iterator* start_vertex, int* tindex, int tid, 
	std::queue<Work_chunk>& work_queue)
{	
	std::vector<Face_iterator> region;
	std::vector<std::pair<Face_iterator, int>> cavity;
	Vertex_location lc;
	Face_iterator fc;
	Point tpoint;
	int li;


	while(!work_queue.empty()){
		Work_chunk wk = work_queue.front();
		work_queue.pop();
		tpoint = p[wk.current];

		// LOCATION
	#ifndef NDEBUG
		start_locate = std::chrono::high_resolution_clock::now();
	#endif
		fc = parallel_locate(tpoint, start_vertex[tid], lc, li);
	#ifndef NDEBUG
		end_locate = std::chrono::high_resolution_clock::now();
		elapsed_locate += (end_locate - start_locate);
	#endif
		if(fc == nullptr){
			work_queue.push(wk);
		#ifndef NDEBUG
			locate_retreats++;
		#endif
		}
		else{
			// CONFLICT REGION
		#ifndef NDEBUG
			start_conflict = std::chrono::high_resolution_clock::now();
		#endif
			bool conflict_done = parallel_conflict_region(fc, tpoint, region, cavity, lc, li);
		#ifndef NDEBUG
			end_conflict = std::chrono::high_resolution_clock::now();
			elapsed_conflict += (end_conflict - start_conflict);
		#endif
			if(!conflict_done){
				region.clear();
				cavity.clear();
				work_queue.push(wk);
			#ifndef NDEBUG
				conflict_retreats++;
			#endif
			}
			else{
				// UPDATE
			#ifndef NDEBUG
				start_update = std::chrono::high_resolution_clock::now();
			#endif
				parallel_insert_in_hole(fc, tpoint, lc, li, region, cavity);
			#ifndef NDEBUG
				end_update = std::chrono::high_resolution_clock::now();
				elapsed_update += (end_conflict - start_conflict);
			#endif
				wk.current++;
				if(!wk.is_done())
					work_queue.push(wk);
			}
		}
	}
}

Face_iterator Delaunay_triangulation::parallel_locate(Point& p, Vertex_iterator start_vertex, Vertex_location& location, int& li)
{
	// Using a vertex as beginning because faces can be deleted on update phase.
	if(start_vertex == nullptr)
		start_vertex = ++vertices().begin();

	Face_iterator fc = start_vertex->incident_face();
	if(!fc->try_lock())
		return nullptr;

	// Location should start with a finite face.
	if(is_infinite(fc)){
		for(int i = 0; i <= dimension(); ++i){
			if(!is_infinite(fc->neighbor(i))){
				Face_iterator to_unlock = fc;
				fc = fc->neighbor(i);
				if(!fc->try_lock()){
					to_unlock->unlock();
					return nullptr;
				}
				to_unlock->vertex(i)->unlock();
				break;
			}
		}
	}

	while(true){
		Point& pa = fc->vertex(0)->point();
		Point& pb = fc->vertex(1)->point();
		Point& pc = fc->vertex(2)->point();

		double o1 = orientation_test(pa, pb, p);
		double o2 = orientation_test(pb, pc, p);
		double o3 = orientation_test(pc, pa, p);

		start_vertex = fc->vertex(0);

		Face_iterator adj = nullptr;
		if(o1 < 0){
			adj = fc->neighbor(2);
			fc->vertex(2)->unlock();
			if(!adj->try_lock()){
				adj->unlock();
				return nullptr;
			}
			if(is_infinite(adj)){
				location = OUTSIDE_CONVEX_HULL;
				return adj;
			}
			fc = adj;
		}
		else if(o2 < 0){
			adj = fc->neighbor(0);
			fc->vertex(0)->unlock();
			if(!adj->try_lock()){
				adj->unlock();
				return nullptr;
			}
			if(is_infinite(adj)){
				location = OUTSIDE_CONVEX_HULL;
				return adj;
			}
			fc = adj;
		}
		else if(o3 < 0){
			adj = fc->neighbor(1);
			fc->vertex(1)->unlock();
			if(!adj->try_lock()){
				adj->unlock();
				return nullptr;
			}
			if(is_infinite(adj)){
				location = OUTSIDE_CONVEX_HULL;
				return adj;
			}
			fc = adj;
		}
		else{
			if(o1 > 0 && o2 > 0 && o3 > 0)
				location = ON_INTERIOR;
			else if((o1 == 0 && o2 == 0) || (o2 == 0 && o3 == 0) || (o1 == 0 && o3 == 0))
				location = ON_VERTEX;
			else{
				location = ON_EDGE;
				if(o1 == 0)
					li = 2;
				else if(o2 == 0)
					li = 0;
				else
					li = 1;
			}

			return fc;
		}
	}
	// Should never reach this line
	return nullptr;
}

bool Delaunay_triangulation::parallel_conflict_region(Face_iterator fc, Point& pp, std::vector<Face_iterator>& region, std::vector<std::pair<Face_iterator, int>>& cavity, 
     	Vertex_location lc, int li)
{
	region.clear();
	cavity.clear();

	std::stack<Face_iterator> stack;
	Face_iterator ff;
	Face_iterator fn;
	Face_iterator key = nullptr;
	// if a point lies on the convex hull boundary and infinite face need to be added to cavity
	bool boundary_edge = false;
	bool done = true;
	std::vector<Face_iterator> locked_faces;

	// Check if fc is already locked, if not try to lock it.
	if(!fc->try_lock())
		done = false;
	else
		locked_faces.push_back(fc);

	if(lc == ON_EDGE){
		if(is_infinite(fc->neighbor(li))){
			boundary_edge = true;
			key = fc->neighbor(li);
		}
	}
	
	std::map<Face_iterator, bool> face_visited;
	stack.push(fc);
	face_visited[fc] = true;
	region.push_back(fc);
	while(!stack.empty() && done){
		ff = stack.top();
		stack.pop();
		face_visited[ff] = true;

		for (int i = 0; i < dimension() + 1; ++i){
			fn = ff->neighbor(i);
			if(!fn->try_lock()){
				done = false;
				break;
			}
			else
				locked_faces.push_back(fn);

			if(!is_infinite(fn)){
				Point& p0 = fn->vertex(0)->point();
       		     	   	Point& p1 = fn->vertex(1)->point();
       		       		Point& p2 = fn->vertex(2)->point();
				
				if(incircle_test(p0, p1, p2, pp) > 0){
					if(!face_visited[fn]){
						region.push_back(fn);
						stack.push(fn);
					}
				}
				else{
					std::pair<Face_iterator, int> pair = std::make_pair(ff, i);
					if(!(ff->neighbor(i))->try_lock()){
						done = false;
						break;
					}
					else
						locked_faces.push_back(ff->neighbor(i));

					cavity.push_back(pair);
				}
			}
			else{ // infinite face, then boundary edge found --> ??
				if(fn != key){
	 				int pos = fn->index(infinite_vertex());
					Point& p0 = fn->vertex(ccw(pos))->point();
					Point& p1 = fn->vertex(cw(pos))->point();
	
					if(orientation_test(p0, p1, pp) > 0){
						if(!face_visited[fn]){
							region.push_back(fn);
							stack.push(fn);
						}
					}
					else{
						std::pair<Face_iterator, int> pair = std::make_pair(ff, i);
						if(!(ff->neighbor(i))->try_lock()){
							done = false;
							break;
						}
						else
							locked_faces.push_back(ff->neighbor(i));

						cavity.push_back(pair);
					}
				}
			}
		}
	}

	if(boundary_edge && done){
		Face_iterator ff = fc->neighbor(li);
		region.push_back(ff);
		int i = ff->index(infinite_vertex());

		if(!(ff->neighbor(cw(i)))->try_lock())
			done = false;
		else
			locked_faces.push_back(ff->neighbor(cw(i)));

		if(!(ff->neighbor(ccw(i)))->try_lock())
			done = false;
		else
			locked_faces.push_back(ff->neighbor(ccw(i)));

		std::pair<Face_iterator, int> pair1 = std::make_pair(ff, cw(i));
		std::pair<Face_iterator, int> pair2 = std::make_pair(ff, ccw(i));
		cavity.push_back(pair1); cavity.push_back(pair2);
	}

	if(!done){
		// Unlock faces.
		for(unsigned int i = 0; i < locked_faces.size(); ++i)
			locked_faces[i]->unlock();
		return false;
	}

	// In case of region done.
	return true;
}

void Delaunay_triangulation::parallel_insert_in_hole(Face_iterator fc, Point& p, Vertex_location lc, int li, std::vector<Face_iterator> region,
	std::vector<std::pair<Face_iterator, int>> cavity)
{
	Vertex_iterator new_vertex = nullptr;
	critical_mutex.lock();
	new_vertex = vertices().insert(Vertex(p));
	critical_mutex.unlock();
	std::vector<Face_iterator> faces_created;

	for(unsigned int i = 0; i < cavity.size(); ++i){
		Face_iterator old_face = std::get<0>(cavity[i]);
		int old_int = std::get<1>(cavity[i]);

		// create new face and set vertices according to old face and new vertex to be inserted
		Face_iterator new_face = nullptr;
		critical_mutex.lock();
		new_face = faces().insert(Face());
		critical_mutex.unlock();
		new_face->vertex(old_face->vertex(cw(old_int)), cw(old_int));
		new_face->vertex(old_face->vertex(ccw(old_int)), ccw(old_int));
		new_face->vertex(new_vertex, old_int);

		//if(!new_face->try_lock())
		//	std::cerr << "ERROR: could not lock new face. (parallel insert)\n";

		// setting adjacencies for neighbors outside cavity
		Face_iterator old_neighbor = old_face->neighbor(old_int);
		old_neighbor->neighbor(new_face, old_neighbor->index(old_face));
		new_face->neighbor(old_neighbor, old_int);

		Vertex_iterator v = new_face->vertex(ccw(old_int));
		v->incident_face(new_face);
		if(i == 0)
			new_vertex->incident_face(new_face);
		faces_created.push_back(new_face);
	}

	// adjacencies among new faces
	for(unsigned int i = 0; i < faces_created.size(); ++i){
		Face_iterator fc = faces_created[i];
		int j = fc->index(new_vertex);
		Vertex_iterator vi = fc->vertex(cw(j));
		Face_iterator ff = vi->incident_face();

		if(fc == ff)
			std::cerr << "ERROR: Identical faces when setting neighbor pointers.\n";

		// indices for setting neighbors
		int m, n;
		if(fc->vertex(cw(j)) == vi)	m = ccw(j);
		else				m = cw(j);

		int k = ff->index(new_vertex);
		if(ff->vertex(cw(k)) == vi)	n = ccw(k);
		else				n = cw(k);

		fc->neighbor(ff, m);
		ff->neighbor(fc, n);
	}

	// Unlocking cavity.
	for(unsigned int i = 0; i < cavity.size(); ++i){
		Face_iterator n = std::get<0>(cavity[i]);
		int j = std::get<1>(cavity[i]);

		(n->neighbor(j))->unlock();
	}

	// removing old faces from the conflict region
	for(unsigned int i = 0; i < region.size(); ++i){
		Face_iterator to_remove = region[i];
		critical_mutex.lock();
		faces().erase(to_remove);
		critical_mutex.unlock();
		if(key_face() == to_remove)
			key_face(nullptr);
	}
	
	//return new_vertex;
}
