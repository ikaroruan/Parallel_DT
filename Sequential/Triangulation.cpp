#include "Triangulation.h"

//		 ==========================================================================		//
//                     CONSTRUCTOR, DESTRUCTOR AND PRIVATES VARIABLE ACCESS/MODIFIERS                   //
//		 ==========================================================================		//

Triangulation::Triangulation()
{
	_infinite_vertex_iterator = vertices().insert(_infinite_vertex); 
	_key_face = nullptr;
	_dimension = -1;
}

Triangulation::~Triangulation()
{

}

int Triangulation::number_of_vertices()
{
	return vertices().size();
}


int Triangulation::dimension() const
{
	return _dimension;
}

void Triangulation::dimension(int i)
{
	_dimension = i;
}

Vertex_iterator Triangulation::infinite_vertex()
{
	return _infinite_vertex_iterator;
}

Vertex_iterator Triangulation::infinite_vertex_iterator()
{
	return _infinite_vertex_iterator;
}

Vertices_container& Triangulation::vertices()
{
	return _vertices;
}

Faces_container& Triangulation::faces()
{
	return _faces;
}

void Triangulation::key_face(Face_iterator fc)
{
	_key_face = fc;
}

Face_iterator Triangulation::key_face()
{
	return _key_face;
}

//					  ===================== 					 //
//                                              UTILITIES                                                //
//					  ===================== 					 //

int Triangulation::ccw(int i)
{
	if(i < 0 || i > 2){
		std::cerr << "ERROR: wrong index for ccw function.\n";
		return -1; // reporting an error
	}
	if(i == 2)
		return 0;
	return i+1;
}

int Triangulation::cw(int i)
{
	if(i < 0 || i > 2){
		std::cerr << "ERROR: wrong index for cw function.\n";
		return -1;
	}
	if(i == 0)
		return 2;
	return i-1;
}

int Triangulation::number_of_faces()
{
	return faces().size();
}

bool Triangulation::is_infinite(Face_iterator fc)
{
	if(fc->vertex(0) == infinite_vertex_iterator() || fc->vertex(1) == infinite_vertex_iterator() ||
		fc->vertex(2) == infinite_vertex_iterator()){
		return true;}
	return false;
}

bool Triangulation::is_equal(Point& p, Point& q)
{
	return ((p.get_x() == q.get_x()) && (p.get_y() == q.get_y()));
}

Vertex_iterator Triangulation::random_vertex()
{
	boost::random::mt19937 rgn;
	boost::random::uniform_int_distribution<> range(0, vertices().size() - 1);
	
	Vertex_iterator vi = vertices().begin();
	do{
		int current = 0;
		int stop = range(rgn);
		while(current < stop){
			++vi;
			++current;
		}
	}while(vi == infinite_vertex_iterator());
	return vi;
}

int Triangulation::mirror_index(Face_iterator fc, int i)
{
	return ccw(fc->neighbor(i)->index(fc->vertex(ccw(i))));
}

// Tests whether q is collinear between p and r
bool Triangulation::collinear_between(Point& p, Point& q, Point& r)
{
	// Precodition that p and r are supposed to be collinear
	double px = p.get_x(); double py = p.get_y();
	double qx = q.get_x(); double qy = q.get_y();
	double rx = r.get_x(); double ry = r.get_y();	

	double d1, d2, dd;

	if(px == qx){
		d1 = std::abs(qy - py);
		d2 = std::abs(ry - qy);
		dd = std::abs(ry - py);
	}else{
		d1 = std::abs(qx - px);
		d2 = std::abs(rx - qx);
		dd = std::abs(rx - px);
	}

	return ((dd > d1) && (dd > d2));
}

// Check if dimension goes from 2 to 1
bool Triangulation::dim_goes_down(Vertex_iterator v)
{
	// Check whether all finite faces have the given vertex to be removed.
	// If one face does not contain, then dimension remains the same.
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		if(!is_infinite(it)){
			if( ! it->contains(v))
				return false;
		}
	}

	// If all finite faces have v, then checking if all vertices except v are collinear is needed.
	Face_iterator fi = (faces().begin());
	while(!fi->contains(v) || is_infinite(fi)) 
		++fi;
	Face_iterator fc = fi;
	int i = fc->index(v);
	Vertex_iterator va = fc->vertex(ccw(i)); Vertex_iterator vb = fc->vertex(cw(i));
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		if(!is_infinite(it)){
			i = it->index(v);
			Vertex_iterator vc = it->vertex(ccw(i));
			Vertex_iterator vd = it->vertex(cw(i));
			

			if(!vertices_are_collinear(va, vb, vc))
				return false;
			if(!vertices_are_collinear(va, vb, vd))
				return false;
		}
	}

	return true;
}

void Triangulation::flip(Face_iterator fc, int i)
{
	Face_iterator ff = fc->neighbor(i);
	int ni = mirror_index(fc, i);
	Vertex_iterator v_cw = fc->vertex(cw(i));
	Vertex_iterator v_ccw = fc->vertex(ccw(i));

	Face_iterator nc = fc->neighbor(ccw(i));
	int nci = mirror_index(fc, ccw(i));
	Face_iterator nf = ff->neighbor(ccw(ni));
	int nfi = mirror_index(ff, ccw(ni));

	fc->vertex(ff->vertex(ni), cw(i));
	ff->vertex(fc->vertex(i), cw(ni));

	// maintaining neighborhood
	fc->neighbor(nf, i); nf->neighbor(fc, nfi);
	ff->neighbor(nc, ni); nc->neighbor(ff, nci);
	ff->neighbor(fc, ccw(ni)); fc->neighbor(ff, ccw(i));

	// check vertices incident faces
	if(v_cw->incident_face() == fc)
		v_cw->incident_face(ff);
	if(v_ccw->incident_face() == ff)
		v_ccw->incident_face(fc);
}

Vertex_iterator Triangulation::is_vertex(Face_iterator fc, Point p)
{
	for(int i = 0; i < 3; ++i){
		double vx, vy, px, py;
		Point vp = fc->vertex(i)->point();
		vx = vp.get_x(); vy = vp.get_y();
		px = p.get_x(); py = p.get_y();

		if(px == vx && py == vy)
			return fc->vertex(i);
	}

	std::cerr << "ERROR: face does not contain the given vertex. (is_vertex function)\n";
	return nullptr;
}	

//			    ===================================================				 //
//                                GEOMETRIC PREDICATES AND POINT LOCATION 	                         //
//			    ===================================================				 //

// Returns value > 0  if they are in counterclockwise order, value < 0 if clockwise order 
// and value == 0 if they are collinear

double Triangulation::orientation_test(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc)
{
	return orientation_test(va->point(), vb->point(), vc->point());
}

double Triangulation::orientation_test(Point& a, Point& b, Point& c)
{
	double pa[2]; double pb[2]; double pc[2];

	pa[0] = a.get_x(); pa[1] = a.get_y();
	pb[0] = b.get_x(); pb[1] = b.get_y();
	pc[0] = c.get_x(); pc[1] = c.get_y();

	return orient2d(pa, pb, pc);
	
}

bool Triangulation::vertices_are_collinear(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc)
{	
	if(orientation_test(va, vb, vc) == 0)
		return true;
	return false;
}

// In this case, key_face is the last finite face returned by locate function. Otherwise, it should be nullptr.
Face_iterator Triangulation::locate(Point& p, Vertex_location& location, int& li)
{	
	return locate(p, key_face(), location, li);
}

Face_iterator Triangulation::locate(Point& p, Face_iterator fc, Vertex_location& location, int& li)
{
	li = -1; // used for indicating edge location for dimension => 2
	if(dimension() < 0){ // only the infinite vertex exists
		location = OUTSIDE_CONVEX_HULL;
		return nullptr;
	}
	if(dimension() == 0){ // triangulation has another vertex rather than the infinite vertex
		Vertex_iterator v0 = ++Vertices_begin(); // first vertex excepting infinite vertex
		if(p.get_x() == v0->point().get_x() && p.get_y() == v0->point().get_y())
			location = ON_VERTEX;
		else
			location = OUTSIDE_CONVEX_HULL;
		return faces().begin();
	}
	if(dimension() == 1){ // there are at least two vertices and the infinite vertex
		Face_iterator fn = faces().begin();
		int vi = fn->index(infinite_vertex());
		Point& p0 = fn->vertex(ccw(vi))->point(); 	
		Point& p1 = fn->vertex(cw(vi))->point();
		Face_iterator fmin = nullptr; Vertex_iterator vmin = nullptr; double dmin = -1;
		Face_iterator fmax = nullptr; Vertex_iterator vmax = nullptr; double dmax = -1;
		
		bool on_edge = (orientation_test(p0, p1, p) == 0);
		
		Face_circulator* fh = new Face_circulator(infinite_vertex());
		Face_iterator done = fh->current_face();
		do{
			int i = fh->current_face()->index(fh->current_vertex()); // current vertex == infinite vertex
			Vertex_iterator u = fh->current_face()->vertex(ccw(i));
			Vertex_iterator w = fh->current_face()->vertex(cw(i));
			
			if(is_equal(u->point(), p)){
				location = ON_VERTEX;
				li = ccw(i);
				fn = fh->current_face();
				delete fh;
				return fn;
			}
			if(is_equal(w->point(), p)){
				location = ON_VERTEX;
				li = cw(i);
				fn = fh->current_face();
				delete fh;
				return fn;
			}
			
			if(collinear_between(u->point(), p, w->point())){
				if(on_edge){
					location = ON_EDGE;
					li = -1;
					fn = fh->current_face();
					delete fh;
					return fn;
				}
			}

			double d1, d2;

			// FIXME: STRANGE NOTATION.
			if(u->point().get_x() == w->point().get_x()){
				d1 = std::abs( u->point().get_y() - p.get_y() );
				d2 = std::abs( w->point().get_y() - p.get_y() );
			}else{
				d1 = std::abs( u->point().get_x() - p.get_x() );
				d2 = std::abs( w->point().get_x() - p.get_x() );
			}
				
			if(dmin == -1){ // First face. Set dmin and dmax for comparison.
				if(d1 < d2){
					dmin = d1; fmin = fh->current_face(); vmin = u;
					dmax = d2; fmax = fh->current_face(); vmax = w;
				}
				else{
					dmin = d2; fmin = fh->current_face(); vmin = w;
					dmax = d1; fmax = fh->current_face(); vmax = u;
				}
			}
 				
			// Updating dmin.
			if(d1 < dmin){
				dmin = d1; fmin = fh->current_face(); vmin = u;
			}
			if(d2 < dmin){
				dmin = d2; fmin = fh->current_face(); vmin = w;
			}

			// Updating dmax
			if(d1 > dmax){
				dmax = d1; fmax = fh->current_face(); vmax = u;
			}
			if(d2 > dmax){
				dmax = d2; fmax = fh->current_face(); vmax = w;
			}
			fh->next();
		}while(fh->current_face() != done);
		delete fh;

		if(on_edge){
			location = ON_EDGE;
			li = fmin->index(vmin);
			//key_face = fmin;
			return fmin;
		}
		else{
			location = OUTSIDE_CONVEX_HULL;
			li = fmax->index(vmax);
			//key_face = fmax;
			return fmax;
		}
	
		return nullptr;
	}

	// if dimension => 2
	// Finite face neighboring the infinite one
	if(fc == nullptr){ // fc should be a finite face for testing orientation
		//std::cout << "Got here.\n";
		fc = faces().begin();
		if(is_infinite(fc)){
			if(!is_infinite((fc)->neighbor(0)))
				fc = (fc)->neighbor(0);
			else if(!is_infinite((fc)->neighbor(1)))
				fc = (fc)->neighbor(1);
			else
				fc = (fc)->neighbor(2);

		}
		if(is_infinite(fc))
			std::cerr << "\nERROR: Face location is infinite.\n";
			
	}

	while(true){
		Point& pa = fc->vertex(0)->point();
		Point& pb = fc->vertex(1)->point();
		Point& pc = fc->vertex(2)->point();

		double o1 = orientation_test(pa, pb, p);
		double o2 = orientation_test(pb, pc, p);
		double o3 = orientation_test(pc, pa, p);

		Face_iterator adj = nullptr;
		//Face_iterator fc_temp = nullptr;
		if(o1 < 0){
			adj = fc->neighbor(2);
			if(is_infinite(adj)){
				location = OUTSIDE_CONVEX_HULL;
				return adj;
			}
			fc = adj;
		}
		else if(o2 < 0){
			adj = fc->neighbor(0);
			if(is_infinite(adj)){
				location = OUTSIDE_CONVEX_HULL;
				return adj;
			}
			fc = adj;
		}
		else if(o3 < 0){
			adj = fc->neighbor(1);
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

			key_face(fc);
			return fc;
		}
	}
	// Should never reach this line
	return nullptr;
}


//				     ==============================		     			 //
//                                         INPUTS AND OUTPUTS                                            //
//				     ============================== 					 //

void Triangulation::show_triangulation()
{
	show_triangulation(std::string());
}

// Based on showme *.ele file
// Do not show the infinite vertex on the visualization
void Triangulation::show_triangulation(std::string filename)
{
	std::ofstream node;
	std::ofstream ele;
	if(filename.size() <= 0){
		node.open("out.node");
		ele.open("out.ele");
	}
	else{
		std::string nodename = filename.substr(0, filename.find_last_of("."));
		std::string elename = nodename + "_out.ele";
		nodename = nodename + "_out.node";

		node.open(nodename.c_str());
		ele.open(elename.c_str());

	}
	std::map<Vertex_iterator, int> vertex_position;

	// *.node file
	node << vertices().size() - 1 << "\t2\t0\t0";
	int num = 1;
	for(Vertex_iterator it = Vertices_begin(); it != Vertices_end(); ++it){
		if(it != infinite_vertex()){
			node << "\n" << num << "\t" << it->point().get_x()
				<< "\t" << it->point().get_y();
			vertex_position[it] = num;
			num++;
		}
	}

	// *.ele file
	int count = 0;
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		if(!is_infinite(it))
			count++;
	}

	num = 1;
	ele << count << "\t3\t0";
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		if(!is_infinite(it)){
			ele << "\n" << num 
				<< "\t" << vertex_position[ it->vertex(0) ]
				<< "\t" << vertex_position[ it->vertex(1) ] 
				<< "\t" << vertex_position[ it->vertex(2) ];
			num++;
		}
	}
}

// Generate a triangulation output with vertices, faces and neighbors adjacencies (based on CGAL output)
void Triangulation::output_triangulation()
{
	output_triangulation(std::string());
}

void Triangulation::output_triangulation(std::string filename)
{
	std::ofstream output;
	if(filename.size() <= 0){
		output.open("out.tri");
	}
	else{
		std::string nodename = filename.substr(0, filename.find_last_of("."));
		nodename = nodename + "_out.tri";

		output.open(nodename.c_str());
	}

	int dim = dimension(); int count = 0;
	std::map<Vertex_iterator, int> vertices_order;
	std::map<Face_iterator, int> faces_order;

	output << vertices().size() << " " << faces().size() << " " << dim;

	vertices_order[infinite_vertex()] = count++;
	for(Vertex_iterator it = Vertices_begin(); it != Vertices_end(); ++it){
		if(it != infinite_vertex()){
			vertices_order[it] = count++;
			output << "\n" << it->point().get_x() << " " << it->point().get_y();
		}
	}
	
	output << "\n";

	count = 0;
	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		faces_order[it] = count++;
		output << "\n";
		for(int i = 0; i < 3; ++i)
			output << vertices_order[it->vertex(i)] << " ";
	}

	output << "\n";

	for(Face_iterator it = Faces_begin(); it != Faces_end(); ++it){
		output << "\n";
		for(int i = 0; i < 3; ++i)
			output << faces_order[it->neighbor(i)] << " ";
	}
}

// Reads a triangulation from a CGAL output
Face_iterator Triangulation::create_triangulation(std::istream& in)
{
	int n, m, dim;
	in >> n >> m >> dim;
	
	Vertex_iterator vertices_vector[n];
	Face_iterator face[m];

	dimension(dim);
	
	vertices_vector[0] = infinite_vertex();

	double x, y;
	for(int i = 1; i < n; ++i){
		in >> x >> y;
		Vertex vi(Point(x,y));
		Vertex_iterator vit = vertices().insert(vi);
		vertices_vector[i] = vit;
	}

	int a, b, c;
	for(int i = 0; i < m; ++i){
		in >> a >> b >> c;
			// FIXME: How about this notation.
			Face fc(vertices_vector[a], vertices_vector[b], vertices_vector[c]);
			Face_iterator fcit = faces().insert(fc);
			face[i] = fcit;
			
			if(vertices_vector[a]->incident_face() == nullptr)
				vertices_vector[a]->incident_face(fcit);
			if(vertices_vector[b]->incident_face() == nullptr)
				vertices_vector[b]->incident_face(fcit);
			if(vertices_vector[c]->incident_face() == nullptr)
				vertices_vector[c]->incident_face(fcit);
	}

	for(int i = 0; i < m; ++i){
		in >> a >> b >> c;
			face[i]->neighbor(face[a], 0);
			face[i]->neighbor(face[b], 1);
			face[i]->neighbor(face[c], 2);
	}	
	
	// returning 'random' face
	boost::random::mt19937 rgn;
	boost::random::uniform_int_distribution<> range(0, m - 1);

	Face_iterator ff ;
	do{
		ff = face[range(rgn)];
	}while(is_infinite(ff));

	return ff;
}

