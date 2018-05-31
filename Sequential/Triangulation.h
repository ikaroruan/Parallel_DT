#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <list>
#include <queue>
#include <stack>
#include <vector>
#include <map> 
#include <fstream>
#include <utility>
#include <iostream>

#include "Vertex.h"
#include "Face.h"
#include "Face_circulator.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <CGAL/spatial_sort.h>
#include <CGAL/Compact_container.h>

extern "C"{
	// get exact predicates from a C code
	#include "predicates.h"
}

// Geometric traits for CGAL spatial sort
struct Less_x{
	bool operator()(const Point& p, const Point& q) const{
		return (p.get_x() < q.get_x());
	}
};

struct Less_y{
	bool operator()(const Point& p, const Point& q) const{
		return (p.get_y() < q.get_y());
	}
};

struct Spatial_sort_traits{
	typedef Point Point_2;
	typedef Less_x Less_x_2;
	typedef Less_y Less_y_2;

	Less_x_2 less_x_2_object() const{
		return Less_x_2();
	}

	Less_y_2 less_y_2_object() const{
		return Less_y_2();
	}
};

typedef enum Vertex_location{
	ON_VERTEX, ON_EDGE, ON_INTERIOR, OUTSIDE_CONVEX_HULL
} Vertex_location;

/*typedef CGAL::Compact_container<Face>	Faces_container;
typedef CGAL::Compact_container<Vertex> Vertices_container;
typedef Faces_container::iterator	Faces_iterator;
typedef Vertices_container::iterator	Vertex_iterator;*/

class Triangulation
{
	public:

	// CONSTRUCTOR AND DESTRUCTOR
	Triangulation();
	~Triangulation();

	// ITERATORS
	Vertex_iterator Vertices_begin(){return _vertices.begin();}
	Vertex_iterator Vertices_end(){return _vertices.end();}
	Face_iterator  Faces_begin(){return _faces.begin();}
	Face_iterator  Faces_end(){return _faces.end();}

	// PRIVATE VARIABLES ACCESS AND MODIFIERS
	int  dimension() const;
	void dimension(int i);
	int  number_of_vertices();
	int  number_of_faces();
	void key_face(Face_iterator fc);
	Face_iterator key_face();
	Vertex_iterator infinite_vertex();
	Vertex_iterator infinite_vertex_iterator();
	Vertices_container& vertices();
	Faces_container& faces();
	
	// UTILITIES
	int ccw(int i);
	int cw(int i);
	void flip(Face_iterator fc, int i);
	int mirror_index(Face_iterator fc, int i);
	bool dim_goes_down(Vertex_iterator v);
	bool is_infinite(Face_iterator fc);
	bool collinear_between(Point& p, Point& q, Point& r);
	bool vertices_are_collinear(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc);
	bool is_equal(Point& p, Point& q);
	Vertex_iterator random_vertex();
	Vertex_iterator is_vertex(Face_iterator fc, Point p);

	// GEOMETRIC PREDICATES AND POINT LOCATION
	double orientation_test(Point& a, Point& b, Point& c);
	double orientation_test(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc);
	Face_iterator locate(Point& p, Face_iterator fc, Vertex_location& location, int& li);
	Face_iterator locate(Point& p, Vertex_location& location, int& li);

	// TRIANGULATION INPUT AND OUTPUTS
	Face_iterator create_triangulation(std::istream& in);
	void show_triangulation();
	void show_triangulation(std::string filename);
	void output_triangulation();
	void output_triangulation(std::string filename);

	private:

	// PRIVATE MEMBER VARIABLES
	Vertices_container _vertices;
	Faces_container _faces;
	int _dimension;
	Vertex _infinite_vertex;
	Vertex_iterator _infinite_vertex_iterator;
	Face_iterator _key_face;
	
};

#endif // TRIANGULATION_H
