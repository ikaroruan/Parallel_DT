#ifndef VERTEX_H
#define VERTEX_H

#include "Point.h"
#include <CGAL/Compact_container.h>

class Vertex;
class Face;

typedef CGAL::Compact_container<Face>	Faces_container;
typedef CGAL::Compact_container<Vertex>	Vertices_container;
typedef Faces_container::iterator	Face_iterator;
typedef Vertices_container::iterator 	Vertex_iterator;

class Vertex : public CGAL::Compact_container_base
{
	public:
	Vertex() : _point(Point()), _incident_face(nullptr) {}
	Vertex(const Point& p) : _point(p), _incident_face(nullptr) {}
	Vertex(Face_iterator f) : _point(Point()), _incident_face(f) {}
	Vertex(const Point& p, Face_iterator f) : _point(p), _incident_face(f) {}
	Vertex(const Vertex& v) : _point(v._point), _incident_face(v._incident_face) {}
	Vertex& operator=(const Vertex& v);
	bool operator==(const Vertex& v);
	bool operator!=(const Vertex& v);

	void point(Point& p);
	void incident_face(Face_iterator f);
	Point& point();
	Face_iterator incident_face();

	//void* for_compact_container() const {return _p;}
	//void *& for_compact_container() {return _p;}

	private:
	Point _point;
	Face_iterator _incident_face;
	//void* _p;
};

#endif
