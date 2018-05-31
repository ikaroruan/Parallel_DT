#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <mutex>
#include "Lock.h"
#include "Point.h"

#include <CGAL/Compact_container.h>

class Face;
class Vertex;

typedef CGAL::Compact_container<Face>	Faces_container;
typedef CGAL::Compact_container<Vertex>	Vertices_container;
typedef CGAL::Compact_container<Face>::iterator	Face_iterator;
typedef CGAL::Compact_container<Vertex>::iterator Vertex_iterator;

class Vertex : public CGAL::Compact_container_base
{
	public:

	Vertex();
	Vertex(Point& p);
	Vertex(const Vertex& v) : _point(v._point), _in_face(v._in_face) {}
	Vertex& operator=(const Vertex& v);
	bool operator==(const Vertex& v);
	bool operator!=(const Vertex& v);

	void point(Point& p);
	Point& point();
	Face_iterator incident_face();
	void incident_face(Face_iterator fc);

	bool lock();
	bool try_lock();
	bool unlock();
	bool is_locked();
	
	private:
	Point _point;
	Face_iterator _in_face;
	Lock _lck;
};

#endif
