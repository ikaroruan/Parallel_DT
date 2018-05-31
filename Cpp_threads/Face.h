#ifndef FACE_H
#define FACE_H

#include "Vertex.h"
#include <mutex>
#include <utility>
#include <iostream>
extern "C"{
	#include "predicates.h"
}

class Face : public CGAL::Compact_container_base
{
	public: 
	Face();
	Face(Vertex_iterator p0, Vertex_iterator p1, Vertex_iterator p2);
	Face(Face_iterator f);
	~Face();

	Face_iterator neighbor(int i);
	Vertex_iterator vertex(int i);
	void neighbor(Face_iterator neighbor, int i);
	void vertex(Vertex_iterator v, int i);
	int index(Face_iterator fc);
	int index(Vertex_iterator vi);
	bool contains(Vertex_iterator v);
	bool is_neighbor(Face_iterator fc);
	bool is_vertex(Vertex_iterator v);
	void reorient();

	bool lock();
	bool try_lock();
	void unlock();
	bool any_vertex_locked();

	private:
	Face_iterator _neighbor[3]; 
	Vertex_iterator _vertices[3];
};

#endif
