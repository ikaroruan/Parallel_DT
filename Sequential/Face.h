#ifndef FACE_H
#define FACE_H

#include <iostream>
#include "Vertex.h"

class Face : public CGAL::Compact_container_base
{
	public:
	Face(){}
	Face(Vertex_iterator v0, Vertex_iterator v1, Vertex_iterator v2);
	Face(Face_iterator f0, Face_iterator f1, Face_iterator f2);
	Face(Vertex_iterator v0, Vertex_iterator v1, Vertex_iterator v2,
	     Face_iterator f0, Face_iterator f1, Face_iterator f2);
	Face(const Face& f);

	Face& operator=(const Face& f);
	bool operator==(const Face& f);
	bool operator!=(const Face& f);

	void vertex(Vertex_iterator v, int i);
	void neighbor(Face_iterator neighbor, int i);
	Vertex_iterator vertex(int i);
	Face_iterator neighbor(int i);
	int index(Vertex_iterator v);
	int index(Face_iterator fc);
	bool contains(Vertex_iterator v);
	bool is_vertex(Vertex_iterator v);
	bool is_neighbor(Face_iterator fc);
	void reorient();

	//void * for_compact_container() const {return _p;}
	//void *& for_compact_container() {return _p;}
	
	private:
	Vertex_iterator _vertices[3];
	Face_iterator _neighbor[3];
	//void* _p;
};

#endif
