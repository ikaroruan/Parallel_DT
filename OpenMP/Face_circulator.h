#ifndef FACE_CIRCULATOR_H
#define FACE_CIRCULATOR_H

#include <iostream>
#include "Vertex.h"
#include "Face.h"

class Face_circulator
{
	public:

	Face_circulator();
	Face_circulator(Vertex_iterator v);
	Face_circulator(Vertex_iterator v, Face& fc);
	Face_circulator(Vertex_iterator v, Face_iterator fc);

	Face_iterator next();
	Face_iterator previous();
	Face_iterator current_face();
	Vertex_iterator current_vertex();
	
	private:

	Face_iterator _current_face;
	Vertex_iterator _current_vertex;
	void current_vertex(Vertex_iterator v);
	void current_face(Face_iterator fc);
	int ccw(int i);
	int cw(int i);
};

#endif
