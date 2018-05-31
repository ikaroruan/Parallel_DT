#include "Face_circulator.h"

Face_circulator::Face_circulator()
{
	_current_face = nullptr;
	_current_face = nullptr;
}

Face_circulator::Face_circulator(Vertex_iterator v)
{
	_current_vertex = v;
	_current_face = v->incident_face();
}

Face_circulator::Face_circulator(Vertex_iterator v, Face& fc)
{
	Face_iterator pfc = &fc;
	_current_vertex = v;
	if(pfc->contains(v))
		_current_face = pfc;
	else
		std::cerr << "ERROR: face does not contain the vertex. (Face_circulator)\n";
}

Face_circulator::Face_circulator(Vertex_iterator v, Face_iterator fc)
{
	_current_vertex = v;
	if(fc->contains(v))
		_current_face = fc;
	else
		std::cerr << "ERROR: face does not contain the vertex. (Face_circulator)\n";
}

Face_iterator Face_circulator::current_face()
{
	return _current_face;
}

void Face_circulator::current_face(Face_iterator fc)
{
	_current_face = fc;
}

Vertex_iterator Face_circulator::current_vertex()
{
	return _current_vertex;
}

void Face_circulator::current_vertex(Vertex_iterator v)
{
	_current_vertex = v;
}

int Face_circulator::ccw(int i)
{
	if(i < 0 || i > 2){
		std::cerr << "ERROR: wrong index for ccw function.\n";
		return -1; // reporting an error
	}
	if(i == 2)
		return 0;
	return i+1;
}

int Face_circulator::cw(int i)
{
	if(i < 0 || i > 2){
		std::cerr << "ERROR: wrong index for cw function.\n";
		return -1;
	}
	if(i == 0)
		return 2;
	return i-1;
}

Face_iterator Face_circulator::next()
{
	if(current_vertex() == nullptr){
		std::cerr << "ERROR: current vertex not initialized.\n";
		return nullptr;
	}

	if(current_face() == nullptr){
		std::cerr << "ERROR: current face not set or not initialized (Face_circulator). \n";
		return nullptr;
	}

	int pos = current_face()->index(current_vertex());
	Face_iterator next = current_face()->neighbor(ccw(pos));
	current_face(next);
	return current_face();
}

Face_iterator Face_circulator::previous()
{
	if(current_vertex() == nullptr){
		std::cerr << "ERROR: current vertex not initialized.\n";
		return nullptr;
	}

	if(current_face() == nullptr){
		std::cerr << "ERROR: current face not set or not initialized (Face_circulator). \n";
		return nullptr;
	}

	int pos = current_face()->index(current_vertex());
	Face_iterator previous = current_face()->neighbor(cw(pos));
	current_face(previous);
	return current_face();
}
