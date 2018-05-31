#include "Face.h"

Face::Face()
{
	for(int i = 0; i < 3; ++i){
		_neighbor[i] = nullptr;
		_vertices[i] = nullptr;
	}
}

Face::Face(Vertex_iterator v0, Vertex_iterator v1, Vertex_iterator v2)
{
	_vertices[0] = v0;
	_vertices[1] = v1;
	_vertices[2] = v2;

	for(int i = 0; i < 3; ++i){
		_neighbor[i] = nullptr;
	}
}

Face::Face(Face_iterator f)
{
	_vertices[0] = f->vertex(0);
	_vertices[1] = f->vertex(1);
	_vertices[2] = f->vertex(2);

	_neighbor[0] = f->neighbor(0);
	_neighbor[1] = f->neighbor(1);
	_neighbor[2] = f->neighbor(2);
}

Face::~Face()
{

}

Face_iterator Face::neighbor(int i)
{
	return _neighbor[i];
}

Vertex_iterator Face::vertex(int i)
{
	return _vertices[i];
}

void Face::neighbor(Face_iterator neighbor, int i)
{
	_neighbor[i] = neighbor;
}

void Face::vertex(Vertex_iterator v, int i)
{
	_vertices[i] = v;
}

int Face::index(Face_iterator fc)
{
	if(fc == neighbor(0))
		return 0;
	if(fc == neighbor(1))
		return 1;
	if(fc == neighbor(2))
		return 2;
	// if fc is not a neighbor
	std::cerr << "ERROR: face does not have this neighbor. (index function)\n";
	return -1;
}

int Face::index(Vertex_iterator vi)
{
	if(vi == vertex(0))
		return 0;
	if(vi == vertex(1))
		return 1;
	if(vi == vertex(2))
		return 2;
	// face does not contain vertex
	std::cerr << "ERROR: face does not have the vertex. (index function)\n";
	return -1;
}

bool Face::contains(Vertex_iterator v)
{
	if(v == vertex(0) || v == vertex(1) || v == vertex(2))
		return true;
	return false;
}

bool Face::is_neighbor(Face_iterator fc)
{
	if(fc == neighbor(0) || fc == neighbor(1) || fc == neighbor(2))
		return true;
	return false;
}

bool Face::is_vertex(Vertex_iterator v)
{
	if(v == vertex(0) || v == vertex(1) || v == vertex(2))
		return true;
	return false;
}

void Face::reorient()
{
	Vertex_iterator v0 = vertex(0); Vertex_iterator v1 = vertex(1);
	Face_iterator n0 = neighbor(0); Face_iterator n1 = neighbor(1);

	vertex(v0, 1); vertex(v1, 0);
	neighbor(n0, 1); neighbor(n1, 0);
}

bool Face::lock()
{
	if(vertex(0)->lock() &&
	   vertex(1)->lock() &&
	   vertex(2)->lock() )
	   	return true;

	vertex(0)->unlock();
	vertex(1)->unlock();
	vertex(2)->unlock();

	return false;
}

bool Face::try_lock()
{
	if(vertex(0)->try_lock() &&
	   vertex(1)->try_lock() &&
	   vertex(2)->try_lock() )
	   	return true;

	vertex(0)->unlock();
	vertex(1)->unlock();
	vertex(2)->unlock();

	return false;
}

void Face::unlock()
{
	vertex(0)->unlock();
	vertex(1)->unlock();
	vertex(2)->unlock();
}

bool Face::any_vertex_locked()
{
	if(vertex(0)->is_locked() ||
	   vertex(1)->is_locked() ||
	   vertex(2)->is_locked() )
	   	return true;
	
	return false;
}

