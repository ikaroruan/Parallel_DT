#include "Face.h"

Face::Face(Vertex_iterator v0, Vertex_iterator v1, Vertex_iterator v2)
{
	_vertices[0] = v0; _vertices[1] = v1; _vertices[2] = v2;
}

Face::Face(Face_iterator f0, Face_iterator f1, Face_iterator f2)
{
	_neighbor[0] = f0; _neighbor[1] = f1; _neighbor[2] = f2;
}

Face::Face(Vertex_iterator v0, Vertex_iterator v1, Vertex_iterator v2,
	   Face_iterator f0, Face_iterator f1, Face_iterator f2)
{
	_vertices[0] = v0; _vertices[1] = v1; _vertices[2] = v2;
	_neighbor[0] = f0; _neighbor[1] = f1; _neighbor[2] = f2;
}

Face::Face(const Face& f)
{
	_vertices[0] = f._vertices[0]; _neighbor[0] = f._neighbor[0];
	_vertices[1] = f._vertices[1]; _neighbor[1] = f._neighbor[1];
	_vertices[2] = f._vertices[2]; _neighbor[2] = f._neighbor[2];
}

Face& Face::operator=(const Face& f)
{
	_vertices[0] = f._vertices[0]; _neighbor[0] = f._neighbor[0];
	_vertices[1] = f._vertices[1]; _neighbor[1] = f._neighbor[1];
	_vertices[2] = f._vertices[2]; _neighbor[2] = f._neighbor[2];

	return *this;
}

bool Face::operator==(const Face& f)
{
	bool vertices = (_vertices[0] == f._vertices[0] && 
			 _vertices[1] == f._vertices[1] &&
			 _vertices[2] == f._vertices[2] );

	bool neighbor = (_neighbor[0] == f._neighbor[0] &&
			 _neighbor[1] == f._neighbor[1] &&
			 _neighbor[2] == f._neighbor[2] );
	
	return (vertices && neighbor);
}

bool Face::operator!=(const Face& f)
{
	return !(*this == f);
}

void Face::vertex(Vertex_iterator v, int i)
{
	_vertices[i] = v;
}

void Face::neighbor(Face_iterator neighbor, int i)
{
	_neighbor[i] = neighbor;
}

Vertex_iterator Face::vertex(int i)
{
	return _vertices[i];
}

Face_iterator Face::neighbor(int i)
{
	return _neighbor[i];
}

int Face::index(Face_iterator fc)
{
	if(fc == neighbor(0))
		return 0;
	if(fc == neighbor(1))
		return 1;
	if(fc == neighbor(2))
		return 2;
	// If fc is not a neighbor.
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
	// Face does not contain vertex.
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
