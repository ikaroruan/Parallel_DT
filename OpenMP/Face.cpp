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

bool Face::lock()
{
	bool done = false;
	#pragma omp critical (Face_lock)
	{
		if(vertex(0)->lock() &&
		   vertex(1)->lock() &&
		   vertex(2)->lock() )
			done = true;

		if(!done){
			vertex(0)->unlock();
			vertex(1)->unlock();
			vertex(2)->unlock();
		}
	}

	return done;
}

bool Face::lock(int* priority)
{
	bool done = false;
	#pragma omp critical (Face_lock_priority)
	{
		if(vertex(0)->lock(priority) &&
		   vertex(1)->lock(priority) &&
		   vertex(2)->lock(priority) )
			done = true;
	
		if(!done){
			vertex(0)->unlock();
			vertex(1)->unlock();
			vertex(2)->unlock();
		}
	}

	return done;
}

void Face::unlock()
{
	#pragma omp critical (Face_unlock)
	{
		vertex(0)->unlock();
		vertex(1)->unlock();
		vertex(2)->unlock();
	}
}

bool Face::is_locked()
{
	if(vertex(0)->is_locked() && 
	   vertex(1)->is_locked() &&
	   vertex(2)->is_locked() )
		return true;
	
	// Else
	return false;
}

bool Face::any_vertex_locked()
{
	if(vertex(0)->is_locked() || 
	   vertex(1)->is_locked() ||
	   vertex(2)->is_locked() )
		return true;
	
	// Else
	return false;
}

bool Face::is_locked_by_owner(int i)
{
	if(vertex(0)->is_locked_by_owner(i) &&
	   vertex(1)->is_locked_by_owner(i) &&
	   vertex(2)->is_locked_by_owner(i) )
	   	return true;

	return false;
}

/*bool Face::double_locked()
{
	// As -1 is owner number when not locked.
	int o0 = -2; int o1 = -3; int o2 = -4;
	#pragma omp critical (locked)
	{
		if(vertex(0)->is_locked())
			o0 = vertex(0)->lock_owner();
		if(vertex(1)->is_locked())
			o1 = vertex(1)->lock_owner();
		if(vertex(2)->is_locked())
			o2 = vertex(2)->lock_owner();
	}

	if(o1 == o2 || o2 == o3 || o3 == o1)
		return true;

	return false;
}

bool Face::double_lock()
{
	if(double_locked())
		return false;
	#pragma omp critical (lock)
	{
		int count = 0;
		int index = 0;
		while(count != 2 || index != 2){
			if(vertex(index) != infinite_vertex() && !vertex(index)->is_locked()){
				vertex(index)->lock();
				count++;
			}
			index++;
		}
	}

	if(count != 2){
		for(int i = 0; i < 2; ++i)
			vertex(i)->unlock();
		return false;
	}
	return true;
}*/
