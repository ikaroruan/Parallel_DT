#include "Vertex.h"
#include <iostream>

Vertex::Vertex()
{
	_in_face = nullptr;
}

Vertex::Vertex(Point& p)
{
	_point = p;
	_in_face = nullptr;
}

void Vertex::point(Point& p)
{
	_point = p;
}

Vertex& Vertex::operator=(const Vertex& v)
{
	_point = v._point;
	_in_face = v._in_face;
	return *this;
}

bool Vertex::operator==(const Vertex& v)
{
	return (_point == v._point && _in_face == v._in_face);
}

bool  Vertex::operator!=(const Vertex& v)
{
	return !(*this == v);
}

Point& Vertex::point()
{
	return _point;
}

Face_iterator Vertex::incident_face()
{
	return _in_face;
}

void Vertex::incident_face(Face_iterator fc)
{
	_in_face = fc;
}	

bool Vertex::lock()
{
	return _lck.lock();
}

bool Vertex::try_lock()
{
	return _lck.try_lock();
}

bool Vertex::unlock()
{
	return _lck.unlock();
}

bool Vertex::is_locked()
{
	return _lck.locked();
}
