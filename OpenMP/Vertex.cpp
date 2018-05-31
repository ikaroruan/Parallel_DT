#include "Vertex.h"

void Vertex::point(Point& p)
{
	_point = p;
}

void Vertex::incident_face(Face_iterator f)
{
	_incident_face = f;
}

Point& Vertex::point()
{
	return _point;
}

Face_iterator Vertex::incident_face()
{
	return _incident_face;
}

Vertex& Vertex::operator=(const Vertex& v)
{
	_point = v._point;
	_incident_face = v._incident_face;
	return *this;
}

bool Vertex::operator==(const Vertex& v)
{
	return (_point == v._point && _incident_face == v._incident_face);
}

bool Vertex::operator!=(const Vertex& v)
{
	return !(*this == v);
}

bool Vertex::lock()
{
	return _lk.lock();
}

bool Vertex::lock(int* priority)
{
	return _lk.lock(priority);
}

bool Vertex::unlock()
{
	return _lk.unlock();
}

bool Vertex::is_locked()
{
	return _lk.locked();
}

bool Vertex::is_locked_by_owner(int i)
{
	return _lk.locked_by_owner(i);
}

int Vertex::lock_owner()
{
	return _lk.owner();
}
