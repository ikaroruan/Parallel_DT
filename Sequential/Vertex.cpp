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
