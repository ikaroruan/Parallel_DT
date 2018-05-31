#include "Point.h"

Point& Point::operator=(const Point& p)
{
	_x = p._x;
	_y = p._y;
	return *this;
}

bool Point::operator==(const Point& p)
{
	return (_x == p._x && _y == p._y);
}

bool Point::operator!=(const Point& p)
{
	return !(*this == p);
}

void Point::set_x(double x)
{
	_x = x;
}

void Point::set_y(double y)
{
	_y = y;
}

double Point::get_x() const
{
	return _x;
}

double Point::get_y() const
{
	return _y;
}
