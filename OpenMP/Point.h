#ifndef POINT_H
#define POINT_H

class Point
{
	public:
	
	Point() : _x(0), _y(0) {}
	Point(const Point& p) : _x(p._x), _y(p._y) {}
	Point(const double x, const double y) : _x(x), _y(y){}
	Point& operator=(const Point& p);
	bool operator==(const Point& p);
	bool operator!=(const Point& p);

	void set_x(double x);
	void set_y(double y);
	double get_x() const;
	double get_y() const;

	private:
	double _x;
	double _y;
};

#endif
