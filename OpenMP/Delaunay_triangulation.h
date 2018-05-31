#ifndef DELAUNAY_TRIANGULATION_H
#define DELAUNAY_TRIANGULATION_H

#include <iostream>
#include <omp.h>
#include <chrono>
#include "Triangulation.h"

#define NUMBER_OF_THREADS 4

class Delaunay_triangulation : public Triangulation
{
	public:
	
	// INSERTION
	Vertex_iterator insert(Point& p);
	Vertex_iterator insert_first(Point& p);
	Vertex_iterator insert_second(Point& p);
	Vertex_iterator insert_in_face(Face_iterator fc, Point& p, Vertex_location lc, int li);
	Vertex_iterator insert_in_hole(Face_iterator fc, Point& p, int li, Vertex_location lc);
	Vertex_iterator insert_on_edge_2(Face_iterator fc, Point& p, Vertex_location lc, int li);
	Vertex_iterator insert_dimension_2(Point& p);
	Vertex_iterator insert_in_single_face(Face_iterator fc, Point& p);
	Vertex_iterator insert_outside_convex_hull_2(Face_iterator fc, Point& p, Vertex_location lc, int li);
	Vertex_iterator insert_dimension_1(Point& p);
	Vertex_iterator insert_on_edge_1(Face_iterator fc, Point& p, Vertex_location lt, int li);
	Vertex_iterator insert_outside_convex_hull_1(Point& p, Face_iterator fc, int li);
	void 	insert(std::vector<Point>::iterator start, std::vector<Point>::iterator end);

	// REMOVAL
	//void remove(Point& p);
	void remove(Vertex_iterator v);
	void remove_dim_0(Vertex_iterator v);
	void remove_dim_1(Vertex_iterator v);
	void remove_dim_2(Vertex_iterator v);
	void remove_dim_decrease(Vertex_iterator v);

	// GEOMETRIC PREDICATES
	double incircle_test(Vertex_iterator va, Vertex_iterator vb, Vertex_iterator vc, Vertex_iterator vd);
	double incircle_test(Point& a, Point& b, Point& c, Point& d);

	// UTILITIES
	void conflict_region(Face_iterator fc, Point& pp, std::vector<Face_iterator>& region, 
		std::vector<std::pair<Face_iterator, int>>& cavity, Vertex_location lc, int li);
	void flip_check(std::queue<std::pair<Face_iterator, int>>& queue);

	// PARALLEL
	void parallel_insert(std::vector<Point>& p);
	Face_iterator parallel_locate(Point& p, Vertex_iterator start_vertex, Vertex_location& location, int& li, int* priority);
	int parallel_conflict_region(Face_iterator fc, Point& pp, std::vector<Face_iterator>& region, 
		std::vector<std::pair<Face_iterator, int>>& cavity, Vertex_location lc, int li, int* priority);
	bool parallel_try_insert(Face_iterator fc, Point& p, Vertex_location lc, int li, 
		std::vector<Face_iterator> region, std::vector<std::pair<Face_iterator, int>> cavity);
	void parallel_insert_in_hole(Face_iterator fc, Point& p, Vertex_location lc, int li, 
		std::vector<Face_iterator> region, std::vector<std::pair<Face_iterator, int>> cavity);

#ifndef NDEBUG
	std::chrono::duration<double> elapsed_locate;
	std::chrono::duration<double> elapsed_conflict;
	std::chrono::duration<double> elapsed_update;
	std::chrono::duration<double> elapsed_parallel;

	int locate_retreats = 0;
	int conflict_retreats = 0;
#endif

};

#endif // DELAUNAY_TRIANGULATION_H
