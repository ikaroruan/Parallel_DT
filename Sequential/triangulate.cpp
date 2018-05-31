#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "Delaunay_triangulation.h"

int main(int argc, char** argv)
{
	/*
	Delaunay_triangulation t;
	Point p0(0.25, 0.25);
	Point p1(1, 0);
	Point p2(2, 1);
	Point p3(1.75, 1.75);
	Point p4(1, 2);
	Point p5(0, 1);
	Point p6(1, 1);
	Point p7(3, 3);
	Point p8(2.05, 0);

	t.insert(p0);
	t.insert(p1);
	t.insert(p2);
	t.insert(p3);
	t.insert(p4);
	t.insert(p5);
	t.insert(p6);
	Vertex_iterator to_delete = t.insert(p7);
	t.insert(p8);

	t.remove(to_delete);

	t.show_triangulation();
	t.output_triangulation();

	return 0;
	*/
	
	if(argc < 2){
		std::cout << "Not enough arguments.\n";
		return 0;
	}
	std::string filename = argv[1];
	std::ifstream input(filename);
	Delaunay_triangulation t;
	int size; int aux1; int aux2; int aux3;
	double x; double y;

	input >> size >> aux1 >> aux2 >> aux3;

	std::cout << "Building triangulation...\n";
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;
	std::chrono::duration<double> elapsed;
	std::vector<Vertex_iterator> to_delete;
	std::vector<Point> points;

	for(int i = 0; i < size; ++i){
		input >> aux1 >> x >> y;
		Point p(x, y);
		points.push_back(p);
	}
	start = std::chrono::high_resolution_clock::now();
	t.insert(points.begin(), points.end());
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Triangulation done.\n";
	std::cout << "Number of vertices: " << t.number_of_vertices() << "\n";
	std::cout << "Number of Faces: " << t.number_of_faces() << "\n";
	std::cout << "Elapsed time = " << elapsed.count() << "s.\n\n";
	
	//Point p(2.5, 2.5);
	//Vertex_iterator del = t.insert(p);
	//t.remove(del);

	t.show_triangulation(filename);
	t.output_triangulation(filename);
	return 0;

	/*
	Delaunay_triangulation t;

	Point p(0, 1);
	Vertex_iterator pit1 = t.insert(p);
	Point p2(1, 1);
	Vertex_iterator pit2 = t.insert(p2);
	Point p3(0.25, 1);
	Vertex_iterator pit3 = t.insert(p3);
	Point p4(0.75, 1);
	Vertex_iterator pit4 = t.insert(p4);

	Point p5(0.5, 2);
	Vertex_iterator pit5 = t.insert(p5);
	Point p6(0.5, 0);
	Vertex_iterator pit6 = t.insert(p6);
	//t.remove(pit2);
	//t.remove(pit4);
	//t.remove(pit1);
	//t.remove(pit3);
	//t.remove(pit5);
	t.remove(pit5);
	t.remove(pit6);

	Point p7(1, 0);
	Vertex_iterator pit7 = t.insert(p7);
	Point p8(0.5, 0.75);
	Vertex_iterator pit8 = t.insert(p8);

	t.remove(pit7);
	t.remove(pit8);
	t.remove(pit1);
	Point p9(1, 2);
	t.insert(p9);

	t.remove(pit3);
	t.remove(pit4);
	Point p10(0, 1.5);
	t.insert(p10);

	t.output_triangulation();
	t.show_triangulation();

	return 0;
	*/
}
