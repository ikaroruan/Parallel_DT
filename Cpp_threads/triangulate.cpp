#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "Delaunay_triangulation.h"

int main(int argc, char** argv)
{

	if(argc < 2){
		std::cerr << "Sorry. No argument has been passed.\n";
		return 0;
	}
	else{
		std::string filename = argv[1];
		std::ifstream input(filename);
		Delaunay_triangulation t;

		int size, aux1, aux2, aux3;
		double x, y;
		std::vector<Point> vector;
		
		input >> size >> aux1 >> aux2 >> aux3;
		for(int i = 0; i < size; ++i){
			input >> aux1 >> x >> y;
			Point p(x, y);
			vector.push_back(p);
		}

		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point end;
		std::chrono::duration<double> elapsed;

		std::cout << "Starting parallel insertion.\n";
		start = std::chrono::high_resolution_clock::now();
		t.parallel_insert(vector);
		end = std::chrono::high_resolution_clock::now();
		elapsed = end - start;
		std::cout << "Insertion done.\n";
		std::cout << "Number of vertices: " << t.number_of_vertices() << "\n";
		std::cout << "Number of faces: " << t.number_of_faces() << "\n";
		std::cout << "Elapsed time = " << elapsed.count() << "\n";

		t.show_triangulation(filename);
		t.output_triangulation(filename);
	}

	return 0;
}
