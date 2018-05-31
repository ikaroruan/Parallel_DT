#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "Delaunay_triangulation.h"

void print_location(Vertex_location location)
{
	if(location == ON_VERTEX)
		std::cout << "ON_VERTEX";
	else if(location == ON_EDGE)
		std::cout << "ON_EDGE";
	else if(location == ON_INTERIOR)
		std::cout << "ON_INTERIOR";
	else
		std::cout << "OUTSIDE_CONVEX_HULL";
}

int main(int argc, char** argv)
{
	/*
	Delaunay_triangulation t;
	std::ifstream input("box_out.tri");
	t.create_triangulation(input);
	
	std::vector<Point> points;
	Point p0(0.75, 0.5); points.push_back(p0);
	Point p1(1.5, 0.5); points.push_back(p1);
	Point p2(4, 4); points.push_back(p2);
	Point p3(1.5, 1.5); points.push_back(p3);
	int nthreads = 4;
	Vertex_location lc[nthreads];
	std::vector<Face_iterator> region[nthreads];
	std::vector<std::pair<Face_iterator, int>> cavity[nthreads];
	Vertex_iterator start_vertex = nullptr;// ++(t.vertices().begin());

	std::cout << "Begin of parallel region.\n";
	#pragma omp parallel private(start_vertex) shared(lc, points, region, cavity) num_threads(nthreads)
	{
		std::queue<Point> queue;
		int tid = omp_get_thread_num();
		int li; int* priority = nullptr;
		Face_iterator fc = t.parallel_locate(points[tid], start_vertex, lc[tid], li, priority);
		if(fc == nullptr)
			queue.push(points[tid]);
		else{
			int num = t.parallel_conflict_region(fc, points[tid], region[tid], cavity[tid], lc[tid], li, priority);
			if(num == 0)
				queue.push(points[tid]);
		}

		while(!queue.empty()){
			Point p = queue.front();
			queue.pop();
			fc = t.parallel_locate(p, start_vertex, lc[tid], li, priority);
			if(fc == nullptr){
				queue.push(points[tid]);
			}
			else{
				int num = t.parallel_conflict_region(fc, p, region[tid], cavity[tid], lc[tid], li, priority);
				if(num == 0)
					queue.push(points[tid]);
			}
		}
	} // END OF PARALLEL REGION
	std::cout << "End of parallel region.\n";

	std::cout << "==============  RESULTS  ================\n";
	std::cout << "Conflic region:\n";
	for(int tid = 0; tid < nthreads; ++tid){
		for(unsigned int i = 0; i < region[tid].size(); ++i){
			Face_iterator ff = (region[tid]).at(i);
			std::cout << "ff" << i << ":\n";
			if(ff->vertex(0) == t.infinite_vertex()) std::cout << "v0: Infinite.\n";
			else std::cout << "v0: x = " << ff->vertex(0)->point().get_x() << "   y = " << ff->vertex(0)->point().get_y() << "\n";
			if(ff->vertex(1) == t.infinite_vertex()) std::cout << "v1: Infinite.\n";
			else std::cout << "v1: x = " << ff->vertex(1)->point().get_x() << "   y = " << ff->vertex(1)->point().get_y() << "\n";
			if(ff->vertex(2) == t.infinite_vertex()) std::cout << "v2: Infinite.\n\n";
			else std::cout << "v2: x = " << ff->vertex(2)->point().get_x() << "   y = " << ff->vertex(2)->point().get_y() << "\n\n";
		}
		std::cout << "-------------------------------\n\n";
	}

	std::cout << "Location:\n";
	for(int i = 0; i < nthreads; ++i){
		std::cout << "Point " << i << ": location = ";
		print_location(lc[i]);
		std::cout << std::endl;
	}

	t.show_triangulation();

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
		//Vertex_iterator it = t.insert(p);
		//to_delete.push_back(it);
	}
	start = std::chrono::high_resolution_clock::now();
	t.parallel_insert(points);
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

	/*std::cout << "Deleting vertices...\n";
	start = std::chrono::high_resolution_clock::now();
	for(unsigned int i = 0; i < to_delete.size(); ++i){
		t.remove(to_delete[i]);
	}
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "All vertices deleted.\n";
	std::cout << "Elapsed time = " << elapsed.count() << "s.\n";
	

	t.show_triangulation(filename);
	t.output_triangulation(filename);
	return 0;
	*/

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
