/*
 * Utility.h
 *
 *  Created on: 2018
 *      Author: mattia
 */

#ifndef UTILITY_H_
#define UTILITY_H_

using namespace std;
#include <omp.h>
#include <xmmintrin.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <cmath>
#include "mytimer.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "mersenneTwister.h"
#include <unordered_set>
#include <CGAL/Line_2.h>
//#include <CGAL/Vector_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Simple_cartesian<double> K1;
//typedef CGAL::Cartesian_converter<K,K1> K_to_K1;
//typedef CGAL::Cartesian_converter<K1,K> K1_to_K;

typedef K::Point_2 Point;
//typedef K::Vector_2 Vector_2;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Line_2<K> Line_2;

typedef CGAL::Segment_2<K> Segment_2;
//typedef CGAL::Triangle_2<K> Triangle_2;
typedef K::Intersect_2 Intersect_2;
typedef std::vector<Point> PointSet;

typedef double areaType;
//

class Utility{

public:


	static void read_instance(std::string,PointSet&);
	static void read_snapshot(std::string,PointSet&);
	static int test_insertion_between(std::vector<Segment_2> &,Point&,Point&,Point&);

	static void write_solution(std::string,PointSet&,double);
	static void write_snapshot(std::string,PointSet&);

	static areaType compute_triangle_area(Point&, Point&, Point&);
	static areaType compute_polygon_area(std::vector<Point>&);

	static double string2double(std::string&);
	static int string2int(std::string&);

	static void genetic(PointSet&,std::string&,bool,int,int,int,int,int,std::string&,int);

	static Polygon_2 pointset2polygon(PointSet&);
	static PointSet polygon2pointset(Polygon_2&);

	static void add_in_position(PointSet&,int,Point);
	static void add_in_position(std::vector<Segment_2>&,int,Segment_2);

	static areaType default_value;
	static areaType bestcost;
	static int inner_block_type;

	static void init_population(int);
	static void evaluate(std::tuple<std::vector<int>,areaType,int>&);
	static void crossover(int&,int&,std::vector<int>&,std::vector<int>&,bool);
	static void mutation(int);

	static int selection(int,bool);
	static int MersenneShuffler(int);
	static int eras;
	static int selection_perc;
	static int mutation_prob;
	static bool minimization;
	static int phenotyper;
	static long  unsigned int wasted;
	static long  unsigned int backtracked;

	static PointSet original_pointset;
	static std::vector<int> original_permutation;
//	static std::vector<areaType> areas;

	static PointSet best_pointset;
	static std::vector<std::tuple<areaType,areaType,areaType>> evolution;
	static std::vector<std::tuple<std::vector<int>,areaType,int>> population;
	static std::vector<int> pointset2permutation(PointSet&);

	static PointSet permutation2pointset(std::vector<int>&);
	static Point index2point(int);
	static int point2index(Point);
    // Check if three points make a right turn using cross product
    static bool right_turn(const Point &p1, const Point &p2, const Point &p3)
    {
        return ((p3.x()-p1.x()) * (p2.y()-p1.y()) - (p3.y()-p1.y())*(p2.x()-p1.x())) > 0;
    }
	static void print_point_set(PointSet&);
	static MersenneTwister randomGenerator;
	static areaType genotype2phenotype(std::vector<int>&,int &);
	static areaType genotype2phenotype(std::vector<int>&,PointSet&ps,int &);
//	static areaType stddev(std::vector<areaType>&);
//	static areaType avg(std::vector<areaType>&);
	static bool checkIfFIle(std::string);

	static std::tuple<std::vector<int>,PointSet,areaType> split_method(bool);
	static std::tuple<std::vector<int>,PointSet,areaType> concentric_hulls(bool);
	static std::tuple<std::vector<int>,PointSet,areaType> polar();

	struct greaterPS{
		bool operator()(std::tuple<std::vector<int>,areaType,int> const &a, std::tuple<std::vector<int>,areaType,int> const &b) const {
			return get<1>(a) > get<1>(b);

		}
	};
	struct smallerPS{

		bool operator()(
				std::tuple<std::vector<int>,areaType,int> const &a,
				std::tuple<std::vector<int>,areaType,int> const &b) const {

			return get<1>(a) < get<1>(b);
		}
	};
	struct XY_Sorter{

		bool operator()(
				Point const &a,
				Point const &b) const {

			return CGAL::to_double(a.x()) < CGAL::to_double(b.x())
			|| (CGAL::to_double(a.x()) == CGAL::to_double(b.x()) && CGAL::to_double(a.y()) < CGAL::to_double(b.y()));
		}
	};
	struct SIN_SORT{

			bool operator()(
					std::tuple<double,double,Point> const &a,
					std::tuple<double,double,Point> const &b) const {

				return std::get<0>(a)<std::get<0>(b);
			}
		};

	struct COS_SORT{

			bool operator()(
					std::tuple<double,double,Point> const &a,
					std::tuple<double,double,Point> const &b) const {
				return std::get<1>(a)<std::get<1>(b);
			}
		};

	struct YX_Sorter{

			bool operator()(
					Point const &a,
					Point const &b) const {

				return CGAL::to_double(a.y()) < CGAL::to_double(b.y())
				|| (CGAL::to_double(a.y()) == CGAL::to_double(b.y()) && CGAL::to_double(a.x()) < CGAL::to_double(b.x()));
			}
		};



};


#endif

