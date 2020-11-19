/*
* Utility.cpp
*
*  Created on: 2017
*      Author: mattia
*/


using namespace std;

#include "Utility.h"
#include "progressBar.h"
#include <CGAL/convex_hull_2.h>
areaType Utility::default_value = 0;
areaType Utility::bestcost = 0;

PointSet Utility::original_pointset;
std::vector<int> Utility::original_permutation;
PointSet Utility::best_pointset;
std::vector<std::tuple<std::vector<int>,areaType,int>> Utility::population;
std::vector<std::tuple<areaType,areaType,areaType>> Utility::evolution;
//std::vector<areaType> Utility::areas;

bool Utility::minimization = false;
int Utility::mutation_prob = 0;
int Utility::eras = 0;

int Utility::selection_perc = 0;
int Utility::phenotyper = 0;
long unsigned int Utility::wasted = 0;
long unsigned int Utility::backtracked = 0;

int Utility::inner_block_type = 0; // DEFAULT CGAL SIMPLE GENERATOR
MersenneTwister Utility::randomGenerator;

namespace {

// Return the number of threads that would be executed in parallel regions
int get_max_threads() {
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}

int get_current_procs(){
#ifdef _OPENMP
	return omp_get_num_procs();
#else
	return 1;
#endif
}

// Set the number of threads that would be executed in parallel regions
void set_num_threads(int num_threads) {
#ifdef _OPENMP
	omp_set_num_threads(num_threads);
#else
	if (num_threads != 1) {
		throw std::runtime_error("compile with -fopenmp");
	}
#endif
}

// Return my thread ID
int get_thread_id() {
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

template<typename T>
struct parallel_vector {
	parallel_vector(size_t size_limit)
	: v(get_max_threads(), std::vector<T>(size_limit)),
	  n(get_max_threads(), 0) {}

	void push_back(const T &x) {
		int id = get_thread_id();
		v[id][n[id]++] = x;
	}

	void clear() {
		for(int i = 0;i<get_max_threads();i++)
			n[i] = 0;
	}

	std::vector<std::vector<T> > v;
	std::vector<size_t> n;
};
}  // namespace

void Utility::read_instance(std::string source,PointSet& poly){
	ifstream ifs;
	ifs.open(source, std::ifstream::in);

	if(!ifs.is_open())
	    throw runtime_error("Error opening file");


	std::string line;



	while(std::getline(ifs,line)){


		if(line.find("#")!=std::string::npos)
			continue;

		std::vector<string> result;
		boost::split(result, line, boost::is_any_of("\t "),boost::token_compress_on);
		poly.push_back(Point(string2double(result[1]),string2double(result[2])));
	}
	std::cout<<"POINTSET_SIZE: "<<poly.size()<<"\n";

	ifs.close();
}

void Utility::write_snapshot(std::string outfile,PointSet& p){


	ofstream ofs;
	std::string namefile = Utility::minimization==true ? outfile+".minresume" : outfile + ".maxresume";

	ofs.open(namefile, std::ifstream::out);

	if(!ofs.is_open())
	    throw runtime_error("Error opening file");


	ofs<<p.size()<<" "<<Utility::population.size()<<" "<<Utility::eras<<" "<<Utility::selection_perc<<" "<<Utility::mutation_prob<<" "<<Utility::minimization<<"\n";//" "<<Utility::discarded_pointsets.size()<<

	for(size_t t=0;t<Utility::population.size();t++){
		for(size_t j=0;j<get<0>(Utility::population[t]).size()-1;j++)
			ofs<<get<0>(Utility::population[t])[j]<<" ";
		ofs<<get<0>(Utility::population[t])[get<0>(Utility::population[t]).size()-1]<<" ";
		ofs<<get<2>(Utility::population[t])<<"\n";

	}


	ofs.close();

}


void Utility::write_solution(std::string outfile,PointSet& input,double elaps){

	if(!pointset2polygon(Utility::best_pointset).is_simple()){
		throw std::runtime_error("NON SIMPLE SOLUTION");
	}

	{
		PointSet test = input;
		std::sort(test.begin(),test.end(),XY_Sorter());
		auto it = std::unique(test.begin(),test.end());
		bool wasUnique = (it == test.end() );
		if(!wasUnique)
			throw std::runtime_error("REPETITIONS IN SOLUTION");
	}

	ofstream ofs;
	std::string namefile = Utility::minimization==true ? outfile+".minsolution" : outfile + ".maxsolution";

	ofs.open(namefile, std::ofstream::out);
	assert(input.size()==Utility::best_pointset.size());

	if(!ofs.is_open())
	    throw runtime_error("Error opening file");

	ofs<<"# file: "<<outfile<<"\n";
	ofs<<"# area: "<<std::fixed<<std::setprecision(2)<<Utility::bestcost<<"\n";
	ofs<<"# index\n";

	for(size_t t=0;t<input.size();t++){
		PointSet::iterator it = std::find(input.begin(), input.end(), Utility::best_pointset[t]);
		int index = std::distance(input.begin(), it);
		ofs<<index<<"\n";
	}
	ofs.close();

	std::time_t rawtime;
	std::tm* timeinfo;
	char buffer [80];

	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);
	std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
	std::string timestring = buffer;

	namefile = Utility::minimization==true ? outfile + "."+timestring + ".mindetails" : outfile + "."+timestring + ".maxdetails";
	ofs.open(namefile, std::ofstream::out);

	assert(input.size()==Utility::best_pointset.size());

	if(!ofs.is_open())
		throw runtime_error("Error opening file");

	std::tuple<std::vector<int>,PointSet,areaType> split1 = Utility::split_method(true);
	std::tuple<std::vector<int>,PointSet,areaType> split2 = Utility::split_method(false);
	std::tuple<std::vector<int>,PointSet,areaType> conce1 = Utility::concentric_hulls(true);
	std::tuple<std::vector<int>,PointSet,areaType> conce2 = Utility::concentric_hulls(false);
	std::tuple<std::vector<int>,PointSet,areaType> polar_ = Utility::polar();
	PointSet convex_points;
	CGAL::convex_hull_2(Utility::original_pointset.begin(),Utility::original_pointset.end(), std::back_inserter(convex_points));
	Polygon_2 CH = pointset2polygon(convex_points);

	ofs<<"file;time;best_area;ch;split1;split2;concentr1;concentr2;polar;popul;eras;sel;mut;pheno;minmax;backtracked;wasted\n";
	ofs<<outfile<<";"<<std::fixed<<std::setprecision(2)
	   <<elaps<<";"
	   <<Utility::bestcost<<";"
	   <<fabs(CGAL::to_double(CH.area()))<<";"
	   <<get<2>(split1)<<";"
	   <<get<2>(split2)<<";"
	   <<get<2>(conce1)<<";"
	   <<get<2>(conce2)<<";"
	   <<get<2>(polar_)<<";"
	   <<Utility::population.size()<<";"
	   <<Utility::eras<<";"
	   <<Utility::selection_perc<<";"
	   <<Utility::mutation_prob<<";"
	   <<Utility::phenotyper<<";"
	   <<Utility::minimization<<";"
	   <<Utility::backtracked<<";"
	   <<Utility::wasted<<std::endl;

	for(size_t t=0;t<Utility::evolution.size();t++){
		if(std::get<0>(Utility::evolution[t])!=Utility::default_value)
			ofs<<std::fixed<<std::setprecision(2)<<std::get<0>(Utility::evolution[t])<<";"<<std::get<1>(Utility::evolution[t])<<";"<<std::get<2>(Utility::evolution[t])<<std::endl;
		else
			ofs<<std::fixed<<std::setprecision(2)<<-1<<";"<<std::get<1>(Utility::evolution[t])<<";"<<std::get<2>(Utility::evolution[t])<<std::endl;

	}
	ofs<<"[";
	for(size_t t=0;t<Utility::best_pointset.size();t++){
		ofs<<"("<<Utility::best_pointset[t].x()<<","<<Utility::best_pointset[t].y()<<")";
		if(t<Utility::best_pointset.size()-1)
			ofs<<",";
	}
	ofs<<"]\n";

	ofs.close();

}


std::tuple<std::vector<int>,PointSet,areaType> Utility::split_method(bool horizontal){




	Point endpoint_one = horizontal ? *CGAL::left_vertex_2(Utility::original_pointset.begin(),Utility::original_pointset.end()) :
			*CGAL::bottom_vertex_2(Utility::original_pointset.begin(),Utility::original_pointset.end());
	Point endpoint_two = horizontal ? *CGAL::right_vertex_2(Utility::original_pointset.begin(),Utility::original_pointset.end()) :
			*CGAL::top_vertex_2(Utility::original_pointset.begin(),Utility::original_pointset.end());


	Line_2 line = Line_2(endpoint_one,endpoint_two);

	PointSet sorted_version = Utility::original_pointset;

	if(horizontal)
		std::sort(sorted_version.begin(),sorted_version.end(),XY_Sorter());
	else
		std::sort(sorted_version.begin(),sorted_version.end(),YX_Sorter());

	assert(sorted_version[0]==endpoint_one);
	assert(sorted_version[sorted_version.size()-1]==endpoint_two);

	Polygon_2 output;
	output.push_back(endpoint_one);
	assert(line.has_on(endpoint_one));
	assert(line.has_on_boundary(endpoint_one));
	assert(line.has_on(endpoint_two));
	assert(line.has_on_boundary(endpoint_two));

	for(size_t t=1;t<sorted_version.size()-1;t++){
		Point p = sorted_version[t];
		assert(!line.has_on(p));
		assert(!line.has_on_boundary(p));
		if(line.has_on_negative_side(p))
			output.push_back(p);
	}
	output.push_back(endpoint_two);

	for(size_t t=sorted_version.size()-2;t>=1;t--){
		Point p = sorted_version[t];
		assert(!line.has_on(p));
		assert(!line.has_on_boundary(p));
		if(line.has_on_positive_side(p))
			output.push_back(p);
	}


	assert(output.is_simple());
	PointSet ps = polygon2pointset(output);
	std::tuple< std::vector<int>,PointSet,areaType > individual = std::make_tuple(
			pointset2permutation(ps),ps,
			fabs(CGAL::to_double(output.area()))
			);
	return individual;


}
std::tuple<std::vector<int>,PointSet,areaType> Utility::polar(){

	std::vector<std::tuple<double,double,Point>> polar_version_original;

	PointSet convex_points;



	CGAL::convex_hull_2(Utility::original_pointset.begin(),Utility::original_pointset.end(), std::back_inserter(convex_points));

	Point endpoint_one = *CGAL::left_vertex_2(convex_points.begin(),convex_points.end());
	Point endpoint_two = *CGAL::right_vertex_2(convex_points.begin(),convex_points.end());

	Point origin = CGAL::midpoint(endpoint_one,endpoint_two);

	PointSet & ref =  Utility::original_pointset;
	std::vector<std::tuple<double,double,Point>> & polarized_version =  polar_version_original;

	for(size_t t=0;t<ref.size();t++){
		Point considered = ref[t];
		Point projection_x = Point(ref[t].x(),origin.y());
		Point projection_y = Point(origin.x(),ref[t].y());

		Segment_2 s_x = Segment_2(origin,projection_x);
		Segment_2 s_y = Segment_2(origin,projection_y);

		Segment_2 L = Segment_2(origin,considered);

		double cosx  = sqrt(CGAL::to_double(s_x.squared_length()));
		double sinx  = sqrt(CGAL::to_double(s_y.squared_length()));


		if(CGAL::to_double(considered.x()) < CGAL::to_double(origin.x()))
			cosx*=-1;
		if(CGAL::to_double(considered.y()) < CGAL::to_double(origin.y()))
			sinx*=-1;
		cosx/=sqrt(CGAL::to_double(L.squared_length()));
		sinx/=sqrt(CGAL::to_double(L.squared_length()));

		assert(sinx<=1 && sinx>=-1);
		assert(cosx<=1 && cosx>=-1);
		polarized_version.push_back(std::make_tuple(std::atan2(sinx,cosx),sqrt(CGAL::to_double(L.squared_length())),ref[t]));

	}

	std::sort(polar_version_original.begin(),polar_version_original.end());


	Polygon_2 output;

	for(size_t t=0;t<polar_version_original.size();t++){
		output.push_back(std::get<2>(polar_version_original[t]));
	}

	assert(output.is_simple());
	PointSet ps = polygon2pointset(output);

	std::tuple< std::vector<int>,PointSet,areaType > individual = std::make_tuple(
					pointset2permutation(ps),ps,
					fabs(CGAL::to_double(output.area()))
			);
			return individual;


}
std::tuple<std::vector<int>,PointSet,areaType> Utility::concentric_hulls(bool horizontal){

	PointSet convex_points;

	CGAL::convex_hull_2(Utility::original_pointset.begin(),Utility::original_pointset.end(), std::back_inserter(convex_points));

	assert(convex_points.size()<=Utility::original_pointset.size());

//	std::cout << convex_points.size() << " points on the convex hull" << std::endl;
//	print_point_set(convex_points);
	Point endpoint_one = horizontal ? *CGAL::left_vertex_2(convex_points.begin(),convex_points.end()) :
			*CGAL::bottom_vertex_2(convex_points.begin(),convex_points.end());

	Point endpoint_two = horizontal ? *CGAL::right_vertex_2(convex_points.begin(),convex_points.end()) :
			*CGAL::top_vertex_2(convex_points.begin(),convex_points.end());

	Line_2 line = Line_2(endpoint_one,endpoint_two);

	PointSet sorted_convex_version = convex_points;
	PointSet sorted_version = Utility::original_pointset;

	if(horizontal){
		std::sort(sorted_convex_version.begin(),sorted_convex_version.end(),XY_Sorter());
		std::sort(sorted_version.begin(),sorted_version.end(),XY_Sorter());
	}
	else{
		std::sort(sorted_convex_version.begin(),sorted_convex_version.end(),YX_Sorter());
		std::sort(sorted_version.begin(),sorted_version.end(),YX_Sorter());

	}

	assert(sorted_convex_version[0]==endpoint_one);
	assert(sorted_convex_version[sorted_convex_version.size()-1]==endpoint_two);

	PointSet output_1,output_2;
	output_1.push_back(endpoint_one);
	output_2.push_back(endpoint_one);

	assert(line.has_on(endpoint_one));
	assert(line.has_on_boundary(endpoint_one));
	assert(line.has_on(endpoint_two));
	assert(line.has_on_boundary(endpoint_two));


	for(size_t t=1;t<sorted_convex_version.size()-1;t++){
		Point p = sorted_convex_version[t];
		assert(!line.has_on(p));
		assert(!line.has_on_boundary(p));

		if(line.has_on_negative_side(p))
			output_1.push_back(p);

		if(line.has_on_positive_side(p))
			output_2.push_back(p);
	}

	output_1.push_back(endpoint_two);
	output_2.push_back(endpoint_two);

	for(size_t t=sorted_version.size()-2;t>=1;t--){
		Point p = sorted_version[t];
		assert(!line.has_on(p));
		assert(!line.has_on_boundary(p));


		if(line.has_on_positive_side(p) || std::find(output_1.begin(),output_1.end(),p)==output_1.end())
			output_1.push_back(p);
		if(line.has_on_negative_side(p) || std::find(output_2.begin(),output_2.end(),p)==output_2.end())
			output_2.push_back(p);
	}

	assert(output_1.size()==output_2.size());
	assert(output_1.size()==Utility::original_pointset.size());

	std::tuple< std::vector<int>,PointSet,areaType > individual;
	Polygon_2 out_1=pointset2polygon(output_1),out_2=pointset2polygon(output_2);

	assert(out_1.is_simple());
	assert(out_2.is_simple());


	if(Utility::minimization){
		if(fabs(CGAL::to_double(out_1.area()))<fabs(CGAL::to_double(out_2.area()))){

			individual = std::make_tuple(
					pointset2permutation(output_1),output_1,
					fabs(CGAL::to_double(out_1.area()))
			);
			return individual;

		}
		else{
			individual = std::make_tuple(
					pointset2permutation(output_2),output_2,
					fabs(CGAL::to_double(out_2.area()))
			);
			return individual;
		}
	}
	else{
		if(fabs(CGAL::to_double(out_1.area()))>fabs(CGAL::to_double(out_2.area()))){
			individual = std::make_tuple(
					pointset2permutation(output_1),output_1,
					fabs(CGAL::to_double(out_1.area()))
			);
			return individual;

		}
		else{
			individual = std::make_tuple(
					pointset2permutation(output_2),output_2,
					fabs(CGAL::to_double(out_2.area()))
			);
			return individual;
		}
	}

}

void Utility::genetic(PointSet& poly,std::string& output, bool min_problem, int popul, int epochs, int mut_prob,
		int selection_percentage,int resume_from_file,std::string& resume_file,int phen){

	mytimer tms;
	Utility::minimization = min_problem;
	Utility::wasted = 0;
	Utility::backtracked = 0;
	Utility::phenotyper = phen;
	assert(phen>=0 && phen<=3);
	Utility::eras = epochs;
	Utility::default_value = Utility::minimization==true ? std::numeric_limits<areaType>::max() : 0.0;
	Utility::best_pointset.clear();
	Utility::original_pointset.clear();


	for(size_t t=0;t<poly.size();t++){
		Utility::original_pointset.push_back(poly[t]);
		Utility::original_permutation.push_back(t);
	}
	std::cout<<"[";
	for(size_t t=0;t<Utility::original_pointset.size();t++){
		std::cout<<"("<<Utility::original_pointset[t].x()<<","<<Utility::original_pointset[t].y()<<")";
		if(t<Utility::original_pointset.size()-1)
			std::cout<<",";
	}
	std::cout<<"]\n";
#ifndef NDEBUG


	{
		PointSet test = Utility::original_pointset;
		std::sort(test.begin(),test.end(),XY_Sorter());
		auto it = std::unique(test.begin(),test.end());
		bool wasUnique = (it == test.end() );
		assert(wasUnique);
	}


	{
		std::vector<int> test = Utility::original_permutation;
		std::sort(test.begin(), test.end());
		auto it = std::unique(test.begin(),test.end());
		bool wasUnique = (it == test.end() );
		assert(wasUnique);

	}

#endif
	int r_n=-1,r_popul = -1, r_eras=-1, r_mut=-1,r_sel=-1,r_minimization=-1;

	bool found_existing = Utility::minimization ? Utility::checkIfFIle(output+".minresume") : Utility::checkIfFIle(output+".maxresume");

	if(resume_from_file==1 || found_existing){


		ifstream ifs;

		if(resume_from_file==1)
			ifs.open(resume_file, std::ifstream::in);

		else{
			if(Utility::minimization)
				ifs.open(output+".minresume", std::ifstream::in);
			else
				ifs.open(output+".maxresume", std::ifstream::in);
		}




		if(!ifs.is_open())
		    throw runtime_error("Error Opening Resume File");

		std::cout<<"Reading Resume Info within File: " << resume_file << std::endl<<std::flush;


		std::string line;
		std::getline(ifs,line);
		std::vector<std::string> resulting_vector_strings;
		boost::split(resulting_vector_strings, line, boost::is_any_of("\t "),boost::token_compress_on);
		assert(resulting_vector_strings.size()==6);

		r_n = string2int(resulting_vector_strings[0]);
		r_popul = string2int(resulting_vector_strings[1]);
		r_eras = string2int(resulting_vector_strings[2]);
		r_sel = string2int(resulting_vector_strings[3]);
		r_mut = string2int(resulting_vector_strings[4]);
		r_minimization = string2int(resulting_vector_strings[5]);


		if(r_n!=poly.size())
			throw new std::runtime_error("RESUMING WRONG INDIVIDUALS IN INSTANCE");
		if(r_popul!=popul)
			throw new std::runtime_error("RESUMING WRONG POPULATION IN INSTANCE");
		if(r_eras!=Utility::eras)
			std::cout<<"WARNING ON ERAS\n";
		if(r_mut!=mut_prob)
			std::cout<<"WARNING ON MUTATION PROB\n";
		if(r_sel!=selection_percentage)
			std::cout<<"WARNING ON SELECTION PERC\n";
		if(r_minimization!=Utility::minimization)
			throw new std::runtime_error("RESUMING MAX FROM MIN INFO OR VICEVERSA");

		int read_population = 0;
		int read_evaluation = 0;
		Utility::population.clear();
		while(true){
			line.clear();
			std::getline(ifs,line);
			if(ifs.eof())
				break;

			std::vector<int> point_indices;
			resulting_vector_strings.clear();
			boost::split(resulting_vector_strings, line, boost::is_any_of("\t "),boost::token_compress_on);
			assert(resulting_vector_strings.size()==r_n+1);
			for(int i=0;i<r_n;i++)
				point_indices.push_back(string2int(resulting_vector_strings[i]));
			assert(point_indices.size()==r_n);
			Utility::population.push_back(std::make_tuple(point_indices, Utility::default_value,string2int(resulting_vector_strings[r_n])));
			++read_population;




		}
		assert(read_population==r_popul);
		assert(Utility::population.size()==r_popul);
		assert(ifs.eof());
		ifs.close();
		Utility::mutation_prob = mut_prob;
		Utility::selection_perc = selection_percentage;
	}
	else{
		Utility::mutation_prob = mut_prob;
		Utility::selection_perc = selection_percentage;
		std::cout<<"GENERATING INITIAL POPULATION..."<<std::flush;
		Utility::init_population(popul);
		std::cout<<"DONE!\n"<<std::flush;
	}

	Utility::bestcost = Utility::default_value;

	if(resume_from_file==1 || found_existing){

#pragma omp parallel for
		for(size_t t=0;t<Utility::population.size();t++){
			evaluate(Utility::population[t]);
		}
#pragma omp flush

		if(Utility::minimization)
			std::sort(Utility::population.begin(),Utility::population.end(),smallerPS());
		else
			std::sort(Utility::population.begin(),Utility::population.end(),greaterPS());

#ifndef NDEBUG
		if(Utility::minimization && Utility::bestcost<get<1>(Utility::population[0])){
			std::cout<<"PROB: "<<Utility::bestcost<<" "<<get<1>(Utility::population[0])<<"\n"<<std::flush;
			throw std::runtime_error("NON MONOTONICAL SOL");
		}
		if(!Utility::minimization && Utility::bestcost>get<1>(Utility::population[0])){
			std::cout<<"PROB: "<<Utility::bestcost<<" "<<get<1>(Utility::population[0])<<"\n"<<std::flush;
			throw std::runtime_error("NON MONOTONICAL SOL");
		}
#endif
		Utility::bestcost = get<1>(Utility::population[0]);
	}


	std::tuple<std::vector<int>,PointSet,areaType> ind1 = Utility::split_method(true);
	std::tuple<std::vector<int>,PointSet,areaType> ind2 = Utility::split_method(false);

	std::cout<<"STARTING BEST COST: "<<Utility::bestcost<<"\n"<<std::flush;
	std::cout<<"COSTS FOUND BY SPLIT METHOD: "<<get<2>(ind1)<<" "<<get<2>(ind2)<<"\n"<<std::flush;

	Utility::population.pop_back();
	Utility::population.insert(Utility::population.begin(),std::make_tuple(get<0>(ind1), get<2>(ind1),0));
	Utility::population.pop_back();
	Utility::population.insert(Utility::population.begin(),std::make_tuple(get<0>(ind2), get<2>(ind2),1));

	ind1 = Utility::concentric_hulls(true);
	ind2 = Utility::concentric_hulls(false);

	std::cout<<"COSTS FOUND BY HULL METHOD: "<<get<2>(ind1)<<" "<<get<2>(ind2)<<"\n"<<std::flush;

	Utility::population.pop_back();
	Utility::population.insert(Utility::population.begin(),std::make_tuple(get<0>(ind1), get<2>(ind1),2));
	Utility::population.pop_back();
	Utility::population.insert(Utility::population.begin(),std::make_tuple(get<0>(ind2), get<2>(ind2),3));

	ind1 = Utility::polar();
	std::cout<<"COST FOUND BY POLAR METHOD: "<<get<2>(ind1)<<"\n"<<std::flush;

	Utility::population.pop_back();
	Utility::population.insert(Utility::population.begin(),std::make_tuple(get<0>(ind1), get<2>(ind1),4));



	/*
	 * NOTE
	 *
	 * I want to modify the coordinates of a point, or the endpoints of a segment, but the CGAL kernel does not seem to support that. Why?
	 * Our kernel concept is representation-independent.
	 * There is no assumption on how a Point_3, for example, is represented.
	 * In particular, it need not be represented by Cartesian coordinates.
	 * Thus writing code that relies on being able to change these coordinates may lead to inefficient code.
	 * Our basic library algorithms are designed around the predicates that operate on the
	 * geometric objects and do not require this mutability. Non-modifiability also makes the
	 * underlying handle-rep scheme used (by default) to allocate, copy, etc., CGAL objects somewhat easier.
	 * We recognize, however, that users may want this flexibility and thus are working on providing
	 * mutable kernel objects.
	 * The only way to do this currently is to construct a new point and assign it to the old one : use p = Point_2(p.x(), p.y()+1); as if you would like to do p.y() += 1;.
	 */

	if(Utility::minimization)
		std::cout << "Performing genetic minimization algorithm for "<<Utility::eras<<" eras with population size "<<Utility::population.size()
		<<" selection perc "<<Utility::selection_perc<<" and mutation prob "<<Utility::mutation_prob<<"\n"<<std::flush;
	else
		std::cout << "Performing genetic maximization algorithm for " <<Utility::eras<<" eras with population size "<<Utility::population.size()
		<<" selection perc "<<Utility::selection_perc<<" and mutation prob "<<Utility::mutation_prob<<"\n"<<std::flush;


	mytimer t_era;
	for(size_t iterations=0;iterations<Utility::eras;iterations++){


		t_era.restart();
		assert(Utility::population.size()==popul);

		int threshold = Utility::selection(iterations,false);
//
		assert(Utility::population.size()==popul);
		assert(threshold<=Utility::population.size());
		threshold = threshold<=1 ? 2 : threshold;


#pragma omp parallel for

		for(int counter=threshold;counter<popul;counter+=2){
			//PARENTS AMONG BEST SOLUTIONS ONLY

			#ifdef TRACE
				std::cout<<"CROSSOVER -- ITERATIONS LEFT: "<<Utility::population.size()-counter<<"\n"<<std::flush;
			#endif

			int parent_id_one =  Utility::MersenneShuffler(threshold);
			int parent_id_two = Utility::MersenneShuffler(threshold);

			while(parent_id_two==parent_id_one)
				parent_id_two = Utility::MersenneShuffler(threshold);

			assert(parent_id_two>=0 && parent_id_one>=0);

			std::vector<int>& parent_one =  get<0>(Utility::population[parent_id_one]);
			std::vector<int>& parent_two =  get<0>(Utility::population[parent_id_two]);


			int split = 1 + (Utility::MersenneShuffler(parent_one.size()-1));
			assert(split>0);
			assert(parent_one.size()==parent_two.size());

			Utility::crossover(split,counter,parent_one,parent_two,true);
			Utility::crossover(split,counter,parent_one,parent_two,false);


		}
#pragma omp flush

		assert(Utility::population.size()==popul);

		Utility::mutation(threshold);
		std::cout<<"TIME PER ERA: "<<t_era.elapsed()<<"\n";

	}
	Utility::selection(Utility::eras,false);
	double elaps = tms.elapsed();
	Utility::write_solution(output,poly,elaps);
	Utility::write_snapshot(output,poly);
	std::cout<<"...done!\n"<<std::flush;
	if(Utility::minimization){
		std::cout<<"MIN_AREA: "<<Utility::bestcost<<"\n";
		std::cout<<"MIN AREA POINTSET: ";
	}

	else{
		std::cout<<"MAX AREA: "<<Utility::bestcost<<"\n";
		std::cout<<"MAX AREA POINTSET: ";
	}
	std::cout<<"[";
	for(size_t t=0;t<Utility::best_pointset.size();t++){
		std::cout<<"("<<Utility::best_pointset[t].x()<<","<<Utility::best_pointset[t].y()<<")";
		if(t<Utility::best_pointset.size()-1)
			std::cout<<",";
	}
	std::cout<<"]\n";
	std::cout<<"BACKTRACKS: "<<Utility::backtracked<<"\n";
	std::cout<<"WASTED EXECUTIONS: "<<Utility::wasted<<"\n";

	std::cout<<"TIME: "<<elaps<<"\n";


}
void Utility::init_population(int size_p){

	Utility::population.clear();


	ProgressStream generator(size_p);
	generator.label() << "Generating " <<size_p<< " individuals";

	for(size_t t=0;t<size_p;t++){
		std::vector<int> individual = Utility::original_permutation;
		std::random_shuffle(individual.begin(), individual.end(), MersenneShuffler);
		Utility::population.push_back(std::make_tuple(individual,Utility::default_value,5));
		++generator;
	}

	assert(Utility::population.size()==size_p);
}

int Utility::MersenneShuffler(int i){
	assert(i>0);
	return Utility::randomGenerator.getRandomInteger()%i;
}

int Utility::selection(int current_era, bool resize){

	double fraction = (double)Utility::selection_perc/(double)1000;
	assert(fraction>0 && fraction<=1);
	assert(round(Utility::population.size()*fraction)>0 && round(Utility::population.size()*fraction)<=Utility::population.size());


#pragma omp parallel for
	for(size_t t=0;t<Utility::population.size();t++){
		evaluate(Utility::population[t]);
	}
#pragma omp flush



	if(Utility::minimization)
		std::sort(Utility::population.begin(),Utility::population.end(),smallerPS());
	else
		std::sort(Utility::population.begin(),Utility::population.end(),greaterPS());

	assert(Utility::selection_perc>0 && Utility::selection_perc<=1000);

#ifndef NDEBUG
	if(Utility::minimization && Utility::bestcost<get<1>(Utility::population[0])){
		std::cout<<"PROB: "<<Utility::bestcost<<" "<<get<1>(Utility::population[0])<<"\n"<<std::flush;
		throw std::runtime_error("NON MONOTONICAL SOL");
	}
	if(!Utility::minimization && Utility::bestcost>get<1>(Utility::population[0])){
		std::cout<<"PROB: "<<Utility::bestcost<<" "<<get<1>(Utility::population[0])<<"\n"<<std::flush;
		throw std::runtime_error("NON MONOTONICAL SOL");
	}
#endif

	int valid = 0;
	double aver = 0;
	double stdv = 0;
	for(size_t t=0;t<Utility::population.size();t++){
//		Utility::areas[t]=get<1>(Utility::population[t]);
		if(get<1>(Utility::population[t])!=Utility::default_value){
			valid++;
			aver+=CGAL::to_double(get<1>(Utility::population[t]));
		}

	}
	aver/=valid;
	for(size_t t=0;t<Utility::population.size();t++){
		if(get<1>(Utility::population[t])!=Utility::default_value)
		stdv+=CGAL::to_double((get<1>(Utility::population[t])-aver)*(get<1>(Utility::population[t])-aver));
	}
	stdv /= valid;
	stdv = sqrt(stdv);

	if(Utility::minimization){
		if(get<1>(Utility::population[0])<Utility::bestcost)
			std::cout<<"MIN UPDATE BEST COST FROM: "<<Utility::bestcost<<" TO: "
			<<std::fixed<<std::setprecision(2)<<get<1>(Utility::population[0])<<" DURING ERA #: "<<current_era
			<<" AVG SOLUTION: "<<aver<<" STDDEV: "<<stdv<<" VALID: "<<valid<<"\n"<<std::flush;
		else
			std::cout<<"MIN NO UPDATE BEST COST DURING ERA #: "<<current_era
			<<std::fixed<<std::setprecision(2)<<" AVG SOLUTION: "<<aver<<" STDDEV: "<<stdv<<" VALID: "<<valid<<"\n"<<std::flush;

	}
	else{
		if(get<1>(Utility::population[0])>Utility::bestcost)
			std::cout<<"MAX UPDATE BEST COST FROM: "<<Utility::bestcost<<" TO: "
			<<std::fixed<<std::setprecision(2)<<get<1>(Utility::population[0])<<" DURING ERA #: "<<current_era
			<<" AVG SOLUTION: "<<aver<<" STDDEV: "<<stdv<<" VALID: "<<valid<<"\n"<<std::flush;
		else
			std::cout<<"MAX NO UPDATE BEST COST DURING ERA #: "<<current_era
			<<std::fixed<<std::setprecision(2)<<" AVG SOLUTION: "<<aver<<" STDDEV: "<<stdv<<" VALID: "<<valid<<"\n"<<std::flush;

	}

	Utility::evolution.push_back(std::make_tuple(Utility::bestcost,aver,stdv));

	Utility::bestcost = get<1>(Utility::population[0]);
	Utility::best_pointset.clear();
	genotype2phenotype(get<0>(Utility::population[0]),Utility::best_pointset,get<2>(Utility::population[0]));

//#ifndef NDEBUG
//	if(Utility::minimization){
//		for(size_t t=0;t<Utility::population.size();t++)
//			if(Utility::bestcost>Utility::population[t].second)
//				throw std::runtime_error("WRONG COMPARATOR MIN");
//
//	}
//	else{
//		for(size_t t=0;t<Utility::population.size();t++)
//			if(Utility::bestcost<Utility::population[t].second)
//				throw std::runtime_error("WRONG COMPARATOR MAX");
//	}
#ifdef TRACE

	std::cout<<"CURRENT BEST POINTSET: ";
	print_point_set(Utility::best_pointset);
#endif





	if(!resize)
		return round(Utility::population.size()*fraction);

	//DROP BAD SOLUTIONS

	Utility::population.resize(round(Utility::population.size()*fraction));
	return round(Utility::population.size()*fraction);

}

void Utility::print_point_set(PointSet& PS){
	std::cout<<"["<<std::flush;
	for(size_t i=0;i<PS.size();i++)
		if(i<PS.size()-1)
			std::cout<<"("<<PS[i].x()<<","<<PS[i].y()<<"),"<<std::flush;
		else
			std::cout<<"("<<PS[i].x()<<","<<PS[i].y()<<")"<<std::flush;
	std::cout<<"]\n"<<std::flush;

}

void Utility::add_in_position(PointSet& p,int index,Point pt){
	p.insert(p.begin()+index+1, pt);
}

void Utility::add_in_position(std::vector<Segment_2> & s,int index,Segment_2 seg){
	s.insert(s.begin()+index+1, seg);
}

areaType Utility::genotype2phenotype(std::vector<int>& perm,int & type_of_ind){

	PointSet ps;

	return genotype2phenotype(perm,ps,type_of_ind);

}
areaType Utility::compute_triangle_area(Point& a, Point& b, Point& c){

	return fabs(0.5 * CGAL::to_double( a.x()*(b.y()-c.y()) + b.x()*(c.y()-a.y()) +  c.x()*(a.y()-b.y())));

}
areaType Utility::compute_polygon_area(std::vector<Point>& P){

	areaType contrib = 0;
	for(size_t t=0;t<P.size();t++){
		contrib += P[t].x()*P[(t+1)%P.size()].y() - P[t].y()*P[(t+1)%P.size()].x();
	}
	return fabs(0.5 * contrib);

}



areaType Utility::genotype2phenotype(std::vector<int>& perm,PointSet &pointset,int & type_of_ind){

 /*
  * i primi tre punti sono un triangolo oppure sono collineari, aggiunge un punto al poligono semplice dato da una sequenza di punti
  * il punto viene aggiunto provando dalla fine, poi in penultima posizione, etc., finché la sequenza non è un poligono
  */


	if(type_of_ind==0){
		std::tuple<std::vector<int>,PointSet,areaType> ind = Utility::split_method(true);
		return get<2>(ind);
	}

	if(type_of_ind==1){
		std::tuple<std::vector<int>,PointSet,areaType> ind = Utility::split_method(false);
		return get<2>(ind);
	}
	if(type_of_ind==2){
		std::tuple<std::vector<int>,PointSet,areaType> ind = Utility::concentric_hulls(true);
		return get<2>(ind);
	}
	if(type_of_ind==3){
		std::tuple<std::vector<int>,PointSet,areaType> ind = Utility::concentric_hulls(false);
		return get<2>(ind);
	}
	if(type_of_ind==4){
		std::tuple<std::vector<int>,PointSet,areaType> ind = Utility::polar();
		return get<2>(ind);
	}

	assert(type_of_ind==5);
	assert(pointset.size()==0);

	std::vector<Segment_2> current_segments;

	pointset.push_back(Utility::index2point(perm[0]));
	pointset.push_back(Utility::index2point(perm[1]));
	pointset.push_back(Utility::index2point(perm[2]));






#ifndef NDEBUG
	Polygon_2 dbg_poly = pointset2polygon(pointset);
	assert(CGAL::collinear(pointset[0],pointset[1],pointset[2]) || dbg_poly.is_simple());
#endif

	int offset = 0;

	if(CGAL::collinear(pointset[offset%perm.size()],pointset[(offset+1)%perm.size()],pointset[(offset+2)%perm.size()])){

		std::sort(pointset.begin(),pointset.end(),XY_Sorter());
		assert(pointset.size()==3);
		for(size_t t=0;t<pointset.size()-1;t++)
			current_segments.push_back(Segment_2(pointset[t],pointset[t+1]));
		while(true){

			assert(offset<perm.size());
			pointset.push_back(Utility::index2point(perm[offset+3]));
			if(!CGAL::collinear(pointset[offset+1],pointset[offset+2],pointset[offset+3])){
				current_segments.push_back(Segment_2(pointset[offset+2],pointset[offset+3]));
				current_segments.push_back(Segment_2(pointset[offset+3],pointset[0]));
				offset++;

				break;
			}
			else{
				std::sort(pointset.begin(),pointset.end(),XY_Sorter());
				current_segments.clear();
				for(size_t t=0;t<pointset.size()-1;t++)
					current_segments.push_back(Segment_2(pointset[t],pointset[t+1]));

				offset++;


			}
		}




	}
	else{
		current_segments.push_back(Segment_2(pointset[0],pointset[1]));
		current_segments.push_back(Segment_2(pointset[1],pointset[2]));
		current_segments.push_back(Segment_2(pointset[2],pointset[0]));

	}

	assert(current_segments.size()==pointset.size());

#ifndef NDEBUG
	dbg_poly = pointset2polygon(pointset);
	assert(dbg_poly.is_simple());
#endif

	unsigned long int backtracks = 0;

	for(int iterations=3+offset;iterations<perm.size();){

#ifdef TRACE
		std::cout<<"Status:\n";
		std::cout<<"PointSet: "<<std::endl;
		print_point_set(pointset);
		std::cout<<"SegmentSet: "<<std::endl;
		for(size_t t=0;t<current_segments.size();t++)
			std::cout<<current_segments[t].source()<<" -> "<<current_segments[t].target()<<"\n";

#endif
		assert(pointset.size()==iterations);

		Point r = index2point(perm[iterations]);

		int best_triangle = -1;

		Polygon_2 current_polygon = pointset2polygon(pointset);

#ifdef TRACE
		std::cout<<"-- TRYING INSERTION OF POINT: "<<r<<"\n"<<std::flush;
#endif
		assert(current_polygon.is_simple());


		if(current_polygon.has_on_bounded_side(r)){
			// point is in the interior
			//Seek for largest triangle if minimization and viceversa


			areaType target_area  = Utility::minimization ? 0.0 : std::numeric_limits<areaType>::max();

			for(int j=0;j<iterations;j++){

				int outcome = Utility::test_insertion_between(current_segments,r,pointset[j],pointset[(j+1)%pointset.size()]);
				if(outcome!=0){
					if(outcome==-1)
						j++;
					continue;
				}




				areaType triangle_area = compute_triangle_area(pointset[j],r,pointset[(j+1)%pointset.size()]);

				assert(triangle_area>0);
				assert(triangle_area<std::numeric_limits<areaType>::max());



				if(Utility::minimization){
					if(triangle_area>target_area){
						target_area=std::max(triangle_area, target_area);
						best_triangle=j;
						assert(current_segments[best_triangle].source()==pointset[best_triangle]);
						assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);
						continue;
					}
				}
				else{
					if(triangle_area<target_area){
						target_area=std::min(triangle_area, target_area);
						best_triangle=j;
						assert(current_segments[best_triangle].source()==pointset[best_triangle]);
						assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);
						continue;

					}

				}
			}

		}
		else if(current_polygon.has_on_unbounded_side(r)){

			//Seek for smaller triangle if minimization and viceversa
			areaType target_area = Utility::minimization ? std::numeric_limits<areaType>::max() : 0;



			for(int j=0;j<iterations;j++){

				int outcome = Utility::test_insertion_between(current_segments,r,pointset[j],pointset[(j+1)%pointset.size()]);
				if(outcome!=0){
					if(outcome==-1)
						j++;
					continue;
				}

				areaType triangle_area = compute_triangle_area(pointset[j],r,pointset[(j+1)%pointset.size()]);

				assert(triangle_area>0);
				assert(triangle_area<std::numeric_limits<areaType>::max());

				if(Utility::minimization){
					if(triangle_area<target_area){
						target_area=std::min(triangle_area, target_area);
						best_triangle=j;
						assert(current_segments[best_triangle].source()==pointset[best_triangle]);
						assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);
						continue;
					}
				}

				else{
					if(triangle_area>target_area){
						target_area=std::max(triangle_area, target_area);
						best_triangle=j;
						assert(current_segments[best_triangle].source()==pointset[best_triangle]);
						assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);
						continue;
					}
				}

			}

		}
		else{
			assert(current_polygon.has_on_boundary(r));
			for(int k=0;k<current_segments.size();k++){
				if(current_segments[k].has_on(r)){
					best_triangle=k;
					assert(current_segments[best_triangle].source()==pointset[best_triangle]);
					assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);
					break;
				}
			}
			assert(best_triangle>=0);

		}


		if(best_triangle==-1){

//
//			std::cout<<"[";
//			for(size_t t=0;t<pointset.size();t++){
//				std::cout<<"("<<pointset[t].x()<<","<<pointset[t].y()<<")";
//				if(t<pointset.size()-1)
//					std::cout<<",";
//			}
//			std::cout<<"]\n";
//
//			std::cout<<"[";
//			for(size_t t=0;t<perm.size();t++){
//				std::cout<<"("<<Utility::index2point(perm[t]).x()<<","<<Utility::index2point(perm[t]).y()<<")";
//				if(t<perm.size()-1)
//					std::cout<<",";
//			}
//			std::cout<<"]\n";

			backtracks++;
			Utility::backtracked++;

#ifdef TRACE
			std::cout<<"LINEAR BACKTRACKING: "<<backtracks<<" ITERATIONS: "<<iterations<<"\n"<<std::flush;

#endif
			if(backtracks>=std::pow(perm.size(), 1)){
				Utility::wasted++;
				return Utility::default_value;
			}

			int num_removals = Utility::phenotyper == 0 ? 1 : 1+Utility::randomGenerator.getRandomInteger()%(iterations - 1);

			if(!((Utility::phenotyper!=0 && num_removals>=1) ||  (Utility::phenotyper==0 && num_removals==1)))
				throw std::runtime_error("something is wrong here");

			while(num_removals>0){

#ifdef TRACE
				std::cout<<"REMOVING POINT: "<<Utility::index2point(perm[position_to_remove])<<"\n";
#endif

//
//				if(iterations<3+offset)
//					throw std::runtime_error("Something is wrong here --");


				int k=0;
				for(;k<current_segments.size();k++){
					if(current_segments[k].target()==Utility::index2point(perm[iterations-1]))
						break;
				}

				assert(std::find(pointset.begin(),pointset.end(),Utility::index2point(perm[iterations-1]))!=pointset.end());
				assert(pointset[k+1]==Utility::index2point(perm[iterations-1]));
				pointset.erase(pointset.begin()+k+1);

				assert(std::find(pointset.begin(),pointset.end(),Utility::index2point(perm[iterations]))==pointset.end());
				assert(std::find(pointset.begin(),pointset.end(),Utility::index2point(perm[iterations-1]))==pointset.end());


				assert(current_segments[k].target()==current_segments[k+1].source());

				assert(k<current_segments.size());

#ifdef TRACE
				std::cout<<"SEGMENTS BEFORE: "<<std::endl<<std::flush;
				for(size_t t=0;t<current_segments.size();t++)
					std::cout<<current_segments[t].source()<<" -> "<<current_segments[t].target()<<"\n";
#endif
				Segment_2 replacement = Segment_2(current_segments[k].source(),current_segments[k+1].target());

				current_segments.erase(current_segments.begin()+k);
				current_segments.erase(current_segments.begin()+k);

				add_in_position(current_segments, k-1, replacement);

#ifdef TRACE
				std::cout<<"SEGMENTS AFTER: "<<std::endl<<std::flush;
				for(size_t t=0;t<current_segments.size();t++)
					std::cout<<current_segments[t].source()<<" -> "<<current_segments[t].target()<<"\n";
#endif

				std::swap(perm[iterations],perm[iterations-1]);

				iterations--;
				num_removals--;
			}
			continue;
		}
		else{

			assert(best_triangle>=0);

#ifdef TRACE
			std::cout<<"INSERTION SUCCEEDED\n";

#endif
			assert(current_segments[best_triangle].source()==pointset[best_triangle]);
			assert(current_segments[best_triangle].target()==pointset[(best_triangle+1)%pointset.size()]);

			current_segments.erase(current_segments.begin()+best_triangle);
			add_in_position(current_segments, best_triangle-1, Segment_2(pointset[best_triangle],r));
			add_in_position(current_segments, best_triangle, Segment_2(r,pointset[(best_triangle+1)%pointset.size()]));

			Utility::add_in_position(pointset,best_triangle, r);

#ifndef NDEBUG
			dbg_poly = pointset2polygon(pointset);
			assert(dbg_poly.is_simple());
#endif
			iterations++;
		}


	}

	assert(pointset.size()==perm.size());

#ifndef NDEBUG
	Polygon_2 polygon = pointset2polygon(pointset);
	if(!polygon.is_simple())
		throw std::runtime_error("NON SIMPLE INDIVIDUAL");
	assert(Utility::compute_polygon_area(pointset)>0);
#endif

	return Utility::compute_polygon_area(pointset);
}
void Utility::evaluate(std::tuple<std::vector<int>,areaType,int>& individual){
	get<1>(individual) = genotype2phenotype(get<0>(individual),get<2>(individual));
}

void Utility::crossover(int& split, int & index, std::vector<int>& parent_one,std::vector<int>& parent_two,bool first_born){


	if(index>=Utility::population.size()){
		assert(!first_born);
		return;
	}

	std::vector<int>& child =  get<0>(Utility::population[index]);
	child.clear();

	if(first_born){
		//FINO A SPLIT DAL PARENT UNO
		std::copy(parent_one.begin(),parent_one.begin()+split,std::back_inserter(child));

		//DA SPLIT IN POI DAL PARENT DUE QUELLO CHE NON STA NELLA PRIMA PARTE
		std::copy_if (parent_two.begin()+split,parent_two.end(), std::back_inserter(child),
					  [child](int i){return std::find(child.begin(),child.end(),i)==child.end();} );
		//PARTE RIMANENTE FILLATA CON ELEMENTI MANCANTI DI PARENT DUE
		std::copy_if (parent_two.begin(),parent_two.begin()+split, std::back_inserter(child),
					  [child](int i){return std::find(child.begin(),child.end(),i)==child.end();} );


	}

	else{
		std::copy(parent_two.begin(),parent_two.begin()+split,std::back_inserter(child));

		std::copy_if (parent_one.begin()+split, parent_one.end(), std::back_inserter(child),
				[child](int i){return std::find(child.begin(),child.end(),i)==child.end();} );

		std::copy_if (parent_one.begin(), parent_one.begin()+split, std::back_inserter(child),
				[child](int i){return std::find(child.begin(),child.end(),i)==child.end();} );

	}

#ifndef NDEBUG
	assert(child.size()==parent_two.size());

	{
		std::vector<int> test = child;
		std::sort(test.begin(), test.end());
		auto it = std::unique(test.begin(),test.end());
		bool wasUnique = (it == test.end() );
		assert(wasUnique);

	}

#endif



}

void Utility::mutation(int renewing_threshold){

	assert(Utility::mutation_prob>0 && Utility::mutation_prob<=1000);

	//MUTATE BAD SOLUTIONS ONLY - KEEP THE OTHERS UNCHANGED
#ifndef NDEBUG
	double fraction = (double)Utility::selection_perc/(double)1000;
	assert(fraction>0 && fraction<=1);
	int best_indices = round(Utility::population.size()*fraction);
	best_indices = best_indices<=1 ? 2 : best_indices;

	assert(renewing_threshold==best_indices);
	assert(renewing_threshold>1);
#endif

	for(size_t t=renewing_threshold;t<Utility::population.size();t++){

		assert(get<0>(Utility::population[t]).size()>0);

		if(1+MersenneShuffler(1000)<=Utility::mutation_prob){

			int counter = 1+round((double)get<0>(Utility::population[t]).size()/100);

//			assert(counter>=1);
//			while(counter>0){
				int to_swap_one = MersenneShuffler(get<0>(Utility::population[t]).size());
				int to_swap_two = MersenneShuffler(get<0>(Utility::population[t]).size());

				while(to_swap_one==to_swap_two)
					to_swap_two = MersenneShuffler(get<0>(Utility::population[t]).size());

				assert(to_swap_one>=0 && to_swap_two>=0);

				std::swap(get<0>(Utility::population[t])[to_swap_one], get<0>(Utility::population[t])[to_swap_two]);
				counter--;
//			}
		}

	}


}

PointSet Utility::polygon2pointset(Polygon_2& p){
	PointSet toreturn;
	for(size_t t=0;t<p.size();t++)
		toreturn.push_back(p[t]);
	return toreturn;

}
Polygon_2 Utility::pointset2polygon(PointSet& p){
	Polygon_2 toreturn;
	for(size_t t=0;t<p.size();t++)
		toreturn.push_back(p[t]);
	return toreturn;

}

Point Utility::index2point(int index){
	return Utility::original_pointset[index];

}
int Utility::point2index(Point p){
	for(size_t t=0;t<Utility::original_pointset.size();t++)
		if(Utility::original_pointset[t]==p)
			return t;
	throw std::runtime_error("not found point2index");
}

PointSet Utility::permutation2pointset(std::vector<int>& perm){
	PointSet to_build;
	assert(perm.size()<=Utility::original_pointset.size());
	for(size_t t=0;t<perm.size();t++)
		to_build.push_back(index2point(perm[t]));
	assert(to_build.size()<=Utility::original_pointset.size());
	return to_build;

}
std::vector<int> Utility::pointset2permutation(PointSet& ps){
	std::vector<int> to_build;
	assert(ps.size()<=Utility::original_pointset.size());
	for(size_t t=0;t<ps.size();t++)
		to_build.push_back(point2index(ps[t]));
	assert(to_build.size()<=Utility::original_pointset.size());
	return to_build;

}
int Utility::test_insertion_between(std::vector<Segment_2> & P,Point & point_to_add,
		Point & first_polygon_point,Point & second_polygon_point){


	Segment_2 first_to_test = Segment_2(point_to_add,first_polygon_point);
	Segment_2 second_to_test = Segment_2(point_to_add,second_polygon_point);

	assert(point_to_add!=first_polygon_point);
	assert(first_polygon_point!=second_polygon_point);

	for(int t=0;t<P.size();t++){


		CGAL::Object first_result = CGAL::intersection(first_to_test,P[t]);
		if (const Point * ipoint = CGAL::object_cast<Point>(&first_result)){
			// handle the point intersection case
			assert(*ipoint!=point_to_add);
			if((*ipoint).x()!=first_polygon_point.x() || (*ipoint).y()!=first_polygon_point.y()){
#ifdef TRACE
				std::cout<<"INTERSECTION OF FIRST SEGMENT "<<"("<<point_to_add<<","<<first_polygon_point<<")\n";
				std::cout<<"WITH SEGMENT "<<"("<<P[t].source()<<","<<P[t].target()<<")\n";
				std::cout<<"INTERSECTION POINT: "<<(*ipoint).x()<<" "<<(*ipoint).y()<<"\n";

#endif
				return 1;


			}
		}
		else if (const Segment_2 * iseg = CGAL::object_cast<Segment_2>(&first_result)){
			// handle the segment intersection case
#ifdef TRACE

			std::cout<<"INTERSECTION OF FIRST SEGMENT-SEG\n";
			std::cout<<"("<<point_to_add<<","<<first_polygon_point<<"),"<<"("<<P[t].source()<<","<<P[t].target()<<")\n";
#endif
			return 1;
		}


		CGAL::Object second_result = CGAL::intersection(second_to_test,P[t]);
		if (const Point * ipoint = CGAL::object_cast<Point>(&second_result)){
			// handle the point intersection case
			assert(*ipoint!=point_to_add);
			if((*ipoint).x()!=second_polygon_point.x() || (*ipoint).y()!=second_polygon_point.y()){

#ifdef TRACE
				std::cout<<"INTERSECTION OF SECOND SEGMENT "<<"("<<point_to_add<<","<<second_polygon_point<<")\n";
				std::cout<<"WITH SEGMENT "<<"("<<P[t].source()<<","<<P[t].target()<<")\n";
				std::cout<<"INTERSECTION POINT: "<<(*ipoint).x()<<" "<<(*ipoint).y()<<"\n";
#endif
				return -1;
			}
		}
		else if (const Segment_2 * iseg = CGAL::object_cast<Segment_2>(&second_result)){
			// handle the segment intersection case
#ifdef TRACE
				std::cout<<"INTERSECTION OF SECOND SEGMENT-SEG\n";
				std::cout<<"("<<point_to_add<<","<<second_polygon_point<<"),"<<"("<<P[t].source()<<","<<P[t].target()<<")\n";
#endif
			return -1;
		}


	}
	return 0;
}

double Utility::string2double(std::string& s){
	return atof(s.c_str());
};

int Utility::string2int(std::string& s){
	return atoi(s.c_str());
};
//
//areaType Utility::stddev(std::vector<areaType>& v){
//  double average = Utility::avg(v);
//  double variation = 0;
//  double nonnull = 0;
//
//  for(size_t t = 0; t < v.size();t++)
//	  if(v[t]!=Utility::default_value){
//		  variation+=((v[t]-average)*(v[t]-average));
//			nonnull++;
//
//	  }
//
//  variation/=nonnull;
//  return sqrt(variation);
//
//}
//areaType Utility::avg(std::vector<areaType>& v){
//	areaType value = 0;
//	double nonnull = 0;
//	for(size_t t = 0; t < v.size();t++){
//		if(v[t]!=Utility::default_value){
//			value+=v[t];
//			nonnull++;
//		}
//	}
//	if(nonnull==0)
//		return Utility::default_value;
//	value/=nonnull;
//	return value;
//}
bool Utility::checkIfFIle(std::string filePath)
{
	try {
		// Create a Path object from given path string
		 boost::filesystem::path pathObj(filePath);
		// Check if path exists and is of a regular file
		if (boost::filesystem::exists(pathObj) && boost::filesystem::is_regular_file(pathObj))
			return true;
	}
	catch (boost::filesystem::filesystem_error & e)
	{
		std::cerr << e.what() << std::endl;
	}
	return false;
}

