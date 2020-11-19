#include "Utility.h"
using namespace std;

int main(int argc, char** argv) {

	//declare supported options
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	("input,i", po::value<std::string>(), "Input File Path")
	("resume_from_file,r", po::value<int>(), "Resume From File [0: no 1: yes]")
	("resume_file_path,f", po::value<std::string>(), "Resume File Path")
	("min_max,x", po::value<int>(), "[0: minimization, 1: maximization]")
	("phenotyper,g", po::value<int>(), "[0: rando, 1: jea 2: jea determ backtracking 3: jea random backtracking]")
	("population_size,p", po::value<int>(), "Pool Size")
	("eras,e", po::value<int>(), "Number of Eras Before Convergence")
	("mutation_probability,m", po::value<int>(), "Probability of Mutation (in permille)")
	("selection_percentage,c", po::value<int>(), "Percentage of Selection (in permille)")

	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);


	std::string input= "",resume_file_path = "";
	int min_max = -1;
	int resume_from_file = -1;
	int population_size = -1;
	int eras = -1;
	int mutation_probability = -1;
	int selection_percentage = -1;
	int phenotyper = -1;

	if (vm.empty()){
		std::cout << desc << "\n";
		throw std::runtime_error("Empty options");
	}

	if(vm.count("input"))
		input = vm["input"].as<std::string>();
	if(vm.count("resume_file_path"))
		resume_file_path = vm["resume_file_path"].as<std::string>();
	if(vm.count("min_max"))
		min_max = vm["min_max"].as<int>();
	if(vm.count("resume_from_file"))
		resume_from_file = vm["resume_from_file"].as<int>();
	if(vm.count("population_size"))
		population_size = vm["population_size"].as<int>();
	if(vm.count("eras"))
		eras = vm["eras"].as<int>();
	if(vm.count("mutation_probability"))
		mutation_probability = vm["mutation_probability"].as<int>();
	if(vm.count("selection_percentage"))
		selection_percentage = vm["selection_percentage"].as<int>();
	if(vm.count("phenotyper"))
		phenotyper = vm["phenotyper"].as<int>();
	if(input=="")
		throw runtime_error("empty input path");

	if(resume_from_file!=0 && resume_from_file!=1)
		throw runtime_error("wrong resume_from_file selection");

	if(resume_file_path=="" && resume_from_file==1)
		throw runtime_error("empty resume_file_path");

	if(min_max!=0 && min_max!=1)
		throw runtime_error("wrong min_max selection");
	if(phenotyper>3 || phenotyper<0)
		throw runtime_error("wrong phenotyper");
	if(population_size<2)
		throw runtime_error("wrong population_size selection");

	if(eras<1)
		throw runtime_error("wrong eras selection");

	if(mutation_probability<=0 || mutation_probability>1000)
		throw runtime_error("wrong mutation_probability selection");

	if(selection_percentage<=0 || selection_percentage>1000)
		throw runtime_error("wrong selection_percentage selection");

	PointSet poly;
	Utility::read_instance(input, poly);
	size_t lastindex = input.find_last_of(".");
	std::string output = input.substr(0, lastindex);






	if(min_max==0){
		std::cout<<"Attacking Min via Genetic Algorithm -- Instance "<<output<<"\n"<<std::flush;
		Utility::genetic(poly,output,true,population_size,eras,mutation_probability,selection_percentage,resume_from_file,resume_file_path,phenotyper);

	}
	else{
		assert(min_max==1);
		std::cout<<"Attacking Max via Genetic Algorithm -- Instance "<<output<<"\n"<<std::flush;
		Utility::genetic(poly,output,false,population_size,eras,mutation_probability,selection_percentage,resume_from_file,resume_file_path,phenotyper);

	}






	return EXIT_SUCCESS;
};


