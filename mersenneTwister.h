#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>

class MersenneTwister
{
public:
    MersenneTwister():generator(time(NULL)),distribution(0, std::numeric_limits<int>::max()),die(generator,distribution)
    {}

    int getRandomInteger()
    {
        return die();
    }

private:
    boost::mt19937 generator;
    boost::random::uniform_int_distribution<> distribution;
    boost::variate_generator<boost::mt19937&, boost::random::uniform_int_distribution<>  > die;
};



#endif //MERSENNETWISTER_H
