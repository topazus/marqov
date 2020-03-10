#ifndef RNDWRAPPER_H
#define RNDWRAPPER_H
#include <random>

// Class which essentially sets up a random number generator containing both integer and double distribution
// mainly this is for convenience as only one object has to be passed to functions

// For performance critical parts of the code, initiate one "global" instance of this class
// and pass it through your functions by reference
// For non-performance critical tasks, one may initiate a seperate instance locally

class RND 
{
	private:
		#ifdef RANDOM_64_BIT	
			std::mt19937_64 gen;
		#else
			std::mt19937 gen;
		#endif

		std::uniform_real_distribution<double> dD;
		std::uniform_int_distribution<> dI;
	public:
		bool integer_distribution_set = false;
		void seed(int s) { gen.seed(s); }

		// constructor, initializes double distribution
		RND(double dmin, double dmax) : dD(dmin, dmax) {}

		// set integer distribtion limits
		void set_integer_range(int imax) 
		{ 
			dI = std::uniform_int_distribution<>(0, imax-1);
			integer_distribution_set = true;
		}
		double d() { return dD(gen); }
		double i() { return dI(gen); }
};

// when creating an RND object (by calling the constructor) the range of doubles is set:
// RND myrnd(0,1);

// the range of the integer generator is not set then, but can be done with
// myrnd.set_integer_range(L);

// now double from [0,1) can be created by calling myrnd.d()
// and integers from [0,L-1] by calling myrnd.i()

#endif
