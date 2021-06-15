#ifndef PROB_DENSITY_H
#define PROB_DENSITY_H

/**
 * A helper class for getting random numbers from a gaussian.
 */
class GaussianPDF
{
	private:
		std::mt19937 gen;
		std::normal_distribution<double> d;

	public:
		GaussianPDF() : gen(std::random_device{}()), d(std::normal_distribution<>{0,1}) {} 
		GaussianPDF(double mu, double sigma) : gen(std::random_device{}()), d(std::normal_distribution<double>{mu,sigma}) {} 

		double draw() {return(d(gen));}
};



/**
 * A helper class for getting random numbers from a bimodal distribution.
 */
class BimodalPDF
{
	private:
		std::mt19937 gen;
		std::discrete_distribution<int> d;
	public:
		BimodalPDF() : gen(std::random_device{}()), d(std::discrete_distribution<int>{1,0,1}) {}
		BimodalPDF(double weight) : gen(std::random_device{}()), d(std::discrete_distribution<int>{weight,0,1-weight}) {}

		int draw() {return(d(gen)-1);}
};

#endif
