#ifndef OBSPARTS_H
#define OBSPARTS_H

#include <complex>
#include <string>
#include <vector>


// Scalar Magnetization Observable
// considers only the first component of the state vector
class ScalarMagnetization
{
	public:
		std::string name, desc;

		ScalarMagnetization(std::string name, std::string description) : name(name), desc(description) {}
		ScalarMagnetization(std::string name) : name(name), desc("scalar magnetization") {}
		ScalarMagnetization() : name("mx"), desc("scalar magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			double retval = 0.0;
			for (int i=0; i<grid.size(); i++) retval += statespace[i][0];
			return std::abs(retval)/double(grid.size());
		}
};

// Magnetization
// Euclidean norm of the sum of state vectors
class Magnetization
{
	public:
		std::string name, desc;

		Magnetization(std::string name, std::string description) : name(name), desc(description) {}
		Magnetization(std::string name) : name(name), desc("magnetization") {}
		Magnetization() : name("m"), desc("magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			std::vector<double> mag(SymD, 0);

			for (int i=0; i<grid.size(); i++) 
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j]; 
				}
			}

			double retval = 0;
			for (int j=0; j<SymD; j++) retval += mag[j]*mag[j];

			return sqrt(retval)/double(grid.size());
		}
};

// Vector Magnetization
// sum of every component of the state vectors
class VectorMagnetization
{
	public:
		std::string name, desc;

		VectorMagnetization(std::string name, std::string description) : name(name), desc(description) {}
		VectorMagnetization(std::string name) : name(name), desc("vector magnetization") {}
		VectorMagnetization() : name("mvec"), desc("vector magnetization") {}

		template <class StateSpace, class Grid>
		std::vector<double> measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			std::vector<double> mag(SymD, 0);

			for (int i=0; i<grid.size(); i++) 
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j]; 
				}
			}
			return mag;
		}
};

// Scalar Magnetization Fourier Component
// considers only the first component of the state vector
class ScalarMagFTComp
{
	public:
		int dir;
		std::string name, desc;

		ScalarMagFTComp(int dir, std::string name, std::string description) : dir(dir), name(name), desc(description) {}
		ScalarMagFTComp(int dir, std::string name) : dir(dir), name(name), desc("Fourier Component of Magnetization") {}
		ScalarMagFTComp(int dir) : dir(dir), name("scalarmagft"+std::to_string(dir)), desc("Fourier Component of Magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			std::complex<double> retval = 0.0;
			std::complex<double> jj(0,1);

			for (int i=0; i<grid.size(); i++)
			{
				const double x = grid.crds(i)[dir];
				retval += double(statespace[i][0]) * std::exp(2.0*M_PI*x*jj);
			}

			return std::pow(std::abs(retval/double(grid.size())), 2);
		}

};





class MagFTComp
{
	public:
		int dir;
		std::string name, desc;
		MagFTComp(int dir, std::string name, std::string description) : dir(dir), name(name), desc(description) {}
		MagFTComp(int dir, std::string name) : dir(dir), name(name), desc("Fourier Component of Magnetization") {}
		MagFTComp(int dir) : dir(dir), name("magft"+std::to_string(dir)), desc("Fourier Component of Magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			const int N = grid.size();

			std::vector<std::complex<double>> magFTcomp(SymD,0);
			const std::complex<double> jj(0,1);

			for (int i=0; i<N; i++)
			{
				double x = grid.crds(i)[dir];
				for (int j=0; j<SymD; j++) magFTcomp[j] += double(statespace[i][j]) * std::exp(2.0*M_PI*x*jj);
			}

			// normalize
			for (int j=0; j<SymD; j++) magFTcomp[j] /= double(N);

			// dot product of complex vector 
			double retval = 0;
			for (int j=0; j<SymD; j++) retval += std::norm(magFTcomp[j]);

			return retval;
		}
};







// Interaction Energy Observable
template <class Hamiltonian>
class InteractionEnergy
{
	private:
		Hamiltonian& ham;

	public:
		InteractionEnergy (Hamiltonian& ham) : ham(ham), name("eint")  {};

		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			double ene = 0.0;

			// interaction part
			for (int a=0; a<ham.interactions.size(); a++)
			{

				double enepart = 0.0;
				for (int idx=0; idx<N; idx++)
				{
					auto nbrs = grid.nbrs(a, idx);
					auto self = ham.interactions[a]->get(statespace[idx]);

					for (std::size_t i = 0; i < nbrs.size(); ++i)
					{
						// index of the neighbour
						auto nbridx = nbrs[i];

						// configuration of the neighbour
						auto nbr = ham.interactions[a]->get(statespace[nbridx]);

						enepart += dot(self,nbr);

					}
				}

				ene += ham.interactions[a]->J * enepart;
			}
			return ene/double(N);
		}
};



// Self Energy Observable
template <class Hamiltonian>
class SelfEnergy
{
	private:
		Hamiltonian& ham;

	public:
		SelfEnergy (Hamiltonian& ham) : ham(ham), name("eself")  {};

		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			double ene = 0.0;

			// self energy part
			for (int b=0; b < ham.onsite.size(); b++)
			{
				double enepart = 0.0;

				for (int idx=0; idx<N; idx++)
				{
					enepart += ham.onsite[b]->get(statespace[idx]);
				}
				ene += ham.onsite[b]->h * enepart;
			}
			return ene/double(N);
		}
};


// Flex Energy Observable
template <class Hamiltonian>
class FlexEnergy
{
	private:
		Hamiltonian& ham;

	public:
		FlexEnergy (Hamiltonian& ham) : ham(ham), name("eflex")  {};

		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			for (int c=0; c<ham.Ngamma; c++)
			{
//				enepart = ham.multisite[c]->energy();
//				ene += ham.onsite[c]->k * enepart;
			}
			return 0;
		}
};






// Full Energy Observable
template <class Hamiltonian>
class Energy
{
	private:
		Hamiltonian& ham;
		InteractionEnergy<Hamiltonian> eint;
		SelfEnergy<Hamiltonian> eself;
		FlexEnergy<Hamiltonian> eflex;

	public:
		Energy (Hamiltonian& ham) : ham(ham), name("e"), eint(ham), eself(ham), eflex(ham)  {};


		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			return eint.measure(statespace,grid) + eself.measure(statespace,grid) + eflex.measure(statespace,grid);
		}
};

#endif
