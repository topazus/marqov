#ifndef ISING_H
#define ISING_H
#include <array>
#include <tuple>
#include <string>
#include <functional>


template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
	Ising_interaction()
	{
		this->J = -1;	// ferro
// 		this->J = +1; 	// antiferro
	}
	StateVector operator() (StateVector& phi) {return phi;};
};

template <class StateVector, class RNG>
class Ising_Initializer
{
public:
    Ising_Initializer()   {}
    Ising_Initializer(RNG&) {}
    StateVector newsv(const StateVector& svold) {StateVector retval(svold); retval[0]=-retval[0];return retval;};
};

class IsingMag
{
public:
    std::string name;
template <class StateSpace, class Grid>
double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			double mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += statespace[i][0];
			}

			return mag/double(N);
		}
		IsingMag() : name("m") {}
};

class IsingMagSquare
{
public:
    std::string name;
template <class StateSpace, class Grid>
int measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			float mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += pow(statespace[i][0],2);
			}
			return mag/double(N);
		}
		IsingMagSquare() : name("m2") {}
};

template <typename SpinType = int>
class Ising
{
public:
    constexpr static double beta = 1/2.26918;
    constexpr static int SymD = 1;
    typedef std::array<SpinType, SymD> StateVector;
    template <typename RNG>
    using MetroInitializer = Ising_Initializer<StateVector, RNG>;
    static constexpr uint Nalpha = 1;
    static constexpr uint Nbeta = 0;
    static constexpr uint Ngamma = 0;
    // requires pointers
    Interaction<StateVector>* interactions[Nalpha];
    OnSite<StateVector, int>* onsite[Nbeta];
    MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
    Ising()
    {
        interactions[0] = new Ising_interaction<StateVector>();
    }
    IsingMag       obs_m;
    IsingMagSquare obs_m2;
    auto getobs()
    {
        return std::make_tuple(obs_m, obs_m2);
    }

    StateVector createnewsv(const StateVector& osv) 
	 {
        StateVector retval(osv);
        retval[0] = -retval[0];
        return  retval;
	}
};
#endif
