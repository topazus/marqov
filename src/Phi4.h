#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include "vectorhelpers.h"
#include "hamiltonianparts.h"


// ------------------------------ OBSERVABLES ---------------------------

class Phi4Mag
{
     public:
          std::string name;
          template <class StateSpace, class Grid>
          double measure(const StateSpace& statespace, const Grid& grid)
          {
               constexpr static int SymD = 3;     // improve me
               const     int N = grid.size();

               std::vector<double> mag(SymD,0) ;

               for (int i=0; i<N; i++)
               {
                    for (int j=0; j<SymD; j++)
                    {
                         mag[j] += statespace[i][j];
                    }
               }

               double retval = 0;
               for (int j=0; j<SymD; j++)
               {
                    retval += mag[j]*mag[j];
               }
               return sqrt(retval)/double(N);
          }
          Phi4Mag() : name("m") {}
};


// ----------------------------------------------------------------------


template <class StateVector, class RNG>
class Phi4_Initializer
{
	public:
		// provide the spin dimension as a compile-time constant expression
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		Phi4_Initializer()   {}
		Phi4_Initializer(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector& osv) 
		{
            double amp = 0.5;
            double r = 2*rng.d() - 1.0;
            double oldlen = std::sqrt(dot(osv, osv));
            double newlen = oldlen + amp*r;
            auto newdir = rnddir<RNG, double, SymD>(rng);
            for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
                newdir[i] *= newlen;
			return newdir;
		};

	private:
		RNG& rng;
};

template <class StateVector>
class Phi4_interaction : public Interaction<StateVector>
{
	public:
		Phi4_interaction()
		{
	 		this->J = -1;	
		}
		StateVector operator() (const StateVector& phi) {return phi;};
};

template <class StateVector>
class Phi4_onsitesquare : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitesquare(double mass, double beta)
		{
	 		this->h = mass/beta;
		}
		double operator() (const StateVector& phi) {return dot(phi,phi);};
};

template <class StateVector>
class Phi4_onsitefour : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitefour(double lambda, double beta)
		{
	 		this->h = lambda/beta;
		}
		double operator() (const StateVector& phi) {return pow(dot(phi,phi)-1.0, 2);};
};



// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class Phi4
{
	public:
		double beta, lambda, mass;

		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;

		template <typename RNG>
		using MetroInitializer =  Phi4_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta  = 2;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta]; //Todo: External fields not yet supported
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		Phi4(double beta, double lambda, double mass) : lambda(lambda), mass(mass), beta(beta)
		{
			interactions[0] = new Phi4_interaction<StateVector>();
			onsite[0]       = new Phi4_onsitesquare<StateVector>(mass, beta);
			onsite[1]       = new Phi4_onsitefour<StateVector>(lambda, beta);
		}

		Phi4Mag obs_m;
		auto getobs() { return std::make_tuple(obs_m);}


		// --- Wolff ---

		template <class A>
		inline auto wolff_coupling(StateVector& sv1, StateVector& sv2, const A a) const 
		{
			return dot(sv1, a) * dot(sv2, a);
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a) const
		{
			const double dotp = dot(sv, a);
			for (int i=0; i<SymD; i++) sv[i] -= 2*dotp*a[i];
		}

		
};
#endif
