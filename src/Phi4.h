#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include "vectorhelpers.h"



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
		StateVector newsv(StateVector& osv) 
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
		StateVector operator() (StateVector& phi) {return phi;};
};

template <class StateVector>
class Phi4_onsitesquare : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitesquare(double beta)
		{
	 		this->h = 1;
		}
		double operator() (StateVector& phi) {return dot(phi,phi);};
};

template <class StateVector>
class Phi4_onsitefour : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitefour(double lambda, double beta)
		{
	 		this->h = lambda/beta;
		}
		double operator() (StateVector& phi) {return pow(dot(phi,phi)-1.0, 2);};
};

template <typename SpinType, typename MyFPType>
class Phi4
{
	public:
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;

		constexpr static MyFPType beta = 1.0/0.5;
		constexpr static MyFPType lambda = 1.0;
		
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

		Phi4()
		{
            interactions[0] = new Phi4_interaction<StateVector>();
            onsite[0] = new Phi4_onsitesquare<StateVector>(beta);
            onsite[1] = new Phi4_onsitefour<StateVector>(1.0,beta);
		}
		
		StateVector createnewsv(const StateVector& osv) 
		{
		}
};
#endif
