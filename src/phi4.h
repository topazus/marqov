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
            double mag = 0.2;
            double r = 2*rng.d() - 1.0;
            double oldlen = std::sqrt(dot(osv, osv));
            double newlen = oldlen + mag*r;
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
	 		this->J = -1;	// -1 ferro, +1 antiferro
		}
		StateVector operator() (StateVector& phi) {return phi;};
};

template <class StateVector>
class Phi4_onsitesquare : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitesquare()
		{
	 		this->h = 1;	// -1 ferro, +1 antiferro
		}
		double operator() (StateVector& phi) {return dot(phi,phi);};
};

template <class StateVector>
class Phi4_onsitefour : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitefour()
		{
	 		this->h =1;	// -1 ferro, +1 antiferro
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
		constexpr static MyFPType beta = 1.0/0.00001;
		
		template <typename RNG>
		using MetroInitializer =  Phi4_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 1;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];//FIXME!!!External fields
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		Phi4()
		{
            interactions[0] = new Phi4_interaction<StateVector>();
            onsite[0] = new Phi4_onsitesquare<StateVector>();
            onsite[1] = new Phi4_onsitefour<StateVector>();
		}
		
		StateVector createnewsv(const StateVector& osv) 
		{
		}
};
#endif
