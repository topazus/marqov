#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include <vector>
#include "../vectorhelpers.h"
#include "../hamparts.h"
#include "../obsparts.h"
#include "termcollection.h" 

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

			/* Francesco version (component-wise)
			const int comp =	rng.integer(3);
			double amp = 0.5;
			double r = rng.real(-1.0, 1.0);
			double oldval = osv[comp];
			double newval = oldval + amp*r;
			
			auto nsv = osv;
			nsv[comp] = newval;
			return nsv;
			*/



			double amp = 0.5;
			auto newdir = rnddir<RNG, double, SymD>(rng);
			auto nsv = osv + mult(amp,newdir);
			return nsv;
		};

	private:
		RNG& rng;
};




template <class StateVector>
class onsite_fourth_minus_one : public OnSite<StateVector, double> 
{
	public:
		onsite_fourth_minus_one(double constant)
		{
	 		this->h = constant;
		}
		double get (const StateVector& phi) {return pow(dot(phi,phi)-1.0, 2);};
};



// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class Phi4
{
	public:
		//  ----  Parameters  ----

		double beta, lambda, mass;
		const double J = -1;
		constexpr static int SymD = 3;
		const std::string name = "Hello, my name is Phi";



		//  ---- Definitions  -----

		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer =  Phi4_Initializer<StateVector, RNG>; 



		//  ----  Hamiltonian terms  ----

		standard_interaction<StateVector> INT;
		onsite_quadratic<StateVector> ONS1;
		onsite_fourth_minus_one<StateVector> ONS2;

		std::vector<Interaction<StateVector>*>                 interactions;
		std::vector<OnSite<StateVector, FPType>*>              onsite; 
		std::array <FlexTerm<StateVector*,  StateVector>*, 0>  multisite;

		Phi4(double beta, double lambda, double mass) : beta(beta), 
												lambda(lambda), 
												mass(mass), 
												name("Phi4"), 
												obs_fx(0), 
												obs_fy(1), 
												obs_fz(2), 
												INT(J),
												ONS1(mass/beta),
												ONS2(lambda/beta)
		{
			interactions.push_back(&INT);
			onsite.push_back(&ONS1);
			onsite.push_back(&ONS2);
		}



		//  ----  Observables  ----

		Magnetization obs_m;
		MagFTComp obs_fx;
		MagFTComp obs_fy;
		MagFTComp obs_fz;
        decltype(std::make_tuple(obs_m, obs_fx, obs_fy, obs_fz)) observables = {std::make_tuple(obs_m, obs_fx, obs_fy, obs_fz)};



		// Provide names for the parameters

		std::string paramname(int i)
		{
			std::string name;
			switch (i)
			{
				case (0): name = "beta"; break;
				case (1): name = "lambda"; break;
				case (2): name = "mass"; break;
				default: break;
			}
			return name;
		}


		//  ----  Wolff  ----

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
