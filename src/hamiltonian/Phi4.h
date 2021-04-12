#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include <vector>
#include "../vectorhelpers.h"
#include "../hamparts.h"
#include "../obsparts.h"


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


		  	/* old (and wrong)

            	double amp = 0.05;
            	double r = rng.real(-1.0, 1.0);
            	double oldlen = std::sqrt(dot(osv, osv));
            	double newlen = oldlen + amp*r;
            	auto newdir = rnddir<RNG, double, SymD>(rng);
            	for(std::size_t i = 0; i < std::tuple_size<StateVector>::value; ++i)
            	    newdir[i] *= newlen;
		  	   return newdir;
		  	*/
		};

	private:
		RNG& rng;
};

template <class StateVector>
class Phi4_interaction
{
	public:
		Phi4_interaction(){}
		StateVector get (const StateVector& phi) {return phi;};
        static constexpr double J = -1;
};

template <class StateVector>
class Phi4_onsitesquare : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitesquare(double mass, double beta)
		{
	 		this->h = mass/beta;
		}
		double get (const StateVector& phi) {return dot(phi,phi);};
};

template <class StateVector>
class Phi4_onsitefour : public OnSite<StateVector, double> 
{
	public:
		Phi4_onsitefour(double lambda, double beta)
		{
	 		this->h = lambda/beta;
		}
		double get (const StateVector& phi) {return pow(dot(phi,phi)-1.0, 2);};
};



// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class Phi4
{
	public:
		double beta, lambda, mass;

		constexpr static int SymD = 3;
		const std::string name;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;

		template <typename RNG>
		using MetroInitializer =  Phi4_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		// requires pointers
        std::array<Phi4_interaction<StateVector>*, 1> interactions = {new Phi4_interaction<StateVector>()};
		std::vector<OnSite<StateVector, FPType>*> onsite; //Todo: External fields not yet supported
		std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;

		Phi4(double beta, double lambda, double mass) : beta(beta), lambda(lambda), mass(mass), name("Phi4"), obs_fx(0), obs_fy(1), obs_fz(2)
		{
            onsite.push_back(new Phi4_onsitesquare<StateVector>(mass, beta));
            onsite.push_back(new Phi4_onsitefour<StateVector>(lambda, beta));
		}

		Magnetization obs_m;
		MagFTComp obs_fx;
		MagFTComp obs_fy;
		MagFTComp obs_fz;
		auto getobs() { return std::make_tuple(obs_m, obs_fx, obs_fy, obs_fz);}


		// provide names for the parameters
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
		
		~Phi4() {delete onsite[1]; delete onsite[0]; delete interactions[0];}
};
#endif
