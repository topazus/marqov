#ifndef HEISENBERG_H
#define HEISENBERG_H
#include <array>
#include <cmath>
template <class RND, typename FPType, int SymD> 
auto rnddir(RND& rn) -> std::array<FPType, SymD>
{

	// Spherical coordinates according to:
	// https://sites.math.washington.edu/~morrow/335_12/sphericalCoords.pdf

    std::array<FPType, SymD> retval;
    retval[0] = 1;//beginning of recursion
    if(SymD > 1) //poor mans constexpr if
    {
        //set up that auxiliary data that might be erased later
        FPType angles[SymD-1];
        for(int i = 0; i < SymD - 1; ++i)
            angles[i] = 2.0*rn.d();
        //recursion for the polar(?) angles
        for(int j = 1; j < SymD - 1; ++j)
            retval[j] = retval[j-1] * angles[j-1]*(2.0 - angles[j-1]);
        for(int j = 0; j < SymD - 2; ++j)
            retval[j] = sqrt(retval[j]) * (1.0-angles[j]);
        retval[SymD - 2] = sqrt(retval[SymD - 2]);
        //fix up the azimuthal angle
        //note the dependency on retval[SymD-2] here. Hence the order is important
        retval[SymD - 1] = retval[SymD - 2] * cos(M_PI*angles[SymD - 2]);
        retval[SymD - 2] = retval[SymD - 2] * sin(M_PI*angles[SymD - 2]);
    }
    return retval;
}

template <class StateVector, class RNG>
class Heisenberg_Initializer
{
public:
    constexpr static int SymD = std::tuple_size<StateVector>::value;
    Heisenberg_Initializer()   {}
    Heisenberg_Initializer(RNG& rn) : rng(rn) {}
    StateVector newsv(StateVector&) {
         return rnddir<RNG, double, SymD>(rng);
    };
private:
    RNG& rng;
};

template <class StateVector>
class Heisenberg_interaction : public Interaction<StateVector> 
{
public:
	Heisenberg_interaction()
	{
// 		this->J = -1;	// ferro
		this->J = +1; 	// antiferro
	}
	StateVector operator() (StateVector& phi) {return phi;};
};


template <typename SpinType, typename MyFPType>
class Heisenberg
{
public:
    constexpr static int SymD = 3;
    typedef MyFPType FPType;
    typedef std::array<SpinType, SymD> StateVector;
    constexpr static MyFPType beta = 2;
    
    template <typename RNG>
    using MetroInitializer =  Heisenberg_Initializer<StateVector, RNG>;//C++11
    
    static constexpr uint Nalpha = 1;
    static constexpr uint Nbeta = 0;
    static constexpr uint Ngamma = 0;
    // requires pointers
    Interaction<StateVector>* interactions;
    OnSite<StateVector>* onsite[Nbeta];
    MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
    Heisenberg()
    {
        interactions = new Heisenberg_interaction<StateVector>();
    }

    StateVector createnewsv(const StateVector& osv) 
    {
	}
};
#endif
