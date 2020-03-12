#ifndef ISING_H
#define ISING_H
#include <array>
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
    Interaction<StateVector>* interactions;
    OnSite<StateVector>* onsite[Nbeta];
    MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
    Ising()
    {
        interactions = new Ising_interaction<StateVector>();
    }

    StateVector createnewsv(const StateVector& osv) 
	 {
        StateVector retval(osv);
        retval[0] = -retval[0];
        return  retval;
	}
};
#endif
