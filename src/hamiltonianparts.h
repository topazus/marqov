#ifndef HAMILTONIANPARTS_H
#define HAMILTONIANPARTS_H
template <class StateVector>
class Interaction
{
	public:
	    double J;
	    virtual StateVector operator() (const StateVector& phi_i) = 0;
};

template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
	    CouplingType h;
	    virtual CouplingType operator() (const StateVector& phi) = 0;
	private:
};

template <class StateSpace, class StateVector>
class MultiSite
{
	public: 
	    double k;
	    virtual double operator() (const StateVector& sv, int svpos, StateSpace s) = 0;
	private:
};
#endif
