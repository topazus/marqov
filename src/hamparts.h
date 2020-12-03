#ifndef HAMILTONIANPARTS_H
#define HAMILTONIANPARTS_H
template <class StateVector>
class Interaction
{
	public:
	    double J;
//	    virtual StateVector operator() (const StateVector& phi_i) = 0;
	    virtual StateVector get(const StateVector& phi_i) = 0;
        virtual ~Interaction() {};
};

template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
	    CouplingType h;
//	    virtual CouplingType operator() (const StateVector& phi) = 0;
	    virtual CouplingType get(const StateVector& phi) = 0;
};

template <class StateSpace, class StateVector>
class MultiSite
{
	public: 
	    double k;
//	    virtual double operator() (const StateVector& sv, int svpos, StateSpace s) = 0;
	    virtual double get(const StateVector& sv, int svpos, StateSpace s) = 0;

		template <class Lattice>
          double diff (const int rsite,
                         const StateVector& svold,
                         const StateVector& svnew,
                         std::vector<int>& nbrs,
                         StateSpace& s,
                         Lattice& grid)
          {

			return 0;
}

};
#endif
