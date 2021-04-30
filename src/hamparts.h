#ifndef HAMILTONIANPARTS_H
#define HAMILTONIANPARTS_H

template <class StateVector>
class Interaction
{
	public:
		double J;
		virtual StateVector get(const StateVector& phi_i) = 0;
		virtual ~Interaction() {};
};

template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
		CouplingType h;
		virtual CouplingType get(const StateVector& phi) = 0;
		virtual ~Onsite(){};
};


template <class StateSpace, class StateVector>
class FlexTerm
{
	public:
		double k;
		virtual double get(const StateVector& sv, int svpos, StateSpace s) = 0;
		virtual ~FlexTerm() {};
		template <class Lattice>
		double diff (const int rsite,
					const StateVector& svold,
					const StateVector& svnew,
					std::vector<int>& nbrs,
					StateSpace& s,
					Lattice& grid) {return 0;}

};


#endif
