#ifndef TERMCOLLECTION_H
#define TERMCOLLECTION_H

#include "../hamparts.h"

template <class StateVector, typename CouplingType = double>
class standard_interaction : public Interaction<StateVector>
{
	public:
		const CouplingType& J;
		standard_interaction(const CouplingType &J) : J(J) {}
		StateVector get (const StateVector& phi) {return phi;};

		// todo: rename J to something more generic, like "constant"
		
};


template <class StateVector, typename CouplingType = double>
class onsite_quadratic : public OnSite<StateVector, CouplingType> 
{
	public:
		onsite_quadratic(CouplingType constant)
		{
			this->h = constant;
		}
		double get (const StateVector& phi) {return dot(phi,phi);};
};



#endif
