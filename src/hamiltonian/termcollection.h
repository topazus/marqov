#ifndef TERMCOLLECTION_H
#define TERMCOLLECTION_H

//template <class StateVector, class ConstType>
template <class StateVector>
class standard_interaction : public Interaction<StateVector>
{
	public:
//		const double& constant;
//		standard_interaction(const double &constant) : constant(constant) {}
//		StateVector get (const StateVector& phi) {return phi;};

		const double& J;
		standard_interaction(const double &J) : J(J) {}
		StateVector get (const StateVector& phi) {return phi;};
		
};


template <class StateVector>
class onsite_quadratic : public OnSite<StateVector, double> 
{
	public:
		onsite_quadratic(double constant)
		{
			this->h = constant;
		}
		double get (const StateVector& phi) {return dot(phi,phi);};
};



#endif
