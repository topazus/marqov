#ifndef TERMCOLLECTION_H
#define TERMCOLLECTION_H

template <class StateVector, typename ConstantType = double>
class standard_interaction : public Interaction<StateVector>
{
	public:
//		const ConstantType& constant;
//		standard_interaction(const ConstantType &constant) : constant(constant) {}
//		StateVector get (const StateVector& phi) {return phi;};

		const ConstantType& J;
		standard_interaction(const ConstantType &J) : J(J) {}
		StateVector get (const StateVector& phi) {return phi;};
		
};


template <class StateVector, typename ConstantType = double>
class onsite_quadratic : public OnSite<StateVector, double> 
{
	public:
		onsite_quadratic(ConstantType constant)
		{
			this->h = constant;
		}
		double get (const StateVector& phi) {return dot(phi,phi);};
};



#endif
