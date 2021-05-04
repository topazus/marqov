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



#endif
