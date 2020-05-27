#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <vector>
#include <type_traits>
#include <cmath>

//A helper to decide in the Metropolis code whether a lattice provides the getbond function
template<class L, class=void> struct has_bonds : std::false_type {};
template<class Lattice> struct has_bonds<Lattice, type_sink_t< decltype( std::declval<Lattice>().getbnds(std::declval<int>(), std::declval<int>(), std::declval<int>()) ) > > : std::true_type {};

template <class Lattice, class NbrType>
auto callbonds_helper(const Lattice& grid, int a, int rsite, int i, NbrType nbr, std::true_type)
{
    auto cpl = grid.getbnds(a, rsite, i);
    return mult(cpl, nbr);
}

template <class Lattice, class NbrType>
auto callbonds_helper(const Lattice& grid, int a, int rsite, int i, NbrType nbr, std::false_type)
{
    return nbr;
}

template <class Lattice, class NbrType>
auto callbonds(const Lattice& grid, int a, int rsite, int i, NbrType nbr)
{
    return callbonds_helper(grid, a, rsite, i, nbr, typename has_bonds<Lattice>::type());
}

template<class A, class B>
struct Promote_Array
{
    //Get Element types
    typedef decltype(dot(std::declval<A>(), std::declval<A>())) AElemType;
    typedef decltype(dot(std::declval<B>(), std::declval<B>())) BElemType;
    // get result type of addition
    typedef typename std::common_type<AElemType, BElemType>::type CommonType;
    // return A if the common type is A, else B
    typedef typename std::conditional<std::is_same<AElemType, CommonType>::value
    , A, B>::type CommonArray;
};

template <class Hamiltonian, class Lattice>
struct Metropolis
{
    template <class StateSpace, class M, class RNG>
    static int move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta, int rsite);
};

template <class Hamiltonian, class Lattice>
template <class StateSpace, class M, class RNG>
int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta, int rsite)
{
    typedef typename Hamiltonian::StateVector StateVector;
    	// old state vector at rsite
	StateVector& svold = statespace[rsite];
	// propose new configuration
	StateVector svnew = metro.newsv(svold);
	
	// interaction part
	double interactionenergydiff = 0;
	for(typename std::remove_cv<decltype(ham.Nalpha)>::type a = 0; a < ham.Nalpha; ++a)
	{
		auto nbrs = grid.getnbrs(a, rsite);
		typedef decltype(ham.interactions[a]->operator()(statespace[0])) InteractionType;
		typedef decltype(MARQOV::callbonds<Lattice>(grid, a, rsite, 0, ham.interactions[a]->operator()(statespace[0]))) BondType;
        
		typename MARQOV::Promote_Array<InteractionType, BondType>::CommonArray averagevector = {0};
		// sum over neighbours
		for (std::size_t i = 0; i < nbrs.size(); ++i)
		{
			auto idx = nbrs[i];
			auto nbr = ham.interactions[a]->operator()(statespace[idx]);
			averagevector = averagevector + MARQOV::callbonds<Lattice>(grid, a, rsite, i, nbr);
		}
		interactionenergydiff += ham.interactions[a]->J * (dot(svnew - svold, averagevector));
	}

    // onsite energy part
    double onsiteenergydiff = 0;
    for (typename std::remove_cv<decltype(ham.Nbeta)>::type b = 0; b < ham.Nbeta; ++b)
    {
       // compute the difference
       auto diff = ham.onsite[b]->operator()(svnew) - ham.onsite[b]->operator()(svold);
       // multiply the constant
       onsiteenergydiff += dot(ham.onsite[b]->h, diff);
    }
    // multi-site energy
    double multisiteenergyold = 0;
    double multisiteenergynew = 0;
    for (typename std::remove_cv<decltype(ham.Ngamma)>::type g = 0; g < ham.Ngamma; ++g)
    {
        multisiteenergynew += ham.multisite[g]->operator()(svnew, rsite, statespace);//FIXME: think about this...
        multisiteenergyold += ham.multisite[g]->operator()(svold, rsite, statespace);
        //forgot k_gamma
    }
    
    // sum up energy differences
    double dE 	= interactionenergydiff + onsiteenergydiff + (multisiteenergynew - multisiteenergyold);

	// improve me: what about models with discrete statevectors where the acceptance probability should be
	// looked up in tables? -> specialized Metropolis routine for this case??

    int retval = 0;
    if ( dE <= 0 )
    {
        svold = svnew;
        retval = 1;
    }
    else if (rng.d() < exp(-beta*dE))
    {
        svold = svnew;
        retval = 1;
    }
    return retval;
}

// Single Metropolis update step statevectors on a lattice
// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
template <class Grid, class Hamiltonian, template<class> class RefType>
inline int Marqov<Grid, Hamiltonian, RefType>::metropolisstep(int rsite)
{
    return Metropolis<Hamiltonian, Grid>::move(this->ham, this->grid, statespace, this->metro, this->rng, this->beta, rsite);
}


/*
// filtered Metropolis prototype ....

// takes a function which takes a StateVector and returns a reduced StateVector

template <class Grid, class Hamiltonian>
template <typename callable1, typename callable2>
inline int Marqov<Grid, Hamiltonian>::metropolisstep(int rsite, callable1 filter_ref, callable2 filter_cpy, int comp)
{
	// old state vector at rsite
	StateVector&     svold = statespace[rsite];
	redStateVector& rsvold = filter_ref(svold, comp);

	// propose new configuration
	redStateVector  rsvnew = metro.newsv(rsvold);
	
	// interaction part
	double interactionenergydiff = 0;
	for(int a = 0; a < ham.Nalpha; ++a)
	{
		// extract neighbours
		auto nbrs = grid.getnbrs(a, rsite);
		double averagevector = {0};
		
		// sum over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			// neighbour index
			auto idx = nbrs[i];
			// full neighbour
			auto nbr  = ham.interactions[a]->operator()(statespace[idx]);
			// reduced neighbour
			auto rnbr = filter_cpy(nbr, comp);
			
			// coupling in the embedded model
			double cpl = ham.metro_coupling(svold, nbr, comp);

			// compute weighted sum of neighbours
			averagevector = averagevector + mult(cpl,rnbr);
		}

		interactionenergydiff += ham.interactions[a]->J * (dot(rsvnew-rsvold, averagevector));
	}

    // (...)

    double dE 	= interactionenergydiff; // + ... + ...

    int retval = 0;
    if ( dE <= 0 )
    {
        rsvold = rsvnew;
        retval = 1;
    }
    else if (rng.d() < exp(-beta*dE))
    {
        rsvold = rsvnew;
        retval = 1;
    }

    return retval;
}

*/

#endif
