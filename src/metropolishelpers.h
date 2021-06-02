/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef METROPOLISHELPERS_H
#define METROPOLISHELPERS_H
#include <utility>
#include <type_traits>

namespace MARQOV
{
	/**
     * Promote_Array utility class.
     *
     * This function takes two types and tries to figure out, which one is wider.
     * @tparam A the first type.
     * @tparam B the other type.
     */
	template<class A, class B>
	struct Promote_Array
	{
		//Get Element types
		typedef decltype(dot(std::declval<A>(), std::declval<A>())) AElemType;
		typedef decltype(dot(std::declval<B>(), std::declval<B>())) BElemType;
		// get result type of addition
		typedef typename std::common_type<AElemType, BElemType>::type CommonType;
		// return A if the common type is A, else B
		typedef typename std::conditional<std::is_same<AElemType, CommonType>::value, A, B>::type CommonArray;
	};




	/** Helper functions to extract Hamiltonian on-site terms 
	*
	* @tparam Grid the type of the Lattice
	* @tparam Hamiltonian the type of the Hamiltonian
	* @param grid the lattice
	* @param ham the Hamiltonian
	* @param idx the site to be considered
	* @param std::false_type marks overload 
	*
	* @return std::vector<int>
	*/
	template <class Grid, class Hamiltonian>
	std::vector<int> get_terms(const Grid& grid, const Hamiltonian& ham, int idx, std::false_type)
	{
		return arange(0, ham.onsite.size());
	}

	/** Helper functions to extract Hamiltonian on-site terms 
	*
	* @tparam Grid the type of the Lattice
	* @tparam Hamiltonian the type of the Hamiltonian
	* @param grid the lattice
	* @param ham the Hamiltonian
	* @param idx the site to be considered
	* @param std::true_type marks overload 
	*
	* @return std::vector<int>
	*/
	template <class Grid, class Hamiltonian>
	std::vector<int> get_terms(const Grid& grid, const Hamiltonian& ham, int idx, std::true_type)
	{
		return grid.termselector(idx);
	}




 
    /** A helper to check whether the lattice provides a bnds method.
     * @tparam L The lattice that we are querying.
     */
	template<class Lattice, class=void>
	struct has_bnds : std::false_type {};

	template<class Lattice>
	struct has_bnds<Lattice, MARQOV::detail::type_sink_t<
		decltype( std::declval<Lattice>().bnds(std::declval<int>(), std::declval<int>()) )
		>> : std::true_type {};




	/** A helper to check whether the lattice provides a termselector method.
     * @tparam Grid the Grid which we check.
     */
	template<class Grid , class = void>
	struct has_trms : std::false_type {};

	template<class Grid>
	struct has_trms<Grid, MARQOV::detail::type_sink_t<
		decltype(std::declval<Grid>().termselector(std::declval<int>()))
		>> : std::true_type {};




	/** A helper to check whether the lattice provides a nbrs function.
     * @tparam Grid the Grid which we check.
     */
	template<class Grid, class = void>
	struct has_nbrs : std::false_type {};

	template<class Grid>
	struct has_nbrs<Grid, MARQOV::detail::type_sink_t<
		decltype( std::declval<Grid>().nbrs(std::declval<int>(), std::declval<int>()) )
		>> : std::true_type {};



	/** A helper to check whether the lattice provides a flexnbrs function.
     * @tparam Grid the Grid which we check.
     */
	template<class Grid, class = void>
	struct has_flexnbrs : std::false_type {};

	template<class Grid>
	struct has_flexnbrs<Grid, MARQOV::detail::type_sink_t<
		decltype( std::declval<Grid>().flexnbrs(std::declval<int>(), std::declval<int>()) )
		>> : std::true_type {};





	template <class Grid>
	auto getnbrs_helper(const Grid& grid, int fam, int idx, std::false_type)
	{
		cout << "nbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}

	template <class Grid>
	auto getnbrs_helper(const Grid& grid, int fam, int idx, std::true_type)
	{
		return grid.nbrs(fam,idx);
	}

	template <class Grid>
	auto getnbrs(const Grid& grid, int fam, int idx)
	{
		return getnbrs_helper<Grid>(grid, fam, idx, has_nbrs<Grid>{});
	}






	template <class Grid>
	auto getflexnbrs_helper(const Grid& grid, int fam, int idx, std::false_type)
	{
		cout << "flexnbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}

	template <class Grid>
	auto getflexnbrs_helper(const Grid& grid, int fam, int idx, std::true_type)
	{
		return grid.flexnbrs(fam,idx);
	}

	/** A helper to detect if the lattice has a getflexnbrs function.
     *
     * @tparam Grid the type of the lattice.
     */
	template <class Grid>
	auto getflexnbrs(const Grid& grid, int fam, int idx)
	{
		return getflexnbrs_helper<Grid>(grid, fam, idx, has_flexnbrs<Grid>{});
	}



	/** has_bnds utility struct
     *
     * This struct checks whether the lattice provides a getbnds method.
     * @tparam Grid the lattice that we query.
     */
//	template<class Grid, class = void>
//	struct has_bnds : std::false_type {};
//
//	template<class Grid>
//	struct has_bnds<Grid, MARQOV::detail::type_sink_t<
//		decltype( std::declval<Grid>().bnds(std::declval<int>(), std::declval<int>()) )


	template <class Grid>
	auto getbnds_helper(const Grid& grid, int fam, int idx, std::false_type)
	{
		cout << "error" << endl;
		exit(0);
		return std::vector<int>{1,1,1,1};
	}

	template <class Grid>
	auto getbnds_helper(const Grid& grid, int fam, int idx, std::true_type)
	{
		return grid.bnds(fam,idx);
	}

	/** A helper to detect if the lattice has a getbnds function.
     *
     * @tparam Grid the type of the lattice.
     */
	template <class Grid>
	auto getbnds(const Grid& grid, int fam, int idx)
	{
		return getbnds_helper<Grid>(grid, fam, idx, has_bnds<Grid>{});
	}


	
	
	
	
	
	
/**
 * Helper function for the Metropolis algorithm, sums over local neighbourhood of spin
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has function grid.nbrs
 * @param std::true_type lattice has function grid.bnds
 *
 * @return sum over neighbours after applying the interaction term and weighted by respective bond strengths
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::true_type, 
											  std::true_type) 
{
	typedef typename Hamiltonian::StateVector StateVector;

	// gather neighbours and bonds
	const auto nbrs = getnbrs<Lattice>(grid, a, rsite);
	const auto bnds = getbnds<Lattice>(grid, a, rsite);

	StateVector neighbourhood = {0};

	// sum over neighbours
	for (std::size_t i=0; i<nbrs.size(); ++i)
	{
		// index of the neighbour
		auto idx = nbrs[i];
		                
		// configuration of the neighbour after applying the interaction
		auto nbr = ham.interactions[a]->get(statespace[idx]);
		                
		// sum up contributions from neighbourbood
		neighbourhood = neighbourhood + mult(bnds[i], nbr);
	}

	return neighbourhood;
}


/**
 * Helper function for the Metropolis algorithm, sums over local neighbourhood of spin
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has function grid.nbrs
 * @param std::false_type lattice has no function grid.bnds
 *
 * @return sum over neighbours after applying the interaction term
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::true_type, 
											  std::false_type) 
{
	typedef typename Hamiltonian::StateVector StateVector;

	// gather neighbours
	const auto nbrs = getnbrs<Lattice>(grid, a, rsite);

	StateVector neighbourhood = {0};

	// sum over neighbours
	for (std::size_t i=0; i<nbrs.size(); ++i)
	{
		// index of the neighbour
		auto idx = nbrs[i];
		                
		// configuration of the neighbour after applying the interaction
		auto nbr = ham.interactions[a]->get(statespace[idx]);
		                
		// sum up contributions from neighbourbood
		neighbourhood = neighbourhood + nbr;
	}

	return neighbourhood;
}


/**
 * Helper function for the Metropolis algorithm, sums over local neighbourhood of spin.
 * this overloading covers the case when the lattice has no function "nbrs".
 * in this case the general Metropolis algorithm cannot be used hence the user is given an error message
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has no function grid.nbrs
 * @param std::false_type lattice has no function grid.bnds
 *
 * @return does not apply 
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::false_type, 
											  std::false_type) 
{
	cout << "[MARQOV] Error: The lattice does not provide the following function: nbrs" << endl;
	cout << "In order to use the general Metropolis algorithm, this function must be implemented" << endl;
	cout << "Alternatively, you can write you own specialization of the Metropolis algorithm" << endl;
	exit(0);
}







    /**
     * Is_Container utility struct.
     *
     * Helpers to determine if the interactions are container-like.
     * By that we mean whether it has a .size() method and an array access
     * operator.
     * See metropolis.h for an example.
     * @see metropolis.h
     * @tparam Cont The Container that we query.
     */
    template <class Cont, class = void, class = void>
    struct Is_Container
    {
        static constexpr bool value = false;
    };

    template <class Cont>
    struct Is_Container<Cont,
    MARQOV::detail::type_sink_t<decltype(std::declval<Cont>().size())>,
    MARQOV::detail::type_sink_t<decltype(std::declval<Cont>().operator[](std::declval<std::size_t>()))>
    >
    {
        static constexpr bool value = true;
    };






	// check of the existence of class of member with name "interactions" at compile time
	// can be of any type!
	template <typename T, typename = int>
	struct HasInteractions : std::false_type { };

	template <typename T>
	struct HasInteractions <T, decltype((void) T::interactions, 0)> : std::true_type { };

	template <class Grid, class Hamiltonian, class StateSpace, class StateVector>
	double compute_interactionenergydiff(Grid& grid, Hamiltonian& ham, StateSpace& statespace, StateVector& svnew, StateVector& svold, int rsite, std::true_type)
	{
		static_assert(MARQOV::Is_Container<decltype(std::declval<Hamiltonian>().interactions)>::value,
			"[MARQOV::Metropolis] COMPILATION FAILED: multisite terms are not a container.");
//		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename std::remove_cv<decltype(ham.interactions.size())>::type InteractionSizeType;
		typedef typename has_nbrs<Grid>::type HasNbrs;
		typedef typename has_bnds<Grid>::type HasBnds;

		double interactionenergydiff = 0;
		for (InteractionSizeType a=0; a<ham.interactions.size(); a++)
		{
			// sum up neighbourhood contributions
			auto nbrhood = nbrhoodloop<Grid,Hamiltonian,StateSpace>(grid, ham, statespace, a, rsite, HasNbrs(), HasBnds());

			// compute interaction energy difference
			interactionenergydiff += ham.interactions[a]->J * (dot(svnew-svold, nbrhood));
		}
		return interactionenergydiff;
	}

	template <class Grid, class Hamiltonian, class StateSpace, class StateVector>
	double compute_interactionenergydiff(Grid& grid, Hamiltonian& ham, StateSpace& statespace, StateVector& svnew, StateVector& svold, int rsite, std::false_type)
	{
		return 0;
	}









	// check of the existence of class of member with name "onsite" at compile time
	// can be of any type!
	template <typename T, typename = int>
	struct HasOnsite : std::false_type { };

	template <typename T>
	struct HasOnsite <T, decltype((void) T::onsite, 0)> : std::true_type { };

	template <class Grid, class Hamiltonian, class StateVector>
	double compute_onsiteenergydiff(Grid& grid, Hamiltonian& ham, StateVector& svnew, StateVector& svold, int rsite, std::true_type)
	{
		static_assert(MARQOV::Is_Container<decltype(std::declval<Hamiltonian>().onsite)>::value,
			"[MARQOV::Metropolis] COMPILATION FAILED: onsite terms are not a container.");
//		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename std::remove_cv<decltype(ham.onsite.size())>::type OnSiteSizeType;


		// find Hamiltonian terms associated to "rsite"
		typedef typename has_trms<Grid>::type HasTrms;

		double onsiteenergydiff = 0;


		auto terms = get_terms<Grid,Hamiltonian>(grid, ham, rsite, HasTrms());

		for (OnSiteSizeType b=0; b<terms.size(); b++)
		{
			// select on-site term
			const int tidx = terms[b];

			// compute the difference
			auto diff = ham.onsite[tidx]->get(svnew) - ham.onsite[tidx]->get(svold);

			// multiply the constant
			onsiteenergydiff += dot(ham.onsite[tidx]->h, diff);
		}

		return onsiteenergydiff;
	}


	template <class Grid, class Hamiltonian, class StateVector>
	double compute_onsiteenergydiff(Grid& grid, Hamiltonian& ham, StateVector& svnew, StateVector& svold, int rsite, std::false_type)
	{
		return 0;
	}






	// check of the existence of class of member with name "multisite" at compile time
	// can be of any type!
	template <typename T, typename = int>
	struct HasFlexTerms : std::false_type { };

	template <typename T>
	struct HasFlexTerms <T, decltype((void) T::multisite, 0)> : std::true_type { };

	template <class Grid, class Hamiltonian, class StateSpace, class StateVector>
	double compute_flexenergydiff(Grid& grid, Hamiltonian& ham, StateSpace& statespace, StateVector& svnew, StateVector& svold, int rsite, std::true_type)
	{
		static_assert(MARQOV::Is_Container<decltype(std::declval<Hamiltonian>().multisite)>::value,
			"[MARQOV::Metropolis] COMPILATION FAILED: multisite terms are not a container.");
//		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename std::remove_cv<decltype(ham.multisite.size())>::type FlexSizeType;

		double flexenergydiff = 0;
		for (FlexSizeType c=0; c<ham.multisite.size(); c++)
		{
			auto nbrs = getflexnbrs<Grid>(grid, c, rsite);
			auto diff = ham.multisite[c]->diff(rsite, svold, svnew, nbrs, statespace);
			flexenergydiff += dot(ham.multisite[c]->k, diff); 
		}
		return flexenergydiff;
	}

	template <class Grid, class Hamiltonian, class StateSpace, class StateVector>
	double compute_flexenergydiff(Grid& grid, Hamiltonian& ham, StateSpace& statespace, StateVector& svnew, StateVector& svold, int rsite, std::false_type)
	{
		return 0;
	}

};

#endif
