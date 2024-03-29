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

#ifndef OBSPARTS_H
#define OBSPARTS_H

#include <complex>
#include <string>
#include <vector>
#include <array>
#include "../../libmarqov/metropolishelpers.h"

/**
 * Magnetization
 * Euclidean norm of the sum of state vectors
 */
class Magnetization
{
	public:
	std::string name{"m"}; ///< The name of the observable
        std::string desc{"magnetization"}; ///< A helpful description that will be used in the HDF5 output files.

        /** Construct a magnetization object
         * 
         * @param name How you want to call it.
         * @param description A helpful description for this observable.
         */
		Magnetization(std::string name, std::string description) : name(name), desc(description) {}
		Magnetization(std::string name) : name(name), desc("magnetization") {}
		Magnetization() = default;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			std::vector<double> mag(SymD, 0);

			for (decltype(grid.size()) i = 0; i < grid.size(); i++) 
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j]; 
				}
			}

			double retval = 0;
			for (int j=0; j<SymD; j++) retval += mag[j]*mag[j];

			return sqrt(retval)/double(grid.size());
		}
};

/**
 * Vector Magnetization.
 * Sum of every component of the state vectors
 */
class VectorMagnetization
{
	public:
		std::string name, desc;

		VectorMagnetization(std::string name, std::string description) : name(name), desc(description) {}
		VectorMagnetization(std::string name) : name(name), desc("vector magnetization") {}
		VectorMagnetization() : name("mvec"), desc("vector magnetization") {}

		template <class StateSpace, class Grid>
		std::vector<double> measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			std::vector<double> mag(SymD, 0);

			for (int i=0; i<grid.size(); i++) 
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j]; 
				}
			}
			return mag;
		}
};

template <typename T, typename = void>
struct Is_tuple_like : std::false_type {};

template <typename T>
struct Is_tuple_like<T, MARQOV::detail::type_sink<std::integral_constant<int, std::tuple_size<T>::value> > > : std::true_type {};

/**
 * MagFTComp
 * 
 */
class MagFTComp
{
	public:
		int dir;
		std::string name;
		std::string desc{"Fourier Component of Magnetization"};
		MagFTComp(int dir, std::string name, std::string description) : dir(dir), name(name), desc(description) {}
		MagFTComp(int dir, std::string name) : dir(dir), name(name) {}
		MagFTComp(int dir) : dir(dir), name("magft"+std::to_string(dir)) {}

		template <class StateSpace, class Grid>
		auto measure(const StateSpace& statespace, const Grid& grid) -> typename std::enable_if<Is_tuple_like<typename StateSpace::value_type>::value, double>::type
		{
			constexpr static int SymD = std::tuple_size<typename StateSpace::value_type>::value;
			const int N = grid.size();

			std::array<std::complex<double>, SymD> magFTcomp{};

			for (int i=0; i<N; i++)
			{
				double x = grid.crds(i)[dir];
                auto expi = std::complex<double>(std::cos(2.0*M_PI*x), std::sin(2.0*M_PI*x));
				for (int j=0; j < SymD; j++) magFTcomp[j] += double(statespace[i][j]) * expi;
			}

			// normalize
			for (int j=0; j < SymD; j++) magFTcomp[j] /= double(N);

			// dot product of complex vector 
			double retval = 0;
			for (int j=0; j < SymD; j++) retval += std::norm(magFTcomp[j]);

			return retval;
		}
		
		template <class StateSpace, class Grid>
		auto measure(const StateSpace& statespace, const Grid& grid) -> typename std::enable_if<!Is_tuple_like<typename StateSpace::value_type>::value, double>::type
		{
			const/*expr static*/ int SymD = statespace[0].size();
			const int N = grid.size();

			std::vector<std::complex<double>> magFTcomp(SymD,0);

			for (int i=0; i<N; i++)
			{
				double x = grid.crds(i)[dir];
                auto expi = std::complex<double>(std::cos(2.0*M_PI*x), std::sin(2.0*M_PI*x));
				for (int j=0; j<SymD; j++) magFTcomp[j] += double(statespace[i][j]) * expi;
			}

			// normalize
			for (int j=0; j<SymD; j++) magFTcomp[j] /= double(N);

			// dot product of complex vector 
			double retval = 0;
			for (int j=0; j<SymD; j++) retval += std::norm(magFTcomp[j]);

			return retval;
		}
};

/**
 * Interaction Energy Observable
 * @tparam Hamiltonian We require the Hamiltonian to be able to calculate energies.
 */
template <class Hamiltonian>
class InteractionEnergy
{
	private:
		const Hamiltonian& ham;

	public:
		InteractionEnergy (const Hamiltonian& ham) : ham(ham), name("eint")  {};
		InteractionEnergy () = delete;

		std::string name;

		template <class StateSpace, class Grid>
		double measure_helper(const StateSpace& statespace, const Grid& grid, std::true_type)
		{
			const int N = grid.size();
			double ene = 0.0;

			// interaction part
			for (decltype(ham.interactions.size()) a = 0; a < ham.interactions.size(); a++)
			{

				double enepart = 0.0;
				for (int idx=0; idx<N; idx++)
				{
					auto nbrs = grid.nbrs(a, idx);
					auto self = ham.interactions[a]->get(statespace[idx]);

					for (std::size_t i = 0; i < nbrs.size(); ++i)
					{
						// index of the neighbour
						auto nbridx = nbrs[i];

						// configuration of the neighbour
						auto nbr = ham.interactions[a]->get(statespace[nbridx]);

						enepart += dot(self,nbr);

					}
				}

				ene += ham.interactions[a]->J * enepart;
			}
			return ene/double(2*N); // account for double counting
		}

		template <class StateSpace, class Grid>
		constexpr static double measure_helper(const StateSpace& statespace, const Grid& grid, std::false_type)
		{
			return 0;
		}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			typedef typename MARQOV::HasInteractions<Hamiltonian>::type HasInteract;
			return measure_helper(statespace, grid, HasInteract());
		}
};

/**
 * Self Energy Observable
 * @tparam Hamiltonian We require the Hamiltonian to be able to calculate energies.
 */
template <class Hamiltonian>
class SelfEnergy
{
	private:
		const Hamiltonian& ham;

	public:
		SelfEnergy (const Hamiltonian& ham) : ham(ham), name("eself")  {};
		SelfEnergy () = delete;

		std::string name;


		template <class StateSpace, class Grid>
		double measure_helper(const StateSpace& statespace, const Grid& grid, std::true_type hasonsite)
		{
			const int N = grid.size();
			double ene = 0.0;

			// self energy part
			for (decltype(ham.onsite.size()) b=0; b < ham.onsite.size(); b++)
			{
				double enepart = 0.0;

				for (int idx=0; idx<N; idx++)
				{
					enepart += ham.onsite[b]->get(statespace[idx]);
				}
				ene += ham.onsite[b]->h * enepart;
			}
			return ene/double(N);
		}

		template <class StateSpace, class Grid>
		constexpr static double measure_helper(const StateSpace& statespace, const Grid& grid, std::false_type hasonsite)
		{
			return 0;
		}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			typedef typename MARQOV::HasOnsite<Hamiltonian>::type HasOns;
			return measure_helper(statespace, grid, HasOns());
		}
};

/**
 * Flex Energy Observable
 * @tparam Hamiltonian We require the Hamiltonian to be able to calculate energies.
 */
template <class Hamiltonian>
class FlexEnergy
{
	private:
		const Hamiltonian& ham;

	public:
		FlexEnergy (const Hamiltonian& ham) : ham(ham), name("eflex")  {};
		FlexEnergy () = delete;

		std::string name;


		template <class StateSpace, class Grid>
		double measure_helper(const StateSpace& statespace, const Grid& grid, std::true_type)
		{
			double retval = 0;
			for (decltype(ham.multisite.size()) c = 0; c < ham.multisite.size(); c++)
			{
				double enepart = ham.multisite[c]->energy(statespace, grid, c);
				retval += ham.multisite[c]->k * enepart;
			}
			return retval;
		}

		template <class StateSpace, class Grid>
		constexpr static double measure_helper(const StateSpace& statespace, const Grid& grid, std::false_type)
		{
			return 0;
		}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			typedef typename MARQOV::HasFlexTerms<Hamiltonian>::type HasFlex;
			return measure_helper(statespace,grid,HasFlex());
		}
};

/**
 * Full Energy Observable
 * @tparam Hamiltonian We require the Hamiltonian to be able to calculate energies.
 */
template <class Hamiltonian>
class Energy
{
	private:
		const Hamiltonian& ham;
		InteractionEnergy<Hamiltonian> eint;
		SelfEnergy<Hamiltonian> eself;
		FlexEnergy<Hamiltonian> eflex;

	public:
		Energy (const Hamiltonian& ham) : ham(ham), eint(ham), eself(ham), eflex(ham), name("e") {}
		Energy() = delete;
		Energy(const Energy&) = default;
		~Energy() = default;

		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			return eint.measure(statespace,grid) + eself.measure(statespace,grid) + eflex.measure(statespace,grid);
		}
};
#endif
