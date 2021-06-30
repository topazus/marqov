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
#include "../../libmarqov/metropolishelpers.h"

/**
 * Magnetization
 * Euclidean norm of the sum of state vectors
 */
class Magnetization
{
	public:
		std::string name; ///< The name of the observable
        std::string desc; ///< A helpful description that will be used in the HDF5 output files.

        /** Construct a magnetization object
         * 
         * @param name How you want to call it.
         * @param description A helpful description for this observable.
         */
		Magnetization(std::string name, std::string description) : name(name), desc(description) {}
		Magnetization(std::string name) : name(name), desc("magnetization") {}
		Magnetization() : name("m"), desc("magnetization") {}

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


/**
 * MagFTComp
 * 
 */
class MagFTComp
{
	public:
		int dir;
		std::string name, desc;
		MagFTComp(int dir, std::string name, std::string description) : dir(dir), name(name), desc(description) {}
		MagFTComp(int dir, std::string name) : dir(dir), name(name), desc("Fourier Component of Magnetization") {}
		MagFTComp(int dir) : dir(dir), name("magft"+std::to_string(dir)), desc("Fourier Component of Magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const/*expr static*/ int SymD = statespace[0].size();
			const int N = grid.size();

			std::vector<std::complex<double>> magFTcomp(SymD,0);
			const std::complex<double> jj(0,1);

			for (int i=0; i<N; i++)
			{
				double x = grid.crds(i)[dir];
				for (int j=0; j<SymD; j++) magFTcomp[j] += double(statespace[i][j]) * std::exp(2.0*M_PI*x*jj);
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
		Hamiltonian& ham;

	public:
		InteractionEnergy (Hamiltonian& ham) : ham(ham), name("eint")  {};

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
			return ene/double(N);
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
		Hamiltonian& ham;

	public:
		SelfEnergy (Hamiltonian& ham) : ham(ham), name("eself")  {};

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
		Hamiltonian& ham;

	public:
		FlexEnergy (Hamiltonian& ham) : ham(ham), name("eflex")  {};

		std::string name;


		template <class StateSpace, class Grid>
		double measure_helper(const StateSpace& statespace, const Grid& grid, std::true_type)
		{
			double retval = 0;
			for (decltype(ham.multisite.size()) c = 0; c < ham.multisite.size(); c++)
			{
				double enepart = ham.multisite[c]->template energy<Grid>(statespace, grid, c);
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
			typedef typename MARQOV::HasOnsite<Hamiltonian>::type HasOns;
			return measure_helper(statespace,grid,HasOns());
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
		Hamiltonian& ham;
		InteractionEnergy<Hamiltonian> eint;
		SelfEnergy<Hamiltonian> eself;
		FlexEnergy<Hamiltonian> eflex;

	public:
		Energy (Hamiltonian& ham) : ham(ham), eint(ham), eself(ham), eflex(ham), name("e") {}

		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			return eint.measure(statespace,grid) + eself.measure(statespace,grid) + eflex.measure(statespace,grid);
		}
};
#endif
