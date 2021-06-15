#ifndef REGULAR_HYPERCUBIC_H
#define REGULAR_HYPERCUBIC_H

#include "util/points.h" 
#include "util/distance.h"
#include "util/regular_lattice.h"

/**
 * The Regular Hypercubic lattice class.
 * This class provides routines for encapsulating all neighbour
 * and coordinate relations for n-dimensional hypercubic lattices.
 */
class RegularHypercubic
{
	private:
		RegularLattice lattice;

	public:
		int len, dim, npoints;

		RegularHypercubic(int len, int dim) : lattice(len, dim), len(len), dim(dim), npoints(pow(len, dim)) {}
		RegularHypercubic(const RegularHypercubic&) = default;
		RegularHypercubic(RegularHypercubic&&) = default;
		const RegularHypercubic& operator=(const RegularHypercubic&) = delete;

		// implement nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		// implement flexnbrs
		std::vector<int> flexnbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		std::size_t size() const {return npoints;}
};

#endif
