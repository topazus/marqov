#ifndef REGULAR_HYPERCUBIC_H
#define REGULAR_HYPERCUBIC_H

#include "util/points.h" 
#include "util/distance.h"
#include "util/regular_lattice.h"
#include "../libmarqov/cachecontainer.h"

/** @class RegularHypercubic
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
		RegularHypercubic(RegularHypercubic&&) noexcept = default;
		RegularHypercubic& operator=(const RegularHypercubic&) = delete;
        RegularHypercubic& operator=(const RegularHypercubic&&) = delete;
        RegularHypercubic() = delete;
		// implement nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		/** Get the real space coordinates of a site.
		*
		* @param i the index of the site.
		* @return the real space coordinates the point.
		*/
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		// implement flexnbrs
		std::vector<int> flexnbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		/** Return the number of sites in this lattice.
         * 
         * @return the number of sites in this lattice.
         */
		std::size_t size() const noexcept {return npoints;}
};

/** This function provides the necessary overload for MARQOV::Core
 *  such that we can dump useful info into the HDF5 file.
 * 
 * @param h5loc the group of the lattice.
 * @param l the lattice that we are dumping info about.
 */
void writelat(H5::Group& h5loc, const RegularHypercubic& l)
{
	dumpscalartoH5(h5loc, "name", std::string("RegularHypercubic"));
	dumpscalartoH5(h5loc, "dim", l.dim);
	dumpscalartoH5(h5loc, "npoints", l.npoints);
}
#endif
