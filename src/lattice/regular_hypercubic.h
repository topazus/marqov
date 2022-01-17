#ifndef REGULAR_HYPERCUBIC_H
#define REGULAR_HYPERCUBIC_H

#include <array>
#include "util/points.h" 
#include "util/distance.h"
#include "util/regular_lattice.h"
#include "../libmarqov/cachecontainer.h"

/** @class RegularHypercubic
 * The Regular Hypercubic lattice class.
 * This class provides routines for encapsulating all neighbour
 * and coordinate relations for n-dimensional hypercubic lattices.
 */
template <int dimT>
class FixedDimRegularHypercubic
{
public:
    int len{0};
    constexpr static int dim{dimT};
    int npoints{0};
    constexpr FixedDimRegularHypercubic(int len) : len(len), npoints(1)
    {
        for(int i = 0; i < dimT; ++i)
        {
            pows[i] = npoints;
            npoints *= len;
        }
    }
		// implement nbrs. FIXME: alpha
		constexpr std::array<int,2*dimT> nbrs(const int alpha, const int i) const
		{
			std::array<int, 2*dimT> temp{0};
			//calculate neighbours for site i
			for(int j = 0; j < dim; ++j)
			{
				int pl = pows[j]*len;

				//positive additions
				int c = i + pows[j];

				//test positive additions for PBCs
				if (c >= (c/pl)*pl)
				{
					temp[2*j] = (i/pl)*pl + c % pl;
				}
				else
				{
					temp[2*j] = c;
				}
				
				//negative additions
				c = i - pows[j];

				//test negative additions for PBCs
				if (c < (i/pl)*pl)
				{
					temp[2*j+1] = (i/pl)*pl + (c + pl) % pl;
				}
				else
				{
					temp[2*j+1] = c;
				}
			}
			return temp;
		}

		/** Get the real space coordinates of a site.
		*
		* @param k the index of the site.
		* @return the real space coordinates the point.
		*/
		constexpr std::array<double,dimT> crds(int k) const
		{
            std::array<double, dimT> retval{0};
            double mya = static_cast<double>(1)/len;
            for (int i=0; i < dimT; i++)
            {
                int index = (k/pows[i])%len;
                k -= index * pows[i];
                retval[i] = index * mya;
                retval[i] += 0.5 * mya;
            }
            return retval;
		}

		// implement flexnbrs
		constexpr std::array<int,2*dimT> flexnbrs(const int alpha, const int i) const
		{
            return this->nbrs(alpha, i);
		}

		/** Return the number of sites in this lattice.
         * 
         * @return the number of sites in this lattice.
         */
		constexpr std::size_t size() const noexcept {return npoints;}
private:
    std::array<int, dim> pows{0};
};

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
