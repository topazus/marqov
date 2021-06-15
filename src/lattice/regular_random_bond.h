#ifndef REGULAR_RANDOM_BOND_H
#define REGULAR_RANDOM_BOND_H

#include "util/disorderbase.h"
#include "util/prob_density.h"



/**
 * A lattice with bond disorder
 * @tparam PDFType The type of disorder distribution to use.
 */
template <class PDFType>
class RegularRandomBond :  public DisorderType
{
	private:
		RegularLattice lattice;
		PDFType PDF;
	
	public:
		int len, dim, npoints;

		using bond_type = decltype(PDF.draw());
		std::vector<std::vector<bond_type>> bonds;
		
		RegularRandomBond(int len, int dim) : lattice(len,dim), len(len), dim(dim), npoints(pow(len,dim))
		{
			// construct bonds
			bonds.resize(lattice.size());
			for (decltype(lattice.size()) i=0; i<lattice.size(); i++)
			{
				for (decltype(lattice[i].size()) j = 0; j < lattice[i].size(); j++) // why does lattice[i].size even work?
				{
					bonds[i].push_back(PDF.draw());
				}
			}

			// "symmetrize" bonds
			for (decltype(lattice.size()) i = 0; i < lattice.size(); i++)
			{
				auto lnbrs = lattice.nbrs(1,i);

				for (decltype(lattice[i].size()) j = 0; j < lattice[i].size(); j++)
				{
					// find i in bonds[lnbr] and replace its value by bonds[i][j]
					auto lnbr = lnbrs[j];
					auto nbrs_temp = lattice.nbrs(1, lnbr);

					auto it  = std::find(nbrs_temp.begin(), nbrs_temp.end(), i);
					auto idx = std::distance(nbrs_temp.begin(), it);

					bonds[lnbr][idx] = bonds[i][j];
				}
			}
		}

		// override nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement bnds
		std::vector<bond_type> bnds(const int alpha, const int i) const
		{
			return bonds[i];
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		std::size_t size() const {return npoints;}

};

#endif
