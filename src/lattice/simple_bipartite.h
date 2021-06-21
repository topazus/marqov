#ifndef SIMPLE_BIPARTITE_H
#define SIMPLE_BIPARTITE_H

/**
 * Simple Bipartite Lattice
 */
class SimpleBipartite
{
	private:
		RegularLattice lattice;

	public:
		int len, dim, npoints;

		SimpleBipartite(int len, int dim) : lattice(len,dim), len(len), dim(dim), npoints(pow(len,dim))
		{
			if (len%2 != 0) cout << "ERROR: linear lattice size must be even!" << endl;
		}


		inline int identify(int i) const // is this correct?
		{
			auto index = IndexOf(i, dim, len);
			
			int quersumme = 0;
			for (decltype(index.size()) j = 0; j < index.size(); j++) quersumme += index[j];

			if (quersumme%2 == 0) return 0;
			else return 1;
		}


		inline std::vector<int> termselector(int rsite) const
		{
			return {this->identify(rsite)};
		}


		// override nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		std::size_t size() const {return npoints;}
};

#endif