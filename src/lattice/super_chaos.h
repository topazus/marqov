#include <vector>
/**
 * Super Chaos Lattice
 * @tparam PointCloud the point cloud that we use
 * @tparam bond_type 
 */
template <class PointCloud, typename bond_type>
class SuperChaos : public DisorderType
{
private:
	int symD;
public:
	std::vector<std::vector<std::vector<bond_type>>> bonds;

	SuperChaos(const PointCloud& cloud) : DisorderType(cloud.size())
	{
		// prepare neighbour array
		const int npoints = cloud.size;
		this->neighbours.resize(npoints);
		this->bonds.resize(npoints);
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> disreal(0.0, 1.0);
        std::uniform_int_distribution<> disint(0, npoints-1);

		// draw bonds
		for (int i=0; i<2*npoints; i++) // make me variable
		{
			const int j = disint(gen);
			const int k = disint(gen);

			auto jcoordinates = cloud.crds(j);
			auto kcoordinates = cloud.crds(k);

			if (i!=k && distancePBSQ_nD(jcoordinates,kcoordinates) < 0.1) // make me 1/L dependent
			// improve me: make sure bond does not already exist!
			{
				this->neighbours[j].push_back(k);
				this->neighbours[k].push_back(j);
			}
		}

		/* under construction
		
		// compute weights
		for (int k=0; k<npoints; k++)
		{
			std::vector<bond_type> temp;
			for (int m=0; m<this->nbrs[k].size(); m++)
			{
				bond_type subtemp;
				for (int n=0; n<sizeof(bond_type); n++)
				{
					bond_type[n] =disreal(gen);
				}
				temp.push_back(subtemp);
			}
			bonds[k] = temp;
		}
				*/
	}

	// override bnds
	std::vector<bond_type> bnds(const int i)
	{
		return bonds[i];
	}
};
