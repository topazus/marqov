#ifndef OBSPARTS_H
#define OBSPARTS_H


class ScalarMagnetization
{
	public:
	std::string name, desc;

	ScalarMagnetization(std::string name, std::string description) : name(name), desc(description) {}
	ScalarMagnetization(std::string name) : name(name), desc("scalar magnetization") {}
	ScalarMagnetization() : name("m"), desc("scalar magnetization") {}

	template <class StateSpace, class Grid>
	double measure(const StateSpace& statespace, const Grid& grid)
	{
		double retval = 0.0;
		for (int i=0; i<grid.size(); i++) retval += statespace[i][0];
		return std::abs(retval)/double(grid.size());
	}
};




template <class Hamiltonian>
class Energy
{
	private:
		Hamiltonian& ham;

	public:
		Energy (Hamiltonian& ham) : ham(ham), name("e")  {};


		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			double ene = 0.0;

			// interaction part
			for (int a=0; a<ham.Nalpha; a++)
			{

				double enepart = 0.0;
				for (int idx=0; idx<N; idx++)
				{
					auto nbrs = grid.getnbrs(a, idx);
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

			// self energy part
			for (int b=0; b<ham.Nbeta; b++)
			{
				double enepart = 0.0;

				for (int idx=0; idx<N; idx++)
				{
					enepart += ham.onsite[b]->get(statespace[idx]);
				}
				ene += ham.onsite[b]->h * enepart;
			}

			// multi-site terms
			for (int c=0; c<ham.Ngamma; c++)
			{
//				enepart = ham.multisite[c]->energy();
//				ene += ham.onsite[c]->k * enepart;
			}



			return ene/double(N);
		}
};

#endif
