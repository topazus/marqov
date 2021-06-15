
/**
 * Erdos-Renyj Graph
 */
class ErdosRenyi : public DisorderType
{
	private:
		int p;
	public:
		ErdosRenyi(int npoints, double p) : DisorderType(npoints), p(p)
		{
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, 1.0);
			// prepare neighbour array
			this->neighbours.resize(npoints);
	
			// the actual implementation
			// of the Erdos-Renyj Graph
			for (int i=0; i<npoints; i++)
			{
				for (int j=i+1; j<npoints; j++)
				{
					if (dis(gen) < this->p)
					{
						this->neighbours[i].push_back(j);
						this->neighbours[j].push_back(i);
					}
				}
			}
		}
};
