#ifndef IO_H
#define IO_H


// export coordinates and neighbour relations

template <class Grid>
void save_geometry(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

	const int dim = grid.dim;

	for (int i=0; i<grid.size(); i++)
	{
		auto nbrs = grid.nbrs(0,i);
		auto crds = grid.crds(i);

		os << nbrs.size() + crds.size();

		for (int j=0; j<crds.size(); j++) os << "\t" << crds[j];
		for (int j=0; j<nbrs.size(); j++) os << "\t" << nbrs[j];

		os << std::endl;
	}
}



// export coordinates and neighbour relations and bond strenghts (only scalars!)

template <class Grid>
void save_geometry_deluxe(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

	const int dim = grid.dim;

	for (int i=0; i<grid.size(); i++)
	{
		auto nbrs = grid.nbrs(0,i);
		auto crds = grid.crds(i);
		auto bnds = grid.getbnds(0,i);

		os << nbrs.size() + crds.size() + bnds.size();

		for (int j=0; j<crds.size(); j++) os << "\t" << crds[j];
		for (int j=0; j<nbrs.size(); j++) os << "\t" << nbrs[j];
		for (int j=0; j<bnds.size(); j++) os << "\t" << bnds[j];

		os << std::endl;
	}
}

#endif
