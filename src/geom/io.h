#ifndef IO_H
#define IO_H

template <class Grid>
void save_geometry(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

	const int dim = grid.dim;

	for (int i=0; i<grid.size(); i++)
	{
		auto nbrs = grid.getnbrs(0,i);
		auto crds = grid.getcrds(i);

		os << nbrs.size() + crds.size();

		for (int j=0; j<crds.size(); j++) os << "\t" << crds[j];
		for (int j=0; j<nbrs.size(); j++) os << "\t" << nbrs[j];

		os << std::endl;
	}
}


#endif
