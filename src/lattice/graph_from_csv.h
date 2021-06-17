#ifndef GRAPHFROMCSV_H 
#define GRAPHFROMCSV_H

#include "util/io.h"

// read minimal graph from csv file
class GraphFromCSV
{
    public:
    int npoints;

    std::vector<std::vector<int>> neighbours;

    GraphFromCSV(std::string filename) 
    {
        npoints = import_geometry(filename, neighbours);
    }
            
    std::vector<int> nbrs(const int a, const int i) const
    {
        return neighbours[i];
    }

    std::size_t size() const {return npoints;}
};


// read graph with coordinates from csv file
class GraphFromCSVwithCoords
{
    public:
    int npoints;

    std::vector<std::vector<int>> neighbours;
    std::vector<std::vector<double>> grid;

    GraphFromCSV(std::string filename, ncoords) 
    {
        npoints = import_geometry(filename, grid, neighbours, ncoords);
    }
            
    std::vector<int> nbrs(const int a, const int i) const
    {
        return neighbours[i];
    }
    
    std::vector<double> crds(const int i) const
    {
        return grid[i];
    }

    std::size_t size() const {return npoints;}
};




#endif
