/*
MIT License

Copyright (C) 2018 - 2020  Manuel Schrauth, Jefferson S.E. Portela, Florian Goth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef NEIGHBOURCLASS_H
#define NEIGHBOURCLASS_H
#include <map>
#include <vector>
#include <algorithm>

/** This function prepares the input data gridx, gridy and n
 * so that n is stored in a format usable by the neighbour class.
 * This task can not be done by the neighbour class since it only encodes
 * the neighbourhood relations but not the coordinates. Nevertheless the coordinates need to 
 * be changed, too.
 * Again: After exection of this function gridx, gridy, and n have changed data!
 *  low amount of neighbours can be found at the beginning
 *  large amounts at the end of n.
 *  For a particular point the neighbours are sorted in ascending order.
 * @param gridx The array of the x coordinats
 * @param gridy the array of the y components
 * @param n the vector containing the neighbour relations that needs to be transformed.
 * */
void prepare_for_neighbours_class(std::vector<double>& gridx, std::vector<double>& gridy, 
                                  std::vector< std::vector<int32_t> >& n)
{
    std::vector<uint> mymap;
    for(const std::vector<int32_t>& v : n)
    {
        mymap.push_back(v.size());
    }
    auto ret = std::minmax_element(mymap.begin(), mymap.end());
    int minneighbours = *ret.first;
    int classes = 1 + *ret.second - minneighbours;
    std::vector<int>* arr = new std::vector<int>[classes];
    for(uint i = 0; i < n.size(); ++i)
    {
        arr[n[i].size() - minneighbours].push_back(i);
    }
    //arr is now an array where we can find the indices of neighbours with a particular number of neighbours
    std::vector<int> newarr;
    for(int i = 0; i < classes; ++i)
        for(uint j = 0; j < arr[i].size(); ++j)
            newarr.push_back(arr[i][j]);
        //newarr is the mapping which new index corresponds to which old index
    std::vector<int> invnewarr(newarr.size());
    for(std::size_t i = 0; i < newarr.size(); ++i)
        invnewarr[newarr[i]] = i;
    //invnewarr holds the info which old index corresponds to which new index
    // now we should have all data to map the indices to the new format.
    delete [] arr;
    std::vector<double> gridxnew(n.size());
    std::vector<double> gridynew(n.size());
    std::vector< std::vector<int32_t> > nnew(n.size());
    for (uint i = 0; i < n.size(); ++i)
    {
        gridxnew[i] = gridx[newarr[i]];
        gridynew[i] = gridy[newarr[i]];
        const std::vector<int32_t>& oldpoint = n[newarr[i]];
        nnew[i].resize(oldpoint.size());
        for(uint j = 0; j < oldpoint.size(); ++j)
            nnew[i][j] = invnewarr[oldpoint[j]];
        //impose some more structure by sorting the neighbour indices.
        std::sort(nnew[i].begin(), nnew[i].end());
    }
    gridx = gridxnew;
    gridy = gridynew;
    n = nnew;
    //The following commented out block can be used to dump the resulting structure
    //     std::ofstream os("test.dat");
    //     for(uint i = 0; i < n.size(); ++i)
    //     {
    //         os<<2+n[i].size()<<" "<<gridx[i]<<" "<<gridy[i]<<" ";
    //         for(uint j = 0; j < n[i].size(); ++j)
    //             os<<n[i][j]<<" ";
    //         os<<std::endl;
    //     }
    //     os.close();
}

/**
* This class encapsulates the edges of the graph.
* It exploits that usually the number of neighbours that the vertices have are small integers 
* whose distribution is not very wide. In the case of lattices with constant coordination numbers for
* example there is no variability and all vertices have the same number of neighbours.
* The number of neighbours is called a "class" of a vertex.
* To obtain a simple linear data structure we have sorted all vertices accordng to the amount of neighbours they have.
* Then every lookup into the array is a twostage process:
* 1.) Determine class of a vertex via a lookup into the bookkeeping array.
      Now we have an offset into the data array where all vertices of a particular class are stored.
* 2.) find the position of the neighbour of the vertex from the index of the vertex.
*/
template <typename IntType = int32_t>
class Neighbours
{
public:
    /** This function returns the position of a particular point in the neighbour array.
     * can be found
     * @param c The class of the point
     * @param pos The index of the point
     * */
    inline std::size_t mypos(int c, const IntType& pos) const
    {
        int res = bookkeeping[c].arr_idx_beg + bookkeeping[c].nbr_class*(pos - bookkeeping[c].nbr_idx_beg);
        return res;
    }
    /** This function queries to which neighour class an index belongs.
     * Hence it can be used to query the number of neighbours.
     * @param pos The position of the index that we ask
     * @return The number of neighbours of point pos.
     * */
    inline int get_class(const std::size_t& pos) const
    {
        int c = 0;
        while (((c) < nr_of_starts) && (pos >= bookkeeping[c].nbr_idx_beg )) ++c;
        return c-1;
    }
//     inline IntType& operator[](const IntType& pos) const
//     {
//         return data[pos*nbsize];
//     }

    inline auto nbrs(int fam, int pos) const
    {
        //Determine neighbour class (i.e. has it three or four or five neighbours)
        int c = get_class(pos);
        //determine the offset of the relevant data into the array
        std::size_t seg = mypos(c, pos);
        return std::vector<IntType>(data + seg, data + seg + c);
    }

    /** This function gives a particular neighbour of a site at a particular position.
     * @param pos The index of the position.
     * @param idx The index of the neighbour
     * @return The index that the neighbour has in this array.
    * */    
    inline IntType& get_neighbour(const IntType& pos, const IntType& idx) const
    {
        //Determine neighbour class (i.e. has it three or four or five neighbours)
        int c = get_class(pos);
        //determine the offset of the relevant data into the array
        std::size_t seg = mypos(c, pos);
        return data[seg + idx];
    }

    inline IntType getRandomNeighbour(std::size_t idx, uint16_t number_of_nn, const IntType& pos, int& rndnumber) const
//    inline IntType getRandomNeighbour(std::size_t idx, uint16_t number_of_nn, const IntType& pos, int rndnumber) const
    {
//	   /*
        switch (number_of_nn)
        {
            case 1:
                break;
            case 2:
                idx += (rndnumber)%2;
                break;
            case 3:
                idx += (rndnumber)%3;
                break;
            case 4:
                idx += (rndnumber)%4;
                break;
            default:
                idx += rndnumber%number_of_nn;
                break;
        }
//	   */
//	   idx += rndnumber;
        return data[idx];
    }

    inline Neighbours(const std::vector<std::vector<IntType> >& n) : Neighbours(&n) {}
    inline Neighbours(const std::vector<std::vector<IntType> >* n)
    {
        /*A helper structure to store some properties that we need for the construction of the data array*/
        struct Props {
            uint occurences; ///< Used to store how often a particular amount of neighbours occured
            IntType min; ///< used to store the minimum neighbour that occured for a particular class
            IntType max; ///< used to store the maximum neighbour that occured for a particular class
            uint lowest_idx; ///< used to find out where the first occurence of a particular class was
            Props () : occurences(0), min(-1), max(0), lowest_idx(-1) {}
        };
        std::map<uint16_t, Props> mymap;
        for (uint i = 0; i < n->size(); ++i )
        {
            uint myneighbours_count = (*n)[i].size();
            Props& prop = mymap[myneighbours_count];//if not found the [] operator calls the default copy constructor.
            ++prop.occurences;
            //we know that the neighbours are sorted in ascending order, therefore we easily find min and max
            const std::vector<IntType>& point = (*n)[i];
            if (point.size() > 0)
            {
                if (point[0] < prop.min) prop.min = point[0];
                if (point.back() > prop.max) prop.max = point.back();
            }
            prop.lowest_idx = std::min(i, prop.lowest_idx);
        }
        //calculate required amount of memory elements
        std::size_t memsize = 0;
        for (auto it = mymap.cbegin(); it != mymap.cend(); ++it)
            memsize += it->first*it->second.occurences;
        //layout: first we store all with one neighbour, then all with two and so on...
        nr_of_starts = mymap.size();
        int ret = posix_memalign((void**)&bookkeeping, 64, nr_of_starts*sizeof(typename Neighbours::Info));//request cache line aligned memory
        if(ret != 0) throw("Can't allocate memory required for bookkeeping!!");
        // fill our array with the number of starts
        auto it = mymap.cbegin();
        bookkeeping[0].nbr_class = it->first;
        bookkeeping[0].arr_idx_beg = 0;
        bookkeeping[0].nbr_idx_beg = it->second.lowest_idx;
        for (int16_t i = 1; i < nr_of_starts; ++i)
        {
            auto prev = it;
            ++it;
            bookkeeping[i].nbr_class = it->first;
            bookkeeping[i].arr_idx_beg = prev->second.occurences*prev->first + bookkeeping[i-1].arr_idx_beg;
            bookkeeping[i].nbr_idx_beg = it->second.lowest_idx;
        }
        ret = posix_memalign((void**)&data, 4096, memsize*sizeof(IntType));//request page-aligned memory
        if(ret != 0) throw("Can't allocate memory for neighbour array!!");
        if (memsize < 1E6)
        {
//            std::cout<<"Allocated "<<memsize*sizeof(IntType)/1024<<" kB of memory."<<std::endl;
        }
        else
        {
//            std::cout<<"Allocated "<<memsize*sizeof(IntType)/1024/1024<<" MB of memory."<<std::endl;
        }
        //copy data to our new array
        int cur_pos = 0;
        for (uint i = 0; i < n->size(); ++i)
        {
            int nbrs = (*n)[i].size();
            for(int j = 0; j < nbrs; ++j)
                data[cur_pos++]  = static_cast<IntType>((*n)[i][j]);
        }
        nr_of_points = n->size();
    }
    inline ~Neighbours() {
        if (data != NULL) free(data);
        if (bookkeeping != NULL) free(bookkeeping);
    }
    std::size_t size() const {return nr_of_points;}

    uint16_t neighbours_in_class(uint c) const {return bookkeeping[c].nbr_class;}

    /*This friend declaration enables a function from readwrite.h to write to this data structure directly. This is needed for a more memory efficient I/O.*/ 
    friend void import_geometry_directly(const int N, std::vector<double>& gridx, std::vector<double>& gridy, Neighbours<int32_t>& outn, const std::string path );

    inline Neighbours() : data(NULL), bookkeeping(NULL), nr_of_starts(0), nr_of_points(0) {}

    /**
     * A move assignment operator to avoid temporary copies.
     * @param n the temporary object.
     * @return An object with the state of n.
	*/
    inline Neighbours& operator=(Neighbours&& n)
    {
        if (this != &n)
        {
            //tidy up our own data
            free(data);
            free(bookkeeping);
            //copy over state from other temporary object
            data = n.data;
            bookkeeping = n.bookkeeping;
            nr_of_starts = n.nr_of_starts;
            nr_of_points = n.nr_of_points;
            //invalidate the OTHER class
            n.nr_of_starts = 0;
            n.nr_of_points = 0;
            n.bookkeeping = NULL;
            n.data = NULL;
        }
        return *this;
    }
    /**
     * A move constructor to avoid temporary copies.
     * @param n the temporary object.
     * @return An object with the state of n.
	*/
    inline Neighbours(Neighbours&& n)
    {
        //copy over state from other temporary object
        data = n.data;
        bookkeeping = n.bookkeeping;
        nr_of_starts = n.nr_of_starts;
        nr_of_points = n.nr_of_points;
        //invalidate the OTHER class
        n.nr_of_starts = 0;
        n.nr_of_points = 0;
        n.bookkeeping = NULL;
        n.data = NULL;
    }
private:
    struct Info {
        std::size_t nbr_class;///< How many neighbours do we have?
        std::size_t arr_idx_beg;///< where does this class start in the neighbour array?
        std::size_t nbr_idx_beg;///< At which neighbour index does this start?
    };
    IntType* data; ///< The data array start stores the actual neighbours
    Info* bookkeeping;///< properties of the data array. The array is sorted in ascending order.
    int16_t nr_of_starts; ///< The length of the bookkeeping array
    std::size_t nr_of_points; ///< with this we can easily access the total amount of neighbours we have
};
#endif
