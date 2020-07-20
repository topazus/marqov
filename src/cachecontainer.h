#ifndef CACHECONTAINER_H
#define CACHECONTAINER_H
/*
MIT License

Copyright (c) 2020 Florian Goth
fgoth@physik.uni-wuerzburg.de

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

#include <vector>
#include <string>
#include <type_traits>
#include <utility>
#include <tuple>
#include <H5Cpp.h>
#include <H5File.h>

/*some predefined HDF5 helpers --------------------------------*/

template <typename T>
class H5Mapper;

/** First we have the POD Type H5 Types
 */
template <>
class H5Mapper<double>
{
    public:
        static constexpr double fillval = 0;
        static constexpr int bytecount = sizeof(double);
        static constexpr int rank = 1;
        static auto H5Type(){return H5::PredType::NATIVE_DOUBLE;}
};

template <>
class H5Mapper<int>
{
    public:
        static constexpr int fillval = 0;
        static constexpr int bytecount = sizeof(int);
        static constexpr int rank = 1;
        static auto H5Type(){return H5::PredType::NATIVE_INT;}
};

/** This maps a 1D vector/array like structure to a custom HDF5 datatype.
 */
template <typename Tp>
class H5Mapper
{
    public:
        static constexpr auto fillval = H5Mapper<typename Tp::value_type>::fillval;
        static constexpr int rank = 1;
        static constexpr int bytecount = std::tuple_size<Tp>::value*H5Mapper<typename Tp::value_type>::bytecount;
        static auto H5Type(){
            hsize_t dims[1] = {std::tuple_size<Tp>::value};
            return H5::ArrayType(H5Mapper<typename Tp::value_type>::H5Type(), rank, dims);
        }
};

/** A helper structure to encapsulate the arguments of a single CacheContainer.
*  The rationale was that it was ambiguous to use tuples of tuples. tuple_cat()
*  would flatten all tuples.
*/
struct CacheContainerArgs
{
    CacheContainerArgs(H5::Group& file, const std::string& n, std::size_t cs=4194304) :
    hfile(file), obsname(n), cachesize(cs) {}
    H5::Group& hfile;
    const std::string& obsname;
    std::size_t cachesize;
};

/** The CacheContainer.
*  It is associated with a dataspace in an already open HDF5 file. It performs caching
*  so that not every new data point leads to I/O. C++ stack unwinding takes care of proper
*  tidy up.
*/
template <class T, class Cont = std::vector<T> >
class CacheContainer
{
public:
    /** @param hfile the HDF5 File that we use to dump the data
    *   @param name the name of the data set in HDF5
    *   @param cs the memory size in Bytes to use for caching. Will be rounded to integers of datatypes
    */
    CacheContainer(H5::Group& hfile, const std::string& name, std::size_t cachesize=4194304) : cachepos(0)
    {
            cachemaxelems = cachesize/sizeof(T);
            constexpr int rank = H5Mapper<T>::rank;
            std::array<hsize_t, rank> maxdims, chunk_dims;
            hsize_t fdims[rank] = {0};
            maxdims.fill(H5S_UNLIMITED);

            H5::DataSpace mspace1(rank, fdims, maxdims.data());
            H5::DSetCreatPropList cparms;
            auto fv = H5Mapper<T>::fillval;
            
            chunk_dims.fill(4096*1024/H5Mapper<T>::bytecount);//4MB chunking
            cparms.setChunk( rank, chunk_dims.data() );
            cparms.setDeflate(9);//Best (1-9) compression
            cparms.setFillValue(  H5Mapper<T>::H5Type(), &fv);
            dataset = hfile.createDataSet(name, H5Mapper<T>::H5Type(), mspace1, cparms);
            dssize = 0;
            cont.resize(cachemaxelems);//allocate space for 1024 entries
    }
    /** A helper constructor that forwards to the main constructor
    * @param args the Argument helper structure
    */
    CacheContainer(CacheContainerArgs args) : CacheContainer<T, Cont>(args.hfile, args.obsname, args.cachesize) {}
    /* Destructor that takes care of flushing the cache.
    */
    ~CacheContainer()
    {
        this->writecache();
    }
    /** pushes data into the cache
    * @param data the element that we write
    */
    void push(const T& data) {
            cont[cachepos] = data;
            cachepos = cachepos + 1;
            if (cachepos >= cachemaxelems)
            {
                this->writecache();
                cachepos = 0;
            }
    }
    /** Convenience function for pushing data.
     * @param data the element that we write
     */
    CacheContainer& operator<<(const T& t) {this->push(t); return *this;}
private:
        /**
        * Writes out the entire current cache of the observable
        */
        void writecache ()
        {
            constexpr int rank = H5Mapper<T>::rank;
            std::array<hsize_t, rank> dims, start, count;
            dims.fill(cachepos);
            H5::DataSpace mspace(rank, dims.data(), NULL);
            start.fill(dssize);

            dssize += cachepos;
            dataset.extend(&dssize);
            auto filespace = dataset.getSpace();
            count.fill(cachepos);
            filespace.selectHyperslab(H5S_SELECT_SET, count.data(), start.data());
            dataset.write(cont.data(), H5Mapper<T>::H5Type(), mspace, filespace);
        }
    H5::DataSet dataset; ///< The HDF5 dataset
    hsize_t dssize; //< the current dataset size
    std::size_t cachemaxelems; ///< How many elements can the cache hold
    std::size_t cachepos; ///< the current position of the cache
    Cont cont;///< the container where the data is held until it is flushed
};

#endif
