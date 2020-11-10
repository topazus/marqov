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

/** First we have the POD Type H5 Types
 */
template <class T>
struct H5MapperBase;

template <> struct H5MapperBase<double> {
    static auto H5Type(){return H5::PredType::NATIVE_DOUBLE;}
};

template <> struct H5MapperBase<float> {
    static auto H5Type(){return H5::PredType::NATIVE_FLOAT;}
};

template <> struct H5MapperBase<int16_t> {
    static auto H5Type(){return H5::PredType::NATIVE_INT16;}
};

template <> struct H5MapperBase<int32_t> {
    static auto H5Type(){return H5::PredType::NATIVE_INT32;}
};

template <> struct H5MapperBase<int64_t> {
    static auto H5Type(){return H5::PredType::NATIVE_INT64;}
};

template <> struct H5MapperBase<long unsigned int> {
    static auto H5Type(){return H5::PredType::NATIVE_ULONG;}
};

/** This maps a 1D vector/array like structure to a custom HDF5 datatype.
 */
template <typename T, class Enable = void>
class H5Mapper
{
    public:
        static constexpr auto fillval = H5Mapper<typename T::value_type>::fillval;
        static constexpr int rank = 1;
        static constexpr int bytecount = std::tuple_size<T>::value*H5Mapper<typename T::value_type>::bytecount;
        static auto H5Type(){
            hsize_t dims[1] = {std::tuple_size<T>::value};
            return H5::ArrayType(H5Mapper<typename T::value_type>::H5Type(), rank, dims);
        }
};

template <typename T>
class H5Mapper<T, typename std::enable_if<std::is_scalar<T>::value>::type> : public H5MapperBase<T>
{
    public:
        static constexpr T fillval = 0;
        static constexpr int bytecount = sizeof(T);
        static constexpr int rank = 1;
};

//Helper to dump a single key-value pair into its own scalar data space
template <typename T>
inline void dumpscalartoH5(H5::Group& h5loc, std::string key, const T& s)
{
    
    H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
    H5::DataSet dset(h5loc.createDataSet(key.c_str(), H5Mapper<T>::H5Type(), dspace));
    dset.write(&s, H5Mapper<T>::H5Type());
}
inline void dumpscalartoH5(H5::Group& h5loc, std::string key, std::string value)
{
    H5::StrType strdatatype(H5::PredType::C_S1, value.size());
    H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
    H5::DataSet dset(h5loc.createDataSet(key.c_str(), strdatatype, dspace));
    dset.write(value.c_str(), strdatatype);
}

/** A helper structure to encapsulate the arguments of a single CacheContainer.
*  The rationale was that it was ambiguous to use tuples of tuples. tuple_cat()
*  would flatten all tuples.
*/
struct CacheContainerArgs
{
    CacheContainerArgs(H5::Group& file, const std::string& n, std::string d = std::string(), std::size_t cs=4194304) :
    hfile(file), obsname(n), desc(d), cachesize(cs) {}
    H5::Group& hfile;
    const std::string& obsname;
    const std::string desc;
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
    *   @param desc A description of the data
    *   @param cs the memory size in Bytes to use for caching. Will be rounded to integers of datatypes
    */
    CacheContainer(H5::Group& hfile, const std::string& name, std::string desc = std::string(), std::size_t cachesize=4194304) : cachepos(0)
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
            if(!desc.empty())
                dataset.setComment(desc.c_str());
            dssize = 0;
            cont.resize(cachemaxelems);//allocate space for 1024 entries
    }
    /** A helper constructor that forwards to the main constructor
    * @param args the Argument helper structure
    */
    CacheContainer(CacheContainerArgs args) : CacheContainer<T, Cont>(args.hfile, args.obsname, args.desc, args.cachesize) {}
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

/** A specialization of the CacheContainer for vector
*  It is associated with a dataspace in an already open HDF5 file. It performs caching
*  so that not every new data point leads to I/O. C++ stack unwinding takes care of proper
*  tidy up.
*/
template <class T, class Cont>
class CacheContainer<std::vector<T>, Cont>
{
public:
    /** @param hf the HDF5 File that we use to dump the data
    *   @param n the name of the data set in HDF5
    *   @param d A description of the data
    *   @param cs the memory size in Bytes to use for caching. Will be rounded to integers of datatypes
    */
    CacheContainer(H5::Group& hf, const std::string& n, std::string d = std::string(), std::size_t cs=4194304) : dssize(0), cachepos(0), hfile(hf), unused(true), cachesize(cs), name(n), desc(d)
    {}
    /** A helper constructor that forwards to the main constructor
    * @param args the Argument helper structure
    */
    CacheContainer(CacheContainerArgs args) : CacheContainer<std::vector<T>, Cont>(args.hfile, args.obsname, args.desc, args.cachesize) {}
    /* Destructor that takes care of flushing the cache.
    */
    ~CacheContainer()
    {
        this->writecache();
    }
    /** pushes data into the cache
    * @param data the element that we write
    */
    void push(const std::vector<T>& data) {
        if (unused) //only a the very first push event do we know the size of the array
        {
            initdataspace(data);
            unused = false;
        }
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
    CacheContainer& operator<<(const std::vector<T>& t) {this->push(t); return *this;}
private:
    /** Init the HDF data space and our internal data structures.
     */
        void initdataspace(const std::vector<T>& t)
        {
            cachemaxelems = cachesize/sizeof(T)/t.size();
            constexpr int rank = H5Mapper<T>::rank;
            std::array<hsize_t, rank> maxdims, chunk_dims;
            hsize_t fdims[rank] = {0};
            maxdims.fill(H5S_UNLIMITED);

            H5::DataSpace mspace1(rank, fdims, maxdims.data());
            H5::DSetCreatPropList cparms;
            auto fv = H5Mapper<T>::fillval;
            
            chunk_dims.fill(4096*1024/H5Mapper<T>::bytecount/t.size());//4MB chunking
            cparms.setChunk( rank, chunk_dims.data() );
            cparms.setDeflate(9);//Best (1-9) compression
            cparms.setFillValue(H5Mapper<T>::H5Type(), &fv);
            hsize_t dims[1] = {t.size()};
            auto arrtype = H5::ArrayType(H5Mapper<T>::H5Type(), rank, dims);

            dataset = hfile.createDataSet(name, arrtype, mspace1, cparms);
            if(!desc.empty())
                dataset.setComment(desc.c_str());
            dssize = 0;
            cont.resize(cachemaxelems);//allocate space for 1024 entries
        }
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
    H5::Group& hfile;
    H5::DataSet dataset; ///< The HDF5 dataset
    hsize_t dssize; //< the current dataset size
    std::size_t cachemaxelems; ///< How many elements can the cache hold
    std::size_t cachepos; ///< the current position of the cache
    Cont cont;///< the container where the data is held until it is flushed
    bool unused;
    uint cachesize;
    std::string name;
    std::string desc;
};


#endif
