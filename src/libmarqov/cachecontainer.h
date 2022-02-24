#ifndef CACHECONTAINER_H
#define CACHECONTAINER_H
/*
MIT License

Copyright (c) 2020-2021 Florian Goth
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
#include <stdexcept>
#include <complex>
#include <H5Cpp.h>
#include <H5File.h>

/* Some predefined HDF5 helpers --------------------------------*/

/** A Class that maps basic, POD C++ Types to HDF5 types.
 */
template <class T>
struct H5MapperBase;

/** A specialization for double.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<double> {
    /**
     * @returns the HDF5 Type of double
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_DOUBLE;}
};

/** A specialization for float.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<float> {
    /**
     * @returns the HDF5 type of float
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_FLOAT;}
};

/** A specialization for int16_t.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<int16_t> {
    /**
     * @returns the HDF5 Type of int16_t
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_INT16;}
};

/** A specialization for int32_t.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<int32_t> {
    /**
     * @returns the HDF5 type of int32_t
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_INT32;}
};

/** A specialization for int64_t.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<int64_t> {
    /**
     * @returns the HDF5 type of int64_t
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_INT64;}
};

/** A specialization for long unsigned int.
 * @see H5MapperBase
 */
template <> struct H5MapperBase<long unsigned int> {
    /**
     * @returns the HDF5 type of long unsigned int.
     */
    static H5::PredType H5Type(){return H5::PredType::NATIVE_ULONG;}
};

/** This maps C++ Types to HDF5 Types.
 * 
 * This is the generic function, valid for tuples, arrays, all 
 * things that have a known length at compile time.
 * 
 * @tparam T The type whose properties we try to infer.
 * @tparam Enable Internal; for distinguishing Scalars from vectors.
 */
template <typename T, class Enable = void>
class H5Mapper
{
    public:
        static constexpr auto fillval = H5Mapper<typename T::value_type>::fillval;///< The initializer value of HDF5. We reuse that from the scalar base.
        static constexpr int rank = 1;///< This template assumes vectors.
        static constexpr int bytecount = std::tuple_size<T>::value*H5Mapper<typename T::value_type>::bytecount;///< The size of a compile time-array
        /** Generate the proper HDF5 Array Type.
         * @returns The Corresponding HDF5 Type of the array.
         */
        static H5::ArrayType H5Type()
        {
            hsize_t dims[1] = {std::tuple_size<T>::value};
            return H5::ArrayType(H5Mapper<typename T::value_type>::H5Type(), rank, dims);
        }
};

/**
 * The partial class specialization of the H5Mapper for scalar data types.
 * 
 * @see H5Mapper
 */
template <typename T>
class H5Mapper<T, typename std::enable_if<std::is_scalar<T>::value>::type> : public H5MapperBase<T>
{
    public:
        static constexpr T fillval = 0; ///< All scalar types get initialized to 0.
        static constexpr int bytecount = sizeof(T);///< Their bytecount on the current platform.
        static constexpr int rank = 1; ///< Scalars are vectors(rank=1) of length 1.
};

/**
 * This specialization enables the support for complex datatypes in MARQOV.
 * 
 * @see H5Mapper
 */
template <typename FPType>
class H5Mapper<std::complex<FPType>, void>
{
    public:
        static constexpr std::complex<FPType> fillval{FPType(0)}; ///< All complex types get initialized to 0+I*0
        static constexpr std::size_t bytecount = sizeof(std::complex<FPType>);///< Their bytecount on the current platform.
        static constexpr int rank = 1; ///< Scalars are vectors(rank=1) of length 1.
        /** Generate the proper HDF5 Array Type.
         * @returns The Corresponding HDF5 Type of the array.
         */
        static auto H5Type()
        {
            std::complex<FPType> dummy;
            H5::CompType cmplx(bytecount);
            cmplx.insertMember("r", 0, H5MapperBase<FPType>::H5Type());
            cmplx.insertMember("i", sizeof(FPType), H5MapperBase<FPType>::H5Type());
            return cmplx;
        }
};

/** Dumps a Scalar into an HDF5 file/group.
 * 
 * Since it's often needed for config related things, this is a small helper
 * function to dump a single scalar into its own scalar data space in a given
 * group.
 * @tparam T the type of the scalar.
 * @param h5loc the HDF5 group / location of the variable.
 * @param key how to name the value in the HDF5 File.
 * @param val the value of the scalar. 
 */
template <typename T>
inline void dumpscalartoH5(H5::Group& h5loc, std::string key, const T& val)
{
    H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
    H5::DataSet dset(h5loc.createDataSet(key.c_str(), H5Mapper<T>::H5Type(), dspace));
    dset.write(&val, H5Mapper<T>::H5Type());
}

/** Dumps a string into an HDF5 file.
 * 
 * Since it's often needed for config related things, this is a small
 * helper that dumps a string into a given group. Strings are interpreted in 
 * HDF5 as arrays of characters.
 * @see dumpscalartoH5
 * @param h5loc the HDF5 group / location of the variable.
 * @param key how to name the value in the HDF5 File.
 * @param value the value of the scalar. 
 */
inline void dumpscalartoH5(H5::Group& h5loc, std::string key, std::string value)
{
    H5::StrType strdatatype(H5::PredType::C_S1, value.size());
    H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
    H5::DataSet dset(h5loc.createDataSet(key.c_str(), strdatatype, dspace));
    dset.write(value.c_str(), strdatatype);
}

/** A helper structure to encapsulate the arguments of a single CacheContainer.
 * 
 *  The rationale was that it was ambiguous to use tuples of tuples. tuple_cat()
 *  would flatten all tuples.
 * @see CacheContainer::CacheContainer
 */
struct CacheContainerArgs
{
    /** Constructor
     * 
     * @param file The group where we write data.
     * @param oname The name of the observable.
     * @param d An extended description of the observable.
     * @param cs The size of the memory that the cache should use.
     */
    CacheContainerArgs(H5::Group& file, const std::string& oname, std::string d = std::string(), std::size_t cs=4194304) :
    hfile(file), obsname(oname), desc(d), cachesize(cs) {}
    H5::Group& hfile; ///< The group where we write data.
    const std::string& obsname; ///< The name of the observable.
    const std::string desc; ///< An extended description of the observable.
    std::size_t cachesize; ///< The size of the memory that the cache should use.
};

/** 
 * The CacheContainer.
 * 
 *  It is associated with a dataspace in an already open HDF5 file. It performs caching
 *  so that not every new data point leads to I/O. C++ stack unwinding takes care of proper
 *  tidy up.
 * @tparam T the type that we put into the container and into the cache
 * @tparam Cont the Container that we use for intermediate storage.
 */
template <class T, class Cont = std::vector<T> >
class CacheContainer
{
public:
    /** construct a cache container.
     * 
     * @param hfile the HDF5 file that we use to dump the data.
     * @param name the name of the data set in HDF5.
     * @param desc A description of the data.
     * @param cachesize the memory size in bytes to use for caching. Will be rounded to integers of datatypes.
     */
    CacheContainer(H5::Group& hfile, const std::string& name, std::string desc = std::string(), std::size_t cachesize=4194304) : cachepos(0)
    {
            if(name.empty()) throw std::runtime_error("[MARQOV::CacheContainer] ERROR: Name of time series not specified!");
            cachemaxelems = cachesize/sizeof(T);
            constexpr int rank = H5Mapper<T>::rank;
            std::array<hsize_t, rank> maxdims, chunk_dims;
            hsize_t fdims[rank] = {};
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
    /** A helper constructor that forwards to the main constructor.
     * 
     * @see CacheContainerArgs
     * @param args the Argument helper structure of type CacheContainerArgs
     */
    CacheContainer(CacheContainerArgs args) : CacheContainer<T, Cont>(args.hfile, args.obsname, args.desc, args.cachesize) {}
    /** Destructor that takes care of flushing the cache.
     */
    ~CacheContainer()
    {
        this->writecache();
    }
    /** Pushes data into the cache.
     * 
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
     * 
     * @param t The element that we write.
     * @returns The modified cache container.
     */
    CacheContainer& operator<<(const T& t) {this->push(t); return *this;}
private:
        /**
        * Writes out the entire current cache of the observable.
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
    H5::DataSet dataset; ///< The HDF5 dataset.
    hsize_t dssize; ///< the current dataset size.
    std::size_t cachemaxelems; ///< How many elements can the cache hold.
    std::size_t cachepos; ///< the current position of the cache.
    Cont cont;///< the container where the data is held until it is flushed.
};

/** A specialization of the CacheContainer for vector.
 * During the caching operation we copy the data to an intermediate linearized version of the time series.
 *  It is associated with a dataspace in an already open HDF5 file. It performs caching
 *  so that not every new data point leads to I/O. C++ stack unwinding takes care of proper
 *  tidy up.
 * @see CacheContainer
 * @tparam T the type of the elements of the vector that we put into the container and into the cache.
 * @tparam Cont the Container that we use for intermediate storage.
 */
template <class T, class Cont>
class CacheContainer<std::vector<T>, Cont>
{
public:
    /** Construct a container for C++ vectors.
     *
     *  @param hf the HDF5 File that we use to dump the data
     *  @param n the name of the data set in HDF5
     *  @param d A description of the data
     *  @param cs the memory size in bytes to use for caching. Will be rounded to integers of datatypes
     */
    CacheContainer(H5::Group& hf, const std::string& n, std::string d = std::string(), std::size_t cs=4194304) : hfile(hf), dssize(0), cachepos(0), unused(true), cachesize(cs), name(n), desc(d)
    {
        if(n.empty()) throw std::runtime_error("[MARQOV::CacheContainer] ERROR: Name of time series not specified!");
    }
    /** A helper constructor that forwards to the main constructor.
     * 
     * @param args the Argument helper structure
     */
    CacheContainer(CacheContainerArgs args) : CacheContainer<std::vector<T>, Cont>(args.hfile, args.obsname, args.desc, args.cachesize) {}
    /** Destructor that takes care of flushing the cache.
     */
    ~CacheContainer()
    {
        if(!unused)
            this->writecache();
    }
    /** Pushes data into the cache.
     * 
     * @param data the element that we write
     */
    void push(const std::vector<T>& data) {
        if (unused) //only at the very first push event do we know the size of the array
        {
            initdataspace(data);
            unused = false;
        }
        for(uint i = 0; i < data.size(); ++i)
        cont[cachepos + i] = data[i];
        cachepos += data.size();
        if (cachepos >= cachemaxelems)
        {
            this->writecache();
            cachepos = 0;
        }
    }
    /** Convenience function for pushing data.
     * @param t the element that we write.
     * @returns the modified cache container.
     */
    CacheContainer& operator<<(const std::vector<T>& t) {this->push(t); return *this;}
private:
    /** Init the HDF data space and our internal data structures.
     * 
     * @see push
     * @param t the vector that we use for initializing the data space. Note that this defines the elements of the data space!
     */
    void initdataspace(const std::vector<T>& t)
    {
        cachemaxelems = (cachesize/sizeof(T)/t.size())*t.size();//properly determine the cache size for vector observables that are non-commensurate with the default cache size
        if (cachemaxelems == 0) cachemaxelems = sizeof(T)*t.size();//If a single measurement of the observable does not fit into the cache, just make it big enough.
        constexpr int rank = H5Mapper<T>::rank;
        std::array<hsize_t, rank> maxdims, chunk_dims;
        hsize_t fdims[rank] = {};
        maxdims.fill(H5S_UNLIMITED);
        
        H5::DataSpace mspace1(rank, fdims, maxdims.data());
        H5::DSetCreatPropList cparms;

        chunk_dims.fill(4096*1024/H5Mapper<T>::bytecount/t.size() + 1);//4MB chunking, irrespective of the length of the vector
        cparms.setChunk(rank, chunk_dims.data() );
        cparms.setDeflate(9);//Best (1-9) compression

        hsize_t dims[1] = {t.size()};
        arrtype = H5::ArrayType(H5Mapper<T>::H5Type(), rank, dims);
        if(t.size() < 8192)//FIXME: This seems to be the practical limit on my Debian Buster system.
        {
            auto fv = std::vector<T>(t.size(), T(H5Mapper<T>::fillval));//The T in front triggers an instantiation...
            cparms.setFillValue(arrtype, fv.data());
        }

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
        arrtype.getArrayDims(&count[0]);//temporarily abuse the count array
        uint nelems = cachepos/count[0]; // hardcoded 1D array
        dims.fill(nelems);
        H5::DataSpace mspace(rank, dims.data(), NULL);
        start.fill(dssize);

        dssize += nelems;
        dataset.extend(&dssize);
        auto filespace = dataset.getSpace();
        count.fill(nelems);
        filespace.selectHyperslab(H5S_SELECT_SET, count.data(), start.data());
        dataset.write(cont.data(), arrtype, mspace, filespace);
    }
    H5::Group& hfile; ///< The HDF5 Group to which we store the time series.
    H5::DataSet dataset; ///< The HDF5 dataset.
    H5::ArrayType arrtype; ///< This is used to store the actual type that gets determined for this vector.
    hsize_t dssize; ///< the current dataset size.
    std::size_t cachemaxelems; ///< How many elements can the cache hold.
    std::size_t cachepos; ///< the current position of the cache.
    std::vector<T> cont; ///< the container where the data is held until it is flushed.
    bool unused; ///< a logical value that becomes true after the first write of data.
    uint cachesize; ///< The amount of cache that we use.
    std::string name; ///< The name of the time series.
    std::string desc; ///< The description of the time series.
};

#endif
