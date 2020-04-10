#ifndef CACHECONTAINER_H
#define CACHECONTAINER_H


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

template <>
class H5Mapper<double>
{
	public:
		static constexpr double fillval = 0;
		static constexpr int rank = 1;
		static auto H5Type(){return H5::PredType::NATIVE_DOUBLE;}
};

template <>
class H5Mapper<int>
{
	public:
		static constexpr int fillval = 0;
		static constexpr int rank = 1;
		static auto H5Type(){return H5::PredType::NATIVE_INT;}
};

template <typename Tp>
class H5Mapper
{
	public:
		static constexpr int fillval = H5Mapper<typename Tp::value_type>::fillval;
		static constexpr int rank = std::tuple_size<Tp>::value;
		static auto H5Type(){return H5Mapper<typename Tp::value_type>::H5Type();}
};

/** A helper structure to encapsulate the arguments of a single CacheContainer.
 *  The rationale was that it was ambiguous to use tuples of tuples. tuple_cat
 *  would flatten all tuples
 */
struct CacheContainerArgs
{
    CacheContainerArgs(H5::H5File& file, const std::string& n, std::size_t cs=4194304) :
    hfile(file), obsname(n), cachesize(cs) {}
    H5::H5File& hfile;
    const std::string& obsname;
    std::size_t cachesize;
};

/** The CacheContainer.
 *  It is associated with a dataspace in an already open HDF5 file. It performs caching
 *  so that not every new data point leads to I/O. C++ stack unwinding takes care of proper
 *  tidy up
 */
template <class T, class Cont = std::vector<T> >
class CacheContainer
{
public:
    /** @param hfile the HDF5 File that we use to dump the data
     *  @param name the name of the data set in HDF5
     *  @param cs the memory size in Bytes to use for caching. Will be rounded to integers of datatypes
     */
    CacheContainer(H5::H5File& hfile, const std::string& name, std::size_t cachesize=4194304) : cachepos(0)
    {
        cachemaxelems = cachesize/sizeof(T);
			constexpr int rank = H5Mapper<T>::rank;
			hsize_t fdims[rank] = {0}; // dim sizes of ds (on disk)
			hsize_t maxdims[rank] = {H5S_UNLIMITED};

			H5::DataSpace mspace1(rank, fdims, maxdims);
			H5::DSetCreatPropList cparms;
			auto fv = H5Mapper<T>::fillval;

			hsize_t chunk_dims[1] = {4096*1024/sizeof(T)};//4MB chunking
			cparms.setChunk( rank, chunk_dims );
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
    /* Convenience function for pushing data
     */
    CacheContainer& operator<<(const T& t) {this->push(t); return *this;}
private:
		/**
         * Writes out the entire current cache of the observable
         */
		void writecache ()
        {
			constexpr int rank = H5Mapper<T>::rank;
			hsize_t dims[rank] = {cachepos};
			H5::DataSpace mspace(rank, dims, NULL);
			hsize_t start[rank] = {dssize};
			dssize += cachepos;
			dataset.extend(&dssize);
			auto filespace = dataset.getSpace();
			hsize_t count[rank] = {cachepos};
			filespace.selectHyperslab(H5S_SELECT_SET, count, start);
			dataset.write(cont.data(), H5Mapper<T>::H5Type(), mspace, filespace);
        }
    H5::DataSet dataset; ///< The HDF5 dataset
    hsize_t dssize; //< the current dataset size
    std::size_t cachemaxelems; ///< How many elements can the cache hold
    std::size_t cachepos; ///< the current position of the cache
    Cont cont;///< the container where the data is held until it is flushed
};

#endif