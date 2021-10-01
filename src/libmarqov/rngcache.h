#ifndef RNGCACHE_H
#define RNGCACHE_H
/* MIT License

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

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <random>
#include <utility>

/** A helper to map the C++ RNG classes to portable names.
 * 
 * @tparam T the used RNG.
 */
template <typename T>
struct RNGName
{
    std::string name = "RNG NAME NOT SPECIFIED!!!!!!!!!!!"; ///< name of the RNG available as a plain string.
};

/** Specialization for ranlux48.
 */
template <>
struct RNGName<std::ranlux48>
{
    std::string name = "ranlux48"; ///< name of the RNG available as a plain string.
};

/** Specialization for knuth_b.
 */
template <>
struct RNGName<std::knuth_b>
{
    std::string name = "knuth_b"; ///< name of the RNG available as a plain string.
};

/** Specialization for ranlux48_base.
 */
template <>
struct RNGName<std::ranlux48_base>
{
    std::string name = "ranlux48_base"; ///< name of the RNG available as a plain string.
};

/** Specialization for ranlux24.
 */
template <>
struct RNGName<std::ranlux24_base>
{
    std::string name = "ranlux24_base"; ///< name of the RNG available as a plain string.
};

/** Specialization for the 64bit version of the Mersenne-Twister.
 */
template <>
struct RNGName<std::mt19937_64>
{
    std::string name = "mt19937_64"; ///< name of the RNG available as a plain string.
};

/** Specialization for the 64bit version of the Mersenne-Twister.
 */
template <>
struct RNGName<std::mt19937>
{
    std::string name = "mt19937"; ///< name of the RNG available as a plain string.
};

/** Specialization for the LCG RAND.
 */
template <>
struct RNGName<std::minstd_rand>
{
    std::string name = "minstd_rand"; ///< name of the RNG available as a plain string.
};

/** Specialization for the LCG RAND0.
 */
template <>
struct RNGName<std::minstd_rand0>
{
    std::string name = "minstd_rand0"; ///< name of the RNG available as a plain string.
};

/** A cache for generated random values.
 * 
 * This implements a cache for random values.
 * If the cache is empty we refill it with new random values from a user
 * defined RNG. Also RNGName how to get a proper name for I/O.
 * @tparam RNG The RNG that the user wishes to use.
 */

template <class RNG>
class RNGCache
{
public:
    template <class ...Args>
    /** Constructor call.
     * 
     * We pass on all arguments to the RNG.
     * @param args the parameters of the RNG.
     */
    RNGCache(Args&&... args) : data(NULL), pos(0), rng(std::forward<Args>(args)...)
    {
        int err = posix_memalign( (void**)(&data), pagesize, nrpages*pagesize );//in C++17 we can replace that with aligned_alloc
        if (err != 0)
            throw std::runtime_error("[RNGCache] error while allocating memory");
        fillcache();
    }
    RNGCache(RNGCache&& other) noexcept : data(std::exchange(other.data, nullptr)), pos(std::move(other.pos)), rng(std::move(other.rng))
    {}
    RNGCache() = delete;
    RNGCache(RNGCache& other) = delete;
    RNGCache& operator=(RNGCache&&) = delete;
    RNGCache& operator=(RNGCache&) = delete;
    /** Frees the memory of the cache.
     */
    ~RNGCache()
    {
        free(data);
    }
    /** Get a random uniform integer from [0, range].
     * 
     * By a cmake switch the debiasing of integers can be enabled.
     * <https://arxiv.org/abs/1805.10941>
     * @param range
     * @return a random integer
     */
    inline uint64_t integer(uint64_t range) noexcept
    {
#ifdef DEBIASINTEGERS
        //From here: https://github.com/imneme/bounded-rands/blob/master/bounded64.cpp
        // which takes inspiration from here: https://arxiv.org/abs/1805.10941
        uint64_t x = integer();
        __uint128_t m = __uint128_t(x) * __uint128_t(range);
        uint64_t l = uint64_t(m);
        if (l < range) {
            uint64_t t = -range;
            if (t >= range) {
                t -= range;
                if (t >= range)
                    t %= range;
            }
            while (l < t) {
                x = integer();
                m = __uint128_t(x) * __uint128_t(range);
                l = uint64_t(m);
            }
        }
        return m >> 64;
#else
        return integer()%range;
#endif
    }
    /** Get a random uniform Integer from [0, max()].
     * 
     * @return a random integer
     */
    inline auto integer() noexcept
    {
        if(pos >= nelems)
        {
            pos = 0;
            fillcache();
        }
        return data[pos++];
    }
    /** Get a random uniform float in the range (min, max).
     * 
     *  If no arguments are specified this gives a double from [0,1)
     * @param max the maximum floating point number to return
     * @param min the minimum floating point number to return
     * @return a random double from the interval [min, max)
     */
    inline auto real(double max = 1.0, double min = 0.0) noexcept
    {
        return min + (max-min)*(double(integer())/double(RNG::max()));
    }
    /** The maximum integer that we support.
     */
    static constexpr auto max() noexcept{return RNG::max();}
    /** Fill the cache.
     * 
     * Can also be used to flush the cache,
     * i.e reset the cache after initally setting the state of the RNG.
     */
    void fillcache() noexcept
    {
        for(int i = 0; i < nelems; ++i)
            data[i] = rng(); //We follow the C++11 convention that operator() advances the state of the RNG
    }
    
    /** Dump the internal state of the RNG in a manner that it can be fully constructed from it.
     * 
     * We assume that the RNG supports the same operations as those from the STL.
     * @return a vector of unsigned 64bit Integers that contain the state.
     */
    std::vector<u_int64_t> dumpstate()
    {
        std::stringstream rngstate;
        rngstate<<rng;//peculiar to the C++11 RNGs
        std::vector<u_int64_t> retval;
        u_int64_t t;
        while (rngstate>>t) 
        {
            retval.push_back(t);
        }
        return retval;
    }

    /** Set the state of the RNG.
     * 
     * This sets the internal state of the RNG.
     * The C++11 STL RNGs dump their state as a sequence of unsigned 64bit integers.
     * We assume that we get this state as a vector.
     * 
     * @param vec A vector containing a sequence of ints for the internal state of the RNG.
     */
    void setstate(std::vector<u_int64_t>& vec)
    {
        std::string rngstring;
        for(decltype(vec.size()) i = 0; i < vec.size(); ++i)
        {
           rngstring+=std::to_string(vec[i])+" ";
        }
        std::istringstream ss(rngstring);
        ss>>rng;
    }
private:
    static constexpr int pagesize = 4096;///< The pagesize of the OS. A common value.
    static constexpr int nrpages = 2; ///< How many memory pages we use.
    typedef decltype(std::declval<RNG>().operator()()) result_type;///< the type that the RNG returns.
    static constexpr int nelems = pagesize*nrpages/sizeof(result_type); ///< How many elements we store in the cache.
    result_type* data{nullptr}; ///< The memory we use for the cache.
    int pos{0}; ///< at which integer are we currently.
    RNG rng{}; ///< the RNG.
};

#endif
