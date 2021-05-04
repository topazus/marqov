#ifndef RNGCACHE_H
#define RNGCACHE_H
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

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <random>

/**
 * A helper to map the C++ RNG classes to portable names
 * @tparam T the used RNG.
 */
template <typename T>
struct RNGName;

template <>
struct RNGName<std::ranlux48_base>
{
    std::string name;
    RNGName() : name("ranlux48_base") {}
};

template <>
struct RNGName<std::ranlux24_base>
{
    std::string name;
    RNGName() : name("ranlux24_base") {}
};

template <>
struct RNGName<std::mt19937_64>
{
    std::string name;
    RNGName() : name("mt19937_64") {}
};

template <>
struct RNGName<std::mt19937>
{
    std::string name;
    RNGName() : name("mt19937") {}
};

template <>
struct RNGName<std::minstd_rand>
{
    std::string name;
    RNGName() : name("minstd_rand") {}
};

template <>
struct RNGName<std::minstd_rand0>
{
    std::string name;
    RNGName() : name("minstd_rand0") {}
};

/**
 * A Cache for generated Random values.
 * This implements a cache for random values.
 * If the cache is empty we refill it with new random values from a user
 * defined RNG.
 * @tparam RNG The RNG that the user wishes to use.
 */

template <class RNG>
class RNGCache
{
public:
    template <class ...Args>
    /** Constructor call. We pass on all arguments to the RNG.
     */
    RNGCache(Args&&... args) : data(NULL), pos(0), rng(std::forward<Args>(args)...)
    {
        int err = posix_memalign( (void**)(&data), pagesize, nrpages*pagesize );//in C++17 we can replace that with aligned_alloc
        if (err != 0)
            throw std::runtime_error("[RNGCache] error while allocating memory");
        fillcache();
    }
    RNGCache(RNGCache&& other) = default;//FIXME: think about it...
    ~RNGCache()
    {
        free(data);
    }
    /** Get a random uniform Integer from [0, range]
     * @param range
     * @return a random integer
     */
    inline auto integer(int range) noexcept
    {
        return integer()%range;
    }
    /** Get a random uniform Integer from [0, max()]
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
    static constexpr auto max() {return RNG::max();}
    /** Fill the cache. Can also be used to flush the cache,
     * i.e reset the cache after initally setting the state of the RNG.
     */
    void fillcache() noexcept
    {
        for(int i = 0; i < nelems; ++i)
            data[i] = rng(); //We follow the C++11 convention that operator() advances the state of the RNG
    }
    /** Dump the internal state of the RNG in a manner that it can be fully constructed from it.
     * @return a vector of 64bit Integers that contain the state.
     */
    std::vector<int64_t> dumpstate()
    {
        typedef typename RNG::result_type IntType;
        std::stringstream rngstate;
        rngstate<<rng;//peculiar to the C++11 RNGs
        std::vector<int64_t> retval;
        int64_t t;
        while (rngstate>>t) retval.push_back(t);
        return retval;
    }
    void setstate(std::vector<int64_t>& vec)
    {
        std::string rngstring;
        for(int i = 0; i < vec.size(); ++i)
        {
           rngstring+=std::to_string(vec[i])+" ";
        }
        std::istringstream ss(rngstring);
        ss>>rng;
    }
private:
    static constexpr int pagesize = 4096;
    static constexpr int nrpages = 2;
    typedef decltype(std::declval<RNG>().operator()()) result_type;
    static constexpr int nelems = pagesize*nrpages/sizeof(result_type);
    result_type* data;
    int pos;
    RNG rng;
};

#endif
