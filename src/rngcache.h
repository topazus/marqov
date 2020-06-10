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
#include <stdexcept>

template <class RNG>
class RNGCache
{
public:
    template <class ...Args>
    RNGCache(Args&&... args) : data(NULL), pos(0), rng(std::forward<Args>(args)...)
    {
        int err = posix_memalign( (void**)(&data), pagesize, nrpages*pagesize );//in C++17 we can replace that with aligned_alloc
        if (err != 0)
            throw std::runtime_error("[RNGCache] error while allocating memory");
        fillcache();
    }
    ~RNGCache()
    {
        free(data);
    }
    auto number() noexcept
    {
        if(pos >= nelems)
        {
            pos = 0;
            fillcache();
        }
        return data[pos++];
    }
    auto real(double max = 1.0, double min = 0.0) noexcept
    {
        return min + number()*(max-min)/RNG::max();
    }
    void dump()
    {//FIXME!!
    }
    static constexpr auto max() {return RNG::max();}
private:
    void fillcache() noexcept
    {
        for(int i = 0; i < nelems; ++i)
            data[i] = rng(); //We follow C++11 convention that operator() advances the state of the RNG
    }
    static constexpr int pagesize = 4096;
    static constexpr int nrpages = 2;
    typedef decltype(std::declval<RNG>().operator()()) RNGValType;
    static constexpr int nelems = pagesize*nrpages/8;
    RNGValType* data;
    int pos;
    RNG rng;
};

#endif
