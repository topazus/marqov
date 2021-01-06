/***************************************************************************
 *   Copyright (C) 2008 - 2017 by Florian Goth   *
 *   fgoth@physik.uni-wuerzburg.de   *
 *                                                                         *
 *   All rights reserved.                                                  *
 *                                                                         *
 *   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: *
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. *
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. *
 *     * Neither the name of Florian Goth nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. *
 *                                                                         *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR *
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR *
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, *
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,   *
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF *
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          *
 ***************************************************************************/
#ifndef PRNG_H
#define PRNG_H

#include <limits>

template <class Int, class Float>
class Prng
{
public:
    typedef typename Float::FPType FPType;
    /** generates a random number in [0...1]
    @return a random number in the range [0...1]
    */
    inline FPType rndfloat();
    /**
    generates a float in the specified range
    @param rng the range in which to generate the float
    @return A float in the specified range
    */
    inline FPType rndfloat(FPType rng);
    /**A random integer as returned by the prng. No guarantees are made.
    @return a random integer
    */
    inline int rndint();
    /**A random integer in the range [0...rng)
    @param rng The maximum Integer to return
    */
    inline unsigned int rndint(unsigned int rng);
    /**
    A random Integer in the range specified by the template Parameter rng
    */
    template <unsigned int rng> inline unsigned int rndInteger();
    /**
    Constructor with two seed values.
    @param xi the seed value for the Integer PRNG
    @param xf the seed value for the floating point PRNG
    */
    inline Prng(unsigned int xi, int xf) : fp_prng(xf), int_prng(xi) {}
    inline void reinit(unsigned int xi, FPType xf);
    Float fp_prng;///< have an instance of PRNG used for generating floating point numbers
    Int int_prng;///< an instance of the PRNG used for generating integers
private:
};

template <class Int, class Float>
template <unsigned int rng>
unsigned int Prng<Int, Float>::rndInteger()
{
    return int_prng.rnd() % rng;//note that we are restricting to the least significant bits. they are often not that random
}

template <class Int, class Float>
unsigned int Prng<Int, Float>::rndint(unsigned int rng)
{
    return int_prng.rnd() % rng;//note that we are restricting to the least significant bits. they are often not that random
}

template <class Int, class Float>
void Prng<Int, Float>::reinit(unsigned int xi, FPType xf)
{
int_prng.reinit(xi);
fp_prng.reinit(xf);
}

template <class Int, class Float>
typename Prng<Int, Float>::FPType Prng<Int, Float>::rndfloat(FPType rng)
{
    return fp_prng.rnd() * rng;
}
template <class Int, class Float>
typename Prng<Int, Float>::FPType Prng<Int, Float>::rndfloat()
{
    return fp_prng.rnd();
}

template <class Int, class Float>
int Prng<Int, Float>::rndint()
{
    return int_prng.rnd();
}

template <class IntegerPRNG, typename FPTypeParam, unsigned int MAX = std::numeric_limits<unsigned int>::max() >
class ToFloat : public IntegerPRNG
{
//Helper structure to convert any PRNG that outputs Integer to floats in the range [0...1]  . The range of the
//IntegerPRNG is specified using the MAX template parameter
public:
    typedef FPTypeParam FPType;
    inline FPType rnd();
    inline ToFloat(int x) : IntegerPRNG(x) {}
    inline void reinit(int x);
private:
};

template <class IntegerPRNG, typename FPTypeParam, unsigned int MAX>
void ToFloat<IntegerPRNG, FPTypeParam, MAX>::reinit(int x)
{
IntegerPRNG::reinit(x);
}

template <class IntegerPRNG, typename FPTypeParam, unsigned int MAX>
FPTypeParam ToFloat<IntegerPRNG, FPTypeParam, MAX>::rnd()
{
    return static_cast<float> (IntegerPRNG::rnd()) / static_cast<float>(MAX);
}

template<unsigned int a, unsigned int c, unsigned int m>
class LinearCongruential
{
    //A Linear Congruential Generator
public:
    inline unsigned int rnd();
    LinearCongruential(unsigned int seed) : x(seed) {}
    inline void reinit(unsigned int);
private:
    unsigned int x;///< the value that is stored in between the generation of a random number
};

template<unsigned int a, unsigned int c, unsigned int m>
void LinearCongruential<a,c,m>::reinit(unsigned int s)
{
x = s;
}

template<unsigned int a, unsigned int c, unsigned int m>
unsigned int LinearCongruential<a,c,m>::rnd()
{
    return (x =(a * x + c) % m);// the familiar scheme for generating random numbers with a linear congruential one
}

typedef LinearCongruential<65539, 0, 2147483647u> RANDU;//A typedef for the famous RANDU PRNG
typedef LinearCongruential<62089911, 4349, 2147483647u> LCG;


#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)
#endif
#if (GCC_VERSION >= 40300) || (__cplusplus > 199711L)

#if (__cplusplus > 199711L)
#include <random>
class MersenneTwister
{
public:
inline std::mt19937::result_type rnd()
{
return rng();
}
inline MersenneTwister(unsigned int seed) : rng(seed) {}
inline void reinit(unsigned int x) {rng.seed(x);}
private:
std::mt19937 rng;
};
#else
#include <tr1/random>
class MersenneTwister
{
public:
inline std::tr1::mt19937::result_type rnd()
{
return rng();
}
inline MersenneTwister(unsigned int seed) : rng(seed) {}
inline void reinit(unsigned int x) {rng.seed(x);}
private:
std::tr1::mt19937 rng;
};
#endif
#endif

/**
This class encapsulates the RNG by Robert M. Ziff. for further info see his paper or Haye Hinrichsens lecture notes
*/
class Ziff
{
public:
    inline Ziff(unsigned int);
    inline void reinit(unsigned int);
    inline unsigned int rnd();
private:
    enum  //various constants. we store them in an enum
    {
        Random_A = 471,
        Random_B = 1586,
        Random_C = 6988,
        Random_D = 9689,
        Random_M = 16383,
    };
    const unsigned int Random_Max;///< the Maximum number that can be generated
    int Random_nd;
    int Random_ra[Random_M+1];///<the array of numbers the RNG is working on
};

void Ziff::reinit(unsigned int seed)
{
    Random_nd = 0;//forgotten in Hinrichsens Code ?
// we initialize via a "good" LCG random number generator
    LCG lcg(seed);
    for (unsigned int k = 0; k <= Random_M; ++k)//some warmup
        Random_ra[k] = lcg.rnd();
    //some warmup
    for (unsigned int k = 0; k < 100; ++k) rnd();
}

unsigned int Ziff::rnd()
{
    ++Random_nd;
    return (Random_ra[Random_nd & Random_M] =
                Random_ra[(Random_nd - Random_A) & Random_M] ^
                Random_ra[(Random_nd - Random_B) & Random_M] ^
                Random_ra[(Random_nd - Random_C) & Random_M] ^
                Random_ra[(Random_nd - Random_D) & Random_M]);
}

Ziff::Ziff(unsigned int seed) : Random_Max(2147483648U)
{
    reinit(seed);
    return;
}

/**
the RNG by Marsaglia and Zaman. Chi^2 show it being a bit bad... 
*/
class Marsaglia_Zaman
{
public:
    inline unsigned int rnd();
    inline Marsaglia_Zaman(unsigned int seed);
    inline void reinit(unsigned int);
private:
    inline unsigned int c();
    unsigned int cp0;
    unsigned int cp1;
    unsigned int r[4];
};

void Marsaglia_Zaman::reinit(unsigned int seed)
{
// we initialize via a "good" LCG random number generator
    LCG lcg(seed);
    r[0] = lcg.rnd();
    r[1] = lcg.rnd();
    r[2] = lcg.rnd();
    r[3] = lcg.rnd();
    cp0 = lcg.rnd();
    cp1 = 0;
//some warmup
    for (unsigned int k = 0; k < 100; ++k) rnd();
}

Marsaglia_Zaman::Marsaglia_Zaman(unsigned int seed)
{
    reinit(seed);
}

unsigned int Marsaglia_Zaman::rnd()
{
    r[3] = r[2];
    r[2] = r[1];
    r[1] = r[0];
    r[0] = (r[2] - r[3] - c()) % (4294967278ul);
    return (69069 * r[1] + 1013904243);//Assume implicit modulo operation
}

unsigned int Marsaglia_Zaman::c()
{
    cp0 = cp1;
    cp1 = ((r[2] - r[3] - cp0) > 0 ? 1 : 0);
    return cp0;
}

template<unsigned int m>
class Coveyou
{
    //Coveyous RNG
    //m should be a power of two?
    //be careful with this...
public:
    inline unsigned int rnd();
    Coveyou(unsigned int seed)
    {
       reinit(seed);
    }
    inline void reinit(unsigned int);
private:
    unsigned int x;
};

template<unsigned int m>
inline void Coveyou<m>::reinit(unsigned int seed)
{
        while ((seed % 4) != 2) ++seed;
        x = seed;
}

template<unsigned int m>
inline unsigned int Coveyou<m>::rnd()
{
    return (x = (x * (x + 1)) % m);
}

/**
this class encapsulates the inverse Congruential Random Number Generator
*/
template<unsigned int a, unsigned int c, unsigned int m>
class InverseCongruential
{
    //The Inverse Congruential Generator
public:
    inline unsigned int rnd();
    InverseCongruential(unsigned int seed) : x(seed) {}
    inline void reinit(unsigned int seed) {x = seed;}
private:
    unsigned int x;
};

template <typename T>
class Tripel
{
    //Helper Structure for the inverse Modulo Operation
public:
    inline T& operator[](unsigned int k)
    {
        return p[k];
    }
    inline Tripel<T>& operator=(const Tripel<T>& rhs)
    {
        if (&rhs != this)
        {
            p[0] = rhs.p[0];
            p[1] = rhs.p[1];
            p[2] = rhs.p[2];
        }
        return *this;
    }
    Tripel(T a1, T a2, T a3)
    {
        p[0]  = a1;
        p[1]  = a2;
        p[2]  = a3;
    }
private:
    T p[3];
};

template <typename T>
inline T inverseModulo(const T a, const T n)
{
    Tripel<T> g0(n,1,0), g1(a,0,1), t(0,0,0);
    while (g1[0] != 0)
    {
        T y = g0[0] / g1[0];
        t[0] = g0[0] % g1[0];
        t[1] = g0[1] - y * g1[1];
        t[2] = g0[2] - y * g1[2];
        g0 = g1;
        g1 = t;
    }
    T retval = g0[2] % n;
    if (retval < 0) retval = (retval + n) % n;
    return retval;
}

template<unsigned int a, unsigned int c, unsigned int m>
inline unsigned int InverseCongruential<a,c,m>::rnd()
{
    unsigned int Xq = (x == 0 ? 0 : inverseModulo<long int>(x,m) );
    return (x = (a * Xq + c) % m);
}

#endif
