#include <vector>
#include "metropolis.h"

using namespace std;

/** This test compares the acceptance rates of the optimized Metropolis rule in MARQOV
 * to a plain implementationl
 */

int acc = 0;

    template <class RNGType>
    inline bool update_accepted(double dE, double beta, RNGCache<RNGType>& rng)
    {
        bool accept = false;
        if ( dE <= 0 )
		{
            accept = true;
		}
		else
        {
//             std::cout<<"here"<<std::endl;
                if (rng.real() < std::exp(-beta*dE) /*1.0 - beta*dE + beta*dE*beta*dE/2 -beta*dE*beta*dE*beta*dE/6*/) 
                {
                    accept = true;
                }
        }
		return accept;
    }

int main()
{
    constexpr int N = 100000;
    std::vector<double> dE(N);
    //disable one of the two scopes to measure the behaviour of one implementation.
    
    {
    RNGCache<mt19937_64> rngcache{0};///< The caching RNG
    for(auto& e: dE)
        e = 10.0*rngcache.real();
    for(int i = 0; i < 500; ++i)
    {
    for(int i = 0; i < N; ++i)
        if(MARQOV::update_accepted(dE[i], 0.44, rngcache)) acc++;
    }
    }
    int acc2 = acc;
//     std::cout<<static_cast<double>(acc)/N/500<<std::endl;
    acc = 0;
    
        {
    RNGCache<mt19937_64> rngcache{0};///< The caching RNG
    for(auto& e: dE)
        e = 10.0*rngcache.real();
    for(int i = 0; i < 500; ++i)
    {
    for(int i = 0; i < N; ++i)
        if(update_accepted(dE[i], 0.44, rngcache)) acc++;
    }
    }
//     std::cout<<static_cast<double>(acc)/N/500<<std::endl;
// std::cout<<acc2-acc<<std::endl;
    return acc2 - acc;
}
