#include <vector>
#include "metropolis.h"

using namespace std;

/** This test compares the acceptance rates of the optimized Metropolis rule in MARQOV
 * to a plain implementation.
 * The defines below control wether the current MARQOV Implementation (NEW)
 * or the baseline (OLD) is executed.
 * In the test we execute both and we determine if the same number of moves is accepted.
 */

#define OLD
#define NEW

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
    int acc = 0;
    constexpr int N = 100000;
    std::vector<double> dE(N);
    //disable one of the two scopes to measure the behaviour of one implementation.
#ifdef NEW
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
#endif
    int acc2 = acc;
//    std::cout<<static_cast<double>(acc)/N/500<<std::endl;
    acc = 0;
    
#ifdef OLD
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
#endif
//     std::cout<<static_cast<double>(acc)/N/500<<std::endl;
// std::cout<<acc2-acc<<std::endl;
    return acc2 - acc;
}
