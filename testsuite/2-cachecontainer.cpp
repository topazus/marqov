#include <vector>
#include <iostream>

#include "cachecontainer.h"

using namespace std;

/** This test tests proper behaviour of the cache for vector valued time
 * series if the the time series is bigger than the default cache size.
 */

int main()
{
    auto obs = std::vector<double>(std::size_t{600000}, 2.0);
    H5::H5File file("test.h5", H5F_ACC_TRUNC);
    H5::Group step(file.createGroup("obs"));
    
    CacheContainer<decltype(obs)> obscache(step, "name", "desc");
    for(int i = 0; i < 5; ++i)
        obscache<<obs;
}
