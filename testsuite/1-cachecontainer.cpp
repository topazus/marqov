#include "cachecontainer.h"
#include <vector>

using namespace std;

/** This test tests proper behaviour of the cache for vector valued time series.
 */

int main()
{
    auto obs = std::vector<double>{1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    H5::H5File file("test.h5", H5F_ACC_TRUNC);
    H5::Group step(file.createGroup("obs"));
    
    CacheContainer<decltype(obs)> obscache(step, "name", "desc");
    for(int i = 0; i < 500000; ++i)
        obscache<<obs;
}
