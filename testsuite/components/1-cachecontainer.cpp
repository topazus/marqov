#include <string>
#include "cachecontainer.h"

using namespace std;

int main()
{
    {
    std::string filename("test.h5");
    auto flag = H5F_ACC_TRUNC;
    H5::H5File retval(filename, flag);
    CacheContainer<double> testcont(retval, "testobs", "a test description");
    CacheContainer<int> testint(retval, "testint", "a test int description");
    for (uint i = 0; i < 2000000; ++i)
    {
        testcont<<42.0;
        testint<<42;
    }
    }
}
