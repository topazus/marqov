#include <string>
#include "cachecontainer.h"

using namespace std;

int main()
{
    {
    std::string filename("test.h5");
    auto flag = H5F_ACC_TRUNC;
    H5::H5File file(filename, flag);
    H5::Group obs1(file.createGroup("obs1"));
    H5::Group obs2(file.createGroup("obs2"));

    CacheContainer<double> testcont(obs1, "testobs", "a test description");
    CacheContainer<int> testint(obs2, "testint", "a test int description");
    for (uint i = 0; i < 2000000; ++i)
    {
        testcont<<42.0;
        testint<<42;
    }
    }
}
