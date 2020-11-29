#include "cachecontainer.h"
#include <chrono>
#include <string>
#include <unistd.h> // provides usleep

//the following macro seems to be necessary for dumping string defines...
#define XSTR(x) STR(x)
#define STR(x) #x

namespace MARQOV
{
    void dumpEnvironmenttoHDF5Group(H5::Group& h5loc)
    {
        h5loc.setComment("Here we have various parameters of the host system.");
        dumpscalartoH5(h5loc, std::string("Code"), std::string("MARQOV"));
        dumpscalartoH5(h5loc, std::string("Version"), std::string("branch: ") + std::string(XSTR(GIT_BRANCH)) + ", SHA1: " +std::string( XSTR(GIT_SHA1) ));
        dumpscalartoH5(h5loc, std::string("Website"), std::string("marqov.physik.uni-wuerzburg.de"));
        dumpscalartoH5(h5loc, std::string("E-Mail"), std::string("marqov@physik.uni-wuerzburg.de"));
        
        //Acquire hostname
        char hostname[HOST_NAME_MAX];
        int err = gethostname(hostname, HOST_NAME_MAX);
        if (err == 0) // We were able to retrieve a hostname, hence we dump it
        {
            dumpscalartoH5(h5loc, std::string("hostname"), std::string(hostname));
        }
        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        std::time_t start_time = std::chrono::system_clock::to_time_t(now);
        char timestamp[100];
        struct tm buf;
        buf = *(std::localtime(&start_time));
        std::strftime(timestamp, sizeof(timestamp), "%A %Y-%m-%d %H:%M:%S", &buf);
        dumpscalartoH5(h5loc, std::string("startingdate"), std::string(timestamp));
    }
};
