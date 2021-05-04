/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
