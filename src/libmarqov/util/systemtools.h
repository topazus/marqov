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

#include <sys/stat.h>
#include <string>
#include <exception>
#include <cerrno>
#include <fstream>

/** Check if file exists
 *
 * This checks wether the file can be opened
 * @param fn the filename
 * @return true if it can be opened, else false
 */
inline bool fileexists(const std::string& fn)
{
    return std::ifstream(fn).good();
}

/** Create Directory.
 * 
 * This creates a directory. If the path exists nothing happens.
 * @param path the path of the directory that should be created.
 * @throws std::runtime_error If an unrecoverable error occurs. If the path exists nothing happens.
 */
inline void makeDir(const std::string path)
{
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status != 0)
    {
        if (errno != EEXIST)
            throw std::runtime_error(std::string("[MARQOV] Failed to create folder ") + path);
    }
}
