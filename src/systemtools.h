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

#include <limits>
//#include <exception>

// make directory
void makeDir(std::string path)
{
	std::string command = "mkdir -p " + path;
	if (system( command.c_str() ) != 0)
	{
		cout << "\nError: Failed to create folder " << path << endl;
//		throw exception();
	}
}


// cross-platform code to get current date/time
// format is YYYY-MM-DD HH:mm:ss

const std::string currentDateTime() 
{
	time_t     now = time(0);
	tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
	
	return buf;
}


