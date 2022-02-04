/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2022, The MARQOV Project
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

#ifndef QUERY_THREADS_H
#define QUERY_THREADS_H

#include <string>
#include "registry.h"

int number_of_threads_per_node(RegistryDB& registry, const std::string& configfile)
{
	int nthreads = 0;
	try
	{
		nthreads = registry.template Get<int>(configfile, "General", "threads_per_node" );
	}
	catch (const Registry_Key_not_found_Exception&)
	{
		std::cout<<"threads_per_node not set -> automatic"<<std::endl;
	}
	return nthreads;
}
#endif
