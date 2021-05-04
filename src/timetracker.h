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

#ifndef TIMETRACKER_H
#define TIMETRACKER_H

#include <chrono>
#include <vector>
#include <utility>
#include <unordered_map>
#include <exception>
#include <unistd.h> // provides usleep, only for testing purposes

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::seconds sec;
typedef std::chrono::milliseconds msec;
typedef std::chrono::microseconds musec;

typedef musec timeformat; // internal time resolution
typedef msec printformat; // time format used for output


namespace marqovtime
{

	/** 
	 * The marqovclock
	 */
	class marqovclock
	{
		public: 
			std::string name;
			std::chrono::system_clock::time_point starttime;
			std::chrono::system_clock::time_point inittime;
			decltype(timeformat::zero()) dur = timeformat::zero();
			marqovclock(std::string name) : name(name), inittime(Time::now()) {}
	
	};
	

	/**
	 * The time tracker
	 */
	class timetracker
	{
		public:
	
			marqovclock wallclock = marqovclock("wallclock"); // the reference 
			std::vector<marqovclock> clocks;	// list of all clocks
			std::unordered_map<std::string, int> clockmap; // recreate the functionality of a Python dictionary
			std::string active_clock = "None"; // there can only be one active clock at a time
	
			timetracker()
			{
				wallclock.starttime = Time::now();
			}
	
			void add_clock(std::string name)
			{
				auto mapentry = std::pair<std::string,int>(name,clocks.size());
				clockmap.insert(mapentry);
				clocks.push_back(marqovclock(name));
			}
	
			void status(bool verbose=true)
			{
				const auto now = Time::now();
				double sum = 0;
	
				cout << endl << endl;
				for (auto& x: clocks) 
				{
					auto dur_print = std::chrono::duration_cast<printformat>(x.dur).count();
					sum = sum + dur_print;
	
					if (x.name != active_clock) // list all clocks
					{
						std::cout << x.name << ": " << dur_print;
					}
					else
					{
						auto diff = std::chrono::duration_cast<printformat>(now-x.starttime).count();
						sum = sum + diff;
						std::cout << x.name << ": " << dur_print+diff << "\t (active)";
					}
					std::cout << endl;
				}
	
				
				if (verbose) // print und wallclock
				{
					std::cout << "-----" << endl << "sum: " << sum << endl;
					wallclock.dur = std::chrono::duration_cast<timeformat>(now-wallclock.inittime);
					auto dur_print = std::chrono::duration_cast<printformat>(wallclock.dur).count();
					std::cout << "wallclock: " << dur_print << endl;
				}
				cout << endl << endl;
			}
	
	
			void run(std::string name)
			{
				if (active_clock != "None") 
				{
					throw std::invalid_argument("There is already a clock running; This may result in unwanted behaviour; use the switch function...."); // todo: catch this exception
				}
				active_clock = name;
				clocks[clockmap[active_clock]].starttime = Time::now();
			}
	
	
			// stop the active clock
			void stop()
			{
				auto now = Time::now();
				auto clockidx = clockmap[active_clock];
				auto& previous_clock = clocks[clockidx];
				previous_clock.dur += std::chrono::duration_cast<timeformat>(now-previous_clock.starttime);
				active_clock = "None";
			}
	
	
			void switch_clock(std::string target)
			{
				auto now = Time::now();
	
	
				auto clockidx = clockmap[active_clock];
				auto& previous_clock = clocks[clockidx];
	
				previous_clock.dur += std::chrono::duration_cast<timeformat>(now-previous_clock.starttime);
	
				auto& target_clock = clocks[clockmap[target]];
				target_clock.starttime = now; 
	
				active_clock = target;
			}
				
	};

}

#endif
