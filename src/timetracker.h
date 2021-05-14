#ifndef TIMETRACKER_H
#define TIMETRACKER_H

#include <chrono>
#include <utility>
#include <unordered_map>
#include <exception>
#include <unistd.h> // provides usleep, only for testing purposes

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::seconds sec;
typedef std::chrono::milliseconds msec;
typedef std::chrono::microseconds musec;

typedef musec timeformat; ///< internal time resolution, do not change


namespace marqovtime
{
	/** A clock which can be used to measure certain elements of the code
     */
	class marqovclock
	{
		public: 
			std::string tag, name;
			std::chrono::system_clock::time_point starttime;
			std::chrono::system_clock::time_point inittime;
			decltype(std::chrono::duration_cast<timeformat>(Time::now()-Time::now())) dur = timeformat::zero(); 

			/** Construct a clock object
			*
			* @param tag short name for this clock (should only be very few characters)
			* @param name longer clock name or description
			*/
			marqovclock(std::string tag, std::string name) : tag(tag), name(name), inittime(Time::now()) {}
			marqovclock(std::string name) : name(name), inittime(Time::now()) { tag = name.at(0);}
	};
	
	
	
	/** The timetracker class which manages the clocks
	*/
	class timetracker
	{
		public:
	
			marqovclock wallclock = marqovclock("wallclock"); // the reference 
			std::vector<marqovclock> clocks;	// list of all clocks
			std::unordered_map<std::string, int> clockmap; // recreate the functionality of a Python dictionary
			std::string active_clock = "None"; // there can only be one active clock at a time
	
			/* Construct a timetracker object
			*/
			timetracker()
			{
				wallclock.starttime = Time::now();
			}
	

			/** Add clock to the timetracker
			*
			* @param tag short name for this clock (should only be very few characters)
			* @param name longer clock name or description
			*/
			void add_clock(std::string name)
			{
				auto mapentry = std::pair<std::string,int>(name,clocks.size());
				clockmap.insert(mapentry);
				clocks.push_back(marqovclock(name));
			}
			void add_clock(std::string tag, std::string name)
			{
				auto mapentry = std::pair<std::string,int>(name,clocks.size());
				clockmap.insert(mapentry);
				clocks.push_back(marqovclock(tag, name));
			}



			/** Print current durations of all clocks
			*
			* @param minimal if true (default), timing results are printed as a one-liner
			* @param verbose if true (default), the sum of all clock timings is printed as well
			*/
			void status(bool minimal=true, bool verbose=true)
			{
				const auto now = Time::now();

				// sum up contributions
				double sum = 0;
				for (auto& x: clocks) sum = sum + std::chrono::duration_cast<timeformat>(x.dur).count();

				// find suitable time format 
				double print_mult = 1.0;
				std::string print_unit = "mus";
				std::cout << std::setprecision(1) << std::fixed;

				const double kilo = 1000;

				if (sum > 10*kilo) // milliseconds
				{
					print_mult = 1./kilo;
					print_unit = "ms";
				}
				if (sum > 10*kilo*kilo) // seconds
				{
					print_mult = 1./kilo/kilo;
					print_unit = "s";
				}
				if (sum > 10*kilo*kilo*60) // minutes
				{
					print_mult = 1./kilo/kilo/60.;
					print_unit = "min";
				}
				if (sum > 10*kilo*kilo*60*60) // hours
				{
					print_mult = 1./kilo/kilo/3600.;
					print_unit = "h";
				}


				// print timing results as a one-liner using only the clock tags
				if (minimal)
				{
					std::cout << "Timing results (in " << print_unit << "): ";
					for (auto& c: clocks)
					{
						auto dur_print = std::chrono::duration_cast<timeformat>(c.dur).count();
						std::cout << c.tag << " " << dur_print*print_mult << " | ";
					}
				}
				// print timing results in a more verbose way using the full clock names
				else
				{
					std::cout << "Timing results:" << endl;
					for (auto& x: clocks) 
					{
						auto dur_print = std::chrono::duration_cast<timeformat>(x.dur).count();
						sum = sum + dur_print;
	
						if (x.name != active_clock) // list all clocks
						{
							std::cout << "  " << x.name << "\t" << dur_print*print_mult << " " << print_unit << endl;
						}
						else
						{
							auto diff = std::chrono::duration_cast<timeformat>(now-x.starttime).count();
							sum = sum + diff;
							std::cout << x.name << " " << dur_print+diff << "\t (active)";
						}
					}
				}
	
				
				if (verbose) // print wallclock
				{
					wallclock.dur = std::chrono::duration_cast<timeformat>(now-wallclock.inittime);
					auto dur_print = std::chrono::duration_cast<timeformat>(wallclock.dur).count();

					if (minimal) std::cout << "wall " << dur_print*print_mult << std::endl;
					else std::cout << "  wallclock     " << dur_print*print_mult << " " << print_unit << endl << endl;
				}
		}
	
	
			
		/** Run a clock 
		* only one clock can be running simultaneously
		*
		* @param name specify name of the clock the be run
		*/
		void run(std::string name)
		{
			if (active_clock != "None") 
			{
				throw std::invalid_argument("There is already a clock running; This may result in unwanted behaviour; use the switch function...."); // todo: catch this exception
			}
			active_clock = name;
			clocks[clockmap[active_clock]].starttime = Time::now();
		}
	
	
		/** Stop the currently active clock */
		void stop()
		{
			auto now = Time::now();
			auto clockidx = clockmap[active_clock];
			auto& previous_clock = clocks[clockidx];
			previous_clock.dur += std::chrono::duration_cast<timeformat>(now-previous_clock.starttime);
			active_clock = "None";
		}
	
	
		/** Switch clocks, that menas stop current clock and start another one
		*
		* @param target name of the clock to be switched to
		*/
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
