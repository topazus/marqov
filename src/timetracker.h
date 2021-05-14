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

typedef musec timeformat; // internal time resolution


namespace marqovtime
{


	class marqovclock
	{
		public: 
			std::string tag, name;
			std::chrono::system_clock::time_point starttime;
			std::chrono::system_clock::time_point inittime;
			decltype(std::chrono::duration_cast<timeformat>(Time::now()-Time::now())) dur = timeformat::zero(); 
			marqovclock(std::string tag, std::string name) : tag(tag), name(name), inittime(Time::now()) {}
			marqovclock(std::string name) : name(name), inittime(Time::now()) { tag = name.at(0);}
	};
	
	
	
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



			void status(bool minimal=false, bool verbose=true)
			{
				const auto now = Time::now();

				// sum up contributions
				double sum = 0;
				for (auto& x: clocks) sum = sum + std::chrono::duration_cast<timeformat>(x.dur).count();


				std::cout << std::setprecision(1) << std::fixed;

				// find suitable time format 
				double print_mult = 1.0;
				std::string print_unit = "mus";

				const double kilo = 1000;

				if (sum > 10*kilo)
				{
					print_mult = 1./kilo;
					print_unit = "ms";
				}
				if (sum > 10*kilo*kilo)
				{
					print_mult = 1./kilo/kilo;
					print_unit = "s";
				}
				if (sum > 10*kilo*kilo*60)
				{
					print_mult = 1./kilo/kilo/60.;
					print_unit = "min";
				}




				if (minimal)
				{

				}
				else
				{
	
					for (auto& x: clocks) 
					{
						auto dur_print = std::chrono::duration_cast<timeformat>(x.dur).count();
						sum = sum + dur_print;
	
						if (x.name != active_clock) // list all clocks
						{
							std::cout << "  " << x.name << "\t" << dur_print*print_mult << print_unit << endl;
						}
						else
						{
							auto diff = std::chrono::duration_cast<timeformat>(now-x.starttime).count();
							sum = sum + diff;
							std::cout << x.name << " " << dur_print+diff << "\t (active)";
						}
					}
	
				
					if (verbose) // print wallclock
					{
//						wallclock.dur = std::chrono::duration_cast<timeformat>(now-wallclock.inittime);
//						auto dur_print = std::chrono::duration_cast<timeformat>(wallclock.dur).count();
//						std::cout << "  wallclock     " << dur_print << " msec" << endl;
					}
					cout << endl << endl;
				}
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
