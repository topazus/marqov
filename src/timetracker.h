#ifndef TIMETRACKER_H
#define TIMETRACKER_H

#include <chrono>
#include <utility>
#include <unordered_map>
#include <unistd.h> // provides usleep, only for testing purposes

typedef std::chrono::milliseconds msec;

class marqovclock
{
	public: 
		std::string name;
		std::chrono::system_clock::time_point starttime, inittime;
		double duration_msec = 0;

		marqovclock(std::string name) : name(name), inittime(std::chrono::system_clock::now())
		{}

};



class timetracker
{
	public:

		std::vector<marqovclock> clocks;
		std::unordered_map<std::string, int> clockmap; // recreate the functionality of a Python dictionary

		std::string active_clock;

		timetracker(){}

		void add_clock(std::string name)
		{
			auto mapentry = std::pair<std::string,int>(name,clocks.size());
			clockmap.insert(mapentry);
			clocks.push_back(marqovclock(name));
		}

		void status()
		{
			cout << endl;
			for (auto& x: clocks) 
			{
				std::cout << x.name << ": " << x.duration_msec;
				if (x.name == active_clock) cout << "\t (active)";
				std::cout << endl;
			}
		}

		void run(std::string name)
		{
			active_clock = name;
			clocks[clockmap[active_clock]].starttime = std::chrono::system_clock::now();
		}

		void switch_clock(std::string target)
		{
			std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

			auto clockidx = clockmap[active_clock];
			auto& previous_clock = clocks[clockidx];

			previous_clock.duration_msec += std::chrono::duration_cast<msec>(now-previous_clock.starttime).count();

			auto& target_clock = clocks[clockmap[target]];
			target_clock.starttime = now; 

			active_clock = target;
		}
			
};

#endif
