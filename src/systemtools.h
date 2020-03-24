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
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
	
	return buf;
}


