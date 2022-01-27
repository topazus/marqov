/* MIT License

Copyright (c) 2022 Florian Goth
fgoth@physik.uni-wuerzburg.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#ifndef FLOG_IMPL_H
#define FLOG_IMPL_H
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "flog_sinks.h"

/** The singleton that does the access into the log objects.
 * Adapted from one on stackoverflow
 */
class FLog
{
public:
    /** Call to get an actual instance of the object
     */
    static FLog& getInstance()
    {
        static FLog instance;
        return instance;
    }
    /**
     * Write to the log buffers.
     */
    template <int level>
    void write (const std::stringstream& msg) {if (level <= loglevel) {for (int i = 0; i < loggers.size(); ++i) loggers[i]->write(msg);}}
    
    /** Set the runtime loglevel
     * @param level the runtime loglevel
     * @return a reference to the object, for chaining calls.
     */
    FLog& setlevel(int l) {loglevel = l; return *this;}
    /** Delete all registered loggers
     * @return a reference to the object, for chaining calls.
     */
    FLog& clearLoggers() {for (int i = 0; i < loggers.size(); ++i) delete loggers[i]; loggers.clear();return *this;}
    /** Create and Append a logger.
     * 
     * We create the logger ourselves and hence are able to manage the lifetime of
     * all loggers.
     * 
     * @tparam T the type of the logger
     * @tparam Args The Argument Types for the constructor call of the logger
     * 
     * @param args the arguments for the constructor call
     * @return a reference to the object, for chaining calls.
     */
    template <typename T, typename... Args>
    FLog& appendLogger(Args&&... args){loggers.emplace_back(new T(std::forward<Args>(args)...));return *this;}
    /** Drop the last added logger
     * 
     * @return a reference to the object, for chaining calls.
     */
    FLog& popLogger(){delete loggers.back(); loggers.pop_back();return *this;}
    
    /** The number of registered loggers
     * 
     * @return the number of registered loggers
     */
    std::size_t nrofloggers() const noexcept {return loggers.size();}
    FLog(const FLog&) = delete;
    FLog(FLog&&) = delete;
    FLog& operator=(const FLog& ) = delete;
    FLog& operator=(FLog&& ) = delete;
private:
    /** Construct a logger in the default state:
     * A file backed log file that writes to flog.txt
     * and a debug level of DEBUG
     */
    FLog() {
         loggers.push_back(new FLogWriter_file("flog.txt"));
    } //I can still call my own constructor and destructor...*/
    ~FLog(){clearLoggers();}
    std::vector<FLogWriterBase*> loggers; ///< the storage wit hthe polymorphic pointers to all loggers
    int loglevel{DEBUG}; ///< the current loglevel
};

template <int c>
struct ChannelName
{
    static std::string getname();
};

template <>
struct ChannelName<DEBUGVERBOSE>
{
    static std::string getname() {return std::string("DEBUGVERBOSE");}
};

template <>
struct ChannelName<DEBUG>
{
    static std::string getname() {return std::string("DEBUG");}
};

template <>
struct ChannelName<RELEASEVERBOSE>
{
    static std::string getname() {return std::string("RELEASEVERBOSE");}
};

template <>
struct ChannelName<RELEASE>
{
    static std::string getname() {return std::string("RELEASE");}
};

/** @class FLogImpl Log statement collector
 * The function collects a log statement and flushes it out at the end of its lifetime.
 * We rely on constant value propagation to optimize a couple of redundant statements away.
 * @tparam Channel Different channels can be utilized using integers.
 */
template <int Channel=DEBUG>
class FLogImpl
{
public:
    /** Construct an object that collects a log statement.
     */
    FLogImpl() {
        msg << "[" + ChannelName<Channel>::getname() + "] ";
    }
    /** The destructor.
     * This flushes out the data to the actual FLog instance.
     */
    ~FLogImpl()
    {
        if (FLOG_EXTLEVEL >= Channel)
            FLog::getInstance().write<Channel>(msg);
    }
    /**Overload for unary stream operators.
     * 
     * @param osmanip an unary stream manipulator like std::endl.
     * @return an object with changed internal state
     */
    FLogImpl& operator<<(std::ostream&(*const osmanip) (std::ostream&))
    {
        if (FLOG_EXTLEVEL >= Channel)
            msg<<(*osmanip);
        return *this;
    }
    /** input stream operator for character arrays
     */
    FLogImpl& operator<<(const char cs[])
    {
        if (FLOG_EXTLEVEL >= Channel)
            msg<<cs;
        return *this;
    }
    /** input stream operator specialization for strings
     */
    FLogImpl& operator<<(const std::string& s)
    {
        if (FLOG_EXTLEVEL >= Channel)
            msg<<s;
        return *this;
    }
    /** input stream operator
     */
    template <typename T>
    FLogImpl& operator<<(const T& arg)
    {
        if (FLOG_EXTLEVEL >= Channel)
            msg<<arg;
        return *this;
    }
private:
    std::stringstream msg;///< This buffer collects the message until the end of the lifetime of this object    
};

/** A free function to initialize the logger to an empty state an initialize it with a file based logger.
 * @param level the runtime loglevel.
 * @param filename the filename of the log file
 */
void FLogInit(int level, const std::string& filename)
{
    FLog::getInstance().clearLoggers().setlevel(level).appendLogger<FLogWriter_file>(filename);
}

void FLogInit(int level, const char* filename)
{
   FLogInit(level, std::string(filename));
}

void FLogClear(int level)
{
    FLog::getInstance().clearLoggers().setlevel(level);
}

/** Function to register a logger
 * 
 * @tparam T the type of logger.
 * @tparam Args the Argument Types for the constructor call of the logger
 * 
 * @param args the arguments for the construction call of the logger.
 */
template <typename T, typename... Args>
void FLogappendLogger(Args&&... args)
{
    FLog::getInstance().appendLogger<T>(std::forward<Args>(args)...);
}

#endif
