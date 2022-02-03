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

#ifndef MLOG_H
#define MLOG_H

//use the respective defines of FLog to configure it's behaviour
#define FLOG_EXTLEVEL MLOG_EXTLEVEL
#define FLOG_DISABLE_VERBOSE_FUN
#include "flog/flog.h"
#include <string>
#include <chrono>

class MLogAppender : public FLogWriterBase
{
public:
    void write(const std::stringstream& s) {file<<s.rdbuf();}
    void reopen(const std::string& fn) {if(file.is_open()) file.close(); file.open(fn);}
    MLogAppender(const std::string& fn) : file(fn, std::ios::out | std::ios::app) {}
private:
    std::ofstream file;
};

class MLogState {
public:
    int level{DEBUG};
    std::string fn;
    MLogState(int l, const std::string& f) : level(l), fn(f) {reset();}
    MLogState() = delete ;
    void reset()
    {
        FLogClear(level);
        FLogappendLogger<MLogAppender>(fn);
    }
};

#define MLOGTIME "( "<< std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count()<<" )"

#define MLOGDEBUGVERBOSE FLOGDEBUGVERBOSE<<MLOGTIME
#define MLOGDEBUG FLOGDEBUG<<MLOGTIME
#define MLOGRELEASEVERBOSE FLOGRELEASEVERBOSE<<MLOGTIME
#define MLOGVERBOSE FLOGVERBOSE<<MLOGTIME
#define MLOGRELEASE FLOGRELEASE<<MLOGTIME

#endif
