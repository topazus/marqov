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
#ifndef FLOG_SINKS_H
#define FLOG_SINKS_H
#include <iostream>
#include <fstream>
#include <sstream>

class FLogWriterBase
{
public:
    virtual void write(const std::stringstream&) = 0;
    virtual ~FLogWriterBase() = default;
};

class FLogWriter_cout : public FLogWriterBase
{
public:
    void write(const std::stringstream& s) {std::cout<<s.rdbuf();}
};

class FLogWriter_file : public FLogWriterBase
{
public:
    void write(const std::stringstream& s) {file<<s.rdbuf();}
    FLogWriter_file(const std::string& fn) : file(fn) {}
protected:
    std::ofstream file{"flogdebug.txt"};
};

#endif
