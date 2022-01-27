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

#ifndef FLOG_MACROS_H
#define FLOG_MACROS_H
#include <string>

#define DEBUGVERBOSE 3
#define DEBUG 2
#define RELEASEVERBOSE 1
#define RELEASE 0

#ifndef FLOG_EXTLEVEL
    #define FLOG_EXTLEVEL RELEASE
#endif

#if defined(__GNUC__) || defined(__clang__)
#define FLOG_FUN __PRETTY_FUNCTION__
#else
#define FLOG_FUN __func__
#endif

#define FLOGPREFIX "("<<std::string(__FILE__)<<", "<<std::string(FLOG_FUN)<<", "<<std::to_string(__LINE__)<<std::string(") ")

#define FLOGDEBUGVERBOSE (FLogImpl<3>()<<FLOGPREFIX)
#define FLOGDEBUG (FLogImpl<2>()<<FLOGPREFIX)
#define FLOGRELEASEVERBOSE (FLogImpl<1>()<<FLOGPREFIX)
#define FLOGVERBOSE FLOGRELEASEVERBOSE
#define FLOGRELEASE (FLogImpl<0>()<<FLOGPREFIX)

#define DP(x) FLOGDEBUG<<#x<<" = "<<(x)<<'\n'
#endif
