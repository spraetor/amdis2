#include <stdarg.h>		// va_list, va_start, ...
#include <stdio.h>		// vsprintf
#include <iostream>
#include <fstream>
#include <sstream>		// std::stringstream
#include <ostream>		// std::ofstream
#include <set>

#include "Log.h"
#include "AMDiS_base.h"

#ifdef HAVE_PARALLEL_PETSC
#include "petsc.h"		// PetscError
#endif

namespace AMDiS
{
  const char* funcName = NULL;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  bool Msg::outputMainRank = true;
#endif

  ThreadPrivate<const char*> Msg::oldFuncName(NULL);
  std::ostream* Msg::out = NULL;
  std::ostream* Msg::error = NULL;
  int Msg::msgInfo = 10;
  bool Msg::msgWait = true;

  // some macros used internally
#if SUPPRESS_OUTPUT
#define PRINT_LINE(stream, line)
#else
#if HAVE_PARALLEL_DOMAIN_AMDIS
#define PRINT_LINE(stream, line) stream << "[" << MPI::COMM_WORLD.Get_rank() << "] " << line
#else
#define PRINT_LINE(stream, line) stream << line
#endif
#endif

  void Msg::wait(bool w)
  {
    FUNCNAME("Msg::wait()");

    if (w)
    {
      std::string line;
      MSG("wait for <enter> ...");
      std::cin >> line;
    }
  }
  

  void Msg::print_funcname(const char* funcName)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    if (outputMainRank && MPI::COMM_WORLD.Get_rank() != 0)
      return;
#endif

    if (!out)
      out = &std::cout;

    if (funcName &&  oldFuncName.get() != funcName)
    {
      PRINT_LINE((*out), funcName << ":" << std::endl);
    }
    else if (!funcName)
    {
      PRINT_LINE((*out), "*unknown function*" << std::endl);
    }
    PRINT_LINE((*out), "               ");

    oldFuncName.set(funcName);
  }


  void Msg::print_error_funcname(const char* funcName, const char* file, int line)
  {
    static int old_line = -1;

    if (!error)
      error = &std::cerr;

    std::stringstream oss;

    if (funcName && oldFuncName.get() != funcName)
    {
      oss << funcName << ": ";
    }
    else if (!funcName)
    {
      if (line-old_line > 5)
        oss << "*unknown function*";
    }

    if (oldFuncName.get() != funcName)
    {
      oss << "ERROR in " << file << ", line " << line << std::endl;;
      oldFuncName.set(funcName);
    }
    else if (line - old_line > 5)
      oss << "ERROR in " << file << ", line " << line << "\n" << std::endl;

    PRINT_LINE((*error), oss.str());
    old_line = line;
  }


  void Msg::print_error_exit(const char* format, ...)
  {
    va_list arg;
    char buff[255];

    if (!error)
      error = &std::cerr;

    va_start(arg, format);
    vsprintf(buff, format, arg);
    PRINT_LINE((*error), buff);
    va_end(arg);
#if defined HAVE_PARALLEL_DOMAIN_AMDIS && !defined HAVE_PARALLEL_MTL4 && (DEBUG == 0 || defined NDEBUG)
#if (PETSC_VERSION_MINOR >= 5)
    PetscError(MPI_COMM_WORLD, __LINE__, "Msg::print_error_exit", "Global.cc", 1, PETSC_ERROR_INITIAL, buff);
#else
    PetscError(MPI_COMM_WORLD, __LINE__, "Msg::print_error_exit", "Global.cc", "src/", 1, PETSC_ERROR_INITIAL, buff);
#endif
#else
    throw std::runtime_error(buff);
#endif
    exit(1);
  }


  void Msg::print_error(const char* format, ...)
  {
    va_list arg;
    char buff[255];


    if (!error)
      error = &std::cerr;

    va_start(arg, format);
    vsprintf(buff, format, arg);
    PRINT_LINE((*error), buff);
    va_end(arg);
  }


  void Msg::print_warn_funcname(const char* funcName,
                                const char* file,
                                int line)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#if (DEBUG == 0)
    if (outputMainRank && MPI::COMM_WORLD.Get_rank() != 0)
      return;
#endif
#endif

    static int old_line = -1;

    if (!out)
      out = &std::cout;

    std::stringstream oss;

    if (funcName  &&  oldFuncName.get() != funcName)
    {
      oss << funcName << ": ";
    }
    else if (!funcName)
    {
      oss << "*unknown function*";
    }

    if (oldFuncName.get() != funcName)
    {
      oss << "WARNING in " << file << ", line " << line << std::endl;
      oldFuncName.set(funcName);
    }
    else if (line - old_line > 5)
    {
      oss << "WARNING in " << file << ", line " << line << std::endl;
    }

    if (oss.str() != "")
      PRINT_LINE((*out), oss.str());

    old_line = line;
  }


  void Msg::print_warn(const char* format, ...)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#if (DEBUG == 0)
    if (outputMainRank && MPI::COMM_WORLD.Get_rank() != 0)
      return;
#endif
#endif
    va_list arg;
    char buff[255];

    if (!out)
      out = &std::cout;

    va_start(arg, format);
    vsprintf(buff, format, arg);
    PRINT_LINE((*out), buff);
    va_end(arg);
  }


  void Msg::print(const char* format, ...)
  {
#ifndef SUPPRESS_OUTPUT
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    if (outputMainRank && MPI::COMM_WORLD.Get_rank() != 0)
      return;
#endif

    va_list arg;
    char buff[255];

    if (!out)
      out = &std::cout;

    va_start(arg, format);
    vsprintf(buff, format, arg);
    (*out) << buff;
    va_end(arg);
#endif
  }



  void waitSec(int seconds)
  {
    clock_t endwait = clock () + seconds * CLOCKS_PER_SEC;
    while (clock() < endwait) {}
  }


  void processMemUsage(double& vm_usage, double& resident_set, bool inMegaByte)
  {
#ifndef _WIN32
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage     = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss;

    // in case x86-64 is configured to use 2MB pages
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;

    if (inMegaByte)
    {
      vm_usage /= 1024.0;
      resident_set /= 1024.0;
    }
#else
    ERROR("Function not available under MS Windows\n");
#endif
  }


} // end namespace AMDiS
