/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include <stdarg.h>
#include <stdio.h>
#include <sstream>

#include "Global.h"
#include "Initfile.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#ifdef HAVE_PARALLEL_PETSC
  #include "petsc.h"
#endif
namespace AMDiS {

  const char *funcName = NULL;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  bool Msg::outputMainRank = true;
#endif

  ThreadPrivate<const char *> Msg::oldFuncName(NULL);
  std::ostream* Msg::out = NULL;
  std::ostream* Msg::error = NULL;
  int Global::dimOfWorld = 0;
  std::vector<std::vector<int> > Global::geoIndexTable;
  int Msg::msgInfo = 10;
  bool Msg::msgWait = true;

  Element *Global::referenceElement[4] = 
    { NULL, 
      new Line(NULL), 
      new Triangle(NULL), 
      new Tetrahedron(NULL) 
    };
    

  void Msg::wait(bool w)
  {
    FUNCNAME("Msg::wait()");

    if (w) {
      char line;
      MSG("wait for <enter> ...");
      std::cin >> line;
//       char* result = fgets(line, 9, stdin);
    }
  }


  void Msg::change_out(std::ostream  *fp)
  {
    FUNCNAME("Msg::change_out()");
  
    if (fp) {
      if (out && *out != std::cout && *out != std::cerr) {
	dynamic_cast< std::ofstream*>(out)->close();
	delete out;
      }

      out = fp;
    } else {
      ERROR("file pointer is pointer to nil;\n");
      ERROR("use previous stream for errors furthermore\n");
    } 
  }


  void Msg::change_error_out(std::ofstream *fp)
  {
    FUNCNAME("Msg::change_error_out()");

    if (fp) {
      if (error && *error != std::cout && *error != std::cerr) {
	dynamic_cast< std::ofstream*>(error)->close();
	delete error;
      }
      
      error = fp;
    } else {
      ERROR("file pointer is pointer to nil;\n");
      ERROR("use previous stream for errors furthermore\n");
    }
  }


  void Msg::open_error_file(const char *filename, OPENMODE type)
  {
    FUNCNAME("Msg::open_error_file()");
    std::ofstream *fp;

    if (filename && (fp = new std::ofstream(filename, type))) {
      if (error && *error != std::cout && *error != std::cerr) {
	dynamic_cast< std::ofstream*>(error)->close();
	delete error;
      }

      error = fp;
    } else {
      if (filename)
	ERROR("can not open %s;\n", filename);
      else
	ERROR("no filename specified;\n");
      ERROR("use previous stream for errors furthermore\n");
    }
  }


  void Msg::print_funcname(const char *funcName)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    if (outputMainRank && MPI::COMM_WORLD.Get_rank() != 0)
      return;
#endif

    if (!out) 
      out = &std::cout;

    if (funcName &&  oldFuncName.get() != funcName) {
      PRINT_LINE((*out), funcName << ":" << std::endl);
    } else if (!funcName) {
      PRINT_LINE((*out), "*unknown function*" << std::endl);
    }
    PRINT_LINE((*out), "               ");

    oldFuncName.set(funcName);
  }


  void Msg::print_error_funcname(const char *funcName, const char *file, int line)
  {
    static int old_line = -1;

    if (!error) 
      error = &std::cerr;

    std::stringstream oss;

    if (funcName && oldFuncName.get() != funcName) {
      oss << funcName << ": ";
    } else if (!funcName) {
      if (line-old_line > 5) 
	oss << "*unknown function*";
    }

    if (oldFuncName.get() != funcName) {
      oss << "ERROR in " << file << ", line " << line << std::endl;;
      oldFuncName.set(funcName);
    } else if (line - old_line > 5)
      oss << "ERROR in " << file << ", line " << line << "\n" << std::endl;

    PRINT_LINE((*error), oss.str());
    old_line = line;
  }


  void Msg::print_error_exit(const char *format, ...)
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


  void Msg::print_error(const char *format, ...)
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


  void Msg::print_warn_funcname(const char *funcName,
				const char *file, 
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

    if (funcName  &&  oldFuncName.get() != funcName) {
      oss << funcName << ": ";
    } else if (!funcName) {
      oss << "*unknown function*";
    }

    if (oldFuncName.get() != funcName) {
      oss << "WARNING in " << file << ", line " << line << std::endl;
      oldFuncName.set(funcName);
    } else if (line - old_line > 5) {
      oss << "WARNING in " << file << ", line " << line << std::endl;
    }

    if (oss.str() != "") 
      PRINT_LINE((*out), oss.str());
    
    old_line = line;
  }


  void Msg::print_warn(const char *format, ...)
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


  void Msg::print(const char *format, ...)
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


  void Global::init()
  {
    int d = -1;

    // get dimension
    TEST_EXIT(Parameters::initialized())("Parameters not initialized!\n");
    Parameters::get("dimension of world",d,0);
    TEST_EXIT(d > 0)("Cannot initialize dimension!\n");
    TEST_EXIT((d == 1) || (d == 2) || (d == 3))("Invalid world dimension %d!\n",d);

    // set dimension
    dimOfWorld = d;

    // prepare geoIndex-Table
    int geoTableSize = abs(static_cast<int>(MINPART)) + MAXPART + 1;
    geoIndexTable.resize(4);
    for (int i = 0; i < 4; i++) {
      geoIndexTable[i].resize(geoTableSize);
      for (int j = 0; j < geoTableSize; j++)
	geoIndexTable[i][j] = 0;      
    }

    geoIndexTable[0][PARTS - MINPART] = 1;
    geoIndexTable[0][VERTEX - MINPART] = 1;
    geoIndexTable[0][EDGE - MINPART] = 0;
    geoIndexTable[0][FACE - MINPART] = 0;
    geoIndexTable[0][WORLD - MINPART] = dimOfWorld;

    for (int i = 1; i < 4; i++) {
      geoIndexTable[i][CENTER - MINPART] = referenceElement[i]->getGeo(CENTER);
      geoIndexTable[i][VERTEX - MINPART] = referenceElement[i]->getGeo(VERTEX);
      geoIndexTable[i][EDGE - MINPART] = referenceElement[i]->getGeo(EDGE);
      geoIndexTable[i][FACE - MINPART] = referenceElement[i]->getGeo(FACE);
      geoIndexTable[i][DIMEN - MINPART] = referenceElement[i]->getGeo(DIMEN);
      geoIndexTable[i][PARTS - MINPART] = referenceElement[i]->getGeo(PARTS);
      geoIndexTable[i][NEIGH - MINPART] = referenceElement[i]->getGeo(NEIGH);
      geoIndexTable[i][WORLD - MINPART] = dimOfWorld;
      geoIndexTable[i][BOUNDARY - MINPART] = referenceElement[i]->getGeo(BOUNDARY);
      geoIndexTable[i][PROJECTION - MINPART] = referenceElement[i]->getGeo(PROJECTION);
    }

    // set msgWait
    Msg::setMsgWait(!(Parameters::getMsgWait() == 0));
  }


  void Global::clear()
  {
    delete referenceElement[1];
    delete referenceElement[2];
    delete referenceElement[3];
  }


  int fac(int i)
  {
    if (i <= 1) 
      return 1;
    else 
      return i * fac(i - 1);
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

    if (inMegaByte) {
      vm_usage /= 1024.0;
      resident_set /= 1024.0;
    }
#else
    ERROR("Function not available under MS Windows\n");
#endif
  }


}
