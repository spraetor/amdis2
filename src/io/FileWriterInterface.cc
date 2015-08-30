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

/** \file FileWriterInterface.cc */


#include "FileWriterInterface.h"
#include "AdaptInfo.h"
#include "Initfile.h"

#include <boost/filesystem.hpp>

namespace AMDiS
{

  bool FileWriterInterface::doWriteTimestep(AdaptInfo* adaptInfo, bool force)
  {
    if (!force)
    {
      if (timeModulo > 0.0)
      {
        if (lastWriteTime > adaptInfo->getStartTime()
            && adaptInfo->getTime() < lastWriteTime + timeModulo)
          return false;
      }
      else
      {
        if (adaptInfo->getTimestepNumber() % tsModulo != 0)
          return false;
      }
    }

    lastWriteTime = adaptInfo->getTime();
    return true;
  }


  void FileWriterInterface::readParameters(std::string name)
  {
    Parameters::get(name + "->filename", filename);
    Parameters::get(name + "->append index", appendIndex);
    Parameters::get(name + "->index length", indexLength);
    Parameters::get(name + "->index decimals", indexDecimals);
    Parameters::get(name + "->write every i-th timestep", tsModulo);
    Parameters::get(name + "->write after timestep", timeModulo);

    Parameters::get(name + "->ParaView create subdirectory", createSubDir);
    if (createSubDir < 0)
      Parameters::get(name + "->create subdirectory", createSubDir);
  }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  void FileWriterInterface::getFilename(AdaptInfo* adaptInfo, std::string& fn, std::string& paraFilename, std::string& postfix)
#else
  void FileWriterInterface::getFilename(AdaptInfo* adaptInfo, std::string& fn)
#endif
  {
    fn = filename;

    if (createSubDir > 0)
    {
      using namespace boost::filesystem;
      path vtu_path = fn;
      path data_basedir("data");
      path vtu_filename = vtu_path.filename();
      vtu_path.remove_filename() /= data_basedir;
      try
      {
        create_directory(vtu_path);
        vtu_path /= vtu_filename;
        fn = vtu_path.string();
      }
      catch (...) {}
    }

#if HAVE_PARALLEL_DOMAIN_AMDIS
    paraFilename = filename;
    fn += "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-";
    postfix = "";
#endif

    if (appendIndex)
    {
      TEST_EXIT(indexLength <= 99)("index lenght > 99\n");
      TEST_EXIT(indexDecimals <= 97)("index decimals > 97\n");
      TEST_EXIT(indexDecimals < indexLength)("index length <= index decimals\n");

      char formatStr[9];
      char timeStr[20];

      sprintf(formatStr, "%%0%d.%df", indexLength, indexDecimals);
      sprintf(timeStr, formatStr, adaptInfo ? adaptInfo->getTime() : 0.0);

      fn += timeStr;
#if HAVE_PARALLEL_DOMAIN_AMDIS
      paraFilename += timeStr;
      postfix += timeStr;
#endif
    }
  }
}
