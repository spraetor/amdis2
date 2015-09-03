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


#include <boost/filesystem.hpp>

#include "FileWriter.h"
#include "Initfile.h"
#include "ValueWriter.h"
#include "MacroWriter.h"
#include "VtkWriter.h"
#include "VtkVectorWriter.h"
#include "PngWriter.h"
#include "PovrayWriter.h"
#include "DofWriter.h"
#include "ArhWriter.h"
#include "Arh2Writer.h"
#include "Arh3Writer.h"
#include "FiniteElemSpace.h"
#include "AdaptInfo.h"
#include "Flag.h"
#include "ElInfo.h"
#include "Mesh.h"
#include "DataCollector.h"

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

using namespace std;

namespace AMDiS
{
  namespace detail
  {
    template<>
    FileWriter<double>::FileWriter(std::string name_,
                                   Mesh* mesh_,
                                   SystemVector* vecs)
      : name(name_),
        mesh(mesh_)
    {
      initialize();

      /*
      * Removed by Siqi. not sure.
      * for (int i = 0; i < static_cast<int>(vecs->getSize()); i++)
      * TEST_EXIT(vecs->getDOFVector(0)->getFeSpace() == vecs->getDOFVector(i)->getFeSpace())
      * 	("All FeSpace have to be equal!\n");
      */
      feSpace = vecs->getDOFVector(0)->getFeSpace();
      solutionVecs.resize(vecs->getSize());
      for (int i = 0; i < static_cast<int>(vecs->getSize()); i++)
        solutionVecs[i] = vecs->getDOFVector(i);

      for (size_t i = 0; i < solutionVecs.size(); i++)
        solutionNames.push_back(solutionVecs[i]->getName());
    }


    template<>
    void FileWriter<double>::writeFiles(AdaptInfo* adaptInfo,
                                        bool force,
                                        int level,
                                        Flag flag,
                                        bool (*writeElem)(ElInfo*))
    {
      FUNCNAME("FileWriter<T>::writeFiles()");
      using namespace ::AMDiS::io;

      if (!super::doWriteTimestep(adaptInfo, force))
        return;

      //-----------------by Siqi---------------------//
      if (writeAMDiSFormat || writePeriodicFormat || writeParaViewFormat
          || writeParaViewVectorFormat || writeParaViewAnimation
          || writeDofFormat || writeArh1 || writePovrayFormat)
      {
        for (int i = 0; i < static_cast<int>(solutionVecs.size()); i++)
          TEST_EXIT(solutionVecs[0]->getFeSpace() == solutionVecs[i]->getFeSpace())
          ("All FeSpaces have to be equal!\n");
      }

      // Containers, which store the data to be written;
      std::vector<DataCollector<>*> dataCollectors(solutionVecs.size());

      if (writeElem)
      {
        for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
          dataCollectors[i] = new DataCollector<>(feSpace, solutionVecs[i],
                                                  level, flag, writeElem);
      }
      else
      {
        for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
          dataCollectors[i] = new DataCollector<>(feSpace, solutionVecs[i],
                                                  traverseLevel,
                                                  flag | traverseFlag,
                                                  writeElement);
      }

      std::string fn, fn_;
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      std::string paraFilename, postfix;
      super::getFilename(adaptInfo, fn, paraFilename, postfix);
      postfix += paraviewFileExt;
      fn_ = paraFilename;
#else
      super::getFilename(adaptInfo, fn);
      fn_ = fn;
#endif


      if (writeAMDiSFormat)
      {
        MacroWriter::writeMacro(dataCollectors[0],
                                const_cast<char*>((fn +  amdisMeshExt).c_str()),
                                adaptInfo ? adaptInfo->getTime() : 0.0);
        MSG("macro file written to %s\n", (fn + amdisMeshExt).c_str());

        ValueWriter::writeValues(dataCollectors[0],
                                 (fn + amdisDataExt).c_str(),
                                 adaptInfo ? adaptInfo->getTime() : 0.0);
        MSG("value file written to %s\n", (fn + amdisDataExt).c_str());
      }

      if (writePeriodicFormat)
      {
        MacroWriter::writePeriodicFile(dataCollectors[0],
                                       (fn + periodicFileExt).c_str());
        MSG("periodic file written to %s\n", (fn + periodicFileExt).c_str());
      }

      if (writeParaViewFormat)
      {
        std::string vtu_file = fn + paraviewFileExt;
        VtkWriter::Aux vtkWriter(&dataCollectors,
                                 solutionNames,
                                 VtkWriter::Vtuformat(paraViewMode), (paraViewPrecision == 1), writeParaViewVectorFormat);
        vtkWriter.writeFile(vtu_file);

#if HAVE_PARALLEL_DOMAIN_AMDIS
        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          // 	  vector<string> componentNames;
          // 	  for (unsigned int i = 0; i < dataCollectors.size(); i++)
          // 	    componentNames.push_back(dataCollectors[i]->getValues()->getName());
          VtkWriter::detail::writeParallelFile(paraFilename + paraviewParallelFileExt,
                                               MPI::COMM_WORLD.Get_size(),
                                               filename,
                                               postfix,
                                               solutionNames,
                                               VtkWriter::Vtuformat(paraViewMode),
                                               (paraViewPrecision == 1),
                                               writeParaViewVectorFormat,
                                               createSubDir > 0);
        }
#endif

        MSG("ParaView file written to %s\n", (fn + paraviewFileExt).c_str());
      }

      // write vtu-vector files
      if (writeParaViewVectorFormat && !writeParaViewFormat)
      {
        VtkVectorWriter::writeFile(solutionVecs, fn_ + paraviewFileExt, true, writeAs3dVector);
        MSG("ParaView file written to %s\n", (fn_ + paraviewFileExt).c_str());
      }

      if (writeParaViewAnimation)
      {
        std::string pvd_file = fn_ + paraviewFileExt;
#if HAVE_PARALLEL_DOMAIN_AMDIS
        pvd_file = fn_ + paraviewParallelFileExt;
        if (MPI::COMM_WORLD.Get_rank() == 0)
#endif
        {
          VtkWriter::detail::updateAnimationFile(adaptInfo,
                                                 pvd_file,
                                                 &paraviewAnimationFrames,
                                                 filename + ".pvd");
        }
      }


      if (writeDofFormat)
      {
        DofWriter::writeFile(solutionVecs, fn + ".dof");
      }

      // write Arh files
      if (writeArh2)
        Arh2Writer::writeFile(solutionVecs, fn_ + ".arh");
      else if (writeArh1)
        ArhWriter::write(fn_ + ".arh", feSpace->getMesh(), solutionVecs);
      else if (writeArh3 || writeArhFormat)
        Arh3Writer::writeFile(solutionVecs, fn_ + ".arh");


#ifdef HAVE_PNG
      if (writePngFormat)
      {
        PngWriter pngWriter(dataCollectors[0]);
        pngWriter.writeFile(fn + ".png", pngType);

        MSG("PNG image file written to %s\n", (fn + ".png").c_str());
      }
#endif

      if (writePovrayFormat)
      {
        PovrayWriter povrayWriter(dataCollectors[0]);
        povrayWriter.writeFile(fn + ".pov");

        MSG("Povray script written to %s\n", (fn + ".pov").c_str());
      }


      for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
        delete dataCollectors[i];
    }

    template<>
    string FileWriter<double>::getParaViewFilename(AdaptInfo* adaptInfo) const
    {
      string ret(filename);
      if (appendIndex)
      {
        TEST_EXIT(indexLength <= 99)("index lenght > 99\n");
        TEST_EXIT(indexDecimals <= 97)("index decimals > 97\n");
        TEST_EXIT(indexDecimals < indexLength)("index length <= index decimals\n");

        char formatStr[9];
        char timeStr[20];

        sprintf(formatStr, "%%0%d.%df", indexLength, indexDecimals);
        sprintf(timeStr, formatStr, adaptInfo ? adaptInfo->getTime() : 0.0);

        ret += timeStr;
      }
      return ret;
    }
  } // end namespace detail
} // end namespace AMDiS
