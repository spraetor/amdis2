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

/** \file FileWriter.hh */

#ifndef AMDIS_FILEWRITER_HH
#define AMDIS_FILEWRITER_HH

#include "boost/lexical_cast.hpp"
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
#include "FiniteElemSpace.h"
#include "AdaptInfo.h"
#include "Flag.h"
#include "ElInfo.h"
#include "Mesh.h"
#include "SystemVector.h"
#include "DataCollector.h"

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

namespace AMDiS
{
  namespace detail
  {
      
    template<typename T>
    FileWriter<T>::FileWriter(std::string str, Mesh *m, DOFVector<T> *vec)
      : name(str),
	mesh(m)
    {
      initialize();

      feSpace = vec->getFeSpace();

      solutionVecs.resize(1);
      solutionVecs[0] = vec;
      
      solutionNames.push_back(vec->getName());
    }


    template<typename T>
    FileWriter<T>::FileWriter(std::string name_,
			  Mesh *mesh_, 
			  std::vector<DOFVector<T>*> vecs,
			  std::vector< std::string > componentNames)
      : name(name_),
	mesh(mesh_),
	solutionNames(componentNames)
    {
      initialize();

  /** removed by Siqi. Not sure.
  *  for (int i = 0; i < static_cast<int>(vecs.size()); i++)
  *    TEST_EXIT(vecs[0]->getFeSpace() == vecs[i]->getFeSpace())
  *      ("All FeSpace have to be equal!\n");
  */
      feSpace = vecs[0]->getFeSpace();
      solutionVecs = vecs;
      
      if (solutionNames.size() < vecs.size()) {
	for (size_t i = solutionNames.size(); i < vecs.size(); i++)
	  solutionNames.push_back(vecs[i]->getName());
      }
    }


    template<typename T>
    FileWriter<T>::FileWriter(std::string name_,
			  Mesh *mesh_,
			  SystemVector *vecs)
      : name(name_),
	mesh(mesh_)
    {
      ERROR_EXIT("SystemVector contains DOFVectors of type double, so the FileWriter<not double> can not be used!\n");
    }


    template<typename T>
    FileWriter<T>::~FileWriter()
    {
      // Do not forget to delete temporal solution vector, if there have been
      // some created in the constructor.
      if (nTmpSolutions > 0)
	for (int i = 0; i < nTmpSolutions; i++)
	  delete solutionVecs[i]; 
    }
    

    template<typename T>
    void FileWriter<T>::initialize()
    {
      amdisMeshExt = ".mesh";
      amdisDataExt = ".dat";
      paraviewFileExt = ".vtu";
      paraviewParallelFileExt = ".pvtu";
      periodicFileExt = ".per";
      writeAMDiSFormat = 0;
      writeParaViewFormat = 0;
      paraViewMode = 0;
      paraViewPrecision = 0;
      writeParaViewVectorFormat = 0;
      writeAs3dVector = false;
      writeParaViewAnimation = 0;
      writePeriodicFormat = 0;
      writePngFormat = 0;
      writePovrayFormat = 0;
      writeDofFormat = 0;
      writeArhFormat = 0;
      writeArh1 = 0;
      writeArh2 = 0;
      writeArh3 = 0;
      pngType = 0;
      nTmpSolutions = 0;
      paraviewAnimationFrames.resize(0),
      compression = NONE;

      readParameters(name);
    }


    template<typename T>
    void FileWriter<T>::readParameters(std::string name)
    {
      super::readParameters(name);
      
      Parameters::get(name + "->AMDiS format", writeAMDiSFormat);
      Parameters::get(name + "->AMDiS mesh ext", amdisMeshExt);
      Parameters::get(name + "->AMDiS data ext", amdisDataExt);
      Parameters::get(name + "->ParaView format", writeParaViewFormat);
      Parameters::get(name + "->ParaView mode", paraViewMode);
      Parameters::get(name + "->ParaView precision", paraViewPrecision);
      Parameters::get(name + "->ParaView vector format", writeParaViewVectorFormat);
      Parameters::get(name + "->write vector as 3d vector", writeAs3dVector);
      Parameters::get(name + "->ParaView animation", writeParaViewAnimation);
      Parameters::get(name + "->ParaView ext", paraviewFileExt);    
      Parameters::get(name + "->Periodic format", writePeriodicFormat);
      Parameters::get(name + "->Periodic ext", periodicFileExt);
      Parameters::get(name + "->PNG format", writePngFormat);
      Parameters::get(name + "->PNG type", pngType);

      Parameters::get(name + "->Povray format", writePovrayFormat);
      Parameters::get(name + "->Povray template", povrayTemplate);
      Parameters::get(name + "->Povray camera location", povrayCameraLocation);
      Parameters::get(name + "->Povray camera look_at", povrayCameraLookAt);

      Parameters::get(name + "->DOF format", writeDofFormat);
      Parameters::get(name + "->ARH format", writeArhFormat);
      Parameters::get(name + "->ARH1 format", writeArh1);
      Parameters::get(name + "->ARH2 format", writeArh2);
      Parameters::get(name + "->ARH3 format", writeArh3);

      std::string compressionStr = "";
      Parameters::get(name + "->compression", compressionStr);
      if (compressionStr == "gzip" || compressionStr == "gz") {
	compression = GZIP;
      } else if (compressionStr == "bzip2" || compressionStr == "bz2") {
	compression = BZIP2;
      }
    }


    template<typename T>
    void FileWriter<T>::writeFiles(AdaptInfo *adaptInfo,
					    bool force,
					    int level,
					    Flag flag,
					    bool (*writeElem)(ElInfo*))
    {
      FUNCNAME("FileWriter<T>::writeFiles()");

      if (!super::doWriteTimestep(adaptInfo, force))
	return;
      
      // Containers, which store the data to be written;
      std::vector<DataCollector<T>*> dataCollectors(solutionVecs.size());
      
      if (writeElem) {
	for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
	  dataCollectors[i] = new DataCollector<T>(feSpace, solutionVecs[i], 
						  level, flag, writeElem);
      } else {
	for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
	  dataCollectors[i] = new DataCollector<T>(feSpace, solutionVecs[i], 
						  traverseLevel, 
						  flag | traverseFlag, 
						  writeElement);
      }
      
     std::string fn;
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      std::string paraFilename, postfix;
      super::getFilename(adaptInfo, fn, paraFilename, postfix);
      postfix += paraviewFileExt;
#else
      super::getFilename(adaptInfo, fn);
#endif

      if (writeParaViewVectorFormat) {
	io::VtkVectorWriter::Aux<T> vtkVectorWriter(&dataCollectors, writeAs3dVector);
#ifdef HAVE_COMPRESSION
	vtkVectorWriter.setCompression(compression);
#endif
	vtkVectorWriter.writeFile(fn + paraviewFileExt);

  #if HAVE_PARALLEL_DOMAIN_AMDIS
	if (MPI::COMM_WORLD.Get_rank() == 0)
	  vtkVectorWriter.writeParallelFile(paraFilename + paraviewParallelFileExt,
					    MPI::COMM_WORLD.Get_size(),
					    filename, postfix);
  #endif
	  
	MSG("ParaView file written to %s\n", (fn + paraviewFileExt).c_str());
      }
      
      if (writeParaViewAnimation) {
  #if HAVE_PARALLEL_DOMAIN_AMDIS
	if (MPI::COMM_WORLD.Get_rank() == 0) {
	  io::VtkWriter::detail::updateAnimationFile(adaptInfo,
					paraFilename + paraviewParallelFileExt,
					&paraviewAnimationFrames,
					filename + ".pvd");

	}
  #else
	io::VtkWriter::detail::updateAnimationFile(adaptInfo,
				      fn + paraviewFileExt,
				      &paraviewAnimationFrames,
				      filename + ".pvd");
  #endif
      }
      
      for (int i = 0; i < static_cast<int>(dataCollectors.size()); i++)
	delete dataCollectors[i];
    }
    
  } // end namespace detail
} // end namespace AMDiS

#endif

