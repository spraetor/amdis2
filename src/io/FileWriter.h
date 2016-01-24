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



/** \file FileWriter.h */

/** \defgroup Output Output module
 * @{ <img src="output.png"> @}
 */

#ifndef AMDIS_FILEWRITER_H
#define AMDIS_FILEWRITER_H

#include <vector>
#include <string>
#include "AMDiS_fwd.h"
#include "FileWriterInterface.h"
#include "Mesh.h"
#include "DataCollector.h"
#include "FileCompression.h"

namespace AMDiS
{

  namespace detail
  {

    /**
    * \ingroup Output
    *
    * \brief
    * Base class of FileWriterScal and FileWriterVec. Manages the file
    * output of solution vectors.
    */
    template<typename T>
    class FileWriter : public ::AMDiS::FileWriterInterface
    {
      typedef ::AMDiS::FileWriterInterface super;

    public:
      /// Constructor for a filewriter for one data component.
      FileWriter(std::string name, Mesh* mesh, DOFVector<T>* vec);

      /// Constructor for a filewriter with more than one data component.
      FileWriter(std::string name,
                 Mesh* mesh,
                 std::vector<DOFVector<T>*> vecs,
                 std::vector<std::string> componentNames = std::vector<std::string>()
                );

      /// Constructor for a filewriter with more than one data component.
      FileWriter(std::string name,
                 Mesh* mesh,
                 SystemVector* vecs);

      /// Destructor
      virtual ~FileWriter();

      /// Implementation of FileWriterInterface::writeFiles().
      virtual void writeFiles(AdaptInfo& adaptInfo, bool force,
                              int level = -1,
                              Flag traverseFlag = Mesh::CALL_LEAF_EL,
                              bool (*writeElem)(ElInfo*) = NULL) override;

      std::vector<std::pair<double, std::string>>& getParaviewAnimationFrames()
      {
        return paraviewAnimationFrames;
      }

      bool getWriteParaViewFormat() const
      {
        return writeParaViewFormat;
      }

      std::string getParaViewFilename(AdaptInfo& info) const;

      const std::vector<std::string>& getSolutionNames() const
      {
        return solutionNames;
      }

    protected:
      /// Initialization of the filewriter.
      void initialize();

      /// Reads all file writer dependend parameters from the init file.
      virtual void readParameters(std::string name) override;

      /// Name of the writer.
      std::string name;

      /// AMDiS mesh-file extension.
      std::string amdisMeshExt;

      /// AMDiS solution-file extension.
      std::string amdisDataExt;

      /// VTK file extension.
      std::string paraviewFileExt;

      /// Parallel VTK file extension.
      std::string paraviewParallelFileExt;

      /// Periodic file extension.
      std::string periodicFileExt;

      /// 0: Don't write AMDiS files; 1: Write AMDiS files.
      int writeAMDiSFormat;

      /// 0: Don't write ParaView files; 1: Write ParaView files.
      int writeParaViewFormat;

      /// 0: ASCII mode; 1: Appended mode; 2:Appended_compressed mode.
      int paraViewMode;

      /// 0: FLOAT32 precision; 1: FLOAT64 precision. Only works in appended and appended_compressed mode.
      int paraViewPrecision;

      /// 0: Don't write ParaView std::vector files; 1: Write ParaView std::vector files.
      int writeParaViewVectorFormat;

      /// 1: extend number of component to 3, so that paraview can display the std::vector as worldstd::vector
      bool writeAs3dVector;

      /// 0: Don't write ParaView animation file; 1: Write ParaView animation file.
      int writeParaViewAnimation;

      /// 0: Don't write periodic files; 1: Write periodic files.
      int writePeriodicFormat;

      /// 0: Don't write png files; 1: Write png image files.
      int writePngFormat;

      /// 0: Gray color picture; 1: RGB picture.
      int pngType;

      /// 0: Don't write Povray scripts; 1: Write Povray scripts
      int writePovrayFormat;

      /// 0: Don't write DOF files; 1: Write DOF files
      int writeDofFormat;

      /// if write latest ARH files
      int writeArhFormat;

      /// write Arh2, prior to writeArhFormat
      int writeArh1;

      /// write Arh2 version 2.1,  prior to writeArhFormat
      int writeArh2;

      /// write Arh2 version 3.0,  prior to writeArhFormat
      int writeArh3;

      /// camera position for povray script files
      std::string povrayCameraLocation;

      /// orientation for camera in povray script files
      std::string povrayCameraLookAt;

      /// name of the template file that will be prepended to all created *.pov files
      std::string povrayTemplate;


      /// Stores a set of std::pairs of timepoint and filename to write a ParaView
      /// animation file.
      std::vector<std::pair<double, std::string>> paraviewAnimationFrames;

      ///
      int timestepNumber;

      /// Mesh used for output.
      Mesh* mesh;

      /// fespace used for output.
      const FiniteElemSpace* feSpace;

      /// Pointers to the std::vectors which store the solution.
      std::vector<DOFVector<T>*> solutionVecs;

      /// Names of the DOFVectors
      std::vector<std::string> solutionNames;

      /** \brief
      * Stores the number of temporal solutions std::vectors, which have been created
      * in the constructor. If this number is greater than zero, the std::vectors
      * stored in solutionVecs_ must be deleted in the destructor.
      */
      int nTmpSolutions;

      /** \brief
      * Defines if, and with what kind of compression, the file should be compressed
      * during writing.
      */
      FileCompression compression;
    };

    template<>
    FileWriter<double>::FileWriter(std::string name_,
                                   Mesh* mesh_,
                                   SystemVector* vecs);

    template<>
    void FileWriter<double>::writeFiles(AdaptInfo& adaptInfo,
                                        bool force,
                                        int level,
                                        Flag flag,
                                        bool (*writeElem)(ElInfo*));

  } // end namespace detail

  typedef detail::FileWriter<double> FileWriter;
  typedef detail::FileWriter<WorldVector<double>> FileVectorWriter;
}

#include "FileWriter.hh"

#endif // AMDIS_FILEWRITER_H
