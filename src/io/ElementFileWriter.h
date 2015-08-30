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



/** \file ElementFileWriter.h */

#ifndef ELEMENTFILEWRITER_H
#define ELEMENTFILEWRITER_H

#include "FileWriter.h"
#include "FiniteElemSpace.h"
#include "MatrixVector.h"
#include "Mesh.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief
    * Filewriter that make it possible to create mesh files, where the values
    * are not defined on DOFs, but instead are defined on the elements.
    *
    * It can be necessary to visualize data defined on elements, e.g. residual
    * error that is defined on elements and not on DOFs. This class takes as
    * input a mesh and a std::map, which defines for each element index a (double)
    * value, and outputs a AMDiS/VTK mesh.
    */
    class ElementFileWriter : public FileWriterInterface
    {
    public:
      /// Constructor.
      ElementFileWriter(std::string name,
                        Mesh* mesh,
                        std::map<int, double>& vec);

      ElementFileWriter(std::string name,
                        Mesh* mesh,
                        std::map<int, std::vector<double>>& vecs);

      /// Implementation of FileWriterInterface::writeFiles().
      void writeFiles(AdaptInfo* adaptInfo, bool force,
                      int level = -1,
                      Flag traverseFlag = Mesh::CALL_LEAF_EL,
                      bool (*writeElem)(ElInfo*) = NULL);

      /// Simple writing procedure for one element std::map.
      static void writeFile(std::map<int, double>& vec,
                            Mesh* mesh,
                            std::string filename,
                            std::string postfix = ".vtu",
                            int level = -1,
                            bool writeAsVector = false);

      static void writeFile(std::map<int, std::vector<double>>& vecs,
                            Mesh* mesh,
                            std::string filename,
                            std::string postfix = ".vtu",
                            int level = -1,
                            bool writeAsVector = true);

    protected:
      /// Writes element data in AMDiS format (1 file !).
      void writeMeshDatValues(std::string filename, double time);

      /// Writes element data in VTK format.
      void writeVtkValues(std::string filename, std::string postfix,
                          int level = -1, bool writeAsVector = false);

      /// Writes a world coordinate to a given file.
      template<typename T>
      void writeCoord(T& file, WorldVector<double> coord)
      {
        for (int i = 0; i < Global::getGeo(WORLD); i++)
          file << " " << std::scientific << coord[i];
        for (int i = Global::getGeo(WORLD); i < 3; i++)
          file << " " << "0.0";

        file << "\n";
      }

    protected:
      /// Name.
      std::string name;

      /// Used filename prefix.
      std::string filename;

      /// AMDiS mesh-data-file extension.
      std::string amdisMeshDatExt;

      /// VTK file extension.
      std::string vtkExt;
      std::string pvdExt;

      /// 0: Don't write AMDiS files.
      /// 1: Write AMDiS files.
      int writeAMDiSFormat;

      /// 0: Don't write VTK files.
      /// 1: Write VTK files.
      int writeVtkFormat;
      int writeVtkVectorFormat;

      bool writeAs3dVector;

      int writeParaViewAnimation;

      /// 0: Don't append time index to filename prefix.
      /// 1: Append time index to filename prefix.
      int appendIndex;

      /// Total length of appended time index.
      int indexLength;

      /// Number of decimals in time index.
      int indexDecimals;

      /// Timestep modulo: write only every tsModulo-th timestep!
      int tsModulo;

      ///
      int timestepNumber;

      /// Mesh used for output.
      Mesh* mesh;

      /// Vector that stores the solution.
      std::map<int, double>* vec;
      std::map<int, std::vector<double>>* vecs;

      std::vector<std::pair<double, std::string>> paraViewAnimationFrames;
    };

  } // end namespace io
} // end namespace AMDiS


#endif  // ELEMENTFILEWRITER_H
