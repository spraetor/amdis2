#pragma once

#include <vector>
#include <string>

#include "AMDiS_fwd.hpp"
#include "io/FileWriter.hpp"

namespace AMDiS
{
  namespace io
  {
    /** \ingroup Output
     * \brief Writer that generates file readable by GNUPlot
     **/
    class GNUPlotWriter : public FileWriterInterface
    {
    public:
      ///
      GNUPlotWriter(std::string filename,
                    const FiniteElemSpace* feSpace,
                    std::vector<DOFVector<double>*>& dofVectors);

      ///
      virtual ~GNUPlotWriter() {}

      ///
      virtual void writeFiles(AdaptInfo& adaptInfo, bool force,
                              int level = -1,
                              Flag traverseFlag = Mesh::CALL_LEAF_EL,
                              bool (*writeElem)(ElInfo*) = NULL);


      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      static void writeFile(Container& /*vec*/, std::string /*filename*/,
                            AdaptInfo& /*adaptInfo*/)
      {
        ERROR_EXIT("GNUPlotWriter not implemented for this container type!\n");
      }


      /// writes a vector of \ref DOFVector to a file.
      static void writeFile(std::vector<DOFVector<double>*>& dofVectors,
                            std::string filename,
                            AdaptInfo& adaptInfo);


      /// writes a \ref DOFVector to a file. Using container pointer.
      static void writeFile(DOFVector<double>* dofVector,
                            std::string filename,
                            AdaptInfo& adaptInfo)
      {
        std::vector<DOFVector<double>*> dofVectors;
        dofVectors.push_back(dofVector);
        writeFile(dofVectors, filename, adaptInfo);
      }


      /// writes a \ref DOFVector to a file. Using container reference.
      static void writeFile(DOFVector<double>& dofVector,
                            std::string filename,
                            AdaptInfo& adaptInfo)
      {
        writeFile(&dofVector, filename, adaptInfo);
      }


      /// writes a \ref SystemVector to a file. Using container pointer.
      static void writeFile(SystemVector* vecs,
                            std::string filename,
                            AdaptInfo& adaptInfo)
      {
        std::vector<DOFVector<double>*> dofVectors;
        for (int i = 0; i < vecs->getSize(); i++)
          dofVectors.push_back(vecs->getDOFVector(i));
        writeFile(dofVectors, filename, adaptInfo);
      }


      /// writes a \ref SystemVector to a file. Using container reference.
      static void writeFile(SystemVector& vecs,
                            std::string filename,
                            AdaptInfo& adaptInfo)
      {
        writeFile(&vecs, filename, adaptInfo);
      }

    protected:
      /// Contains the mesh
      const FiniteElemSpace* feSpace_;

      /// vector of dof vectors to write
      std::vector<DOFVector<double>*> dofVectors_;

      /// file name
      std::string filename_;
    };

  } // end namespace io
} // end namespace AMdiS
