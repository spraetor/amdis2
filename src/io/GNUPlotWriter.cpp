#include "io/GNUPlotWriter.hpp"

#include "AdaptInfo.hpp"
#include "DOFVector.hpp"
#include "ElInfo.hpp"
#include "FiniteElemSpace.hpp"
#include "Global.hpp"
#include "Mesh.hpp"
#include "Traverse.hpp"

using namespace std;

namespace AMDiS
{
  namespace io
  {

    GNUPlotWriter::GNUPlotWriter(std::string fn,
                                 const FiniteElemSpace* feSpace,
                                 std::vector<DOFVector<double>*>& dofVectors)
      : feSpace_(feSpace),
        dofVectors_(dofVectors),
        filename_(filename)
    {
      filename = fn;
    }


    void GNUPlotWriter::writeFiles(AdaptInfo& adaptInfo, bool /*force*/,
                                   int, Flag, bool (*)(ElInfo*))
    {
      DOFVector<WorldVector<double>> coords(feSpace_, "coords");
      Mesh* mesh = feSpace_->getMesh();

      double time = adaptInfo.getTime();

      mesh->getDofIndexCoords(coords);

      std::ofstream file(filename.c_str(), ios_base::out);
      TEST_EXIT(file.is_open())("could not open file %s for writing\n", filename_.c_str());

      file.precision(10);
      file.setf(std::ios_base::scientific);

      file << "# line format: time x y z val1[x,y,z] val2[x,y,z] ...\n";

      DOFVector<WorldVector<double>>::Iterator coordsIt(&coords, USED_DOFS);

      int index, numVecs = static_cast<int>(dofVectors_.size());
      for (coordsIt.reset(); !coordsIt.end(); ++coordsIt)
      {
        index = coordsIt.getDOFIndex();
        file << time << " ";
        for (int i = 0; i < Global::getGeo(WORLD); i++)
          file << (*coordsIt)[i] << " ";
        for (int i = 0; i < numVecs; i++)
          file << (*(dofVectors_[i]))[index] << " ";
        file << "\n";
      }
    }


    void GNUPlotWriter::writeFile(std::vector<DOFVector<double>*>& dofVectors,
                                  std::string filename,
                                  AdaptInfo& adaptInfo)
    {
      GNUPlotWriter gnuplotWriter(filename, dofVectors[0]->getFeSpace(), dofVectors);
      gnuplotWriter.writeFiles(adaptInfo, true);
    }

  } // end namespace io
} // end namespace AMDiS
