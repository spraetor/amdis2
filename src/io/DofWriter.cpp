#include "io/DofWriter.hpp"

#include "DOFVector.hpp"
#include "BasisFunction.hpp"

using namespace std;

namespace AMDiS
{
  namespace io
  {

    namespace DofWriter
    {

      void writeFile(std::vector<DOFVector<double>*>& vec, std::string filename)
      {
        const FiniteElemSpace* feSpace = vec[0]->getFeSpace();
        Mesh* mesh = feSpace->getMesh();

        DOFVector<WorldVector<double>> coordDof(feSpace, "tmp");
        const BasisFunction* basFcts = feSpace->getBasisFcts();
        int nBasFcts = basFcts->getNumber();
        std::vector<DegreeOfFreedom> dofVec(nBasFcts);

        TraverseStack stack;
        ElInfo* elInfo =
          stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
        while (elInfo)
        {
          basFcts->getLocalIndices(elInfo->getElement(),
                                   feSpace->getAdmin(), dofVec);
          for (int i = 0; i < nBasFcts; i++)
          {
            DimVec<double>* baryCoords = basFcts->getCoords(i);
            elInfo->coordToWorld(*baryCoords, coordDof[dofVec[i]]);
          }

          elInfo = stack.traverseNext(elInfo);
        }

        std::ofstream outfile;
        outfile.open(filename.c_str());
        outfile.setf(std::ios::scientific, std:: ios::floatfield);
        outfile.precision(16);

        DOFIterator<WorldVector<double>> it(&coordDof, USED_DOFS);
        for (it.reset(); !it.end(); ++it)
        {
          outfile << it.getDOFIndex() << " ";

          for (int i = 0; i < Global::getGeo(WORLD); i++)
            outfile << (*it)[i] << " ";

          for (unsigned int i = 0; i < vec.size(); i++)
            outfile << (*(vec[i]))[it.getDOFIndex()] << " ";
          outfile << endl;
        }

        outfile.close();
      }

    } // end namespace DofWriter
  }
} // end namespace io, AMDiS
