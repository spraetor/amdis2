#include <map>
#include <utility>

#include "ElementDofIterator.hpp"
#include "ProblemStatDbg.hpp"
#include "SystemVector.hpp"

namespace AMDiS
{
  using namespace std;

  void ProblemStatDbg::writeDbgMatrix(string filename)
  {
    FUNCNAME("ProblemStatDbg::writeDbMatrix()");

    using mtl::tag::major;
    using mtl::tag::nz;
    using mtl::begin;
    using mtl::end;
    namespace traits = mtl::traits;
    typedef DOFMatrix::base_matrix_type Matrix;
    typedef traits::range_generator<major, Matrix>::type cursor_type;
    typedef traits::range_generator<nz, cursor_type>::type icursor_type;

    // Create map with all world coordinates of all DOFs.
    DofToCoord dofToCoords;
    createDofToCoordMap(dofToCoords);

    // Start writing output file.
    ofstream out(filename.c_str());
    out.precision(15);

    // First, we write the number of DOFs the file contains.
    out << dofToCoords.size() << endl;

    int dimOfWorld = Global::getGeo(WORLD);
    // And now write all the DOF's coords. The first line contains the dof index, all
    // the other lines contain one component of the coordinates.
    for (DofToCoord::iterator it = dofToCoords.begin(); it != dofToCoords.end(); ++it)
    {
      out << it->first << endl;
      for (int j = 0; j < dimOfWorld; j++)
        out << (it->second)[j] << endl;
    }

    // Write the matrices.
    for (int i = 0; i < nComponents; i++)
    {
      for (int j = 0; j < nComponents; j++)
      {
        DOFMatrix* dofmatrix = (*systemMatrix)[i][j];

        if (!dofmatrix)
          continue;

        Matrix& matrix = dofmatrix->getBaseMatrix();
        int nnz = matrix.nnz();
        int testNnz = 0;

        // Write to file, how many entries the matrix conatins.
        out << nnz << endl;

        traits::row<Matrix>::type row(matrix);
        traits::col<Matrix>::type col(matrix);
        traits::const_value<Matrix>::type value(matrix);

        for (cursor_type cursor = begin<major>(matrix), cend = end<major>(matrix); cursor != cend; ++cursor)
        {
          for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); icursor != icend; ++icursor)
          {
            out << row(*icursor) << " "
                << col(*icursor) << " "
                << value(*icursor) << endl;
            testNnz++;
          }
        }

        // Just to check, if nnz() is correct.
        TEST_EXIT(testNnz == nnz)("NNZ does not fit together!\n");
      }
    }


    for (int i = 0; i < nComponents; i++)
    {
      out << rhs->getDOFVector(i)->getUsedSize() << endl;
      DOFIterator<double> it(rhs->getDOFVector(i), USED_DOFS);
      int c = 0;
      for (it.reset(); !it.end(); ++it)
      {
        out << c << " " << *it << endl;
        c++;
      }

      TEST_EXIT(c == rhs->getDOFVector(i)->getUsedSize())
      ("RHS NNZ does not fit together!\n");
    }


    for (int i = 0; i < nComponents; i++)
    {
      out << solution->getDOFVector(i)->getUsedSize() << endl;
      DOFIterator<double> it(solution->getDOFVector(i), USED_DOFS);
      int c = 0;
      for (it.reset(); !it.end(); ++it)
      {
        out << c << " " << *it << endl;
        c++;
      }

      TEST_EXIT(c == rhs->getDOFVector(i)->getUsedSize())
      ("RHS NNZ does not fit together!\n");
    }

    out.close();
  }


  void ProblemStatDbg::readAndCompareDbgMatrix(vector<string> filenames)
  {
    using mtl::tag::major;
    using mtl::tag::nz;
    using mtl::begin;
    using mtl::end;
    namespace traits = mtl::traits;
    //     typedef DOFMatrix::base_matrix_type Matrix;
    //     typedef traits::range_generator<major, Matrix>::type cursor_type;
    //     typedef traits::range_generator<nz, cursor_type>::type icursor_type;

    // Create a map from coords of all DOFs, to the DOF indices in this problem.
    CoordToDof coordToDof;
    createCoordToDofMap(coordToDof);

    int dimOfWorld = Global::getGeo(WORLD);
    vector<vector<map<pair<DegreeOfFreedom, DegreeOfFreedom>, double>>> nnzValues(nComponents);
    vector<map<DegreeOfFreedom, double>> rhsValues(nComponents);
    vector<map<DegreeOfFreedom, double>> solValues(nComponents);
    for (int i = 0; i < nComponents; i++)
      nnzValues[i].resize(nComponents);

    // Stores to each dof index of this problem a map from rank indices (of each rank
    // that also has this dof) to the corresponding local dof index.
    map<DegreeOfFreedom, vector<DegreeOfFreedom>> dofMapHereToFiles;

    for (vector<string>::iterator fileIt = filenames.begin();
         fileIt != filenames.end(); ++fileIt)
    {

      // Open file and read the number of DOFs the file contains.
      ifstream in(fileIt->c_str());
      int nReadDofs;
      in >> nReadDofs;

      // Is used to map the dof indices in the files to the global coordinates.
      DofToCoord dofToCoord;

      // Read all DOF indices and their world coordinates.
      for (int i = 0; i < nReadDofs; i++)
      {
        DegreeOfFreedom dof;
        WorldVector<double> coords;

        in >> dof;
        for (int j = 0; j < dimOfWorld; j++)
          in >> coords[j];
        dofToCoord[dof] = coords;
      }

      map<DegreeOfFreedom, DegreeOfFreedom> dofMapFileToHere;

      for (DofToCoord::iterator it = dofToCoord.begin(); it != dofToCoord.end(); ++it)
      {
        DegreeOfFreedom dofIndexInFile = it->first;
        WorldVector<double>& dofCoords = it->second;

        if (coordToDof.find(dofCoords) == coordToDof.end())
        {
          cout << "Cannot find dof index for coords: " << endl;
          cout << dofCoords;
          exit(0);
        }

        DegreeOfFreedom dofIndexHere = coordToDof[dofCoords];
        dofMapHereToFiles[dofIndexHere].push_back(dofIndexInFile);
        dofMapFileToHere[dofIndexInFile] = dofIndexHere;
      }

      // Now we traverse all matrices and check them.
      for (int i = 0; i < nComponents; i++)
      {
        for (int j = 0; j < nComponents; j++)
        {
          DOFMatrix* dofmatrix = (*systemMatrix)[i][j];

          if (!dofmatrix)
            continue;

          int readNnz;
          in >> readNnz;

          // Read each entry in file and check it with the corresponding matrix value
          // in this problem.
          for (int k = 0; k < readNnz; k++)
          {
            DegreeOfFreedom row, col;
            double value;

            in >> row;
            in >> col;
            in >> value;

            if (dofMapFileToHere.count(row) == 0)
            {
              cout << "Cannot find row index for: " << row << endl;
              exit(0);
            }

            if (dofMapFileToHere.count(col) == 0)
            {
              cout << "Cannot find col index for: " << col << endl;
              exit(0);
            }

            // Get dof indices for row and col of this problem matrix.
            DegreeOfFreedom rowHere = dofMapFileToHere[row];
            DegreeOfFreedom colHere = dofMapFileToHere[col];
            pair<DegreeOfFreedom, DegreeOfFreedom> rowcol = make_pair(rowHere, colHere);

            if (nnzValues[i][j].count(rowcol) == 0)
              nnzValues[i][j][rowcol] = value;
            else
              nnzValues[i][j][rowcol] += value;
          }
        }
      }

      for (int i = 0; i < nComponents; i++)
      {
        int readNnz;
        in >> readNnz;
        for (int k = 0; k < readNnz; k++)
        {
          int row;
          double value;

          in >> row;
          in >> value;

          WorldVector<double> rowCoords = dofToCoord[row];

          if (coordToDof.find(rowCoords) == coordToDof.end())
          {
            cout << "Cannot find row index for coords: " << endl;
            cout << rowCoords;
            exit(0);
          }

          DegreeOfFreedom rowHere = coordToDof[rowCoords];

          if (rhsValues[i].count(rowHere) == 0)
            rhsValues[i][rowHere] = value;
          else
            rhsValues[i][rowHere] += value;
        }
      }

      for (int i = 0; i < nComponents; i++)
      {
        int readNnz;
        in >> readNnz;
        for (int k = 0; k < readNnz; k++)
        {
          int row;
          double value;

          in >> row;
          in >> value;

          WorldVector<double> rowCoords = dofToCoord[row];

          if (coordToDof.find(rowCoords) == coordToDof.end())
          {
            cout << "Cannot find row index for coords: " << endl;
            cout << rowCoords;
            exit(0);
          }

          DegreeOfFreedom rowHere = coordToDof[rowCoords];

          if (solValues[i].count(rowHere) == 0)
          {
            solValues[i][rowHere] = value;
          }
          else
          {
            double diff = fabs(solValues[i][rowHere] - value);
            if (diff > 1e-8)
            {
              cout << "DIFFERENT values in solution vector!" << endl;
              exit(0);
            }
          }
        }
      }

      in.close();

      cout << "File read!" << endl;
    }

    cout << "Now checking ..." << endl;

    cout.precision(15);

    for (int i = 0; i < nComponents; i++)
    {
      for (int j = 0; j < nComponents; j++)
      {
        DOFMatrix* dofmatrix = (*systemMatrix)[i][j];

        if (!dofmatrix)
          continue;

        for (map<pair<DegreeOfFreedom, DegreeOfFreedom>, double>::iterator nnzIt =
               nnzValues[i][j].begin(); nnzIt != nnzValues[i][j].end(); ++nnzIt)
        {

          DegreeOfFreedom row = nnzIt->first.first;
          DegreeOfFreedom col = nnzIt->first.second;
          double value = nnzIt->second;

          double valueHere = (dofmatrix->getBaseMatrix())[row][col];
          double diff = fabs(valueHere - value);

          // And so we check the difference between the value read from file and the
          // corresponding value in the matrix.
          if (diff > 1e-8)
          {
            cout << "Wrong value in component matrix " << i
                 << "/" << j << ": " << endl
                 << "  DOF in this matrix[" << row << "][" << col << "] = "
                 << valueHere << endl
                 << "  DOF in other matrix[" << dofMapHereToFiles[row][0] << "]["
                 << dofMapHereToFiles[col][0] << "] = "
                 << value << endl;

            exit(0);
          }
        }
      }
    }

    for (int i = 0; i < nComponents; i++)
    {
      for (map<DegreeOfFreedom, double>::iterator rhsIt = rhsValues[i].begin();
           rhsIt != rhsValues[i].end(); ++rhsIt)
      {

        DegreeOfFreedom row = rhsIt->first;
        double value = rhsIt->second;

        double valueHere = (*(rhs->getDOFVector(i)))[row];
        double diff = fabs(valueHere - value);

        if (diff > 1e-8)
        {
          cout << "WRONG value in rhs[" << i << "]!" << endl
               << "  DOF in other rhs[" << row << "] = " << value << endl
               << "  DOF in this rhs[" << row << "] = " << valueHere << endl;

          exit(0);
        }
      }
    }

    for (int i = 0; i < nComponents; i++)
    {
      for (map<DegreeOfFreedom, double>::iterator solIt = solValues[i].begin();
           solIt != solValues[i].end(); ++solIt)
      {

        DegreeOfFreedom row = solIt->first;
        double value = solIt->second;

        double valueHere = (*(solution->getDOFVector(i)))[row];
        double diff = fabs(valueHere - value);

        if (diff > 1e-8)
        {
          cout << "WRONG value in sol[" << i << "]!" << endl
               << "  DOF in other sol[" << row << "] = " << value << endl
               << "  DOF in this sol[" << row << "] = " << valueHere << endl;

          exit(0);
        }
      }
    }

    cout << "FINISHED COMPARING!" << endl;

    exit(0);
  }


  void ProblemStatDbg::createDofToCoordMap(DofToCoord& dofMap)
  {
    const BasisFunction* basisFcts = componentSpaces[0]->getBasisFcts();
    WorldVector<double> coords;
    ElementDofIterator elDofIter(componentSpaces[0], true);

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(componentMeshes[0], -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo)
    {
      elDofIter.reset(elInfo->getElement());
      int i = 0;
      do
      {
        DimVec<double>* baryCoord = basisFcts->getCoords(i);
        elInfo->coordToWorld(*baryCoord, coords);
        i++;
        dofMap[elDofIter.getDof()] = coords;
      }
      while (elDofIter.next());

      elInfo = stack.traverseNext(elInfo);
    }
  }


  void ProblemStatDbg::createCoordToDofMap(CoordToDof& dofMap)
  {
    const BasisFunction* basisFcts = componentSpaces[0]->getBasisFcts();
    WorldVector<double> coords;
    ElementDofIterator elDofIter(componentSpaces[0], true);

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(componentMeshes[0], -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo)
    {
      elDofIter.reset(elInfo->getElement());
      int i = 0;
      do
      {
        DimVec<double>* baryCoord = basisFcts->getCoords(i);
        elInfo->coordToWorld(*baryCoord, coords);
        i++;
        dofMap[coords] = elDofIter.getDof();
      }
      while (elDofIter.next());

      elInfo = stack.traverseNext(elInfo);
    }
  }

} // end namespace AMDiS
