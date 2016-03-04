#include <fstream>

#include "FixVec.hpp"
#include "DOFVector.hpp"
#include "BasisFunction.hpp"
#include "Mesh.hpp"

namespace AMDiS
{
  namespace io
  {

    namespace GridWriter
    {

      template<typename T>
      void writeGrid(const WorldVector<double>* p,
                     int* numPoints,
                     double* dist,
                     DOFVector<T>* vec,
                     const char* filename,
                     int outFilePrecision)
      {
        FUNCNAME("writeGrid()");

        TEST_EXIT(vec)("no dof vector\n");
        TEST_EXIT(filename)("no filename\n");
        TEST_EXIT(vec->getFeSpace()->getBasisFcts())("no basis fcts\n");

        std::ofstream outFile(filename);
        outFile.precision(outFilePrecision);
        outFile.setf(std::ios_base::scientific);
        const Mesh* mesh = vec->getFeSpace()->getMesh();
        int dim = mesh->getDim();

        TEST_EXIT(dim == Global::getGeo(WORLD))("not for DIM != DIM_OF_WORLD\n");

        WorldVector<double>* basis = new WorldVector<double>[dim];
        double* lengthBasis = new double[dim];
        WorldVector<double>* step = new WorldVector<double>[3];
        for (int i = 0; i < 3; i++)
          step[i].set(0.0);

        for (int i = 0; i < dim; i++)
        {
          TEST_EXIT(numPoints[i] > 0)("numPoints < 1\n");
          TEST_EXIT(dist[i] >= 0)("dist < 0\n");
          lengthBasis[i] = 0;
        }

        WorldVector<double> curCoord;
        DimVec<double> bary(dim);
        Element* elp;
        const BasisFunction* basFcts = vec->getFeSpace()->getBasisFcts();
        DenseVector<double> uhLoc(basFcts->getNumber());

        // get basis of grid
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            basis[i][j] = p[i + 1][j] ;
            lengthBasis[i] += basis[i][j] * basis[i][j];
          }
          lengthBasis[i] = std::sqrt(lengthBasis[i]);
        }

        // norm basis, get steps
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            basis[i][j] /= lengthBasis[i];
            step[i][j] = basis[i][j] * dist[j];
          }
        }

        /* write grid points */
        int localNumPoints[3] = {1, 1, 1};
        for (int i = 0; i < dim; i++)
          localNumPoints[i] = numPoints[i];

        // Warning "Coords not in mesh domain" should be printed at most one time.
        bool warning = false;

        for (int i = 0; i < localNumPoints[0]; i++)
        {
          for (int j = 0; j < localNumPoints[1]; j++)     // pseudo-loop for dim < 2
          {
            for (int k = 0; k < localNumPoints[2]; k++)   // pseudo-loop for dim < 3
            {
              // clac current coords
              for (int l = 0; l < dim; l++)
              {
                curCoord[l] = p[0][l]
                              + (i * step[0][l])
                              + (j * step[1][l])     // j = 0 for dim > 1
                              + (k * step[2][l]);    // k = 0 for dim > 2
              }

              int inside = (const_cast<Mesh*>(mesh))->findElementAtPoint(curCoord,
                           &elp,
                           bary,
                           NULL, NULL, NULL);

              // write coords
              for (int l = 0; l < dim; l++)
                outFile << curCoord[l] << " ";

              if (!inside)
              {
                if (!warning)
                {
                  WARNING("Coords not in mesh domain\n");
                  warning = true;
                }
                // write value
                outFile << "0.0" << std::endl;
              }
              else
              {
                // get value at coords
                vec->getLocalVector(elp, uhLoc);
                double value = basFcts->evalUh(bary, uhLoc);

                // write value
                outFile << value << std::endl;
              }
            }
            if (localNumPoints[2] > 1)
              outFile << std::endl;
          }
          if (localNumPoints[1] > 1)
            outFile << std::endl;
        }

        delete [] basis;
        delete [] lengthBasis;
        delete [] step;
      }


      template<typename T>
      void writeGrid(const WorldVector<double>* p,
                     int* numPoints,
                     double* dist,
                     std::vector<DOFVector<T>*> vec,
                     const char* filename,
                     int outFilePrecision)
      {
        FUNCNAME("writeGrid()");

        TEST_EXIT(vec.size() > 0)("no dof vector\n");
        TEST_EXIT(filename)("no filename\n");
        TEST_EXIT(vec[0]->getFeSpace()->getBasisFcts())("no basis fcts\n");

        std::ofstream outFile(filename);
        outFile.precision(outFilePrecision);
        outFile.setf(std::ios_base::scientific);
        const Mesh* mesh = vec[0]->getFeSpace()->getMesh();
        int dim = mesh->getDim();

        TEST_EXIT(dim == Global::getGeo(WORLD))("not for DIM != DIM_OF_WORLD\n");

        WorldVector<double>* basis = new WorldVector<double>[dim];
        double* lengthBasis = new double[dim];
        WorldVector<double>* step = new WorldVector<double>[3];
        for (int i = 0; i < 3; i++)
          step[i].set(0.0);

        for (int i = 0; i < dim; i++)
        {
          TEST_EXIT(numPoints[i] > 0)("numPoints < 1\n");
          TEST_EXIT(dist[i] >= 0)("dist < 0\n");
          lengthBasis[i] = 0;
        }

        WorldVector<double> curCoord;
        DimVec<double> bary(dim);
        Element* elp;
        const BasisFunction* basFcts = vec[0]->getFeSpace()->getBasisFcts();
        DenseVector<double> uhLoc(basFcts->getNumber());

        // get basis of grid
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            basis[i][j] = p[i + 1][j] ;
            lengthBasis[i] += basis[i][j] * basis[i][j];
          }
          lengthBasis[i] = std::sqrt(lengthBasis[i]);
        }

        // norm basis, get steps
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            basis[i][j] /= lengthBasis[i];
            step[i][j] = basis[i][j] * dist[j];
          }
        }

        /* write grid points */
        int localNumPoints[3] = {1, 1, 1};
        for (int i = 0; i < dim; i++)
          localNumPoints[i] = numPoints[i];

        // Warning "Coords not in mesh domain" should be printed at most one time.
        bool warning = false;

        for (int i = 0; i < localNumPoints[0]; i++)
        {
          for (int j = 0; j < localNumPoints[1]; j++)     // pseudo-loop for dim < 2
          {
            for (int k = 0; k < localNumPoints[2]; k++)   // pseudo-loop for dim < 3
            {
              // clac current coords
              for (int l = 0; l < dim; l++)
              {
                curCoord[l] = p[0][l]
                              + (i * step[0][l])
                              + (j * step[1][l])     // j = 0 for dim > 1
                              + (k * step[2][l]);    // k = 0 for dim > 2
              }

              int inside = (const_cast<Mesh*>(mesh))->findElementAtPoint(curCoord,
                           &elp,
                           bary,
                           NULL, NULL, NULL);

              // write coords
              for (int l = 0; l < dim; l++)
                outFile << curCoord[l] << " ";

              if (!inside)
              {
                if (!warning)
                {
                  WARNING("Coords not in mesh domain\n");
                  warning = true;
                }
                // write value
                outFile << "0.0" << std::endl;
              }
              else
              {
                // get value at coords
                double value;
                for (int rr = 0; rr < vec.size(); rr++)
                {
                  vec[rr]->getLocalVector(elp, uhLoc);
                  value = basFcts->evalUh(bary, uhLoc);
                  // write value
                  outFile << value <<"\t";
                }
                outFile << std::endl;
              }
            }
            if (localNumPoints[2] > 1)
              outFile << std::endl;
          }
          if (localNumPoints[1] > 1)
            outFile << std::endl;
        }

        delete [] basis;
        delete [] lengthBasis;
        delete [] step;
      }

    } // end namespace GridWriter
  }
} // end namespace io, AMDiS
