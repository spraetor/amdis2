#include "ElInfo2d.h"
#include "BasisFunction.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "FiniteElemSpace.h"
#include "Flag.h"
#include "MacroElement.h"
#include "Mesh.h"
#include "Global.h"
#include "FixVec.h"
#include "DOFVector.h"
#include "MatrixVectorOperations.h"

namespace AMDiS
{

  double ElInfo2d::mat_d1_left_val[3][3] = {{0.0, 1.0, 0.5},
    {0.0, 0.0, 0.5},
    {1.0, 0.0, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d1_left(mat_d1_left_val);


  double ElInfo2d::mat_d1_right_val[3][3] = {{0.0, 0.0, 0.5},
    {1.0, 0.0, 0.5},
    {0.0, 1.0, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d1_right(mat_d1_right_val);



  double ElInfo2d::mat_d2_left_val[6][6] = {{0.0, 1.0, 0.0, 0.375, -0.125, 0.0},
    {0.0, 0.0, 0.0, -0.125, -0.125, 0.0},
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.5, 1.0},
    {0.0, 0.0, 1.0, 0.75, 0.25, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d2_left(mat_d2_left_val);

  double ElInfo2d::mat_d2_right_val[6][6] = {{0.0, 0.0, 0.0, -0.125, -0.125, 0.0},
    {1.0, 0.0, 0.0, -0.125, 0.375, 0.0},
    {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.5, 0.0, 1.0},
    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.25, 0.75, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d2_right(mat_d2_right_val);



  double ElInfo2d::mat_d3_left_val[10][10] = {{0.0,  1.0, -6.25e-02,  3.125e-01,  0.0,  0.0,  6.25e-02,  0.0,  0.0, -6.25e-02},
    {0.0,  0.0, -6.25e-02,  6.25e-02,  0.0,  0.0,  6.25e-02,  0.0,  0.0, 6.25e-02},
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.5e-01,  0.0,  0.0, -0.125},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.5e-01,  0.0,  1.0, 0.375},
    {0.0,  0.0,  5.625e-01,  9.375e-01,  1.0,  0.0, -6.25e-02,  0.0,  0.0, 1.875e-01},
    {0.0,  0.0,  5.625e-01, -3.125e-01,  0.0,  0.0, -6.25e-02,  0.0,  0.0, -1.875e-01},
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0, 7.5e-01}
  };
  mtl::dense2D<double> ElInfo2d::mat_d3_left(mat_d3_left_val);

  double ElInfo2d::mat_d3_right_val[10][10] = {{0.0,  0.0, -6.25e-02,  6.25e-02,  0.0,  0.0,  6.25e-02,  0.0,  0.0, 6.25e-02},
    {1.0,  0.0, -6.25e-02,  6.25e-02,  0.0,  0.0,  3.125e-01,  0.0,  0.0, -6.25e-02},
    {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0, -2.5e-01,  0.0,  0.0,  0.0,  1.0,  0.0, 0.375},
    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0, -2.5e-01,  0.0,  0.0,  0.0,  0.0,  0.0, -0.125},
    {0.0,  0.0,  5.625e-01, -6.25e-02,  0.0,  0.0, -3.125e-01,  0.0,  0.0, -1.875e-01},
    {0.0,  0.0,  5.625e-01, -6.25e-02,  0.0,  1.0,  9.375e-01,  0.0,  0.0, 1.875e-01},
    {0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 7.5e-01}
  };
  mtl::dense2D<double> ElInfo2d::mat_d3_right(mat_d3_right_val);



  double ElInfo2d::mat_d4_left_val[15][15] = {{0.0,  1.0,  0.0,  2.734375e-01,  0.0, -3.906250e-02,  2.343750e-02,  0.0, -3.906250e-02,  0.0,  0.0,  0.0,  2.343750e-02, -3.906250e-02, 0.0},
    {0.0,  0.0,  0.0, -3.906250e-02,  0.0,  2.343750e-02,  2.343750e-02,  0.0, -3.906250e-02,  0.0,  0.0,  0.0, -3.906250e-02, -3.906250e-02, 0.0},
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -6.25e-02,  0.0,  1.875e-01,  0.0,  0.0,  0.0,  0.125,  6.25e-02, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.375,  0.0,  0.0,  0.0, -0.125,  0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.375,  0.0,  1.0,  0.0,  0.375,  0.0, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -6.25e-02,  0.0,  1.875e-01,  0.0,  0.0,  1.0, -0.125,  3.125e-01, 0.0},
    {0.0,  0.0,  0.0,  1.093750e+00,  1.0,  4.687500e-01, -9.375e-02,  0.0,  3.125e-02,  0.0,  0.0,  0.0, -3.125e-02,  1.562500e-01, 0.0},
    {0.0,  0.0,  1.0, -5.468750e-01,  0.0,  7.031250e-01,  1.406250e-01,  0.0,  1.562500e-02,  0.0,  0.0,  0.0, -4.687500e-02, -2.343750e-01, 0.0},
    {0.0,  0.0,  0.0,  2.187500e-01,  0.0, -1.562500e-01, -9.375e-02,  0.0,  3.125e-02,  0.0,  0.0,  0.0,  9.375e-02,  1.562500e-01, 0.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  5.625e-01,  0.0, -1.875e-01,  0.0,  0.0,  0.0,  0.375,  9.375e-01, 1.0},
    {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  5.625e-01,  0.0, -1.875e-01,  0.0,  0.0,  0.0, -0.375, -3.125e-01, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 7.5e-01, 0.0, 0.0, 0.0, 7.5e-01, 0.0, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d4_left(mat_d4_left_val);

  double ElInfo2d::mat_d4_right_val[15][15] = {{0.0,  0.0,  0.0, -3.906250e-02,  0.0,  2.343750e-02,  2.343750e-02,  0.0, -3.906250e-02,  0.0,  0.0,  0.0, -3.906250e-02, -3.906250e-02, 0.0},
    {1.0,  0.0,  0.0, -3.906250e-02,  0.0,  2.343750e-02, -3.906250e-02,  0.0,  2.734375e-01,  0.0,  0.0,  0.0, -3.906250e-02,  2.343750e-02, 0.0},
    {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0,  1.875e-01,  0.0, -6.25e-02,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  3.125e-01, -0.125, 0.0},
    {0.0,  0.0,  0.0, -0.375,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.375, 0.0},
    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0,  0.0,  0.0, -0.375,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.125, 0.0},
    {0.0,  0.0,  0.0,  1.875e-01,  0.0, -6.25e-02,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  6.25e-02,  0.125, 0.0},
    {0.0,  0.0,  0.0,  3.125e-02,  0.0, -9.375e-02, -1.562500e-01,  0.0,  2.187500e-01,  0.0,  0.0,  0.0,  1.562500e-01,  9.375e-02, 0.0},
    {0.0,  0.0,  1.0,  1.562500e-02,  0.0,  1.406250e-01,  7.031250e-01,  0.0, -5.468750e-01,  0.0,  0.0,  0.0, -2.343750e-01, -4.687500e-02, 0.0},
    {0.0,  0.0,  0.0,  3.125e-02,  0.0, -9.375e-02,  4.687500e-01,  1.0,  1.093750e+00,  0.0,  0.0,  0.0,  1.562500e-01, -3.125e-02, 0.0},
    {0.0,  0.0,  0.0, -1.875e-01,  0.0,  5.625e-01,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -3.125e-01, -0.375, 0.0},
    {0.0,  0.0,  0.0, -1.875e-01,  0.0,  5.625e-01,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  9.375e-01,  0.375, 1.0},
    {0.0, 0.0, 0.0, 7.5e-01, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.5e-01, 0.0}
  };
  mtl::dense2D<double> ElInfo2d::mat_d4_right(mat_d4_right_val);


  ElInfo2d::ElInfo2d(Mesh* aMesh)
    : ElInfo(aMesh)
  {}


  ElInfo2d::~ElInfo2d()
  {}


  void ElInfo2d::fillMacroInfo(const MacroElement* mel)
  {
    FUNCNAME("ElInfo::fillMacroInfo()");

    macroElement = const_cast<MacroElement*>(mel);
    element = const_cast<Element*>(mel->getElement());
    parent = NULL;
    level = 0;

    if (fillFlag.isSet(Mesh::FILL_COORDS) ||
        fillFlag.isSet(Mesh::FILL_DET)    ||
        fillFlag.isSet(Mesh::FILL_GRD_LAMBDA))
    {

      int vertices = mesh->getGeo(VERTEX);
      for (int i = 0; i < vertices; i++)
        coord[i] = mel->getCoord(i);
    }

    int neighbours = mesh->getGeo(NEIGH);

    if (fillFlag.isSet(Mesh::FILL_OPP_COORDS) ||
        fillFlag.isSet(Mesh::FILL_NEIGH))
    {

      bool fill_opp_coords = (fillFlag.isSet(Mesh::FILL_OPP_COORDS));

      for (int i = 0; i < neighbours; i++)
      {
        MacroElement* macroNeighbour = mel->getNeighbour(i);

        if (macroNeighbour)
        {
          neighbour[i] = macroNeighbour->getElement();
          Element* nb = const_cast<Element*>(neighbour[i]);

          int edgeNo = oppVertex[i] = mel->getOppVertex(i);


          if (nb->getFirstChild() && edgeNo != 2)
          {
            /*
            * Search for the next neighbour. In many cases, the neighbour element
            * may be refinemed in a way, such that there is no new vertex on the
            * common edge. This situation is shown in the following picture:
            *
            *               /|\
            *              / | \
            *             /  |  \
            *            /\  |   \
            *           /  \ |    \
            *          /    \|     \
            *          -------------
            *
            *            nb     el
            *
            * Note that we know (because of the last if statement), that the
            * neighbour element has children and the common edge is not the
            * refinement edge, which has always the number 2, of our element.
            */

            if (edgeNo == 0)
            {
              /*
              * The situation is as follows:
              *
              *          -------
              *          \    /|\
              *           \  / | \
              *            \/  |  \
              *             \  |   \
              *              \ |    \
              *               \|     \
              *                -------
              *
              *            nb     el
                    *
              * That means, the edge 0 of the same level neighbour is the common
              * edge, i.e., the direct neighbour is the second child of the same
              * level neighbour.
              */
              nb = neighbour[i] = nb->getSecondChild();
            }
            else
            {
              // The situation is as shown in the picture above. So the next
              // neighbour is the first child of the same level neighbour element.
              nb = neighbour[i] = nb->getFirstChild();
            }

            // In both cases the opp vertex number is 2, as one can see in the
            // pictures above.
            oppVertex[i] = 2;

            if (fill_opp_coords)
            {
              if (nb->isNewCoordSet())
              {
                oppCoord[i] = *(nb->getNewCoord());
              }
              else
              {
                // In both cases, that are shown in the pictures above, the opp
                // vertex of the neighbour edge is the midpoint of the vertex 0
                // and vertex 1 of the same level neighbour element.
                oppCoord[i] = (macroNeighbour->getCoord(0) +
                               macroNeighbour->getCoord(1)) * 0.5;
              }

              switch (i)
              {
              case 0:
                // The common edge is the edge 0 of this element.

                switch (edgeNo)
                {
                case 1:
                  neighbourCoord[i][0] = macroNeighbour->getCoord(2);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(0);
                  break;
                case 0:
                  neighbourCoord[i][0] = macroNeighbour->getCoord(1);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(2);
                  break;
                default:
                  ERROR_EXIT("Should not happen!\n");
                }

                neighbourCoord[i][2] = oppCoord[i];
                break;

              case 1:
                // The common edge is the edge 1 of this element.
                switch (edgeNo)
                {
                case 0:
                  neighbourCoord[i][0] = macroNeighbour->getCoord(1);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(2);
                  break;
                case 1:
                  neighbourCoord[i][0] = macroNeighbour->getCoord(2);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(0);
                  break;
                default:
                  ERROR_EXIT("Should not happen!\n");
                }

                neighbourCoord[i][2] = oppCoord[i];
                break;

              case 2:
                if (*(macroNeighbour->getElement()->getDof(2)) == *(element->getDof(0)))
                {
                  neighbourCoord[i][0] = macroNeighbour->getCoord(2);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(1);
                }
                else if (*(macroNeighbour->getElement()->getDof(2)) == *(element->getDof(1)))
                {
                  neighbourCoord[i][0] = macroNeighbour->getCoord(0);
                  neighbourCoord[i][1] = macroNeighbour->getCoord(2);
                }
                else
                {
                  ERROR_EXIT("Should not happen! Non-conforming AMDiS-mesh? Periodic mesh with corrected index-circles?\n");
                }

                // I've deleted here some code, be I think that this case is not
                // possible. If an error occurs in this line, please check AMDiS
                // revision <= 476 at the same position.
                //		ERROR_EXIT("Should not happen!\n");

                break;

              default:
                std::cout << "------------- Error --------------" << std::endl;
                std::cout << "  Neighbour counter = " << i << "\n";
                std::cout << "  Element index     = " << element->getIndex() << "\n\n";
                for (int j = 0; j < neighbours; j++)
                {
                  if (mel->getNeighbour(j))
                  {
                    std::cout << "  Neighbour " << j << ": "
                              << mel->getNeighbour(j)->getElement()->getIndex()
                              << std::endl;
                  }
                  else
                  {
                    std::cout << "  Neighbour " << j << ": not existing" << std::endl;
                  }
                  std::cout << "  OppVertex " << j << ": "
                            << static_cast<int>(mel->getOppVertex(j))
                            << std::endl << std::endl;
                }
                ERROR_EXIT("should not happen!\n");
                break;
              }
            }
          }
          else
          {

            // In this case, we know that the common edge is the refinement edge.
            // This makes everything much more simpler, because we know that the
            // next neighbour is equal to the samel level neighbour. If the same
            // level neighbour would be refinement, also this element must to be
            // refinement, because they share the refinement edge.

            if (fill_opp_coords)
            {
              oppCoord[i] = macroNeighbour->getCoord(edgeNo);
              neighbourCoord[i] = macroNeighbour->getCoord();
            }
          }
        }
        else
        {
          neighbour[i] = NULL;
        }
      }
    }

    if (fillFlag.isSet(Mesh::FILL_BOUND))
    {
      for (int i = 0; i < element->getGeo(BOUNDARY); i++)
        boundary[i] = mel->getBoundary(i);
      for (int i = 0; i < element->getGeo(PROJECTION); i++)
        projection[i] = mel->getProjection(i);
    }
  }


  /****************************************************************************/
  /*   fill ElInfo structure for one child of an element   		    */
  /****************************************************************************/

  void ElInfo2d::fillElInfo(int ichild, const ElInfo* elInfoOld)
  {
    FUNCNAME("ElInfo::fillElInfo()");

    Element* elem = elInfoOld->getElement();
    Flag fill_flag = elInfoOld->getFillFlag();

    TEST_EXIT_DBG(elem->getFirstChild())("no children?\n");
    element = const_cast<Element*>((ichild == 0) ?
                                   elem->getFirstChild() :
                                   elem->getSecondChild());
    TEST_EXIT_DBG(element)("missing child %d?\n", ichild);

    macroElement  = elInfoOld->getMacroElement();
    fillFlag = fill_flag;
    parent = elem;
    level = elInfoOld->getLevel() + 1;
    iChild = ichild;

    if (fillFlag.isSet(Mesh::FILL_COORDS) ||
        fillFlag.isSet(Mesh::FILL_DET)    ||
        fillFlag.isSet(Mesh::FILL_GRD_LAMBDA))
    {

      if (elem->isNewCoordSet())
        coord[2] = *(elem->getNewCoord());
      else
        coord[2].setMidpoint(elInfoOld->getCoord(0), elInfoOld->getCoord(1));

      if (ichild == 0)
      {
        coord[0] = elInfoOld->getCoord(2);
        coord[1] = elInfoOld->getCoord(0);
      }
      else
      {
        coord[0] = elInfoOld->getCoord(1);
        coord[1] = elInfoOld->getCoord(2);
      }
    }

    bool fill_opp_coords = (fill_flag.isSet(Mesh::FILL_OPP_COORDS));

    if (fill_flag.isSet(Mesh::FILL_NEIGH) || fill_opp_coords)
    {
      if (ichild == 0)
      {
        // Calculation of the neighbour 2, its oppCoords and the
        // cooresponding oppVertex.

        neighbour[2] = elInfoOld->getNeighbour(1);
        oppVertex[2] = elInfoOld->getOppVertex(1);

        if (neighbour[2] && fill_opp_coords)
        {
          oppCoord[2] = elInfoOld->getOppCoord(1);
          neighbourCoord[2] = elInfoOld->getNeighbourCoord(1);
        }


        // Calculation of the neighbour 1, its oppCoords and the
        // cooresponding oppVertex.

        if (elem->getFirstChild()  &&
            elem->getSecondChild()->getFirstChild()  &&
            elem->getSecondChild()->getFirstChild())
        {

          neighbour[1] = elem->getSecondChild()->getSecondChild();
          oppVertex[1] = 2;

          if (fill_opp_coords)
          {
            if (elem->getSecondChild()->isNewCoordSet())
              oppCoord[1] = *(elem->getSecondChild()->getNewCoord());
            else
              oppCoord[1].setMidpoint(elInfoOld->getCoord(1), elInfoOld->getCoord(2));

            neighbourCoord[1][0] = coord[0];
            neighbourCoord[1][1] = coord[2];
            neighbourCoord[1][2] = oppCoord[1];
          }
        }
        else
        {
          neighbour[1] = elem->getSecondChild();
          oppVertex[1] = 0;

          if (fill_opp_coords)
          {
            oppCoord[1] = elInfoOld->getCoord(1);

            neighbourCoord[1][0] = elInfoOld->getCoord(1);
            neighbourCoord[1][1] = elInfoOld->getCoord(2);
            neighbourCoord[1][2] = coord[2];
          }
        }


        // Calculation of the neighbour 0, its oppCoords and the
        // cooresponding oppVertex.

        Element* nb = elInfoOld->getNeighbour(2);
        if (nb)
        {
          TEST_EXIT_DBG(elInfoOld->getOppVertex(2) == 2)
          ("Fill child %d of element %d (mel %d): Invalid neighbour %d!\n",
           ichild,
           elInfoOld->getElement()->getIndex(),
           elInfoOld->getMacroElement()->getIndex(),
           nb->getIndex());

          TEST_EXIT_DBG(nb->getFirstChild())
          ("Missing first child in element %d!\n", nb->getIndex());
          TEST_EXIT_DBG(nb->getSecondChild())
          ("Missing second child in element %d!\n", nb->getIndex());

          nb = nb->getSecondChild();

          if (nb->getFirstChild())
          {
            oppVertex[0] = 2;

            if (fill_opp_coords)
            {
              if (nb->isNewCoordSet())
              {
                oppCoord[0] = *(nb->getNewCoord());
              }
              else
              {
                oppCoord[0].setMidpoint(elInfoOld->getNeighbourCoord(2)[1],
                                        elInfoOld->getNeighbourCoord(2)[2]);
              }

              neighbourCoord[0][0].setMidpoint(elInfoOld->getNeighbourCoord(2)[0],
                                               elInfoOld->getNeighbourCoord(2)[1]);
              neighbourCoord[0][1] = elInfoOld->getNeighbourCoord(2)[1];
              neighbourCoord[0][2] = oppCoord[0];
            }

            nb = nb->getFirstChild();
          }
          else
          {
            oppVertex[0] = 1;

            if (fill_opp_coords)
            {
              oppCoord[0] = elInfoOld->getOppCoord(2);

              neighbourCoord[0][0] = elInfoOld->getNeighbourCoord(2)[0];
              neighbourCoord[0][1] = elInfoOld->getNeighbourCoord(2)[2];
              neighbourCoord[0][2].setMidpoint(elInfoOld->getNeighbourCoord(2)[0],
                                               elInfoOld->getNeighbourCoord(2)[1]);
            }
          }
        }

        neighbour[0] = nb;
      }
      else       /* ichild == 1 */
      {
        // Calculation of the neighbour 2, its oppCoords and the
        // cooresponding oppVertex.

        neighbour[2] = elInfoOld->getNeighbour(0);
        oppVertex[2] = elInfoOld->getOppVertex(0);

        if (neighbour[2] && fill_opp_coords)
        {
          oppCoord[2] = elInfoOld->getOppCoord(0);
          neighbourCoord[2] = elInfoOld->getNeighbourCoord(0);
        }


        // Calculation of the neighbour 0, its oppCoords and the
        // cooresponding oppVertex.

        if (elem->getFirstChild()->getFirstChild())
        {
          neighbour[0] = elem->getFirstChild()->getFirstChild();
          oppVertex[0] = 2;

          if (fill_opp_coords)
          {
            if (elem->getFirstChild()->isNewCoordSet())
            {
              oppCoord[0] = *(elem->getFirstChild()->getNewCoord());
            }
            else
            {
              oppCoord[0].setMidpoint(elInfoOld->getCoord(0),
                                      elInfoOld->getCoord(2));
            }

            neighbourCoord[0][0] = coord[2];
            neighbourCoord[0][1] = coord[1];
            neighbourCoord[0][2] = oppCoord[0];
          }
        }
        else
        {
          neighbour[0] = elem->getFirstChild();
          oppVertex[0] = 1;

          if (fill_opp_coords)
          {
            oppCoord[0] = elInfoOld->getCoord(0);

            neighbourCoord[0][0] = elInfoOld->getCoord(2);
            neighbourCoord[0][1] = elInfoOld->getCoord(0);
            neighbourCoord[0][2] = coord[2];
          }
        }

        // Calculation of the neighbour 1, its oppCoords and the
        // cooresponding oppVertex.

        Element* nb = elInfoOld->getNeighbour(2);
        if (nb)
        {
          TEST(elInfoOld->getOppVertex(2) == 2)("invalid neighbour\n");
          TEST((nb = nb->getFirstChild()))("missing child?\n");

          if (nb->getFirstChild())
          {
            oppVertex[1] = 2;

            if (fill_opp_coords)
            {
              if (nb->isNewCoordSet())
              {
                oppCoord[1] = *(nb->getNewCoord());
              }
              else
              {
                oppCoord[1].setMidpoint(elInfoOld->getNeighbourCoord(2)[0],
                                        elInfoOld->getNeighbourCoord(2)[2]);
              }

              neighbourCoord[1][0] = elInfoOld->getNeighbourCoord(2)[0];
              neighbourCoord[1][1].setMidpoint(elInfoOld->getNeighbourCoord(2)[0],
                                               elInfoOld->getNeighbourCoord(2)[1]);
              neighbourCoord[1][2] = oppCoord[1];
            }
            nb = nb->getSecondChild();

          }
          else
          {
            oppVertex[1] = 0;

            if (fill_opp_coords)
            {
              oppCoord[1] = elInfoOld->getOppCoord(2);

              neighbourCoord[1][0] = elInfoOld->getNeighbourCoord(2)[2];
              neighbourCoord[1][1] = elInfoOld->getNeighbourCoord(2)[0];
              neighbourCoord[1][2].setMidpoint(elInfoOld->getNeighbourCoord(2)[0],
                                               elInfoOld->getNeighbourCoord(2)[1]);
            }
          }
        }
        neighbour[1] = nb;
      } // if (ichild == 0) {} else
    } // if (fill_flag.isSet(Mesh::FILL_NEIGH) || fillFlag.isSet(Mesh::FILL_OPP_COORDS))


    if (fill_flag.isSet(Mesh::FILL_BOUND))
    {
      if (elInfoOld->getBoundary(2))
        boundary[5] = elInfoOld->getBoundary(2);
      else
        boundary[5] = INTERIOR;

      if (ichild == 0)
      {
        boundary[3] = elInfoOld->getBoundary(5);
        boundary[4] = elInfoOld->getBoundary(3);
        boundary[0] = elInfoOld->getBoundary(2);
        boundary[1] = INTERIOR;
        boundary[2] = elInfoOld->getBoundary(1);
      }
      else
      {
        boundary[3] = elInfoOld->getBoundary(4);
        boundary[4] = elInfoOld->getBoundary(5);
        boundary[0] = INTERIOR;
        boundary[1] = elInfoOld->getBoundary(2);
        boundary[2] = elInfoOld->getBoundary(0);
      }

      if (elInfoOld->getProjection(0) &&
          elInfoOld->getProjection(0)->getType() == VOLUME_PROJECTION)
      {

        projection[0] = elInfoOld->getProjection(0);
      }
      else     // boundary projection
      {
        if (ichild == 0)
        {
          projection[0] = elInfoOld->getProjection(2);
          projection[1] = NULL;
          projection[2] = elInfoOld->getProjection(1);
        }
        else
        {
          projection[0] = NULL;
          projection[1] = elInfoOld->getProjection(2);
          projection[2] = elInfoOld->getProjection(0);
        }
      }
    }
  }


  double ElInfo2d::calcGrdLambda(DimVec<WorldVector<double>>& grd)
  {
    FUNCNAME("ElInfo2d::calcGrdLambda()");

    testFlag(Mesh::FILL_COORDS);

    double adet = 0.0;
    int dim = mesh->getDim();

    e1 = coord[1] - coord[0];
    e2 = coord[2] - coord[0];

    if (dimOfWorld == 2)
    {
      double sdet = e1[0] * e2[1] - e1[1] * e2[0];
      adet = math::abs(sdet);

      if (adet < 1.0E-25)
      {
        MSG("abs(det) = %f\n", adet);
        for (int i = 0; i <= dim; i++)
          grd[i].set(0.0);
      }
      else
      {
        double det1 = 1.0 / sdet;

        grd[1][0] = e2[1] * det1;  // a11: (a_ij) = A^{-T}
        grd[1][1] = -e2[0] * det1; // a21
        grd[2][0] = -e1[1] * det1; // a12
        grd[2][1] = e1[0] * det1;  // a22
        grd[0][0] = -grd[1][0] - grd[2][0];
        grd[0][1] = -grd[1][1] - grd[2][1];
      }
    }
    else
    {
      normal = cross(e1, e2);

      adet = norm(normal);

      if (adet < 1.0E-15)
      {
        MSG("abs(det) = %lf\n", adet);
        for (int i = 0; i <= dim; i++)
          for (int j = 0; j < dimOfWorld; j++)
            grd[i][j] = 0.0;
      }
      else
      {
        grd[1] = cross(e2, normal);
        grd[2] = cross(normal, e1);

        double adet2 = 1.0 / (adet * adet);

        for (int i = 0; i < dimOfWorld; i++)
        {
          grd[1][i] *= adet2;
          grd[2][i] *= adet2;
        }

        grd[0][0] = -grd[1][0] - grd[2][0];
        grd[0][1] = -grd[1][1] - grd[2][1];
        grd[0][2] = -grd[1][2] - grd[2][2];
      }
    }

    return adet;
  }


  int ElInfo2d::worldToCoord(const WorldVector<double>& xy,
                             DimVec<double>& lambda) const
  {
    FUNCNAME("ElInfo::worldToCoord()");

    DimVec<WorldVector<double>> edge(mesh->getDim());
    WorldVector<double> x;
    static DimVec<double> vec(mesh->getDim());

    int dim = mesh->getDim();

    for (int j = 0; j < dimOfWorld; j++)
    {
      double x0 = coord[dim][j];
      x[j] = xy[j] - x0;
      for (int i = 0; i < dim; i++)
        edge[i][j] = coord[i][j] - x0;
    }

    double det  = edge[0][0] * edge[1][1] - edge[0][1] * edge[1][0];
    double det0 =       x[0] * edge[1][1] -       x[1] * edge[1][0];
    double det1 = edge[0][0] * x[1]       - edge[0][1] * x[0];

    if (math::abs(det) < DBL_TOL)
    {
      ERROR("det = %le; abort\n", det);
      for (int i = 0; i <= dim; i++)
        lambda[i] = 1.0 / dim;
      return 0;
    }


    lambda[0] = det0 / det;
    lambda[1] = det1 / det;
    lambda[2] = 1.0 - lambda[0] - lambda[1];

    int k = -1;
    double lmin = 0.0;
    for (int i = 0; i <= dim; i++)
    {
      if (lambda[i] < -1.e-5)
      {
        if (lambda[i] < lmin)
        {
          k = i;
          lmin = lambda[i];
        }
      }
    }

    return k;
  }


  double ElInfo2d::getNormal(int side, WorldVector<double>& normal) const
  {
    FUNCNAME_DBG("ElInfo::getNormal()");

    int i0 = (side + 1) % 3;
    int i1 = (side + 2) % 3;

    if (dimOfWorld == 2)
    {
      normal[0] = coord[i1][1] - coord[i0][1];
      normal[1] = coord[i0][0] - coord[i1][0];
    }
    else     // dow == 3
    {
      WorldVector<double> e0, e1, e2, elementNormal;

      e0 = coord[i1];
      e0 -= coord[i0];
      e1 = coord[i1];
      e1 -= coord[side];
      e2 = coord[i0];
      e2 -= coord[side];

      elementNormal = cross(e1, e2);
      normal = cross(elementNormal, e0);
    }

    double detn = norm(normal);

    TEST_EXIT_DBG(detn > 1.e-30)("det = 0 on face %d\n", side);

    normal *= 1.0 / detn;

    return detn;
  }


  /****************************************************************************/
  /*  calculate the normal of the element for dim of world = dim + 1          */
  /*  return the absulute value of the determinant from the                   */
  /*  transformation to the reference element                                 */
  /****************************************************************************/
  double ElInfo2d::getElementNormal(WorldVector<double>& elementNormal) const
  {
    FUNCNAME_DBG("ElInfo::getElementNormal()");

    TEST_EXIT_DBG(dimOfWorld == 3)
    (" element normal only well defined for  DIM_OF_WORLD = DIM + 1 !!");

    WorldVector<double> e0 = coord[1] - coord[0];
    WorldVector<double> e1 = coord[2] - coord[0];

    elementNormal = cross(e0, e1);

    double detn = norm(elementNormal);

    TEST_EXIT_DBG(detn > 1.e-30)("det = 0");

    elementNormal *= 1.0 / detn;

    return detn;
  }


  mtl::dense2D<double>& ElInfo2d::getSubElemCoordsMat(int degree) const
  {
    FUNCNAME("ElInfo2d::getSubElemCoordsMat()");

    using namespace mtl;

    if (subElemMatrices[degree].count(std::make_pair(refinementPathLength, refinementPath)) == 0)
    {
      switch (degree)
      {
      case 1:
      {
        dense2D<double> mat(3, 3), tmpMat(3, 3);
        mat = 1;

        for (int i = 0; i < refinementPathLength; i++)
        {
          if (refinementPath & (1 << i))
          {
            tmpMat = mat * mat_d1_right;
            mat = tmpMat;
          }
          else
          {
            tmpMat = mat * mat_d1_left;
            mat = tmpMat;
          }
        }

        subElemMatrices[1][std::make_pair(refinementPathLength, refinementPath)] = mat;
      }
      break;
      case 2:
      {
        dense2D<double> mat(6, 6), tmpMat(6, 6);
        mat = 1;

        for (int i = 0; i < refinementPathLength; i++)
        {
          if (refinementPath & (1 << i))
          {
            tmpMat = mat * mat_d2_right;
            mat = tmpMat;
          }
          else
          {
            tmpMat = mat * mat_d2_left;
            mat = tmpMat;
          }
        }

        subElemMatrices[2][std::make_pair(refinementPathLength, refinementPath)] = mat;
      }
      break;
      case 3:
      {
        dense2D<double> mat(10, 10), tmpMat(10, 10);
        mat = 1;

        for (int i = 0; i < refinementPathLength; i++)
        {
          if (refinementPath & (1 << i))
          {
            tmpMat = mat * mat_d3_right;
            mat = tmpMat;
          }
          else
          {
            tmpMat = mat * mat_d3_left;
            mat = tmpMat;
          }
        }

        subElemMatrices[3][std::make_pair(refinementPathLength, refinementPath)] = mat;
      }
      break;
      case 4:
      {
        dense2D<double> mat(15, 15), tmpMat(15, 15);
        mat = 1;

        for (int i = 0; i < refinementPathLength; i++)
        {
          if (refinementPath & (1 << i))
          {
            tmpMat = mat * mat_d4_right;
            mat = tmpMat;
          }
          else
          {
            tmpMat = mat * mat_d4_left;
            mat = tmpMat;
          }
        }

        subElemMatrices[4][std::make_pair(refinementPathLength, refinementPath)] = mat;
      }
      break;

      default:
        ERROR_EXIT("Not supported for basis function degree: %d\n", degree);
      }
    }

    return subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)];
  }

} // end namespace AMDiS
