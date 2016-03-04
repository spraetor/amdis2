#include "est/Estimator.hpp"

#include "Traverse.hpp"
#include "Initfile.hpp"
// #include "DualTraverse.hpp"
#include "DOFVector.hpp"

namespace AMDiS
{

  Estimator::Estimator(std::string name_, int r)
    : name(name_),
      row(r)
  {
    Parameters::get(name + "->error norm", norm);
  }


  double Estimator::estimate(double ts)
  {
    //     FUNCNAME("Estimator::estimate()");

    bool dualTraverse = false;

    mesh = uh[row == -1 ? 0 : row]->getFeSpace()->getMesh();
    auxMesh = NULL;

    init(ts);

    if (!dualTraverse)
      singleMeshTraverse();
    else
      dualMeshTraverse();

    exit();

    return est_sum;
  }


  void Estimator::singleMeshTraverse()
  {
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
    while (elInfo)
    {
      estimateElement(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }
  }


  void Estimator::dualMeshTraverse()
  {
    ERROR_EXIT("Not supportet any more!\n");
//     DualTraverse dualTraverse;
//     DualElInfo dualElInfo;
//
//     bool cont = dualTraverse.traverseFirst(mesh, auxMesh, -1, -1,
//                                            traverseFlag, traverseFlag,
//                                            dualElInfo);
//     while (cont)
//     {
//       estimateElement(dualElInfo.rowElInfo, &dualElInfo);
//       cont = dualTraverse.traverseNext(dualElInfo);
//     }
  }
}
