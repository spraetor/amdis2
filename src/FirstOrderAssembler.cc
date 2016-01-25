#include <vector>

// #include <boost/numeric/mtl/mtl.hpp>

#include "Assembler.h"
#include "FirstOrderAssembler.h"
#include "Operator.h"
#include "QPsiPhi.h"
#include "FiniteElemSpace.h"
#include "Quadrature.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS
{
  ThreadPrivate<vector<SubAssembler*>>
                                    FirstOrderAssembler::optimizedSubAssemblersGrdPhi;
  ThreadPrivate<vector<SubAssembler*>>
                                    FirstOrderAssembler::optimizedSubAssemblersGrdPsi;

  ThreadPrivate<vector<SubAssembler*>>
                                    FirstOrderAssembler::standardSubAssemblersGrdPhi;
  ThreadPrivate<vector<SubAssembler*>>
                                    FirstOrderAssembler::standardSubAssemblersGrdPsi;


  FirstOrderAssembler::FirstOrderAssembler(Operator* op,
      Assembler* assembler,
      Quadrature* quad,
      bool optimized,
      FirstOrderType type)
    : SubAssembler(op, assembler, quad, 1, optimized, type)
  {
    FUNCNAME_DBG("FirstOrderAssmebler::FirstOrderAssembler()");

    TEST_EXIT_DBG(dim > 0)("Should not happen!\n");

    Lb.resize(1);
    Lb[0].change_dim(dim + 1);
  }


  FirstOrderAssembler* FirstOrderAssembler::getSubAssembler(Operator* op,
      Assembler* assembler,
      Quadrature* quad,
      FirstOrderType type,
      bool optimized)
  {
    vector<SubAssembler*>& subAssemblers =
      optimized
      ? (type == GRD_PSI
         ? optimizedSubAssemblersGrdPsi.get()
         : optimizedSubAssemblersGrdPhi.get())
      : (type == GRD_PSI
         ? standardSubAssemblersGrdPsi.get()
         : standardSubAssemblersGrdPhi.get());

    vector<OperatorTerm*> opTerms =
      (type == GRD_PSI) ? op->getFirstOrderGrdPsi() : op->getFirstOrderGrdPhi();

    // check if a assembler is needed at all
    if (opTerms.empty())
      return NULL;

    sort(opTerms.begin(), opTerms.end());

    // check if a new assembler is needed
    for (SubAssembler* subAssembler : subAssemblers)
    {
      vector<OperatorTerm*> assTerms = *(subAssembler->getTerms());

      sort(assTerms.begin(), assTerms.end());

      if (opTerms == assTerms && subAssembler->getQuadrature() == quad)
        return dynamic_cast<FirstOrderAssembler*>(subAssembler);
    }

    // check if all terms are pw_const
    bool pwConst = std::all_of(opTerms.begin(), opTerms.end(),
                               [](OperatorTerm* term)
    {
      return term->isPWConst();
    });

    // create new assembler
    FirstOrderAssembler* newAssembler;
    if (!optimized)
    {
      newAssembler = (type == GRD_PSI)
                     ? dynamic_cast<FirstOrderAssembler*>(new Stand10(op, assembler, quad))
                     : dynamic_cast<FirstOrderAssembler*>(new Stand01(op, assembler, quad));
    }
    else
    {
      if (pwConst)
      {
        newAssembler = (type == GRD_PSI)
                       ? dynamic_cast<FirstOrderAssembler*>(new Pre10(op, assembler, quad))
                       : dynamic_cast<FirstOrderAssembler*>(new Pre01(op, assembler, quad));
      }
      else
      {
        newAssembler = (type == GRD_PSI)
                       ? dynamic_cast<FirstOrderAssembler*>(new Quad10(op, assembler, quad))
                       : dynamic_cast<FirstOrderAssembler*>(new Quad01(op, assembler, quad));
      }
    }

    subAssemblers.push_back(newAssembler);
    return newAssembler;
  }


  Stand10::Stand10(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, false, GRD_PSI)
  {
    name = "standard first order assembler";

    psi = rowFeSpace->getBasisFcts();
    phi = colFeSpace->getBasisFcts();
  }


  void Stand10::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    DenseVector<double> grdPsi(dim + 1, 0.0);
    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);
    vector<double> phival(nCol);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }
    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq] *= elInfo->getDet();

      for (int i = 0; i < nCol; i++)
        phival[i] = (*(phi->getPhi(i)))(quadrature->getLambda(iq));

      for (int i = 0; i < nRow; i++)
      {
        (*(psi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPsi);
        for (int j = 0; j < nCol; j++)
          mat[i][j] += quadrature->getWeight(iq) * phival[j] * dot(Lb[iq], grdPsi);
      }
    }
  }


  void Stand10::calculateElementVectorImpl(const ElInfo* elInfo, DenseVector<double>& vec)
  {
    DenseVector<double> grdPsi(dim + 1, 0.0);
    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq] *= elInfo->getDet();

      for (int i = 0; i < nRow; i++)
      {
        (*(psi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPsi);
        vec[i] += quadrature->getWeight(iq) * dot(Lb[iq], grdPsi);
      }
    }
  }


  Quad10::Quad10(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, true, GRD_PSI)
  {
    name = "fast quadrature first order assembler";
  }


  void Quad10::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    if (firstCall)
    {
      const BasisFunction* basFcts = rowFeSpace->getBasisFcts();
      psiFast = updateFastQuadrature(psiFast, basFcts, INIT_GRD_PHI);
      basFcts = colFeSpace->getBasisFcts();
      phiFast = updateFastQuadrature(phiFast, basFcts, INIT_PHI);
      firstCall = false;
    }

    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);
    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    const mtl::dense2D<double>& phi = phiFast->getPhi();

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq] *= elInfo->getDet();

      const vector<DenseVector<double>>& grdPsi = psiFast->getGradient(iq);
      double factor = quadrature->getWeight(iq);

      for (int i = 0; i < nRow; i++)
      {
        for (int j = 0; j < nCol; j++)
          mat[i][j] += factor * phi[iq][j] * dot(Lb[iq], grdPsi[i]);
      }
    }
  }


  void Quad10::calculateElementVectorImpl(const ElInfo* elInfo, DenseVector<double>& vec)
  {
    if (firstCall)
    {
      const BasisFunction* basFcts = rowFeSpace->getBasisFcts();
      psiFast = updateFastQuadrature(psiFast, basFcts, INIT_GRD_PHI);
      basFcts = colFeSpace->getBasisFcts();
      phiFast = updateFastQuadrature(phiFast, basFcts, INIT_PHI);
      firstCall = false;
    }

    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);
    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq] *= elInfo->getDet();
      const vector<DenseVector<double>>& grdPsi = psiFast->getGradient(iq);

      for (int i = 0; i < nRow; i++)
        vec[i] += quadrature->getWeight(iq) * dot(Lb[iq], grdPsi[i]);
    }
  }


  Pre10::Pre10(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, true, GRD_PSI)
  {
    name = "precalculated first order assembler";
  }


  void Pre10::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    if (firstCall)
    {
      q10 = Q10PsiPhi::provideQ10PsiPhi(rowFeSpace->getBasisFcts(),
                                        colFeSpace->getBasisFcts(),
                                        quadrature);
      q1 = Q1Psi::provideQ1Psi(rowFeSpace->getBasisFcts(), quadrature);
      firstCall = false;
    }

    const int** nEntries = q10->getNumberEntries();
    // Do not need do resize Lb, because it's size is always at least one.
    Lb[0] = 0.0;

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    Lb[0] *= elInfo->getDet();

    for (int i = 0; i < nRow; i++)
    {
      for (int j = 0; j < nCol; j++)
      {
        const int* k = q10->getKVec(i, j);
        const double* values = q10->getValVec(i, j);
        double val = 0.0;
        for (int m = 0; m < nEntries[i][j]; m++)
          val += values[m] * Lb[0][k[m]];
        mat[i][j] += val;
      }
    }
  }


  Stand01::Stand01(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, false, GRD_PHI)
  {
    FUNCNAME_DBG("Stand01::Stand01()");

    TEST_EXIT_DBG(dim > 0)("Should not happen!\n");

    grdPhi.resize(nCol);
    for (int i = 0; i < nCol; i++)
      grdPhi[i].change_dim(dim + 1);

    psi = rowFeSpace->getBasisFcts();
    phi = colFeSpace->getBasisFcts();
  }


  void Stand01::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);
    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq] *= elInfo->getDet();

      for (int i = 0; i < nCol; i++)
        (*(phi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPhi[i]);

      for (int i = 0; i < nRow; i++)
      {
        double psival = (*(psi->getPhi(i)))(quadrature->getLambda(iq));
        for (int j = 0; j < nCol; j++)
          mat[i][j] += quadrature->getWeight(iq) * psival * dot(Lb[iq], grdPhi[j]);
      }
    }
  }


  Quad01::Quad01(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, true, GRD_PHI)
  {
    name = "fast quadrature first order assembler";
  }


  void Quad01::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    if (firstCall)
    {
      const BasisFunction* basFcts = rowFeSpace->getBasisFcts();
      psiFast = updateFastQuadrature(psiFast, basFcts, INIT_PHI);
      basFcts = colFeSpace->getBasisFcts();
      phiFast = updateFastQuadrature(phiFast, basFcts, INIT_GRD_PHI);
      firstCall = false;
    }

    int nPoints = quadrature->getNumPoints();
    Lb.resize(nPoints);
    for (int iq = 0; iq < nPoints; iq++)
    {
      Lb[iq].change_dim(dim + 1);
      Lb[iq] = 0.0;
    }

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    const mtl::dense2D<double>& psi = psiFast->getPhi();

    mtl::dense2D<double> facMat(nPoints, nCol);
    for (int iq = 0; iq < nPoints; iq++)
    {
      double weight = quadrature->getWeight(iq);
      DenseVector<double>& Lb_iq = Lb[iq];
      const vector<DenseVector<double>>& grdPsi = phiFast->getGradient(iq);

      Lb_iq *= elInfo->getDet();

      for (int i = 0; i < nCol; i++)
        facMat[iq][i] = weight * dot(Lb_iq, grdPsi[i]);
    }

    mat += trans(psi) * facMat;
  }


  Pre01::Pre01(Operator* op, Assembler* assembler, Quadrature* quad)
    : FirstOrderAssembler(op, assembler, quad, true, GRD_PHI)
  {
    name = "precalculated first order assembler";
  }


  void Pre01::calculateElementMatrixImpl(const ElInfo* elInfo, ElementMatrix& mat)
  {
    if (firstCall)
    {
      q01 = Q01PsiPhi::provideQ01PsiPhi(rowFeSpace->getBasisFcts(),
                                        colFeSpace->getBasisFcts(),
                                        quadrature);
      q1 = Q1Psi::provideQ1Psi(rowFeSpace->getBasisFcts(), quadrature);
      firstCall = false;
    }

    const int** nEntries = q01->getNumberEntries();
    // Do not need to resize Lb, because it's size is always at least one!
    Lb[0] = 0.0;

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    Lb[0] *= elInfo->getDet();

    for (int i = 0; i < nRow; i++)
    {
      for (int j = 0; j < nCol; j++)
      {
        const int* l = q01->getLVec(i, j);
        const double* values = q01->getValVec(i, j);
        double val = 0.0;
        for (int m = 0; m < nEntries[i][j]; m++)
          val += values[m] * Lb[0][l[m]];
        mat[i][j] += val;
      }
    }
  }


  void Pre10::calculateElementVectorImpl(const ElInfo* elInfo, DenseVector<double>& vec)
  {
    if (firstCall)
    {
      q10 = Q10PsiPhi::provideQ10PsiPhi(rowFeSpace->getBasisFcts(),
                                        colFeSpace->getBasisFcts(),
                                        quadrature);
      q1 = Q1Psi::provideQ1Psi(rowFeSpace->getBasisFcts(), quadrature);
      firstCall = false;
    }

    const int* nEntries = q1->getNumberEntries();
    // Do not need to resize Lb, because it's size is always at least one!
    Lb[0] = 0.0;

    for (OperatorTerm* term : terms)
      static_cast<FirstOrderTerm*>(term)->getLb(elInfo, Lb);

    Lb[0] *= elInfo->getDet();

    for (int i = 0; i < nRow; i++)
    {
      const int* k = q1->getKVec(i);
      const double* values = q1->getValVec(i);
      double val = 0.0;
      for (int m = 0; m < nEntries[i]; m++)
        val += values[m] * Lb[0][k[m]];
      vec[i] += val;
    }
  }

} // end namespace AMDiS
