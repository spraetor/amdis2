#include <vector>

#include <boost/numeric/mtl/mtl.hpp>

#include "Assembler.h"
#include "ZeroOrderAssembler.h"
#include "Operator.h"
#include "QPsiPhi.h"
#include "FiniteElemSpace.h"
#include "Quadrature.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS 
{
  ThreadPrivate<vector<SubAssembler*> > 
  ZeroOrderAssembler::optimizedSubAssemblers;
  
  ThreadPrivate<vector<SubAssembler*> >
  ZeroOrderAssembler::standardSubAssemblers;


  ZeroOrderAssembler::ZeroOrderAssembler(Operator* op,
                              					 Assembler* assembler,
                              					 Quadrature* quad,
                              					 bool optimized)
    : SubAssembler(op, assembler, quad, 0, optimized)
  {}


  ZeroOrderAssembler* ZeroOrderAssembler::getSubAssembler(Operator* op,
                                          							  Assembler* assembler,
                                          							  Quadrature* quad,
                                          							  bool optimized)
  {
    // check if an assembler is needed at all
    if (op->zeroOrder.empty())
      return NULL;   

    vector<SubAssembler*>& subAssemblers =
      optimized ? optimizedSubAssemblers.get() : standardSubAssemblers.get();

    vector<OperatorTerm*> opTerms = op->zeroOrder;
    sort(opTerms.begin(), opTerms.end());

    // check if a new assembler is needed
    if (quad) {
      for (SubAssembler* subAssembler : subAssemblers) {
      	vector<OperatorTerm*> assTerms = *(subAssembler->getTerms());
      	sort(assTerms.begin(), assTerms.end());
      
      	if (opTerms == assTerms && subAssembler->getQuadrature() == quad)
      	  return dynamic_cast<ZeroOrderAssembler*>(subAssembler);
      }
    }
 
    // check if all terms are pw_const
    bool pwConst = std::all_of(op->zeroOrder.begin(), op->zeroOrder.end(), 
                               [](OperatorTerm* term){ return term->isPWConst(); });

    // create new assembler
    ZeroOrderAssembler *newAssembler;
    if (!optimized)
      newAssembler = new StandardZOA(op, assembler, quad);
    else {
      if (pwConst)
        newAssembler = new PrecalcZOA(op, assembler, quad);
      else
     	newAssembler = new FastQuadZOA(op, assembler, quad);      
    }

    subAssemblers.push_back(newAssembler);

    return newAssembler;
  }


  StandardZOA::StandardZOA(Operator *op, Assembler *assembler, Quadrature *quad)
    : ZeroOrderAssembler(op, assembler, quad, false)      
  {
    name = "standard zero order assembler";
  }


  void StandardZOA::calculateElementMatrixImpl(const ElInfo *elInfo, 
				                                       ElementMatrix& mat)
  {
    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    int nPoints = quadrature->getNumPoints();
    DenseVector<double> c(nPoints, 0.0);
    DenseVector<double> phival(nCol);

    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, nPoints, c);

    // TODO: use the case that psi == phi and phival == psival
    if (symmetric) {
      TEST_EXIT_DBG(nCol == nRow)("nCol != nRow, but symmetric assembling!\n");
            
      for (int iq = 0; iq < nPoints; iq++) {
      	c[iq] *= elInfo->getDet();
      
      	// calculate phi at QPs only once!
      	for (int i = 0; i < nCol; i++)
      	  phival[i] = (*(phi->getPhi(i)))(quadrature->getLambda(iq));
      
      	for (int i = 0; i < nRow; i++) {
      	  double psival = (*(psi->getPhi(i)))(quadrature->getLambda(iq));
      	  mat[i][i] += quadrature->getWeight(iq) * c[iq] * psival * phival[i];
      	  for (int j = i + 1; j < nCol; j++) {
      	    double val = quadrature->getWeight(iq) * c[iq] * psival * phival[j];
      	    mat[i][j] += val;
      	    mat[j][i] += val;
      	  }
      	}
      }
    } else {      //  non symmetric assembling 
      for (int iq = 0; iq < nPoints; iq++) {
      	c[iq] *= elInfo->getDet();
      
      	// calculate phi at QPs only once!
      	for (int i = 0; i < nCol; i++)
      	  phival[i] = (*(phi->getPhi(i)))(quadrature->getLambda(iq));
      
      	for (int i = 0; i < nRow; i++) {
      	  double psival = (*(psi->getPhi(i)))(quadrature->getLambda(iq));
      	  for (int j = 0; j < nCol; j++)
      	    mat[i][j] += quadrature->getWeight(iq) * c[iq] * psival * phival[j];	  
      	}
      }
    }
  }

  // TODO: add combination of elementMatrix and elementVector assembling
  void StandardZOA::calculateElementVectorImpl(const ElInfo *elInfo, DenseVector<double>& vec)
  {
    int nPoints = quadrature->getNumPoints();
    DenseVector<double> c(nPoints, 0.0);

    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, nPoints, c);

    for (int iq = 0; iq < nPoints; iq++) {
      c[iq] *= elInfo->getDet();

      for (int i = 0; i < nRow; i++) {
      	double psi = (*(rowFeSpace->getBasisFcts()->getPhi(i)))
      	  (quadrature->getLambda(iq));
      	vec[i] += quadrature->getWeight(iq) * c[iq] * psi;
      }
    }
  }


  FastQuadZOA::FastQuadZOA(Operator *op, Assembler *assembler, Quadrature *quad)
    : ZeroOrderAssembler(op, assembler, quad, true)
  {
    name = "fast quadrature zero order assembler";
  }


  void FastQuadZOA::calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat)
  {
    int nPoints = quadrature->getNumPoints();

    if (firstCall) {
      const BasisFunction *basFcts = rowFeSpace->getBasisFcts();
      psiFast = updateFastQuadrature(psiFast, basFcts, INIT_PHI);
      basFcts = colFeSpace->getBasisFcts();
      phiFast = updateFastQuadrature(phiFast, basFcts, INIT_PHI);
      firstCall = false;      
    }

    if (num_rows(c) < static_cast<size_t>(nPoints))
      c.change_dim(nPoints);
    set_to_zero(c);

    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, nPoints, c);

    const DenseMatrix<double>& psi = psiFast->getPhi();
    const DenseMatrix<double>& phi = phiFast->getPhi();

    if (symmetric) {
      TEST_EXIT_DBG(nCol == nRow)("nCol != nRow, but symmetric assembling!\n");

      for (int iq = 0; iq < nPoints; iq++) {
      	c[iq] *= elInfo->getDet() * quadrature->getWeight(iq);
      
      	for (int i = 0; i < nRow; i++) {
      	  mat[i][i] += c[iq] * psi[iq][i] * phi[iq][i];
      	  for (int j = i + 1; j < nCol; j++) {
      	    double val = c[iq] * psi[iq][i] * phi[iq][j];
      	    mat[i][j] += val;
      	    mat[j][i] += val;
      	  }
      	}
      }
    } else {      /*  non symmetric assembling   */
      for (int iq = 0; iq < nPoints; iq++)
        c[iq] *= elInfo->getDet() * quadrature->getWeight(iq);

      for (int i = 0; i < nRow; i++)
      	for (int j = 0; j < nCol; j++) {
      	  double v = 0.0;
      	  for (int iq = 0; iq < nPoints; iq++)
      	    v += c[iq] * psi[iq][i] * phi[iq][j];
      	  mat[i][j] += v;
      	}
    }
  }


  void FastQuadZOA::calculateElementVectorImpl(const ElInfo *elInfo, DenseVector<double>& vec)
  {
    int nPoints = quadrature->getNumPoints();

    if (firstCall) {
      const BasisFunction *basFcts = rowFeSpace->getBasisFcts();
      psiFast = updateFastQuadrature(psiFast, basFcts, INIT_PHI);
      basFcts = colFeSpace->getBasisFcts();
      phiFast = updateFastQuadrature(phiFast, basFcts, INIT_PHI);
      firstCall = false;      
    }

    DenseVector<double> c(nPoints, 0.0);
    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, nPoints, c);

    const DenseMatrix<double>& psi = psiFast->getPhi();

    for (int iq = 0; iq < nPoints; iq++) {
      c[iq] *= elInfo->getDet();

      for (int i = 0; i < nRow; i++)
        vec[i] += quadrature->getWeight(iq) * c[iq] * psi[iq][i];
    }
  }


  PrecalcZOA::PrecalcZOA(Operator *op, Assembler *assembler, Quadrature *quad) 
    : ZeroOrderAssembler(op, assembler, quad, true)
  {
    name = "precalculated zero order assembler";
  }


  void PrecalcZOA::calculateElementMatrixImpl(const ElInfo *elInfo, 
					                                    ElementMatrix& mat)
  {
    if (firstCall) {
      q00 = Q00PsiPhi::provideQ00PsiPhi(rowFeSpace->getBasisFcts(), 
					colFeSpace->getBasisFcts(), 
					quadrature);
      q0 = Q0Psi::provideQ0Psi(rowFeSpace->getBasisFcts(), quadrature);
      firstCall = false;      
    }

    DenseVector<double> c(1, 0.0);
    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, 1, c);    

    c[0] *= elInfo->getDet();
    
    if (symmetric) {
      for (int i = 0; i < nRow; i++) {
      	mat[i][i] += c[0] * q00->getValue(i,i);
      	for (int j = i + 1; j < nCol; j++) {
      	  double val = c[0] * q00->getValue(i, j);
      	  mat[i][j] += val;
      	  mat[j][i] += val;
      	}
      }
    } else {
      for (int i = 0; i < nRow; i++)
      	for (int j = 0; j < nCol; j++)
      	  mat[i][j] += c[0] * q00->getValue(i, j);
    }
  }


  void PrecalcZOA::calculateElementVectorImpl(const ElInfo *elInfo, 
					                                    DenseVector<double>& vec)
  {
    if (firstCall) {
      q00 = Q00PsiPhi::provideQ00PsiPhi(rowFeSpace->getBasisFcts(), 
                                        colFeSpace->getBasisFcts(), 
                                        quadrature);
      q0 = Q0Psi::provideQ0Psi(rowFeSpace->getBasisFcts(), quadrature);	
      firstCall = false;      
    }

    DenseVector<double> c(1, 0.0);    
    for (auto& term : terms)
      (static_cast<ZeroOrderTerm*>(term))->getC(elInfo, 1, c);

    c[0] *= elInfo->getDet();

    for (int i = 0; i < nRow; i++)
      vec[i] += c[0] * q0->getValue(i);
  }

} // end namespace AMDiS
