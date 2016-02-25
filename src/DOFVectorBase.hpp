/** \file DOFVector.h */

#pragma once

// std c++ headers
#include <vector>
#include <map>

// AMDiS headers
#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"
#include "Boundary.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "DOFIndexed.hpp"

namespace AMDiS
{
  template <class T>
  class DOFVectorBase : public DOFIndexed<T>
  {
  public:
    DOFVectorBase()
      : feSpace(NULL),
        elementVector(3),
        boundaryManager(NULL),
        nBasFcts(0)
    {}

    DOFVectorBase(const FiniteElemSpace* f, std::string n);

    virtual ~DOFVectorBase();

    /// For the given element, this function returns an array of all DOFs of
    /// this DOFVector that are defined on this element.
    virtual void getLocalVector(const Element* el,
                                DenseVector<T>& localVec) const;

    /// Evaluates the DOF vector at a set of quadrature points defined on the
    /// given element.
    void getVecAtQPs( const ElInfo* elInfo,
                      const Quadrature* quad,
                      const FastQuadrature* quadFast,
                      DenseVector<T>& vecAtQPs) const;

    /// Evaluates the gradient of a DOF vector at a set of quadrature points
    /// defined on the given element.
    void getGrdAtQPs( const ElInfo* elInfo,
                      const Quadrature* quad,
                      const FastQuadrature* quadFast,
                      DenseVector<Gradient_t<T>>& grdAtQPs) const;

    /// Evaluates the comp'th component of the derivative of a DOF vector at a
    /// set of quadrature points defined on the given element.
    void getDerivativeAtQPs(const ElInfo* elInfo,
                            const Quadrature* quad,
                            const FastQuadrature* quadFast,
                            int comp,
                            DenseVector<T>& derivativeAtQPs) const;

    /// Evaluates the jacobian of a DOF vector at a set of quadrature points
    ///  defined on the given element.
    void getD2AtQPs(const ElInfo* elInfo,
                    const Quadrature* quad,
                    const FastQuadrature* quadFast,
                    DenseVector<typename D2Type<T>::type>& d2AtQPs) const;

    /// Returns the FE space the DOF vector is defined on.
    const FiniteElemSpace* getFeSpace() const
    {
      return feSpace;
    }

    /// Assembles the element vector for the given ellement and adds the
    /// element matrix to the current DOF vector.
    void assemble(T factor, ElInfo* elInfo,
                  const BoundaryType* bound,
                  Operator* op = NULL);

    void addElementVector(T sign,
                          const DenseVector<double>& elVec,
                          const BoundaryType* bound,
                          ElInfo* elInfo,
                          bool add = true);

    ///
    void assembleOperator(Operator& op);

    /// That function must be called after the matrix assembling has been
    /// finished. This makes it possible to start some cleanup or matrix
    /// data compressing procedures.
    void finishAssembling();

    void addOperator(Operator* op,
                     double* factor = NULL,
                     double* estFactor = NULL)
    {
      operators.push_back(op);
      operatorFactor.push_back(factor);
      operatorEstFactor.push_back(estFactor);
    }

    std::vector<double*>::iterator getOperatorFactorBegin()
    {
      return operatorFactor.begin();
    }

    std::vector<double*>::iterator getOperatorFactorEnd()
    {
      return operatorFactor.end();
    }

    std::vector<double*>::iterator getOperatorEstFactorBegin()
    {
      return operatorEstFactor.begin();
    }

    std::vector<double*>::iterator getOperatorEstFactorEnd()
    {
      return operatorEstFactor.end();
    }

    std::vector<Operator*>::iterator getOperatorsBegin()
    {
      return operators.begin();
    }

    std::vector<Operator*>::iterator getOperatorsEnd()
    {
      return operators.end();
    }

    Flag getAssembleFlag();

    /// Evaluates \f[ u_h(x(\lambda)) = \sum_{i=0}^{m-1} vec[ind[i]] *
    /// \varphi^i(\lambda) \f] where \f$ \varphi^i \f$ is the i-th basis
    /// function, \f$ x(\lambda) \f$ are the world coordinates of lambda
    /// and \f$ m \f$ is the number of basis functions
    T evalUh(const DimVec<double>& lambda, DegreeOfFreedom* ind);

    std::vector<Operator*>& getOperators()
    {
      return operators;
    }

    std::vector<double*>& getOperatorFactor()
    {
      return operatorFactor;
    }

    std::vector<double*>& getOperatorEstFactor()
    {
      return operatorEstFactor;
    }

    /// Returns \ref name
    std::string getName() const
    {
      return name;
    }

    void setName(std::string n)
    {
      name = n;
    }

    BoundaryManager* getBoundaryManager() const
    {
      return boundaryManager;
    }

    void setBoundaryManager(BoundaryManager* bm)
    {
      boundaryManager = bm;
    }

    void setDirichletDofValue(DegreeOfFreedom dof, T value)
    {
      dirichletDofValues[dof] = value;
    }

    std::map<DegreeOfFreedom, T>& getDirichletValues()
    {
      return dirichletDofValues;
    }

  protected:
    ///
    const FiniteElemSpace* feSpace;

    ///
    std::string name;

    ///
    DenseVector<double> elementVector;

    ///
    std::vector<Operator*> operators;

    ///
    std::vector<double*> operatorFactor;

    ///
    std::vector<double*> operatorEstFactor;

    ///
    BoundaryManager* boundaryManager;

    /// Number of basis functions of the used finite element space.
    int nBasFcts;

    /// Dimension of the mesh this DOFVectorBase belongs to
    int dim;

    std::map<DegreeOfFreedom, T> dirichletDofValues;
  };


} // end namespace AMDiS
