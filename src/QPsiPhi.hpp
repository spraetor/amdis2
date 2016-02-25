#pragma once

#include <list>

#include "AMDiS_fwd.hpp"
#include "DOFAdmin.hpp"

namespace AMDiS
{
  /** \ingroup Integration
   * \brief
   * Calculates element stiffness matrices by preevaluated integrals over the
   * the reference element (integral of the product of the derivatives of psi
   * and phi).
   */
  class Q11PsiPhi
  {
  protected:
    /// Constructor
    Q11PsiPhi(BasisFunction const* psi,
              BasisFunction const* phi,
              Quadrature const* q);

  public:
    /// Destructor
    ~Q11PsiPhi();

    /// Returns a Q11PsiPhi object.
    static Q11PsiPhi const* provideQ11PsiPhi(BasisFunction const*,
					     BasisFunction const*,
					     Quadrature const*);

    /// Compares two Q11PsiPhi objects.
    bool operator==(Q11PsiPhi const&) const;

    /// Compares two Q11PsiPhi objects.
    bool operator!=(Q11PsiPhi const& q11pp) const
    {
      return !(operator==(q11pp));
    }

    /// Returns \ref values[i][j][k]
    double getValue(unsigned int i,
                    unsigned int j,
                    unsigned int v) const
    {
      if (values && values[i] && values[i][j] &&
          (static_cast<int>(v) < nrEntries[i][j]))
        return values[i][j][v];

      return 0.0;
    }

    /// Returns \ref nrEntries[i][j]
    int getNumberEntries(unsigned int i, unsigned int j) const
    {
      if (nrEntries && nrEntries[i])
        return nrEntries[i][j];

      return 0;
    }

    /// Returns \ref nrEntries
    const int** getNumberEntries() const
    {
      return const_cast<const int**>(nrEntries);
    }

    /// Returns \ref k[i1][i2][i3]
    int getK(unsigned int i1,
             unsigned int i2,
             unsigned int i3) const
    {
      if (k && k[i1] && k[i1][i2] && (static_cast<int>(i3) < nrEntries[i1][i2]))
        return k[i1][i2][i3];

      return 0;
    }

    /// Returns \ref l[i][j][v]
    int getL(unsigned int i, unsigned int j, unsigned int v) const
    {
      if (l && l[i] && l[i][j] && (static_cast<int>(v) < nrEntries[i][j]))
        return l[i][j][v];

      return 0;
    }

    /// Returns \values[i][j]
    const double* getValVec(unsigned int i, unsigned int j) const
    {
      if (values && values[i] && values[i][j])
        return values[i][j];

      return NULL;
    }

    /// Returns \ref k[i][j]
    const int* getKVec(unsigned int i, unsigned int j) const
    {
      if (k && k[i] && k[i][j])
        return k[i][j];

      return NULL;
    }

    /// Returns \ref l[i][j]
    const int* getLVec(unsigned int i, unsigned int j) const
    {
      if (l && l[i] && l[i][j])
        return l[i][j];

      return NULL;
    }

  protected:
    /// List of pointers to all Q11PsiPhi objects.
    static std::list<Q11PsiPhi*> preList;

    /// Pointer to the first set of basis functions.
    BasisFunction const* psi;

    /// Pointer to the second set of basis functions.
    BasisFunction const* phi;

    /// Pointer to the Quadrature which is used for the integration.
    Quadrature const* quadrature;

    /** \brief
     * Matrix of size psi->getNumber() * phi->getNumber() storing the count of
     * non zero integrals; nrEntries[i][j] is the count of non zero values of
     * \f$ \hat{Q}_{ij,kl}^{11} \f$
     * (0 <= k,l <= DIM) for the pair (psi[i], phi[j]), 0 <= i <
     * psi->getNumber(), 0 <= j < phi->getNumber()
     */
    int** nrEntries;

    /** \brief
     * tensor storing the non zero integrals; values[i][j] is a vector of length
     * \ref nrEntries[i][j] storing the non zero values for the pair (psi[i],
     * phi[j])
     */
    double*** values;

    /// Tensor storing the indices k of the non zero integrals.
    int*** k;

    /// Tensor storing the indices l of the non zero integrals.
    int*** l;

    /// Number of non zero entries.
    int allEntries;

    /// Pointer to an array in values.
    double* val_vec;

    /// Pointer to an array in k.
    int* k_vec;

    /// Pointer to an array in l.
    int* l_vec;
  };


  /** \ingroup Integration
   * \brief
   * Calculates element stiffness matrices by preevaluated integrals over the
   * the reference element (integral of the product of the derivative of psi
   * and phi).
   */
  class Q10PsiPhi
  {
  protected:
    /// Constructor
    Q10PsiPhi(BasisFunction const* psi,
              BasisFunction const* phi,
              Quadrature const* q);

  public:
    /// Destructor
    ~Q10PsiPhi();

    /// Returns a Q10PsiPhi object.
    static Q10PsiPhi const* provideQ10PsiPhi(BasisFunction const*,
					     BasisFunction const*,
					     Quadrature const*);

    /// Compares two Q10PsiPhi objects.
    bool operator==(Q10PsiPhi const&) const;

    /// Compares two Q10PsiPhi objects.
    bool operator!=(Q10PsiPhi const& q10pp) const
    {
      return !(operator==(q10pp));
    }

    /// Returns \ref values[i][j][k]
    double getValue(unsigned int i,
                    unsigned int j,
                    unsigned int v) const
    {
      if (values && values[i] && values[i][j] &&
          (static_cast<int>(v) < nrEntries[i][j]))
        return values[i][j][v];

      return 0.0;
    }

    /// Returns \ref nrEntries[i][j]
    int getNumberEntries(unsigned int i, unsigned int j) const
    {
      if (nrEntries && nrEntries[i])
        return nrEntries[i][j];

      return 0;
    }

    /// Returns \ref nrEntries
    const int** getNumberEntries() const
    {
      return const_cast<const int**>(nrEntries);
    }

    /// Returns \ref k[i1][i2][i3]
    int getK(unsigned int i1,
             unsigned int i2,
             unsigned int i3) const
    {
      if (k && k[i1] && k[i1][i2] && (static_cast<int>(i3) < nrEntries[i1][i2]))
        return k[i1][i2][i3];

      return 0;
    }

    /// Returns \values[i][j]
    const double* getValVec(unsigned int i, unsigned int j) const
    {
      if (values && values[i] && values[i][j])
        return values[i][j];

      return NULL;
    }

    /// Returns \ref k[i][j]
    const int* getKVec(unsigned int i, unsigned int j) const
    {
      if (k && k[i] && k[i][j])
        return k[i][j];

      return NULL;
    }

  protected:
    /// List of pointers to all Q11PsiPhi objects
    static std::list<Q10PsiPhi*> preList;

    /// Pointer to the first set of basis functions
    BasisFunction const* psi;

    /// Pointer to the second set of basis functions
    BasisFunction const* phi;

    /// Pointer to the Quadrature which is used for the integration
    Quadrature const* quadrature;

    /** \brief
     * Matrix of size psi->getNumber() * phi->getNumber() storing the count of
     * non zero integrals; nrEntries[i][j] is the count of non zero values of
     * \f$ \hat{Q}_{ij,kl}^{11} \f$
     * (0 <= k,l <= DIM) for the pair (psi[i], phi[j]), 0 <= i <
     * psi->getNumber(), 0 <= j < phi->getNumber()
     */
    int** nrEntries;

    /** \brief
     * tensor storing the non zero integrals; values[i][j] is a vector of length
     * \ref nrEntries[i][j] storing the non zero values for the pair (psi[i],
     * phi[j])
     */
    double*** values;

    /// Tensor storing the indices k of the non zero integrals;
    int*** k;

    /// Number of all non zero entries.
    int allEntries;

    /// Pointer to an array in values.
    double* val_vec;

    /// Pointer to an array in k.
    int* k_vec;
  };


  /** \ingroup Integration
   * \brief
   * Calculates element stiffness matrices by preevaluated integrals over the
   * the reference element (integral of the product of psi and the derivative
   * of phi).
   */
  class Q01PsiPhi
  {
  protected:
    /// Constructor
    Q01PsiPhi(BasisFunction const* psi,
              BasisFunction const* phi,
              Quadrature const* q);

  public:
    /// Destructor
    ~Q01PsiPhi();

    /// Returns a Q01PsiPhi object.
    static Q01PsiPhi const* provideQ01PsiPhi(BasisFunction const*,
					     BasisFunction const*,
					     Quadrature const*);

    /// Compares two Q01PsiPhi objects.
    bool operator==(Q01PsiPhi const&) const;

    /// Compares two Q01PsiPhi objects.
    bool operator!=(Q01PsiPhi const& q01pp) const
    {
      return !(operator==(q01pp));
    }

    /// Returns \ref values[i][j][k]
    double getValue(unsigned int i,
                    unsigned int j,
                    unsigned int v) const
    {
      if (values && values[i] && values[i][j] &&
          (static_cast<int>(v) < nrEntries[i][j]))
        return values[i][j][v];

      return 0.0;
    }

    /// Returns \ref nrEntries[i][j]
    int getNumberEntries(unsigned int i, unsigned int j) const
    {
      if (nrEntries && nrEntries[i])
        return nrEntries[i][j];

      return 0;
    }

    /// Returns \ref nrEntries
    const int** getNumberEntries() const
    {
      return const_cast<const int**>(nrEntries);
    }

    /// Returns \ref k[i1][i2][i3]
    inline int getK(unsigned int i1,
                    unsigned int i2,
                    unsigned int i3) const;

    /// Returns \values[i][j]
    const double* getValVec(unsigned int i, unsigned int j) const
    {
      if (values && values[i] && values[i][j])
        return values[i][j];

      return NULL;
    }

    /// Returns \ref k[i][j]
    const int* getLVec(unsigned int i, unsigned int j) const
    {
      if (l && l[i] && l[i][j])
        return l[i][j];

      return NULL;
    }

    /// Returns \ref k[i][j][v]
    int getL(unsigned int i,
             unsigned int j,
             unsigned int v) const
    {
      if (l && l[i] && l[i][j] && (static_cast<int>(v) < nrEntries[i][j]))
        return l[i][j][v];

      return 0;
    }

  protected:
    /// List of pointers to all Q11PsiPhi objects
    static std::list<Q01PsiPhi*> preList;

    /// Pointer to the first set of basis functions
    BasisFunction const* psi;

    /// Pointer to the second set of basis functions
    BasisFunction const* phi;

    /// Pointer to the Quadrature which is used for the integration
    Quadrature const* quadrature;

    /** \brief
     * Matrix of size psi->getNumber() * phi->getNumber() storing the count of
     * non zero integrals; nrEntries[i][j] is the count of non zero values of
     * \f$ \hat{Q}_{ij,kl}^{11} \f$
     * (0 <= k,l <= DIM) for the pair (psi[i], phi[j]), 0 <= i <
     * psi->getNumber(), 0 <= j < phi->getNumber()
     */
    int** nrEntries;

    /** \brief
     * tensor storing the non zero integrals; values[i][j] is a vector of length
     * \ref nrEntries[i][j] storing the non zero values for the pair (psi[i],
     * phi[j])
     */
    double*** values;

    /// Tensor storing the indices l of the non zero integrals;
    int*** l;

    /// Number of all non zero entries.
    int allEntries;

    /// Pointer to an array in values
    double* val_vec;

    /// Pointer to an array in l.
    int* l_vec;
  };


  /** \ingroup Integration
   * \brief
   * Calculates element stiffness matrices by preevaluated integrals over the
   * the reference element (integral of the product of psi and phi).
   */
  class Q00PsiPhi
  {
    /// List of pointers to all Q11PsiPhi objects
    static std::list<Q00PsiPhi*> preList;

    /// Pointer to the first set of basis functions
    BasisFunction const* psi;

    /// Pointer to the second set of basis functions
    BasisFunction const* phi;

    /// Pointer to the Quadrature which is used for the integration
    Quadrature const* quadrature;

    /** \brief
     * Matrix storing the integrals
     * \f[ values[i][j] = \hat{Q}_{ij}^{00} = \hat{Q}(\overline{\psi}^i
     * \overline{\phi}^j) \f]
     * for the pair (psi[i], phi[j]), 0 <= i <= psi->getNumber(),
     * 0 <= j <= phi->getNumber()
     */
    double** values;

  protected:
    /// Constructor
    Q00PsiPhi(BasisFunction const*, BasisFunction const*, Quadrature const*);

  public:
    /// Destructor
    ~Q00PsiPhi();

    /// Returns a Q00PsiPhi object.
    static Q00PsiPhi const* provideQ00PsiPhi(BasisFunction const*,
					     BasisFunction const*,
					     Quadrature const*);

    /// Compares two Q00PsiPhi objects.
    bool operator==(Q00PsiPhi const& q00pp) const
    {
      return (q00pp.psi == psi && q00pp.phi == phi && q00pp.quadrature == quadrature);
    }

    /// Compares two Q00PsiPhi objects.
    bool operator!=(Q00PsiPhi const& q00pp) const
    {
      return !(operator==(q00pp));
    }

    /// Returns \ref values[i][j]
    double getValue(unsigned int i, unsigned int j) const;

    /// Returns \ref values[i]
    const double* getValVec(unsigned int i) const;
  };



  inline double Q00PsiPhi::getValue(unsigned int i,unsigned  int j) const
  {
    if ((values)&&(values[i])) return values[i][j];
    return 0.;
  }

  inline const double* Q00PsiPhi::getValVec(unsigned int i) const
  {
    if ((values)&&(values[i])) return values[i];
    return NULL;
  }


  /** \ingroup Integration
   * \brief
   * Integral of psi.
   */
  class Q0Psi
  {
  protected:
    /// List of pointers to all Q0Psi objects
    static std::list<Q0Psi*> preList;

    /// Pointer to the first set of basis functions
    BasisFunction const* psi;

    /// Pointer to the Quadrature which is used for the integration
    Quadrature const* quadrature;

    /// Vector storing the integrals
    double* values;

  protected:
    /// Constructor
    Q0Psi(BasisFunction const*, Quadrature const*);

  public:
    /// Destructor
    ~Q0Psi();

    /// Returns a Q0Psi object.
    static Q0Psi* provideQ0Psi(BasisFunction const*, Quadrature const*);

    /// Compares two Q0Psi objects.
    bool operator==(Q0Psi const& q0p) const
    {
      return (q0p.psi == psi && q0p.quadrature == quadrature);
    }

    /// Compares two Q0Psi objects.
    bool operator!=(Q0Psi const& q0p) const
    {
      return !(operator==(q0p));
    }

    /// Returns \ref value
    double getValue(int i) const
    {
      return values[i];
    }

    /// Returns \ref values[i]
    const double* getValVec() const
    {
      return values;
    }
  };


  /** \ingroup Integration
   * \brief
   * Integral of the derivative of psi.
   */
  class Q1Psi
  {
  protected:
    /// Constructor
    Q1Psi(BasisFunction const* psi, Quadrature const* q);


  public:
    /// Destructor
    ~Q1Psi();

    /// Returns a Q1Psi object.
    static Q1Psi const* provideQ1Psi(BasisFunction const*,
                                     Quadrature const*);

    /// Compares two Q1Psi objects.
    bool operator==(Q1Psi const&) const;

    /// Compares two Q1Psi objects.
    bool operator!=(Q1Psi const& q1p) const
    {
      return !(operator==(q1p));
    }

    /// Returns \ref values[i][j]
    double getValue(unsigned int i,
                    unsigned int j) const
    {
      if (values && values[i] && (static_cast<int>(j) < nrEntries[i]))
        return values[i][j];

      return 0.0;
    }

    /// Returns \ref nrEntries[i]
    int getNumberEntries(unsigned int i) const
    {
      return (nrEntries ? nrEntries[i] : 0);
    }

    /// Returns \ref nrEntries
    int const* getNumberEntries() const
    {
      return const_cast<const int*>(nrEntries);
    }

    /// Returns \ref k[i1][i2]
    int getK(unsigned int i1, unsigned int i2) const
    {
      if (k && k[i1] && (static_cast<int>(i2) < nrEntries[i1]))
        return k[i1][i2];

      return 0;
    }

    /// Returns \ref k[i]
    int const* getKVec(unsigned int i) const
    {
      if (k && k[i])
        return k[i];

      return NULL;
    }

    /// Returns \values[i]
    double const* getValVec(unsigned int i) const
    {
      if (values && values[i])
        return values[i];

      return NULL;
    }

  protected:
    /// List of pointers to all Q1Psi objects
    static std::list<Q1Psi*> preList;

    /// Pointer to the first set of basis functions
    BasisFunction const* psi;

    /// Pointer to the Quadrature which is used for the integration
    Quadrature const* quadrature;

    /** \brief
     * Array of size psi->getNumber() storing the count of
     * non zero integrals; nrEntries[i] is the count of non zero values of
     * \f$ \hat{Q}_{i,kl}^{11} \f$
     * (0 <= k <= DIM) for psi[i], 0 <= i < psi->getNumber().
     */
    int* nrEntries;

    /// Number of all non zero entries.
    int allEntries;

    /// Tensor storing the non zero integrals; values[i] is a vector of length
    /// \ref nrEntries[i] storing the non zero values for psi[i].
    double** values;

    /// Matrix storing the indices k of the non zero integrals;
    int** k;

    /// Pointer to an array in values.
    double* val_vec;

    /// Pointer to an array in k.
    int* k_vec;
  };

} // end namespace AMDiS
