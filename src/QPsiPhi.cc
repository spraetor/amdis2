#include <algorithm>
#include <functional>

#include "QPsiPhi.h"
#include "BasisFunction.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Quadrature.h"
#include "Mesh.h"
#include "FixVec.h"
#include "DOFMatrix.h"
#include "DOFVector.h"
#include "Global.h"

namespace AMDiS
{
  constexpr double TOO_SMALL = 1.e-15;

  std::list<Q11PsiPhi*> Q11PsiPhi::preList;
  std::list<Q01PsiPhi*> Q01PsiPhi::preList;
  std::list<Q00PsiPhi*> Q00PsiPhi::preList;
  std::list<Q10PsiPhi*> Q10PsiPhi::preList;

  std::list<Q0Psi*> Q0Psi::preList;
  std::list<Q1Psi*> Q1Psi::preList;

  /****************************************************************************/
  /* information about precomputed integrals of basis functions on the        */
  /* standard element;                                                        */
  /****************************************************************************/

  /****************************************************************************/
  /* second order term:                                                       */
  /****************************************************************************/

  Q11PsiPhi::Q11PsiPhi(const BasisFunction* ps, const BasisFunction* ph,
                       const Quadrature* quadrat)
    : psi(ps),
      phi(ph),
      quadrature(quadrat),
      nrEntries(NULL),
      values(NULL),
      k(NULL),
      l(NULL)
  {
    FUNCNAME_DBG("Q11PsiPhi::Q11PsiPhi()");

    const FastQuadrature* q_phi, *q_psi;
    int j, lk, ll, n, iq, all_entries, n_psi, n_phi;
    double val, grdi, grdj;
    int d = ps->getDim();

    if (!psi)
      psi = phi;
    if (!phi)
      phi = psi;

    if (!quadrature)
      quadrature = Quadrature::provideQuadrature(d, psi->getDegree() +
                   phi->getDegree() - 2);

    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_GRD_PHI);
    n_phi = phi->getNumber();
    q_phi = FastQuadrature::provideFastQuadrature(phi, *quadrature, INIT_GRD_PHI);

    nrEntries = new int* [n_psi];
    values = new double** [n_psi];
    k = new int** [n_psi];
    l = new int** [n_psi];
    for (int i = 0; i <n_psi; i++)
    {
      nrEntries[i] = new int[n_phi];
      values[i] = new double*[n_phi];
      k[i] = new int* [n_phi];
      l[i] = new int* [n_phi];
    }


    //****************************************************************************
    //*  compute first the number of all non zero entries                        *
    //****************************************************************************

    int numPoints = quadrature->getNumPoints();

    all_entries = 0;
    for (int i = 0; i < n_psi; i++)
    {
      for (j = 0; j < n_phi; j++)
      {
        for (lk = 0; lk < d+1; lk++)
        {
          for (ll =  0; ll < d+1; ll++)
          {
            for (val = iq = 0; iq < numPoints; iq++)
            {
              grdi = q_psi->getGradient(iq,i,lk);
              grdj = q_phi->getGradient(iq,j,ll);

              val += quadrature->getWeight(iq)*grdi*grdj;
            }
            if (math::abs(val) > TOO_SMALL)
              all_entries++;
          }
        }
      }
    }

    //****************************************************************************
    //* now, access memory for all information                                   *
    //****************************************************************************
    allEntries = all_entries;

    val_vec = new double[all_entries];
    k_vec = new int[all_entries];
    l_vec = new int[all_entries];

    //***************************************************************************
    // and now, fill information                                                *
    //***************************************************************************

    for (int i = 0; i < n_psi; i++)
    {
      for (j = 0; j < n_phi; j++)
      {
        values[i][j] = val_vec;
        k[i][j]     = k_vec;
        l[i][j]     = l_vec;

        for (n = lk = 0; lk < d+1; lk++)
        {
          for (ll =  0; ll < d+1; ll++)
          {
            for (val = iq = 0; iq < numPoints; iq++)
            {
              grdi = q_psi->getGradient(iq,i,lk);
              grdj = q_phi->getGradient(iq,j,ll);

              val += quadrature->getWeight(iq)*grdi*grdj;
            }
            if (math::abs(val) > TOO_SMALL)
            {
              TEST_EXIT_DBG(all_entries > 0)("now more entries found than counted before\n");
              all_entries--;
              n++;
              *val_vec = val;
              val_vec++;
              *k_vec = lk;
              k_vec++;
              *l_vec = ll;
              l_vec++;
            }
          }
        }
        nrEntries[i][j] = n;
      }
    }
  }

  Q11PsiPhi::~Q11PsiPhi()
  {
    if (nrEntries)
    {
      for (int i = 0; i < psi->getNumber(); i++)
      {
        delete [] nrEntries[i];
        delete [] values[i];
        delete [] k[i];
        delete [] l[i];
      }
      delete [] nrEntries;
      delete [] values;
      delete [] k;
      delete [] l;
      delete [] val_vec;
      delete [] k_vec;
      delete [] l_vec;
    }
  }


  bool Q11PsiPhi::operator==(const Q11PsiPhi& q11pp) const
  {
    return (q11pp.psi == psi && q11pp.phi == phi && q11pp.quadrature == quadrature);
  }


  const Q11PsiPhi* Q11PsiPhi::provideQ11PsiPhi(const BasisFunction* ps,
      const BasisFunction* ph,
      const Quadrature*    quadrat)
  {
    std::list<Q11PsiPhi*>::iterator list;

    if (!ps && !ph) return NULL;

    if (!ps)  ps = ph;
    if (!ph)  ph = ps;

    if (!quadrat)
      quadrat = Quadrature::provideQuadrature(ps->getDim(),
                                              ps->getDegree() +
                                              ph->getDegree() - 2);

    compareQPsiPhi<Q11PsiPhi> comp(ps,ph,quadrat);

    //***************************************************************************
    //  look for an existing entry in the list                                  *
    //***************************************************************************

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q11PsiPhi* newQ11PsiPhi = new Q11PsiPhi(ps,ph,quadrat);
      preList.push_back(newQ11PsiPhi);
      return newQ11PsiPhi;
    }
    else
    {
      return *list;
    }
  }

  /****************************************************************************/
  /* first order term                                                         */
  /****************************************************************************/

  Q10PsiPhi::Q10PsiPhi(const BasisFunction* ps, const BasisFunction* ph,
                       const Quadrature* quadrat):psi(ps),phi(ph),quadrature(quadrat)
  {
    FUNCNAME_DBG("Q10PsiPhi::Q10PsiPhi");
    const FastQuadrature*      q_phi, *q_psi;
    int                       i, j, lk, n=0, iq, all_entries, n_psi, n_phi;
    double                      val, psij, grdi;
    int                       d=quadrature->getDim();

    int numPoints = quadrature->getNumPoints();

    if (!psi && !phi)
    {
      nrEntries=NULL;
      k=NULL;
      values=NULL;
    }

    if (!psi)  psi = phi;
    if (!phi)  phi = psi;

    if (!quadrature)
      quadrature = Quadrature::provideQuadrature(psi->getDim(),
                   psi->getDegree()+phi->getDegree()-1);

    /****************************************************************************/
    /*  create a new one                                                        */
    /****************************************************************************/
    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_PHI);
    n_phi = phi->getNumber();
    q_phi = FastQuadrature::provideFastQuadrature(phi, *quadrature, INIT_GRD_PHI);

    nrEntries = new int* [n_psi];
    values = new double** [n_psi];
    k = new int** [n_psi];
    for (int i = 0; i < n_psi; i++)
    {
      nrEntries[i] = new int[n_phi];
      values[i] = new double*[n_phi];
      k[i] = new int* [n_phi];
    }

    /****************************************************************************/
    /*  compute first the number of all non zero entries                        */
    /****************************************************************************/

    all_entries = 0;
    for(i = 0; i < n_psi; i++)
    {
      for(j = 0; j < n_phi; j++)
      {
        for(lk =  0; lk < d+1; lk++)
        {
          for(val = iq = 0; iq < numPoints; iq++)
          {
            psij = q_psi->getPhi(iq,j);
            grdi = q_phi->getGradient(iq,i,lk);

            val += quadrature->getWeight(iq)*psij*grdi;
          }
          if (math::abs(val) > TOO_SMALL)
            all_entries++;
        }
      }
    }

    /****************************************************************************/
    /* now, access memory for all information                                   */
    /****************************************************************************/

    val_vec = new double[all_entries];
    k_vec = new int[all_entries];

    /****************************************************************************/
    /* and now, fill information                                                */
    /****************************************************************************/

    for(i = 0; i < n_psi; i++)
    {
      for(j = 0; j < n_phi; j++)
      {
        values[i][j] = val_vec;
        k[i][j]     = k_vec;

        for(n = lk =  0; lk < d+1; lk++)
        {
          for(val = iq = 0; iq < numPoints; iq++)
          {
            psij = q_psi->getPhi(iq,j);
            grdi = q_phi->getGradient(iq,i,lk);

            val += quadrature->getWeight(iq)*psij*grdi;
          }
          if (math::abs(val) > TOO_SMALL)
          {
            TEST_EXIT_DBG(all_entries > 0)("now more entries found than counted before\n");
            all_entries--;
            n++;
            *val_vec++ = val;
            *k_vec++ = lk;
          }
        }
        nrEntries[i][j] = n;
      }
    }
  }

  Q10PsiPhi::~Q10PsiPhi()
  {
    if (nrEntries)
    {
      for (int i = 0; i < psi->getNumber(); i++)
      {
        delete [] nrEntries[i];
        delete [] values[i];
        delete [] k[i];
      }
      delete [] nrEntries;
      delete [] values;
      delete [] k;

      delete [] val_vec;
      delete [] k_vec;
    }
  }

  bool Q10PsiPhi::operator==(const Q10PsiPhi& q10pp) const
  {
    return (q10pp.psi == psi && q10pp.phi == phi && q10pp.quadrature == quadrature);
  }


  const Q10PsiPhi* Q10PsiPhi::provideQ10PsiPhi(const BasisFunction* ps,
      const BasisFunction* ph,
      const Quadrature* quadrat)
  {
    std::list<Q10PsiPhi*>::iterator list;

    if (!ps && !ph) return NULL;

    if (!ps)  ps = ph;
    if (!ph)  ph = ps;

    if (!quadrat)  quadrat =
        Quadrature::provideQuadrature(ps->getDim(),
                                      ps->getDegree() +
                                      ph->getDegree() - 1);

    compareQPsiPhi<Q10PsiPhi> comp(ps,ph,quadrat);

    /****************************************************************************/
    /*  look for an existing entry in the list                                  */
    /****************************************************************************/

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q10PsiPhi* newQ10PsiPhi = new Q10PsiPhi(ps,ph,quadrat);
      preList.push_back(newQ10PsiPhi);
      return newQ10PsiPhi;
    }
    else
    {
      return *list;
    }
  }


  Q01PsiPhi::Q01PsiPhi(const BasisFunction* ps, const BasisFunction* ph,
                       const Quadrature* quadrat)
    : psi(ps),phi(ph),quadrature(quadrat)
  {
    FUNCNAME_DBG("Q01PsiPhi::Q01PsiPhi");
    const FastQuadrature*      q_phi, *q_psi;
    int                       i, j, ll, n=0, iq, all_entries, n_psi, n_phi;
    double                      val, grdj, psii;
    int                       d=quadrature->getDim();

    int numPoints = quadrature->getNumPoints();

    if (!psi && !phi)
    {
      nrEntries=NULL;
      l=NULL;
      values=NULL;
    }

    if (!psi)  psi = phi;
    if (!phi)  phi = psi;

    if (!quadrature)
      quadrature = Quadrature::provideQuadrature(psi->getDim(),
                   psi->getDegree() +
                   phi->getDegree() - 1);

    /****************************************************************************/
    /*  create a new one                                                        */
    /****************************************************************************/
    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_GRD_PHI);
    n_phi = phi->getNumber();
    q_phi = FastQuadrature::provideFastQuadrature(phi, *quadrature, INIT_PHI);

    nrEntries = new int* [n_psi];
    values = new double** [n_psi];
    l = new int** [n_psi];
    for (int i = 0; i < n_psi; i++)
    {
      nrEntries[i] = new int[n_phi];
      values[i] = new double*[n_phi];
      l[i] = new int* [n_phi];
    }


    /****************************************************************************/
    /*  compute first the number of all non zero entries                        */
    /****************************************************************************/

    all_entries = 0;
    for(i = 0; i < n_psi; i++)
    {
      for(j = 0; j < n_phi; j++)
      {
        for(ll =  0; ll < d+1; ll++)
        {
          for(val = iq = 0; iq < numPoints; iq++)
          {
            grdj = q_phi->getGradient(iq,j,ll);
            psii = q_psi->getPhi(iq,i);

            val += quadrature->getWeight(iq)*grdj*psii;
          }
          if(math::abs(val) > TOO_SMALL)
            all_entries++;
        }
      }
    }

    /****************************************************************************/
    /* now, access memory for all information                                   */
    /****************************************************************************/

    val_vec = new double[all_entries];
    l_vec = new int[all_entries];

    /****************************************************************************/
    /* and now, fill information                                                */
    /****************************************************************************/

    for(i = 0; i < n_psi; i++)
    {
      for(j = 0; j < n_phi; j++)
      {
        values[i][j] = val_vec;
        l[i][j]     = l_vec;

        for(n = ll =  0; ll < d+1; ll++)
        {
          for(val = iq = 0; iq < numPoints; iq++)
          {
            grdj = q_phi->getGradient(iq,j,ll);
            psii = q_psi->getPhi(iq,i);

            val += quadrature->getWeight(iq)*grdj*psii;
          }
          if (math::abs(val) > TOO_SMALL)
          {
            TEST_EXIT_DBG(all_entries > 0)("now more entries found than counted before\n");
            all_entries--;
            n++;
            *val_vec++ = val;
            *l_vec++ = ll;
          }
        }
        nrEntries[i][j] = n;
      }
    }
  }

  Q01PsiPhi::~Q01PsiPhi()
  {
    if (nrEntries)
    {
      for (int i = 0; i < psi->getNumber(); i++)
      {
        delete [] nrEntries[i];
        delete [] values[i];
        delete [] l[i];
      }
      delete [] nrEntries;
      delete [] values;
      delete [] l;

      delete [] val_vec;
      delete [] l_vec;
    }
  }

  bool Q01PsiPhi::operator==(const Q01PsiPhi& q01pp) const
  {
    return (q01pp.psi == psi && q01pp.phi == phi && q01pp.quadrature == quadrature);
  }


  const Q01PsiPhi* Q01PsiPhi::provideQ01PsiPhi(const BasisFunction* ps,
      const BasisFunction* ph,
      const Quadrature*    quadrat)
  {
    std::list<Q01PsiPhi*>::iterator list;

    if (!ps && !ph) return NULL;

    if (!ps)  ps = ph;
    if (!ph)  ph = ps;

    if (!quadrat) quadrat = Quadrature::provideQuadrature(ps->getDim(),
                              ps->getDegree() +
                              ph->getDegree() - 1);

    compareQPsiPhi<Q01PsiPhi> comp(ps,ph,quadrat);

    /****************************************************************************/
    /*  look for an existing entry in the list                                  */
    /****************************************************************************/

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q01PsiPhi* newQ01PsiPhi = new Q01PsiPhi(ps,ph,quadrat);
      preList.push_back(newQ01PsiPhi);
      return newQ01PsiPhi;
    }
    else
    {
      return *list;
    }
  }

  /****************************************************************************/
  /*  first order term:                                                       */
  /****************************************************************************/


  Q00PsiPhi::Q00PsiPhi(const BasisFunction* ps,
                       const BasisFunction* ph,
                       const Quadrature* quadrat)
    : psi(ps),
      phi(ph),
      quadrature(quadrat)
  {

    const FastQuadrature*      q_phi, *q_psi;
    int                       i, j,iq, n_psi, n_phi;
    double                      val;

    int numPoints = quadrature->getNumPoints();

    if (!psi && !phi)
    {
      values=NULL;
    }

    if (!psi)  psi = phi;
    if (!phi)  phi = psi;

    if (!quadrature)
      quadrature = Quadrature::provideQuadrature(psi->getDim(),
                   psi->getDegree() +
                   phi->getDegree());

    /****************************************************************************/
    /*  create a new one                                                        */
    /****************************************************************************/
    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_PHI);
    n_phi = phi->getNumber();
    q_phi = FastQuadrature::provideFastQuadrature(phi, *quadrature, INIT_PHI);

    values = new double*[n_psi];
    for (int i = 0; i < n_psi; i++)
      values[i] = new double[n_phi];


    /****************************************************************************/
    /*  compute first the number of all non zero entries                        */
    /****************************************************************************/

    for(i = 0; i < n_psi; i++)
    {
      for(j = 0; j < n_phi; j++)
      {
        for(val = iq = 0; iq < numPoints; iq++)
          val += quadrature->getWeight(iq)*q_psi->getPhi(iq,i)*q_phi->getPhi(iq,j);

        if (math::abs(val) < TOO_SMALL)
          values[i][j]=0.0;
        else
          values[i][j]=val;
      }
    }
  }

  Q00PsiPhi::~Q00PsiPhi()
  {
    for (int i = 0; i < psi->getNumber(); i++)
      delete [] values[i];

    delete [] values;
  }


  Q00PsiPhi* Q00PsiPhi::provideQ00PsiPhi(const BasisFunction* ps,
                                         const BasisFunction* ph,
                                         const Quadrature*    quadrat)
  {
    std::list<Q00PsiPhi*>::iterator list;

    if (!ps && !ph) return NULL;

    if (!ps)  ps = ph;
    if (!ph)  ph = ps;

    if (!quadrat)
      quadrat = Quadrature::provideQuadrature(ps->getDim(),
                                              ps->getDegree()+ph->getDegree());

    compareQPsiPhi<Q00PsiPhi> comp(ps,ph,quadrat);

    /****************************************************************************/
    /*  look for an existing entry in the list                                  */
    /****************************************************************************/

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q00PsiPhi* newQ00PsiPhi = new Q00PsiPhi(ps,ph,quadrat);
      preList.push_back(newQ00PsiPhi);
      return newQ00PsiPhi;
    }
    else
    {
      return *list;
    }
  }


  Q0Psi::Q0Psi(const BasisFunction* ps, const Quadrature* quadrat)
    : psi(ps), quadrature(quadrat)
  {
    const FastQuadrature* q_psi;

    int iq, n_psi;
    double val, psii = 0.0;
    int numPoints = quadrature->getNumPoints();

    if (!psi)
      values = NULL;

    if (!quadrature)
      quadrature = Quadrature::provideQuadrature(psi->getDim(), 2*psi->getDegree());

    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_PHI);

    values = new double[n_psi];

    for (int i = 0; i < n_psi; i++)
    {
      for (val = iq = 0; iq < numPoints; iq++)
      {
        psii = q_psi->getPhi(iq,i);
        val += quadrature->getWeight(iq)*psii;
      }

      if (math::abs(val) < TOO_SMALL)
        values[i]=0.0;
      else
        values[i]=val;
    }
  }

  Q0Psi::~Q0Psi()
  {
    delete [] values;
  }

  Q0Psi* Q0Psi::provideQ0Psi(const BasisFunction* ps, const Quadrature* quadrat)
  {
    std::list<Q0Psi*>::iterator list;

    if (!ps) return NULL;
    if (!quadrat)  quadrat = Quadrature::provideQuadrature(ps->getDim(),
                               2*ps->getDegree());

    compareQPsi<Q0Psi> comp(ps, quadrat);

    /****************************************************************************/
    /*  look for an existing entry in the list                                  */
    /****************************************************************************/

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q0Psi* newQ0Psi = new Q0Psi(ps, quadrat);
      preList.push_back(newQ0Psi);
      return newQ0Psi;
    }
    else
    {
      return *list;
    }
  }


  Q1Psi::Q1Psi(const BasisFunction* ps, const Quadrature* quadrat)
    : psi(ps), quadrature(quadrat),nrEntries(NULL),values(NULL),k(NULL)
  {
    FUNCNAME_DBG("Q1Psi::Q1Psi");
    const FastQuadrature* q_psi;
    int                       i, lk, n, iq, all_entries, n_psi;
    double                    val, grdi;

    if (!quadrature)  quadrature = Quadrature::provideQuadrature(psi->getDim(),
                                     2*psi->getDegree() - 1);

    n_psi = psi->getNumber();
    q_psi = FastQuadrature::provideFastQuadrature(psi, *quadrature, INIT_GRD_PHI);

    nrEntries = new int[n_psi];
    values = new double*[n_psi];
    k = new int* [n_psi];

    //****************************************************************************
    //*  compute first the number of all non zero entries                        *
    //****************************************************************************

    int numPoints = quadrature->getNumPoints();
    int d = psi->getDim();

    all_entries = 0;
    for(i = 0; i < n_psi; i++)
    {
      for(lk = 0; lk < d+1; lk++)
      {
        for(val = iq = 0; iq < numPoints; iq++)
        {
          grdi = q_psi->getGradient(iq,i,lk);
          //psii = q_psi->getPhi(iq, i);
          val += quadrature->getWeight(iq)*grdi;
        }
        if (math::abs(val) > TOO_SMALL)
          all_entries++;
      }
    }



    //****************************************************************************
    //* now, access memory for all information                                   *
    //****************************************************************************
    allEntries = all_entries;

    val_vec = new double[all_entries];
    k_vec = new int[all_entries];

    //***************************************************************************
    // and now, fill information                                                *
    //***************************************************************************

    for(i = 0; i < n_psi; i++)
    {
      values[i] = val_vec;
      k[i] = k_vec;

      for(n = lk = 0; lk < d+1; lk++)
      {
        for(val = iq = 0; iq < numPoints; iq++)
        {
          grdi = q_psi->getGradient(iq,i,lk);
          //psii = q_psi->getPhi(iq, i);
          val += quadrature->getWeight(iq)*grdi;
        }
        if (math::abs(val) > TOO_SMALL)
        {
          TEST_EXIT_DBG(all_entries > 0)("now more entries found than counted before\n");
          all_entries--;
          n++;
          *val_vec = val;
          val_vec++;
          *k_vec = lk;
          k_vec++;
        }
        nrEntries[i] = n;
      }
    }
  }

  Q1Psi::~Q1Psi()
  {
    if (nrEntries)
    {
      delete [] nrEntries;
      delete [] values;
      delete [] k;
      delete [] val_vec;
      delete [] k_vec;
    }
  }


  bool Q1Psi::operator==(const Q1Psi& q1p) const
  {
    return (q1p.psi == psi && q1p.quadrature == quadrature);
  }


  const Q1Psi* Q1Psi::provideQ1Psi(const BasisFunction* ps,
                                   const Quadrature*    quadrat)
  {
    std::list<Q1Psi*>::iterator list;

    if (!ps) return NULL;

    if (!quadrat)
      quadrat = Quadrature::provideQuadrature(ps->getDim(),
                                              2 * ps->getDegree() - 1);

    compareQPsi<Q1Psi> comp(ps, quadrat);

    //***************************************************************************
    //  look for an existing entry in the list                                  *
    //***************************************************************************

    list = find_if(preList.begin(),
                   preList.end(),
                   comp);

    if (list==preList.end())
    {
      Q1Psi* newQ1Psi = new Q1Psi(ps, quadrat);
      preList.push_back(newQ1Psi);
      return newQ1Psi;
    }
    else
    {
      return *list;
    }
  }

} // end namespace AMDiS
