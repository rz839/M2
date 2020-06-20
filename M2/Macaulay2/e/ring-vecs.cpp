// Copyright 2003  Michael E. Stillman

#include "ring.hpp"
#include "text-io.hpp"
#include <vector>
#include "matrix.hpp"
#include "geovec.hpp"
#include "ringmap.hpp"
#include "poly.hpp"
//  Notes: ring_elem's are treated as immutable objects: they are not changed,
//  and
// the fact that one cannot change is used throughout.

vec Ring::tensor(const FreeModule *F, vec v, const FreeModule *G, vec w) const
{
  vecHeap H(F);
  for (; v != NULL; v = v->next)
    {
      vec w1 = component_shift(v->comp * M2::bugfix::rank(G), w);
      mult_vec_to(w1, v->coeff, false);
      H.add(w1);
    }
  return H.value();
}


//////////////////////////////////////////////
//  Homogeniety and the grading //////////////
//////////////////////////////////////////////

bool Ring::vec_is_homogeneous(const FreeModule *F, const vec f) const
{
  if (!this->is_graded()) return false;
  if (f == NULL) return true;
  int *d = degree_monoid()->make_one();
  int *e = degree_monoid()->make_one();
  bool result = multi_degree(f->coeff, d);
  if (result)
    {
      degree_monoid()->mult(d, F->degree(f->comp), d);
      for (vecterm *t = f->next; (t != NULL) && result; t = t->next)
        {
          bool ishom = multi_degree(t->coeff, e);
          result = result && ishom;
          if (result)
            {
              degree_monoid()->mult(e, F->degree(t->comp), e);
              if (0 != degree_monoid()->compare(d, e)) result = false;
            }
        }
    }
  degree_monoid()->remove(d);
  degree_monoid()->remove(e);
  return result;
}

//////////////////////////////////////////////
//  Divisibility checks               ////////
//                                    ////////
//////////////////////////////////////////////

bool static check_nterm_multiples(const PolyRing *R,
                                  ring_elem f1,  // in R
                                  ring_elem g1,  // in R
                                  ring_elem c,   // in flat coeffs of R
                                  ring_elem d)   // in flat coeffs of R
{
  Nterm *f;
  Nterm *g;
  const Monoid *M = R->getMonoid();
  const Ring *K = R->getCoefficients();
  for (f = f1, g = g1; f != 0 && g != 0; f = f->next, g = g->next)
    {
      if (M->compare(f->monom, g->monom) != 0) return false;
      ring_elem c1 = K->mult(c, g->coeff);
      ring_elem d1 = K->mult(d, f->coeff);
      int isequal = K->is_equal(c1, d1);
      if (!isequal) return false;
    }
  if (f == NULL && g == NULL) return true;
  return false;
}

bool Ring::vec_is_scalar_multiple(vec f, vec g) const
// is df = cg, some scalars c,d?
// These scalars are over the very bottom base field/ZZ.
{
  if (f == NULL) return true;
  if (g == NULL) return true;
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0) return true;
  const PolyRing *PR1 = PR->getNumeratorRing();
#ifdef DEVELOPMENT
#warning "use numerator only"
#endif
  if (f->comp != g->comp) return false;
  Nterm *f1 = f->coeff;
  Nterm *g1 = g->coeff;
  ring_elem c = f1->coeff;
  ring_elem d = g1->coeff;
  vec p, q;
  for (p = f, q = g; p != NULL && q != NULL; p = p->next, q = q->next)
    {
      if (p->comp != q->comp) return 0;
      if (!check_nterm_multiples(PR1, p->coeff, q->coeff, c, d)) return false;
    }
  if (q == NULL && p == NULL) return true;
  return false;
}

vec Ring::vec_remove_monomial_factors(vec f, bool make_squarefree_only) const
{
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0) return copy_vec(f);
  if (f == 0) return 0;

  int *exp = newarray_atomic(int, PR->n_vars());

  Nterm *t = f->coeff;
  PR->getMonoid()->to_expvector(t->monom, exp);  // Get the process started

  for (vec a = f; a != NULL; a = a->next) monomial_divisor(a->coeff, exp);

  if (make_squarefree_only)
    // Now divide each term by exp[i]-1, if exp[i] >= 2
    for (int i = 0; i < PR->n_vars(); i++)
      if (exp[i] >= 1) exp[i]--;

  vec result = vec_divide_by_expvector(exp, f);

  deletearray(exp);
  return result;
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
