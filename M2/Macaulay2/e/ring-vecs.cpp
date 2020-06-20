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

vec Ring::vec_diff(vec v, int rankFw, vec w, int use_coeff) const
// rankFw is the rank of the free module corresponding to w.
{
  vec result = NULL;
  for (; v != NULL; v = v->next)
    for (vecterm *p = w; p != NULL; p = p->next)
      {
        ring_elem a = diff(v->coeff, p->coeff, use_coeff);
        if (is_zero(a))
          {
            remove(a);
            continue;
          }
        vecterm *t = new_vec();
        t->comp = rankFw * v->comp + p->comp;
        t->coeff = a;
        t->next = result;
        result = t;
      }
  vec_sort(result);
  return result;
}

int Ring::vec_in_subring(int nslots, const vec v) const
{
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0 || v == NULL) return true;
  const Monoid *M = PR->getMonoid();
  for (vec w = v; w != NULL; w = w->next)
    if (!M->in_subring(nslots, PR->lead_flat_monomial(w->coeff))) return false;
  return true;
}

void Ring::vec_degree_of_var(int n, const vec v, int &lo, int &hi) const
{
  if (v == NULL)
    {
      ERROR("attempting to find degree of zero vector");
      return;
    }
  degree_of_var(n, v->coeff, lo, hi);
  for (vec w = v->next; w != 0; w = w->next)
    {
      int lo1, hi1;
      degree_of_var(n, w->coeff, lo1, hi1);
      if (lo1 < lo) lo = lo1;
      if (hi1 > hi) hi = hi1;
    }
}

vec Ring::vec_divide_by_var(int n, int d, const vec v) const
{
  vecterm head;
  vecterm *result = &head;
  for (vec w = v; w != 0; w = w->next)
    {
      ring_elem a = divide_by_var(n, d, w->coeff);
      if (!is_zero(a))
        {
          vec t = make_vec(w->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = 0;
  return head.next;
}

vec Ring::vec_divide_by_expvector(const int *exp, const vec v) const
{
  vecterm head;
  vecterm *result = &head;
  for (vec w = v; w != 0; w = w->next)
    {
      ring_elem a = divide_by_expvector(exp, w->coeff);
      if (!is_zero(a))
        {
          vec t = make_vec(w->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = 0;
  return head.next;
}

//////////////////////////////////////////////
//  Homogeniety and the grading //////////////
//////////////////////////////////////////////

bool Ring::vec_multi_degree(const FreeModule *F, const vec f, int *degf) const
// Returns true if the element is homogeneous
// Sets degf to be the highest degree found (actually, the join of the
//   degree vectors occuring).
{
  int *degv;
  degree_monoid()->one(degf);
  if (f == NULL) return true;
  bool result = multi_degree(f->coeff, degf);
  degree_monoid()->mult(degf, F->degree(f->comp), degf);
  degv = degree_monoid()->make_one();

  for (vec v = f->next; v != 0; v = v->next)
    {
      bool ishom = multi_degree(v->coeff, degv);
      result = result && ishom;
      degree_monoid()->mult(degv, F->degree(v->comp), degv);

      if (0 != degree_monoid()->compare(degf, degv))
        {
          result = false;
          degree_monoid()->lcm(degf, degv, degf);
        }
    }
  degree_monoid()->remove(degv);
  return result;
}

void Ring::vec_degree(const FreeModule *F, const vec f, int *degf) const
{
  vec_multi_degree(F, f, degf);
}

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

void Ring::vec_degree_weights(const FreeModule *F,
                              const vec f,
                              M2_arrayint wts,
                              int &lo,
                              int &hi) const
{
  vecterm *t = f;
  if (t == NULL)
    {
      lo = hi = 0;
      return;
    }
  degree_weights(t->coeff, wts, lo, hi);
  lo += F->primary_degree(t->comp);
  hi += F->primary_degree(t->comp);
  for (t = t->next; t != NULL; t = t->next)
    {
      int lo1, hi1;
      degree_weights(t->coeff, wts, lo1, hi1);
      lo1 += F->primary_degree(t->comp);
      hi1 += F->primary_degree(t->comp);
      if (hi1 > hi) hi = hi1;
      if (lo1 < lo) lo = lo1;
    }
}

vec Ring::vec_homogenize(const FreeModule *F,
                         const vec f,
                         int v,
                         int d,
                         M2_arrayint wts) const
// Any terms which can't be homogenized are silently set to 0
{
  vecterm head;
  vecterm *result = &head;
  assert(wts->array[v] != 0);
  // If an error occurs, then return 0, and set ERROR

  for (vec w = f; w != 0; w = w->next)
    {
      int e = F->primary_degree(w->comp);
      ring_elem a = homogenize(w->coeff, v, d - e, wts);
      if (!is_zero(a))
        {
          result->next = make_vec(w->comp, a);
          result = result->next;
        }
    }
  result->next = 0;
  return head.next;
}

vec Ring::vec_homogenize(const FreeModule *F,
                         const vec f,
                         int v,
                         M2_arrayint wts) const
{
  vecterm *result = NULL;
  if (f == NULL) return result;
  int lo, hi;
  vec_degree_weights(F, f, wts, lo, hi);
  assert(wts->array[v] != 0);
  int d = (wts->array[v] > 0 ? hi : lo);
  return vec_homogenize(F, f, v, d, wts);
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

ring_elem Ring::vec_content(vec f) const
{
  if (f == 0) return zero();
  ring_elem c = content(f->coeff);
  for (vec t = f->next; t != 0; t = t->next) lower_content(c, t->coeff);
  return c;
}

vec Ring::vec_divide_by_given_content(vec f, ring_elem c) const
{
  if (f == 0) return 0;
  vecterm head;
  vec result = &head;
  for (const vecterm *p = f; p != 0; p = p->next)
    {
      vec w = new_vec();
      result->next = w;
      result = w;
      w->comp = p->comp;
      w->coeff = divide_by_given_content(p->coeff, c);
    }
  result->next = 0;
  return head.next;
}

vec Ring::vec_divide_by_content(vec f) const
{
  if (f == 0) return 0;
  ring_elem c = vec_content(f);
  return vec_divide_by_given_content(f, c);
}

ring_elem Ring::vec_split_off_content(vec f, vec &result) const
{
  ring_elem c = vec_content(f);
  if (f == 0)
    result = 0;
  else
    result = vec_divide_by_given_content(f, c);
  return c;
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
