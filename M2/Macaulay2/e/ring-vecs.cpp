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

namespace M2::bugfix {
bool check_nterm_multiples(const PolyRing *R,
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
  return !f && !g;
}
}  // M2::bugfix

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
