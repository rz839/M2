//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "Ring.h"
#include "base/bugfix.h"

template <typename D>
const Monoid* RingBase<D>::degree_monoid() const
{
  return crtp()->impl_degree_monoid();
}

//template <typename D>
//ring_elem RingBase<D>::impl_mult(ring_elem f, ring_elem g)
//{
//  return ring_elem(0);
//}

template <typename D>
ring_elem RingBase<D>::impl_power(ring_elem f, int n) const
{
  ring_elem base = crtp()->copy(f);
  ring_elem res = crtp()->from_long(1);

  if (n<0)
  {
    base = crtp()->invert(base);  // could throw
    n = -n;
  }

  while (n)
  {
    if (n&1)
      res = crtp()->mult(res, base);
    n >>= 1;
    if (!n)
      break;
    base = crtp()->mult(base, base);
  }
  return res;
}

template <typename D>
ring_elem RingBase<D>::impl_power(const ring_elem gg, mpz_srcptr m) const
{
  ring_elem ff = gg;
  int cmp = mpz_sgn(m);
  if (cmp == 0) return one();
  mpz_t n;
  mpz_init_set(n, m);
  if (cmp < 0)
    {
      mpz_neg(n, n);
      ff = invert(ff);
      if (is_zero(ff))
        {
          ERROR(
              "either element not invertible, or no method available to "
              "compute its inverse");
          return ff;
        }
    }
  ring_elem prod = from_long(1);
  ring_elem base = copy(ff);
  ring_elem tmp;

  for (;;)
    {

      if (mpz_tstbit(n, 0u))
        {
          tmp = mult(prod, base);
          prod = tmp;
        }
      mpz_tdiv_q_2exp(n, n, 1);
      if (mpz_sgn(n) == 0)
        {
          mpz_clear(n);
          return prod;
        }
      else
        {
          tmp = mult(base, base);
          base = tmp;
        }
    }
}

template <typename D>
bool RingBase<D>::is_field() const
{
  return m_isfield == 1;
}

template <typename D>
bool RingBase<D>::declare_field()
{
  if (m_isfield >= 0)
  {
    m_isfield = 1;
    return true;
  }
  else
  {
    ERROR("attempting to declare a ring with known non-units to be a field");
    return false;
  }
}

template <typename D>
ring_elem RingBase<D>::get_non_unit() const
{
  if (m_isfield >= 0) return zero();
  return crtp()->copy(m_non_unit);
}

template <typename D>
void RingBase<D>::set_non_unit(ring_elem non_unit) const
{
  if (m_isfield == 1)  // i.e. declared to be a field
    ERROR("a non unit was found in a ring declared to be a field");
  m_isfield = -1;
  m_non_unit = non_unit;
}

//template <typename D>
//vec RingBase<D>::new_vec() const
//{
//  return 0;
//}

template <typename D>
ring_elem RingBase<D>::preferred_associate(ring_elem f) const
{
  // Here we assume that 'this' is a field:
  if (crtp()->is_zero(f))
    return crtp()->from_long(1);
  return crtp()->invert(f);
}

template <typename D>
bool RingBase<D>::lower_associate_divisor(ring_elem &f, const ring_elem g) const
{
  if (crtp()->is_zero(f))
  {
    f = g;
    return !crtp()->is_zero(f);
  }
  return true;
}

template <typename D>
void RingBase<D>::lower_content(ring_elem &result, ring_elem g) const
{
  // default implementation
  // The default implementation here ASSUMES that result and g are in the same
  // ring!
  if (crtp()->is_zero(result)) result = g;
}

template <typename D>
ring_elem RingBase<D>::content(ring_elem f) const
// default implementation
{
  return f;
}

template <typename D>
ring_elem RingBase<D>::content(ring_elem f, ring_elem g) const
// default implementation
{
  lower_content(f, g);
  return f;
}

template <typename D>
ring_elem RingBase<D>::divide_by_given_content(ring_elem f, ring_elem c) const
// default implementation
{
  // The default implementation here ASSUMES that f and c are in the same ring!
  return crtp()->divide(f, c);
}

template <typename D>
ring_elem RingBase<D>::divide_by_content(ring_elem f) const
{
  ring_elem c = content(f);
  return divide_by_given_content(f, c);
}

template <typename D>
ring_elem RingBase<D>::split_off_content(ring_elem f, ring_elem &result) const
{
  ring_elem c = content(f);
  result = divide_by_given_content(f, c);
  return c;
}

template <typename D>
bool RingBase<D>::is_equal(const vecterm *a, const vecterm *b) const
{
  for (;; a = a->next, b = b->next)
  {
    if (!a) return b == nullptr;
    if (!b) return false;
    if (a->comp != b->comp) return false;
    if (!this->is_equal(a->coeff, b->coeff)) return false;
  }
}

template <typename D>
ring_elem RingBase<D>::divide(const ring_elem f, const ring_elem g) const
{
  return crtp()->mult(f, crtp()->invert(g));
}

template <typename D>
int RingBase<D>::index_of_var(const ring_elem a) const { return -1; }

template <typename D>
M2_arrayint RingBase<D>::support(const ring_elem a) const
{
  M2_arrayint result = M2_makearrayint(0);
  return result;
}

template <typename D>
ring_elem RingBase<D>::divide_by_var(int n, int d, const ring_elem a) const
{
  if (d == 0) return a;
  return from_long(0);
}

template <typename D>
ring_elem RingBase<D>::homogenize(const ring_elem f, int, int deg, M2_arrayint) const
{
  if (deg != 0) ERROR("homogenize: no homogenization exists");
  return f;
}

template <typename D>
void RingBase<D>::degree(const ring_elem, int *d) const { degree_monoid()->one(d); }

template <typename D>
bool RingBase<D>::multi_degree(const ring_elem f, int *d) const
// returns true iff f is homogeneous
{
  degree_monoid()->one(d);
  return true;
}

template <typename D>
ring_elem RingBase<D>::remainder(const ring_elem f, const ring_elem g) const
{
  if (crtp()->is_zero(g)) return f;
  return crtp()->zero();
}

template <typename D>
ring_elem RingBase<D>::quotient(const ring_elem f, const ring_elem g) const
{
  if (crtp()->is_zero(g)) return g;
  return crtp()->divide(f, g);
}

template <typename D>
ring_elem RingBase<D>::remainderAndQuotient(const ring_elem f,
                                     const ring_elem g,
                                     ring_elem &quot) const
{
  if (crtp()->is_zero(g))
  {
    quot = g;  // zero
    return f;
  }
  quot = crtp()->divide(f, g);
  return crtp()->zero();
}

template <typename D>
void RingBase<D>::vec_increase_maxnorm(gmp_RRmutable norm, const vec v) const
// If any real number appearing in f has larger absolute value than norm,
// replace norm.
// Default for rings not over RRR or CCC is to do nothing.
{
  for (const vecterm *p = v; p != 0; p = p->next)
    increase_maxnorm(norm, p->coeff);
}

template <typename D>
vec RingBase<D>::vec_zeroize_tiny(gmp_RR epsilon, const vec v) const
{
  vecterm head;
  vec result = &head;
  for (const vecterm *p = v; p != 0; p = p->next)
    {
      ring_elem a = zeroize_tiny(epsilon, p->coeff);
      if (!is_zero(a))
        {
          vec w = new_vec();
          result->next = w;
          result = w;
          w->comp = p->comp;
          w->coeff = a;
        }
    }
  result->next = 0;
  return head.next;
}

template <typename D>
bool RingBase<D>::get_entry(const vecterm *v, int r, ring_elem &result) const
{
  for (const vecterm *p = v; p != 0; p = p->next)
    if (p->comp < r)
      break;
    else if (p->comp == r)
      {
        result = p->coeff;
        return true;
      }
  return false;
}

template <typename D>
ring_elem RingBase<D>::get_entry(vec v, int r) const
{
  while (v != NULL)
    {
      if (v->comp == r) return v->coeff;
      if (v->comp < r) return from_long(0);
      v = v->next;
    }
  return from_long(0);
}

template <typename D>
int RingBase<D>::n_nonzero_terms(const vecterm *v) const
{
  int result = 0;
  for (; v != NULL; v = v->next) result++;
  return result;
}

template <typename D>
vec RingBase<D>::negate_vec(vec v) const
{
  vecterm result;
  vecterm *b = &result;
  for (vecterm *a = v; a != NULL; a = a->next)
    {
      b->next = make_vec(a->comp, negate(a->coeff));
      b = b->next;
    }
  b->next = NULL;
  return result.next;
}

template <typename D>
vec RingBase<D>::add_vec(vec v, vec w) const
{
  vec f = copy_vec(v);
  vec g = copy_vec(w);
  add_vec_to(f, g);
  return f;
}

template <typename D>
vec RingBase<D>::subtract_vec(vec v, vec w) const
{
  vec f = copy_vec(v);
  vec g = negate_vec(w);
  add_vec_to(f, g);
  return f;
}

template <typename D>
vec RingBase<D>::mult_vec(int n, vec v) const
{
  ring_elem f = from_long(n);
  vec result = mult_vec(f, v);
  return result;
}

template <typename D>
vec RingBase<D>::mult_vec(const ring_elem f, const vec w) const
{
  if (is_zero(f)) return NULL;
  vecterm head;
  vec result = &head;
  for (vec v = w; v != 0; v = v->next)
    {
      ring_elem a = mult(f, v->coeff);
      if (!is_zero(a))
        {
          vec t = make_vec(v->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = NULL;
  return head.next;
}

template <typename D>
vec RingBase<D>::rightmult_vec(const vec w, const ring_elem f) const
{
  if (is_zero(f)) return NULL;
  vecterm head;
  vec result = &head;
  for (vec v = w; v != 0; v = v->next)
    {
      ring_elem a = mult(v->coeff, f);
      if (!is_zero(a))
        {
          vec t = make_vec(v->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = NULL;
  return head.next;
}

template <typename D>
vec RingBase<D>::sub_vector(const vecterm *v, M2_arrayint r) const
{
  if (v == 0) return 0;
  // Largest component which occurs in v occurs first.
  VECTOR(int) trans(v->comp + 1);
  for (int i = 0; i < v->comp; i++) trans.push_back(-1);

  for (unsigned j = 0; j < r->len; j++)
    if (r->array[j] >= 0 && r->array[j] <= v->comp) trans[r->array[j]] = j;

  vecterm head;
  vecterm *result = &head;
  for (; v != NULL; v = v->next)
    if (trans[v->comp] != -1)
      {
        result->next = new_vec();
        result = result->next;
        result->next = 0;
        result->coeff = v->coeff;
        result->comp = trans[v->comp];
      }
  result->next = NULL;
  result = head.next;

  vec_sort(result);
  return result;
}

template <typename D>
vec RingBase<D>::component_shift(int n, vec v) const
{
  vecterm head;
  vec result = &head;
  for (const vecterm *p = v; p != 0; p = p->next)
    {
      vec w = new_vec();
      result->next = w;
      result = w;
      w->comp = p->comp + n;
      w->coeff = p->coeff;  // copy is not done
    }
  result->next = 0;
  return head.next;
}

template <typename D>
vec RingBase<D>::tensor_shift(int n, int m, vec v) const
{
  vecterm head;
  vecterm *result = &head;

  for (; v != NULL; v = v->next)
    {
      vec w = new_vec();
      result->next = w;
      result = w;
      w->comp = n * v->comp + m;
      w->coeff = v->coeff;  // copy is not done
    }
  result->next = NULL;
  return head.next;
}

template <typename D>
void RingBase<D>::vec_text_out(buffer &o,
                        const vecterm *v,
                        bool p_one,
                        bool p_plus,
                        bool p_parens) const
{
  if (v == NULL)
    {
      o << "0";
      return;
    }

  p_one = false;
  for (const vecterm *t = v; t != NULL; t = t->next)
    {
      this->elem_text_out(o, t->coeff, p_one, p_plus, p_parens);
      o << "<" << t->comp << ">";
      p_plus = true;
    }
}

///////////////////////////////////////
// Routines which modify a vec ////////
///////////////////////////////////////

template <typename D>
void RingBase<D>::mult_vec_to(vec &v, const ring_elem r, bool opposite_mult) const
{
  if (this->is_zero(r))
    {
      remove_vec(v);
      v = 0;
      return;
    }
  vecterm head;
  head.next = v;
  vec p = &head;
  while (p->next != 0)
    {
      // old version: this->mult_to(p->next->coeff, a);
      ring_elem c;
      if (opposite_mult)
        c = this->mult(p->next->coeff, r);
      else
        c = this->mult(r, p->next->coeff);
      p->next->coeff = c;
      if (this->is_zero(p->next->coeff))
        {
          vec tmp = p->next;
          p->next = tmp->next;
          remove_vec_node(tmp);
        }
      else
        p = p->next;
    }
  v = head.next;
}

template <typename D>
void RingBase<D>::mult_row(vec &v, const ring_elem r, int i, bool opposite_mult) const
{
  vecterm head;
  head.next = v;
  for (vec p = &head; p->next != 0; p = p->next)
    if (p->next->comp < i)
      break;
    else if (p->next->comp == i)
      {
        ring_elem c;
        if (opposite_mult)
          c = mult(p->next->coeff, r);
        else
          c = mult(r, p->next->coeff);
        p->next->coeff = c;
        if (this->is_zero(p->next->coeff))
          {
            vec tmp = p->next;
            p->next = tmp->next;
            remove_vec_node(tmp);
          }
        break;
      }
  v = head.next;
}

template <typename D>
void RingBase<D>::divide_vec_to(vec &v, const ring_elem a) const
{
  if (this->is_zero(a))
    {
      remove_vec(v);
      v = 0;
    }
  vecterm head;
  head.next = v;
  vec p = &head;
  while (p->next != 0)
    {
      // old version: this->mult_to(p->next->coeff, a);
      ring_elem c =
          this->divide(p->next->coeff, a);  // exact or quotient?? MES MES
      p->next->coeff = c;
      if (this->is_zero(p->next->coeff))
        {
          vec tmp = p->next;
          p->next = tmp->next;
          remove_vec_node(tmp);
        }
      else
        p = p->next;
    }
  v = head.next;
}

template <typename D>
void RingBase<D>::divide_row(vec &v, int r, const ring_elem a) const
{
  vecterm head;
  head.next = v;
  for (vec p = &head; p->next != 0; p = p->next)
    if (p->next->comp < r)
      break;
    else if (p->next->comp == r)
      {
        ring_elem c =
            this->divide(p->next->coeff, a);  // exact or quotient?? MES MES
        p->next->coeff = c;
        if (this->is_zero(p->next->coeff))
          {
            vec tmp = p->next;
            p->next = tmp->next;
            remove_vec_node(tmp);
          }
        break;
      }
}

template <typename D>
void RingBase<D>::negate_vec_to(vec &v) const
{
  vec w = v;
  while (w != NULL)
    {
      negate_to(w->coeff);
      w = w->next;
    }
}

template <typename D>
void RingBase<D>::subtract_vec_to(vec &v, vec &w) const
{
  negate_vec_to(w);
  add_vec_to(v, w);
}

template <typename D>
ring_elem RingBase<D>::dot_product(const vecterm *v, const vecterm *w) const
{
  ring_elem result = this->from_long(0);
  while (true)
    {
      if (v == 0) return result;
      if (w == 0) return result;
      if (v->comp > w->comp)
        v = v->next;
      else if (v->comp < w->comp)
        w = w->next;
      else
        {
          ring_elem a = this->mult(v->coeff, w->coeff);
          result = this->add(result, a);
          v = v->next;
          w = w->next;
        }
    }
}

template <typename D>
void RingBase<D>::set_entry(vec &v, int r, ring_elem a) const
{
  vec p;
  bool iszero = this->is_zero(a);
  vecterm head;
  head.next = v;
  for (p = &head; p->next != 0; p = p->next)
    if (p->next->comp <= r) break;

  if (p->next == 0 || p->next->comp < r)
    {
      if (iszero) return;
      vec w = new_vec();
      w->next = p->next;
      w->comp = r;
      w->coeff = a;
      p->next = w;
    }
  else if (p->next->comp == r)
    {
      if (iszero)
        {
          // delete node
          vec tmp = p->next;
          p->next = tmp->next;
          remove_vec_node(tmp);
        }
      else
        p->next->coeff = a;
    }
  v = head.next;
}

template <typename D>
vec RingBase<D>::vec_lead_term(int nparts, const FreeModule *F, vec v) const
{
  // May be over-ridden by subclasses.  In particular, by polynomial classes.
  if (v == 0) return 0;
  return make_vec(v->comp, v->coeff);
}

template <typename D>
ring_elem RingBase<D>::vec_content(vec f) const
{
  if (f == 0) return zero();
  ring_elem c = content(f->coeff);
  for (vec t = f->next; t != 0; t = t->next) lower_content(c, t->coeff);
  return c;
}

template <typename D>
vec RingBase<D>::vec_divide_by_given_content(vec f, ring_elem c) const
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

template <typename D>
vec RingBase<D>::vec_divide_by_content(vec f) const
{
  if (f == 0) return 0;
  ring_elem c = vec_content(f);
  return vec_divide_by_given_content(f, c);
}

template <typename D>
ring_elem RingBase<D>::vec_split_off_content(vec f, vec &result) const
{
  ring_elem c = vec_content(f);
  if (f == 0)
    result = 0;
  else
    result = vec_divide_by_given_content(f, c);
  return c;
}

template <typename D>
vec RingBase<D>::vec_homogenize(const FreeModule *F,
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
      int e = M2::bugfix::primary_degree(F, w->comp);
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

template <typename D>
vec RingBase<D>::vec_homogenize(const FreeModule *F,
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

template <typename D>
void RingBase<D>::vec_degree_weights(const FreeModule *F,
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
  lo += M2::bugfix::primary_degree(F, t->comp);
  hi += M2::bugfix::primary_degree(F, t->comp);
  for (t = t->next; t != NULL; t = t->next)
    {
      int lo1, hi1;
      degree_weights(t->coeff, wts, lo1, hi1);
      lo1 += M2::bugfix::primary_degree(F, t->comp);
      hi1 += M2::bugfix::primary_degree(F, t->comp);
      if (hi1 > hi) hi = hi1;
      if (lo1 < lo) lo = lo1;
    }
}

template <typename D>
void RingBase<D>::vec_degree(const FreeModule *F, const vec f, int *degf) const
{
  vec_multi_degree(F, f, degf);
}


// returns true iff f is homogeneous
template <typename D>
bool RingBase<D>::vec_multi_degree(const FreeModule *F, const vec f, int *degf) const
// Returns true if the element is homogeneous
// Sets degf to be the highest degree found (actually, the join of the
//   degree vectors occuring).
{
  int *degv;
  degree_monoid()->one(degf);
  if (f == NULL) return true;
  bool result = multi_degree(f->coeff, degf);
  degree_monoid()->mult(degf, M2::bugfix::degree(F, f->comp), degf);
  degv = degree_monoid()->make_one();

  for (vec v = f->next; v != 0; v = v->next)
    {
      bool ishom = multi_degree(v->coeff, degv);
      result = result && ishom;
      degree_monoid()->mult(degv, M2::bugfix::degree(F, v->comp), degv);

      if (0 != degree_monoid()->compare(degf, degv))
        {
          result = false;
          degree_monoid()->lcm(degf, degv, degf);
        }
    }
  degree_monoid()->remove(degv);
  return result;
}

template <typename D>
bool RingBase<D>::vec_is_scalar_multiple(vec f, vec g) const
// is df = cg, some scalars c,d?
// These scalars are over the very bottom base field/ZZ.
{
  if (f == NULL) return true;
  if (g == NULL) return true;
  const PolynomialRing *PR = cast_to_PolynomialRing();
  return PR == nullptr;
}

template <typename D>
vec RingBase<D>::vec_remove_monomial_factors(vec f, bool make_squarefree_only) const
{
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0) return copy_vec(f);
  return 0;
}
