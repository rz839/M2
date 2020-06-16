//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "Ring.h"
#include "RingEnum.h"

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
