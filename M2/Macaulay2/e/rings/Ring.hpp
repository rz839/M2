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
