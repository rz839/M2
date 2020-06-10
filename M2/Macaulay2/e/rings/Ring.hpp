//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "Ring.h"

//template <typename D>
//ring_elem RingBase<D>::impl_mult(ring_elem f, ring_elem g)
//{
//  return ring_elem(0);
//}

template <typename D>
ring_elem RingBase<D>::impl_power(ring_elem f, int n) const
{
  if (n<0)
    throw std::runtime_error("elements cannot be inverted in this ring");

  ring_elem base = crtp()->copy(f);
  ring_elem res = crtp()->from_long(1);

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
