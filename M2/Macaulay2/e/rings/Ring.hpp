//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "Ring.h"
#include "RingEnum.h"

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
