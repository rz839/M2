//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "../ringelem.hpp"

class IRing
{
  /* arithmetic
   */
  virtual ring_elem power(ring_elem f, int n) const = 0;
//  virtual ring_elem mult(ring_elem f, ring_elem g) const = 0;
  virtual ring_elem invert(ring_elem f) const = 0;
};

template <typename Derived>
class RingBase : public virtual IRing,
                 public MutableEngineObject  // TODO(RZ): make sure that's OK
{
protected:
  const Derived* crtp() const { return static_cast<const Derived*>(this); }

public:
  ring_elem power(ring_elem f, int n) const override { return crtp()->impl_power(f, n); }
//  ring_elem mult(ring_elem f, ring_elem g) override { return crtp()->impl_mult(f, g); }
  ring_elem invert(ring_elem f) const override { throw std::runtime_error("inverse not supported in this ring"); }

protected:
  [[maybe_unused]] ring_elem impl_power(ring_elem f, int n) const;
// [[maybe_unused]] ring_elem impl_mult(ring_elem f, ring_elem g) const;
};

// #include "Ring.hpp"  // safe but confuses CLion
