//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "../ringelem.hpp"
#include <type_traits>

static constexpr bool RZ_CRTP = false;

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

protected:
  long m_char{0};
  const PolynomialRing *m_degree_ring;

  ring_elem m_zeroV;          /// zero element in the ring
  ring_elem m_oneV;
  ring_elem m_minus_oneV;

public:
  long characteristic() const { return m_char; }
  const PolynomialRing* get_degree_ring() const { return m_degree_ring; }
  ring_elem one() const { return m_oneV; }
  ring_elem minus_one() const { return m_minus_oneV; }
  ring_elem zero() const { return m_zeroV; }


 public:
  ring_elem power(ring_elem f, int n) const override { return crtp()->impl_power(f, n); }
//  ring_elem mult(ring_elem f, ring_elem g) override { return crtp()->impl_mult(f, g); }
  ring_elem invert(ring_elem f) const override { throw std::runtime_error("inverse not supported in this ring"); }

protected:
  [[maybe_unused]] ring_elem impl_power(ring_elem f, int n) const;
// [[maybe_unused]] ring_elem impl_mult(ring_elem f, ring_elem g) const;
};

// #include "Ring.hpp"  // safe but confuses CLion
