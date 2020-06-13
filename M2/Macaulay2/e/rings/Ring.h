//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "../ringelem.hpp"
#include <type_traits>

static constexpr bool RZ_CRTP = false;
class ARing;
class CoefficientRingR;

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

  M2_arrayint m_heft_vector{nullptr};
  // This vector, if NULL, and if there are any variables in the ring imply that
  // the heft vector should be taken as the default: the first degree should be
  // used
  // If non-NULL, this should dot to a positive value for every variable in the
  // ring.
  // Question: does this include coefficient variables in the ring?

  ring_elem m_zeroV;          /// zero element in the ring
  ring_elem m_oneV;
  ring_elem m_minus_oneV;

  mutable ring_elem m_non_unit;
  mutable int m_isfield;  // 1: means yes, or declared yes.
  // 0: means not declared to be so.
  // -1: means a non unit was found, and that _non_unit contains
  //    a non-zero non-unit.
  // If a non unit is found, then _isfield is set to -1.

  const ARing *getARing() const { return m_AR; }
  const ARing *m_AR;

 public:
  long characteristic() const { return m_char; }
  const PolynomialRing* get_degree_ring() const { return m_degree_ring; }
  ring_elem one() const { return m_oneV; }
  ring_elem minus_one() const { return m_minus_oneV; }
  ring_elem zero() const { return m_zeroV; }

  bool is_field() const;
  bool declare_field();  // if false is returned, then an ERROR has been set.
  ring_elem get_non_unit() const;
  void set_non_unit(ring_elem zero_div) const;

  M2_arrayint get_heft_vector() const
  {
    return m_heft_vector;
  }  // This CAN BE NULL

 public:
  ring_elem power(ring_elem f, int n) const override { return crtp()->impl_power(f, n); }
//  ring_elem mult(ring_elem f, ring_elem g) override { return crtp()->impl_mult(f, g); }
  ring_elem invert(ring_elem f) const override { throw std::runtime_error("inverse not supported in this ring"); }

protected:
  [[maybe_unused]] ring_elem impl_power(ring_elem f, int n) const;
// [[maybe_unused]] ring_elem impl_mult(ring_elem f, ring_elem g) const;
};

// #include "Ring.hpp"  // safe but confuses CLion
