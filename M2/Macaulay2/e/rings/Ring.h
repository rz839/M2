//
// Created by Radoslav Zlatev on 6/10/20.
//

#pragma once

#include "../ringelem.hpp"
#include "RingEnum.h"
#include <type_traits>
#include <exception>

static constexpr bool RZ_CRTP = false;
class ARing;
class CoefficientRingR;

class IRing
{
/*******************************************************************************
 *   TRAITS
 ******************************************************************************/
public:
  virtual bool has_trait(const M2::RingTrait trait) const
  {
    switch (trait)
      {
        case M2::RingTrait::BASIC_RING:
        case M2::RingTrait::COMMUTATIVE:
        case M2::RingTrait::GRADED:
          return true;
        default:
          return false;
      }
  }
  virtual M2::RingTypeId get_typeid() const { return M2::RingTypeId::RING; }

  virtual M2::RingID ringID() const { return M2::ring_old; }
  virtual bool is_basic_ring() const
  {
    return true;
  }  // The default is to be a basic ring.
  virtual bool isFinitePrimeField() const { return false; }
  virtual bool is_ZZ() const { return false; }
  virtual bool is_QQ() const { return false; }
  virtual bool is_RRR() const { return false; }
  virtual bool is_CCC() const { return false; }
  virtual bool is_fraction_field() const { return false; }
  virtual bool is_fraction_poly_ring() const { return false; }
  // returns true if this ring has fractions.  This includes
  // polynomial rings over QQ, polynomial rings over fraction fields,
  // fraction rings, and localizations.
  // If this returns true, then 'get_denominator_ring()' returns non-NULL value.
  //

  virtual bool is_poly_ring() const { return false; }
  // Returns true if this is a polynomial ring, possibly with fractions
  // and possibly with quotient ideal, and possibly with non-commutative
  // multiplication.  Equivalent to (cast_to_PolynomialRing() != 0).

  virtual bool is_commutative_ring() const { return true; }
  // Returns true iff this is a commutative ring.

  virtual bool is_quotient_ring() const { return false; }
  // Returns true if this is a polynomial ring, (possibly with fractions),
  // with a quotient ideal.  This could be a non-commutative ring
  // with skew-commutative, Weyl algebra, or other multiplication.

  virtual bool is_weyl_algebra() const { return false; }
  // Returns true if this is a polynomial ring (possibly with quotient)
  // (possibly with ZZ fractions, or other commutative fractions)
  // but with Weyl algebra multiplication on some of the variables.

  virtual bool is_skew_commutative_ring() const { return false; }
  // Returns true if this is a polynomial ring (possibly with quotient)
  // (possibly with ZZ fractions, or other commutative fractions)
  // but with some variables anti-commuting.

  virtual bool is_solvable_algebra() const { return false; }
  virtual bool is_graded() const { return true; }
  // Is this ring graded, with the given grading?
  // ZZ, QQ, ZZ/p, GF, RR, ... are all graded.
  // polynomial rings are graded
  // Weyl algebras can be graded or not
  // quotient polynomial rings can be graded or not.

  typedef enum { COEFF_ZZ, COEFF_QQ, COEFF_BASIC } CoefficientType;
  virtual CoefficientType coefficient_type() const { return COEFF_BASIC; }
  // What the ultimate coefficient type is.  ZZ, QQ, finite fields return these
  // three values.  Fraction fields return their ultimate value, as do poly
  // rings.

  virtual bool has_associate_divisors() const { return true; }
  // There are only a few rings which do not have such divisors: frac rings
  //   over quotients of poly rings.

/*******************************************************************************
 *   ARITHMETIC
 ******************************************************************************/

  // For some finite fields, if a = (getGenerator())^r, return r.
  // If it is not implemented for this ring, an exception is thrown
  // If a is zero, then r is set to -1.
  virtual long discreteLog(const ring_elem &a) const
  {
    throw exc::engine_error("cannot compute discrete logarithm in this ring");
  }
public:
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
  const Monoid * degree_monoid() const;
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

  /** @name Vector Methods *****************************************************
   *
   * These are methods for vector operations. Move to an adaptor or a decorator
   * class later.
   */
  vec new_vec() const { return new vecterm; }
};

#include "Ring.hpp"  // safe but confuses CLion
