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
class Tower;
class PolyRing;
class PolyRingFlat;
class SchurRing;
class WeylAlgebra;
class SolvableAlgebra;
class RRR;
class CCC;

#define M2_BAD_LEVEL throw exc::engine_error("not implemented at this level");

class IRing
{
public:
  virtual unsigned int computeHashValue(const ring_elem a) const = 0;

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
  virtual bool is_basic_ring() const { return true; }
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

  /*****************************************************************************
   *  DOWN-CASTING - to be removed
   ****************************************************************************/

  virtual const Tower *cast_to_Tower() const { return nullptr; }
  virtual Tower *cast_to_Tower() { return nullptr; }
  virtual const PolynomialRing *cast_to_PolynomialRing() const { return nullptr; }
  virtual PolynomialRing *cast_to_PolynomialRing() { return nullptr; }
  virtual const PolyRing *cast_to_PolyRing() const { return nullptr; }
  virtual PolyRing *cast_to_PolyRing() { return nullptr; }
  virtual const PolyRingFlat *cast_to_PolyRingFlat() const { return nullptr; }
  virtual PolyRingFlat *cast_to_PolyRingFlat() { return nullptr; }

  virtual const SchurRing *cast_to_SchurRing() const { return nullptr; }
  virtual SchurRing *cast_to_SchurRing() { return nullptr; }
  virtual const SolvableAlgebra *cast_to_SolvableAlgebra() const { return nullptr; }
  virtual SolvableAlgebra *cast_to_SolvableAlgebra() { return nullptr; }
  virtual const WeylAlgebra *cast_to_WeylAlgebra() const { return nullptr; }

  virtual RRR *cast_to_RRR() { return nullptr; }
  virtual const RRR *cast_to_RRR() const { return nullptr; }
  virtual CCC *cast_to_CCC() { return nullptr; }
  virtual const CCC *cast_to_CCC() const { return nullptr; }

/*******************************************************************************
 *   OTHER
 ******************************************************************************/
public:
  // Galois Field routines.  These three routines only return non-NULL values
  // if this was created as a Galois field, isom to A = kk[b]/(f(b)), kk = prime
  // field of char p.

  // For some finite fields, if a = (getGenerator())^r, return r.
  // If it is not implemented for this ring, an exception is thrown
  // If a is zero, then r is set to -1.
  virtual long discreteLog(const ring_elem &a) const
  {
    throw exc::engine_error("cannot compute discrete logarithm in this ring");
  }

  // Returns NULL if not a GF.  Returns f(b) in the ring kk[b].  (Variable name
  // might be different)
  virtual const RingElement *getMinimalPolynomial() const { return 0; }
  // Returns NULL if not a GF.  Returns an element of 'this', whose powers give
  // all non-zero elements
  // of the field.
  virtual const RingElement *getGenerator() const
  {
    throw exc::engine_error("not implemented for this ring");
  }

  // Returns the element in the polynomial ring A corresponding to the element
  // a.
  // Returns NULL if not a GF field.
  // Essentially the same as 'lift', except that more information, not readily
  // available, is needed
  // for that call.
  virtual const RingElement *getRepresentation(const ring_elem &a) const
  {
    return 0;
  }

  virtual void text_out(buffer &o) const = 0;

/*******************************************************************************
 *   ARITHMETIC
 ******************************************************************************/

public:
  virtual std::pair<bool, long> coerceToLongInteger(ring_elem a) const = 0;

  virtual ring_elem from_long(long n) const = 0;
  virtual ring_elem from_int(mpz_srcptr n) const = 0;

  /**
   * The next batch of functions perform element conversions.
   * If successful, true is returned and `result` holds the value;
   * otherwise, false is returned and `result` is unchanged.
   */
  virtual bool from_rational(const mpq_srcptr q, ring_elem &result) const = 0;
  virtual bool from_BigReal(gmp_RR a, ring_elem &result) const = 0;
  virtual bool from_BigComplex(gmp_CC z, ring_elem &result) const = 0;
  virtual bool from_double(double a, ring_elem &result) const = 0;
  virtual bool from_complex_double(double re, double im, ring_elem &result) const = 0;

  virtual ring_elem var(int v) const = 0;

  /**
   * Returns an invertible element c of the same ring such that c*f is the
   * preferred associate of the element f.
   * WARNING: The default implementation is for a field.
   */
  virtual ring_elem preferred_associate(ring_elem f) const = 0;

  /**
  // Replaces f with the unit c such that (fx+g)//c is the preferred associate
  //   of fx+g, in the ring A[x], where A is 'this'.
  // Returns false if f will never be changed after this
  // (This happens over ZZ if f is non-zero (therefore 1 or -1, over a finite
  // filed if f != 0,
  // but over QQ will never happen)
  // WARNING: The default implementation is for a field.
   */
  virtual bool lower_associate_divisor(ring_elem &f, ring_elem g) const = 0;

  /**
   * content functions
   */
  virtual void lower_content(ring_elem &c, ring_elem g) const = 0;
  virtual ring_elem content(ring_elem f) const = 0;
  virtual ring_elem content(ring_elem f, ring_elem g) const = 0;
  virtual ring_elem divide_by_given_content(ring_elem f, ring_elem c) const = 0;
  virtual ring_elem divide_by_content(ring_elem f) const = 0;
  virtual ring_elem split_off_content(ring_elem f, ring_elem &result) const = 0;

  virtual bool promote(const Ring *R, const ring_elem f, ring_elem &result) const = 0;
  virtual bool lift(const Ring *R, const ring_elem f, ring_elem &result) const = 0;

  virtual bool is_unit(const ring_elem f) const = 0;
  virtual bool is_zero(const ring_elem f) const = 0;

  virtual bool is_equal(const ring_elem f, const ring_elem g) const = 0;
  virtual bool is_equal(const vecterm *a, const vecterm *b) const = 0;
  virtual int compare_elems(const ring_elem f, const ring_elem g) const = 0;  /// {-1,0,1} for f<=>g, resp.

  virtual ring_elem copy(const ring_elem f) const = 0;
  virtual void remove(ring_elem &f) const = 0;

  virtual ring_elem negate(const ring_elem f) const = 0;
  virtual ring_elem add(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem subtract(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem mult(const ring_elem f, const ring_elem g) const = 0;

  virtual void negate_to(ring_elem &f) const = 0;
  virtual void add_to(ring_elem &f, ring_elem &g) const = 0;
  virtual void subtract_to(ring_elem &f, ring_elem &g) const = 0;
  virtual void mult_to(ring_elem &f, const ring_elem g) const = 0;

  virtual ring_elem power(ring_elem f, int n) const = 0;
  virtual ring_elem power(const ring_elem f, mpz_srcptr n) const = 0;
  virtual ring_elem invert(const ring_elem f) const = 0;
  virtual ring_elem divide(const ring_elem f, const ring_elem g) const = 0;

  // The default version is for a field:
  //   f % 0 is f, otherwise f % g is 0.
  //   f // 0 is 0, otherwise f // g is f/g
  // These three routines: remainder, quotient and remainderAndQuotient
  // satisfy these properties:
  // If r = remainder(f,g), q = quotient(f,g), then
  // (1) f = q*g + r
  // (2) If f is in ideal(g), then r = 0.
  // (3) If g is invertible, then r = 0, and q = f * g^(-1).
  // (4) If the ring is ZZ, then the remainder is "balanced": -[g/2] < r <=
  // [g/2]
  // remainderAndQuotient combines remainder and quotient into one routine.
  virtual ring_elem remainder(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem quotient(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem remainderAndQuotient(const ring_elem f,
                                         const ring_elem g,
                                         ring_elem &quot) const = 0;

  // Constructs elements x and y in the ring s.t. ax + by = 0.  This syzygy is
  // chosen as simply as possible.  For example, over QQ, x is chosen
  // to be positive.  The routine must handle the case when a=0, but can
  // ignore the case when b=0... (Really?)
  virtual void syzygy(const ring_elem a,
                      const ring_elem b,
                      ring_elem &x,
                      ring_elem &y) const = 0;

  virtual void elem_text_out(buffer &o,
                             const ring_elem f,
                             bool p_one = true,
                             bool p_plus = false,
                             bool p_parens = false) const = 0;

/*******************************************************************************
 * POLYNOMIAL ROUTINES
 ******************************************************************************/

  virtual ring_elem eval(const RingMap *map,
                         const ring_elem f,
                         int first_var) const = 0;

  virtual int index_of_var(const ring_elem a) const = 0;
  virtual M2_arrayint support(const ring_elem a) const = 0;

  virtual void monomial_divisor(const ring_elem a, int *exp) const = 0;
  virtual ring_elem diff(ring_elem a, ring_elem b, int use_coeff) const = 0;
  virtual bool in_subring(int nslots, const ring_elem a) const = 0;
  virtual void degree_of_var(int n, const ring_elem a, int &lo, int &hi) const = 0;
  virtual ring_elem divide_by_var(int n, int d, const ring_elem a) const = 0;
  virtual ring_elem divide_by_expvector(const int *exp, const ring_elem a) const = 0;

  virtual ring_elem homogenize(const ring_elem f, int v, int deg, M2_arrayint wts) const = 0;
  virtual ring_elem homogenize(const ring_elem f, int v, M2_arrayint wts) const = 0;

  // Routines expecting a grading.  The default implementation
  // is that the only degree is 0.
  virtual bool is_homogeneous(const ring_elem f) const = 0;
  virtual void degree(const ring_elem f, int *d) const = 0;
  virtual bool multi_degree(const ring_elem f, int *d) const = 0;
  // returns true iff f is homogeneous
  virtual void degree_weights(const ring_elem f,
                              M2_arrayint wts,
                              int &lo,
                              int &hi) const = 0;

  // antipode: for non skew commuting poly rings, this is the identity.
  // Otherwise, this changes the signs of the monomials, implementing the
  // (anti)-isomorphism
  // of the ring with its opposite ring.
  virtual ring_elem antipode(ring_elem f) const { return f; }

/*******************************************************************************
 *   OTHER
 ******************************************************************************/

  virtual ring_elem random() const
  {
    throw exc::engine_error("random scalar elements for this ring are not implemented");
  }

  // These next three routines are only overridden by RRR,CCC,polynomial rings,
  // and quotient rings
  virtual unsigned long get_precision() const { return 0; }  /// if the ring is not over RRR or CCC, returns 0.
  virtual ring_elem zeroize_tiny(gmp_RR epsilon, const ring_elem f) const { return f; } /// Default is to return f itself.
  // If any real number appearing in f has larger absolute value than norm,
  // replace norm.
  // Default for rings not over RRR or CCC is to do nothing.
  // If any real number appearing in f has larger absolute value than norm,
  // replace norm.
  virtual void increase_maxnorm(gmp_RRmutable norm, const ring_elem f) const {};
  // Default is to return f itself.

  /* Polynomial routines.  These all set an error if the ring is not
     a polynomial ring.  OR, they will be moved to polyring.hpp  */

  virtual vec vec_diff(vec v, int rankFw, vec w, int use_coeff) const { M2_BAD_LEVEL }
  virtual int vec_in_subring(int n, const vec v) const { M2_BAD_LEVEL }
  virtual void vec_degree_of_var(int n, const vec v, int &lo, int &hi) const { M2_BAD_LEVEL }
  virtual vec vec_divide_by_var(int n, int d, const vec v) const { M2_BAD_LEVEL }
  virtual vec vec_divide_by_expvector(const int *exp, const vec v) const { M2_BAD_LEVEL }

  virtual ring_elem vec_content(vec f) const = 0;
  virtual vec vec_divide_by_given_content(vec f, ring_elem c) const = 0;
  virtual vec vec_divide_by_content(vec f) const = 0;
  virtual ring_elem vec_split_off_content(vec f, vec &result) const = 0;

  virtual vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     int deg,
                     M2_arrayint wts) const = 0;
  virtual vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     M2_arrayint wts) const = 0;
  virtual void vec_degree_weights(const FreeModule *F,
                          const vec f,
                          M2_arrayint wts,
                          int &lo,
                          int &hi) const = 0;
  virtual void vec_degree(const FreeModule *F, const vec f, int *d) const = 0;
  virtual bool vec_multi_degree(const FreeModule *F, const vec f, int *degf) const = 0;
};

template <typename Derived>
class RingBase : public virtual IRing,
                 public MutableEngineObject  // TODO(RZ): make sure that's OK
{
protected:
  const Derived* crtp() const { return static_cast<const Derived*>(this); }
  template <typename ReturnT>
  auto not_impl() const -> ReturnT
  {
    throw exc::engine_error("crtp method not implemented at this level");
  }

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

  M2_arrayint get_heft_vector() const { return m_heft_vector; }  // can be NULL

public:
  std::pair<bool, long> coerceToLongInteger(ring_elem a) const override { return {}; }
  ring_elem from_long(long n) const override { return crtp()->template not_impl<ring_elem>(); }
  ring_elem from_int(mpz_srcptr n) const override { return crtp()->template not_impl<ring_elem>(); }

  bool from_rational(const mpq_srcptr q, ring_elem &result) const override { return crtp()->template not_impl<bool>(); }
  bool from_BigReal(gmp_RR a, ring_elem &result) const override { result=from_long(0); return false; }
  bool from_BigComplex(gmp_CC z, ring_elem &result) const override { result=from_long(0); return false; }
  bool from_double(double a, ring_elem &result) const override { result=from_long(0); return false; }
  bool from_complex_double(double re, double im, ring_elem &result) const override { result=from_long(0); return false; }

  ring_elem var(int v) const override { return zero(); }

  ring_elem preferred_associate(ring_elem f) const override;
  bool lower_associate_divisor(ring_elem &f, ring_elem g) const override;

  void lower_content(ring_elem &c, ring_elem g) const override;
  ring_elem content(ring_elem f) const override;
  ring_elem content(ring_elem f, ring_elem g) const override;
  ring_elem divide_by_given_content(ring_elem f, ring_elem c) const override;
  ring_elem divide_by_content(ring_elem f) const override;
  ring_elem split_off_content(ring_elem f, ring_elem &result) const override;

  bool promote(const Ring *R, const ring_elem f, ring_elem &result) const override { return crtp()->template not_impl<bool>(); }
  bool lift(const Ring *R, const ring_elem f, ring_elem &result) const override { return crtp()->template not_impl<bool>(); }

  bool is_unit(const ring_elem f) const override { return crtp()->template not_impl<bool>(); }
  bool is_zero(const ring_elem f) const override { return crtp()->template not_impl<bool>(); }

  bool is_equal(const ring_elem f, const ring_elem g) const override { return crtp()->template not_impl<bool>(); }
  bool is_equal(const vecterm *a, const vecterm *b) const override;
  int compare_elems(const ring_elem f, const ring_elem g) const override { return crtp()->template not_impl<bool>(); }

  ring_elem copy(const ring_elem f) const override { throw exc::engine_error("crtp method not implemented on this level"); }
  void remove(ring_elem &f) const override { throw exc::engine_error("crtp method not implemented on this level"); }

  ring_elem negate(const ring_elem f) const override { throw exc::engine_error("crtp method not implemented on this level"); }
  ring_elem add(const ring_elem f, const ring_elem g) const override { throw exc::engine_error("crtp method not implemented on this level"); }
  ring_elem subtract(const ring_elem f, const ring_elem g) const override { throw exc::engine_error("crtp method not implemented on this level"); }
  ring_elem mult(const ring_elem f, const ring_elem g) const override { throw exc::engine_error("crtp method not implemented on this level"); }

  void negate_to(ring_elem &f) const override { f=crtp()->negate(f); }
  void add_to(ring_elem &f, ring_elem &g) const override { f=crtp()->add(f,g); }
  void subtract_to(ring_elem &f, ring_elem &g) const override { f=crtp()->subtract(f,g); }
  void mult_to(ring_elem &f, const ring_elem g) const override { f=crtp()->mult(f,g); }

  ring_elem power(const ring_elem f, int n) const override { return crtp()->impl_power(f,n); }
  ring_elem power(const ring_elem f, mpz_srcptr n) const override { return crtp()->impl_power(f,n); }
  ring_elem invert(const ring_elem f) const override { throw exc::engine_error("crtp method not implemented on this level"); }
  ring_elem divide(const ring_elem f, const ring_elem g) const override;

  ring_elem remainder(const ring_elem f, const ring_elem g) const override;
  ring_elem quotient(const ring_elem f, const ring_elem g) const override;
  ring_elem remainderAndQuotient(const ring_elem f, const ring_elem g, ring_elem &quot) const override;

  void syzygy(const ring_elem a, const ring_elem b, ring_elem &x, ring_elem &y) const override {}

  int index_of_var(const ring_elem a) const override;
  M2_arrayint support(const ring_elem a) const override;

  void monomial_divisor(const ring_elem a, int *exp) const override {}
  ring_elem diff(ring_elem a, ring_elem b, int use_coeff) const override { return crtp()->mult(a, b); }
  bool in_subring(int nslots, const ring_elem a) const override { return true; }
  void degree_of_var(int n, const ring_elem a, int &lo, int &hi) const override { hi = lo = 0; }
  ring_elem divide_by_var(int n, int d, const ring_elem a) const override;
  ring_elem divide_by_expvector(const int *exp, const ring_elem a) const override { return a; }

  ring_elem homogenize(const ring_elem f, int v, int deg, M2_arrayint wts) const override;
  ring_elem homogenize(const ring_elem f, int v, M2_arrayint wts) const override { return f; }

  // Routines expecting a grading.  The default implementation
  // is that the only degree is 0.
  bool is_homogeneous(const ring_elem f) const override { return true; }
  void degree(const ring_elem f, int *d) const override;
  bool multi_degree(const ring_elem f, int *d) const override;
  // returns true iff f is homogeneous
  void degree_weights(const ring_elem f,
                              M2_arrayint wts,
                              int &lo,
                              int &hi) const override { hi=lo=0; }


 protected:
  [[maybe_unused]] ring_elem impl_power(ring_elem f, int n) const;
  [[maybe_unused]] ring_elem impl_power(ring_elem f, mpz_srcptr n) const;

//  [[maybe_unused]] ring_elem impl_mult(ring_elem f, ring_elem g) const;

  /** @name Vector Methods *****************************************************
   *
   * These are methods for vector operations. Move to an adaptor or a decorator
   * class later.
   */
protected:
  vec new_vec() const { return new vecterm; }
  void remove_vec_node(vec n) const;

public:
  vec vec_zeroize_tiny(gmp_RR epsilon, const vec f) const;
  // If any real number appearing in f has larger absolute value than norm,
  // replace norm.
  // Default for rings not over RRR or CCC is to do nothing.
  // Default for rings not over RRR or CCC is to do nothing.
  void vec_increase_maxnorm(gmp_RRmutable norm, const vec f) const;


  void vec_sort(vecterm *&f) const;

  int compare_vecs(vec v, vec w) const;

  vec e_sub_i(int r) const;
  vec make_vec(int r, ring_elem a) const;
  vec make_vec_from_array(int len, Nterm **array) const;  /// takes ownership of the Nterm's!!

  vec copy_vec(const vecterm *v) const;
  void remove_vec(vec v) const;

  void add_vec_to(vec &v, vec &w) const;       /// v <- v+w, w is set to 0.

  bool get_entry(const vecterm *v, int r, ring_elem &result) const;
  ring_elem get_entry(vec v, int r) const;
  vec sub_vector(const vecterm *v, M2_arrayint r) const;
  int n_nonzero_terms(const vecterm *v) const;
  void vec_text_out(buffer &o,
                    const vecterm *v,
                    bool p_one = true,
                    bool p_plus = false,
                    bool p_parens = false) const;


  virtual vec vec_lead_term(int nparts, const FreeModule *F, vec v) const;

  vec negate_vec(vec v) const;
  vec add_vec(vec v, vec w) const;
  vec subtract_vec(vec v, vec w) const;
  vec mult_vec(int n, vec v) const;
  vec mult_vec(const ring_elem f, const vec w) const;
  vec rightmult_vec(const vec w, const ring_elem f) const;

  void set_entry(vec &v, int i, ring_elem r) const;
  void mult_vec_to(vec &v,
                   const ring_elem r,
                   bool opposite_mult) const;  // multiplies v <- r * v or v * r
  void mult_row(vec &v, const ring_elem r, int i, bool opposite_mult) const;
  void negate_vec_to(vec &v) const;            // v <- -v.
  void subtract_vec_to(vec &v, vec &w) const;  // v <- v-w, w is set to 0.

  vec component_shift(int n, vec v) const;

  vec tensor_shift(int n, int m, vec v) const;

  void divide_vec_to(vec &v, const ring_elem a) const;
  void divide_row(vec &v, int r, const ring_elem a) const;
  ring_elem dot_product(const vecterm *v, const vecterm *w) const;

  ring_elem vec_content(vec f) const override;
  vec vec_divide_by_given_content(vec f, ring_elem c) const override;
  vec vec_divide_by_content(vec f) const override;
  ring_elem vec_split_off_content(vec f, vec &result) const override;

  vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     int deg,
                     M2_arrayint wts) const override;
  vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     M2_arrayint wts) const override;
  void vec_degree_weights(const FreeModule *F,
                          const vec f,
                          M2_arrayint wts,
                          int &lo,
                          int &hi) const override;
  void vec_degree(const FreeModule *F, const vec f, int *d) const override;
  bool vec_multi_degree(const FreeModule *F, const vec f, int *degf) const override;
};

#include "Ring.hpp"  // safe but confuses CLion
#include "RingVec.hpp"
