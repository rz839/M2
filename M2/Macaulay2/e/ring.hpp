// Copyright 1995 Michael E. Stillman

#ifndef _ring_hh_
#define _ring_hh_

#include "hash.hpp"
#include "error.h"
#include "ringelem.hpp"
#include "monoid.hpp"
#include "aring.hpp"
#include "exceptions.hpp"
///// Ring Hierarchy ///////////////////////////////////

/***************
 * CRTP includes
 ***************/
#include "rings/Ring.hpp"

class RingZZ;
class RRR;
class CCC;
class Tower;
class PolynomialRing;
class PolyRing;
class PolyRingFlat;
class SkewPolynomialRing;
class SchurRing;
class SchurRing2;
class SchurSnRing;
class WeylAlgebra;
class SolvableAlgebra;

class FreeModule;
class RingMap;

class gbvectorHeap;
class gbvector;
class buffer;

class SumCollector;

class ARing;
class MutableMatrix;

class CoefficientRingR;

/**
    @ingroup rings

    @brief xxx
    xxx
    xxx
*/
class Ring : public RingBase<Ring>
{
public:
  using RingBase<Ring>::power;
  using Base = RingBase<Ring>;

protected:
  void initialize_ring(long charac,
                       const PolynomialRing *DR = nullptr,
                       const M2_arrayint heft_vec = nullptr);
  Ring() = default;

public:
  virtual ~Ring() = default;
  const Monoid *impl_degree_monoid() const;

  // ---------------------------------------------------------------------------
  //   Other Methods (used for specific rings)
  // ---------------------------------------------------------------------------

  virtual MutableMatrix *makeMutableMatrix(size_t nrows,
                                           size_t ncols,
                                           bool dense) const
  {
    return 0;
  }

  virtual FreeModule *make_FreeModule() const;
  virtual FreeModule *make_Schreyer_FreeModule() const;
  virtual FreeModule *make_FreeModule(int n) const;

  virtual SumCollector *make_SumCollector() const;

  // ---------------------------------------------------------------------------
  //   Ring Arithmetic and Conversion
  // ---------------------------------------------------------------------------
  virtual std::pair<bool, long> coerceToLongInteger(ring_elem a) const;

  virtual ring_elem from_long(long n) const = 0;
  virtual ring_elem from_int(mpz_srcptr n) const = 0;

  // from_rational: if the rational q cannot be placed into this ring, false is
  // returned, and result is not touched.
  virtual bool from_rational(const mpq_srcptr q, ring_elem &result) const = 0;

  // The default version calls from_long(0) and returns false.
  virtual bool from_BigReal(gmp_RR a, ring_elem &result) const;
  // The default version calls from_long(0) and returns false.
  virtual bool from_BigComplex(gmp_CC z, ring_elem &result) const;
  // Returns false if this ring cannot coerce a double to an element in this
  // ring
  virtual bool from_double(double a, ring_elem &result) const;
  // Returns false if this ring cannot coerce a complex double (re+im*ii) to an
  // element in this ring
  virtual bool from_complex_double(double re,
                                   double im,
                                   ring_elem &result) const;

  virtual ring_elem var(int v) const;

  virtual ring_elem preferred_associate(ring_elem f) const;
  // Returns an invertible element c of the same ring such that c*f is the
  // preferred associate of the element f.
  // WARNING: The default implementation is for a field.

  virtual bool lower_associate_divisor(ring_elem &f, ring_elem g) const;
  // Replaces f with the unit c such that (fx+g)//c is the preferred associate
  //   of fx+g, in the ring A[x], where A is 'this'.
  // Returns false if f will never be changed after this
  // (This happens over ZZ if f is non-zero (therefore 1 or -1, over a finite
  // filed if f != 0,
  // but over QQ will never happen)
  // WARNING: The default implementation is for a field.

  virtual bool promote(const Ring *R,
                       const ring_elem f,
                       ring_elem &result) const = 0;
  virtual bool lift(const Ring *R,
                    const ring_elem f,
                    ring_elem &result) const = 0;

  virtual bool is_unit(const ring_elem f) const = 0;
  virtual bool is_zero(const ring_elem f) const = 0;
  virtual bool is_equal(const ring_elem f, const ring_elem g) const = 0;

  virtual int compare_elems(const ring_elem f, const ring_elem g) const = 0;
  // Returns -1 if f < g, 0 if f and g are either equal, or considered equal,
  // and 1 if f > g.
  // The specific definition is ring dependent, but for numeric rings it
  // defaults to
  // the standard order.

  virtual ring_elem copy(const ring_elem f) const = 0;
  virtual void remove(ring_elem &f) const = 0;

  void negate_to(ring_elem &f) const;
  void add_to(ring_elem &f, ring_elem &g) const;
  void subtract_to(ring_elem &f, ring_elem &g) const;
  void mult_to(ring_elem &f, const ring_elem g) const;

  virtual ring_elem negate(const ring_elem f) const = 0;
  virtual ring_elem add(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem subtract(const ring_elem f, const ring_elem g) const = 0;
  virtual ring_elem mult(const ring_elem f, const ring_elem g) const = 0;

  virtual ring_elem power(const ring_elem f, mpz_srcptr n) const;
//  virtual ring_elem power(const ring_elem f, int n) const;  // TODO(RZ): remove - crtp'ed
  // These two power routines can be used for n >= 0.

//  virtual ring_elem invert(const ring_elem f) const = 0;  // TODO(RZ): remove - crtp'ed
  virtual ring_elem divide(const ring_elem f, const ring_elem g) const = 0;

  virtual ring_elem remainder(const ring_elem f, const ring_elem g) const;
  virtual ring_elem quotient(const ring_elem f, const ring_elem g) const;
  virtual ring_elem remainderAndQuotient(const ring_elem f,
                                         const ring_elem g,
                                         ring_elem &quot) const;
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

  virtual void syzygy(const ring_elem a,
                      const ring_elem b,
                      ring_elem &x,
                      ring_elem &y) const = 0;
  // Constructs elements x and y in the ring s.t. ax + by = 0.  This syzygy is
  // chosen as simply as possible.  For example, over QQ, x is chosen
  // to be positive.  The routine must handle the case when a=0, but can
  // ignore the case when b=0... (Really?)

  virtual ring_elem random() const;

  virtual void elem_text_out(buffer &o,
                             const ring_elem f,
                             bool p_one = true,
                             bool p_plus = false,
                             bool p_parens = false) const = 0;

  virtual ring_elem eval(const RingMap *map,
                         const ring_elem f,
                         int first_var) const = 0;

  // Polynomial routines
  // The default implementation is for non-polynomial rings
  virtual int index_of_var(const ring_elem a) const;
  virtual M2_arrayint support(const ring_elem a) const;

  virtual void monomial_divisor(const ring_elem a, int *exp) const;
  virtual ring_elem diff(ring_elem a, ring_elem b, int use_coeff) const;
  virtual bool in_subring(int nslots, const ring_elem a) const;
  virtual void degree_of_var(int n, const ring_elem a, int &lo, int &hi) const;
  virtual ring_elem divide_by_var(int n, int d, const ring_elem a) const;
  virtual ring_elem divide_by_expvector(const int *exp,
                                        const ring_elem a) const;

  virtual ring_elem homogenize(const ring_elem f,
                               int v,
                               int deg,
                               M2_arrayint wts) const;
  virtual ring_elem homogenize(const ring_elem f, int v, M2_arrayint wts) const;

  // Routines expecting a grading.  The default implementation
  // is that the only degree is 0.
  virtual bool is_homogeneous(const ring_elem f) const;
  virtual void degree(const ring_elem f, int *d) const;
  virtual bool multi_degree(const ring_elem f, int *d) const;
  // returns true iff f is homogeneous
  virtual void degree_weights(const ring_elem f,
                              M2_arrayint wts,
                              int &lo,
                              int &hi) const;

  // antipode: for non skew commuting poly rings, this is the identity.
  // Otherwise, this changes the signs of the monomials, implementing the
  // (anti)-isomorphism
  // of the ring with its opposite ring.
  virtual ring_elem antipode(ring_elem f) const { return f; }
  //////////////////////////////////////////
  // Cleaning real and complex numbers /////
  //////////////////////////////////////////
  virtual unsigned long get_precision()
      const;  // if the ring is not over RRR or CCC, returns 0.
  virtual ring_elem zeroize_tiny(gmp_RR epsilon, const ring_elem f) const;
  // Default is to return f itself.
  virtual void increase_maxnorm(gmp_RRmutable norm, const ring_elem f) const;
  // If any real number appearing in f has larger absolute value than norm,
  // replace norm.
  // Default for rings not over RRR or CCC is to do nothing.
  vec vec_zeroize_tiny(gmp_RR epsilon, const vec f) const;
  // Default is to return f itself.
  void vec_increase_maxnorm(gmp_RRmutable norm, const vec f) const;
  // If any real number appearing in f has larger absolute value than norm,
  // replace norm.
  // Default for rings not over RRR or CCC is to do nothing.

  //////////////////////////////////////////
  /// vector operations ////////////////////
  //////////////////////////////////////////
  // These routines all act on linked lists
  // of vecterm's, sorted by descending component.
  // We always assume that ringelem's are immutable:
  // The same value might be shared in several vecterms.
  //
  // These routines are implemented in ring-vec.cpp
  //////////////////////////////////////////
 protected:
//  vec new_vec() const;
  void remove_vec_node(vec n) const;

 public:
  void vec_sort(vecterm *&f) const;

  int compare_vecs(vec v, vec w) const;

  vec e_sub_i(int r) const;
  vec make_vec(int r, ring_elem a) const;
  vec make_vec_from_array(int len, Nterm **array)
      const;  // takes ownership of the Nterm's!!

  vec copy_vec(const vecterm *v) const;
  void remove_vec(vec v) const;

  bool is_equal(const vecterm *a, const vecterm *b) const;
  bool get_entry(const vecterm *v, int r, ring_elem &result) const;
  ring_elem get_entry(vec v, int r) const;
  vec sub_vector(const vecterm *v, M2_arrayint r) const;
  int n_nonzero_terms(const vecterm *v) const;
  void vec_text_out(buffer &o,
                    const vecterm *v,
                    bool p_one = true,
                    bool p_plus = false,
                    bool p_parens = false) const;
  vec vec_eval(const RingMap *map, const FreeModule *F, const vec v) const;

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
  void add_vec_to(vec &v, vec &w) const;       // v <- v+w, w is set to 0.
  void subtract_vec_to(vec &v, vec &w) const;  // v <- v-w, w is set to 0.

  vec mult_vec_matrix(const Matrix *m, vec v, bool opposite_mult) const;

  vec component_shift(int n, vec v) const;

  vec tensor_shift(int n, int m, vec v) const;

  vec tensor(const FreeModule *F, vec v, const FreeModule *G, vec w) const;

  void divide_vec_to(vec &v, const ring_elem a) const;
  void divide_row(vec &v, int r, const ring_elem a) const;
  ring_elem dot_product(const vecterm *v, const vecterm *w) const;

  /* Polynomial routines.  These all set an error if the ring is not
     a polynomial ring.  OR, they will be moved to polyring.hpp  */
  vec vec_diff(vec v, int rankFw, vec w, int use_coeff) const;
  int vec_in_subring(int n, const vec v) const;
  void vec_degree_of_var(int n, const vec v, int &lo, int &hi) const;
  vec vec_divide_by_var(int n, int d, const vec v) const;
  vec vec_divide_by_expvector(const int *exp, const vec v) const;

  // Some divisibility routines
  bool vec_is_scalar_multiple(vec f, vec g)
      const;  // is cf = dg, some scalars c,d? (not both zero).
  vec vec_remove_monomial_factors(vec f, bool make_squarefree_only) const;

  bool vec_multi_degree(const FreeModule *F, const vec f, int *degf) const;
  // returns true iff f is homogeneous

  void vec_degree(const FreeModule *F, const vec f, int *d) const;
  void vec_degree_weights(const FreeModule *F,
                          const vec f,
                          M2_arrayint wts,
                          int &lo,
                          int &hi) const;
  bool vec_is_homogeneous(const FreeModule *F, const vec f) const;
  vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     int deg,
                     M2_arrayint wts) const;
  vec vec_homogenize(const FreeModule *F,
                     const vec f,
                     int v,
                     M2_arrayint wts) const;

  // content of vectors and ring elements, default implementation is for basic
  // fields
  virtual void lower_content(ring_elem &c, ring_elem g)
      const;  // c is a content elem, g is in ring
  virtual ring_elem content(ring_elem f) const;
  virtual ring_elem content(ring_elem f, ring_elem g) const;
  virtual ring_elem divide_by_given_content(ring_elem f, ring_elem c) const;

  ring_elem divide_by_content(ring_elem f) const;
  ring_elem split_off_content(ring_elem f, ring_elem &result) const;
  ring_elem vec_content(vec f) const;
  vec vec_divide_by_given_content(vec f, ring_elem c) const;
  vec vec_divide_by_content(vec f) const;
  ring_elem vec_split_off_content(vec f, vec &result) const;
};

class SumCollector : public our_new_delete
{
 public:
  SumCollector() {}
  virtual ~SumCollector() {}
  virtual void add(ring_elem f) = 0;
  virtual ring_elem getValue() = 0;
};

#define ZERO_RINGELEM (ring_elem(static_cast<Nterm *>(0)))

#include "ZZ.hpp"
extern RingZZ *globalZZ;
extern RingZZ *makeIntegerRing();

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e ring.o "
// indent-tabs-mode: nil
// End:
