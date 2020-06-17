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




  virtual ring_elem random() const
  {
    throw exc::engine_error("random scalar elements for this ring are not implemented");
  }


  // Polynomial routines
  // The default implementation is for non-polynomial rings


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
