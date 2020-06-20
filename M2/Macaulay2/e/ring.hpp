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
  // this is poorly written, as so much else;
  // you only really call this from other derived (general) rings further down
  // the inheritance, so move it to RingBase<D>
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

  virtual MutableMatrix *makeMutableMatrix(size_t nrows, size_t ncols, bool dense) const { return nullptr; }
  virtual FreeModule *make_FreeModule() const;
  virtual FreeModule *make_Schreyer_FreeModule() const;
  virtual FreeModule *make_FreeModule(int n) const;
  virtual SumCollector *make_SumCollector() const;


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








 public:

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
