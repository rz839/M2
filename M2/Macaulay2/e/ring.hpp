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
  vec tensor(const FreeModule *F, vec v, const FreeModule *G, vec w) const;

  /* Polynomial routines.  These all set an error if the ring is not
     a polynomial ring.  OR, they will be moved to polyring.hpp  */

  // Some divisibility routines
  bool vec_is_scalar_multiple(vec f, vec g)
      const;  // is cf = dg, some scalars c,d? (not both zero).
  vec vec_remove_monomial_factors(vec f, bool make_squarefree_only) const;

  bool vec_multi_degree(const FreeModule *F, const vec f, int *degf) const;
  // returns true iff f is homogeneous

  void vec_degree(const FreeModule *F, const vec f, int *d) const;
  bool vec_is_homogeneous(const FreeModule *F, const vec f) const;
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
