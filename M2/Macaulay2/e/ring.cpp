// Copyright 1995 Michael E. Stillman

#include "ring.hpp"
#include "aring-RRR.hpp"
#include "aring-CCC.hpp"
#include "monoid.hpp"
#include "poly.hpp"
#include "base/bugfix.h"

RingZZ *makeIntegerRing() { return new RingZZ; }

const Monoid *Ring::impl_degree_monoid() const { return get_degree_ring()->getMonoid(); }

void Ring::initialize_ring(long P0,
                           const PolynomialRing *DR,
                           const M2_arrayint heft_vec)
{
  // Remember: if this is a poly ring, the ring is K[M].
  // If this is a basic routine, K = this, M = trivial monoid.
  // If this is a frac field, K = R, M = trivial monoid.
  Base::m_char = P0;
  if (DR == 0)
    Base::m_degree_ring = PolyRing::get_trivial_poly_ring();
  else
    Base::m_degree_ring = DR;
  Base::m_heft_vector = heft_vec;

  Base::m_zeroV = ZERO_RINGELEM;
  Base::m_oneV = ZERO_RINGELEM;
  Base::m_minus_oneV = ZERO_RINGELEM;

  Base::m_non_unit = ZERO_RINGELEM;
  Base::m_isfield = 0;
}


///////////////////////////////////
// SumCollector: default version //
///////////////////////////////////
class SumCollectorDefault : public SumCollector
{
  const Ring *R;
  ring_elem result;

 public:
  SumCollectorDefault(const Ring *R0) : R(R0), result(R->zero()) {}
  virtual ~SumCollectorDefault() {}
  virtual void add(ring_elem f) { R->add_to(result, f); }
  virtual ring_elem getValue()
  {
    ring_elem val = result;
    result = R->zero();
    return val;
  }
};

SumCollector *Ring::make_SumCollector() const
{
  return new SumCollectorDefault(this);
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
