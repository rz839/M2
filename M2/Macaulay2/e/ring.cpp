// Copyright 1995 Michael E. Stillman

#include "ring.hpp"
#include "aring-RRR.hpp"
#include "aring-CCC.hpp"
#include "monoid.hpp"
#include "poly.hpp"

#include "freemod.hpp"
#include "coeffrings.hpp"

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

FreeModule *Ring::make_FreeModule() const
{
  return new FreeModule(this, 0, false);
}

FreeModule *Ring::make_Schreyer_FreeModule() const
{
  return new FreeModule(this, 0, true);
}

FreeModule *Ring::make_FreeModule(int n) const
{
  return new FreeModule(this, n, false);
}

ring_elem Ring::power(const ring_elem gg, mpz_srcptr m) const
{
  ring_elem ff = gg;
  int cmp = mpz_sgn(m);
  if (cmp == 0) return one();
  mpz_t n;
  mpz_init_set(n, m);
  if (cmp < 0)
    {
      mpz_neg(n, n);
      ff = invert(ff);
      if (is_zero(ff))
        {
          ERROR(
              "either element not invertible, or no method available to "
              "compute its inverse");
          return ff;
        }
    }
  ring_elem prod = from_long(1);
  ring_elem base = copy(ff);
  ring_elem tmp;

  for (;;)
    {
      if (RingZZ::mod_ui(n, 2) == 1)
        {
          tmp = mult(prod, base);
          prod = tmp;
        }
      mpz_tdiv_q_2exp(n, n, 1);
      if (mpz_sgn(n) == 0)
        {
          mpz_clear(n);
          return prod;
        }
      else
        {
          tmp = mult(base, base);
          base = tmp;
        }
    }
}

//ring_elem Ring::power(const ring_elem gg, int n) const
//{
//  ring_elem ff = gg;
//  if (n == 0) return one();
//  if (n < 0)
//    {
//      n = -n;
//      ff = invert(ff);
//      if (is_zero(ff))
//        {
//          ERROR("negative power of noninvertible element requested");
//          return ff;
//        }
//    }
//
//  // The exponent 'n' should be > 0 here.
//  ring_elem prod = from_long(1);
//  ring_elem base = copy(ff);
//  ring_elem tmp;
//
//  for (;;)
//    {
//      if ((n % 2) != 0)
//        {
//          tmp = mult(prod, base);
//          prod = tmp;
//        }
//      n >>= 1;
//      if (n == 0)
//        {
//          return prod;
//        }
//      else
//        {
//          tmp = mult(base, base);
//          base = tmp;
//        }
//    }
//}

ring_elem Ring::remainder(const ring_elem f, const ring_elem g) const
{
  if (is_zero(g)) return f;
  return zero();
}

ring_elem Ring::quotient(const ring_elem f, const ring_elem g) const
{
  if (is_zero(g)) return g;
  return divide(f, g);
}

ring_elem Ring::remainderAndQuotient(const ring_elem f,
                                     const ring_elem g,
                                     ring_elem &quot) const
{
  if (is_zero(g))
    {
      quot = g;  // zero
      return f;
    }
  quot = divide(f, g);
  return zero();
}

void Ring::monomial_divisor(const ring_elem a, int *exp) const
{
  // Do nothing
}

ring_elem Ring::diff(ring_elem a, ring_elem b, int use_coeff) const
{
  return mult(a, b);
}

bool Ring::in_subring(int nslots, const ring_elem a) const { return true; }
void Ring::degree_of_var(int n, const ring_elem a, int &lo, int &hi) const
{
  lo = 0;
  hi = 0;
}

ring_elem Ring::divide_by_var(int n, int d, const ring_elem a) const
{
  if (d == 0) return a;
  return from_long(0);
}

ring_elem Ring::divide_by_expvector(const int *exp, const ring_elem a) const
{
  return a;
}

ring_elem Ring::homogenize(const ring_elem f, int, int deg, M2_arrayint) const
{
  if (deg != 0) ERROR("homogenize: no homogenization exists");
  return f;
}

ring_elem Ring::homogenize(const ring_elem f, int, M2_arrayint) const
{
  return f;
}

bool Ring::is_homogeneous(const ring_elem) const { return true; }
void Ring::degree(const ring_elem, int *d) const { degree_monoid()->one(d); }
bool Ring::multi_degree(const ring_elem f, int *d) const
// returns true iff f is homogeneous
{
  degree_monoid()->one(d);
  return true;
}

void Ring::degree_weights(const ring_elem, M2_arrayint, int &lo, int &hi) const
{
  lo = hi = 0;
}
int Ring::index_of_var(const ring_elem a) const { return -1; }
M2_arrayint Ring::support(const ring_elem a) const
{
  M2_arrayint result = M2_makearrayint(0);
  return result;
}

// These next three routines are only overridden by RRR,CCC,polynomial rings,
// and quotient rings
unsigned long Ring::get_precision() const { return 0; }
ring_elem Ring::zeroize_tiny(gmp_RR epsilon, const ring_elem f) const
// Default is to return f itself.
{
  return f;
}

void Ring::increase_maxnorm(gmp_RRmutable norm, const ring_elem f) const
// If any real number appearing in f has larger absolute value than norm,
// replace norm.
{
  // Default for rings not over RRR or CCC is to do nothing.
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
