#include "util.hpp"
#include "polyring.hpp"
#include "ring.hpp"
#include "monoid.hpp"
#include "qring.hpp"
#include "polyquotient.hpp"
#include "matrix.hpp"
#include "matrix-con.hpp"
#include "geopoly.hpp"

PolynomialRing::~PolynomialRing() {}
void PolynomialRing::setQuotientInfo(QRingInfo *qinfo0)
{
  qinfo_ = qinfo0;
  const PolyRing *numerR = getNumeratorRing();  // might be 'this'

  for (int i = 0; i < n_quotients(); i++)
    {
      if (!numerR->is_homogeneous(quotient_element(i)))
        {
          setIsGraded(false);
          break;
        }
    }

  overZZ_ = (coeff_type_ == Ring::COEFF_ZZ);
}

void PolynomialRing::initialize_PolynomialRing(const Ring *K,
                                               const Monoid *M,
                                               const PolyRing *numeratorR,
                                               const PolynomialRing *ambientR,
                                               const Ring *denomR)
{
  nvars_ = M->n_vars();
  K_ = K;
  M_ = M;
  numerR_ = numeratorR;
  ambientR_ = ambientR;
  denomR_ = denomR;

  exp_size = EXPONENT_BYTE_SIZE(nvars_);

  if (K->is_QQ() || (K == globalZZ && denomR != 0))
    coeff_type_ = Ring::COEFF_QQ;
  else if (K == globalZZ && denomR == 0)
    coeff_type_ = Ring::COEFF_ZZ;
  else
    coeff_type_ = Ring::COEFF_BASIC;

  is_weyl_ = false;
  is_solvable_ = false;
  is_skew_ = false;
  overZZ_ = false;
  qinfo_ = new QRingInfo;
  is_ZZ_quotient_ = false;
  ZZ_quotient_value_ = ZERO_RINGELEM;

  if (numeratorR != this)
    {
      // We must set the non-commutative settings ourselves at this time
      if (numeratorR->cast_to_WeylAlgebra() != 0)
//      if (dynamic_cast<const WeylAlgebra*>(numeratorR))
        is_weyl_ = true;
      else if (numeratorR->cast_to_SolvableAlgebra() != 0)
//      else if (dynamic_cast<const SolvableAlgebra*>(numeratorR))
        is_solvable_ = true;
      else if (numeratorR->is_skew_commutative())
        {
          is_skew_ = true;
          skew_ = numeratorR->getSkewInfo();
        }
    }

  poly_size_ = 0;  // The callee needs to set this later
  gb_ring_ = 0;    // The callee needs to set this later

  // Also: callee should call setIsGraded, and set oneV, minus_oneV, zeroV
}

PolynomialRing *PolynomialRing::create_quotient(const PolynomialRing *R,
                                                VECTOR(Nterm *) & elems)
// Grabs 'elems'.  Each element of 'elems' should be in the ring R.
// They should also form a GB.
{
  // Here are the cases:
  // (1) R is a polynomial ring over a basic field
  // (2) R is a polynomial ring over ZZ
  // (3) R is a polynomial ring over QQ

  // case (1), (2): PolyRingQuotient
  // case (3): PolyQQ

  PolynomialRing *result = NULL;
  Ring::CoefficientType coeff_type = R->coefficient_type();

  QRingInfo *qrinfo = NULL;
  switch (coeff_type)
    {
      case COEFF_BASIC:
        qrinfo = new QRingInfo_field_basic(R->getNumeratorRing(), elems);
        result = new PolyRingQuotient;
        break;
      case COEFF_QQ:
        qrinfo = new QRingInfo_field_QQ(R->getNumeratorRing(), elems);
        result = new PolyRingQuotient;
        break;
      case COEFF_ZZ:
        QRingInfo_ZZ *qrinfoZZ = new QRingInfo_ZZ(R->getNumeratorRing(), elems);
        qrinfo = qrinfoZZ;
        result = new PolyRingQuotient;
        result->is_ZZ_quotient_ = qrinfoZZ->is_ZZ_quotient();
        result->ZZ_quotient_value_ = qrinfoZZ->ZZ_quotient_value();
        break;
    }

  result->initialize_ring(
      R->characteristic(), R->get_degree_ring(), R->get_heft_vector());

  result->initialize_PolynomialRing(R->getCoefficients(),
                                    R->getMonoid(),
                                    R->getNumeratorRing(),
                                    R->getAmbientRing(),
                                    R->getDenominatorRing());

  result->gb_ring_ = R->get_gb_ring();
  result->setQuotientInfo(qrinfo);  // Also sets graded-ness

  result->Base::m_zeroV = result->from_long(0);
  result->Base::m_oneV = result->from_long(1);
  result->Base::m_minus_oneV = result->from_long(-1);

  return result;
}

PolynomialRing *PolynomialRing::create_quotient(const PolynomialRing *R,
                                                const Matrix *M)
{
  if (M->get_ring() != R)
    {
      ERROR("quotient elements not in the expected polynomial ring");
      return 0;
    }
  VECTOR(Nterm *) elems;

  for (int i = 0; i < M->n_cols(); i++)
    {
      Nterm *f = R->numerator(M->elem(0, i));
      elems.push_back(f);
    }

  for (int i = 0; i < R->n_quotients(); i++)
    elems.push_back(R->quotient_element(i));

  return create_quotient(R->getAmbientRing(), elems);
}

PolynomialRing *PolynomialRing::create_quotient(const PolynomialRing *R,
                                                const PolynomialRing *B)
// R should be an ambient poly ring
// B should have: ambient of B is the logical coeff ring of R
//   i.e. R = A[x], B = A/I
// return A[x]/I.
{
  VECTOR(Nterm *) elems;

  for (int i = 0; i < B->n_quotients(); i++)
    {
      ring_elem f;
      R->promote(B->getNumeratorRing(), B->quotient_element(i), f);
      elems.push_back(f);
    }
  return create_quotient(R, elems);
}

Matrix *PolynomialRing::getPresentation() const
{
  const PolynomialRing *R = getAmbientRing();

  MatrixConstructor mat(R->make_FreeModule(1), 0);
  for (int i = 0; i < n_quotients(); i++)
    // NEED: to make this into a fraction, if R has fractions.
    mat.append(R->make_vec(0, quotient_element(i)));
  return mat.to_matrix();
}

class SumCollectorPolyHeap : public SumCollector
{
  polyheap H;

 public:
  SumCollectorPolyHeap(const PolynomialRing *R0) : H(R0) {}
  ~SumCollectorPolyHeap() {}
  virtual void add(ring_elem f) { H.add(f); }
  virtual ring_elem getValue() { return H.value(); }
};

SumCollector *PolynomialRing::make_SumCollector() const
{
  return new SumCollectorPolyHeap(this);
}

unsigned int PolynomialRing::computeHashValue(const ring_elem a) const
{
  unsigned int hash = 0;
  unsigned int seed1 = 103;
  unsigned int seed2 = 347654;
  for (const Nterm *t = a.poly_val; t != 0; t = t->next)
    {
      unsigned int hash1 = getCoefficientRing()->computeHashValue(t->coeff);
      unsigned int hash2 = getMonoid()->computeHashValue(t->monom);
      hash += seed1 * hash1 + seed2 * hash2;
      seed1 += 463633;
      seed2 += 7858565;
    }
  return hash;
}

vec PolynomialRing::vec_diff(vec v, int rankFw, vec w, int use_coeff) const
// rankFw is the rank of the free module corresponding to w.
{
  vec result = NULL;
  for (; v != NULL; v = v->next)
    for (vecterm *p = w; p != NULL; p = p->next)
      {
        ring_elem a = diff(v->coeff, p->coeff, use_coeff);
        if (is_zero(a))
          {
            remove(a);
            continue;
          }
        vecterm *t = new_vec();
        t->comp = rankFw * v->comp + p->comp;
        t->coeff = a;
        t->next = result;
        result = t;
      }
  vec_sort(result);
  return result;
}

int PolynomialRing::vec_in_subring(int nslots, const vec v) const
{
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0 || v == NULL) return true;
  const Monoid *M = PR->getMonoid();
  for (vec w = v; w != NULL; w = w->next)
    if (!M->in_subring(nslots, PR->lead_flat_monomial(w->coeff))) return false;
  return true;
}

void PolynomialRing::vec_degree_of_var(int n, const vec v, int &lo, int &hi) const
{
  if (v == NULL)
    {
      ERROR("attempting to find degree of zero vector");
      return;
    }
  degree_of_var(n, v->coeff, lo, hi);
  for (vec w = v->next; w != 0; w = w->next)
    {
      int lo1, hi1;
      degree_of_var(n, w->coeff, lo1, hi1);
      if (lo1 < lo) lo = lo1;
      if (hi1 > hi) hi = hi1;
    }
}

vec PolynomialRing::vec_divide_by_var(int n, int d, const vec v) const
{
  vecterm head;
  vecterm *result = &head;
  for (vec w = v; w != 0; w = w->next)
    {
      ring_elem a = divide_by_var(n, d, w->coeff);
      if (!is_zero(a))
        {
          vec t = make_vec(w->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = 0;
  return head.next;
}

vec PolynomialRing::vec_divide_by_expvector(const int *exp, const vec v) const
{
  vecterm head;
  vecterm *result = &head;
  for (vec w = v; w != 0; w = w->next)
    {
      ring_elem a = divide_by_expvector(exp, w->coeff);
      if (!is_zero(a))
        {
          vec t = make_vec(w->comp, a);
          result->next = t;
          result = t;
        }
    }
  result->next = 0;
  return head.next;
}

bool PolynomialRing::vec_is_scalar_multiple(vec f, vec g) const
// is df = cg, some scalars c,d?
// These scalars are over the very bottom base field/ZZ.
{
  if (f == NULL) return true;
  if (g == NULL) return true;
  const PolynomialRing *PR = cast_to_PolynomialRing();
  if (PR == 0) return true;
  const PolyRing *PR1 = PR->getNumeratorRing();
#ifdef DEVELOPMENT
#warning "use numerator only"
#endif
  if (f->comp != g->comp) return false;
  Nterm *f1 = f->coeff;
  Nterm *g1 = g->coeff;
  ring_elem c = f1->coeff;
  ring_elem d = g1->coeff;
  vec p, q;
  for (p = f, q = g; p != NULL && q != NULL; p = p->next, q = q->next)
    {
      if (p->comp != q->comp) return 0;
      if (!M2::bugfix::check_nterm_multiples(PR1, p->coeff, q->coeff, c, d)) return false;
    }
  return !p && !q;
}


// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
