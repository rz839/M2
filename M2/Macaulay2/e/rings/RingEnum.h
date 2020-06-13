//
// Created by Radoslav Zlatev on 6/11/20.
//

#pragma once

namespace M2 {

/**
 * Ring traits.
 */
enum class RingTrait : uint8_t
{
  BASIC_RING,
  FINITE_PRIME_FIELD,
  GALOIS_FIELD,
  FRACTION_FIELD,
  FRACTION_POLYNRING,
  POLYRING,

  ZZ,
  QQ,
  RR,
  CC,
  WEYL_ALGEBRA,
  SOLVABLE_ALGEBRA,

  COMMUTATIVE,
  SKEWCOMMUTATIVE,
  GRADED,
  QUOTIENT,
  FIELD_DECLARED,       /// declared to be a field
  FIELD_FOUND,          /// the class may or may not be field, but this instance was found to be one

  ASSOC_DIVISORS,
  GCD,
};

/**
 * RingTypeId identifies the ring class, both generically and concretely.
 * That is, a ring could be both RingZZ and RingZZFlint.
 * This is to be used across promotion and lifting, as well as in the runtime.
 */
enum class RingTypeId : uint8_t {
  // generic rings
  RING,
  SCHURRING,  // SchurRing2

  // concrete rings
  FLINT_ZZ,
};

}  // namespace M2
