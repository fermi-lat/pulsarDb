/** \file AbsoluteTime
    \brief Declaration of AbsoluteTime class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_AbsoluteTime_h
#define pulsarDb_AbsoluteTime_h

#include <iosfwd>

#include "pulsarDb/Duration.h"

namespace pulsarDb {

  class CanonicalTime;

  /** \class AbsoluteTime
      \brief An absolute moment in time.
  */
  class AbsoluteTime {
    public:
      virtual ~AbsoluteTime() {}

      /** \brief Subtract the given time from this one, to determine the duration of the interval between them.
          \param t The other time.
      */
      virtual Duration operator -(const AbsoluteTime & t) const = 0;

      /** \brief Add the given time increment to this one.
          \param d The duration (elapsed time, relative time) being added.
      */
      virtual AbsoluteTime & operator +=(const Duration & d) = 0;

      /** \brief Convert this time to a CanonicalTime representation of it.
          \param t The target time system time.
      */
      virtual void to(CanonicalTime & t) const = 0;

      virtual void write(std::ostream & os) const = 0;

      virtual AbsoluteTime * clone() const = 0;

      bool operator ==(const AbsoluteTime & t) const { return *this - t == Duration(0., UnitSec); }
      bool operator !=(const AbsoluteTime & t) const { return *this - t != Duration(0., UnitSec); }
      bool operator <(const AbsoluteTime & t) const { return *this - t < Duration(0., UnitSec); }
      bool operator >(const AbsoluteTime & t) const { return *this - t > Duration(0., UnitSec); }
      bool operator <=(const AbsoluteTime & t) const {
Duration d = operator -(t);
Duration zero(0., UnitSec);
bool r =  d <= zero;
return r;
return *this - t <= Duration(0., UnitSec);
}
      bool operator >=(const AbsoluteTime & t) const { return *this - t >= Duration(0., UnitSec); }
  };

  inline std::ostream & operator <<(std::ostream & os, const AbsoluteTime & t) {
    t.write(os);
    return os;
  }
}

#endif
