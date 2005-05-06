/** \file Duration
    \brief Declaration of Duration class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_Duration_h
#define pulsarDb_Duration_h

#include <cmath>
#include <iostream>

#include "pulsarDb/TimeConstants.h"

namespace pulsarDb {

  /// \enum TimeUnit_e Enumerated type representing different time units.
  enum TimeUnit_e {
    UnitDay,
    UnitSec
  };

  /** \class Duration
      \brief Class representing a measureable length of time, i.e. elapsed time,
             or a difference between two absolute times. Note that the first AbsoluteTime object
             used determines the time system in which the duration is valid!
  */
  class Duration {
    public:
      /** \brief Construct a duration object
          \param t The amount of time.
          \param unit The unit in which time t is measured.
      */
      Duration(double t, TimeUnit_e unit): m_t(t), m_unit(unit) {}

      /// \brief Return the current value of this time in days.
      double day() const;

      /// \brief Return the current value of this time in seconds.
      double sec() const;

      Duration op(double (*func)(double)) const;

      Duration & operator -=(const Duration & dur);

      Duration operator -(const Duration & dur) const;

      bool operator !=(const Duration & dur) const { return sec() != dur.sec(); }
      bool operator ==(const Duration & dur) const { return sec() == dur.sec(); }
      bool operator <(const Duration & dur) const { return sec() < dur.sec(); }
      bool operator <=(const Duration & dur) const { return sec() <= dur.sec(); }
      bool operator >(const Duration & dur) const { return sec() > dur.sec(); }
      bool operator >=(const Duration & dur) const { return sec() >= dur.sec(); }

      void write(std::ostream & os) const;

    private:
      double m_t;
      TimeUnit_e m_unit;
  };

  inline double Duration::day() const {
    switch (m_unit) {
      case UnitDay:
        break; // No modifcation necessary.
      case UnitSec:
        return m_t * DayPerSec();
    }
    return m_t;
  }

  inline double Duration::sec() const {
    switch (m_unit) {
      case UnitDay:
        return m_t * SecPerDay();
      case UnitSec:
        break; // No modifcation necessary.
    }
    return m_t;
  }

  inline Duration Duration::op(double (*func)(double)) const { return Duration(func(m_t), m_unit); }

  inline Duration & Duration::operator -=(const Duration & dur) {
    switch (m_unit) {
      case UnitDay:
        m_t -= dur.day();
        break;
      case UnitSec:
        m_t -= dur.sec();
        break;
    }
    return *this;
  }

  inline Duration Duration::operator -(const Duration & dur) const {
    Duration result(*this);
    result -= dur;
    return result;
  }

  inline void Duration::write(std::ostream & os) const {
    os << m_t << (m_unit == UnitDay ? " day" : " second");
    if (m_t == 1.) os << "s";
  }

  inline std::ostream & operator <<(std::ostream & os, const Duration & dur) {
    dur.write(os);
    return os;
  }

}

#endif
