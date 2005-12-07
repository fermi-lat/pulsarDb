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
      Duration(long double t, TimeUnit_e unit): m_day(), m_sec() {
        if (unit == UnitDay) {
          m_day = std::floor(t);
          m_sec = (t - m_day) * SecPerDay();
        } else if (unit == UnitSec) {
          m_day = std::floor(t * DayPerSec());
          m_sec = t - (m_day * SecPerDay());
        }
      }

      /// \brief Return the current value of this time in days.
      long double day() const;

      /// \brief Return the current value of this time in seconds.
      long double sec() const;

      Duration op(long double (*func)(long double)) const;

      Duration operator +(const Duration & dur) const { Duration result = *this; result += dur; return result; }

      Duration operator +=(const Duration & dur) { m_day += dur.m_day; m_sec += dur.m_sec; return *this; }

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
      long double m_day;
      long double m_sec;
  };

  inline long double Duration::day() const {
    return m_day + m_sec * DayPerSec();
  }

  inline long double Duration::sec() const {
    return m_sec + m_day * SecPerDay();
  }

  inline Duration Duration::op(long double (*func)(long double)) const { return Duration(func(day()), UnitDay); }

  inline Duration & Duration::operator -=(const Duration & dur) {
    m_day -= dur.m_day;
    m_sec -= dur.m_sec;
    return *this;
  }

  inline Duration Duration::operator -(const Duration & dur) const {
    Duration result(*this);
    result -= dur;
    return result;
  }

  inline void Duration::write(std::ostream & os) const {
    double dday = day();
    os << dday << " day";
    if (dday == 1.) os << "s";
  }

  inline std::ostream & operator <<(std::ostream & os, const Duration & dur) {
    dur.write(os);
    return os;
  }

}

#endif
