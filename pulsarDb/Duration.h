/** \file Duration
    \brief Declaration of Duration class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_Duration_h
#define pulsarDb_Duration_h

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "pulsarDb/TimeConstants.h"

#include "st_stream/Stream.h"

namespace pulsarDb {

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
      Duration(long day, double sec): m_time(add(time_type(day, 0.), splitSec(sec))) {}

      /// \brief Return the current value of this time in days.
      long double day() const;

      /// \brief Return the current value of this time in seconds.
      long double sec() const;

      Duration op(long double (*func)(long double)) const;

      Duration operator +(const Duration & dur) const;

      Duration & operator +=(const Duration & dur);

      Duration & operator -=(const Duration & dur);

      Duration operator -(const Duration & dur) const;

      Duration operator -() const;

      bool operator !=(const Duration & dur) const
        { return m_time.first != dur.m_time.first && m_time.second != dur.m_time.second; }
      bool operator ==(const Duration & dur) const
        { return m_time.first == dur.m_time.first && m_time.second == dur.m_time.second; }
      bool operator <(const Duration & dur) const
        { return m_time.first != dur.m_time.first ? m_time.first < dur.m_time.first : m_time.second < dur.m_time.second; }
      bool operator <=(const Duration & dur) const
        { return m_time.first != dur.m_time.first ? m_time.first < dur.m_time.first : m_time.second <= dur.m_time.second; }
      bool operator >(const Duration & dur) const
        { return m_time.first != dur.m_time.first ? m_time.first > dur.m_time.first : m_time.second > dur.m_time.second; }
      bool operator >=(const Duration & dur) const
        { return m_time.first != dur.m_time.first ? m_time.first > dur.m_time.first : m_time.second >= dur.m_time.second; }

      void write(std::ostream & os) const;

      void write(st_stream::OStream & os) const;

    private:
      typedef std::pair<long, double> time_type;

      Duration(const time_type & new_time): m_time(new_time) {}

#if 0
      // This is commented out because there are problems converting days into days + seconds
      // without losing precision. Not needed with current scheme because constructors work
      // directly from long values for days.
      /** \brief Convert any number of days into days & seconds in range [0, 86400).
          \param day Input number of days.
      */
      time_type splitDay(double day) const {
        double offset;
        if (0. > day) {
          offset = -.5;
        } else {
          offset = +.5;
        }
        double int_part;
        double frac_part = std::modf(day, &int_part);
        int num_digit_all = std::numeric_limits<double>::digits10;
        int num_digit_int = int_part == 0 ? 0 : int(floor(log10(fabs(int_part))) + 0.5) + 1;
        int num_digit_frac = num_digit_all - num_digit_int;
        double factor = floor(exp(num_digit_frac * log(10.0)));
        frac_part = floor(frac_part * factor) / factor;
        return time_type(long(int_part + offset), frac_part * SecPerDay());

        long del = long(std::floor(day) + offset);
        return time_type(del, (day - del) * SecPerDay());
      }
#endif

      /** \brief Convert any number of seconds into days & seconds in range [0, 86400).
          \param sec Input number of seconds.
      */
      time_type splitSec(double sec) const {
        double offset;
        if (0. > sec) {
          offset = -.5;
        } else if (SecPerDay() <= sec) {
          offset = +.5;
        } else {
          return time_type(0, sec);
        }
        double double_day = std::floor(sec * DayPerSec()) + offset;
        long day;
        if (std::numeric_limits<long>::min() <= double_day && double_day <= std::numeric_limits<long>::max()) {
          day = long(double_day);
        } else if (std::numeric_limits<long>::min() > double_day) {
          std::ostringstream os;
          os << "Time " << double_day << " is too small to be expressed as a long integer";
          throw std::runtime_error(os.str());
        } else { // std::numeric_limits<long>::max() < double_day
          std::ostringstream os;
          os << "Time " << double_day << " is too large to be expressed as a long integer";
          throw std::runtime_error(os.str());
        }
        return time_type(day, sec - day * SecPerDay());
      }

      /** \brief Add two times which are represented by long day and double second fields. Seconds
            part of the result is guaranteed to be in the range [0., SecPerDay())
          \param t1 The first time being added.
          \param t2 The second time being added.
      */
      time_type add(time_type t1, time_type t2) const {
        // Day portions are added no matter what.
        long day = t1.first + t2.first;
    
        // Sum the two seconds portions.
        double sec = t1.second + t2.second;
    
        // Check for overflow.
        if (sec >= SecPerDay()) {
          ++day;
          // Do not reuse sum variable from above, in order to preserve maximum precision.
          sec = (t1.second - SecPerDay()) + t2.second;
        }
        return time_type(day, sec);
      }

      /** \brief Multiply by -1 a time represented by long day and double second fields. Seconds
            part of the result is guaranteed to be in the range [0., SecPerDay())
          \param t1 The first time being negated.
      */
      time_type negate(time_type t1) const {
        return time_type(-t1.first - 1, SecPerDay() - t1.second);
      }

      time_type m_time;
  };

  inline long double Duration::day() const {
    return m_time.first + m_time.second * DayPerSec();
  }

  inline long double Duration::sec() const {
    return m_time.second + m_time.first * SecPerDay();
  }

  inline Duration Duration::operator +(const Duration & dur) const {
    return Duration(add(m_time, dur.m_time));
  }

  inline Duration & Duration::operator +=(const Duration & dur) {
    m_time = add(m_time, dur.m_time);
    return *this;
  }

  inline Duration Duration::op(long double (*func)(long double)) const { return Duration(0, func(sec())); }

  inline Duration & Duration::operator -=(const Duration & dur) {
    m_time = add(m_time, negate(dur.m_time));
    return *this;
  }

  inline Duration Duration::operator -(const Duration & dur) const {
    return Duration(add(m_time, negate(dur.m_time)));
  }

  inline Duration Duration::operator -() const {
    return Duration(negate(m_time));
  }

  inline void Duration::write(std::ostream & os) const {
    os << m_time.first << " day";
    if (m_time.first != 1) os << "s";
    std::streamsize prec = os.precision(std::numeric_limits<double>::digits10);
    os << ", " << m_time.second << " second";
    if (m_time.second != 1.) os << "s";
    os.precision(prec);
  }

  inline void Duration::write(st_stream::OStream & os) const {
    os << m_time.first << " day";
    if (m_time.first != 1) os << "s";
    std::streamsize prec = os.precision(std::numeric_limits<double>::digits10);
    os << ", " << m_time.second << " second";
    if (m_time.second != 1.) os << "s";
    os.precision(prec);
  }

  inline std::ostream & operator <<(std::ostream & os, const Duration & dur) {
    dur.write(os);
    return os;
  }

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const Duration & dur) {
    dur.write(os);
    return os;
  }

}

#endif
