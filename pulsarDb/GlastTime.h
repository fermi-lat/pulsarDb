/** \file GlastTime.h
    \brief Declaration of GLAST mission elapsed time classes.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_GlastTime_h
#define pulsarDb_GlastTime_h

#include <iostream>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/Duration.h"
#include "pulsarDb/CanonicalTime.h"

#include "st_stream/Stream.h"

namespace pulsarDb {

  /** \class GlastTime
      \brief An absolute moment in time represented in one of a limited set of independent time systems
             which have well-established relationships with each other.
  */
  template <typename ACanonicalTime>
  class GlastTime : public AbsoluteTime {
    public:
      GlastTime(long double elapsed = 0.): m_elapsed(0, elapsed) {}

      GlastTime(const AbsoluteTime & t);

      /** \brief Subtract the given time from this one, to determine the duration of the interval between them.
          \param t The other time.
      */
      virtual Duration operator -(const AbsoluteTime & t) const;

      /** \brief Add the given time increment to this one.
          \param d The duration (elapsed time, relative time) being added.
      */
      virtual AbsoluteTime & operator +=(const Duration & d) { m_elapsed += d; return *this; }

      /** \brief Subtract the given duration from this time.
          \param d The duration (elapsed time, relative time) being added.
      */
      virtual AbsoluteTime & operator -=(const Duration & d) { m_elapsed -= d; return *this; }

      /** \brief Convert this time to a CanonicalTime representation of it.
          \param t The target time system time.
      */
      virtual void to(CanonicalTime & t) const;

      inline void from(const AbsoluteTime & t);

      virtual void write(std::ostream & os) const;

      virtual void write(st_stream::OStream & os) const;

      virtual AbsoluteTime * clone() const { return new GlastTime(*this); }

      long double elapsed() const { return m_elapsed.sec(); }

      virtual long double value() const { return m_elapsed.sec(); }

    private:
      Duration m_elapsed;
  };

  // TODO: This should be a member of GlastTime? Or a constant? Or in some other glast constant class?
  static long s_mjdref = 51910;

  template <typename ACanonicalTime>
  inline GlastTime<ACanonicalTime>::GlastTime(const AbsoluteTime & t): m_elapsed(0, 0.) { from(t); }

  template <typename ACanonicalTime>
  inline Duration GlastTime<ACanonicalTime>::operator -(const AbsoluteTime & t) const {
    return m_elapsed - GlastTime(t).m_elapsed;
  }

  /** \brief Convert this time to a CanonicalTime representation of it.
      \param t The target time system time.
  */
  template <typename ACanonicalTime>
  inline void GlastTime<ACanonicalTime>::to(CanonicalTime & t) const {
    // Convert this time to the canonical time system implied by this GlastTime, and use that to convert to the target time.
    ACanonicalTime my_time(s_mjdref);
    my_time += m_elapsed;
    my_time.to(t);
  }

  template <typename ACanonicalTime>
  inline void GlastTime<ACanonicalTime>::from(const AbsoluteTime & t) {
    // See if other time happens to be already in GlastTime.
    const GlastTime * glast_t = dynamic_cast<const GlastTime *>(&t);
    if (0 != glast_t) {
      m_elapsed = glast_t->m_elapsed;
    } else {
      // Convert the other time to the canonical time system implied by this GlastTime, and assign from that.
      ACanonicalTime other_t(t);
      m_elapsed = other_t.getMjd() - Duration(s_mjdref, 0.);
    }
  }

  template <typename ACanonicalTime>
  inline void GlastTime<ACanonicalTime>::write(std::ostream & os) const {
    static const ACanonicalTime canonical_time(0.);
    os << m_elapsed.sec() << " seconds GLAST MET " << canonical_time.timeSystemName();
  }

  template <typename ACanonicalTime>
  inline void GlastTime<ACanonicalTime>::write(st_stream::OStream & os) const {
    static const ACanonicalTime canonical_time(0.);
    os << m_elapsed.sec() << " seconds GLAST MET " << canonical_time.timeSystemName();
  }

  typedef GlastTime<TdbTime> GlastTdbTime;

  typedef GlastTime<TtTime> GlastTtTime;

}

#endif
