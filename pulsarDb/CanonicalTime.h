/** \file CanonicalTime
    \brief Declaration of CanonicalTime class, and some related subclasses.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_CanonicalTime_h
#define pulsarDb_CanonicalTime_h

#include <stdexcept>
#include <string>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/Duration.h"

namespace pulsarDb {

  class TaiTime;

  class TtTime;

  /** \class CanonicalTime
      \brief An absolute moment in time represented in one of a limited set of independent time systems
             which have well-established relationships with each other.
  */
  class CanonicalTime : public AbsoluteTime {
    public:
      /** \brief Assign this object's time to another CanonicalTime object.
          \param t The destination time object.
      */
      virtual void to(CanonicalTime & t) const;

      virtual void write(std::ostream & os) const;

      /** \brief Assign this object's time to a TaiTime object.
          \param t The destination time object.
      */
      virtual void toSystem(TaiTime & t) const = 0;

      /** \brief Assign this object's time to a TtTime object.
          \param t The destination time object.
      */
      virtual void toSystem(TtTime & t) const = 0;

      /// \brief Return a string identifying the name of time system used by this object.
      virtual const char * timeSystemName() const = 0;

      virtual long double mjd() const = 0;

      virtual double value() const { return mjd(); }

  };

  /** \class TaiTime
      \brief An absolute moment in time represented in the TAI system.
  */
  class TaiTime : public CanonicalTime {
    public:
      /** \brief Create a TaiTime object from the given MJD time.
          \param mjd The time expressed in MJD.
      */
      TaiTime(long double mjd): m_mjd(mjd) {}

      /** \brief Create a TaiTime object from the given AbsoluteTime object.
          \param t The absolute time object.
      */
      TaiTime(const AbsoluteTime & t): m_mjd(0.) { t.to(*this); }

      /** \brief Determine the elapsed time difference between the given time object and this one.
          \param t The absolute time object.
      */
      virtual Duration operator -(const AbsoluteTime & t) const;

      /** \brief Add the given time duration to this time.
          \param d The duration to add.
      */
      virtual AbsoluteTime & operator +=(const Duration & d);

      virtual AbsoluteTime * clone() const { return new TaiTime(*this); }

      virtual void toSystem(TaiTime & t) const { t.m_mjd = m_mjd; }

      virtual void toSystem(TtTime & t) const;

      virtual const char * timeSystemName() const { return "TAI"; }

      virtual long double mjd() const { return m_mjd; }

      void setMjd(long double mjd) { m_mjd = mjd; }

    private:
      long double m_mjd;
  };

  class TtTime : public CanonicalTime {
    public:
      TtTime(long double mjd): m_mjd(mjd) {}

      TtTime(const AbsoluteTime & t): m_mjd(0.) { t.to(*this); }

      virtual Duration operator -(const AbsoluteTime & t) const;

      virtual AbsoluteTime & operator +=(const Duration & d);

      virtual AbsoluteTime * clone() const { return new TtTime(*this); }

      virtual void toSystem(TaiTime & t) const { t.setMjd(m_mjd + TaiMinusTtDay()); }

      virtual void toSystem(TtTime & t) const { t.m_mjd = m_mjd; }

      virtual const char * timeSystemName() const { return "TT"; }

      virtual long double mjd() const { return m_mjd; }

      void setMjd(long double mjd) { m_mjd = mjd; }

    private:
      long double m_mjd;
  };

  // TODO: This is obviously not right!!!!!
  typedef TtTime TdbTime;

  inline void CanonicalTime::to(CanonicalTime & t) const {
    TaiTime * tai_t = dynamic_cast<TaiTime *>(&t);
    if (0 != tai_t) {
      toSystem(*tai_t);
      return;
    }

    TtTime * tt_t = dynamic_cast<TtTime *>(&t);
    if (0 != tt_t) {
      toSystem(*tt_t);
      return;
    }

    throw std::logic_error(std::string("Cannot convert from time system ") + timeSystemName() + " to " + t.timeSystemName());
  }

  inline void CanonicalTime::write(std::ostream & os) const {
    os << mjd() << " MJD " << timeSystemName();
  }

  inline Duration TaiTime::operator -(const AbsoluteTime & t) const {
    TaiTime my_sys_time(t);
    return Duration(m_mjd - my_sys_time.m_mjd, UnitDay);
  }

  inline AbsoluteTime & TaiTime::operator +=(const Duration & d) {
    m_mjd += d.day();
    return *this;
  }

  inline void TaiTime::toSystem(TtTime & t) const { t.setMjd(m_mjd + TtMinusTaiDay()); }

  inline Duration TtTime::operator -(const AbsoluteTime & t) const {
    TtTime my_sys_time(t);
    return Duration(m_mjd - my_sys_time.m_mjd, UnitDay);
  }

  inline AbsoluteTime & TtTime::operator +=(const Duration & d) {
    m_mjd += d.day();
    return *this;
  }

}

#endif
