/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <iostream>
#include <vector>

#include "pulsarDb/AbsoluteTime.h"
#include "st_stream/Stream.h"

namespace pulsarDb {

  /** \class PulsarEph
      \brief Class representing a single pulsar ephemeris. Warning: f0, f1, f2 depend on the time system.
      While AbsoluteTime objects are interchangeable, these other values are not!
  */
  class PulsarEph {
    public:
      PulsarEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until, const AbsoluteTime & epoch):
        m_since(valid_since.clone()), m_until(valid_until.clone()), m_epoch(epoch.clone()) {}

      virtual ~PulsarEph() { delete m_epoch; delete m_until; delete m_since; }

      virtual const AbsoluteTime & valid_since() const { return *m_since; }
      virtual const AbsoluteTime & valid_until() const { return *m_until; }
      virtual const AbsoluteTime & epoch() const { return *m_epoch; }
      virtual long double phi0() const = 0;
      virtual long double f0() const = 0;
      virtual long double f1() const = 0;
      virtual long double f2() const = 0;
      virtual PulsarEph * clone() const = 0;

    protected:
      AbsoluteTime * m_since;
      AbsoluteTime * m_until;
      AbsoluteTime * m_epoch;
      PulsarEph(const PulsarEph & eph): m_since(eph.m_since->clone()), m_until(eph.m_until->clone()),
        m_epoch(eph.m_epoch->clone()) {}
      PulsarEph & operator =(const PulsarEph & eph) {
        delete m_since; m_since = eph.m_since->clone();
        delete m_until; m_until = eph.m_until->clone();
        delete m_epoch; m_epoch = eph.m_epoch->clone();
        return *this;
      }
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph);

  /** \class FrequencyEph
      \brief Class representing a single pulsar ephemeris.
  */
  class FrequencyEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param epoch The epoch (time origin).
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      FrequencyEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until, const AbsoluteTime & epoch,
        long double phi0, long double f0, long double f1, long double f2): PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0), m_f0(f0),
        m_f1(f1), m_f2(f2) {}

      FrequencyEph(const FrequencyEph & eph): PulsarEph(eph), m_phi0(eph.m_phi0), m_f0(eph.m_f0), m_f1(eph.m_f1), m_f2(eph.m_f2) {}

      virtual ~FrequencyEph() {}

      virtual long double phi0() const { return m_phi0; }
      virtual long double f0() const { return m_f0; }
      virtual long double f1() const { return m_f1; }
      virtual long double f2() const { return m_f2; }
      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

    private:
      long double m_phi0;
      long double m_f0;
      long double m_f1;
      long double m_f2;
  };

  /** \class PeriodEph
      \brief Class representing a single pulsar ephemeris.
  */
  class PeriodEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param epoch The epoch (time origin).
          \param p0 The period at the epoch (time origin).
          \param p1 The first time derivative of the period at the epoch (time origin).
          \param p2 The second time derivative of the period at the epoch (time origin).
      */
      PeriodEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until, const AbsoluteTime & epoch, long double phi0,
        long double p0, long double p1, long double p2): PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0),
        m_p0(p0), m_p1(p1), m_p2(p2) {}

      PeriodEph(const PeriodEph & eph): PulsarEph(eph), m_phi0(eph.m_phi0), m_p0(eph.m_p0), m_p1(eph.m_p1), m_p2(eph.m_p2) {}

      virtual ~PeriodEph() {}

      virtual long double phi0() const { return m_phi0; }
      virtual long double f0() const { return 1. / m_p0; }
      virtual long double f1() const { return - m_p1 / (m_p0 * m_p0); }
      virtual long double f2() const { long double p0sq = m_p0 * m_p0; return 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq; }
      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

    private:
      long double m_phi0;
      long double m_p0;
      long double m_p1;
      long double m_p2;
  };

  typedef std::vector<PulsarEph *> PulsarEphCont;
}

#endif
