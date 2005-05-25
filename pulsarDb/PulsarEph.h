/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <vector>

#include "pulsarDb/AbsoluteTime.h"

namespace pulsarDb {

  class TimingModel;

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
      virtual double phi0() const = 0;
      virtual double f0() const = 0;
      virtual double f1() const = 0;
      virtual double f2() const = 0;
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

  /** \class DatabaseEph
      \brief Class representing a single pulsar ephemeris.
  */
  class DatabaseEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch The epoch (time origin).
          \param toa The pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      DatabaseEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until,
        const AbsoluteTime & epoch, const AbsoluteTime & toa, double f0, double f1, double f2);

      DatabaseEph(const DatabaseEph & eph): PulsarEph(eph), m_toa(eph.m_toa->clone()),
        m_f0(eph.m_f0), m_f1(eph.m_f1), m_f2(eph.m_f2) {}

      virtual ~DatabaseEph();

      DatabaseEph & operator =(const DatabaseEph & eph) {
        PulsarEph::operator =(eph);
        delete m_toa; m_toa = eph.m_toa->clone();
        m_f0 = eph.m_f0;
        m_f1 = eph.m_f1;
        m_f2 = eph.m_f2;
        return *this;
      }

      virtual double phi0() const;
      virtual double f0() const { return m_f0; }
      virtual double f1() const { return m_f1; }
      virtual double f2() const { return m_f2; }

      const AbsoluteTime & toa() const { return *m_toa; }
      virtual PulsarEph * clone() const { return new DatabaseEph(*this); }

    private:
      AbsoluteTime * m_toa;
      double m_f0;
      double m_f1;
      double m_f2;
  };

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
        double phi0, double f0, double f1, double f2): PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0), m_f0(f0),
        m_f1(f1), m_f2(f2) {}

      FrequencyEph(const FrequencyEph & eph): PulsarEph(eph), m_phi0(eph.m_phi0), m_f0(eph.m_f0), m_f1(eph.m_f1), m_f2(eph.m_f2) {}

      virtual ~FrequencyEph() {}

      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return m_f0; }
      virtual double f1() const { return m_f1; }
      virtual double f2() const { return m_f2; }
      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

    private:
      double m_phi0;
      double m_f0;
      double m_f1;
      double m_f2;
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
      PeriodEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until, const AbsoluteTime & epoch, double phi0,
        double p0, double p1, double p2): PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0),
        m_p0(p0), m_p1(p1), m_p2(p2) {}

      PeriodEph(const PeriodEph & eph): PulsarEph(eph), m_phi0(eph.m_phi0), m_p0(eph.m_p0), m_p1(eph.m_p1), m_p2(eph.m_p2) {}

      virtual ~PeriodEph() {}

      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return 1. / m_p0; }
      virtual double f1() const { return - m_p1 / (m_p0 * m_p0); }
      virtual double f2() const { double p0sq = m_p0 * m_p0; return 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq; }
      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

    private:
      double m_phi0;
      double m_p0;
      double m_p1;
      double m_p2;
  };

  typedef std::vector<PulsarEph *> PulsarEphCont;
}

#endif
