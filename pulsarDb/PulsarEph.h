/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "st_stream/Stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/TimeSystem.h"

namespace pulsarDb {

  /** \class PulsarEph
      \brief Class representing a single pulsar ephemeris. Warning: f0, f1, f2 depend on the time system.
      While AbsoluteTime objects are interchangeable, these other values are not!
  */
  class PulsarEph {
    public:
      PulsarEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch):
        m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_unit_time(0, 1.) {}

      virtual ~PulsarEph() {}

      virtual const timeSystem::AbsoluteTime & valid_since() const { return m_since; }
      virtual const timeSystem::AbsoluteTime & valid_until() const { return m_until; }
      virtual const timeSystem::AbsoluteTime & epoch() const { return m_epoch; }
      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }
      virtual double ra() const = 0;
      virtual double dec() const = 0;
      virtual double phi0() const = 0;
      virtual double f0() const = 0;
      virtual double f1() const = 0;
      virtual double f2() const = 0;
      virtual PulsarEph * clone() const = 0;

      /** \brief Compute the spin phase of the given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
          \param phase_offset Phase value to be added to the computed pulse phase.
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      /** \brief Compute the pulse frequency at a given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
          \param derivative_order Order of time derivative of frequency to compute. Set 0 (zero) to obtain
                 frequency, 1 (one) the first time derivative of frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const = 0;

      /** \brief Compute the Right Ascension and Declination at a given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
      */
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const = 0;

    protected:
      virtual double dt(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
      virtual double dt(const timeSystem::AbsoluteTime & at) const { return dt(at, m_epoch); }

      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_epoch;
      // TODO: Remove the unit time, or else move it into subclasses, depending on how the design shakes out.
      // For now it is harmless but unneccessary.
      timeSystem::Duration m_unit_time;
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
      FrequencyEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch, double ra, double dec,
        double phi0, double f0, double f1, double f2): PulsarEph(time_system_name, valid_since, valid_until, epoch),
          m_ra(ra), m_dec(dec), m_phi0(phi0), m_f0(f0), m_f1(f1), m_f2(f2) {}

      FrequencyEph(const FrequencyEph & eph): PulsarEph(eph), m_ra(eph.m_ra), m_dec(eph.m_dec), m_phi0(eph.m_phi0),
        m_f0(eph.m_f0), m_f1(eph.m_f1), m_f2(eph.m_f2) {}

      virtual ~FrequencyEph() {}

      virtual double ra() const { return m_ra; }
      virtual double dec() const { return m_dec; }
      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return m_f0; }
      virtual double f1() const { return m_f1; }
      virtual double f2() const { return m_f2; }
      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      // TODO: Write a test code for this method.
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      // TODO: Write a test code for this method.
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

    private:
      double m_ra;
      double m_dec;
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
      PeriodEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch, double ra, double dec,
        double phi0, double p0, double p1, double p2): PulsarEph(time_system_name, valid_since, valid_until, epoch),
        m_ra(ra), m_dec(dec), m_phi0(phi0), m_p0(p0), m_p1(p1), m_p2(p2) {}

      PeriodEph(const PeriodEph & eph): PulsarEph(eph), m_ra(eph.m_ra), m_dec(eph.m_dec), m_phi0(eph.m_phi0),
        m_p0(eph.m_p0), m_p1(eph.m_p1), m_p2(eph.m_p2) {}

      virtual ~PeriodEph() {}

      virtual double ra() const { return m_ra; }
      virtual double dec() const { return m_dec; }
      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return 1. / m_p0; }
      virtual double f1() const { return - m_p1 / (m_p0 * m_p0); }
      virtual double f2() const { double p0sq = m_p0 * m_p0; return 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq; }
      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      // TODO: Write a test code for this method.
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      // TODO: Write a test code for this method.
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

    private:
      double m_ra;
      double m_dec;
      double m_phi0;
      double m_p0;
      double m_p1;
      double m_p2;
  };

  typedef std::vector<PulsarEph *> PulsarEphCont;
}

#endif
