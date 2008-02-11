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
      \brief Base class representing a single pulsar ephemeris.
  */
  class PulsarEph {
    public:
      virtual ~PulsarEph() {}

      virtual const timeSystem::TimeSystem & getSystem() const = 0;
      virtual const timeSystem::AbsoluteTime & getValidSince() const = 0;
      virtual const timeSystem::AbsoluteTime & getValidUntil() const = 0;
      virtual const timeSystem::AbsoluteTime & getEpoch() const = 0;

      /** \brief Compute the spin phase of the given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
          \param phase_offset Phase value to be added to the computed pulse phase.
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute the pulse frequency at a given time in the time system given to this object
                 upon its construction. Call getSytem method to obtain the time system to interpret
                 the return value of this method. The unit of frequency and its derivatives must be
                 based on seconds, i.e., Hz (sE-1), Hz/s (sE-2), and so on.
                 Note: validity of the ephemeris (valid since and valid until) are not checked. 
          \param ev_time Time of the event.
          \param derivative_order Order of time derivative of frequency to compute. Set 0 (zero) to obtain
                 frequency, 1 (one) the first time derivative of frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const = 0;

      /** \brief Compute the Right Ascension and Declination at a given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
      */
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /** \brief Create a copy of this object.
      */
      virtual PulsarEph * clone() const = 0;

      /** \brief Output text expression of this PulsarEph to a given output stream.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const = 0;

    protected:
      PulsarEph() {};

      /** \brief Compute the number of revolutions that the pulsar will make between the ephemeris epoch and the given time,
                 including a fractional part of it.
                 Note: The fractional part of this will be used as a pulse phase computed by calcPulsePhase method.
          \param ev_time Time of the event.
      */
      virtual double calcCycleCount(const timeSystem::AbsoluteTime & ev_time) const = 0;
  };

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) { return eph.write(os); }

  /** \class FrequencyEph
      \brief Class representing a single pulsar ephemeris expressed with three frequency coefficients.
             Note: f0, f1, f2 depend on the time system, while AbsoluteTime objects don't.
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
        double phi0, double f0, double f1, double f2): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_unit_time(0, 1.), m_ra(ra), m_dec(dec), m_phi0(phi0),
        m_f0(f0), m_f1(f1), m_f2(f2) {}

      virtual ~FrequencyEph() {}

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_epoch; }

      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      virtual double calcCycleCount(const timeSystem::AbsoluteTime & ev_time) const;

    private:
      // TODO: Remove these methods.
      double ra() const { return m_ra; }
      double dec() const { return m_dec; }
      double phi0() const { return m_phi0; }
      double f0() const { return m_f0; }
      double f1() const { return m_f1; }
      double f2() const { return m_f2; }

    private:
      // TODO: Rename dt method to calcElapsedSecond? Or computeElapsedSecond like in PulsarToolApp?
      virtual double dt(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
      virtual double dt(const timeSystem::AbsoluteTime & at) const { return dt(at, m_epoch); }

      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_epoch;
      // TODO: Remove the unit time, or else move it into subclasses, depending on how the design shakes out.
      // For now it is harmless but unneccessary.
      timeSystem::Duration m_unit_time;
      double m_ra;
      double m_dec;
      double m_phi0;
      double m_f0;
      double m_f1;
      double m_f2;
  };

  /** \class PeriodEph
      \brief Class representing a single pulsar ephemeris expressed with three period coefficients.
             Note: p0, p1, p2 depend on the time system, while AbsoluteTime objects don't.
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
        double phi0, double p0, double p1, double p2): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_unit_time(0, 1.), m_ra(ra), m_dec(dec), m_phi0(phi0),
        m_p0(p0), m_p1(p1), m_p2(p2) {}

      virtual ~PeriodEph() {}

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_epoch; }

      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      virtual double calcCycleCount(const timeSystem::AbsoluteTime & ev_time) const;

    private:
      // TODO: Remove these methods?
      double ra() const { return m_ra; }
      double dec() const { return m_dec; }
      double phi0() const { return m_phi0; }
      double f0() const { return 1. / m_p0; }
      double f1() const { return - m_p1 / (m_p0 * m_p0); }
      double f2() const { double p0sq = m_p0 * m_p0; return 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq; }

    private:
      // TODO: Rename dt method to calcElapsedSecond? Or computeElapsedSecond like in PulsarToolApp?
      virtual double dt(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
      virtual double dt(const timeSystem::AbsoluteTime & at) const { return dt(at, m_epoch); }

      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_epoch;
      // TODO: Remove the unit time, or else move it into subclasses, depending on how the design shakes out.
      // For now it is harmless but unneccessary.
      timeSystem::Duration m_unit_time;
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
