/** \file TimingModel.h
    \brief Encapsulation of standard timing model.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC
*/
#ifndef pulsarDb_TimingModel_h
#define pulsarDb_TimingModel_h

#include <cmath>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeSystem.h"

namespace pulsarDb {

  class TimingModel {
    public:
      virtual ~TimingModel() {}

      /** \brief Compute frequency and its derivatives at a given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      virtual FrequencyEph calcEphemeris(const PulsarEph & eph, const timeSystem::AbsoluteTime & ev_time) const {
        double dt = eph.dt(ev_time);
        long double f0 = eph.f0() + eph.f1() * dt + eph.f2()/2.0 * dt * dt;
        long double f1 = eph.f1() + eph.f2() * dt;
        long double phi0 = calcPulsePhase(eph, ev_time);
        return FrequencyEph(eph.getSystem().getName(), ev_time, ev_time, ev_time, phi0, f0, f1, eph.f2());
      }

      /** \brief Correct event time to account for pdot cancellation. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      virtual void cancelPdot(const PulsarEph & eph, timeSystem::AbsoluteTime & ev_time) const {
        double dt = eph.dt(ev_time);
        long double dt_squared = dt * dt;
        timeSystem::Duration corrected_dt(0, eph.f1()/eph.f0()/2.0 * dt_squared + eph.f2()/eph.f0()/6.0 * dt * dt_squared);
        ev_time = ev_time + timeSystem::ElapsedTime(eph.getSystem().getName(), corrected_dt);
      }

      /** \brief Compute the spin phase of the given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      virtual long double calcPulsePhase(const PulsarEph & eph, const timeSystem::AbsoluteTime & ev_time) const {
        double dt = eph.dt(ev_time);
        long double dt_squared = dt * dt;
        double int_part = 0.;
        long double phase =
          std::modf(eph.phi0() + eph.f0() * dt + eph.f1()/2.0 * dt_squared + eph.f2()/6.0 * dt * dt_squared, &int_part);
        if (phase < 0.) ++phase;
        return phase;
      }

      /** \brief Compute the orbital phase of the given time.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      virtual long double calcOrbitalPhase(const OrbitalEph & eph, const timeSystem::AbsoluteTime & ev_time) const;

      /** \brief
          \param eph The ephemeris.
          \param ev_time Before this method executes, this is the time the photon was emitted at the pulsar.
          After this method executes, this is the time the photon was detected by the instrument.
      */
      virtual void modulateBinary(const OrbitalEph & eph, timeSystem::AbsoluteTime & emission_time) const;

      /** \brief
          \param eph The ephemeris.
          \param ev_time Before this method executes, this is the time the photon was detected by the instrument.
          After this method executes, this is the time the photon was emitted at the pulsar.
      */
      virtual void demodulateBinary(const OrbitalEph & eph, timeSystem::AbsoluteTime & arrival_time) const;

      virtual TimingModel * clone() const;

    protected:
      timeSystem::ElapsedTime calcOrbitalDelay(const OrbitalEph & eph, const timeSystem::AbsoluteTime & emission_time) const;
  };

}

#endif
