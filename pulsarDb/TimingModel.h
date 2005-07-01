/** \file TimingModel.h
    \brief Encapsulation of standard timing model.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC
*/
#ifndef pulsarDb_TimingModel_h
#define pulsarDb_TimingModel_h

#include <cmath>

#include "pulsarDb/Duration.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace pulsarDb {

  class TimingModel {
    public:
      virtual ~TimingModel() {}

      /** \brief Compute frequency and its derivatives at a given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      virtual FrequencyEph calcEphemeris(const PulsarEph & eph, const AbsoluteTime & ev_time) const {
        double dt = (ev_time - eph.epoch()).sec();
        double f0 = eph.f0() + eph.f1() * dt + eph.f2()/2.0 * dt * dt;
        double f1 = eph.f1() + eph.f2() * dt;
        double phi0 = calcPhase(eph, ev_time);
        return FrequencyEph(eph.valid_since(), eph.valid_until(), eph.epoch(), phi0, f0, f1, eph.f2());
      }

      /** \brief Correct event time to account for pdot cancellation. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      // TODO: rename correctPdot cancelPdot?
      virtual void correctPdot(const PulsarEph & eph, AbsoluteTime & ev_time) const {
        double dt = (ev_time - eph.epoch()).sec();
        double dt_squared = dt * dt;
        ev_time += Duration(eph.f1()/eph.f0()/2.0 * dt_squared + eph.f2()/eph.f0()/6.0 * dt * dt_squared, UnitSec);
      }

      /** \brief Compute the phase of the given time. Note: validity of the
                 ephemeris (valid since and valid until) are not checked.
          \param eph The ephemeris.
          \param ev_time Time of the event.
      */
      // TODO: create calcOrbitalPhase.
      // TODO: related to calcOrbitalPhase, rename calcPhase calcPulsePhase?
      virtual double calcPhase(const PulsarEph & eph, const AbsoluteTime & ev_time) const {
        double dt = (ev_time - eph.epoch()).sec();
        double int_part; // ignored, needed for modf.
        double dt_squared = dt * dt;
        double phase = modf(eph.phi0() + eph.f0() * dt + eph.f1()/2.0 * dt_squared + eph.f2()/6.0 * dt * dt_squared, &int_part);
        if (phase < 0.) ++phase;
        return phase;
      }

      /** \brief
          \param eph The ephemeris.
          \param ev_time Before this method executes, this is the time the photon was emitted at the pulsar.
          After this method executes, this is the time the photon was detected by the instrument.
      */
      virtual void modulateBinary(const OrbitalEph & eph, AbsoluteTime & emission_time) const;

      /** \brief
          \param eph The ephemeris.
          \param ev_time Before this method executes, this is the time the photon was detected by the instrument.
          After this method executes, this is the time the photon was emitted at the pulsar.
      */
      virtual void demodulateBinary(const OrbitalEph & eph, AbsoluteTime & arrival_time) const;

      virtual TimingModel * clone() const;

    protected:
      Duration calcOrbitalDelay(const OrbitalEph & eph, const AbsoluteTime & emission_time) const;
  };

}

#endif
