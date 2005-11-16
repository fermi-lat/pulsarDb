/** \file TimingModel.cxx
    \brief Implementation of TimingModel class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC
*/
#include <cmath>
#include <stdexcept>

#include <iostream>

#include "pulsarDb/CanonicalTime.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

namespace {

//======================================================================
// copied from atErrors.h
//----------------------------------------------------------------------
/* Error Code of atFunctions. */
#define NOT_CONVERGED 10        /*equation was not solved*/
//----------------------------------------------------------------------
// copied from atKepler.c (& modified)
//----------------------------------------------------------------------
//#include "atFunctions.h"
//#include "atError.h"
//#include <math.h>

/*
 * solve Kepler equation (KEPLER)  g + e sin E = E
 */
int atKepler(
        long double g,        /* input: mean anomaly */
        long double eccent,        /* input: eccentricity */
        long double *e)        /* output: eccentric anomaly */
{
    static long double eps = 1e-6;
    static int imax = 50;

    int i;
    static long double error, deltae, d__1;

    *e = g;
    if (g == 0.) return 0;

    for (i=0; i<imax; i++) {
        deltae = (g - *e + eccent * sin(*e)) / (1. - eccent * cos(*e));
        *e += deltae;
        error = (d__1 = deltae / *e, std::fabs(d__1));
        if (error < eps) return 0;
    }
    return NOT_CONVERGED;
}

}

namespace pulsarDb {

  long double TimingModel::calcOrbitalPhase(const OrbitalEph & eph, const AbsoluteTime & ev_time) const {
    // compute elapsed time from epoch of periastron in seconds
    long double delta_second = (ev_time - TdbTime(eph[T0])).sec();

    // Compute the time difference as a fraction of the period.
    long double delta_period = delta_second / eph[PB];

    // Compute the complete phase.
    long double phase = delta_period * (1. - delta_period * eph[PBDOT] / 2.0);

    // Express phase as a value between 0. and 1.
    long double int_part;
    phase = std::modf(phase, &int_part);
    if (phase < 0.) ++phase;
    return phase;
  }

  void TimingModel::modulateBinary(const OrbitalEph & eph, AbsoluteTime & ev_time) const {
    ev_time += calcOrbitalDelay(eph, ev_time);
  }

  void TimingModel::demodulateBinary(const OrbitalEph & eph, AbsoluteTime & ev_time) const {
    static const int s_max_iteration = 1000;
    const long double epsilon = 10.0e-9; // 10 nanoseconds.

    // Save target arrival time (ev_time) in orig_time.
    // TODO: Do not use TdbTime explicitly.  Generalize it.
    TdbTime orig_time(ev_time);

    // Initial guess of orbital delay.
    Duration delay = calcOrbitalDelay(eph, ev_time);

    // Iterative approximation of demodulated time.
    int ii;
    for (ii=0; ii<s_max_iteration; ii++) {

      // Compute next candidate of demodulated time by:
      // ev_time = orig_time - delay, i.e.:
      ev_time = orig_time;
      ev_time -= delay;

      // Compute orbital delay at ev_time
      delay = calcOrbitalDelay(eph, ev_time);

      // Compare time difference between candidate demodulated time
      // (ev_time) and target arrival time (orig_time) with the
      // estimated orbital delay based on the binary model (delay).
      if (std::fabs(((orig_time - ev_time) - delay).sec()) < epsilon) break;
    }

    // Check for non-convergence.
    if (ii == s_max_iteration) throw std::runtime_error("Binary demodulation did not converge.");

  }

  TimingModel * TimingModel::clone() const {
    return new TimingModel(*this);
  }

  Duration TimingModel::calcOrbitalDelay(const OrbitalEph & eph, const AbsoluteTime & ev_time) const {
    // compute elapsed time from epoch of periastron in seconds
    long double delta_second = (ev_time - TdbTime(eph[T0])).sec();

    // calculate mean anomaly
    long double delta_period = delta_second / eph[PB];
    long double mean_anomaly = OrbitalEph::s_two_pi * delta_period
      * (1. - delta_period * eph[PBDOT] / 2.0);

    // solve Kepler's equasion
    long double eccen = eph[ECC] + eph[ECCDOT] * delta_second; // eccentricity
    long double eccen_anomaly = 0.0; // eccentric anomaly
    int status = atKepler(mean_anomaly, eccen, &eccen_anomaly);

    // atKepler not converged
    if (0 != status) {
      throw std::runtime_error("atKepler did not converge.");
    }

    // convert eccentric anomaly to true anomaly
    long double true_anomaly = 2.0 * atan(sqrt((1.0+eccen)/(1.0-eccen))
                                     * tan(eccen_anomaly/2.0));
    true_anomaly += OrbitalEph::s_two_pi * floor((eccen_anomaly - true_anomaly)/ OrbitalEph::s_two_pi);
    while ((true_anomaly - eccen_anomaly) > OrbitalEph::s_one_pi) true_anomaly -= OrbitalEph::s_two_pi;
    while ((eccen_anomaly - true_anomaly) > OrbitalEph::s_one_pi) true_anomaly += OrbitalEph::s_two_pi;

    // compute periastron longitude
    long double omega = eph[OM]
      + eph[OMDOT] * true_anomaly * eph[PB] / OrbitalEph::s_two_pi;

    // compute projected semimajor axis
    long double semiax = eph[A1] + eph[XDOT] * delta_second;

    // compute time delays due to orbital motion
    long double roemer_frac = sin(omega) * (cos(eccen_anomaly) - eccen)
      + sqrt(1.0 - eccen*eccen) * cos(omega) * sin(eccen_anomaly);
    long double roemer = semiax * roemer_frac;
    long double einstein = eph[GAMMA] * sin(eccen_anomaly);
    long double shapiro = - 2.0 * eph[SHAPIRO_R]
      * log(1.0 - eccen*cos(eccen_anomaly) - eph[SHAPIRO_S]*roemer_frac);

    // return total delay
    return Duration(roemer + einstein + shapiro, UnitSec);
  }

}
