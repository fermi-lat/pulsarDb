/** \file TimingModel.cxx
    \brief Implementation of TimingModel class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC
*/
#include <cmath>
#include <stdexcept>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "timeSystem/AbsoluteTime.h"

using namespace timeSystem;

#if 0
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
        double g,        /* input: mean anomaly */
        double eccent,        /* input: eccentricity */
        double *e)        /* output: eccentric anomaly */
{
    static double eps = 1e-6;
    static int imax = 50;

    int i;
    static double error, deltae, d__1;

    *e = g;
    if (g == 0.) return 0;

    for (i=0; i<imax; i++) {
        deltae = (g - *e + eccent * std::sin(*e)) / (1. - eccent * std::cos(*e));
        *e += deltae;
        error = (d__1 = deltae / *e, std::fabs(d__1));
        if (error < eps) return 0;
    }
    return NOT_CONVERGED;
}

}
#endif

namespace pulsarDb {

  double TimingModel::calcOrbitalPhase(const OrbitalEph & eph, const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    return eph.calcOrbitalPhase(ev_time, phase_offset);
#if 0
    // compute elapsed time from epoch of periastron in seconds
    double delta_second = eph.dt(ev_time);

    // Compute the time difference as a fraction of the period.
    double delta_period = delta_second / eph[PB];

    // Compute the complete phase.
    double phase = delta_period * (1. - delta_period * eph[PBDOT] / 2.0);

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    double int_part; // ignored, needed for modf.
    phase = std::modf(phase_offset + phase, &int_part);
    if (phase < 0.) ++phase;
    return phase;
#endif
  }

  void TimingModel::modulateBinary(const OrbitalEph & eph, timeSystem::AbsoluteTime & ev_time) const {
    eph.modulateBinary(ev_time);
#if 0
    ev_time += calcOrbitalDelay(eph, ev_time);
#endif
  }

  void TimingModel::demodulateBinary(const OrbitalEph & eph, timeSystem::AbsoluteTime & ev_time) const {
    eph.demodulateBinary(ev_time);
#if 0
    static const int s_max_iteration = 100;
    const ElapsedTime epsilon(eph.getSystem().getName(), Duration(0, 10.e-9)); // 10 nanoseconds.

    // Save target arrival time (ev_time) in orig_time.
    AbsoluteTime orig_time = ev_time;

    // Initial guess of orbital delay.
    ElapsedTime delay = calcOrbitalDelay(eph, ev_time);

    // Iterative approximation of demodulated time.
    int ii;
    for (ii=0; ii<s_max_iteration; ii++) {

      // Compute next candidate of demodulated time.
      ev_time = orig_time - delay;

      // Compute orbital delay at ev_time.
      delay = calcOrbitalDelay(eph, ev_time);

      // Compare time difference between candidate demodulated time
      // (ev_time) and target arrival time (orig_time) with the
      // estimated orbital delay based on the binary model (delay).
      if (orig_time.equivalentTo(ev_time + delay, epsilon)) break;
    }

    // Check for non-convergence.
    if (ii == s_max_iteration) throw std::runtime_error("Binary demodulation did not converge.");
#endif
  }

  TimingModel * TimingModel::clone() const {
    return new TimingModel(*this);
  }

  timeSystem::ElapsedTime TimingModel::calcOrbitalDelay(const OrbitalEph & eph, const timeSystem::AbsoluteTime & ev_time) const {
    return eph.calcOrbitalDelay(ev_time);
#if 0
    // compute elapsed time from epoch of periastron in seconds
    double delta_second = eph.dt(ev_time);

    // calculate mean anomaly
    double delta_period = delta_second / eph[PB];
    double mean_anomaly = OrbitalEph::s_two_pi * delta_period
      * (1. - delta_period * eph[PBDOT] / 2.0);

    // solve Kepler's equasion
    double eccen = eph[ECC] + eph[ECCDOT] * delta_second; // eccentricity
    double eccen_anomaly = 0.0; // eccentric anomaly
    int status = atKepler(mean_anomaly, eccen, &eccen_anomaly);

    // atKepler not converged
    if (0 != status) {
      throw std::runtime_error("atKepler did not converge.");
    }

    // convert eccentric anomaly to true anomaly
    double true_anomaly = 2.0 * std::atan(std::sqrt((1.0+eccen)/(1.0-eccen))
        * std::tan(eccen_anomaly/2.0));
    true_anomaly += OrbitalEph::s_two_pi * floor((eccen_anomaly - true_anomaly)/ OrbitalEph::s_two_pi);
    while ((true_anomaly - eccen_anomaly) > OrbitalEph::s_one_pi) true_anomaly -= OrbitalEph::s_two_pi;
    while ((eccen_anomaly - true_anomaly) > OrbitalEph::s_one_pi) true_anomaly += OrbitalEph::s_two_pi;

    // compute periastron longitude
    double omega = eph[OM]
      + eph[OMDOT] * true_anomaly * eph[PB] / OrbitalEph::s_two_pi;

    // compute projected semimajor axis
    double semiax = eph[A1] + eph[XDOT] * delta_second;

    // compute time delays due to orbital motion
    double roemer_frac = std::sin(omega) * (std::cos(eccen_anomaly) - eccen)
      + std::sqrt(1.0 - eccen*eccen) * std::cos(omega) * std::sin(eccen_anomaly);
    double roemer = semiax * roemer_frac;
    double einstein = eph[GAMMA] * std::sin(eccen_anomaly);
    double shapiro = - 2.0 * eph[SHAPIRO_R]
      * std::log(1.0 - eccen*std::cos(eccen_anomaly) - eph[SHAPIRO_S]*roemer_frac);

    // return total delay
    return ElapsedTime(eph.getSystem().getName(), Duration(0, roemer + einstein + shapiro));
#endif
  }
}
