/** \file BtModelEph.cxx
    \brief Implementation of the BtModelEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <limits>
#include <iomanip>
#include <stdexcept>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/BtModelEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/MjdFormat.h"

using namespace timeSystem;

namespace tip {
  class Header;
}

namespace {

//======================================================================
// Copied from atErrors.h.
//----------------------------------------------------------------------
/* Error Code of atFunctions. */
#define NOT_CONVERGED 10        /*equation was not solved*/
//----------------------------------------------------------------------
// Copied from atKepler.c (& modified).
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

namespace pulsarDb {

  const double BtModelEph::s_one_pi = M_PI;
  const double BtModelEph::s_two_pi = 2. * BtModelEph::s_one_pi;
  const double BtModelEph::s_rad_per_deg  = BtModelEph::s_one_pi / 180.;
  const double BtModelEph::s_sec_per_day  = 86400.;
  const double BtModelEph::s_sec_per_year = 365.25 * BtModelEph::s_sec_per_day;
  const double BtModelEph::s_rad_year_per_deg_sec = BtModelEph::s_rad_per_deg / BtModelEph::s_sec_per_year;

  BtModelEph::BtModelEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot,
    double ecc, double ecc_dot, double om, double om_dot, const AbsoluteTime & t0, double gamma):
    OrbitalEph(ElapsedTime(time_system_name, Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem(time_system_name)),
    m_pb(pb), m_pb_dot(pb_dot), m_a1(a1), m_x_dot(x_dot), m_ecc(ecc), m_ecc_dot(ecc_dot), m_om(om * s_rad_per_deg),
    m_om_dot(om_dot * s_rad_year_per_deg_sec), m_t0(t0), m_gamma(gamma) {}

  BtModelEph::BtModelEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    OrbitalEph(ElapsedTime("TDB", Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem("TDB")), m_t0("TDB", 0, 0.) {
    // Get parameters from record.
    // Required fields: PB, A1, ECC, OM, T0.
    // Optional fields: PBDOT, XDOT, ECCDOT, OMDOT, GAMMA.
    read(record, "PB",        m_pb);
    read(record, "PBDOT",     m_pb_dot, 0.);
    read(record, "A1",        m_a1);
    read(record, "XDOT",      m_x_dot, 0.);
    read(record, "ECC",       m_ecc);
    read(record, "ECCDOT",    m_ecc_dot, 0.);
    read(record, "OM",        m_om);
    read(record, "OMDOT",     m_om_dot, 0.);
    double dbl_t0 = 0.;
    read(record, "T0",        dbl_t0);
    read(record, "GAMMA",     m_gamma, 0.);

    // Create an AbsoluteTime object from the value of "T0" column.
    m_t0 = AbsoluteTime("TDB", Mjd1(dbl_t0));

    // Adjust units.
    m_om *= s_rad_per_deg;
    m_om_dot *= s_rad_year_per_deg_sec;
  }

  void BtModelEph::writeModelParameter(st_stream::OStream & os) const {
    os << format("PB",        m_pb,      "s")         << std::endl;
    os << format("PBDOT",     m_pb_dot,  "")          << std::endl;
    os << format("A1",        m_a1,      "lt-s")      << std::endl;
    os << format("XDOT",      m_x_dot,   "lt-s/s")    << std::endl;
    os << format("ECC",       m_ecc,     "")          << std::endl;
    os << format("ECCDOT",    m_ecc_dot, "s**(-1)")   << std::endl;
    os << format("OM",        m_om,      "radians")   << std::endl;
    os << format("OMDOT",     m_om_dot,  "radians/s") << std::endl;
    std::string t0_string;
    try {
      t0_string = m_t0.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      t0_string = m_t0.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("T0",        t0_string, "")   << std::endl;
    os << format("GAMMA",     m_gamma,   "s");
  }

  OrbitalEph * BtModelEph::clone() const { return new BtModelEph(*this); }

  double BtModelEph::calcOrbitalPhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Compute elapsed time from epoch of periastron in seconds.
    double delta_second = calcElapsedSecond(ev_time);

    // Compute orbital period.
    // --- Eq. 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double porb = m_pb + m_pb_dot * delta_second / 2.0;

    // Compute the complete phase.
    // --- The right-hand side of Eq. 2.27 in Blandford & Teuolsky, ApJ 205, 580 (1976), divided by two pi
    double phase = delta_second / porb;

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  ElapsedTime BtModelEph::calcOrbitalDelay(const AbsoluteTime & ev_time) const {
    // Note: Use Eqs. 2.27, 2.30, and 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976), because the forward
    //       conversion (from a pulsar proper time to an infinite-frequency barycenter arrival time) is needed
    //       in this method, and Eqs. 5 and 6 in Taylor & Weisberg, ApJ 345, 434 (1989) form the inverse formula
    //       which computes a pulsar proper time from an infinite-frequency barycenter arrival time.

    // Compute elapsed time from epoch of periastron in seconds.
    double delta_second = calcElapsedSecond(ev_time);

    // Compute orbital period.
    // --- Eq. 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double porb = m_pb + m_pb_dot * delta_second / 2.0;

    // Compute projected semimajor axis.
    // --- Eq. 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double semiax = m_a1 + m_x_dot * delta_second;

    // Compute eccentricity.
    // --- Eq. 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double eccen = m_ecc + m_ecc_dot * delta_second;

    // Compute periastron longitude.
    // --- Eq. 2.38 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double omega = m_om + m_om_dot * delta_second;

    // Calculate mean anomaly.
    // --- The right-hand side of Eq. 2.27 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double mean_anomaly = BtModelEph::s_two_pi * delta_second / porb;

    // Solve Kepler's equasion.
    // --- Eq. 2.27 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double eccen_anomaly = 0.0; // eccentric anomaly
    int status = atKepler(mean_anomaly, eccen, &eccen_anomaly);

    // atKepler not converged.
    if (0 != status) {
      throw std::runtime_error("Could not solve Kepler equation numerically (atKepler did not converge)");
    }

    // Compute time delays due to orbital motion.
    // --- Eqs. 2.30 and 2.31 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    double alpha = semiax * std::sin(omega);
    double beta = semiax * std::cos(omega) * std::sqrt(1. - eccen * eccen);
    double delay = alpha * (std::cos(eccen_anomaly) - eccen) + (beta + m_gamma) * std::sin(eccen_anomaly);

    // Return total delay.
    // --- Eq. 2.30 in Blandford & Teuolsky, ApJ 205, 580 (1976)
    return ElapsedTime(m_system->getName(), Duration(delay, "Sec"));
  }

}
