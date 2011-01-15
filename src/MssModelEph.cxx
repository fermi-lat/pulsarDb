/** \file MssModelEph.cxx
    \brief Implementation of the MssModelEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <limits>
#include <iomanip>
#include <stdexcept>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/MssModelEph.h"

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

  const double MssModelEph::s_one_pi = M_PI;
  const double MssModelEph::s_two_pi = 2. * MssModelEph::s_one_pi;
  const double MssModelEph::s_rad_per_deg  = MssModelEph::s_one_pi / 180.;
  const double MssModelEph::s_sec_per_day  = 86400.;
  const double MssModelEph::s_sec_per_year = 365.25 * MssModelEph::s_sec_per_day;
  const double MssModelEph::s_rad_year_per_deg_sec = MssModelEph::s_rad_per_deg / MssModelEph::s_sec_per_year;
  const double MssModelEph::s_rad_year2_per_deg_sec2 = MssModelEph::s_rad_year_per_deg_sec / MssModelEph::s_sec_per_year;
  const double MssModelEph::s_sec_per_microsec = 1.e-6;

  MssModelEph::MssModelEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot, double x_2dot,
    double ecc, double ecc_dot, double om, double om_dot, double om_2dot, const AbsoluteTime & t0, double delta_r, double delta_theta,
    double gamma, double shapiro_r, double shapiro_s, double aberration_a, double aberration_b):
    OrbitalEph(ElapsedTime(time_system_name, Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem(time_system_name)),
    m_pb(pb), m_pb_dot(pb_dot), m_a1(a1), m_x_dot(x_dot), m_x_2dot(x_2dot), m_ecc(ecc), m_ecc_dot(ecc_dot),
    m_om(om * s_rad_per_deg), m_om_dot(om_dot * s_rad_year_per_deg_sec), m_om_2dot(om_2dot * s_rad_year2_per_deg_sec2),
    m_t0(t0), m_delta_r(delta_r), m_delta_theta(delta_theta), m_gamma(gamma), m_shapiro_r(shapiro_r * s_sec_per_microsec),
    m_shapiro_s(shapiro_s), m_aberration_a(aberration_a), m_aberration_b(aberration_b) {}

  MssModelEph::MssModelEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    OrbitalEph(ElapsedTime("TDB", Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem("TDB")), m_t0("TDB", 0, 0.) {
    // Get parameters from record.
    // Required fields: PB, A1, ECC, OM, T0.
    // Optional fields: PBDOT, XDOT, X2DOT, ECCDOT, OMDOT, OM2DOT, DELTA_R, DELTA_THETA, GAMMA, SHAPIRO_R, SHAPIRO_S,
    //                  ABERRATION_A, ABERRATION_B.
    read(record, "PB",           m_pb);
    read(record, "PBDOT",        m_pb_dot, 0.);
    read(record, "A1",           m_a1);
    read(record, "XDOT",         m_x_dot, 0.);
    read(record, "X2DOT",        m_x_2dot, 0.);
    read(record, "ECC",          m_ecc);
    read(record, "ECCDOT",       m_ecc_dot, 0.);
    read(record, "OM",           m_om);
    read(record, "OMDOT",        m_om_dot, 0.);
    read(record, "OM2DOT",       m_om_2dot, 0.);
    double dbl_t0 = 0.;
    read(record, "T0",           dbl_t0);
    read(record, "GAMMA",        m_gamma, 0.);
    read(record, "DELTA_R",      m_delta_r, 0.);
    read(record, "DELTA_THETA",  m_delta_theta, 0.);
    read(record, "SHAPIRO_R",    m_shapiro_r, 0.);
    read(record, "SHAPIRO_S",    m_shapiro_s, 0.);
    read(record, "ABERRATION_A", m_aberration_a, 0.);
    read(record, "ABERRATION_B", m_aberration_b, 0.);

    // Create an AbsoluteTime object from the value of "T0" column.
    m_t0 = AbsoluteTime("TDB", Mjd1(dbl_t0));

    // Adjust units.
    m_om *= s_rad_per_deg;
    m_om_dot *= s_rad_year_per_deg_sec;
    m_om_2dot *= s_rad_year2_per_deg_sec2;
    m_shapiro_r *= s_sec_per_microsec;
  }

  void MssModelEph::writeModelParameter(st_stream::OStream & os) const {
    os << format("PB",           m_pb,           "s")            << std::endl;
    os << format("PBDOT",        m_pb_dot,       "")             << std::endl;
    os << format("A1",           m_a1,           "lt-s")         << std::endl;
    os << format("XDOT",         m_x_dot,        "lt-s/s")       << std::endl;
    os << format("X2DOT",        m_x_2dot,       "lt-s/s**2")    << std::endl;
    os << format("ECC",          m_ecc,          "")             << std::endl;
    os << format("ECCDOT",       m_ecc_dot,      "s**(-1)")      << std::endl;
    os << format("OM",           m_om,           "radians")      << std::endl;
    os << format("OMDOT",        m_om_dot,       "radians/s")    << std::endl;
    os << format("OM2DOT",       m_om_2dot,      "radians/s**2") << std::endl;
    std::string t0_string;
    try {
      t0_string = m_t0.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      t0_string = m_t0.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("T0",           t0_string,      "")  << std::endl;
    os << format("DELTA_R",      m_delta_r,      "")  << std::endl;
    os << format("DELTA_THETA",  m_delta_theta,  "")  << std::endl;
    os << format("GAMMA",        m_gamma,        "s") << std::endl;
    os << format("SHAPIRO_R",    m_shapiro_r,    "s") << std::endl;
    os << format("SHAPIRO_S",    m_shapiro_s,    "")  << std::endl;
    os << format("ABERRATION_A", m_aberration_a, "s") << std::endl;
    os << format("ABERRATION_B", m_aberration_b, "s") << std::endl;
  }

  OrbitalEph * MssModelEph::clone() const { return new MssModelEph(*this); }

  double MssModelEph::calcOrbitalPhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Compute elapsed time from epoch of periastron in seconds.
    double delta_second = calcElapsedSecond(ev_time);

    // Compute the time difference as a fraction of the period.
    // --- The term that appears twice in Eq. 12 in Taylor & Weisberg, ApJ 345, 434 (1989), divided by two pi.
    double delta_period = delta_second / m_pb;

    // Compute the complete phase.
    // --- The right-hand side of Eq. 12 in Taylor & Weisberg, ApJ 345, 434 (1989), divided by two pi.
    double phase = delta_period * (1. - delta_period * m_pb_dot / 2.0);

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  ElapsedTime MssModelEph::calcOrbitalDelay(const AbsoluteTime & ev_time) const {
    // Compute elapsed time from epoch of periastron in seconds.
    double delta_second = calcElapsedSecond(ev_time);

    // Calculate mean anomaly.
    // --- The right-hand side of Eq. 12 in Taylor & Weisberg, ApJ 345, 434 (1989)
    double delta_period = delta_second / m_pb;
    double mean_anomaly = MssModelEph::s_two_pi * delta_period
      * (1. - delta_period * m_pb_dot / 2.0);

    // Solve Kepler's equasion.
    // --- Eq. 12 in Taylor & Weisberg, ApJ 345, 434 (1989)
    double eccen = m_ecc + m_ecc_dot * delta_second; // eccentricity
    double eccen_anomaly = 0.0; // eccentric anomaly
    int status = atKepler(mean_anomaly, eccen, &eccen_anomaly);

    // atKepler not converged.
    if (0 != status) {
      throw std::runtime_error("Could not solve Kepler equation numerically (atKepler did not converge)");
    }

    // Convert eccentric anomaly to true anomaly.
    // --- Eq. 13 in Taylor & Weisberg, ApJ 345, 434 (1989)
    double true_anomaly = 2.0 * std::atan(std::sqrt((1.0+eccen)/(1.0-eccen))
        * std::tan(eccen_anomaly/2.0));
    true_anomaly += MssModelEph::s_two_pi * floor((eccen_anomaly - true_anomaly)/ MssModelEph::s_two_pi);
    while ((true_anomaly - eccen_anomaly) > MssModelEph::s_one_pi) true_anomaly -= MssModelEph::s_two_pi;
    while ((eccen_anomaly - true_anomaly) > MssModelEph::s_one_pi) true_anomaly += MssModelEph::s_two_pi;

    // Compute projected semimajor axis.
    // --- Eqs. 64 and 73 in Wex, MNRAS 298, 67 (1998)
    double delta_second_squared = delta_second * delta_second;
    double semiax = m_a1 + m_x_dot * true_anomaly * m_pb / MssModelEph::s_two_pi + m_x_2dot / 2.0 * delta_second_squared;

    // Compute periastron longitude.
    // --- Eqs. 65 and 74 in Wex, MNRAS 298, 67 (1998)
    double omega = m_om + m_om_dot * true_anomaly * m_pb / MssModelEph::s_two_pi + m_om_2dot / 2.0 * delta_second_squared;

    // Compute time delays due to orbital motion.
    // --- Eqs. 8, 9, and 10 in Taylor & Weisberg, ApJ 345, 434 (1989)
    double eccen_r = eccen * (1.0 + m_delta_r);
    double eccen_theta = eccen * (1.0 + m_delta_theta);
    double roemer = semiax * (std::sin(omega) * (std::cos(eccen_anomaly) - eccen_r)
      + std::sqrt(1.0 - eccen_theta*eccen_theta) * std::cos(omega) * std::sin(eccen_anomaly));
    double einstein = m_gamma * std::sin(eccen_anomaly);
    double shapiro = - 2.0 * m_shapiro_r * std::log(1.0 - eccen*std::cos(eccen_anomaly)
      - m_shapiro_s * (std::sin(omega) * (std::cos(eccen_anomaly) - eccen)
      + std::sqrt(1.0 - eccen*eccen) * std::cos(omega) * std::sin(eccen_anomaly)));
    double aberration = m_aberration_a * (std::sin(omega + true_anomaly) + eccen * std::sin(omega))
      + m_aberration_b * (std::cos(omega + true_anomaly) + eccen * std::cos(omega));

    // Return total delay.
    // --- Eq. 7 in Taylor & Weisberg, ApJ 345, 434 (1989)
    return ElapsedTime(m_system->getName(), Duration(roemer + einstein + shapiro + aberration, "Sec"));
  }

}
