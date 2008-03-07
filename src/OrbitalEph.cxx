#include <cmath>
#include <iomanip>

#include "pulsarDb/OrbitalEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"

using namespace timeSystem;

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

  // TODO: Avoid duplication of these get methods (another copy is in PulsarEph.h).
  inline bool IsNotANumber(double x) {
#ifdef WIN32
    return 0 != _isnan(x);
#else
    return 0 != std::isnan(x);
#endif
  }

}

namespace pulsarDb {

  OrbitalEph::OrbitalEph(const timeSystem::ElapsedTime & tolerance, int max_iteration): m_tolerance(tolerance),
    m_max_iteration(max_iteration) {}

  void OrbitalEph::modulateBinary(timeSystem::AbsoluteTime & ev_time) const {
    ev_time += calcOrbitalDelay(ev_time);
  }

  void OrbitalEph::demodulateBinary(timeSystem::AbsoluteTime & ev_time) const {
    // Save target arrival time (ev_time) in orig_time.
    AbsoluteTime orig_time = ev_time;

    // Initial guess of orbital delay.
    ElapsedTime delay = calcOrbitalDelay(ev_time);

    // Iterative approximation of demodulated time.
    int ii;
    for (ii=0; ii<m_max_iteration; ii++) {

      // Compute next candidate of demodulated time.
      ev_time = orig_time - delay;

      // Compute orbital delay at ev_time.
      delay = calcOrbitalDelay(ev_time);

      // Compare time difference between candidate demodulated time
      // (ev_time) and target arrival time (orig_time) with the
      // estimated orbital delay based on the binary model (delay).
      if (orig_time.equivalentTo(ev_time + delay, m_tolerance)) break;
    }

    // Check for non-convergence.
    if (ii == m_max_iteration) throw std::runtime_error("Binary demodulation did not converge.");

  }

  const double SimpleDdEph::s_one_pi = M_PI;
  const double SimpleDdEph::s_two_pi = 2. * SimpleDdEph::s_one_pi;
  const double SimpleDdEph::s_rad_per_deg  = SimpleDdEph::s_one_pi / 180.;
  const double SimpleDdEph::s_sec_per_day  = 86400.;
  const double SimpleDdEph::s_sec_per_year = 365. * SimpleDdEph::s_sec_per_day;
  const double SimpleDdEph::s_rad_year_per_deg_sec = SimpleDdEph::s_rad_per_deg / SimpleDdEph::s_sec_per_year;
  const double SimpleDdEph::s_sec_per_microsec = 1.e-6;

  SimpleDdEph::SimpleDdEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot,
    double ecc, double ecc_dot, double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma,
    double shapiro_r, double shapiro_s):
    OrbitalEph(timeSystem::ElapsedTime(time_system_name, timeSystem::Duration(0, 10.e-9)), 100),
    m_system(&timeSystem::TimeSystem::getSystem(time_system_name)), m_pb(pb), m_pb_dot(pb_dot), m_a1(a1), m_x_dot(x_dot),
    m_ecc(ecc), m_ecc_dot(ecc_dot), m_om(om * s_rad_per_deg), m_om_dot(om_dot * s_rad_year_per_deg_sec), m_t0(t0), m_gamma(gamma),
    m_shapiro_r(shapiro_r * s_sec_per_microsec), m_shapiro_s(shapiro_s) {}

  SimpleDdEph::SimpleDdEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    OrbitalEph(timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, 10.e-9)), 100),
    m_system(&timeSystem::TimeSystem::getSystem("TDB")), m_t0("TDB", Duration(0, 0.), Duration(0, 0.)) {
    m_pb = get(record["PB"]);
    m_pb_dot = get(record["PBDOT"]);
    m_a1 = get(record["A1"]);
    m_x_dot = get(record["XDOT"]);
    m_ecc = get(record["ECC"]);
    m_ecc_dot = get(record["ECCDOT"]);
    m_om = get(record["OM"]) * s_rad_per_deg;
    m_om_dot = get(record["OMDOT"]) * s_rad_year_per_deg_sec;
    double dbl_t0 = get(record["T0"]);
    m_gamma = get(record["GAMMA"]);
    m_shapiro_r = get(record["SHAPIRO_R"]) * s_sec_per_microsec;
    m_shapiro_s = get(record["SHAPIRO_S"]);

    // Handle any INDEFs.
    if (IsNotANumber(m_pb) || IsNotANumber(m_a1) || IsNotANumber(m_ecc) || IsNotANumber(m_om) || IsNotANumber(dbl_t0)) {
      throw std::runtime_error("SimpleDdEph: invalid orbital ephemeris passed to the constructor");
    }
    if (IsNotANumber(m_pb_dot)) m_pb_dot = 0.;
    if (IsNotANumber(m_x_dot)) m_x_dot = 0.;
    if (IsNotANumber(m_ecc_dot)) m_ecc_dot = 0.;
    if (IsNotANumber(m_om_dot)) m_om_dot = 0.;
    if (IsNotANumber(m_gamma)) m_gamma = 0.;
    if (IsNotANumber(m_shapiro_r)) m_shapiro_r = 0.;
    if (IsNotANumber(m_shapiro_s)) m_shapiro_s = 0.;

    // Create an AbsoluteTime object from the value of "T0" column.
    m_t0 = timeSystem::AbsoluteTime("TDB", Duration(IntFracPair(dbl_t0), Day), Duration(0, 0.));
  }

  SimpleDdEph::~SimpleDdEph() {}

  double SimpleDdEph::calcElapsedSecond(const timeSystem::AbsoluteTime & at) const {
    return (at - m_t0).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  // TODO: Can this part be unified to PulsarEph::write?
  st_stream::OStream & SimpleDdEph::write(st_stream::OStream & os) const {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14); os << "PB = " << m_pb << std::endl;
    os.prefix().width(14); os << "PBDOT = " << m_pb_dot << std::endl;
    os.prefix().width(14); os << "A1 = " << m_a1 << std::endl;
    os.prefix().width(14); os << "XDOT = " << m_x_dot << std::endl;
    os.prefix().width(14); os << "ECC = " << m_ecc << std::endl;
    os.prefix().width(14); os << "ECCDOT = " << m_ecc_dot << std::endl;
    os.prefix().width(14); os << "OM = " << m_om << std::endl;
    os.prefix().width(14); os << "OMDOT = " << m_om_dot << std::endl;
    MjdRep mjd_rep(m_system->getName(), 0, 0.);
    mjd_rep = m_t0;
    double dbl_t0 = mjd_rep.getValue().getDouble();
    os.prefix().width(14); os << "T0 = " << dbl_t0 << std::endl;
    os.prefix().width(14); os << "GAMMA = " << m_gamma << std::endl;
    os.prefix().width(14); os << "SHAPIRO_R = " << m_shapiro_r << std::endl;
    os.prefix().width(14); os << "SHAPIRO_S = " << m_shapiro_s;
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

  OrbitalEph * SimpleDdEph::clone() const { return new SimpleDdEph(*this); }

  double SimpleDdEph::calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    // compute elapsed time from epoch of periastron in seconds
    double delta_second = calcElapsedSecond(ev_time);

    // Compute the time difference as a fraction of the period.
    double delta_period = delta_second / m_pb;

    // Compute the complete phase.
    double phase = delta_period * (1. - delta_period * m_pb_dot / 2.0);

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    double int_part; // ignored, needed for modf.
    phase = std::modf(phase_offset + phase, &int_part);
    if (phase < 0.) ++phase;
    return phase;
  }

  timeSystem::ElapsedTime SimpleDdEph::calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const {
    // compute elapsed time from epoch of periastron in seconds
    double delta_second = calcElapsedSecond(ev_time);

    // calculate mean anomaly
    double delta_period = delta_second / m_pb;
    double mean_anomaly = SimpleDdEph::s_two_pi * delta_period
      * (1. - delta_period * m_pb_dot / 2.0);

    // solve Kepler's equasion
    double eccen = m_ecc + m_ecc_dot * delta_second; // eccentricity
    double eccen_anomaly = 0.0; // eccentric anomaly
    int status = atKepler(mean_anomaly, eccen, &eccen_anomaly);

    // atKepler not converged
    if (0 != status) {
      throw std::runtime_error("atKepler did not converge.");
    }

    // convert eccentric anomaly to true anomaly
    double true_anomaly = 2.0 * std::atan(std::sqrt((1.0+eccen)/(1.0-eccen))
        * std::tan(eccen_anomaly/2.0));
    true_anomaly += SimpleDdEph::s_two_pi * floor((eccen_anomaly - true_anomaly)/ SimpleDdEph::s_two_pi);
    while ((true_anomaly - eccen_anomaly) > SimpleDdEph::s_one_pi) true_anomaly -= SimpleDdEph::s_two_pi;
    while ((eccen_anomaly - true_anomaly) > SimpleDdEph::s_one_pi) true_anomaly += SimpleDdEph::s_two_pi;

    // compute periastron longitude
    double omega = m_om
      + m_om_dot * true_anomaly * m_pb / SimpleDdEph::s_two_pi;

    // compute projected semimajor axis
    double semiax = m_a1 + m_x_dot * delta_second;

    // compute time delays due to orbital motion
    double roemer_frac = std::sin(omega) * (std::cos(eccen_anomaly) - eccen)
      + std::sqrt(1.0 - eccen*eccen) * std::cos(omega) * std::sin(eccen_anomaly);
    double roemer = semiax * roemer_frac;
    double einstein = m_gamma * std::sin(eccen_anomaly);
    double shapiro = - 2.0 * m_shapiro_r
      * std::log(1.0 - eccen*std::cos(eccen_anomaly) - m_shapiro_s*roemer_frac);

    // return total delay
    return ElapsedTime(m_system->getName(), Duration(0, roemer + einstein + shapiro));
  }

}
