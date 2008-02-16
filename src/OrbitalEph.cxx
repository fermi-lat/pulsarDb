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

  // TODO: Hide this in SimpleDdEph class.
  enum BinaryIndex { PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S, NUMBER_ORBITAL_PAR };

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
    m_system(&timeSystem::TimeSystem::getSystem(time_system_name)), m_par(NUMBER_ORBITAL_PAR, 0.), m_t0(t0) {
    m_par[PB] = pb;
    m_par[PBDOT] = pb_dot;
    m_par[A1] = a1;
    m_par[XDOT] = x_dot;
    m_par[ECC] = ecc;
    m_par[ECCDOT] = ecc_dot;
    m_par[OM] = om;
    m_par[OMDOT] = om_dot;
    MjdRep mjd_rep(time_system_name, 0, 0.);
    mjd_rep = t0;
    m_par[T0] = mjd_rep.getValue().getDouble();
    m_par[GAMMA] = gamma;
    m_par[SHAPIRO_R] = shapiro_r;
    m_par[SHAPIRO_S] = shapiro_s;

    // Adjust units.
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
  }

  SimpleDdEph::SimpleDdEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    OrbitalEph(timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, 10.e-9)), 100),
    m_system(&timeSystem::TimeSystem::getSystem("TDB")), m_par(NUMBER_ORBITAL_PAR, 0.),
    m_t0("TDB", Duration(0, 0.), Duration(0, 0.)) {
    m_par[PB] = get(record["PB"]);
    m_par[PBDOT] = get(record["PBDOT"]);
    m_par[A1] = get(record["A1"]);
    m_par[XDOT] = get(record["XDOT"]);
    m_par[ECC] = get(record["ECC"]);
    m_par[ECCDOT] = get(record["ECCDOT"]);
    m_par[OM] = get(record["OM"]);
    m_par[OMDOT] = get(record["OMDOT"]);
    m_par[T0] = get(record["T0"]);
    m_par[GAMMA] = get(record["GAMMA"]);
    m_par[SHAPIRO_R] = get(record["SHAPIRO_R"]);
    m_par[SHAPIRO_S] = get(record["SHAPIRO_S"]);

    // Handle any INDEFs.
    for (size_t index = 0; index != sizeof(m_par) / sizeof(double); ++index) {
      if (0 != IsNotANumber(m_par[index])) {
        switch (index) {
        case PB:
        case A1:
        case ECC:
        case OM:
        case T0:
          throw std::runtime_error("PulsarDb::getEph(): invalid orbital ephemeris");
          break;
        case PBDOT:
        case XDOT:
        case ECCDOT:
        case OMDOT:
        case GAMMA:
        case SHAPIRO_R:
        case SHAPIRO_S:
        default:
          m_par[index] = 0.;
          break;
        }
      }
    }

    m_t0 = timeSystem::AbsoluteTime("TDB", Duration(IntFracPair(m_par[T0]), Day), Duration(0, 0.));
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
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
    os.prefix().width(14); os << "PB = " << m_par[PB] << std::endl;
    os.prefix().width(14); os << "PBDOT = " << m_par[PBDOT] << std::endl;
    os.prefix().width(14); os << "A1 = " << m_par[A1] << std::endl;
    os.prefix().width(14); os << "XDOT = " << m_par[XDOT] << std::endl;
    os.prefix().width(14); os << "ECC = " << m_par[ECC] << std::endl;
    os.prefix().width(14); os << "ECCDOT = " << m_par[ECCDOT] << std::endl;
    os.prefix().width(14); os << "OM = " << m_par[OM] << std::endl;
    os.prefix().width(14); os << "OMDOT = " << m_par[OMDOT] << std::endl;
    os.prefix().width(14); os << "T0 = " << m_par[T0] << std::endl;
    os.prefix().width(14); os << "GAMMA = " << m_par[GAMMA] << std::endl;
    os.prefix().width(14); os << "SHAPIRO_R = " << m_par[SHAPIRO_R] << std::endl;
    os.prefix().width(14); os << "SHAPIRO_S = " << m_par[SHAPIRO_S];
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

  OrbitalEph * SimpleDdEph::clone() const { return new SimpleDdEph(*this); }

  double SimpleDdEph::calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    // compute elapsed time from epoch of periastron in seconds
    double delta_second = calcElapsedSecond(ev_time);

    // Compute the time difference as a fraction of the period.
    double delta_period = delta_second / m_par[PB];

    // Compute the complete phase.
    double phase = delta_period * (1. - delta_period * m_par[PBDOT] / 2.0);

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
    double delta_period = delta_second / m_par[PB];
    double mean_anomaly = SimpleDdEph::s_two_pi * delta_period
      * (1. - delta_period * m_par[PBDOT] / 2.0);

    // solve Kepler's equasion
    double eccen = m_par[ECC] + m_par[ECCDOT] * delta_second; // eccentricity
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
    double omega = m_par[OM]
      + m_par[OMDOT] * true_anomaly * m_par[PB] / SimpleDdEph::s_two_pi;

    // compute projected semimajor axis
    double semiax = m_par[A1] + m_par[XDOT] * delta_second;

    // compute time delays due to orbital motion
    double roemer_frac = std::sin(omega) * (std::cos(eccen_anomaly) - eccen)
      + std::sqrt(1.0 - eccen*eccen) * std::cos(omega) * std::sin(eccen_anomaly);
    double roemer = semiax * roemer_frac;
    double einstein = m_par[GAMMA] * std::sin(eccen_anomaly);
    double shapiro = - 2.0 * m_par[SHAPIRO_R]
      * std::log(1.0 - eccen*std::cos(eccen_anomaly) - m_par[SHAPIRO_S]*roemer_frac);

    // return total delay
    return ElapsedTime(m_system->getName(), Duration(0, roemer + einstein + shapiro));
  }

}
