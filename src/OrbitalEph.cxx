#include <cmath>
#include <iomanip>

#include "pulsarDb/OrbitalEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeRep.h"

using namespace timeSystem;

namespace pulsarDb {

  const long double OrbitalEph::s_one_pi = M_PI;
  const long double OrbitalEph::s_two_pi = 2.L * OrbitalEph::s_one_pi;
  const long double OrbitalEph::s_rad_per_deg  = OrbitalEph::s_one_pi / 180.0L;
  const long double OrbitalEph::s_sec_per_day  = 86400.0L;
  const long double OrbitalEph::s_sec_per_year = 365.0L * OrbitalEph::s_sec_per_day;
  const long double OrbitalEph::s_rad_year_per_deg_sec = OrbitalEph::s_rad_per_deg / OrbitalEph::s_sec_per_year;
  const long double OrbitalEph::s_sec_per_microsec = 1.e-6L;

  OrbitalEph::OrbitalEph(const std::string & time_system_name, double parameters[NUMBER_ORBITAL_PAR], double unit_time_sec):
    m_system(&timeSystem::TimeSystem::getSystem(time_system_name)), m_par(parameters, parameters + NUMBER_ORBITAL_PAR),
    m_t0(time_system_name, Duration(IntFracPair(parameters[T0]), Day), Duration(0, 0.)), m_unit_time(unit_time_sec) {
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
  }

  OrbitalEph::OrbitalEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot,
    double ecc, double ecc_dot, double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma,
    double shapiro_r, double shapiro_s, double unit_time_sec): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
    m_par(NUMBER_ORBITAL_PAR, 0.), m_t0(t0), m_unit_time(unit_time_sec) {
    m_par[PB] = pb;
    m_par[PBDOT] = pb_dot;
    m_par[A1] = a1;
    m_par[XDOT] = x_dot;
    m_par[ECC] = ecc;
    m_par[ECCDOT] = ecc_dot;
    m_par[OM] = om;
    m_par[OMDOT] = om_dot;
    // TODO: Can the following steps to set m_par[T0] be improved?
    MjdRep mjd_rep(time_system_name, 0, 0.);
    mjd_rep.setAbsoluteTime(t0);
    IntFracPair time_pair = mjd_rep.getValue();
    m_par[T0] = time_pair.getIntegerPart() + time_pair.getFractionalPart();
    m_par[GAMMA] = gamma;
    m_par[SHAPIRO_R] = shapiro_r;
    m_par[SHAPIRO_S] = shapiro_s;

    // Adjust units.
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
  }

  OrbitalEph::~OrbitalEph() {}

  double OrbitalEph::dt(const timeSystem::AbsoluteTime & at) const {
    IntFracPair numerator = (at - m_t0).computeElapsedTime(m_system->getName()).getTime().getValue(Sec);
    return (numerator.getIntegerPart() + numerator.getFractionalPart()) / m_unit_time;
  }

  st_stream::OStream & OrbitalEph::write(st_stream::OStream & os) const {
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
}
