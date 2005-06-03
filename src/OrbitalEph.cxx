#include <cmath>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/CanonicalTime.h"
#include "pulsarDb/OrbitalEph.h"

namespace pulsarDb {

  const double OrbitalEph::s_one_pi = M_PI;
  const double OrbitalEph::s_two_pi = 2. * OrbitalEph::s_one_pi;
  const double OrbitalEph::s_rad_per_deg  = OrbitalEph::s_one_pi / 180.0;
  const double OrbitalEph::s_sec_per_day  = 86400.0;
  const double OrbitalEph::s_sec_per_year = 365.0 * OrbitalEph::s_sec_per_day;
  const double OrbitalEph::s_rad_year_per_deg_sec = OrbitalEph::s_rad_per_deg / OrbitalEph::s_sec_per_year;
  const double OrbitalEph::s_sec_per_microsec = 1.e-6;

  OrbitalEph::OrbitalEph(double parameters[NUMBER_ORBITAL_PAR]): m_par(parameters, parameters + NUMBER_ORBITAL_PAR),
    m_t0(new TdbTime(parameters[T0])) {
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
  }

  OrbitalEph::OrbitalEph(double pb, double pb_dot, double a1, double x_dot, double ecc, double ecc_dot, double om, double om_dot,
    const AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s): m_par(NUMBER_ORBITAL_PAR, 0.), m_t0(t0.clone()) {
    m_par[PB] = pb;
    m_par[PBDOT] = pb_dot;
    m_par[A1] = a1;
    m_par[XDOT] = x_dot;
    m_par[ECC] = ecc;
    m_par[ECCDOT] = ecc_dot;
    m_par[OM] = om;
    m_par[OMDOT] = om_dot;
    m_par[T0] = TdbTime(t0).mjd();
    m_par[GAMMA] = gamma;
    m_par[SHAPIRO_R] = shapiro_r;
    m_par[SHAPIRO_S] = shapiro_s;

    // Adjust units.
    m_par[OM] *= s_rad_per_deg;
    m_par[OMDOT] *= s_rad_year_per_deg_sec;
    m_par[SHAPIRO_R] *= s_sec_per_microsec;
  }

  OrbitalEph::~OrbitalEph() { delete m_t0; }

}
