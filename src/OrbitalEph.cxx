#include <cmath>
#include "pulsarDb/OrbitalEph.h"

namespace pulsarDb {
  const double OrbitalEph::s_one_pi = M_PI;
  const double OrbitalEph::s_two_pi = 2. * OrbitalEph::s_one_pi;
  const double OrbitalEph::s_rad_per_deg  = OrbitalEph::s_one_pi / 180.0;
  const double OrbitalEph::s_sec_per_day  = 86400.0;
  const double OrbitalEph::s_sec_per_year = 365.0 * OrbitalEph::s_sec_per_day;
  const double OrbitalEph::s_rad_year_per_deg_sec = OrbitalEph::s_rad_per_deg / OrbitalEph::s_sec_per_year;
  const double OrbitalEph::s_sec_per_microsec = 1.e-6;
}
