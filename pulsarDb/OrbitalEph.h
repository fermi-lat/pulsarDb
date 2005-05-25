/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include <vector>

namespace pulsarDb {

  enum BinaryIndex { PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S };

  class OrbitalEph {
    public:
      typedef std::vector<double>::size_type size_type;

      static const double s_one_pi;
      static const double s_two_pi;
      static const double s_rad_per_deg;
      static const double s_sec_per_day;
      static const double s_sec_per_year;
      static const double s_rad_year_per_deg_sec;
      static const double s_sec_per_microsec;

      OrbitalEph(double parameters[12]): m_par(parameters, parameters + 12) {
        m_par[OM] *= s_rad_per_deg;
        m_par[OMDOT] *= s_rad_year_per_deg_sec;
        m_par[SHAPIRO_R] *= s_sec_per_microsec;
      }

      virtual ~OrbitalEph() {}

      const double & operator [](size_type index) const { return m_par[index]; }

      virtual OrbitalEph * clone() const { return new OrbitalEph(*this); }

    private:
      std::vector<double> m_par;
  };

}

#endif
