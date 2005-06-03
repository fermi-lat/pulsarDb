/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include <vector>

namespace pulsarDb {

  class AbsoluteTime;

  enum BinaryIndex { PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S, NUMBER_ORBITAL_PAR };

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

      OrbitalEph(double pb, double pb_dot, double a1, double x_dot, double ecc, double ecc_dot, double om, double om_dot,
        const AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s);

      OrbitalEph(double parameters[NUMBER_ORBITAL_PAR]);

      virtual ~OrbitalEph();

      const double & operator [](size_type index) const { return m_par[index]; }

      virtual const AbsoluteTime & t0() const { return *m_t0; }

      virtual OrbitalEph * clone() const { return new OrbitalEph(*this); }

    private:
      std::vector<double> m_par;
      AbsoluteTime * m_t0;
  };

  typedef std::vector<OrbitalEph *> OrbitalEphCont;
}

#endif
