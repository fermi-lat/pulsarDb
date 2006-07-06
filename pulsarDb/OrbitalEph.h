/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include "st_stream/Stream.h"

#include "timeSystem/AbsoluteTime.h"

#include <iostream>
#include <vector>

namespace pulsarDb {

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

      OrbitalEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot, double ecc, double ecc_dot,
        double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s,
        double unit_time_sec = 1.);

      OrbitalEph(const std::string & time_system_name, double parameters[NUMBER_ORBITAL_PAR], double unit_time_sec = 1.);

      virtual ~OrbitalEph();

      const double & operator [](size_type index) const { return m_par[index]; }

      virtual const timeSystem::AbsoluteTime & t0() const { return m_t0; }

      virtual double dt(const timeSystem::AbsoluteTime & at) const;

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual OrbitalEph * clone() const { return new OrbitalEph(*this); }

      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    private:
      const timeSystem::TimeSystem * m_system;
      std::vector<double> m_par;
      timeSystem::AbsoluteTime m_t0;
      double m_unit_time;
  };

  typedef std::vector<OrbitalEph *> OrbitalEphCont;

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const OrbitalEph & eph) { return eph.write(os); }

}

#endif
