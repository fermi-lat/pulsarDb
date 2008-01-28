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

  // TODO: Hide this in SimpleDdEph class.
  enum BinaryIndex { PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S, NUMBER_ORBITAL_PAR };

  class OrbitalEph {
    public:
      void modulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      virtual const timeSystem::AbsoluteTime & t0() const = 0;

      // TODO: Remove OrbitalEph::dt method.
      virtual double dt(const timeSystem::AbsoluteTime & at) const = 0;

      virtual const timeSystem::TimeSystem & getSystem() const = 0;

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const = 0;

      virtual OrbitalEph * clone() const = 0;

      virtual st_stream::OStream & write(st_stream::OStream & os) const = 0;
  };

  typedef std::vector<OrbitalEph *> OrbitalEphCont;

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const OrbitalEph & eph) { return eph.write(os); }

  class SimpleDdEph : public OrbitalEph {
    public:
      typedef std::vector<double>::size_type size_type;

      static const double s_one_pi;
      static const double s_two_pi;
      static const double s_rad_per_deg;
      static const double s_sec_per_day;
      static const double s_sec_per_year;
      static const double s_rad_year_per_deg_sec;
      static const double s_sec_per_microsec;

      SimpleDdEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot, double ecc, double ecc_dot,
        double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s,
        double unit_time_sec = 1.);

      SimpleDdEph(const std::string & time_system_name, double parameters[NUMBER_ORBITAL_PAR], double unit_time_sec = 1.);

      virtual ~SimpleDdEph();

      virtual const timeSystem::AbsoluteTime & t0() const { return m_t0; }

      // TODO: Make OrbitalEph::dt method private.
      virtual double dt(const timeSystem::AbsoluteTime & at) const;

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const;

      virtual OrbitalEph * clone() const;

      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    private:
      const timeSystem::TimeSystem * m_system;
      std::vector<double> m_par;
      timeSystem::AbsoluteTime m_t0;
      double m_unit_time;
  };
}

#endif
