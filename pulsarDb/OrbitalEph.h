/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include "st_stream/Stream.h"
#include "tip/Table.h"

#include "pulsarDb/FormattedEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"

#include <iostream>
#include <vector>

namespace tip {
  class Header;
}

namespace pulsarDb {

  class OrbitalEph: public FormattedEph {
    public:
      void modulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      virtual const timeSystem::AbsoluteTime & t0() const = 0;

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /// \brief Create a copy of this object.
      virtual OrbitalEph * clone() const = 0;

      /// \brief Output text expression of this PulsarEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      OrbitalEph(const timeSystem::ElapsedTime & tolerance, int max_iteration);

      /// \brief Output text expression of subclass-specific parameters of this PulsarEph to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const = 0;

    private:
      timeSystem::ElapsedTime m_tolerance;
      int m_max_iteration;
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
        double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s);

      SimpleDdEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      virtual ~SimpleDdEph();

      virtual const timeSystem::AbsoluteTime & t0() const { return m_t0; }

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const;

      virtual OrbitalEph * clone() const;

    protected:
      /// \brief Output text expression of subclass-specific parameters of this PulsarEph to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const;

      const timeSystem::TimeSystem * m_system;
      double m_pb, m_pb_dot;
      double m_a1, m_x_dot;
      double m_ecc, m_ecc_dot;
      double m_om, m_om_dot;
      timeSystem::AbsoluteTime m_t0;
      double m_gamma;
      double m_shapiro_r, m_shapiro_s;
  };
}

#endif
