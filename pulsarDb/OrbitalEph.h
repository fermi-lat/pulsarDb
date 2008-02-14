/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include "st_stream/Stream.h"
#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"

#include <iostream>
#include <vector>

namespace pulsarDb {

  class OrbitalEph {
    public:
      OrbitalEph(const timeSystem::ElapsedTime & tolerance, int max_iteration);

      void modulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      virtual const timeSystem::AbsoluteTime & t0() const = 0;

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const = 0;

      virtual OrbitalEph * clone() const = 0;

      virtual st_stream::OStream & write(st_stream::OStream & os) const = 0;

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

      SimpleDdEph(const std::string & time_system_name, const tip::Table::ConstRecord & record);

      virtual ~SimpleDdEph();

      virtual const timeSystem::AbsoluteTime & t0() const { return m_t0; }

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const;

      virtual OrbitalEph * clone() const;

      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    private:
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const;

      // TODO: Avoid duplication of these get methods (another copy is in PulsarEph.h).
      /** \brief Helper method to get a value from a cell, returning it as a double, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
      */
      inline double get(const tip::TableCell & cell) const {
        double value = 0.;
        get(cell, value);
        return value;
      }

      // TODO: Avoid duplication of these get methods (another copy is in PulsarEph.h).
      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
          \param value Variable to store the value.
      */
      template <typename T>
      inline void get(const tip::TableCell & cell, T & value) const {
        // WARNING: This will break for a string column.
        if (cell.isNull()) value = 0;
        else cell.get(value);
      }

      const timeSystem::TimeSystem * m_system;
      std::vector<double> m_par;
      timeSystem::AbsoluteTime m_t0;
  };
}

#endif
