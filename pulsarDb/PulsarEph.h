/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <vector>

namespace pulsarDb {

  /** \struct PulsarEph
      \brief Class representing a single pulsar ephemeris.
  */
  struct PulsarEph {
      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch_int The integer portion of the epoch (time origin).
          \param epoch_frac The fractional portion of the epoch (time origin).
          \param toa_geo_int The integer portion of the geocentric pulse arrival time.
          \param toa_geo_frac The fractional portion of the geocentric pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      PulsarEph(double valid_since, double valid_until, long epoch_int, double epoch_frac, long toa_geo_int, double toa_geo_frac,
        double phi0, double f0, double f1, double f2);

      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch The epoch (time origin).
          \param toa_geo The geocentric pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      PulsarEph(double valid_since, double valid_until, long double epoch, long double toa_geo, double phi0, double f0,
        double f1, double f2);

      virtual ~PulsarEph();

      long double m_epoch;
      long double m_toa; // Deprecated. TOA will be used to compute phi0
      double m_since; // Deprecated. Used only to select ephemeris.
      double m_until; // Deprecated. Used only to select ephemeris.
      double m_phi0;
      double m_f0;
      double m_f1;
      double m_f2;
  };

  typedef std::vector<PulsarEph> PulsarEphCont;

}

#endif
