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
          \param t0geo_int The integer portion of the geocentric pulse arrival time.
          \param t0geo_frac The fracitonal portion of the geocentric pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      PulsarEph(double valid_since, double valid_until, long epoch_int, double epoch_frac, long t0geo_int, double t0geo_frac,
        double f0, double f1, double f2);

      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch The epoch (time origin).
          \param t0geo The geocentric pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      PulsarEph(double valid_since, double valid_until, long double epoch, long double t0geo, double f0, double f1, double f2);

      virtual ~PulsarEph();

      long double m_epoch;
      long double m_t0;
      double m_since;
      double m_until;
      double m_f0;
      double m_f1;
      double m_f2;
  };

  typedef std::vector<PulsarEph> PulsarEphCont;

}

#endif
