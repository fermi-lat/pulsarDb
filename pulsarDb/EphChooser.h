/** \file EphChooser.h
    \brief Interface for EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphChooser_h
#define pulsarDb_EphChooser_h

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace pulsarDb {

  class AbsoluteTime;

  /** \class EphChooser
      \brief Abstraction providing method for choosing best ephemeris from a container of candidate.
  */
  class EphChooser {
    public:
      virtual ~EphChooser() {}

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const AbsoluteTime & t) const;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const AbsoluteTime & t) const;
  };

  /** \class SloppyEphChooser
      \brief Abstraction providing method for choosing best ephemeris from a container of candidate, being
      less pedantic than EphChooser.
  */
  class SloppyEphChooser : public EphChooser {
    public:
      /** \brief Choose the ephemeris closest to the given time even if the time is outside its interval of validity.
                 If no ephemeris meets these requirements, an exception is thrown.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const AbsoluteTime & t) const;
  };

}
#endif
