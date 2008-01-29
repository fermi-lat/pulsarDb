/** \file EphChooser.h
    \brief Interface for EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphChooser_h
#define pulsarDb_EphChooser_h

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"

namespace timeSystem {
  class AbsoluteTime;
}

namespace pulsarDb {

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
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const = 0;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const = 0;

      virtual EphChooser * clone() const = 0;

      const PulsarEph & findClosest(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

      const OrbitalEph & findClosest(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

    protected:
      double deltaTime(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
  };

  /** \class StrictEphChooser
      \brief Abstraction providing method for choosing best ephemeris from a container of candidate.
  */
  class StrictEphChooser : public EphChooser {
    public:
      StrictEphChooser(const timeSystem::ElapsedTime & tolerance = timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, 1.e-6)));

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

      virtual EphChooser * clone() const;

    protected:
      timeSystem::ElapsedTime m_tolerance;
  };

  /** \class SloppyEphChooser
      \brief Abstraction providing method for choosing best ephemeris from a container of candidate, being
      less pedantic than EphChooser.
  */
  class SloppyEphChooser : public EphChooser {
    public:
      SloppyEphChooser();

      SloppyEphChooser(const timeSystem::ElapsedTime & tolerance);

      /** \brief Choose the ephemeris closest to the given time even if the time is outside its interval of validity.
                 If no ephemeris meets these requirements, an exception is thrown.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const;

      virtual EphChooser * clone() const;

    protected:
      StrictEphChooser m_strict_chooser;
  };

}
#endif
