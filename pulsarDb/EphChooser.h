/** \file EphChooser.h
    \brief Interface for EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphChooser_h
#define pulsarDb_EphChooser_h

#include "pulsarDb/EphStatus.h"
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
          \param ephemerides Container of pulsar ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const = 0;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param ephemerides Container of orbital ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const = 0;

      /** \brief Examine time coverage by a given set of pulsar ephemerides and ephemeris remarks in the given time interval,
                 and return the result as a list of ephemeris remarks of various types, such as ephemeris gaps, glitches,
                 timing noise, and informative notes. Interpretation of "time coverage" may depend on each subclass.
                 One typical behavior is to define an ephemeris gap as a time interval during which its choose method will
                 throw an exception.
          \param ephemerides Container of pulsar ephemerides whose time coverage is to be examined.
          \param start_time The start time of a time interval of interest.
          \param stop_time The stop time of a time interval of interest.
          \param eph_status Container of ephemeris status that stores the result of examination.
      */
      virtual void examine(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & start_time,
        const timeSystem::AbsoluteTime & stop_time, EphStatusCont & eph_status) const = 0;

      /// \brief Returns an object of the same type as this object.
      virtual EphChooser * clone() const = 0;

    protected:
      /** \brief Find a pulsar ephemeris whose reference epoch is the closest to the given absolute time.
          \param ephemerides Container of pulsar ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      const PulsarEph & findClosest(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Find an orbital ephemeris whose reference epoch is the closest to the given absolute time.
          \param ephemerides Container of orbital ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      const OrbitalEph & findClosest(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Compute a time separation between two absolute times, to be used to compare more than one ephemerides
                 to find a best candidate. Return value of this method is non-negative.
          \param at1 One absolute time.
          \param at2 The other absolute time.
      */
      double measureTimeSeparation(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
  };

  /** \class StrictEphChooser
      \brief Abstraction providing method for choosing best ephemeris from a container of candidate.
  */
  class StrictEphChooser : public EphChooser {
    public:
      StrictEphChooser(const timeSystem::ElapsedTime & tolerance = timeSystem::ElapsedTime("TDB", timeSystem::Duration(1.e-6, "Sec")));

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param ephemerides Container of pulsar ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param ephemerides Container of orbital ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Return ephemeris gaps and ephemeris remarks overlapped with a given time interval.
          \param ephemerides Container of pulsar ephemerides whose time coverage is to be examined.
          \param start_time The start time of a time interval of interest.
          \param stop_time The stop time of a time interval of interest.
          \param eph_status Container of ephemeris status that stores the result of examination. The orginal contents
                 of the container will be removed by this method.
      */
      virtual void examine(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & start_time,
        const timeSystem::AbsoluteTime & stop_time, EphStatusCont & eph_status) const;

      /// \brief Returns an object of the same type as this object.
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
          \param ephemerides Container of pulsar ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const PulsarEph & choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Choose the best ephemeris for the given absolute time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param ephemerides Container of orbital ephemerides from which the closest ephemeris is to be found.
          \param abs_time The time of interest.
      */
      virtual const OrbitalEph & choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Return ephemeris gaps and ephemeris remarks overlapped with a given time interval.
          \param ephemerides Container of pulsar ephemerides whose time coverage is to be examined. Not used by this class.
          \param start_time The start time of a time interval of interest.
          \param stop_time The stop time of a time interval of interest.
          \param eph_status Container of ephemeris status that stores the result of examination. The orginal contents
                 of the container will be removed by this method.
      */
      virtual void examine(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & start_time,
        const timeSystem::AbsoluteTime & stop_time, EphStatusCont & eph_status) const;

      /// \brief Returns an object of the same type as this object.
      virtual EphChooser * clone() const;

    protected:
      StrictEphChooser m_strict_chooser;
  };

}
#endif
