/** \file EphComputer.h
    \brief Interface for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputer_h
#define pulsarDb_EphComputer_h

#include <string>
#include <utility>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace timeSystem {
  class AbsoluteTime;
}

namespace pulsarDb {
  class EphChooser;
  class PdotCanceler;
  class PulsarDb;

  /** \class EphComputer
      \brief Abstraction providing high-level access to all main functions of pulsarDb package.
  */
  class EphComputer {
    public:
      /// \brief Construct an EphComputer object.
      EphComputer();

      /** \brief Construct an EphComputer object, with a specific ephemeris chooser.
          \par chooser Ephemeris chooser for this object to use in choosing the "best" ephemeris for a given time.
      */
      EphComputer(const EphChooser & chooser);

      /// \brief Destruct this EphComputer object.
      ~EphComputer();

      /** \brief Load spin ephemerides, orbital ephemerides, and ephemeris remarks from a given database to an internal storage.
          \param database Pulsar ephemeris database from which ephemeris data are loaded.
      */
      void load(const PulsarDb & database);

      /** \brief Load spin ephemerides from a given database.
          \param database Pulsar ephemeris database from which spin ephemerides are loaded to an internal storage.
      */
      void loadPulsarEph(const PulsarDb & database);

      /** \brief Load a given spin ephemeris to an internal storage.
          \param pulsar_eph Spin ephemeris to be loaded.
      */
      void loadPulsarEph(const PulsarEph & pulsar_eph);

      /** \brief Load orbital ephemerides from a given database.
          \param database Pulsar ephemeris database from which orbital ephemerides are loaded to an internal storage.
      */
      void loadOrbitalEph(const PulsarDb & database);

      /** \brief Load a given orbital ephemeris to an internal storage.
          \param pulsar_eph Orbital ephemeris to be loaded.
      */
      void loadOrbitalEph(const OrbitalEph & orbital_eph);

      /** \brief Load ephemeris remarks from a given database.
          \param database Pulsar ephemeris database from which ephemeris remarks are loaded to an internal storage.
      */
      void loadEphRemark(const PulsarDb & database);

      /** \brief Load a given ephemeris remark to an internal storage.
          \param pulsar_eph Ephemeris remark to be loaded.
      */
      void loadEphRemark(const EphStatus & eph_remark);

      /** \brief Set parameters for the p-dot cancellation to an internal storage.
          \param time_system_name Name of time system in which given parameters are interpreted.
          \param time_origin Origin of time for the pdot-cancellation.
          \param fdot_ratio Vector of frequency derivatives at the time origin divided by a frequency at the time origin,
                 where fdot_ratio[0] stores a ratio of the first time derivative of frequency over the frequency,
                 fdot_ratio[1] stores a ratio of the second time derivative of frequency over the frequency, and so on.
      */
      void setPdotCancelParameter(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
        const std::vector<double> & fdot_ratio);


      /** \brief Compute parameters for the p-dot cancellation from a give spin ephemeris, and set the result to an internal storage.
          \param time_origin Origin of time for the pdot-cancellation.
          \param pulsar_eph Pulsar's spin ephemeris from which ratios of frequency derivatives over a frequency is to be computed.
          \param max_derivative Maximum order of frequency derivative to be computed for use in the p-dot cancellation.
      */
      void setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, const PulsarEph & pulsar_eph, int max_derivative);

      /** \brief Compute parameters for the p-dot cancellation from spin ephemerides stored in the internal storage of
                 spin ephemeris, and set the result to an internal storage.
          \param time_origin Origin of time for the pdot-cancellation.
          \param max_derivative Maximum order of frequency derivative to be computed for use in the p-dot cancellation.
      */
      void setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, int max_derivative);

      /** \brief Apply p-dot cancellation to a given time, and update the argument.
          \param ev_time Absolute time to be corrected for p-dot cancellation.
      */
      void cancelPdot(timeSystem::AbsoluteTime & ev_time) const;

      /** \brief Compute a pulse phase for a given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time for which a pulse phase is to be compted.
          \param phase_offset Value to be added to a pulse phase. This value is added to the computed pulse phase
                 before truncated to a value in range [0, 1).
      */
      double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute a sky position at a given time, and return it. The returned value is a pair of Right Ascension
                 and Declination of the position.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which a sky position is to be computed.
      */
      std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

      /** \brief Compute an orbital phase of a given time, and return it.
          \param ev_time Absolute time for which an orbital phase is to be computed.
          \param phase_offset Value to be added to an orbital phase. This value is added to the computed orbital phase
                 before truncated to a value in range [0, 1).
      */
      double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute an orbital delay for a photon emission time, add it to the time, and set the result to the time.
          \param emission_time Time at which a photo leaves a star in a binary system. This argument is to be updated with
                 the result of the computation.
      */
      void modulateBinary(timeSystem::AbsoluteTime & emission_time) const;

      /** \brief Estimate a photon emission time from a photon arrival time, and set the result to the argument as a return value.
          \param arrival_time Time at which a photon arrives at an observer outside of a binary system. This argument is to be updated
                 with the result of the computation.
      */
      void demodulateBinary(timeSystem::AbsoluteTime & arrival_time) const;

      /// \brief Return the number of spin ephemerides stored in the internal storage.
      PulsarEphCont::size_type getNumPulsarEph() const;

      /// \brief Return the number of orbital ephemerides stored in the internal storage.
      OrbitalEphCont::size_type getNumOrbitalEph() const;

      /// \brief Return the number of ephemeris remarks stored in the internal storage.
      EphStatusCont::size_type getNumEphRemark() const;

      /** \brief Choose the "best" spin ephemeris for a given time, and return it.
          \param ev_time Absolute time for which a spin ephemeris is to be chosen.
      */
      const PulsarEph & choosePulsarEph(const timeSystem::AbsoluteTime & ev_time) const;

      /** \brief Choose the "best" orbital ephemeris for a given time, and return it.
          \param ev_time Absolute time for which an orbital ephemeris is to be chosen.
      */
      const OrbitalEph & chooseOrbitalEph(const timeSystem::AbsoluteTime & ev_time) const;

      /** \brief Examine time coverage by spin ephemerides in the internal storage in the given time interval,
                 and return the result as a list of ephemeris remarks of various types, such as ephemeris gaps
                 and informative notes. Interpretation of "time coverage" may depend on the type of an EphChooser
                 object stored in this object.
          \param start_time The start time of a time interval of interest.
          \param stop_time The stop time of a time interval of interest.
          \param eph_status Container of ephemeris status that stores the result of examination. The original contents
                 of the container will be removed by this method.
      */
      void examinePulsarEph(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
        EphStatusCont & eph_status_cont) const;

      /** \brief Search for ephemeris remarks overlapped with a given time interval, and return the retuls as a list of
                 ephemeris remarks set to the argument of this method.
          \param start_time The start time of a time interval of interest.
          \param stop_time The stop time of a time interval of interest.
          \param eph_status Container of ephemeris status that stores the ephemeris remarks found. The original contents
                 of the container will be removed by this method.
      */
      void getEphRemark(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
        EphStatusCont & eph_status_cont) const;

      /** \brief Summarize time systems of pulsar's spin ephemerides, orbital ephemerides, and pdot cancellation,
                 and return character strings summarizing them for human reading through arguments of this method.
          \param spin_summary Character string to receive a summary of time systems for spin ephemerides. The original contents
                 of the string will be replaced with the summary.
          \param orbital_summary Character string to receive a summary of time systems for orbital ephemerides. The original contents
                 of the string will be replaced with the summary.
          \param pdot_summary Character string to receive a summary of time system for pdot cancellation. The original contents
                 of the string will be replaced with the summary.
      */
      void summarizeTimeSystem(std::string & spin_summary, std::string & orbital_summary, std::string & pdot_summary) const;

    private:
      /// \brief Copy constructors.
      EphComputer(const EphComputer &);
      EphComputer & operator =(const EphComputer &);

      PulsarEphCont m_pulsar_eph_cont;
      OrbitalEphCont m_orbital_eph_cont;
      EphStatusCont m_eph_remark_cont;
      PdotCanceler * m_pdot_canceler;
      EphChooser * m_chooser;
  };

}

#endif
