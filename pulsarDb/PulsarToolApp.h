/** \file PulsarToolApp.h
    \brief Declaration of base class for pulsar tool applications.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarToolApp_h
#define pulsarDb_PulsarToolApp_h

#include <map>
#include <set>
#include <string>
#include <vector>

#include "EphStatus.h"

#include "st_app/StApp.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/EventTimeHandler.h"

#include "tip/Table.h"

namespace st_app {
  class AppParGroupd;
}

namespace st_stream {
  class OStream;
}

namespace timeSystem {
  class EventTimeHandler;
}

namespace tip {
  class Header;
}

namespace pulsarDb {

  class EphChooser;
  class EphComputer;

  /** \class PulsarToolApp
      \brief Base class for pulsar tool application classes.
  */
  class PulsarToolApp : public st_app::StApp {
    public:
      typedef std::vector<timeSystem::EventTimeHandler *> handler_cont_type;

      enum TimeCorrectionMode_e {
        REQUIRED, ALLOWED, SUPPRESSED
      };

      /// \brief Construct a PulsarToolApp object.
      PulsarToolApp();

      /// \brief Destruct this PulsarToolApp object.
      virtual ~PulsarToolApp() throw();

      /** \brief Run the application. This method calls runApp method after cleaning up the internal data members,
                 and clean them up again after the call of runApp method.
      */
      virtual void run();

      /// \brief Run the tool-specific part of the application.
      virtual void runApp() = 0;

      /** \brief Interpret a character string as a time representation in a specified format, create an AbsoluteTime object
                 that represents the given time string, and return it.
          \param time_format Character string that specifies a time format (e.g., "MJD"). If this parameter is set to "FILE",
                 then the same time format is chosen as the opened event file that is listed first.
          \param time_system Name of time system in which the time string is to be interpreted. If this parameter is set to
                 "FILE", then the time system name is read from TIMESYS keyword of a given FITS header.
          \param time_value Character string that represents a time of interest.
          \param parsed_time_format Name of time format that this method determines to use is to be set as a return value.
          \param parsed_time_system Name of time system that this method determines to use is to be set as a return value.
          \param event_time_handler Pointer to an EventTimeHandler object that is used in case that time_format and/or time_system
                 is "FILE". If set to zero (0), this method throws an exception if time_format or time_system is "FILE".
      */
      virtual timeSystem::AbsoluteTime parseTime(const std::string & time_format, const std::string & time_system,
        const std::string & time_value, std::string & parsed_time_format, std::string & parsed_time_system,
        const timeSystem::EventTimeHandler * event_time_handler = 0) const;

      /** \brief Open event file(s).
          \param pars Set of parameters from which event file names, etc., are to be read.
          \param read_only If true, the file(s) are to be opened in a read-only mode, otherwise in a read-write mode.
      */
      void openEventFile(const st_app::AppParGroup & pars, bool read_only = true);

      /** \brief Search for a specified field name in the opened event file(s), and append a new column if missing.
          \param field_name Name of FITS column to be searched for.
          \param field_format Format of the FITS column to be used in case of appending a new column.
      */
      void reserveOutputField(const std::string & field_name, const std::string & field_format);

      /** \brief Define a set of time corrections and give it a name for future reference.
          \param mode_name Name of a time correction mode to define.
          \param tcmode_bary Time correction mode on the arrival time correction.
          \param tcmode_bin Time correction mode on the binary demodulation.
          \param tcmode_pdot Time correction mode on the p-dot cancellation.
      */
      void defineTimeCorrectionMode(const std::string & mode_name, TimeCorrectionMode_e tcmode_bary, TimeCorrectionMode_e tcmode_bin,
        TimeCorrectionMode_e tcmode_pdot);

      /** \brief Select a set of timing modes to use after a call of this method.
          \param mode_name Name of the time correction mode to select.
      */
      void selectTimeCorrectionMode(const std::string & mode_name);

      /** \brief Select a set of timing modes to use after a call of this method. The mode name is read from a given set of parameters.
          \param pars Set of parameters from which a time correction mode is read.
      */
      void selectTimeCorrectionMode(const st_app::AppParGroup & pars);

      /** \brief Initialize ephemeris computations, with the ephemeris style explicitly specified.
          \param pars Set of parameters from which initialization information is to be read.
          \param chooser Ephemeris chooser to be used in the application.
          \param eph_style Character string that specifies ephemeris style such as "FREQ" or "DB".
          \param os Output stream to write the ancestry records of pulsar ephemeris database to.
      */
      void initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser, const std::string & eph_style,
        st_stream::OStream & os);

      /** \brief Initialize ephemeris computations, with the ephemeris style to be read from a given set of parameters.
          \param pars Set of parameters from which initialization information is to be read.
          \param chooser Ephemeris chooser to be used in the application.
          \param os Output stream to write the ancestry records of pulsar ephemeris database to.
      */
      void initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser, st_stream::OStream & os);

      /** \brief Initializa arrival time corrections, with the time origin to be read from a given set of parameters.
          \param pars Set of parameters from which initialization information is to be read.
          \param vary_ra_dec If true, the sky position of a celestial source is to be read from the stored spin ephemerides
                 for each event. Otherwise, one position is to be used for all events.
          \param guess_pdot If true, parameters for p-dot cancellation is to be computed from the stored spin ephemerides.
                 Otherwise, the method reads them from the given set of parameters.
      */
      void initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot);

      /** \brief Initializa arrival time corrections, with the time origin specified by a character string.
          \param pars Set of parameters from which initialization information is to be read.
          \param vary_ra_dec If true, the sky position of a celestial source is to be read from the stored spin ephemerides
                 for each event. Otherwise, one position is to be used for all events.
          \param guess_pdot If true, parameters for p-dot cancellation is to be computed from the stored spin ephemerides.
                 Otherwise, the method reads them from the given set of parameters.
          \param str_origin Character string representing the time origin.
      */
      void initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot, const std::string & str_origin);

      /** \brief Initializa arrival time corrections, with the time origin explicitly given.
          \param pars Set of parameters from which initialization information is to be read.
          \param vary_ra_dec If true, the sky position of a celestial source is to be read from the stored spin ephemerides
                 for each event. Otherwise, one position is to be used for all events.
          \param guess_pdot If true, parameters for p-dot cancellation is to be computed from the stored spin ephemerides.
                 Otherwise, the method reads them from the given set of parameters.
          \param abs_origin AbsoluteTime object representing the time origin.
      */
      void initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot,
        const timeSystem::AbsoluteTime & abs_origin);

      /** \brief Compute the number of elepased seconds since the time origin.
          \param abs_time Absolute time for which the number of elapsed seconds is to be computed.
      */
      double computeElapsedSecond(const timeSystem::AbsoluteTime & abs_time) const;

      /** \brief Create an AbsoluteTime object that represents the time after a given elapsed time since the time origin.
          \param elapsed_time Number of elapsed seconds since the time origin.
      */
      timeSystem::AbsoluteTime computeAbsoluteTime(double elapsed_time) const;

      /// \brief Set the internal event iterator to point to the first event in the opened FITS file.
      void setFirstEvent();

      /// \brief Increment the internal event iterator to point to the next event in the opened FITS file.
      void setNextEvent();

      /// \brief Return if the internal event iterator points to one past the last event in the opened FITS file.
      bool isEndOfEventList() const;

      /// \brief Create an AbsoluteTime object that represents the current event time, and return it.
      timeSystem::AbsoluteTime getEventTime() const;

      /** \brief Write a value to a field of the current FITS row.
          \param field_name Name of the field to be written a given value.
          \param field_value Value of the field to be written.
      */
      template <typename DataType>
      void setFieldValue(const std::string & field_name, const DataType & field_value) {
        tip::TableRecord & record = (*m_event_handler_itor)->getCurrentRecord();
        record[field_name].set(field_value);
      }

      /// \brief Create an AbsoluteTime object that represents the earliest time in GTI extension(s).
      timeSystem::AbsoluteTime getStartTime() const;

      /// \brief Create an AbsoluteTime object that represents the latest time in GTI extension(s).
      timeSystem::AbsoluteTime getStopTime() const;

      /// \brief Return a reference to the internal ephemeris computer for external use.
      EphComputer & getEphComputer() const;

      /** \brief Examine ephemeris coverage and ephemeris remarks relevant for the opened event file(s),
                 and output a summary report of the ephemeris status into a given output stream.
          \param os Output stream to write the summary report to.
          \param code_to_report List of ephemeris status code to be reported.
      */
      void reportEphStatus(st_stream::OStream & os, const std::set<EphStatusCodeType> & code_to_report) const;

      /** \brief Examine ephemeris coverage and ephemeris remarks relevant for a given time interval,
                 and output a summary report of the ephemeris status into a given output stream.
          \param os Output stream to write the summary report to.
          \param start_time Start time of the time interval to examine ephemerides.
          \param stop_time Stop time of the time interval to examine ephemerides.
          \param code_to_report List of ephemeris status code to be reported.
      */
      void reportEphStatus(st_stream::OStream & os, const timeSystem::AbsoluteTime & start_time,
        const timeSystem::AbsoluteTime & stop_time, const std::set<EphStatusCodeType> & code_to_report) const;

      /** \brief Write a list of parameters to the header(s) of the event file(s).
          \param pars Set of parameters to be written to the event file(s).
      */
      void writeParameter(const st_app::AppParGroup & pars);

    private:
      handler_cont_type m_event_handler_cont;
      handler_cont_type m_gti_handler_cont;
      std::string m_time_field;
      std::string m_gti_start_field;
      std::string m_gti_stop_field;
      std::vector<std::pair<std::string, std::string> > m_output_field_cont;
      timeSystem::EventTimeHandler * m_reference_handler;
      EphComputer * m_computer;
      std::map<const std::string, TimeCorrectionMode_e> m_tcmode_dict_bary;
      std::map<const std::string, TimeCorrectionMode_e> m_tcmode_dict_bin;
      std::map<const std::string, TimeCorrectionMode_e> m_tcmode_dict_pdot;
      TimeCorrectionMode_e m_tcmode_bary;
      TimeCorrectionMode_e m_tcmode_bin;
      TimeCorrectionMode_e m_tcmode_pdot;
      bool m_request_bary;
      bool m_demod_bin;
      bool m_cancel_pdot;
      bool m_vary_ra_dec;
      const timeSystem::TimeSystem * m_target_time_system;
      timeSystem::AbsoluteTime m_target_time_origin;
      handler_cont_type::iterator m_event_handler_itor;
      bool m_report_eph_status;

      /** \brief Reset all members of application. This should be called from the subclass's run() method
          to allow multiple runs to work without leaking memory or coupling consecutive runs of the tool
          accidentally.
      */
      virtual void resetApp();

      /** \brief Read a time value from a specified FITS column, create an AbosoluteTime object that represents the time,
                 apply time correctios on it if requested, and return it.
          \param handler EventTimeHandler object to read a time value from.
          \param column_name Name of FITS column from which a time value to be read.
          \param request_time_correction If true, arrival time corrections are to be requested on the time. Otherwise,
                 a time read from the event file is converted into an AbsoluteTime object without arrival time corrections.
      */
      timeSystem::AbsoluteTime readTimeColumn(timeSystem::EventTimeHandler & handler, const std::string & column_name,
        bool request_time_correction) const;

      /** \brief Compute a time boundary (start time or stop time) of the opened event file(s), and return it.
          \param request_start_time If true, start time is computed. Otherwise, stop time is computed.
          \param request_time_correction If true, arrival time corrections are to be requested on the time. Otherwise,
                 a time read from the event file is converted into an AbsoluteTime object without arrival time corrections.
      */
      timeSystem::AbsoluteTime computeTimeBoundary(bool request_start_time, bool request_time_correction) const;

      /// \brief Set up the internal variables when a new FITS table is opened.
      void setupCurrentEventTable();
  };

}

#endif
