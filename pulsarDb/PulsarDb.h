/** \file PulsarDb.h
    \brief Interface for PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarDb_h
#define pulsarDb_PulsarDb_h

#include <string>

#include "pulsarDb/PulsarEph.h"

namespace tip {
  class Table;
}

namespace pulsarDb {

  /** \class PulsarDb
      \brief Abstraction providing access to pulsar ephemerides database.
  */
  class PulsarDb {
    public:
      /** \brief Create a data base access object for the given ephermerides db file.
                 This opens a copy of the file in memory. The version on disk will be
                 unaffected.
          \param in_file The input file name.
      */
      PulsarDb(const std::string & in_file);

      virtual ~PulsarDb() throw();

      /** \brief Filter the current database using the given row filtering expression. The filtering is
                 performed in place.
          \param expression The row filtering expression. Examples: f2 != 0., #row > 50 && #row < 100
      */
      virtual void filter(const std::string & expression);

      /** \brief Select ephemerides whose validation interval [VALID_SINCE, VALID_UNTIL] overlaps the given time range.
          \param t_start The start time of the interval.
          \param t_stop The stop time of the interval. (An exception will be thrown if t_stop < t_start).
      */
      virtual void filterInterval(double t_start, double t_stop);

      /** \brief Select ephemerides whose name matches the input name. Matching is case insensitive.

                 The name will be looked up two ways. It will be directly compared to the PSRNAME field in
                 the SPIN_PARAMETERS table. In addition, synonyms will be found in the ALTERNATIVE_NAMES table
                 and these will then be used to look up ephemerides from the SPIN_PARAMETERS table.
          \param name The name of the pulsar. Examples: Crab, PSR B0531+21, PSR J0534+2200
      */
      virtual void filterName(const std::string & name);

      /** \brief Save the currently selected ephemerides into an output file.

                 All tables from the input file will be copied to the output file, but in each table, records which
                 do not match the current selection of ephemerides will be omitted. The matching is done based on PSRNAME for
                 the ORBITAL_PARAMETERS and ALTERNATIVE_NAMES, and on OBSERVER_CODE for the OBSERVERS table.
                 For example, if none of the ephemerides are for binary pulsars, the output ORBITAL_PARAMETERS
                 table will not include any ephemerides.
          \param out_file The name of the output file.
      */
      virtual void save(const std::string & out_file);

      /** \brief Choose the best ephemeris for the given MJD time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL + 1.)
                 (One is added to the endpoint because the validity expires at the end of that day.)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param mjd The MJD time (TT system).
      */
      virtual const PulsarEph & chooseEph(long double mjd, bool extrapolate = false) const;

      /// \brief Get the currently selected container of ephemerides.
      virtual const PulsarEphCont & getEph() const;

      /// \brief Get the number of currently selected ephemerides.
      virtual int getNumEph() const;

    protected:
      virtual const PulsarEph & chooseValidEph(long double mjd) const;
      virtual const PulsarEph & extrapolateEph(long double mjd) const;

    private:
      std::string m_in_file;
      mutable PulsarEphCont m_ephemerides;
      tip::Table * m_spin_par_table;
      const tip::Table * m_psr_name_table;
  };

}
#endif
