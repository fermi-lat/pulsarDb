/** \file PulsarDb.h
    \brief Interface for PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarDb_h
#define pulsarDb_PulsarDb_h

#include <map>
#include <string>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "tip/FileSummary.h"

namespace tip {
  class Table;
}

namespace pulsarDb {

  class AbsoluteTime;

  /** \class PulsarEphCont
      \brief Abstraction providing access to a container of pulsar ephemerides.
  */
  class PulsarEphCont {
    public:
      typedef std::vector<PulsarEph *> Cont_t;

      virtual ~PulsarEphCont() { clear(); }

      void clear() {
        for (Cont_t::reverse_iterator itor = m_ephemerides.rbegin(); itor != m_ephemerides.rend(); ++itor)
          delete (*itor);
        m_ephemerides.clear();
      }

      void insertEph(const PulsarEph & eph) { m_ephemerides.push_back(eph.clone()); }

      /** \brief Choose the best ephemeris for the given MJD time. Throws an exception if no ephemeris is found.

                 The ephmeris returned contains the time in the half-open interval [VALID_SINCE, VALID_UNTIL)
                 If more than one candidate ephemeris contains the time, the ephemeris with the latest start
                 time is chosen. If more than one candidate has the same start time, the one with the latest
                 stop time is chosen. If more than one candidate has the same start and stop times, the ephemeris
                 which appears last in the table is selected.
          \param t The time of interest.
          \param strict_validity If false, the ephemeris closest to the given time will be returned even if the time
                 is outside its interval of validity. If true, only an ephemeris which contains the time in its interval
                 of validity will be returned. In either case, if no ephemeris meets the requirements, an exception
                 is thrown.
      */
      virtual const PulsarEph & chooseEph(const AbsoluteTime & t, bool strict_validity = true) const;

    protected:
      virtual const PulsarEph & chooseFromAllEph(const AbsoluteTime & t) const;
      virtual const PulsarEph & chooseFromValidEph(const AbsoluteTime & t) const;

      Cont_t m_ephemerides;
  };

  /** \class PulsarDb
      \brief Abstraction providing access to pulsar ephemerides database.
  */
  class PulsarDb {
    public:
      typedef std::map<std::string, tip::Table *> TableCont;

      /** \brief Create a data base access object for the given ephermerides db file.
                 This opens a copy of the file in memory. The version on disk will be
                 unaffected.
          \param in_file The input file name.
          \param edit_in_place Changes will affect input file.
      */
      PulsarDb(const std::string & in_file, bool edit_in_place = false);

      virtual ~PulsarDb();

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

      /// \brief Get the currently selected container of ephemerides.
      virtual void getEph(PulsarEphCont & cont) const;

      /// \brief Get the number of currently selected ephemerides.
      virtual int getNumEph() const;

    private:
      TimingModel m_model;
      std::string m_in_file;
      tip::FileSummary m_summary;
      TableCont m_table;
      tip::Table * m_spin_par_table;
      const tip::Table * m_psr_name_table;
  };

}
#endif
