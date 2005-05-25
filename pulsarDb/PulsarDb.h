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

#include "tip/FileSummary.h"

namespace tip {
  class Extension;
  class Table;
}

namespace pulsarDb {

  class AbsoluteTime;

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
          \param tpl_file The name of the template file, only used if the output file does not already exist.
      */
      virtual void save(const std::string & out_file, const std::string & tpl_file) const;

      void updateKeywords(tip::Extension & ext) const;

      /// \brief Get the currently selected container of ephemerides.
      virtual void getEph(PulsarEphCont & cont) const;

      /// \brief Get the number of currently selected ephemerides.
      virtual int getNumEph() const;

    protected:
      PulsarDb();
      virtual void loadTables(bool edit_in_place);

      std::string m_in_file;
      tip::FileSummary m_summary;
      TableCont m_table;
      tip::Table * m_spin_par_table;
      const tip::Table * m_psr_name_table;
  };

}
#endif
