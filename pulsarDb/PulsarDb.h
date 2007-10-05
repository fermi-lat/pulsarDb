/** \file PulsarDb.h
    \brief Interface for PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarDb_h
#define pulsarDb_PulsarDb_h

#include <map>
#include <string>
#include <set>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "tip/FileSummary.h"
#include "tip/Table.h"

namespace tip {
  class Extension;
  class Table;
}

namespace pulsarDb {

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

      /** \brief Select ephemerides whose SOLAR_SYSTEM_EPHEMERIS value matches the input solar_eph. Matching is case insensitive.
          \param solar_eph The name of the solar system ephemeris. Examples: JPL DE405, JPL DE200
      */
      virtual void filterSolarEph(const std::string & solar_eph);

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

      /// \brief Get the currently selected container of spin (pulsar) ephemerides.
      virtual void getEph(PulsarEphCont & cont, const TimingModel & model = TimingModel()) const;

      /// \brief Get the currently selected container of orbital ephemerides.
      virtual void getEph(OrbitalEphCont & cont) const;

      /// \brief Get the number of currently selected ephemerides.
      virtual int getNumEph(bool spin_table = true) const;

    protected:
      PulsarDb();
      /** \brief Load data tables from the associated pulsar ephemerides database file.
          \param edit_in_place If true any filterings will affect the input file. If false, filtering
          will take place on a copy of the file in memory only.
      */
      virtual void loadTables(bool edit_in_place);

      /** \brief Clean up all extensions based on current set of selected spin and orbital ephemerides. All
          information in the OBSERVERS and ALTERNATIVE_NAMES extension which is not associated with a pulsar
          contained in the SPIN_PARAMETERS or ORBITAL_PARAMETERS extensions will be removed.
      */
      virtual void clean();

      /** \brief Creates a filtering expression from an input field name and a set of accepted values.
          \param field_name The name of the field on which to filter.
          \param values The set of allowed values.
      */
      virtual std::string createFilter(const std::string & field_name, const std::set<std::string> & values) const;

      /** \brief Helper method to get a value from a cell, returning it as a double, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
      */
      double get(const tip::TableCell & cell) const;

      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
          \param value Variable to store the value.
      */
      template <typename T>
      void get(const tip::TableCell & cell, T & value) const;

      std::string m_in_file;
      tip::FileSummary m_summary;
      TableCont m_table;
      tip::Table * m_spin_par_table;
      tip::Table * m_orbital_par_table;
      tip::Table * m_obs_code_table;
      tip::Table * m_psr_name_table;
  };

  inline double PulsarDb::get(const tip::TableCell & cell) const {
    double value = 0.;
    get(cell, value);
    return value;
  }

  template <typename T>
  inline void PulsarDb::get(const tip::TableCell & cell, T & value) const {
    // WARNING: This will break for a string column.
    if (cell.isNull()) value = 0;
    else cell.get(value);
  }

}
#endif
