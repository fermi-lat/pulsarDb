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

#include "tip/FileSummary.h"
#include "tip/Table.h"

namespace tip {
  class Extension;
  class Header;
}

namespace pulsarDb {

  /** \brief Base class to define the interface of ephemeris factory classes.
             Note: Template parameter, EPHTYPE, must be either PulsarEph or OrbitalEph.
  */
  template <typename EPHTYPE>
  class IEphFactory {
    public:
      virtual EPHTYPE * create(const tip::Table::ConstRecord & record, const tip::Header & header) const = 0;
  };

  /** \brief Concrete class to create a PulsarEph object or an OrbitalEph object.
             Note: The first emplate parameter, EPHTYPE, must be either PulsarEph or OrbitalEph.
                   The second template parameter, EPHSTYLE, must be one of their subclasses.
  */
  template <typename EPHTYPE, typename EPHSTYLE>
  class EphFactory: public IEphFactory<EPHTYPE> {
    public:
      virtual EPHTYPE * create(const tip::Table::ConstRecord & record, const tip::Header & header) const {
        return new EPHSTYLE(record, header);
      }

      static EphFactory & getFactory() {
        static EphFactory s_factory;
        return s_factory;
      }

    private:
      EphFactory() {};
  };

  /** \class PulsarDb
      \brief Abstraction providing access to pulsar ephemerides database.
  */
  class PulsarDb {
    public:
      typedef std::vector<tip::Table *> TableCont;

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
      virtual void getEph(PulsarEphCont & cont) const;

      /// \brief Get the currently selected container of orbital ephemerides.
      virtual void getEph(OrbitalEphCont & cont) const;

      /// \brief Get the number of currently selected ephemerides.
      virtual int getNumEph(bool spin_table = true) const;

      /** \brief Associate a keyword with a PulsarEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a PulsarEph class.
                 Note: Template parameter, EPHSTYLE, must be one of PulsarEph subclasses.
      */
      // TODO: Write a test code for this method.
      template <typename EPHSTYLE>
      void registerPulsarEph(const std::string & eph_style) {
        m_spin_factory_cont[eph_style] = &EphFactory<PulsarEph, EPHSTYLE>::getFactory();
      }

      /** \brief Associate a keyword with a OrbitalEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a OrbitalEph class.
                 Note: Template parameter, EPHSTYLE, must be one of OrbitalEph subclasses.
      */
      // TODO: Write a test code for this method.
      template <typename EPHSTYLE>
      void registerOrbitalEph(const std::string & eph_style) {
        m_orbital_factory_cont[eph_style] = &EphFactory<OrbitalEph, EPHSTYLE>::getFactory();
      }

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

      std::string m_in_file;
      tip::FileSummary m_summary;
      TableCont m_all_table;
      TableCont m_spin_par_table;
      TableCont m_orbital_par_table;
      TableCont m_obs_code_table;
      TableCont m_psr_name_table;

    private:
      typedef std::map<std::string, IEphFactory<PulsarEph> *> spin_factory_cont_type;
      typedef std::map<std::string, IEphFactory<OrbitalEph> *> orbital_factory_cont_type;

      spin_factory_cont_type m_spin_factory_cont;
      orbital_factory_cont_type m_orbital_factory_cont;
  };
}
#endif
