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

#include "pulsarDb/EphStatus.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

#include "tip/FileSummary.h"
#include "tip/Table.h"
#include "tip/TipException.h"
#include "tip/TipFile.h"

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
      virtual ~IEphFactory() {}

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

      /** \brief Create a data base object for the given FITS template file.
          \param tpl_file The name of the FITS template file used to create the file in memory.
          \param default_spin_ext The extension number of SPIN_PARAMETER extension in a given template file (tpl_file),
                 to which spin ephemerides in pulsar database file in the original format (w/o EPHSTYLE) will be loaded.
                 If zero (0) is given, an exception will be thrown for spin-ephemeris extensions in the original format.
                 Give one (1) for the first HDU after the primary HDU.
          \param default_orbital_ext Same as above, but for ORBITAL_PARAMETER extensions.
      */
      PulsarDb(const std::string & tpl_file, TableCont::size_type default_spin_ext = 0, TableCont::size_type default_orbital_ext = 0);

      virtual ~PulsarDb();

      /** \brief Load ephemerides and related information from the given ephemerides database file.
                 This copies the file in memory, and the version on disk will be unaffected.
          \param in_file The name of the input ephemerides database file.
      */
      virtual void load(const std::string & in_file);

      /** \brief Filter the current database using the given row filtering expression. The filtering is
                 performed in place.
          \param expression The row filtering expression. Examples: f2 != 0., #row > 50 && #row < 100
          \param filter_spin If true, spin parameters will be filtered by the given expression.
          \param filter_orbital If true, orbital parameters will be filtered by the given expression.
          \param filter_remark If true, ephemeris remarks will be filtered by the given expression.
      */
      virtual void filterExpression(const std::string & expression, bool filter_spin = true, bool filter_orbital = false,
        bool filter_remark = false);

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
          \param creator The character string to be assigned to a value of CREATOR header keyword.
          \param author The character string to be assigned to a value of AUTHOR header keyword.
          \param clobber If true, it overwrites the output file even if it already exists.  If no and the output file
                 already exists, it throws an exception.
      */
      virtual void save(const std::string & out_file, const std::string & creator, const std::string & author,
        bool clobber = false) const;

      /** \brief Get the currently selected spin (pulsar) ephemerides.
          \param cont The container to fill the currently selected spin ephemerides in it. The previous contents of
                 the container will be removed by this method.
      */
      virtual void getEph(PulsarEphCont & cont) const { getEphBody(m_spin_par_table, m_spin_factory_cont, cont); }

      /** \brief Get the currently selected orbital ephemerides.
          \param cont The container to fill the currently selected orbital ephemerides in it. The spin ephemeris objects
                 that are stored in this container before calling this method will be destroyed and removed from the container.
      */
      virtual void getEph(OrbitalEphCont & cont) const { getEphBody(m_orbital_par_table, m_orbital_factory_cont, cont); }

      /** \brief Get the number of currently selected ephemerides.
          \param cont The container to fill the currently selected orbital ephemerides in it. The orbital ephemeris objects
                 that are stored in this container before calling this method will be destroyed and removed from the container.
      */
      virtual int getNumEph(bool spin_table = true) const;

      /** \brief Associate a keyword with a PulsarEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a PulsarEph class.
                 Note: Template parameter, EPHSTYLE, must be one of PulsarEph subclasses.
      */
      template <typename EPHSTYLE>
      void registerPulsarEph(const std::string & eph_style) {
        m_spin_factory_cont[eph_style] = &EphFactory<PulsarEph, EPHSTYLE>::getFactory();
      }

      /** \brief Associate a keyword with a OrbitalEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a OrbitalEph class.
                 Note: Template parameter, EPHSTYLE, must be one of OrbitalEph subclasses.
      */
      template <typename EPHSTYLE>
      void registerOrbitalEph(const std::string & eph_style) {
        m_orbital_factory_cont[eph_style] = &EphFactory<OrbitalEph, EPHSTYLE>::getFactory();
      }

      /** \brief Get the remarks stored in this pulsar ephemerides database.
          \param cont The container to fill the remarks in it. The EphStatus objects that are stored in this container
                 before calling this method will be destroyed and removed from the container.
      */
      virtual void getRemark(EphStatusCont & cont) const;

    private:
      typedef std::vector<std::string> ParsedLine;

      /** \brief Helper method for the constructor, to check and pick up an appropriate table for
          a default spin or orbital extension to support the original ephemerides file format.
          \param ext_number The extension number. Same as constructor arguments, default_spin/orbital_ext.
          \param is_spin_table If true, it looks for a spin parameter table. If false, it looks for an
          orbital parameter table.
      */
      tip::Table * findDefaultTable(TableCont::size_type ext_number, bool is_spin_table) const;

      /** \brief Load ephemerides and related information from the given FITS file.
          \param in_file The name of the input FITS file.
      */
      virtual void loadFits(const std::string & in_file);

      /** \brief Helper method to find an extension that contains all required keyword-value pairs,
          and update its header if necessary.
          \param required_keyword List of keyword-value pairs that must be found in an extension to find.
      */
      virtual tip::Table * updateMatchingHeader(const tip::Header::KeySeq_t & required_keyword);

      /** \brief Load ephemerides and related information from the given TEXT file.
          \param in_file The name of the input TEXT file.
      */
      virtual void loadText(const std::string & in_file);

      /** \brief Helper method to parse a text line for loadText method.
          \param line The line to parse (input).
          \param parsed_line The parsed line (output).
      */
      virtual void parseLine(const char * line, ParsedLine & parsed_line);

      /** \brief Helper method to strip out leading and trailing white spaces.
          \param string_value The string object from which white spaces are to be stripped.
      */
      virtual void stripWhiteSpace(std::string & string_value);

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

      /** \brief Helper method for getEph methods for PulsarEphCont and OrbitalEphCont.
          \param table_cont The container of tip tables, from which ephemerides are extracted.
          \param factory_cont The container of ephemeris factories to be used to create appropriate ephemeris objects.
          \param eph_cont The ephemeris container to be filled with appropriate ephemeris objects.
          Note: the original contents of this object will be removed from the container.
      */
      template <typename FactoryCont, typename EphCont>
      void getEphBody(const TableCont & table_cont, const FactoryCont & factory_cont, EphCont & eph_cont) const;

      std::string m_tpl_file;
      tip::TipFile m_tip_file;
      TableCont m_all_table;
      TableCont m_spin_par_table;
      TableCont m_orbital_par_table;
      TableCont m_eph_remark_table;
      TableCont m_obs_code_table;
      TableCont m_psr_name_table;
      tip::Table * m_default_spin_par_table;
      tip::Table * m_default_orbital_par_table;
      std::map<std::string, IEphFactory<PulsarEph> *> m_spin_factory_cont;
      std::map<std::string, IEphFactory<OrbitalEph> *> m_orbital_factory_cont;
  };

  template <typename FactoryCont, typename EphCont>
  void PulsarDb::getEphBody(const TableCont & table_cont, const FactoryCont & factory_cont, EphCont & eph_cont) const {
    // Empty container then refill it.
    for (typename EphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
    eph_cont.clear();

    // Reserve space for all ephemerides.
    int num_record = 0;
    for (TableCont::const_iterator itor = table_cont.begin(); itor != table_cont.end(); ++itor) {
      num_record += (*itor)->getNumRecords();
    }
    eph_cont.reserve(num_record);

    for (TableCont::const_iterator table_itor = table_cont.begin(); table_itor != table_cont.end(); ++table_itor) {
      const tip::Table & table = **table_itor;
      const tip::Header & header(table.getHeader());

      // Get the extension name.
      std::string ext_name;
      header["EXTNAME"].get(ext_name);

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const tip::TipException &) {
        // Note: EPHSTYLE must exist in SPIN_PARAMETERS and ORBITAL_PARAMETERS extensions, and it is enforced
        //       in the constructor of this class. Not finding EPHSTYLE here suggests inconsistency in methods
        //       of this class, or a bug most likely.
        throw std::logic_error("EPHSTYLE header keyword is missing in " + ext_name + " extension");
      }

      // Use a registered subclass of PulsarEph or OrbitalEph whichever appropriate, if EPHSTYLE keyword exists.
      typename FactoryCont::mapped_type factory(0);
      typename FactoryCont::const_iterator factory_itor = factory_cont.find(eph_style);
      if (factory_itor != factory_cont.end()) {
        factory = factory_itor->second;
      } else {
        throw std::runtime_error("Unknown ephemeris style for " + ext_name + " extension: EPHSTYLE = " + eph_style);
      }

      // Iterate over current selection.
      for (tip::Table::ConstIterator record_itor = table.begin(); record_itor != table.end(); ++record_itor) {
        // For convenience, get record from iterator.
        tip::Table::ConstRecord & record(*record_itor);

        // Add the ephemeris to the container.
        eph_cont.push_back(factory->create(record, header));
      }
    }
  }
}
#endif
