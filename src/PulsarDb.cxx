/** \file PulsarDb.cxx
    \brief Implementation of the PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/

#include <cfloat>
#include <cmath>
#include <ctime>
#include <cctype>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "facilities/commonUtilities.h"

#include "pulsarDb/PulsarDb.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace timeSystem;
using namespace tip;

namespace pulsarDb {

  PulsarDb::PulsarDb(const std::string & in_file, bool edit_in_place): m_in_file(in_file), m_summary(), m_table(),
    m_spin_par_table(0), m_orbital_par_table(0), m_obs_code_table(0), m_psr_name_table(0), m_spin_factory_cont(),
    m_orbital_factory_cont() {
    loadTables(edit_in_place);
  }

  PulsarDb::~PulsarDb() {
    // Clear out orbital ephemeris factories.
    for (orbital_factory_cont_type::reverse_iterator itor = m_orbital_factory_cont.rbegin();
      itor != m_orbital_factory_cont.rend(); ++itor) {
      delete itor->second;
    }
    m_orbital_factory_cont.clear();

    // Clear out pulsar ephemeris factories.
    for (spin_factory_cont_type::reverse_iterator itor = m_spin_factory_cont.rbegin();
      itor != m_spin_factory_cont.rend(); ++itor) {
      delete itor->second;
    }
    m_spin_factory_cont.clear();

    // TODO: Change Tip so it has reverse_iterator, then use it to delete in the correct order.
    for (FileSummary::const_iterator itor = m_summary.begin(); itor != m_summary.end(); ++itor) {
      TableCont::iterator table_itor = m_table.find(itor->getExtId());
      if (m_table.end() != table_itor) delete table_itor->second;
    }
  }

  void PulsarDb::filter(const std::string & expression) {
    m_spin_par_table->filterRows(expression);
  }

  void PulsarDb::filterInterval(double t_start, double t_stop) {
    std::ostringstream os;

    if (t_start > t_stop) {
      os << "PulsarDb::filterInterval was passed interval with start > stop: [" << t_start << ", " << t_stop << "]" << std::endl;
      throw std::runtime_error(os.str());
    }

    // Condition for disjoint intervals is:
    // t_stop < VALID_SINCE || VALID_UNTIL < t_start
    // Condition for overlapping is the negation:
    // !(t_stop < VALID_SINCE || VALID_UNTIL < t_start)
    // This is logically equivalent to:
    // t_stop > VALID_SINCE && VALID_UNTIL > t_start
    os << t_stop << ">VALID_SINCE&&VALID_UNTIL>" << t_start;
    m_spin_par_table->filterRows(os.str());
    clean();
  }

  void PulsarDb::filterName(const std::string & name) {
    std::string up_name = name;

    // Convert to uppercase.
    for (std::string::iterator s_itor = up_name.begin(); s_itor != up_name.end(); ++s_itor) *s_itor = toupper(*s_itor);

    // Check for the (limited) wildcard.
    if (up_name == "ALL" || up_name == "ANY") return;

    std::set<std::string> name_matches;
    // Look up the given name in the extension containing alternate names.
    for (Table::Iterator itor = m_psr_name_table->begin(); itor != m_psr_name_table->end(); ++itor) {
      // See if there is an alternate name.
      std::string alt_name;
      (*itor)["ALTNAME"].get(alt_name);

      // Convert to uppercase.
      for (std::string::iterator s_itor = alt_name.begin(); s_itor != alt_name.end(); ++s_itor) *s_itor = toupper(*s_itor);

      if (alt_name == up_name) {
        // Found a match.
        std::string real_name;
        (*itor)["PSRNAME"].get(real_name);
        name_matches.insert(real_name);
      }
    }

    // If the all-upper case name is different from the original name, add it.
    if (up_name != name) {
      name_matches.insert(up_name);
    }

    // Finally (just to be sure) add the original name given.
    name_matches.insert(name);

    std::string expression = createFilter("PSRNAME", name_matches);

    // OK, do the filter for real on the spin parameters and orbital parameters tables.
    m_spin_par_table->filterRows(expression);
    m_orbital_par_table->filterRows(expression);
    clean();
  }

  void PulsarDb::filterSolarEph(const std::string & solar_eph) {
    // Convert to uppercase.
    std::string solar_eph_uc(solar_eph);
    for (std::string::iterator itor = solar_eph_uc.begin(); itor != solar_eph_uc.end(); ++itor) *itor = toupper(*itor);

    // Create a list of solar ephemeris names to match.
    std::set<std::string> name_matches;
    name_matches.insert(solar_eph_uc);

    // If the all-upper case name is different from the original name, add it.
    if (solar_eph_uc != solar_eph) {
      name_matches.insert(solar_eph_uc);
    }

    // Finally (just to be sure) add the original name given.
    name_matches.insert(solar_eph);

    std::string expression = createFilter("SOLAR_SYSTEM_EPHEMERIS", name_matches);

    // OK, do the filter for real on the spin parameters and orbital parameters tables.
    m_spin_par_table->filterRows(expression);
    m_orbital_par_table->filterRows(expression);
    clean();
  }

  void PulsarDb::save(const std::string & out_file, const std::string & tpl_file) const {
    // Get first extension in input file.
    FileSummary::const_iterator ext_itor = m_summary.begin();

    // If input file was empty, throw exception.
    if (m_summary.end() == ext_itor) throw std::runtime_error("Input ephemeris file is empty");

    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Find data directory for this app.
    std::string data_dir = facilities::commonUtilities::getDataPath("pulsarDb");

    if (!file_svc.fileExists(out_file)) {
      // Create output file.
      file_svc.createFile(out_file, tpl_file);
    }

    // Update keywords only in primary extension.
    std::auto_ptr<Extension> primary(file_svc.editExtension(out_file, ext_itor->getExtId()));

    // Update standard keywords.
    updateKeywords(*primary);

    // Loop over all extensions in input file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    for (++ext_itor; ext_itor != m_summary.end(); ++ext_itor) {
      TableCont::const_iterator table_itor = m_table.find(ext_itor->getExtId());
      if (m_table.end() == table_itor)
        throw std::logic_error("Could not find extension \"" + ext_itor->getExtId() + "\" in file \"" + m_in_file + "\"");

      // Input table pointer.
      const Table * in_table = table_itor->second;

      // Open output table.
      std::auto_ptr<Table> out_table(file_svc.editTable(out_file, ext_itor->getExtId()));

      // Update standard keywords.
      updateKeywords(*out_table);

      // Start at beginning of both tables.
      Table::ConstIterator in_itor = in_table->begin();
      Table::Iterator out_itor = out_table->end();

      // TODO: Resize output to match input. Requires tip to handle random access iterators.
//      out_table->setNumRecords(in_table->getNumRecords() + out_table->getNumRecords());

      // Copy all rows in table.
      for (; in_itor != in_table->end(); ++in_itor, ++out_itor) *out_itor = *in_itor;
    }
  }

  void PulsarDb::updateKeywords(tip::Extension & ext) const {
    // Update necessary keywords.
    tip::Header::KeyValCont_t keywords;
    tip::Header & header(ext.getHeader());

    keywords.push_back(tip::Header::KeyValPair_t("DATE", header.formatTime(time(0))));
    keywords.push_back(tip::Header::KeyValPair_t("CREATOR", "gtpulsardb"));

    header.update(keywords);
  }

  void PulsarDb::clean() {
    // Set of pulsar names and observer codes.
    std::set<std::string> pulsar_names;
    std::set<std::string> observer_codes;

    // Compose a set of all pulsar names and observer codes in spin extension. This is used below
    // to filter the other extensions in the file.
    for (Table::Iterator in_itor = m_spin_par_table->begin(); in_itor != m_spin_par_table->end(); ++in_itor) {
      std::string tmp;

      (*in_itor)["PSRNAME"].get(tmp);
      pulsar_names.insert(tmp);

      (*in_itor)["OBSERVER_CODE"].get(tmp);
      observer_codes.insert(tmp);
    }

    // Compose a set of all pulsar names and observer codes in orbital extension. This is used below
    // to filter the other extensions in the file.
    for (Table::Iterator in_itor = m_orbital_par_table->begin(); in_itor != m_orbital_par_table->end(); ++in_itor) {
      std::string tmp;

      (*in_itor)["PSRNAME"].get(tmp);
      pulsar_names.insert(tmp);

      (*in_itor)["OBSERVER_CODE"].get(tmp);
      observer_codes.insert(tmp);
    }

    // Filter observer codes based on codes found in the spin/orbital extensions.
    std::string expression = createFilter("OBSERVER_CODE", observer_codes);
    m_obs_code_table->filterRows(expression);

    // Filter alternative names based on names found in the spin/orbital extensions.
    expression = createFilter("PSRNAME", pulsar_names);
    m_psr_name_table->filterRows(expression);

  }

  std::string PulsarDb::createFilter(const std::string & field_name, const std::set<std::string> & values) const {
    std::string expression;
    for (std::set<std::string>::const_iterator itor = values.begin(); itor != values.end(); ++itor) {
      expression += "(" + field_name + "==\"" + *itor + "\")||";
    }

    // Trim off the last 2 characters ||
    return expression.substr(0, expression.size() - 2);
  }

  void PulsarDb::getEph(PulsarEphCont & cont) const {
    // Refill this container.
    cont.clear();

    cont.reserve(m_spin_par_table->getNumRecords());
    const tip::Header & header(m_spin_par_table->getHeader());

    // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
    const IEphFactory<PulsarEph> * factory(0);
    std::string eph_style;
    try {
      header["EPHSTYLE"].get(eph_style);
    } catch (const tip::TipException &) {
      // Use FrequencyEph if EPHSTYLE keyword is missing.
      static const EphFactory<PulsarEph, FrequencyEph> s_default_spin_factory;
      factory = &s_default_spin_factory;
    }

    // Use a registered subclass of OrbitalEph, if EPHSTYLE keyword exists.
    if (0 == factory) {
      spin_factory_cont_type::const_iterator itor = m_spin_factory_cont.find(eph_style);
      if (itor != m_spin_factory_cont.end()) {
        factory = itor->second;
      } else {
        throw std::runtime_error("Unknown pulsar ephemeris style: EPHSTYLE = " + eph_style);
      }
    }

    // Iterate over current selection.
    for (Table::Iterator itor = m_spin_par_table->begin(); itor != m_spin_par_table->end(); ++itor) {
      // For convenience, get record from iterator.
      Table::ConstRecord & record(*itor);

      // Add the ephemeris to the container.
      cont.push_back(factory->create(record, header));
    }
  }

  void PulsarDb::getEph(OrbitalEphCont & cont) const {
    // Empty container then refill it.
    cont.clear();

    if (0 != m_orbital_par_table) {
      cont.reserve(m_orbital_par_table->getNumRecords());
      const tip::Header & header(m_orbital_par_table->getHeader());

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      const IEphFactory<OrbitalEph> * factory(0);
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const tip::TipException &) {
        // Use SimpleDdEph if EPHSTYLE keyword is missing.
        static const EphFactory<OrbitalEph, SimpleDdEph> s_default_orbital_factory;
        factory = &s_default_orbital_factory;
      }

      // Use a registered subclass of OrbitalEph, if EPHSTYLE keyword exists.
      if (0 == factory) {
        orbital_factory_cont_type::const_iterator itor = m_orbital_factory_cont.find(eph_style);
        if (itor != m_orbital_factory_cont.end()) {
          factory = itor->second;
        } else {
          throw std::runtime_error("Unknown orbital ephemeris style: EPHSTYLE = " + eph_style);
        }
      }

      for (Table::Iterator itor = m_orbital_par_table->begin(); itor != m_orbital_par_table->end(); ++itor) {
        // For convenience, get record from iterator.
        Table::ConstRecord & record(*itor);

        // Add the ephemeris to the container.
        cont.push_back(factory->create(record, header));
      }
    }
  }

  int PulsarDb::getNumEph(bool spin_table) const {
    int num_eph = 0;

    // Get the number of entry in the requested table.
    if (spin_table) {
      num_eph = m_spin_par_table->getNumRecords();
    } else {
      num_eph = m_orbital_par_table->getNumRecords();
    }

    // Return the number of entry in the requested table.
    return num_eph;
  }

  PulsarDb::PulsarDb(): m_in_file(), m_summary(), m_table(), m_spin_par_table(0), m_orbital_par_table(0), m_obs_code_table(0),
    m_psr_name_table(0) {
  }

  void PulsarDb::loadTables(bool edit_in_place) {
    // Alias for file service singleton.
    IFileSvc & file_svc(IFileSvc::instance());

    // Get summary of extensions in input file.
    file_svc.getFileSummary(m_in_file, m_summary);

    // Loop over remaining extensions in output file.
    FileSummary::const_iterator ext_itor = m_summary.begin();

    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    for (++ext_itor; ext_itor != m_summary.end(); ++ext_itor) {
      Table * table = 0;
      std::string ext_name = ext_itor->getExtId();
      if (edit_in_place) {
        table = file_svc.editTable(m_in_file, ext_name);
      } else {
        table = file_svc.editTable(m_in_file, ext_name, "#row>0");
      }
      m_table.insert(std::make_pair(ext_name, table));
    }

    // Find spin parameter table.
    TableCont::iterator itor = m_table.find("SPIN_PARAMETERS");
    if (m_table.end() == itor) throw std::runtime_error("Could not find SPIN_PARAMETERS table in file " + m_in_file);
    m_spin_par_table = itor->second;
    
    // Find orbital parameter table.
    itor = m_table.find("ORBITAL_PARAMETERS");
    if (m_table.end() == itor) throw std::runtime_error("Could not find ORBITAL_PARAMETERS table in file " + m_in_file);
    m_orbital_par_table = itor->second;
    
    // Find observer codes.
    itor = m_table.find("OBSERVERS");
    if (m_table.end() == itor) throw std::runtime_error("Could not find OBSERVERS table in file " + m_in_file);
    m_obs_code_table = itor->second;
    
    // Find pulsar name table.
    itor = m_table.find("ALTERNATIVE_NAMES");
    if (m_table.end() == itor) throw std::runtime_error("Could not find ALTERNATIVE_NAMES table in file " + m_in_file);
    m_psr_name_table = itor->second;
  }

}
