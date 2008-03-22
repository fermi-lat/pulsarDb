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
#include <fstream>
#include <list>
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

  PulsarDb::PulsarDb(const std::string & tpl_file, TableCont::size_type default_spin_ext, TableCont::size_type default_orbital_ext):
    m_tpl_file(tpl_file), m_tip_file(), m_all_table(), m_spin_par_table(), m_orbital_par_table(), m_obs_code_table(), m_psr_name_table(),
    m_default_spin_par_table(0), m_default_orbital_par_table(0), m_spin_factory_cont(), m_orbital_factory_cont() {
    // Create a FITS file in memory that holds ephemeris data.
    m_tip_file = IFileSvc::instance().createMemFile("pulsardb.fits", m_tpl_file);

    // Get summary of extensions in the memory FITS file.
    FileSummary file_summary;
    IFileSvc::instance().getFileSummary(m_tip_file.getName(), file_summary);

    // If memory FITS file is empty, throw exception.
    FileSummary::const_iterator ext_itor = file_summary.begin();
    if (file_summary.end() == ext_itor) throw std::runtime_error("Given FITS template produces an empty FITS file");

    // Loop over remaining extensions in output file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    int ext_number = 1;
    for (++ext_itor; ext_itor != file_summary.end(); ++ext_itor, ++ext_number) {
      std::string ext_name = ext_itor->getExtId();

      // Get a table and its extension name.
      Table * table = 0;
      std::ostringstream oss;
      oss << ext_number;
      table = m_tip_file.editTable(oss.str());

      // Check EPHSTYLE header keyword in SPIN_PARAMETERS and ORBITAL_PARAMETERS extensions.
      if ("SPIN_PARAMETERS" == ext_name || "ORBITAL_PARAMETERS" == ext_name) {
        Header & header(table->getHeader());
        std::string eph_style;
        try {
          header["EPHSTYLE"].get(eph_style);
        } catch (const TipException &) {
          throw std::runtime_error("Could not find header keyword \"EPHSTYLE\" in extension \"" + ext_name + "\" in FITS template \""
            + tpl_file + "\"");
        }
      }

      // Add this table to an appropriate list of tables.
      m_all_table.push_back(table);
      if ("SPIN_PARAMETERS" == ext_name) {
        m_spin_par_table.push_back(table);
      } else if ("ORBITAL_PARAMETERS" == ext_name) {
        m_orbital_par_table.push_back(table);
      } else if ("OBSERVERS" == ext_name) {
        m_obs_code_table.push_back(table);
      } else if ("ALTERNATIVE_NAMES" == ext_name) {
        m_psr_name_table.push_back(table);
      } else {
        // Do nothing for other extensions.
        ;
      }
    }

    // Set default spin and orbital tables.
    m_default_spin_par_table = findDefaultTable(default_spin_ext, true);
    m_default_orbital_par_table = findDefaultTable(default_orbital_ext, false);
  }

  Table * PulsarDb::findDefaultTable(TableCont::size_type ext_number, bool is_spin_table) const {
    Table * return_table(0);
     std::string table_type = (is_spin_table ? "spin" : "orbital");

    if (ext_number >= m_all_table.size()) { 
      // Throw an exception for extension number out of bound.
      std::ostringstream os;
      os << "Default extension number for " << table_type << " parameters out of bound: " << ext_number << std::endl;
      throw std::runtime_error(os.str());

    } else if (ext_number == 0) {
      // Return zero to indicate no default.
      return_table = 0;

    } else {
      // Set a designated table to the return value.
      return_table = m_all_table[ext_number - 1];

      // Check whether the chosen table is really a table of the requested type (spin or orbital).
      const TableCont & table_cont = (is_spin_table ? m_spin_par_table : m_orbital_par_table);
      const std::set<const Table *> table_set(table_cont.begin(), table_cont.end());
      if (table_set.find(return_table) == table_set.end()) {
        std::ostringstream os;
        os << "Default " << table_type << " parameter extension is set to an extension not for " << table_type << " parameters";
        throw std::runtime_error(os.str());
      }
    }
 
    return return_table;
 }

  void PulsarDb::load(const std::string & in_file) {
    // Determine a type of the input file, FITS or TEXT.
    bool is_fits = true;
    try {
      // Try opening a primary extension of a FITS file.
      std::auto_ptr<const Extension> ext(IFileSvc::instance().readExtension(in_file, "0"));
    } catch (const TipException &) {
      is_fits = false;
    }

    // Call an appropriate loader method.
    if (is_fits) loadFits(in_file);
    else loadText(in_file);

    // Refresh the tip table objects.
    // TODO: Do we really need to do this?
    for (TableCont::iterator itor = m_all_table.begin(); itor != m_all_table.end(); ++itor) {
      (*itor)->filterRows("#row>0");
    }
  }

  PulsarDb::~PulsarDb() {
    // Clear out orbital ephemeris factories.
    m_orbital_factory_cont.clear();

    // Clear out pulsar ephemeris factories.
    m_spin_factory_cont.clear();

    // Delete all tip::Table objects.
    for (TableCont::reverse_iterator itor = m_all_table.rbegin(); itor != m_all_table.rend(); ++itor) delete *itor;
  }

  void PulsarDb::filter(const std::string & expression) {
    for (TableCont::iterator itor = m_spin_par_table.begin(); itor != m_spin_par_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
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
    for (TableCont::iterator itor = m_spin_par_table.begin(); itor != m_spin_par_table.end(); ++itor) {
      (*itor)->filterRows(os.str());
    }
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
    for (TableCont::iterator table_itor = m_psr_name_table.begin(); table_itor != m_psr_name_table.end(); ++table_itor) {
      const Table & name_table = **table_itor;
      for (Table::ConstIterator record_itor = name_table.begin(); record_itor != name_table.end(); ++record_itor) {
        Table::ConstRecord & record(*record_itor);
        // See if there is an alternate name.
        std::string alt_name;
        record["ALTNAME"].get(alt_name);

        // Convert to uppercase.
        for (std::string::iterator s_itor = alt_name.begin(); s_itor != alt_name.end(); ++s_itor) *s_itor = toupper(*s_itor);

        if (alt_name == up_name) {
          // Found a match.
          std::string real_name;
          record["PSRNAME"].get(real_name);
          name_matches.insert(real_name);
        }
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
    for (TableCont::iterator itor = m_spin_par_table.begin(); itor != m_spin_par_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
    for (TableCont::iterator itor = m_orbital_par_table.begin(); itor != m_orbital_par_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
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
    for (TableCont::iterator itor = m_spin_par_table.begin(); itor != m_spin_par_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
    for (TableCont::iterator itor = m_orbital_par_table.begin(); itor != m_orbital_par_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
    clean();
  }

  void PulsarDb::save(const std::string & out_file, const std::string & creator, bool clobber) const {
    // Copy the entire contents of the memory FITS file to an output file.
    m_tip_file.copyFile(out_file, clobber);

    // Loop over all extensions in output file and update header keywords.
    FileSummary file_summary;
    IFileSvc::instance().getFileSummary(out_file, file_summary);
    int ext_number = 0;
    for (FileSummary::const_iterator ext_itor = file_summary.begin(); ext_itor != file_summary.end();
       ++ext_itor, ++ext_number) {
      // Open this extension by extension number.
      std::ostringstream oss;
      oss << ext_number;
      std::auto_ptr<Extension> ext(IFileSvc::instance().editExtension(out_file, oss.str()));
      Header & header(ext->getHeader());

      // Update header keywords.
      Header::KeyValCont_t keywords;
      keywords.push_back(Header::KeyValPair_t("DATE", header.formatTime(time(0))));
      keywords.push_back(Header::KeyValPair_t("CREATOR", creator));
      header.update(keywords);
    }
  }

  void PulsarDb::clean() {
    // Set of pulsar names and observer codes.
    std::set<std::string> pulsar_names;
    std::set<std::string> observer_codes;

    // Compose a set of all pulsar names and observer codes in spin extension. This is used below
    // to filter the other extensions in the file.
    for (TableCont::iterator table_itor = m_spin_par_table.begin(); table_itor != m_spin_par_table.end(); ++table_itor) {
      const Table & spin_table = **table_itor;
      for (Table::ConstIterator in_itor = spin_table.begin(); in_itor != spin_table.end(); ++in_itor) {
        std::string tmp;

        (*in_itor)["PSRNAME"].get(tmp);
        pulsar_names.insert(tmp);

        (*in_itor)["OBSERVER_CODE"].get(tmp);
        observer_codes.insert(tmp);
      }
    }

    // Compose a set of all pulsar names and observer codes in orbital extension. This is used below
    // to filter the other extensions in the file.
    for (TableCont::iterator table_itor = m_orbital_par_table.begin(); table_itor != m_orbital_par_table.end(); ++table_itor) {
      const Table & orbital_table = **table_itor;
      for (Table::ConstIterator in_itor = orbital_table.begin(); in_itor != orbital_table.end(); ++in_itor) {
        std::string tmp;

        (*in_itor)["PSRNAME"].get(tmp);
        pulsar_names.insert(tmp);

        (*in_itor)["OBSERVER_CODE"].get(tmp);
        observer_codes.insert(tmp);
      }
    }

    // Filter observer codes based on codes found in the spin/orbital extensions.
    std::string expression = createFilter("OBSERVER_CODE", observer_codes);
    for (TableCont::iterator itor = m_obs_code_table.begin(); itor != m_obs_code_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }

    // Filter alternative names based on names found in the spin/orbital extensions.
    expression = createFilter("PSRNAME", pulsar_names);
    for (TableCont::iterator itor = m_psr_name_table.begin(); itor != m_psr_name_table.end(); ++itor) {
      (*itor)->filterRows(expression);
    }
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

    // Reserve space for all ephemerides.
    int num_record = 0;
    for (TableCont::const_iterator itor = m_spin_par_table.begin(); itor != m_spin_par_table.end(); ++itor) {
      num_record += (*itor)->getNumRecords();
    }
    cont.reserve(num_record);

    for (TableCont::const_iterator table_itor = m_spin_par_table.begin(); table_itor != m_spin_par_table.end(); ++table_itor) {
      const Table & spin_table = **table_itor;
      const Header & header(spin_table.getHeader());

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const TipException &) {
        // Note: EPHSTYLE must exist in SPIN_PARAMETERS extension, and it is enforced in the constructor of this class.
        //       Not finding EPHSTYLE here suggests inconsistency in methods of this class.
        throw std::logic_error("EPHSTYLE header keyword is missing in SPIN_PARAMETERS extension");
      }

      // Use a registered subclass of PulsarEph, if EPHSTYLE keyword exists.
      const IEphFactory<PulsarEph> * factory(0);
      spin_factory_cont_type::const_iterator factory_itor = m_spin_factory_cont.find(eph_style);
      if (factory_itor != m_spin_factory_cont.end()) {
        factory = factory_itor->second;
      } else {
        throw std::runtime_error("Unknown pulsar ephemeris style: EPHSTYLE = " + eph_style);
      }

      // Iterate over current selection.
      for (Table::ConstIterator record_itor = spin_table.begin(); record_itor != spin_table.end(); ++record_itor) {
        // For convenience, get record from iterator.
        Table::ConstRecord & record(*record_itor);

        // Add the ephemeris to the container.
        cont.push_back(factory->create(record, header));
      }
    }
  }

  void PulsarDb::getEph(OrbitalEphCont & cont) const {
    // Empty container then refill it.
    cont.clear();

    for (TableCont::const_iterator table_itor = m_orbital_par_table.begin(); table_itor != m_orbital_par_table.end(); ++table_itor) {
      const Table & orbital_table = **table_itor;
      cont.reserve(orbital_table.getNumRecords());
      const Header & header(orbital_table.getHeader());

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const TipException &) {
        // Note: EPHSTYLE must exist in ORBITAL_PARAMETERS extension, and it is enforced in the constructor of this class.
        //       Not finding EPHSTYLE here suggests inconsistency in methods of this class.
        throw std::logic_error("EPHSTYLE header keyword is missing in ORBITAL_PARAMETERS extension");
      }

      // Use a registered subclass of OrbitalEph, if EPHSTYLE keyword exists.
      const IEphFactory<OrbitalEph> * factory(0);
      orbital_factory_cont_type::const_iterator factory_itor = m_orbital_factory_cont.find(eph_style);
      if (factory_itor != m_orbital_factory_cont.end()) {
        factory = factory_itor->second;
      } else {
        throw std::runtime_error("Unknown orbital ephemeris style: EPHSTYLE = " + eph_style);
      }

      for (Table::ConstIterator record_itor = orbital_table.begin(); record_itor != orbital_table.end(); ++record_itor) {
        // For convenience, get record from iterator.
        Table::ConstRecord & record(*record_itor);

        // Add the ephemeris to the container.
        cont.push_back(factory->create(record, header));
      }
    }
  }

  int PulsarDb::getNumEph(bool spin_table) const {
    int num_eph = 0;

    // Get the right table container.
    const TableCont & table_cont = (spin_table ? m_spin_par_table : m_orbital_par_table);

    // Get the number of entry in the requested table.
    for (TableCont::const_iterator itor = table_cont.begin(); itor != table_cont.end(); ++itor) {
      num_eph += (*itor)->getNumRecords();
    }

    // Return the number of entry in the requested table.
    return num_eph;
  }

  void PulsarDb::loadFits(const std::string & in_file) {
    // Get summary of extensions in input file.
    FileSummary in_summary;
    IFileSvc::instance().getFileSummary(in_file, in_summary);

    // Get first extension in input file.
    FileSummary::const_iterator ext_itor = in_summary.begin();

    // If input file was empty, do nothing.
    if (in_summary.end() == ext_itor) return;

    // Loop over all extensions in input file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    int ext_number = 1;
    for (++ext_itor; ext_itor != in_summary.end(); ++ext_itor, ++ext_number) {
      const std::string ext_name = ext_itor->getExtId();

      // Open input table by extension number, in order to work with multiple extensions of the same extension name.
      std::ostringstream oss;
      oss << ext_number;
      std::auto_ptr<const Table> in_table(IFileSvc::instance().readTable(in_file, oss.str()));

      // Collect header keyword to require.
      const Header & in_header(in_table->getHeader());
      Header::KeySeq_t required_keyword;
      bool ephstyle_found = false;
      for (Header::ConstIterator key_itor = in_header.begin(); key_itor != in_header.end(); ++key_itor) {
        // Ignore some header keywords: HISTORY, COMMENT, a blank name, CHECKSUM, DATASUM, DATE, CREATOR, NAXIS2.
        std::string keyword_name(key_itor->getName());
        if (!keyword_name.empty()) {
          for (std::string::iterator str_itor = keyword_name.begin(); str_itor != keyword_name.end(); ++str_itor) {
            *str_itor = std::toupper(*str_itor);
          }
          if (keyword_name != "HISTORY" && keyword_name != "COMMENT" && keyword_name != "CHECKSUM" && keyword_name != "DATASUM"
            && keyword_name != "DATE" && keyword_name != "CREATOR" && keyword_name != "NAXIS2") required_keyword.push_back(*key_itor);

          // Use EPHSTYLE header keyword as an indicator of the current format.
          if (keyword_name == "EPHSTYLE") ephstyle_found = true;
        }
      }

      // Find the right tip::Table object.
      Table * target_table(0);
      if ("SPIN_PARAMETERS" == ext_name && !ephstyle_found) {
        // Treat the input as a spin parameter table in the original format.
        target_table = m_default_spin_par_table;

      } else if ("ORBITAL_PARAMETERS" == ext_name && !ephstyle_found) {
        // Treat the input as an orbital parameter table in the original format.
        target_table = m_default_orbital_par_table;

      } else {
        // Treat the input as any table in the current format.
        target_table = updateMatchingHeader(required_keyword);
      }

      if (target_table == 0) {
        // Throw an exception if no extension is found to load to.
        throw std::runtime_error("Could not find an extension to load the contents of extension \"" + ext_name + "\" in file \""
          + in_file + "\"");
      }

      // Copy all rows in the table if a matching extension is found, start at beginning of both tables.
      Table::ConstIterator in_itor = in_table->begin();
      Table::Iterator target_itor = target_table->end();

      // TODO: Resize output to match input. Requires tip to handle random access iterators.
      // target_table->setNumRecords(in_table->getNumRecords() + target_table->getNumRecords());

      // Copy all rows in table.
      for (; in_itor != in_table->end(); ++in_itor, ++target_itor) *target_itor = *in_itor;
    }
  }

  Table * PulsarDb::updateMatchingHeader(const Header::KeySeq_t & required_keyword) {
    // Prepare for the matching and updating.
    Table * matched_table(0);
    Header::KeySeq_t::size_type num_required_keyword(required_keyword.size());
    typedef std::list<std::pair<tip::Header::Iterator, tip::Header::ConstIterator> > MatchedPairCont;
    MatchedPairCont updater_pair;

    // Walk through all tables.
    for (TableCont::iterator table_itor = m_all_table.begin(); table_itor != m_all_table.end(); ++table_itor) {
      Table & table(**table_itor);
      Header & header(table.getHeader());

      // Prepare containers of matched pairs.
      MatchedPairCont tight_match; // Contains pairs whose keyword names and values match.
      MatchedPairCont loose_match; // Contains pairs whose keyword names match and at least one of their values is empty.

      // Find matching keyword-value pairs in the header, by case-insensitive string comparison.
      for (Header::Iterator header_key_itor = header.begin(); header_key_itor != header.end(); ++header_key_itor) {
        // Get keyword name and make it case-insensitive.
        std::string name_header(header_key_itor->getName());
        for (std::string::iterator str_itor = name_header.begin(); str_itor != name_header.end(); ++str_itor) {
          *str_itor = std::toupper(*str_itor);
        }

        // Check all required records.
        for (Header::ConstIterator required_key_itor = required_keyword.begin(); required_key_itor != required_keyword.end();
          ++required_key_itor) {
          // Get keyword name and make it case-insensitive.
          std::string name_required(required_key_itor->getName());
          for (std::string::iterator str_itor = name_required.begin(); str_itor != name_required.end(); ++str_itor) {
            *str_itor = std::toupper(*str_itor);
          }

          // Compare two records.
          if (name_header == name_required) {
            // Mark as a loosely matched pair if one of them is empty.
            if (header_key_itor->empty() || required_key_itor->empty()) {
              loose_match.push_back(std::make_pair(header_key_itor, required_key_itor));
            } else {
              // Get keyword values and make them case-insensitive.
              std::string string_value_header(header_key_itor->getValue());
              for (std::string::iterator str_itor = string_value_header.begin(); str_itor != string_value_header.end(); ++str_itor) {
                *str_itor = std::toupper(*str_itor);
              }
              std::string string_value_required(required_key_itor->getValue());
              for (std::string::iterator str_itor = string_value_required.begin(); str_itor != string_value_required.end(); ++str_itor) {
                *str_itor = std::toupper(*str_itor);
              }

              // Compare keyword values as string.
              bool value_identical = false;
              if (string_value_header == string_value_required) {
                value_identical = true;
              } else {
                // Get keyword values as double.
                double double_nan = std::numeric_limits<double>::quiet_NaN();
                double double_value_header(double_nan);
                header_key_itor->getValue(double_value_header);
                double double_value_required(double_nan);
                header_key_itor->getValue(double_value_required);
                if (double_value_header == double_value_required && double_value_header != double_nan) value_identical = true;
              }
              if (value_identical) tight_match.push_back(std::make_pair(header_key_itor, required_key_itor));
            }
          }
        }
      }

      // Append loosely matched pairs to tightly matched pairs.
      MatchedPairCont matched_pair(tight_match);
      matched_pair.insert(matched_pair.end(), loose_match.begin(), loose_match.end());

      // Prepare for containers for match evaluation.
      std::set<Header::Iterator> header_key_taken;
      std::set<Header::ConstIterator> required_key_taken;

      // Evaluate goodness of matching; first check tightly matched pairs.
      updater_pair.clear();
      for (MatchedPairCont::const_iterator pair_itor = matched_pair.begin(); pair_itor != matched_pair.end(); ++pair_itor) {
        const Header::Iterator & header_key(pair_itor->first);
        const Header::ConstIterator & required_key(pair_itor->second);
        if (header_key_taken.find(header_key) == header_key_taken.end()
          && required_key_taken.find(required_key) == required_key_taken.end()) {
          header_key_taken.insert(header_key);
          required_key_taken.insert(required_key);

          // Collect matched pairs which requires to update this header.
          if (header_key->empty() && !required_key->empty()) updater_pair.push_back(*pair_itor);
        }
      }

      // Judge completeness of matchup.
      if (required_key_taken.size() == num_required_keyword) {
        if (matched_table) throw std::runtime_error("More than one extension contain the given list of header keywords.");
        else matched_table = &table;
      }
    }

    // Update the matched header.
    if (matched_table) {
      for (MatchedPairCont::iterator pair_itor = updater_pair.begin(); pair_itor != updater_pair.end(); ++pair_itor) {
        Header::Iterator & header_key(pair_itor->first);
        Header::ConstIterator & required_key(pair_itor->second);
        header_key->setValue(required_key->getValue());
      }
    }

    // Return the matched table, if only one match found.
    return matched_table;
  }

  void PulsarDb::loadText(const std::string & in_file) {
    // Open input text file.
    std::ifstream in_table(in_file.c_str());

    // Make sure it worked.
    if (!in_table) throw std::runtime_error("Could not open file " + in_file);

    // Create a text buffer.
    static const size_t s_line_size = 2048;
    char buf[s_line_size];

    Table * target_table(0);
    ParsedLine parsed_line;
    ParsedLine field_name;
    Header::KeySeq_t header_keyword;

    // Read table until extension name is found.
    std::string ext_name;
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

      if (parsed_line.size() > 1) throw std::runtime_error(std::string("Expected name of an extension on line ") + buf);
      else if (parsed_line.size() == 1) {
        ext_name = *parsed_line.begin();
        header_keyword.push_back(KeyRecord("EXTNAME", ext_name, "name of this binary table extension"));
        break;
      }
    }
      
    std::map<std::string, size_t> found_map;

    // Read table until names of fields are found.
    bool keyword_found = false;
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

      // Ignore comment lines, etc.
      if (parsed_line.size() == 0) continue;

      // Check whether it is a header keyword line, by looking for the equal sign ("=") in the line.
      std::string str_buf(buf);
      std::string::size_type pos_equal = str_buf.find('=');

      if (pos_equal != std::string::npos) {
        // Parse out keyword name.
        std::string str_name = str_buf.substr(0, pos_equal);
        stripWhiteSpace(str_name);
        str_buf.erase(0, pos_equal + 1);

        // Parse out keyword value.
        std::string::size_type pos_slash = str_buf.find('/');
        std::string str_value = str_buf.substr(0, pos_slash);
        stripWhiteSpace(str_value);
        if (pos_slash == std::string::npos) str_buf.clear();
        else str_buf.erase(0, pos_slash + 1);

        // Parse out keyword comment.
        std::string str_comment = str_buf;
        stripWhiteSpace(str_comment);

        // Add this header keyword, in order to require for an appropriate table for this text file.
        KeyRecord key_record(str_name, str_value, str_comment);
        header_keyword.push_back(key_record);
        keyword_found = true;

      } else {
        // Find an appropriate table and update its header.
        if ("SPIN_PARAMETERS" == ext_name && !keyword_found) {
          // Treat the input as a spin parameter table in the original format.
          target_table = m_default_spin_par_table;

        } else if ("ORBITAL_PARAMETERS" == ext_name && !keyword_found) {
          // Treat the input as an orbital parameter table in the original format.
          target_table = m_default_orbital_par_table;

        } else {
          // Treat the input as any table in the current format.
          target_table = updateMatchingHeader(header_keyword);
        }

        // Throw an exception if no extension is found to load to.
        if (target_table == 0) {
          throw std::runtime_error("Could not find an extension to load the contents of extension \"" + ext_name + "\" in file \""
            + in_file + "\"");
        }

        // Get list of fields in the table.
        Table::FieldCont field = target_table->getValidFields();
        field_name.assign(field.begin(), field.end());

        // Convert to lowercase, because field names are all lowercase.
        for (ParsedLine::iterator itor = parsed_line.begin(); itor != parsed_line.end(); ++itor) {
          for (std::string::iterator s_itor = itor->begin(); s_itor != itor->end(); ++s_itor) *s_itor = tolower(*s_itor);
        }

        // Populate input_column number as a function of output column name.
        for (ParsedLine::iterator itor = field_name.begin(); itor != field_name.end(); ++itor) {
          std::vector<std::string>::iterator found_itor = std::find(parsed_line.begin(), parsed_line.end(), *itor);
          if (parsed_line.end() != found_itor) {
            found_map[*itor] = found_itor - parsed_line.begin();
          }
        }
        break;
      }
    }
  
    if (found_map.empty())
      throw std::runtime_error("File " + in_file + " does not have any columns in common with output extension " + ext_name);

    // Populate output table starting at the beginning.
    Table::Iterator target_itor = target_table->begin();

    // Read the rest of input table to populate fields.
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

      // Ignore comment lines, etc.
      if (parsed_line.size() == 0) continue;

      // Make sure the correct number of fields are found.
      if (parsed_line.size() != found_map.size()) {
        std::ostringstream os;
        os << "Line " << buf << " has " << parsed_line.size() << " token(s), not " << found_map.size() << ", as expected";
        throw std::runtime_error(os.str());
      }
  
      // Populate the parsed line by iterating over output columns.
      for (ParsedLine::iterator in_itor = field_name.begin(); in_itor != field_name.end(); ++in_itor) {
        // See if this column was found amongst the input table columns.
        std::map<std::string, size_t>::iterator found_itor = found_map.find(*in_itor);
        if (found_map.end() != found_itor) {
          // Found column, so look up the input table column number from the output column number.
          (*target_itor)[*in_itor].set(parsed_line[found_itor->second]);
        } else {
          // TODO: Make this write an INDEF in the field (requires further changes to tip).
          // (*target_itor)[*in_itor].set([found_itor->second]);
        }
      }
      ++target_itor;
    }
  }

  void PulsarDb::parseLine(const char * line, ParsedLine & parsed_line) {
    parsed_line.clear();

    const char * begin = line;

    // Skip leading whitespace.
    while(isspace(*begin)) ++begin;

    // Handle the case of a completely blank line.
    if ('\0' == *begin) return;

    // Check for comment and return if it is one.
    if ('#' == *begin) return;

    // At this point there is at least one character, hence one field.
    parsed_line.resize(1);

    // Start with the first field.
    std::vector<std::string>::iterator current_field = parsed_line.begin();

    // Start with no quotes.
    bool in_quote = false;
    const char * current = begin;
    while ('\0' != *current) {
      // Check for quote, skip them but toggle the flag.
      if ('"' == *current) {
        in_quote = !in_quote;
        ++current;
        continue;
      }

      // Unquoted space is taken as the end of this field.
      if (!in_quote && isspace(*current)) {
        // Add a new field at the end of the container.
        current_field = parsed_line.insert(parsed_line.end(), "");

        // Skip spaces after this field.
        while (isspace(*current)) ++current;

        // Start the next field.
        begin = current;
        continue;
      }

      // Add the character to the current field.
      current_field->push_back(*current);

      // Continue to next character.
      ++current;
    }

    if (in_quote) throw std::runtime_error(std::string("Line ") + line + " contains an unbalanced quote");
  }

  void PulsarDb::stripWhiteSpace(std::string & string_value) {
    // Remove leading white spaces.
    for (std::string::iterator itor = string_value.begin(); itor != string_value.end() && std::isspace(*itor); ++itor) {
      string_value.erase(itor);
    }

    // Remove trailing white spaces.
    for (std::string::reverse_iterator r_itor = string_value.rbegin(); r_itor != string_value.rend() && std::isspace(*r_itor);
      ++r_itor) {
      std::string::iterator itor = r_itor.base();
      --itor;
      string_value.erase(itor);
    }
  }
}
