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

  PulsarDb::PulsarDb(const std::string & tpl_file): m_tpl_file(tpl_file), m_tip_file(), m_summary(), m_all_table(),
    m_spin_par_table(), m_orbital_par_table(), m_obs_code_table(), m_psr_name_table(), m_spin_factory_cont(),
    m_orbital_factory_cont() {
    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Create a FITS file in memory that holds ephemeris data.
    m_tip_file = file_svc.createMemFile("pulsardb.fits", m_tpl_file);

    // Get summary of extensions in the memory FITS file.
    file_svc.getFileSummary(m_tip_file.getName(), m_summary);

    // If memory FITS file is empty, throw exception.
    FileSummary::const_iterator ext_itor = m_summary.begin();
    if (m_summary.end() == ext_itor) throw std::runtime_error("Given FITS template produces an empty FITS file");

    // Loop over remaining extensions in output file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    int ext_number = 1;
    for (++ext_itor; ext_itor != m_summary.end(); ++ext_itor, ++ext_number) {
      // Get a table and its extension name.
      Table * table = 0;
      std::ostringstream oss;
      oss << ext_number;
      table = file_svc.editTable(m_tip_file.getName(), oss.str());

      // Add this table to an appropriate list of tables.
      std::string ext_name = ext_itor->getExtId();
      if ("SPIN_PARAMETERS" == ext_name) {
        m_spin_par_table.push_back(table);
        m_all_table.push_back(table);
      } else if ("ORBITAL_PARAMETERS" == ext_name) {
        m_orbital_par_table.push_back(table);
        m_all_table.push_back(table);
      } else if ("OBSERVERS" == ext_name) {
        m_obs_code_table.push_back(table);
        m_all_table.push_back(table);
      } else if ("ALTERNATIVE_NAMES" == ext_name) {
        m_psr_name_table.push_back(table);
        m_all_table.push_back(table);
      } else {
        // Do nothing for other extensions.
        ;
      }
    }
  }

  void PulsarDb::load(const std::string & in_file) {
    // Create a temporary PulsarDb object in the same format as this object.
    std::auto_ptr<PulsarDb> pulsar_db(new PulsarDb(m_tpl_file));

    // First try opening a tip file containing the ephemerides.
    bool try_text = false;
    try {
      pulsar_db->loadFits(in_file);
    } catch (const tip::TipException &) {
      // Discard the temporary PulsarDb object because it might be corrupted.
      pulsar_db.reset(0);

      // Hope that maybe it is a text file.
      try_text = true;
    }

    // Next try opening it as a text file if there was a problem.
    if (try_text) {
      pulsar_db.reset(new PulsarDb(m_tpl_file));
      try {
        pulsar_db->loadText(in_file);
      } catch (const std::exception &) {
        throw std::runtime_error("Could not open file " + in_file + " as either a FITS file or a text data file");
      }
    }

    // Merge the temporary PulsarDb object into this object.
    loadFits(pulsar_db->m_tip_file.getName());

    // Refresh the tip table objects.
    // TODO: Do we really need to do this?
    for (TableCont::iterator itor = m_all_table.begin(); itor != m_all_table.end(); ++itor) {
      (*itor)->filterRows("#row>0");
    }
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

  void PulsarDb::save(const std::string & out_file, bool clobber) const {
    // Get first extension in input file.
    FileSummary::const_iterator ext_itor = m_summary.begin();

    // If input file was empty, throw exception.
    if (m_summary.end() == ext_itor) throw std::runtime_error("Input ephemeris file is empty");

    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Find data directory for this app.
    std::string data_dir = facilities::commonUtilities::getDataPath("pulsarDb");

    // Create output file using the same FITS template as this object, respecting clobber.
    file_svc.createFile(out_file, m_tpl_file, clobber);

    // Update keywords only in primary extension.
    std::auto_ptr<Extension> primary(file_svc.editExtension(out_file, ext_itor->getExtId()));

    // Update standard keywords.
    updateKeywords(*primary);

    // Loop over all extensions in input file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    // TODO: Revise this section for new design of D4 FITS definition.
    for (++ext_itor; ext_itor != m_summary.end(); ++ext_itor) {
      const std::string ext_name = ext_itor->getExtId();

      // Find the right tip::Table object.
      const TableCont * table_cont(0);
      const Table * in_table(0);
      if ("SPIN_PARAMETERS" == ext_name) table_cont = &m_spin_par_table;
      else if ("ORBITAL_PARAMETERS" == ext_name) table_cont = &m_orbital_par_table;
      else if ("OBSERVERS" == ext_name) table_cont = &m_obs_code_table;
      else if ("ALTERNATIVE_NAMES" == ext_name) table_cont = &m_psr_name_table;
      else throw std::runtime_error("Unknown extension name: " + ext_name);
      if (table_cont->empty()) {
        throw std::logic_error("Extension \"" + ext_name + "\" does not exist in pulsar ephemerides file");
      } else {
        in_table = *(table_cont->begin());
      }

      // Open output table.
      std::auto_ptr<Table> out_table(file_svc.editTable(out_file, ext_name));

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
      const tip::Header & header(spin_table.getHeader());

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      const IEphFactory<PulsarEph> * factory(0);
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const tip::TipException &) {
        // Use FrequencyEph if EPHSTYLE keyword is missing.
        factory = &EphFactory<PulsarEph, FrequencyEph>::getFactory();
      }

      // Use a registered subclass of OrbitalEph, if EPHSTYLE keyword exists.
      if (0 == factory) {
        spin_factory_cont_type::const_iterator factory_itor = m_spin_factory_cont.find(eph_style);
        if (factory_itor != m_spin_factory_cont.end()) {
          factory = factory_itor->second;
        } else {
          throw std::runtime_error("Unknown pulsar ephemeris style: EPHSTYLE = " + eph_style);
        }
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
      const tip::Header & header(orbital_table.getHeader());

      // Try to read EPHSTYLE keyword to select a proper ephemeris factory.
      const IEphFactory<OrbitalEph> * factory(0);
      std::string eph_style;
      try {
        header["EPHSTYLE"].get(eph_style);
      } catch (const tip::TipException &) {
        // Use SimpleDdEph if EPHSTYLE keyword is missing.
        factory = &EphFactory<OrbitalEph, SimpleDdEph>::getFactory();
      }

      // Use a registered subclass of OrbitalEph, if EPHSTYLE keyword exists.
      if (0 == factory) {
        orbital_factory_cont_type::const_iterator factory_itor = m_orbital_factory_cont.find(eph_style);
        if (factory_itor != m_orbital_factory_cont.end()) {
          factory = factory_itor->second;
        } else {
          throw std::runtime_error("Unknown orbital ephemeris style: EPHSTYLE = " + eph_style);
        }
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
    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Get summary of extensions in input file.
    tip::FileSummary in_summary;
    file_svc.getFileSummary(in_file, in_summary);

    // Get first extension in input file.
    FileSummary::const_iterator ext_itor = in_summary.begin();

    // If input file was empty, do nothing.
    if (in_summary.end() == ext_itor) return;

    // Loop over all extensions in input file.
    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    int ext_number = 1;
    for (++ext_itor; ext_itor != in_summary.end(); ++ext_itor, ++ext_number) {
      const std::string ext_name = ext_itor->getExtId();

      // Find the right tip::Table object.
      // TODO: Generalize this part for new D4 FITS definition with multiple ephemerides models.
      TableCont * target_table_cont(0);
      Table * target_table(0);
      if ("SPIN_PARAMETERS" == ext_name) target_table_cont = &m_spin_par_table;
      else if ("ORBITAL_PARAMETERS" == ext_name) target_table_cont = &m_orbital_par_table;
      else if ("OBSERVERS" == ext_name) target_table_cont = &m_obs_code_table;
      else if ("ALTERNATIVE_NAMES" == ext_name) target_table_cont = &m_psr_name_table;
      else throw std::runtime_error("Unknown extension name: " + ext_name);
      if (target_table_cont->empty()) {
        throw std::logic_error("Could not find extension \"" + ext_name + "\" in file \"" + in_file + "\"");
      } else {
        target_table = *(target_table_cont->begin());
      }

      // Open input table by extension number, in order to work with multiple extensions of the same extension name.
      std::ostringstream oss;
      oss << ext_number;
      std::auto_ptr<const Table> in_table(file_svc.readTable(in_file, oss.str()));

      // Start at beginning of both tables.
      Table::ConstIterator in_itor = in_table->begin();
      Table::Iterator target_itor = target_table->end();

      // TODO: Resize output to match input. Requires tip to handle random access iterators.
      // target_table->setNumRecords(in_table->getNumRecords() + target_table->getNumRecords());

      // Copy all rows in table.
      for (; in_itor != in_table->end(); ++in_itor, ++target_itor) *target_itor = *in_itor;
    }
  }

  void PulsarDb::loadText(const std::string & in_file) {
    // Open input text file.
    std::ifstream in_table(in_file.c_str());

    // Make sure it worked.
    if (!in_table) throw std::runtime_error("Could not open file " + in_file);

    // Create a text buffer.
    static const size_t s_line_size = 2048;
    char buf[s_line_size];

    std::auto_ptr<Table> out_table(0);
    ParsedLine parsed_line;
    ParsedLine field_name;

    // Read table until extension name is found.
    std::string extension;
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

      if (parsed_line.size() > 1) throw std::runtime_error(std::string("Expected name of an extension on line ") + buf);
      else if (parsed_line.size() == 1) {
        extension = *parsed_line.begin();
        out_table.reset(m_tip_file.editTable(extension));
        // Get list of fields in the table.
        Table::FieldCont field = out_table->getValidFields();
        field_name.assign(field.begin(), field.end());
        break;
      }
    }
      
    std::map<std::string, size_t> found_map;

    // Read table until names of fields are found.
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

      if (parsed_line.size() > 0) {
        for (ParsedLine::iterator itor = parsed_line.begin(); itor != parsed_line.end(); ++itor) {
          // Convert to lowercase, because field names are all lowercase.
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
      throw std::runtime_error("File " + in_file + " does not have any columns in common with output extension " + extension);

    // Populate output table starting at the beginning.
    Table::Iterator out_itor = out_table->begin();

    // Read the rest of input table to populate fields.
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, parsed_line);

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
          (*out_itor)[*in_itor].set(parsed_line[found_itor->second]);
        } else {
          // TODO: Make this write an INDEF in the field (requires further changes to tip).
          // (*out_itor)[*in_itor].set([found_itor->second]);
        }
      }
      ++out_itor;
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
}
