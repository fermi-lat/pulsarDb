/** \file TextPulsarDb.cxx
    \brief Implementation of the TextPulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <fstream>
#include <stdexcept>

#include "pulsarDb/TextPulsarDb.h"

#include "st_facilities/Env.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace tip;

namespace pulsarDb {

  TextPulsarDb::TextPulsarDb(const std::string & in_file): PulsarDb() {
    // Find data directory for this app.
    std::string data_dir = st_facilities::Env::getDataDir("pulsarDb");

    // Find template file.
    std::string tpl_file = st_facilities::Env::appendFileName(data_dir, "PulsarEph.tpl");

    // Alias for file service singleton.
    IFileSvc & file_svc(IFileSvc::instance());

    // Name input file after temporary memory file.
    m_in_file = "mem://text_pulsardb.fits";

    // Create a temporary FITS file in memory.
    file_svc.createFile(m_in_file, tpl_file);

    // Read the input text file and store it in the temporary FITS file in memory.
    readTextFile(in_file);

    // Read the temporary FITS file copy into the input table.
    loadTables(true);
  }

  void TextPulsarDb::readTextFile(const std::string & in_file) {
    // Open input text file.
    std::ifstream in_table(in_file.c_str());

    // Make sure it worked.
    if (!in_table) throw std::runtime_error("Could not open file " + in_file);

    char buf[s_line_size];

    std::auto_ptr<Table> out_table(0);
    Row row;
    Row field_name;

    // Read table until extension name is found.
    std::string extension;
    while (in_table) {
      // Get next line.
      in_table.getline(buf, s_line_size);

      // Parse it into fields.
      parseLine(buf, row);

      if (row.size() > 1) throw std::runtime_error(std::string("Expected name of an extension on line ") + buf);
      else if (row.size() == 1) {
        extension = *row.begin();
        out_table.reset(IFileSvc::instance().editTable(m_in_file, extension));
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
      parseLine(buf, row);

      if (row.size() > 0) {
        for (Row::iterator itor = row.begin(); itor != row.end(); ++itor) {
          // Convert to lowercase, because field names are all lowercase.
          for (std::string::iterator s_itor = itor->begin(); s_itor != itor->end(); ++s_itor) *s_itor = tolower(*s_itor);
        }

        // Populate input_column number as a function of output column name.
        for (Row::iterator itor = field_name.begin(); itor != field_name.end(); ++itor) {
          std::vector<std::string>::iterator found_itor = std::find(row.begin(), row.end(), *itor);
          if (row.end() != found_itor) {
            found_map[*itor] = found_itor - row.begin();
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
      parseLine(buf, row);

      if (row.size() == 0) continue;

      // Make sure the correct number of fields are found.
      if (row.size() != found_map.size()) {
        std::ostringstream os;
        os << "Line " << buf << " has " << row.size() << " row, not " << found_map.size() << ", as expected";
        throw std::runtime_error(os.str());
      }
  
// TODO: handle conversion errors
      // Populate the row by iterating over output columns.
      for (Row::iterator in_itor = field_name.begin(); in_itor != field_name.end(); ++in_itor) {
        // See if this column was found amongst the input table columns.
        std::map<std::string, size_t>::iterator found_itor = found_map.find(*in_itor);
        if (found_map.end() != found_itor) {
          // Found column, so look up the input table column number from the output column number.
          (*out_itor)[*in_itor].set(row[found_itor->second]);
        } else {
          // TODO: Make this write an INDEF in the field.
          // (*out_itor)[*in_itor].set([found_itor->second]);
        }
      }
      ++out_itor;
    }
  }

  void TextPulsarDb::parseLine(const char * line, Row & row) {
    row.clear();

    const char * begin = line;

    // Skip leading whitespace.
    while(isspace(*begin)) ++begin;

    // Handle the case of a completely blank line.
    if ('\0' == *begin) return;

    // Check for comment and return if it is one.
    if ('#' == *begin) return;

    // At this point there is at least one character, hence one field.
    row.resize(1);

    // Start with the first field.
    std::vector<std::string>::iterator current_field = row.begin();

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
        current_field = row.insert(row.end());

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
