/** \file PulsarDb.cxx
    \brief Implementation of the PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <ctime>
#include <cctype>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "pulsarDb/CanonicalTime.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/TimingModel.h"

#include "st_facilities/Env.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace tip;

namespace pulsarDb {

  PulsarDb::PulsarDb(const std::string & in_file, bool edit_in_place): m_in_file(in_file), m_summary(), m_table(),
    m_spin_par_table(0), m_orbital_par_table(0), m_psr_name_table(0) {
    loadTables(edit_in_place);
  }

  PulsarDb::~PulsarDb() {
    // TODO: CHANge Tip so it has reverse_iterator, then use it to delete in the correct order.
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
    filter(os.str());
  }

  void PulsarDb::filterName(const std::string & name) {
    std::string expression;
    std::string up_name = name;

    // Convert to uppercase.
    for (std::string::iterator s_itor = up_name.begin(); s_itor != up_name.end(); ++s_itor) *s_itor = toupper(*s_itor);

    // Check for the (limited) wildcard.
    if (up_name == "ALL" || up_name == "ANY") return;

    // Look up the given name in the extension containing alternate names.
    for (Table::ConstIterator itor = m_psr_name_table->begin(); itor != m_psr_name_table->end(); ++itor) {
      // See if there is an alternate name.
      std::string alt_name;
      (*itor)["ALTNAME"].get(alt_name);

      // Convert to uppercase.
      for (std::string::iterator s_itor = alt_name.begin(); s_itor != alt_name.end(); ++s_itor) *s_itor = toupper(*s_itor);

      if (alt_name == up_name) {
        // Found a match.
        std::string real_name;
        (*itor)["PSRNAME"].get(real_name);
        expression = "(PSRNAME==\"" + real_name + "\")||";
      }
    }

    // If the all-upper case name is different from the original name, add it.
    if (up_name != name) {
      expression += "(PSRNAME==\"" + up_name + "\")||";
    }

    // Finally (just to be sure) add the original name given.
    expression += "(PSRNAME==\"" + name + "\")";

    // OK, do the filter for real on the spin parameters table.
    filter(expression);
  }

  void PulsarDb::save(const std::string & out_file, const std::string & tpl_file) const {
    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Find data directory for this app.
    std::string data_dir = st_facilities::Env::getDataDir("pulsarDb");

    if (!file_svc.fileExists(out_file)) {
      // Create output file.
      file_svc.createFile(out_file, tpl_file);
    }

    // Loop over all extensions in input file.
    FileSummary::const_iterator ext_itor = m_summary.begin();

    // TODO: should this be an error?
    // If input file was empty, just return.
    if (m_summary.end() == ext_itor) return;

    // Update keywords only in primary extension.
    std::auto_ptr<Extension> primary(file_svc.editExtension(out_file, ext_itor->getExtId()));

    // Update standard keywords.
    updateKeywords(*primary);

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

#if 0
// TODO: This was originally a save() method which removed unneeded information. If we want such
// a "pruning" capability, it could be modified to prune the input file. Then it would work similarly
// to the filter methods, which also affect the input file.
  void PulsarDb::prune() {
    // Get alias to tip's file service, for brevity.
    IFileSvc & file_svc(IFileSvc::instance());

    // Find data directory for this app.
    std::string data_dir = st_facilities::Env::getDataDir("pulsarDb");

    // Find template file.
    std::string tpl_file = st_facilities::Env::appendFileName(data_dir, "PulsarEph.tpl");

    // Set of pulsar names and observer codes.
    std::set<std::string> pulsar_names;
    std::set<std::string> observer_codes;

    // Compose a set of all pulsar names and observer codes. This is used below
    // to filter the other extensions in the file.
    for (Table::Iterator in_itor = m_spin_par_table->begin(); in_itor != m_spin_par_table->end(); ++in_itor) {
      std::string tmp;

      (*in_itor)["PSRNAME"].get(tmp);
      pulsar_names.insert(tmp);

      (*in_itor)["OBSERVER_CODE"].get(tmp);
      observer_codes.insert(tmp);
    }

    if (!append || !file_svc.fileExists(out_file)) {
      // Create output file.
      file_svc.createFile(out_file, tpl_file);
    }

    // Loop over remaining extensions in output file.
    FileSummary::const_iterator ext_itor = m_summary.begin();

    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    for (++ext_itor; ext_itor != m_summary.end(); ++ext_itor) {
      TableCont::const_iterator table_itor = m_table.find(ext_itor->getExtId());
      if (m_table.end() == table_itor)
        throw std::logic_error("Could not find extension \"" + ext_itor->getExtId() + "\" in file \"" + m_in_file + "\"");

      // Input table pointer.
      const Table * in_table = table_itor->second;

      // Open output table.
      std::auto_ptr<Table> out_table(file_svc.editTable(out_file, ext_itor->getExtId()));

      // Start at beginning of both tables.
      Table::ConstIterator in_itor = in_table->begin();
      Table::Iterator out_itor = out_table->end();

      if (ext_itor->getExtId() == "ORBITAL_PARAMETERS" || ext_itor->getExtId() == "ALTERNATIVE_NAMES") {
        // Copy only records which match the selected pulsar name.
        for (; in_itor != in_table->end(); ++in_itor) {
          std::string match;

          // Get the pulsar name.
          (*in_itor)["PSRNAME"].get(match);

          // Only copy the row if the name is in the set of names.
          if (pulsar_names.end() != pulsar_names.find(match)) {
            *out_itor = *in_itor;
            ++out_itor;
          }
        }
      } else if (ext_itor->getExtId() == "OBSERVERS") {
        // Copy only records which match the selected observer codes.
        for (; in_itor != in_table->end(); ++in_itor) {
          std::string match;
          (*in_itor)["OBSERVER_CODE"].get(match);
          if (observer_codes.end() != observer_codes.find(match)) {
            *out_itor = *in_itor;
            ++out_itor;
          }
        }
      } else {
        // Resize output to match input.
//        out_table->setNumRecords(in_table->getNumRecords());

        // Copy all rows in table.
        for (; in_itor != in_table->end(); ++in_itor, ++out_itor) *out_itor = *in_itor;
      }
    }
  }
#endif

  void PulsarDb::getEph(PulsarEphCont & cont) const {
    // Refill this container.
    cont.clear();

    // Determine which columns hold toa info.
    std::string toa_int_col = "TOABARY_INT";
    std::string toa_frac_col = "TOABARY_FRAC";

    cont.reserve(m_spin_par_table->getNumRecords());

    // Iterate over current selection.
// TODO: why doesn't ConstIterator work here???
    for (Table::Iterator itor = m_spin_par_table->begin(); itor != m_spin_par_table->end(); ++itor) {
      // For convenience, get record from iterator.
      Table::ConstRecord & r(*itor);

      // Epoch and toa are split into int and frac parts.
      long epoch_int = 0;
      double epoch_frac = 0.;
      long toa_int = 0;
      double toa_frac = 0.;

      // Read the separate parts from the file.
      r["EPOCH_INT"].get(epoch_int);
      r["EPOCH_FRAC"].get(epoch_frac);
      r[toa_int_col].get(toa_int);
      r[toa_frac_col].get(toa_frac);

      // Combine separate parts of epoch and toa to get long double values.
      TdbTime epoch = (long double)(epoch_int) + epoch_frac;
      long double toa = (long double)(toa_int) + toa_frac;

      // TODO confirm that db uses tdb!
      TdbTime valid_since = r["VALID_SINCE"].get();

      // One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
      // whereas the valid_until argument to DatabaseEph is the absolute cutoff.
      TdbTime valid_until = r["VALID_UNTIL"].get() + 1.L;

      // Add the ephemeris to the container.
      cont.push_back(DatabaseEph(valid_since, valid_until, epoch, TdbTime(toa),
        r["F0"].get(), r["F1"].get(), r["F2"].get()).clone());
    }
  }

  void PulsarDb::getEph(OrbitalEphCont & cont) const {
    // Empty container then refill it.
    cont.clear();

    if (0 != m_orbital_par_table) {
      cont.reserve(m_orbital_par_table->getNumRecords());

      for (Table::ConstIterator itor = m_orbital_par_table->begin(); itor != m_orbital_par_table->end(); ++itor) {
        // For convenience, get record from iterator.
        Table::ConstRecord & r(*itor);
        
        double par[] = {
          r["PB"].get(),
          r["PBDOT"].get(),
          r["A1"].get(),
          r["XDOT"].get(),
          r["ECC"].get(),
          r["ECCDOT"].get(),
          r["OM"].get(),
          r["OMDOT"].get(),
          r["T0"].get(),
          r["GAMMA"].get(),
          r["SHAPIRO_R"].get(),
          r["SHAPIRO_S"].get()
        };

        // Handle any INDEFs.
        for (size_t index = 0; index != sizeof(par) / sizeof(double); ++index) {
          if (0 != isnan(par[index])) {
            switch (index) {
              case PB:
              case A1:
              case ECC:
              case OM:
              case T0:
                throw std::runtime_error("PulsarDb::getEph(): invalid orbital ephemeris");
                break;
              case PBDOT:
              case XDOT:
              case ECCDOT:
              case OMDOT:
              case GAMMA:
              case SHAPIRO_R:
              case SHAPIRO_S:
              default:
                par[index] = 0.;
                break;
            }
          }
        }
        cont.push_back(new OrbitalEph(par));
      }
    }
  }

  int PulsarDb::getNumEph() const {
    return m_spin_par_table->getNumRecords();
  }

  PulsarDb::PulsarDb(): m_in_file(), m_summary(), m_table(), m_spin_par_table(0), m_orbital_par_table(0), m_psr_name_table(0) {
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
    
    // Find pulsar name table.
    itor = m_table.find("ALTERNATIVE_NAMES");
    if (m_table.end() == itor) throw std::runtime_error("Could not find ALTERNATIVE_NAMES table in file " + m_in_file);
    m_psr_name_table = itor->second;
  }

}
