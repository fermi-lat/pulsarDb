/** \file PulsarDb.cxx
    \brief Implementation of the PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cctype>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/TimingModel.h"

#include "st_facilities/Env.h"

#include "tip/FileSummary.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace tip;

namespace pulsarDb {

  PulsarDb::PulsarDb(const std::string & in_file, bool bary_toa): m_in_file(in_file), m_spin_par_table(0), m_psr_name_table(0),
    m_bary_toa(bary_toa) {
    m_spin_par_table = IFileSvc::instance().editTable(in_file, "SPIN_PARAMETERS", "#row>0");
    m_psr_name_table = IFileSvc::instance().readTable(in_file, "ALTERNATIVE_NAMES");
  }

  PulsarDb::~PulsarDb() { delete m_psr_name_table; delete m_spin_par_table; }

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

  void PulsarDb::save(const std::string & out_file) {
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

    // Create output file.
    file_svc.createFile(out_file, tpl_file);

    // Get summary of extensions in output file.
    FileSummary summary;
    file_svc.getFileSummary(out_file, summary);

    // Loop over remaining extensions in output file.
    FileSummary::const_iterator ext_itor = summary.begin();

    // Skip the primary by incrementing the iterator in the first clause of the for loop.
    for (++ext_itor; ext_itor != summary.end(); ++ext_itor) {
      // Input table pointer.
      const Table * in_table = 0;

      // Autopointer which is used iff a new table is read.
      std::auto_ptr<const Table> auto_table(0);
      if (ext_itor->getExtId() == "SPIN_PARAMETERS") {
        // Just reuse the open, filtered, fully prepared table.
        in_table = m_spin_par_table;
      } else {
        // Read the extension from the input file.
        auto_table.reset(file_svc.readTable(m_in_file, ext_itor->getExtId()));
        in_table = auto_table.get();
      }

      // Open output table.
      std::auto_ptr<Table> out_table(file_svc.editTable(out_file, ext_itor->getExtId()));

      // Start at beginning of both tables.
      Table::ConstIterator in_itor = in_table->begin();
      Table::Iterator out_itor = out_table->begin();

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
        out_table->setNumRecords(in_table->getNumRecords());

        // Copy all rows in table.
        for (; in_itor != in_table->end(); ++in_itor, ++out_itor) *out_itor = *in_itor;
      }
    }
  }

  void PulsarDb::getEph(PulsarEphCont & cont) const {
    // Refill this container.
    cont.clear();

    // Determine which columns hold toa info.
    std::string toa_int_col = m_bary_toa ? "TOABARY_INT" : "TOAGEO_INT";
    std::string toa_frac_col = m_bary_toa ? "TOABARY_FRAC" : "TOAGEO_FRAC";

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
      long double epoch = (long double)(epoch_int) + epoch_frac;
      long double toa = (long double)(toa_int) + toa_frac;

      // One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
      // whereas the valid_until argument to DatabaseEph is the absolute cutoff.
      long double valid_until = r["VALID_UNTIL"].get() + 1.L;

      // Add the ephemeris to the container.
      cont.insertEph(DatabaseEph(m_model, r["VALID_SINCE"].get(), valid_until, epoch, toa, r["F0"].get(), r["F1"].get(),
        r["F2"].get()));
    }
  }

  int PulsarDb::getNumEph() const {
    return m_spin_par_table->getNumRecords();
  }

  const PulsarEph & PulsarEphCont::chooseEph(long double mjd, bool strict_validity) const {
    if (m_ephemerides.empty()) throw std::runtime_error("PulsarEphCont::chooseEph found no candidate ephemerides");

    if (!strict_validity) {
      try {
        return chooseFromValidEph(mjd);
      } catch (const std::exception &) {
        return chooseFromAllEph(mjd);
      }
    }

    return chooseFromValidEph(mjd);
  }

  const PulsarEph & PulsarEphCont::chooseFromValidEph(long double mjd) const {
    Cont_t::const_iterator candidate = m_ephemerides.end();

    for (Cont_t::const_iterator itor = m_ephemerides.begin(); itor != m_ephemerides.end(); ++itor) {
      // See if this ephemeris contains the mjd.
      if ((*itor)->valid_since() <= mjd && mjd < (*itor)->valid_until()) {

        // See if this is the first candidate, which is automatically accepted.
        if (m_ephemerides.end() == candidate) {
          candidate = itor;
        // Otherwise, prefer the eph which starts later.
        } else if ((*itor)->valid_since() > (*candidate)->valid_since()) {
          candidate = itor;
        } else if ((*itor)->valid_since() == (*candidate)->valid_since()) {
          // The two start at the same time, so break the tie based on which one is valid longer.
          // Note that in a tie here, the one selected is the one appearing last in the sequence.
          if ((*itor)->valid_until() >= (*candidate)->valid_until())
            candidate = itor;
        }
      }
    }

    // If no candidate was found, throw an exception.
    if (m_ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "PulsarEphCont::chooseFromValidEph could not find an ephemeris for time " << mjd;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const PulsarEph & PulsarEphCont::chooseFromAllEph(long double mjd) const {
    Cont_t::const_iterator candidate = m_ephemerides.begin();

    double diff = std::min(fabs(mjd - (*candidate)->valid_since()), fabs(mjd - (*candidate)->valid_until()));
    
    for (Cont_t::const_iterator itor = m_ephemerides.begin(); itor != m_ephemerides.end(); ++itor) {
      double new_diff = std::min(fabs(mjd - (*itor)->valid_since()), fabs(mjd - (*itor)->valid_until()));
      if (new_diff < diff) {
        candidate = itor;
        diff = new_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (m_ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "PulsarEphCont::chooseFromAllEph could not find an ephemeris for time " << mjd;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

}
