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

#include "st_facilities/Env.h"

#include "tip/FileSummary.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace tip;

namespace pulsarDb {

  PulsarDb::PulsarDb(const std::string & in_file): m_in_file(in_file), m_spin_par_table(0), m_psr_name_table(0) {
    m_spin_par_table = IFileSvc::instance().editTable(in_file, "SPIN_PARAMETERS", "#row>0");
    m_psr_name_table = IFileSvc::instance().readTable(in_file, "ALTERNATIVE_NAMES");
  }

  PulsarDb::~PulsarDb() throw() { delete m_psr_name_table; delete m_spin_par_table; }

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

  const PulsarEph & PulsarDb::chooseEph(long double mjd, bool extrapolate) const {
    // Make sure ephemerides are already loaded.
    getEph();
    
    if (m_ephemerides.empty()) throw std::runtime_error("PulsarDb::chooseEph found no candidate ephemerides");

    if (extrapolate) {
      try {
        return chooseValidEph(mjd);
      } catch (const std::exception &) {
        return extrapolateEph(mjd);
      }
    }

    return chooseValidEph(mjd);
  }

  const PulsarEphCont & PulsarDb::getEph() const {
    // Refill this container.
    m_ephemerides.clear();

    // Iterate over current selection.
// TODO: why doesn't ConstIterator work here???
    for (Table::Iterator itor = m_spin_par_table->begin(); itor != m_spin_par_table->end(); ++itor) {
      // For convenience, get record from iterator.
      Table::ConstRecord & r(*itor);

      // Integral values must be looked up less conveniently.
      long epoch_int = 0;
      long t0_int = 0;
      r["EPOCH_INT"].get(epoch_int);
      r["T0GEO_INT"].get(t0_int);

      // Add the ephemeris to the container.
      m_ephemerides.push_back(PulsarEph(r["VALID_SINCE"].get(), r["VALID_UNTIL"].get(), epoch_int, r["EPOCH_FRAC"].get(),
        t0_int, r["T0GEO_FRAC"].get(), r["F0"].get(), r["F1"].get(), r["F2"].get()));
    }

    return m_ephemerides;
  }

  int PulsarDb::getNumEph() const {
    return m_spin_par_table->getNumRecords();
  }

  const PulsarEph & PulsarDb::chooseValidEph(long double mjd) const {
    PulsarEphCont::const_iterator candidate = m_ephemerides.end();

    for (PulsarEphCont::const_iterator itor = m_ephemerides.begin(); itor != m_ephemerides.end(); ++itor) {
      // See if this ephemeris contains the mjd.
      if (itor->m_since <= mjd && mjd < itor->m_until + 1.) {

        // See if this is the first candidate, which is automatically accepted.
        if (m_ephemerides.end() == candidate) {
          candidate = itor;
        // Otherwise, prefer the eph which starts later.
        } else if (itor->m_since > candidate->m_since) {
          candidate = itor;
        } else if (itor->m_since == candidate->m_since) {
          // The two start at the same time, so break the tie based on which one is valid longer.
          // Note that in a tie here, the one selected is the one appearing last in the sequence.
          if (itor->m_until >= candidate->m_until)
            candidate = itor;
        }
      }
    }

    // If no candidate was found, throw an exception.
    if (m_ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "PulsarDb::chooseValidEph could not find an ephemeris for time " << mjd;
      throw std::runtime_error(os.str());
    }
    return *candidate;
  }

  // Note this method assumes there is at least one candidate ephemeris.
  const PulsarEph & PulsarDb::extrapolateEph(long double mjd) const {
    PulsarEphCont::const_iterator candidate = m_ephemerides.begin();

    double diff = std::min(fabs(mjd - candidate->m_since), fabs(mjd - candidate->m_until));
    
    for (PulsarEphCont::const_iterator itor = m_ephemerides.begin(); itor != m_ephemerides.end(); ++itor) {
      double new_diff = std::min(fabs(mjd - itor->m_since), fabs(mjd - itor->m_until));
      if (new_diff < diff) {
        candidate = itor;
        diff = new_diff;
      }
    }

    return *candidate;
  }

}
