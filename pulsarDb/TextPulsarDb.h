/** \file TextPulsarDb.h
    \brief Interface for TextPulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_TextPulsarDb_h
#define pulsarDb_TextPulsarDb_h

#include <map>
#include <string>
#include <vector>

#include "pulsarDb/PulsarDb.h"

namespace pulsarDb {

  /** \class TextPulsarDb
      \brief Abstraction providing access to pulsar ephemerides database from a text file.
  */
  class TextPulsarDb : public PulsarDb {
    public:
      typedef std::vector<std::string> Row;
      /** \brief Create a data base access object for the given ephermerides db file.
                 This opens a copy of the file in memory. The version on disk will be
                 unaffected.
          \param in_file The input file name.
      */
      TextPulsarDb(const std::string & in_file);

    protected:
      static const size_t s_line_size = 2048;

      virtual void readTextFile(const std::string & in_file);

      virtual void parseLine(const char * line, Row & row);
  };

}
#endif
