/** \file EphComputerApp.h
    \brief Interface for EphComputerApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputerApp_h
#define pulsarDb_EphComputerApp_h

#include "pulsarDb/PulsarToolApp.h"

#include "st_app/StApp.h"

#include "st_stream/StreamFormatter.h"

namespace pulsarDb {

  class EphComputer;

  /** \class EphComputerApp
      \brief Main application class for pulsar database access.
  */
  class EphComputerApp: public pulsarDb::PulsarToolApp {
    public:
      EphComputerApp();
      virtual ~EphComputerApp() throw();
      virtual void run();

    private:
      st_stream::StreamFormatter m_os;
  };

}
#endif
