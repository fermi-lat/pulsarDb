/** \file EphComputerApp.h
    \brief Interface for EphComputerApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputerApp_h
#define pulsarDb_EphComputerApp_h

#include "pulsarDb/PulsarToolApp.h"

#include "st_stream/StreamFormatter.h"

namespace pulsarDb {

  class EphComputer;

  /** \class EphComputerApp
      \brief Main application class for ephemeris computations using pulsar database.
  */
  class EphComputerApp: public pulsarDb::PulsarToolApp {
    public:
      /// \brief Construct an EphComputerApp object.
      EphComputerApp();

      /// \brief Destruct this EphComputerApp object.
      virtual ~EphComputerApp() throw();

      /// \brief Run the application.
      virtual void runApp();

    private:
      st_stream::StreamFormatter m_os;
  };

}
#endif
