##### Library ######
add_library(
  pulsarDb STATIC
  src/BtModelEph.cxx
  src/Ell1ModelEph.cxx
  src/EphChooser.cxx
  src/EphComputer.cxx
  src/EphComputerApp.cxx
  src/EphStatus.cxx
  src/FrequencyEph.cxx
  src/HighPrecisionEph.cxx
  src/MssModelEph.cxx
  src/OrbitalEph.cxx
  src/PdotCanceler.cxx
  src/PeriodEph.cxx
  src/PulsarDb.cxx
  src/PulsarDbApp.cxx
  src/PulsarEph.cxx
  src/PulsarToolApp.cxx
  src/SimpleDdEph.cxx
)
target_link_libraries(
  pulsarDb
  PUBLIC st_app st_stream tip timeSystem
  PRIVATE st_facilities
)
target_include_directories(
  pulsarDb PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

###### Executables ######
add_executable(gtephem src/gtephem/gtephem.cxx)
add_executable(gtpulsardb src/gtpulsardb/gtpulsardb.cxx)

target_link_libraries(gtephem PRIVATE pulsarDb)
target_link_libraries(gtpulsardb PRIVATE pulsarDb)

###### Tests ######
add_executable(test_pulsarDb src/test/test_pulsarDb.cxx)
target_link_libraries(test_pulsarDb PRIVATE pulsarDb)


###############################################################
# Installation
###############################################################
install(DIRECTORY pulsarDb DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})
install(DIRECTORY data/outref DESTINATION ${FERMI_INSTALL_DATADIR}/pulsarDb)
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_REFDATADIR}/pulsarDb)

install(
  TARGETS pulsarDb gtephem gtpulsardb test_pulsarDb
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
