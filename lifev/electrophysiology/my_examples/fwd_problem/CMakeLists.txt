
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  fwd_problem
  SOURCES main.cpp 
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_IonicMinimalModel_pacingECG_data
  CREATE_SYMLINK
  SOURCE_FILES data MonodomainSolverParamList.xml SolverParamList2.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

