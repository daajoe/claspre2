file(GLOB LIBEXST_TEST_SRC *.cpp)
include_directories(${CMAKE_SOURCE_DIR}/lib_exst)
add_executable(exst_test ${LIBEXST_TEST_SRC})
target_link_libraries(exst_test lib_exst cppunit htd)

set_target_properties(exst_test PROPERTIES OUTPUT_NAME exst_test)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
