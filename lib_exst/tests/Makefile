SHELL = /bin/sh
CXX = g++
CXXFLAGS = -std=c++11 -W -Wall
INCLUDE = -I../ -I.
ALL_CXXFLAGS = $(INCLUDE) $(CXXFLAGS)
TARGET = clasp_test
CPPUNIT_INCLUDE:=/usr/local/include/cppunit
CPPUNIT_LIB:= /usr/local/lib/libcppunit.a
HTD_INCLUDE:=/usr/local/include/htd
HTD_LIB:= /usr/local/lib/libhtd.so
exst_lib:=../lib/exst.a
sources:=$(wildcard *.cpp)
objects:=$(patsubst %.cpp,%.o,$(sources))
deps:=$(patsubst %.cpp,%.dep,$(sources))

$(TARGET):	$(exst_lib) $(objects) $(deps)
	$(CXX) $(ALL_CXXFLAGS) $(objects) $(exst_lib) $(CPPUNIT_LIB) $(HTD_LIB) -o $(TARGET)

ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),distclean)
    -include $(deps)
  endif
endif

%.o:	%.cpp
	$(CXX) $(ALL_CXXFLAGS) -I $(CPPUNIT_INCLUDE) -I $(HTD_INCLUDE) -c $< -o $@

%.dep: %.cpp
	$(CXX) -std=c++11 $(INCLUDE) -I $(CPPUNIT_INCLUDE) -I $(HTD_INCLUDE) -MM $< > $@

$(exst_lib):
	$(MAKE)  -C ../src

.PHONY:	clean distclean
clean:
	@$(RM) $(TARGET)
	@$(RM) $(objects)
	@$(RM) $(deps)

distclean: clean


