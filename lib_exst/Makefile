SHELL        ?= /bin/sh
AR           ?= ar
CXX          ?= g++
ROOT_DIR     ?= .
TARGET_DIR   ?= $(ROOT_DIR)/lib
BUILD_DIR    ?= $(ROOT_DIR)/lib
CXXFLAGS     ?= -O3 -DNDEBUG
FLAGS        ?=
ifneq ($(FLAGS),)
-include $(FLAGS)
endif
INCLUDE      := -I $(ROOT_DIR)
ALL_CXXFLAGS := -std=c++11 -W -Wall $(INCLUDE) $(CXXFLAGS)
TARGET       := $(TARGET_DIR)/exst.a
SOURCES      := $(patsubst $(ROOT_DIR)/src/%.cpp,%.cpp,$(wildcard $(ROOT_DIR)/src/*.cpp))
DEPS         := $(patsubst %.cpp,$(BUILD_DIR)/%.dep, $(SOURCES))
OBJECTS      := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(SOURCES))
vpath %.cpp $(ROOT_DIR)/src
vpath %.h   $(ROOT_DIR)/exst
HTD_INCLUDE:=/usr/local/include/htd
HTD_LIB:= /usr/local/lib/libhtd.so

$(TARGET):	$(FLAGS) $(DEPS) $(OBJECTS)
	$(AR)  cvrs $(TARGET) $(OBJECTS)


ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

$(BUILD_DIR)/%.dep: %.cpp $(FLAGS)
	$(CXX) $(ALL_CXXFLAGS) -MT $(BUILD_DIR)/$*.o  -I $(HTD_INCLUDE) -MM $< > $@

$(BUILD_DIR)/%.o: %.cpp $(FLAGS)
	$(CXX) $(ALL_CXXFLAGS)  -I $(HTD_INCLUDE) -c $< -o $@


.PHONY:	clean

ifeq ($(findstring Windows,$(OS)),)
clean:
	@$(RM) $(TARGET)
	@$(RM) $(OBJECTS)
	@$(RM) $(DEPS)
else
RM := del /F /Q
clean:
	@$(RM) $(subst /,\,$(TARGET))  2>nul
	@$(RM) $(subst /,\,$(OBJECTS)) 2>nul
	@$(RM) $(subst /,\,$(DEPS))    2>nul
endif
