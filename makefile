CXX       := mpicxx
CXXFLAGS  := -std=c++14 -Ofast -Iinclude

SRC_DIR   := src
INC_DIR   := include
OBJ_DIR   := build
BIN_DIR   := bin

SRCS      := $(wildcard $(SRC_DIR)/*.cpp)
HDRS      := $(wildcard $(INC_DIR)/*.hpp)
OBJS      := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# filter files with main() function for inclusion in library
LIBO      := $(filter-out $(OBJ_DIR)/driver.o,$(filter-out $(OBJ_DIR)/single_rhs.o,$(OBJS)))

# get lapack linker flags
UNAME     := $(shell uname -n -o -s)
ifeq ($(findstring Darwin,$(UNAME)),Darwin)
  CXXFLAGS := $(CXXFLAGS) -framework Accelerate
endif
ifeq ($(findstring derecho,$(UNAME)),derecho)
  CXXFLAGS := $(CXXFLAGS) -qmkl
endif

all: directories $(BIN_DIR)/driver $(BIN_DIR)/single_rhs $(BIN_DIR)/mylib.so

$(BIN_DIR)/driver: $(filter-out $(OBJ_DIR)/single_rhs.o,$(OBJS))
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/single_rhs: $(filter-out $(OBJ_DIR)/driver.o,$(OBJS))
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR)/mylib.so: $(LIBO)
	$(CXX) $(CXXFLAGS) -shared -fPIC -o $@ $^

.PHONY: clean
clean:
	rm -r $(OBJ_DIR)
	rm $(BIN_DIR)/driver $(BIN_DIR)/single_rhs $(BIN_DIR)/mylib.so

.PHONY: directories
directories: $(BIN_DIR) $(OBJ_DIR)
	mkdir -p $(BIN_DIR)
	mkdir -p $(OBJ_DIR)
