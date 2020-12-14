#template taken from https://github.com/TheNetAdmin/Makefile-Templates

# tool macros
CC := g++
CCFLAGS := -fopenmp
DBGFLAGS := -g
CCOBJFLAGS := $(CCFLAGS) -c  -fPIC

# path macros
BIN_PATH := build/bin
OBJ_PATH := build/obj
SRC_PATH := src/cpp
DBG_PATH := build/debug

# lib macros
EIGEN_LIB = /usr/include/eigen3

# compile macros
TARGET_NAME := OrthOpt
ifeq ($(OS),Windows_NT)
	TARGET_NAME := $(addsuffix .exe,$(TARGET_NAME))
endif
TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

# src files & obj files
SRC := $(foreach x, $(SRC_PATH), $(wildcard $(addprefix $(x)/*,.c*)))
OBJ := $(addprefix $(OBJ_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC)))))
OBJ_DEBUG := $(addprefix $(DBG_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC)))))

# clean files list
DISTCLEAN_LIST := $(OBJ) \
                  $(OBJ_DEBUG)
CLEAN_LIST := $(TARGET) \
			  $(TARGET_DEBUG) \
			  $(DISTCLEAN_LIST)

# default rule
default: makedir all

# non-phony targets
$(TARGET): $(OBJ)
	$(CC) $(CCFLAGS) -o $@ $(OBJ)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c*
	$(CC) $(CCOBJFLAGS) -I$(EIGEN_LIB) -o $@ $<

$(DBG_PATH)/%.o: $(SRC_PATH)/%.c*
	$(CC) $(CCOBJFLAGS) -I$(EIGEN_LIB) $(DBGFLAGS) -o $@ $<

$(TARGET_DEBUG): $(OBJ_DEBUG)
	$(CC) $(CCFLAGS) $(DBGFLAGS) $(OBJ_DEBUG) -o $@

# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(BIN_PATH) $(OBJ_PATH) $(DBG_PATH)

.PHONY: all
all: $(TARGET)

.PHONY: debug
debug: $(TARGET_DEBUG)

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(CLEAN_LIST)

.PHONY: distclean
distclean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(DISTCLEAN_LIST)
