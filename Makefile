# ============================================================================ #
#   MAKEFILE FOR THIN FILM WITH FEEDBACK CONTROL                               #
# ============================================================================ #
UNAME := $(shell uname)

# ================ #
#   Definitiions   #
# ================ #
# Basilsik-C compiler (building and linking)
CC=qcc
LD=$(CC)

# C flags
CFLAGS=-O3

# required libraries (avoid openmp on macOS since it is slow)
LDFLAGS=-fopenmp

# file/folder names
EXE=film
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out
DUMP_DIR=./dump
PLT_DIR=./plots
SRC=$(wildcard $(SRC_DIR)/*.c)
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.c=.o)))


# =============== #
#   Build Rules   #
# =============== #
# Default build target
film: directories $(OBJ)
	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
	$(LD) $(CFLAGS) $(LDFLAGS) $(IFLAGS) -o $(EXE) $(OBJ)

# Build rule for binaries
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(CC) $(CFLAGS) $(LDFLAGS) $(IFLAGS) -c -o $@ $<

# Force rebuild of all files
.PHONY: all
all: clean film

# Create required directories
.PHONY: directories
directories:
	@printf "`tput bold``tput setaf 3`Creating output directories`tput sgr0`\n"
	mkdir -p $(OBJ_DIR) $(OUT_DIR) $(PLT_DIR) $(DUMP_DIR)

# Purge build files and executable
.PHONY: clean
clean:
	@printf "`tput bold``tput setaf 1`Cleaning`tput sgr0`\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE)

# Purge build files, executable, and all outputs
.PHONY: deepclean
deepclean:
	@printf "`tput bold``tput setaf 1`Deep cleaning`tput sgr0`\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE) $(OUT_DIR)/* $(PLT_DIR)/* $(DUMP_DIR)/*

# generate a pure C source file
.PHONY: source
source:
	@printf "`tput bold``tput setaf 5`Building source`tput sgr0`\n"
	cd $(SRC_DIR) && \
	$(CC) -source $(CFLAGS) -DPARALLEL $(LDFLAGS) $(IFLAGS) film.c
