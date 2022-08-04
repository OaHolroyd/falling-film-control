# ============================================================================ #
#   MAKEFILE FOR THIN FILM WITH FEEDBACK CONTROL                               #
# ============================================================================ #

# ================ #
#   Definitiions   #
# ================ #
# Compilers
BC=qcc # Basilisk
CC=gcc # C99
LD=$(CC) # linker

# C flags
BFLAGS=-O3 -Wall
CFLAGS=-O3 -Wall -Wextra -pedantic -Wno-unused-parameter -Wshadow \
       -Waggregate-return -Wbad-function-cast -Wcast-align -Wcast-qual \
       -Wfloat-equal -Wformat=2 -Wlogical-op -Wmissing-include-dirs \
       -Wnested-externs -Wpointer-arith -Wconversion -Wno-sign-conversion \
       -Wredundant-decls -Wsequence-point -Wstrict-prototypes -Wswitch -Wundef \
       -Wunused-but-set-parameter -Wwrite-strings
# CFLAGS=-O0 -g -fbounds-check -fsanitize=address -fsanitize=bounds -fsanitize=bounds-strict

# required libraries
LDFLAGS=-fopenmp -lm -llapacke

# file/folder names
EXE=film
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out
DUMP_DIR=./dump
PLT_DIR=./plots

SRC=$(filter-out $(SRC_DIR)/_$(EXE).c $(SRC_DIR)/$(EXE).c, $(wildcard $(SRC_DIR)/*.c)) $(SRC_DIR)/_$(EXE).c
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.c=.o)))


# =============== #
#   Build Rules   #
# =============== #
# Default build target
film: directories source link
	rm $(SRC_DIR)/_$(EXE).c

# link to make executable
.PHONY: link
link: $(OBJ)
	@printf "\033[1;32mLinking\033[0m\n"
	$(LD) $(CFLAGS) -o $(EXE) $(OBJ) $(LDFLAGS) $(IFLAGS)

# Build rule for basilisk binary
$(OBJ_DIR)/_film.o: $(SRC_DIR)/_film.c
	@printf "\033[1;36mBuilding %s\033[0m\n" $@
	$(CC) $(BFLAGS) -c -o $@ $< $(LDFLAGS) $(IFLAGS)

# Build rule for binaries
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@printf "\033[1;36mBuilding %s\033[0m\n" $@
	$(CC) $(CFLAGS) -c -o $@ $< $(LDFLAGS) $(IFLAGS)

# Force rebuild of all files
.PHONY: all
all: clean film

# Create required directories
.PHONY: directories
directories:
	@printf "\033[1;33mCreating output directories\033[0m\n"
	mkdir -p $(OBJ_DIR) $(OUT_DIR) $(PLT_DIR) $(DUMP_DIR)

# Purge build files and executable
.PHONY: clean
clean:
	@printf "\033[1;31mCleaning\033[0m\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE) $(SRC_DIR)/_$(EXE).c

# Purge build files, executable, and all outputs
.PHONY: deepclean
deepclean:
	@printf "\033[1;31mDeep cleaning\033[0m\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE) $(OUT_DIR)/* $(PLT_DIR)/* $(DUMP_DIR)/*

# generate a pure C source file from film.c
.PHONY: source
source: $(SRC_DIR)/$(EXE).c
	@printf "\033[1;35mBuilding source\033[0m\n"
	cd $(SRC_DIR) && \
	$(BC) -source $(BFLAGS) -DPARALLEL $(LDFLAGS) $(IFLAGS) $(EXE).c
