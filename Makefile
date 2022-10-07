# ============================================================================ #
#   MAKEFILE FOR THIN FILM WITH FEEDBACK CONTROL                               #
# ============================================================================ #

# =============== #
#   Definitions   #
# =============== #
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
EXE_BE=film-benney
EXE_WR=film-wr
EXE_NS=film-ns
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out
DUMP_DIR=./dump
PLT_DIR=./plots

SRC=$(filter-out $(wildcard $(SRC_DIR)/_*.c) $(SRC_DIR)/$(EXE).c, $(wildcard $(SRC_DIR)/*.c))
SRC_BE=$(filter-out $(SRC_DIR)/$(EXE_WR).c $(SRC_DIR)/$(EXE_NS).c, $(SRC))
SRC_WR=$(filter-out $(SRC_DIR)/$(EXE_NS).c $(SRC_DIR)/$(EXE_BE).c, $(SRC))
SRC_NS=$(filter-out $(SRC_DIR)/$(EXE_BE).c $(SRC_DIR)/$(EXE_WR).c $(SRC_DIR)/$(EXE_NS).c, $(SRC)) $(SRC_DIR)/_$(EXE_NS).c

OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.c=.o)))
OBJ_BE=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_BE:.c=.o)))
OBJ_WR=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_WR:.c=.o)))
OBJ_NS=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_NS:.c=.o)))


# =============== #
#   Build Rules   #
# =============== #
# default is build all executables
.PHONY: films
films: $(EXE_BE) $(EXE_WR) $(EXE_NS)

# force rebuild of all files
.PHONY: all
all: clean films

# build Benney executable
$(EXE_BE): directories $(OBJ_BE)
	@printf "\033[1;32mLinking (Benney)\033[0m\n"
	$(LD) $(CFLAGS) -o $(EXE_BE) $(OBJ_BE) $(LDFLAGS) $(IFLAGS)

# build weighted residuals executable
$(EXE_WR): directories $(OBJ_WR)
	@printf "\033[1;32mLinking (weighted residuals)\033[0m\n"
	$(LD) $(CFLAGS) -o $(EXE_WR) $(OBJ_WR) $(LDFLAGS) $(IFLAGS)

# build Navier-Stokes executable
$(EXE_NS): directories source $(OBJ_NS)
	@printf "\033[1;32mLinking (Navier-Stokes)\033[0m\n"
	$(LD) $(CFLAGS) -o $(EXE_NS) $(OBJ_NS) $(LDFLAGS) $(IFLAGS)
	rm $(SRC_DIR)/_$(EXE_NS).c

# build basilisk binary
$(OBJ_DIR)/_$(EXE_NS).o: $(SRC_DIR)/_$(EXE_NS).c
	@printf "\033[1;36mBuilding %s\033[0m\n" $@
	$(CC) $(BFLAGS) -c -o $@ $< $(LDFLAGS) $(IFLAGS)

# build C binary
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@printf "\033[1;36mBuilding %s\033[0m\n" $@
	$(CC) $(CFLAGS) -c -o $@ $< $(LDFLAGS) $(IFLAGS)

# create required directories
.PHONY: directories
directories:
	@printf "\033[1;33mCreating output directories\033[0m\n"
	mkdir -p $(OBJ_DIR) $(OUT_DIR) $(PLT_DIR) $(DUMP_DIR)

# delete build files and executable
.PHONY: clean
clean:
	@printf "\033[1;31mCleaning\033[0m\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE_BE) ./$(EXE_WR) ./$(EXE_NS) $(SRC_DIR)/_$(EXE_NS).c

# delete build files, executable, and all outputs
.PHONY: deepclean
deepclean: clean
	@printf "\033[1;31mDeep cleaning\033[0m\n"
	rm -rf $(OUT_DIR)/* $(PLT_DIR)/* $(DUMP_DIR)/*

# generate a pure C source file from film.c
.PHONY: source
source: $(SRC_DIR)/$(EXE_NS).c
	@printf "\033[1;35mBuilding source\033[0m\n"
	cd $(SRC_DIR) && \
	$(BC) -source $(BFLAGS) -DPARALLEL $(LDFLAGS) $(IFLAGS) $(EXE_NS).c
