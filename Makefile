# ============================================================================ #
#   MAKEFILE FOR THIN FILM WITH FEEDBACK CONTROL                               #
# ============================================================================ #

# ================ #
#   Definitiions   #
# ================ #
# Basilsik-C compiler (building and linking)
CC=qcc
LD=$(CC)

# C flags
CFLAGS=-O3 #-Wall -Wextra

# required libraries
IFLAGS=
LDFLAGS=-lm

# file/folder names
EXE=film
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out
PLT_DIR=./plots
SRC=$(wildcard $(SRC_DIR)/*.c)
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.c=.o)))


# =============== #
#   Build Rules   #
# =============== #
# Default build target
film: directories $(OBJ)
	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
	$(LD) $(CFLAGS) -o $(EXE) $(OBJ) $(LDFLAGS)

# Build rule for binaries
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(CC) $(CFLAGS) $(IFLAGS) -c -o $@ $<

# Force rebuild of all files
.PHONY: all
all: clean film

# Create required directories
.PHONY: directories
directories:
	@printf "`tput bold``tput setaf 3`Creating output directories`tput sgr0`\n"
	mkdir -p $(OBJ_DIR) $(OUT_DIR) $(PLT_DIR)

# Purge build files and executable
.PHONY: clean
clean:
	@printf "`tput bold``tput setaf 1`Cleaning`tput sgr0`\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE)

# Purge build files, executable, and all outputs
.PHONY: deepclean
deepclean:
	@printf "`tput bold``tput setaf 1`Deep cleaning`tput sgr0`\n"
	rm -rf $(OBJ_DIR)/*.o ./$(EXE) $(OUT_DIR)/* $(PLT_DIR)/* *.gif *.png
