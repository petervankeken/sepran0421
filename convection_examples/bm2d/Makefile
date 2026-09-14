# Find all source files, create a list of corresponding object files
SRCS=$(wildcard *.f90)
OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard *mod.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))
#
# # Compiler/Linker settings
PROGRAM = bm2d
PRG_OBJ = $(PROGRAM).o
#
# # make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	compile $<
#
# Linker
$(PROGRAM) : $(OBJS)
	ccompile csecond.c
	seplink $(PROGRAM)
	cp $(PROGRAM) $(SEPBINDIR)/bm2d

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

#
# # Blocks and allocations depends on shared
# mod_blocks.o mod_allocations.o : mod_shared.o
#
