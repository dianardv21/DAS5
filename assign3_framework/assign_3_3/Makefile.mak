PROGNAME = assign3_3
TARNAME = assign3_3.tgz

RUNARGS = 1000 1000 # i_max t_max, increase this when testing on the DAS4!
NODES = 1 # How many DAS4 nodes
PROCESSESPERNODE = 1 # How many processes to spawn on one machine.

PRUNARGS = -v -np $(NODES) -$(PROCESSESPERNODE) \
		   -sge-script $$PRUN_ETC/prun-openmpi
CC = mpicc

WARNFLAGS = -Wall -Werror-implicit-function-declaration -Wshadow \
		  -Wstrict-prototypes -pedantic-errors
CFLAGS = -std=c99 -ggdb -O2 $(WARNFLAGS) -D_POSIX_C_SOURCE=200112
LFLAGS = -lm -lrt

# Do some substitution to get a list of .o files from the given .c files.
OBJFILES = $(patsubst %.c,%.o,$(SRCFILES))

.PHONY: all run runlocal plot clean dist todo

all: $(PROGNAME)

$(PROGNAME): $(OBJFILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

run: $(PROGNAME)
	prun $(PRUNARGS) $(PROGNAME) $(RUNARGS)

runlocal: $(PROGNAME)
	mpirun -n $(PROCESSESPERNODE) $(PROGNAME) $(RUNARGS)

