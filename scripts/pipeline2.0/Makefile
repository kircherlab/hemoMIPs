SHEDSKIN_LIBDIR=/fast/users/pkleiner_c/.local/lib/python2.7/site-packages/shedskin/lib
CC=g++
CCFLAGS=-O2 -march=native -Wno-deprecated $(CPPFLAGS) -I. -I${SHEDSKIN_LIBDIR} -g -fPIC -D__SS_BIND -I/fast/groups/ag_kircher/hemophilia/scripts/gc-7.6.4/bin/include/ -I/fast/users/pkleiner_c/miniconda2/include/python2.7 -I/fast/users/pkleiner_c/miniconda2/include/python2.7  
LFLAGS=-L/fast/groups/ag_kircher/hemophilia/scripts/gc-7.6.4/bin/lib/ -lgc -lpcre $(LDFLAGS) -shared -Xlinker -export-dynamic -lpthread -ldl -lutil -lm -lpython2.7 -lutil

CPPFILES=/fast/groups/ag_kircher/hemophilia/scripts/pipeline2.0/MergeTrimReads.cpp \
	${SHEDSKIN_LIBDIR}/sys.cpp \
	${SHEDSKIN_LIBDIR}/string.cpp \
	${SHEDSKIN_LIBDIR}/stat.cpp \
	${SHEDSKIN_LIBDIR}/re.cpp \
	${SHEDSKIN_LIBDIR}/os/path.cpp \
	${SHEDSKIN_LIBDIR}/os/__init__.cpp \
	${SHEDSKIN_LIBDIR}/math.cpp \
	${SHEDSKIN_LIBDIR}/builtin.cpp

HPPFILES=/fast/groups/ag_kircher/hemophilia/scripts/pipeline2.0/MergeTrimReads.hpp \
	${SHEDSKIN_LIBDIR}/sys.hpp \
	${SHEDSKIN_LIBDIR}/string.hpp \
	${SHEDSKIN_LIBDIR}/stat.hpp \
	${SHEDSKIN_LIBDIR}/re.hpp \
	${SHEDSKIN_LIBDIR}/os/path.hpp \
	${SHEDSKIN_LIBDIR}/os/__init__.hpp \
	${SHEDSKIN_LIBDIR}/math.hpp \
	${SHEDSKIN_LIBDIR}/builtin.hpp

all:	MergeTrimReads.so

MergeTrimReads.so:	$(CPPFILES) $(HPPFILES)
	$(CC)  $(CCFLAGS) $(CPPFILES) $(LFLAGS) -o MergeTrimReads.so

clean:
	rm -f MergeTrimReads.so

.PHONY: all clean

