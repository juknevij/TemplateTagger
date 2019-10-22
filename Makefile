CC=g++
CFLAGS= -c -O3 -funroll-loops -Wall -pedantic `$(FASTJETLOCATION)/fastjet-config --cxxflags` 
LDFLAGS=
COMMON_SOURCES=
SOURCE1=short_example.cc
SOURCE2=01-basic.cc
SOURCE3=02-tagger.cc
SOURCE4=03-templdef.cc
SOURCE5=04-boosted_top.cc
SOURCE6=build_template.cc


COMMON_OBJECTS=$(COMMON_SOURCES:.cc=.o)
OBJECT1=$(SOURCE1:.cc=.o)
OBJECT2=$(SOURCE2:.cc=.o)
OBJECT3=$(SOURCE3:.cc=.o)
OBJECT4=$(SOURCE4:.cc=.o)
OBJECT5=$(SOURCE5:.cc=.o)
OBJECT6=$(SOURCE6:.cc=.o)

EXECUTABLE1=short_example
EXECUTABLE2=01-basic
EXECUTABLE3=02-tagger
EXECUTABLE4=03-templdef
EXECUTABLE5=04-boosted_top
EXECUTABLE6=build_template


FASTJETLOCATION=/home/juknevij/MonteCarlos/fj3/bin
LDLIBS= `$(FASTJETLOCATION)/fastjet-config  --libs` 


.PHONY: all clean short templ basic tagger templdef boostedtop

all: short_example 01-basic 02-tagger 03-templdef 04-boosted_top build_template

short: ${EXECUTABLE1}

basic: ${EXECUTABLE2}

tagger: ${EXECUTABLE3}

templdef: ${EXECUTABLE4}

boostedtop: ${EXECUTABLE5}

templ: ${EXECUTABLE6}


fj: CFLAGS += `$(FASTJETLOCATION)/fastjet-config --cxxflags` 

$(EXECUTABLE1): $(COMMON_OBJECTS) $(OBJECT1)
	$(CC) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(EXECUTABLE2): $(COMMON_OBJECTS) $(OBJECT2)
	$(CC) $(LDFLAGS) $(LDLIBS) $^ -o $@ 

$(EXECUTABLE3): $(COMMON_OBJECTS) $(OBJECT3)
	$(CC) $(LDFLAGS) $(LDLIBS)  $^ -o $@ 
	
$(EXECUTABLE4): $(COMMON_OBJECTS) $(OBJECT4)
	$(CC) $(LDFLAGS) $(LDLIBS)  $^ -o $@ 
	
$(EXECUTABLE5): $(COMMON_OBJECTS) $(OBJECT5)
	$(CC) $(LDFLAGS) $(LDLIBS)  $^ -o $@ 
	
$(EXECUTABLE6): $(COMMON_OBJECTS) $(OBJECT6)
	$(CC) $(LDFLAGS) $(LDLIBS)  $^ -o $@ 


	
.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o short_example 01-basic 02-tagger 03-templdef 04-boosted_top build_template *~
