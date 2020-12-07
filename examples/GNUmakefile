CXXFLAGS += $(shell pkg-config --cflags libpdb-redo)
LIBS += $(shell pkg-config --libs libpdb-redo)

all: example

OBJECTS = example.o example_rsrc.o

example.o: example.cpp
	$(CXX) -c -std=c++17 -o $@ $? $(CXXFLAGS)

bond-info.bin:
	ln -s ../rsrc/bond-info.bin.gz $@

example_rsrc.o: bond-info.bin
	mrc -o $@ $^

example: $(OBJECTS)
	$(CXX) -std=c++17 -o $@ $(OBJECTS) $(CXXFLAGS) $(LIBS)
