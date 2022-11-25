# C++ compiler configuration
CXX             = g++
CXXFLAGS        = -O3

OBJS=		strtokenizer.o dataset.o utils.o model.o lda.o
MAIN=		lda

all:	$(OBJS)
	$(CXX) -g -o $(MAIN) $(OBJS) ${LIB} ${CXXFLAGS}

strtokenizer.o:	strtokenizer.h strtokenizer.cpp
	$(CXX) -c -o strtokenizer.o strtokenizer.cpp $(CXXFLAGS)

dataset.o:	dataset.h dataset.cpp
	$(CXX) -c -o dataset.o dataset.cpp $(CXXFLAGS)

utils.o:	utils.h utils.cpp
	$(CXX) -c -o utils.o utils.cpp $(CXXFLAGS)

model.o:	model.h model.cpp
	$(CXX) -c -o model.o model.cpp $(CXXFLAGS)

lda.o: lda.cpp
	$(CXX) -c -o lda.o lda.cpp $(CXXFLAGS)

test:


clean:
	rm $(OBJS)
	rm $(MAIN)
