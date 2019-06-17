CXXFLAGS += -g -Wall -O3 --std=c++11

OBJS = calculateDxy calculatePolymorphism listPolyDivSites nonOverlappingWindows softmaskFromHardmask sitePatterns

all: $(OBJS)

#test: $(OBJS)
#	./test_scripts.sh

clean:
	rm $(OBJS)
