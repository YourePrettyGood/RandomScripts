CXXFLAGS += -g -Wall -O3 --std=c++11

OBJS = calculateDxy calculatePolymorphism listPolyDivSites nonOverlappingWindows softmaskFromHardmask sitePatterns

all: $(OBJS)

clean:
	rm $(OBJS)
