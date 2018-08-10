CXXFLAGS += -g -Wall -O3 --std=c++11

OBJS = calculateDxy calculatePolymorphism listPolyDivSites nonOverlappingWindows sitePatterns

.PHONY: all,clean

all: $(OBJS)

clean:
	rm $(OBJS)
