CXXFLAGS += -g -Wall -O3 --std=c++11

#OBJS = calculateDxy calculatePolymorphism listPolyDivSites nonOverlappingWindows softmaskFromHardmask sitePatterns
#Exclude sitePatterns for the time being, as it sometimes has a compiling error
OBJS = calculateDxy calculatePolymorphism listPolyDivSites nonOverlappingWindows softmaskFromHardmask

all: $(OBJS)

clean:
	rm $(OBJS)
