CC      = mpicxx
CFLAGS  = -std=c++14 -O3 -framework Accelerate
SRCS    = $(wildcard *.cpp)
HDRS    = $(wildcard *.hpp)
OBJS    = $(subst .cpp,.o,$(SRCS))

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@ $(LDLIBS)

driver: $(filter-out single_rhs.o,$(OBJS))
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

single_rhs: $(filter-out driver.o,$(OBJS))
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

all: driver single_rhs
