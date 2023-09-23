CC      = mpicxx
CFLAGS  = -std=c++14 -O3
SRCS    = $(wildcard *.cpp)
HDRS    = $(wildcard *.hpp)
OBJS    = $(subst .cpp,.o,$(SRCS))

UNAME   = $(shell uname -n)
ifeq ($(UNAME),cygnari-mbp.local)
  LAPACK_FLAG = -framework Accelerate
endif
ifeq ($(findstring derecho,$(UNAME)),derecho)
  LAPACK_FLAG = -qmkl
endif

all: driver single_rhs

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(LAPACK_FLAG) $(LDFLAGS) -c $< -o $@ $(LDLIBS)

driver: $(filter-out single_rhs.o,$(OBJS))
	$(CC) $(CFLAGS) $(LAPACK_FLAG) $(LDFLAGS) $^ -o $@ $(LDLIBS)

single_rhs: $(filter-out driver.o,$(OBJS))
	$(CC) $(CFLAGS) $(LAPACK_FLAG) $(LDFLAGS) $^ -o $@ $(LDLIBS)
