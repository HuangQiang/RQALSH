SRCS=pri_queue.cpp util.cpp random.cpp block_file.cpp b_node.cpp b_tree.cpp \
	rqalsh.cpp afn.cpp main.cpp
OBJS=$(SRCS:.cpp=.o)

CXX?=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: $(OBJS)
	$(CXX) -o rqalsh $(OBJS)

pri_queue.o: pri_queue.h

util.o: util.h

random.o: random.h

block_file.o: block_file.h

b_node.o: b_node.h

b_tree.o: b_tree.h

rqalsh.o: rqalsh.h

afn.o: afn.h

main.o:

clean:
	-rm $(OBJS) rqalsh
