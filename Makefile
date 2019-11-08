SRCS=util.cc random.cc pri_queue.cc block_file.cc \
	b_node.cc b_tree.cc qab_node.cc qab_tree.cc \
	qdafn.cc drusilla_select.cc rqalsh.cc rqalsh_star.cc \
	afn.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o rqalsh ${OBJS}

util.o: util.h

random.o: random.h

pri_queue.o: pri_queue.h

block_file.o: block_file.h

b_node.o: b_node.h

b_tree.o: b_tree.h

qab_node.o: qab_node.h

qab_tree.o: qab_tree.h

qdafn.o: qdafn.h

drusilla_select.o: drusilla_select.h

rqalsh.o: rqalsh.h

rqalsh_star.o: rqalsh_star.h

afn.o: afn.h

main.o:

clean:
	-rm ${OBJS}
