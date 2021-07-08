# ------------------------------------------------------------------------------
#  Compile with C++ 11
# ------------------------------------------------------------------------------
SRCS=util.cc random.cc pri_queue.cc block_file.cc b_node.cc b_tree.cc \
	qab_node.cc qab_tree.cc qdafn.cc drusilla_select.cc rqalsh.cc \
	rqalsh_star.cc afn.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o rqalsh ${OBJS}

clean:
	-rm ${OBJS} rqalsh
