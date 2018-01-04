#!/bin/bash
make
rm *.o


# Parameters of RQALSH:
#    -alg   (integer)   options of algorithms (0 - 3)
#    -n     (integer)   cardinality of the dataset
#    -qn    (integer)   number of queries
#    -d     (integer)   dimensionality of the dataset
#    -B     (integer)   page size
#    -beta  (real)      the percentage of false positive
#    -delta (real)      error probability
#    -c     (real)      approximation ratio (c > 1)
#    -ds    (string)    file path of the dataset
#    -qs    (string)    file path of the query set
#    -ts    (string)    file path of the ground truth set
#    -df    (string)    data folder to store new format of data
#    -of    (string)    output folder to store the results of rqalsh
#
# The options of algorithms (-alg) are:
#    0 - Ground-Truth
#        Parameters: -alg 0 -n -qn -d -ds -qs -ts -of
#
#    1 - Indexing
#        Parameters: -alg 1 -n -d -B -beta -delta -c -ds -df -of
#
#    2 - RQALSH
#        Parameters: -alg 2 -qn -d -qs -ts -df -of
#
#    3 - Linear Scan
#        Parameters: -alg 3 -n -qn -d -B -qs -ts -df -of
#
# NOTE: Each parameter is required to be separated by one space

n=30159
d=5408
B=65536
dname=P53

qn=1000
beta=100
delta=0.49
c=2.0
dPath=./data/${dname}/${dname}
dFolder=./data/${dname}/
oFolder=./results/${dname}/


### Ground Truth ###
./rqalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.fn2.0 -of ${oFolder}


### Indexing ###
./rqalsh -alg 1 -n ${n} -d ${d} -B ${B} -beta ${beta} -delta ${delta} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}


### k-FN search ###
./rqalsh -alg 2 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}ID.fn2.0 -df ${dFolder} -of ${oFolder}


### Brute-Force Search ###
./rqalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -qs ${dPath}.q -ts ${dPath}ID.fn2.0 -df ${dFolder} -of ${oFolder}






