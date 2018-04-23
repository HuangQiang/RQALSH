#!/bin/bash
make
rm *.o

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=Sift
n=999000
qn=1000
d=128
B=4096
beta=100
delta=0.49
c=2.0

dPath=./data/${dname}/${dname}
dFolder=./data/${dname}/
oFolder=./results${c}/${dname}/

# ------------------------------------------------------------------------------
#  Ground Truth 
# ------------------------------------------------------------------------------
./rqalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
  -ts ${dPath}.fn2.0

# ------------------------------------------------------------------------------
#  RQALSH_Star
# ------------------------------------------------------------------------------
L_list=(2 3 4 5 6 10 15 20 30 50 60 75 100 150) 
M_list=(150 100 75 60 50 30 20 15 10 6 5 4 3 2)
length=`expr ${#L_list[*]} - 1`

for j in $(seq 0 ${length})
do 
  L=${L_list[j]}
  M=${M_list[j]}

  ./rqalsh -alg 1 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -beta ${beta} \
    -delta ${delta} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

  ./rqalsh -alg 2 -qn ${qn} -d ${d} -L ${L} -M ${M} -qs ${dPath}.q \
    -ts ${dPath}.fn2.0 -df ${dFolder} -of ${oFolder}
done

# ------------------------------------------------------------------------------
#  RQALSH
# ------------------------------------------------------------------------------
./rqalsh -alg 3 -n ${n} -d ${d} -B ${B} -beta ${beta} -delta ${delta} -c ${c} \
  -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./rqalsh -alg 4 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn2.0 \
  -df ${dFolder} -of ${oFolder}

# ------------------------------------------------------------------------------
#  Linear Scan
# ------------------------------------------------------------------------------
./rqalsh -alg 5 -n ${n} -qn ${qn} -d ${d} -B ${B} -qs ${dPath}.q \
  -ts ${dPath}.fn2.0 -df ${dFolder} -of ${oFolder}
