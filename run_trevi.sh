#!/bin/bash
make
rm *.o

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=Trevi
n=99900
qn=1000
d=4096
B=65536
beta=100
delta=0.49
c=2.0

dPath=../../data/${dname}/${dname}
dFolder=../../data/${dname}/
oFolder=../../results${c}/${dname}/

# ------------------------------------------------------------------------------
#  Ground Truth 
# ------------------------------------------------------------------------------
./rqalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
  -ts ${dPath}.fn2.0

# ------------------------------------------------------------------------------
#  RQALSH_Star
# ------------------------------------------------------------------------------
L_list=(2 3 5 6 10 15) 
M_list=(15 10 6 5 3 2)
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
