#!/bin/bash
make
make clean

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=Mnist
n=59000
qn=1000
d=50
B=4096
c=2.0

dPath=data/${dname}/${dname}
rPath=results${c}/${dname}
dFolder=data/${dname}/

# ------------------------------------------------------------------------------
#  Ground Truth 
# ------------------------------------------------------------------------------
./rqalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
    -ts ${dPath}.fn${c}

# ------------------------------------------------------------------------------
#  RQALSH*
# ------------------------------------------------------------------------------
beta=100
delta=0.49
L_list=(1500 1000 750 600 500 300 200 150 100 75 60 50 40 30 20 15 10 6 5 4 3 2) 
M_list=(2 3 4 5 6 10 15 20 30 40 50 60 75 100 150 200 300 500 600 750 1000 1500)
length=`expr ${#L_list[*]} - 1`

for j in $(seq 0 ${length})
do 
    L=${L_list[j]}
    M=${M_list[j]}
    oFolder=${rPath}/rqalsh_star/${L}_${M}/

    ./rqalsh -alg 1 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -beta ${beta} \
        -delta ${delta} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

    ./rqalsh -alg 2 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
        -df ${dFolder} -of ${oFolder}
done

# ------------------------------------------------------------------------------
#  RQALSH
# ------------------------------------------------------------------------------
beta=100
delta=0.49
oFolder=${rPath}/rqalsh/

./rqalsh -alg 3 -n ${n} -d ${d} -B ${B} -beta ${beta} -delta ${delta} -c ${c} \
    -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./rqalsh -alg 4 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
    -df ${dFolder} -of ${oFolder}

# ------------------------------------------------------------------------------
#  Drusilla_Select
# ------------------------------------------------------------------------------
L_list=(16) 
M_list=(18)
length=`expr ${#L_list[*]} - 1`

for j in $(seq 0 ${length})
do 
    L=${L_list[j]}
    M=${M_list[j]}
    oFolder=${rPath}/drusilla_select/${L}_${M}/

    ./rqalsh -alg 5 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -ds ${dPath}.ds \
        -df ${dFolder} -of ${oFolder}

    ./rqalsh -alg 6 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
        -df ${dFolder} -of ${oFolder}
done

# ------------------------------------------------------------------------------
#  QDAFN (Guarantee mode)
# ------------------------------------------------------------------------------
L=0
M=0
oFolder=${rPath}/qdafn/guarantee/

./rqalsh -alg 7 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -c ${c} -ds ${dPath}.ds \
    -df ${dFolder} -of ${oFolder}

./rqalsh -alg 8 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
    -df ${dFolder} -of ${oFolder}

# ------------------------------------------------------------------------------
#  QDAFN (Heuristic mode)
# ------------------------------------------------------------------------------
proj=10
for ((i=2; i<=10; i=i+1))
do
    proj=$(($proj + 10))
    cand=$((273 - $proj))	
    for ((j=1; j<=5; j=j+1))
    do
        oFolder=${rPath}/qdafn/heuristic/${proj}_${j}/

        ./rqalsh -alg 7 -n ${n} -d ${d} -B ${B} -L ${proj} -M ${cand} -c ${c} \
            -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}
      
        ./rqalsh -alg 8 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
            -df ${dFolder} -of ${oFolder}
    done
done

# ------------------------------------------------------------------------------
#  Linear Scan
# ------------------------------------------------------------------------------
oFolder=${rPath}/

./rqalsh -alg 9 -n ${n} -qn ${qn} -d ${d} -B ${B} -qs ${dPath}.q \
    -ts ${dPath}.fn${c} -df ${dFolder} -of ${oFolder}
