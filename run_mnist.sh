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
dFolder=data/${dname}/
oFolder=results${c}/${dname}/

# # ------------------------------------------------------------------------------
# #  Ground Truth 
# # ------------------------------------------------------------------------------
# ./rqalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
#   -ts ${dPath}.fn${c}

# # ------------------------------------------------------------------------------
# #  RQALSH_Star
# # ------------------------------------------------------------------------------
# beta=100
# delta=0.49

# L_list=(1500 1000 750 600 500 300 150 100) 
# M_list=(2 3 4 5 6 10 20 30)
# length=`expr ${#L_list[*]} - 1`

# for j in $(seq 0 ${length})
# do 
#   L=${L_list[j]}
#   M=${M_list[j]}

#   ./rqalsh -alg 1 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -beta ${beta} \
#     -delta ${delta} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

#   ./rqalsh -alg 2 -qn ${qn} -d ${d} -L ${L} -M ${M} -qs ${dPath}.q \
#     -ts ${dPath}.fn${c} -df ${dFolder} -of ${oFolder}
# done

# # ------------------------------------------------------------------------------
# #  RQALSH
# # ------------------------------------------------------------------------------
# beta=100
# delta=0.49

# ./rqalsh -alg 3 -n ${n} -d ${d} -B ${B} -beta ${beta} -delta ${delta} -c ${c} \
#   -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

# ./rqalsh -alg 4 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
#   -df ${dFolder} -of ${oFolder}

# # ------------------------------------------------------------------------------
# #  Drusilla_Select
# # ------------------------------------------------------------------------------
# L_list=(16) 
# M_list=(18)
# length=`expr ${#L_list[*]} - 1`

# for j in $(seq 0 ${length})
# do 
#   L=${L_list[j]}
#   M=${M_list[j]}
#   outFolder=${oFolder}drusilla/${L}-${M}/

#   ./rqalsh -alg 5 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -ds ${dPath}.ds \
#     -df ${dFolder} -of ${outFolder}

#   ./rqalsh -alg 6 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
#     -df ${dFolder} -of ${outFolder}
# done

# ------------------------------------------------------------------------------
#  QDAFN (Guarantee mode)
# ------------------------------------------------------------------------------
L=0
M=0
outFolder=${oFolder}qdafn/guarantee/

./rqalsh -alg 7 -n ${n} -d ${d} -B ${B} -L ${L} -M ${M} -c ${c} -ds ${dPath}.ds \
  -df ${dFolder} -of ${outFolder}

./rqalsh -alg 8 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
  -df ${dFolder} -of ${outFolder}

# # ------------------------------------------------------------------------------
# #  QDAFN (Heuristic mode)
# # ------------------------------------------------------------------------------
# proj=10
# for ((i=2; i<=10; i=i+1))
# do
#   proj=$(($proj + 10))
#   cand=$((273 - $proj))	
#   for ((j=1; j<=5; j=j+1))
#   do
# 	  outFolder=${oFolder}qdafn/heuristic/${proj}/${j}/

# 	  ./rqalsh -alg 7 -n ${n} -d ${d} -B ${B} -L ${proj} -M ${cand} -c ${c} \
#       -ds ${dPath}.ds -df ${dFolder} -of ${outFolder}
	
#     ./rqalsh -alg 8 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.fn${c} \
#       -df ${dFolder} -of ${outFolder}
#   done
# done

# # ------------------------------------------------------------------------------
# #  Linear Scan
# # ------------------------------------------------------------------------------
# ./rqalsh -alg 5 -n ${n} -qn ${qn} -d ${d} -B ${B} -qs ${dPath}.q \
#   -ts ${dPath}.fn${c} -df ${dFolder} -of ${oFolder}
