#!/bin/bash
make
rm *.o

# Parameters of QALSH:
#    -alg  (integer)   options of algorithms (0 - 3)
#    -d    (integer)   dimensionality of the dataset
#    -n    (integer)   cardinality of the dataset
#    -qn   (integer)   number of queries
#    -B    (integer)   page size
#    -c    (real)      approximation ratio (c > 1)
#    -ds   (string)    file path of the dataset
#    -qs   (string)    file path of the query set
#    -ts   (string)    file path of the ground truth set
#    -df   (string)    data folder to store new format of data
#    -of   (string)    output folder to store info of rqalsh
#
# The options of algorithms (-alg) are:
#    0 - Ground-Truth
#        Parameters: -alg 0 -n -qn -d -ds -qs -ts -of
#
#    1 - Indexing
#        Parameters: -alg 1 -n -d -B -c -ds -df -of
#
#    2 - RQALSH
#        Parameters: -alg 2 -qn -d -qs -ts -df -of
#
#    3 - Linear Scan
#        Parameters: -alg 3 -n -qn -d -B -qs -ts -df -of
#
# NOTE: Each parameter is required to be separated by one space


dataPath=./data
outPath=.


for ratio in 1.5 2.5 3.0
do	
	###### Sift ######
	./rqalsh -alg 1 -n 999000 -d 128 -B 4096 -beta 100 -delta 0.49 -c ${ratio} -ds ${dataPath}/Sift/Sift.ds -df ${dataPath}/Sift/ -of ${outPath}/result${ratio}/rqalsh/Sift/
	
	./rqalsh -alg 2 -qn 1000 -d 128 -qs ${dataPath}/Sift/Sift.q -ts ${dataPath}/Sift/SiftID.fn2.0 -df ${dataPath}/Sift/ -of ${outPath}/result${ratio}/rqalsh/Sift/
	
	
	###### Gist ######
	./rqalsh -alg 1 -n 999000 -d 960 -B 16384 -beta 100 -delta 0.49 -c ${ratio} -ds ${dataPath}/Gist/Gist.ds -df ${dataPath}/Gist/ -of ${outPath}/result${ratio}/rqalsh/Gist/
	
	./rqalsh -alg 2 -qn 1000 -d 960 -qs ${dataPath}/Gist/Gist.q -ts ${dataPath}/Gist/GistID.fn2.0 -df ${dataPath}/Gist/ -of ${outPath}/result${ratio}/rqalsh/Gist/


	###### Trevi ######
	./rqalsh -alg 1 -n 99900 -d 4096 -B 65536 -beta 100 -delta 0.49 -c ${ratio} -ds ${dataPath}/Trevi/Trevi.ds -df ${dataPath}/Trevi/ -of ${outPath}/result${ratio}/rqalsh/Trevi/
	
	./rqalsh -alg 2 -qn 1000 -d 4096 -qs ${dataPath}/Trevi/Trevi.q -ts ${dataPath}/Trevi/TreviID.fn2.0 -df ${dataPath}/Trevi/ -of ${outPath}/result${ratio}/rqalsh/Trevi/


	###### P53 ######
	./rqalsh -alg 1 -n 30159 -d 5408 -B 65536 -beta 100 -delta 0.49 -c ${ratio} -ds ${dataPath}/P53/P53.ds -df ${dataPath}/P53/ -of ${outPath}/result${ratio}/rqalsh/P53/
	
	./rqalsh -alg 2 -qn 1000 -d 5408 -qs ${dataPath}/P53/P53.q -ts ${dataPath}/P53/P53ID.fn2.0 -df ${dataPath}/P53/ -of ${outPath}/result${ratio}/rqalsh/P53/
done



for beta in 10 50 200 500
do	
	###### Sift ######
	./rqalsh -alg 1 -n 999000 -d 128 -B 4096 -beta ${beta} -delta 0.49 -c 2.0 -ds ${dataPath}/Sift/Sift.ds -df ${dataPath}/Sift/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Sift/
	
	./rqalsh -alg 2 -qn 1000 -d 128 -qs ${dataPath}/Sift/Sift.q -ts ${dataPath}/Sift/SiftID.fn2.0 -df ${dataPath}/Sift/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Sift/
	
	
	###### Gist ######
	./rqalsh -alg 1 -n 999000 -d 960 -B 16384 -beta ${beta} -delta 0.49 -c 2.0 -ds ${dataPath}/Gist/Gist.ds -df ${dataPath}/Gist/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Gist/
	
	./rqalsh -alg 2 -qn 1000 -d 960 -qs ${dataPath}/Gist/Gist.q -ts ${dataPath}/Gist/GistID.fn2.0 -df ${dataPath}/Gist/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Gist/


	###### Trevi ######
	./rqalsh -alg 1 -n 99900 -d 4096 -B 65536 -beta ${beta} -delta 0.49 -c 2.0 -ds ${dataPath}/Trevi/Trevi.ds -df ${dataPath}/Trevi/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Trevi/
	
	./rqalsh -alg 2 -qn 1000 -d 4096 -qs ${dataPath}/Trevi/Trevi.q -ts ${dataPath}/Trevi/TreviID.fn2.0 -df ${dataPath}/Trevi/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/Trevi/


	###### P53 ######
	./rqalsh -alg 1 -n 30159 -d 5408 -B 65536 -beta ${beta} -delta 0.49 -c 2.0 -ds ${dataPath}/P53/P53.ds -df ${dataPath}/P53/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/P53/
	
	./rqalsh -alg 2 -qn 1000 -d 5408 -qs ${dataPath}/P53/P53.q -ts ${dataPath}/P53/P53ID.fn2.0 -df ${dataPath}/P53/ -of ${outPath}/result2.0/rqalsh_beta=${beta}/P53/
done



for delta in 0.1 0.2 0.3 0.4
do	
	###### Sift ######
	./rqalsh -alg 1 -n 999000 -d 128 -B 4096 -beta 100 -delta ${delta} -c 2.0 -ds ${dataPath}/Sift/Sift.ds -df ${dataPath}/Sift/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Sift/
	
	./rqalsh -alg 2 -qn 1000 -d 128 -qs ${dataPath}/Sift/Sift.q -ts ${dataPath}/Sift/SiftID.fn2.0 -df ${dataPath}/Sift/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Sift/
	
	
	###### Gist ######
	./rqalsh -alg 1 -n 999000 -d 960 -B 16384 -beta 100 -delta ${delta} -c 2.0 -ds ${dataPath}/Gist/Gist.ds -df ${dataPath}/Gist/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Gist/
	
	./rqalsh -alg 2 -qn 1000 -d 960 -qs ${dataPath}/Gist/Gist.q -ts ${dataPath}/Gist/GistID.fn2.0 -df ${dataPath}/Gist/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Gist/


	###### Trevi ######
	./rqalsh -alg 1 -n 99900 -d 4096 -B 65536 -beta 100 -delta ${delta} -c 2.0 -ds ${dataPath}/Trevi/Trevi.ds -df ${dataPath}/Trevi/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Trevi/
	
	./rqalsh -alg 2 -qn 1000 -d 4096 -qs ${dataPath}/Trevi/Trevi.q -ts ${dataPath}/Trevi/TreviID.fn2.0 -df ${dataPath}/Trevi/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/Trevi/


	###### P53 ######
	./rqalsh -alg 1 -n 30159 -d 5408 -B 65536 -beta 100 -delta ${delta} -c 2.0 -ds ${dataPath}/P53/P53.ds -df ${dataPath}/P53/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/P53/
	
	./rqalsh -alg 2 -qn 1000 -d 5408 -qs ${dataPath}/P53/P53.q -ts ${dataPath}/P53/P53ID.fn2.0 -df ${dataPath}/P53/ -of ${outPath}/result2.0/rqalsh_delta=${delta}/P53/
done
