# RQALSH: Reverse Query-Aware Locality-Sensitive Hashing Scheme

Version: 1.0.1

Release date: 01-01-2018

Last Modified: 23-04-2018


Introduction
--------

This package is written in the C++ programming language. It provides two 
external LSH schemes RQALSH and RQALSH* for c-Approximate Furthest Neighbor 
(or simply c-AFN) search under Euclidean distance.

How to use this Package?
--------

We provide a Makefile and a script (i.e., run_mnist.sh) as a running example 
for comipling and running this package. Before start running this package, 
please ensure the input format of the dataset and query set is correct. We 
provide a sample dataset and query set (i.e., Mnist) for your reference.

We also provide the scripts (i.e., run_sift.sh, run_gist.sh, run_trevi.sh,
run_p53.sh, and run_para.sh) for the users who would like to reproduce our 
results presented in ICDE 2017 and TKDE 2017. The datasets Sift, Gist, Trevi, 
and P53 we used can be downloaded from the following links:

* Sift: https://drive.google.com/open?id=1tgcUU9X61TehVa_Klj5skVdYRoYZ7CgX

* Gist: https://drive.google.com/open?id=1fvUTGUbYgg8oaGNbZbAMLnfmxoU8UDhh

* Trevi: https://drive.google.com/open?id=1XSiiQ6D1zoxGXULl3sHxsjPO8JCM-md1

* P53: https://drive.google.com/open?id=1hjGvcq29WsgHpGoz0vCdCYAUR453aY29


Author
--------

* **Qiang Huang**

  Smart Systems Institute, National University of Singapore (NUS),
  
  Singapore, 119613 
  
  huangq2011@gmail.com, huangq25@mail2.sysu.edu.cn
  
  https://sites.google.com/site/qianghuang2017/


Relevant Papers
--------

The paper for the package of RQALSH has been published in ICDE 2017 and TKDE 2017, 
which are displayed as follows:

* **Qiang Huang, Jianlin Feng, Qiong Fang. Reverse Query-Aware Locality-Sensitive 
Hashing for High-Dimensional Furthest Neighbor Search. IEEE International Conference 
on Data Engineering (ICDE), 167 - 170, 2017.**

* **Qiang Huang, Jianlin Feng, Qiong Fang, Wilfred Ng. Two Efficient Hashing Schemes 
for High-Dimensional Furthest Neighbor Search. IEEE Transactions on Knowledge and 
Data Engineering (TKDE) 29(12), 2772 - 2785, 2017.**

If you use the package for publications, please cite the papers above.
