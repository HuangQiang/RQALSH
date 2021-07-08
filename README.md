## RQALSH: Reverse Query-Aware Locality-Sensitive Hashing

<p align="center"><img src="imgs/rqalsh.jpg" alt="drawing" width="450"/></p>

### Introduction

This package provides two external LSH schemes RQALSH and RQALSH<sup>*</sup> for high-dimensional ```c-Approximate Furthest Neighbor (c-AFN)``` search from the following two papers:

```bash
Qiang Huang, Jianlin Feng, Qiong Fang. Reverse Query-Aware Locality-Sensitive
Hashing for High-Dimensional Furthest Neighbor Search. 2017 IEEE 33rd International 
Conference on Data Engineering (ICDE), pages 167-170, 2017.

Qiang Huang, Jianlin Feng, Qiong Fang, Wilfred Ng. Two Efficient Hashing Schemes for 
High-Dimensional Furthest Neighbor Search. IEEE Transactions on Knowledge and Data 
Engineering (TKDE), 29(12): 2772–2785, 2017.
```

### Compilation

The package requires ```g++``` with ```c++11``` support. To download and compile the code, type:

```bash
$ git clone https://github.com/HuangQiang/RQALSH.git
$ cd RQALSH
$ make
```

### Datasets

We use four real-life datasets [Sift](https://drive.google.com/open?id=1tgcUU9X61TehVa_Klj5skVdYRoYZ7CgX), [Gist](https://drive.google.com/open?id=1fvUTGUbYgg8oaGNbZbAMLnfmxoU8UDhh), [Trevi](https://drive.google.com/open?id=1XSiiQ6D1zoxGXULl3sHxsjPO8JCM-md1), and [P53](https://drive.google.com/open?id=1hjGvcq29WsgHpGoz0vCdCYAUR453aY29) for comparison. We randomly remove 1,000 data objects from each dataset and use them as queries. The statistics of datasets and queries are summarized in the following table:

| Datasets | #Objects  | #Queries | Dimensionality | Page Size | Domain Size | Data Size |
| -------- | --------- | -------- | -------------- | --------- | ----------- | --------- |
| Sift     | 1,000,000 | 1000     | 128            | 4 KB      | [0, 218]    | 337.8 MB  |
| Gist     | 1,000,000 | 1000     | 960            | 16 KB     | [0, 14,772] | 4.0 GB    |
| Trevi    | 100,900   | 1000     | 4,096          | 64 KB     | [0, 255]    | 1.5 GB    |
| P53      | 31,159    | 1000     | 5,408          | 64 KB     | [0, 10,000] | 833.7 MB  |

### Run Experiments

```bash
Usage: qalsh [OPTIONS]

This package supports 10 options to evaluate the performance of RQALSH, RQALSH*,
QDAFN, QDAFN*, Drusilla_Select, and Linear_Scan for c-AFN search. The parameters
are introduced as follows.

  -alg    integer    options of algorithms (0 - 9)
  -n      integer    cardinality of dataset
  -d      integer    dimensionality of dataset and query set
  -qn     integer    number of queries
  -B      integer    page size
  -L      integer    number of projections for RQALSH*, QDAFN*, Drusilla_Select
  -M      integer    number of candidates  for RQALSH*, QDAFN*, Drusilla_Select
  -beta   integer    number of false positives
  -delta  float      error probability
  -c      float      approximation ratio for c-AFN search (c > 1)
  -ds     string     address of data  set
  -qs     string     address of query set
  -ts     string     address of truth set
  -df     string     data folder to store new format of data
  -of     string     output folder to store output results
```

We provide the scripts to repeat experiments reported in TKDE 2017. A quick example is shown as follows (run RQALSH<sup>*</sup> and RQALSH on ```Mnist```):

```bash
# RQALSH*
./rqalsh -alg 1 -n 59000 -d 50 -B 4096 -L 1000 -M 3 -beta 100 -delta 0.49 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results/rqalsh_star/
./rqalsh -alg 2 -qn 1000 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.fn2.0 -df data/Mnist/ -of results/rqalsh_star/

# RQALSH
./rqalsh -alg 3 -n 59000 -d 50 -B 4096 -beta 100 -delta 0.49 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results/rqalsh/
./rqalsh -alg 4 -qn 1000 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.fn2.0 -df data/Mnist/ -of results/rqalsh/
```

If you would like to get more information to run other algorithms, please check the scripts in the package. When you run the package, please ensure that the path for the dataset, query set, and truth set is correct. Since the package will automatically create folder for the output path, please keep the path as short as possible.

### Related Publications

If you use this package for publications, please cite the papers as follows.

```bib
@inproceedings{huang2017reverse,
    title={Reverse query-aware locality-sensitive hashing for high-dimensional furthest neighbor search}
    author={Huang, Qiang and Feng, Jianlin and Fang, Qiong},
    booktitle={2017 IEEE 33rd International Conference on Data Engineering (ICDE)},
    pages={167--170},
    year={2017},
    organization={IEEE}
}

@article{huang2017two,
    title={Two efficient hashing schemes for high-dimensional furthest neighbor search}
    author={Huang, Qiang and Feng, Jianlin and Fang, Qiong and Ng, Wilfred},
    booktitle={IEEE Transactions on Knowledge and Data Engineering},
    volumn={29},
    number={12},
    pages={2772--2785},
    year={2017},
    organization={IEEE}
}
```
