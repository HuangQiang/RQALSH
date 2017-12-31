// -----------------------------------------------------------------------------
// Copyright (c) 2013-2018 Sun Yat-Sen University (SYSU). All Rights Reserved.
//
// SYSU grants permission to use, copy, modify, and distribute this software
// and its documentation for NON-COMMERCIAL purposes and without fee, provided 
// that this copyright notice appears in all copies.
//
// SYSU provides this software "as is," without representations or warranties
// of any kind, either expressed or implied, including but not limited to the
// implied warranties of merchantability, fitness for a particular purpose, 
// and noninfringement. SYSU shall not be liable for any damages arising from
// any use of this software.
//
// Authors: Qiang Huang  (huangq2011@gmail.com, huangq25@mail2.sysu.edu.cn)
//
// Created on:       2016-01-12
// Last Modified on: 2017-09-01
// Version 1.0.0
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <errno.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

// -----------------------------------------------------------------------------
//  For Linux directory
// -----------------------------------------------------------------------------
#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#define  LINUX_

// -----------------------------------------------------------------------------
//  For Windows directory
// -----------------------------------------------------------------------------
//#include <direct.h>
//#include <io.h>


#include "def.h"
#include "pri_queue.h"
#include "util.h"
#include "random.h"
#include "block_file.h"
#include "b_node.h"
#include "b_tree.h"
#include "rqalsh.h"
#include "afn.h"

using namespace std;

