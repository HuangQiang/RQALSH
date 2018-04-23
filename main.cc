#include "headers.h"

// -----------------------------------------------------------------------------
void usage() 						// display usage of the package
{
	printf("\n"
		"--------------------------------------------------------------------\n"
		" Usage of the Package for External c-k-AFN Search:                  \n"
		"--------------------------------------------------------------------\n"
		"    -alg   (integer)   options of algorithms\n"
		"    -n     (integer)   number of data  objects\n"
		"    -qn    (integer)   number of query objects\n"
		"    -d     (integer)   dimensionality\n"
		"    -B     (integer)   page size\n"
		"    -L	    (integer)   number of projection\n"
		"    -M     (integer)   number of candidates\n"
		"    -beta  (integer)   false positive percentage\n"
		"    -delta (real)      error probability\n"
		"    -c     (real)      approximation ratio (c > 1)\n"
		"    -ds    (string)    address of data  set\n"
		"    -qs    (string)    address of query set\n"
		"    -ts    (string)    address of truth set\n"
		"    -df    (string)    data folder to store new format of data\n"
		"    -of    (string)    output folder to store output results\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" The Options of Algorithms (-alg) are:                              \n"
		"--------------------------------------------------------------------\n"
		"    0 - Ground-Truth\n"
		"        Params: -alg 0 -n -qn -d -ds -qs -ts\n"
		"\n"
		"    1 - Indexing of RQALSH*\n"
		"        Params: -alg 1 -n -d -B -L -M -beta -delta -c -ds -df -of\n"
		"\n"
		"    2 - c-k-AFN Search of RQALSH*\n"
		"        Params: -alg 2 -qn -d -L -M -qs -ts -df -of\n"
		"\n"
		"    3 - Indexing of RQALSH\n"
		"        Params: -alg 3 -n -d -B -beta -delta -c -ds -df -of\n"
		"\n"
		"    4 - c-k-AFN Search of RQALSH\n"
		"        Params: -alg 4 -qn -d -qs -ts -df -of\n"
		"\n"
		"    5 - k-FN Search of Linear Scan\n"
		"        Params: -alg 5 -n -qn -d -B -qs -ts -df -of\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" Author: Qiang HUANG  (huangq2011@gmail.com)                        \n"
		"--------------------------------------------------------------------\n"
		"\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char** args)
{
	srand((unsigned)time(NULL));	// set the random seed
	// usage();

	int   alg   = -1;				// option of algorithm
	int   n     = -1;				// cardinality
	int   qn    = -1;				// query number
	int   d     = -1;				// dimensionality
	int   B     = -1;				// page size
	int   L     = -1;				// number of projection
	int   M     = -1;				// number of candidates
	int   beta  = -1;			    // false positive percentage
	float delta = -1.0f;			// error probability
	float ratio = -1.0f;			// approximation ratio
	char  data_set[200];			// address of data  set
	char  query_set[200];			// address of query set
	char  truth_set[200];			// address of truth set
	char  data_folder[200];			// data folder
	char  output_folder[200];		// output folder

	bool failed = false;
	int  cnt = 1;

	while (cnt < nargs && !failed) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg           = %d\n", alg);
			if (alg < 0 || alg > 5) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n             = %d\n", n);
			if (n <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn            = %d\n", qn);
			if (qn <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d             = %d\n", d);
			if (d <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-B") == 0) {
			B = atoi(args[++cnt]);
			printf("B             = %d\n", B);
			if (B <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-L") == 0) {
			L = atoi(args[++cnt]);
			printf("L             = %d\n", L);
			if (L <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-M") == 0) {
			M = atoi(args[++cnt]);
			printf("M             = %d\n", M);
			if (M <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-beta") == 0) {
			beta = atoi(args[++cnt]);
			printf("beta          = %d\n", beta);
			if (beta <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-delta") == 0) {
			delta = (float)atof(args[++cnt]);
			printf("delta         = %.2f\n", delta);
			if (delta < 0.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-c") == 0) {
			ratio = (float) atof(args[++cnt]);
			printf("c             = %.1f\n", ratio);
			if (ratio <= 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("data_set      = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query_set     = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth_set     = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-df") == 0) {
			strncpy(data_folder, args[++cnt], sizeof(data_folder));
			printf("data_folder   = %s\n", data_folder);

			int len = (int) strlen(data_folder);
			if (data_folder[len - 1] != '/') {
				data_folder[len] = '/';
				data_folder[len + 1] = '\0';
			}
			create_dir(data_folder);
		}
		else if (strcmp(args[cnt], "-of") == 0) {
			strncpy(output_folder, args[++cnt], sizeof(output_folder));
			printf("output_folder = %s\n", output_folder);

			int len = (int) strlen(output_folder);
			if (output_folder[len - 1] != '/') {
				output_folder[len] = '/';
				output_folder[len + 1] = '\0';
			}
			create_dir(output_folder);
		}
		else {
			failed = true;
			break;
		}
		cnt++;
	}
	printf("\n");

	switch (alg) {
	case 0:
		ground_truth(n, qn, d, data_set, query_set, truth_set);
		break;
	case 1:
		indexing_of_rqalsh_star(n, d, B, L, M, beta, delta, ratio,
			data_set, data_folder, output_folder);
		break;
	case 2:
		kfn_of_rqalsh_star(qn, d, L, M, query_set, truth_set, 
			data_folder, output_folder);
		break;
	case 3:
		indexing_of_rqalsh(n, d, B, beta, delta, ratio, data_set, 
			data_folder, output_folder);
		break;
	case 4:
		kfn_of_rqalsh(qn, d, query_set, truth_set, data_folder,
			output_folder);
		break;
	case 5:
		linear_scan(n, qn, d, B, query_set, truth_set, data_folder, 
			output_folder);
		break;
	default:
		printf("Parameters Error!\n");
		usage();
		break;
	}

	return 0;
}


