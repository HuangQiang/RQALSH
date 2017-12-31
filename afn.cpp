#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// output the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	char* data_set,						// address of data set
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* output_folder)				// folder to store info of rqalsh
{
	clock_t start_time = (clock_t) -1;
	clock_t end_time   = (clock_t) -1;

	int i, j;
	FILE* fp = NULL;
	
	int maxk = MAXK;
	MaxK_List* list = new MaxK_List(MAXK);
	g_memory += (SIZEFLOAT + SIZEINT) * maxk;

	// -------------------------------------------------------------------------
	//  Read data set and query set
	// -------------------------------------------------------------------------
	start_time = clock();
	g_memory += SIZEFLOAT * (n + qn) * d;
	if (check_mem()) return 1;

	float** data = new float*[n];
	for (i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data)) {
		error("Reading Dataset Error!\n", true);
	}

	float** query = new float*[qn];
	for (i = 0; i < qn; i++) query[i] = new float[d];
	if (read_set(qn, d, query_set, query) == 1) {
		error("Reading Query Set Error!\n", true);
	}
	end_time = clock();
	printf("Read Dataset and Query Set: %.6f Seconds\n\n", 
		((float) end_time - start_time) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  output ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	start_time = clock();
	fp = fopen(truth_set, "w");		// open output file
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	fprintf(fp, "%d %d\n", qn, maxk);
	for (i = 0; i < qn; i++) {
		list->reset();
		for (j = 0; j < n; j++) {	// find k-fn points of query
			float dist = calc_l2_dist(data[j], query[i], d);
			int id = j + 1;
			list->insert(dist, id);
		}

		fprintf(fp, "%d", i + 1);	// output Lp dist of k-fn points
		for (j = 0; j < maxk; j++) {
			fprintf(fp, " %f", list->ith_largest_key(j));
		}
		for (j = 0; j < maxk; j++) {
			fprintf(fp, " %d", list->ith_largest_id(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);						// close output file
	end_time = clock();
	printf("Generate Ground Truth: %.6f Seconds\n\n", 
		((float) end_time - start_time) / CLOCKS_PER_SEC);
		
	// -------------------------------------------------------------------------
	//  Calculate relative contrast for each query
	// -------------------------------------------------------------------------
	start_time = clock();
	char fname[200];
	strcpy(fname, output_folder);	// generate output file
	strcat(fname, "rc.out");
	
	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;
	}
	calc_relative_constrast(n, qn, d, fp, data, query);
	fclose(fp);						// close output file
	end_time = clock();
	printf("Calculate Relative Contrast: %.6f Seconds\n\n", 
		((float) end_time - start_time) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (data != NULL) {				// release <data>
		for (i = 0; i < n; i++) {
			delete[] data[i]; data[i] = NULL;
		}
		delete[] data; data = NULL;
		g_memory -= SIZEFLOAT * n * d;
	}
	if (query != NULL) {			// release <query>
		for (i = 0; i < qn; i++) {
			delete[] query[i]; query[i] = NULL;
		}
		delete[] query; query = NULL;
		g_memory -= SIZEFLOAT * qn * d;
	}
	if (list != NULL) {				// release <list>
		delete list; list = NULL;
		g_memory -= (SIZEFLOAT + SIZEINT) * maxk;
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int indexing(						// build hash tables for the dataset
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   beta,
	float delta,
	float ratio,						// approximation ratio
	char* data_set,						// address of data set
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// folder to store info of rqalsh
{
	clock_t start_time = (clock_t) -1;
	clock_t end_time   = (clock_t) -1;

	// -------------------------------------------------------------------------
	//  Read data set
	// -------------------------------------------------------------------------
	start_time = clock();
	g_memory += SIZEFLOAT * n * d;
	if (check_mem()) return 1;

	float** data = new float*[n];
	for (int i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data) == 1) {
		error("Reading Dataset Error!\n", true);
	}
	end_time = clock();
	printf("Read Dataset: %.6f Seconds\n\n", 
		((float) end_time - start_time) / CLOCKS_PER_SEC);


	// -------------------------------------------------------------------------
	//  Write the data set in new format to disk
	// -------------------------------------------------------------------------
	start_time = clock();
	write_data_new_form(n, d, B, data, data_folder);
	end_time = clock();
	printf("Write Dataset in New Format: %.6f Seconds\n\n", 
		((float) end_time - start_time) / CLOCKS_PER_SEC);
	
	// -------------------------------------------------------------------------
	//  Indexing
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, output_folder);
	strcat(fname, "index.out");

	FILE* fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;					// fail to return
	}

	start_time = clock();
	RQALSH* lsh = new RQALSH();
	lsh->init(n, d, B, beta, delta, ratio, output_folder);
	lsh->bulkload(data);
	end_time = clock();

	float indexing_time = ((float) end_time - start_time) / CLOCKS_PER_SEC;
	printf("\nIndexing Time: %.6f seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time: %.6f seconds\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (data != NULL) {
		for (int i = 0; i < n; i++) {
			delete[] data[i]; data[i] = NULL;
		}
		delete[] data; data = NULL;
		g_memory -= SIZEFLOAT * n * d;
	}
	if (lsh != NULL) {
		delete lsh; lsh = NULL;
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh_afn(						// c-k-AFN search
	int   qn,							// number of query points
	int   d,							// dimensionality
	char* query_set,					// path of query set
	char* truth_set,					// groundtrue file
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// output folder
{
	clock_t start_time = (clock_t) -1;
	clock_t end_time   = (clock_t) -1;

	int maxk = MAXK;
	int i, j;
	FILE* fp = NULL;				// file pointer

	// -------------------------------------------------------------------------
	//  Read query set
	// -------------------------------------------------------------------------
	g_memory += (SIZEFLOAT * qn * d);
	float** query = new float*[qn];
	for (i = 0; i < qn; i++) {
		query[i] = new float[d];
	}

	if (read_set(qn, d, query_set, query)) {
		error("Reading Query Set Error!\n", true);
	}

	// -------------------------------------------------------------------------
	//  Read the ground truth file
	// -------------------------------------------------------------------------
	g_memory += ((SIZEFLOAT + SIZEINT) * qn * maxk);
	float** R = new float*[qn];
	int**   ID = new int*[qn];
	for (i = 0; i < qn; i++) {
		R[i]  = new float[maxk];
		ID[i] = new int[maxk];
		for (j = 0; j < maxk; j++) {
			R[i][j]  = -1.0;
			ID[i][j] = -1;
		}
	}

	fp = fopen(truth_set, "r");		// open ground truth file
	if (!fp) {
		printf("Could not open %s.\n", truth_set);
		return 1;
	}

	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (int i = 0; i < qn; i++) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < maxk; j++) {
			fscanf(fp, "%f", &R[i][j]);
		}
		for (j = 0; j < maxk; j++) {
			fscanf(fp, "%d", &ID[i][j]);
		}
	}
	fclose(fp);						// close groundtrue file

	// -------------------------------------------------------------------------
	//  c-k-Approximate Furthest Neighbor (c-k-AFN) search via rqalsh
	// -------------------------------------------------------------------------
	int kFNs[] = {1, 2, 5, 10};
	int max_round = 4;
	int top_k = 0;

	float all_time    = -1.0f;		// average running time
	float all_ratio   = -1.0f;		// average ratio
	float this_ratio  = -1.0f;

	float all_recall  = -1.0f;		// average recall
	float this_recall = -1.0f;
	int   this_io = 0;				// average page IOs
	int   all_io  = 0;

	RQALSH* lsh = new RQALSH();		// restore RQALSH
	if (lsh->restore(output_folder)) {
		error("Could not restore rqalsh\n", true);
	}

	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "rqalsh.out");

	fp = fopen(output_set, "w");	// open output file
	if (!fp) {
		printf("Could not create the output file.\n");
		return 1;
	}

	printf("RQALSH for c-k-AFN Search: \n");
	printf("    Top-k\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		start_time = clock();
		top_k = kFNs[num];
		MaxK_List* list = new MaxK_List(top_k);

		all_ratio  = 0.0f;
		all_recall = 0.0f;
		all_io     = 0;
		for (i = 0; i < qn; i++) {
			list->reset();
			this_io = lsh->kfn(top_k, query[i], data_folder, list);
			all_io += this_io;

			this_ratio = 0.0f;
			for (j = 0; j < top_k; j++) {
				this_ratio += (R[i][j] / list->ith_largest_key(j));
			}
			this_ratio /= top_k;
			all_ratio += this_ratio;

			this_recall = calc_recall(top_k, list, ID[i]);
			all_recall += this_recall;
			/*
			printf("        No = %2d, Top-k = %d, IO = %4d, Ratio = %0.4f, "
				"Recall = %.2f%%\n", i + 1, top_k, this_io, this_ratio, 
				this_recall * 100.0f);
			*/
		}
		delete list; list = NULL;
		end_time = clock();
		all_time = ((float) end_time - start_time) / CLOCKS_PER_SEC;

		all_ratio  = all_ratio / qn;
		all_time   = (all_time * 1000.0f) / qn;
		all_recall = (all_recall * 100.0f) / qn;
		all_io     = (int) ceil((float) all_io / (float) qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", 
			top_k, all_ratio, all_io, all_time, all_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", 
			top_k, all_ratio, all_io, all_time, all_recall);
	}
	printf("\n");
	fclose(fp);						// close output file

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (query != NULL) {			// release <query>
		for (i = 0; i < qn; i++) {
			delete[] query[i]; query[i] = NULL;
		}
		delete[] query; query = NULL;
		g_memory -= SIZEFLOAT * qn * d;
	}
	if (lsh != NULL) {				// release <lsh>
		delete lsh; lsh = NULL;
	}
	if (R != NULL || ID != NULL) {	// release <R> and <ID>
		for (i = 0; i < qn; i++) {
			delete[] R[i]; R[i] = NULL;
			delete[] ID[i]; ID[i] = NULL;
		}
		delete[] R; R = NULL;
		delete[] ID; ID = NULL;
		g_memory -= ((SIZEFLOAT + SIZEINT) * qn * maxk);
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   B,							// page size
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// output folder
{
	clock_t start_time = (clock_t) -1;
	clock_t end_time   = (clock_t) -1;

	int i, j;
	int maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  Open the output file, and read the ground true results
	// -------------------------------------------------------------------------
	g_memory += ((SIZEFLOAT + SIZEINT) * qn * maxk);
	float** R = new float*[qn];
	int**  ID = new int*[qn];
	for (i = 0; i < qn; i++) {
		R[i]  = new float[maxk];
		ID[i] = new int[maxk];
		for (j = 0; j < maxk; j++) {
			R[i][j]  = -1.0F;
			ID[i][j] = -1;
		}
	}

	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}
									// read top-k furthest distance
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (int i = 0; i < qn; i++) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < maxk; j++) {
			fscanf(fp, "%f", &R[i][j]);
		}
		for (j = 0; j < maxk; j++) {
			fscanf(fp, "%d", &ID[i][j]);
		}
	}
	fclose(fp);						// close ground true file

	// -------------------------------------------------------------------------
	//  Calc the number of data object in one page and the number of data file.
	//  <num> is the number of data in one data file
	//  <total_file> is the total number of data file
	// -------------------------------------------------------------------------
	int num = (int) floor((float) B / (d * SIZEFLOAT));
	int total_file = (int) ceil((float) n / num);
	if (total_file == 0) return 1;

	// -------------------------------------------------------------------------
	//  Brute-Force Linear Scan method for Furthest Neighbor Search
	//  For each query, we limit that we can ONLY read one page of data.
	// -------------------------------------------------------------------------
	g_memory += (SIZEFLOAT * d * 2 + SIZECHAR * (600 + B));
	float* data     = new float[d];	// one data object
	float* query    = new float[d];	// one query object

	char* buffer    = new char[B];	// every time one can read one page
	char* fname     = new char[200];// file name for data
	char* data_path = new char[200];// data path
	char* out_set   = new char[200];// output file

	strcpy(out_set, output_folder);	// generate output file
	strcat(out_set, "linear.out");

	fp = fopen(out_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", out_set);
		return 1;
	}

	int kFNs[] = {1, 2, 5, 10};
	int max_round = 4;

	float all_time    = -1.0f;
	float this_ratio  = -1.0f;
	float all_ratio   = -1.0f;
	float this_recall = -1.0f;
	float all_recall  = -1.0f;

	int   top_k = -1;
	int   count = -1;
	float dist  = -1.0f;
									// generate the data path
	strcpy(data_path, data_folder);
	strcat(data_path, "data/");

	printf("Linear Scan Search:\n");
	printf("    Top-k\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < max_round; round++) {
		top_k = kFNs[round];
		all_ratio = 0.0f;
		all_recall = 0.0f;

		start_time = clock();
		MaxK_List* list = new MaxK_List(top_k);
		FILE* qfp = fopen(query_set, "r");
		if (!qfp) error("Could not open the query set.\n", true);

		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  Step 1: read a query from disk
			// -----------------------------------------------------------------
			fscanf(qfp, "%d", &j);
			for (j = 0; j < d; j++) {
				fscanf(qfp, " %f", &query[j]);
			}

			// -----------------------------------------------------------------
			//  Step 2: find k-fn results for the query
			// -----------------------------------------------------------------
			list->reset();
			int id = 0;
			for (j = 0; j < total_file; j++) {
				// -------------------------------------------------------------
				//  Step 2.1: get the file name of current data page
				// -------------------------------------------------------------
				get_data_filename(j, data_path, fname);

				// -------------------------------------------------------------
				//  Step 2.2: read one page of data into buffer
				// -------------------------------------------------------------
				if (read_buffer_from_page(B, fname, buffer) == 1) {
					error("error to read a data page", true);
				}

				// -------------------------------------------------------------
				//  Step 2.3: find the k-fn results in this page. NOTE: the 
				//  number of data in the last page may be less than <num>
				// -------------------------------------------------------------
				if (j < total_file - 1) count = num;
				else count = n - num * (total_file - 1);

				for (int z = 0; z < count; z++) {
					read_data_from_buffer(z, d, data, buffer);
					dist = calc_l2_dist(data, query, d);
					id = id + 1;
					list->insert(dist, id);
				}
			}

			this_ratio = 0.0f;
			for (j = 0; j < top_k; j++) {
				this_ratio += (R[i][j] / list->ith_largest_key(j));
			}
			this_ratio /= top_k;
			all_ratio += this_ratio;

			this_recall = calc_recall(top_k, list, ID[i]);
			all_recall += this_recall;
		}
		// ---------------------------------------------------------------------
		//  Step 3: output the results of top-k furthest neighbor points
		// ---------------------------------------------------------------------
		delete list; list = NULL;
		fclose(qfp);
		end_time   = clock();
		all_time   = ((float) end_time - start_time) / CLOCKS_PER_SEC;
		all_time   = (all_time * 1000.0f) / qn;
		all_ratio  = all_ratio / qn;
		all_recall = (all_recall * 100.0f) / qn;
									// output results
		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", 
			top_k, all_ratio, total_file, all_time, all_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", 
			top_k, all_ratio, total_file, all_time, all_recall);
	}
	printf("\n");
	fclose(fp);						// close output file

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (R != NULL || ID != NULL) {
		for (i = 0; i < qn; i++) {
			delete[] R[i];  R[i]  = NULL;
			delete[] ID[i]; ID[i] = NULL;
		}
		delete[] R;  R  = NULL;
		delete[] ID; ID = NULL;
		g_memory -= ((SIZEFLOAT + SIZEINT) * qn * maxk);
	}
	if (buffer != NULL || data != NULL || query != NULL) {
		delete[] buffer; buffer = NULL;
		delete[] data; data = NULL;
		delete[] query; query = NULL;
		g_memory -= (SIZEFLOAT * d * 2 + SIZECHAR * B);
	}
	if (fname != NULL || data_path != NULL || out_set != NULL) {
		delete[] fname; fname = NULL;
		delete[] data_path; data_path = NULL;
		delete[] out_set; out_set = NULL;
		g_memory -= (SIZECHAR * 600);
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}
