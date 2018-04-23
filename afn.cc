#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set and query set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data and Query: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MaxK_List *list = new MaxK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {
			float dist = calc_l2_dist(d, data[j], query[i]);
			list->insert(dist, j + 1);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float truth_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_rqalsh_star(		// indexing of RQALSH_Star
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read dataset
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  write the data set in new format to disk
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	write_data_new_form(n, d, B, (const float **) data, data_folder);

	gettimeofday(&end_time, NULL);
	float write_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Write Dataset in New Format: %f Seconds\n\n", write_file_time);

	// -------------------------------------------------------------------------
	//  indexing of RQALSH_Star
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	char index_path[200];
	sprintf(index_path, "%srqalsh_star_L=%d_M=%d/", output_folder, L, M);

	RQALSH_Star *lsh = new RQALSH_Star();
	lsh->build(n, d, B, L, M, beta, delta, ratio, index_path, (const float**) data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  write indexing time to disk
	// -------------------------------------------------------------------------
	char fname[200];
	sprintf(fname, "%srqalsh_star_L=%d_M=%d.index", output_folder, L, M);

	FILE *fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;
	}
	fprintf(fp, "Indexing Time: %f seconds\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_rqalsh_star(				// c-k-AFN search of RQALSH_Star
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection
	int   M,							// number of candidates
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  load RQALSH_Star
	// -------------------------------------------------------------------------
	char index_path[200];
	sprintf(index_path, "%srqalsh_star_L=%d_M=%d/", output_folder, L, M);

	RQALSH_Star *lsh = new RQALSH_Star();
	if (lsh->load(index_path)) {
		printf("Could not load RQALSH_Star\n");
		exit(1);
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search by RQALSH_Star
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%srqalsh_star_L=%d_M=%d.out", output_folder, L, M);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	int   io_cost = -1;

	printf("c-k-AFN Search by RQALSH_Star: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; ++num) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[num];
		MaxK_List *list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		io_cost = 0;
		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			io_cost += lsh->kfn(top_k, (const float *) query[i], 
				data_folder, list);
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;
		io_cost       = (int) ceil((float) io_cost / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_rqalsh(				// indexing of rqalsh
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read dataset
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  write the data set in new format to disk
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	write_data_new_form(n, d, B, (const float **) data, data_folder);

	gettimeofday(&end_time, NULL);
	float write_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Write Dataset in New Format: %f Seconds\n\n", write_file_time);

	// -------------------------------------------------------------------------
	//  indexing of RQALSH
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "rqalsh/");

	RQALSH *lsh = new RQALSH();
	lsh->build(n, d, B, beta, delta, ratio, index_path, (const float **) data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  write indexing time to disk
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, output_folder);
	strcat(fname, "rqalsh.index");

	FILE *fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;
	}
	fprintf(fp, "Indexing Time: %f seconds\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_rqalsh(					// c-k-AFN search of RQALSH
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  load RQALSH
	// -------------------------------------------------------------------------
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "rqalsh/");

	RQALSH *lsh = new RQALSH();
	if (lsh->load(index_path)) {
		printf("Could not load RQALSH\n");
		exit(1);
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search by RQALSH
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "rqalsh.out");
	
	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	int   io_cost = -1;

	printf("c-k-AFN Search by RQALSH: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; ++num) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[num];
		MaxK_List *list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		io_cost = 0;
		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			io_cost += lsh->kfn(top_k, (const float *) query[i], 
				data_folder, list);
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;
		io_cost       = (int) ceil((float) io_cost / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  k-FN search by Linear Scan
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "linear.out");

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	long long io_cost = -1;

	printf("k-FN Search by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < maxRound; round++) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[round];
		MaxK_List *list = new MaxK_List(top_k);
		
		overall_ratio = 0.0f;
		recall = 0.0f;
		io_cost = 0;
		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			io_cost += linear(n, d, B, top_k, (const float*) query[i], 
				data_folder, list);
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;
		io_cost       = (int) ceil((float) io_cost / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f%%\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;
	
	return 0;
}
