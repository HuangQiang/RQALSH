#include <algorithm>
#include <cstring>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qdafn.h"
#include "drusilla_select.h"
#include "rqalsh.h"
#include "rqalsh_star.h"
#include "afn.h"


// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "linear.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  k-FN search by Linear Scan
	// -------------------------------------------------------------------------
	printf("Top-k FN Search by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);
		
		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += linear(n, d, B, top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;
		g_io      = (int) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f%%\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_rqalsh_star(		// indexing of RQALSH*
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	sprintf(output_set, "%srqalsh_star.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing of RQALSH*
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	RQALSH_Star *lsh = new RQALSH_Star();
	lsh->build(n, d, B, L, M, beta, delta, ratio, data, output_folder);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);
	fclose(fp);

	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_rqalsh_star(				// c-k-AFN search of RQALSH*
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	sprintf(output_set, "%srqalsh_star.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  load RQALSH*
	// -------------------------------------------------------------------------
	RQALSH_Star *lsh = new RQALSH_Star();
	if (lsh->load(output_folder)) return 1;
	lsh->display();

	// -------------------------------------------------------------------------
	//  c-k-AFN search by RQALSH*
	// -------------------------------------------------------------------------
	printf("Top-k FN Search by RQALSH*: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += lsh->kfn(top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;
		g_io      = (int) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_rqalsh(				// indexing of RQALSH
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "rqalsh.out");
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing of RQALSH
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "indices/");

	RQALSH *lsh = new RQALSH();
	lsh->build(n, d, B, beta, delta, ratio, data, index_path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);
	fclose(fp);

	delete lsh; lsh = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_rqalsh(					// c-k-AFN search of RQALSH
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "rqalsh.out");
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  load RQALSH
	// -------------------------------------------------------------------------
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "indices/");

	RQALSH *lsh = new RQALSH();
	if (lsh->load(index_path)) return 1;
	lsh->display();

	// -------------------------------------------------------------------------
	//  c-k-AFN search by RQALSH
	// -------------------------------------------------------------------------
	printf("Top-k FN Search by RQALSH: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += lsh->kfn(top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;
		g_io      = (int) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_drusilla_select(	// indexing of Drusilla_Select
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   L,							// number of projection
	int   M,							// number of candidates
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "drusilla.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  Indexing of Drusilla_Select
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Drusilla_Select* drusilla = new Drusilla_Select();
	drusilla->build(n, d, L, M, B, data, output_folder);
	drusilla->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);
	fclose(fp);

	delete drusilla; drusilla = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_drusilla_select(			// c-k-AFN via Drusilla_Select
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "drusilla.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  load index of Drusilla_Select
	// -------------------------------------------------------------------------
	Drusilla_Select *drusilla = new Drusilla_Select();
	if (drusilla->load(output_folder)) return 1;
	drusilla->display();

	// -------------------------------------------------------------------------
	//  c-k-AFN search via Drusilla_Select
	// -------------------------------------------------------------------------
	printf("Top-k FN Search by Drusilla_Select: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += drusilla->search(query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;
		g_io      = (int) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete drusilla; drusilla = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_qdafn(				// indexing of QDAFN
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	int   L,							// number of projections
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qdafn.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  Indexing of QDAFN
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char index_path[200];
	sprintf(index_path, "%sindices/", output_folder);

	QDAFN* qdafn = new QDAFN();
	qdafn->build(n, d, B, L, M, ratio, data, index_path);
	qdafn->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete qdafn; qdafn = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kfn_of_qdafn(					// c-k-AFN via QDAFN
	int   qn,							// number of query points
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qdafn.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  load index of QDAFN
	// -------------------------------------------------------------------------
	char index_path[200];
	sprintf(index_path, "%sindices/", output_folder);

	QDAFN *qdafn = new QDAFN();
	if (qdafn->load(index_path)) return 1;
	qdafn->display();

	// -------------------------------------------------------------------------
	//  c-k-AFN search via QDAFN
	// -------------------------------------------------------------------------
	printf("Top-k FN Search by QDAFN: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += qdafn->search(top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;
		g_io      = (int) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete qdafn; qdafn = NULL;

	return 0;
}
