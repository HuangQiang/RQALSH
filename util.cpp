#include "headers.h"


// -----------------------------------------------------------------------------
//  Global variables
// -----------------------------------------------------------------------------
long g_memory = 0;


// -----------------------------------------------------------------------------
//  Utility functions
// -----------------------------------------------------------------------------
int compfloats(						// compare two real values 
	float v1,							// 1st value (type float)
	float v2)							// 2nd value (type float)
{
	if (v1 - v2 < -FLOATZERO) return -1;
	else if (v1 - v2 > FLOATZERO) return 1;
	else return 0;
}

// -----------------------------------------------------------------------------
void error(							// display an error message
	char* msg,							// message
	bool is_exit)						// whether exit
{
	printf(msg);
	if (is_exit) exit(1);
}

// -----------------------------------------------------------------------------
int check_mem()						// check memory is enough
{
	if (g_memory > MAXMEMORY) {		// default maximum memory is 1 GB
		printf("I am going to need around %.2f MB memory.\n", 
			(float) g_memory / (1024.0f * 1024.0f));

		printf("Can you afford it (y/n)?\n");
		char c = getchar();			// ask for more memory
		getchar();
		if (c != 'y' && c != 'Y') {
			return 1;				// fail to return
		}
	}
	return 0;
}

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char* data_path)					// input path
{
#ifdef LINUX_						// create directory under Linux
	int len = (int) strlen(data_path);
	for (int i = 0; i < len; i++) {
		if (data_path[i] == '/') {
			char ch = data_path[i + 1];
			data_path[i + 1] = '\0';
									// check whether the directory exists
			int ret = access(data_path, F_OK);
			if (ret != 0) {			// create the directory
				ret = mkdir(data_path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", data_path);
					error("", true);
				}
			}
			data_path[i + 1] = ch;
		}
	}
#else								// create directory under Windows
	int len = (int) strlen(data_path);
	for (int i = 0; i < len; i++) {
		if (data_path[i] == '/') {
			char ch = data_path[i + 1];
			data_path[i + 1] = '\0';
									// check whether the directory exists
			int ret = _access(data_path, 0);
			if (ret != 0) {			// create the directory
				ret = _mkdir(data_path);
				if (ret != 0) {
					printf("Could not create directory %s\n", data_path);
					error("", true);
				}
			}
			data_path[i + 1] = ch;
		}
	}
#endif
}

// -----------------------------------------------------------------------------
float calc_l2_dist(				// calc L_2 distance (type float)
	float* vec1,						// 1st point
	float* vec2,						// 2nd point
	int    dim)							// dimension
{
	float diff = 0.0F;
	float ret  = 0.0F;

	for (int i = 0; i < dim; i++) {
		diff = vec1[i] - vec2[i];
		ret += diff * diff;
	}
	return sqrt(ret);
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall between two ID list
	int n,								// length of ID list
	MaxK_List* list,					// ID list of k-nn method
	int* ID)							// ground truth ID list
{
	int count = 0;
	for (int i = 0; i < n; i++) {
		int obj_id = list->ith_largest_id(i);

		for (int j = 0; j < n; j++) {
			if (obj_id == ID[j]) {
				count++;
				break;
			}
		}
	}
	float recall = (float) count / (float) n;
	return recall;
}

// -----------------------------------------------------------------------------
void calc_relative_constrast(		// calc relative constrast
	int n,								// number of data objects
	int qn,								// number of queries
	int d,								// dimensionality
	FILE* fp,							// file pointer
	float** data,						// data objects
	float** query)						// queries
{
	int i = 0;
	int j = 0;

	g_memory += (SIZEFLOAT * qn * 3);
	float* avg_dist = new float[qn];
	float* min_dist = new float[qn];
	float* max_dist = new float[qn];
	for (i = 0; i < qn; i++) {
		avg_dist[i] = 0.0;
		min_dist[i] = MAX_FLT;
		max_dist[i] = MIN_FLT;
	}

	float dist   = -1.0F;
	float rc_min = -1.0F;
	float rc_max = -1.0F;

	//printf("ID\tMin\t\tAvg\t\tMax\t\tRC_NN\tRC_FN\n");
	fprintf(fp, "ID\tMin\tAvg\tMax\tRC_NN\tRC_FN\n");
	for (i = 0; i < qn; i++) {
		avg_dist[i] = 0.0;
		for (j = 0; j < n; j++) {	// find nn object of query
			dist = calc_l2_dist(data[j], query[i], d);

			avg_dist[i] += dist;
			if (min_dist[i] > dist) min_dist[i] = dist;
			if (max_dist[i] < dist) max_dist[i] = dist;
		}
		avg_dist[i] /= (float) n;
		rc_min = avg_dist[i] / min_dist[i];
		rc_max = max_dist[i] / avg_dist[i];

		//printf("%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t%.2f\n", 
		//	i+1, min_dist[i], avg_dist[i], max_dist[i], rc_min, rc_max);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\n", 
			i+1, min_dist[i], avg_dist[i], max_dist[i], rc_min, rc_max);
	}

	float min = 0.0F, max = 0.0F, avg = 0.0F;
	for (i = 0; i < qn; i++) {
		min += min_dist[i]; max += max_dist[i]; avg += avg_dist[i];
	}
	min /= qn; max /= qn; avg /= qn;

	rc_min = avg / min;
	rc_max = max / avg;
	//printf("Total\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t%.2f\n", 
	//	min, avg, max, rc_min, rc_max);
	fprintf(fp, "Avg\t%f\t%f\t%f\t%f\t%f\n", 
		min, avg, max, rc_min, rc_max);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (avg_dist != NULL && min_dist != NULL && max_dist != NULL) {
		delete[] avg_dist; avg_dist = NULL;
		delete[] min_dist; min_dist = NULL;
		delete[] max_dist; max_dist = NULL;
		g_memory -= (SIZEFLOAT * qn * 3);
	}
}	

// -----------------------------------------------------------------------------
//  Read the original dataset (or query set) from disk
// -----------------------------------------------------------------------------
int read_set(						// read (data or query) set from disk
	int n,								// number of data points
	int d,								// dimensionality
	char* set,							// address of dataset
	float** points)						// data or queries (return)
{
	int i = 0;
	int j = 0;
	FILE* fp = NULL;

	fp = fopen(set, "r");			// open data file
	if (!fp) {
		printf("I could not open %s.\n", set);
		return 1;
	}

	i = 0;
	while (!feof(fp) && i < n) {	// read data file
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; j++) {
			fscanf(fp, " %f", &points[i][j]);
		}
		fscanf(fp, "\n");

		i++;
	}
	if (!feof(fp) && i == n) {		// check the size of set
		error("The size of set is larger than you input\n", false);
	}
	else if (feof(fp) && i < n) {
		printf("Set the size of dataset to be %d. ", i);
		error("And try again\n", true);
	}
	fclose(fp);						// close data file
	return 0;
}

// -----------------------------------------------------------------------------
//  Write the data set in new format to the disk
// -----------------------------------------------------------------------------
int write_data_new_form(			// write dataset with new format
	int n,								// cardinality
	int d,								// dimensionality
	int B,								// page size
	float** data,						// data set
	char* output_path)					// output path
{
									// number of data in one data file
	int num = (int) floor((float) B / (d * SIZEFLOAT));
									// total number of data file
	int total_file = (int) ceil((float) n / num);
	if (total_file == 0) {
		return 1;					// fail to return
	}

	// -------------------------------------------------------------------------
	//  Check whether the directory exists. If the directory does not exist, we 
	//  create the directory for each folder.
	// -------------------------------------------------------------------------
	char* data_path = new char[200];
	g_memory += SIZECHAR * 200;

	strcpy(data_path, output_path);
	strcat(data_path, "data/");
	create_dir(data_path);

	// -------------------------------------------------------------------------
	//  Write data of qalsh
	// -------------------------------------------------------------------------
	char* fname = new char[200];
	char* buffer = new char[B];		// allocate one page size
	for (int i = 0; i < B; i++) {
		buffer[i] = 0;
	}
	g_memory += SIZECHAR * (200 + B);

	int left  = 0;
	int right = 0;
	for (int i = 0; i < total_file; i++) {
									// get file name of data
		get_data_filename(i, data_path, fname);

		left = i * num;
		right = left + num;
		if (right > n) right = n;	// write data to buffer
		write_data_to_buffer(d, left, right, data, buffer);

									// write one page of data to disk
		if (write_buffer_to_page(B, fname, buffer) == 1) {
			error("write_data_new_form error to write a page", true);
		}
	}

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (buffer != NULL) {
		delete[] buffer; buffer = NULL;
		g_memory -= SIZECHAR * 200;
	}
	if (data_path != NULL || fname != NULL) {
		delete[] data_path; data_path = NULL;
		delete[] fname; fname = NULL;
		g_memory -= SIZECHAR * (200 + B);
	}
	return 0;
}

// -----------------------------------------------------------------------------
void get_data_filename(				// get file name of data
	int file_id,						// data file id
	char* data_path,					// path to store data in new format
	char* fname)						// file name of data (return)
{
	char c[20];

	strcpy(fname, data_path);
	sprintf(c, "%d", file_id);
	strcat(fname, c);
	strcat(fname, ".data");
}

// -----------------------------------------------------------------------------
void write_data_to_buffer(			// write data to buffer
	int d,								// dimensionality
	int left,							// left data id
	int right,							// right data id
	float** data,						// data set
	char* buffer)						// buffer to store data
{
	int c = 0;
	for (int i = left; i < right; i++) {
		for (int j = 0; j < d; j++) {
			memcpy(&buffer[c], &data[i][j], SIZEFLOAT);
			c += SIZEFLOAT;
		}
	}
}

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write buffer to one page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer)						// buffer to store data
{
	if (fname == NULL || buffer == NULL) return 1;

	FILE* fp = NULL;
	fp = fopen(fname, "wb");		// open data file to write
	if (!fp) {
		printf("I could not create %s.\n", fname);
		return 1;					// fail to return
	}

	fwrite(buffer, B, 1, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
//  Read data in new format from disk
// -----------------------------------------------------------------------------
int read_data(						// read data from page
	int id,								// index of data
	int d,								// dimensionality
	int B,								// page size
	float* data,						// real data (return)
	char* output_path)					// output path
{
	// -------------------------------------------------------------------------
	//  Get file name of data
	// -------------------------------------------------------------------------
	char* fname     = new char[200];
	char* data_path = new char[200];
	g_memory += SIZECHAR * 400;

	strcpy(data_path, output_path);
	strcat(data_path, "data/");
									// number of data in one data file
	int num = (int) floor((float) B / (d * SIZEFLOAT));
									// data file id
	int file_id = (int) floor((float) id / num);

	get_data_filename(file_id, data_path, fname);

	// -------------------------------------------------------------------------
	//  Read buffer (one page of data) in new format from disk
	// -------------------------------------------------------------------------
	char* buffer = new char[B];		// allocate one page size
	g_memory += SIZECHAR * B;
	for (int i = 0; i < B; i++) {
		buffer[i] = 0;
	}
	if (read_buffer_from_page(B, fname, buffer) == 1) {
		error("read_data() error to read a page", true);
	}

	// -------------------------------------------------------------------------
	//  Read data from buffer
	// -------------------------------------------------------------------------
	int index = id % num;
	read_data_from_buffer(index, d, data, buffer);

	if (buffer != NULL) {
		delete[] buffer; buffer = NULL;
		g_memory -= SIZECHAR * B;
	}
	if (data_path != NULL || fname != NULL) {
		delete[] data_path; data_path = NULL;
		delete[] fname; fname = NULL;
		g_memory -= SIZECHAR * 400;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read buffer from page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer)						// buffer to store data
{
	if (fname == NULL || buffer == NULL) return 1;

	FILE* fp = NULL;
	fp = fopen(fname, "rb");
	if (!fp) {
		printf("read_buffer_from_page could not open %s.\n", fname);
		return 1;					// fail to return
	}

	fread(buffer, B, 1, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
void read_data_from_buffer(			// read data from buffer
	int index,							// index of data in buffer
	int d,								// dimensionality
	float* data,						// data set
	char* buffer)						// buffer to store data
{
	int c = index * d * SIZEFLOAT;
	for (int i = 0; i < d; i++) {
		memcpy(&data[i], &buffer[c], SIZEFLOAT);
		c += SIZEFLOAT;
	}
}

