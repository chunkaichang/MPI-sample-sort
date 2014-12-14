#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <float.h>
#include <sys/time.h>
#include <limits.h>
#include <stdbool.h>
#define SEED 100
#define OUTPUT 0
#define CHECK 1
//Sequential sampleSort.  
//Assume size is a multiple of nbuckets*nbuckets

double get_clock() {
   struct timeval tv; int ok;
   ok = gettimeofday(&tv, (void *) 0);
   if (ok<0) { printf("gettimeofday error");  }
   return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6); 
}

int compare(const void *num1, const void *num2) {
	unsigned long long* n1 = (unsigned long long*)num1;
	unsigned long long* n2 = (unsigned long long*)num2;
	return (*n1 > *n2) - (*n1 < *n2);
}

int main(int argc, char *argv[]) {
	int i,j,size,nbuckets,count;
	double t1,t2;
	unsigned long long *splitters,*elmnts,*sample,**buckets; 
	//unsigned long long **glob_buckets;
	unsigned long long *my_elmnts,*my_sample;
	int my_rank,num_proc;
	unsigned long long check, check_before,glob_check;
	int *bucket_sizes;
	//int *glob_bucket_sizes, *recv_count, *displ;
	bool checkMax, glob_checkMax;

	if (argc != 2) {
		fprintf(stderr,
			"Wrong number of arguments.\nUsage: %s N\n",
			argv[0]);
		return -1;
	}
	
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	nbuckets = num_proc;
	size = atoi(argv[1]);

	// Allocate memory on each processor
	if(my_rank == 0)
	{
		elmnts = (unsigned long long*)malloc(sizeof(unsigned long long)*size);
		//glob_checkMax = 1;
		//glob_bucket_sizes = (int*)malloc(sizeof(int)*nbuckets);
		
		//recv_count = (int*)malloc(sizeof(int)*nbuckets);

		//displ = (int*)malloc(sizeof(int)*nbuckets);
		
		//glob_buckets = (unsigned long long**)malloc(sizeof(unsigned long long*)*nbuckets);
		/*
		for(i=0;i<nbuckets;i++) {
			glob_bucket_sizes[i] = 0;
			recv_count[i] = 0;
			displ[i] = 0;
			glob_buckets[i] = (unsigned long long*)malloc(sizeof(unsigned long long)*2*size/nbuckets);
		}	
		*/					                
	}
	
	my_elmnts = (unsigned long long*)malloc(sizeof(unsigned long long)*(size/nbuckets));
	my_sample = (unsigned long long*)malloc(sizeof(unsigned long long)*(nbuckets-1));
	
	splitters = (unsigned long long*)malloc(sizeof(unsigned long long)*nbuckets);
	sample = (unsigned long long*)malloc(sizeof(unsigned long long)*nbuckets*(nbuckets-1));
	buckets = (unsigned long long**)malloc(sizeof(unsigned long long*)*nbuckets);
	
	//the size of each bucket is guaranteed to be less than
	//2*size/nbuckets becuase of the way we choose the sample
	for(i=0;i<nbuckets;i++) {
		buckets[i] = (unsigned long long*)malloc(sizeof(unsigned long long)*2*size/nbuckets);
	}
	bucket_sizes = (int*)malloc(sizeof(int)*nbuckets);
	for(i=0;i<nbuckets;i++) {
		bucket_sizes[i] = 0;
	}

	//Rank 0 fills elmnts with random numbers and scatter elmnts to others
	if(my_rank == 0){
		srand(SEED);
		for(i=0;i<size;i++) {
			elmnts[i] = rand()%100;
			//printf("%d ",elmnts[i]);
		}
		//printf("\n");
	}	
	MPI_Scatter(elmnts, size/nbuckets, MPI_UNSIGNED_LONG_LONG,
				my_elmnts, size/nbuckets, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	
	/*
	for(i=0;i<size/nbuckets;i++)
		printf("%d ",my_elmnts[i]);
	
	printf("\n");
	*/
	#if CHECK
	if(my_rank == 0){
	check_before = 0;
	for(i=0;i<size;i++) {
		check_before ^= elmnts[i];
	}
	//printf("Check before:%llu",check_before);
	}
	#endif
	
	MPI_Bcast(&check_before, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	t1 = get_clock();

	//select the samples
	
	qsort(my_elmnts,size/nbuckets,sizeof(unsigned long long),compare);
	
	for(j=0;j<nbuckets-1;j++) {
		my_sample[j] = my_elmnts[size/nbuckets/nbuckets*(j+1)];
	}
	
	//send samples to everyone	
	MPI_Allgather(my_sample, nbuckets-1, MPI_UNSIGNED_LONG_LONG, sample, nbuckets-1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

	//select the splitters
	qsort(sample,nbuckets*(nbuckets-1),sizeof(unsigned long long),compare);
	for(i=1;i<nbuckets;i++) {
		splitters[i-1] = sample[i*(nbuckets-1)];
	}
	splitters[nbuckets-1] = ULLONG_MAX;
	/*
	if(my_rank == 0)
	{
		for(i=0;i<nbuckets*(nbuckets-1);i++)
			printf("Samples: %d", sample[i]);
		printf("\n");	
		for(i=0;i<nbuckets-1;i++)
			printf("Splitter:%d \n",splitters[i]);
	}
	*/
	//put into buckets based on splitters
	for(i=0;i<size/nbuckets;i++) {
		j = 0;
		while(j < nbuckets) {
			if(my_elmnts[i]<splitters[j]) {
				buckets[j][bucket_sizes[j]] = my_elmnts[i];
				bucket_sizes[j]++;
				j = nbuckets;
			}
			j++;
		}
	}

	//Sort each buckets
	for(i=0;i<nbuckets;i++) {
		qsort(buckets[i], bucket_sizes[i], sizeof(unsigned long long), compare);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(my_rank == 0){
	t2 = get_clock();
	printf("Time: %lf\n",(t2-t1));
	}
	/*
	if(my_rank == 1)
	{
		for(i=0;i<nbuckets;i++)
		{
			printf("Bucket:%d\n",i);
			for(j=0;j<bucket_sizes[i];j++)
				printf(" %llu",buckets[i][j]);
			printf("\n");
		}
	}
	*/
	//rank0 gather recv_count and displacement
	//if(my_rank == 0)
	//	printf("Rank 0 bucket1 size: %d\n",bucket_sizes[1]);
	/*
	for(i=0;i<nbuckets;i++){
		MPI_Gather(&bucket_sizes[i], 1, MPI_INT, recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(my_rank == 0){
			//printf("Bucket %d recv_count: ",i);
			for(j=1;j<nbuckets;j++){
				//printf("%d ",recv_count[j]);
				
				displ[j] = displ[j-1] + recv_count[j-1];
				
			}	
			//printf("\n");	
		}
		// send bucket[i] to rank 0
		MPI_Gatherv(buckets[i], bucket_sizes[i], MPI_UNSIGNED_LONG_LONG,
				glob_buckets[i], recv_count, displ, MPI_UNSIGNED_LONG_LONG, 0,  MPI_COMM_WORLD);

		if(my_rank == 0){
			for(j=0;j<nbuckets;j++){
				recv_count[j] = 0;
				displ[j] = 0;
			}
		}
	}	
	
	//reduce bucket_sizes
	MPI_Reduce(bucket_sizes, glob_bucket_sizes, nbuckets, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


	//Rank 0 sorts each bucket
	if(my_rank == 0){
		for(i=0;i<nbuckets;i++) {
			qsort(glob_buckets[i], glob_bucket_sizes[i], sizeof(unsigned long long), compare);
		}
	
	t2 = get_clock();
	printf("Time: %lf\n",(t2-t1));
	}
	*/
	
	#if CHECK
	//if(my_rank == 0){
	count = 0;
	for(i=0;i<nbuckets;i++) {
		//printf("Bucket%d\n",i);
		for(j=0;j<bucket_sizes[i];j++) {
			check ^= buckets[i][j];
			//printf("%llu ",glob_buckets[i][j]);
		}
		count += bucket_sizes[i];
		//printf("\n");
	}
	
	
	//reduce check
	MPI_Reduce(&check, &glob_check, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR, 0, MPI_COMM_WORLD);
	
	if(my_rank == 0){
		printf("The bitwise xor is %llu\n", glob_check^check_before);
	}
	
	checkMax = true;
	for(i=0;i<nbuckets-1;i++) {
		if(buckets[i][bucket_sizes[i]-1] > buckets[i+1][0]) {
			checkMax = false;
		}
	}

	MPI_Reduce(&checkMax, &glob_checkMax, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	
	if(my_rank == 0){
	printf("The max of each bucket is not greater than the min of the next:	%s\n",
		glob_checkMax ? "true" : "false");
	}
	
	//}	
	#endif
	#if OUTPUT
	if(my_rank == 0){
	count = 0;
	for(i=0;i<nbuckets;i++) {
		for(j=0;j<glob_bucket_sizes[i];j++) {
			elmnts[count+j] = buckets[i][j];
		}
		count += bucket_sizes[i];
	}
	for(i=0;i<size;i++) {
		printf("%llu\n",elmnts[i]);
	}
	}
	#endif

	if(my_rank == 0){
		free(elmnts);
		/*
		free(glob_bucket_sizes);
		free(recv_count);
		free(displ);
		for(i=0;i<nbuckets;i++) {
			free(glob_buckets[i]);
	        }
	        free(glob_buckets);
		*/
	}
	free(my_elmnts);
	free(my_sample);
	free(splitters);
	free(sample);
	free(bucket_sizes);
	for(i=0;i<nbuckets;i++) {
		free(buckets[i]);
	}
	free(buckets);
	
	MPI_Finalize();
	return 0;
}
