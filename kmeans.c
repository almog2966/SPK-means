#include <stdio.h>
#include <stdlib.h>
#include <math.h>





double* kmeans_calc(int k, int n_rows, int n_col, int max_iter, double e, double* dataPoints, double *centroid)
{
    double *clusters;
    int *counter;
    int iter,i,j,r,c,min_centroid_index,first;
    double min_dist,curr_dist,norm_dist;
    clusters = (double *)malloc(k * n_col * sizeof(double));
    counter = (int *)calloc(k,sizeof(int));
    iter = 0;
    if(clusters == NULL || counter == NULL)
    {
        printf("An Error Has Occurred");
        free(clusters);
        free(counter);
        return(NULL);
    }
    do{
        /*Reset Data Before iteration*/
        for(i=0; i<k;i++)
        {
            for(j=0; j<n_col; j++)
            {
                clusters[i*n_col +j] = 0;
            }
            counter[i] = 0;
        }

        for(i=0; i<n_rows; i++)
        {   
            /*the i row of dataPoints*/
            first = 0;
            min_dist = 0;
            min_centroid_index = 0;
            for(r=0; r<k; r++)
            {
                /*the R row of centroids*/
                curr_dist = 0;
                /*calculating the min distance*/
                for(c=0; c<n_col; c++)
                {
                    /*the C column*/
                    curr_dist += ((dataPoints[i*n_col+c] - centroid[r*n_col+c])*(dataPoints[i*n_col+c] - centroid[r*n_col+c]));
                }
                if(first == 0 || curr_dist<min_dist)
                {
                    first = 1;
                    min_dist = curr_dist;
                    min_centroid_index = r;
                }
            }
            counter[min_centroid_index]++;
            /*updating the cluster*/
            for(c=0; c<n_col; c++)
            {
                clusters[min_centroid_index*n_col + c] += dataPoints[i*n_col+c];
            }

        }
        /*update cluster to to new centroids*/
        for(i=0; i<k; i++)
        {
            for(j=0; j<n_col; j++)
            {
                clusters[i*n_col+j] = clusters[i*n_col+j]/counter[i];
            }
        }
        /*calculating normal distacne*/
        norm_dist = 0;
        for(i = 0; i<k*n_col; i++)
        {
            norm_dist += ((centroid[i] - clusters[i])*(centroid[i] - clusters[i]));
        }
        norm_dist = sqrt(norm_dist);
        for(i=0; i<k*n_col; i++)
        {
            centroid[i] = clusters[i];
        }
        iter++;
    } while (iter<max_iter && norm_dist >= e);

    

    free(clusters);
    free(counter);

    return centroid;

    


}
