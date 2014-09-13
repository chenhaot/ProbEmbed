#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LogisticEmbed_common.h"
#include "LandmarkTable.h"

void Init_LMTable(LMTABLE* plmt, double** X, int k, int d, int lm_num, int nn_num)
{
    int i;
    int j;
    int s;
    int t;
    plmt -> lm_num = lm_num;
    plmt -> nn_num = nn_num;
    plmt -> k = k;

    if(lm_num < 0)
    {
	lm_num = k;
	plmt -> lm_num = k;
    }


    if(lm_num > k)
    {
	printf("Error: number of the landmarks is greater than the total number of songs.\n");
	exit(1);
    }

    if(nn_num > lm_num)
    {
	printf("Error: number of nearest landmarks is greater thant the total number of landmarks.\n");
	exit(1);
    }

    if(lm_num == k)
    {
	printf("..........................!!!!!!!!!!\n");
	(plmt -> landmarks_idx) = (int*)malloc(lm_num * sizeof(int));
	for(i = 0; i < lm_num; i++)
	    (plmt -> landmarks_idx)[i] = i;
    }
    else
    {
	srand(time(NULL));
	(plmt -> landmarks_idx) = (int*)malloc(lm_num * sizeof(int));
	for(i = 0; i < lm_num; i++)
	    (plmt -> landmarks_idx[i]) = -1;
	j = 0;
	while(j < lm_num)
	{
	    int temp_idx = rand() % k;
	    if(!exist_in_list(temp_idx, plmt -> landmarks_idx, lm_num))
	    {
		(plmt -> landmarks_idx)[j] = temp_idx;
		j++;
	    }
	}
    }


    double** landmarks = zerosarray(lm_num, d); 
    for(i = 0; i < lm_num; i++)
	Veccopy(X[(plmt -> landmarks_idx)[i]], landmarks[i], d);

    plmt -> song_lm_hash = (int**)malloc(k * sizeof(int*));

    plmt -> lm_song_hash = (int**)malloc(lm_num * sizeof(int*));
    for(i = 0; i < lm_num; i++)
	(plmt -> lm_song_hash)[i] = (int*)malloc(k * sizeof(int));
    (plmt -> lm_song_length) = (int*)malloc(lm_num * sizeof(int));
    for(i = 0; i < lm_num; i++)
	(plmt -> lm_song_length)[i] = 0;

    for(i = 0; i < k; i++)
    {
	plmt -> song_lm_hash[i] = (int*)malloc(nn_num * sizeof(int));
	for(j = 0; j < nn_num; j++)
	    (plmt -> song_lm_hash[i])[j] = -1;

	double* temp_near_dist = (double*)malloc(nn_num * sizeof(double)); 
	for(j = 0; j < nn_num; j++)
	    temp_near_dist[j] = 10000000.0; 

	for(j = 0; j < lm_num; j++)
	{
	    double temp_dist = vec_sq_dist(X[i], landmarks[j], d); 
	    for(t = 0; t < nn_num; t++)
	    {
		if(temp_dist < temp_near_dist[t])
		{
		    for(s = nn_num - 1; s > t; s--)
		    {
			temp_near_dist[s] = temp_near_dist[s - 1];
			(plmt -> song_lm_hash[i])[s] = (plmt -> song_lm_hash[i])[s - 1];
		    }
		    temp_near_dist[t] = temp_dist;
		    (plmt -> song_lm_hash)[i][t] = j;
		    break;
		}
	    }

	}
	free(temp_near_dist);

	for(j = 0; j < nn_num; j++)
	{
	    int temp_idx = (plmt -> song_lm_hash[i])[j];
	    (plmt -> lm_song_hash)[temp_idx][plmt -> lm_song_length[temp_idx]] = i;
	    (plmt -> lm_song_length[temp_idx]) += 1;
	}
    }

    for(i = 0; i < lm_num; i++)
	(plmt -> lm_song_hash)[i] = (int*)realloc((plmt -> lm_song_hash)[i], (plmt -> lm_song_length[i]) * sizeof(int));

    Array2Dfree(landmarks, lm_num, d);
}

void Init_LMTable_kmeans(LMTABLE* plmt, double** X, int k, int d, int lm_num)
{
    int i;
    int j;
    plmt -> lm_num = lm_num;
    plmt -> nn_num = 1;
    plmt -> k = k;
    plmt -> landmarks_idx = NULL;

    int* assignment =  kmeans(X, k, d, lm_num, 1e-5, 0);
    plmt -> song_lm_hash = (int**)malloc(k * sizeof(int*));
    for(i = 0; i < k; i++)
    {
	plmt -> song_lm_hash[i] = (int*)malloc(sizeof(int));
	plmt -> song_lm_hash[i][0] = assignment[i];
    }

    plmt -> lm_song_length = (int*)calloc(lm_num, sizeof(int));
    plmt -> lm_song_hash = (int**)malloc(lm_num * sizeof(int*));
    for(i = 0; i < lm_num; i++)
	(plmt -> lm_song_hash)[i] = (int*)malloc(k * sizeof(int));
    for(i = 0; i < k; i++)
    {
	int current_lm = assignment[i];
	(plmt -> lm_song_hash)[current_lm][(plmt -> lm_song_length)[current_lm]] = i;
	(plmt -> lm_song_length)[current_lm] += 1;
    }
    for(i = 0; i < lm_num; i++)
	(plmt -> lm_song_hash)[i] = (int*)realloc((plmt -> lm_song_hash)[i], (plmt -> lm_song_length)[i] * sizeof(int));
    free(assignment);
}


void Free_LMTable(LMTABLE* plmt)
{
    int i;
    for(i = 0; i < (plmt -> k); i++)
	free(plmt -> song_lm_hash[i]);
    free(plmt -> song_lm_hash);
    for(i = 0; i < (plmt -> lm_num); i++)
	free(plmt -> lm_song_hash[i]);
    free(plmt -> lm_song_hash);
    free(plmt -> lm_song_length);
    if((plmt -> landmarks_idx) != NULL)
	free(plmt -> landmarks_idx);
}

int* get_neaby_songs_by_landmark(LMTABLE lmt, int num_landmark, int* num_candidates, double** X, double* x, int k, int d)
{
    int i;
    int j;
    int t;
    int s;

    if(num_landmark > lmt.lm_num)
    {
	printf("Error the required number of nearest landmarks is greater than the total number of landmarks.\n");
	exit(1);
    }

    double* temp_near_dist = (double*)malloc(num_landmark * sizeof(double)); 
    for(j = 0; j < num_landmark; j++)
	temp_near_dist[j] = 10000000.0; 
    int* temp_landmark_idx = (int*)malloc(num_landmark * sizeof(int));
    for(j = 0; j < num_landmark; j++)
	temp_landmark_idx[j] = -1;

    for(j = 0; j < lmt.lm_num; j++)
    {
	double temp_dist = vec_sq_dist(X[lmt.landmarks_idx[j]], x, d); 
	for(t = 0; t < num_landmark; t++)
	{
	    if(temp_dist < temp_near_dist[t])
	    {
		for(s = num_landmark - 1; s > t; s--)
		{
		    temp_near_dist[s] = temp_near_dist[s - 1];
		    temp_landmark_idx[s] = temp_landmark_idx[s - 1];
		}
		temp_near_dist[t] = temp_dist;
		temp_landmark_idx[t] = j;
		break;
	    }
	}

    }


    *num_candidates = 0;
    for(i = 0; i < num_landmark; i++)
	(*num_candidates) += lmt.lm_song_length[temp_landmark_idx[i]];

    int* nearby_songs = (int*)malloc((*num_candidates) * sizeof(int));
    
    t = 0;
    for(i = 0; i < num_landmark; i++)
    {
	for(j = 0; j < lmt.lm_song_length[temp_landmark_idx[i]]; j++)
	{
	    nearby_songs[t] = lmt.lm_song_hash[temp_landmark_idx[i]][j];
	    t++;
	}
    }
    free(temp_near_dist);
    free(temp_landmark_idx);
    return nearby_songs;
}

int* get_neaby_songs_by_landmark_with_lowerbound(LMTABLE lmt, int lower_bound, int* num_candidates, double** X, double* x, int k, int d)
{
    int i;
    int j;
    int t;
    int s;

    if(lower_bound > k)
    {
	printf("Error: the required number of songs that can be transitioned to is greater thant the total number of songs\n");
	exit(1);
    }

    double* temp_near_dist = (double*)malloc(lmt.lm_num * sizeof(double)); 
    for(j = 0; j < lmt.lm_num; j++)
	temp_near_dist[j] = 10000000.0; 
    int* temp_landmark_idx = (int*)malloc(lmt.lm_num * sizeof(int));
    for(j = 0; j < lmt.lm_num; j++)
	temp_landmark_idx[j] = -1;

    for(j = 0; j < lmt.lm_num; j++)
    {
	double temp_dist = vec_sq_dist(X[lmt.landmarks_idx[j]], x, d); 
	for(t = 0; t < lmt.lm_num; t++)
	{
	    if(temp_dist < temp_near_dist[t])
	    {
		for(s = lmt.lm_num - 1; s > t; s--)
		{
		    temp_near_dist[s] = temp_near_dist[s - 1];
		    temp_landmark_idx[s] = temp_landmark_idx[s - 1];
		}
		temp_near_dist[t] = temp_dist;
		temp_landmark_idx[t] = j;
		break;
	    }
	}

    }


    *num_candidates = 0;
    for(i = 0; i < lmt.lm_num; i++)
    {
	(*num_candidates) += lmt.lm_song_length[temp_landmark_idx[i]];
	if((*num_candidates) >= lower_bound)
	{
	    s = i;
	    break;
	}
    }

    int* nearby_songs = (int*)malloc((*num_candidates) * sizeof(int));
    
    t = 0;
    for(i = 0; i <= s; i++)
    {
	for(j = 0; j < lmt.lm_song_length[temp_landmark_idx[i]]; j++)
	{
	    nearby_songs[t] = lmt.lm_song_hash[temp_landmark_idx[i]][j];
	    t++;
	}
    }
    free(temp_near_dist);
    free(temp_landmark_idx);
    return nearby_songs;
}
