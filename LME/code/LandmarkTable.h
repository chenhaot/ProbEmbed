#ifndef LANDMARK_TABLE
#define LANDMARK_TABLE


typedef struct landmarktable
{
    int k;
    int lm_num; //number of landmarks
    int nn_num; //number of nearest neighbor landmarks
    int** song_lm_hash;
    int** lm_song_hash;
    int* lm_song_length;
    int* landmarks_idx;
}
LMTABLE;



void Init_LMTable(LMTABLE* plmt, double** X, int k, int d, int lm_num, int nn_num);
void Init_LMTable_kmeans(LMTABLE* plmt, double** X, int k, int d, int lm_num);
int* get_neaby_songs_by_landmark(LMTABLE lmt, int num_landmark, int* num_candidates, double** X, double* x, int k, int d);
int* get_neaby_songs_by_landmark_with_lowerbound(LMTABLE lmt, int lower_bound, int* num_candidates, double** X, double* x, int k, int d);
void Free_LMTable(LMTABLE* plmt);
#endif
