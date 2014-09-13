#ifndef LOGISTIC_COMMON
#define LOGISTIC_COMMON


typedef struct PlaylistData
{
    int num_songs;
    int num_playlists;
    int num_appearance;
    int* all_ids;
    int* id_counts;
    int* playlists_length;
    int** playlists;
}
PDATA;

typedef struct TagData
{
    int num_songs;
    int num_tags;
    int* num_tags_for_song;
    int** tags;
}
TDATA;

typedef struct Paras
{
    //Normalize the 2-norm of each embedding to 1 or not
    int do_normalization;
    int method;
    int d;
    double ita;
    double eps;
    double lambda;
    double nu_multiplier;
    int random_init;
    int fast_collection;
    double alpha;
    double beta;
    int least_iter;
    int radius;
    int num_points;
    int stoc_grad;
    int allow_self_transition;
    int regularization_type;
    int output_distr; //Whether output the log likelihood distribution files agains song and transition
    int grid_heuristic;
    int regeneration_interval;
    int bias_enabled;
    int num_llhood_track;
    char tagfile[200];
    int tag_regularizer;
    int hessian;
    int landmark_heuristic;
    int num_landmark;
    double lowerbound_ratio;
    int reboot_enabled;
    int landmark_burnin_iter;
    int use_hash_TTable;
}
PARAS;

typedef struct NeighborNode
{
    int id;
    struct NeighborNode* pnext;
}
NNODE;

typedef struct NeighborList
{
    int length;
    NNODE* pheader;
}
NLIST;


PDATA read_playlists_data(char* filename);
void free_playlists_data(PDATA pd);
void print_playlists_data(PDATA pd);
TDATA read_tag_data(char* filename);
void free_tag_data(TDATA td);
void print_tag_data(TDATA td);
void Array2Dcopy(double** src, double** dest, int m, int n);
void Veccopy(double* src, double* dest, int length);
void Array2Dfree(double** X, int m, int n);
double** zerosarray(int m, int n);
double** randarray(int m, int n, double c);
double* randvec(int m, double c);
void Array2Dexclcopy(double** src, double** dest, int m, int n, int excl_idx);
double innerprod(double* vec1, double* vec2, int length);
void exp_on_vec(double* vec, int length);
void inverse_on_vec(double* vec, int length);   //ADITH: To support inversion
void log_on_vec(double* vec, int length);
void sum_along_direct(double** X, double* vec, int m, int n, int direct);
double mat_norm_diff(double** X1, double** X2, int m, int n);
void add_vec(double* org_vec, double* update_vec, int length, double coeff);
void add_mat(double** org_mat, double** update_mat, int m, int n, double coeff);
void tile_vec(double* vec, double** mat, int m, int n, int direct);
void mat_mult(double** src1, double** src2, double** dest, int m, int n);
void mat_neg_euc(double** src1, double** src2, double** dest, int m, int n);
double sum_vec(double* vec, int length);
void scale_vec(double* vec, int length, double scale);
void scale_mat(double** X, int m, int n, double scale);
double frob_norm(double** X, int m, int n);
double vec_norm(double* vec, int length);
double avg_norm(double** X, int m, int n);
double Cal_log_likelihood(double** X, PDATA pd, int d);
void write_embedding_to_file(double** X, int m, int n, char* filename, double* bias_terms);
int split_double_line(char* str, double* parray);
char* remove_eof(char* str);
double** read_matrix_file(char* filename, int* m, int* n, double** bias_terms);
void Cal_exp_affinity_mat(double** X, double** A, int m, int n);
void Vecexccopy(double* src, double* dest, int length, int exc);
void Update_affinity_mat(double**X, double** A, int idx, int m, int n);
void insert_in_Nlist(NLIST* nl, int to_insert);
void free_Nlist(NLIST* nl);
//void free_Nnode(NNODE* pnn);
NLIST** find_neighbors(PDATA pd);
void print_nlist(NLIST* p);
void Array2Dcopy_nn(double** src, double** dest, int d, NLIST* p);

double vec_sq_dist(double* x1, double* x2, int length);
void Array2Dexclcopy_minus(double** src, double** dest, int m, int n, int excl_idx, double* bigx);
void vec_mat_sum(double* vec, double vec_c, double** mat, double mat_c, double** dest, int m, int n);
void vec_scalar_sum(double* vec, double sc, int length);

int* merge_two_lists(int* total_length, int* list1, int length1, int* list2, int length2);
int exist_in_list(int target, int* list, int length);
int exist_in_nonneg_list(int target, int* list, int length);

void calculate_realX(double** X, double** realX, TDATA td, int k, int m, int d, int num_points);

int* kmeans(double** X, int m, int n, int k, double eps, int verb);

int* read_hash(char* hash_filename, int* num_idx);
void calculate_realX_with_hash(double** X, double** realX, TDATA td, int k, int m, int d, int num_points, int k_train, int* hash);

int* get_test_ids(int k, int k_train, int* train_ids);
void int_list_copy(int* src, int* dest, int length);

#endif
