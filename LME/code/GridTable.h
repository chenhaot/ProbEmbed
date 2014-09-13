#ifndef GRID_TABLE
#define GRID_TABLE

typedef struct index
{
    int idx;
    struct index* pnext;
}
INDEX;

typedef struct gridcoordinate
{
    INDEX* first_idx_struct;
    int* coor;
}
GRIDCOORDINATE;

typedef struct gridtable
{
    int d;
    int num_slot;
    int num_grids_per_axis;
    GRIDCOORDINATE* gridarray;
    double* min_vals;
    double* max_vals;
    double* intervals;
}
GRIDTABLE;

void Init_Grid_Table(GRIDTABLE* gt, double** X, int k, int d, int num_grids_per_axis);
void Free_Grid_Table(GRIDTABLE* gt);
int is_null_coordinate(GRIDCOORDINATE* gc);
void calc_grid_coordinate(GRIDTABLE* gt, double* x, int* grid_coor);
int coordinate_hash(int* grid_coor, int d, int num_slot);
int insert_into_grid_table(GRIDTABLE* gt, int* grid_coor, int idx);
int slot_not_used(GRIDTABLE* gt, int hash_val);
int same_coordinate(int* coor1, int* coor2, int d);
int* get_indices(GRIDTABLE* gt, int* grid_coor, int* num_indices);
int* get_nearby_indices(GRIDTABLE* gt, double* x, int* num_indices, int grid_radius);
int is_in_grid(GRIDTABLE* gt, int* coor);
int* append_lists(int* total_length, int* list1, int length1, int* list2, int length2);
void print_list(int* list, int length);
void convert_num_to_coor(int num, int n, int d, int* coor);

#endif
