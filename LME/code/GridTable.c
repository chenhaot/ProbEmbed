#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GridTable.h"

void Init_Grid_Table(GRIDTABLE* gt, double** X, int k, int d, int num_grids_per_axis)
{
    if(k < 1)
    {
	printf("Empty sample matrix.\n");
	exit(1);
    }

    int i;
    int j;

    gt -> d = d;

    (gt -> min_vals) = (double*)malloc(gt -> d * sizeof(double));
    (gt -> max_vals) = (double*)malloc(gt -> d * sizeof(double));
    (gt -> intervals) = (double*)malloc(gt -> d * sizeof(double));

    //find extreeme values of each columns
    for(i = 0; i < gt -> d; i++)
    {
	(gt -> min_vals)[i] = X[0][i];
	(gt -> max_vals)[i] = X[0][i];
    }
    for(i = 0; i < k; i++)
    {
	for(j = 0; j < gt -> d; j++)
	{
	    (gt -> min_vals[j]) = X[i][j] < (gt -> min_vals[j]) ? X[i][j] : (gt -> min_vals)[j];
	    (gt -> max_vals[j]) = X[i][j] > (gt -> max_vals[j]) ? X[i][j] : (gt -> max_vals[j]);
	}
    }

    for(i = 0; i < gt -> d; i++)
	gt -> intervals[i] = ((gt -> max_vals)[i] - (gt -> min_vals)[i]) / ((double)num_grids_per_axis);

    gt -> num_slot = k;
    gt -> num_grids_per_axis = num_grids_per_axis;

    gt -> gridarray = (GRIDCOORDINATE*)malloc(gt -> num_slot * sizeof(GRIDCOORDINATE));
    for(i = 0; i < gt -> num_slot; i++)
	(gt -> gridarray)[i].first_idx_struct = NULL;

    int* temp_coor = (int*) malloc((gt -> d) * sizeof(int));
    for(i = 0; i < k; i++)
    {
	calc_grid_coordinate(gt, X[i], temp_coor);
	insert_into_grid_table(gt, temp_coor, i);
    }

    free(temp_coor);

}

void Free_Grid_Table(GRIDTABLE* gt)
{
    int i;
    for(i = 0; i < gt -> num_slot; i++)
    {
	if(!slot_not_used(gt, i))
	{
	    //printf("Free nonempty slot.\n");
	    INDEX* p = (gt -> gridarray)[i].first_idx_struct;
	    INDEX* tempp; 
	    while(p != NULL)
	    {
		tempp = p;
		p = p -> pnext;
		free(tempp);
	    }

	    free((gt -> gridarray)[i].coor);
	}
    }
    free(gt -> gridarray);
    free(gt -> min_vals);
    free(gt -> max_vals);
    free(gt -> intervals);
}


void calc_grid_coordinate(GRIDTABLE* gt, double* x, int* grid_coor)
{
    int i;
    for(i = 0; i < gt -> d; i++)
    {
	grid_coor[i] = (int)floor((x[i] - (gt -> min_vals[i])) / (gt -> intervals[i]));
	if(grid_coor[i] >= (gt -> num_grids_per_axis))
	    grid_coor[i] = (gt -> num_grids_per_axis - 1);
	if(grid_coor[i] < 0)
	    grid_coor[i] = 0;
	if(grid_coor[i] < 0 || grid_coor[i] >= (gt -> num_grids_per_axis))
	{
	    printf("Error: coordinate %d exceeds the range %d.\n", grid_coor[i], gt -> num_grids_per_axis);
	    printf("Float value %f.\n", (x[i] - (gt -> min_vals[i])) / (gt -> intervals[i]));
	    printf("Interval value %f.\n", gt -> intervals[i]);
	    printf("Value %f.\n", (x[i] - (gt -> min_vals[i])));
	    exit(1);
	}
    }
}

int coordinate_hash(int* grid_coor, int d, int num_slot)
{
    int i;
    int r = 1;
    for(i = 0; i < d; i++)
	if(!grid_coor[i] == 0)
	    r *= grid_coor[i];
    return (r % num_slot);
}

int insert_into_grid_table(GRIDTABLE* gt, int* grid_coor, int idx)
{
    int d = gt -> d;
    int num_slot = gt -> num_slot;
    int i;
    int start_hash_val = coordinate_hash(grid_coor, d, num_slot);
    int hash_val = start_hash_val;

    //Find the right slot to insert
    while((!slot_not_used(gt, hash_val)) && (!same_coordinate(grid_coor, (gt -> gridarray)[hash_val].coor, d)))
    {
	//printf("Collision detected.\n");
	hash_val = (hash_val + 1) % num_slot;
	if(hash_val == start_hash_val)
	{
	    printf("Error: the table is already full.\n");
	    exit(1);
	}
    }

    GRIDCOORDINATE* tempp = (gt -> gridarray) + hash_val;
    if(slot_not_used(gt, hash_val))
    {
	tempp -> first_idx_struct = (INDEX*)malloc(sizeof(INDEX));
	tempp -> first_idx_struct -> idx = idx;
	tempp -> first_idx_struct -> pnext = NULL;
	tempp -> coor = (int*)malloc(d * sizeof(int));
	for(i = 0; i < d; i++)
	    (tempp -> coor)[i] = grid_coor[i];
    }
    else
    {
	//printf("Same slot\n");
	INDEX* indexp = tempp -> first_idx_struct;
	while(indexp -> pnext != NULL)
	    indexp = indexp -> pnext;
	indexp -> pnext = (INDEX*)malloc(sizeof(INDEX));
	indexp = indexp -> pnext;
	indexp -> idx = idx;
	indexp -> pnext = NULL;
    }
}

int slot_not_used(GRIDTABLE* gt, int hash_val)
{
    return((gt -> gridarray)[hash_val].first_idx_struct == NULL);
}

int same_coordinate(int* coor1, int* coor2, int d)
{
    int i;
    for(i = 0; i < d; i++)
	if(coor1[i] != coor2[i])
	    return 0;
    return 1;
}

int* get_indices(GRIDTABLE* gt, int* grid_coor, int* num_indices)
{
    int d = gt -> d;
    int num_slot = gt -> num_slot;
    int start_hash_val = coordinate_hash(grid_coor, d, num_slot);
    int hash_val = start_hash_val;
    while(1)
    {
	if(slot_not_used(gt, hash_val))
	{
	    *num_indices = 0;
	    return NULL;
	}
	else if(!same_coordinate(grid_coor, (gt -> gridarray)[hash_val].coor, d))
	{
	    hash_val = (hash_val + 1) % num_slot;
	    if(hash_val == start_hash_val)
	    {
		*num_indices = 0;
		return NULL;
	    }
	}
	else
	    break;
    }
    int* indices_array;
    *num_indices = 0;
    INDEX* p = (gt -> gridarray)[hash_val].first_idx_struct;
    while(p != NULL)
    {
	(*num_indices) += 1;
	p = p -> pnext;
    }
    indices_array = (int*)malloc((*num_indices) * sizeof(int));
    int i;
    p = (gt -> gridarray)[hash_val].first_idx_struct;
    for(i = 0; i < (*num_indices); i++)
    {
	indices_array[i] = p -> idx;
	p = p -> pnext;
    }
    return indices_array;
}

int* get_nearby_indices(GRIDTABLE* gt, double* x, int* num_indices, int grid_radius)
{
    int i;
    if(grid_radius < 0)
    {
	printf("The grid radius cannot be negative.\n");
	exit(1);
    }
    int* center_grid_coor = (int*)malloc(gt -> d * sizeof(int));
    calc_grid_coordinate(gt, x, center_grid_coor);
    int* p = get_indices(gt, center_grid_coor, num_indices);
    if(grid_radius == 0)
    {
	free(center_grid_coor);
	return p;
    }
    else
    {
	int* working_grid_coor = (int*)malloc(gt -> d * sizeof(int));
	for(i = 0; i < gt -> d; i++)
	    working_grid_coor[i] = center_grid_coor[i] - grid_radius;
	int num;
	for(num = 0; num < (int)pow(2 * grid_radius + 1, gt -> d); num++)
	{
	    convert_num_to_coor(num, 2 * grid_radius + 1, gt -> d, working_grid_coor);
	    for(i = 0; i < gt -> d; i++)
		working_grid_coor[i] += (center_grid_coor[i] - grid_radius);
	    if(is_in_grid(gt, working_grid_coor) && !same_coordinate(working_grid_coor, center_grid_coor, gt -> d))
	    {
		int add_num;
		int* add_p = get_indices(gt, working_grid_coor, &add_num);
		p = append_lists(num_indices, p, *num_indices, add_p, add_num);
		free(add_p);
	    }
	}
	free(center_grid_coor);
	free(working_grid_coor);
	return p;
    }
}

int is_in_grid(GRIDTABLE* gt, int* coor)
{
    int i;
    for(i = 0; i < gt -> d; i++)
    {
	if(coor[i] >= (gt -> num_grids_per_axis) || coor[i] < 0)
	    return 0;
    }
    return 1;
}

int* append_lists(int* total_length, int* list1, int length1, int* list2, int length2)
{
    int i;
    *total_length = length1 + length2;
    list1 = (int*)realloc(list1, (*total_length) * sizeof(int));
    for(i = length1; i < (*total_length); i++)
	list1[i] = list2[i - length1];
    return list1;
}


void print_list(int* list, int length)
{
    int i;
    for(i = 0; i < length; i++)
    {
	printf("%d", list[i]);
	if(i != length - 1)
	    putchar(' ');
	else
	    putchar('\n');
    }
}


void convert_num_to_coor(int num, int n, int d, int* coor)
{
    int i;
    for(i = 0; i < d; i++)
    {
	coor[i] = num / ((int)pow(n, (d - 1 - i)));
	num = num % ((int)pow(n, (d - 1 - i)));
    }
}
