#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "PairHashTable.h"


PHASH* create_empty_hash(int l)
{
    PHASH* ph = (PHASH*)malloc(sizeof(PHASH));
    ph -> length = l;
    ph -> num_used = 0;
    ph -> p = (HELEM*) malloc(l * sizeof(HELEM));
    int i;
    for(i = 0; i < l; i++)
    {
	((ph -> p) + i) -> key.fr = -1;
	((ph -> p) + i) -> key.to = -1;
	((ph -> p) + i) -> val = 0.0;
    }
    return ph;
}

void free_hash(PHASH* ph)
{
    free(ph -> p);
    free(ph);
}

int is_null_entry(TPAIR tp)
{
    return(tp.fr == -1 && tp.to == -1);
}

int pair_equal(TPAIR tp1, TPAIR tp2)
{
    return (tp1.fr == tp2.fr && tp1.to == tp2.to);
}

long hash_fun(TPAIR tp)
{
    return ((fabs(tp.fr) + 1) ) * ((fabs(tp.to) + 1) );
}

//Return the idx of the element, -1 if not exist
int exist_in_hash(PHASH* ph, TPAIR tp)
{
    int start_idx = hash_fun(tp) % (ph -> length);
    int idx = start_idx;
    TPAIR temp_pair;
    while(1)
    {
	temp_pair = (ph -> p)[idx].key;
	if(is_null_entry(temp_pair))
	    return -1;
	else if(pair_equal(temp_pair, tp))
	    return idx;
	else
	{
	    idx = (idx + 1) % (ph -> length);
	    if(idx == start_idx)
		return -1;
	}
    }

}

void add_entry(PHASH* ph, HELEM elem)
{
    if(ph -> num_used == ph -> length)
    {
	printf("The hash table is already full.\n");
	exit(1);
    }
    int start_idx = hash_fun(elem.key) % (ph -> length);
    int idx = start_idx; 
    while(!is_null_entry((ph -> p)[idx].key))
	idx = (idx + 1) % (ph -> length);
    (ph -> p)[idx] = elem;
    (ph -> num_used) += 1;

}

double retrieve_value(PHASH* ph, TPAIR tp)
{
    int start_idx = hash_fun(tp) % (ph -> length);
    int idx = start_idx;
    while(!pair_equal((ph -> p)[idx].key, tp))
	idx = (idx + 1) % (ph -> length);
    return (ph -> p)[idx].val;
}

void update_with(PHASH* ph, int idx, double add_value)
{
    ((ph -> p) + idx) -> val += add_value;
}

double retrieve_value_with_idx(PHASH* ph, int idx)
{
    return ((ph -> p) + idx) -> val;
}
