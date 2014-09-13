#ifndef PAIR_HASH_TABLE
#define PAIR_HASH_TABLE

//Transition pair
typedef struct
{
    int fr;
    int to;
}
TPAIR;

//element for hash table
typedef struct
{
    TPAIR key;
    double val;
}
HELEM;

typedef struct
{
    int length; //Total number of available elements
    int num_used;
    HELEM* p;
}
PHASH;

PHASH* create_empty_hash(int l);
void free_hash(PHASH* ph);
int is_null_entry(TPAIR tp);
int pair_equal(TPAIR tp1, TPAIR tp2);
long hash_fun(TPAIR tp);
int exist_in_hash(PHASH* ph, TPAIR tp);
void add_entry(PHASH* ph, HELEM elem);
double retrieve_value(PHASH* ph, TPAIR tp);
void update_with(PHASH* ph, int idx, double add_value);
double retrieve_value_with_idx(PHASH* ph, int idx);
#endif
