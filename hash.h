#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "uthash.h"

struct my_struct{
  double * ar;
  int id;
  UT_hash_handle hh;
};

struct my_struct *vnormals;


double* search(double key0, double key1, double key2);
void insert(double key0, double key1, double key2, double data0, double data1, double data2);
void print_hash();
