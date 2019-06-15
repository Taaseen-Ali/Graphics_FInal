#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "uthash.h"

struct my_struct{
  double * key;
  double * normals;
  UT_hash_handle hh;
};

struct my_struct *vnormals;


double* search(double key0, double key1, double key2);
void insert(double p0, double p1, double p2, double n0, double n1, double n2);
void print_hash();
