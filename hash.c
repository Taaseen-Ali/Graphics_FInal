#include "hash.h" 

struct my_struct * vnormals = NULL;

double* search(double key0, double key1, double key2) {
  double * search_key = malloc(sizeof(double)*3);
  search_key[0] = key0;
  search_key[1] = key1;
  search_key[2] = key2;

  struct my_struct * result;
  HASH_FIND(hh, vnormals, search_key, sizeof(double)*3, result);
  return result->normals;
}

void insert(double p0, double p1, double p2, double n0, double n1, double n2){
  struct my_struct * point;
  point = malloc(sizeof(struct my_struct));

  point->key = malloc(sizeof(double)*3);
  point->key[0]=p0;
  point->key[1]=p1;
  point->key[2]=p2;

  point->normals = malloc(sizeof(double)*3);
  point->normals[0]=n0;
  point->normals[1]=n1;
  point->normals[2]=n2;

  HASH_ADD_KEYPTR(hh, vnormals, point->key,
		  sizeof(double)*3, point);
}
