#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "vector.h"
#include "matrix.h"
#include "wannorb.h"
#include "wanndata.h"
#include "rotate_ham.h"

void init_vector(vector * v, double x, double y, double z) {
    v->x=x;
    v->y=y;
    v->z=z;
}

int translate_match(vector * rv, vector x1, vector x2) {
  int ii;
  vector tmp;
  rv->x=(int)(rint(x1.x-x2.x));
  rv->y=(int)(rint(x1.y-x2.y));
  rv->z=(int)(rint(x1.z-x2.z));
  tmp=vector_add(x2, (*rv));
  if (distance(x1, tmp, NULL)<eps5)
    return 1;
  else
    return 0;
}

vector vector_scale(double a, vector v) {
  vector r;
  r.x=a*v.x;
  r.y=a*v.y;
  r.z=a*v.z;

  return r;
}

vector vector_multiply(vector v1, int * n) {
  vector r;
  r.x=v1.x*n[0];
  r.y=v1.y*n[1];
  r.z=v1.z*n[2];

  return r;
}

vector vector_sub(vector v1, vector v2) {
  vector vr;
  vr.x=v1.x-v2.x;
  vr.y=v1.y-v2.y;
  vr.z=v1.z-v2.z;
  return vr;
}

vector vector_add(vector v1, vector v2) {
  vector vr;
  vr.x=v1.x+v2.x;
  vr.y=v1.y+v2.y;
  vr.z=v1.z+v2.z;
  return vr;
}

int equal(vector v1, vector v2) {
  if ((fabs(v1.x-v2.x)<eps5) &&
      (fabs(v1.y-v2.y)<eps5) &&
      (fabs(v1.z-v2.z)<eps5))
    return 1;
  else
    return 0;
}

int equale(vector v1, vector v2, double epsdiff) {
  if ((fabs(v1.x-v2.x)<epsdiff) &&
      (fabs(v1.y-v2.y)<epsdiff) &&
      (fabs(v1.z-v2.z)<epsdiff))
    return 1;
  else
    return 0;
}

vector matrix_dot(vector * Tmat, vector v) {    //=v*T
  vector r;
  r.x=v.x*Tmat[0].x+v.y*Tmat[1].x+v.z*Tmat[2].x;
  r.y=v.x*Tmat[0].y+v.y*Tmat[1].y+v.z*Tmat[2].y;
  r.z=v.x*Tmat[0].z+v.y*Tmat[1].z+v.z*Tmat[2].z;
  return r;
}

double distance(vector v1, vector v2, vector * Tmat) {
  double r;
  vector x1, x2;
  if (Tmat!=NULL) {
    x1=matrix_dot(Tmat, v1);
    x2=matrix_dot(Tmat, v2);
  }
  else {
    x1=v1;
    x2=v2;
  }

  r=sqrt((x1.x-x2.x)*(x1.x-x2.x)+(x1.y-x2.y)*(x1.y-x2.y)+(x1.z-x2.z)*(x1.z-x2.z));
  return r;
}

vector cross_product(vector v1, vector v2) {
  vector r;
  r.x=v1.y*v2.z-v1.z*v2.y;
  r.y=v1.z*v2.x-v1.x*v2.z;
  r.z=v1.x*v2.y-v1.y*v2.x;

  return r;
}

double dot_product(vector v1, vector v2) {
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

double volume_product(vector v1, vector v2, vector v3) {
  double vol;
  vol=v1.x*v2.y*v3.z+v1.y*v2.z*v3.x+v1.z*v2.x*v3.y-
     (v1.x*v2.z*v3.y+v1.y*v2.x*v3.z+v1.z*v2.y*v3.x);
  return vol;
}

int isenclosed(vector v1, vector v2) {
  return ( fabs(v1.x)<=v2.x &&
           fabs(v1.y)<=v2.y &&
           fabs(v1.z)<=v2.z );
}

int find_vector(vector v, vector * list, int nlist) {
  int ii;
  for (ii=0; ii<nlist; ii++) {
    if (equal(v, list[ii]))
      return ii;
  }
  return -1;
}

vector vector_rotate(vector in, double  rotation[3][3]) {
  vector out;
  vector symm[3];
  int i;

  for(i=0;i<3;i++){
      init_vector(symm+i, rotation[i][0], rotation[i][1], rotation[i][2]);
  }

  out.x=dot_product(symm[0], in);
  out.y=dot_product(symm[1], in);
  out.z=dot_product(symm[2], in);
  
  //out = vector_scale(1.0/dot_product(out,out), out); //normalization
  return out;
}

vector vector_Rotate(vector in, vector * Rot){
  vector out;
  out.x=dot_product(Rot[0], in);
  out.y=dot_product(Rot[1], in);
  out.z=dot_product(Rot[2], in);
  return out;
}


double vector_norm(vector v)
{
  return sqrt(dot_product(v,v));
}

vector vector_normalization(vector v)
{
    vector out;
    out = vector_scale(1.0/vector_norm(v), v);
    return out;
}

vector vector_round(vector v){
  vector out;
  out.x = round(v.x);
  out.y = round(v.y);
  out.z = round(v.z);
  return out;
}

vector array2vector(double * in){
    vector v;
    v.x=in[0];
    v.y=in[1];
    v.z=in[2];
    return v;

}

int kpt_equivalent(vector kpt1, vector kpt2, double lattice[3][3]){
    vector dk;
    double rlatt[3][3];
    matrix3x3_inverse(rlatt, lattice);
    dk = vector_sub(kpt1,kpt2);
    //dk = vector_rotate(dk , rlatt);
    if( fabs(dk.x - round(dk.x)) > eps5) return 0;
    if( fabs(dk.y - round(dk.y)) > eps5) return 0;
    if( fabs(dk.z - round(dk.z)) > eps5) return 0;
    return 1;
}


int vec_comp(vector v1, vector v2, double tolerance){
    if(fabs(v1.x - v2.x) > tolerance) 
        return (int) ((v1.x - v2.x)/tolerance);
    else if (fabs(v1.y - v2.y) > tolerance)
        return (int) ((v1.y - v2.y)/tolerance);
    else if (fabs(v1.z - v2.z) > tolerance)
        return (int) ((v1.z - v2.z)/tolerance);
    else
        return 0;
}

void vec_llist_init(vec_llist ** p2head){
    * p2head = NULL;
}

void vec_llist_add(vec_llist ** p2head, vector val){
    // add to end of linked list
    vec_llist * curr = *p2head;
    vec_llist * new_node;

    new_node       = (vec_llist *) malloc(sizeof(vec_llist));
    new_node->val  = val;
    new_node->next = NULL;
    if(*p2head == NULL){
        *p2head = new_node;
        return;
    }
    else{
        while(curr->next != NULL){
            curr = curr->next;
        }
        curr->next = new_node;
    }
}

int vec_llist_add_inorder(vec_llist ** p2head, vector val, int flag_force){
    // add to list with last->val < add->val < next->val
    // if flag_force != 1 and if val is already in the llist, return -1 and do not add it.
    vec_llist * new_node;
    vec_llist * last = NULL;
    vec_llist * curr = *p2head;
    while(curr != NULL && vec_comp(curr->val, val, eps6) < 0 ){
        last = curr;
        curr = curr->next;
    }
    if(flag_force == 1 || curr == NULL || vec_comp(curr->val, val, eps6) != 0){
        new_node       = (vec_llist *) malloc(sizeof(vec_llist));
        new_node->val  = val;
        new_node->next = curr;
        if(last == NULL){
            *p2head = new_node;
        } else {
            last->next = new_node;
        }
        return 1;
    }
    else {
        return -1;
    }
}

vector vec_llist_del(vec_llist ** p2head, vector val, int * nerr){
    // delete an element and return it.
    vector ret;
    vec_llist * last = NULL;
    vec_llist * next = NULL;
    vec_llist * curr = *p2head;
    *nerr = 1; // nerr==1 for no error
    init_vector(&ret, 0.5, 0.5, 0.5);
    while(curr != NULL && vec_comp(curr->val, val, eps6) != 0){
        last = curr;
        curr = curr->next;
    }
    if(*p2head == NULL){
        *nerr = -1;
        return ret;
    } else if(curr == NULL){
        *nerr = -1;
        return ret;
    } else if( last == NULL){
        ret  = curr->val;
        next = curr->next;
        free(curr);
        *p2head = next;
        return ret;
    } else {
        ret  = curr->val;
        next = curr->next;
        free(curr);
        last->next = next;
        return ret;
    }
    *nerr = 0;
    return ret;
} 

void vec_llist_free(vec_llist ** p2head){
    vec_llist * curr = *p2head;
    vec_llist * next;
    while(curr != NULL){
        next = curr->next;
        free(curr);
        curr = next;
    }
}

vector vec_llist_pop(vec_llist ** p2head, int * nerr){
    // pop an element and return it.
    // head = * p2head
    vector ret;
    vec_llist * next=NULL;
    *nerr = 1;   //no err

    init_vector(&ret, 0.5, 0.5, 0.5);
    if(*p2head != NULL){
        next = (*p2head)->next;
        ret  = (*p2head)->val; 
        free(*p2head);
        *p2head = next;
    }
    else{
        *nerr = -1;  // llist is empty.
    }
    return ret;
}

int vec_llist_find(vec_llist ** p2head, vector val){
    vec_llist * curr = *p2head;
    int ival = 0;

    while(curr != NULL && vec_comp(curr->val, val, eps6) != 0){
        ival++;
        curr = curr->next;
    }
    if(curr == NULL)
        ival = -1;
    return ival;
}
