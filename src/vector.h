#ifdef VECTOR_H
#else
#define VECTOR_H
typedef struct __vector {
  double x;
  double y;
  double z;
} vector;

typedef struct __symm{
    vector rot[3];
    vector trans;
} symm;

typedef struct __vec_llist{
    vector val;
    struct __vec_llist * next;
} vec_llist;

void init_vector(vector * v, double x, double y, double z);
int equal(vector v1, vector v2);
int equale(vector v1, vector v2, double epsdiff);
double distance(vector v1, vector v2, vector * Tmat);
vector cross_product(vector v1, vector v2);
double dot_product(vector v1, vector v2);
double volume_product(vector v1, vector v2, vector v3);
vector vector_scale(double a, vector v);
vector vector_multiply(vector v1, int * n);
vector vector_add(vector v1, vector v2);
vector vector_sub(vector v1, vector v2);    //v1-v2
int translate_match(vector *rv, vector x1, vector x2);
int isenclosed(vector v1, vector v2);
int find_vector(vector v, vector * list, int nlist);
vector vector_rotate(vector in, double  rotation[3][3]);
vector vector_Rotate(vector in, vector * Rot);

double vector_norm(vector v);
vector vector_normalization(vector v);
vector vector_round(vector v);
vector array2vector(double * in);
int kpt_equivalent(vector kpt1, vector kpt2, double lattice[3][3]);

int vec_comp(vector v1, vector v2, double tolerance);

// init
void vec_llist_init(vec_llist ** p2head);

// add to end of linked list
void vec_llist_add(vec_llist ** p2head, vector val);    

// add to linked list with last->val.x <= val.x && last->val.y <= val.y && last->val.z <= val.z
int vec_llist_add_inorder(vec_llist ** p2head, vector val, int flag_force);

// delete an element and return it.
vector vec_llist_del(vec_llist ** p2head, vector val, int * nerr);

// free up whole list
void vec_llist_free(vec_llist ** p2head);

// pop an element and return it.
vector vec_llist_pop(vec_llist ** p2head, int * nerr);

// find the location in linked list, if not found return -1
int vec_llist_find(vec_llist ** p2head, vector val);

#endif
