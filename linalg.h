#ifndef __LINALG_H
#define __LINALG_H

typedef struct
{
    int rows;
    int cols;
    double** data;
} Matrix;

/* memory management */
Matrix* new_matrix(int, int);
void del_matrix(Matrix*);

/* misc */
void print_matrix(Matrix*);

/* unary operations */
void rref(Matrix*);
void transpose(Matrix*);
Matrix* subset(Matrix*, int, int, int, int);
Matrix* inverse(Matrix*);

/* binary operations */
Matrix* add(Matrix*, Matrix*);
Matrix* subtract(Matrix*, Matrix*);
Matrix* multiply(Matrix*, Matrix*);
Matrix* row_bind(Matrix*, Matrix*);
Matrix* col_bind(Matrix*, Matrix*);

#endif
