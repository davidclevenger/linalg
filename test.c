#include <stdio.h>
#include "linalg.h"

int main()
{
    Matrix* matrix = new_matrix(2, 1);
    Matrix* b = new_matrix(2, 2);

    matrix->data[0][0] = 1.0;
    matrix->data[1][0] = 1.0;

    b->data[0][0] = 2.0;
    b->data[0][1] = 2.0;
    b->data[1][0] = 2.0;
    b->data[1][1] = 2.0;

    print_matrix(b);
    Matrix* x = multiply(b, matrix);
    print_matrix(x);

    del_matrix(matrix);
    del_matrix(b);
    del_matrix(x);
}
