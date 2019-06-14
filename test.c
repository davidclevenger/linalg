#include <stdio.h>
#include "linalg.h"

int main()
{
    Matrix* b = new_matrix(2, 2);

    b->data[0][0] = 1.0;
    b->data[0][1] = 2.0;
    b->data[1][0] = 4.0;
    b->data[1][1] = 2.0;

    print_matrix(b);
    inverse(b);

    del_matrix(b);
}
