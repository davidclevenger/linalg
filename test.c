#include <stdio.h>
#include "linalg.h"

int main()
{
    Matrix* b = new_matrix(2, 2);

    b->data[0][0] = 3.0;
    b->data[0][1] = 2.0;
    b->data[1][0] = 4.0;
    b->data[1][1] = 1.0;

    Matrix* inv = inverse(b);
	print_matrix(inv);
	del_matrix(inv);
    del_matrix(b);
}
