#include <stdlib.h>
#include <stdio.h>
#include "linalg.h"

/* private prototypes */

void __value_swap(double*, double*);
void __row_swap(Matrix*, int, int);
void __col_swap(Matrix*, int, int);

void __error_dimension(void);
void __error_parameter(const char*);

/* private definitions */

void
__value_swap(double* lhs, double* rhs)
{
    double tmp = *lhs;
    *lhs = *rhs;
    *rhs = tmp;
}

void
__row_swap(Matrix* matrix, int r1_index, int r2_index)
{   
    int i;
    double save;

    if( r1_index >= matrix->rows
        || r2_index >= matrix->rows )
    {
        __error_parameter("index out of bounds");
        return; // return error_parameters
    }

    /* swap values one by one */
    for(i = 0; i < matrix->cols; i++)
    {
        save = matrix->data[r1_index][i];
        matrix->data[r1_index][i] = matrix->data[r2_index][i];
        matrix->data[r2_index][i] = save;
    }
}

void
__col_swap(Matrix* matrix, int c1_index, int c2_index)
{
    int i;
    double save;

    if( c1_index >= matrix->cols
        || c2_index >= matrix->cols )
    {
        __error_parameter("__col_swap: index out of bounds\n");
        return;
    }

    /* swap values one by one */
    for(i = 0; i < matrix->rows; i++)
    {
        save = matrix->data[i][c1_index];
        matrix->data[i][c1_index] = matrix->data[i][c2_index];
        matrix->data[i][c2_index] = save;
    }
}

void
__error_dimension(void)
{
    fprintf(stderr, "incorrect dimensions.\n"); 
}

void
__error_parameter(const char* msg)
{
    if( msg )
    {
        fprintf(stderr, "%s\n", msg);
    }
    else
    {
        fprintf(stderr, "parameter error.\n");
    }
}

/* public definitions */

Matrix*
new_matrix(int rows, int cols)
{
    int i;
    Matrix* rv = (Matrix*) calloc(sizeof(Matrix), 1);
    rv->rows = rows;
    rv->cols = cols;
    rv->data = (double**) calloc(sizeof(double*) * rv->rows, 1);
    for(i = 0; i < rv->rows; i++)
    {
        rv->data[i] = (double*) calloc(sizeof(double) * rv->cols, 1);
    }
    
    return rv; 
}

void
del_matrix(Matrix* matrix)
{
    int i;
    for(i = 0; i < matrix->rows; i++)
    {
        free(matrix->data[i]);
    }

    if( matrix->data )
    {
        free(matrix->data);
    }

    if( matrix )
    {
        free(matrix);
    }
}

void
print_matrix(Matrix* matrix)
{
    int i, j;

    for(i = 0; i < matrix->rows; i++)
    {
        if( i == 0 )
        {
            printf("[");
        }
        else
        {
            printf(" ");
        }

        /* ROW begin */
        printf("[");
        for(j = 0; j < matrix->cols; j++)
        {   
			if( matrix->data[i][j] > 0 )
			{
				printf(" ");
			}

            if( j < matrix->cols-1 )
            {
                printf("%.4f ", matrix->data[i][j]);
            }
            else
            {
                printf("%.4f", matrix->data[i][j]);
            }
        }
        printf("]");
        /* ROW end */

        if( i < matrix->rows-1 )
        {
            printf("\n");
        }
        else
        {
            printf("]\n");
        }
    }
}

void
rref(Matrix* matrix)
{
	int i, j, r, brk;
	double f, t;

	if( matrix->rows <= matrix->cols )
	{
		brk = matrix->rows;
	}
	else
	{
		brk = matrix->cols;
	}
	

	for(r = 0; r < brk; r++)
	{
		f = matrix->data[r][r];
	
		for(j = 0; j < matrix->cols; j++)
		{
			matrix->data[r][j] /= f;
		}

		for(i = 0; i < matrix->rows; i++)
		{
			if( i == r ) continue;

			t = matrix->data[i][r] / matrix->data[r][r];

			for(j = 0; j < matrix->cols; j++)
			{
				matrix->data[i][j] -= t * matrix->data[r][j];
			}
		}
	}
}

void
transpose(Matrix* matrix)
{
    int i, j;

    Matrix* transpose = new_matrix(matrix->cols, matrix->rows);

    for(i = 0; i < matrix->rows; i++)
    {
        for(j = 0; j < matrix->cols; j++)
        {
            transpose->data[j][i] = matrix->data[i][j];
        }
    }
    
    del_matrix(matrix);

    //TODO this might leak
    matrix = transpose;
}

Matrix*
subset(Matrix* matrix,
       int row_start_index,
       int col_start_index,
       int num_rows,
       int num_cols)
{
	Matrix* out;
	int i, j;

	if( row_start_index < 0
		|| col_start_index < 0
		|| row_start_index >= matrix->rows
		|| col_start_index >= matrix->cols
		|| num_rows == 0
		|| num_cols == 0
		|| row_start_index + num_rows > matrix->rows
		|| col_start_index + num_cols > matrix->cols )
	{
		__error_parameter("subset: incompatible dimensions.");
		return NULL;
	}

	out = new_matrix(num_rows, num_cols);

	for(i = row_start_index; i < row_start_index + num_rows; i++)
	{
		for(j = col_start_index; j < col_start_index + num_cols; j++)
		{
			out->data[i-row_start_index][j-col_start_index] = matrix->data[i][j];
		}
	}

	return out;
}

Matrix* inverse(Matrix* matrix)
{
    int i, j;
    /* check if square */
    if( matrix->rows != matrix->cols )
    {
        __error_dimension();
        return NULL;
    }

    Matrix* inverted = new_matrix(matrix->rows, matrix->cols);
    
    /* build identity matrix */
    for(i = 0; i < inverted->rows; i++)
    {
        for(j = 0; j < inverted->cols; j++)
        {
            if( i == j )
            {
                inverted->data[i][j] = 1.0;
            }
            else
            {
                inverted->data[i][j] = 0.0;
            }
        }
    }

    Matrix* joined = col_bind(matrix, inverted);

	/*
	* [A|I(n)] -> [I(n)|A^-1]
	*
	*/
    rref(joined);
    Matrix* out = subset(joined,
                         0,
                         inverted->cols,
                         inverted->rows,
                         inverted->cols);

    /* free intermediate matrices */
    del_matrix(inverted);
    del_matrix(joined);

    return out;
}

Matrix*
add(Matrix* lhs, Matrix* rhs)
{
    Matrix* rv;
    int i, j;

    if( lhs->rows != rhs->rows
        || lhs->cols != rhs->cols )
    {
        __error_dimension();
    }

    rv = new_matrix(lhs->rows, lhs->cols);

    for(i = 0; i < rv->rows; i++)
    {
        for(j = 0; j < rv->cols; j++)
        {
            rv->data[i][j] = lhs->data[i][j] + rhs->data[i][j];
        }
    }
    
    return rv;
}

Matrix*
subtract(Matrix* lhs, Matrix* rhs)
{
    Matrix* rv;
    int i, j;

    if( lhs->rows != rhs->rows
        || lhs->cols != rhs->cols )
    {
        __error_dimension();
    }

    rv = new_matrix(lhs->rows, lhs->cols);

    for(i = 0; i < rv->rows; i++)
    {
        for(j = 0; j < rv->cols; j++)
        {
            rv->data[i][j] = lhs->data[i][j] - rhs->data[i][j];
        }
    }
    
    return rv;
}

Matrix*
multiply(Matrix* lhs, Matrix* rhs)
{
    Matrix* rv;
    int i, j, x;

    if( lhs->cols != rhs->rows )
    {
       __error_dimension();
    }

    rv = new_matrix(lhs->rows, rhs->cols);

    for(i = 0; i < rv->rows; i++)
    {
        for(j = 0; j < rv->cols; j++)
        {
            rv->data[i][j] = 0;

            for(x = 0; x < lhs->cols; x++)
            {
                rv->data[i][j] += lhs->data[i][x] * rhs->data[x][j];
            }
        }
    }
    
    return rv;
}

Matrix*
row_bind(Matrix* top, Matrix* bottom)
{   
    int i, j, bound;

    if( top->cols != bottom->cols )
    {
        __error_dimension();
        return NULL;
    }

    Matrix* out = new_matrix(top->rows + bottom->rows, top->cols);

    for(i = 0; i < top->rows; i++)
    {
        for(j = 0; j < top->cols; j++)
        {
            out->data[i][j] = top->data[i][j];
        }
    }

    // TODO check boundaries
	bound = i + bottom->rows;
    for(; i < bound; i++)
    {
        for(j = 0; j < bottom->cols; j++)
        {
            out->data[i][j] = bottom->data[i][j];
        }
    }

    return out;
}

Matrix*
col_bind(Matrix* left, Matrix* right)
{
    int i, j, bound;

    if( left->cols != right->cols )
    {
        __error_dimension();
        return NULL;
    }

    Matrix* out = new_matrix(left->rows, left->cols + right->cols);

    for(i = 0; i < left->rows; i++)
    {
        for(j = 0; j < left->cols; j++)
        {
            out->data[i][j] = left->data[i][j];
        }
    }

    // TODO check boundaries

	bound = left->cols + right->cols;
    for(i = 0; i < right->rows; i++)
    {
        for(j = left->cols; j < bound; j++)
        {
            out->data[i][j] = right->data[i][j-(right->cols)];
        }
    }
    
    return out;
}

