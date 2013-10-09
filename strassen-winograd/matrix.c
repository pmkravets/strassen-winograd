//
//  matrix.c
//  strassen-winograd
//
//  Created by Pavel Kravets on 03.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"

void zero_matrix(unsigned int rows, unsigned int cols, matrix* m)
{
    m->matrix = (double*) malloc(rows * cols * sizeof(double));
    m->rowNum = rows;
    m->colNum = cols;
    m->offset = 0;
    m->strand = 0;
    m->skip = 0;
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            set_element(m, i, j, .0);
        }
    }
}

void random_matrix(unsigned int rows, unsigned int cols, matrix* m)
{
    srand((unsigned int) time(NULL));
    m->matrix = (double*) malloc(rows * cols * sizeof(double));
    m->rowNum = rows;
    m->colNum = cols;
    m->offset = 0;
    m->strand = 0;
    m->skip = 0;
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            set_element(m, i, j, (double) rand() / RAND_MAX);
        }
    }
}

void eye(unsigned int size, matrix* m)
{
    m->matrix = (double*) malloc(size * size * sizeof(double));
    m->rowNum = size;
    m->colNum = size;
    m->offset = 0;
    m->strand = 0;
    m->skip = 0;
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            if (i==j)
            {
                set_element(m, i, j, 1.);
            }
            else
            {
                set_element(m, i, j, .0);
            }
        }
    }
}

void print_matrix(matrix* m)
{
    int rows = m->rowNum;
    int cols = m->colNum;
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            printf("%.3f ", element(m, i, j));
        }
        printf("\n");
    }

}

void add_matrices(matrix* m1, matrix* m2, matrix* res)
{
    if (m1->rowNum != m2->rowNum || m1->colNum != m2->colNum)
    {
        printf("Matrices can not be added");
        return;
    }
    if (!res->matrix) {
        res->matrix = (double*) malloc(m1->rowNum * m1->colNum * sizeof(double));
        res->rowNum = m1->rowNum;
        res->colNum = m2->colNum;
        res->offset = 0;
        res->strand = 0;
        res->skip = 0;
    }
    for (int i=0; i<m1->rowNum; i++)
    {
        for (int j=0; j<m1->colNum; j++)
        {
            set_element(res, i, j, element(m1, i, j) + element(m2, i, j));
        }
    }
}

void subtract_matrices(matrix* m1, matrix* m2, matrix* res)
{
    if (m1->rowNum != m2->rowNum || m1->colNum != m2->colNum)
    {
        printf("Matrices can not be subtracted");
        return;
    }
    if (!res->matrix) {
        res->matrix = (double*) malloc(m1->rowNum * m1->colNum * sizeof(double));
        res->rowNum = m1->rowNum;
        res->colNum = m2->colNum;
        res->offset = 0;
        res->strand = 0;
        res->skip = 0;
    }
    for (int i=0; i<m1->rowNum; i++)
    {
        for (int j=0; j<m1->colNum; j++)
        {
            set_element(res, i, j, element(m1, i, j) - element(m2, i, j));
        }
    }
}

double element(matrix* m, unsigned int i, unsigned int j)
{
    return m->matrix[m->offset + i*(m->colNum + m->skip) + j];
}

void set_element(matrix* m, unsigned int i, unsigned int j, double element)
{
    m->matrix[m->offset + i*(m->colNum + m->skip) + j] = element;
}

void add_to_element(matrix* m, unsigned int i, unsigned int j, double element)
{
    m->matrix[m->offset + i*(m->colNum + m->skip) + j] += element;
}

void submatrix(matrix* m, unsigned int rowNum, unsigned int colNum, unsigned int offset, unsigned int strand, unsigned int skip, matrix* res)
{
    res->matrix = m->matrix;
    res->rowNum = rowNum;
    res->colNum = colNum;
    res->offset = offset + m->offset;
    res->strand = strand;
    res->skip = skip + m->skip;
}

void copy_matrix(matrix* src, matrix* dest)
{
    unsigned int row, col;
    row = src->rowNum < dest->rowNum ? src->rowNum : dest->rowNum;
    col = src->colNum < dest->colNum ? src->colNum : dest->colNum;
    for (int i = 0; i<row; i++)
    {
        for (int j = 0; j<col; j++)
        {
            set_element(dest, i, j, element(src, i, j));
        }
    }
    
}

matrix* alloc_matrix()
{
    matrix* m = malloc(sizeof(matrix));
    m->matrix = NULL;
    m->offset = 0;
    m->skip = 0;
    m->strand = 0;
    return m;
}

void free_matrix(matrix* m)
{
    free(m->matrix);
    free(m);
}


