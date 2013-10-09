//
//  strassen.c
//  strassen-winograd
//
//  Created by Pavel Kravets on 06.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "strassen.h"
#include "naive.h"


int next_power(int n)
{
    n--;
    n |= n >> 1;   // Divide by 2^k for consecutive doublings of k up to 32,
    n |= n >> 2;   // and then or the results.
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;           // The result is a number of 1 bits equal to the number
    // of bits in the original number, plus 1. That's the
    // next highest power of 2.
    return n;
}

void pad_matrix(matrix* old, matrix* new)
{
    int maxDim = old->rowNum > old->colNum ? old->colNum : old->colNum;
    int newDim = next_power(maxDim);
    zero_matrix(newDim, newDim, new);
    for (int i=0; i<old->rowNum; i++)
    {
        for (int j=0; j<old->colNum; j++)
        {
            set_element(new, i, j, element(old, i, j));
        }
    }
}

void strip_zeros(matrix* m, matrix* res, unsigned int row, unsigned int col)
{
    zero_matrix(row, col, res);
    copy_matrix(m, res);
    
}

void partition_matrix(matrix* m , matrix* a11, matrix* a12, matrix* a21, matrix* a22)
{
    int cols = m->colNum / 2;
    int rows = m->rowNum / 2;
    int strand, skip;
    strand = skip = cols;
    submatrix(m, rows, cols, 0, strand, skip, a11);
    submatrix(m, rows, cols, cols, strand, skip, a12);
    submatrix(m, rows, cols, cols*m->colNum + cols*m->skip, strand, skip, a21);
    submatrix(m, rows, cols, cols*m->colNum+cols + cols*m->skip, strand, skip, a22);
}

void strassen_helper(matrix* m1, matrix* m2, matrix* res, unsigned int depth, unsigned int maxDepth, unsigned int naiveSize)
{
    if (depth>=maxDepth || m1->rowNum==naiveSize || m1->rowNum==2)
    {
        naive_mult(m1, m2, res);
        return;
    }
    
    matrix* a11 = malloc(sizeof(matrix));
    matrix* a12 = malloc(sizeof(matrix));
    matrix* a21 = malloc(sizeof(matrix));
    matrix* a22 = malloc(sizeof(matrix));
    
    partition_matrix(m1, a11, a12, a21, a22);
    
    matrix* b11 = malloc(sizeof(matrix));
    matrix* b12 = malloc(sizeof(matrix));
    matrix* b21 = malloc(sizeof(matrix));
    matrix* b22 = malloc(sizeof(matrix));
    
    partition_matrix(m2, b11, b12, b21, b22);
    
    zero_matrix(m1->rowNum, m1->colNum, res);
    
    matrix* res11 = malloc(sizeof(matrix));
    matrix* res12 = malloc(sizeof(matrix));
    matrix* res21 = malloc(sizeof(matrix));
    matrix* res22 = malloc(sizeof(matrix));
    
    partition_matrix(res, res11, res12, res21, res22);
    
    matrix* temp1 = alloc_matrix();
    add_matrices(a11, a22, temp1);
    
    matrix* temp2 = alloc_matrix();
    add_matrices(b11, b22, temp2);
    
    matrix* p1 = alloc_matrix();
    strassen_helper(temp1, temp2, p1, depth+1, maxDepth, naiveSize);
    
    add_matrices(a21, a22, temp1);
    
    matrix* p2 = alloc_matrix();
    strassen_helper(temp1, b11, p2, depth+1, maxDepth, naiveSize);
    
    subtract_matrices(b12, b22, temp1);
    
    matrix* p3 = alloc_matrix();
    strassen_helper(a11, temp1, p3, depth+1, maxDepth, naiveSize);
    
    subtract_matrices(b21, b11, temp1);
    
    matrix* p4 = alloc_matrix();
    strassen_helper(a22, temp1, p4, depth+1, maxDepth, naiveSize);
    
    add_matrices(a11, a12, temp1);
    
    matrix* p5 = alloc_matrix();
    strassen_helper(temp1, b22, p5, depth+1, maxDepth, naiveSize);
    
    subtract_matrices(a21, a11, temp1);
    add_matrices(b11, b12, temp2);
    
    matrix* p6 = alloc_matrix();
    strassen_helper(temp1, temp2, p6, depth+1, maxDepth, naiveSize);
    
    subtract_matrices(a12, a22, temp1);
    add_matrices(b21, b22, temp2);
    
    matrix* p7 = alloc_matrix();
    strassen_helper(temp1, temp2, p7, depth+1, maxDepth, naiveSize);
    
    add_matrices(p1, p4, res11);
    subtract_matrices(res11, p5, res11);
    add_matrices(res11, p7, res11);
    
    add_matrices(p3, p5, res12);
    
    add_matrices(p2, p4, res21);
    
    subtract_matrices(p1, p2, res22);
    add_matrices(res22, p3, res22);
    add_matrices(res22, p6, res22);

    free_matrix(temp1);
    free_matrix(temp2);
    free_matrix(p1);
    free_matrix(p2);
    free_matrix(p3);
    free_matrix(p4);
    free_matrix(p5);
    free_matrix(p6);
    free_matrix(p7); 
}

void strassen_mult(matrix* orig1, matrix* orig2, matrix* res, unsigned int maxDepth, unsigned int naiveSize)
{
    if (orig1->colNum != orig2->rowNum)
    {
        printf("Error. Matrices can not be multiplicated");
        return;
    }
    if (orig1->rowNum == orig1->colNum && orig1->rowNum==orig2->colNum && next_power(orig1->rowNum)==orig1->rowNum)
    {
        strassen_helper(orig1, orig2, res, 0, maxDepth, naiveSize);
    }
    else
    {
        matrix* m1 = malloc(sizeof(matrix));
        matrix* m2  = malloc(sizeof(matrix));
        pad_matrix(orig1, m1);
        pad_matrix(orig2, m2);
        matrix* r = alloc_matrix();
        strassen_helper(m1, m2, r, 0, maxDepth, naiveSize);
        strip_zeros(r, res, orig1->rowNum, orig2->colNum);
    }
    
}
