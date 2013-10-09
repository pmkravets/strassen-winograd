//
//  winograd.c
//  strassen-winograd
//
//  Created by Pavel Kravets on 04.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "winograd.h"

void winograd_preprocess(matrix* m1, matrix* m2, double* row, double* col)
{
    int d = m1->colNum/2;
    
    for (int i=0; i<m1->rowNum; i++)
    {
        row[i] = 0;
        for (int j=0; j<d; j++)
        {
            row[i] += element(m1, i, 2*j) * element(m1, i, 2*j+1);
        }
    }
    
    for (int i=0; i<m2->colNum; i++)
    {
        col[i] = 0;
        for (int j=0; j<d; j++)
        {
            col[i] += element(m2, 2*j, i) * element(m2, 2*j+1, i);
        }
    }
    
}

void winograd_mult(matrix* m1, matrix* m2, matrix* res)
{
    if (m1->colNum != m2->rowNum)
    {
        printf("Error. Matrices can not be multiplicated");
        return;
    }
    res->matrix = (double*) malloc(m1->rowNum * m2->colNum * sizeof(double));
    double* row = (double*) malloc(m1->rowNum * sizeof(double));
    double* col = (double*) malloc(m2->colNum * sizeof(double));
    
    int d = m1->colNum/2;
    
    res->rowNum = m1->rowNum;
    res->colNum = m2->colNum;
    
    winograd_preprocess(m1, m2, row, col);
    
    for (int i=0; i<res->rowNum; i++)
    {
        for (int j=0; j<res->colNum; j++)
        {
            set_element(res, i, j, -row[i]-col[j]);
            for (int k=0; k<d; k++)
            {
                add_to_element(res, i, j, (element(m1, i, 2*k)+element(m2, 2*k+1, j)) * (element(m1, i, 2*k+1)+element(m2, 2*k, j)));
            }
            if (m1->colNum % 2 !=0)
            {
                add_to_element(res, i, j, element(m1, i, m1->colNum-1)*element(m2, m2->rowNum-1, j));
            }
        }
    }
    
}
