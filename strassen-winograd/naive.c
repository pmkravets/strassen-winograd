//
//  naive.c
//  strassen-winograd
//
//  Created by Pavel Kravets on 07.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "naive.h"

void naive_mult(matrix* m1, matrix* m2, matrix* res)
{
    if (m1->colNum != m2->rowNum)
    {
        printf("Error. Matrices can not be multiplicated");
        return;
    }
    
    int rows = res->rowNum = m1->rowNum;
    int cols = res->colNum = m2->colNum;
    
    res->matrix = (double*) malloc(rows * cols * sizeof(double));
    
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            double sum = 0.;
            for (int k=0; k<m1->colNum; k++)
            {
                sum += element(m1, i, k) * element(m2, k, j);
            }
            set_element(res, i, j, sum);
        }
    }
}
