//
//  winograd.h
//  strassen-winograd
//
//  Created by Pavel Kravets on 04.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#ifndef strassen_winograd_winograd_h
#define strassen_winograd_winograd_h

#include "matrix.h"

void winograd_mult(matrix* m1, matrix* m2, matrix* res);

#endif
