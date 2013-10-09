//
//  strassen.h
//  strassen-winograd
//
//  Created by Pavel Kravets on 06.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#ifndef strassen_winograd_strassen_h
#define strassen_winograd_strassen_h

#include "matrix.h"

void strassen_mult(matrix* m1, matrix* m2, matrix* res, unsigned int maxDepth, unsigned int naiveSize);

#endif
