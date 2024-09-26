#pragma once

#include "xcore.h"

struct ArrayI {
    E_Int count;
    E_Int *ptr;
};

void ArrayI_free(ArrayI *arr);