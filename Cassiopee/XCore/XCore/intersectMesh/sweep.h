#pragma once

#include <vector>

#include "queue.h"
#include "status.h"
#include "dcel.h"
#include "segment.h"

void sweep(Queue &Q, Status &T, std::vector<Segment *> &S,
    std::vector<Vertex *> &I, std::vector<Hedge *> &H);
