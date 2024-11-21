#pragma once

#include "xcore.h"

struct Box3;
struct FaceSort;
struct ArrayI;

#define MAX_FACES_PER_POINT 8

struct Point {
    E_Float x, y, z;
};

struct PointFaces {
    E_Int count;
    E_Int ptr[MAX_FACES_PER_POINT];
};

bool Point_in_tri(E_Float px, E_Float py, E_Float pz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz);

bool Point_in_Box3D(const Point *p, const Box3 *box);

bool Point_in_FaceSort(const Point *p, const FaceSort *f);

void PointFaces_extract_by_threshold
(
    const PointFaces *sploc, E_Int spcount,
    const E_Int *skin, E_Int mcount,
    const E_Int threshold,
    ArrayI *faces
);

void points_write(const char *fname, const std::vector<Point> &P);

void point_write(E_Float px, E_Float py, E_Float pz);

void point_write(const Point *p);