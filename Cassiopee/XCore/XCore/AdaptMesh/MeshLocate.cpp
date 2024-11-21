#include "Mesh.h"
#include "Point.h"

bool Mesh_point_in_tri(const Mesh *M, const Point *p, E_Int tid)
{
    const E_Int *face = Mesh_get_face(M, tid);
    E_Int A = face[0], B = face[2], C = face[4];
    return Point_in_tri(p->x, p->y, p->z,
                        M->X[A], M->Y[A], M->Z[A],
                        M->X[B], M->Y[B], M->Z[B],
                        M->X[C], M->Y[C], M->Z[C]);
}

bool Mesh_point_in_quad(const Mesh *M, const Point *p, E_Int qid)
{
    // TODO(Imad): maybe compute face centers once in pre-pass
    // Star the quad into 4 triangles
    E_Float O[3] = {0.0, 0.0, 0.0};
    const E_Int *face = Mesh_get_face(M, qid);
    E_Int A = face[0], B = face[2], C = face[4], D = face[6];
    O[0] = (M->X[A] + M->X[B] + M->X[C] + M->X[D]) * 0.25;
    O[1] = (M->Y[A] + M->Y[B] + M->Y[C] + M->Y[D]) * 0.25;
    O[2] = (M->Z[A] + M->Z[B] + M->Z[C] + M->Z[D]) * 0.25;

    bool hit = false;

    // First triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[A], M->Y[A], M->Z[A],
                       M->X[B], M->Y[B], M->Z[B]);
    if (hit) return true;

    // Second triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[B], M->Y[B], M->Z[B],
                       M->X[C], M->Y[C], M->Z[C]);
    if (hit) return true;


    // Third triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[C], M->Y[C], M->Z[C],
                       M->X[D], M->Y[D], M->Z[D]);
    if (hit) return true;

    // Fourth triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[D], M->Y[D], M->Z[D],
                       M->X[A], M->Y[A], M->Z[A]);
    if (hit) return true;

    return false;
}

bool Mesh_point_in_face(const Mesh *M, const Point *p, E_Int fid)
{
    if (M->ftype[fid] == QUAD) return Mesh_point_in_quad(M, p, fid);
    assert(M->ftype[fid] == TRI);
    return Mesh_point_in_tri(M, p, fid);
}