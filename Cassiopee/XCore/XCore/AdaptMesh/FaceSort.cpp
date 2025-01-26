#include "FaceSort.h"
#include "Mesh.h"

void FaceSort_compute_data(const Mesh *M, FaceSort *mfaces, E_Int mcount)
{
    for (E_Int i = 0; i < mcount; i++) {
        FaceSort *face = &mfaces[i];

        E_Int tid = face->fid;
        assert(M->ftype[tid] == TRI);

        E_Int *tri = Mesh_get_face(M, tid);
        E_Int A = tri[0], B = tri[2], C = tri[4];

        face->UX = (M->X[B] - M->X[A]);
        face->UY = (M->Y[B] - M->Y[A]);
        face->UZ = (M->Z[B] - M->Z[A]);

        face->VX = (M->X[C] - M->X[A]);
        face->VY = (M->Y[C] - M->Y[A]);
        face->VZ = (M->Z[C] - M->Z[A]);

        face->UU = face->UX*face->UX + face->UY*face->UY + face->UZ*face->UZ;
        face->VV = face->VX*face->VX + face->VY*face->VY + face->VZ*face->VZ;
        face->UV = face->UX*face->VX + face->UY*face->VY + face->UZ*face->VZ;

        face->inv_denom = face->UU*face->VV - face->UV*face->UV;

        assert(face->inv_denom != 0.0);

        face->inv_denom = 1.0 / face->inv_denom;

        face->fc[0] = (M->X[A] + M->X[B] + M->X[C]) / 3.0;
        face->fc[1] = (M->Y[A] + M->Y[B] + M->Y[C]) / 3.0;
        face->fc[2] = (M->Z[A] + M->Z[B] + M->Z[C]) / 3.0;

        // Store A
        face->xa = M->X[A];
        face->ya = M->Y[A];
        face->za = M->Z[A];
    }
}