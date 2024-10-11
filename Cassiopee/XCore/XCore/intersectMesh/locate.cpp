#include "mesh.h"
#include "triangle.h"
#include "primitives.h"

std::vector<PointLoc> IMesh::locate(const Smesh &Sf)
{
    std::vector<PointLoc> ploc(Sf.np);

    for (E_Int pid = 0; pid < Sf.np; pid++) {

        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];

        E_Int I = floor((x - xmin) / HX);
        E_Int J = floor((y - ymin) / HY);
        E_Int K = floor((z - zmin) / HZ);
        E_Int voxel = get_voxel(I, J, K);

        const auto &pf = bin_faces[voxel];

        bool found = false;

        // w, u, v
        auto &loc = ploc[pid];
        assert(loc.fid == -1);
        assert(loc.v_idx == -1);
        assert(loc.e_idx == -1);
        assert(loc.u == -1);
        assert(loc.v == -1);
        assert(loc.w == -1);

        for (size_t i = 0; i < pf.size() && !found; i++) {
            E_Int fid = pf[i];
            const auto &pn = F[fid];
            E_Int a = pn[0], b = pn[1], c = pn[2];
            E_Float u, v, w;
            if (Triangle::is_point_inside(
                x, y, z,
                X[a], Y[a], Z[a],
                X[b], Y[b], Z[b],
                X[c], Y[c], Z[c],
                u, v, w)) {
                
                found = true;

                loc.fid = fid;
                loc.u = u;
                loc.v = v;
                loc.w = w;

                if      (Sign(1-u, 1e-3) == 0) loc.v_idx = 0;
                else if (Sign(1-v, 1e-3) == 0) loc.v_idx = 1;
                else if (Sign(1-w, 1e-3) == 0) loc.v_idx = 2;
                else if (Sign(w, 1e-3) == 0) loc.e_idx = 0;
                else if (Sign(u, 1e-3) == 0) loc.e_idx = 1;
                else if (Sign(v, 1e-3) == 0) loc.e_idx = 2;

                break;
            }
        }

        assert(found);
    }

    return ploc;
}
