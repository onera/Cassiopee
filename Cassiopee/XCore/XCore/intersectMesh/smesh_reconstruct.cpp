#include "smesh.h"
#include "mesh.h"
#include "io.h"

void Smesh::get_leaves(E_Int fid, std::vector<E_Int> &leaves) const
{
    leaves.push_back(fid);
    auto it = fchildren.find(fid);
    if (it == fchildren.end()) return;

    const auto &children_lists = it->second;
    for (const auto &children_list : children_lists) {
        for (E_Int child : children_list)
            get_leaves(child, leaves);
    }
}

void Smesh::reconstruct(IMesh &M)
{
    // POINTS

    E_Int NP = M.np + (np - np_before_adapt);
    
    M.X.resize(NP);
    M.Y.resize(NP);
    M.Z.resize(NP);
    
    NP = M.np;

    std::vector<Point> points;

    for (E_Int pid = np_before_adapt; pid < np; pid++) {
        assert(l2gp.find(pid) == l2gp.end());
        l2gp[pid] = NP;
        g2lp[NP] = pid;

        M.X[NP] = X[pid];
        M.Y[NP] = Y[pid];
        M.Z[NP] = Z[pid];

        points.push_back({M.X[NP], M.Y[NP], M.Z[NP]});

        NP++;
    }

    // Isolate the smesh cells

    std::set<E_Int> conf_cells;

    for (E_Int fid = 0; fid < nf_before_adapt; fid++) {
        E_Int own = M.owner[l2gf.at(fid)];
        assert(M.neigh[l2gf.at(fid)] == -1);
        conf_cells.insert(own);
    }

    // Isolate the faces to conformize
    std::set<E_Int> conf_faces;

    for (E_Int cid : conf_cells) {
        const auto &pf = M.C[cid];
        for (E_Int fid : pf) {
            if (g2lf.find(fid) == g2lf.end()) conf_faces.insert(fid);
        }
    }

    // Delete from owners the smesh faces

    for (E_Int cid : conf_cells) {
        const auto &pf = M.C[cid];
        std::vector<E_Int> PF;
        for (E_Int fid : pf) {
            if (g2lf.find(fid) == g2lf.end()) PF.push_back(fid);
        }
        M.C[cid] = PF;
    }

    // FACES

    E_Int NF = M.nf + (nf - nf_before_adapt);

    M.F.resize(NF);
    M.owner.resize(NF);
    M.neigh.resize(NF, -1);

    NF = M.nf;

    for (E_Int fid = nf_before_adapt; fid < nf; fid++) {
        assert(l2gf.find(fid) == l2gf.end());
        l2gf[fid] = NF;
        assert(g2lf.find(NF) == g2lf.end());
        g2lf[NF] = fid;
        NF++;
    }

    // Update old faces and add new faces

    for (E_Int fid = 0; fid < nf; fid++) {
        std::vector<E_Int> PN(Fc[fid]);
        for (auto &pid : PN) pid = l2gp.at(pid);

        E_Int gfid = l2gf.at(fid);
        assert(M.neigh[gfid] == -1);

        M.F[gfid] = PN;
    }


    for (E_Int fid : conf_faces) {
        auto it = g2lf.find(fid);
        assert (it == g2lf.end());

        const auto &pn = M.F[fid];

        std::vector<E_Int> new_pn;

        for (size_t i = 0; i < pn.size(); i++) {
            E_Int gp = pn[i];
            E_Int gq = pn[(i+1)%pn.size()];

            auto itp = g2lp.find(gp);
            auto itq = g2lp.find(gq);

            if (itp == g2lp.end() || itq == g2lp.end()) {
                new_pn.push_back(gp);
                continue;
            }

            E_Int lp = itp->second;
            E_Int lq = itq->second;

            std::vector<E_Int> local;
            get_edge_centers(lp, lq, local);
            for (auto p : local) new_pn.push_back(l2gp.at(p));
        }

        M.F[fid] = new_pn;
    }

    // Add skin faces leaves to owners

    for (E_Int fid = 0; fid < nf_before_adapt; fid++) {
        std::vector<E_Int> leaves;
        get_leaves(fid, leaves);
    
        // Parent element to update
        E_Int gfid = l2gf.at(fid);
        assert(M.owner[gfid] != -1);
        assert(M.neigh[gfid] == -1);
        E_Int own = M.owner[gfid];
        assert(conf_cells.find(own) != conf_cells.end());
        auto &pf = M.C[own];

        // Add the children
        for (E_Int leaf : leaves) {
            E_Int gleaf = l2gf.at(leaf);
            pf.push_back(gleaf);

            M.owner[gleaf] = own;
            M.neigh[gleaf] = -1;
        }
    }

    /*
    {
        if (M.nc > 24) {
            auto pf = M.C[24];
            std::map<E_Int, E_Int> local_map;
            int npts = 0;
            int i = 0;
            for (E_Int fid : pf) {
                const auto &pn = M.F[fid];
                printf("face%d: ", i);
                for (E_Int p : pn) {
                    auto it = local_map.find(p);
                    if (it == local_map.end()) {
                        local_map[p] = npts;
                        printf("%d ", npts);

                        if (npts == 10) {
                            point_write("10.im", M.X[p], M.Y[p], M.Z[p]);
                        }

                        if (npts == 11) {
                            point_write("11.im", M.X[p], M.Y[p], M.Z[p]);
                        }

                        if (npts == 12) {
                            point_write("12.im", M.X[p], M.Y[p], M.Z[p]);
                        }

                        if (npts == 13) {
                            point_write("13.im", M.X[p], M.Y[p], M.Z[p]);
                        }

                        if (npts == 14) {
                            point_write("14.im", M.X[p], M.Y[p], M.Z[p]);
                        }

                        npts++;
                    } else {
                        printf("%d ", it->second);
                    }
                }
                puts("");
                i++;
            }

        }
    }
    */

    M.np = NP;
    M.nf = NF;
}