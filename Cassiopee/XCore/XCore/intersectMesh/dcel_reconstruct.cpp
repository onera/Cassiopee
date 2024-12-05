/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "dcel.h"
#include "smesh.h"

void Dcel::reconstruct(Smesh &Mf, int color) const
{
    std::map<E_Int, std::vector<E_Int>> ofid_to_ifids;

    std::vector<E_Int> faces_to_keep;
    
    E_Int nif = 0;

    for (size_t fid = 0; fid < F.size(); fid++) {
        Face *f = F[fid];
        E_Int ofid = f->oids[color];
        if (ofid == -1) continue;

        ofid_to_ifids[ofid].push_back(fid);

        faces_to_keep.push_back(fid);

        nif++;
    }

    //write_ngon("keep.im", faces_to_keep);

    Mf.F.resize(faces_to_keep.size());

    // POINTS

    std::map<E_Int, E_Int> new_pids;
    E_Int NP = Mf.np;

    for (size_t i = 0; i < faces_to_keep.size(); i++) {
        E_Int fid = faces_to_keep[i];
        //Face *f = F[fid];
        const auto &vertices = Fv[fid];
        for (size_t j = 0; j < vertices.size(); j++) {
            Vertex *v = vertices[j];
            if (v->oids[color] != -1) continue;
            E_Int vid = v->id;
            auto it = new_pids.find(vid);
            if (it == new_pids.end()) {
                new_pids[vid] = NP;
                NP++;
            }
        }
    }

    E_Int NF = Mf.nf;

    for (const auto &fdat : ofid_to_ifids) {
        E_Int parent = fdat.first;
        const auto &children = fdat.second;

        // In Mf, replace parent with first child
        const auto &vertices = Fv[children[0]];
        std::vector<E_Int> PN(vertices.size());
        for (size_t i = 0; i < vertices.size(); i++) {
            Vertex *v = vertices[i];
            E_Int oid = v->oids[color];
            if (oid != -1) {
                PN[i] = oid;
            } else {
                PN[i] = new_pids.at(v->id);
            }
        }
        Mf.F[parent] = PN;

        // Construct fchildren for volume mesh reconstruction
        std::vector<E_Int> fchildren(children.size()-1);

        // Add the rest of the children
        for (size_t i = 1; i < children.size(); i++) {
            //Face *f = F[children[i]];
            const auto &vertices = Fv[children[i]];
            std::vector<E_Int> PN(vertices.size());
            for (size_t j = 0; j < vertices.size(); j++) {
                Vertex *v = vertices[j];
                E_Int oid = v->oids[color];
                if (oid != -1) {
                    PN[j] = oid;
                } else {
                    PN[j] = new_pids.at(v->id);
                }
            }
            Mf.F[NF] = PN;

            fchildren[i-1] = NF;

            NF++;
        }

        Mf.fchildren[parent].push_back(fchildren);
    }

    Mf.Fc = Mf.F;

    Mf.X.resize(NP);
    Mf.Y.resize(NP);
    Mf.Z.resize(NP);

    for (const auto &pdat : new_pids) {
        E_Int vid = pdat.first;
        E_Int pid = pdat.second;
        Mf.X[pid] = V[vid]->x;
        Mf.Y[pid] = V[vid]->y;
        Mf.Z[pid] = V[vid]->z;
    }

    Mf.nf = NF;
    Mf.np = NP;

    Mf.clear_conformal_data();
    Mf.make_edges();
    Mf.make_point_faces();
    Mf.make_point_edges();

    // Construct edge centers for volume mesh reconstruction

    const auto &vc = vcenter[color];

    for (const auto &vdat : vc) {
        const std::pair<Vertex *, Vertex *> &e = vdat.first;
        const Vertex *c = vdat.second;
        const Vertex *p = e.first;
        const Vertex *q = e.second;

        E_Int pid =  p->oids[color] != -1 ? p->oids[color] : new_pids.at(p->id);
        E_Int qid =  q->oids[color] != -1 ? q->oids[color] : new_pids.at(q->id);

        Mf.ecenter[{pid, qid}] = new_pids.at(c->id);
    }

    puts("constructed ecenters");
}