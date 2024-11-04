#include "dcel.h"
#include "smesh.h"

void Dcel::reconstruct_smesh(Smesh &Mf, int color, bool check_Euler) const
{
    std::map<E_Int, std::vector<Face *>> oid_to_faces;

    std::map<Vertex *, E_Int> new_pids;

    E_Int NP = Mf.np;

    E_Int face_incr = 0;

    // Which faces to keep?
    for (Face *f : F) {
        if (f->oids[color] == -1) continue;

        oid_to_faces[f->oids[color]].push_back(f);

        // Which points to keep?

        auto vertices = get_face_vertices(f);

        for (Vertex *v : vertices) {
            auto it = new_pids.find(v);
            if (it == new_pids.end()) {
                // Was this point before the intersection?
                if (v->oids[color] != -1) {
                    // Yes, keep its original id
                    new_pids[v] = v->oids[color];
                } else {
                    // No, add it
                    new_pids[v] = NP++;
                }
            }
        }

        face_incr++;
    }

    face_incr -= oid_to_faces.size();

    // Replace the old faces with the new intersection faces

    Mf.F.resize(Mf.nf+face_incr);
    Mf.Fc.resize(Mf.nf+face_incr);
    E_Int NF = Mf.nf;
    for (const auto &fdat : oid_to_faces) {
        E_Int fid = fdat.first;
        const auto &faces = fdat.second;

        std::vector<E_Int> new_pn;

        // Replace fid with the first face and add the rest of the faces
        for (size_t i = 0; i < faces.size(); i++) {
            new_pn.clear();
            auto vertices = get_face_vertices(faces[i]);
            for (Vertex *v : vertices) new_pn.push_back(new_pids[v]);
            if (i > 0) {
                Mf.F[NF] = new_pn;
                Mf.Fc[NF] = new_pn;
                NF++;
            } else {
                Mf.F[fid] = new_pn;
                Mf.Fc[fid] = new_pn;
            }
        }
    }


    // Resize the points
    Mf.X.resize(NP, EFLOATMAX);
    Mf.Y.resize(NP, EFLOATMAX);
    Mf.Z.resize(NP, EFLOATMAX);

    for (const auto &vdat : new_pids) {
        Vertex *v = vdat.first;
        E_Int new_pid = vdat.second;
        if (v->oids[color] == -1) {
            assert(Mf.X[new_pid] == EFLOATMAX);
            assert(Mf.Y[new_pid] == EFLOATMAX);
            assert(Mf.Z[new_pid] == EFLOATMAX);
            Mf.X[new_pid] = v->x;
            Mf.Y[new_pid] = v->y;
            Mf.Z[new_pid] = v->z;
        }
    }


    Mf.nf = NF;
    Mf.np = NP;
}


Smesh Dcel::export_smesh(bool check_Euler) const
{
    Smesh smesh;
    smesh.check_Euler = check_Euler;

    smesh.np = V.size();
    smesh.X.resize(smesh.np);
    smesh.Y.resize(smesh.np);
    smesh.Z.resize(smesh.np);

    for (Vertex *v : V) {
        smesh.X[v->id] = v->x;
        smesh.Y[v->id] = v->y;
        smesh.Z[v->id] = v->z;
    }

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        Hedge *REP = f->rep;
        Cycle *c = REP->cycle;
        if (c->inout != Cycle::OUTER) continue;
        Hedge *h = REP;

        std::vector<E_Int> pn;

        do {
            Vertex *v = h->orig;
            pn.push_back(v->id);
            h = h->next;
        } while (h != REP);

        smesh.F.push_back(pn);
    }

    smesh.nf = (E_Int)smesh.F.size();

    smesh.Fc = smesh.F;

    smesh.make_edges();
    
    return smesh;
}

std::vector<Dcel::Cycle *> Dcel::extract_cycles_of_indices(
    const std::vector<E_Int> &indices) const
{
    std::vector<Cycle *> ret;
    ret.reserve(indices.size());

    for (E_Int index : indices) {
        ret.push_back(C[index]);
    }

    return ret;
}

std::vector<E_Int> Dcel::extract_indices_of_type(int type) const
{
    std::vector<E_Int> ret;

    for (size_t i = 0; i < C.size(); i++) {
        if (C[i]->inout == type) ret.push_back(i);
    }

    return ret;
}
