#include "dcel.h"
#include "smesh.h"

/*
void Dcel::reconstruct(Smesh &Mf, int color, bool check_Euler) const
{
    Smesh ret;
    ret.check_Euler = check_Euler;
    auto &new_F = ret.F;
    std::map<Vertex *, E_Int> new_pids;
    ret.np = ret.nf = 0;

    std::vector<Face *> fids;

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        Hedge *h = f->rep;
        if (h->color != color) continue;
        Hedge *w = h->next;
        while (w != h) {
            if (w->color != h->color) break;
            w = w->next;
        }
        if (w == h) {
            fids.push_back(f);
        }
    }

    char fname[128] = {0};
    sprintf(fname, "single_color_%d.im", color);
    write_ngon(fname, fids);

    for (Face *f : F) {
        if (f->oids[color] == -1) continue;

        std::vector<E_Int> pn;

        std::vector<Vertex *> vertices = get_face_vertices(f);
        for (Vertex *v : vertices) {
            auto it = new_pids.find(v);
            if (it == new_pids.end()) {
                new_pids[v] = ret.np;
                pn.push_back(ret.np);
                ret.np++;
            } else {
                pn.push_back(it->second);
            }
        }

        new_F.push_back(pn);
        ret.nf++;
    }

    auto &new_X = ret.X;
    auto &new_Y = ret.Y;
    auto &new_Z = ret.Z;

    new_X.resize(ret.np), new_Y.resize(ret.np), new_Z.resize(ret.np);
    for (const auto &vdat : new_pids) {
        new_X[vdat.second] = vdat.first->x;
        new_Y[vdat.second] = vdat.first->y;
        new_Z[vdat.second] = vdat.first->z;
    }

    ret.Fc = ret.F;

    ret.make_edges();
}
*/

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
