#include "mesh.h"

struct EdgeNode {
    E_Int p, q;
    E_Int i, posi;
    mutable E_Int j, posj;

    EdgeNode(E_Int P, E_Int Q)
    {
        p = std::min(P, Q);
        q = std::max(P, Q);
        i = posi = j = posj = -1;
    }

    bool operator<(const EdgeNode &e) const
    {
        return (p < e.p) || (p == e.p && q < e.q);
    }
};

void IMesh::make_skin_graph()
{
    auto &xadj = sgraph.xadj;
    auto &fpts = sgraph.fpts;
    auto &fadj = sgraph.fadj;

    xadj.clear();
    fpts.clear();
    fadj.clear();
    
    xadj.resize(skin.size()+1);
    xadj[0] = 0;

    for (size_t i = 0; i < skin.size(); i++) {
        E_Int fid = skin[i];
        xadj[i+1] = xadj[i] + F[fid].size();
    }

    fpts.resize(xadj[skin.size()]);
    fadj.resize(xadj[skin.size()], -1);

    for (size_t i = 0; i < skin.size(); i++) {
        E_Int fid = skin[i];
        const auto &pn = F[fid];
        E_Int start = xadj[i];
        E_Int end = xadj[i+1];
        for (E_Int j = 0; j < end-start; j++)
            fpts[start+j] = pn[j];
    }

    std::set<EdgeNode> edge_set;

    for (size_t i = 0; i < skin.size(); i++) {
        E_Int start = xadj[i];
        E_Int np = xadj[i+1] - start;
        const E_Int *pn = &fpts[start];
        for (E_Int j = 0; j < np; j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%np];
            EdgeNode e(p, q);
            auto it = edge_set.find(e);
            if (it == edge_set.end()) {
                e.i = i;
                e.posi = j;
                edge_set.insert(e);
            } else {
                assert(it->i != -1);
                assert(it->posi != -1);
                assert(it->j == -1);
                assert(it->posj == -1);
                it->j = i;
                it->posj = j;
            }
        }
    }

    for (const auto &e : edge_set) {
        E_Int pi = xadj[e.i] + e.posi;
        assert(fadj[pi] == -1);
        fadj[pi] = e.j;

        E_Int pj = xadj[e.j] + e.posj;
        assert(fadj[pj] == -1);
        fadj[pj] = e.i;
    }
}

