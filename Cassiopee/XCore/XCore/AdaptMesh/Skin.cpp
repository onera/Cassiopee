#include "Skin.h"
#include "common/mem.h"
#include "Mesh.h"
#include "DynMesh.h"

#include <stack>
#include <set>

void SkinGraph_free(SkinGraph *skin_graph)
{
    skin_graph->nf = 0;
    XFREE(skin_graph->skin);
    XFREE(skin_graph->xadj);
    XFREE(skin_graph->fpts);
    XFREE(skin_graph->fnei);
}

struct EdgeNode {
    E_Int m_p, m_q;
    E_Int m_i, m_posi;
    mutable E_Int m_j, m_posj;

    EdgeNode(E_Int p, E_Int q)
    {
        m_p = std::min(p, q);
        m_q = std::max(p, q);
        m_i = -1;
        m_posi = -1;
        m_j = -1;
        m_posj = -1;
    }

    bool operator<(const EdgeNode &e) const
    {
        return (m_p < e.m_p) || (m_p == e.m_p && m_q < e.m_q);
    }
};

void SkinGraph_make_skin_neighbours(SkinGraph *skin_graph)
{
    E_Int nf = skin_graph->nf;
    const E_Int *xadj = skin_graph->xadj;
    const E_Int *fpts = skin_graph->fpts;

    skin_graph->fnei = (E_Int *)XMALLOC(xadj[nf] * sizeof(E_Int));
    memset(skin_graph->fnei, -1, xadj[nf] * sizeof(E_Int));
    E_Int *fnei = skin_graph->fnei; 

    std::set<EdgeNode> edge_set;

    for (E_Int i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int np = xadj[i+1] - start;
        const E_Int *pn = &fpts[start];
        for (E_Int j = 0; j < np; j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%np];
            EdgeNode e(p, q);
            auto it = edge_set.find(e);
            if (it != edge_set.end()) {
                assert(it->m_i    != -1);
                assert(it->m_posi != -1);
                assert(it->m_j    == -1);
                assert(it->m_posj == -1);
                it->m_j = i; 
                it->m_posj = j;
            } else {
                e.m_i = i;
                e.m_posi = j;
                edge_set.insert(e);
            }
        } 
    }

    for (const auto &e : edge_set) {
        E_Int pi = xadj[e.m_i] + e.m_posi;
        assert(fnei[pi] == -1);
        fnei[pi] = e.m_j;

        E_Int pj = xadj[e.m_j] + e.m_posj;
        assert(fnei[pj] == -1);
        fnei[pj] = e.m_i;
    }
}

void SkinGraph_smooth_ref_data(const SkinGraph *skin_graph, E_Int *fdat,
    const Mesh *M)
{
    E_Int nf = skin_graph->nf;

    std::stack<E_Int> stk;

    for (E_Int i = 0; i < nf; i++) {
        if (fdat[i] > 0)
            stk.push(i);
    }

    const E_Int *xadj = skin_graph->xadj;
    const E_Int *fnei = skin_graph->fnei;

    while (!stk.empty()) {

        E_Int fid = stk.top();
        stk.pop();

        E_Int start = xadj[fid];
        E_Int nneis = xadj[fid+1] - start;
        const E_Int *neis = &fnei[start];

        for (E_Int i = 0; i < nneis; i++) {
            E_Int nei = neis[i];
            E_Int incr_nei = fdat[nei] + M->flevel[skin_graph->skin[nei]];
            E_Int incr_fid = fdat[fid] + M->flevel[skin_graph->skin[fid]];
            E_Int diff = abs(incr_nei - incr_fid);
            if (diff <= 1) continue;
            E_Int idx_to_modify = incr_fid > incr_nei ? nei : fid;
            fdat[idx_to_modify] += diff-1;
            stk.push(idx_to_modify);
        }
    }
}
