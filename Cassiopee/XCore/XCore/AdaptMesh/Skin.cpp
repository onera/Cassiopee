#include "Skin.h"
#include "common/mem.h"

void SkinGraph_free(SkinGraph *skin_graph)
{
    skin_graph->nf = 0;
    XFREE(skin_graph->skin);
    XFREE(skin_graph->xadj);
    XFREE(skin_graph->fpts);
    XFREE(skin_graph->fnei);
}

struct EdgeNode {
    E_Int p, q;
    E_Int i, j;
    E_Int posi, posj;
    EdgeNode *next;
};

static
EdgeNode *make_edge_node(E_Int p, E_Int q, E_Int i, E_Int posi)
{
    EdgeNode *node = (EdgeNode *)XMALLOC(sizeof(EdgeNode));
    node->p = p < q ? p : q;
    node->q = p < q ? q : p;
    node->i = i;
    node->posi = posi;
    node->j = -1;
    node->posj = -1;
    node->next = NULL;
    return node;
}

static
EdgeNode *find_edge_node(EdgeNode **ht, E_Int hsize, E_Int p, E_Int q)
{
    E_Int p_ = p < q ? p : q;
    E_Int q_ = p < q ? q : p;
    E_Int bucket = p_ % hsize;
    EdgeNode *current = ht[bucket];

    while (current) {
        E_Int P = current->p, Q = current->q;
        if (P == p_ && Q == q_) {
            return current;
        }
        current = current->next;
    }

    return NULL;
}

static
void insert_edge_node(EdgeNode *node, const E_Int hsize, EdgeNode **ht)
{
    assert(node->p < node->q);
    E_Int bucket = node->p % hsize;
    EdgeNode *current = ht[bucket];

    if (current) {
        EdgeNode *tmp = current;
        ht[bucket] = node;
        node->next = tmp;
    } else {
        ht[bucket] = node;
    }
}

void SkinGraph_make_skin_neighbours(SkinGraph *skin_graph)
{
    E_Int nf = skin_graph->nf;
    const E_Int *xadj = skin_graph->xadj;
    const E_Int *fpts = skin_graph->fpts;

    skin_graph->fnei = (E_Int *)XMALLOC(xadj[nf] * sizeof(E_Int));
    E_Int *fnei = skin_graph->fnei;
    memset(fnei, -1, xadj[nf] * sizeof(E_Int));

    EdgeNode **ht = (EdgeNode **)XMALLOC(nf * sizeof(EdgeNode *));
    memset(ht, 0, nf * sizeof(EdgeNode *));

    for (E_Int i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int np = xadj[i+1] - start;
        const E_Int *pn = &fpts[start];
        for (E_Int j = 0; j < np; j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%np];
            EdgeNode *node = find_edge_node(ht, nf, p, q);
            if (node) {
                assert(node->i    != -1);
                assert(node->posi != -1);
                assert(node->j    == -1);
                assert(node->posj == -1);
                node->j = i; 
                node->posj = j;
            } else {
                node = make_edge_node(p, q, i, j);
                insert_edge_node(node, nf, ht);
            }
        } 
    }

    for (E_Int i = 0; i < nf; i++) {
        EdgeNode *node = ht[i];
        while (node) {
            E_Int pi = xadj[node->i] + node->posi;
            assert(fnei[pi] == -1);
            fnei[pi] = node->j;
            
            E_Int pj = xadj[node->j] + node->posj;
            assert(fnei[pj] == -1);
            fnei[pj] = node->i;

            node = node->next;
        }
    }
}