#include "Skin.h"
#include "common/mem.h"

void SkinGraph_free(SkinGraph *skin_graph)
{
    XFREE(skin_graph->skin);
    XFREE(skin_graph->xadj);
    XFREE(skin_graph->fpts);
    XFREE(skin_graph->fnei);
}
