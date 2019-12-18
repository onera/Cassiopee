#ifndef __CEDRE_IO_MPI_IO_NODE_RANK_H__
#define __CEDRE_IO_MPI_IO_NODE_RANK_H__

#include "pdm.h"
#include "pdm_mpi.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Définition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Retourne le rang du procesus courant dans le noeud
 *
 * parameters :
 *   comm            <-- Communicateur MPI
 *
 * return :
 *   Rank
 *
 *----------------------------------------------------------------------------*/

int
PDM_io_mpi_node_rank
(
PDM_MPI_Comm comm
);

/*----------------------------------------------------------------------------
 * Retourne le nom du noeud
 *
 * parameters :
 *   hostname_ptr       <-> Pointeur sur la chaine contenant le nom
 *   hostname_length    <-> Longueur de la chaine
 *
 * return :
 *   Code d'erreur
 *
 *----------------------------------------------------------------------------*/

void
PDM_io_get_hostname
(
char  **hostname_ptr,
size_t *hostname_length
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CEDRE_IO_MPI_IO_NODE_RANK_H__ */
