/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/* #if !defined(_XOPEN_SOURCE) || !defined(_BSD_SOURCE) */
/* #define _XOPEN_SOURCE 500 */
/* #endif */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <stddef.h>
#include <stdint.h>

#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Définitions des macro locales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Default hostname size
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Max
 *----------------------------------------------------------------------------*/

#define _CEDRE_IO_MAX(a,b)        (a) > (b) ? (a) : (b)

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*============================================================================
 * Variables globales
 *============================================================================*/

static const int  local_hostname_default_len = 64;

/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction de hachage adler32
 *
 * parameters :
 *   buf              <-- Buffer
 *   buflength        <-- Taille du buffer
 * return
 *   Cle de hachage
 *----------------------------------------------------------------------------*/

static uint32_t
 _adler32
(
const void *buf,
size_t buflength
)
{

  const uint8_t * buffer = (const uint8_t *)buf;

  uint32_t s1 = 1;
  uint32_t s2 = 0;

  for (size_t n = 0; n < buflength; n++) {
    s1 = (s1 + buffer[n]) % 65521;
    s2 = (s2 + s1) % 65521;
  }

  return (s2 << 16) | s1;
}


/*----------------------------------------------------------------------------
 * Rang du processus courant sur son noeud
 *
 * parameters :
 *   comm             <-- Communicateur MPI global
 * return
 *   Rank
 *----------------------------------------------------------------------------*/

static int
_node_rank_by_hash
(
PDM_MPI_Comm comm
)
{

  char * hostname = NULL;
  size_t hostnameLength = 0;

// TODO: Use MPI_Get_processor_name instead of PDM_io_get_hostname
//  int MPI_Get_processor_name(char *name, int *resultlen)
  PDM_io_get_hostname(&hostname, &hostnameLength);

  /* Détermine la clé de hachage */

  uint32_t checkSum = _adler32(hostname, hostnameLength);

  if (0 == 1)
    PDM_printf( "-- hostname adlerCode %s %ud\n", hostname, checkSum);

  int commRank = -1;
  PDM_MPI_Comm_rank(comm, &commRank);

  PDM_MPI_Comm nodeComm = PDM_MPI_COMM_NULL;

  /* Détermine le communicateur local au noeud */

  PDM_MPI_Comm_split(comm, checkSum, commRank, &nodeComm);

  /* Détermine le rang dans communicateur du noeud */

  int nodeRank;
  PDM_MPI_Comm_rank(nodeComm, &nodeRank);

  int nodeSize;
  PDM_MPI_Comm_size(nodeComm, &nodeSize);

  /* Vérifie si 2 noms de noeud ne rentrent pas en collision :
     Si c'est le cas, les processus de ces 2 noeuds sont
     regroupés dans un même communicateur */

  int nSend = (int) hostnameLength;

  PDM_MPI_Allreduce((void *) &hostnameLength, (void *) &nSend, 1,
                PDM_MPI_INT, PDM_MPI_MAX, nodeComm);

  nSend += 1; /* \0 */

  char *send = (char *)malloc(sizeof(char) * nSend);

  for (int i = 0; i < nSend; i++)
    send[i] = 0x00;

  strncpy(send, hostname, nSend);

  char * recv = (char *)malloc(sizeof(char) * nSend * nodeSize);
  PDM_MPI_Allgather((void * )send, nSend, PDM_MPI_CHAR, (void *)recv, nSend, PDM_MPI_CHAR, nodeComm);

  char *neighbor = recv;
  int localNodeRank = 0;

  for (int i = 0; i < nodeSize; ++i) {
    if (strcmp(send, neighbor) == 0) {
      if (i < nodeRank) {
        ++localNodeRank;
      }
      else {
        break;
      }
    }
    neighbor += nSend;
  }

  int isSameNode = (nodeRank != localNodeRank);
  int isSameNodeMax;

  PDM_MPI_Allreduce((void *) &isSameNode, (void *) &isSameNodeMax, 1,
                PDM_MPI_INT, PDM_MPI_MAX, nodeComm);

  if (isSameNodeMax == 1) {

    /* Traitement des collisions */

    neighbor = recv;
    char *hostInComm = (char *) malloc(sizeof(char) * nSend * nodeSize);
    strncpy(hostInComm, recv, nSend);
    int nHostInComm = 1;

    for (int j = 1; j < nodeSize; j++) {
      char *ptHostInComm = hostInComm;
      for (int i = 0; i < nHostInComm; i++) {
        if (strcmp(ptHostInComm, neighbor) != 0) {
          strncpy(hostInComm + nHostInComm * nSend, neighbor, nSend);
          nHostInComm += 1;
          break;
        }
        ptHostInComm += nSend;
      }
      neighbor += nSend;
    }

    int tag = 0;
    for (int i = 0; i < nHostInComm; i++) {
      char *ptHostInComm = hostInComm;
      if (strcmp(ptHostInComm, send) == 0) {
        tag = i;
        break;
      }
      ptHostInComm += nSend;
    }

    PDM_MPI_Comm nodeComm2;
    PDM_MPI_Comm_split(nodeComm, tag, nodeRank, &nodeComm2);

    int nodeRank2;
    PDM_MPI_Comm_rank(nodeComm2, &nodeRank2);

    nodeRank = nodeRank2;

    PDM_MPI_Comm_free(&nodeComm2);
    free(hostInComm);

  }

  // Clean up.

  free(send);
  send = NULL;

  free(recv);
  recv = NULL;

  PDM_MPI_Comm_free(&nodeComm);

  free(hostname);
  hostname = NULL;

  // Affichage

  if (0 == 1) {

    int masterNodeRank = 0;
    int nbIORank = 0;
    if (nodeRank == 0)
      masterNodeRank = 1;


    PDM_MPI_Allreduce((void *) &masterNodeRank, (void *) &nbIORank, 1,
                  PDM_MPI_INT, PDM_MPI_SUM, comm);

    int commSize;
    PDM_MPI_Comm_size(comm, &commSize);
    if (commRank == 0)
      PDM_printf("Part of ranks for parallel IO : %d/%d\n", nbIORank, commSize);

  }

  return nodeRank;
}

/*============================================================================
 * Definition des fonctions publiques
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
)
{

  int nodeRank = -1;

  nodeRank = _node_rank_by_hash(comm);

  char * hostname;
  size_t hostnameLength;

  PDM_io_get_hostname(&hostname, &hostnameLength);
  free(hostname);
  hostname = NULL;

  return nodeRank;
}

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
)
{

  char * hostname = NULL;
  size_t nHostname = 0;

  int  allocateMore = 0;

  *hostname_ptr = NULL;
  *hostname_length = 0;

  do {
    nHostname += local_hostname_default_len;

    hostname = (char *)malloc(sizeof(char) * nHostname);

    if (hostname == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "allocating %lu bytes of memory for hostname failed: %s.\n",
              sizeof(char) * nHostname, strerror(errno));
      abort();
    }

    int error;

    error = gethostname(hostname, nHostname);

    if (error == -1) {
      if (errno == ENAMETOOLONG) {
        allocateMore = 1;
        free(hostname); hostname = NULL;
      }
      else {
        free(hostname);
        hostname = NULL;

        PDM_error(__FILE__, __LINE__, 0, "gethostname failed with error %d: %s.\n", errno, strerror(errno));
        abort();
      }

    }
    else {
      allocateMore = 0;
    }

  } while (allocateMore == 1);

  hostname[nHostname - 1] = 0x00;

  *hostname_length = strlen(hostname);
  *hostname_ptr = hostname;
}

#undef _CEDRE_IO_MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
