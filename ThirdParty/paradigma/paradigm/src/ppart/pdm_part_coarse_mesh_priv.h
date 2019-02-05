#ifndef PDM_PART_COARSE_MESH_PRIV_H
#define	PDM_PART_COARSE_MESH_PRIV_H

/* 
 * File:   pdm_part_coarse_mesh_priv.h
 * Author: jmagnene
 *
 * Created on July 12, 2016, 10:54 AM
 */

#include "pdm_part_priv.h"
#include "pdm_timer.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_mpi.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 * 
 * _part_t defines a coarse partition 
 *
 */

typedef struct  {
  
  _part_t      *part;        //Coarse mesh
  
  int           nCoarseCellWanted;    /*!< Number Cell wanted for agglomeration     */ 
  
  // int *cellWeight;    /*!< Integer weight for graoh partitionning  */ 
  // int *faceWeight;    /*!< Number Cell wanted for agglomeration     */ 
  
  int *coarseCellCellIdx;    //Array of indexes of the connected partitions (size : nCoarseCell + 1)
  
  int *coarseCellCell;       //Partitioning array (size : coarseCellCellIdx[nCoarseCell]) 
  
  int *coarseFaceGroupToFineFaceGroup; //Coarse face group - fine face group connectivity (size = faceGroupIdx[nFaceGroup])
  
  int *coarseFaceToFineFace; //Coarse face - fine face connectivity (size = nCoarseFace)

  int *coarseVtxToFineVtx;   //Coarse vertex - fine vertex connectivity (size = nCoarseVtx)

  void *specific_data;   /*!< Specific data  */

  /* Array specific to anisotropic agglomeration */
  /* int *agglomerationLines;       */
  /* int *agglomerationLinesIdx;    */
  /* int  agglomerationLinesIdx_size;   */
  /* int *isOnFineBnd;   */
  
  /* /\* Array specific to anisotropic agglomeration if Initialise from a finer grid *\/ */
  /* int *agglomerationLinesInit;       */
  /* int *agglomerationLinesInitIdx;    */
  /* int  agglomerationLinesInitIdx_size;   */
  /* int *isOnFineBndInit;   */
  
  
} _coarse_part_t;


/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 * 
 * _coarse_mesh_t defines a coarse mesh 
 *
 */

typedef struct  {

  /* Partitions */          
    
  int nPart;        /*!< Number of partitions to define
                                      on this process */ 
  
  int method;       /*!< Partitioning method */
  int nTPart;       /*!< Total number of partitions */
  int nFaceGroup;   /*!< Number of boundaries */
  
  int have_cellTag;
  int have_faceTag;
  int have_vtxTag;
  int have_cellWeight;
  int have_faceWeight;
  int have_faceGroup;

  void *specific_data;
  
  //TIMER
  
  PDM_timer_t *timer;             /*!< Timer */ 

  double times_elapsed [18];          /*!< Elapsed times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu[18];             /*!< CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_u[18];           /*!< User CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_s[18];          /*!< Systeme CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  
  /* Communicator */
  
  PDM_MPI_Comm  comm;   /*!< Communicator */
    
  _part_t        **part_ini;               /*!< Partition: fine mesh                            */
  
  _coarse_part_t **part_res;               //Coarse mesh
  


} _coarse_mesh_t;
    

/**
 * \struct PDM_part_renum_fct_t
 *
 * \brief  Function pointer used to define a coarse mesh method
 *
 */

typedef void (*PDM_coarse_mesh_fct_t) (_coarse_mesh_t  *cm,
                                       const int       ipart,
                                       int             *nCoarseCellComputed,
                                       int             *cellCellIdx,
                                       int             *cellCell,
                                       int             *cellPart);

/**
 * \struct _coarse_mesh_method_t
 * \brief coarse mesh method
 *
 */

typedef struct _renum_method_t {

  char                  *name; /*!< Name of method */
  PDM_coarse_mesh_fct_t   fct;  /*!< Renumbering function */

} _coarse_mesh_method_t;


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Return an initialized coarse part object
 *
 */

static inline _coarse_part_t * 
_coarse_part_create
(
void
 )
{
  _coarse_part_t *cp = (_coarse_part_t *) malloc(sizeof(_coarse_part_t));
  cp->part = _part_create();
  
  cp->coarseCellCell = NULL;
  
  cp->coarseCellCellIdx = NULL;
  
  cp->coarseFaceGroupToFineFaceGroup = NULL;
  
  cp->coarseFaceToFineFace = NULL; 

  cp->coarseVtxToFineVtx = NULL;

  cp->specific_data = NULL;
  
  return cp;
  
}


/*============================================================================
 * Public function prototypes
 *============================================================================*/



/**
 *
 * \brief Add a new coarse mesh method
 *
 * \param [in]      name          Mesh entity to renumber
 * \param [in]      fct           Function
 *
 */

int
PDM_coarse_mesh_method_add
(
 const char                 *name,     /*!< Name          */
 PDM_coarse_mesh_fct_t       fct       /*!< Function      */
 );

  
/**
 *
 * \brief Get index of a coarse mesh method from it's name
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_idx_get_cf, PDM_COARSE_MESH_METHOD_IDX_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 );

int
PDM_coarse_mesh_method_idx_get
(
const char *name
 );


/**
 *
 * \brief Get name of a coarse mesh method from it's index
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_name_get_cf, PDM_COARSE_MESH_METHOD_NAME_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 );

char *
PDM_coarse_mesh_method_name_get
(
const int id
 );
  

/**
 *
 * \brief Get the number of coarse mesh method
 *
 * \return Number of methods
 *
 */

void
PROCF (pdm_coarse_mesh_method_n_get, PDM_COARSE_MESH_METHOD_N_GET)
(
 int  *n_method
 );
  
int
PDM_coarse_mesh_method_n_get
(
void
);
  

/**
 *
 * \brief Purge coarse mesh methods catalog
 *
 */

void
PDM_coarse_mesh_method_purge
(
void
 );


/**
 *
 * \brief Load local coarse mesh methods
 *
 */

void
PDM_coarse_mesh_method_load_local
(
void
 );


/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

_coarse_mesh_t *
PDM_part_coarse_mesh_get_from_id
(
 int  cmId
 );


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_COARSE_MESH_PRIV_H */
