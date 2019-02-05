#ifndef __PDM_GLOBAL_POINT_MEAN_H__
#define __PDM_GLOBAL_POINT_MEAN_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute a global mean   
 *
 * \param [in]   n_part       Number of local partitions 
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_global_mean_create
(
 const int n_part,
 const PDM_MPI_Comm comm
);

void
PROCF (pdm_global_mean_create, PDM_GLOBAL_MEAN_CREATE)
(
 const int *n_part,
 const PDM_MPI_Fint *fcomm,
       int *id
);


/**
 *
 * \brief Set absolute number   
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_point      Number of points in the partition
 * \param [in]   numabs       Absolute number of points
 *
 */

void
PDM_global_mean_set
(
 const int          id,
 const int          i_part,
 const int          n_point,
 const PDM_g_num_t *numabs
);

void
PROCF (pdm_global_mean_set, PDM_GLOBAL_MEAN_SET)
(
 const int         *id,
 const int         *i_part,
 const int         *n_point,
 const PDM_g_num_t *numabs
);


/**
 *
 * \brief Free a global point mean structure   
 *
 * \param [in]   id           Identifier
 *
 * \return     Identifier    
 */

void
PDM_global_mean_free
(
 const int          id
);

void
PROCF (pdm_global_mean_free, PDM_GLOBAL_MEAN_FREE)
(
 const int         *id
);


/**
 *
 * \brief Set local field and it associated weight    
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part                Current partition
 * \param [in]   stride                Stride of the field
 * \param [in]   local_field           Local value of field
 * \param [in]   local_weight          Local weight used to compute the mean value
 * \param [in]   global_mean_field_ptr Pointer where global mean field
 *                                     will be stored after computing
 */

void
PDM_global_mean_field_set
(
 const int          id,
 const int          i_part,
 const int          stride, 
 const double      *local_field,
 const double      *local_weight,
 double            *global_mean_field_ptr
);

void
PROCF (pdm_global_mean_field_set, PDM_GLOBAL_MEAN_FIELD_SET)
(
 const int         *id,
 const int         *i_part,
 const int         *stride, 
 const double      *local_field,
 const double      *local_weight,
 double            *global_mean_field_ptr
);


/**
 *
 * \brief Compute the global average field   
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_global_mean_field_compute
(
 const int          id
);

void
PROCF (pdm_global_mean_field_compute, PDM_GLOBAL_MEAN_FIELD_COMPUTE)
(
 const int         *id
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GLOBAL_POINT_MEAN_H__ */
