#ifndef __PDM_DHASH_TAB_H__
#define	__PDM_DHASH_TAB_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_dhash_table_key_t
 * \brief 2 Types of keys                
 *
 */
 
typedef enum {

  PDM_DHASH_TABLE_KEY_INT  = 0,  /*!< Integer key  */
  PDM_DHASH_TABLE_KEY_LONG = 1,  /*!< Long key */

} PDM_dhash_table_key_t;

/**
 * \struct PDM_dhash_tab_t
 * \brief Distributed hash Table 
 * 
 *  PDM_dhash_tab_t define a distributed hasn table structure
 *
 */

typedef struct _dhash_tab_t PDM_dhash_tab_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 * \brief Return an initialized distributed hash table
 *
 * This function returns an initialized \ref PDM_hash_tab_t structure
 *
 * \param [in]  tKey   \ref PDM_HASH_TAB_KEY_INT or \ref PDM_HASH_TAB_KEY_LONG
 * \param [in]  keyMax key max 
 *
 * \return      A new initialized \ref PDM_hash_tab_t structure
 *
 */

PDM_dhash_tab_t *
PDM_dhash_tab_create
(
const PDM_dhash_table_key_t tKey,
const PDM_MPI_Comm              comm,
void                        *keyMax
);


/**
 * \brief Store keys in continuous block  
 *
 * This function stores keys in continuous block and distributes them on their
 * associated processes. 
 *
 * \param [in]  hash_table    Hash table
 *
 */

void
PDM_dhash_block_store
(
PDM_hash_tab_t *hash_table
);


/**
 * \brief Store keys in continuous block  
 *
 * This function returns the current block properties 
 *
 * \param [in]  hash_table    Hash table
 * \param [out] key_min       Min key in the current block 
 * \param [out] key_max       Max key in the current block
 *
 */

void
PDM_dhash_block_properties_get
(
PDM_hash_tab_t *hash_table,
void           *key_min,
void           *key_max
);


/**
 * \brief Add a new data for a key
 *
 * This function adds a new data for a key
 *
 * \param [in]  hash_table    Hash table
 * \param [in]  key           key 
 * \param [in]  data          data
 *
 */

void
PDM_dhash_tab_data_add
(
PDM_hash_tab_t *hash_table,
void           *key,
void           *data
);


/**
 * \brief Get number of data associated to a key
 *
 * This function gets the number of data associated to a key
 *
 * \param [in]  hash_table    Hash table
 * \param [in]  key           key 
 *
 * \return Number of data 
 *    
 */

int
PDM_dhash_tab_n_data_get
(
PDM_hash_tab_t *hash_table,
void           *key
);


/**
 * \brief Get data associated to a key
 *
 * This function gets data associated to a key
 *
 * \param [in]  hash_table    Hash table
 * \param [in]  key           key 
 *
 * \return data 
 *    
 */

void **
PDM_dhash_tab_data_get
(
PDM_hash_tab_t *hash_table,
void           *key
);


/**
 * \brief Free a \ref PDM_dhash_table_t object
 *
 * This function frees an \ref PDM_dhash_table_t object
 *
 * \param [in]  hash_table    Hash table to free
 *
 * \return      NULL
 *
 */

PDM_dhash_tab_t *
PDM_dhash_table_free
(
PDM_dhash_tab_t *hash_table
);


#ifdef	__cplusplus
}
#endif

#endif	/* __PDM_DHASH_TAB_H__ */

