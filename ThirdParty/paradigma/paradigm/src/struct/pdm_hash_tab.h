#ifndef __PDM_HASH_TAB_H__
#define	__PDM_HASH_TAB_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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
 * \enum PDM_hash_table_key_t
 * \brief 2 Types of keys
 *
 */

typedef enum {

  PDM_HASH_TAB_KEY_INT  = 0,  /*!< Integer key  */
  PDM_HASH_TAB_KEY_LONG = 1,  /*!< Long key */

} PDM_hash_tab_key_t;

/**
 * \struct PDM_hash_tab_t
 * \brief Hash Table
 *
 *  PDM_hash_tab_t define a hasn table structure
 *
 */

typedef struct _hash_tab_t PDM_hash_tab_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an initialized hash table
 *
 * This function returns an initialized \ref PDM_hash_tab_t structure
 *
 * \param [in]  tKey   \ref PDM_HASH_TAB_KEY_INT or \ref PDM_HASH_TAB_KEY_LONG
 * \param [in]  keyMax  key max
 *
 * \return      A new initialized \ref PDM_hash_tab_t structure
 *
 */

PDM_hash_tab_t *
PDM_hash_tab_create
(
const PDM_hash_tab_key_t  tKey,
void                       *keyMax
);


/**
 * \brief Add a new data for a key
 *
 * This function adds a new data for a key
 *
 * \param [in]  ht        Hash table
 * \param [in]  key       key
 * \param [in]  data      data
 *
 */

void
PDM_hash_tab_data_add
(
PDM_hash_tab_t *ht,
void           *key,
void           *data
);


/**
 * \brief Free data for a key
 *
 *
 * \param [in]  ht        Hash table
 * \param [in]  key       key
 * \param [in]  data      data
 *
 */

void
PDM_hash_tab_data_free
(
PDM_hash_tab_t *ht,
void           *key
);


/**
 * \brief Get number of data associated to a key
 *
 * This function gets the number of data associated to a key
 *
 * \param [in]  ht       Hash table
 * \param [in]  key      key
 *
 * \return Number of data
 *
 */

int
PDM_hash_tab_n_data_get
(
PDM_hash_tab_t *ht,
void           *key
);


/**
 * \brief Get data associated to a key
 *
 * This function gets data associated to a key
 *
 * \param [in]  ht       Hash table
 * \param [in]  key      key
 *
 * \return data
 *
 */

void **
PDM_hash_tab_data_get
(
PDM_hash_tab_t *ht,
void           *key
);


/**
 * \brief Purge a \ref PDM_hash_table_t object
 *
 * This function empties an \ref PDM_hash_table_t object
 *
 * \param [in]  hash_table    Hash table to purge
 * \param [in]  remove_data   \ref PDM_FALSE or \ref PDM_TRUE
 *
 *
 */

void
PDM_hash_tab_purge
(
PDM_hash_tab_t *ht,
PDM_bool_t remove_data

);


/**
 * \brief Free a \ref PDM_hash_table_t object
 *
 * This function frees an \ref PDM_hash_table_t object
 *
 * \param [in]  ht    Hash table to free
 *
 * \return      NULL
 *
 */

PDM_hash_tab_t *
PDM_hash_tab_free
(
PDM_hash_tab_t *ht
);


/**
 * \brief Get maximum key of hash table
 *
 * This function returns the maximum key
 *
 * \param [in]  ht        Hash table
 *
 * \return max Key
 *
 */

void *
PDM_hash_tab_keyMax_get
(
PDM_hash_tab_t *ht
);


/**
 * \brief Get key type of hash table
 *
 * This function returns the key type
 *
 * \param [in]  ht        Hash table
 *
 * \return key type
 *
 */

PDM_hash_tab_key_t
PDM_hash_tab_keyType_get
(
PDM_hash_tab_t *ht
);


#ifdef	__cplusplus
}
#endif

#endif	/* __PDM_HASH_TAB_H__ */

