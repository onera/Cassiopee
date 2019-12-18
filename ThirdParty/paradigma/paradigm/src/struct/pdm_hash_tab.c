
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_hash_tab.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_gnum.h"

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
 * \struct _hash_tab_t
 * \brief  Hash table
 *
 * \ref _hash_tab_t defines a hash table structure
 *
 */

typedef struct {

  PDM_hash_tab_key_t   tKey;    /*!< Type of key */
  int                 *nDataKey;/*!< Number of data for each key */
  PDM_g_num_t           keyMax;  /*!< Key max */
  void              ***data;    /*!< Data */
  int                 *mDataKey;/*!< Max data for each key */
  int                  n_key_info; /*!< Number keys with information */
  int                  l_key_info; /*!< Size of \ref key_info*/
  PDM_g_num_t         *key_info; /*!< list of Keys with info */

} _hash_tab_t;

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
void                     *keyMax
)
{
  _hash_tab_t *ht = malloc (sizeof(_hash_tab_t));
  const int nDataDefault = 4;

  ht->tKey = tKey;

  if (tKey == PDM_HASH_TAB_KEY_INT) {
    ht->keyMax = (PDM_g_num_t) *((int *) keyMax);
  }
  else if (tKey == PDM_HASH_TAB_KEY_LONG) {
    ht->keyMax = *((PDM_g_num_t *) keyMax);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_hash_tab_create : Unknown hey type");
    abort();
  }

  ht->l_key_info = PDM_MAX (ht->keyMax/10, 2);
  ht->key_info = malloc (sizeof(PDM_g_num_t)*ht->l_key_info);
  ht->n_key_info = 0;

  ht->data = malloc (sizeof(void **) * ht->keyMax);
  ht->nDataKey = malloc (sizeof(int) * ht->keyMax);
  ht->mDataKey = malloc (sizeof(int) * ht->keyMax);
  for (int i = 0; i < ht->keyMax; i++) {
    ht->nDataKey[i] = 0;
    ht->mDataKey[i] = nDataDefault;
    ht->data[i] = malloc (sizeof(void *) * nDataDefault);
  }

  return (PDM_hash_tab_t *) ht;
}


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
PDM_hash_tab_data_add
(
PDM_hash_tab_t *ht,
void           *key,
void           *data
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  PDM_g_num_t _key = -1;;

  if (_ht->tKey == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  }
  else if (_ht->tKey == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  }
  else {
	  PDM_error(__FILE__, __LINE__, 0, "PDM_hash_tab_data_add error : unknown PDM_hash_tab_key_t\n");
		abort();
	}

  if (_ht->nDataKey[_key] >= _ht->mDataKey[_key]) {
    _ht->mDataKey[_key] *= 2;
    _ht->data[_key] = realloc (_ht->data[_key], sizeof(void *) *
                                                _ht->mDataKey[_key]);
  }

  if (_ht->nDataKey[_key] == 0) {
    if (_ht->n_key_info >= _ht->l_key_info) {
      _ht->l_key_info *= 2;
      _ht->key_info = realloc(_ht->key_info, sizeof(PDM_g_num_t) *  _ht->l_key_info);
    }
    _ht->key_info[_ht->n_key_info] = _key;
    _ht->n_key_info += 1;
  }

  _ht->data[_key][_ht->nDataKey[_key]++] = data;

}


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
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  PDM_g_num_t _key = -1;

  if (_ht->tKey == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  }
  else if (_ht->tKey == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  }
  else {
	  PDM_error(__FILE__, __LINE__, 0, "PDM_hash_tab_data_free error : unknown PDM_hash_tab_key_t\n");
		abort();
	}

  for (int i = 0; i < _ht->nDataKey[_key]; i++) {
    if (_ht->data[_key][i] != NULL) {
      free (_ht->data[_key][i]);
    }
    _ht->data[_key][i] = NULL;
  }

  _ht->nDataKey[_key] = 0;

}


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
PDM_hash_tab_n_data_get
(
PDM_hash_tab_t *ht,
void           *key
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  PDM_g_num_t _key = -1;

  if (_ht->tKey == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  }
  else if (_ht->tKey == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  }
  else {
	  PDM_error(__FILE__, __LINE__, 0, "PDM_hash_tab_data_get error : unknown PDM_hash_tab_key_t\n");
		abort();
	}
  return _ht->nDataKey[_key];
}


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
PDM_hash_tab_data_get
(
PDM_hash_tab_t *ht,
void           *key
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  PDM_g_num_t _key = -1;

  if (_ht->tKey == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  }
  else if (_ht->tKey == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  }
  else {
	  PDM_error(__FILE__, __LINE__, 0, "PDM_hash_tab_data_get error : unknown PDM_hash_tab_key_t\n");
		abort();
	}
  return _ht->data[_key];
}


/**
 * \brief Free a \ref PDM_hash_table_t object
 *
 * This function frees an \ref PDM_hash_table_t object
 *
 * \param [in]  hash_table    Hash table to free
 *
 * \return      NULL
 *
 */

PDM_hash_tab_t *
PDM_hash_tab_free
(
PDM_hash_tab_t *ht
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;

  for (int i = 0; i < _ht->keyMax; i++) {
    free (_ht->data[i]);
  }
  free (_ht->data);
  free (_ht->nDataKey);
  free (_ht->mDataKey);
  free (_ht->key_info);

  free (_ht);

  return NULL;
}


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
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;

  if (remove_data) {

    for (int i = 0; i < _ht->n_key_info; i++) {
      PDM_g_num_t _key = _ht->key_info[i];
      PDM_hash_tab_data_free (ht, &_key);
    }
  }

  int *_nDataKey = _ht->nDataKey;

  for (int i = 0; i < _ht->n_key_info; i++) {
    PDM_g_num_t _key = _ht->key_info[i];
    _nDataKey[_key] = 0;
  }

  _ht->n_key_info = 0;
}


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
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  return &(_ht->keyMax);
}


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
)
{
  _hash_tab_t *_ht = (_hash_tab_t *) ht;
  return _ht->tKey;
}


#ifdef	__cplusplus
}
#endif

