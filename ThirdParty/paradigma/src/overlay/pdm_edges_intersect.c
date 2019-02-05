/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_line.h"
#include "pdm_plane.h"
#include "pdm_edges_intersect.h"
#include "pdm_edges_intersect_priv.h"
#include "pdm_hash_tab.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_mpi.h"
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 *
 * \brief Create a new \ref _edges_intersect_res_t object
 *
 * \param [in]   nGEdgeA    Global number mesh A edge
 * \param [in]   nGEdgeB    Global number mesh B edge
 * \param [in]   nNewPointA Number of new point for mesh A
 * \param [in]   nNewPointB Number of new point for mesh B
 *
 * \return      A new \ref _edges_intersect_res_t     
 */

static _edges_intersect_res_t *
_edges_intersect_res_create
(
const PDM_g_num_t   nGEdgeA,
const PDM_g_num_t   nGEdgeB,
const int          nNewPointsA,
const int          nNewPointsB        
)
{
  _edges_intersect_res_t *newInter = malloc (sizeof(_edges_intersect_res_t));

  newInter->nGEdgeA     = nGEdgeA;
  newInter->nGEdgeB     = nGEdgeB;
  newInter->originEdgeA = -1;
  newInter->originEdgeB = -1;
  newInter->tIntersect  = PDM_LINE_INTERSECT_UNDEF;

  newInter->nNewPointsA = nNewPointsA;
  newInter->uA          = malloc (sizeof(double) * nNewPointsA);
  newInter->coordsA     = malloc (sizeof(double) * 3 * nNewPointsA);
  newInter->linkA       = malloc (sizeof(PDM_g_num_t) * nNewPointsA);
  newInter->gNumA       = malloc (sizeof(PDM_g_num_t) * nNewPointsA);
  newInter->oNewPointsA = malloc (sizeof(PDM_edges_intersect_point_t) * nNewPointsA);

  newInter->nNewPointsB = nNewPointsB;
  newInter->uB          = malloc (sizeof(double) * nNewPointsB);
  newInter->coordsB     = malloc (sizeof(double) * 3 * nNewPointsB);
  newInter->linkB       = malloc (sizeof(PDM_g_num_t) * nNewPointsB);
  newInter->gNumB       = malloc (sizeof(PDM_g_num_t) * nNewPointsB);
  newInter->oNewPointsB = malloc (sizeof(PDM_edges_intersect_point_t) * nNewPointsB);
  
  return newInter;
}


/**
 *
 * \brief Free a \ref _edges_intersect_res_t object
 *
 * \param [in]   eir  Edges intersection results object
 *
 * \return NULL     
 */

static _edges_intersect_res_t *
_edges_intersect_res_free (_edges_intersect_res_t *eir)
{
  if (eir->linkA != NULL) {
    free (eir->gNumA);
    free (eir->linkA);
    free (eir->uA);
    free (eir->coordsA);
    free (eir->oNewPointsA);
    eir->gNumA      = NULL;
    eir->linkA      = NULL;
    eir->uA           = NULL;
    eir->oNewPointsA  = NULL;
    eir->coordsA      = NULL;
  }
  eir->nNewPointsB      = 0;
  if (eir->linkB != NULL) {
    free (eir->gNumB);
    free (eir->linkB);
    free (eir->uB);
    free (eir->coordsB);
    free (eir->oNewPointsB);
    eir->gNumB      = NULL;
    eir->linkB      = NULL;
    eir->uB           = NULL;
    eir->oNewPointsB  = NULL;
    eir->coordsB      = NULL;
  }
  free (eir);
  return NULL;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a new \ref PDM_edges_intersect_t object
 *
 * \param [in]   maxGNEdgeA    Max global number of edges in mesh A
 * \param [in]   maxGNEdgeN    Max global number of edges in mesh B
 * \param [in]   vtxCarLengthTol Absolute tolerance for characteristic length 
 * \param [in]   sMSGComm        size of mpicomm
 *
 * \return      A new \ref PDM_edges_intersect_t     
 */

PDM_edges_intersect_t *
PDM_edges_intersect_create 
(
const PDM_g_num_t maxGNEdgeA,
const PDM_g_num_t maxGNEdgeB,
const double     vtxCarLengthTol,
const PDM_MPI_Comm   comm        
)
{
  int sMSGComm;
  PDM_MPI_Comm_size (comm, &sMSGComm);
  
  _edges_intersect_t *ei = malloc( sizeof(_edges_intersect_t));
  
  PDM_g_num_t _keyMax = (maxGNEdgeA + maxGNEdgeB) % sMSGComm;
  int keyMax = (int) _keyMax;
  ei->ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMax);
  
  PDM_g_num_t _keyMaxA = maxGNEdgeA % sMSGComm;
  int keyMaxA = (int) _keyMaxA;
  ei->htA = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMaxA);

  PDM_g_num_t _keyMaxB = maxGNEdgeB % sMSGComm;
  int keyMaxB = (int) _keyMaxB;
  ei->htB = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMaxB);
  
  ei->maxGNEdgeA = maxGNEdgeA;
  ei->maxGNEdgeB = maxGNEdgeB;
  ei->vtxCarLengthTol = vtxCarLengthTol;
  ei->comm = comm;
  ei->sMSGComm = sMSGComm;

  return (PDM_edges_intersect_t *) ei;
}


/**
 *
 * \brief Get result of the intersection 
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   get_t            Type of key to return data
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [out]  n_intersect      Number of intersections
 *
 * \return    Result of the intersection     
 */

PDM_edges_intersect_res_t **
PDM_edges_intersect_get 
(
PDM_edges_intersect_t  *ei,
PDM_edges_get_t         get_t,        
const PDM_g_num_t       nGEdgeA,
const PDM_g_num_t       nGEdgeB,
int                    *n_intersect        
)
{ 
  
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;
  
  if (get_t == PDM_EDGES_GET_FROM_AB) {
  
    PDM_hash_tab_t *ht = _ei->ht;
    PDM_g_num_t _key = (nGEdgeA + nGEdgeB) % _ei->sMSGComm;
    int key = (int) _key;

    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas = 
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeA == data->nGEdgeB) {
        *n_intersect = 1;
        _edges_intersect_res_t **_datas_edge = 
            malloc(sizeof(_edges_intersect_res_t *));
        _datas_edge[0] = data;
        return (PDM_edges_intersect_res_t **) _datas_edge;
      }
    }
    
    PDM_error(__FILE__, __LINE__, 0,"Error PDM_edges_intersect_get : "
                   "Requested intersection not found "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n", nGEdgeA, nGEdgeB);
    abort();
  }

  else if (get_t == PDM_EDGES_GET_FROM_A) {
  
    PDM_hash_tab_t *ht = _ei->htA;
    PDM_g_num_t _key = nGEdgeA % _ei->sMSGComm;
    int key = (int) _key;

    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas = 
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    int _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeA == nGEdgeA) {
        _nData += 1;
      }
    }
    
    _edges_intersect_res_t **_datas_edge = 
            malloc(sizeof(_edges_intersect_res_t *) * _nData);

    _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeA == nGEdgeA) {
        _datas_edge[_nData++] = data;
      }
    }
    *n_intersect = _nData;
    return (PDM_edges_intersect_res_t **) _datas_edge;
  }

  else if (get_t == PDM_EDGES_GET_FROM_B) {
  
    PDM_hash_tab_t *ht = _ei->htB;
    PDM_g_num_t _key = nGEdgeB % _ei->sMSGComm;
    int key = (int) _key;
    
    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas = 
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    int _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeB == nGEdgeB) {
        _nData += 1;
      }
    }
    
    _edges_intersect_res_t **_datas_edge = 
            malloc(sizeof(_edges_intersect_res_t *) * _nData);

    _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeB == nGEdgeB) {
        _datas_edge[_nData++] = data;
      }
    }
    *n_intersect = _nData;
    return (PDM_edges_intersect_res_t **) _datas_edge;
  }
  
  return NULL;
}


/**
 *
 * \brief Perform an intersection between a meshA edge and a meshB edge 
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGVtxA           Global number of edgeA vertices
 * \param [in]   charLgthVtxA[2]  Characteristic length of edgeA vertices     
 * \param [in]   coordsVtxA[6]    Coordinates of edgeA vertices
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [in]   nGVtxB           Global number of edgeB vertices
 * \param [in]   charLgthVtxB[2]  Characteristic length of edgeB vertices     
 * \param [in]   coordsVtxB[6]    Coordinates of edgeB vertices
 *
 * \return    Result of the intersection     
 */

PDM_edges_intersect_res_t *
PDM_edges_intersect_add 
(
PDM_edges_intersect_t       *ei,
const PDM_g_num_t             nGEdgeA,
const PDM_g_num_t             nGVtxA[2],
const double                 charLgthVtxA[2],
const double                 coordsVtxA[6],
const PDM_g_num_t             nGEdgeB,
const PDM_g_num_t             nGVtxB[2],
const double                 charLgthVtxB[2],
const double                 coordsVtxB[6]
)
{ 
  
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;
  PDM_hash_tab_t *ht = _ei->ht;
  PDM_hash_tab_t *htA = _ei->htA;
  PDM_hash_tab_t *htB = _ei->htB;
  
  const double minMin = 1e-18;
    
  double _charLgthVtxA[2] = {PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxA[0], minMin),
                             PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxA[1], minMin)};
  
  double _charLgthVtxB[2] = {PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxB[0], minMin),
                             PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxB[1], minMin)}; 
  
  PDM_g_num_t _key  = (nGEdgeA + nGEdgeB) % _ei->sMSGComm;
  PDM_g_num_t _keyA = nGEdgeA % _ei->sMSGComm;
  PDM_g_num_t _keyB = nGEdgeB % _ei->sMSGComm;
  
  int key = (int) _key;
  int keyA = (int) _keyA;
  int keyB = (int) _keyB;
  
  const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
  _edges_intersect_res_t ** datas = 
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                               (void *) &key);
  /**********************************************
   * Check if intersection is already preformed *
   **********************************************/
  
  for (int i = 0; i < nData; i++) {
    _edges_intersect_res_t *data = datas[i];
    if (data->nGEdgeA == data->nGEdgeB) {
      return (PDM_edges_intersect_res_t *) data;
    }
  }

  double u1;
  double v1;
  
  /******************************
   * Perform a new intersection *
   ******************************/

  /*
   * Line-line intersection
   */

  PDM_line_intersect_t tIntersect =
          PDM_line_intersection (coordsVtxA, &(coordsVtxA[3]),
                                 coordsVtxB, &(coordsVtxB[3]),
                                 &u1, &v1);
  
  bool isInitialOnLine = false;
  if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {
    isInitialOnLine = true;
  }
  
  /*
   * Initialization
   */

  _edges_intersect_res_t *newInter = NULL;
  
  double vA[3];
  double vB[3];
    
  for (int i = 0; i < 3; i++) {
    vA[i] = coordsVtxA[3 + i] - coordsVtxA[i]; 
    vB[i] = coordsVtxB[3 + i] - coordsVtxB[i]; 
  }    
 
  double vA_norm = PDM_MODULE (vA); 
  double vB_norm = PDM_MODULE (vB); 
  
  double nVA[3];
  double nVB[3];
  
  for (int i = 0; i < 3; i++) {
    nVB[i] = vB[i] / vB_norm;
    nVA[i] = vA[i] / vA_norm;
  }
  
  double dVtxA[2] = {PDM_ABS (u1 * vA_norm), 
                     PDM_ABS ((1. - u1) * vA_norm)}; 
  
  double dVtxB[2] = {PDM_ABS (v1 * vB_norm), 
                     PDM_ABS ((1. - v1) * vB_norm)};


  bool isSetted = false;
  
//  double coordInter[3] = {coordsVtxA[0] + u1 * vA[0],
//                          coordsVtxA[1] + u1 * vA[1],
//                          coordsVtxA[2] + u1 * vA[2]};   

    
//  double A1Inter[3] = {coordInter[0] - coordsVtxA[0],
//                       coordInter[1] - coordsVtxA[1],
//                       coordInter[2] - coordsVtxA[2]};
//
//  double A2Inter[3] = {coordInter[0] - coordsVtxA[3],
//                       coordInter[1] - coordsVtxA[4],
//                       coordInter[2] - coordsVtxA[5]};
//
//  double B1Inter[3] = {coordInter[0] - coordsVtxB[0],
//                       coordInter[1] - coordsVtxB[1],
//                       coordInter[2] - coordsVtxB[2]};
//
//  double B2Inter[3] = {coordInter[0] - coordsVtxB[3],
//                       coordInter[1] - coordsVtxB[4],
//                       coordInter[2] - coordsVtxB[5]};

//  double dA1Inter = PDM_MODULE (A1Inter);
//  double dA2Inter = PDM_MODULE (A2Inter);
//  double dB1Inter = PDM_MODULE (B1Inter);
//  double dB2Inter = PDM_MODULE (B2Inter);

  /*
   * Resolution of inconsistencies
   */

  if ((tIntersect == PDM_LINE_INTERSECT_NO) ||
      (tIntersect == PDM_LINE_INTERSECT_YES)) {
    
    /* 
     * Check if A vertex and B vertex are the same
     */

    int isA1B1 = (dVtxA[0] < _charLgthVtxA[0]) &&
                 (dVtxB[0] < _charLgthVtxB[0]);
    
    int isA1B2 = (dVtxA[0] < _charLgthVtxA[0]) &&
                 (dVtxB[1] < _charLgthVtxB[1]);

    int isA2B1 = (dVtxA[1] < _charLgthVtxA[1]) &&
                 (dVtxB[0] < _charLgthVtxB[0]);
    
    int isA2B2 = (dVtxA[1] < _charLgthVtxA[1]) &&
                 (dVtxB[1] < _charLgthVtxB[1]);
    
    /* 
     * Check if A vertex is on B Edge and vice versa
     */
    
    double closestA1EdgeB[3];
    double tA1EdgeB;
    double d2A1EdgeB = PDM_line_distance (coordsVtxA,
                                          coordsVtxB,
                                          coordsVtxB + 3,
                                          &tA1EdgeB,
                                          closestA1EdgeB);
    
    double closestA2EdgeB[3];
    double tA2EdgeB;
    double d2A2EdgeB = PDM_line_distance (coordsVtxA + 3,
                                          coordsVtxB,
                                          coordsVtxB + 3,
                                          &tA2EdgeB,
                                          closestA2EdgeB);
    
    double closestB1EdgeA[3];
    double tB1EdgeA;
    double d2B1EdgeA = PDM_line_distance (coordsVtxB,
                                          coordsVtxA,
                                          coordsVtxA + 3,
                                          &tB1EdgeA,
                                          closestB1EdgeA);
    
    double closestB2EdgeA[3];
    double tB2EdgeA;
    double d2B2EdgeA = PDM_line_distance (coordsVtxB + 3,
                                          coordsVtxA,
                                          coordsVtxA + 3,
                                          &tB2EdgeA,
                                          closestB2EdgeA);
     
    double _deltaEdgA = _charLgthVtxA[1] - _charLgthVtxA[0];
    double _deltaEdgB = _charLgthVtxB[1] - _charLgthVtxB[0];
    
    double A1closestEdgeB[3] = {closestA1EdgeB[0] - coordsVtxA[0],
                                closestA1EdgeB[1] - coordsVtxA[1],
                                closestA1EdgeB[2] - coordsVtxA[2]};
    double dA1closestEdgeB = PDM_MODULE (A1closestEdgeB);
    
    double A2closestEdgeB[3] = {closestA2EdgeB[0] - coordsVtxA[3+0],
                                closestA2EdgeB[1] - coordsVtxA[3+1],
                                closestA2EdgeB[2] - coordsVtxA[3+2]};
    double dA2closestEdgeB = PDM_MODULE (A2closestEdgeB);


    double B1closestEdgeA[3] = {closestB1EdgeA[0] - coordsVtxB[0],
                                closestB1EdgeA[1] - coordsVtxB[1],
                                closestB1EdgeA[2] - coordsVtxB[2]};
    double dB1closestEdgeA = PDM_MODULE (B1closestEdgeA);
    
    double B2closestEdgeA[3] = {closestB2EdgeA[0] - coordsVtxB[3+0],
                                closestB2EdgeA[1] - coordsVtxB[3+1],
                                closestB2EdgeA[2] - coordsVtxB[3+2]};
    double dB2closestEdgeA = PDM_MODULE (B2closestEdgeA);
     
    int isA1OnEdgeB = (d2A1EdgeB < ((_charLgthVtxB[0] + tA1EdgeB * _deltaEdgB) *
                                    (_charLgthVtxB[0] + tA1EdgeB * _deltaEdgB))) && 
                                    (dA1closestEdgeB < _charLgthVtxA[0]);
     
    int isA2OnEdgeB = (d2A2EdgeB < ((_charLgthVtxB[0] + tA2EdgeB * _deltaEdgB) *
                                    (_charLgthVtxB[0] + tA2EdgeB * _deltaEdgB))) && 
                                    (dA2closestEdgeB < _charLgthVtxA[1]);
     
    int isB1OnEdgeA = (d2B1EdgeA < ((_charLgthVtxA[0] + tB1EdgeA * _deltaEdgA) *
                                    (_charLgthVtxA[0] + tB1EdgeA * _deltaEdgA))) && 
                                    (dB1closestEdgeA < _charLgthVtxB[0]);
     
    int isB2OnEdgeA = (d2B2EdgeA < ((_charLgthVtxA[0] + tB2EdgeA * _deltaEdgA) *
                                    (_charLgthVtxA[0] + tB2EdgeA * _deltaEdgA))) && 
                                    (dB2closestEdgeA < _charLgthVtxB[1]);
    /* 
     * Check if A Edge and B are the same line
     */

    if ((isA1B1 && isA2B2) || (isA2B1 && isA1B2)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE; 
    }

    if ((isA1OnEdgeB && isA2OnEdgeB) || (isB1OnEdgeA && isB2OnEdgeA)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE; 
    }

    if ((isA1B1 && isA2OnEdgeB) || (isA1B1 && isB2OnEdgeA) ||
        (isA2B1 && isA1OnEdgeB) || (isA2B1 && isB2OnEdgeA) ||
        (isA1B2 && isA2OnEdgeB) || (isA1B2 && isB1OnEdgeA) ||
        (isA2B2 && isA1OnEdgeB) || (isA2B2 && isB1OnEdgeA)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE; 
    }

    if (tIntersect != PDM_LINE_INTERSECT_ON_LINE) {
      
      /* 
       * Define intersection for same vertex inconsistencies
       */

      if (isA1B1 || isA1B2 || isA2B1 || isA2B2) {

        int nNewPointsA = 1;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);


        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        
        
        int iA = 0;
        int iB = 0;
                
        if (isA1B1) {

          iA = 0;
          iB = 0;
          
        }

        else if (isA1B2) {

          iA = 0;
          iB = 1;

        }

        else if (isA2B1) {

          iA = 1;
          iB = 0;

        }

        else if (isA2B2) {

          iA = 1;
          iB = 1;

        }

        if (newInter->originEdgeA == nGVtxA[1]) {
          iA = (iA + 1) % 2; 
        }  

        if (newInter->originEdgeB == nGVtxB[1]) {
          iB = (iB + 1) % 2; 
        }  

        newInter->uA[0] = (double) iA;
        newInter->uB[0] = (double) iB;

        newInter->linkA[0] = nGVtxB[iB];
        newInter->linkB[0] = nGVtxA[iA];
        newInter->gNumA[0] = 0;
        newInter->gNumB[0] = 0;

        for (int i = 0; i < 3; i++) {
          double _comp = (coordsVtxB[3*iB+i] + coordsVtxA[3*iA+i]) / 2;
          newInter->coordsA[i] = _comp;
          newInter->coordsB[i] = _comp;   
        }

        isSetted = true;

      }    

      /* 
       * Define intersection for vertex on edge inconsistencies
       */

      else if (isA1OnEdgeB) {
        int nNewPointsA = 0;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;
        
        double B1ClosestA1[3] = {closestA1EdgeB[0] - coordsVtxB[0],
                                 closestA1EdgeB[1] - coordsVtxB[1],
                                 closestA1EdgeB[2] - coordsVtxB[2]};
        
        double  mB1ClosestA1 = PDM_MODULE (B1ClosestA1);

        for (int i = 0; i < 3; i++) {
          B1ClosestA1[i] = B1ClosestA1[i] / mB1ClosestA1;
        }

        newInter->uB[0] = PDM_DOT_PRODUCT (B1ClosestA1, nVB);
        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
 
        for (int i = 0; i < 3; i++) {
          newInter->coordsB[i] = closestA1EdgeB[i];
        }

        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumB[0] = 0;

        isSetted = true;

      }

      else if (isA2OnEdgeB) {

        int nNewPointsA = 0;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;
        
        double B1ClosestA2[3] = {closestA2EdgeB[0] - coordsVtxB[0],
                                 closestA2EdgeB[1] - coordsVtxB[1],
                                 closestA2EdgeB[2] - coordsVtxB[2]};
        
        double  mB1ClosestA2 = PDM_MODULE (B1ClosestA2);

        for (int i = 0; i < 3; i++) {
          B1ClosestA2[i] = B1ClosestA2[i] / mB1ClosestA2;
        }

        newInter->uB[0] = PDM_DOT_PRODUCT (B1ClosestA2, nVB);
        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
 
        for (int i = 0; i < 3; i++) {
          newInter->coordsB[i] = closestA2EdgeB[i];
        }

        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumB[0] = 0;

        isSetted = true;

      }

      else if (isB1OnEdgeA) {

        int nNewPointsA = 1;
        int nNewPointsB = 0;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

        double A1ClosestB1[3] = {closestB1EdgeA[0] - coordsVtxA[0],
                                 closestB1EdgeA[1] - coordsVtxA[1],
                                 closestB1EdgeA[2] - coordsVtxA[2]};
        
        double  mA1ClosestB1 = PDM_MODULE (A1ClosestB1);

        for (int i = 0; i < 3; i++) {
          A1ClosestB1[i] = A1ClosestB1[i] / mA1ClosestB1;
        }

        newInter->uA[0] = PDM_DOT_PRODUCT (A1ClosestB1, nVA);
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }
 
        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = closestB1EdgeA[i];
        }

        newInter->linkA[0] = nGVtxB[0];
        newInter->gNumA[0] = 0;

        isSetted = true;

      }

      else if (isB2OnEdgeA) {


        int nNewPointsA = 1;
        int nNewPointsB = 0;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

        double A1ClosestB2[3] = {closestB2EdgeA[0] - coordsVtxA[0],
                                 closestB2EdgeA[1] - coordsVtxA[1],
                                 closestB2EdgeA[2] - coordsVtxA[2]};
        
        double  mA1ClosestB2 = PDM_MODULE (A1ClosestB2);

        for (int i = 0; i < 3; i++) {
          A1ClosestB2[i] = A1ClosestB2[i] / mA1ClosestB2;
        }

        newInter->uA[0] = PDM_DOT_PRODUCT (A1ClosestB2, nVA);
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }
 
        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = closestB2EdgeA[i];
        }

        newInter->linkA[0] = nGVtxB[1];
        newInter->gNumA[0] = 0;

        isSetted = true;

      }
    }  
  } 
 
  
  if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {
  
    if (isInitialOnLine) {
      int isA1B1 = (dVtxA[0] < _charLgthVtxA[0]) &&
                   (dVtxB[0] < _charLgthVtxB[0]);

      if (!isA1B1) {
  
        double A1B1[3] = {coordsVtxB[0] - coordsVtxA[0],
                          coordsVtxB[1] - coordsVtxA[1],
                          coordsVtxB[2] - coordsVtxA[2]};
        
        double mA1B1 = PDM_MODULE (A1B1);
        
        for (int i = 0; i < 3; i++) {
          A1B1[i] = A1B1[i] / mA1B1;
        }
        
        int isSameLine = (PDM_ABS (PDM_DOT_PRODUCT (A1B1, nVA)) > 0.99);
               
        if (!isSameLine) {
          tIntersect = PDM_LINE_INTERSECT_NO;
          newInter = _edges_intersect_res_create (nGEdgeA, 
                                                  nGEdgeB, 
                                                  0,
                                                  0);
          newInter->tIntersect = tIntersect;
          newInter->originEdgeA = nGVtxA[0];
          newInter->originEdgeB = nGVtxB[0];

          isSetted = true;

        }
      }
    }
    
    if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {
      
      double A1B1[3] = {coordsVtxB[0] - coordsVtxA[0],
                        coordsVtxB[1] - coordsVtxA[1],
                        coordsVtxB[2] - coordsVtxA[2]};
      double mA1B1 = PDM_MODULE (A1B1);

      double A1B2[3] = {coordsVtxB[3+0] - coordsVtxA[0],
                        coordsVtxB[3+1] - coordsVtxA[1],
                        coordsVtxB[3+2] - coordsVtxA[2]};
      double mA1B2 = PDM_MODULE (A1B2);

      double B1A2[3] = {coordsVtxA[3+0] - coordsVtxB[0],
                        coordsVtxA[3+1] - coordsVtxB[1],
                        coordsVtxA[3+2] - coordsVtxB[2]};
      double mB1A2 = PDM_MODULE (B1A2);
        
      for (int i = 0; i < 3; i++) {
        A1B1[i] = A1B1[i] / mA1B1;
        A1B2[i] = A1B2[i] / mA1B2;
        B1A2[i] = B1A2[i] / mB1A2;
      }
      
      int isA1B1 = (dVtxA[0] < _charLgthVtxA[0]) &&
                   (dVtxB[0] < _charLgthVtxB[0]);
    
      int isA1B2 = (dVtxA[0] < _charLgthVtxA[0]) &&
                   (dVtxB[1] < _charLgthVtxB[1]);

      int isA2B1 = (dVtxA[1] < _charLgthVtxA[1]) &&
                   (dVtxB[0] < _charLgthVtxB[0]);
    
      int isA2B2 = (dVtxA[1] < _charLgthVtxA[1]) &&
                   (dVtxB[1] < _charLgthVtxB[1]);

      int nNewPointsA = 0;
      int nNewPointsB = 0;

/*
      double *coordsNewPointsA = NULL;             
      double *coordsNewPointsB = NULL;
              
      PDM_g_num_t *gNNewPointsA = NULL;                 
      PDM_g_num_t *gNNewPointsB = NULL;                 
*/
     
      double closestB2EdgeA[3];
      double tB2EdgeA;
      
      PDM_line_distance (coordsVtxB + 3,
                         coordsVtxA,
                         coordsVtxA + 3,
                        &tB2EdgeA,
                         closestB2EdgeA);
          
      double A1ClosestB2[3] = {closestB2EdgeA[0] - coordsVtxA[0],
                               closestB2EdgeA[1] - coordsVtxA[1],
                               closestB2EdgeA[2] - coordsVtxA[2]};
       
      double mA1ClosestB2 = PDM_MODULE (A1ClosestB2);

      for (int i = 0; i < 3; i++) {
        A1ClosestB2[i] = A1ClosestB2[i] / mA1ClosestB2;
      }
    
      double closestA2EdgeB[3];
      double tA2EdgeB;
      PDM_line_distance (coordsVtxA + 3,
                         coordsVtxB,
                         coordsVtxB + 3,
                         &tA2EdgeB,
                         closestA2EdgeB);

      double B1ClosestA2[3] = {closestA2EdgeB[0] - coordsVtxB[0],
                               closestA2EdgeB[1] - coordsVtxB[1],
                               closestA2EdgeB[2] - coordsVtxB[2]};

      double  mB1ClosestA2 = PDM_MODULE (B1ClosestA2);

      for (int i = 0; i < 3; i++) {
        B1ClosestA2[i] = B1ClosestA2[i] / mB1ClosestA2;
      }
     
      double closestB1EdgeA[3];
      double tB1EdgeA;
      PDM_line_distance (coordsVtxB,
                         coordsVtxA,
                         coordsVtxA + 3,
                        &tB1EdgeA,
                         closestB1EdgeA);
          
      double A1ClosestB1[3] = {closestB1EdgeA[0] - coordsVtxA[0],
                               closestB1EdgeA[1] - coordsVtxA[1],
                               closestB1EdgeA[2] - coordsVtxA[2]};
       
      double mA1ClosestB1 = PDM_MODULE (A1ClosestB1);

      for (int i = 0; i < 3; i++) {
        A1ClosestB1[i] = A1ClosestB1[i] / mA1ClosestB1;
      }

      
      double closestA1EdgeB[3];
      double tA1EdgeB;
      PDM_line_distance (coordsVtxA,
                         coordsVtxB,
                         coordsVtxB + 3,
                        &tA1EdgeB,
                         closestA1EdgeB);

      double B1ClosestA1[3] = {closestA1EdgeB[0] - coordsVtxB[0],
                               closestA1EdgeB[1] - coordsVtxB[1],
                               closestA1EdgeB[2] - coordsVtxB[2]};

      double  mB1ClosestA1 = PDM_MODULE (B1ClosestA1);

      for (int i = 0; i < 3; i++) {
        B1ClosestA1[i] = B1ClosestA1[i] / mB1ClosestA1;
      }
      
      
      if ((isA1B1 && isA2B2) || (isA2B1 && isA1B2)) {
        
        nNewPointsA = 2;
        nNewPointsB = 2;

        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
                
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        
        if (isA1B1 && isA2B2) {
          newInter->uA[0] = 0.;
          newInter->uA[1] = 1.;
          newInter->uB[0] = 0.;
          newInter->uB[1] = 1.;
          
          newInter->linkA[0] = nGVtxB[0];
          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[0] = 0;
          newInter->gNumA[1] = 0;
       
          newInter->linkB[0] = nGVtxA[0];
          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[0] = 0;
          newInter->gNumB[1] = 0;
          
          double _coords1[3];
          double _coords2[3];
          
          for (int i = 0; i < 3; i++) {
            _coords1[i] = (coordsVtxB[i] + coordsVtxA[i]) / 2;
            _coords2[i] = (coordsVtxB[3+i] + coordsVtxA[3+i]) / 2;
          }

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[i]     = _coords1[i]; 
            newInter->coordsB[i]     = _coords1[i]; 
            newInter->coordsA[3 + i] = _coords2[i]; 
            newInter->coordsB[3 + i] = _coords2[i]; 
          }
          
        }

        else {
          newInter->uA[0] = 1.;
          newInter->uA[1] = 0.;
          newInter->uB[0] = 1.;
          newInter->uB[1] = 0.;
          
          newInter->linkA[0] = nGVtxB[1];
          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[0] = 0;
          newInter->gNumA[1] = 0;
        
          newInter->linkB[0] = nGVtxA[1];
          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[0] = 0;
          newInter->gNumB[1] = 0;

          double _coords1[3];
          double _coords2[3];
          
          for (int i = 0; i < 3; i++) {
            _coords1[i] = (coordsVtxB[3+i] + coordsVtxA[  i]) / 2;
            _coords2[i] = (coordsVtxB[  i] + coordsVtxA[3+i]) / 2;
          }

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[i]     = _coords1[i]; 
            newInter->coordsB[i]     = _coords2[i]; 
            newInter->coordsA[3 + i] = _coords2[i]; 
            newInter->coordsB[3 + i] = _coords1[i]; 
          }

        }
        
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
          newInter->uA[1] = 1 - newInter->uA[1];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
          newInter->uB[1] = 1 - newInter->uB[1];
        }
        
      }
      
      else if (isA1B1) {
        
        nNewPointsA = 1;
        nNewPointsB = 1;
        
        double u = PDM_DOT_PRODUCT (A1ClosestB2, nVA);
        double v = PDM_DOT_PRODUCT (B1ClosestA2, nVB);
        
        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }
        
        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->linkA[0] = nGVtxB[0];
        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumA[0] = 0;
        newInter->gNumB[0] = 0;
 
        newInter->uA[0] = 0.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
                
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                
        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[i] + coordsVtxA[i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i]; 
          newInter->coordsB[i] = _coords[i]; 
        }
        
        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestB2EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA2EdgeB[i];
          }
        }
      }

      else if (isA2B2) {

        nNewPointsA = 1;
        nNewPointsB = 1;
        
        double u = PDM_DOT_PRODUCT (A1ClosestB1, nVA);
        double v = PDM_DOT_PRODUCT (B1ClosestA1, nVB);
        
        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }
        
        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->linkA[0] = nGVtxB[1];
        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumA[0] = 0;
        newInter->gNumB[0] = 0;
 
        newInter->uA[0] = 1.;
        newInter->uB[0] = 1.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
                
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                
        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[3+i] + coordsVtxA[3+i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i]; 
          newInter->coordsB[i] = _coords[i]; 
        }
        
        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestB1EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA1EdgeB[i];
          }
        }

      }

      else if (isA1B2) {

        nNewPointsA = 1;
        nNewPointsB = 1;
        
        double u = PDM_DOT_PRODUCT (A1ClosestB1, nVA);
        double v = PDM_DOT_PRODUCT (B1ClosestA2, nVB);
        
        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }
        
        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->linkA[0] = nGVtxB[1];
        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumA[0] = 0;
        newInter->gNumB[0] = 0;
 
        newInter->uA[0] = 1.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
                
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                
        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[3+i] + coordsVtxA[i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i]; 
          newInter->coordsB[i] = _coords[i]; 
        }
        
        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestB1EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA2EdgeB[i];
          }
        }

      }
      
      else if (isA2B1) {

        nNewPointsA = 1;
        nNewPointsB = 1;
        
        double u = PDM_DOT_PRODUCT (A1ClosestB2, nVA);
        double v = PDM_DOT_PRODUCT (B1ClosestA1, nVB);
        
        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }
        
        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        newInter->linkA[0] = nGVtxB[0];
        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumA[0] = 0;
        newInter->gNumB[0] = 0;
 
        newInter->uA[0] = 1.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }
                
        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                
        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[i] + coordsVtxA[3+i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i]; 
          newInter->coordsB[i] = _coords[i]; 
        }
        
        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestB2EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA1EdgeB[i];
          }
        }

      }

      else {
        nNewPointsA = 0;
        nNewPointsB = 0;

        double u[2] = {PDM_DOT_PRODUCT (A1ClosestB1, nVA),
                       PDM_DOT_PRODUCT (A1ClosestB2, nVA)};
        double v[2] = {PDM_DOT_PRODUCT (B1ClosestA1, nVB), 
                       PDM_DOT_PRODUCT (B1ClosestA2, nVB)};

        if ((u[0] <= 1.) && (u[0] >= 0.)) {
          nNewPointsA += 1;
        }

        if ((u[1] <= 1.) && (u[1] >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v[0] <= 1.) && (v[0] >= 0.)) {
          nNewPointsB += 1;
        }

        if ((v[1] <= 1.) && (v[1] >= 0.)) {
          nNewPointsB += 1;
        }
                
        newInter = _edges_intersect_res_create (nGEdgeA, 
                                                nGEdgeB, 
                                                nNewPointsA,
                                                nNewPointsB);
        
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        
        
        nNewPointsA = 0;
        nNewPointsB = 0;
        
        if ((u[0] <= 1.) && (u[0] >= 0.)) {
          newInter->oNewPointsA[nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[nNewPointsA] = u[0];
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[nNewPointsA] = 1 - newInter->uA[nNewPointsA];
          }

          newInter->linkA[nNewPointsA] = nGVtxB[0];
          newInter->gNumA[nNewPointsA] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsA + i] =  closestB1EdgeA[i];
          }
          nNewPointsA += 1;
        }

        if ((u[1] <= 1.) && (u[1] >= 0.)) {
          newInter->oNewPointsA[nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[nNewPointsA] = u[1];
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[nNewPointsA] = 1 - newInter->uA[nNewPointsA];
          }

          newInter->linkA[nNewPointsA] = nGVtxB[1];
          newInter->gNumA[nNewPointsA] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsA + i] =  closestB2EdgeA[i];
          }
          nNewPointsA += 1;
        }
        
        if ((v[0] <= 1.) && (v[0] >= 0.)) {
          newInter->oNewPointsB[nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[nNewPointsB] = v[0];
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uA[nNewPointsB] = 1 - newInter->uA[nNewPointsB];
          }

          newInter->linkB[nNewPointsB] = nGVtxA[0];
          newInter->gNumB[nNewPointsB] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsB + i] =  closestA1EdgeB[i];
          }
          nNewPointsB += 1;
        }

        if ((v[1] <= 1.) && (v[1] >= 0.)) {
          newInter->oNewPointsB[nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[nNewPointsB] = v[1];
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uA[nNewPointsB] = 1 - newInter->uA[nNewPointsB];
          }

          newInter->linkB[nNewPointsB] = nGVtxA[1];
          newInter->gNumB[nNewPointsB] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsB + i] =  closestA2EdgeB[i];
          }
          nNewPointsB += 1;
        }

      }
  
      isSetted = true;

    }
    
  }
  
  /*
   * Storage intersection
   */

  if (!isSetted && (tIntersect == PDM_LINE_INTERSECT_YES)) {

    int nNewPointsA = 1;
    int nNewPointsB = 1;

    newInter = _edges_intersect_res_create (nGEdgeA, 
                                            nGEdgeB, 
                                            nNewPointsA,
                                            nNewPointsB);

    newInter->tIntersect = tIntersect;
    newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
    newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);

    newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_NEW;
    newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_NEW;

    newInter->linkA[0] = -1;
    newInter->linkB[0] = -1;
    newInter->gNumA[0] = 0;
    newInter->gNumB[0] = 0;

    newInter->uA[0] = u1;
    newInter->uB[0] = v1;

    if (newInter->originEdgeA == nGVtxA[1]) {
      newInter->uA[0] = 1 - newInter->uA[0];
    }

    if (newInter->originEdgeB == nGVtxB[1]) {
      newInter->uB[0] = 1 - newInter->uB[0];
    }

    for (int i = 0; i < 3; i++) {
      newInter->coordsA[i] =  coordsVtxA[0] + u1 * vA[0];
    }

    for (int i = 0; i < 3; i++) {
      newInter->coordsB[i] =  coordsVtxB[0] + u1 * vB[0];
    }
        
    isSetted = true;
  }

  PDM_hash_tab_data_add (ht, (void *) &key, newInter);
  PDM_hash_tab_data_add (htA, (void *) &keyA, newInter);
  PDM_hash_tab_data_add (htB, (void *) &keyB, newInter);
  
  return (PDM_edges_intersect_res_t *) newInter;
}


/**
 *
 * \brief Free \ref PDM_edges_intersect_t object
 *
 * \param [in]   ei               Current edges intersection pointer
 *
 * \return     NULL     
 */

PDM_edges_intersect_t *
PDM_edges_intersect_free 
(
PDM_edges_intersect_t *ei        
)
{
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;
  
  int *keyMax;
  
  keyMax = (int *) PDM_hash_tab_keyMax_get (_ei->ht);

  for (int i = 0; i < *keyMax; i++) {
    int nData = PDM_hash_tab_n_data_get (_ei->ht, &i);
    _edges_intersect_res_t **eir = 
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (_ei->ht, &i);

    for (int j = 0; j < nData; j++) {
      _edges_intersect_res_free (eir[j]);
    }
    
    free (eir);
  }

  PDM_hash_tab_free (_ei->ht);
  PDM_hash_tab_free (_ei->htA);
  PDM_hash_tab_free (_ei->htB);
  
  free (_ei);
  
  return NULL;
}


/**
 *
 * \brief Get data intersection
 *
 * \param [in]  eir             Current edges intersection result pointer
 * \param [in]  mesh            Origin mesh \ref PDM_edges_intersect_MESHA or
 *                                          \ref PDM_edges_intersect_MESHB
 * \param [out] nGEdge          Global number of meshA edge
 * \param [out] originEdge      Global number of edge origin
 * \param [out] tIntersect      Intersection type
 * \param [out] nNewPoints      Number of intersection points
 * \param [out] oNewPoints      Origin of intersection points                               
 * \param [out] link            Linked vertex in linked mesh
 * \param [out] gNum            Global number in overlay mesh  
 * \param [out] coords          Coordinates of intersection point  
 * \param [out] u               Parameter of the intersections in edges 
 *                              parametric coordinates
 *
 */

void
PDM_edges_intersect_res_data_get 
(
PDM_edges_intersect_res_t   *eir,        
PDM_edges_intersect_mesh_t   mesh,
PDM_g_num_t                  *nGEdge,
PDM_g_num_t                  *originEdge, 
PDM_line_intersect_t        *tIntersect,
int                         *nNewPoints,
PDM_edges_intersect_point_t **oNewPoints,                               
PDM_g_num_t                  **link,
PDM_g_num_t                  **gNum,
double                      **coords, 
double                      **u 
)
{
  _edges_intersect_res_t *_eir = (_edges_intersect_res_t *) eir;
  if (mesh == PDM_EDGES_INTERSECT_MESHA) {
    *nGEdge          = _eir->nGEdgeA;
    *originEdge      = _eir->originEdgeA;
    *tIntersect      = _eir->tIntersect;
    *nNewPoints      = _eir->nNewPointsA;
    *oNewPoints      = _eir->oNewPointsA;                               
    *link            = _eir->linkA;
    *gNum            = _eir->gNumA;
    *coords          = _eir->coordsA;
    *u               = _eir->uA; 
  }
  if (mesh == PDM_EDGES_INTERSECT_MESHB) {
    *nGEdge          = _eir->nGEdgeB;
    *originEdge      = _eir->originEdgeB;
    *tIntersect      = _eir->tIntersect;
    *nNewPoints      = _eir->nNewPointsB;
    *oNewPoints      = _eir->oNewPointsB;                               
    *link            = _eir->linkB;
    *gNum            = _eir->gNumB;
    *coords          = _eir->coordsB;
    *u               = _eir->uB; 
  }
}


/**
 *
 * \brief Perform edges intersection from two polygons
 *
 * \param [in]    intersect            Edges intersection management
 * \param [in]    nVtxA                Number of polygon A vertices
 * \param [in]    faceToEdgeA          Polygon A face to edge connectivity
 * \param [in]    faceToVtxA           Polygon A face to vertex connectivity
 * \param [in]    faceVtxCooA          Polygon A vertex coordinates
 * \param [in]    faceVtxEpsA          Polygon A vertex characteristic length 
 * \param [in]    gNumB                Polygon B global number
 * \param [in]    nVtxB                Number of polygon B vertices
 * \param [in]    faceToEdgeB          Polygon B face to edge connectivity
 * \param [in]    faceToVtxB           Polygon B face to vertex connectivity
 * \param [in]    faceVtxCooB          Polygon B vertex coordinates
 * \param [in]    faceVtxEpsB          Polygon B vertex characteristic length 
 *
 */

void
PDM_edges_intersect_poly_add 
(
PDM_edges_intersect_t  *ei,
const int               nVtxA,
PDM_g_num_t             *faceToEdgeA,
PDM_g_num_t             *faceToVtxA,
double                 *faceVtxCooA,
double                 *faceVtxEpsA,
const int               nVtxB,
PDM_g_num_t             *faceToEdgeB,
PDM_g_num_t             *faceToVtxB,
double                 *faceVtxCooB,
double                 *faceVtxEpsB
)
{
  PDM_g_num_t *_faceToEdgeA = faceToEdgeA; 
  PDM_g_num_t *_faceToVtxA  = faceToVtxA; 
  double     *_faceVtxCooA = faceVtxCooA;
  double     *_faceVtxEpsA = faceVtxEpsA;
  
  PDM_g_num_t *_faceToEdgeB = faceToEdgeB; 
  PDM_g_num_t *_faceToVtxB  = faceToVtxB; 
  double     *_faceVtxCooB = faceVtxCooB;
  double     *_faceVtxEpsB = faceVtxEpsB;

  /*
   * Compute Normal
   */
  
  double nA[3];
  PDM_plane_normal (nVtxA, faceVtxCooA, nA);

  double nB[3];
  PDM_plane_normal (nVtxB, faceVtxCooB, nB);

  double dot = PDM_DOT_PRODUCT (nA, nB);
  
  bool revert = false;
  
  if (dot < 0) {
    revert = true;
  }
    
  /*
   * Reorient if necessary
   */

  if (revert) {
  
    _faceToEdgeB = malloc (sizeof(PDM_g_num_t) * nVtxB); 
    _faceToVtxB  = malloc (sizeof(PDM_g_num_t) * nVtxB); 
    _faceVtxCooB = malloc (sizeof(double) * 3 * nVtxB);
    _faceVtxEpsB = malloc (sizeof(double) * nVtxB);
    
    int j = nVtxB - 1;
    for (int i = 0; i < nVtxB; i++) {
      _faceToEdgeB[i] = faceToEdgeB[j];
      _faceToVtxB[i] = faceToVtxB[j];
      for (int k = 0; k < 3; k++) {
        _faceVtxCooB[3*i+k] = faceVtxCooB[3*j+k];
      }
      _faceVtxEpsB[i] = faceVtxEpsB[j];
      j += -1;
    }
    
  }
  
  /*
   *   Compute Edges Intersection :
   *   - First step : Remove case with vertex located on 2 two different edges
   *   - Second step : Compute other intersection
   */
  

  // 
  // FIXME: faire un test pour verifier que le resultat de 2 intersections d'une 
  // arete de A de B donne toujours des resultats differents
  // Faire l'inverse : Detection  des aretes a pb .
  // Echange MPI des aretes a pb broad cast apres un premier clipping
  // Mise a jour des clipping concernes
  //

  int *vtxAOnEdgeB = malloc(sizeof(int) * nVtxA);
  _edges_intersect_res_t **vtxAOnEdgeBEir = malloc(sizeof(_edges_intersect_res_t *) * nVtxA);
  for (int i = 0; i < nVtxA; i++) {
    vtxAOnEdgeBEir[i] = NULL;
  }
  
  int *vtxBOnEdgeA = malloc(sizeof(int) * nVtxB);
  _edges_intersect_res_t **vtxBOnEdgeAEir = malloc(sizeof(_edges_intersect_res_t *) * nVtxB);
  for (int i = 0; i < nVtxB; i++) {
    vtxBOnEdgeAEir[i] = NULL;
  }
  
  for (int i = 0; i < nVtxA; i++) {

    int inext = (i + 1) % nVtxA;
    
    PDM_g_num_t nGVtxA[2]   = {_faceToVtxA[i], _faceToVtxA[inext]};
    double charLgthVtxA[2] = {_faceVtxEpsA[i], _faceVtxEpsA[inext]};
    double coordsVtxA[6]   = {_faceVtxCooA[3*i],
                              _faceVtxCooA[3*i+1],
                              _faceVtxCooA[3*i+2],
                              _faceVtxCooA[3*inext],
                              _faceVtxCooA[3*inext+1],
                              _faceVtxCooA[3*inext+2]};

    for (int j = 0; j < nVtxB; j++) {
    
      int jnext = (j + 1) % nVtxB;

      PDM_g_num_t nGVtxB[2]   = {_faceToVtxB[j], _faceToVtxB[jnext]};
      double charLgthVtxB[2] = {_faceVtxEpsB[j], _faceVtxEpsB[jnext]};
      double coordsVtxB[6]   = {_faceVtxCooB[3*j],
                                _faceVtxCooB[3*j+1],
                                _faceVtxCooB[3*j+2],
                                _faceVtxCooB[3*jnext],
                                _faceVtxCooB[3*jnext+1],
                                _faceVtxCooB[3*jnext+2]};

      /*
       * Perform intersection 
       */

      PDM_edges_intersect_res_t *eir = 
                PDM_edges_intersect_add (ei,
                                         _faceToEdgeA[i],
                                         nGVtxA,
                                         charLgthVtxA,
                                         coordsVtxA,
                                         _faceToEdgeB[j],
                                         nGVtxB,
                                         charLgthVtxB,
                                         coordsVtxB);
      /*
       * Get intersection properties
       */
      
      PDM_line_intersect_t         tIntersect;

      PDM_g_num_t                   nGEdgeA;
      PDM_g_num_t                   originEdgeA;
      int                          nNewPointsA;
      PDM_edges_intersect_point_t *oNewPointsA;                               
      PDM_g_num_t                  *linkA;
      PDM_g_num_t                  *gNumA;
      double                      *coordsA; 
      double                      *uA; 


      PDM_edges_intersect_res_data_get (eir,        
                                        PDM_EDGES_INTERSECT_MESHA,
                                        &nGEdgeA,
                                        &originEdgeA,
                                        &tIntersect,
                                        &nNewPointsA,
                                        &oNewPointsA,                               
                                        &linkA,
                                        &gNumA,
                                        &coordsA,
                                        &uA);

      PDM_g_num_t                   nGEdgeB;
      PDM_g_num_t                   originEdgeB;
      int                          nNewPointsB;
      PDM_edges_intersect_point_t *oNewPointsB;                               
      PDM_g_num_t                  *linkB;
      PDM_g_num_t                  *gNumB;
      double                      *coordsB; 
      double                      *uB; 

      PDM_edges_intersect_res_data_get (eir,        
                                        PDM_EDGES_INTERSECT_MESHB,
                                        &nGEdgeB,
                                        &originEdgeB,
                                        &tIntersect,
                                        &nNewPointsB,
                                        &oNewPointsB,                               
                                        &linkB,
                                        &gNumB,
                                        &coordsB,
                                        &uB);
      
      _edges_intersect_res_t *_eir = (_edges_intersect_res_t *) eir;

      /*
       * Remove inconsistencies : 
       * Check if a vertex is not on two different edges
       *   - case 1 : B vertex on two 2 different A edges 
       *   - case 2 : A vertex on two 2 different B edges 
       */
      
      for (int k = 0; k < nNewPointsA; k++) {
        if (oNewPointsA[k] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
          int ind = j;
          if (linkA[k] != nGVtxB[0]) {
            ind = jnext;
            assert(linkA[k] != nGVtxB[1]);
          }
          
          /*
           * If B vertex is already on an other A edge :
           * Update intersections to merge vertex case 
           */
          
          if (vtxBOnEdgeA[ind] != -1) {
            int i1 = vtxBOnEdgeA[ind];
            int i1next = (i1 + 1) % nVtxA;

            /*
             * Look for common vertex
             */
            
            PDM_g_num_t prenGVtxA[2] = {_faceToVtxA[i1], _faceToVtxA[i1next]};

            int common_isom = -1;

            for (int k1 = 0; k1 < 2; k1++) {
              for (int k2 = 0; k2 < 2; k2++) {
                if (nGVtxA[k1] == prenGVtxA[k2]) {
                  common_isom = k1;
                  break;
                }
              }
            }

            if (common_isom == -1) {
              PDM_error(__FILE__, __LINE__, 0, "Probleme simplication sommet proche de 2 aretes :\n"
                      "les 2 aretes n'ont pas de sommet commun\n");
              abort();
            }
            
            double coords_new[3];
            for (int k1 = 0; k1 < 3; k1++) {
              coords_new[k1] = (_faceVtxCooB[3*ind+k1] + coordsVtxA[3*common_isom+k1]) / 2;
            }
              
            /*
             * Update preEir and eir to merge vertex case
             *   - From A : switch to PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB
             *   - From B : Add a new modified A vertex
             *   - Move vertices to half-distance. 
             */

            /* From A */

            
            _edges_intersect_res_t *preEir = vtxBOnEdgeAEir[ind];
            
            for (int k1 = 0; k1 < preEir->nNewPointsA; k1++) {
              if ((preEir->oNewPointsA[k1] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) &&
                  (preEir->linkA[k1] == linkA[k])) {    
                preEir->oNewPointsA[k1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                preEir->linkA[k1]     = preEir->linkA[k1];
                preEir->gNumA[k1]     = preEir->gNumA[k1];
                if (preEir->originEdgeA == nGVtxA[common_isom]) {
                  preEir->uA[k1] = 0;
                }
                else {
                  preEir->uA[k1] = 1.;                  
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  preEir->coordsA[3*k1+k2] = coords_new[k2];
                }
              }
            }
            
            _eir->oNewPointsA[k] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            _eir->linkA[k]     = _eir->linkA[k];
            _eir->gNumA[k]     = 0;
            if (_eir->originEdgeA == nGVtxA[common_isom]) {
              _eir->uA[k] = 0;
            }
            else {
              _eir->uA[k] = 1.;                  
            }
            for (int k2 = 0; k2 < 3; k2++) {
              _eir->coordsA[3*k+k2] = coords_new[k2];
            }
   
            /* From B */

            preEir->oNewPointsB = realloc(preEir->oNewPointsB, 
                                          (preEir->nNewPointsB + 1) * sizeof(PDM_edges_intersect_point_t));
            preEir->linkB = realloc(preEir->linkB, sizeof(PDM_g_num_t) * (preEir->nNewPointsB + 1));
            preEir->gNumB = realloc(preEir->gNumB, sizeof(PDM_g_num_t) * (preEir->nNewPointsB + 1));
            preEir->uB  = realloc(preEir->uB, sizeof(double) * (preEir->nNewPointsB + 1));
            
            preEir->oNewPointsB[preEir->nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            preEir->linkB[preEir->nNewPointsB] = nGVtxA[common_isom];
            preEir->gNumB[preEir->nNewPointsB] = 0;
            if (preEir->originEdgeB == linkA[k]) {
              preEir->uB[preEir->nNewPointsB] = 0;
            }
            else {
              preEir->uB[preEir->nNewPointsB] = 1.;                  
            }
            for (int k2 = 0; k2 < 3; k2++) {
              preEir->coordsB[3*preEir->nNewPointsB+k2] = coords_new[k2];
            }
            preEir->nNewPointsB += 1 ;
            

            _eir->oNewPointsB = realloc(_eir->oNewPointsB, 
                                          (_eir->nNewPointsB + 1) * sizeof(PDM_edges_intersect_point_t));
            _eir->linkB = realloc(_eir->linkB, sizeof(PDM_g_num_t) * (_eir->nNewPointsB + 1));
            _eir->gNumB = realloc(_eir->gNumB, sizeof(PDM_g_num_t) * (_eir->nNewPointsB + 1));
            _eir->uB  = realloc(_eir->uB, sizeof(double) * (_eir->nNewPointsB + 1));
            
            _eir->oNewPointsB[_eir->nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            _eir->linkB[_eir->nNewPointsB] = nGVtxA[common_isom];
            _eir->gNumB[_eir->nNewPointsB] = 0;
            if (_eir->originEdgeB == linkA[k]) {
              _eir->uB[_eir->nNewPointsB] = 0;
            }
            else {
              _eir->uB[_eir->nNewPointsB] = 1.;                  
            }
            for (int k2 = 0; k2 < 3; k2++) {
              _eir->coordsB[3*_eir->nNewPointsB+k2] = coords_new[k2];
            }
            _eir->nNewPointsB += 1 ;
            
          }

          /*
           * If point B is not on an other edge : tag it
           */

          else {
            vtxBOnEdgeA[ind]    = i;
            vtxBOnEdgeAEir[ind] = _eir;
          }            
        }       
      }

      for (int k = 0; k < nNewPointsB; k++) {
        if (oNewPointsB[k] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {
          int ind = i;
          if (linkB[k] != nGVtxA[0]) {
            ind = inext;
            assert(linkB[k] != nGVtxA[1]);
          }
          
          /*
           * If B vertex is already on an other A edge :
           * Update intersections to merge vertex case 
           */
          
          if (vtxAOnEdgeB[ind] != -1) {
            int j1 = vtxAOnEdgeB[ind];
            int j1next = (j1 + 1) % nVtxB;

            /*
             * Look for common vertex
             */
            
            PDM_g_num_t prenGVtxB[2] = {_faceToVtxB[j1], _faceToVtxB[j1next]};

            int common_isom = -1;

            for (int k1 = 0; k1 < 2; k1++) {
              for (int k2 = 0; k2 < 2; k2++) {
                if (nGVtxB[k1] == prenGVtxB[k2]) {
                  common_isom = k1;
                  break;
                }
              }
            }

            if (common_isom == -1) {
              PDM_error(__FILE__, __LINE__, 0, "Probleme simplication sommet proche de 2 aretes :\n"
                      "les 2 aretes n'ont pas de sommet commun\n");
              abort();
            }

            double coords_new[3];
            for (int k1 = 0; k1 < 3; k1++) {
              coords_new[k1] = (_faceVtxCooA[3*ind+k1] + coordsVtxB[3*common_isom+k1]) / 2;
            }

            /*
             * Update preEir and eir to merge vertex case
             *   - From B : switch to PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB
             *   - From A : Add a new modified B vertex
             *   - Move vertices to half-distance. 
             */

            /* From B */
            
            _edges_intersect_res_t *preEir = vtxAOnEdgeBEir[ind];
            
            for (int k1 = 0; k1 < preEir->nNewPointsB; k1++) {
              if ((preEir->oNewPointsB[k1] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) &&
                  (preEir->linkB[k1] == linkB[k])) {    
                preEir->oNewPointsB[k1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                preEir->linkB[k1]     = preEir->linkB[k1];
                preEir->gNumB[k1]     = preEir->gNumB[k1];
                if (preEir->originEdgeB == nGVtxB[common_isom]) {
                  preEir->uB[k1] = 0;
                }
                else {
                  preEir->uB[k1] = 1.;                  
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  preEir->coordsB[3*k1+k2] = coords_new[k2];
                }
              }
            }
            
            _eir->oNewPointsB[k] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            _eir->linkB[k]     = _eir->linkB[k];
            _eir->gNumB[k]     = _eir->gNumB[k];
            if (_eir->originEdgeB == nGVtxB[common_isom]) {
              _eir->uB[k] = 0;
            }
            else {
              _eir->uB[k] = 1.;                  
            }
            for (int k2 = 0; k2 < 3; k2++) {
              preEir->coordsB[3*k+k2] = coords_new[k2];
            }

            /* From A */

            preEir->oNewPointsA = realloc(preEir->oNewPointsA, 
                                          (preEir->nNewPointsA + 1) * sizeof(PDM_edges_intersect_point_t));
            preEir->linkA = realloc(preEir->linkA, sizeof(PDM_g_num_t) * (preEir->nNewPointsA + 1));
            preEir->gNumA = realloc(preEir->gNumA, sizeof(PDM_g_num_t) * (preEir->nNewPointsA + 1));
            preEir->uA  = realloc(preEir->uA, sizeof(double) * (preEir->nNewPointsA + 1));
            
            preEir->oNewPointsA[preEir->nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            preEir->linkA[preEir->nNewPointsA] = nGVtxB[common_isom];
            preEir->gNumA[preEir->nNewPointsA] = 0;
            if (preEir->originEdgeA == linkB[k]) {
              preEir->uA[preEir->nNewPointsA] = 0;
            }
            else {
              preEir->uA[preEir->nNewPointsA] = 1.;                  
            }
            preEir->nNewPointsA += 1 ;

            _eir->oNewPointsA = realloc(_eir->oNewPointsA, 
                                          (_eir->nNewPointsA + 1) * sizeof(PDM_edges_intersect_point_t));
            _eir->linkA = realloc(_eir->linkA, sizeof(PDM_g_num_t) * (_eir->nNewPointsA + 1));
            _eir->gNumA = realloc(_eir->gNumA, sizeof(PDM_g_num_t) * (_eir->nNewPointsA + 1));
            _eir->uA  = realloc(_eir->uA, sizeof(double) * (_eir->nNewPointsA + 1));
            
            _eir->oNewPointsA[_eir->nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
            _eir->linkA[_eir->nNewPointsA] = nGVtxB[common_isom];
            _eir->gNumA[_eir->nNewPointsA] = 0;
            if (_eir->originEdgeA == linkB[k]) {
              _eir->uA[_eir->nNewPointsA] = 0;
            }
            else {
              _eir->uA[_eir->nNewPointsA] = 1.;                  
            }
            _eir->nNewPointsA += 1 ;
            
          }

          /*
           * If point B is not on an other edge : tag it
           */

          else {
            vtxAOnEdgeB[ind]    = j;
            vtxAOnEdgeBEir[ind] = _eir;
          }            
        }       
      } 
    }
  }

  /*
   * Free local memory
   */
    
  if (_faceToEdgeA != faceToEdgeA) {
    free (_faceToEdgeA);
  }
  
  if (_faceToVtxA  == faceToVtxA) {
    free (_faceToVtxA);
  }
    
  if (_faceVtxCooA == faceVtxCooA) {
    free (_faceVtxCooA);
  }
    
  if (_faceVtxEpsA == faceVtxEpsA) {
    free (_faceVtxEpsA);
  }

  if (_faceToEdgeB != faceToEdgeB) {
    free (_faceToEdgeB);
  }
  
  if (_faceToVtxB == faceToVtxB) {
    free (_faceToVtxB);
  }
    
  if (_faceVtxCooB == faceVtxCooB) {
    free (_faceVtxCooB);
  }
    
  if (_faceVtxEpsB == faceVtxEpsB) {
    free (_faceVtxEpsB);
  }
}


/**
 *
 * \brief Remove inconsistencies between processes
 *
 * \param [in]   ei           Current edges intersection pointer
 * \param [in]   nAbsVtxA     Absolute number of vertices in initial A mesh
 * \param [in]   nAbsVtxB     Absolute number of vertices in initial B mesh
 * \param [out]  nAbsNewVtxA  Absolute number of vertices in A mesh after intersections
 * \param [out]  nAbsNewVtxB  Absolute number of vertices in B mesh after intersections
 *
 */

void
PDM_edges_intersect_synchronize 
(
PDM_edges_intersect_t       *ei,
PDM_g_num_t             nAbsVtxA,
PDM_g_num_t             nAbsVtxB,
PDM_g_num_t            *nAbsNewVtxA,
PDM_g_num_t            *nAbsNewVtxB
)
{
  
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;

  int sMSGComm;
  PDM_MPI_Comm_size (_ei->comm, &sMSGComm);
  int lastRank = sMSGComm - 1; 

  int myRank;
  PDM_MPI_Comm_rank (_ei->comm, &myRank);
  
  /*
   * Part to block hash table
   */
  
  PDM_hash_tab_t *ht = _ei->ht;
  
  int keyMax  = * ((int *) PDM_hash_tab_keyMax_get (ht));
  int *nDataKey = malloc (sizeof(int) * keyMax);

  for (int key = 0; key < keyMax; key++) {
    nDataKey[key] = 0;
  }
  
  int nProcData = 0;
  for (int key = 0; key < keyMax; key++) {
    
    int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    nDataKey[key] = nData;
    nProcData += nData;

/*

    _edges_intersect_res_t ** datas = 
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                               (void *) &key);   
*/
    
  }
    
  PDM_g_num_t *keys        = malloc (sizeof(PDM_g_num_t) * nProcData);

  int        *tIntersects = malloc (sizeof(int) * nProcData); 

  PDM_g_num_t *gNumEdgeA   = malloc (sizeof(PDM_g_num_t) * nProcData);
  int        *nNewPointsA = malloc (sizeof(int) * nProcData); 
  PDM_edges_intersect_point_t *oNewPointsA = 
          malloc (sizeof(PDM_edges_intersect_point_t) * 2 * nProcData);
  PDM_g_num_t *connectPointA = malloc (sizeof(PDM_g_num_t) * 2 * nProcData);
  double *uPointA = malloc (sizeof(double) * 2 * nProcData);
  double *coordsPointA = malloc (sizeof(double) * 6 * nProcData);

  int        *nNewPointsB = malloc (sizeof(int) * nProcData); 
  PDM_edges_intersect_point_t *oNewPointsB = 
          malloc (sizeof(PDM_edges_intersect_point_t) * 2 * nProcData);
  PDM_g_num_t *connectPointB = malloc (sizeof(PDM_g_num_t) * 2 * nProcData);
  double *uPointB = malloc (sizeof(double) * 2 * nProcData);
  double *coordsPointB = malloc (sizeof(double) * 6 * nProcData);
  
  nProcData = 0;
  for (int key = 0; key < keyMax; key++) {
  
    int idxA = 0;
    int idxB = 0;

    _edges_intersect_res_t ** datas = 
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                               (void *) &key);   

    int nData = PDM_hash_tab_n_data_get (ht, &key);
    
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *_data = datas[i];
      
      keys[nProcData]        = _data->nGEdgeA + _data->nGEdgeB;
      tIntersects[nProcData] = _data->tIntersect;

      gNumEdgeA[nProcData  ] = _data->nGEdgeA;
      nNewPointsA[nProcData] = _data->nNewPointsA;
      for (int j = 0; j < nNewPointsA[nProcData]; j++) {
        oNewPointsA[idxA] = _data->oNewPointsA[j];
        
        PDM_g_num_t gNum = _data->linkA[j];
        for (int k = 0; k < 3; k++) { 
          coordsPointA[3 * idxA + k] = _data->coordsA[3*j+k];
        }
        connectPointA[idxA] = gNum;
        uPointA[idxA] = _data->uA[j];
        idxA += 1;
      }
     
      nNewPointsB[nProcData] = _data->nNewPointsB;
      for (int j = 0; j < nNewPointsB[nProcData]; j++) {
        oNewPointsB[idxB] = _data->oNewPointsB[j];
        
        PDM_g_num_t gNum = _data->linkB[j];
        for (int k = 0; k < 3; k++) { 
          coordsPointB[3 * idxB + k] = _data->coordsB[3*j+k];
        }
        connectPointB[idxB] = gNum;
        uPointB[idxB] = _data->uB[j];
        idxB += 1;
      }

      nProcData++;
    }
  }
  
  /*
   * TODO : - A single exchange to optimize communications !
   *        - Add a revert param to PDM_part_to_block_exch :
   *          To perform block to part communication 
   *          (less expensive than PDM_block_to_part))
   */
  
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCK_POST_MERGE,
                                                     1.,
                                                     (PDM_g_num_t **) &keys,
                                                     &nProcData,
                                                     1,
                                                     _ei->comm);
 
  PDM_g_num_t *block_gnum = (PDM_g_num_t *) PDM_part_to_block_block_gnum_get (ptb);
  int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);

  /*
   * A info
   */ 

  int *strideOne   = malloc (sizeof(PDM_g_num_t) * nProcData);
  for (int i = 0; i < nProcData; i++) {
    strideOne[i] = 1;
  }

  int *b_tIntersects = NULL;
  int *b_strideOne = NULL;
  
  PDM_part_to_block_exch (ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR,
                         1,
                         &strideOne,
                         (void **) &tIntersects,
                         &b_strideOne,
                         (void **) &b_tIntersects);
  
  free (b_strideOne);

  PDM_g_num_t *b_gNumEdgeA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         1,
                         &strideOne,
                         (void **) &gNumEdgeA,
                         &b_strideOne,
                         (void **) &b_gNumEdgeA);
 

  PDM_edges_intersect_point_t *b_oNewPointsA  = NULL;
  int                         *b_nNewPointsA;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsA,
                         (void **)&oNewPointsA,
                         &b_nNewPointsA,
                         (void **)&b_oNewPointsA);
  
  free(b_nNewPointsA);

  PDM_g_num_t *b_connectPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsA,
                         (void **)&connectPointA,
                         &b_nNewPointsA,
                         (void **)&b_connectPointA);

  free(b_nNewPointsA);


  for (int k = 0; k < nProcData; k++) {
    nNewPointsA[k] = 2 * nNewPointsA[k]; 
  }

  double *b_uPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsA,
                         (void **)&uPointA,
                         &b_nNewPointsA,
                         (void **)&b_uPointA);

  free(b_nNewPointsA);

  for (int k = 0; k < nProcData; k++) {
    nNewPointsA[k] = 3 * nNewPointsA[k]/2; 
  }
  
  double *b_coordsPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsA,
                         (void **)&coordsPointA,
                         &b_nNewPointsA,
                         (void **)&b_coordsPointA);

  for (int k = 0; k < nProcData; k++) {
    nNewPointsA[k] = nNewPointsA[k] / 3;
  }

  for (int k = 0; k < n_elt_block; k++) {
    b_nNewPointsA[k] = b_nNewPointsA[k] / 3;
  }
  
  /*
   * B info
   */ 
  
  PDM_edges_intersect_point_t *b_oNewPointsB  = NULL;
  int                         *b_nNewPointsB;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsB,
                         (void **)&oNewPointsB,
                         &b_nNewPointsB,
                         (void **)&b_oNewPointsB);
  
  free(b_nNewPointsB);

  PDM_g_num_t *b_connectPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsB,
                         (void **)&connectPointB,
                         &b_nNewPointsB,
                         (void **)&b_connectPointB);

  free(b_nNewPointsB);

  for (int k = 0; k < nProcData; k++) {
    nNewPointsB[k] = 2 * nNewPointsB[k]; 
  }

  double *b_uPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsB,
                         (void **)&uPointB,
                         &b_nNewPointsB,
                         (void **)&b_uPointB);

  free(b_nNewPointsB);


  for (int k = 0; k < nProcData; k++) {
    nNewPointsB[k] = 3 * nNewPointsB[k]/2; 
  }
  
  double *b_coordsPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &nNewPointsB,
                         (void **)&coordsPointB,
                         &b_nNewPointsB,
                         (void **)&b_coordsPointB);

  for (int k = 0; k < nProcData; k++) {
    nNewPointsB[k] = nNewPointsB[k] / 3;
  }

  for (int k = 0; k < n_elt_block; k++) {
    b_nNewPointsB[k] = b_nNewPointsB[k] / 3;
  }
    
  free (tIntersects); 

  free (gNumEdgeA);
  
  free (nNewPointsA); 
  free (oNewPointsA);
  free (connectPointA);
  free (coordsPointA);
  free (uPointA);

  free (nNewPointsB); 
  free (oNewPointsB);
  free (connectPointB);
  free (coordsPointB);
  free (uPointB);

  /*
   * Remove inconsistencies
   *   - Fill true intersection properties
   */

  int *b_strideOne_idx = malloc (sizeof(int) * (n_elt_block + 1));
  b_strideOne_idx[0] = 0;
  for (int i = 0; i < n_elt_block; i++) {
    b_strideOne_idx[i+1] = b_strideOne_idx[i] + b_strideOne[i];  
  }
  
  free (b_strideOne);
  
  int *tag = malloc(sizeof(int) * b_strideOne_idx[n_elt_block]);
  for (int i = 0; i < b_strideOne_idx[n_elt_block]; i++) {
    tag[i] = 0;
  }

  int *b_nNewPointsA_idx = malloc(sizeof(int) * (n_elt_block + 1));
  b_nNewPointsA_idx[0] = 0;
  for (int i = 0; i < n_elt_block; i++) {
    b_nNewPointsA_idx[i+1] = b_nNewPointsA_idx[i] + b_nNewPointsA[i]; 
  }
  
  int *b_nNewPointsB_idx = malloc(sizeof(int) * (n_elt_block + 1));
  b_nNewPointsB_idx[0] = 0;
  for (int i = 0; i < n_elt_block; i++) {
    b_nNewPointsB_idx[i+1] = b_nNewPointsB_idx[i] + b_nNewPointsB[i]; 
  }

  int *b_strideOne_idx_true = malloc(sizeof(int) * (n_elt_block + 1));

  int *b_tIntersects_true = malloc(sizeof(int) * b_strideOne_idx[n_elt_block]);
  PDM_g_num_t *b_gNumEdgeA_true = malloc(sizeof(PDM_g_num_t) * b_strideOne_idx[n_elt_block]);
  
  int *b_nNewPointsA_true = malloc(sizeof(int) * b_strideOne_idx[n_elt_block]);
  PDM_edges_intersect_point_t *b_oNewPointsA_true = 
          malloc(sizeof(PDM_edges_intersect_point_t) * 2 * b_strideOne_idx[n_elt_block]);
  PDM_g_num_t *b_connectPointA_true = 
          malloc(sizeof(PDM_g_num_t) * 2 * b_strideOne_idx[n_elt_block]);
  double *b_uPointA_true = 
          malloc(sizeof(double) * 2 * b_strideOne_idx[n_elt_block]);
  double *b_coordsPointA_true = 
          malloc(sizeof(double) * 6 * b_strideOne_idx[n_elt_block]);
  
  int *b_nNewPointsB_true = malloc(sizeof(int) * b_strideOne_idx[n_elt_block]);
  PDM_edges_intersect_point_t *b_oNewPointsB_true = 
          malloc(sizeof(PDM_edges_intersect_point_t) * 2 * b_strideOne_idx[n_elt_block]);
  PDM_g_num_t *b_connectPointB_true = 
          malloc(sizeof(PDM_g_num_t) * 2 * b_strideOne_idx[n_elt_block]);
  double *b_uPointB_true = malloc(sizeof(double) * 2 * b_strideOne_idx[n_elt_block]);
  double *b_coordsPointB_true = malloc(sizeof(double) * 6 * b_strideOne_idx[n_elt_block]);
    
  int idx_true   = 0;
  int idx_newPtA = 0;
  int idx_newPtB = 0;

  for (int i = 0; i < n_elt_block; i++) {
    b_strideOne_idx_true[i] = idx_true; 
    PDM_g_num_t gNum = block_gnum[i];
    for (int j = b_strideOne_idx[i]; j < b_strideOne_idx[i+1]; j++) {
      if (tag[j] == 0) {
        tag[j] = 1;
        
        for (int k = j+1; k < b_strideOne_idx[i+1]; k++) {
          if ((tag[k] == 0) && (b_gNumEdgeA[j] == b_gNumEdgeA[k])) {
            tag[j] += 1;
            if (tag[j] > 2) {
              PDM_error(__FILE__, __LINE__, 0, "Sortie en erreur : intersection calculée + de 2 fois\n");
              abort();
            }
            tag[k] = -1;
            
            b_tIntersects_true[idx_true] = PDM_MAX (b_tIntersects[j], b_tIntersects[k]);
            b_gNumEdgeA_true[idx_true] = b_gNumEdgeA[j];            

            /*
             * Merge new points for A
             */
            
            b_nNewPointsA_true[idx_true] = PDM_MAX (b_nNewPointsA[j], b_nNewPointsA[k]);
            if ((b_nNewPointsA[k] == 2) && (b_nNewPointsA[j] == 2)) {
              int link[2] = {-1, -1};
              for (int k1 = 0; k1 < b_nNewPointsA[j]; k1++) {
                for (int k2 = 0; k2 < b_nNewPointsA[k]; k2++) {
                  if (b_connectPointA[b_nNewPointsA_idx[j]+k1] ==
                      b_connectPointA[b_nNewPointsA_idx[k]+k2]) {    
                    link[k1] = k2;
                    break;
                  }
                }
              }
              
              assert ((link[0] != -1) && (link[1] != -1));
              b_oNewPointsA_true[idx_newPtA] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j]], 
                                                        b_oNewPointsA[b_nNewPointsA_idx[k] + link[0]]);
              int i_true1 = b_nNewPointsA_idx[j];
              if (b_oNewPointsA_true[idx_newPtA] == b_oNewPointsA[b_nNewPointsA_idx[k] + link[0]]) {
                i_true1 = b_nNewPointsA_idx[k] + link[0];
              }

              b_oNewPointsA_true[idx_newPtA + 1] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j] + 1], 
                                                            b_oNewPointsA[b_nNewPointsA_idx[k] + link[1]]);
              int i_true2 = b_nNewPointsA_idx[j] + 1;
              if (b_oNewPointsA_true[idx_newPtA + 1] == b_oNewPointsA[b_nNewPointsA_idx[k] + link[1]]) {
                i_true2 = b_nNewPointsA_idx[k] + link[1];
              }
              
              b_connectPointA_true[idx_newPtA    ] = b_connectPointA[i_true1];
              b_connectPointA_true[idx_newPtA + 1] = b_connectPointA[i_true2];
              
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true1 + k1];
                b_coordsPointA_true[3*(idx_newPtA+1) + k1] = b_coordsPointA[3*i_true2 + k1];
              }
              
              for (int k1 = 0; k1 < 2; k1++) {
                b_uPointA_true[2*idx_newPtA + k1] = b_uPointA[2*i_true1 + k1];
                b_uPointA_true[2*(idx_newPtA+1) + k1] = b_uPointA[2*i_true2 + k1];
              }

              idx_newPtA += 2;
              
            }
            
            else if ((b_nNewPointsA[k] == 1) && (b_nNewPointsA[j] == 1)) {
              b_oNewPointsA_true[idx_newPtA] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j]], 
                                                        b_oNewPointsA[b_nNewPointsA_idx[k]]);
              int i_true1 = b_nNewPointsA_idx[j];
              if (b_oNewPointsA_true[idx_newPtA] == b_oNewPointsA[b_nNewPointsA_idx[k]]) {
                i_true1 = b_nNewPointsA_idx[k];
              }

              b_connectPointA_true[idx_newPtA    ] = b_connectPointA[i_true1];
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true1 + k1];
              }
              
              for (int k1 = 0; k1 < 2; k1++) {
                b_uPointA_true[2*idx_newPtA + k1] = b_uPointA[2*i_true1 + k1];
              }
              idx_newPtA += 1;
            }

            else if (b_nNewPointsA[k] != b_nNewPointsA[j]) {
              int idx_sup;
              int idx_inf;
              if (b_nNewPointsA[j] > b_nNewPointsA[k]) {
                idx_sup = j;
                idx_inf = k;
              }
              else {
                idx_sup = k;
                idx_inf = j;
              }

              if (b_nNewPointsA[idx_inf] == 0) {
                
                int i_true = b_nNewPointsA_idx[idx_sup];
                b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true]; 

                b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true];
                for (int k1 = 0; k1 < 3; k1++) {
                  b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true + k1];
                }

                for (int k1 = 0; k1 < 2; k1++) {
                  b_uPointA_true[2*idx_newPtA + k1] = b_uPointA[2*i_true + k1];
                }
                idx_newPtA += 1;
                
              }

              else {
                
                assert (b_nNewPointsA[idx_inf] == 1);
                assert (b_nNewPointsA[idx_sup] == 2);

                /*
                 * if "idx_inf" is "Vertex A on vertex B" : 
                 *    - Look for the related point in "idx_sup"
                 *    - copy "idx_inf"
                 *    - copy the second "idx_sup" vertex
                 * Otherwise :
                 *    - copy "idx_sup" vertices
                 *
                 */
                
                int i_true[2] = {-1, -1};
                
                if (b_oNewPointsA[b_nNewPointsA_idx[idx_inf]] 
                    == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                                   
                  i_true[0] = b_nNewPointsA_idx[idx_inf]; 
                  if (b_connectPointA[b_nNewPointsA_idx[idx_sup]] == gNum) {
                    i_true[1] = b_nNewPointsA_idx[idx_sup];
                  }
                  else if (b_connectPointA[b_nNewPointsA_idx[idx_sup] + 1] == gNum) {
                    i_true[1] = b_nNewPointsA_idx[idx_sup] + 1;
                  }
                  assert (i_true[1] != -1);
                  
                }
                
                else {
                  i_true[0] = b_nNewPointsA_idx[idx_sup]; 
                  i_true[1] = b_nNewPointsA_idx[idx_sup] + 1; 
                }
                
                for (int k2 = 0; k2 < 2; k2++) {
                  b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true[k2]]; 

                  b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true[k2]];
                  for (int k1 = 0; k1 < 3; k1++) {
                    b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3 * i_true[k2] + k1];
                  }

                  for (int k1 = 0; k1 < 2; k1++) {
                    b_uPointA_true[2*idx_newPtA + k1] = b_uPointA[2 * i_true[k2] + k1];
                  }
                  idx_newPtA += 1;
                }
              }
            }

            else {
              b_nNewPointsA_true[idx_true] = 0;
            }

            /*
             * Merge new points for B
             */
            
            b_nNewPointsB_true[idx_true] = PDM_MAX (b_nNewPointsB[j], b_nNewPointsB[k]);
            if ((b_nNewPointsB[k] == 2) && (b_nNewPointsB[j] == 2)) {
              int link[2] = {-1, -1};
              for (int k1 = 0; k1 < b_nNewPointsB[j]; k1++) {
                for (int k2 = 0; k2 < b_nNewPointsB[k]; k2++) {
                  if (b_connectPointB[b_nNewPointsB_idx[j]+k1] ==
                      b_connectPointB[b_nNewPointsB_idx[k]+k2]) {    
                    link[k1] = k2;
                    break;
                  }
                }
              }
              
              assert ((link[0] != -1) && (link[1] != -1));
              b_oNewPointsB_true[idx_newPtB] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j]], 
                                                        b_oNewPointsB[b_nNewPointsB_idx[k] + link[0]]);
              int i_true1 = b_nNewPointsB_idx[j];
              if (b_oNewPointsB_true[idx_newPtB] == b_oNewPointsB[b_nNewPointsB_idx[k] + link[0]]) {
                i_true1 = b_nNewPointsB_idx[k] + link[0];
              }

              b_oNewPointsB_true[idx_newPtB + 1] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j] + 1], 
                                                            b_oNewPointsB[b_nNewPointsB_idx[k] + link[1]]);
              int i_true2 = b_nNewPointsB_idx[j] + 1;
              if (b_oNewPointsB_true[idx_newPtB + 1] == b_oNewPointsB[b_nNewPointsB_idx[k] + link[1]]) {
                i_true2 = b_nNewPointsB_idx[k] + link[1];
              }
              
              b_connectPointB_true[idx_newPtB    ] = b_connectPointB[i_true1];
              b_connectPointB_true[idx_newPtB + 1] = b_connectPointB[i_true2];
              
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true1 + k1];
                b_coordsPointB_true[3*(idx_newPtB+1) + k1] = b_coordsPointB[3*i_true2 + k1];
              }
              
              for (int k1 = 0; k1 < 2; k1++) {
                b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2*i_true1 + k1];
                b_uPointB_true[2*(idx_newPtB+1) + k1] = b_uPointB[2*i_true2 + k1];
              }

              idx_newPtB += 2;
              
            }
            
            else if ((b_nNewPointsB[k] == 1) && (b_nNewPointsB[j] == 1)) {
              b_oNewPointsB_true[idx_newPtB] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j]], 
                                                        b_oNewPointsB[b_nNewPointsB_idx[k]]);
              int i_true1 = b_nNewPointsB_idx[j];
              if (b_oNewPointsB_true[idx_newPtB] == b_oNewPointsB[b_nNewPointsB_idx[k]]) {
                i_true1 = b_nNewPointsB_idx[k];
              }

              b_connectPointB_true[idx_newPtB    ] = b_connectPointB[i_true1];
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true1 + k1];
              }
              
              for (int k1 = 0; k1 < 2; k1++) {
                b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2*i_true1 + k1];
              }
              idx_newPtB += 1;
            }

            else if (b_nNewPointsB[k] != b_nNewPointsB[j]) {
              int idx_sup;
              int idx_inf;
              if (b_nNewPointsB[j] > b_nNewPointsB[k]) {
                idx_sup = j;
                idx_inf = k;
              }
              else {
                idx_sup = k;
                idx_inf = j;
              }

              if (b_nNewPointsB[idx_inf] == 0) {
                
                int i_true = b_nNewPointsB_idx[idx_sup];
                b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true]; 

                b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true];
                for (int k1 = 0; k1 < 3; k1++) {
                  b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true + k1];
                }

                for (int k1 = 0; k1 < 2; k1++) {
                  b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2*i_true + k1];
                }
                idx_newPtB += 1;
                
              }

              else {
                
                assert (b_nNewPointsB[idx_inf] == 1);
                assert (b_nNewPointsB[idx_sup] == 2);

                /*
                 * if "idx_inf" is "Vertex A on vertex B" : 
                 *    - Look for the related point in "idx_sup"
                 *    - copy "idx_inf"
                 *    - copy the second "idx_sup" vertex
                 * Otherwise :
                 *    - copy "idx_sup" vertices
                 *
                 */
                
                int i_true[2] = {-1, -1};
                
                if (b_oNewPointsB[b_nNewPointsB_idx[idx_inf]] 
                    == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                                   
                  i_true[0] = b_nNewPointsB_idx[idx_inf]; 
                  if (b_connectPointB[b_nNewPointsB_idx[idx_sup]] == gNum) {
                    i_true[1] = b_nNewPointsB_idx[idx_sup];
                  }
                  else if (b_connectPointB[b_nNewPointsB_idx[idx_sup] + 1] == gNum) {
                    i_true[1] = b_nNewPointsB_idx[idx_sup] + 1;
                  }
                  assert (i_true[1] != -1);
                  
                }
                
                else {
                  i_true[0] = b_nNewPointsB_idx[idx_sup]; 
                  i_true[1] = b_nNewPointsB_idx[idx_sup] + 1; 
                }
                
                for (int k2 = 0; k2 < 2; k2++) {
                  b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true[k2]]; 

                  b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true[k2]];
                  for (int k1 = 0; k1 < 3; k1++) {
                    b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3 * i_true[k2] + k1];
                  }

                  for (int k1 = 0; k1 < 2; k1++) {
                    b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2 * i_true[k2] + k1];
                  }
                  idx_newPtB += 1;
                }
              }
            }

            else {
              b_nNewPointsB_true[idx_true] = 0;
            }
            
            idx_true += 1;
            
          }
        }
        
        /*
         * Copy intersection if only one computation
         */  
        
        if (tag[j] == 1) {
          
          b_tIntersects_true[idx_true] = b_tIntersects[j];
          b_gNumEdgeA_true[idx_true] = b_gNumEdgeA[j];
          
          b_nNewPointsA_true[idx_true] = nNewPointsA[j];

          for (int k = 0; k < nNewPointsA[j]; k++) {
            int i_true = b_nNewPointsA_idx[j] + k ;
            b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true]; 

            b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true];
            for (int k1 = 0; k1 < 3; k1++) {
              b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true + k1];
            }

            for (int k1 = 0; k1 < 2; k1++) {
              b_uPointA_true[2*idx_newPtA + k1] = b_uPointA[2*i_true + k1];
            }
            idx_newPtA += 1;
          }
          
          b_nNewPointsB_true[idx_true] = nNewPointsB[j];

          for (int k = 0; k < nNewPointsB[j]; k++) {
            int i_true = b_nNewPointsB_idx[j] + k ;
            b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true]; 

            b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true];
            for (int k1 = 0; k1 < 3; k1++) {
              b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true + k1];
            }

            for (int k1 = 0; k1 < 2; k1++) {
              b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2*i_true + k1];
            }
            idx_newPtB += 1;
          }

          idx_true += 1;
          
        }
      }
    }
  }
  
  b_strideOne_idx_true[n_elt_block] = idx_true;
  
  free (tag);
  free (b_nNewPointsB_idx);
  free (b_nNewPointsA_idx);
  
  free (b_connectPointA);
  free (b_connectPointB);
  free (b_coordsPointA);  
  free (b_coordsPointB);
  free (b_gNumEdgeA);
  free (b_nNewPointsA);
  free (b_nNewPointsB);
  free (b_nNewPointsA_idx);
  free (b_nNewPointsB_idx);
  free (b_oNewPointsA);
  free (b_oNewPointsB);
  free (b_strideOne);
  free (b_strideOne_idx);
  free (b_tIntersects);
    
  /*
   * Perform absolute number of new points
   *  - grace a b_oNewPointsA_true et b_oNewPointsA_true
   *  - verifier qu'on a le meme nombre de nouveaux points de chaque côté ! 
   *  - mpi_all_reduce ou equivalent
   *    
   */
  
  PDM_g_num_t nNewPtsA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      nNewPtsA += 1;
    }
  }
  
  PDM_g_num_t beg_nNewPtsA = 0;
  PDM_MPI_Request request1;
  PDM_MPI_Iscan(&nNewPtsA, &beg_nNewPtsA, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);
  
  PDM_g_num_t nNewPtsB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      nNewPtsB += 1;
    }
  }

  PDM_g_num_t beg_nNewPtsB = 0;
  PDM_MPI_Request request2;
  PDM_MPI_Iscan(&nNewPtsB, &beg_nNewPtsB, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request2);
  
  int nPtsFromBForA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromBForA += 1;
    }
  }
  
  PDM_MPI_Wait (&request1);

  PDM_g_num_t *b_gNumPointsA_true = malloc(sizeof(PDM_g_num_t) * idx_newPtA);  
  PDM_g_num_t end_nNewPtsA = beg_nNewPtsA + nAbsVtxA + 1; 
  beg_nNewPtsA += -nNewPtsA + 1 + nAbsVtxA;

  PDM_g_num_t beg_nPtsFromBForA = 0;
  
  PDM_MPI_Igather (&end_nNewPtsA, 1, PDM__PDM_MPI_G_NUM,
              &beg_nPtsFromBForA, 1, PDM__PDM_MPI_G_NUM,
              lastRank,
              _ei->comm,
              &request1);
          
  for (int i = 0; i < idx_newPtA; i++) {
    if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      b_gNumPointsA_true[i] = beg_nNewPtsA + i; 
    }
  }
  
  int nPtsFromAForB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromAForB += 1;
    }
  }

  PDM_MPI_Wait (&request2);
  
  PDM_g_num_t *b_gNumPointsB_true = malloc(sizeof(PDM_g_num_t) * idx_newPtB);  
  PDM_g_num_t end_nNewPtsB = beg_nNewPtsB + nAbsVtxB + 1; 
  beg_nNewPtsB += -nNewPtsB + 1 + nAbsVtxB;

  PDM_g_num_t beg_nPtsFromAForB = 0;
  
  PDM_MPI_Igather (&end_nNewPtsB, 1, PDM__PDM_MPI_G_NUM,
              &beg_nPtsFromAForB, 1, PDM__PDM_MPI_G_NUM,
              lastRank,
              _ei->comm,
              &request2);

  for (int i = 0; i < idx_newPtB; i++) {
    if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      b_gNumPointsB_true[i] = beg_nNewPtsB + i; 
    }
  }
  
  /* 
   * New points from B for A
   *    - Compute absolute number
   *    - Synchronize coordinates  
   */

  PDM_MPI_Wait (&request1);
  
  double *b_cNewPointsA_true_pack = malloc (sizeof(double) * 3 * nPtsFromBForA);
  PDM_edges_intersect_point_t *b_oNewPointsA_true_pack = 
          malloc (sizeof(PDM_edges_intersect_point_t) * nPtsFromBForA);
  PDM_g_num_t *b_lNewPointsA_true_pack = malloc (sizeof(PDM_g_num_t) * nPtsFromBForA);

  nPtsFromBForA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromBForA += 1;
      for (int k = 0; k < 3; k++) {
        b_cNewPointsA_true_pack[3*nPtsFromBForA + k] = b_coordsPointA_true[3*i+k];
      }
      b_oNewPointsA_true_pack[nPtsFromBForA] = b_oNewPointsA_true[i];
      b_lNewPointsA_true_pack[nPtsFromBForA] = b_connectPointA_true[i];
      nPtsFromBForA += 1;
    }
  }

  PDM_part_to_block_t *ptbBForA = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          (PDM_g_num_t **) &b_lNewPointsA_true_pack, // Max des points selectionnes 
                                                          &nPtsFromBForA,
                                                          1,
                                                          _ei->comm);
 
  int n_BForA_gnum = PDM_part_to_block_n_elt_block_get (ptbBForA);
  
  int *b_stride_packA = malloc (sizeof(int) * nPtsFromBForA);
  for (int i = 0; i < nPtsFromBForA; i++) {
    b_stride_packA[i] = 3;
  }
  
  int *b_b_stride_packA = NULL;
  double *b_b_cNewPointsA_true_pack = NULL;
  PDM_part_to_block_exch (ptbBForA,
                        sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &b_stride_packA,
                         (void **)&b_cNewPointsA_true_pack,
                         &b_b_stride_packA,
                         (void **)&b_b_cNewPointsA_true_pack);
  
  free (b_b_stride_packA);
  for (int i = 0; i < nPtsFromBForA; i++) {
    b_stride_packA[i] = 1;
  }

  PDM_edges_intersect_point_t *b_b_oNewPointsA_true_pack = NULL;
  PDM_part_to_block_exch (ptbBForA,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR,
                         1,
                         &b_stride_packA,
                         (void **)&b_oNewPointsA_true_pack,
                         &b_b_stride_packA,
                         (void **)&b_b_oNewPointsA_true_pack);
  
  free (b_stride_packA);
  
  /* Synchronize coordinates */
  
  PDM_g_num_t n_active_BForA_gnum = 0;
  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_stride_packA[i];
    int end = b_b_stride_packA[i+1];  
    int cpt = 0;
    int state = 0;
    PDM_edges_intersect_point_t eip;
    for (int j = beg; j < end; j++) {
      n_active_BForA_gnum += 1;
      if (state == 0) {
        state = 1;
        cpt   = 1;
        eip   = b_b_oNewPointsA_true_pack[j];
      }
      else {
        if (eip != b_b_oNewPointsA_true_pack[j]) {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersection : Inconsistencies : A B point is on a A vertex and on a A edge");
          abort();
        }
        cpt += 1;
        
        for (int k = 0; k < 3; k++) {      
          b_b_cNewPointsA_true_pack[3 * beg + k] += b_b_cNewPointsA_true_pack[3 * j + k];
        }
      }
    }  
    for (int k = 0; k < 3; k++) {      
      b_b_cNewPointsA_true_pack[3 * beg + k] = b_b_cNewPointsA_true_pack[3 * beg + k] / cpt;
    }
  }
  
  /* 
   * Compress 
   */
  
  int idx = 0;
  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_stride_packA[i];
    int end = b_b_stride_packA[i+1];
    if (end > beg) {
      for (int k = 0; k < 3; k++) {      
        b_b_cNewPointsA_true_pack[idx++] = b_b_cNewPointsA_true_pack[3*beg + k];
      }
    }
  }
  
  /* 
   * Define and copy Numabs 
   */
  
  PDM_g_num_t *b_b_gNumVtxFromBForA = malloc (sizeof(PDM_g_num_t) * b_b_stride_packA[n_BForA_gnum]);
  
  PDM_g_num_t beg_n_BForA_gnum;
  
  PDM_MPI_Iscan (&n_active_BForA_gnum, &beg_n_BForA_gnum, 1, 
             PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);
  
  PDM_MPI_Wait (&request1);
  
  beg_n_BForA_gnum += -n_active_BForA_gnum + beg_nPtsFromBForA + 1; 

  idx = 0;
  
  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_stride_packA[i];
    int end = b_b_stride_packA[i+1];
    if (end > beg) {
      b_b_gNumVtxFromBForA[idx++] = beg_n_BForA_gnum++;
    }
  }
  
  beg_n_BForA_gnum += -1;
          
  PDM_MPI_Gather (&beg_n_BForA_gnum, 1, PDM__PDM_MPI_G_NUM,
              &nAbsNewVtxA, 1, PDM__PDM_MPI_G_NUM,
              lastRank,
              _ei->comm);
  
  /* 
   * New points from A for B  :
   *   - Compute absolute number
   *   - Synchronize coordinates  
   */

  PDM_MPI_Wait (&request2);
  
  double *b_cNewPointsB_true_pack = malloc (sizeof(double) * 3 * nPtsFromAForB);
  PDM_edges_intersect_point_t *b_oNewPointsB_true_pack = 
          malloc (sizeof(PDM_edges_intersect_point_t) * nPtsFromAForB);
  PDM_g_num_t *b_lNewPointsB_true_pack = malloc (sizeof(PDM_g_num_t) * nPtsFromAForB);

  nPtsFromAForB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if ((b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) ||
        (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromAForB += 1;
      for (int k = 0; k < 3; k++) {
        b_cNewPointsB_true_pack[3*nPtsFromAForB + k] = b_coordsPointB_true[3*i+k];
      }
      b_oNewPointsB_true_pack[nPtsFromAForB] = b_oNewPointsB_true[i];
      b_lNewPointsB_true_pack[nPtsFromAForB] = b_connectPointB_true[i];
      nPtsFromAForB += 1;
    }
  }

  PDM_part_to_block_t *ptbAForB = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          (PDM_g_num_t **) &b_lNewPointsB_true_pack,  
                                                          &nPtsFromAForB,
                                                          1,
                                                          _ei->comm);
 
  int n_AForB_gnum = PDM_part_to_block_n_elt_block_get (ptbAForB);
  
  int *b_stride_packB = malloc (sizeof(int) * nPtsFromAForB);
  for (int i = 0; i < nPtsFromAForB; i++) {
    b_stride_packB[i] = 3;
  }
  
  int *b_b_stride_packB = NULL;
  double *b_b_cNewPointsB_true_pack = NULL;
  PDM_part_to_block_exch (ptbAForB,
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         1,
                         &b_stride_packB,
                         (void **)&b_cNewPointsB_true_pack,
                         &b_b_stride_packB,
                         (void **)&b_b_cNewPointsB_true_pack);
  
  free (b_b_stride_packB);
  for (int i = 0; i < nPtsFromAForB; i++) {
    b_stride_packB[i] = 1;
  }

  PDM_edges_intersect_point_t *b_b_oNewPointsB_true_pack = NULL;
  PDM_part_to_block_exch (ptbAForB,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR,
                         1,
                         &b_stride_packB,
                         (void **)&b_oNewPointsB_true_pack,
                         &b_b_stride_packB,
                         (void **)&b_b_oNewPointsB_true_pack); 
  free (b_stride_packB);
  
  /* 
   * Synchronize coordinates 
   */
  
  PDM_g_num_t n_active_AForB_gnum = 0;
  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_stride_packB[i];
    int end = b_b_stride_packB[i+1];  
    int cpt = 0;
    int state = 0;
    PDM_edges_intersect_point_t eip;
    for (int j = beg; j < end; j++) {
      n_active_AForB_gnum += 1;
      if (state == 0) {
        state = 1;
        cpt   = 1;
        eip   = b_b_oNewPointsB_true_pack[j];
      }
      else {
        if (eip != b_b_oNewPointsB_true_pack[j]) {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersection : Inconsistencies : A B point is on a A vertex and on a A edge");
          abort();
        }
        cpt += 1;
        
        for (int k = 0; k < 3; k++) {      
          b_b_cNewPointsB_true_pack[3 * beg + k] += b_b_cNewPointsB_true_pack[3 * j + k];
        }
      }
    }  
    for (int k = 0; k < 3; k++) {      
      b_b_cNewPointsB_true_pack[3 * beg + k] = b_b_cNewPointsB_true_pack[3 * beg + k] / cpt;
    }
  }
  
  /* 
   * Compress 
   */
  
  idx = 0;
  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_stride_packB[i];
    int end = b_b_stride_packB[i+1];
    if (end > beg) {
      for (int k = 0; k < 3; k++) {      
        b_b_cNewPointsB_true_pack[idx++] = b_b_cNewPointsB_true_pack[3*beg + k];
      }
    }
  }
  
  /* 
   * Define and copy Numabs 
   */
  
  PDM_g_num_t *b_b_gNumVtxFromAForB = malloc (sizeof(PDM_g_num_t) * b_b_stride_packB[n_AForB_gnum]);
  
  PDM_g_num_t beg_n_AForB_gnum;
  
  PDM_MPI_Iscan (&n_active_AForB_gnum, &beg_n_AForB_gnum, 1, 
             PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);
  
  PDM_MPI_Wait (&request1);
  
  beg_n_AForB_gnum += -n_active_AForB_gnum + beg_nPtsFromAForB + 1; 

  idx = 0;
  
  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_stride_packB[i];
    int end = b_b_stride_packB[i+1];
    if (end > beg) {
      b_b_gNumVtxFromAForB[idx++] = beg_n_AForB_gnum++;
    }
  } 

  beg_n_AForB_gnum += -1;
          
  PDM_MPI_Gather (&beg_n_AForB_gnum, 1, PDM__PDM_MPI_G_NUM,
              &nAbsNewVtxB, 1, PDM__PDM_MPI_G_NUM,
              lastRank,
              _ei->comm);
  
  /*
   * block to part the true intersection and absolute number
   *   - For A
   *   - For B
   */

  PDM_g_num_t    *blockDistribIdxA = PDM_part_to_block_distrib_index_get (ptbBForA);
  
  PDM_block_to_part_t *pbtBForA = PDM_block_to_part_create (blockDistribIdxA,
                                                           &b_lNewPointsA_true_pack,  
                                                           &nPtsFromBForA,
                                                            1,
                                                            _ei->comm);

  int    *part_strideA = NULL;
  double *b_cNewPointsA_true_merge = NULL;
  PDM_g_num_t  *b_cNewPointsA_true_gnum = NULL;

  PDM_block_to_part_exch (pbtBForA,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_b_stride_packA, 
                          (void *) b_b_gNumVtxFromBForA,
                          &part_strideA,
                          (void **) &b_cNewPointsA_true_gnum);

  for (int i = 0; i < (n_BForA_gnum + 1); i++) {
    b_b_stride_packA[i] *= 3;
  }
  
  free (part_strideA);

  PDM_block_to_part_exch (pbtBForA,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_b_stride_packA, 
                          (void *) b_b_gNumVtxFromBForA,
                          &part_strideA,
                          (void **) &b_cNewPointsA_true_merge);

  
  PDM_g_num_t    *blockDistribIdxB = PDM_part_to_block_distrib_index_get (ptbAForB);
  
  PDM_block_to_part_t *pbtAForB = PDM_block_to_part_create (blockDistribIdxB,
                                                           &b_lNewPointsB_true_pack,  
                                                           &nPtsFromAForB,
                                                            1,
                                                            _ei->comm);

  int    *part_strideB = NULL;
  double *b_cNewPointsB_true_merge = NULL;
  PDM_g_num_t  *b_cNewPointsB_true_gnum = NULL;

  PDM_block_to_part_exch (pbtAForB,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_b_stride_packB, 
                          (void *) b_b_gNumVtxFromAForB,
                          &part_strideB,
                          (void **) &b_cNewPointsB_true_gnum);

  for (int i = 0; i < (n_AForB_gnum + 1); i++) {
    b_b_stride_packB[i] *= 3;
  }

  free (part_strideB);
  
  PDM_block_to_part_exch (pbtAForB,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_b_stride_packB, 
                          (void *) b_b_gNumVtxFromAForB,
                          &part_strideB,
                          (void **) &b_cNewPointsB_true_merge);

  PDM_block_to_part_free (pbtBForA);
  PDM_part_to_block_free (ptbBForA);

  PDM_block_to_part_free (pbtAForB);
  PDM_part_to_block_free (ptbAForB);

  /*
   * Update _edges_intersect_structure (pdm_block_to_part)
   *  transfer :
   *
   * b_strideOne_idx_true;
   * b_tIntersects_true;
   * b_gNumEdgeA_true;
   * b_nNewPointsA_true;
   * b_oNewPointsA_true;
   * b_connectPointA_true;
   * b_uPointA_true;
   * b_coordsPointA_true;
   * b_gNumA_true;   
   * b_nNewPointsB_true;
   * b_oNewPointsB_true;
   * b_connectPointB_true; 
   * b_uPointB_true;
   * b_coordsPointB_true;
   * b_gNumB_true;   
   * 
   */

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb); 
  
  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       &keys,
                                                       &nProcData,
                                                       1,
                                                       _ei->comm);
  
  int *b_strideOne_true = malloc (sizeof(int) * n_elt_block);
  
  for (int i = 0; i < n_elt_block; i++) {
    b_strideOne_true [i] = b_strideOne_idx_true[i+1] - b_strideOne_idx_true[i]; 
  } 
  
  int *r_strideOne_true;
  int *r_tIntersects_true;
  
  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          b_strideOne_true, 
                          (void *) b_tIntersects_true,
                          &r_strideOne_true,
                          (void **) &r_tIntersects_true);
  free (r_strideOne_true);
  
  PDM_g_num_t *r_gNumEdgeA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_strideOne_true, 
                          (void *) b_gNumEdgeA_true,
                          &r_strideOne_true,
                          (void **) &r_gNumEdgeA_true);
  free (r_strideOne_true);
  
  int *r_nNewPointsA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          b_strideOne_true, 
                          (void *) b_nNewPointsA_true,
                          &r_strideOne_true,
                          (void **) &r_nNewPointsA_true);
  free (r_strideOne_true);
  
  int *b_stridePtsADep_true = malloc (sizeof(int) * n_elt_block);             
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsADep_true[i] = 0;   
    for (int j = b_strideOne_idx_true[i]; j < b_strideOne_idx_true[i+1]; j++) {
      b_stridePtsADep_true[i] += b_nNewPointsA_true[j];   
    }
  }   
  
  int *r_stridePtsADep_true = NULL;
  PDM_edges_intersect_point_t *r_oNewPointsA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_edges_intersect_point_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsADep_true, 
                          (void *) b_oNewPointsA_true,
                          &r_stridePtsADep_true,
                          (void **) &b_oNewPointsA_true);

  free (r_stridePtsADep_true);
  PDM_g_num_t *r_connectPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsADep_true, 
                          (void *) b_connectPointA_true,
                          &r_stridePtsADep_true,
                          (void **) &r_connectPointA_true);

  free (r_stridePtsADep_true);
  double *r_uPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_stridePtsADep_true, 
                          (void *) b_uPointA_true,
                          &r_stridePtsADep_true,
                          (void **) &r_uPointA_true);
  
  free (r_stridePtsADep_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsADep_true[i] *= 3;
  }
  double *r_coordsPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_stridePtsADep_true, 
                          (void *) b_coordsPointA_true,
                          &r_stridePtsADep_true,
                          (void **) &r_coordsPointA_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsADep_true[i] = b_stridePtsADep_true[i]/3;
  }

  free (r_stridePtsADep_true);
  PDM_g_num_t *r_gNumPointsA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsADep_true, 
                          (void *) b_gNumPointsA_true,
                          &r_stridePtsADep_true,
                          (void **) &r_gNumPointsA_true);
  
  int *r_nNewPointsB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          b_strideOne_true, 
                          (void *) b_nNewPointsB_true,
                          &r_strideOne_true,
                          (void **) &r_nNewPointsB_true);
  
  int *b_stridePtsBDep_true = malloc (sizeof(int) * n_elt_block);             
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsBDep_true[i] = 0;   
    for (int j = b_strideOne_idx_true[i]; j < b_strideOne_idx_true[i+1]; j++) {
      b_stridePtsBDep_true[i] += b_nNewPointsB_true[j];   
    }
  }   
  
  int *r_stridePtsBDep_true = NULL;
  PDM_edges_intersect_point_t *r_oNewPointsB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_edges_intersect_point_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsBDep_true, 
                          (void *) b_oNewPointsB_true,
                          &r_stridePtsBDep_true,
                          (void **) &b_oNewPointsB_true);

  free (r_stridePtsBDep_true);
  PDM_g_num_t *r_connectPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsBDep_true, 
                          (void *) b_connectPointB_true,
                          &r_stridePtsBDep_true,
                          (void **) &r_connectPointB_true);

  free (r_stridePtsBDep_true);
  double *r_uPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_stridePtsBDep_true, 
                          (void *) b_uPointB_true,
                          &r_stridePtsBDep_true,
                          (void **) &r_uPointB_true);
  
  free (r_stridePtsBDep_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsBDep_true[i] *= 3;
  }
  double *r_coordsPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          b_stridePtsBDep_true, 
                          (void *) b_coordsPointB_true,
                          &r_stridePtsBDep_true,
                          (void **) &r_coordsPointB_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsBDep_true[i] = b_stridePtsBDep_true[i]/3;
  }

  free (r_stridePtsBDep_true);
  PDM_g_num_t *r_gNumPointsB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          b_stridePtsBDep_true, 
                          (void *) b_gNumPointsB_true,
                          &r_stridePtsBDep_true,
                          (void **) &r_gNumPointsB_true);

  /*
   * Cleanup 
   */

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);
  
  free (b_cNewPointsB_true_merge);
  free (b_cNewPointsA_true_merge);
  free (part_strideA);
  free (part_strideB);
  
  free (b_connectPointA_true);
  free (b_connectPointB_true);
  free (b_coordsPointA_true);  
  free (b_coordsPointB_true);
  free (b_gNumEdgeA_true);
  free (b_nNewPointsA_true);
  free (b_nNewPointsB_true);
  free (b_oNewPointsA_true);
  free (b_oNewPointsB_true);
  free (b_strideOne_idx_true);
  free (b_tIntersects_true);
  free (b_gNumPointsA_true);
  free (b_gNumPointsB_true);
 
  free (b_b_gNumVtxFromAForB);
  free (b_b_stride_packB);
  free (b_b_cNewPointsB_true_pack);
  
  free (b_b_gNumVtxFromBForA);
  free (b_b_stride_packA);
  free (b_b_cNewPointsA_true_pack);

  free (b_stridePtsADep_true);

  /*
   * Update edges intersection structure structures
   * 
   *  keys : absolu ou local ?
   *  Comment relier au mieux la structure _edges_intersect_res_t* avec
   *  les tableau r_ 
   * 
   */

  int idxData = 0;
  
  int *r_strideOne_idx_true     = malloc (sizeof(int) * (nProcData + 1));
  int *r_stridePtsADep_idx_true = malloc (sizeof(int) * (nProcData + 1));
  int *r_stridePtsBDep_idx_true = malloc (sizeof(int) * (nProcData + 1));
  
  r_strideOne_idx_true[0]     = 0;
  r_stridePtsADep_idx_true[0] = 0;
  r_stridePtsBDep_idx_true[0] = 0;
  
  for (int i = 0; i < nProcData; i++) {
    r_strideOne_idx_true[i+1] = r_strideOne_idx_true[i] + r_strideOne_true[i]; 
    r_stridePtsADep_idx_true[i+1] = r_stridePtsADep_idx_true[i] + r_stridePtsADep_true[i]; 
    r_stridePtsBDep_idx_true[i+1] = r_stridePtsBDep_idx_true[i] + r_stridePtsBDep_true[i]; 
  }
  
  for (int key = 0; key < keyMax; key++) {
  
    _edges_intersect_res_t ** datas = 
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht, 
                                                               (void *) &key);   

    int nData = PDM_hash_tab_n_data_get (ht, &key);
    
    for (int i = 0; i < nData; i++) {
      
      _edges_intersect_res_t *_data = datas[i];

      int isFound = 0;
      
      int idxA = r_stridePtsADep_idx_true[idxData];
      int idxB = r_stridePtsBDep_idx_true[idxData];
      
      for (int k = r_strideOne_idx_true[idxData]; k < r_strideOne_idx_true[idxData+1]; k++) {
        
        if (_data->nGEdgeA != r_gNumEdgeA_true[k]) {
          idxA += r_nNewPointsA_true[k];
          idxB += r_nNewPointsB_true[k];
          continue;
        }
        
        isFound = 1;
        _data->tIntersect = (PDM_line_intersect_t) r_tIntersects_true[k];       
        const int _nNewPointsA = r_nNewPointsA_true[k];

        if (_data->nNewPointsA != _nNewPointsA) {
          _data->nNewPointsA = _nNewPointsA;        
          _data->uA          = realloc (_data->uA, sizeof(double) * _nNewPointsA);
          _data->coordsA     = realloc (_data->coordsA, sizeof(double) * 3 * _nNewPointsA);
          _data->linkA       = realloc (_data->linkA, sizeof(PDM_g_num_t) * _nNewPointsA);
          _data->gNumA       = realloc (_data->gNumA, sizeof(PDM_g_num_t) * _nNewPointsA);
          _data->oNewPointsA = realloc (_data->oNewPointsA, sizeof(PDM_edges_intersect_point_t) * _nNewPointsA);
        }
        
        _data->nNewPointsA = _nNewPointsA;
        
        for (int k1 = 0; k1 < _nNewPointsA; k1++) {
          _data->oNewPointsA[k1] = r_oNewPointsA_true[idxA+k1];
          _data->gNumA[k1]       = r_gNumPointsA_true[idxA+k1];
          _data->linkA[k1]       = r_connectPointA_true[idxA+k1];
          _data->uA[k1]          = r_uPointA_true[idxA+k1];
          for (int k2 = 0; k2 < 3; k2++) {
            _data->coordsA[3*k1+k2] = r_coordsPointA_true[3*(idxA+k1)+k2];
          }
        }
        
        const int _nNewPointsB = r_nNewPointsB_true[k];

        if (_data->nNewPointsB != _nNewPointsB) {
          _data->nNewPointsB = _nNewPointsB;        
          _data->uB          = realloc (_data->uB, sizeof(double) * _nNewPointsB);
          _data->coordsB     = realloc (_data->coordsB, sizeof(double) * 3 * _nNewPointsB);
          _data->linkB       = realloc (_data->linkB, sizeof(PDM_g_num_t) * _nNewPointsB);
          _data->gNumB       = realloc (_data->gNumB, sizeof(PDM_g_num_t) * _nNewPointsB);
          _data->oNewPointsB = realloc (_data->oNewPointsB, sizeof(PDM_edges_intersect_point_t) * _nNewPointsB);
        }
        
        _data->nNewPointsB = _nNewPointsB;
        
        for (int k1 = 0; k1 < _nNewPointsB; k1++) {
          _data->oNewPointsB[k1] = r_oNewPointsB_true[idxB+k1];
          _data->gNumB[k1]       = r_gNumPointsB_true[idxB+k1];
          _data->linkB[k1]       = r_connectPointB_true[idxB+k1];
          _data->uB[k1]          = r_uPointB_true[idxB+k1];
          for (int k2 = 0; k2 < 3; k2++) {
            _data->coordsA[3*k1+k2] = r_coordsPointA_true[3*(idxB+k1)+k2];
          }
        }

        break;
      }
      
      if (!isFound) {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersections :"
                         "No Data found to update current intersection\n");
        abort();
      }
      idxData += 1;
    }
  }
  
  /*
   * Cleanup  
   */

  free (r_strideOne_true);
  free (r_tIntersects_true);
  free (r_gNumEdgeA_true);

  free (r_stridePtsADep_true);
  free (r_nNewPointsA_true);
  free (r_oNewPointsA_true);
  free (r_connectPointA_true);
  free (r_uPointA_true);
  free (r_coordsPointA_true);
  free (r_gNumPointsA_true);

  free (r_stridePtsBDep_true);
  free (r_nNewPointsB_true);
  free (r_oNewPointsB_true);
  free (r_connectPointB_true); 
  free (r_uPointB_true);
  free (r_coordsPointB_true);
  free (r_gNumPointsB_true);  
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

