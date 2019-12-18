/*============================================================================
 * Encapsulation de MPI
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/


/*============================================================================
 * Definition des variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Indirection sur le code d'erreur
 *----------------------------------------------------------------------------*/

static const int mpi_err[] = {

  MPI_SUCCESS,
  MPI_ERR_BUFFER,
  MPI_ERR_COUNT,
  MPI_ERR_TYPE,
  MPI_ERR_TAG,
  MPI_ERR_COMM,
  MPI_ERR_RANK,
  MPI_ERR_ROOT,
  MPI_ERR_TRUNCATE,
  MPI_ERR_GROUP,
  MPI_ERR_OP,
  MPI_ERR_REQUEST,
  MPI_ERR_TOPOLOGY,
  MPI_ERR_DIMS,
  MPI_ERR_ARG,
  MPI_ERR_UNKNOWN,
  MPI_ERR_OTHER,
  MPI_ERR_INTERN,
  MPI_ERR_IN_STATUS,
  MPI_ERR_PENDING,
  MPI_MAX_ERROR_STRING,



  MPI_ERR_ACCESS,
  MPI_ERR_AMODE,
  MPI_ERR_BAD_FILE,
  MPI_ERR_CONVERSION,
  MPI_ERR_DUP_DATAREP,
  MPI_ERR_FILE_EXISTS,
  MPI_ERR_FILE_IN_USE,
  MPI_ERR_FILE,
  MPI_ERR_INFO_KEY,
  MPI_ERR_INFO_NOKEY,
  MPI_ERR_INFO_VALUE,
  MPI_ERR_IO,
  MPI_ERR_NO_MEM,
  MPI_ERR_NOT_SAME,
  MPI_ERR_NO_SPACE,
  MPI_ERR_NO_SUCH_FILE,
  MPI_ERR_QUOTA,
  MPI_ERR_READ_ONLY,
  MPI_ERR_UNSUPPORTED_DATAREP,
  MPI_ERR_UNSUPPORTED_OPERATION,
  MPI_ERR_WIN,
  MPI_ERR_LASTCODE,
  MPI_ERR_ASSERT,
  MPI_ERR_BASE,
  MPI_ERR_DISP,
  MPI_ERR_KEYVAL,
  MPI_ERR_LOCKTYPE,
  MPI_ERR_RMA_CONFLICT,
  MPI_ERR_RMA_SYNC,
  MPI_ERR_SIZE

};

/*----------------------------------------------------------------------------
 * Indirection sur le mode du fichier
 *----------------------------------------------------------------------------*/

static const int mpi_file_mode[] = {

MPI_MODE_CREATE,
MPI_MODE_RDONLY,
MPI_MODE_WRONLY,
MPI_MODE_RDWR,
MPI_MODE_DELETE_ON_CLOSE,
MPI_MODE_UNIQUE_OPEN,
MPI_MODE_EXCL,
MPI_MODE_APPEND,
MPI_MODE_SEQUENTIAL,
MPI_DISPLACEMENT_CURRENT,
MPI_SEEK_SET,
MPI_SEEK_CUR,
MPI_SEEK_END,
MPI_MODE_WRONLY | MPI_MODE_APPEND,
MPI_MODE_WRONLY | MPI_MODE_CREATE

};

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Datatype -> MPI_Datatype
 *----------------------------------------------------------------------------*/

static const MPI_Datatype mpi_datatype_cste[] = {

  MPI_BYTE,
  MPI_PACKED,
  MPI_CHAR,
  MPI_SHORT,
  MPI_INT,
  MPI_LONG,
  MPI_FLOAT,
  MPI_DOUBLE,
  MPI_LONG_DOUBLE,
  MPI_UNSIGNED_CHAR,
  MPI_UNSIGNED_SHORT,
  MPI_UNSIGNED_LONG,
  MPI_UNSIGNED,
  MPI_FLOAT_INT,
  MPI_DOUBLE_INT,
  MPI_LONG_DOUBLE_INT,
  MPI_LONG_INT,
  MPI_SHORT_INT,
  MPI_2INT,
  MPI_CHARACTER,
  MPI_INTEGER,
  MPI_REAL,
  MPI_DOUBLE_PRECISION,
  MPI_DATATYPE_NULL,
  MPI_INT8_T,
  MPI_INT16_T,
  MPI_INT32_T,
  MPI_INT64_T,
  MPI_UINT8_T,
  MPI_UINT16_T,
  MPI_UINT32_T,
  MPI_UINT64_T
};

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Op -> MPI_Op
 *----------------------------------------------------------------------------*/

static const MPI_Op mpi_op[] = {

  MPI_MAX,
  MPI_MIN,
  MPI_SUM,
  MPI_OP_NULL

};



/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_File ->constantes MPI_File
 *----------------------------------------------------------------------------*/

static const MPI_File mpi_file_cste[] = {

  MPI_FILE_NULL

};

/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_Comm ->constantes MPI_Comm
 *----------------------------------------------------------------------------*/

static const MPI_Comm mpi_comm_cste[] = {

  MPI_COMM_NULL,
  MPI_COMM_WORLD

};

/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_Request ->constantes MPI_Request
 *----------------------------------------------------------------------------*/

static const MPI_Request mpi_request_cste[] = {

  MPI_REQUEST_NULL,

};

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_File -> MPI_File
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_File **mpi_file   = NULL; /* Tableau de stockage */
static int       l_mpi_file = 0;     /* Taille du tableau */
static int       n_mpi_file = 0;     /* Nombre de fichiers stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Comm -> MPI_Comm
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_Comm **mpi_comm   = NULL; /* Tableau de stockage */
static int       l_mpi_comm = 0;     /* Taille du tableau */
static int       n_mpi_comm = 0;     /* Nombre de communicateurs stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Request -> MPI_Request
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_Request **mpi_request = NULL; /* Tableau de stockage */
static int       l_mpi_request = 0;   /* Taille du tableau */
static int       n_mpi_request = 0;   /* Nombre de communicateurs stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Datatype -> MPI_Datatype
 * stockage dans un tableau des types utilisateurs
 *----------------------------------------------------------------------------*/

static MPI_Datatype **mpi_datatype   = NULL; /* Tableau de stockage */
static int           l_mpi_datatype = 0;     /* Taille du tableau */
static int           n_mpi_datatype = 0;     /* Nombre de communicateurs stockes */

/*============================================================================
 * Defintion des fonctions pprivees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * mpi_err -> pdm_mpi_err
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

static int _mpi_2_pdm_mpi_err(int code_mpi)
{
  int code = PDM_MPI_ERR_OTHER;
  switch(code_mpi) {
  case MPI_SUCCESS:
    code = PDM_MPI_SUCCESS;
  break;
  case MPI_ERR_BUFFER:
    code = PDM_MPI_ERR_BUFFER;
      break;
  case MPI_ERR_COUNT:
    code = PDM_MPI_ERR_COUNT;
      break;
  case MPI_ERR_TYPE:
    code = PDM_MPI_ERR_TYPE;
    break;
  case MPI_ERR_TAG:
    code = PDM_MPI_ERR_TAG;
    break;
  case MPI_ERR_COMM:
    code = PDM_MPI_ERR_COMM;
    break;
  case MPI_ERR_RANK:
    code = PDM_MPI_ERR_RANK;
    break;
  case MPI_ERR_ROOT:
    code = PDM_MPI_ERR_ROOT;
    break;
  case MPI_ERR_TRUNCATE:
    code = PDM_MPI_ERR_TRUNCATE;
    break;
  case MPI_ERR_GROUP:
    code = PDM_MPI_ERR_GROUP;
    break;
  case MPI_ERR_OP:
    code = PDM_MPI_ERR_OP;
    break;
  case MPI_ERR_REQUEST:
    code = PDM_MPI_ERR_REQUEST;
    break;
  case MPI_ERR_TOPOLOGY:
    code = PDM_MPI_ERR_TOPOLOGY;
    break;
  case MPI_ERR_DIMS:
    code = PDM_MPI_ERR_DIMS;
    break;
  case MPI_ERR_ARG:
    code = PDM_MPI_ERR_ARG;
    break;
  case MPI_ERR_UNKNOWN:
    code = PDM_MPI_ERR_UNKNOWN;
    break;
  case MPI_ERR_OTHER:
    code = PDM_MPI_ERR_OTHER;
    break;
  case MPI_ERR_INTERN:
    code = PDM_MPI_ERR_INTERN;
    break;
  case MPI_ERR_IN_STATUS:
    code = PDM_MPI_ERR_IN_STATUS;
    break;
  case MPI_ERR_PENDING:
    code = PDM_MPI_ERR_PENDING;
    break;



  case MPI_ERR_ACCESS:
    code = PDM_MPI_ERR_ACCESS;
    break;
  case MPI_ERR_AMODE:
    code = PDM_MPI_ERR_AMODE;
    break;
  case MPI_ERR_BAD_FILE:
    code = PDM_MPI_ERR_BAD_FILE;
    break;
  case MPI_ERR_CONVERSION:
    code = PDM_MPI_ERR_CONVERSION;
    break;
  case MPI_ERR_DUP_DATAREP:
    code = PDM_MPI_ERR_DUP_DATAREP;
    break;
  case MPI_ERR_FILE_EXISTS:
    code = PDM_MPI_ERR_FILE_EXISTS;
    break;
  case MPI_ERR_FILE_IN_USE:
    code = PDM_MPI_ERR_FILE_IN_USE;
    break;
  case MPI_ERR_FILE:
    code = PDM_MPI_ERR_FILE;
    break;
  case MPI_ERR_INFO_KEY:
    code = PDM_MPI_ERR_INFO_KEY;
    break;
  case MPI_ERR_INFO_NOKEY:
    code = PDM_MPI_ERR_INFO_NOKEY;
    break;
  case MPI_ERR_INFO_VALUE:
    code = PDM_MPI_ERR_INFO_VALUE;
    break;
  case MPI_ERR_IO:
    code = PDM_MPI_ERR_IO;
    break;
  case MPI_ERR_NO_MEM:
    code = PDM_MPI_ERR_NO_MEM;
    break;
  case MPI_ERR_NOT_SAME:
    code = PDM_MPI_ERR_NOT_SAME;
    break;
  case MPI_ERR_NO_SPACE:
    code = PDM_MPI_ERR_NO_SPACE;
    break;
  case MPI_ERR_NO_SUCH_FILE:
    code = PDM_MPI_ERR_NO_SUCH_FILE;
    break;
  case MPI_ERR_QUOTA:
    code = PDM_MPI_ERR_QUOTA;
    break;
  case MPI_ERR_READ_ONLY:
    code = PDM_MPI_ERR_READ_ONLY;
    break;
  case MPI_ERR_UNSUPPORTED_DATAREP:
    code = PDM_MPI_ERR_UNSUPPORTED_DATAREP;
    break;
  case MPI_ERR_UNSUPPORTED_OPERATION:
    code = PDM_MPI_ERR_UNSUPPORTED_OPERATION;
    break;
  case MPI_ERR_WIN:
    code = PDM_MPI_ERR_WIN;
    break;
  case MPI_ERR_LASTCODE:
    code = PDM_MPI_ERR_LASTCODE;
    break;



  case MPI_ERR_ASSERT:
    code = PDM_MPI_ERR_ASSERT;
    break;
  case MPI_ERR_BASE:
    code = PDM_MPI_ERR_BASE;
    break;
  case MPI_ERR_DISP:
    code = PDM_MPI_ERR_DISP;
    break;
  case MPI_ERR_KEYVAL:
    code = PDM_MPI_ERR_KEYVAL;
    break;
  case MPI_ERR_LOCKTYPE:
    code = PDM_MPI_ERR_LOCKTYPE;
    break;
  case MPI_ERR_RMA_CONFLICT:
    code = PDM_MPI_ERR_RMA_CONFLICT;
    break;
  case MPI_ERR_RMA_SYNC:
    code = PDM_MPI_ERR_RMA_SYNC;
    break;
  case MPI_ERR_SIZE:
    code = PDM_MPI_ERR_SIZE;

  }
  return code;
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

static MPI_Comm _pdm_mpi_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_comm < 0)
    return mpi_comm_cste[-pdm_mpi_comm - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_comm < l_mpi_comm)
      return *(mpi_comm[pdm_mpi_comm]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_comm :"
            " pdm_mpi_comm '%d' non valide\n", pdm_mpi_comm);
      abort();
    }
  }
}

/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_comm
 *
 * MPI_Comm -> PDM_MPI_Comm
 *----------------------------------------------------------------------------*/

static PDM_MPI_Comm _mpi_2_pdm_mpi_comm(MPI_Comm comm)
{

  /* Traitement des communicateurs predefinis */

  if (comm == MPI_COMM_NULL)
    return PDM_MPI_COMM_NULL;

  else if (comm == MPI_COMM_WORLD)
    return PDM_MPI_COMM_WORLD;

  /* Traitement des communicateurs utilisateurs */

  else {

    /* Recherche du communicateur MSG correspondant au communicateur MPI */


    if (mpi_comm != NULL) {
      for (int i = 0; i < l_mpi_comm; i++)
        if (mpi_comm[i] != NULL)
          if (*(mpi_comm[i]) == comm)
            return (PDM_MPI_Comm) i;
    }

    /* Si non trouve cree un nouveau communicateur MSG */

    if (mpi_comm == NULL) {
      l_mpi_comm = 4;
      mpi_comm = (MPI_Comm **) malloc(sizeof(MPI_Comm *) * l_mpi_comm);
      for (int i = 0; i < l_mpi_comm; i++)
        mpi_comm[i] = NULL;
    }

    if (l_mpi_comm <= n_mpi_comm) {
      int  p_l_mpi_comm = l_mpi_comm;
      l_mpi_comm = 2 * l_mpi_comm;
      mpi_comm = (MPI_Comm **) realloc((void*) mpi_comm,
                                       l_mpi_comm *
                                       sizeof(MPI_Comm *));
      for (int i = p_l_mpi_comm; i < l_mpi_comm; i++)
        mpi_comm[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_comm[i] != NULL)
      i++;

    mpi_comm[i] = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    *(mpi_comm[i]) = comm;
    n_mpi_comm += 1;

    return (PDM_MPI_Comm) i;
  }
}


/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_request
 *
 * PDM_MPI_Request -> MPI_Request
 *----------------------------------------------------------------------------*/

static MPI_Request _pdm_mpi_2_mpi_request(PDM_MPI_Request pdm_mpi_request)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_request < 0)
    return mpi_request_cste[-pdm_mpi_request - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_request < l_mpi_request)
      return *(mpi_request[pdm_mpi_request]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_request :"
            " pdm_mpi_request '%d' non valide\n", pdm_mpi_request);
      abort();
    }
  }
}


/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_request
 *
 * MPI_Request -> PDM_MPI_Request
 *----------------------------------------------------------------------------*/

static PDM_MPI_Request _mpi_2_pdm_mpi_request(MPI_Request request)
{

  /* Traitement des communicateurs predefinis */

  if (request == MPI_REQUEST_NULL)
    return PDM_MPI_REQUEST_NULL;

  /* Traitement des communicateurs utilisateurs */

  else {

    /* Recherche du communicateur MSG correspondant au communicateur MPI */

    if (mpi_request != NULL) {
      for (int i = 0; i < l_mpi_request; i++)
        if (mpi_request[i] != NULL)
          if (*(mpi_request[i]) == request)
            return (PDM_MPI_Request) i;
    }

    /* Si non trouve cree un nouveau communicateur MSG */

    if (mpi_request == NULL) {
      l_mpi_request = 4;
      mpi_request = (MPI_Request **) malloc(sizeof(MPI_Request *) * l_mpi_request);
      for (int i = 0; i < l_mpi_request; i++)
        mpi_request[i] = NULL;
    }

    if (l_mpi_request <= n_mpi_request) {
      int  p_l_mpi_request = l_mpi_request;
      l_mpi_request = 2 * l_mpi_request;
      mpi_request = (MPI_Request **) realloc((void*) mpi_request,
                                             l_mpi_request *
                                             sizeof(MPI_Request *));
      for (int i = p_l_mpi_request; i < l_mpi_request; i++)
        mpi_request[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_request[i] != NULL)
      i++;

    mpi_request[i] = (MPI_Request *) malloc(sizeof(MPI_Request));
    *(mpi_request[i]) = request;
    n_mpi_request += 1;

    return (PDM_MPI_Request) i;
  }
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_datatype
 *
 * PDM_MPI_Datatype -> MPI_Datatype
 *----------------------------------------------------------------------------*/

static MPI_Datatype _pdm_mpi_2_mpi_datatype(PDM_MPI_Datatype pdm_mpi_datatype)
{

  /* Traitement des MPI_Datatype connus  */

  if (pdm_mpi_datatype < 0)
    return mpi_datatype_cste[-pdm_mpi_datatype - 1];

  /* Traitement des MPI_Datatype utilisateurs  */

  else {
    if (pdm_mpi_datatype < l_mpi_datatype)
      return *(mpi_datatype[pdm_mpi_datatype]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_datatype :"
            " pdm_mpi_datatype '%d' non valide\n", pdm_mpi_datatype);
      abort();
    }
  }
}


/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_data
 *
 * MPI_Datatype -> PDM_MPI_Datatype
 *----------------------------------------------------------------------------*/

static PDM_MPI_Datatype _mpi_2_pdm_mpi_datatype(MPI_Datatype datatype)
{

  /* Traitement des communicateurs predefinis */

  if (datatype == MPI_BYTE)
    return PDM_MPI_BYTE;
  else if (datatype == MPI_PACKED)
    return  PDM_MPI_PACKED;
  else if (datatype == MPI_CHAR)
    return  PDM_MPI_CHAR;
  else if (datatype == MPI_SHORT)
    return  PDM_MPI_SHORT;
  else if (datatype == MPI_INT)
    return  PDM_MPI_INT;
  else if (datatype == MPI_LONG)
    return  PDM_MPI_LONG;
  else if (datatype == MPI_FLOAT)
    return  PDM_MPI_FLOAT;
  else if (datatype == MPI_DOUBLE)
    return  PDM_MPI_DOUBLE;
  else if (datatype == MPI_LONG_DOUBLE)
    return  PDM_MPI_LONG_DOUBLE;
  else if (datatype == MPI_UNSIGNED_CHAR)
    return  PDM_MPI_UNSIGNED_CHAR;
  else if (datatype == MPI_UNSIGNED_SHORT)
    return  PDM_MPI_UNSIGNED_SHORT;
  else if (datatype == MPI_UNSIGNED_LONG)
    return  PDM_MPI_UNSIGNED_LONG;
  else if (datatype == MPI_UNSIGNED)
    return PDM_MPI_UNSIGNED;
  else if (datatype == MPI_FLOAT_INT)
    return  PDM_MPI_FLOAT_INT;
  else if (datatype == MPI_DOUBLE_INT)
    return  PDM_MPI_DOUBLE_INT;
  else if (datatype == MPI_LONG_DOUBLE_INT)
    return PDM_MPI_LONG_DOUBLE_INT;
  else if (datatype == MPI_LONG_INT)
    return PDM_MPI_LONG_INT;
  else if (datatype == MPI_SHORT_INT)
    return PDM_MPI_SHORT_INT;
  else if (datatype == MPI_2INT)
    return PDM_MPI_2INT;
  else if (datatype == MPI_CHARACTER)
    return PDM_MPI_CHARACTER;
  else if (datatype == MPI_INTEGER)
    return PDM_MPI_INTEGER;
  else if (datatype == MPI_REAL)
    return PDM_MPI_REAL;
  else if (datatype == MPI_DOUBLE_PRECISION)
    return PDM_MPI_DOUBLE_PRECISION;
  else if (datatype == MPI_DATATYPE_NULL)
    return PDM_MPI_DATATYPE_NULL;
  else if (datatype == MPI_INT8_T)
    return PDM_MPI_INT8_T;
  else if (datatype == MPI_INT16_T)
    return PDM_MPI_INT16_T;
  else if (datatype == MPI_INT32_T)
    return PDM_MPI_INT32_T;
  else if (datatype == MPI_INT64_T)
    return PDM_MPI_INT64_T;
  else if (datatype == MPI_UINT8_T)
    return PDM_MPI_UINT8_T;
  else if (datatype == MPI_UINT16_T)
    return PDM_MPI_UINT16_T;
  else if (datatype == MPI_UINT32_T)
    return PDM_MPI_UINT32_T;
  else if (datatype == MPI_UINT64_T)
    return PDM_MPI_UINT64_T;

  /* Traitement des communicateurs utilisateurs */

  else {

    /* Recherche du datatype MSG correspondant au datatype MPI */

    if (mpi_datatype != NULL) {
      for (int i = 0; i < l_mpi_datatype; i++)
        if (mpi_datatype[i] != NULL)
          if (*(mpi_datatype[i]) == datatype)
            return (PDM_MPI_Datatype) i;
    }

    /* Si non trouve cree un nouveau datatype MSG */

    if (mpi_datatype == NULL) {
      l_mpi_datatype = 4;
      mpi_datatype = (MPI_Datatype **)
        malloc(sizeof(MPI_Datatype *) * l_mpi_datatype);
      for (int i = 0; i < l_mpi_datatype; i++)
        mpi_datatype[i] = NULL;
    }

    if (l_mpi_datatype <= n_mpi_datatype) {
      int  p_l_mpi_datatype = l_mpi_datatype;
      l_mpi_datatype = 2 * l_mpi_datatype;
      mpi_datatype = (MPI_Datatype **) realloc((void*) mpi_datatype,
                                       l_mpi_datatype *
                                       sizeof(MPI_Datatype *));
      for (int i = p_l_mpi_datatype; i < l_mpi_datatype; i++)
        mpi_datatype[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_datatype[i] != NULL)
      i++;

    mpi_datatype[i] = (MPI_Datatype *) malloc(sizeof(MPI_Datatype));
    *(mpi_datatype[i]) = datatype;
    n_mpi_datatype += 1;

    return (PDM_MPI_Datatype) i;
  }
}



/*----------------------------------------------------------------------------
 * PDM_MPI_File_Create
 *
 * PDM_MPI_File -> MPI_File
 *----------------------------------------------------------------------------*/

static PDM_MPI_File _pdm_mpi_file_create()
{

  /* Si non trouve, on cree un nouveau fichier MSG */

  if (mpi_file == NULL) {
    l_mpi_file = 4;
      mpi_file = (MPI_File **) malloc(sizeof(MPI_File *) * l_mpi_file);
      for (int i = 0; i < l_mpi_file; i++)
        mpi_file[i] = NULL;
  }

  if (l_mpi_file <= n_mpi_file) {
    int  p_l_mpi_file = l_mpi_file;
    l_mpi_file = 2 * l_mpi_file;
    mpi_file = (MPI_File **) realloc((void*) mpi_file,
                                     l_mpi_file *
                                     sizeof(MPI_File *));
    for (int i = p_l_mpi_file; i < l_mpi_file; i++)
      mpi_file[i] = NULL;
  }

  /* Recherche de la premiere place libre pour stocker le fichier */

  int i = 0;
  while (mpi_file[i] != NULL)
    i++;

  mpi_file[i] = (MPI_File *) malloc(sizeof(MPI_File));
  n_mpi_file += 1;
  return i;
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_file
 *
 * PDM_MPI_File -> MPI_File
 *----------------------------------------------------------------------------*/

static MPI_File _pdm_mpi_2_mpi_file(PDM_MPI_File pdm_mpi_file)
{

  /* Traitement des MPI_File connus  */

  if (pdm_mpi_file < 0)
    return mpi_file_cste[-pdm_mpi_file - 1];

  /* Traitement des MPI_File utilisateurs  */

  else {
    if (pdm_mpi_file < l_mpi_file)
      return *(mpi_file[pdm_mpi_file]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_file :"
              " pdm_mpi_file '%d' non valide\n", pdm_mpi_file);
      abort();
    }
  }
}

/*============================================================================
 * Defintion des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Init(int *argc, char ***argv)
{
  return MPI_Init(argc, argv);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Finalize (void)
{
  if (mpi_file != NULL) {
    for (int i = 0; i < l_mpi_file; i++) {
      if (mpi_file[i] != NULL) {
        MPI_File_close(mpi_file[i]);
        mpi_file[i] = NULL;
      }
    }
    free(mpi_file);
    l_mpi_file = 0;
    n_mpi_file = 0;
  }
  if (mpi_comm != NULL) {
    for (int i = 0; i < l_mpi_comm; i++) {
      if (mpi_comm[i] != NULL) {
        MPI_Comm_free(mpi_comm[i]);
        free (mpi_comm[i]);
        mpi_comm[i] = NULL;
      }
    }
    free(mpi_comm);
    l_mpi_comm = 0;
    n_mpi_comm = 0;
  }

  if (mpi_request != NULL) {
    for (int i = 0; i < l_mpi_request; i++) {
      if (mpi_request[i] != NULL) {
        MPI_Request_free(mpi_request[i]);
        mpi_request[i] = NULL;
      }
    }
    free(mpi_request);
    l_mpi_request = 0;
    n_mpi_request = 0;
  }

  if (mpi_datatype != NULL) {
    for (int i = 0; i < l_mpi_datatype; i++) {
      if (mpi_datatype[i] != NULL) {
        MPI_Type_free(mpi_datatype[i]);
        mpi_datatype[i] = NULL;
      }
    }
    free(mpi_datatype);
    l_mpi_datatype = 0;
    n_mpi_datatype = 0;
  }
  return MPI_Finalize();
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_comm < 0)
    return (void *) &mpi_comm_cste[-pdm_mpi_comm - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_comm < l_mpi_comm)
      return (void *) mpi_comm[pdm_mpi_comm];
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_comm :"
            " pdm_mpi_comm '%d' non valide\n", pdm_mpi_comm);
      abort();
    }
  }
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_free_mpi_comm(void *pt_mpi_comm)
{

  MPI_Comm *comm = (MPI_Comm *) pt_mpi_comm;
  MPI_Comm_free (comm);
  return NULL;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_mpi_2_pdm_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *pt_mpi_comm)
{

  MPI_Comm _mpi_comm = *((MPI_Comm *) pt_mpi_comm);
  return _mpi_2_pdm_mpi_comm(_mpi_comm);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_open(PDM_MPI_Comm comm, char *filename, int amode, PDM_MPI_File *fh)
{

  *fh = _pdm_mpi_file_create();

  char *hints = getenv("PDM_IO_HINTS");

  MPI_Info hints_mpi = MPI_INFO_NULL;

  if (hints != NULL) {

    MPI_Info_create (&hints_mpi);

    char *cp_hints = malloc (sizeof(char *) * (strlen(hints) + 1));
    char *name = malloc (sizeof(char *) * (strlen(hints) + 1));
    char *value = malloc (sizeof(char *) * (strlen(hints) + 1));
    strcpy (cp_hints, hints);

    char *pch;
    char *str2 = cp_hints;

    do {
      pch = strtok (str2,"=");
      str2 = NULL;
      if (pch != NULL) {
        strcpy(name, pch);
        pch = strtok (str2, ":");
        if (pch == NULL) {
          PDM_printf ("Error PDM_MPI_File_open : No value for hint \"%s\"."
                  " Check \"PDM_IO_HINTS\" environment variable\n", name);
          exit(1);
        }
        else {
          strcpy(value, pch);
          MPI_Info_set (hints_mpi, name, value);
          PDM_printf ("MPI/IO hint \"%s\" = \"%s\"\n", name, value);
        }
      }
    } while (pch != NULL);

    free (cp_hints);
    free (name);
    free (value);

  }

  int code = MPI_File_open(_pdm_mpi_2_mpi_comm(comm),
                           filename,
                           mpi_file_mode[amode],
                           hints_mpi,
                           mpi_file[*fh]);

  if (hints != NULL) {

    MPI_Info_free(&hints_mpi);

  }

  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_close(PDM_MPI_File *fh)
{
  int code =  MPI_File_close(mpi_file[*fh]);

  free(mpi_file[*fh]);

  mpi_file[*fh] = NULL;
  n_mpi_file -= 1;

  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_seek(PDM_MPI_File fh, PDM_MPI_Offset offset, int whence)
{
  int code = MPI_File_seek(_pdm_mpi_2_mpi_file(fh),
                           (MPI_Offset) offset,
                           whence);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_size(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_size(_pdm_mpi_2_mpi_file(fh),
                           (MPI_Offset*) &_tmp_offset);
  *offset = _tmp_offset;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_position(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_position(_pdm_mpi_2_mpi_file(fh),
                                    &_tmp_offset);
  *offset = _tmp_offset;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_set_view(PDM_MPI_File fh, PDM_MPI_Offset disp, PDM_MPI_Datatype etype,
	              PDM_MPI_Datatype filetype, char *datarep)
{
  int code = MPI_File_set_view(_pdm_mpi_2_mpi_file(fh),
                               (MPI_Offset) disp,
                               _pdm_mpi_2_mpi_datatype(etype),
                               _pdm_mpi_2_mpi_datatype(filetype),
                               datarep,
                               MPI_INFO_NULL);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_view(PDM_MPI_File fh, PDM_MPI_Offset *disp,
                      PDM_MPI_Datatype *etype, PDM_MPI_Datatype *filetype, char *datarep)
{

  MPI_Datatype mpi_etype;
  MPI_Datatype mpi_filetype;
  MPI_Offset _disp = (MPI_Offset) *disp;

  int code = MPI_File_get_view(_pdm_mpi_2_mpi_file(fh),
                               &_disp,
                               &mpi_etype,
                               &mpi_filetype,
                               datarep);

  *etype    = _mpi_2_pdm_mpi_datatype(mpi_etype);
  *filetype = _mpi_2_pdm_mpi_datatype(mpi_filetype);
  *disp     = (PDM_MPI_Offset) _disp;

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                     int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at(_pdm_mpi_2_mpi_file(fh),
                              (MPI_Offset) offset,
                              buf,
                              count,
                              _pdm_mpi_2_mpi_datatype(datatype),
                              &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at_all(_pdm_mpi_2_mpi_file(fh),
                                  (MPI_Offset) offset,
                                  buf,
                                  count,
                                  _pdm_mpi_2_mpi_datatype(datatype),
                                  &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                      int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at(_pdm_mpi_2_mpi_file(fh),
                               _offset,
                               buf,
                               count,
                               _pdm_mpi_2_mpi_datatype(datatype),
                               &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at_all(_pdm_mpi_2_mpi_file(fh),
                                   (MPI_Offset) _offset,
                                   buf,
                                   count,
                                   _pdm_mpi_2_mpi_datatype(datatype),
                                   &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read(PDM_MPI_File fh, void *buf, int count,
                  PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code =  MPI_File_read(_pdm_mpi_2_mpi_file(fh),
                            buf,
                            count,
                            _pdm_mpi_2_mpi_datatype(datatype),
                            &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_all(PDM_MPI_File fh, void *buf, int count,
                      PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code = MPI_File_read_all(_pdm_mpi_2_mpi_file(fh),
                                buf,
                                count,
                                _pdm_mpi_2_mpi_datatype(datatype),
                                &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write(PDM_MPI_File fh, void *buf, int count,
                   PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code =  MPI_File_write(_pdm_mpi_2_mpi_file(fh),
                             buf,
                             count,
                             _pdm_mpi_2_mpi_datatype(datatype),
                             &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_all(PDM_MPI_File fh, void *buf, int count,
                       PDM_MPI_Datatype datatype, int *n_octet_lus)

{

  MPI_Status status;

  int code =  MPI_File_write_all(_pdm_mpi_2_mpi_file(fh),
                                 buf,
                                 count,
                                 _pdm_mpi_2_mpi_datatype(datatype),
                                 &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                        recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                        root, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Igather (wrapping de la fonction MPI_Igather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Igather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                        recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                        root, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gatherv(sendbuf,
                         sendcount,
                         _pdm_mpi_2_mpi_datatype(sendtype),
                         recvbuf,
                         recvcounts,
                         displs,
                         _pdm_mpi_2_mpi_datatype(recvtype),
                         root,
                         _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
             int tag, PDM_MPI_Comm comm)
{
  int code =  MPI_Recv(buf, count, _pdm_mpi_2_mpi_datatype(datatype), source,
                       tag, _pdm_mpi_2_mpi_comm(comm), MPI_STATUS_IGNORE);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
              int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code =  MPI_Irecv(buf, count, _pdm_mpi_2_mpi_datatype(datatype), source,
                       tag, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
             int tag, PDM_MPI_Comm comm)
{
  int code = MPI_Send(buf, count, _pdm_mpi_2_mpi_datatype(datatype), dest,
                      tag, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Issend (wrapping de la fonction MPI_Issend)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Issend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
               PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Issend(buf, count, _pdm_mpi_2_mpi_datatype(datatype), dest,
                        tag, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Wait(PDM_MPI_Request *request)

{
  MPI_Request _request = _pdm_mpi_2_mpi_request(*request);
  int code = MPI_Wait(&_request, MPI_STATUS_IGNORE);

  free(mpi_request[*request]);
  mpi_request[*request] = NULL;
  n_mpi_request += -1;
  *request = PDM_MPI_COMM_NULL;

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *
 *----------------------------------------------------------------------------*/


int PDM_MPI_Type_create_hindexed (int count,
                              const int array_of_blocklengths[],
                              const PDM_MPI_Aint array_of_displacements[],
                              PDM_MPI_Datatype oldtype,
                              PDM_MPI_Datatype *newtype)
{
  MPI_Datatype mpi_newtype;
  MPI_Aint *_array_of_displacements = malloc (sizeof(MPI_Aint) * count);

  for (int i = 0; i < count; i++) {
    _array_of_displacements[i] = array_of_displacements[i];
  }

  int code = MPI_Type_create_hindexed(count,
                               array_of_blocklengths,
                               _array_of_displacements,
                               _pdm_mpi_2_mpi_datatype(oldtype),
                               &mpi_newtype);

  *newtype = _mpi_2_pdm_mpi_datatype(mpi_newtype);
  free (_array_of_displacements);
  return _mpi_2_pdm_mpi_err(code);

}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype)
{
  MPI_Datatype mpi_type = _pdm_mpi_2_mpi_datatype(*datatype);
  int code =  MPI_Type_commit(&mpi_type);
  *datatype = _mpi_2_pdm_mpi_datatype(mpi_type);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype)
{
  MPI_Datatype mpi_type = _pdm_mpi_2_mpi_datatype(*datatype);
  int code = MPI_Type_free(&mpi_type);
  free(mpi_datatype[*datatype]);
  mpi_datatype[*datatype] = NULL;
  *datatype = PDM_MPI_DATATYPE_NULL;
  n_mpi_datatype += -1;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm)
{

  /* Conversion Fortran vers C */

  MPI_Comm _mpi_comm = MPI_Comm_f2c(comm);
  PDM_MPI_Comm c_comm = _mpi_2_pdm_mpi_comm(_mpi_comm);
  return c_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm)
{

  /* Conversion Fortran vers C */

  MPI_Comm _mpi_comm = _pdm_mpi_2_mpi_comm(comm);
  PDM_MPI_Fint f_comm = (PDM_MPI_Fint) MPI_Comm_c2f(_mpi_comm);
  return f_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scatter (wrapping de la fonction MPI_Scatter)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scatter(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                int root, PDM_MPI_Comm comm)
{
  int code = MPI_Scatter(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                         recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                         root, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Barrier(PDM_MPI_Comm comm)
{
  int code =  MPI_Barrier(_pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Bcast (wrapping de la fonction MPI_Bcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Bcast(void *buffer, int count, PDM_MPI_Datatype datatype,
              int root, PDM_MPI_Comm comm)
{
  int code = MPI_Bcast(buffer,
                       count,
                       _pdm_mpi_2_mpi_datatype(datatype),
                       root,
                       _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                  void *recvbuf, int recvcount,
                  PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code =  MPI_Allgather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                            recvbuf, recvcount,
                            _pdm_mpi_2_mpi_datatype(recvtype),
                            _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgatherv (wrapping de la fonction MPI_Allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts,
                   int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Allgatherv(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                            recvbuf, recvcounts, displs,
                            _pdm_mpi_2_mpi_datatype(recvtype), _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce (wrapping de la fonction MPI_Reduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   PDM_MPI_Datatype datatype, PDM_MPI_Op op,
		   int root, PDM_MPI_Comm comm)
{
  int code = MPI_Reduce(sendbuf, recvbuf, count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           mpi_op[op], root,
                           _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allreduce (wrapping de la fonction MPI_Allreduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  int code = MPI_Allreduce(sendbuf, recvbuf, count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           mpi_op[op],
                           _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scan (wrapping de la fonction MPI_Scan)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  int code = MPI_Scan(sendbuf, recvbuf, count,
                      _pdm_mpi_2_mpi_datatype(datatype),
                      mpi_op[op],
                      _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


int PDM_MPI_Iscan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm,
             PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Iscan(sendbuf, recvbuf, count,
                      _pdm_mpi_2_mpi_datatype(datatype),
                      mpi_op[op],
                      _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoall (wrapping de la fonction MPI_Alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Alltoall(sendbuf, sendcount,
                          _pdm_mpi_2_mpi_datatype(sendtype),
                          recvbuf, recvcount,
                          _pdm_mpi_2_mpi_datatype(recvtype),
                          _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Ialltoall(sendbuf, sendcount,
                          _pdm_mpi_2_mpi_datatype(sendtype),
                          recvbuf, recvcount,
                          _pdm_mpi_2_mpi_datatype(recvtype),
                          _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv (wrapping de la fonction MPI_Alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                      PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                      int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Alltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           _pdm_mpi_2_mpi_datatype(sendtype),
                           recvbuf,
                           recvcounts,
                           rdispls,
                           _pdm_mpi_2_mpi_datatype(recvtype),
                           _pdm_mpi_2_mpi_comm(comm));

  return _mpi_2_pdm_mpi_err(code);
}


int PDM_MPI_Alltoallv_l(void *sendbuf, int *sendcounts, size_t *sdispls,
                      PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                      size_t *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_SUCCESS;

  int size;
  MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), &size);

  INT_MAX;
  int large = 0;
  for (int i = 0; i < size; i++) {
    if ((sdispls[i] > INT_MAX) || (rdispls[i] > INT_MAX)) {
      large = 1;
    }
  }

  if (!large) {

    int *_sdispls = malloc(sizeof(int) * size);
    int *_rdispls = malloc(sizeof(int) * size);

    for (int i = 0; i < size; i++) {
      _sdispls[i] = (int) sdispls[i];
      _rdispls[i] = (int) rdispls[i];
    }

    MPI_Alltoallv(sendbuf,
                  sendcounts,
                  _sdispls,
                  _pdm_mpi_2_mpi_datatype(sendtype),
                  recvbuf,
                  recvcounts,
                  _rdispls,
                  _pdm_mpi_2_mpi_datatype(recvtype),
                  _pdm_mpi_2_mpi_comm(comm));

    free (_sdispls);
    free (_rdispls);
  }

  else {

    MPI_Request *request_r = malloc(sizeof(MPI_Request) * size);
    MPI_Request *request_s = malloc(sizeof(MPI_Request) * size);

    int size_sendType;
    MPI_Type_size(_pdm_mpi_2_mpi_datatype(sendtype), &size_sendType);

    int size_recvType;
    MPI_Type_size(_pdm_mpi_2_mpi_datatype(recvtype), &size_recvType);

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recvType);
        code = MPI_Irecv(buf, recvcounts[i], _pdm_mpi_2_mpi_datatype(recvtype), i,
                         0, _pdm_mpi_2_mpi_comm(comm), request_r + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
      if (sendcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_sendType);
        code = MPI_Issend(buf, sendcounts[i], _pdm_mpi_2_mpi_datatype(sendtype), i,
                          0, _pdm_mpi_2_mpi_comm(comm), request_s + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
    }

    if (code != MPI_SUCCESS) {
      return _mpi_2_pdm_mpi_err(code);
    }

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        code = MPI_Wait(request_r + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (code != MPI_SUCCESS) {
      return _mpi_2_pdm_mpi_err(code);
    }

    for (int i = 0; i < size; i++) {
      if (sendcounts[i] != 0) {
        code = MPI_Wait(request_s + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    free (request_r);
    free (request_s);
  }

  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoallv (wrapping de la fonction MPI_Ialltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                       PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                       int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                       PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Ialltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           _pdm_mpi_2_mpi_datatype(sendtype),
                           recvbuf,
                           recvcounts,
                           rdispls,
                           _pdm_mpi_2_mpi_datatype(recvtype),
                           _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request(_mpi_request);

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen)
{
   int code = MPI_Error_string(mpi_err[errorcode], string, resultlen);
   return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank)
{
  int code = MPI_Comm_rank(_pdm_mpi_2_mpi_comm(comm), rank);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size)
{
  int code = MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), size);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_get_max_error_string(void)
{
  return MPI_MAX_ERROR_STRING;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_free(PDM_MPI_Comm *comm)
{
 int code = 0;
  if ((*comm != PDM_MPI_COMM_NULL) || (*comm != PDM_MPI_COMM_WORLD)) {

    MPI_Comm mpi_comm_loc = _pdm_mpi_2_mpi_comm(*comm);
    code = MPI_Comm_free(&mpi_comm_loc);

    free(mpi_comm[*comm]);
    mpi_comm[*comm] = NULL;
    n_mpi_comm += -1;
    return _mpi_2_pdm_mpi_err(code);
  }
  *comm = PDM_MPI_COMM_NULL;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm)
{
  MPI_Comm _newcomm;
  int code = MPI_Comm_split(_pdm_mpi_2_mpi_comm(comm), color, key, &_newcomm);
  *newcomm = _mpi_2_pdm_mpi_comm(_newcomm);
  return _mpi_2_pdm_mpi_err(code);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
