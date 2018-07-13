#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#if defined(_WIN32) && !defined(__NUTC__)
# include <io.h>
# define ACCESS _access
# define UNLINK _unlink
#else
# include <unistd.h>
# define ACCESS access
# define UNLINK unlink
#endif
#include <errno.h>

#define CGIO_MODE_READ   0
#define CGIO_MODE_WRITE  1
#define CGIO_MODE_MODIFY 2

#define CGIO_FILE_NONE   0
#define CGIO_FILE_ADF    1
#define CGIO_FILE_HDF5   2
#define CGIO_FILE_ADF2   3
#define CGIO_FILE_PHDF5  4

/* currently these are the same as for ADF */

#define CGIO_MAX_DATATYPE_LENGTH  2
#define CGIO_MAX_DIMENSIONS      12
#define CGIO_MAX_NAME_LENGTH     32
#define CGIO_MAX_LABEL_LENGTH    32
#define CGIO_MAX_VERSION_LENGTH  32
#define CGIO_MAX_DATE_LENGTH     32
#define CGIO_MAX_ERROR_LENGTH    80
#define CGIO_MAX_LINK_DEPTH     100
#define CGIO_MAX_FILE_LENGTH   1024
#define CGIO_MAX_LINK_LENGTH   4096

/* these are the cgio error codes */

#define CGIO_ERR_NONE          0
#define CGIO_ERR_BAD_CGIO     -1
#define CGIO_ERR_MALLOC       -2
#define CGIO_ERR_FILE_MODE    -3
#define CGIO_ERR_FILE_TYPE    -4
#define CGIO_ERR_NULL_FILE    -5
#define CGIO_ERR_TOO_SMALL    -6
#define CGIO_ERR_NOT_FOUND    -7
#define CGIO_ERR_NULL_PATH    -8
#define CGIO_ERR_NO_MATCH     -9
#define CGIO_ERR_FILE_OPEN   -10
#define CGIO_ERR_READ_ONLY   -11
#define CGIO_ERR_NULL_STRING -12
#define CGIO_ERR_BAD_OPTION  -13
#define CGIO_ERR_FILE_RENAME -14
#define CGIO_ERR_TOO_MANY    -15
#define CGIO_ERR_DIMENSIONS  -16
#define CGIO_ERR_BAD_TYPE    -17
#define CGIO_ERR_NOT_HDF5    -18

#include <stdio.h>

int cgio_check_file (const char *filename, int *file_type)
{
    int n;
    char buf[32];
    FILE *fp;
    static const char *HDF5sig = "\211HDF\r\n\032\n";
    struct stat st;

    if (ACCESS (filename, 0) || stat (filename, &st) ||
        S_IFREG != (st.st_mode & S_IFREG)) {
        return CGIO_ERR_NOT_FOUND;
    }

    *file_type = CGIO_FILE_NONE;
    if (NULL == (fp = fopen (filename, "rb"))) {
        if (errno == EMFILE)
            return CGIO_ERR_TOO_MANY;
        return CGIO_ERR_FILE_OPEN;
    }
    fread (buf, 1, sizeof(buf), fp);
    buf[sizeof(buf)-1] = 0;
    fclose (fp);

    /* check for ADF */
    if (0 == strncmp (&buf[4], "ADF Database Version", 20)) {
        *file_type = CGIO_FILE_ADF;
        return CGIO_ERR_NONE;
    }

    /* check for HDF5 */

    for (n = 0; n < 8; n++) {
        if (buf[n] != HDF5sig[n]) break;
    }
    if (n == 8) {
        *file_type = CGIO_FILE_HDF5;
        return CGIO_ERR_NONE;
    }

    return CGIO_ERR_FILE_TYPE;
}

/* Cette fonction a ete hackee, elle ne leve plus d'erreurs */
int cgio_find_file (const char *parentfile, const char *filename,
    int file_type, int max_path_len, char *pathname)
{
    int n, size, len, type;
    char *p, *s;

    if (filename == NULL || !*filename)
        return CGIO_ERR_NULL_FILE;
    size = max_path_len - 1 - (int)strlen(filename);
    if (size < 0) return CGIO_ERR_TOO_SMALL;

    /* full path */

    if (*filename == '/'
#ifdef _WIN32
        || *filename == '\\' || *(filename+1) == ':'
#endif
    ) {
        if (cgio_check_file(filename, &type) == CGIO_ERR_NONE &&
           (file_type == CGIO_FILE_NONE || file_type == type)) {
            strcpy(pathname, filename);
            return CGIO_ERR_NONE;
        }
        //if (get_error() == CGIO_ERR_TOO_MANY)
        //    return CGIO_ERR_TOO_MANY;
        return CGIO_ERR_NOT_FOUND;
    }

    /* check relative to parent's directory */

    if (parentfile != NULL && *parentfile && (int)strlen(parentfile) < max_path_len-1) {
        strcpy(pathname, parentfile);
        p = strrchr(pathname, '/');
#ifdef _WIN32
        if (p == NULL) p = strrchr(pathname, '\\');
#endif
        if (p != NULL) {
            *++p = 0;
            if ((int)strlen(pathname) <= size) {
              strcpy(p, filename);
              if (cgio_check_file(pathname, &type) == CGIO_ERR_NONE &&
                    (file_type == CGIO_FILE_NONE || file_type == type))
                return CGIO_ERR_NONE;
            }
        }
    }
    
    /* check current directory */

    if (cgio_check_file(filename, &type) == CGIO_ERR_NONE &&
       (file_type == CGIO_FILE_NONE || file_type == type)) {
        strcpy(pathname, filename);
        return CGIO_ERR_NONE;
    }

        /* check file type environment variable */

    if (file_type == CGIO_FILE_ADF || file_type == CGIO_FILE_ADF2)
        p = getenv ("ADF_LINK_PATH");
#ifdef BUILD_HDF5
    else if (file_type == CGIO_FILE_HDF5 ||
             file_type == CGIO_FILE_PHDF5)
        p = getenv ("HDF5_LINK_PATH");
#endif
    else
        p = NULL;
    while (p != NULL && *p) {
#ifdef _WIN32
        if (NULL == (s = strchr (p, ';')))
#else
        if (NULL == (s = strchr (p, ':')))
#endif
            len = (int)strlen(p);
        else
            len = (int)(s++ - p);
        if (len) {
            if (len > size) return CGIO_ERR_TOO_SMALL;
            strncpy (pathname, p, len);
#ifdef _WIN32
            for (n = 0; n < len; n++) {
                if (*p == '\\') *p = '/';
            }
#endif
            p = pathname + len;
            if (*(p-1) != '/')
                *p++ = '/';
            strcpy (p, filename);
            if (cgio_check_file(pathname, &type) == CGIO_ERR_NONE &&
               (file_type == CGIO_FILE_NONE || file_type == type))
                return CGIO_ERR_NONE;
        }
        p = s;
    }

    /* check $CGNS_LINK_PATH environment variable */

    p = getenv ("CGNS_LINK_PATH");
    while (p != NULL && *p) {
#ifdef _WIN32
        if (NULL == (s = strchr (p, ';')))
#else
        if (NULL == (s = strchr (p, ':')))
#endif
            len = (int)strlen(p);
        else
            len = (int)(s++ - p);
        if (len) {
            if (len > size) return CGIO_ERR_TOO_SMALL;
            strncpy (pathname, p, len);
#ifdef _WIN32
            for (n = 0; n < len; n++) {
                if (*p == '\\') *p = '/';
            }
#endif
            p = pathname + len;
            if (*(p-1) != '/')
                *p++ = '/';
            strcpy (p, filename);
            if (cgio_check_file(pathname, &type) == CGIO_ERR_NONE &&
               (file_type == CGIO_FILE_NONE || file_type == type))
                return CGIO_ERR_NONE;
        }
        p = s;
    }


    // Je n'ai pas integre la recherche dans differents paths

    return CGIO_ERR_NOT_FOUND;
}
