/*============================================================================
 * Manage case files associated with the EnSight Gold writer
 *============================================================================*/
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>  /* toupper() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_writer.h"
#include "pdm_writer_ensight_case.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Time set entry structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int           n_time_values;   /* Number of time step values */
  double       *time_value;      /* Time step values */

} PDM_writer_ensight_case_time_t;

/*----------------------------------------------------------------------------
 * Variable entry structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char         *name;            /* Variable name */
  char         *case_line;       /* Line in case file */
  PDM_writer_statut_t   time_dep;        /* time dependant */
  char         *file_name_base;  /* file name base */
  char         *file_name;       /* file name base */

  int           dim;             /* Associated dimension: 0 (constant), 1
                                    1 (scalar), 3 (vector), 6 (symmetrical
                                    tensor), or 9 (asymmetrical tensor) */
  PDM_writer_var_loc_t  loc;  /* variable at nodes, elements, or particles */

} PDM_writer_ensight_case_var_t;

/*----------------------------------------------------------------------------
 * EnSight case file structure
 *----------------------------------------------------------------------------*/

struct _PDM_writer_ensight_case_t {

  char          *name;              /* Case name */
  char          *case_file_name;    /* Case file name */

  char          *file_name_prefix;  /* File name prefix */
  int            dir_name_length;   /* Associated directory name length
                                       (may be 0); index in file_name_prefix
                                       corresponding to the base file name */

  char          *geom_file_name_base; /* Geometry file name */
  char          *geom_file_name;      /* Geometry file name */

  PDM_writer_ensight_case_time_t  *time_set;     /* Time Set entries */

  int                      n_vars;       /* Number of variables */
  PDM_writer_ensight_case_var_t   **var;          /* Variable entries */

  PDM_writer_topologie_t   time_dependency;    /* Mesh time dependency */

} _PDM_writer_ensight_case_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const int _l_max_chaine_ens = 1024;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add a new time step number and value to a time set if necessary:
 * if the corresponding physical time is not present in the structure,
 * the corresponding elements are added.
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   time_step  <-- number of time step to add
 *   time_value <-- associated time value
 *
 * returns:
 *   0 if no time was added, 1 if a new time was added
 *----------------------------------------------------------------------------*/

static int
_add_time(PDM_writer_ensight_case_time_t  *const time_set,
          const double                   time_value)
{
  /*const char time_value_err_string[] =
    "The time value associated with time step <%d> equals <%g>,\n"
    "but time value <%g> has already been associated with this time step.\n";*/

  /* Finally, add a new time step and value if necessary */

  time_set->n_time_values += 1;
  
  time_set->time_value = (double *) realloc (time_set->time_value, 
                                             time_set->n_time_values * sizeof(double));
  
  time_set->time_value[time_set->n_time_values - 1] = time_value;

  return 1;
}

/*----------------------------------------------------------------------------
 * Add a new variable entry
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   name       <-- variable name
 *   dimension  <-- variable dimension (0: constant, 1: scalar, 3: vector,
 *                  6: symmetrical tensor, 9: asymmetrical tensor)
 *   location   <-- variable definition location (nodes, elements, or particles)
 *   time_set   <-- associated time set index
 *----------------------------------------------------------------------------*/

static void
_add_var(PDM_writer_ensight_case_t       *const this_case,
         const char              *const name,
         const int                      dimension,
         const PDM_writer_statut_t              time_dep,
         const PDM_writer_var_loc_t             location)
{
  char line[1024], description[50];
  int i;

  PDM_writer_ensight_case_var_t  *var;

  size_t l = strlen(name);

  this_case->n_vars += 1;
  var = (PDM_writer_ensight_case_var_t *) malloc (sizeof(PDM_writer_ensight_case_var_t));

  var->name = (char *) malloc( (l + 1) * sizeof(char));
  strcpy(var->name, name);

  /* Create description (49 chars max) */

  if (l > 49)
      l = 49;

  strncpy(description, name, l);
  description[l] = '\0';

  /* Some characters not allowed in format, replaced by '_' */

  for (i = 0 ; i < l ; i++) {
    switch (description[i]) {
    case '(':
    case ')':
    case ']':
    case '[':
    case '+':
    case '-':
    case '@':
    case ' ':
    case '\t':
    case '!':
    case '#':
    case '*':
    case '^':
    case '$':
    case '/':
      description[i] = '~';
      break;
    default:
       break;
    }
  }

  /* Assign time set to obtain file index, and dimension and location
     before case line creation */

  var->dim      = dimension;
  var->loc      = location;
  var->time_dep = time_dep;

  /* Create associated case file line, up to file name
     (which may depend on the number of remaining characters,
     so as to avoid going beyond 1024 characters if possible) */

  switch(var->dim) {
  case 0:
    strcpy(line, "constant per case file: ");
    break;
  case 1:
    strcpy(line, "scalar per ");
    break;
  case 3:
    strcpy(line, "vector per ");
    break;
  case 6:
    strcpy(line, "tensor symm per ");
    break;
  case 9:
    strcpy(line, "tensor asym per ");
    break;
  }

  if (var->dim > 0) {
    switch(var->loc) {
    case PDM_WRITER_VAR_SOMMETS:
      strcat(line, "node:    ");
      break;
    case PDM_WRITER_VAR_ELEMENTS:
      strcat(line, "element: ");
      break;
    case PDM_WRITER_VAR_PARTICULES:
      strcat(line, "measured node: ");
      break;
    }
  }

  l = strlen(line); /* At this stage, l = 31 at most */

  if (var->time_dep == 1)
    sprintf(line + l, "%d ", 1);
  else
    strcat(line, "  ");

  l = strlen(line);  /* At this stage, l = 35 at most with
                        time set number < 100 (EnSight maximum:
                        16, apparently only for measured data) */

  sprintf(line + l, "%32s ", description); /* Description max 49 chars,
                                              (32 recommended for compatibility
                                              with other formats)
                                              so 1024 char line limit not
                                              exceeded here */

  for (l = strlen(line) ; l < 61 ; l++)
    line[l] = ' ';
  line[l] = '\0'; /* Line length = 35 + 49 + 1 = 85 max at this stage,
                     usually 31 + 32 + 1 = 64 */

  /* Create (current) file name. */

  size_t prefix_len =   strlen(this_case->file_name_prefix)
               - this_case->dir_name_length + 1;
  
  size_t base_len = strlen(name);

  var->file_name_base = (char *) malloc (
             (  this_case->dir_name_length + prefix_len
                + base_len + 1) * sizeof(char));
  sprintf(var->file_name_base, "%s.", this_case->file_name_prefix);

  strcat(var->file_name_base, name);
  for (size_t i1 = this_case->dir_name_length + prefix_len ;
       i1 < this_case->dir_name_length + prefix_len + base_len ;
       i1++) {
    switch (var->file_name_base[i1]) {
    case '@':
    case ' ':
    case '\t':
    case '!':
    case '#':
    case '*':
    case '^':
    case '$':
    case '/':
      var->file_name_base[i1] = '_';
    default:
      var->file_name_base[i1] = (char) tolower(var->file_name_base[i1]);
    }
  }

  var->file_name = NULL;

  /* Now we may finish associated case file line */

  var->case_line = (char *) malloc (
             (  strlen(line) + strlen(var->file_name_base) + 6
                - this_case->dir_name_length + 1) * sizeof(char));

  strcpy(var->case_line, line);
  strcat(var->case_line, var->file_name_base + this_case->dir_name_length);

  /* Replace current time index by wildcards */

  if (var->time_dep == 1)
    strcat(var->case_line, ".*****");

  /* Finally, associate variable entry in case file */

  if (strlen(var->case_line) > 1024) {
    PDM_printf ("Line of the EnSight case file \"%s\"\n"
            "for variable \"%s\",\n"
            "exceeds 1024 characters, so this file must be edited and variable\n"
            "descriptions or referenced files renamed so as to be readable.\n",
            this_case->case_file_name, name);
  }

  this_case->var = 
    (PDM_writer_ensight_case_var_t **) realloc (this_case->var, this_case->n_vars * sizeof(PDM_writer_ensight_case_var_t *));

  this_case->var[this_case->n_vars - 1] = var;
}

/*----------------------------------------------------------------------------
 * Remove variable entries
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *----------------------------------------------------------------------------*/

static void
_del_vars(PDM_writer_ensight_case_t  *const this_case)
{
  int i;

  for (i = 0 ; i < this_case->n_vars ; i++) {

    PDM_writer_ensight_case_var_t  *var = this_case->var[i];

    free(var->name);
    free(var->case_line);
    free(var->file_name_base);
    if (var->file_name != NULL)
      free(var->file_name);

    free(var);

  }

  free(this_case->var);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Create a new case file structure.
 *
 * parameters:
 *   name            <-- case name
 *   restart         <-- if restart == 1, case file is read
 *   dir_prefix      <-- associated local or absolute directory name
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to new case file structure
 *----------------------------------------------------------------------------*/

PDM_writer_ensight_case_t *
PDM_writer_ensight_case_cree 
(
const char                   *const name,
const int                           restart,
const char                   *const dir_prefix,
const PDM_writer_topologie_t                time_dependency
)
{
  size_t  i, name_len, prefix_len;

  PDM_writer_ensight_case_t   *this_case = NULL;

  /* Create and initialize structure */

  this_case = (PDM_writer_ensight_case_t *) malloc (sizeof(PDM_writer_ensight_case_t));

  /* Initialize base name and partial file names */

  this_case->name = (char *) malloc((strlen(name) + 1) * sizeof(char));
  strcpy(this_case->name, name);
  name_len = strlen(name);

  for (i = 0 ; i < name_len ; i++) {
    if (   (this_case->name[i] == ' ')
           || (this_case->name[i] == '\t'))
      this_case->name[i] = '_';
  }

  if (dir_prefix != NULL)
    prefix_len = strlen(dir_prefix) + 1;
  else
    prefix_len = 0;

  this_case->dir_name_length = (int) prefix_len;

  this_case->case_file_name = (char *) malloc((prefix_len + name_len + 6) * sizeof(char));
  if (dir_prefix != NULL) {
    strcpy(this_case->case_file_name, dir_prefix);
    strcat(this_case->case_file_name, "/");
  }
  else
    this_case->case_file_name[0] = '\0';

  for (i = 0 ; i < name_len ; i++)
    this_case->case_file_name[prefix_len + i] = (char) toupper(name[i]);
  this_case->case_file_name[prefix_len + name_len] = '\0';

  this_case->file_name_prefix = (char *) malloc((strlen(this_case->case_file_name) + 1) *
                                                sizeof(char));
  strcpy(this_case->file_name_prefix, this_case->case_file_name);
  for (i = 0 ; i < name_len ; i++)
    this_case->file_name_prefix[prefix_len + i]
      = (char) tolower((int) this_case->case_file_name[prefix_len + i]);

  strcat(this_case->case_file_name, ".case");

  /* Initialize other members */

  this_case->time_set = NULL;

  this_case->n_vars = 0;
  this_case->var = NULL;

  this_case->time_dependency = time_dependency;

  /* Geometry file name (after time dependency) */

  this_case->geom_file_name_base = NULL;
  this_case->geom_file_name      = NULL;

  char extension[5] = ".geo";
 
  this_case->geom_file_name_base = (char *) malloc((strlen(this_case->file_name_prefix) 
                                                    + strlen(extension) + 1 ) * sizeof(char));

  strcpy(this_case->geom_file_name_base, this_case->file_name_prefix);
  strcat(this_case->geom_file_name_base, extension);

  /* Status information */

  if (restart) {

    FILE *f = fopen(this_case->case_file_name, "r");

    int               geom_timeset = -1;
    int               geom_fileset = -1;
    size_t            ind = -1;
    
    char              ligne[_l_max_chaine_ens];
    char              nom_fic_geo_base[_l_max_chaine_ens];

    int              fmt_ensight      = 0;
    
    char             *retval = NULL;
    
    int          *time_set_num = NULL;
    int          *time_set_n_step = NULL;

    if (f != NULL) {

      /* Vérification du format */
      
      do {
        retval = fgets(ligne, _l_max_chaine_ens, f);
      } while (retval != NULL && strncmp(ligne, "FORMAT", strlen("FORMAT")) != 0);

      if (retval != NULL) {

        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "type:", strlen("type:")) != 0);
        
      }

      if (retval != NULL) {

        for (ind = strlen("type:");
             ligne[ind] != '\0' && (ligne[ind] == ' ' || ligne[ind] == '\t');
             ind++);

        if (strncmp(ligne + ind, "ensight", strlen("ensight")) == 0) {

          fmt_ensight = 1;

          ind += strlen("ensight");
          while (ligne[ind] != '\0' && (ligne[ind] == ' ' || ligne[ind] == '\t'))
            ind++;
  
        }
      }
 
      if (fmt_ensight == 0) {
        PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : "
                "File \"%s\" does not seem to be a valid\n"
                "EnSight 6 or Gold case file.",
                __FILE__, __LINE__, this_case->case_file_name);
      }

      /* Recherche des infos sur le fichier géométrique */

      do {
        retval = fgets(ligne, _l_max_chaine_ens, f);
      } while (retval != NULL && strncmp(ligne, "GEOMETRY", strlen("GEOMETRY")) != 0);

      if (retval != NULL) {

        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "model:", strlen("model:")) != 0);

      }

      if (retval != NULL) {

        /* La rubrique model: contient deux numéros optionnels (numéro de pas de
           temps et de jeux de fichiers), le nom de base du ou des fichiers
           géométriques, et éventuellement l'option "change_coords_only" */

        if (sscanf(ligne, "%*s %d %d %s",
                   &geom_timeset, &geom_fileset, nom_fic_geo_base) != 3) {
          if (sscanf(ligne, "%*s %d %s",
                     &geom_timeset, nom_fic_geo_base) != 2) {
            if (sscanf(ligne, "%*s %s",
                       nom_fic_geo_base) != 1) {
              PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : "
                      "The \"%s\" case file does not seem to\n"
                      "indicate a geometry file",
                      __FILE__, __LINE__,
                      this_case->case_file_name);
              abort();
            }
          }
        }
      }

      if (((time_dependency != PDM_WRITER_TOPO_CONSTANTE) && (geom_timeset == -1)) ||
          ((time_dependency == PDM_WRITER_TOPO_CONSTANTE) && (geom_timeset != -1))) {
        PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : Inconsistency geom time dependency between "
                "The \"%s\" case and the function argument\n",
                __FILE__, __LINE__,
                this_case->case_file_name);
        abort();
      }

      /* Recherche des infos sur les fichiers variables */

      do {
        retval = fgets(ligne, _l_max_chaine_ens, f);
      } while (retval != NULL && strncmp(ligne, "VARIABLE", strlen("VARIABLE")) != 0);

      while (1) {

        if (retval != NULL) {
        
          do {
            retval = fgets(ligne, _l_max_chaine_ens, f);
          } while (retval != NULL && (strncmp(ligne, "constant", strlen("constant")) != 0 &&
                                      strncmp(ligne, "scalar", strlen("scalar")) != 0 &&
                                      strncmp(ligne, "vector", strlen("vector")) != 0 &&
                                      strncmp(ligne, "tensor", strlen("tensor")) != 0 &&
                                      strncmp(ligne, "complex", strlen("complex")) != 0 &&
                                      strncmp(ligne, "TIME", strlen("TIME")) && 0 ));

          if (retval != NULL && strncmp(ligne, "TIME", strlen("TIME")) != 0 )
            break;
          else {

            int index_loc;
            int var_dim;
          
            if (strncmp(ligne, "constant per case file:", strlen("constant per case file:")) == 0) {
              index_loc = strlen("constant per case file:");
              var_dim = PDM_WRITER_VAR_CSTE;
            }
            else if (strncmp(ligne, "scalar per ", strlen("scalar per ")) == 0) {
              index_loc = strlen("scalar per ");
              var_dim = PDM_WRITER_VAR_SCALAIRE;
            }
            else if (strncmp(ligne, "vector per ", strlen("vector per ")) == 0) {
              index_loc = strlen("vector per ");
              var_dim = PDM_WRITER_VAR_VECTEUR;
            }
            else if (strncmp(ligne, "tensor symm per ", strlen("tensor symm per ")) == 0) {
              index_loc = strlen("tensor symm per ");
              var_dim = PDM_WRITER_VAR_TENSEUR_SYM;
            }
            else if (strncmp(ligne, "tensor asymm per ", strlen("tensor asymm per ")) == 0) {
              index_loc = strlen("tensor asymm per ");
              var_dim = PDM_WRITER_VAR_TENSEUR_ASYM;
            }
            else {
              PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : "
                      "The dimension \"%s\" is not yet implemented", 
                      __FILE__, __LINE__,
                      ligne);
              abort();
            }

            char *ligne_ss_dim = ligne + index_loc;
            int index_dim;
            PDM_writer_var_loc_t var_loc;

            if (strncmp(ligne_ss_dim, "node:    ", strlen("node:    ")) == 0) {
              index_dim = strlen("node:    ");
              var_loc = PDM_WRITER_VAR_SOMMETS;
            }
            else if (strncmp(ligne_ss_dim, "element: ", strlen("element: ")) == 0) {
              index_dim = strlen("element: ");
              var_loc = PDM_WRITER_VAR_ELEMENTS;
            }
            else if (strncmp(ligne_ss_dim, "measured node: ", strlen("measured node: ")) == 0) {
              index_dim = strlen("measured node: ");
              var_loc = PDM_WRITER_VAR_PARTICULES;
            }
            else {
              PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : "
                      "The location \"%s\" is not yet implemented", 
                      __FILE__, __LINE__,
                      ligne_ss_dim);
              abort();
            }

            char nom_var[_l_max_chaine_ens];
            char *ligne_ss_loc = ligne_ss_dim + index_dim;
          
            int var_file_set = -1;
            int var_time_set = -1;
            
            if (sscanf(ligne_ss_loc, "%d %d %s", &var_time_set, &var_file_set, nom_var) != 3) {
              if (sscanf(ligne_ss_loc, "%d %s", &var_time_set, nom_var) != 2) {
                if (sscanf(ligne_ss_loc, "%s", nom_var) != 1) {
                  PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d\n",
                          __FILE__, __LINE__);
                  abort();
                }
              }
            }
            
            PDM_writer_statut_t time_dep = PDM_WRITER_OFF;

            if (var_time_set == 1)
              time_dep = PDM_WRITER_ON;

            _add_var(this_case,
                     nom_var,
                     var_dim,
                     time_dep,
                     var_loc);
          }
        }

        else
          break;
       
      }

      /* Recherche des infos sur les rubriques Time : il n'y a qu'u time set*/

      do {
        retval = fgets(ligne, _l_max_chaine_ens, f);
      }  while (retval != NULL && strncmp(ligne, "TIME", strlen("TIME")) != 0);
        
      if (retval != NULL) {

        this_case->time_set = (PDM_writer_ensight_case_time_t *) malloc(sizeof(PDM_writer_ensight_case_time_t));
        this_case->time_set->time_value = NULL;
        this_case->time_set->n_time_values = 0;

        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "time set:", strlen("time set:")) != 0);
              
        if (retval != NULL) {
          if (sscanf(ligne, "%*s %d", time_set_num) != 1) {
            PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d\n",
                    __FILE__, __LINE__);
            abort();
          }
        }
            
        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "number of steps:", strlen("number of steps:")) != 0);

        if (retval != NULL) {
          if (sscanf(ligne, "%*s %d", time_set_n_step) != 1) {
            PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d\n",
                    __FILE__, __LINE__);
            abort();
          }
        }

        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "filename start number: 1", strlen("filename start number: 1")) != 0);

        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "filename increment:    1", strlen("filename increment:    1")) != 0);
              
        do {
          retval = fgets(ligne, _l_max_chaine_ens, f);
        } while (retval != NULL && strncmp(ligne, "time values:", strlen("time values:")) != 0);
        
        if (retval != NULL) {
        
          for (int j = 0; j < *time_set_n_step; j++) {
          
            retval = fgets(ligne, _l_max_chaine_ens, f);
          
            if (retval == NULL) {
            
              PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : "
                      "The dimension \"%s\" is not yet implemented", 
                      __FILE__, __LINE__,
                      ligne);
              abort();
            }
            else {
              double time_value = atof(ligne);
              _add_time(this_case->time_set, time_value); 
            }
          }
          retval = fgets(ligne, _l_max_chaine_ens, f);
        }
      }
    }
  }      

  /* Return new case structure */

  return this_case;
}


 /*----------------------------------------------------------------------------
  * Destroy a case file structure.
  *
  * parameters:
  *   this_case  <-- case structure
  *
  * returns:
  *   NULL pointer
  *----------------------------------------------------------------------------*/

 PDM_writer_ensight_case_t *
 PDM_writer_ensight_case_lib(PDM_writer_ensight_case_t  *this_case)
 {

   /* Free names */

   free(this_case->name);
   free(this_case->case_file_name);
   free(this_case->file_name_prefix);

   free(this_case->geom_file_name_base);
    if (this_case->geom_file_name != NULL)
      free(this_case->geom_file_name);

   /* Free variable entries */

   _del_vars(this_case);

   /* Free time sets */

   if (this_case->time_set->time_value != NULL)
     free(this_case->time_set->time_value);
   free(this_case->time_set);

  /* Free structure and return */

  free(this_case);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

PDM_writer_topologie_t
PDM_writer_ensight_case_geo_time_dep_get(PDM_writer_ensight_case_t  *this_case)
{
  return  this_case->time_dependency;
}


/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight var
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

PDM_writer_statut_t
PDM_writer_ensight_case_var_time_dep_get
(
 PDM_writer_ensight_case_t  *this_case,
 const char         *name
)
{
  int i = 0;
  PDM_writer_ensight_case_var_t  *var = this_case->var[i];
  for (i = 0 ; i < this_case->n_vars ; i++) {
    var = this_case->var[i];
    if (strcmp(var->name, name) == 0) {
      break;
    }
  }

  if (i >= this_case->n_vars) {
    PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : Unknown variable \"%s\":\n\n",
           __FILE__, __LINE__, name);
    abort();
  }
  
  return var->time_dep;
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

char *
PDM_writer_ensight_case_var_file_name_get
(
PDM_writer_ensight_case_t  *this_case, 
const char* name
)
{
  int i = 0;
  PDM_writer_ensight_case_var_t  *var = this_case->var[i];
  for (i = 0 ; i < this_case->n_vars ; i++) {
    var = this_case->var[i];
    if (strcmp(var->name, name) == 0) {
      break;
    }
  }

  if (i >= this_case->n_vars) {
    PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : Unknown variable \"%s\":\n\n",
           __FILE__, __LINE__, name);
    abort();
  }
  
  char *file_name;
  if (var->time_dep == PDM_WRITER_ON) {
    file_name = var->file_name;
  }
  else {
    file_name = var->file_name_base;
  }

  return file_name;
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

char *
PDM_writer_ensight_case_geo_file_name_get
(
PDM_writer_ensight_case_t  *this_case
)
{
  char *name;
  if (this_case->time_dependency != PDM_WRITER_TOPO_CONSTANTE) {
    name = this_case->geom_file_name;
  }
  else {
    name = this_case->geom_file_name_base;
  }
  
  return name;
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/


void
PDM_writer_ensight_case_var_cree
(
PDM_writer_ensight_case_t  *this_case,
const char         *const name,
const PDM_writer_var_dim_t  dimension,
const PDM_writer_statut_t   time_dep,
const PDM_writer_var_loc_t  location
)
{
  _add_var(this_case,
           name,
           (int) dimension,
           time_dep,
           location);
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_case_time_step_add
(
PDM_writer_ensight_case_t  *this_case, 
const double time_value
)
{
  if (this_case->time_set == NULL) {
    this_case->time_set = (PDM_writer_ensight_case_time_t *) malloc(sizeof(PDM_writer_ensight_case_time_t));
    this_case->time_set->time_value = NULL;
    this_case->time_set->n_time_values = 0;
  }

  _add_time(this_case->time_set,
            time_value);
  if (this_case->time_dependency != PDM_WRITER_TOPO_CONSTANTE) {
    
    if (this_case->geom_file_name == NULL) {
      this_case->geom_file_name = (char *) malloc(sizeof(char) * strlen(this_case->geom_file_name_base) + 7);
    } 

    int geom_index = this_case->time_set->n_time_values;
    char extension[7];
    sprintf(extension, ".%05d", geom_index);
    strcpy(this_case->geom_file_name, this_case->geom_file_name_base);
    strcat(this_case->geom_file_name, extension);
  }

  for (int i = 0; i < this_case->n_vars; i++) {

    PDM_writer_ensight_case_var_t  *var = this_case->var[i];
    if (var->time_dep == PDM_WRITER_ON) {
      if (var->file_name == NULL) {
        var->file_name = (char *) malloc(sizeof(char) * strlen(var->file_name_base) + 7);
      } 

      int geom_index = this_case->time_set->n_time_values;
      char extension[7];
      sprintf(extension, ".%05d", geom_index);
      strcpy(var->file_name, var->file_name_base);
      strcat(var->file_name, extension);
    }
  }
}


/*----------------------------------------------------------------------------
 * Write an EnSight Gold case file.
 *
 * This function should only be called by one process in parallel mode.
 *
 * parameters:
 *   this_case  <-- case structure
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_case_write(PDM_writer_ensight_case_t  *const this_case,
                      int                       rank)
{
  int          i, j;
  FILE         *f;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* if (this_case->modified == 0) */
  /*   return; */

  /* this_case->modified = 0; */

  if (rank > 0)
    return;

  /* Open case file (overwrite it if present) */

  f = fopen(this_case->case_file_name, "w");

  if (f == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : Error opening file \"%s\":\n\n",
           __FILE__, __LINE__, this_case->case_file_name);
    abort();
  }

  /* Output FORMAT */

  fprintf(f,
          "FORMAT\n"
          "type: ensight gold\n");

  /* Output geometry */

  fprintf(f,
          "\n"
          "GEOMETRY\n");

  if (this_case->time_dependency == PDM_WRITER_TOPO_CONSTANTE)
    fprintf(f, "model: %s.geo\n",
            this_case->file_name_prefix + this_case->dir_name_length);


  else if (this_case->time_dependency == PDM_WRITER_TOPO_DEFORMABLE)
    fprintf(f, "model: %d %s.geo.*****  change_coords_only\n",
            1,
            this_case->file_name_prefix + this_case->dir_name_length);

  
  else
    fprintf(f, "model: %d %s.geo.*****\n",
            1,
            this_case->file_name_prefix + this_case->dir_name_length);

  /* Output variables */

  if (this_case->n_vars > 0) {

    fprintf(f,
            "\n"
            "VARIABLE\n");

    for (i = 0 ; i < this_case->n_vars ; i++) {
      const PDM_writer_ensight_case_var_t  *var = this_case->var[i];
      fprintf(f, "%s\n", var->case_line);
    }

  }

  if (this_case->time_set != NULL) {
    
    fprintf(f,
            "\n"
            "TIME\n");
    
    const PDM_writer_ensight_case_time_t  *ts = this_case->time_set;
    
    fprintf(f, "time set:              %d\n", 1);
    fprintf(f, "number of steps:       %d\n", ts->n_time_values);
    fprintf(f, "filename start number: 1\n");
    fprintf(f, "filename increment:    1\n");
    fprintf(f, "time values:\n");
    
    for (j = 0 ; j < ts->n_time_values ; j++)
      fprintf(f, "            %g\n", ts->time_value[j]);
  }

  /* Close case file */

  if (fclose(f) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Error in %s line %d : Error closing file \"%s\":\n\n",
           __FILE__, __LINE__, this_case->case_file_name);
    abort();
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
