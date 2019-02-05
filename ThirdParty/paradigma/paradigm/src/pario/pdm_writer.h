#ifndef __PDM_WRITER_H__
#define __PDM_WRITER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Macro permettant de prendre en compte les arguments caches des longueurs
 * de chaine du Fortran (Evite les erreurs de compilation)
 *----------------------------------------------------------------------------*/

#if defined (__uxpv__) 
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Statut
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_OFF,    
  PDM_WRITER_ON

} PDM_writer_statut_t;

/*----------------------------------------------------------------------------
 * Type de topologie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_TOPO_CONSTANTE ,     
  PDM_WRITER_TOPO_DEFORMABLE,     
  PDM_WRITER_TOPO_VARIABLE     

} PDM_writer_topologie_t;

/*----------------------------------------------------------------------------
 * Type d'elements géometriques (It's the same than the type defined into PDM_Mesh_nodal)
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_POINT = PDM_MESH_NODAL_POINT,     
  PDM_WRITER_BAR2 = PDM_MESH_NODAL_BAR2,     
  PDM_WRITER_TRIA3 = PDM_MESH_NODAL_TRIA3,     
  PDM_WRITER_QUAD4 = PDM_MESH_NODAL_QUAD4,     
  PDM_WRITER_POLY_2D = PDM_MESH_NODAL_POLY_2D,     
  PDM_WRITER_TETRA4 = PDM_MESH_NODAL_TETRA4,     
  PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5,     
  PDM_WRITER_PRISM6 = PDM_MESH_NODAL_PRISM6,     
  PDM_WRITER_HEXA8 = PDM_MESH_NODAL_HEXA8,     
  PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D    

} PDM_writer_elt_geom_t;


typedef struct _PDM_writer_t PDM_writer_t;

typedef struct _PDM_writer_geom_t PDM_writer_geom_t;

typedef struct _PDM_writer_var_t PDM_writer_var_t;



/**
 * \struct PDM_writer_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_t instance 
 *
 */

typedef void (*PDM_writer_fct_t) (PDM_writer_t *pw);


/**
 * \struct PDM_writer_geom_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_geom_t instance 
 *
 */

typedef void (*PDM_writer_geom_fct_t) (PDM_writer_geom_t *geom);


/**
 * \struct PDM_writer_var_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_var_t instance 
 *
 */

typedef void (*PDM_writer_var_fct_t) (PDM_writer_var_t *var);


/*----------------------------------------------------------------------------
 * Format du fichie de sortie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_FMT_BIN,        
  PDM_WRITER_FMT_ASCII       

} PDM_writer_fmt_fic_t;

/*----------------------------------------------------------------------------
 * Types des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_CSTE           = 0,
  PDM_WRITER_VAR_SCALAIRE       = 1,
  PDM_WRITER_VAR_VECTEUR        = 3,
  PDM_WRITER_VAR_TENSEUR_SYM    = 6,
  PDM_WRITER_VAR_TENSEUR_ASYM   = 9

} PDM_writer_var_dim_t;

/*----------------------------------------------------------------------------
 * Localisation des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_SOMMETS,
  PDM_WRITER_VAR_ELEMENTS,
  PDM_WRITER_VAR_PARTICULES

} PDM_writer_var_loc_t;


/*=============================================================================
 * Variables globales
 *============================================================================*/

/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet 
 *
 * parameters :
 *   fmt             <-- Format de sortie
 *   fmt_fic         <-- Binary or ASCII
 *   topologie       <-- Indique le maillage est mobile ou non
 *   st_reprise      <-- Complete les sorties des calculs precedents en reprise
 *   rep_sortie      <-- Repertoire de sortie                  
 *   nom_sortie      <-- Nom de la sortie                       
 *   pdm_mpi_com         <-- Communicateur MSG                      
 *   acces           <-- Type d'acces                           
 *   prop_noeuds_actifs <-- Proportion des noeuds actifs dans les acces au fichier
 *                            *  -1 : tous les processus actifs
 *                            *   1 : un processus par noeud
 *                            * 0 < val < 1 : un processus par noeud actif
 *   options         <-- Options complementaires propres au format sous
 *                       la forme ("nom_1 = val_1 : ... : nom_n = val_n")                          
 *
 * return :
 *                   --> Identificateur de l'objet cree
 *
 *----------------------------------------------------------------------------*/

void 
PROCF (pdm_writer_create_cf, PDM_WRITER_CREATE_CF)
(
const char          *fmt,
const int           *l_fmt,
const int           *fmt_fic,
const int           *topologie,
const int           *st_reprise,
const char          *rep_sortie,
const char          *nom_sortie,
const int           *l_rep_sortie,
const int           *l_nom_sortie,
const PDM_MPI_Fint  *pdm_mpi_comm,   
const int           *acces,
const double        *prop_noeuds_actifs,
const char          *options,
const int           *l_options,
int                 *id_cs     
ARGF_SUPP_CHAINE
);

int
PDM_writer_create
(
const char                  *fmt,
const PDM_writer_fmt_fic_t   fmt_fic,   
const PDM_writer_topologie_t topologie,
const PDM_writer_statut_t    st_reprise,
const char                  *rep_sortie,
const char                  *nom_sortie,
const PDM_MPI_Comm           pdm_mpi_comm,
const PDM_io_acces_t         acces,
const double                 prop_noeuds_actifs,
const char                  *options
);  

/*----------------------------------------------------------------------------
 * Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *
 *
 *----------------------------------------------------------------------------*/

void 
PROCF (pdm_writer_free, PDM_WRITER_FREE)
(
int        *id_cs
);

void  
PDM_writer_free
(
const int   id_cs
);

/*----------------------------------------------------------------------------
 * Debut d'increment
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   delta_t         <-- Delta de temps par rapport au dernier increment
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_step_beg, PDM_WRITER_STEP_BEG)
(
int           *id_cs,
double        *physical_time
);

void
PDM_writer_step_beg
(
const int      id_cs,
const double   physical_time
);

/*----------------------------------------------------------------------------
 * Fin d'increment
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_step_end, PDM_WRITER_STEP_END)
(
int          *id_cs
);

void
PDM_writer_step_end
(
const int     id_cs
);

/*----------------------------------------------------------------------------
 * Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
 *
 * parameters :
 *   id_cs            <-- Identificateur de l'objet cs
 *   nom_geom         <-- Nom de l'objet geometrique
 *   st_decoup_poly2d <-- Active le decoupage des polygones 
 *   st_decoup_poly3d <-- Active le decoupage des polyedres
 *
 * return :
 *                   --> Identificateur de l'objet geom dans cs 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_create_cf, PDM_WRITER_GEOM_CREATE_CF)
(
int           *id_cs,
char          *nom_geom,
int           *st_decoup_poly2d,
int           *st_decoup_poly3d,
int           *l_nom_geom,
int           *n_part,
int           *id_geom
ARGF_SUPP_CHAINE
);

int 
PDM_writer_geom_create
(
const int           id_cs,
const char         *nom_geom,
const PDM_writer_statut_t   st_decoup_poly2d,
const PDM_writer_statut_t   st_decoup_poly3d,
const int           n_part
);


/*----------------------------------------------------------------------------
 * Definition des coordonnees de la partition courante          
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part          <-- Indice de partition
 *   n_som           <-- Nombre de sommets de la partition
 *   coords          <-- Coordonnes des sommets            
 *   numabs          <-- Numerotation absolue des sommets     
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_coord_set, PDM_WRITER_GEOM_COORD_SET)
(
int             *id_cs,
int             *id_geom,  
int             *id_part, 
int             *n_som,  
PDM_real_t       *coords,  
PDM_g_num_t       *numabs
);

void
PDM_writer_geom_coord_set
(
const int        id_cs,
const int        id_geom,  
const int        id_part, 
const int        n_som,  
const PDM_real_t *coords,  
const PDM_g_num_t *numabs
);

/*----------------------------------------------------------------------------
 * Definition des coordonnees des sommets de la partition courante a partir
 *          
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Indice de partition
 *   n_som           <-- Nombre de sommets de la partition
 *   n_som_parent    <-- Nombre de sommets parent
 *   numabs          <-- Numerotation absolue des sommets (size = n_som)    
 *   num_parent      <-- Numerotation des sommets dans la numerotation parente (size = n_som)    
 *   coords_parent   <-- Coordonnes des sommets parents (size = 3 * n_som_parent)            
 *   numabs_parent   <-- Numerotation absolue des sommets parents (size = n_som_parent)    
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_coord_from_parent_set, PDM_WRITER_GEOM_COORD_FROM_PARENT_SET)
(
int             *id_cs,
int             *id_geom,  
int             *id_part, 
int             *n_som,  
int             *n_som_parent,  
PDM_g_num_t       *numabs,
int             *num_parent,
PDM_real_t       *coords_parent,  
PDM_g_num_t       *numabs_parent
);

void
PDM_writer_geom_coord_from_parent_set
(
const int        id_cs,
const int        id_geom,  
const int        id_part, 
const int        n_som,  
const int        n_som_parent,  
const PDM_g_num_t *numabs,
const int       *num_parent,
const PDM_real_t *coords_parent,  
const PDM_g_num_t *numabs_parent
);


/*----------------------------------------------------------------------------
 * Ajout d'un bloc d'elements d'un type donne
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   t_elt           <-- Type d'element
 *
 * return :
 *   id_bloc         --> Identificateur du bloc
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_add, PDM_WRITER_GEOM_BLOC_ADD)
(
int   *id_cs,
int   *id_geom,
PDM_writer_statut_t  *st_free_data,  
int   *t_elt,
int   *id_bloc  
); 

int
PDM_writer_geom_bloc_add 
(
const int            id_cs,
const int            id_geom,   
PDM_writer_statut_t          st_free_data,  
const PDM_writer_elt_geom_t  t_elt
); 


/*----------------------------------------------------------------------------
 * Ajout d'un bloc d'elements d'un type donne dans la partition courante
 *
 *  - PDM_writer_POINT :
 *
 *   1 x            
 *
 *  - PDM_writer_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_writer_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_writer_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_writer_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_writer_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_writer_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_writer_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_bloc         <-- Identificateur du bloc
 *   id_part         <-- Indice de partition
 *   t_elt           <-- Type d'element
 *   n_elt           <-- Nombre d'elements dans le bloc 
 *   connec          <-- Table de connectivite des elements
 *   num_part        <-- Numerotation dans la partition    
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_std_set, PDM_WRITER_GEOM_BLOC_STD_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc,  
int           *id_part, 
int           *n_elt,    
PDM_l_num_t      *connec,   
PDM_g_num_t     *numabs
); 

void
PDM_writer_geom_bloc_std_set 
(
const int            id_cs,
const int            id_geom,  
const int            id_bloc,     
const int            id_part, 
const int            n_elt,    
      PDM_l_num_t      *connec,   
      PDM_g_num_t     *numabs
); 
 
/*----------------------------------------------------------------------------
 * Ajout d'un bloc de polygones dans la partition courante
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part          <-- Indice de partition
 *   n_elt           <-- Nombre d'elements dans le bloc 
 *   connec_idx      <-- Index dans la table de connectivite (dim = n_elt+1)
 *   connec          <-- Table de connectivite des elements (dim = connec_idx[n_elt])
 *   numabs          <-- Numerotation absolue des elements
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_poly2d_set, PDM_WRITER_GEOM_BLOC_POLY2D_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc, 
int           *id_part, 
PDM_l_num_t      *n_elt,    
PDM_l_num_t      *connec_idx,   
PDM_l_num_t      *connec,
PDM_g_num_t     *numabs
); 
 
void
PDM_writer_geom_bloc_poly2d_set 
(
const int            id_cs,
const int            id_geom,  
const int            id_bloc, 
const int            id_part, 
const PDM_l_num_t       n_elt,    
      PDM_l_num_t      *connec_idx,   
      PDM_l_num_t      *connec,
      PDM_g_num_t     *numabs
); 

 
/*----------------------------------------------------------------------------
 * Ajout d'un bloc de polyedres dans la partition courante
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part          <-- Indice de partition
 *   n_elt           <-- Nombre d'elements dans le bloc 
 *   n_face          <-- Nombre de faces de chaque element (dim = n_elt)
 *   facsom_idx      <-- Index dans la table de connectivite des faces (dim = n_face_total+1)
 *   facsom          <-- Table de connectivite des faces (dim = facsom_idx[n_face_total}
 *   cellfac_idx     <-- Index dans la table de connectivite des cellules (dim = n_elt+1)
 *   cellfac         <-- Table de connectivite des elements (dim = cellfac_idx[n_elt])
 *   numabs          <-- Numerotation absolue des elements
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_poly3d_set, PDM_WRITER_GEOM_BLOC_POLY3D_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc, 
int           *id_part, 
PDM_l_num_t      *n_elt,    
PDM_l_num_t      *n_face,   
PDM_l_num_t      *facsom_idx,   
PDM_l_num_t      *facsom,
PDM_l_num_t      *cellfac_idx,   
PDM_l_num_t      *cellfac,
PDM_g_num_t     *numabs
); 

void
PDM_writer_geom_bloc_poly3d_set 
(
const int        id_cs,
const int        id_geom,  
const int        id_bloc, 
const int        id_part, 
const PDM_l_num_t   n_elt,    
const PDM_l_num_t   n_face,   
      PDM_l_num_t  *facsom_idx,   
      PDM_l_num_t  *facsom,
      PDM_l_num_t  *cellfac_idx,   
      PDM_l_num_t  *cellfac,
      PDM_g_num_t *numabs
); 


/*----------------------------------------------------------------------------
 *
 * Ajout de cellules 3D decrites en fonctions des faces. Cette fonction
 * détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Identificateur de partition           
 *   n_elt           <-- Nombre de cellules 3D ajoutées         
 *   n_face          <-- Nombre de faces décrites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   cell_face_idx   <-- Index de connectivite cellules -> faces  
 *   cell_face       <-- Connectivite cellules -> faces
 *   numabs          <-- Numerotatio absolue des cellules 
 *   ind_num         --> Indirection vers la nouvelle numerotation des cellules 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_cell3d_cellface_add, PDM_WRITER_GEOM_CELL3D_CELLFACE_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part, 
int         *n_elt,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
); 

void
PDM_writer_geom_cell3d_cellface_add
(
const int    id_cs,
const int    id_geom,
const int    id_part, 
const int    n_elt,
const int    n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
); 


/*----------------------------------------------------------------------------
 *
 * Ajout de cellules 2D decrites en fonctions des faces. Cette fonction
 * détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   n_elt           <-- Nombre de cellules 3D ajoutées         
 *   n_face          <-- Nombre de faces décrites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   cell_face_idx   <-- Index de connectivite cellules -> faces  
 *   cell_face       <-- Connectivite cellules -> faces
 *   numabs          <-- Numerotatio absolue des cellules 
 *   ind_num         --> Indirection vers la nouvelle numerotation des cellules 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_cell2d_cellface_add, PDM_WRITER_GEOM_CELL2D_CELLFACE_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part, 
int         *n_elt,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
); 

void
PDM_writer_geom_cell2d_cellface_add
(
const int          id_cs,
const int          id_geom,
const int          id_part, 
const int          n_elt,
const int          n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
); 


/*----------------------------------------------------------------------------
 *
 * Ajout de faces decrites en fonctions des sommets. Cette fonction
 * détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   n_elt           <-- Nombre de cellules 3D ajoutées         
 *   n_face          <-- Nombre de faces décrites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   numabs          <-- Numerotation absolue des faces    
 *   ind_num         --> Indirection vers la nouvelle numerotation des faces    
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_faces_facesom_add, PDM_WRITER_GEOM_FACES_FACESOM_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_g_num_t   *numabs
); 

void
PDM_writer_geom_faces_facesom_add
(
const int          id_cs,
const int          id_geom,
const int          id_part, 
const int          n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_g_num_t   *numabs
); 


/*----------------------------------------------------------------------------
 * Ecriture du maillage courant                                  
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_write, PDM_WRITER_GEOM_WRITE)
(
int           *id_cs,
int           *id_geom
); 

void
PDM_writer_geom_write
(
const int            id_cs,
const int            id_geom
); 

/*----------------------------------------------------------------------------
 * Liberation des donnees decrivant le maillage courant
 * les indirections sur les numérotation absolues sont conservées
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_data_free, PDM_WRITER_GEOM_DATA_FREE)
(
int           *id_cs,
int           *id_geom
); 

void
PDM_writer_geom_data_free
(
const int      id_cs,
const int      id_geom
); 


/*----------------------------------------------------------------------------
 * Liberation des donnees decrivant le maillage courant
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_free, PDM_WRITER_GEOM_FREE)
(
int           *id_cs,
int           *id_geom
); 

void
PDM_writer_geom_free
(
const int      id_cs,
const int      id_geom
); 

/*----------------------------------------------------------------------------
 * Creation d'une variable                                                     
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   st_dep_temps    <-- Indique si la variable est dependante du temps
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   dim             <-- Dimension de la variable
 *   loc             <-- Localisation de la variable
 *   nom_var         <-- Nom de la variable
 *
 * return :
 *                   --> Identificateur de l'objet variable     
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_create_cf, PDM_WRITER_VAR_CREATE_CF)
(
int         *id_cs,
int         *st_dep_tps,
int         *dim,
int         *loc,
char        *nom_var,   
int         *l_nom_var,   
int         *id_var 
ARGF_SUPP_CHAINE
);

int
PDM_writer_var_create
(
const int          id_cs,
const PDM_writer_statut_t  st_dep_tps,
const PDM_writer_var_dim_t dim,
const PDM_writer_var_loc_t loc,
const char        *nom_var
);

/*----------------------------------------------------------------------------
 * Mapping des noms de variable                                                     
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   public_name     <-- Nom Public de la variable
 *   pivate_name     <-- Nom privé de la variable
 *
 * return :
 *                   --> Identificateur de l'objet variable     
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_name_map_add_cf, PDM_WRITER_NAME_MAP_ADD_CF)
(
int         *id_cs,
char        *public_name,
int         *l_public_name,
char        *private_name,
int         *l_private_name
ARGF_SUPP_CHAINE
);

void
PDM_writer_name_map_add
(
const int   id_cs,
const char *public_name,
const char *private_name
);

/*----------------------------------------------------------------------------
 * Mise a jour des valeurs de la variable
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Identificateur de la partition dans l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise à jour
 *   val             <-- Valeurs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_write, PDM_WRITER_VAR_WRITE)
(
int         *id_cs,
int         *id_var
);

void
PDM_writer_var_write
(
const int        id_cs,
const int        id_var
);


/*----------------------------------------------------------------------------
 * Mise a jour des valeurs de la variable. Attention, les valeurs définies aux
 * elements doivent être définies suivant l'ordre de définition des blocs !
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Identificateur de la partition dans l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise à jour
 *   val             <-- Valeurs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_set, PDM_WRITER_VAR_SET)
(
int         *id_cs,
int         *id_var,
int         *id_geom,
int         *id_part,
PDM_real_t   *val
);

void
PDM_writer_var_set
(
const int        id_cs,
const int        id_var,
const int        id_geom,
const int        id_part,
const PDM_real_t *val
);

/*----------------------------------------------------------------------------
 * Liberation du tableau de donnees des variables
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise à jour
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_data_free, PDM_WRITER_VAR_DATA_FREE)
(
int         *id_cs,
int         *id_var
);

void
PDM_writer_var_data_free
(
const int    id_cs,
const int    id_var
);


/*----------------------------------------------------------------------------
 * Liberation du tableau de donnees des variables
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise à jour
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_free, PDM_WRITER_VAR_FREE)
(
int         *id_cs,
int         *id_var
);

void
PDM_writer_var_free
(
const int    id_cs,
const int    id_var
);


/**
 * \brief Add a writer format
 *
 * Define a new format writer
 *
 * \param [in] name            Name                                                    
 * \param [in] create_fct      Customize \ref PDM_writer_create function for the new format  (or NULL)
 * \param [in] free_fct        Customize \ref PDM_writer_free function for the new format (or NULL)
 * \param [in] beg_step_fct    Customize \ref PDM_writer_step_beg function for the new format (or NULL)
 * \param [in] end_step_fct    Customize \ref PDM_writer_step_end function for the new format (or NULL)
 * \param [in] geom_create_fct Customize \ref PDM_writer_geom_create function for the new format (or NULL)
 * \param [in] geom_write_fct  Customize \ref PDM_writer_geom_write function for the new format
 * \param [in] geom_free_fct   Customize \ref PDM_writer_geom_free function for the new format (or NULL)
 * \param [in] var_create_fct  Customize \ref PDM_writer_var_create function for the new format (or NULL)
 * \param [in] var_write_fct   Customize \ref PDM_writer_var_write function for the new format
 * \param [in] var_free_fct    Customize \ref PDM_writer_var_free function for the new format (or NULL)
 *
 */

void
PDM_writer_fmt_add
(
 const char                  *name,           /*!< Name                                                     */
 const PDM_writer_fct_t      create_fct,      /*!< Customize \ref PDM_writer_create function for the format */
 const PDM_writer_fct_t      free_fct,        /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t      beg_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t      end_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_create_fct, /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_write_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_free_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_create_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_write_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_free_fct     /*!< Customize \ref PDM_writer_free function for the format   */
);


/**
 * \brief Free format
 *
 */

void
PDM_writer_fmt_free
(
 void
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_H__ */
