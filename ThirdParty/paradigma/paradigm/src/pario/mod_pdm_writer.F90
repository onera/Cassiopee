module mod_pdm_writer
  
  use mod_pdm
  
  implicit none

  !
  ! Statut

  integer, parameter :: PDM_WRITER_OFF  = 0
  integer, parameter :: PDM_WRITER_ON   = 1

  !
  ! Type de topologie

  integer, parameter :: PDM_WRITER_TOPO_CONSTANTE  = 0     
  integer, parameter :: PDM_WRITER_TOPO_DEFORMABLE = 1     
  integer, parameter :: PDM_WRITER_TOPO_VARIABLE   = 2   

  !
  ! Type d'elements géometriques

  integer, parameter :: PDM_WRITER_POINT    = 0    
  integer, parameter :: PDM_WRITER_BAR2     = 1     
  integer, parameter :: PDM_WRITER_TRIA3    = 2     
  integer, parameter :: PDM_WRITER_QUAD4    = 3     
  integer, parameter :: PDM_WRITER_POLY_2D  = 4     
  integer, parameter :: PDM_WRITER_TETRA4   = 5     
  integer, parameter :: PDM_WRITER_PYRAMID5 = 6     
  integer, parameter :: PDM_WRITER_PRISM6   = 7     
  integer, parameter :: PDM_WRITER_HEXA8    = 8     
  integer, parameter :: PDM_WRITER_POLY_3D  = 9     

  !
  ! Format de sortie

  integer, parameter ::  PDM_WRITER_FMT_ENSIGHT = 0     

  !
  ! Format du fichier

  integer, parameter ::  PDM_WRITER_FMT_BIN   = 0
  integer, parameter ::  PDM_WRITER_FMT_ASCII = 1

  !
  ! Dimension géométrique de la sortie                      

  integer, parameter ::  PDM_WRITER_DIM_2 = 0     
  integer, parameter ::  PDM_WRITER_DIM_3 = 1     

  !
  ! Dim des variables

  integer, parameter ::  PDM_WRITER_VAR_CSTE         = 0
  integer, parameter ::  PDM_WRITER_VAR_SCALAIRE     = 1
  integer, parameter ::  PDM_WRITER_VAR_VECTEUR      = 3
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_SYM  = 6
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_ASYM = 9

  !
  ! Localisation des variables

  integer, parameter ::  PDM_WRITER_VAR_SOMMETS      = 0
  integer, parameter ::  PDM_WRITER_VAR_ELEMENTS     = 1
  integer, parameter ::  PDM_WRITER_VAR_PARTICULES   = 2

  !
  ! Interfaces publiques
  ! --------------------

  interface pdm_writer_create       ; module procedure &
    pdm_writer_create_
  end interface

  interface pdm_writer_geom_create  ; module procedure &
    pdm_writer_geom_create_
  end interface

  interface pdm_writer_var_create   ; module procedure &
    pdm_writer_var_create_
  end interface
  
  interface pdm_writer_name_map_add   ; module procedure &
    pdm_writer_name_map_add_
  end interface
  
  !
  ! Fonctions Privees
  ! -----------------

  private :: pdm_writer_create_
  private :: pdm_writer_geom_create_
  private :: pdm_writer_var_create_
  private :: pdm_writer_name_map_add_

contains 

  !----------------------------------------------------------------------------
  ! Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet 
  !
  ! parameters :
  !   fmt             <-- Format de sortie
  !   topologie       <-- Indique le maillage est mobile ou non
  !   st_decoup_poly  <-- Active le decoupage des polyedres
  !   st_reprise      <-- Complete les sorties des calcul precedents en reprise
  !   dim_geo         <-- Dimension geometrique             
  !   rep_sortie      <-- Repertoire de sortie                  
  !   nom_sortie      <-- Nom de la sortie                       
  !   msg_com         <-- Communicateur MSG                      
  !
  ! return :
  !                   --> Identificateur de l'objet cree
  !
  !----------------------------------------------------------------------------

  subroutine pdm_writer_create_ ( &
    fmt,                &
    fmt_fic,            &
    topologie,          &
    st_reprise,         &
    rep_sortie,         &
    nom_sortie,         &
    msg_comm,           &
    acces,              &
    prop_noeuds_actifs, &
    options,            &
    id_cs)

    implicit none
    ! Arguments

    character (len = *),          intent(in)    :: fmt
    integer (kind = pdm_l_num_s), intent(in)    :: fmt_fic
    integer (kind = pdm_l_num_s), intent(in)    :: topologie
    integer (kind = pdm_l_num_s), intent(in)    :: st_reprise
    character (len = *),          intent(in)    :: rep_sortie
    character (len = *),          intent(in)    :: nom_sortie
    integer (kind = pdm_l_num_s), intent(in)    :: msg_comm   
    integer (kind = pdm_l_num_s), intent(in)    :: acces      
    real (kind = 8),              intent(in)    :: prop_noeuds_actifs
    character (len = *),          intent(in)    :: options
    integer (kind = pdm_l_num_s), intent(inout) :: id_cs     

    ! Variables locales
    integer :: l_fmt
    integer :: l_rep_sortie
    integer :: l_nom_sortie
    integer :: l_options

    !
    ! Calcul de la longueur des chaines pour la conversion en string C 
    l_fmt = len(fmt)
    l_rep_sortie = len(rep_sortie)
    l_nom_sortie = len(nom_sortie)
    l_options= len(options)

    call pdm_writer_create_cf( &
         fmt,                &
         l_fmt,              &
         fmt_fic,            &
         topologie,          &
         st_reprise,         &
         rep_sortie,         &
         nom_sortie,         &
         l_rep_sortie,       &
         l_nom_sortie,       &
         msg_comm,           &
         acces,              &
         prop_noeuds_actifs, &
         options,            &
         l_options,          &
         id_cs)

  end subroutine pdm_writer_create_

  !----------------------------------------------------------------------------
  ! Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
  ! (codee en C)
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_lib (id_cs)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(inout) :: id_cs     


  !----------------------------------------------------------------------------
  ! Debut d'increment
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   delta_t         <-- Delta de temps par rapport au dernier increment
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_incr_deb (id_cs, delta_t)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     real (kind = 8)             , intent(in)    :: delta_t   

  !----------------------------------------------------------------------------
  ! Fin d'increment
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_incr_fin (id_cs)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     

  !----------------------------------------------------------------------------
  ! Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   nom_geom        <-- Nom de l'objet geometrique
  !
  ! return :
  !                   --> Identificateur de l'objet geom dans cs 
  !
  !----------------------------------------------------------------------------

  subroutine pdm_writer_geom_create_ (id_cs,           &
                            nom_geom,        &
                            st_decoup_poly2d, &
                            st_decoup_poly3d, &
                            n_part, &
                            id_geom)

    implicit none

    !
    ! Arguments

    integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
    character (len = *),       intent(in)    :: nom_geom   
    integer (kind = pdm_l_num_s), intent(inout) :: st_decoup_poly2d
    integer (kind = pdm_l_num_s), intent(inout) :: st_decoup_poly3d
    integer (kind = pdm_l_num_s), intent(inout) :: n_part   
    integer (kind = pdm_l_num_s), intent(inout) :: id_geom   

    !
    ! Variables locales

    integer :: l_nom_geom

    !
    ! Calcul de la longueur des chaines pour la conversion en string C 

    l_nom_geom = len(nom_geom)

    call pdm_writer_geom_create_cf (id_cs,            &
                          nom_geom,         &
                          st_decoup_poly2d, &
                          st_decoup_poly3d, &
                          l_nom_geom,       &
                          n_part,           &
                          id_geom)

  end subroutine pdm_writer_geom_create_

  !----------------------------------------------------------------------------
  ! Definition du nombre de partitions sur le processus courant 
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   n_part          <-- Nombre de partitions sur le processus courant                
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_n_part_set (id_cs, id_geom, n_part)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_part 

  !----------------------------------------------------------------------------
  ! Definition des coordonnees de la partition courante          
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_part          <-- Indice de partition
  !   n_som           <-- Nombre de sommets de la partition
  !   coords          <-- Coordonnes des sommets            
  !   numabs          <-- Numerotation absolue des sommets     
  !
  !---------------------------------------------------------------------------*/
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_coord_set (id_cs, id_geom, n_part)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_part 
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_som
  !     real (kind = 8),              intent(in)    :: coords(*)
  !     integer (kind = pdm_g_num_s), intent(in)    :: numabs(*)

  !----------------------------------------------------------------------------
  ! Ajout d'un bloc d'elements d'un type donne dans la partition courante
  !
  !  - PDM_WRITER_POINT :
  !
  !   1 x            
  !
  !  - PDM_WRITER_BAR2 :
  !
  !   1 x-------x 2
  !
  !  - PDM_WRITER_TRIA3 :   
  !
  !   1 x-------x 3
  !      \     /
  !       \   /
  !        \ /
  !         x 2
  !
  !  - PDM_WRITER_QUAD4 :          
  !
  !      4 x-------x 3
  !       /       /
  !      /       /
  !   1 x-------x2
  !
  !   - PDM_WRITER_TETRA4 :    
  !
  !         x 4
  !        /|\
  !       / | \
  !      /  |  \
  !   1 x- -|- -x 3
  !      \  |  /
  !       \ | /
  !        \|/
  !         x 2
  !
  !   - PDM_WRITER_PYRAMID5 :
  !
  !          5 x
  !           /|\
  !          //| \
  !         // |  \
  !      4 x/--|---x 3
  !       //   |  /
  !      //    | /
  !   1 x-------x 2
  !
  !  - PDM_WRITER_PRSIM6 :
  !
  !   4 x-------x 6
  !     |\     /|
  !     | \   / |
  !   1 x- \-/ -x 3
  !      \ 5x  /
  !       \ | /
  !        \|/
  !         x 2
  !
  !  - PDM_WRITER_HEXA8 :   
  !
  !      8 x-------x 7
  !       /|      /|
  !      / |     / |
  !   5 x-------x6 |
  !     | 4x----|--x 3
  !     | /     | /
  !     |/      |/
  !   1 x-------x 2
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_part          <-- Indice de partition
  !   t_elt           <-- Type d'element
  !   n_elt           <-- Nombre d'elements dans le bloc 
  !   connec          <-- Table de connectivite des elements
  !   numabs          <-- Numerotation absolue des elements
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_bloc_add (id_cs, id_geom, id_part, t_elt, n_elt, connec
  !                              numabs)  
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_part 
  !     integer (kind = pdm_l_num_s), intent(in)    :: t_elt 
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_elt  
  !     integer (kind = pdm_l_num_s), intent(in)    :: connec(*)
  !     integer (kind = pdm_g_num_s), intent(in)   :: numabs(*)
  
  !----------------------------------------------------------------------------
  ! Ajout d'un bloc de polygones dans la partition courante
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_part          <-- Indice de partition
  !   n_elt           <-- Nombre d'elements dans le bloc 
  !   connec_idx      <-- Index dans la table de connectivite (dim = n_elt+1)
  !   connec          <-- Table de connectivite des elements (dim = connec_idx[n_elt])
  !   numabs          <-- Numerotation absolue des elements
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_bloc_poly2d_add (id_cs, id_geom, id_part, n_elt, &
  !                                     connec_idx, connec, numabs)  
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_part 
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_elt  
  !     integer (kind = pdm_l_num_s), intent(in)    :: connec_idx(*)
  !     integer (kind = pdm_l_num_s), intent(in)    :: connec(*)
  !     integer (kind = pdm_g_num_s), intent(in)   :: numabs(*)
  
  !----------------------------------------------------------------------------
  ! Ajout d'un bloc de polyedres dans la partition courante
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_part          <-- Indice de partition
  !   n_elt           <-- Nombre d'elements dans le bloc 
  !   n_face          <-- Nombre de faces de chaque element (dim = n_elt)
  !   facsom_idx      <-- Index dans la table de connectivite des faces (dim = n_face_total+1)
  !   facsom          <-- Table de connectivite des faces (dim = facsom_idx[n_face_total}
  !   cellfac_idx     <-- Index dans la table de connectivite des cellules (dim = n_elt+1)
  !   cellfac         <-- Table de connectivite des elements (dim = cellfac_idx[n_elt])
  !   numabs          <-- Numerotation absolue des elements
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_bloc_poly3d_add (id_cs, id_geom, id_part, n_elt, &
  !                                     n_face, facsom_idx, facsom, &
  !                                     cellfac_idx, cellfac, numabs)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_part 
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_elt  
  !     integer (kind = pdm_l_num_s), intent(in)    :: n_face(*)
  !     integer (kind = pdm_l_num_s), intent(in)    :: facsom_idx(*)
  !     integer (kind = pdm_l_num_s), intent(in)    :: facsom(*)
  !     integer (kind = pdm_l_num_s), intent(in)    :: cellfac_idx(*)
  !     integer (kind = pdm_l_num_s), intent(in)    :: cellfac(*)
  !     integer (kind = pdm_g_num_s), intent(in)   :: numabs(*)
  
  !----------------------------------------------------------------------------
  ! Ecriture du maillage courant                                  
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_write (id_cs, id_geom)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  
  !----------------------------------------------------------------------------
  ! Liberation des donnees decrivant le maillage courant                        
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_geom_donnees_lib (id_cs, id_geom)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  
  !----------------------------------------------------------------------------
  ! Creation d'une variable                                                     
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   st_dep_temps    <-- Indique si la variable est dependante du temps
  !   t_var           <-- Type de variable
  !   nom_var         <-- Nom de la variable
  !
  ! return :
  !                   --> Identificateur de l'objet variable     
  !
  !----------------------------------------------------------------------------

  subroutine pdm_writer_var_create_ (id_cs,           &
                           st_dep_tps,      &
                           dim,             &
                           loc,             &
                           nom_var,         &   
                           id_var) 

    implicit none

    !
    ! Arguments

    integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
    integer (kind = pdm_l_num_s), intent(in)    :: st_dep_tps
    integer (kind = pdm_l_num_s), intent(in)    :: dim     
    integer (kind = pdm_l_num_s), intent(in)    :: loc     
    character (len = *),       intent(in)    :: nom_var
    integer (kind = pdm_l_num_s), intent(inout) :: id_var    

    !
    ! Variables locales

    integer :: l_nom_var

    !
    ! Calcul de la longueur des chaines pour la conversion en string C 

    l_nom_var = len(nom_var)

    call pdm_writer_var_create_cf (id_cs,           &
                         st_dep_tps,      &
                         dim,             &
                         loc,             &
                         nom_var,         &   
                         l_nom_var,       &   
                         id_var) 

  end subroutine pdm_writer_var_create_
  

  !----------------------------------------------------------------------------
  ! Mapping des noms de variable                                                     
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   public_name     <-- Nom Public de la variable
  !   pivate_name     <-- Nom privé de la variable
  !
  ! return :
  !                   --> Identificateur de l'objet variable     
  !
  !----------------------------------------------------------------------------

  subroutine pdm_writer_name_map_add_ ( &
    id_cs,        &
    public_name,  &
    private_name )
    implicit none
    ! Arguments
    integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
    character (len = *),          intent(in)    :: public_name
    character (len = *),          intent(in)    :: private_name

    ! Variables locales
    integer :: l_public_name
    integer :: l_private_name

    ! Calcul de la longueur des chaines pour la conversion en string C 
    l_public_name = len(public_name)
    l_private_name = len(private_name)

    call pdm_writer_name_map_add_cf ( &
         id_cs,           &
         public_name,     &
         l_public_name,   &
         private_name,    &
         l_private_name)

  end subroutine pdm_writer_name_map_add_


  !----------------------------------------------------------------------------
  ! Mise a jour des valeurs de la variable
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_part         <-- Identificateur de la partition dans l'objet geometrique
  !   id_var          <-- Identificateur de la variable mise à jour
  !   val             <-- Valeurs
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_var_ecr (id_cs, id_geom, id_part, id_var, val)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_part 
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_var
  !     real*8,                    intent(in)    :: val(*)    
  
  !----------------------------------------------------------------------------
  ! Liberation du tableau de donnees des variables
  !
  ! parameters :
  !   id_cs           <-- Identificateur de l'objet cs
  !   id_geom         <-- Identificateur de l'objet geometrique
  !   id_var          <-- Identificateur de la variable mise à jour
  !
  !----------------------------------------------------------------------------
  !
  ! Fonction definie en C dans cedre_sorties.c
  !
  ! subroutine pdm_writer_ var_donnees_lib (id_cs, id_geom, id_var)
  !
  !     ! 
  !     ! Arguments
  !
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_cs     
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_geom   
  !     integer (kind = pdm_l_num_s), intent(in)    :: id_var
  
end module mod_pdm_writer
