
module mod_pdm_io
  
  use mod_pdm
  
  implicit none

  !
  ! Parametres
  ! ----------

  !
  ! Types de sufixe                          

  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_auto        = 0 ! Suffixe detemine automatiquement                      
  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_man         = 1 ! Suffixe fourni par l'utilisateur                       

  !
  ! Endianness                                 

  integer (kind = pdm_l_num_s), parameter :: pdm_io_bigendian        = 0 ! Contenu en bigendian
  integer (kind = pdm_l_num_s), parameter :: pdm_io_littleendian     = 1 ! Contenu en little endian
  integer (kind = pdm_l_num_s), parameter :: pdm_io_native           = 3 ! Contenu natif machine
  
  !
  ! Types de données

  integer (kind = pdm_l_num_s), parameter :: pdm_io_t_int    = 0 ! Type de donnee int
  integer (kind = pdm_l_num_s), parameter :: pdm_io_t_long   = 1 ! Type de donnee long
  integer (kind = pdm_l_num_s), parameter :: pdm_io_t_double = 2 ! Type de donnee double
  integer (kind = pdm_l_num_s), parameter :: pdm_io_t_float  = 3 ! Type de donnee float
  integer (kind = pdm_l_num_s), parameter :: pdm_io_t_char   = 4 ! Type de donnee char

  !
  ! Types d'entrees/sorties paralleles

  integer (kind = pdm_l_num_s), parameter :: pdm_io_acces_mpiio_eo   = 0 ! Acces parallele avec MPIIO (explicit offset)
  integer (kind = pdm_l_num_s), parameter :: pdm_io_acces_mpiio_ip   = 1 ! Acces parallele avec MPIIO (individual pointer)
  integer (kind = pdm_l_num_s), parameter :: pdm_io_acces_mpi_simple = 2 ! Acces parallele sans MPI
  integer (kind = pdm_l_num_s), parameter :: pdm_io_acces_seq        = 3 ! Acces 1 fichier par processus

  !
  ! Mode d'acces lecture, ecriture, lecture/ecriture

  integer (kind = pdm_l_num_s), parameter :: pdm_io_mode_lecture  = 0    ! Acces en lecture
  integer (kind = pdm_l_num_s), parameter :: pdm_io_mode_ecriture = 1    ! Acces en ecriture
  integer (kind = pdm_l_num_s), parameter :: pdm_io_mode_ajout    = 2    ! Acces en lecture/ecriture

  !
  ! Indique si le nombre de sous-variables d'une variable est constant ou variable
  ! selon le numero de la donnee

  integer (kind = pdm_l_num_s), parameter :: pdm_io_n_composante_constant = 0  
  integer (kind = pdm_l_num_s), parameter :: pdm_io_n_composante_variable = 1  

  !
  ! Indique si les donnees lues ou écrites sont rangÃ©es par blocs contigus en mÃ©moire
  ! ou respectent une indirection (donnees entrelacees)

  integer(kind = pdm_l_num_s), parameter :: pdm_io_donnees_bloc = 0  
  integer(kind = pdm_l_num_s), parameter :: pdm_io_donnees_entrelacee = 1  

  !
  ! Indique si le fichier contient une entete IOCEDRE

  integer(kind = pdm_l_num_s), parameter :: pdm_io_entete_on  = 0
  integer(kind = pdm_l_num_s), parameter :: pdm_io_entete_off = 1

  !
  ! Indique le format du fichier

  integer(kind = pdm_l_num_s), parameter :: pdm_io_fmt_txt = 0
  integer(kind = pdm_l_num_s), parameter :: pdm_io_fmt_bin = 1

  !
  ! Active ou non le backup d'un fichier preexistant

  integer(kind = pdm_l_num_s), parameter :: pdm_io_backup_on  = 0
  integer(kind = pdm_l_num_s), parameter :: pdm_io_backup_off = 1

  !
  ! Interfaces publiques
  ! --------------------

  interface pdm_io_open ; module procedure &
    pdm_io_open_
  end interface

  interface pdm_io_fmt_donnee_set ; module procedure &
    pdm_io_fmt_donnee_set_
  end interface

  !
  ! Fonctions Privees
  ! -----------------

  private :: pdm_io_open_
  private :: pdm_io_fmt_donnee_set_

contains 

!----------------------------------------------------------------------------
! Ouverture d'un fichier pour acces parallele
!
! parameters :
!    nom             <-- Nom du fichier
!    fmt             <-- Fichier text ou binaire
!    s_entete        <-- Fichier avec ou sans entete 
!    suff_t          <-- Type de suffixe (manuel ou automatique)
!    suff_u          <-- Suffixe (si suffixe manuel)
!    s_backup        <-- Active le backup d'un fichier preexistant en mode ecriture
!    accesio         <-- Type (parallele avec mpiio, parallele sans mpiio,
!                              sequentiel)
!    mode            <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
!    msg_comm        <-- Communicateur lie au fichier
!    unite           --> Unite du fichier
!    ierr            --> Indique si le fichier est de type pdm_io ou non     
!                        Utiliser uniquement pour une ouverture en lecture
!
!----------------------------------------------------------------------------

  subroutine pdm_io_open_(nom, &
                          fmt, &
                          suff_t, &
                          suff_u, &
                          s_backup, &
                          type_io, &
                          mode, &
                          endian, &
                          msg_comm, &
                          prop_noeuds_actifs,  &
                          unite, &
                          ierr)
    implicit none

    !
    ! Arguments

    character (len = *),          intent(in)  :: nom
    integer (kind = pdm_l_num_s), intent(in)  :: fmt       
    integer (kind = pdm_l_num_s), intent(in)  :: suff_t  
    character (len = *),          intent(in)  :: suff_u
    integer (kind = pdm_l_num_s), intent(in)  :: s_backup
    integer (kind = pdm_l_num_s), intent(in)  :: type_io
    integer (kind = pdm_l_num_s), intent(in)  :: mode
    integer (kind = pdm_l_num_s), intent(in)  :: endian       
    integer (kind = pdm_l_num_s), intent(in)  :: msg_comm
    real                        , intent(in)  :: prop_noeuds_actifs
    integer (kind = pdm_l_num_s), intent(out) :: unite
    integer (kind = pdm_l_num_s), intent(out) :: ierr 
    
    !
    ! Variables locales

    integer            :: l_nom ! Longueur de la chaine nom
    integer            :: l_suff_u ! Longueur de la chaine nom

    !
    ! Calcul de la longueur des chaines pour la conversion en string C 

    l_nom = len(nom)
    l_suff_u = len(suff_u)

    !
    ! Appel de la fonction C

    call pdm_io_open_cf(nom, &
                        l_nom, &
                        fmt, &
                        suff_t, &
                        suff_u, &
                        l_suff_u, &
                        s_backup, &
                        type_io, &
                        mode, &
                        endian, &
                        msg_comm, &
                        prop_noeuds_actifs, &
                        unite, &
                        ierr)

  end subroutine pdm_io_open_

!----------------------------------------------------------------------------
! Ouverture d'un fichier pour acces parallele
!
! parameters :
!   nom             <-- Nom du fichier
!   suff_t          <-- Type de suffixe :             
!                           - pdm_io_suff_auto
!                           - pdm_io_suff_man
!   type_io         <-- Type d'acces a choisir entre :
!                           - pdm_io_acces_mpiio_eo
!                           - pdm_io_acces_mpiio_ip
!                           - pdm_io_acces_mpi_simple
!                           - pdm_io_acces_seq
!   mode            <-- Mode d'acces a choisir entre : 
!                           - pdm_io_mode_lecture
!                           - pdm_io_mode_ecriture
!                           - pdm_io_mode_ajout
!   msg_comm        <-- Communicateur lie au fichier
!   unite           --> Unite du fichier
!   ierr            --> Indicateur d'erreur  
!
!----------------------------------------------------------------------------

  subroutine pdm_io_fmt_donnee_set_(unite, n_char_fmt, data_type, fmt)
    implicit none

    !
    ! Arguments

    integer (kind = pdm_l_num_s), intent(in)  :: unite  
    integer (kind = pdm_l_num_s), intent(in)  :: n_char_fmt  
    integer (kind = pdm_l_num_s), intent(in)  :: data_type
    character (len = *),          intent(in)  :: fmt
    
    !
    ! Variables locales

    integer (kind = pdm_l_num_s) :: l_fmt ! Longueur de la chaine nom

    !
    ! Calcul de la longueur des chaines pour la conversion en string C 

    l_fmt = len(fmt)

    !
    ! Appel de la fonction C

    call pdm_io_fmt_donnee_set_cf(unite, n_char_fmt, data_type, fmt, l_fmt)

  end subroutine pdm_io_fmt_donnee_set_

!----------------------------------------------------------------------------
! Lecture globale : Le processus maitre accede seul au fichier et redistribue
! l'information a l'ensemble des processus du communicateur
!
! arguments :
!   unite           <-- Unite du fichier
!   taille_donnee   <-- Taille unitaire de la donnnee en octet
!   n_donnees       <-- Nombre de donnees a lire
!   donnees         --> Donnees lues
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_lecture_globale(unite, &
!                                       taille_donnee, & 
!                                       n_donnees, & 
!                                       donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),          intent(in)  :: unite
!     integer (kind = pdm_l_num_s),          intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),          intent(in)  :: n_donnees
!     "Pas de type : void* en C", dimension(*), intent(out) :: donnees


!----------------------------------------------------------------------------
! Ecriture globale : Le processus maitre ecrit seul sur le fichier
!
! arguments :
!   unite             <-- Unite du fichier
!   taille_donnee     <-- Taille unitaire de la donnnee en octet
!   n_donnees         <-- Nombre de donnees a ecrire
!   donnees           <-- Donnees a ecrire
!  
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_ecriture_global(unite, &
!                                       taille_donnee, & 
!                                       n_donnees, & 
!                                       donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),          intent(in)  :: unite
!     integer (kind = pdm_l_num_s),          intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),          intent(in)  :: n_donnees
!     "Pas de type : void* en C", dimension(*), intent(int) :: donnees


!----------------------------------------------------------------------------
! Initialise une phase d'Ã©criture parallÃ¨le de tableaux de donnÃ©es associÃ©es
! aux numÃ©ros de variables CEDRE
! Chaque tableau a ses propres caractÃ©ristiques :
!         - taille de donnÃ©es
!         - nombre de donnÃ©e
!         - indirection (numÃ©rotation absolue)
! 
!  arguments :
!    unite             <-- Unite du fichier
!    rangement         <-- Type de rangement A choisir entre :
!                              - pdm_io_donnees_bloc
!                              - pdm_io_donnees_entrelacee
!    num_var_pdm_max <-- NumÃ©ro max de variable CEDRE
!    n_partition_local <-- Nombre de partitions locales
! 
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!   subroutine pdm_io_tab_ecr_debut(unite, &
!                                     t_rangement, &
!                                     num_var_pdm_max, &
!                                     n_partition_local)
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),      intent(in)  :: unite
!     integer,                              intent(in)  :: t_rangement
!     integer (kind = pdm_l_num_s),      intent(in)  :: num_var_pdm_max
!     integer (kind = pdm_l_num_s),      intent(in)  :: n_partition_local


!----------------------------------------------------------------------------
! Ajoute une partie des donnees dans un tableau associÃ©s Ã  une variable
! CEDRE
!
! arguments :
!   num_var_cedre         <-- NumÃ©ro de variable CEDRE 
!   num_indirection_cedre <-- NumÃ©ro d'indirection CEDRE         
!   i_part                <-- indice de partition 
!   t_n_composantes       <-- Type de tailles composantes A choisir entre :
!                               - pdm_io_n_composante_constant
!                               - pdm_io_n_composante_variable
!   n_composantes         <-- Nombre de composantes pour chaque donnee
!   taille_donnee         <-- Taille unitaire de la donnnee
!   n_donnees             <-- Nombre de donnees a lire
!   indirection           <-- Indirection de redistribition des donnees si
!                             les donnÃ©es sonr entrelacÃ©es et adresse de dÃ©but
!                             de bloc pour des donnÃ©es dÃ©finies par bloc
!   donnees               <-- Donnees a Ã©crire
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!  subroutine pdm_io_tab_ecr_ajout_donnees(num_var_cedre, &         
!                                            num_indirection_cedre, &         
!                                            i_part, &         
!                                            t_n_composantes, &         
!                                            n_composantes, &
!                                            taille_donnee, &
!                                            n_donnees, &
!                                            indirection, &
!                                            donnees)
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s),           intent(in)  :: num_var_cedre
!     integer (kind = pdm_l_num_s),           intent(in)  :: num_indirection_cedre
!     integer,                                   intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     integer (kind = pdm_g_num_s), dimension(*), intent(in)  :: indirection
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Finalise une phase d'Ã©criture parallÃ¨le de tableaux de donnÃ©es associÃ©es
! aux numÃ©ros de variables CEDRE. Cette fonction dÃ©clenche rÃ©ellement
! les Ã©critures
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!
!  subroutine pdm_io_tab_ecr_fin


!----------------------------------------------------------------------------
! Initialise une phase d'Ã©criture parallÃ¨le de tableaux de donnÃ©es associÃ©es
! aux numÃ©ros de variables CEDRE
! Chaque tableau a ses propres caractÃ©ristiques :
!         - taille de donnÃ©es
!         - nombre de donnÃ©e
!         - indirection (numÃ©rotation absolue)
! 
!  arguments :
!    unite             <-- Unite du fichier
!    rangement         <-- Type de rangement A choisir entre :
!                              - pdm_io_donnees_bloc
!                              - pdm_io_donnees_entrelacee
!    num_var_pdm_max <-- NumÃ©ro max de variable CEDRE
!    n_partition_local <-- Nombre de partitions locales
! 
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!   subroutine pdm_io_tab_ecr_debut(unite, &
!                                     t_rangement, &
!                                     num_var_pdm_max, &
!                                     n_partition_local)
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),      intent(in)  :: unite
!     integer,                              intent(in)  :: t_rangement
!     integer (kind = pdm_l_num_s),      intent(in)  :: num_var_pdm_max
!     integer (kind = pdm_l_num_s),      intent(in)  :: n_partition_local


!----------------------------------------------------------------------------
! Ajoute une partie des donnees dans un tableau associÃ©s Ã  une variable
! CEDRE
!
! arguments :
!   num_var_cedre         <-- NumÃ©ro de variable CEDRE 
!   num_indirection_cedre <-- NumÃ©ro d'indirection CEDRE         
!   i_part                <-- indice de partition 
!   t_n_composantes       <-- Type de tailles composantes A choisir entre :
!                               - pdm_io_n_composante_constant
!                               - pdm_io_n_composante_variable
!   n_composantes         <-- Nombre de composantes pour chaque donnee
!   taille_donnee         <-- Taille unitaire de la donnnee
!   n_donnees             <-- Nombre de donnees a lire
!   indirection           <-- Indirection de redistribition des donnees si
!                             les donnÃ©es sonr entrelacÃ©es et adresse de dÃ©but
!                             de bloc pour des donnÃ©es dÃ©finies par bloc
!   donnees               <-- Tableau de stockage des donnÃ©e lues
!                             (les donnÃ©es sont accessibles aprÃ¨s appel
!                              de la fonction  pdm_io_tab_lec_fin)
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!  subroutine pdm_io_tab_lec_ajout_donnees(num_var_cedre, &         
!                                            num_indirection_cedre, &         
!                                            i_part, &         
!                                            t_n_composantes, &         
!                                            n_composantes, &
!                                            taille_donnee, &
!                                            n_donnees, &
!                                            indirection, &
!                                            donnees)
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s),           intent(in)  :: num_var_cedre
!     integer (kind = pdm_l_num_s),           intent(in)  :: num_indirection_cedre
!     integer,                                   intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     integer (kind = pdm_g_num_s), dimension(*), intent(in)  :: indirection
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Finalise une phase delecture parallÃ¨le de tableaux de donnÃ©es associÃ©es
! aux numÃ©ros de variables CEDRE. Cette fonction dÃ©clenche rÃ©ellement
! les Ã©critures
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io_tab.c
!
!
!  subroutine pdm_io_tab_lec_fin


!----------------------------------------------------------------------------
! Lecture parallele de donnees entrelacees
!
! arguments :
!   unite           <-- Unite du fichier
!   t_n_composantes <-- Indique si le nombre de composantes est constant
!                       A choisir entre :
!                             - pdm_io_n_composante_constant
!                             - pdm_io_n_composante_variable
!   n_composantes   <-- Nombre de composantes (tableau ou scalaire suivant 
!                                              la valeur de t_n_composantes)
!   taille_donnee   <-- Taille unitaire de la donnnee en octet
!   n_donnees       <-- Nombre de donnees a lire
!   indirection     <-- Indirection (Numerotation absolue)
!   donnees         --> Donnees lues
!  
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_lec_par_entrelacee(unite, &
!                                          t_n_composantes, &
!                                          n_composantes, &
!                                          taille_donnee, & 
!                                          n_donnees, & 
!                                          indirection, &
!                                          donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite
!     integer,                                   intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     integer (kind = pdm_g_num_s), dimension(*), intent(in)  :: indirection
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Lecture parallele de blocs continus de donnees
!
! arguments :
!   unite           <-- Unite du fichier
!   t_n_composantes <-- Indique si le nombre de composantes est constant
!                       A choisir entre :
!                             - pdm_io_n_composante_constant
!                             - pdm_io_n_composante_variable
!   n_composantes   <-- Nombre de composantes (tableau ou scalaire suivant 
!                                              la valeur de t_n_composantes)
!   taille_donnee   <-- Taille unitaire de la donnnee en octet
!   n_donnees       <-- Nombre de donnees a ecrire
!   debut_bloc      <-- Index du debut du bloc
!   donnees         --> Donnees lues
!  
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_lec_par_bloc(unite, &
!                                    t_n_composantes, &
!                                    n_composantes, &
!                                    taille_donnee, & 
!                                    n_donnees, & 
!                                    debut_bloc, &
!                                    donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite
!     integer                                    intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_g_num_s),          intent(in)  :: debut_bloc
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Ecriture parallele de donnees entrelacees
!
! arguments :
!   unite           <-- Unite du fichier
!   t_n_composantes <-- Indique si le nombre de composantes est constant
!                       A choisir entre :
!                             - pdm_io_n_composante_constant
!                             - pdm_io_n_composante_variable
!   n_composantes   <-- Nombre de composantes (tableau ou scalaire suivant 
!                                              la valeur de t_n_composantes)
!   taille_donnee   <-- Taille unitaire de la donnnee en octet
!   n_donnees       <-- Nombre de donnees a lire
!   indirection     <-- Indirection (Numerotation absolue)
!   donnees         <-- Donnees a ecrire
!  
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_ecr_par_entrelacee(unite, &
!                                          t_n_composantes, &
!                                          n_composantes, &
!                                          taille_donnee, & 
!                                          n_donnees, & 
!                                          indirection, &
!                                          donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite
!     integer,                                   intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     integer (kind = pdm_g_num_s), dimension(*), intent(in)  :: indirection
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Ecriture parallele de blocs continus de donnees
!
! arguments :
!   unite           <-- Unite du fichier
!   t_n_composantes <-- Indique si le nombre de composantes est constant
!                       A choisir entre :
!                             - pdm_io_n_composante_constant
!                             - pdm_io_n_composante_variable
!   n_composantes   <-- Nombre de composantes (tableau ou scalaire suivant 
!                                              la valeur de t_n_composantes)
!   taille_donnee   <-- Taille unitaire de la donnnee en octet
!   debut_bloc      <-- Index du debut du bloc
!   n_donnees       <-- Nombre de donnees a ecrire
!   donnees         <-- Donnees a ecrire
!  
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_ecr_par_bloc(unite, &
!                                    t_n_composantes, &
!                                    n_composantes, &
!                                    taille_donnee, & 
!                                    n_donnees, & 
!                                    debut_bloc, &
!                                    donnees)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite
!     integer                                    intent(in)  :: t_n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_composantes
!     integer (kind = pdm_l_num_s),           intent(in)  :: taille_donnee
!     integer (kind = pdm_g_num_s),          intent(in)  :: debut_bloc
!     integer (kind = pdm_l_num_s),           intent(in)  :: n_donnees
!     "Pas de type : void* en C" , dimension(*), intent(out) :: donnees 


!----------------------------------------------------------------------------
! Fermeture du fichier
!
! arguments :
!   unite           <-- Unite du fichier
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_close(unite)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite


!----------------------------------------------------------------------------
! Destruction de la structure fichier liee a l'unite
!
! arguments :
!   unite           <-- Unite du fichier
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_detruit(unite)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite


!----------------------------------------------------------------------------
! Affiche des informations sur le fichier
!
! arguments :
!   unite           <-- Unite du fichier
!
!----------------------------------------------------------------------------*/
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_dump(unite)
!
!     ! 
!     ! Arguments
! 
!     integer (kind = pdm_l_num_s),           intent(in)  :: unite


!----------------------------------------------------------------------------
! Retourne le temps cumule d'acces aux fichiers
!
! arguments :
!   unite           <-- Unite du fichier
!   t_cpu           --> Temps CPU
!   t_elapsed       --> Temps elapsed
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_timer_fichier (unite, &
!                                          t_cpu, &
!                                          t_elapsed)
!
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s), intent(in)  :: unite
!     double precision,                intent(out) :: t_cpu
!     double precision,                intent(out) :: t_elapsed


!----------------------------------------------------------------------------
! Retourne le temps cumule de distribution des donnees
!
! arguments :
!   unite           <-- Unite du fichier
!   t_cpu           --> Temps CPU
!   t_elapsed       --> Temps elapsed
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_timer_ditrib (unite, &
!                                         t_cpu, &
!                                         t_elapsed)
!
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s), intent(in)  :: unite
!     double precision,                intent(out) :: t_cpu
!     double precision,                intent(out) :: t_elapsed


!----------------------------------------------------------------------------
! Retourne le temps cumule total
!
! arguments :
!   unite           <-- Unite du fichier
!   t_cpu           --> Temps CPU
!   t_elapsed       --> Temps elapsed
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_timer_total (unite, &
!                                        t_cpu, &
!                                        t_elapsed)
!
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s), intent(in)  :: unite
!     double precision,                intent(out) :: t_cpu
!     double precision,                intent(out) :: t_elapsed


!----------------------------------------------------------------------------
! Retourne le temps cumule consacrÃ© au swap endian
!
! arguments :
!   unite           <-- Unite du fichier
!   t_cpu           --> Temps CPU
!   t_elapsed       --> Temps elapsed
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_timer_swap_endian(unite, &
!                                             t_cpu, &
!                                             t_elapsed)
!
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s), intent(in)  :: unite
!     double precision,                intent(out) :: t_cpu
!     double precision,                intent(out) :: t_elapsed


!----------------------------------------------------------------------------
! Retourne le numero de version de pdm_io ayant permis 
! de generer le fichier
!
! arguments :
!   unite  <-- Unite du fichier
!   majeur  --> Numero de version majeur
!   mineur  --> Numero de version mineur
!   release --> Numero de version release
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_version(unite, &
!                                   majeur, &
!                                   mineur, &
!                                   release)
!     ! 
!     ! Arguments
!
!     integer (kind = pdm_l_num_s), intent(in)  :: unite
!     integer                        , intent(out) :: majeur
!     integer                        , intent(out) :: mineur
!     integer                        , intent(out) :: release


!----------------------------------------------------------------------------
! Retourne le numero de version courant de pdm_io
!
! arguments :
!   majeur  --> Numero de version majeur
!   mineur  --> Numero de version mineur
!   release --> Numero de version release
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_version_courante(majeur, &
!                                            mineur, &
!                                            release)
!     ! 
!     ! Arguments
!
!     integer                        , intent(out) :: majeur
!     integer                        , intent(out) :: mineur
!     integer                        , intent(out) :: release


!----------------------------------------------------------------------------
! Retourne le numero de version de pdm_io ayant permis 
! de generer le fichier
!
! arguments :
!   unite  <-- Unite du fichier
!   annee   --> Annee
!   mois    --> mois
!   jour    --> jour
!   heure   --> heure
!   minute  --> minute
!   seconde --> seconde
!
!----------------------------------------------------------------------------
!
! Fonction definie en C dans pdm_io.c
!
!   subroutine pdm_io_get_date(unite, &
!                                annee, &
!                                mois, &
!                                jour, &
!                                heure, &
!                                minute, &
!                                seconde)
               
end module mod_pdm_io
