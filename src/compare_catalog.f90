program compare_catalog
    
    use ira_precision
    use ira_pbc, only: pbc_vec
    
  IMPLICIT NONE
  integer                   :: nref_catalog   ! size of the catalog
  integer                   :: nstruc_analyze   ! size of the trajectory
  integer                   :: nbegin, nend, n_list,idx 
  integer                   :: iref, istruc, iat, jat   ! Counters 
  character(len=500)        :: line
  real                      :: dmax, dist, dist_lowest
  integer, allocatable      :: list(:)
  LOGICAL                   :: fail
  CHARACTER(len=60)         :: file_catalog, file_analyze, file_out
  TYPE struct
  INTEGER, ALLOCATABLE  :: typ(:)
  REAL,  ALLOCATABLE    :: coords(:,:)
  REAL, DIMENSION(3,3)  :: lat
  REAL, DIMENSION(3,3)  :: beta 
  INTEGER               :: nat
  REAL                  :: dist_k
  ENDTYPE struct 
  TYPE(struct), ALLOCATABLE :: catalog (:)    ! Catalog of structure
  TYPE(struct), ALLOCATABLE :: analyze (:)    ! Structures to analyze
  TYPE(struct)              :: loc        ! local structure of global structure to analyze  
  
  CALL getarg(1, file_analyze)
  CALL getarg(2, file_out)
 
  dmax = 0.0
  
  ! ... Read catalog structures xyz
  file_catalog = 'catalog.xyz'
  OPEN(10,FILE=TRIM(file_catalog))
  
  READ(10,*) nref_catalog
  ALLOCATE( catalog(nref_catalog))
  
  DO iref= 1,nref_catalog
   READ(10,*) catalog(iref)%nat
   READ(10,*)
   ALLOCATE( catalog(iref)%typ( 1:catalog(iref)%nat )    )
   ALLOCATE( catalog(iref)%coords( 1:3, 1:catalog(iref)%nat )  )
   DO iat = 1, catalog(iref)%nat
    READ(10,*) catalog(iref)%typ(iat),   &
       catalog(iref)%coords(1,iat),&
       catalog(iref)%coords(2,iat),&
       catalog(iref)%coords(3,iat)
   ENDDO
   
   ! ... Assume first atom in is central and recenter
   DO iat = 1, catalog(iref)%nat
    catalog(iref)%coords(:,iat) = catalog(iref)%coords(:,iat) - catalog(iref)%coords(:, 1)
   ENDDO
   
   ! ... Find maximal distance of any atom from the central atom in ref.
   DO iat = 1, catalog(iref)%nat 
    dmax = max( dmax, norm2( catalog(iref)%coords(:,iat) - catalog(iref)%coords(:,1) ))
   ENDDO
   
   ! .. put ref in beta
   DO iat = 2, catalog(iref)%nat
    DO jat = 2, catalog(iref)%nat
     call set_orthonorm_bas( catalog(iref)%coords(:,iat), &
             catalog(iref)%coords(:,jat), &
             catalog(iref)%beta, fail)
     IF( .NOT. fail ) EXIT
    ENDDO
    IF( .NOT. fail ) EXIT
   ENDDO
   catalog(iref)%dist_k = max( norm2( catalog(iref)%coords(:,iat)), norm2( catalog(iref)%coords(:,jat)))
   catalog(iref)%dist_k = catalog(iref)%dist_k*2.5
   
   DO iat = 1, catalog(iref)%nat
    catalog(iref)%coords(:,iat) = matmul( catalog(iref)%beta, catalog(iref)%coords(:,iat))
   ENDDO
   
  ENDDO
  CLOSE(10)
  ! ... increase dmax slightly, just to be sure
  dmax = dmax * 1.4
  
  ! ... Read the structure that needs to be analyzed xyz
  OPEN(11,FILE=TRIM(file_analyze))
  READ(11,*) nstruc_analyze
  ALLOCATE( analyze(nstruc_analyze))
  DO istruc= 1,nstruc_analyze
   PRINT*, istruc
   READ(11,*) analyze(istruc)%nat
   READ(11, '(a500)') line
   nbegin = scan(line,'"')
   line = line(nbegin+1:)
   nend=scan(line,'"')
   line = line(:nend-1)
   READ(line,*) analyze(istruc)%lat

   ALLOCATE( analyze(istruc)%typ(1:analyze(istruc)%nat) )
   ALLOCATE( analyze(istruc)%coords(1:3,1:analyze(istruc)%nat) )
   DO iat = 1, analyze(istruc)%nat
    READ(11,*) analyze(istruc)%typ(iat),   &
       analyze(istruc)%coords(1,iat),&
       analyze(istruc)%coords(2,iat),&
       analyze(istruc)%coords(3,iat)
   ENDDO
  
  ENDDO
  CLOSE(11)
  PRINT*, "Fin de la lecture"
  
  ! ... Loop over sites of analyze struc, get loc and compare it to catalog to have smallest distance
  OPEN(12,FILE=TRIM(file_out))
  WRITE(12,*) 'Number of frames: ', nstruc_analyze
  WRITE(12,*) 'Number of refs:   ', nref_catalog
  DO istruc = 1, nstruc_analyze
   ALLOCATE( list(1:analyze(istruc)%nat) )
   PRINT*, istruc
   WRITE(12,'(a8,i10)') 'Frame=', istruc
   WRITE(12,'(a8,a5,a5,a10,a11)') 'Atom', 'Type', 'Ref', 'Distance', 'Neighbours'
   DO iat = 1, analyze(istruc)%nat
    dist_lowest= 999.0
    !
    !print*, iat
    ! ... Gget local structure around atom iat
    IF (analyze(istruc)%typ(iat) /= 999) then
     CALL set_atoms_dist( analyze(istruc)%nat, analyze(istruc)%coords,&
             analyze(istruc)%lat, iat, dmax, n_list, list )
     ALLOCATE( loc%typ(1:n_list+1))
     ALLOCATE( loc%coords(1:3,1:n_list+1))
    !print*, 'b'
     CALL get_atoms( analyze(istruc)%nat, analyze(istruc)%typ, analyze(istruc)%coords, analyze(istruc)%lat, &
       iat, n_list, list, loc%nat, loc%typ, loc%coords )
     ! ... Compare this local structure to the ones in catalog
     DO iref=1,nref_catalog
      !print*, 'c'
      CALL compare_Q( catalog(iref)%nat, catalog(iref)%typ, catalog(iref)%coords,&
          loc%nat, loc%typ, loc%coords, catalog(iref)%dist_k, dist )
      !print*, 'd'
      IF( dist .lt. dist_lowest) THEN
       dist_lowest = dist
       idx   = iref
      ENDIF
     ENDDO
     DEALLOCATE( loc%typ, loc%coords )
     
      IF (dist_lowest <= 5) THEN
       WRITE(12,'(i8,i5,i5,f10.7)',ADVANCE='NO') iat, analyze(istruc)%typ(iat), idx, dist_lowest
       DO jat = 2, n_list
        IF (analyze(istruc)%typ(list(jat)) /= 999 .AND. (list(jat)+42)/43 /= (list(jat-1)+42)/43) THEN
         WRITE(12,'(i8)',ADVANCE='NO') list(jat)
        ENDIF
       ENDDO
       WRITE(12,'(a)')
      ENDIF
    ENDIF 
   ENDDO
   DEALLOCATE( list )
  ENDDO
  WRITE(12,*) "Fin du fichier"
  CLOSE(12)
  PRINT*, "Fin de l'analyse"
  
  ! ... Final Cleaning
  DO istruc = 1, nstruc_analyze
   DEALLOCATE( analyze(istruc)%coords, analyze(istruc)%typ )
  ENDDO
  DO iref=1,nref_catalog
   DEALLOCATE( catalog(iref)%coords, catalog(iref)%typ )
  ENDDO
  DEALLOCATE ( catalog )
  
end program compare_catalog


subroutine set_atoms_dist( nat, coords, lat, isite, rcut, n_list, list )
  !! set list of indices based on distance
  use ira_precision
  use ira_pbc, only: pbc_vec
  implicit none
  integer,       intent(in) :: nat
  real, dimension(3,nat),  intent(in) :: coords
  real, dimension(3,3),  intent(in) :: lat
  integer,       intent(in) :: isite
  real,        intent(in) :: rcut
  integer,       intent(out) :: n_list
  integer, dimension(nat), intent(out) :: list
  !!
  integer :: i
  real, dimension(3) :: rij
  real :: dist
  !!
  !! initialize to 0
  !!
  list(:) = 0
  n_list = 0
  !!
  !! find norm of each vector
  !!

  do i = 1, nat
   if(i .eq. isite) cycle
   rij = coords(:,i) - coords(:,isite)
   !call cart_to_crist(rij,lat)
   !call periodic(rij)
   !call crist_to_cart(rij,lat)
   call pbc_vec(rij,lat)
   dist = sqrt(dot_product(rij,rij))
   if( dist .le. rcut ) then
    n_list = n_list + 1
    list(n_list) = i
   endif
  end do

end subroutine set_atoms_dist


subroutine get_atoms( nat, typ, coords, lat, isite, n_list, array, nat_loc, typ_loc, coords_loc )
  !! get atoms from list into *_loc arrays
  !! Attention: nat_loc is computed before entering in this routine, and
  !! the arrays typ_loc and coords_loc are not allocated here!
  !!
  !! input array is integer list of indices to be included
  use ira_precision
  use ira_pbc, only: pbc_vec
  implicit none
  integer,       intent(in) :: nat
  integer, dimension(nat),   intent(in) :: typ
  real, dimension(3,nat),  intent(in) :: coords
  real, dimension(3,3),    intent(in) :: lat
  integer,       intent(in) :: isite
  integer,       intent(in) :: n_list
  integer, dimension(n_list),   intent(in) :: array
  integer,        intent(out) :: nat_loc
  integer, dimension(n_list+1), intent(out) :: typ_loc
  real, dimension(3,n_list+1),  intent(out) :: coords_loc

  integer :: i

  nat_loc = n_list + 1
  !!
  !! put isite on first index
  coords_loc(:,1) = coords(:,isite)
  typ_loc(1) = typ(isite)
  !!
  !! get others from list
  coords_loc(:,2:) = coords(:, array(1:n_list) )
  typ_loc(2:) = typ( array(1:n_list) )
  !!
  !! recenter, periodic
  do i = 1, nat_loc
   coords_loc(:,i) = coords_loc(:,i) - coords(:,isite)
   !call cart_to_crist(coords_loc(:,i),lat)
   !call periodic(coords_loc(:,i))
   !call crist_to_cart(coords_loc(:,i),lat)
   call pbc_vec(coords_loc(:,i),lat)
  end do

end subroutine get_atoms


subroutine compare_Q( nat_ref, typ_ref_in, coords_ref_in, &
        nat_loc, typ_loc_in, coords_loc_in, &
        dist_k, hd_out )
  !!
  !! calculate hd, based on input structures.
  !! coords_ref is already in beta, so need to find
  !! gamma and then do svd.
  !! Assume central atom is first in list, and both strucs already centered.
  implicit none
  integer, intent(in) :: nat_ref
  integer, dimension(nat_ref), intent(in) :: typ_ref_in
  real, dimension(3,nat_ref), intent(in) :: coords_ref_in
  integer, intent(in) :: nat_loc
  integer, dimension(nat_loc), intent(in) :: typ_loc_in
  real, dimension(3,nat_loc), intent(in) :: coords_loc_in
  real, intent(in) :: dist_k
  real, intent(out) :: hd_out

  real :: some_thr
  real, dimension(3,3) :: gamma, svd_r
  integer, dimension(3) :: gamma_idx
  real, dimension(3) :: svd_t, rdum
  real :: dd
  integer :: i
  integer :: ierr
  integer, dimension(nat_loc) :: found
  real, dimension(nat_loc) :: dists
  !! local copies
  integer, dimension(nat_ref) :: typ_ref
  integer, dimension(nat_loc) :: typ_loc
  real, dimension(3,nat_ref) :: coords_ref
  real, dimension(3,nat_loc) :: coords_loc

  typ_ref = typ_ref_in
  coords_ref = coords_ref_in
  typ_loc = typ_loc_in
  coords_loc = coords_loc_in

  !!
  !! The next two lines make all atomic types identical. If you want to
  !! ignore the atomic types uncomment the next two lines. This will make it
  !! possible to compare the geometry of two chemically different structures.
  !!
  ! typ_ref(:) = 1
  ! typ_loc(:) = 1

  some_thr = 9999.9
  !print*, 'aa'

  call get_gamma_m( nat_ref, typ_ref, coords_ref, &
   nat_loc, typ_loc, coords_loc, &
   dist_k, gamma, gamma_idx, dd, some_thr)
  !!
  !! apply gamma
  !!
  !print*, 'ab'
  do i = 1, nat_ref
   coords_ref(:,i) = matmul(transpose(gamma),coords_ref(:,i))
  end do
  !!
  !! find permutation
  !!
  !print*, 'ac'
  call cshda( nat_ref, typ_ref, coords_ref, &
   nat_loc, typ_loc, coords_loc, &
   some_thr, found(1:nat_loc), dists(1:nat_loc) )

  !print*, 'ad'
  !print*, found
  !! apply permutation
  if (found(1) .ne. 0) then
   typ_loc(:) = typ_loc(found(:))
   coords_loc(:,:) = coords_loc(:,found(:))
   !!
   !! call svd
   !!
   !print*, 'ae'
   call svdrot_m( nat_ref, typ_ref, coords_ref, &
    nat_ref, typ_loc(1:nat_ref), coords_loc(:,1:nat_ref), &
    svd_r, svd_t, ierr )

   !print*, 'af'
   !! apply svd
   do i = 1, nat_ref
    coords_ref(:,i) = matmul( svd_r, coords_ref(:,i) ) + svd_t
   end do


   !!
   !! measure hd
   !!
   !print*, 'ag'
   dd = 0.0
   do i = 1, nat_ref
    rdum = coords_ref(:,i) - coords_loc(:,i)
    dd = max( norm2(rdum), dd )
   end do

   hd_out = dd

  endif

  ! write(*,*) nat_loc + nat_ref
  ! write(*,*) dd
  ! do j = 1, nat_ref
  !  write(*,*) 1, coords_ref(:,j)
  ! end do
  ! do j = 1, nat_loc
  !  write(*,*) 2, coords_loc(:,j)
  ! end do

end subroutine compare_Q


