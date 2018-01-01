!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
! Module containing routines for fft with a custom energy cutoff
!--------------------------------------------------------------------
!
MODULE fft_custom

  USE kinds, ONLY: DP
  USE parallel_include
  
  USE fft_types, ONLY: fft_type_descriptor
  USE fft_ggen, ONLY: fft_set_nl, fft_set_nlm
  
  IMPLICIT NONE

  TYPE fft_cus
  
     ! ... data structure containing information about "custom" fft grid:
     ! ... G-vectors and the like - FIXME: to be aligned with QE style

     REAL(kind=DP) :: ecutt
     ! Custom cutoff (rydberg)
     REAL(kind=DP) :: gcutmt
     INTEGER :: ngmt_g
     REAL(kind=DP), DIMENSION(:), POINTER :: ggt
     REAL(kind=DP), DIMENSION(:,:),POINTER :: gt
     INTEGER :: gstart_t
     INTEGER :: npwt
     LOGICAL :: initialized = .FALSE.
     
  END TYPE fft_cus

!--------------------------------------------------------------------
CONTAINS
!=----------------------------------------------------------------------------=!

  !-----------------------------------------------------------------------
  SUBROUTINE ggenx( ngm, g, comm, dfftt, fc )
  !-----------------------------------------------------------------------
    !
    ! Initialize g-vectors for custom grid, in exactly the same ordering
    ! as for the dense and smooth grids
    !
    ! FIXME: Should be merged with similar case for dfftp/dffts
    !
    !--------------------------------------------------------------------
    !
    USE kinds,              ONLY : DP
    USE cell_base,          ONLY : at, tpiba2
    USE control_flags,      ONLY : gamma_only
    USE constants,          ONLY : eps8
    USE mp,                 ONLY : mp_max, mp_sum
    
    IMPLICIT NONE
    ! number of input G-vectors
    INTEGER, INTENT(IN) :: ngm
    ! G-vectors in FFT grid
    REAL(dp), INTENT(IN) :: g(3,ngm)
    ! communicator of the group on which g-vecs are distributed
    INTEGER, INTENT(IN) :: comm
    TYPE (fft_type_descriptor), INTENT(INOUT) :: dfftt 
    TYPE(fft_cus), INTENT(INOUT) :: fc
    !
    INTEGER,  DIMENSION(:), ALLOCATABLE :: mill(:,:)
    INTEGER :: ngmt, n1, n2, n3, i
    !
    ngmt = dfftt%ngm
    !
    !  calculate fc%ngmt_g, sum over all processors
    !
    fc%ngmt_g = ngmt
    CALL mp_sum( fc%ngmt_g, comm )
    !
    !  allocate arrays
    !
    ALLOCATE( fc%ggt(ngmt) )
    ALLOCATE( fc%gt (3, ngmt) )
    !  
    ! fc%npwt = number of PW in sphere of radius ecutwfc (useful for Gamma)
    !
    fc%npwt=0
    !
    DO i = 1, ngmt
       !
       fc%gt(:,i) = g(:,i)
       !
       ! compute fc%npwt
       !
       fc%ggt(i) = SUM(fc%gt (1:3,i)**2)
       IF ( fc%ggt(i) <= fc%ecutt / tpiba2) fc%npwt = fc%npwt + 1
       !
    END DO
    !
    !     determine first nonzero g vector
    !
    IF (fc%ggt(1).LE.eps8) THEN
       fc%gstart_t=2
    ELSE
       fc%gstart_t=1
    ENDIF
    !
    ALLOCATE( mill(3,ngmt) )
    !
    CALL fft_set_nl ( dfftt, at, g, mill  )
    IF ( gamma_only) CALL fft_set_nlm (dfftt, mill)
    !
    !     Miller indices no longer needed
    !
    DEALLOCATE ( mill )
    !
  END SUBROUTINE ggenx
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE ggent(ngm_, comm, dfftt, fc)
    !--------------------------------------------------------------------
    !
    ! Initialize g-vectors for custom grid
    !
    ! FIXME: Should be merged with ggen
    !
    USE kinds,              ONLY : DP
    USE cell_base,          ONLY : at, bg, tpiba2
    USE control_flags,      ONLY : gamma_only
    USE constants,          ONLY : eps8
    USE mp,                 ONLY: mp_max, mp_sum
    
    IMPLICIT NONE
    
    TYPE (fft_type_descriptor), INTENT(INOUT) :: dfftt 
    TYPE(fft_cus), INTENT(INOUT) :: fc
    INTEGER, INTENT(IN) :: ngm_
    INTEGER, INTENT(IN) :: comm  ! communicator of the group over which
                                 ! g-vectors are distributed
    !
    INTEGER,  DIMENSION(:), ALLOCATABLE :: mill(:,:)
    INTEGER :: ngmt, ngmx, n1, n2, n3, n1s, n2s, n3s
    REAL(DP) ::  t (3), tt, swap
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing all g vectors, on all processors: replicated data
    INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors:
    !     replicated data
    INTEGER, ALLOCATABLE :: igsrt(:)
    !
    INTEGER :: m1, m2, mc
    INTEGER :: i, j, k, ipol, ng, igl, iswap, indsw, ni, nj, nk
    !    
    ngmt = ngm_
    !
    !  calculate sum over all processors
    !
    fc%ngmt_g = ngm_
    CALL mp_sum( fc%ngmt_g, comm )
    !
    !  allocate arrays - only those that are always kept until the end
    !
    ALLOCATE( fc%ggt(ngmt) )
    ALLOCATE( fc%gt (3, ngmt) )
    !
    ALLOCATE( mill_g( 3, fc%ngmt_g ) )
    ALLOCATE( mill_unsorted( 3, fc%ngmt_g ) )
    ALLOCATE( igsrt( fc%ngmt_g ) )
    ALLOCATE( g2sort_g( fc%ngmt_g ) )
   
    g2sort_g(:) = 1.0d20
    !
    ! save present value of ngm in ngmx variable
    !
    ngmx = ngmt
    !
    ngmt = 0
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftt%nr1-1)/2
    nj = (dfftt%nr2-1)/2
    nk = (dfftt%nr3-1)/2
    !
    iloop: DO i = -ni, ni
       !
       ! gamma-only: exclude space with x < 0
       !
       IF ( gamma_only .AND. i < 0) CYCLE iloop
       jloop: DO j = -nj, nj
          !
          ! gamma-only: exclude plane with x = 0, y < 0
          !
          IF ( gamma_only .AND. i == 0 .AND. j < 0) CYCLE jloop
          kloop: DO k = -nk, nk
             !
             ! gamma-only: exclude line with x = 0, y = 0, z < 0
             !
             IF ( gamma_only .AND. i == 0 .AND. j == 0 .AND. k < 0) CYCLE kloop
             t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
             tt = SUM(t(:)**2)
             IF (tt <= fc%gcutmt) THEN
                ngmt = ngmt + 1
                IF (ngmt > fc%ngmt_g) CALL errore ('ggent', 'too many g-vectors', ngmt)
                mill_unsorted( :, ngmt ) = (/ i,j,k /)
                IF ( tt > eps8 ) THEN
                   g2sort_g(ngmt) = tt
                ELSE
                   g2sort_g(ngmt) = 0.d0
                ENDIF
             ENDIF
          ENDDO kloop
       ENDDO jloop
    ENDDO iloop
    
    IF (ngmt  /= fc%ngmt_g ) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(ngmt - fc%ngmt_g))

    igsrt(1) = 0
    CALL hpsort_eps( fc%ngmt_g, g2sort_g, igsrt, eps8 )
    mill_g(1,:) = mill_unsorted(1,igsrt(:))
    mill_g(2,:) = mill_unsorted(2,igsrt(:))
    mill_g(3,:) = mill_unsorted(3,igsrt(:))
    DEALLOCATE( g2sort_g, igsrt, mill_unsorted )
    ngmt = 0
    
    ngloop: DO ng = 1, fc%ngmt_g

       i = mill_g(1, ng)
       j = mill_g(2, ng)
       k = mill_g(3, ng)
       
       IF ( dfftt%lpara ) THEN
          m1 = MOD (i, dfftt%nr1) + 1
          IF (m1 < 1) m1 = m1 + dfftt%nr1
          m2 = MOD (j, dfftt%nr2) + 1
          IF (m2 < 1) m2 = m2 + dfftt%nr2
          mc = m1 + (m2 - 1) * dfftt%nr1x
          IF ( dfftt%isind ( mc ) == 0) CYCLE ngloop
       END IF
       
       ngmt = ngmt + 1
       
       !  To map local (ngmt) and global (ng) G-vector index: 
       !     fc%ig_l2gt( ngmt ) = ng
       !  The global G-vector arrangement depends on the number of processors
       !
       
       fc%gt (1:3, ngmt) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
       fc%ggt (ngmt) = SUM(fc%gt (1:3, ngmt)**2)
       
       IF (ngmt > ngmx) CALL errore ('ggent', 'too many g-vectors', ngmt)
    ENDDO ngloop

    DEALLOCATE( mill_g )

    IF (ngmt /= ngmx) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(ngmt - ngmx))
    !
    !     determine first nonzero g vector
    !
    IF (fc%ggt(1).LE.eps8) THEN
       fc%gstart_t=2
    ELSE
       fc%gstart_t=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    ALLOCATE( mill(3,ngmt) )
    !
    CALL fft_set_nl ( dfftt, at, fc%gt, mill  )
    IF ( gamma_only) CALL fft_set_nlm (dfftt, mill)
    !
    !     Miller indices no longer needed
    !
    DEALLOCATE ( mill )
    !
    ! set npwt - it should eventually be calculated somewhere else with 
    ! n_plane_waves() but it is good enough for gamma_only case

    IF(gamma_only) THEN
       fc%npwt=0
       DO ng = 1, ngmt
          tt = (fc%gt (1, ng) ) **2 + (fc%gt (2, ng) ) **2 + (fc%gt&
               & (3, ng) ) **2
          IF (tt <= fc%ecutt / tpiba2) THEN
             !
             ! here if |k+G|^2 <= Ecut increase the number of G
             !  inside the sphere
             !
             fc%npwt = fc%npwt + 1
          ENDIF
       ENDDO
    ENDIF

    RETURN
    !    
  END SUBROUTINE ggent
  
  SUBROUTINE deallocate_fft_custom(fc)
    !this subroutine deallocates all the fft custom stuff
    USE fft_types, ONLY : fft_type_deallocate
    
    IMPLICIT NONE

    TYPE(fft_cus) :: fc

    IF(.NOT. fc%initialized) RETURN

    IF ( ASSOCIATED (fc%gt)  )  DEALLOCATE(fc%gt)
    IF ( ASSOCIATED (fc%ggt) )  DEALLOCATE(fc%ggt)
    fc%initialized=.FALSE.

    RETURN

  END SUBROUTINE deallocate_fft_custom
  !----------------------------------------------------------------------------  
END MODULE fft_custom
