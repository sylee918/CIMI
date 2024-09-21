! Modification history:
! March 14, 2002
!   * Rename subroutine Circle to Cirular to avoid having the same name as the
!     subroutine CIRCLE in t96_01.f.
! February 20, 2008
!   * Add dealloacte arrays just before the end of subroutine Gmresm to 
!     avoid memory problem.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   Potential solver package:
!
!
    SUBROUTINE Gmresm (bndloc, c, v, idim, jdim, jj1, jj2, ncoeff)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: idim, jdim, jj1, jj2, ncoeff
    REAL, INTENT (IN) :: c (ncoeff,idim,jdim), bndloc (jdim)
    REAL, INTENT (IN OUT) :: v (idim,jdim)
!
!
!     In this subroutine we will solve difference equations on a part of
!     the RCM grid. The part of the grid where V is found is defined as
!     jj1 <= j <= jj2,  MINVAL(min_j) <= i <= idim, where jj1 and jj2 are
!     RCM variables from common block NDICES. On this part of the grid,
!     we will renumber the grid points and re-write the RCM difference
!     equations as A*X=B, where X is the vector of unknowns, A is a matrix
!     and B is a vector both formed from C1-C5 coefficients.
!
!     Since A is sparse, it is stored in the CRS format (compressed-row).
!     Linear system A*X=B is solved iteratively by GMRES(m) algorithm
!     (Generalized Minimized Residuals with restarts).
!
!
!     Local variables:
!       nmax  : max # of grid points to treat
!       M     : max # of Krylov vectors (when exceeded, goes into restart)
!       itermx: max # of restarts allowed
!       ilow  : min value of min_j over all J
!       ni    : I-size (# of pts) of the square part of grid treated here.
!       nj    : J-size of the square part of grid treated
!       nij   : total number of grid points over which method loops.
!       nnz   : the total number of nonzero elements of matrix A
!       nzmax : max # of non-zero coefficients of A.
!       amatrx: non-zero elements of A ordered into a linear sequence
!       colind: vector holding column numbers for each element of AMATRX
!       ROWPTR: vector with locations of 1st element of A in each row in AMATRX
!       DIAPTR: vector with locations of diagonal elements of A in AMATRX
!       PIVOTS: vector with inverses of diagonal of the preconditioner matrix
!       TOL   : relative error. Empirically, 1e-5 is equiv. to
!               sum of resids < 1 Volt
!       H     : Hessenberg matrix holding dot products
!       B     : holds right-hand side vector of linear system to be solved
!       CS, SN: vectors holding Givens rotations parameters
!       X0    : vector holding initial approx. on entry and solution on exit
!       RESID : vector of residuals (b-A*x)
!       AX, W : work vectors
!       Y     : vector holding coefficients of the solution expanded in the
!               Krylov vectors
!       S     : initially unit vector, then rotated by Givens rotations
!
!       This algorithm uses preconditioning based on incomplete LU
!       factorization. Namely, it uses D-ILU type.
!

    INTEGER, PARAMETER :: M = 300, itermx = 300
    INTEGER            :: nzmax, nij, nnz
    INTEGER, ALLOCATABLE   :: row_ptr (:), i_column (:), diag_ptr(:)
    REAL   , ALLOCATABLE   :: a_mtrx(:), b_mtrx(:), pivots(:), &
                                      x0(:), x(:,:), resid(:), w (:), ax (:)
    REAL , POINTER :: window (:)
    REAL , TARGET  :: h (m+1,m), s (m+1)
    REAL           :: sn (m), cs (m), y (m)
    INTEGER :: nj, ni, krow, ii, jj, i, jkryl, irstrt
!
    INTEGER :: imin
    REAL, PARAMETER :: tol_gmres = 1.0E-5, machine_tiny = TINY (1.0)
!
      REAL :: bnorm, rnorm, relerr
!
!     REAL :: error_s, error_1
!
!  1. Arrange difference equations into a matrix form A*x=b. This call returns
!     A as (AMATRX, COLIND, ROWPTR), and also B (RHS) and X0:
!
      CALL Gmresm_define_matrix ( )
!
!
! 2.  Compute preconditioner:
!
      CALL Gmresm_compute_DILU ()
!
!
!
!  3. Now begin GMRES(m) algorithm to solve A*x=b with X0 as init. approx.
!
      bnorm = SQRT ( DOT_PRODUCT (b_mtrx, b_mtrx))
      IF (ABS(bnorm) < machine_tiny) bnorm = 1.0
!
!
!     Begin GMRES(m) iterations (restarts):
!
      Restart_loop: DO irstrt = 1, itermx
!
!
!        ... Compute the norm of initial residuals:
!
         ax     = b_mtrx - Gmresm_Mtrx_times_vect (X0)
         resid  = Gmresm_Msolve (ax )
         rnorm  = SQRT (DOT_PRODUCT (resid, resid ))
         relerr = rnorm / bnorm
! Stan Sazykin's correction on 16 April 2009
      !  IF (relerr < TOL_gmres ) RETURN     ! V already holds solution
         IF (relerr < tol_gmres ) THEN
            IF (irstrt == 1) RETURN ! V already holds solution
            EXIT restart_loop
         ENDIF
!
!        .. Set 1st Krylov vector to R/||R||:
!
         x (:, 1) = resid / rnorm
!
!        .. Set up unit vector E1 of length RNORM:
!
         s (1) = rnorm
         s (2:m+1) = 0.0
!
!        .. Loop to generate orthonormal vectors in Krylov subspace:
!
         iterate_loop: DO jkryl = 1, M
!
!            ... Compute A*X(Jkryl) and solve M*w=A*X(kryl) for w:
!
             ax = Gmresm_Mtrx_times_vect (x (:,jkryl) )
             w  = Gmresm_Msolve ( ax )
!
!            ... Form J-th column of H-matrix and X (Jkryl+1)
!                (modified Gramm-Schmidt process):
!
            DO i = 1, jkryl
               H (i,jkryl)   = DOT_PRODUCT ( w , x (:,i) )
               w             = w  - h (i,jkryl) * x (:,i)
            END DO
            h (jkryl+1,jkryl)  = SQRT (DOT_PRODUCT (w, w))
            x (:, jkryl+1)     = w / h(jkryl+1,jkryl)
!
!
!           .. Update QR-factorization of H. For that, 
!           .... first, apply 1, ..., (Jkryl-1)th rotations
!                to the new (Jkryl-th) column of H:
!
            DO i = 1, Jkryl-1
               window => h (i:i+1, jkryl)
               window = Gmresm_Rotate_vector ( window, cs (i), sn (i) )
            END DO
!
!           .... second, compute the Jkryl-th rotation that
!                will zero H (jkryl+1,jkryl):
!
            window => h (jkryl:jkryl+1, jkryl)
            CALL Gmresm_Get_rotation ( window, cs (jkryl), sn (jkryl) )
!
!           .... third, apply Jkryl-th rotation to Jkryl-th column of H
!                and to S (rhs):
!
            window => h (jkryl:jkryl+1, jkryl)
            window = Gmresm_Rotate_vector ( window, cs (jkryl), sn (jkryl) )
            h (jkryl+1,jkryl) = 0.0
!
            window => s (jkryl : jkryl+1)
            window = Gmresm_Rotate_vector ( window, cs (jkryl), sn (jkryl) )
!
!
!           .. Approximate the norm of current residual:
!
            relerr = ABS (s (jkryl+1)) / bnorm
!
            IF (relerr < TOL_gmres) THEN
!
!              .. Residual is small, compute solution, exit:
!
               y(1:Jkryl) = Gmresm_Solve_upper_triang (A = h, B_RHS = s, N = Jkryl)
               DO i = 1, Jkryl
                  x0  = x0 + y(i)* X(:,i)
               END DO
               EXIT restart_loop
!
            END IF
!
         END DO iterate_loop
!
!
!        We got here because after a maximum number of Krylov vectors
!        was reached, approximated norm of residual was not small enough.
!        However, need to compute approx solution and check the actual norm
!        of residual (because the approx. norm may not be accurate due to
!        round offs):
!
         y(1:m) = Gmresm_Solve_upper_triang (A = h, B_RHS = s, N = m)
         DO i = 1, m
            x0 = x0 + y(i)* X(:,i)
         END DO
! 
         resid = b_mtrx - Gmresm_Mtrx_times_vect ( x0 )
         rnorm   = SQRT (DOT_PRODUCT (resid, resid))
         relerr = rnorm / bnorm
!
!        .. If the actual norm of residual is indeed small, exit:
! 
         IF (relerr < TOL_gmres) EXIT restart_loop
!
!        .. If not, continue by restarting...
!
      END DO restart_loop
!
!
!     Finished GMRES(m) loop. We get here either because
!     the solution was found, or because maximum number of
!     iterations was exceeded. Check for this:
!     
      IF (relerr >= TOL_gmres) THEN
         STOP 'convergence in GMRES(m) not achieved, stopping'
      END IF
!
!
!  4. Solution was found. The final step is to decode solution and put it
!     back into V array.
!
      DO jj = jj1, jj1+nj-1
         DO ii = imin, idim
            krow = ni*(jj-jj1) + (ii-imin+1)
            IF (ii >= CEILING(bndloc(jj))) v (ii,jj) = X0 (krow)
         END DO
      END DO
!
!     CALL Circle ( v, idim, jdim, jj1 )
      CALL Circular ( v, idim, jdim, jj1 )       ! MCF, March 14, 2002
!
!
! ************* Residual check ********************
!  error_s = 0.
!  DO jj = jj1, jj2
!     DO ii = min_j(jj), idim -1
!        error_1 = v (ii,jj) - c(1,ii,jj)*v(ii+1,jj) &
!                  -c(2,ii,jj)*v(ii-1,jj)&
!                  -c(3,ii,jj)*v(ii,jj+1)&
!                  -c(4,ii,jj)*v(ii,jj-1) - c(5,ii,jj)
!        error_s = error_s + ABS (error_1)
!     END DO
!     ii = idim
!     error_1 = v(ii,jj) - c(2,ii,jj)*v(ii-1,jj)&
!               -c(3,ii,jj)*v(ii,jj+1)-c(4,ii,jj)*v(ii,jj-1)&
!               -c(5,ii,jj)
!     error_s = error_s + ABS (error_1)
!  END DO
!  WRITE (*, &
!  &'(A11,ES9.2, 2X, A7,ES9.2, 2X, A14,ES9.2, 2X, A5,I3, 2X,&
!  & A2,I3)')  &
!        'SUM(resid)=', error_s, &
!        'RNORM2=', rnorm, &
!        'RNORM2/BNORM2=',rnorm/bnorm, &
!        'ITER=', irstrt,  &
!        'J=', jkryl
!
      IF (ALLOCATED (a_mtrx)) DEALLOCATE (a_mtrx)       ! Added by MCF on
      IF (ALLOCATED (b_mtrx)) DEALLOCATE (b_mtrx)       ! February 20, 2008
      IF (ALLOCATED (i_column)) DEALLOCATE (i_column)   !
      IF (ALLOCATED (diag_ptr)) DEALLOCATE (diag_ptr)   !
      IF (ALLOCATED (row_ptr)) DEALLOCATE (row_ptr)     !
      IF (ALLOCATED (pivots)) DEALLOCATE (pivots)       !
      IF (ALLOCATED (ax)) DEALLOCATE (ax)               !
      IF (ALLOCATED (w)) DEALLOCATE (w)                 !
      IF (ALLOCATED (resid)) DEALLOCATE (resid)         !
      IF (ALLOCATED (x0)) DEALLOCATE (x0)               !
      IF (ALLOCATED (x)) DEALLOCATE (x)                 !
!
      RETURN
!
!
    CONTAINS  !----------------------------------------------------
!
!
!
      SUBROUTINE Gmresm_Define_matrix ( )
      IMPLICIT NONE
!_____________________________________________________________________________
!     Subroutine returns matrix A stored in 3 vectors AMATRX, COLIND, ROWPTR
!     also NNZ is the number of non-zero elements of A, and DIAPTR vector
!     holds locations of the diagonal elements of A in AMATRX.
!
!     This subroutine will compute:
!     -- nij, size of smallest rect. grid area enclosing the modeling region,
!     -- nzmax, upper limit on the number of non-zero coeffs of matrix A,
!     -- nnz, actual number of non-zero elements of A,
!     -- b_mtrx(1:nij), right-hand side of system of linear equations,
!     -- a_mtrx(1:nnz)--matrix A of linear system encoded in CRS,
!     -- row_ptr(1:nij+1), i_column(1:nnz), diag_ptr (1:nij), encoding of A in CRS,
!
!     -- also allocate pivots (1:nzmax) for the pre-conditioner.
!_____________________________________________________________________________
!
!     Local variables:
!
      INTEGER :: i, j, L, krow
!
!
!     In this subroutine we take the coefficients of the RCM difference
!     equations approximating the MI-coupling PDE, and reformulate these
!     equations as to cast them into a linear system A*X=B, where A is
!     an NxN square matrix, X is the unknown vector whose elements are
!     unknown values of the potential on grid points V(i,j), and B is the
!     right-hand-side vector.
!     Apparently, such reformulation requires: (1) to number all grid points
!     in the modeling region sequentially into a 1-dimensional sequence, and
!     then (2) to form A from c1-c4 and B from c5 RCM coefficients.
!     As A is going to be sparse, an additional task is to store (encode) A
!     in the Compressed-Row-Storage format for using in the potential solver.
!
!     Matrix A is stored in one-dim REAL array AMATRX and two INTEGER 1-dim
!     arrays COLIND and ROWPTR. We simply go along each row of A starting with
!     the 1st row, then 2nd, etc, and for each non-zero element a(p,q), we
!     write AMATRX(L)=a(p,q), COLIND(L)=q, and L-index numbers those non-zero
!     elements sequentially. ROWPTR(p) has the L-index of where p-th row
!     starts in AMATRX.
!
!  1. Numbering grid points into a 1-dim. sequence.
!     Imagine the RCM grid as extending vertically in I from I=1 (highest lat,
!     top) to I=IDIM (lowest lat., bottom) and horizontally in J from J=jj1
!     (noon, left) to J=JJ2 (last point before noon, right). If MIN_J(j) gives
!     the first I-point inside the modeling region, then we will consider all
!     grid points (i,j) such that ilow <= i <= idim, jj1 <= j <= jj2, where
!     ilow is MINVAL(min_j(:)). This will result in inclusion of some points
!     that are outside the modeling region, but we will define the difference
!     equations for them such that they don't matter.
!
!     The rectangular region of the grid we treat has the size NI by NJ, with
!     total of NIJ points in it:
!
      imin = MINVAL (CEILING(bndloc))
      nj  = jj2 - jj1 + 1
      ni  = idim - imin + 1
      nij = ni*nj
      nzmax = nij * ncoeff
      CALL Gmresm_Make_storage ( )
!
!
!     Number all points and store only non-zero elements of A. Order grid
!     points along J-lines from ILOW to IDIM, occasionally including points
!     outside the modeling region. Each point (i,j)
!     has number KROW (so that coefficients of the difference equation on that
!     point are on the krow-th row of A).
!
      L = 0
      DO j = jj1, jj2
      DO i = imin, idim
!
         krow          = ni * (j - jj1) + (i-imin + 1)
         b_mtrx (krow) = c (5,i,j)   ! this will be RHS vector
         X0 (krow)     = v (i,j) ! initial approximation taken from prev solution
         row_ptr(krow) = L + 1
!
         IF (i < CEILING(bndloc(j))) THEN
!
!           we are outside the main modeling region. In this case, the value
!           of the potential at this point is irrelevant, but we need to define
!           coefficients so as to proceed with this point as efficiently as possible.
!           Therefore, make V(i,j) a solution; that is, add a difference equation
!           V (i,j) = V (i,j)_on_input; the row of matrix A is then:
!
!           . . . . . .   1 . . . . . . . . . .
!           and RHS is just V(i,j)
!
            L               = L + 1
            a_mtrx (L)      = 1.0
            i_column (L)    = krow
            b_mtrx (krow)   = v (i,j)
            diag_ptr (krow) = L
!
         ELSE IF (j == jj1) THEN   ! first NI rows of A:
!
!           we are on the left (J=JJ1) side of the grid, and j-1 neighbors
!           wrap around from the other side (j=jj2)--periodic boundary
!           condition in J. So for these first NI rows of the matrix A,
!           c4 is the last coefficient. The first NI matrix rows look like:
!
!            1 -c1   0 .............-c3........ . . .      -c4...........0
!          -c2  1  -c1 ............. 0  -c3.... . . .       0  -c4.......0
!            0  -c2  1 -c1.................-c3. . . .       .......-c4....0
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!                   -c2 1 ....................... -c3  . .  ..........  -c4
!
            IF (i /= imin) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (2,i,j)
               i_column (L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = 1.0
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= idim) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  krow + NI
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            i_column (L) =  nij - NI + krow
!
         ELSE IF (j == jj2) THEN  ! last NI rows of A:
!
!           we are on the right (J=JJ2) side of the grid, periodic boundary
!           condition means c3-neighbors will come from J=JJ1, so that c3
!           coefficients will be first, c4-last on each row of matrix A. The
!           last NI rows of matrix A look like:
!
!          -c3  0 ..... . . .     -c4............  1  -c1 .................
!           0  -c3..... . . .      0 -c4......... -c2  1  -c1 .............
!           ......-c3.. . . .      .....-c4......  0  -c2  1  -c1 .........
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .. .
!           ...............-c3. . . ............-c4 .................-c2  1
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  i - imin + 1
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            I_column (L) =  krow - NI
!
            IF (i /= imin) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (2,i,j)
               i_column (L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = 1.0
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= idim) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
         ELSE          ! the rest of A:
!
!           we are neither on the left (j=jj1) nor on the right (j=jj2) side
!           of the grid, so don't need to worry about periodic boundary
!           conditions in J. RCM difference equation looks like:
!           V(i,j)=c(1,i,j)*V(i+1,j)+c(2,i,j)*V(i-1,j)+c(3,i,j)*V(i,j+1)+
!                  c(4,i,j)*V(i,j-1)+c(5,i,j);
!           except that if i=idim, then there is no C(1) term, and
!           if i=min_j(j), there is no c(2) term.
!           For such a grid point, the row of A looks like:
!
!           . . . . . -c4 . . . . . . -c2  1  -c1 . . . . . . -c3 . . . .
!
!           except that there is no c2 if i=ilow and no c1 if i=idim.
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            i_column (L) =  krow - ni
!
            IF (i /= imin) THEN
               L           =  L + 1
               a_mtrx(L)   = -c (2,i,j)
               i_column(L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = 1.0
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= idim) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  krow + ni
!
         END IF
      END DO
      END DO
!
      nnz = L ! number of non-zero elements of A
      row_ptr (nij+1) = nnz + 1    ! by definition of CRS format.
!
      RETURN
      END SUBROUTINE Gmresm_Define_matrix
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Make_storage ( )
      IMPLICIT NONE
      IF (ALLOCATED (a_mtrx)) DEALLOCATE (a_mtrx)
      IF (ALLOCATED (b_mtrx)) DEALLOCATE (b_mtrx)
      IF (ALLOCATED (i_column)) DEALLOCATE (i_column)
      IF (ALLOCATED (diag_ptr)) DEALLOCATE (diag_ptr)
      IF (ALLOCATED (row_ptr)) DEALLOCATE (row_ptr)
      IF (ALLOCATED (pivots)) DEALLOCATE (pivots)
      IF (ALLOCATED (ax)) DEALLOCATE (ax)
      IF (ALLOCATED (w)) DEALLOCATE (w)
      IF (ALLOCATED (resid)) DEALLOCATE (resid)
      IF (ALLOCATED (x0)) DEALLOCATE (x0)
      IF (ALLOCATED (x)) DEALLOCATE (x)
      ALLOCATE ( a_mtrx (nzmax), b_mtrx (nij), pivots (nij), &
                 i_column (nzmax), row_ptr (nij+1), diag_ptr (nij), &
                 ax (nij), w (nij), resid (nij), x0(nij), x (nij,m+1) )
      RETURN
      END SUBROUTINE Gmresm_Make_storage
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Compute_DILU ()
      IMPLICIT NONE
!_____________________________________________________________________________
!     Compute the preconditioner M. Matrix A is split as A = L_a + D_a + U_a
!     (strictly-lower triangular, diagonal and strictly-upper triangular).
!     Then M = L * U = (D + L_a) * D^(-1) * (D + U_a), so only need to find
!     and store D (one diagonal). In fact, PIVOTS holds inverses of D since
!     will divide by them later. D-ILU preconditioner M is kept in PIVOTS.
!     All structures are accessed from the host subroutine.
!
!     This subroutine only modifies (computes) PIVOTS (1:nij).
!_____________________________________________________________________________
!
      INTEGER  :: irow, jcol, krow
      REAL     :: element
      LOGICAL         :: found
!
      pivots  = a_mtrx (diag_ptr )
      DO irow = 1, nij
         pivots (irow) = 1.0 / pivots (irow)
         DO jcol = diag_ptr(irow)+1, row_ptr(irow+1)-1
            found = .FALSE.
            DO krow = row_ptr (i_column (jcol)), diag_ptr (i_column (jcol)) - 1
               IF (i_column (krow) == irow) THEN
                  found = .TRUE.
                  element = a_mtrx (krow)
               END IF
            END DO
            IF (found) THEN
               pivots (i_column (jcol)) = &
                  pivots (i_column (jcol)) - element*pivots(irow)*a_mtrx(jcol)
            END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE Gmresm_Compute_DILU
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Mtrx_times_vect (x)
      IMPLICIT NONE
      REAL , INTENT (IN)         :: x (:)
      REAL , DIMENSION (SIZE(x)) :: Gmresm_Mtrx_times_vect
!____________________________________________________________________________
!     subroutine to form matrix-vector product. Matrix A of size NxN
!     is assumed to be sparse with NNZ non-zero elements and is stored
!     in the compressed-row (CRS) format.
!     We compute y = A*x, where y and x are both vectors of length NNZ.
!     On entry, X holds x, and on exit Y is the result.
!____________________________________________________________________________
!
      INTEGER  :: i, j
!
      DO i = 1, SIZE (x)
         Gmresm_Mtrx_times_vect (i) = 0.0
         DO j = row_ptr (i), row_ptr (i+1)-1
            Gmresm_Mtrx_times_vect (i) = &
               Gmresm_Mtrx_times_vect (i) + a_mtrx(j) * x (i_column(j) )
         END DO
      END DO
      RETURN
      END FUNCTION Gmresm_Mtrx_times_vect
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Get_rotation ( vector_in, cos_theta, sin_theta )
      IMPLICIT NONE
      REAL , INTENT (IN)  :: vector_in (2)
      REAL , INTENT (OUT) :: cos_theta, sin_theta
!_____________________________________________________________________________
!     Compute a Givens (plane) rotation that will act on 2 elements of a
!     vector, A and B, and will zero B. Returns cosine and sine of THETA, the
!     angle of rotation. The transformation is
!
!        A_prime = A * COS(theta) - B * SIN(theta)
!        B_prime = A * CIN(theta) - B * COS(theta)
!
!     In matrix-vector terms,
!        X = (...... A ..... B .....)^T,
!        T =
!        X_prime = T * X = (..... A_prime ...... 0 .....)^T,
!        only 2 elements of X are changed by the rotation.
!_____________________________________________________________________________
!
      REAL :: temp
!
      IF ( ABS(vector_in (2)) < machine_tiny ) THEN
         cos_theta = 1.0
         sin_theta = 0.0
      ELSE IF ( ABS ( vector_in(2) ) > ABS ( vector_in (1) ) ) THEN
         temp = -vector_in (1) / vector_in (2)
         sin_theta = 1.0 / SQRT( 1.0 + temp**2 )
         cos_theta = temp * sin_theta
      ELSE
         temp = -vector_in (2) / vector_in (1)
         cos_theta = 1.0 / SQRT( 1.0 + temp**2 )
         sin_theta = temp * cos_theta
      END IF
!
      RETURN
      END SUBROUTINE Gmresm_Get_rotation
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Rotate_vector ( vec_in, cos_theta, sin_theta )
      IMPLICIT NONE
      REAL , INTENT (IN) :: vec_in (2), cos_theta, sin_theta
      REAL               :: Gmresm_Rotate_vector (2)
!___________________________________________________________________________
!     Apply a plane (Givens) rotation with cos_theta, sin_theta) to a vector.
!     Rotation acts on only two elements of the vector, X and Y.
!___________________________________________________________________________
!
      REAL :: temp
!
      temp                = cos_theta * vec_in (1) - sin_theta * vec_in (2)
      Gmresm_Rotate_vector(2)    = sin_theta * vec_in (1) + cos_theta * vec_in (2)
      Gmresm_Rotate_vector(1)    = temp
!
      RETURN
      END FUNCTION Gmresm_Rotate_vector
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Solve_upper_triang (a, b_rhs, n)
      IMPLICIT NONE
      INTEGER , INTENT (IN) :: n
      REAL , INTENT (IN) :: a (:,:), b_rhs (:)
      REAL  :: Gmresm_Solve_upper_triang (n)
!
!     Given an upper triangular matrix A and right-hand side vector B(n),
!     solves linear system A*x=b. A is a Hessenberg matrix (nmax+1 by nmax),
!     but only n by n section is used.
!
      INTEGER :: j
!
      IF (UBOUND(a,DIM=1) /= UBOUND (a,DIM=2) + 1) STOP 'PROBLEM 1 IN SOLVETR'
      IF (n > UBOUND (a,DIM = 2) .OR. n < 1) STOP 'PROBLEM 2 IN SOLVETR'
!
      Gmresm_Solve_upper_triang = b_rhs(1:n)
      DO j = N, 1, -1
         Gmresm_Solve_upper_triang (j)     = Gmresm_Solve_upper_triang (j) / a (j,j)
         Gmresm_Solve_upper_triang (1:j-1) = Gmresm_Solve_upper_triang (1:j-1) - &
                                      Gmresm_Solve_upper_triang (j) * a(1:j-1,j)
      END DO
      RETURN
      END FUNCTION Gmresm_Solve_upper_triang
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Msolve (x)
      IMPLICIT NONE
      REAL , INTENT (IN) :: x (:)
      REAL , DIMENSION (SIZE(x)) :: Gmresm_Msolve
!_____________________________________________________________________________
!     This subroutine solves the system L *  U * y = x, where
!     M = L * U = (D + L_a) * (I + D^(-1) * U_a) is the D-ILU preconditioner.
!     Matrices L_a and U_a are strictly lower and strictly upper triangular,
!     so that A = D_a + L_a + U_a, and D comes from incomplete LU factorization
!     when computing preconditioner. Solution proceeds in the regular way by
!     forward- and then back-substition (solving L*z=x, then U*y=z).
!     A is in the compressed-row-storage format (A_MTRX, I_COLUMN, ROW_PTR).
!     Diagonal matrix D (in fact, its inverse) is stored in the PIVOTS:
!     PIVOTS(i)=1/D(i,i), and DIAPTR vector holds locations of d_i_i in amatrx.
!     Since A_MTRX and PIVOTS do not change in the potential solver once is
!     has been called, we access them by host association from GMRESM. Only
!     vector X changes from invocation to invocaton of MSOLVE, and we pass it
!     as an argument.
!**** NOTE: book by Barrett et al ("templates ...") has an error in the back-
!           substitution algorithm (p.73). Here I do it correctly.
!_____________________________________________________________________________
!
      INTEGER :: i, j
      REAL    :: tmp
!
      DO i = 1, SIZE (x)
         tmp = 0.0
         DO j = row_ptr (i), diag_ptr (i) - 1
            tmp = tmp + a_mtrx (j) * Gmresm_msolve (i_column (j))
         END DO
         Gmresm_Msolve (i) = pivots (i) * (x(i) - tmp)
      END DO
      DO i = SIZE(x), 1, -1
         tmp = 0.0
         DO j = diag_ptr (i) + 1, row_ptr (i+1)-1
            tmp = tmp + a_mtrx(j) * Gmresm_Msolve (i_column(j))
         END DO
         Gmresm_Msolve (i) = Gmresm_Msolve (i) - pivots (i) * tmp
      END DO
      RETURN
      END FUNCTION Gmresm_Msolve
!
!
    END SUBROUTINE Gmresm
