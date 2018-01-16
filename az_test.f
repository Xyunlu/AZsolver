C***********************************************************************
C
       program main
C
C---------------------------------------------------------------
       use mymodule

       implicit none
       integer n, nrow
C..    orgin matrix
       integer :: nnz
       integer, dimension(:), allocatable :: na, ia
       double precision, dimension(:), allocatable :: a

c       integer, dimension(:), allocatable :: update, update_index
c       integer, dimension(:), allocatable :: external, extern_index

       integer, dimension(:), allocatable :: update
       integer, dimension(:), allocatable :: bindx
       double precision, dimension(:), allocatable :: val, b, x

       integer, dimension(:), allocatable :: nupdate

       integer myid, nproc
       integer N_update, nextern, ierror
       integer i

       include "az_aztecf.h"
       include "mpif.h"

       n = 3
       nrow = n*n

       allocate(na(nrow+1))
       na(1) = 1
       do i=1, nrow
         call create_mat_na(i, na, n)
       enddo
       nnz = na(nrow+1) - 1
       allocate(ia(nnz), a(nnz))
       do i=1, nrow
         call create_mat_val(i, n, na, ia, a)
       enddo

       call My_INIT(myid, nproc)

       N_update = nrow/nproc
       if( myid .lt. mod(nrow, nproc) ) N_update = N_update+1
c       print *,'myid, N_update=', myid, N_update

       allocate(nupdate(0:nproc))
       if( myid .eq. 0) then
          nupdate(myid+1) = N_update
          do i=1, nproc-1
            call My_recvai(0,i,nupdate(i+1),1)
          enddo
          nupdate(0) = 0
          do i=1, nproc
            nupdate(i) = nupdate(i)+nupdate(i-1)
          enddo
          do i=1, nproc-1
            call My_sendai(i,0,nupdate,nproc+1)
          enddo
       else
          call My_sendai(0,myid,N_update,1)
          call My_recvai(myid,0,nupdate,nproc+1)
       endif

       allocate(update(N_update))
       do i=1, N_update
         update(i) = nupdate(myid)+i
       enddo
       call crs2dmsr(myid,N_update,nextern,update,na,ia,a,
     +              bindx,val,b,x)
       print *,'myid, N_update, nextern=', myid,N_update,nextern
       print *,'update:', (update(i),i=1,N_update)

c       call azsolv(N_update, update, update_index,
c     +            external, extern_index,
c     +            bindx, val, b, x)
       call azsolv(N_update, nextern, update,
     +            bindx, val, b, x)

       call MPI_FINALIZE(ierror)
C
       end program

       subroutine create_mat_na(row, na, n)
       implicit none
       integer row, n, na(*)
       integer k

       k = 5
       if( mod((row-1)/n, n) .eq. 0 .or. mod((row-1)/n,n) .eq. n-1) k = k-1
       if( mod((row-1),n) .eq. 0 .or. mod((row-1),n) .eq.  n-1 ) k = k-1
       na(row+1) = na(row) + k

       return
       end

       subroutine create_mat_val(row, n, na, ia, a)
       implicit none
       integer row, n
       integer na(*), ia(*)
       double precision a(*)
       integer i, n0, n1

       n0 = na(row)
       n1 = na(row+1)-1
       if( mod((row-1)/n,n) .ne. 0) then
          a(n0) = -1.D0
          ia(n0) = row-n
          n0 = n0+1
       endif
       if( mod((row-1),n) .ne. 0) then
          a(n0) = -1.D0
          ia(n0) = row-1
          n0 = n0+1
       endif
       a(n0) = 4.D0
       ia(n0) = row
       n0 = n0+1
       if( mod((row-1),n) .ne. n-1) then
          a(n0) = -1.D0
          ia(n0) = row+1
          n0 = n0+1
       endif
       if( mod((row-1)/n,n) .ne. n-1) then
          a(n0) = -1.D0
          ia(n0) = row+n
          n0 = n0+1
       endif

       return
       end

