      subroutine mpi_distribute(num_job,job,message)
      use mpi
      implicit none

      integer my_id,num_procs
      common /MPIpar/ my_id,num_procs

      external job

      character*(*) message

      integer num_job,i,j,k,complete
      integer source,tag,ierr,status(mpi_status_size)


      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes

      i=0
      if (my_id.eq.0) then
        i=1
        j=num_job
      endif
      complete=0
      do while (complete.eq.0)
        if (my_id.ne.0) then
          call MPI_SEND(i,1,mpi_int,0,0,MPI_cOMM_WORLD,ierr)
          call MPI_REcV(i,1,mpi_int,0,0,MPI_cOMM_WORLD,status,ierr)
          if (i.eq.0) then
            complete=1
          else
            call job(i)
          endif
        else
          call MPI_REcV(k,1,mpi_int,MPI_ANY_SOURcE
     .,MPI_ANY_TAG,MPI_cOMM_WORLD,status,ierr)
          tag=status(MPI_TAG)
          source=status(MPI_SOURcE)
          call MPI_SEND(i,1,mpi_int,source,tag,MPI_cOMM_WORLD,ierr)
          if (k.gt.0) then
            j=j-1
          endif
          write(*,*) source,i,j,trim(message)
          if (i.ne.0) i=i+1
          if (i.gt.num_job) i=0
          if (j.eq.0) complete=1
        endif
      enddo

      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer iexpo,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) EXPO_FILE(NMAX_EXPO),filename
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO
      integer ierr

      nchip=0
      open(unit=10,file=EXPO_FILE(iexpo),status='old',iostat=ierr)
      rewind 10
      if (ierr.ne.0) then
        write(*,*) 'EXPO_FILE reading error!!',EXPO_FILE(iexpo)
        stop
      endif
      do while (ierr.ge.0)
        read(10,'(A)',iostat=ierr) filename
        if (ierr.lt.0) cycle
        nchip=nchip+1
        IMAGE_FILE(nchip)=filename
      enddo
      close(10)

      call get_dir(IMAGE_FILE(1),DIR_OUTPUT,2)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_dir(imagefile,dir,level)
      implicit none

      character*(*) imagefile,dir
      integer i,n,level,u

      u=0
      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'/') u=u+1
        if (u.eq.level) exit
      enddo
      if (u.ne.level) stop 'Image_file name is NOT normal !'

      dir=imagefile(1:i-1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PREFIX(imagefile,PREFIX)
      implicit none

      character*(*) imagefile,PREFIX
      integer i,p_dot,p_slash,n

      p_dot=0
      p_slash=0

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'.'.and.p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/'.and.p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0.or.p_slash.eq.0..or.p_dot.le.p_slash+1) then
        write(*,*) 'Image_file name is NOT normal !'
        read(*,*)
      endif

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PREFIX_expo(imagefile,PREFIX)
      implicit none

      character*(*) imagefile,PREFIX
      integer i,p_dot,p_slash,n

      p_dot=0
      p_slash=0

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'_'.and.p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/'.and.p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0.or.p_slash.eq.0..or.p_dot.le.p_slash+1) 
     .stop 'Image_file name is NOT normal !'

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reording_cat_multi(buf, nfeat, ngal, idec, ira)
      implicit none
      integer nfeat, ngal, idec, ira
      real buf(nfeat, ngal)

      integer :: i
      integer, allocatable :: idx(:)
      real,    allocatable :: tmp(:,:)

      allocate(idx(ngal))
      allocate(tmp(nfeat,ngal))

      ! initial idx = [1,2,...,ngal]
      do i = 1, ngal
          idx(i) = i
      enddo

      ! multi-key quicksort： dec, ra
      call quicksort_multi(buf, idx, idec, ira, 1, ngal, nfeat, ngal)

      ! reordering
      do i = 1, ngal
          tmp(:, i) = buf(:, idx(i))
      enddo

      buf(:,:) = tmp(:,:)

      deallocate(idx)
      deallocate(tmp)
      end
C===========================================================
C Quick Sort 
C===========================================================
      RECURSIVE subroutine quicksort_multi
     .(buf, idx, idec, ira, left, right, nfeat, ngal)
      implicit none
      integer nfeat, ngal
      real buf(nfeat, ngal)
      integer idx(ngal), left, right, idec, ira
      integer i, j, tmpi, piv

      if (left >= right) return

      ! pivot chosen as middle element
      piv = idx((left + right) / 2)

      i = left
      j = right

      do
          ! right forward buf(idec, idx(i)) < buf(idec, piv)
          do while ( buf(idec,idx(i)) < buf(idec,piv) .or.
     &            (buf(idec,idx(i)) == buf(idec,piv) .and.
     &             buf(ira, idx(i)) < buf(ira, piv) ) )
              i = i + 1
          enddo

          ! left forward buf(idec, idx(j)) > buf(idec, piv)
          do while ( buf(idec,idx(j)) > buf(idec,piv) .or.
     &            (buf(idec,idx(j)) == buf(idec,piv) .and.
     &             buf(ira, idx(j)) > buf(ira, piv) ) )
              j = j - 1
          enddo

          if (i <= j) then
              tmpi = idx(i)
              idx(i) = idx(j)
              idx(j) = tmpi
              i = i + 1
              j = j - 1
          endif

          if (i > j) exit
      enddo

      call quicksort_multi(buf, idx, idec, ira, left, j, nfeat, ngal)
      call quicksort_multi(buf, idx, idec, ira, i, right, nfeat, ngal)
      end
C===========================================================
C 2D Grid K-D Tree Bisection Subroutines
C===========================================================
      recursive subroutine quicksort_bins(dec, ra, idx, left, right,
     .                                    dim, ntot)
      implicit none
      integer ntot, left, right, dim
      real dec(ntot), ra(ntot)
      integer idx(ntot)
      integer i, j, tmpi, piv
      real piv_dec, piv_ra

      if (left .ge. right) return

      piv = idx((left + right) / 2)
      piv_dec = dec(piv)
      piv_ra = ra(piv)

      i = left
      j = right

      do
         if (dim .eq. 0) then
            do while (dec(idx(i)) .lt. piv_dec .or.
     .               (dec(idx(i)) .eq. piv_dec .and.
     .                ra(idx(i)) .lt. piv_ra))
               i = i + 1
            enddo
            do while (dec(idx(j)) .gt. piv_dec .or.
     .               (dec(idx(j)) .eq. piv_dec .and.
     .                ra(idx(j)) .gt. piv_ra))
               j = j - 1
            enddo
         else
            do while (ra(idx(i)) .lt. piv_ra .or.
     .               (ra(idx(i)) .eq. piv_ra .and.
     .                dec(idx(i)) .lt. piv_dec))
               i = i + 1
            enddo
            do while (ra(idx(j)) .gt. piv_ra .or.
     .               (ra(idx(j)) .eq. piv_ra .and.
     .                dec(idx(j)) .gt. piv_dec))
               j = j - 1
            enddo
         endif

         if (i .le. j) then
            tmpi = idx(i)
            idx(i) = idx(j)
            idx(j) = tmpi
            i = i + 1
            j = j - 1
         endif

         if (i .gt. j) exit
      enddo

      call quicksort_bins(dec, ra, idx, left, j, dim, ntot)
      call quicksort_bins(dec, ra, idx, i, right, dim, ntot)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      recursive subroutine kdtree_partition(dec, ra, weight, idx,
     .                                      left, right, k, dim,
     .                                      ntot, p_start,
     .                                      part_id_all)
      implicit none
      integer ntot, left, right, k, dim, p_start
      real dec(ntot), ra(ntot)
      integer weight(ntot), idx(ntot)
      integer part_id_all(ntot)

      integer k1, k2, i, mid
      double precision W_total, W_target, W_run, W_prev

      if (k .eq. 1) then
         do i = left, right
            part_id_all(idx(i)) = p_start
         enddo
         return
      endif

      if (left .ge. right) then
         do i = left, right
            part_id_all(idx(i)) = p_start
         enddo
         return
      endif

      call quicksort_bins(dec, ra, idx, left, right, dim, ntot)

      W_total = 0.0d0
      do i = left, right
         W_total = W_total + dble(weight(idx(i)))
      enddo

      k1 = k / 2
      k2 = k - k1
      W_target = W_total * dble(k1) / dble(k)

      W_run = 0.0d0
      mid = left
      do i = left, right - 1
         W_prev = W_run
         W_run = W_run + dble(weight(idx(i)))
         if (W_run .ge. W_target) then
            if (i .gt. left) then
               if (abs(W_prev - W_target) .lt. abs(W_run - W_target))
     .         then
                  mid = i - 1
               else
                  mid = i
               endif
            else
               mid = left
            endif
            exit
         endif
         mid = i
      enddo

      call kdtree_partition(dec, ra, weight, idx, left, mid, k1,
     .                      1-dim, ntot, p_start, part_id_all)
      call kdtree_partition(dec, ra, weight, idx, mid+1, right,
     .                      k2, 1-dim, ntot, p_start+k1,
     .                      part_id_all)
      end
