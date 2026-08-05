      subroutine Gen_cat_loop(Dst_DIR)
        use mpi
        implicit none
        include 'para.inc'

      integer ierr, my_id, num_procs
      common /MPIpar/my_id, num_procs

      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      real cat(nmax_per_core,ichi2)
      integer ng
      common /shear_data_pass/ cat,ng

      integer nchip
      character*(strl) Dst_DIR,DIR_OUTPUT,filename
      character*(strl) PREFIX,IMAGE_FILE(NMAX_cHIP)
      character*400 cat_list
      integer ierror, i, u
      integer subnum

      ! Grid variables
      integer NTILES_MAX
      parameter (NTILES_MAX = 6480000)
      integer, allocatable :: tile_counts(:), global_tile_counts(:)
      integer, allocatable :: tile_partition(:)
      integer N_active, i_dec, i_ra, tile_index, i_dec_shifted
      real dra, ddec

      ! Active bin arrays
      real, allocatable :: active_dec(:), active_ra(:)
      integer, allocatable :: active_weight(:), active_idx(:)
      integer, allocatable :: active_part_id(:)

      ! Bisection variables
      integer N_total, K, p_id, target_rank
      integer, allocatable :: target_proc(:), part_id(:)

      ! Redistribution variables
      integer, allocatable :: send_counts(:), send_displs(:)
      integer, allocatable :: recv_counts_all(:), recv_displs(:)
      real, allocatable :: send_buf(:,:), recv_buf(:,:)
      integer, allocatable :: send_part_buf(:), recv_part_buf(:)
      integer, allocatable :: send_counts_real(:), recv_counts_real(:)
      integer, allocatable :: send_displs_real(:), recv_displs_real(:)
      integer, allocatable :: curr_disp(:)
      integer total_recv, local_idx, p
      integer status(mpi_status_size)

      ! Output variables
      real, allocatable :: sort_buf(:,:)
      character*6 c_part
      real p_dec_min, p_dec_max, p_ra_min, p_ra_max
      logical, save :: initialized = .false.

      ! Parallel summary variables
      integer max_parts_pp, j_part
      integer, allocatable :: sum_pid(:), sum_count(:)
      real, allocatable :: sum_decmin(:), sum_decmax(:)
      real, allocatable :: sum_ramin(:), sum_ramax(:)
      integer, allocatable :: all_pid(:), all_count(:)
      real, allocatable :: all_decmin(:), all_decmax(:)
      real, allocatable :: all_ramin(:), all_ramax(:)

      external read_shear_dat
 500  format(I6, 2X, I8, 2X, F9.4, 2X, F9.4, 2X, F9.4, 2X, F9.4)

      ! Step 1: Parallel exposure reading
      call mpi_distribute(N_EXPO,read_shear_dat,'read shear.dat')
      write(*,*) 'Rank', my_id, ' read galaxies:', ng
      call flush(6)

      ! Step 2: Local binning on 0.1 degree grid
      allocate(tile_counts(NTILES_MAX))
      allocate(global_tile_counts(NTILES_MAX))
      tile_counts = 0

      do i = 1, ng
         dra = cat(i, icra)
         if (dra .ge. 360.0) dra = dra - 360.0
         if (dra .lt. 0.0) dra = dra + 360.0

         ddec = cat(i, icdec)
         if (ddec .ge. 90.0) ddec = 89.999
         if (ddec .lt. -90.0) ddec = -90.0

         i_dec = floor(ddec * 10.0)
         i_ra = floor(dra * 10.0)
         if (i_dec .lt. -900) i_dec = -900
         if (i_dec .gt. 899) i_dec = 899
         if (i_ra .lt. 0) i_ra = 0
         if (i_ra .gt. 3599) i_ra = 3599
         tile_index = (i_dec + 900) * 3600 + i_ra + 1
         tile_counts(tile_index) = tile_counts(tile_index) + 1
      enddo

      if (my_id .eq. 0) then
         write(*,*) '[Step 2] Local binning done.'
         call flush(6)
      endif

      ! Step 3: Global sum of grid counts
      call MPI_Allreduce(tile_counts, global_tile_counts, NTILES_MAX,
     .                   MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      deallocate(tile_counts)

      if (my_id .eq. 0) then
         write(*,*) '[Step 3] MPI_Allreduce done.'
         call flush(6)
      endif

      ! Step 4: Extract active bins
      N_active = 0
      do i = 1, NTILES_MAX
         if (global_tile_counts(i) .gt. 0) then
            N_active = N_active + 1
         endif
      enddo

      allocate(active_dec(max(N_active,1)))
      allocate(active_ra(max(N_active,1)))
      allocate(active_weight(max(N_active,1)))
      allocate(active_idx(max(N_active,1)))
      allocate(active_part_id(max(N_active,1)))

      N_active = 0
      do i = 1, NTILES_MAX
         if (global_tile_counts(i) .gt. 0) then
            N_active = N_active + 1
            active_idx(N_active) = N_active
            active_weight(N_active) = global_tile_counts(i)
            
            i_dec_shifted = (i - 1) / 3600
            i_ra = mod(i - 1, 3600)
            active_dec(N_active) = (real(i_dec_shifted) - 900.0) / 10.0
     .                             + 0.05
            active_ra(N_active) = real(i_ra) / 10.0 + 0.05
         endif
      enddo

      ! Step 5: Deterministic parallel K-D tree partition
      N_total = 0
      do i = 1, N_active
         N_total = N_total + active_weight(i)
      enddo

      K = (N_total + n_target_subcat - 1) / n_target_subcat
      if (K .lt. 1) K = 1

      if (my_id .eq. 0) then
         write(*,*) '[Step 5] N_total=', N_total, ' N_active=',
     .              N_active, ' K=', K
         call flush(6)
      endif

      call kdtree_partition(active_dec, active_ra, active_weight,
     .                      active_idx, 1, N_active, K, 0, N_active,
     .                      1, active_part_id)

      if (my_id .eq. 0) then
         write(*,*) '[Step 5] K-D tree partition done.'
         call flush(6)
      endif

      ! Step 6: Build global lookup mapping
      allocate(tile_partition(NTILES_MAX))
      tile_partition = 0

      N_active = 0
      do i = 1, NTILES_MAX
         if (global_tile_counts(i) .gt. 0) then
            N_active = N_active + 1
            tile_partition(i) = active_part_id(N_active)
         endif
      enddo

      deallocate(global_tile_counts)
      deallocate(active_dec, active_ra, active_weight, active_idx)
      deallocate(active_part_id)

      ! Step 7: Map local galaxies to target process and partition
      allocate(target_proc(max(ng,1)))
      allocate(part_id(max(ng,1)))

      do i = 1, ng
         dra = cat(i, icra)
         if (dra .ge. 360.0) dra = dra - 360.0
         if (dra .lt. 0.0) dra = dra + 360.0

         ddec = cat(i, icdec)
         if (ddec .ge. 90.0) ddec = 89.999
         if (ddec .lt. -90.0) ddec = -90.0

         i_dec = floor(ddec * 10.0)
         i_ra = floor(dra * 10.0)
         if (i_dec .lt. -900) i_dec = -900
         if (i_dec .gt. 899) i_dec = 899
         if (i_ra .lt. 0) i_ra = 0
         if (i_ra .gt. 3599) i_ra = 3599
         tile_index = (i_dec + 900) * 3600 + i_ra + 1

         p_id = tile_partition(tile_index)
         part_id(i) = p_id
         target_proc(i) = mod(p_id - 1, num_procs)
      enddo

      deallocate(tile_partition)

      if (my_id .eq. 0) then
         write(*,*) '[Step 7] Galaxy mapping done.'
         call flush(6)
      endif

      ! Step 8: Alltoallv redistribution of galaxy data
      allocate(send_counts(0:num_procs-1))
      allocate(send_displs(0:num_procs-1))
      allocate(recv_counts_all(0:num_procs-1))
      allocate(recv_displs(0:num_procs-1))

      send_counts = 0
      do i = 1, ng
         p = target_proc(i)
         send_counts(p) = send_counts(p) + 1
      enddo

      call MPI_Alltoall(send_counts, 1, MPI_INTEGER,
     .                  recv_counts_all, 1, MPI_INTEGER,
     .                  MPI_COMM_WORLD, ierr)

      send_displs(0) = 0
      recv_displs(0) = 0
      total_recv = recv_counts_all(0)
      do i = 1, num_procs-1
         send_displs(i) = send_displs(i-1) + send_counts(i-1)
         recv_displs(i) = recv_displs(i-1) + recv_counts_all(i-1)
         total_recv = total_recv + recv_counts_all(i)
      enddo

      if (my_id .eq. 0) then
         write(*,*) '[Step 8] Rank 0 total_recv=', total_recv
         call flush(6)
      endif

      allocate(send_buf(ichi2, max(ng,1)))
      allocate(send_part_buf(max(ng,1)))
      allocate(recv_buf(ichi2, max(total_recv,1)))
      allocate(recv_part_buf(max(total_recv,1)))

      allocate(curr_disp(0:num_procs-1))
      curr_disp = send_displs
      do i = 1, ng
         p = target_proc(i)
         local_idx = curr_disp(p) + 1
         do u = 1, ichi2
            send_buf(u, local_idx) = cat(i, u)
         enddo
         send_part_buf(local_idx) = part_id(i)
         curr_disp(p) = curr_disp(p) + 1
      enddo
      deallocate(curr_disp, target_proc, part_id)

      allocate(send_counts_real(0:num_procs-1))
      allocate(recv_counts_real(0:num_procs-1))
      allocate(send_displs_real(0:num_procs-1))
      allocate(recv_displs_real(0:num_procs-1))

      send_counts_real = send_counts * ichi2
      recv_counts_real = recv_counts_all * ichi2
      send_displs_real = send_displs * ichi2
      recv_displs_real = recv_displs * ichi2

      call MPI_Alltoallv(send_buf, send_counts_real,
     .   send_displs_real, MPI_REAL, recv_buf, recv_counts_real,
     .   recv_displs_real, MPI_REAL, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
         write(*,*) '[Step 8] Alltoallv for galaxy data done.'
         call flush(6)
      endif

      call MPI_Alltoallv(send_part_buf, send_counts,
     .   send_displs, MPI_INTEGER, recv_part_buf, recv_counts_all,
     .   recv_displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      deallocate(send_counts, send_displs, recv_counts_all)
      deallocate(recv_displs, send_counts_real, recv_counts_real)
      deallocate(send_displs_real, recv_displs_real, send_buf)
      deallocate(send_part_buf)

      if (my_id .eq. 0) then
         write(*,*) '[Step 8] Alltoallv complete.'
         call flush(6)
      endif

      ! Step 9: Get headers
      if (.not. initialized) then
        if (my_id .eq. 0) then
          call get_image_list(1,IMAGE_FILE,nchip,DIR_OUTPUT)
          call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
          filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//
     .             '_all.cat'
          open(unit=11,file=filename,status='old',action='read'
     .                                              ,iostat=ierror)
          if (ierror.eq.0) then
              read(11,'(A)') cat_list
              close(11)
          else
            write(*,*) 'Rank 0 cannot open header file:',
     .                 trim(filename)
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr) 
          endif
        endif
        call MPI_Bcast(cat_list, 400, MPI_CHARACTER, 0,
     .                 MPI_COMM_WORLD, ierr)
        initialized = .true.
      endif

      ! Step 10: Write sub-catalogs (fully parallel, no barriers)
      max_parts_pp = (K + num_procs - 1) / num_procs

      allocate(sum_pid(max_parts_pp))
      allocate(sum_count(max_parts_pp))
      allocate(sum_decmin(max_parts_pp))
      allocate(sum_decmax(max_parts_pp))
      allocate(sum_ramin(max_parts_pp))
      allocate(sum_ramax(max_parts_pp))
      sum_pid = 0
      sum_count = 0

      if (my_id .eq. 0) then
         write(*,*) '[Step 10] Writing sub-catalogs in parallel...'
         call flush(6)
      endif

      j_part = 0
      do p_id = my_id + 1, K, num_procs
         j_part = j_part + 1
         sum_pid(j_part) = p_id

         subnum = 0
         do i = 1, total_recv
            if (recv_part_buf(i) .eq. p_id) subnum = subnum + 1
         enddo
         sum_count(j_part) = subnum
         if (subnum .le. 0) cycle

         allocate(sort_buf(ichi2, subnum))
         local_idx = 0
         do i = 1, total_recv
            if (recv_part_buf(i) .eq. p_id) then
               local_idx = local_idx + 1
               do u = 1, ichi2
                  sort_buf(u, local_idx) = recv_buf(u, i)
               enddo
            endif
         enddo

         call reording_cat_multi(sort_buf, ichi2, subnum,
     .                           icdec, icra)

         write(c_part, '(I6.6)') p_id
         filename = trim(Dst_DIR)//'/subcat_'//
     .      trim(adjustl(c_part))//'.cat'

         open(unit=20, file=trim(filename), status='replace',
     .        iostat=ierror)
         write(20,'(A)') trim(cat_list)
         do i = 1, subnum
            write(20,*) (sort_buf(u, i), u=1, ichi2)
         enddo
         close(20)

         write(*,*) 'Rank', my_id, ' wrote ', subnum,
     .      ' galaxies to ', trim(filename)
         call flush(6)

         sum_decmin(j_part) = sort_buf(icdec, 1)
         sum_decmax(j_part) = sort_buf(icdec, 1)
         sum_ramin(j_part) = sort_buf(icra, 1)
         sum_ramax(j_part) = sort_buf(icra, 1)
         do i = 2, subnum
            if (sort_buf(icdec, i) .lt. sum_decmin(j_part))
     .         sum_decmin(j_part) = sort_buf(icdec, i)
            if (sort_buf(icdec, i) .gt. sum_decmax(j_part))
     .         sum_decmax(j_part) = sort_buf(icdec, i)
            if (sort_buf(icra, i) .lt. sum_ramin(j_part))
     .         sum_ramin(j_part) = sort_buf(icra, i)
            if (sort_buf(icra, i) .gt. sum_ramax(j_part))
     .         sum_ramax(j_part) = sort_buf(icra, i)
         enddo
         deallocate(sort_buf)
      enddo

      deallocate(recv_buf, recv_part_buf)

c     Gather summary data to Rank 0 via MPI_Gather
      if (my_id .eq. 0) then
         allocate(all_pid(max_parts_pp * num_procs))
         allocate(all_count(max_parts_pp * num_procs))
         allocate(all_decmin(max_parts_pp * num_procs))
         allocate(all_decmax(max_parts_pp * num_procs))
         allocate(all_ramin(max_parts_pp * num_procs))
         allocate(all_ramax(max_parts_pp * num_procs))
      else
         allocate(all_pid(1), all_count(1))
         allocate(all_decmin(1), all_decmax(1))
         allocate(all_ramin(1), all_ramax(1))
      endif

      call MPI_Gather(sum_pid, max_parts_pp, MPI_INTEGER,
     .               all_pid, max_parts_pp, MPI_INTEGER,
     .               0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(sum_count, max_parts_pp, MPI_INTEGER,
     .               all_count, max_parts_pp, MPI_INTEGER,
     .               0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(sum_decmin, max_parts_pp, MPI_REAL,
     .               all_decmin, max_parts_pp, MPI_REAL,
     .               0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(sum_decmax, max_parts_pp, MPI_REAL,
     .               all_decmax, max_parts_pp, MPI_REAL,
     .               0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(sum_ramin, max_parts_pp, MPI_REAL,
     .               all_ramin, max_parts_pp, MPI_REAL,
     .               0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(sum_ramax, max_parts_pp, MPI_REAL,
     .               all_ramax, max_parts_pp, MPI_REAL,
     .               0, MPI_COMM_WORLD, ierr)

c     Rank 0 writes sorted summary file
      if (my_id .eq. 0) then
         filename = './catalog_summary.txt'
         open(unit=21, file=trim(filename), status='replace',
     .        iostat=ierror)
         write(21, '(A)') 'part_id  count  dec_min  dec_max  '//
     .                    'ra_min   ra_max'
         do p_id = 1, K
            target_rank = mod(p_id - 1, num_procs)
            j_part = (p_id - 1) / num_procs + 1
            i = target_rank * max_parts_pp + j_part
            if (all_count(i) .gt. 0) then
               write(21, 500) all_pid(i), all_count(i),
     .                        all_decmin(i), all_decmax(i),
     .                        all_ramin(i), all_ramax(i)
            endif
         enddo
         close(21)
         write(*,*) '[Step 10] All sub-catalogs written. K=', K
         call flush(6)
      endif

      deallocate(sum_pid, sum_count, sum_decmin, sum_decmax)
      deallocate(sum_ramin, sum_ramax)
      deallocate(all_pid, all_count, all_decmin, all_decmax)
      deallocate(all_ramin, all_ramax)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_shear_dat(iexpo)
        implicit none
        include 'para.inc'

      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO
      
      integer iexpo,nchip,ierror,ierr
      integer u,undefined_flag
      integer,save :: ntot
      data ntot /0/

      character*(strl) filename,DIR_OUTPUT
      character*(strl) PREFIX,IMAGE_FILE(NMAX_cHIP)
      real item(ichi2)
      
      real cat(nmax_per_core,ichi2)
      integer ng
      common /shear_data_pass/ cat,ng

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)
      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)

      filename=trim(DIR_OUTPUT)//'/result/'
     .//trim(PREFIX)//'_all.cat'
      open(unit=10,file=filename,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
          write(*,*) filename, ' does not exist!!'
          return
      endif
      ! read the first line
      read(10,*)
      do while (ierror.ge.0)
          read(10,*,iostat=ierror) (item(u),u=1,ichi2)
          if (ierror.lt.0) cycle
          undefined_flag=0
          do u=1,ichi2
            if (item(u).ne.item(u)) then
                ! NaN check（NaN ont equal to itself）
                write(*,*) 'Error / NaN detected!',filename
                undefined_flag = 1
                exit
            endif
            if (abs(item(u)).gt.1.0d300) then
                ! check Inf
                write(*,*) 'Error / Inf detected!',filename
                undefined_flag = 1
                exit
            endif
          enddo         
          if (undefined_flag.eq.1) cycle
          ntot=ntot+1
          if (ntot.gt.nmax_per_core) then 
            write(*,*) 'Error / nmax_per_core is too small!!!'
            stop
          endif
          do u=1,ichi2
            cat(ntot,u)=item(u)
          enddo
      enddo
      close(10)
      ng=ntot

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
