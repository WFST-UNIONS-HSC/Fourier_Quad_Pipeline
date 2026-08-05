      program main
      use mpi
      include 'para.inc'

      integer ierr, my_id, num_procs
      common /MPIpar/my_id, num_procs

      character*(strl) EXPO_LIST,DST_DIR
      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      !----- MPI initialization ---------
      call MPI_Init( ierr )
      call MPI_comm_rank(MPI_cOMM_WORLD, my_id, ierr )
      call MPI_comm_size(MPI_cOMM_WORLD, num_procs, ierr )
      ! my_id = 0,1,...,num_procs-1
      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes

      call getarg(1,EXPO_LIST)
      call getarg(2,DST_DIR)
      ! -----------------------------------------------------------------
      if (my_id.eq.0) call initialize(EXPO_LIST)
      call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_cOMM_WORLD,ierr)
      call MPI_Bcast(EXPO_FILE,NMAX_EXPO*strl,mpi_character,0
     .,MPI_cOMM_WORLD,ierr)

      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      ! -----------------------------------------------------------------
      call Gen_cat_loop(DST_DIR)

      if (my_id.eq.0) then
        write(*,*) 'All done!'
      endif


      !----- MPI finalization -------------------------------------------
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      call MPI_Finalize (ierr)

      stop
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initialize(EXPO_LIST)
      implicit none
      include 'para.inc'

      character*(strl) EXPO_LIST
      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      character*(strl) expo_name
      integer ierror,nchip

      N_EXPO=0

      open(unit=10,file=EXPO_LIST,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) 'EXPO_LIST reading error!!'
        stop
      endif
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) expo_name , nchip
        if (ierror.lt.0) cycle
        N_EXPO=N_EXPO+1
        EXPO_FILE(N_EXPO)=trim(expo_name)
      enddo
      close(10)
      write(*,*) 'Total number of EXPOSURE: ',N_EXPO


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
