      subroutine pre_process(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer i

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),i)

        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT,i)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT,cid)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      character*(*) IMAGE_FILE,DIR_OUTPUT
      integer proc_error
      character*(strl) filename,catfile

      integer nx,ny
      real array(npx,npy),normap(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nxc
      double precision cRPIX(2),cD(2,2),cRVAL(2)

      real sigabc(2,3)
      integer order
      common /sig_pass/ sigabc,order

      real aa,bb,cc

      character*(strl) MASK_FILE
      integer nxx , nyy
      real flat_weight(npx,npy)

      integer cid
      proc_error=0
      do i=1,2
        do j=1,3
          sigabc(i,j)=0.
        enddo
      enddo

      call readimage_para(IMAGE_FILE
     .,nx,ny,npx,npy,array,cRPIX,cD,cRVAL)
      if (array(1,1).lt.(-99990.0)) then
        write(*,*) 'Error / Pre_proc cant find image file!'
        proc_error=1
      endif

      do i=1,nx
        do j=1,ny
          weight(i,j)=1
          if (array(i,j).gt.saturation_thresh) weight(i,j)=0
          normap(i,j)=array(i,j)
        enddo
      enddo

c ==========================================
c Function: Apply the DQ mask before astrometry and defect merging
c ==========================================
      if (proc_error.eq.0) then
        call fq_expo_ccd_product_path(IMAGE_FILE,DIR_OUTPUT,
     .  DIR_DQ,cid,
     .    '.fits',MASK_FILE)
        call readimage(MASK_FILE,nxx,nyy,npx,npy,flat_weight)
        if (flat_weight(1,1).lt.(-99990.0)) then
          write(*,*) 'Error / cant find mask file!'
          proc_error=1
        elseif (nxx.ne.nx .or. nyy.ne.ny) then
          write(*,*) 'Error / wrong size of DQ file!'
          proc_error=1
        endif
        if (proc_error.eq.0) then
          do i=1,nx
            do j=1,ny
              if (abs(flat_weight(i,j)).gt.1e-7) then
                weight(i,j)=0
              endif
            enddo
          enddo
        endif
      endif
c------------------------------------------------------

      nxc=nx/2

      if (ccD_split.eq.2) then

        call set_background(1,nxc,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_background(nxc+1,nx,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nxc,1,ny,npx,npy,normap
     .,aa,bb,cc,proc_error)
        if (proc_error.eq.0) then
          sigabc(1,1)=aa
          sigabc(1,2)=bb
          sigabc(1,3)=cc
        endif
        if (proc_error.eq.0) then
          call set_sig(nxc+1,nx,1,ny,npx,npy,normap
     .,aa,bb,cc,proc_error)
        endif
        if (proc_error.eq.0) then
          sigabc(2,1)=aa
          sigabc(2,2)=bb
          sigabc(2,3)=cc
        endif
      else
        call set_background(1,nx,1,ny,npx,npy,normap        
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nx,1,ny,npx,npy,normap
     .,aa,bb,cc,proc_error)
        if (proc_error.eq.0) then
          sigabc(1,1)=aa
          sigabc(1,2)=bb
          sigabc(1,3)=cc
        endif
      endif


c--------------------------------------------------------------
      call fq_chip_product_path(IMAGE_FILE,DIR_OUTPUT,
     .  DIR_ASTRO_DATA,
     .  '_astro.dat',filename)
      catfile=ASTROMETRY_cAT
      call generate_gaia_file_name(cRVAL,catfile,proc_error)
      call gen_astrometry_data(catfile,nx,ny,npx,npy
     .,normap,weight,cRPIX,cD,cRVAL,filename,proc_error)

      call locate_defects(nx,ny,npx,npy,array,normap
     .,weight,area_max,area_thresh,proc_error)
      call merge_defects(nx,ny,npx,npy,weight,normap
     .,area_max,source_thresh,area_thresh,proc_error)

c------------------------------------------------------------------

      if (proc_error.eq.0) then
c Serialize the final combined saturation, DQ, and detected-defect mask.
        do i=1,nx
          do j=1,ny
            if (weight(i,j).eq.0) normap(i,j)=-1000.
          enddo
        enddo
      else
        do i=1,nx
          do j=1,ny
            normap(i,j)=-1000.
          enddo
        enddo
      endif

      if (proc_error.eq.0) then
        normap(1,1)=-1.
      else
        normap(1,1)=1.
      endif

      do i=1,ccD_split
        do j=1,3
          normap(1+i,j)=sigabc(i,j)
        enddo
      enddo

      call fq_chip_product_path(IMAGE_FILE,DIR_OUTPUT,
     .  DIR_NORM,
     .  '_norm.fits',filename)
      call writeimage_copyhdu(IMAGE_FILE,filename
     .,nx,ny,npx,npy,normap)

c-----------------------------------------------------------
    !   if (proc_error.eq.0) then
    !     write(*,*) 'Status of processing ',trim(IMAGE_FILE)
    !  .,': OK.'
    !   else
    !     write(*,*) 'Status of processing ',trim(IMAGE_FILE)
    !  .,': ERROR!'
    !   endif
    !   if (proc_error.ne.0) 
    !  .write(*,*) 'Error / pre-processing ', trim(IMAGE_FILE)

c      pause

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_background(nx1,nx2,ny1,ny2,npx,npy,image
     .,blocksize,nct,ncx,ierror)
      implicit none

c ==========================================
c Function: Subtract the frozen historical Type-1 background
c Method: Use random middle-third rough fitting and rank-median blocks
c ==========================================
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,blocksize,nct,ncx
      real image(npx,npy)
      integer npp,nfit
      parameter (npp=1000)
      parameter (nfit=5000)
      real pix(npp),arr(npp,3),arr2(npp,3),c(nct),crough(4)
      real mean_sam(nfit,3),sam(nfit),tempo(nfit)
      integer indx(npp),nsam,nsam1,changed,nr
      integer ix,iy,i,j,k,nbx,nby,xmin,xmax,ymin,ymax,nct_min
      real mean,sig,aa,bb,cc,x,y,arr_min,arr_max,ratio
      real ran1,func_val

      if (ierror.eq.1) return
c ==========================================
c Function: Restore the original two-level F77 background workflow
c Method: Use random samples and absolute scaled pixel coordinates
c ==========================================
        ratio=1./max(nx2-nx1,ny2-ny1)

c Perform the historical rough flatten directly in absolute coordinates.
        do i=1,npp
          ix=int(ran1()*(nx2-nx1)+nx1)
          iy=int(ran1()*(ny2-ny1)+ny1)
          pix(i)=image(ix,iy)
          arr2(i,1)=ix*ratio
          arr2(i,2)=iy*ratio
          arr2(i,3)=pix(i)
        enddo
        call sort(npp,npp,pix)

        if (pix(1).eq.pix(npp)) then
          ierror=1
          return
        endif

        arr_min=pix(npp/3)
        arr_max=pix(2*npp/3)
        nr=0
        do i=1,npp
          if (arr2(i,3).ge.arr_min.and.arr2(i,3).le.arr_max) then
            nr=nr+1
            arr(nr,1)=arr2(i,1)
            arr(nr,2)=arr2(i,2)
            arr(nr,3)=arr2(i,3)
          endif
        enddo

        if (nr.lt.npp/10) then
          ierror=1
          return
        endif

        call fit_2D(npp,nr,arr,4,2,crough)

        do i=nx1,nx2
          do j=ny1,ny2
            x=i*ratio
            y=j*ratio
            image(i,j)=image(i,j)-func_val(x,y,4,2,crough)
          enddo
        enddo

c Fit and subtract the legacy fine background from the rough residual.
        nbx=max((nx2-nx1)/blocksize,1)
        nby=max((ny2-ny1)/blocksize,1)

        nsam=0
        do i=1,nbx
          xmin=(i-1)*blocksize+nx1
          xmax=min(xmin+blocksize,nx2)
          do j=1,nby
            ymin=(j-1)*blocksize+ny1
            ymax=min(ymin+blocksize,ny2)

            do k=1,npp
              ix=int(ran1()*(xmax-xmin)+xmin)
              iy=int(ran1()*(ymax-ymin)+ymin)
              arr(k,1)=ix*ratio
              arr(k,2)=iy*ratio
              arr(k,3)=image(ix,iy)
              pix(k)=arr(k,3)
            enddo

            call indexx(npp,npp,pix,indx)
            k=indx(npp/2)
            nsam=nsam+1
            mean_sam(nsam,1)=arr(k,1)
            mean_sam(nsam,2)=arr(k,2)
            mean_sam(nsam,3)=arr(k,3)
          enddo
        enddo

        nct_min=nct*3/2
        nsam1=nsam
        changed=1
        do while (changed.eq.1 .and. nsam1.ge.nct_min)
          call find_slope_2D(nfit,nsam1,mean_sam,aa,bb,cc)

          do i=1,nsam1
            sam(i)=mean_sam(i,3)-aa-bb*mean_sam(i,1)
     .             -cc*mean_sam(i,2)
            tempo(i)=sam(i)
          enddo
          call sort(nsam1,nfit,tempo)
          mean=tempo(nsam1/2)
          sig=0.5*(tempo(nsam1*5/6)-tempo(nsam1/6))
          nsam=0
          changed=0
          do i=1,nsam1
            if (abs(sam(i)-mean).lt.3.*sig) then
              nsam=nsam+1
              mean_sam(nsam,1)=mean_sam(i,1)
              mean_sam(nsam,2)=mean_sam(i,2)
              mean_sam(nsam,3)=mean_sam(i,3)
            else
              changed=1
            endif
          enddo
          nsam1=nsam
        enddo

        if (nsam1.lt.nct_min) then
          write(*,*) 'Background not stable enough!',nsam1
          ierror=1
          return
        endif

        call fit_2D(nfit,nsam1,mean_sam,nct,ncx,c)

        do i=nx1,nx2
          do j=ny1,ny2
            x=i*ratio
            y=j*ratio
            image(i,j)=image(i,j)-func_val(x,y,nct,ncx,c)
          enddo
        enddo
        return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_sig(nx1,nx2,ny1,ny2,npx,npy,image
     .,aa,bb,cc,ierror)
      implicit none

c ==========================================
c Function: Apply the frozen historical Type-1 sigma estimator
c Method: Fit the middle third of 2000 random squared differences
c ==========================================
      integer nx1,nx2,ny1,ny2,npx,npy,ierror
      integer npp_old
      parameter (npp_old=2000)
      integer i,j,ix_old,iy_old,nr_old
      real image(npx,npy),aa,bb,cc
      real pix_old(npp_old),arr_old(npp_old,3)
      real arr2_old(npp_old,3),arr_min_old,arr_max_old,temp_old
      real ran1

      aa=0.
      bb=0.
      cc=0.
      if (ierror.ne.0) return
c ==========================================
c Function: Restore the original random-2000 F77 sigma workflow
c Method: Trim the middle third and fit absolute pixel coordinates
c ==========================================
        do i=1,npp_old
          ix_old=int(ran1()*(nx2-nx1-1)+nx1)
          iy_old=int(ran1()*(ny2-ny1-1)+ny1)
          pix_old(i)=0.5*((image(ix_old,iy_old)
     .               -image(ix_old+1,iy_old))**2
     .               +(image(ix_old,iy_old)
     .               -image(ix_old,iy_old+1))**2)
          arr_old(i,1)=ix_old
          arr_old(i,2)=iy_old
          arr_old(i,3)=pix_old(i)
        enddo
        call sort(npp_old,npp_old,pix_old)

        arr_min_old=pix_old(npp_old/3)
        arr_max_old=pix_old(2*npp_old/3)
        nr_old=0
        do i=1,npp_old
          if (arr_old(i,3).ge.arr_min_old .and.
     .        arr_old(i,3).le.arr_max_old) then
            nr_old=nr_old+1
            arr2_old(nr_old,1)=arr_old(i,1)
            arr2_old(nr_old,2)=arr_old(i,2)
            arr2_old(nr_old,3)=arr_old(i,3)
          endif
        enddo

        if (nr_old.lt.npp_old/10) then
          ierror=1
          return
        endif

        call find_slope_2D(npp_old,nr_old,arr2_old,aa,bb,cc)

        do i=nx1,nx2
          do j=ny1,ny2
            temp_old=aa+bb*i+cc*j
            image(i,j)=image(i,j)/sqrt(0.5*temp_old)
          enddo
        enddo
        return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_defects(nx,ny,npx,npy
     .,weight,normap,area_max,source_thresh,area_thresh
     .,ierror)
      implicit none

      integer nx,ny,npx,npy,area_max,area_thresh,ierror
      real normap(npx,npy),source_thresh
      integer weight(npx,npy),mark(nx,ny)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k


      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).eq.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      if (nb.eq.area_max) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).gt.1
     . .or. (mark(u,v).eq.0 .and. weight(u,v).eq.0)) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=area_max
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine locate_defects(nx,ny,npx,npy,array,normap,weight
     .,area_max,area_thresh,ierror)
      implicit none

      integer nx,ny,npx,npy,ierror,area_max,area_thresh
      real normap(npx,npy),map(nx,ny),array(npx,npy)
      integer weight(npx,npy)
      real diffx(nx,ny),diffy(nx,ny)
      integer i,j,ix,iy
      character filename*100

      integer margin
      parameter (margin=10)
      real defect_halo_thresh
      parameter (defect_halo_thresh=1.)
      integer y_smooth,x_smooth
      parameter (y_smooth=200)
      parameter (x_smooth=100)

      real sig,med,sigx,medx,sigy,medy

      real loga,iden
      external loga,iden

      do j=1,ny
        do i=nx/2-margin,nx/2+margin
          weight(i,j)=0
        enddo
        do i=1,margin
          weight(i,j)=0
          weight(nx+1-i,j)=0
        enddo
      enddo
      do i=1,nx
        do j=1,margin
          weight(i,j)=0
          weight(i,ny+1-j)=0
        enddo
      enddo

      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
c          map(i,j)=loga(normap(i,j),1)
c          map(i,j)=log(normap(i,j)
          map(i,j)=loga(array(i,j),1)
        enddo
      enddo


      call remove_continuous(nx,ny,nx,ny,map,iden,4)
c      filename='map1.fits'
c      call writeimage(filename,nx,ny,nx,ny,map)
c      call get_sig_med(nx,ny,map,sig,med)
c      do i=1,nx
c        do j=1,ny
c          if (abs(map(i,j)).gt.5.*sig) weight(i,j)=0
c        enddo
c      enddo
      do i=1,nx
        ix=mod(i,nx)+1
        do j=1,ny
          iy=mod(j,ny)+1
          diffx(i,j)=map(i,j)-map(ix,j)
          diffy(i,j)=map(i,j)-map(i,iy)
        enddo
      enddo
      call get_sig_med(nx,ny,diffx,sigx,medx)
      do i=1,nx
        do j=1,ny
          if (abs(diffx(i,j)).gt.8.*sigx) weight(i,j)=0
        enddo
      enddo
      call get_sig_med(nx,ny,diffy,sigy,medy)
      do i=1,nx
        do j=1,ny
          if (abs(diffy(i,j)).gt.8.*sigy) weight(i,j)=0
        enddo
      enddo

      call mask_source_regions(nx,ny,npx,npy
     .,weight,normap,area_max,defect_halo_thresh*2.,area_thresh)

      call detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)

      call detect_artificial_stripes(nx,ny,npx,npy,weight
     .,diffx,diffy,sigx,sigy,medx,medy)

      do i=1,nx
        do j=1,ny
          if (weight(i,j).gt.1) weight(i,j)=1
        enddo
      enddo

      call detect_stellar_halo(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)

      call detect_dent(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)



      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_artificial_stripes(nx,ny,npx,npy,weight
     .,diffx,diffy,sigx,sigy,medx,medy)
      implicit none

      integer nx,ny,npx,npy
      integer weight(npx,npy)
      real diffx(nx,ny),diffy(nx,ny)
      real entropy(nx,ny),sig,med,sigx,sigy,medx,medy

      integer i,j

      call get_entropy(nx,ny,diffx,sigx,medx,2,entropy)
      call get_sig_med(nx,ny,entropy,sig,med)
c      call writeimage('entropy1.fits',nx,ny,nx,ny,entropy)
      do i=1,nx
        do j=1,ny
          if (weight(i,j).eq.1 .and. abs(entropy(i,j)-med).gt.10.*sig)
     . weight(i,j)=0
        enddo
      enddo

      call get_entropy(nx,ny,diffy,sigy,medy,2,entropy)
      call get_sig_med(nx,ny,entropy,sig,med)
c      call writeimage('entropy2.fits',nx,ny,nx,ny,entropy)
      do i=1,nx
        do j=1,ny
          if (weight(i,j).eq.1 .and. abs(entropy(i,j)-med).gt.10.*sig)
     . weight(i,j)=0
        enddo
      enddo

c      pause

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mask_source_regions(nx,ny,npx,npy
     .,weight,normap,area_max,source_thresh,area_thresh)
      implicit none

      integer nx,ny,npx,npy,area_max,area_thresh
      real normap(npx,npy),source_thresh
      integer weight(npx,npy),mark(nx,ny)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).eq.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      if (nb.eq.area_max) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=area_max
                weight(buffer(k,1),buffer(k,2))=2
              enddo
            else
              if (nb.ge.area_thresh) then
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                  weight(ix,iy)=2
                enddo
              else
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                enddo
              endif
            endif
          endif
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Detect stripe outliers including partial edge blocks
c Method: Use ceiling block counts and sum/sqrt(nvalid), so full,
c         partial, and DQ-masked blocks share one noise scale.
c ==========================================
      subroutine detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)
      implicit none

      integer nx,ny,npx,npy
      integer y_smooth,x_smooth
      real normap(npx,npy)
      real ymap(nx,(ny+y_smooth-1)/y_smooth)
      real xmap((nx+x_smooth-1)/x_smooth,ny)
      integer weight(npx,npy)
      integer yvalid(nx,(ny+y_smooth-1)/y_smooth)
      integer xvalid((nx+x_smooth-1)/x_smooth,ny)

      integer i,j,startj,endj,jj,numy,numx
      integer starti,endi,ii,nvalid
      real sig,med

      numy=(ny+y_smooth-1)/y_smooth
      do i=1,nx
        do j=1,numy
          ymap(i,j)=0.
          yvalid(i,j)=0
          nvalid=0
          startj=(j-1)*y_smooth+1
          endj=min(j*y_smooth,ny)
          do jj=startj,endj
            if (weight(i,jj).eq.1) then
              ymap(i,j)=ymap(i,j)+normap(i,jj)
              nvalid=nvalid+1
            endif
          enddo
          if (nvalid.gt.0) then
            ymap(i,j)=ymap(i,j)/sqrt(real(nvalid))
            yvalid(i,j)=1
          endif
        enddo
      enddo
      call get_stripe_sig_med(nx*numy,ymap,yvalid,sig,med)
      do i=1,nx
        do j=1,numy
          if (yvalid(i,j).eq.1 .and. sig.gt.0.
     .        .and. abs(ymap(i,j)-med).gt.sig*4.) then
            startj=(j-1)*y_smooth+1
            endj=min(j*y_smooth,ny)
            do jj=startj,endj
              weight(i,jj)=0
            enddo
          endif
        enddo
      enddo

      numx=(nx+x_smooth-1)/x_smooth
      do j=1,ny
        do i=1,numx
          xmap(i,j)=0.
          xvalid(i,j)=0
          nvalid=0
          starti=(i-1)*x_smooth+1
          endi=min(i*x_smooth,nx)
          do ii=starti,endi
            if (weight(ii,j).eq.1) then
              xmap(i,j)=xmap(i,j)+normap(ii,j)
              nvalid=nvalid+1
            endif
          enddo
          if (nvalid.gt.0) then
            xmap(i,j)=xmap(i,j)/sqrt(real(nvalid))
            xvalid(i,j)=1
          endif
        enddo
      enddo
      call get_stripe_sig_med(numx*ny,xmap,xvalid,sig,med)
      do j=1,ny
        do i=1,numx
          if (xvalid(i,j).eq.1 .and. sig.gt.0.
     .        .and. abs(xmap(i,j)-med).gt.sig*4.) then
            starti=(i-1)*x_smooth+1
            endi=min(i*x_smooth,nx)
            do ii=starti,endi
              weight(ii,j)=0
            enddo
          endif
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Estimate stripe median and width from valid blocks
c Method: Preserve the 1000-sample robust estimator while drawing only
c         from blocks that contain at least one usable pixel.
c ==========================================
      subroutine get_stripe_sig_med(np,map,valid,sig,med)
      implicit none

      integer np,valid(np)
      real map(np),sig,med

      integer npp
      parameter (npp=1000)
      integer i,k,nuse
      real values(np),pix(npp),ran1

      nuse=0
      do i=1,np
        if (valid(i).eq.1) then
          nuse=nuse+1
          values(nuse)=map(i)
        endif
      enddo

      if (nuse.eq.0) then
        sig=-1.
        med=0.
        return
      endif

      do k=1,npp
        i=int(ran1()*real(nuse))+1
        pix(k)=values(i)
      enddo

      call sort(npp,npp,pix)
      sig=0.5*(pix(5*npp/6)-pix(npp/6))
      med=pix(npp/2)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_stellar_halo(nx,ny,npx,npy,normap,weight
     .,npmax,defect_halo_thresh)
      implicit none

c The purpose of this code is to mark the parts affected by stellar halos.

c input & output:
      integer nx,ny,npx,npy,npmax
      real normap(npx,npy),defect_halo_thresh
      integer weight(npx,npy)
      real smoothed(nx,ny)

c local variables:
      integer dmark(nx,ny),buffer(npmax,2)
      integer nb,nbb,toobig,i,j,ix,iy,jx,jy,u,v,k1,k2,k


      do i=1,nx
        do j=1,ny
          smoothed(i,j)=normap(i,j)
        enddo
      enddo

      call smooth_image55(nx,ny,smoothed,1)


      do i=1,nx
        do j=1,ny
          if (smoothed(i,j).ge.defect_halo_thresh) then
            dmark(i,j)=1
          else
            dmark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (dmark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            dmark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (dmark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      dmark(u,v)=-1
                      if (nb.eq.npmax) then
                        toobig=1
                        goto 20
                      endif
                    elseif (dmark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=npmax
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_dent(nx,ny,npx,npy,normap,weight
     .,npmax,defect_halo_thresh)
      implicit none

c input & output:
      integer nx,ny,npx,npy,npmax
      real normap(npx,npy),defect_halo_thresh
      integer weight(npx,npy)

c local variables:
      integer dmark(nx,ny),buffer(npmax,2)
      integer nb,nbb,toobig,i,j,ix,iy,jx,jy,u,v,k1,k2,k


      do i=1,nx
        do j=1,ny
          if (normap(i,j).le.-defect_halo_thresh) then
            dmark(i,j)=1
          else
            dmark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (dmark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            dmark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (dmark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      dmark(u,v)=-1
                      if (nb.eq.npmax) then
                        toobig=1
                        goto 20
                      endif
                    elseif (dmark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=npmax
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
