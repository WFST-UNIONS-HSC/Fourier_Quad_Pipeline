      subroutine pre_process(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT,FLAT_FILE
      integer i
      character*2 chip_name

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

20      format(I2.2)

      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),i)
        write(chip_name,20) i
        FLAT_FILE=trim(FLAT_PATH)//'/flat_'
     .//trim(chip_name)//'_weight.fits'

        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT
     .,FLAT_FILE)
c        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT,FLAT_FILE)
c      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      character*(*) IMAGE_FILE,DIR_OUTPUT
      integer proc_error
      character*(strl) PREFIX,filename,catfile

      integer nx,ny
      real array(npx,npy),normap(npx,npy)
      integer weight(npx,npy)
c      real flat(npx,npy)

      integer i,j,u,v,nxc
      double precision cRPIX(2),cD(2,2),cRVAL(2)

      real sigabc(2,3)
      integer order
      common /sig_pass/ sigabc,order

      real aa,bb,cc

      character*(strl) FLAT_FILE
      integer nxx , nyy
      real flat_weight(npx,npy)


      proc_error=0

c      call readimage(FLAT_FILE,nx,ny,npx,npy,flat)
      call readimage_para(IMAGE_FILE                              
     .,nx,ny,npx,npy,array,cRPIX,cD,cRVAL)

      do i=1,nx
        do j=1,ny
          weight(i,j)=1
c          if (include_FLAT.eq.1) then
c            if (flat(i,j).lt.0.5) then
c              weight(i,j)=0
c            else
c              array(i,j)=array(i,j)*flat(i,j)
c            endif
c          endif
          if (array(i,j).gt.saturation_thresh) weight(i,j)=0
          normap(i,j)=array(i,j)
        enddo
      enddo

c------------------------------------------------------

      nxc=nx/2

      if (ccD_split.eq.2) then

        call set_background(1,nxc,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_background(nxc+1,nx,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nxc,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(1,1)=aa
        sigabc(1,2)=bb
        sigabc(1,3)=cc
        call set_sig(nxc+1,nx,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(2,1)=aa
        sigabc(2,2)=bb
        sigabc(2,3)=cc
      else
        call set_background(1,nx,1,ny,npx,npy,normap        
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nx,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(1,1)=aa
        sigabc(1,2)=bb
        sigabc(1,3)=cc
      endif


c--------------------------------------------------------------

      call get_PREFIX(IMAGE_FILE,PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
      if (ASTROMETRY_trivial.eq.1) then
        call gen_astrometry_data_trivial(cRPIX,cD,cRVAL,filename)
      else
        catfile=ASTROMETRY_cAT
        call generate_gaia_file_name(cRVAL,catfile)
        call gen_astrometry_data(catfile,nx,ny,npx,npy   
     .,normap,weight,cRPIX,cD,cRVAL,filename,proc_error)
      endif

      call locate_defects(nx,ny,npx,npy,array,normap
     .,weight,area_max,area_thresh,proc_error)
      call merge_defects(nx,ny,npx,npy,weight,normap
     .,area_max,source_thresh,area_thresh,proc_error)

c------------------------------------------------------------------

      if (proc_error.eq.0) then
c     Super flat field correction
        if (include_FLAT.eq.1) then
          call readimage(FLAT_FILE,nxx,nyy,npx,npy,flat_weight)
          if (nxx.ne.nx .or. nyy.ne.ny) then
            write(*,*) 'Error: wrong size of flat file!'
            proc_error=1
          endif
          do i=1,nx
            do j=1,ny
              if (flat_weight(i,j).lt.(-900.0)) weight(i,j) = 0
              if (weight(i,j).eq.0) normap(i,j)=-1000. 
            enddo
          enddo    
        else
          do i=1,nx
            do j=1,ny
              if (weight(i,j).eq.0) normap(i,j)=-1000.
            enddo
          enddo
        endif
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

      filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_norm.fits'
      call writeimage_copyhdu(IMAGE_FILE,filename
     .,nx,ny,npx,npy,normap)

c-----------------------------------------------------------
      if (proc_error.eq.0) then
        write(*,*) 'Status of processing ',trim(IMAGE_FILE)
     .,': OK.'
      else
        write(*,*) 'Status of processing ',trim(IMAGE_FILE)
     .,': ERROR!'
      endif

c      pause

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_background(nx1,nx2,ny1,ny2,npx,npy,image
     .,blocksize,nct,ncx,ierror)
      implicit none

c The purpose of this code is to achieve a fine adjustment of the background.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,blocksize,nct,ncx
      real image(npx,npy)

c parameters:
      integer npp,nfit
      parameter (npp=1000)
      parameter (nfit=5000)

c local variables:
      real pix(npp),arr(npp,3),c(nct)
      real mean_sam(nfit,3),sam(nfit),tempo(nfit),tmp_fit(nfit,3)
      integer indx(npp),nsam,nsam1,changed
      integer ix,iy,i,j,k,nbx,nby,xmin,xmax,ymin,ymax,nct_min
      real i2,j2,ij,mean,sig,aa,bb,cc,med,x,y
      real bound1,bound2
c uses:
      real ran1,func_val,ratio

      if (ierror.eq.1) return

c a rough flattening of the field first:
      call flatten_chip(nx1,nx2,ny1,ny2,npx,npy,image,4,2,ierror)

c      return

      if (ierror.eq.1) return

      ratio=1./max((nx2-nx1),(ny2-ny1))

c divide the chip into nbx*nby blocks, each with sidelength of "blocksize":
      nbx=max((nx2-nx1)/blocksize,1)
      nby=max((ny2-ny1)/blocksize,1)

c get the median of each block:
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

c remove the blocks that are outliers:

      nct_min=nct*3/2

      nsam1=nsam
      changed=1
      do while (changed.eq.1 .and. nsam1.ge.nct_min)
        call find_slope_2D(nfit,nsam1,mean_sam,aa,bb,cc)

        do i=1,nsam1
          sam(i)=mean_sam(i,3)-aa-bb*mean_sam(i,1)-cc*mean_sam(i,2)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine flatten_chip(nx1,nx2,ny1,ny2,npx,npy,array
     .,nct,ncx,ierror)
      implicit none

c The purpose of this code is to achieve a rough flattening of the background.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,nct,ncx
      real array(npx,npy)

c parameter:
      integer npp
      parameter (npp=1000)

c local variables:
      real pix(npp),arr(npp,3),arr2(npp,3),c(nct)
      real aa,bb,cc,arr_min,arr_max,x,y
      integer ix,iy,i,j,nr
c Uses:
      real ran1,func_val,ratio

      if (ierror.eq.1) return

      ratio=1./max(nx2-nx1,ny2-ny1)

c pick "npp" random positions:
      do i=1,npp
        ix=int(ran1()*(nx2-nx1)+nx1)
        iy=int(ran1()*(ny2-ny1)+ny1)
        pix(i)=array(ix,iy)
        arr2(i,1)=ix*ratio
        arr2(i,2)=iy*ratio
        arr2(i,3)=pix(i)
      enddo
      call sort(npp,npp,pix)

      if (pix(1).eq.pix(npp)) then
        ierror=1
        return
      endif

c remove the pixels of extreme values:
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

      call fit_2D(npp,nr,arr,nct,ncx,c)

      do i=nx1,nx2
        do j=ny1,ny2
          x=i*ratio
          y=j*ratio
          array(i,j)=array(i,j)-func_val(x,y,nct,ncx,c)
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_sig(nx1,nx2,ny1,ny2,npx,npy,image
     .,aa,bb,cc,ierror)
      implicit none

c The purpose of this code is to get the rms of noise as a function of position.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror
      real image(npx,npy)

c parameters:
      integer npp
      parameter (npp=2000)

c local variables:
      real pix(npp),arr(npp,3),arr2(npp,3),c(6)
      integer ix,iy,i,j,nr
      real aa,bb,cc,arr_min,arr_max,temp
c uses:
      real ran1

      if (ierror.eq.1) return

c get a sample of points:

      do i=1,npp
        ix=int(ran1()*(nx2-nx1-1)+nx1)
        iy=int(ran1()*(ny2-ny1-1)+ny1)
        pix(i)=0.5*((image(ix,iy)-image(ix+1,iy))**2
     .+(image(ix,iy)-image(ix,iy+1))**2)
        arr(i,1)=ix
        arr(i,2)=iy
        arr(i,3)=pix(i)
      enddo
      call sort(npp,npp,pix)

      arr_min=pix(npp/3)
      arr_max=pix(2*npp/3)

c remove the points of extreme values:
      nr=0
      do i=1,npp
        if (arr(i,3).ge.arr_min.and.arr(i,3).le.arr_max) then
          nr=nr+1
          arr2(nr,1)=arr(i,1)
          arr2(nr,2)=arr(i,2)
          arr2(nr,3)=arr(i,3)
        endif
      enddo

c fit a 2D linear function with the remaining points:

      if (nr.lt.npp/10) then
        ierror=1
        return
      endif

      call find_slope_2D(npp,nr,arr2,aa,bb,cc)

      do i=nx1,nx2
        do j=ny1,ny2
          temp=aa+bb*i+cc*j
          image(i,j)=image(i,j)/sqrt(0.5*temp)
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      subroutine detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)
      implicit none

      integer nx,ny,npx,npy
      integer y_smooth,x_smooth
      real normap(npx,npy),ymap(nx,ny/y_smooth),xmap(nx/x_smooth,ny)
      integer weight(npx,npy)

      integer np,i,j,startj,jj,numy,numx,starti,ii
      real sig,med

c      character filename*100

      numy=ny/y_smooth
      do i=1,nx
        do j=1,numy
          ymap(i,j)=0.
          startj=(j-1)*y_smooth
          do jj=1,y_smooth
            if (weight(i,startj+jj).eq.1) then
              ymap(i,j)=ymap(i,j)+normap(i,startj+jj)
            endif
          enddo
        enddo
      enddo
      call get_sig_med(nx,numy,ymap,sig,med)
      do i=1,nx
        do j=1,numy
          if (abs(ymap(i,j)-med).gt.sig*4.) then
            startj=(j-1)*y_smooth
            do jj=1,y_smooth
              weight(i,startj+jj)=0
            enddo
          endif
        enddo
      enddo

      numx=nx/x_smooth
      do j=1,ny
        do i=1,numx
          xmap(i,j)=0.
          starti=(i-1)*x_smooth
          do ii=1,x_smooth
            if (weight(starti+ii,j).eq.1) then
              xmap(i,j)=xmap(i,j)+normap(starti+ii,j)
            endif
          enddo
        enddo
      enddo
      call get_sig_med(numx,ny,xmap,sig,med)
      do j=1,ny
        do i=1,numx
          if (abs(xmap(i,j)-med).gt.sig*4.) then
            starti=(i-1)*x_smooth
            do ii=1,x_smooth
              weight(starti+ii,j)=0
            enddo
          endif
        enddo
      enddo


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
