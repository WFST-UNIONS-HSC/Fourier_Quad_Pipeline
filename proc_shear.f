      subroutine proc_shear(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      call expo_shear(nchip,IMAGE_FILE,DIR_OUTPUT)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine expo_shear(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer ichip,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      integer proc_error,nstar,status,nstar_tot
      character*(strl) PREFIX,headname,filename,PREFIX1,PREFIX2
      character*(strl) PREFIX_head
      integer ngal,ierror,nn1,nn2,i,j,u,v,nx,ny,k
      real g1,g2,de,h1,h2
      real gal_p_coll(ngal_max,ns,ns),gal_para(ngal_max,npara)
      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd)
      real gal_p(ns,ns),psf_model(ns,ns),aa(npara)
      double precision ra,dec,gf1,gf2,temp,cos2,sin2,cos4,sin4
      double precision x,y,xx,yy
      integer parity
      double precision local_coe(ns,ns,npl),PSF_coe(ns,ns,npo)
      real ePSF(ns,ns),ePSF_p(ns,ns),psf_FWHM,psfmap(npx,npy)
      real px,py,dstar

      proc_error=0

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      headname=trim(DIR_OUTPUT)//'/astrometry/'//trim(PREFIX)//'.head'

      if (ext_PSF.eq.1) then
        filename=trim(PSF_PATH)//'/PSF.fits'
        call readimage(filename,nx,ny,ns,ns,ePSF)
        call get_power(ns,ns,ePSF,ePSF_p,0)
        do i=1,ns
          do j=1,ns
            local_coe(i,j,1)=ePSF_p(i,j)
          enddo
        enddo
      endif


      do ichip=1,nchip
        proc_error=0
        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        PREFIX1=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        PREFIX2=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)
        if (ext_PSF.ne.1 .and. PSF_type.eq.1) then
          filename=trim(PREFIX1)//'_PSF_coe_local.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          read(10,*) nstar,status
          if (status.eq.-1) then
            proc_error=1
          else
            do i=1,ns
              do j=1,ns
                read(10,*) (local_coe(i,j,k),k=1,npl)
              enddo
            enddo
          endif
          close(10)
        elseif (ext_PSF.ne.1 .and. PSF_type.eq.2) then
          filename=trim(PREFIX1)//'_PSF_local.fits'
          call readimage(filename,nx,ny,npx,npy,psfmap)
          nstar=int(psfmap(step_psf-1,step_psf-1)+0.5)
          if (psfmap(step_psf,step_psf).lt.-1.) proc_error=1
        endif
        if (proc_error.eq.1) then
          ngal=0
          goto 50
        endif
        call read_astrometry_para(headname,ichip,cRPIX,cD,cRVAL,PU,npd
     .,proc_error)
        if (proc_error.eq.1) then
          ngal=0
          goto 50
        endif

        ngal=0
        filename=trim(PREFIX1)//'_source_info.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) filename
          stop 'catalog file error in shear_proc!!'
        endif
        read(10,*)
c        read(10,*) 'ig xc yc sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F'
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,iSNR_F)
          if (ierror.lt.0) cycle
          ngal=ngal+1
          do i=1,iSNR_F
            gal_para(ngal,i)=aa(i)
          enddo
        enddo
        close(10)

        if (ngal.eq.0) then
          goto 50
        endif

        nn1=ns*len_g
        nn2=ns*(int(ngal/len_g)+1)
        filename=trim(PREFIX1)//'_source_p.fits'
        call read_stamps(ngal_max,1,ngal,ns,ns,gal_p_coll
     .,nn1,nn2,filename)

50        filename=trim(PREFIX2)//'_shear.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) 'ig xc yc sigma nstar imax jmax ' 
     .,'half_light_flux half_light_area flag psf_FWHM SNR_F '  
     .,'ra dec gf1 gf2 g1 g2 de h1 h2 cos2 sin2 parity'

        do i=1,ngal
          do u=1,ns
            do v=1,ns
              gal_p(u,v)=gal_p_coll(i,u,v)
            enddo
          enddo
          x=gal_para(i,2)
          y=gal_para(i,3)
          if (ext_PSF.eq.1) then
            call get_PSF_model(ns,1,1,local_coe,x,y,psf_model)
          else
            if (PSF_type.eq.1) then
              call get_PSF_model(ns,npl,nplx,local_coe,x,y,psf_model)
            elseif (PSF_type.eq.2) then
              px=x
              py=y
              call get_PSF_model_very_local(psfmap,px,py
     .,psf_model,dstar)
            endif
          endif
          call get_PSF_area(psf_model,psf_FWHM)

          gal_para(i,iPSF)=psf_FWHM
          gal_para(i,istar)=nstar
          call coordinate_transfer_PU(ra,dec,x,y,1,cRPIX,cD
     .,cRVAL,PU,npd)
          gal_para(i,ira)=ra
          gal_para(i,idec)=dec
          call field_distortion_PU(x,y,npd,PU,cD,cRPIX,gf1,gf2
     .,cos2,sin2,parity)
          gal_para(i,igf1)=gf1
          gal_para(i,igf2)=gf2
          call get_shear(ns,gal_p,psf_model,g1,g2,de,h1,h2)
          gal_para(i,ig1)=g1*cos2+g2*sin2
          gal_para(i,ig2)=g2*cos2-g1*sin2
          gal_para(i,ide)=de
          cos4=cos2*cos2-sin2*sin2
          sin4=2d0*sin2*cos2
          gal_para(i,ih1)=h1*cos4+h2*sin4
          gal_para(i,ih2)=h2*cos4-h1*sin4
          if (parity.eq.-1) then
            gal_para(i,ig2)=-gal_para(i,ig2)
            gal_para(i,ih2)=-gal_para(i,ih2)
          endif
          gal_para(i,icos2)=cos2
          gal_para(i,isin2)=sin2
          gal_para(i,iparity)=parity

          write(10,*) (gal_para(i,j),j=1,iparity)

        enddo
        close(10)
      enddo

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE get_shear(n,gal,psf,g1,g2,de,h1,h2)
      implicit none

      integer n
      real g1,g2,de,h1,h2,gal(n,n),psf(n,n)
      integer i,j,n_2,cc
      real ks,peak,thresh,area,ks_2,kx,kx2,ky,ky2,k2,k,temp,temp1
      real filter_deriv,filter,ff,norm

      real pi,PSFr_ratio
      parameter (pi=3.1415926)
      parameter (PSFr_ratio=0.75)

      peak=psf(1,1)
      do i=1,n
        do j=1,n
          if (psf(i,j).gt.peak) peak=psf(i,j)
        enddo
      enddo
ccccccccccccccccccc
      thresh=exp(-1.)*peak
c     thresh=exp(-0.5)*peak
ccccccccccccccccccc
      area=0.
      do i=1,n
        do j=1,n
          if (psf(i,j).ge.thresh) area=area+1.
        enddo
      enddo

      ks=sqrt(area/pi)
      ks_2=(ks*PSFr_ratio)**(-2)

      thresh=peak*1e-5
      n_2=n/2
      cc=1+n_2

      g1=0.
      g2=0.
      de=0.
      h1=0.
      h2=0.

      do i=1,n
        kx=i-cc
        kx2=kx*kx
        do j=1,n
          ky=j-cc
          ky2=ky*ky
          k2=kx2+ky2
          k=sqrt(k2)
          if (psf(i,j).gt.thresh) then
            ff=k2*ks_2
            temp=exp(-ff)/psf(i,j)              ! 相当于T
            temp1=temp*gal(i,j)                 ! 相当于T*M

            g1=g1-temp1*(kx2-ky2)
            g2=g2-temp1*2.*kx*ky
            de=de+temp1*k2*(2.-ff)              ! N
            h1=h1+temp1*ks_2*(k2*k2-8.*kx2*ky2) ! U
            h2=h2+temp1*ks_2*4.*kx*ky*(kx2-ky2)

          endif
        enddo
      enddo


      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PSF_area(model,FWHM)
      implicit none
      include 'para.inc'

      real model(ns,ns),area,thresh,FWHM,r,beta
      integer i,j

      thresh=exp(-1.) * model(ns/2+1,ns/2+1)

      area=-1e-5
      do i=1,ns
        do j=1,ns
          if (model(i,j).ge.thresh) area=area+1.
        enddo
      enddo
      if (area.le.0.) then
        FWHM=-1.
        return
      endif
      r=sqrt(area/pi)
      beta=ns/(2.*pi)/r
      FWHM=beta*2.*sqrt(2.*log(2.))*0.2628

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc