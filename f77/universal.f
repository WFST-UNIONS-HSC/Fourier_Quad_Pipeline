      subroutine get_chi2(ns,n1,n2,image,model,chi2)
      implicit none

      integer ns,n1,n2
      real image(ns,ns)
      integer i,j
      real model(ns,ns),cc,chi2

      cc=ns/2.+1.

      chi2=0.
      do i=n1,n2
        do j=n1,n2
          chi2=chi2
     .+((image(i,j)-model(i,j))*((i-cc)**2+(j-cc)**2))**2
        enddo
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ana_chi2(n,map1,map2,p)
      implicit none

      integer n
      real map1(n,n),map2(n,n),p,flux
      integer i,j,n1,n2

      n1=n/4
      n2=(n/4)*3
      p=0.
      flux=0.
      do i=n1,n2
        do j=n1,n2
          flux=flux+(map1(i,j)+map2(i,j))*0.5
          p=p+(map1(i,j)-map2(i,j))**2
        enddo
      enddo
      p=p/flux

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fit_func(x,y,nc,nx,n)
      implicit none

      real x,y,fit_func
      integer nc,nx,n
      integer px,py

      px=mod(n-1,nx)
      py=(n-1)/nx

      fit_func=x**px*y**py

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function func_val(x,y,nc,nx,c)
      implicit none

      real x,y,func_val
      integer nc,nx
      real c(nc),fit_func
      integer i

      func_val=0.
      do i=1,nc
        func_val=func_val+fit_func(x,y,nc,nx,i)*c(i)
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fit_2D(np,n,arr,nc,nx,c)
      implicit none

c      fit a 2D function with f(x,y)=\Sigma_{i=1}^{nc}func_i(x,y)c_i

      integer np,n,nc,nx
      real arr(np,3),c(nc)
      real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
      real x,y,f,smax
      real fit_func
      external fit_func
      integer i,j,u,v

      do i=1,nc
        do j=1,nc
          coe(i,j)=0.
        enddo
        vec(i)=0.
      enddo

      do i=1,n
        x=arr(i,1)
        y=arr(i,2)
        f=arr(i,3)

        do j=1,nc
          temp(j)=fit_func(x,y,nc,nx,j)
        enddo

        do u=1,nc
          do v=1,nc
            coe(u,v)=coe(u,v)+temp(u)*temp(v)
          enddo
          vec(u)=vec(u)+f*temp(u)
        enddo

      enddo

      call matrix_inverse(coe,nc,coe_1)

      do j=1,nc
        c(j)=0.
        do i=1,nc
          c(j)=c(j)+coe_1(j,i)*vec(i)
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fit_func_2(x,y,n)
      implicit none

      real x,y,fit_func_2
      integer n,px,py,order,nn

      px=0
      py=-1
      order=-1
      do nn=1,n
        if (py.eq.order) then
          order=order+1
          px=order
          py=0
        else
          px=px-1
          py=py+1
        endif
      enddo

      fit_func_2=x**px*y**py

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function func_val_2(x,y,nc,c)
      implicit none

      real x,y,func_val_2
      integer nc,i
      real c(nc),fit_func_2
      external fit_func_2

      func_val_2=0.
      do i=1,nc
        func_val_2=func_val_2+fit_func_2(x,y,i)*c(i)
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fit_2D_2(n,np,arr,nc,c)
      implicit none

c      fit a 2D function with f(x,y)=\Sigma_{i=1}^{nc}func_i(x,y)c_i

      integer n,np,nc
      real arr(np,3),c(nc)
      real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
      real x,y,f,smax
      real fit_func_2
      external fit_func_2
      integer i,j,u,v

      coe=0.
      vec=0.

      do i=1,n
        x=arr(i,1)
        y=arr(i,2)
        f=arr(i,3)
        do j=1,nc
          temp(j)=fit_func_2(x,y,j)
        enddo
        do u=1,nc
          do v=1,nc
            coe(u,v)=coe(u,v)+temp(u)*temp(v)
          enddo
          vec(u)=vec(u)+f*temp(u)
        enddo
      enddo

      call matrix_inverse(coe,nc,coe_1)

      do j=1,nc
        c(j)=0.
        do i=1,nc
          c(j)=c(j)+coe_1(j,i)*vec(i)
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_slope_2D(np,n,arr,aa,bb,cc)
      implicit none

      integer np,n
      real arr(np,3),aa,bb,cc
      integer i,j
      real c(3,3),vec(3),c_1(3,3)

      do i=1,3
        do j=1,3
          c(i,j)=0.
        enddo
        vec(i)=0.
      enddo

      do i=1,n
        c(1,1)=c(1,1)+1.
        c(1,2)=c(1,2)+arr(i,1)
        c(1,3)=c(1,3)+arr(i,2)
        vec(1)=vec(1)+arr(i,3)

        c(2,1)=c(2,1)+arr(i,1)
        c(2,2)=c(2,2)+arr(i,1)*arr(i,1)
        c(2,3)=c(2,3)+arr(i,2)*arr(i,1)
        vec(2)=vec(2)+arr(i,3)*arr(i,1)

        c(3,1)=c(3,1)+arr(i,2)
        c(3,2)=c(3,2)+arr(i,1)*arr(i,2)
        c(3,3)=c(3,3)+arr(i,2)*arr(i,2)
        vec(3)=vec(3)+arr(i,3)*arr(i,2)
      enddo

      call matrix_inverse(c,3,c_1)

      aa=0.
      bb=0.
      cc=0.
      do i=1,3
        aa=aa+c_1(1,i)*vec(i)
        bb=bb+c_1(2,i)*vec(i)
        cc=cc+c_1(3,i)*vec(i)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matrix_inverse(ma,n,ma_1)
      implicit none

      integer n
      real ma(n,n),ma_1(n,n),d,b(n)
      integer indx(n),i,j

      call ludcmp(ma,n,n,indx,d)

      do i=1,n
        do j=1,n
          if (j.eq.i) then
            b(j)=1.
          else
            b(j)=0.
          endif
        enddo

        call lubksb(ma,n,n,indx,b)

        do j=1,n
          ma_1(j,i)=b(j)
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine generate_gaia_file_name(cRVAL,filename)
      implicit none

      double precision cRVAL(2)
      character filename*(*),c_ra*2,c_dec
      integer ra,dec

10    format(I1.1)
20    format(I2.2)

      dec=min(int(abs(cRVAL(2))/10d0)+1,9)
      write(c_dec,10) dec

      ra=min(int(cRVAL(1)/10d0),35)
      write(c_ra,20) ra

      if (cRVAL(2).ge.0) then
        filename=trim(filename)//'/gaia_p'
      else
        filename=trim(filename)//'/gaia_m'
      endif

      filename=trim(filename)//c_dec
      if (dec.eq.9) then
        filename=trim(filename)//'.cat'
      else
        filename=trim(filename)//'_'//c_ra//'.cat'
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine generate_gal_cat_file_name(cRVAL,filename
     .,sortfile,sortnum)
      implicit none
      include 'para.inc'

      double precision cRVAL(2) , m_ra , m_dec
      character filename*(*)
      character*(strl) sortfile(27)
      integer dec1,dec2
      character c_dec*2
      character c_ra*3
      integer ra1,ra2
      integer ra,dec,sortnum,mra,mdec

10    format(I2.2)
20    format(I3.3)

      m_dec=1.0d0
      ! first judge the dec
      if (abs(cRVAL(2)).lt.30.) then
        m_ra=1.4d0
      elseif (abs(cRVAL(2)).lt.40.) then
        m_ra=1.6d0
      elseif (abs(cRVAL(2)).lt.50.) then
        m_ra=1.8d0
      elseif (abs(cRVAL(2)).lt.60.) then
        m_ra=2.5d0
      else
        m_ra=3.0d0
      endif
      ! --------------dec-------------------------
      dec1=int(floor(cRVAL(2)-m_dec))
      dec2=int(floor(cRVAL(2)+m_dec))
      !--------------- ra -----------------------
      ra1=int(floor(cRVAL(1)-m_ra))
      ra2=int(floor(cRVAL(1)+m_ra))
      !------------------------------------------
      ! write in the file name
      sortnum=0
      do dec=dec1,dec2
        do ra=ra1,ra2
          sortnum=sortnum+1
          mra=ra
          mdec=dec
          if (ra.lt.0) then
            mra=mra+360
          elseif (ra.ge.360) then
            mra=mra-360
          endif 
          write(c_ra,20) mra

          if (dec.ge.0) then
            sortfile(sortnum)=trim(filename)//'/gal_p'
            write(c_dec,10) mdec
          elseif (dec.lt.0) then
            sortfile(sortnum)=trim(filename)//'/gal_m'
            mdec = -mdec - 1
            write(c_dec,10) mdec
          endif
          sortfile(sortnum)=trim(sortfile(sortnum))//c_dec
     .//'_'//c_ra//'.cat'
        enddo
      enddo
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
        if (imagefile(i:i).eq.'.' .and. p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/' .and. p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0 .or. p_slash.eq.0. .or. p_dot.le.p_slash+1) then
        write(*,*) 'Image_file name is NOT normal !'
        read(*,*)
      endif

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PREFIX_expo(imagefile,PREFIX)
      implicit none

      character*(*) imagefile,PREFIX
      integer i,p_dot,p_slash,n

      p_dot=0
      p_slash=0

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'_' .and. p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/' .and. p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0 .or. p_slash.eq.0. .or. p_dot.le.p_slash+1) 
     .stop 'Image_file name is NOT normal !'

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_chip_id(imagefile,id)
      implicit none

      character*(*) imagefile
      character*2 id_letter
      integer i,p1,p2,n,id

c      p1=0
c      p2=0
c      n=len(trim(imagefile))
c      do i=n,1,-1
c        if (imagefile(i:i).eq.'.'.and.p1.eq.0) then
c          p1=i
c          cycle
c        endif
c        if (imagefile(i:i).eq.'_'.and.p2.eq.0) then
c          p2=i
c          cycle
c        endif
c      enddo
c      if (p1.eq.0.or.p2.eq.0) stop 'Image name is NOT normal !'
c      id_letter=imagefile(p2+1:p1-1)
c      read(id_letter,*) id

      call read_ccDNUM(imagefile,id)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_med_sig(np,n,dat,med,sig)
      implicit none

      integer n,np
      real dat(np),med,sig,arr(n)
      integer i

      do i=1,n
        arr(i)=dat(i)
      enddo

      call sort(n,n,arr)

      med=arr(n/2)
      sig=0.5*(arr(n*5/6)-arr(n/6))

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function iden(x,n)
      implicit none

      integer n
      real x,iden

      iden=x

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function loga(x,n)
      implicit none

      integer n
      real x,loga

      if (n.eq.1) then
        if (x.lt.-30.) then
          loga=-log(-2.*x)
        else
          loga=log(x+sqrt(x**2+1.))
        endif
      else
        loga=0.5*(exp(x)-exp(-x))
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fit_linear_2D(n,np,x,xt,coe)
      implicit none

      integer n,np
      double precision x(np,2),xt(np,2),coe(2,3)
      double precision matx(3,3),bvec1(3),bvec2(3),vec(3)
      double precision matx_1(3,3)
      integer i,j,u,v


      do u=1,3
        do v=1,3
          matx(u,v)=0d0
        enddo
        bvec1(u)=0d0
        bvec2(u)=0d0
      enddo

      do i=1,n
c        xt(i,1)=coe(1,1)*x(i,1)+coe(1,2)*x(i,2)+coe(1,3)
c        xt(i,2)=coe(2,1)*x(i,1)+coe(2,2)*x(i,2)+coe(2,3)

        vec(1)=x(i,1)
        vec(2)=x(i,2)
        vec(3)=1d0

        do u=1,3
          do v=1,3
            matx(u,v)=matx(u,v)+vec(u)*vec(v)
          enddo
          bvec1(u)=bvec1(u)+xt(i,1)*vec(u)
          bvec2(u)=bvec2(u)+xt(i,2)*vec(u)
        enddo
      enddo

      call matrix_inverse_doub(matx,3,3,matx_1)

      do i=1,3
        coe(1,i)=0d0
        coe(2,i)=0d0
        do j=1,3
          coe(1,i)=coe(1,i)+matx_1(i,j)*bvec1(j)
          coe(2,i)=coe(2,i)+matx_1(i,j)*bvec2(j)
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matrix_inverse_doub(ma,n,np,ma_1)
      implicit none

      integer n,np
      double precision ma(np,np),ma_1(np,np),d,b(n),norm
      integer indx(n),i,j

      norm=0d0
      do i=1,n
        do j=1,n
          norm=max(norm,abs(ma(i,j)))
        enddo
      enddo

      do i=1,n
        do j=1,n
          ma(i,j)=ma(i,j)/norm
        enddo
      enddo

      call ludcmp_doub(ma,n,np,indx,d)

      do i=1,n
        do j=1,n
          if (j.eq.i) then
            b(j)=1d0
          else
            b(j)=0d0
          endif
        enddo

        call lubksb_doub(ma,n,np,indx,b)

        do j=1,n
          ma_1(j,i)=b(j)
        enddo
      enddo

      do i=1,n
        do j=1,n
          ma_1(i,j)=ma_1(i,j)/norm
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_peak_width(np,n,arr,p,sig,status,direc)
      implicit none

      integer np,n,status,direc
      real arr(np),p,sig,d
      integer i,j,ip,u,v,b1,b2
      integer nb,smooth
      parameter (nb=200)
      parameter (smooth=3)
      real den(0:nb+1),posi(nb),vol1,vol2,tmp,thresh,vol
      integer mark(0:nb+1),bound1,bound2,num,bb(2,2)

      call sort(n,np,arr)

      if (direc.gt.0) n=(n*9)/10

      d=(arr(n)-arr(1))/(nb-1.)

      do i=1,nb
        posi(i)=arr(1)+(i-1.)*d
        den(i)=0.
        mark(i)=0
      enddo
      den(0)=0.
      den(nb+1)=0.
      mark(0)=0
      mark(nb+1)=0

      do i=1,n
        tmp=(arr(i)-arr(1))/d+1
        ip=int(tmp+0.5)
        do j=max(ip-4*smooth,1),min(ip+4*smooth,nb)
          den(j)=den(j)+exp(-0.5*((tmp-j)/smooth)**2)
        enddo
      enddo

      vol1=0
      vol2=0
      bb(2,1)=0
      bb(2,2)=0
      bb(1,1)=0
      bb(1,2)=0

      do i=1,nb
        if (mark(i).gt.0) cycle
        mark(i)=1
        if (den(i).gt.den(i-1).and.den(i).gt.den(i+1)) then
          bound1=i
          bound2=i
          ip=i
          thresh=den(i)*0.5
          do while ((den(bound2+1).gt.thresh
     . .or. den(bound2+1).lt.den(bound2)).and.bound2.lt.nb)
            bound2=bound2+1
            mark(bound2)=2
            if (den(bound2)*0.5.gt.thresh) then
              thresh=den(bound2)*0.5
              ip=bound2
            endif
          enddo
          do while ((den(bound1-1).gt.thresh
     . .or. den(bound1-1).lt.den(bound1)).and.bound1.gt.1
     . .and. mark(bound1-1).ne.2)
            bound1=bound1-1
            if (den(bound1)*0.5.gt.thresh) then
              thresh=den(bound1)*0.5
              ip=bound1
            endif
          enddo

          thresh=thresh*0.5
          b2=ip
          do while (b2.lt.bound2.and.den(b2).gt.thresh)
            b2=b2+1
          enddo
          b1=ip
          do while (b1.gt.bound1.and.den(b1).gt.thresh)
            b1=b1-1
          enddo
          bound1=b1
          bound2=b2

          vol=0.
          do j=bound1,bound2
            vol=vol+den(j)
          enddo

          if (vol.gt.vol1) then
            vol2=vol1
            bb(2,1)=bb(1,1)
            bb(2,2)=bb(1,2)
            vol1=vol
            bb(1,1)=bound1
            bb(1,2)=bound2
          elseif (vol.gt.vol2) then
            vol2=vol
            bb(2,1)=bound1
            bb(2,2)=bound2
          endif
        endif
      enddo

      if (direc.gt.0) then
        if (bb(1,1).gt.bb(2,1)) then
          bound1=bb(1,1)
          bound2=bb(1,2)
          status=1
        else
          bound1=bb(2,1)
          bound2=bb(2,2)
          status=2
        endif
      else
        if (bb(1,1).gt.bb(2,1)) then
          bound1=bb(2,1)
          bound2=bb(2,2)
          status=2
        else
          bound1=bb(1,1)
          bound2=bb(1,2)
          status=1
        endif
      endif

      p=0
      sig=0
      num=0
      do i=1,n
        if (arr(i).lt.posi(bound1).or.arr(i).gt.posi(bound2)) cycle
        p=p+arr(i)
        sig=sig+arr(i)**2
        num=num+1
      enddo
      p=p/num
      sig=max(sqrt(sig/num-p*p),0.02*(posi(bound2)-posi(bound1)))


c      write(*,*) bound1,bound2,num,posi(bound1),posi(bound2)
c      do i=1,nb
c        write(*,*) 'den(posi):',i,posi(i),den(i)
c      enddo
c      write(*,*) bb(1,1),bb(1,2),vol1
c      write(*,*) bb(2,1),bb(2,2),vol2
c      write(*,*) 'final:',p,sig,p-3*sig,p+3*sig,status


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_peak_width_low_side(np,n,arr,p,sig)
      implicit none

      integer np,n
      real arr(np),p,sig,den(n),den2(n),thresh
      integer i,j,ip

      call sort(n,np,arr)

      do i=2,n-1
        den(i)=((arr(i+1)-arr(i-1))/2.)**2
      enddo

      den(1)=(arr(2)-arr(1))**2
      den(n)=(arr(n)-arr(n-1))**2

      den2=0.
      do i=3,n-2
        do j=i-2,i+2
          den2(i)=den2(i)+den(j)
        enddo
        den2(i)=1./sqrt(den2(i)/5.)
      enddo

      den2(1)=1./sqrt((den(1)+den(2)+den(3))/3.)
      den2(2)=1./sqrt((den(1)+den(2)+den(3)+den(4))/4.)
      den2(n)=1./sqrt((den(n)+den(n-1)+den(n-2))/3.)
      den2(n-1)=1./sqrt((den(n)+den(n-1)+den(n-2)+den(n-3))/4.)

      ip=1
      p=den2(1)

      do i=2,n
        if (den2(i).gt.p) then
          p=den2(i)
          ip=i
        endif
      enddo

      thresh=p*0.5
      do i=ip+1,n
        if (den2(i).lt.thresh) exit
      enddo

      p=arr(ip)
      sig=arr(i)-arr(ip)

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
      if (ierr.ne.0) stop 'EXPO_FILE reading error!!'
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
      Function bound_ra(ra_input)
      implicit none
      include 'para.inc'

      double precision ra_input,bound_ra

      if (ra_input.ge.360d0) then
        bound_ra=ra_input-360d0
      elseif (ra_input.lt.0d0) then
        bound_ra=ra_input+360d0
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc