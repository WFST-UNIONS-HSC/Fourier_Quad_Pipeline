      real pi,arc_convert
      parameter (pi=3.1415926)
      parameter (arc_convert=pi/180.)
      real pixel_size
      parameter (pixel_size=0.2628)
c in unit of arcsec

      integer strl,NMAX_EXPO,NMAX_CHIP
      parameter (strl=150)
      parameter (NMAX_EXPO=25000)
      parameter (NMAX_CHIP=62)
	
      integer nmax_per_core
      parameter (nmax_per_core=4500000)

      integer n_target_subcat
      parameter (n_target_subcat=500000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer iSNR_F,igf1,igf2,ira,idec,ih_flux,iPSF
      integer i_imax,ig1,ig2,ide,ih1,ih2,ipeak,iflag
      integer i_jmax,ih_area,istar,isig
      integer icos2,isin2,iparity,ichi2
      integer ccd_num,icdec,icra,izp
      integer imag_i
      integer iflags_ft,iflags_fg,iflags_gold,iext_mash,iid
      parameter (iflags_ft=1)
      parameter (iflags_fg=2)
      parameter (iflags_gold=3)
      parameter (iext_mash=4)
      parameter (icra=5)
      parameter (icdec=6)
      parameter (imag_i=11)
      parameter (ccd_num=1+18)
      parameter (iid=1+ccd_num)
      parameter (izp=17)

      parameter (isig=4+ccd_num)
      parameter (istar=5+ccd_num)
      parameter (ipeak=5+ccd_num)
      parameter (i_imax=6+ccd_num)
      parameter (i_jmax=7+ccd_num)
      parameter (ih_flux=8+ccd_num)
      parameter (ih_area=9+ccd_num)
      parameter (iflag=10+ccd_num)
      parameter (iPSF=11+ccd_num)
      parameter (iSNR_F=12+ccd_num)
      parameter (ira=13+ccd_num)
      parameter (idec=14+ccd_num)
      parameter (igf1=15+ccd_num)
      parameter (igf2=16+ccd_num)
      parameter (ig1=17+ccd_num)
      parameter (ig2=18+ccd_num)
      parameter (ide=19+ccd_num)
      parameter (ih1=20+ccd_num)
      parameter (ih2=21+ccd_num)
      parameter (icos2=22+ccd_num)
      parameter (isin2=23+ccd_num)
      parameter (iparity=24+ccd_num)
      parameter (ichi2=25+ccd_num)
cccccccccccccccccccccccccccccccccccccccccccccc
