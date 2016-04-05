      subroutine getmodel11(zin,agein,optdir,nchar,
     &          hlam,wave,bolo,z,y,age,mags,Lsq,Lavg,mtol,
     &       widths,Fcont,Fline,indexunits,labels,maglabels,
     &       massinit,massnow,imfpar,zsolar,
     &       ysolar,yprimo,imflow,imfupp,hbmorph,pagb,eta,imftype,
     &       status,high)

c     input variables: what age and metallicity do you want?
      real         zin,agein
      character    optdir*25
c     output variables - storage of model results
      integer      nind, status, nchar, fstatus, nmag
      parameter    (nind = 179, nmag=145)
C     primary output variables
      real         hlam(1221),wave(1221), bolo, z, y, age
      real         high(14000)
      real         mags(nmag),Lsq(nmag),Lavg(nmag),mtol(nmag)
      real         widths(nind),Fcont(nind),Fline(nind)
      integer      indexunits(nind)
      character*10 labels(nind)
      character*12 maglabels(nmag)
      character*20 imftype
      real         massinit, massnow
C     other output variables
      real         imfpar, zsolar, ysolar, yprimo, imflow, imfupp
      real         hbmorph, pagb, eta

C     local variables: versions of the primary output variables
      real         hlam2(1221,2,2), bolo2(2,2), high2(14000,2,2)
      real         mags2(nmag,2,2),Lsq2(nmag,2,2)
      real         Lavg2(nmag,2,2),mtol2(nmag,2,2)
      real         Fcont2(nind,2,2),Fline2(nind,2,2)

C     local variables - navigating the model grid points
      real  AA(18),ZZ(9),fehlist(9)
      real ztarg, atarg
      integer   jz, ja
      character*80 infile


C     [Fe/H] -2.23 -1.63 -1.23 -0.93 -0.53 -0.33 0.00 0.37 0.62
      data fehlist / -2.23, -1.63, -1.23, -0.93, -0.53, -0.33, 
     z          0.00, 0.37, 0.62 /
c      data ZZ / .000169, .0005344, .00169, .0053442, .0095036, .0100667,
c     z          .0169,   .030053, .0534425 /
      data AA / 0.0112, 0.0200, 0.0355, 0.0631, 0.1000, 0.1778, 0.3162,
     z          0.5623, 1.0000, 1.7783, 3.1623, 5.6234, 10.0, 11.2202,
     z          12.5893, 14.1254, 15.8490, 17.7828 /

      do i=1,9
         ZZ(i) = 0.0169*10.0**fehlist(i)
      end do

C     check for out-of-grid request. If so, return with error message
      ierr = 1
      if ( zin .lt. .0000995 ) ierr = -1
      if ( zin .gt. .0705 ) ierr = -1
      if ( agein .lt. 0.011 ) ierr = -1
      if ( agein .gt. 18.01 ) ierr = -1
      if ( ierr .lt. 0 ) then
         print*, 'Your request was out of the bounds covered by the',
     z         ' model grid.'
C        returning status = 5  means bad news. Otherwise we leave status alone
         status = 5
         return
      end if

C     set up for interpolation
C     find corners of grid
      call locate (AA,18,agein,ja)
      call locate (ZZ,9,zin,jz)
      if ( ja.eq.0) ja = 1
      if ( ja.eq.18) ja = 17
      if ( jz.eq.0) jz = 1
      if ( jz.eq.9) jz = 8

      do kaloop = ja,ja+1
         do kzloop = jz, jz+1
            atarg = AA(kaloop)
            ztarg = ZZ(kzloop)
            write(infile,'(a23,a11,a2,i3.3,a1,i4.4,a5)') 
     z        '/space/www/dial/export/',optdir(1:11),'/a',
     z        int(10.0*atarg),'z',int(100000.1*ztarg),'.fits'
            ka = kaloop - ja + 1
            kz = kzloop - jz + 1

C           actually read in the selected model data via FITSIO routines
            call rfits2(infile,high2(1,ka,kz),hlam2(1,ka,kz),wave,
     &       bolo2(ka,kz),z,y,age,
     &       mags2(1,ka,kz),Lsq2(1,ka,kz),Lavg2(1,ka,kz),mtol2(1,ka,kz),
     &       widths,Fcont2(1,ka,kz),Fline2(1,ka,kz),
     &       indexunits,labels,massinit,massnow,imfpar,zsolar,ysolar,
     &       yprimo,imflow,
     &       imfupp,hbmorph,pagb,eta,imftype,maglabels,fstatus)

c            print*, 'Debugging info:'
c            print*, 'Labels.',labels(1),labels(101)

         end do
      end do

C     relevant data should be stored in the 2x2 "2" arrays. Interpolate
C     linearly in log AGE and log Z. First Age.
      afrac = (alog10(agein) - alog10(AA(ja)))/
     &               (alog10(AA(ja+1))-alog10(AA(ja)))
      zfrac = (alog10(zin) - alog10(ZZ(jz)))/
     &               (alog10(ZZ(jz+1))-alog10(ZZ(jz)))
      do iz = 1,2
         do j=1,14000
            if (high2(j,1,iz).gt.0.0 .and. high2(j,2,iz).gt.0.0) then
               x1 = (1.0-afrac)*alog10(high2(j,1,iz)) + 
     &                    afrac*alog10(high2(j,2,iz))
               high2(j,1,iz) = 10.0**x1
            else
               high2(j,1,iz) = (1.0-afrac)*high2(j,1,iz) +
     &                               afrac*high2(j,2,iz)
               if (high2(j,1,iz) .lt. 0.0) then
                  print*, 'Oooo, big high-res boo boo !'
                  high2(j,1,iz) = 0.0
               end if
            end if
         end do
         do j=1,1221
            if (hlam2(j,1,iz).gt.0.0 .and. hlam2(j,2,iz).gt.0.0) then
               x1 = (1.0-afrac)*alog10(hlam2(j,1,iz)) + 
     &                    afrac*alog10(hlam2(j,2,iz))
               hlam2(j,1,iz) = 10.0**x1
            else
               hlam2(j,1,iz) = (1.0-afrac)*hlam2(j,1,iz) +
     &                               afrac*hlam2(j,2,iz)
               if (hlam2(j,1,iz) .lt. 0.0) then
                  print*, 'Oooo, big boo boo !'
                  hlam2(j,1,iz) = 0.0
               end if
            end if
         end do
         bolo2(1,iz) = (1.0-afrac)*bolo2(1,iz) + afrac*bolo2(2,iz)
         do j=1,nmag
            mags2(j,1,iz) = (1.0-afrac)*mags2(j,1,iz) + 
     &                      afrac*mags2(j,2,iz)
            Lsq2(j,1,iz) = (1.0-afrac)*Lsq2(j,1,iz) + 
     &                      afrac*Lsq2(j,2,iz)
            Lavg2(j,1,iz) = (1.0-afrac)*Lavg2(j,1,iz) + 
     &                      afrac*Lavg2(j,2,iz)
            mtol2(j,1,iz) = (1.0-afrac)*mtol2(j,1,iz) + 
     &                      afrac*mtol2(j,2,iz)
         end do
         do j=1,nind
            Fcont2(j,1,iz) = (1.0-afrac)*Fcont2(j,1,iz) + 
     &                        afrac*Fcont2(j,2,iz)
            Fline2(j,1,iz) = (1.0-afrac)*Fline2(j,1,iz) + 
     &                        afrac*Fline2(j,2,iz)
         end do
      end do

C------------
C     now, the final interpolation; in log Z
      do j=1,14000
            if (high2(j,1,1).gt.0.0 .and. high2(j,1,2).gt.0.0) then
               x1 = (1.0-zfrac)*alog10(high2(j,1,1)) + 
     &                    zfrac*alog10(high2(j,1,2))
               high(j) = 10.0**x1
            else
               high(j) = (1.0-zfrac)*high2(j,1,1) + 
     &                        zfrac*high2(j,1,2)
            end if
      end do
      do j=1,1221
            if (hlam2(j,1,1).gt.0.0 .and. hlam2(j,1,2).gt.0.0) then
               x1 = (1.0-zfrac)*alog10(hlam2(j,1,1)) + 
     &                    zfrac*alog10(hlam2(j,1,2))
               hlam(j) = 10.0**x1
            else
               hlam(j) = (1.0-zfrac)*hlam2(j,1,1) + 
     &                        zfrac*hlam2(j,1,2)
            end if
      end do
      bolo = (1.0-zfrac)*bolo2(1,1) + zfrac*bolo2(1,2)
      do j=1,nmag
         mags(j) = (1.0-zfrac)*mags2(j,1,1) + 
     &                      zfrac*mags2(j,1,2)
         Lsq(j) = (1.0-zfrac)*Lsq2(j,1,1) + 
     &                      zfrac*Lsq2(j,1,2)
         Lavg(j) = (1.0-zfrac)*Lavg2(j,1,1) + 
     &                      zfrac*Lavg2(j,1,2)
         mtol(j) = (1.0-zfrac)*mtol2(j,1,1) + 
     &                      zfrac*mtol2(j,1,2)
      end do
      do j=1,nind
         Fcont(j) = (1.0-zfrac)*Fcont2(j,1,1) + 
     &                         zfrac*Fcont2(j,1,2)
         Fline(j) = (1.0-zfrac)*Fline2(j,1,1) + 
     &                          zfrac*Fline2(j,1,2)
      end do

C-----make sure everything is ready to return to the main program
      z = zin
      age = agein
      y = 0.228 + 2.7*zin

      return
      end
