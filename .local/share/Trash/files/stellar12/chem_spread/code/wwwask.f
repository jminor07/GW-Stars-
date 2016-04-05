      program wwwask
      integer     histpar
      parameter   ( histpar = 49 )
      real        ages(histpar),metals(histpar),shist(histpar)
      integer     status, npop
      character dirname*80, sline*40
      integer     iwidth
      real        zshift, age

C     interpret user input (stored in the temporary file 'tmp_input')
C     example: "12 0. 100. 
C              "DISP
C              "Wo__vanilla
C            : this means the user wants a 12 Gyr, [Fe/H]=0, pop with
C            : a relative weight of 100. 
C            : Then there may be one or two carriage returns.
C            : DISP is "display results" 
C            : SEND adds the low-res flux
C            : HIGH is low-res flux followed by hi-res flux
C            : 
C            : Finally, the "Wo__vanilla" string is a directory name


      npop = 49
      p = 0.00019
      a = 0.004
      Z = 0.02
      status = 3
      write(dirname,'(a11)') 'Wo__vanilla'

C     print*, 'Please enter the population age in Gyr'      
C     read*, age
C     print*, 'Setting age to ',age,' Gyr.'

      print*, 'Enter 1 for narrow, 2 for normal, 3 for wide,',
     z' 4 for Simple, 5 for Delta'
      read*, iwidth
      print*, 'Enter the shift in [Fe/H] '
      read*, zshift

      if ( zshift.lt.-1 .or. zshift.gt.0.3) then
         print*, ' Zshift is out of range. Stopping.'
         stop
      end if

      if ( iwidth .le. 0 .or. iwidth .ge. 6 ) then
         print*, 'No, that was not one of the choices!'
         stop
      end if

      print*, 'Looping over ages between 8 and 17 Gyr.'
         
         sum = 0.0
         do i=1,npop
            ages(i) = age
C           hardwired abundance scale running from -2.0 to 0.5 in 49 bins
            delta = 2.5/real(npop+1)
            f = -2.0 + delta*real(i)
            metals(i) = f
C           line to include user-requested shift in the abundance scale 
            f = f - zshift
            if( iwidth.eq.1 ) then
            b = (Z * 10.0**(1.5*f-0.0800762535152676) + a) / p
            e = exp( (- (Z*10.0**(1.5*f-0.0800762535152676) + 
     z                                      a)**2 - a**2 )/(2.0*p) )
            c = ( Z*10.0**(1.5*f-0.0800762535152676)) / 0.43429444819
            end if
            if( iwidth.eq.2 ) then
            b = (Z * 10.0**(f-0.0800762535152676) + a) / p
            e = exp( (- (Z * 10.0**(f-0.0800762535152676) +
     z                                       a)**2 - a**2 )/(2.0*p) )
            c = ( Z*10.0**(f-0.0800762535152676 )) / 0.43429444819
            end if
            if( iwidth.eq.3 ) then
            b = (Z * 10.0**(0.667*f-0.0800762535152672) + a) / p
            e = exp( (- (Z*10.0**(0.667*f-0.0800762535152676) + 
     z                                         a)**2 - a**2 )/(2.0*p) )
            c = ( Z*10.0**(0.667*f-0.0800762535152676 )) / 0.43429444819
            end if
            shist(i) = b*c*e
            sum = sum + shist(i)
         end do


C     normalize the integral to one
      do i=1,npop
         shist(i) = shist(i)/sum
      end do

C     start loop over ages
      do iage = 1,10
         age = 7.0 + real(iage)

         do i=1,npop
            ages(i) = age
            print*, i,ages(i),metals(i),shist(i)
         end do

C     Simple model
         if ( iwidth .eq. 4 ) then
            call simplemodel(zshift,histpar,metals,shist)
         end if
C     Delta function - put all in one bin
         if ( iwidth .eq. 5 ) then
            rmax = 100.
            do j=1,npop
               shist(j) = 0.0
               radius = abs( metals(j) - zshift )
               if ( radius .lt. rmax ) then
                  rmax = radius
                  jsuccess = j
                  rmin = radius
               end if
            end do
            shist(jsuccess) = 1.0
         end if

C     if things are OK, call subroutine
C      dirname = '/space/www/dial/export/'//dirname
         nchar = 34
         print*, dirname
         call ask (iwidth,zshift,metals,ages,shist,
     z             npop,status,dirname,nchar)

         print*, 'status = ',status
         do i=1,npop
            print*, i,ages(i),metals(i),shist(i)
         end do

C     end loop over ages
      end do

C     write status flag to alert cgi-perl script
      open(unit=34,file='tmp_flag',status='unknown')
      write(34,'(i2)') status
      close(34)

      stop
      end


      subroutine ask(iwidth,zshift,metalref,ageref,shist,
     z                 npop,status,dirname,nchar)
C     noninteractively 
C     interpolates in a model grid and combines populations.

C                  nind is the number of indices
C                  histpar is the max number of populations to combine
C                  npop is the number of populations to combine
C                  nclrs is the number of colors
      integer      nind,histpar,status,npop,nclrs
      parameter    (nind = 179, histpar = 49, nclrs=145)
      real         metalref(histpar),ageref(histpar),shist(histpar)
C     dirname is the directory for the model grid, nchar how many characters
C     in the directory string.
      character    dirname*80
      integer      nchar,iwidth
      real         zshift

C     primary passed variables
      real         hlam(1221),wave(1221), bolo, z, y, age
      real         mags(nclrs),Lsq(nclrs),Lavg(nclrs),mtol(nclrs)
      real         widths(nind),Fcont(nind),Fline(nind)
      integer      indexunits(nind)
      character*10 labels(nind)
      character*12 maglabels(nclrs)
      character*20 imftype
      real         massinit,massnow
C     other passed variables
      real         imfpar, zsolar, ysolar, yprimo, imflow, imfupp
      real         hbmorph, pagb, eta

C    
      real         zin, agein
C     variables for storing things.
      real         ahlam(1221,histpar), abolo(histpar)
      real         ahigh(14000,histpar),high(14000)
      real         amags(nclrs,histpar),amtol(nclrs,histpar)
      real         aLsq(nclrs,histpar),aLavg(nclrs,histpar)
      real         aFcont(nind,histpar),aFline(nind,histpar)
      real         amassinit(histpar),amassnow(histpar)
      real         aind(nind,histpar)
C     variables for output
      real         index(nind), Lbar(nclrs)

C     other variables
      character*24   ctime,t
      integer        time

      

C     normalize mass histogram to integral one
 22   xsum = 0.0
      print*, ' npop = ',npop
      do j=1,npop
         xsum = xsum + shist(j)
      end do
      do j=1,npop
         shist(j) = shist(j)/xsum
      end do

C     get and store the interpolated population models
      do i=1,npop
         zin = 0.0169*10**(metalref(i))
         agein = ageref(i)

         call getmodel(zin,agein,dirname,nchar,
     &          hlam,wave,bolo,z,y,age,mags,Lsq,Lavg,mtol,
     &       widths,Fcont,Fline,indexunits,labels,maglabels,
     &       massinit,massnow,imfpar,zsolar,
     &       ysolar,yprimo,imflow,imfupp,hbmorph,pagb,eta,imftype,
     &       status,high)

         if (status.eq.5) then
            open(unit=34,
     &        file='tmp_flag',status='unknown')
            write(34,'(i2)') status
            close(34)
            stop
         end if

C        store important results
         do j=1,14000
            ahigh(j,i) = high(j)
         end do
         do j=1,1221
            ahlam(j,i) = hlam(j)
         end do
         abolo(i) = bolo
         amassinit(i) = massinit
         amassnow(i) = massnow
c         print*, ' bolo, massnow = ',bolo,massnow
         do j=1,nclrs
            amags(j,i) = mags(j)
            aLsq(j,i) = Lsq(j)
            aLavg(j,i) = Lavg(j)
            amtol(j,i) = mtol(j)
         end do
c         do j=1,15
c            print*, ' mag, Lsq, Lavg = ',mags(j),Lsq(j),Lavg(j),mtol(j)
c         end do
         do j=1,nind
            aFcont(j,i) = Fcont(j)
            aFline(j,i) = Fline(j)
            if (indexunits(j).eq. 1) then
               aind(j,i) = widths(j)*(1.0-Fline(j)/Fcont(j))
            else
               aind(j,i) = -2.5*alog10(Fline(j)/Fcont(j))
            end if
         end do
C     end of loop over ingredient populations
      end do


C     prepare to sum by clearing arrays
      do j=1,14000
         high(j) = 0.0
      end do
      do j=1,1221
         hlam(j) = 0.0
      end do
      bolo = 0.0
      massinit = 0.0
      massnow = 0.0
      do j=1,nclrs
         mags(j) = 0
         Lsq(j) = 0.0
         Lavg(j) = 0.0
C        set mtol equal to the set of solar absolute mags
         mtol(j) = amags(j,1) + 2.5*alog10(amassnow(1)/amtol(j,1))
c         print*, j,maglabels(j),mtol(j)
      end do
      do j=1,nind
         Fcont(j) = 0.0
         Fline(j) = 0.0
      end do

C     loop over histogram bins and sum
      do iz = 1,npop
C        if there are few stars in a bin, skip it.
         if (shist(iz).lt.0.0001) goto 55
         do j=1,14000
            high(j) = high(j) + shist(iz)*ahigh(j,iz)
         end do
         do j=1,1221
            hlam(j) = hlam(j) + shist(iz)*ahlam(j,iz)
         end do
C        note that we go to non-log units
c         print*, iz,bolo, shist(iz),abolo(iz)
         bolo = bolo + shist(iz)*10**(abolo(iz))
         massinit = massinit + shist(iz)*amassinit(iz)
         massnow = massnow + shist(iz)*amassnow(iz)
         do j=1,nclrs
            mags(j) = mags(j) +shist(iz)*10**(-0.4*amags(j,iz))
            Lsq(j) = Lsq(j) + shist(iz)*aLsq(j,iz)
            Lavg(j) = Lavg(j) + shist(iz)*aLavg(j,iz)
         end do
         do j=1,nind
            Fcont(j) = Fcont(j) + shist(iz)*aFcont(j,iz)
            Fline(j) = Fline(j) + shist(iz)*aFline(j,iz)
         end do
 55      continue
      end do

C     open some output files
      open(unit=45,file='wwwask.output',status='unknown',
     z access='append')
      t = ctime(time())
      write(45,'(a,a24)') ' W94 model interpolation output ',t
      write(45,'(x)')
      write(45,'(a)') 'Composite population comprised of: '
      write(45,'(a)') ' Age (Gyr)    [Fe/H]    weight'
      do i=1,npop
         write(45,'(2f10.3,f10.5)') ageref(i), metalref(i), shist(i)
      end do
      write(45,'(x)')
      write(45,'(x)')

C     form indices and output them to a file
      do k=1,nind
         if (indexunits(k).eq. 1) then
            index(k) = widths(k)*(1.0-Fline(k)/Fcont(k))
         else
            index(k) = -2.5*alog10(Fline(k)/Fcont(k))
         end if
      end do

      open(unit=56,file='wwwask.200',status='unknown',access='append')
      write(56,'(a,79(a10,1x))') ' age width  zshift ',
     z                        (labels(k),k=101,179)
      write(56,'(f7.1,i2,f7.2,79(2x,f8.4,1x))') 
     z       ageref(1),iwidth,zshift,(index(k),k=101,179)
      close(56)

      write(45,'(a)') 'Spectral Indices (Lick system)'
      write(45,'(7(a10,1x))') (labels(k),k=1,7)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=1,7)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=8,14)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=8,14)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=15,21)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=15,21)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=22,25)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=22,25)
      write(45,'(x)')
      write(45,'(a)') 'Spectral Indices (200 km/s system)'
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=101,107)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=101,107)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=108,114)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=108,114)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=115,121)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=115,121)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=122,128)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=122,128)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=129,135)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=129,135)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=136,142)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=136,142)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=143,149)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=143,149)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=150,156)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=150,156)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=157,163)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=157,163)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=164,170)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=164,170)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=171,177)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=171,177)
      write(45,'(x)')
      write(45,'(7(a10,1x))') (labels(k),k=178,179)
      write(45,'(7(2x,f8.4,1x))') (index(k),k=178,179)
      write(45,'(x)')
      write(45,'(x)')
C     form Lbar and output

      do k=1,nclrs
         Lbar(k) = -2.5 * alog10( Lsq(k)/Lavg(k) )
      end do
      write(45,'(a)') ' SBF magnitudes '
      write(45,'(a,a,a)') 
     &      '   Ubar   Bbar   Vbar   Rbar   Ibar',
     &      '   Jbar   Hbar   Kbar   Lbar   Lpbar',
     &       '  Mbar   Kpbar'
      write(45,'(12f7.3)')  Lbar(1),
     z         (Lbar(k),k=3,13)
      write(45,'(x)')
      write(45,'(x)')

C     calculate and output mass-to-light ratios
C     divide the mass by the population luminosity in solar units
      do k=1,nclrs
         mags(k) = -2.5*alog10(mags(k))
         mtol(k) = massnow / (10.0**(0.4*(mtol(k)-mags(k))))
      end do
      write(45,'(a)') ' Mass-to-light ratios '
      write(45,'(a,a,a)') 
     &      '  M/L_U  M/L_B  M/L_V  M/L_Rc M/L_Ic',
     &       ' M/L_J  M/L_H  M/L_K  M/L_L  M/L_Lp',
     &       ' M/L_M  M/L_Kp'
      write(45,'(12f7.2)') mtol(1),(mtol(k),k=3,13)
      write(45,'(x)')
      write(45,'(x)')

C     output magnitudes and colors
      write(45,'(a)') ' Ordinary magnitudes '
      write(45,'(a,a,a)') 
     &      '      U      B      V      Rc     Ic',
     &       '     J      H      K      L      Lp',
     &       '     M      Kp'
      write(45,'(12f7.2)') mags(1),(mags(k),k=3,13)
      write(45,'(x)')
      write(45,'(x)')
      write(45,'(a)') 'Some colors'
      write(45,'(a,a)') 
     &      '     U-V     B-V     V-Rc    V-Ic    V-K',
     &      '     J-K     H-K'
      write(45,'(7f8.3)') mags(1)-mags(4),mags(3)-mags(4),
     &  mags(4)-mags(5),mags(4)-mags(6),mags(4)-mags(9),
     &  mags(7)-mags(9),mags(8)-mags(9)
      write(45,'(x)')
      write(45,'(a,f6.3,a)') 'B-V =',mags(3)-mags(4)-0.05,
     &      ' might be more accurate (see W94).'
      write(45,'(x)')
      write(45,'(x)')
      write(45,'(a,f9.0)') 'Initial population mass   : ',massinit
      write(45,'(a,f9.0)') 'Mass in stars and remnants: ',massnow
      write(45,'(a,f7.3)') 'log Lbol/Lbol_sun         : ',alog10(bolo)
      if (imftype(2:4).eq.'SAL' .or. imftype(2:4).eq.'sal') then
         write(45,'(a,f7.3)') 'IMF power-law slope       : ',imfpar
      else if (imftype(2:4).eq.'MIL'.or.imftype(2:4).eq.'mil') then
         write(45,'(a)') 'We have a Miller-Scalo (1979) piecewise IMF'
      else
         write(45,'(a)') 'I do not recognise this IMF type!'
      end if
      write(45,'(a,f7.3)') 'IMF lower mass cut        : ',imflow
      write(45,'(a,f7.3)') 'IMF upper mass cut        : ',imfupp
      write(45,'(x)')

C     complete listing of magnitudes
      write(45,'(a)') 'Complete Photometry Rundown:'
      write(45,'(x)')
      write(45,'(a)') 'Filter         Abs_Mag   SBF_Mag  M/L ratio'
      do i=1,nclrs
         write(45,'(a12,3f10.3)') maglabels(i),mags(i),Lbar(i),mtol(i)
      end do
      write(45,'(x)')

C     output SED
      if (status.eq.2 .or. status.eq.3) then
      write(45,'(a)')'Spectral energy distribution in wavelength units.'
      write(45,'(x)')
      write(45,'(a)')'Normalized so that the integral over wavelength'
      write(45,'(a)')'in nm is equal to the bolometric luminosity.'
      write(45,'(a)')'Interpolation error can cause up to 1% mismatch.'
      write(45,'(x)')
      write(45,'(a)')'Column 1: wavelength in Angstroms '
      write(45,'(a)')'Column 2: log wavelength (in Angstroms) '
      write(45,'(a)')'Column 3: Flux in wavelength units'
      write(45,'(a)')'Column 4: log of column 3'
      write(45,'(x)')
      write(45,'(x)')
      do j=1,1221
         if ( hlam(j).gt.0.0 ) then
            x1 = alog10(hlam(j))
         else
            x1 = 0.0
         end if
         write(45,'(4e14.5)') 10.0*wave(j),alog10(wave(j))+1.0,
     &        hlam(j),x1
      end do
      end if

C     high-res flux
      if (status.eq.3) then
         write(45,'(x)')
         write(45,'(x)')
         write(45,'(a)') 'Higher-resolution wavelength and flux'
         write(45,'(x)')
         do j=1,14000
            wv = 3000.0 + 0.5*real(j-1)
            write(45,'(f9.3,e14.5)') wv,high(j)
         end do
      end if

      close(45)


      open(unit=46,file='wwwask.hspectrum',status='unknown',
     z access='append')
      do j=1,14000
         wv = 3000.0 + 0.5*real(j-1)
         write(46,'(f9.3,e14.5)') wv,high(j)
      end do
      close(46)

      open(unit=47,file='wwwask.mags',status='unknown',access='append')
      write(47,'(150a9)') (maglabels(k),k=1,nclrs)
      write(47,'(150f9.3)') (mags(k),k=1,nclrs)
      close(47)

      open(unit=48,file='wwwask.lbars',status='unknown',access='append')
      write(48,'(150a9)') (maglabels(k),k=1,nclrs)
      write(48,'(150f9.3)') (Lbar(k),k=1,nclrs)
      close(48)

      return
      end



C----------------------------------------------------------------------------
      subroutine      simplemodel(yield,histpar,metals,shist)
C     subroutine      simplemodel
C                     returns histogram filled with Simple model predictions

      real        yield
      integer     histpar
      real        shist(histpar),metals(histpar)
C     local variables
      real        sum, p, zsol, xZP
      integer     j

C     Z solar      
      zsol = 0.0169
C     p is the yield (yield is in logarithmic form)
      p = zsol*10**yield

C     sample the simple models at 0.1-dex intervals and store result in shist
      do j=1,histpar
C        xZP is the metallicity in units of the yield
         xZP = zsol*(10**metals(j))/p
C        the magic number 2.3025851 is 1/log(e)
C        dM = d[Fe/H] * 1/log(e) * xZP * exp(xZP)
         shist(j) = 0.1*2.3025851*xZP*exp(-xZP)
      end do

C    normalize histogram to one
      sum = 0.0
      do j=1,histpar
         sum = sum + shist(j)
      end do
      do j=1,histpar
         shist(j) = shist(j)/sum
      end do
c      print*, 'For yield ',yield,'   the sample/bin ratio is ',sum

C and return
      return
      end
