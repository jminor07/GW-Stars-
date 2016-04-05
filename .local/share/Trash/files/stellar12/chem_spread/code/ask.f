      program ask
C     Interactively inquires of user which model he would like, and
C     interpolates in a model grid to supply this.

C                  histpar is the number of bins in the [Fe/H] histogram
C                  nind is the number of indices
      integer      nind,histpar
      parameter    (nind = 25, histpar = 28)

C     primary passed variables
      real         hlam(1221),wave(1221), bolo, z, y, age
      real         mags(11),Lsq(11),Lavg(11)
      real         widths(nind),Fcont(nind),Fline(nind)
      integer      indexunits(nind)
      character*10 labels(nind)
      real         mass
C     other passed variables
      real         imfpar, zsolar, ysolar, yprimo, imflow, imfupp
      real         hbmorph, pagb, eta

C    
      real         shist(histpar), zin, agein
      character*25 optdir
C     variables for storing things.
      real         ahlam(1221,histpar), abolo(histpar)
      real         amags(11,histpar)
      real         aLsq(11,histpar),aLavg(11,histpar)
      real         aFcont(nind,histpar),aFline(nind,histpar)
      real         amass(histpar),aind(nind,histpar)
      real         metalref(histpar),ageref(histpar)
C     variables for output
      real         index(nind), Lbar(11)

C     other variables
      character*80   text
      character*24   ctime,t
      integer        time,status

      print*, ' '
      print*, ' '
      print*, 'ASK returns a flux, colors, SBF mags, and indices'
      print*, 'for combinations of single-burst populations. '
      print*, ' '
      print*, 'You can enter up to ',histpar,
     &        ' separate populations, along '
      print*, 'with weights (which the program will normalize so that '
      print*, 'the total population mass is a million Mo). '
      print*, 'Metallicities are restricted to -2 < [Fe/H] < +0.5, and '
      print*, 'ages to between 1 and 18 Gyr. In addition, the age<8, '
      print*, '[Fe/H]<-0.225 corner is not allowed. '
      print*, ' '
      print*, 'Enter trios of numbers separated by spaces: '
      print*, '( age  [Fe/H]  weight ). '
      print*, 'Enter a blank to exit the loop '
      print*, ' ' 

      do i=1,histpar
 11      print*, 'Enter (age [Fe/H] weight) no. ',i
         read(*,'(a80)') text
         if (text.eq.' ') then
            npop =  i-1
            goto 22
         end if
         read(text,*,err=99)  ageref(i),metalref(i),shist(i)
         goto 12
 99      print*, 'I had trouble reading that. Try again.'
         goto 11
 12   end do
      npop = 28

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
c        optdir is 25 characters long, gives subdirectory for model grid
         optdir = 'vanilla'

         call getmodel(zin,agein,optdir,
     &          hlam,wave,bolo,z,y,age,mags,Lsq,Lavg,
     &       widths,Fcont,Fline,indexunits,labels,mass,imfpar,zsolar,
     &       ysolar,yprimo,imflow,imfupp,hbmorph,pagb,eta,status)


C        store important results
         do j=1,1221
            ahlam(j,i) = hlam(j)
         end do
         abolo(i) = bolo
         amass(i) = mass
c         print*, ' bolo, mass = ',bolo,mass
         do j=1,11
            amags(j,i) = mags(j)
            aLsq(j,i) = Lsq(j)
            aLavg(j,i) = Lavg(j)
c            print*, ' mag, Lsq, Lavg = ',mags(j),Lsq(j),Lavg(j)
         end do
         do j=1,nind
            aFcont(j,i) = Fcont(j)
            aFline(j,i) = Fline(j)
            if (indexunits(j).eq. 1) then
               aind(j,i) = widths(j)*(1.0-Fline(j)/Fcont(j))
            else
               aind(j,i) = -2.5*alog10(Fline(j)/Fcont(j))
            end if
c            print*, 'index, Fline, Fcont = ',aind(j,i),
c     &                      Fline(j),Fcont(j)
         end do
C     end of loop over ingredient populations
      end do


C     prepare to sum by clearing arrays
      do j=1,1221
         hlam(j) = 0.0
      end do
      bolo = 0.0
      mass = 0.0
      do j=1,11
         mags(j) = 0
         Lsq(j) = 0.0
         Lavg(j) = 0.0
      end do
      do j=1,nind
         Fcont(j) = 0.0
         Fline(j) = 0.0
      end do

C     loop over histogram bins and sum
      do iz = 1,npop
C        if there are few stars in a bin, skip it.
         if (shist(iz).lt.0.0001) goto 55
         do j=1,1221
            hlam(j) = hlam(j) + shist(iz)*ahlam(j,iz)
         end do
C        note that we go to non-log units
c         print*, iz,bolo, shist(iz),abolo(iz)
         bolo = bolo + shist(iz)*10**(abolo(iz))
         mass = mass + shist(iz)*amass(iz)
         do j=1,11
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
      open(unit=45,file='ask_most.out',status='unknown')
      open(unit=48,file='ask_sed.out',status='unknown')
      t = ctime(time())
      write(45,'(a,a24)') ' W94 model interpolation output ',t
      write(45,'(x)')
      write(45,'(a)') 'Composite population comprised of: '
      write(45,'(a)') ' Age (Gyr)    [Fe/H]    weight'
      write(48,'(a,a24)') ' W94 model interpolation output ',t
      write(48,'(x)')
      write(48,'(a)') 'Composite population comprised of: '
      write(48,'(a)') ' Age (Gyr)    [Fe/H]    weight'
      do i=1,npop
         write(45,'(2f10.3,f10.5)') ageref(i), metalref(i), shist(i)
         write(48,'(2f10.3,f10.5)') ageref(i), metalref(i), shist(i)
      end do
      write(45,'(x)')
      write(45,'(x)')
      write(48,'(x)')
      write(48,'(x)')

c      print*, ' summed mass = ',mass
c      print*, ' Log Lbol = ',alog10(bolo)
C     form indices and output them to a file
      do k=1,nind
         if (indexunits(k).eq. 1) then
            index(k) = widths(k)*(1.0-Fline(k)/Fcont(k))
         else
            index(k) = -2.5*alog10(Fline(k)/Fcont(k))
         end if
      end do
      write(45,'(a)') 'Spectral Indices'
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
      write(45,'(x)')
C     form Lbar and output

      do k=1,11
         Lbar(k) = -2.5 * alog10( Lsq(k)/Lavg(k) )
      end do
      write(45,'(a)') ' SBF magnitudes '
      write(45,'(a,a,a)') 
     &      '   Ubar   Bbar   Vbar   Rbar   Ibar',
     &      '   Jbar   Hbar   Kbar   Lbar   Lpbar',
     &       '  Mbar'
      write(45,'(11f7.3)')  
     z         (Lbar(k),k=1,11)
      write(45,'(x)')
      write(45,'(x)')

C     output magnitudes and colors
      do k=1,11
         mags(k) = -2.5*alog10(mags(k))
      end do
      write(45,'(a)') ' Ordinary magnitudes '
      write(45,'(a,a,a)') 
     &      '       U       B       V       Rc      I',
     &      'c      J       H       K'
      write(45,'(11f8.3)') (mags(k),k=1,8)
      write(45,'(x)')
      write(45,'(a,a)') '       L       Lp','      M'
      write(45,'(11f8.3)') (mags(k),k=9,11)
      write(45,'(x)')
      write(45,'(x)')
      write(45,'(a)') 'Some colors'
      write(45,'(a,a)') 
     &      '     U-V     B-V     V-Rc    V-Ic    V-K',
     &      '     J-K     H-K'
      write(45,'(7f8.3)') mags(1)-mags(3),mags(2)-mags(3),
     &  mags(3)-mags(4),mags(3)-mags(5),mags(3)-mags(8),
     &  mags(6)-mags(8),mags(7)-mags(8)
      
      write(45,'(x)')
      write(45,'(x)')
      write(45,'(a,f9.0)') 'Total population mass: ',mass
      write(45,'(a,f7.3)') 'log Lbol/Lbol_sun    : ',alog10(bolo)
      write(45,'(a,f7.3)') 'IMF power-law slope  : ',imfpar
      write(45,'(a,f7.3)') 'IMF lower mass cut   : ',imflow
      write(45,'(a,f7.3)') 'IMF upper mass cut   : ',imfupp
      write(45,'(x)')


C     output SED
      write(48,'(a)')'Spectral energy distribution in wavelength units.'
      write(48,'(x)')
      write(48,'(a)')'Normalized so that the integral over wavelength'
      write(48,'(a)')'in nm is equal to the bolometric luminosity.'
      write(48,'(a)')'Interpolation error can cause up to 1% mismatch.'
      write(48,'(x)')
      write(48,'(a)')'Column 1: wavelength in Angstroms '
      write(48,'(a)')'Column 2: log wavelength (in Angstroms) '
      write(48,'(a)')'Column 3: Flux in wavelength units'
      write(48,'(a)')'Column 4: log of column 3'
      write(48,'(x)')
      write(48,'(x)')
      do j=1,1221
         if ( hlam(j).gt.0.0 ) then
            x1 = alog10(hlam(j))
         else
            x1 = 0.0
         end if
         write(48,'(4e14.5)') 10.0*wave(j),alog10(wave(j))+1.0,
     &        hlam(j),x1
      end do

C     let us do some checks, shall we? Does integral of flux = Lbol?
c      xsum = 0.0
c      do j=2,1221
c         xsum = xsum + 0.5*(hlam(j)+hlam(j-1))*(wave(j)-wave(j-1))
c      end do
c      print*, ' integral = ',xsum,'   bolo =',bolo

      print*, ' Bye!'

      close(45)
      close(46)
      close(47)
      close(48)

      stop 
      end
