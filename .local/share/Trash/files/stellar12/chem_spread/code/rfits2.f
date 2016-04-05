C     subroutine rfits2
      subroutine rfits2(infile,high,hlam,wave,bolo,z,y,age,mags,
     &   Lsq,Lavg,mtol,
     &   widths,Fcont,Fline,indexunits,labels,massinit,massnow,imfpar,
     &   zsolar,ysolar,yprimo,imflow,imfupp,hbmorph,pagb,eta,imftype,
     &   maglabels,status)

C     a subroutine to read selected parts of a diag-output FITS file
C     Uses the FITSIO routines of William Pence

      integer      nind,nmag
      parameter    (nind=179,nmag=145)
C     input variables
      character*80 infile
C     primary passed variables
      real         high(14000)
      real         hlam(1221),wave(1221), bolo, z, y, age
      real         mags(nmag),Lsq(nmag),Lavg(nmag),mtol(nmag)
      real         widths(nind),Fcont(nind),Fline(nind)
      integer      indexunits(nind)
      character*10 labels(nind)
      character*12 maglabels(nmag)
      character*20 imftype
      real         massinit,massnow
C     other passed variables
      real         imfpar, zsolar, ysolar, yprimo, imflow, imfupp
      real         hbmorph, pagb, eta
C     local variables
      character*80 comment
C     variables for reading in header records
      integer iunit, status, block, htype
C     htype indicates the kind of HDU: 0=primary binary, 1=ASCII table,
C     and 2=binary table. We don't deal with type 2 in this subroutine.
      character*30 errtxt
C     variables for reading in data arrays and header keywords
      logical      simple, extend
      integer      bitpix, naxis, naxes(99), pcount, gcount, group
      integer      fpixel 
      logical      anyf
C     for ascii tables.  "20" is a size large enough for this subroutine
      integer      rowlen, nrows, tfield, tbcol(20)
      character*12  ttype(20), tform(20), tunit(20)
      character*16 extnam


      print*, ' sub rfits2 ok '

      status = 0
      iunit = 15

C     open FITS file
      call ftopen(iunit,infile,0,block,status)
C     error check
      if (status.gt.0) goto 444

C     kk is the HDU number
      kk = 1
C     go to HDU number kk
      call ftmahd(iunit,kk,htype,status)
      if (status .gt. 0 ) goto 444
C     get primary keywords
      call ftghpr(iunit,99,simple,bitpix,naxis,naxes,pcount,
     &                  gcount,extend,status)
C     get some keyword values
      call ftgkey(iunit,'imftype',imftype,comment,status)
      call ftgkye(iunit,'loglbol',bolo,   comment,status)
      call ftgkye(iunit,'z',      z,      comment,status)
      call ftgkye(iunit,'y',      y,      comment,status)
      call ftgkye(iunit,'age',    age,    comment,status)
      call ftgkye(iunit,'massinit',massinit,comment,status)
      call ftgkye(iunit,'massnow',massnow,comment,status)
      call ftgkye(iunit,'imfpar', imfpar, comment,status)
      call ftgkye(iunit,'zsolar', zsolar, comment,status)
      call ftgkye(iunit,'ysolar', ysolar, comment,status)
      call ftgkye(iunit,'yprimo', yprimo, comment,status)
      call ftgkye(iunit,'imflow', imflow, comment,status)
      call ftgkye(iunit,'imfupp', imfupp, comment,status)
      call ftgkye(iunit,'hbmorph',hbmorph,comment,status)
      call ftgkye(iunit,'pagb',   pagb,   comment,status)
      call ftgkye(iunit,'eta',    eta,    comment,status)

C     read first and third rows of the image (wavelength scale and flux)
      group = 0
      fpixel = 1
      call ftgpve(iunit,group,fpixel,naxes(1),0,wave,anyf,status)
      fpixel = 1 + 2*naxes(1)
      call ftgpve(iunit,group,fpixel,naxes(1),0,hlam,anyf,status)
      if (status .gt. 0) goto 444

C        -----------------------------------------------------------
C     go to HDU number kk
      kk = 2
      call ftmahd(iunit,kk,htype,status)
      if (status .gt. 0 ) goto 444

C     read required ascii table header keywords
      call ftghtb(iunit,20,rowlen,nrows,tfield,ttype,
     &                  tbcol,tform,tunit,extnam,status)
            
C     read columns for Absolute magnitudes, Lsq, and Lavg
      call ftgcvs(iunit,1,1,1,nmag,' ',maglabels,anyf,status)
      call ftgcve(iunit,2,1,1,nmag,0,mags,anyf,status)
      call ftgcve(iunit,3,1,1,nmag,0,mtol,anyf,status)
      call ftgcve(iunit,4,1,1,nmag,0,Lsq,anyf,status)
      call ftgcve(iunit,5,1,1,nmag,0,Lavg,anyf,status)

      print*, ' after mags '

C     go to HDU number kk
      kk = 3
      call ftmahd(iunit,kk,htype,status)
      print*, ' debug 0 '
      if (status .gt. 0 ) goto 444
C     read required ascii table header keywords
      print*, ' debug 1 '
      call ftghtb(iunit,20,rowlen,nrows,tfield,ttype,
     &                  tbcol,tform,tunit,extnam,status)
      print*, ' debug 2 '
C     read columns for index name, index width, Continuum flux, line flux,
C     and the 0/1 code for angstroms or magnitude units
      print*, ' before labels '
      call ftgcvs(iunit,1,1,1,nind,' ',labels,anyf,status)
      print*, ' after labels '
      call ftgcve(iunit,2,1,1,nind,0,widths,anyf,status)
      call ftgcve(iunit,4,1,1,nind,0,Fcont,anyf,status)
      call ftgcve(iunit,5,1,1,nind,0,Fline,anyf,status)
      call ftgcvj(iunit,6,1,1,nind,0,indexunits,anyf,status)

      print*, ' before high-res fluxes '

C     go to HDU number kk
C     5 = "image" extension of 0.5-res fluxes
      kk = 5
      call ftmahd(iunit,kk,htype,status)
      if (status .gt. 0 ) goto 444
      call ftghpr(iunit,maxdim,simple,bitpix,naxis,naxes,pcount,
     z            gcount,extend,status)
c      call ftghbn(iunit,maxdim,nrows,tfields,ttype,tform,tunit,
c     z       extname,varidat,status)

      print*, 'maxdim = ',maxdim
      print*, 'naxis = ',naxis
      print*, 'naxes(1) = ',naxes(1)
      print*, 'group = ',group

      call ftgpve(iunit,group,1,14000,1,high,anyf,status)
c      call ftgcve(iunit,2,1,1,14000,0,high,anyf,status)

      if (status .gt. 0) goto 444

C     close files
      call ftclos(iunit,status)


      return

 444  print*, 'Error. Subroutine rfits2 Aborts.'
      call ftgerr(status,errtxt)
      print*, status,':  ',errtxt
      call ftclos(iunit,status)
      close(16)

      return

      end
