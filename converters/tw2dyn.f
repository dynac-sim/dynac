       program tw2dyn
c TRACEWIN to DYNAC converter
c Version 6.0R1 21-Jul-2014
c Added reading TTF factors from  independent file for DTL
       implicit real*8 (a-h,o-z)
       parameter(ncards=59,ntwcrds=12,npar=20,maxcell1=3000)
       DATA VL,XMAT,RPEL,QST/2.99792458E8,938.27231,28.17938E-14,1./
c KLE describes dynac cards, TWKLE describes tracewin cards       
       character*8 KLE(ncards)     
       character*9 TWKLE(ntwcrds),twkley    
       DATA KLE/'GEBEAM','INPUT','RDBEAM','ETAC','DRIFT',
     *          'QUADRUPO','SEXTUPO','QUADSXT','SOLENO','SOQUAD',
     *          'BMAGNET','CAVMC','CAVSC','FIELD','HARM',
     *          'BUNCHER','RFQCL','NEWF','NREF','SCDYNAC',
     *          'SCDYNEL','SCPOS','TILT','TILZ','CHANGREF',
     *          'TOF','REJECT','ZROT','ALINER','ACCEPT',
     *          'EMIT','EMITGR','COMMENT','WRBEAM','ENVEL',
     *          'CHASE','RWFIELD','RANDALI','TWQA','EMIPRT',
     *          'MMODE','RFQPTQ','STRIPPER','STEER','ZONES',
     *          'PROFGR','SECORD','RASYN','FDRIFT','FSOLE',
     *          'EGUN','COMPRES','REFCOG','FPART','QUAELEC',
     *          'QUAFK','CAVNUM','EDFLEC','STOP'/
       DATA TWKLE/'DRIFT','QUAD','EDGE','BEND','DTL_CEL',
     *            'NCELLS','MULTIPOLE','FREQ','GAP','END',
     *            'FIELD_MAP','LATTICE'/
       COMMON/PRMTRS/param(npar)
       character*80 inarg,myarg(10),infiln,oufiln,text,fieldn,inphas
       character*80 ineo,inttf
       character*128 inline
       dimension pha(maxcell1),phb(maxcell1),dtleo(maxcell1)
       dimension dtlt(maxcell1),dtltp(maxcell1),dtltpp(maxcell1)
       logical ffound
       pi=4.*atan(1.)
       in=7
c in2 is LUN for file containing the phases, inphas is name of that file       
       in2=9
       inphas='OptimusPlus_Phases_wField.txt'
c in3 is LUN for file containing the Eo for DTL, ineo is name of that file
       in3=10
       ineo='ESS_DTL_Eo.dat'
c in4 is LUN for file containing the TTF for DTL, inttf is name of that file
       in4=11
       inttf='DTL_TTF.txt'
       iou=16
       WRITE(6,3101)
       DO
         call GET_COMMAND_ARGUMENT(narg,inarg)
         if(LEN_TRIM(inarg).eq.0) exit
         narg=narg+1
         myarg(narg)=TRIM(inarg)
       ENDDO
c       INPUT ARGUMENTS:
c ********************************************************************************
       ffound=.false.
       do i=2,narg
         text=myarg(i)
         if(text(1:1).ne.'-') then
c the input argument is the name of the input file
           if(i.eq.2) then
             infiln=myarg(i)
             write(6,*) 'Input file: ',infiln
             ffound=.true.
           endif
           if(i.eq.3) then
             oufiln=myarg(i)
             write(6,*) 'Output file: ',oufiln
             ffound=.true.
           endif
         else
           if(myarg(i).eq.'-h') then
c     print out of help message
             write(6,*) 'Command format:'
             write(6,*) 'convert file1 file2 [-h]'
             write(6,*) 'where file1 is the TRACEWIN file (input ',
     *                  'to the convert program) '
             write(6,*) 'where file1 is the DYNAC file (output ',
     *                  'from the convert program) '
             write(6,*) 'Optional arguments:'
             write(6,*) '-h will list the argument options (this ',
     *                  'list)'
             stop
           endif
         endif
       enddo
3101   FORMAT('******** TW2DYN V6.0R1 (Beta), 21-Jul-2014 *********')
c
c read the phases
       open(in2,file=inphas,status='unknown')
       nphas=0
       DO
         READ(in2,'(A)',ERR=777,END=800,IOSTAT=N) inline
         nphas=nphas+1
         READ(inline,*) idum,zzz,pha(nphas),phb(nphas)
       ENDDO
800    close(in2)
c
c read the TTF values for the DTL
       open(in4,file=inttf,status='unknown')
       nttf=0
       DO
         READ(in4,'(A)',ERR=789,END=850,IOSTAT=N) inline
         nttf=nttf+1
         READ(inline,*) idum,dtlt(nttf),dtltp(nttf),dtltpp(nttf)
       ENDDO
850    close(in3)
c read the Eo values for the DTL
       open(in3,file=ineo,status='unknown')
       neo=0
       DO
          READ(in3,'(A)',ERR=799,END=900,IOSTAT=N) inline
          neo=neo+1
          READ(inline,*) idum,dtleo(neo)
       ENDDO
900    close(in3)
c
       open(in,file=infiln,status='unknown')
       open(iou,file=oufiln,status='unknown')
c
c initialize some parameters
c counter for number of DTL cells
       ndtlcl=0
c counter for number of electric (cavity) fields       
       ncavf=0
       freq=0.
       ncav=0
       neo=0
       DO
         READ(in,'(A)',ERR=888,END=999,IOSTAT=N) inline
         IF(inline(1:1).eq.';') then
c comment line
           write(iou,'(A)') inline
         ELSE
c key line
           ktw=0
           twkley=''
           DO k=1,ntwcrds
             twkley=twkle(k)
             IF(twkley(1:3).eq.inline(1:3)) then
               ktw=k
             ENDIF
           ENDDO
           if(ktw.eq.1) then
c TRACEWIN DRIFT
             READ(inline(6:128),*) param(1)
             write(iou,'(A5)') KLE(5)             
             write(iou,*) param(1)/10.
           endif           
           if(ktw.eq.2) then
c TRACEWIN QUAD
             READ(inline(5:128),*) param(1),param(2),param(3)
             write(iou,'(A8)') KLE(6)
             param(3)=param(3)/20.
             param(2)=0.1*param(2)*param(3)
             write(iou,*) param(1)/10.,param(2),param(3)
           endif           
           if(ktw.eq.3) then
c TRACEWIN EDGE (TBD)
             READ(inline(5:128),*) param(1),param(2),param(3)
c             write(iou,*) KLE(11)
             param(3)=param(3)/10.
             param(2)=0.1*param(2)*param(3)
             write(iou,*) 'ERROR: EDGE NOT YET IMPLEMENTED'
           endif           
           if(ktw.eq.4) then
c TRACEWIN BEND (TBD)
             READ(inline(5:128),*) param(1),param(2),param(3)
c             write(iou,*) KLE(11)
             param(3)=param(3)/10.
             param(2)=0.1*param(2)*param(3)
             write(iou,*) 'ERROR: BEND NOT YET IMPLEMENTED'
           endif           
           if(ktw.eq.5) then
c TRACEWIN DTL_CEL
             if(freq.eq.0.) then
               write(6,*) 'ERROR: DTL frequency not defined at input'
               STOP
             endif  
             ndtlcl=ndtlcl+1
             READ(inline(8:128),*) (param(i),i=1,14)
c cell length
             cl=param(1)/10.
c length first 1/2 quad
             ql1=param(2)/10.
c length second 1/2 quad
             ql2=param(3)/10.
c cell center
             shift=param(4)/10.
c gradient first 1/2 quad (need to convert to field)
             qb1=param(5)
c gradient second 1/2 quad (need to convert to field)
             qb2=param(6)
c EoTL
             ef=param(7)
c PHIs
             phis=param(8)
c aperture
             apt=param(9)/10.
c abs or rel phase
             ipar=int(param(10))
c BETAr
             br=param(11)
c TTF
             ttf=param(12)
c TTF'
             ttfp=param(13)/(-2.*pi)
c TTF''
             ttfpp=param(14)/(-4.*pi*pi)
c*******************************************             
c now write output data: 1/2*Q1,GAP,1/2*Q2 *
c*******************************************
             if(ef.ne.0.) then
             neo=neo+1
c*case1*start EoTL.ne.0
c               ef=ef/(ttf*cl*10000.)
                ef=dtleo(neo)
c calculate k required for deriving T' and T'' from kT' and kT''               
c               fk=freq*2000000.*pi/(br*vl)
               fk=1.
               ttfp=ttfp/fk
               ttfpp=ttfpp/(fk*fk)
               ttf=dtlt(neo)
               ttfp=dtltp(neo)
               ttfpp=dtltpp(neo)
               write(6,*) 'check',ndtlcl,neo,ef,ttf,ttfp,ttfpp
c first half quad (if field.ne.0)
               if(qb1.ne.0.) then
                 write(iou,'(A8)') KLE(6)
                 qb1=0.1*qb1*apt/2.
                 write(iou,*) ql1,qb1,apt/2.
c drift back over the length of the first half quad
                 write(iou,'(A5)') KLE(5)
                 write(iou,*) -ql1
               else
c drift back over the length of the "g" parameter as listed in TRACEWIN
c                 write(iou,'(A5)') KLE(5)
c                 write(iou,*) -shift
               endif
c gap             
               write(iou,'(A5)') KLE(13)
               write(iou,*) ndtlcl,' 0. ',br,cl,ttf,ttfp,
     *          ' 0. 0. 0. 0.',ef,phis,' 0. ',ttfpp,freq,' 1.'
               if(qb2.ne.0.) then
c drift back over the length of the second half quad
                 write(iou,'(A5)') KLE(5)
                 write(iou,*) -ql2
c second half quad (if field.ne.0)       
                 write(iou,'(A8)') KLE(6)
                 qb2=0.1*qb2*apt/2.
                 write(iou,*) ql2,qb2,apt/2.
               else
c drift over the length of the "g" parameter as listed in TRACEWIN
c                 write(iou,'(A5)') KLE(5)
c                 write(iou,*) shift
               endif
c*case1*end EoTL.ne.0
             else 
c*case2*start EoTL.eq.0 
c first half quad (if field.ne.0)
               write(iou,'(A7,I4,A10)')';CAVSC ',ndtlcl,' has Eo=0.'
               if(qb1.ne.0.) then
                 write(iou,'(A8)') KLE(6)
                 qb1=0.1*qb1*apt/2.
                 write(iou,*) ql1,qb1,apt/2.
c drift back over the length of the first half quad
                 write(iou,'(A5)') KLE(5)
                 write(iou,*) -ql1
                 write(iou,'(A5)') KLE(5)
                 write(iou,*) cl
               else
c drift back over the length of the "g" parameter as listed in TRACEWIN
c                 write(iou,'(A5)') KLE(5)
c                 write(iou,*) -shift
               endif
               if(qb2.ne.0.) then
               write(iou,'(A5)') KLE(5)
               write(iou,*) cl
c drift back over the length of the second half quad
                 write(iou,'(A5)') KLE(5)
                 write(iou,*) -ql2
c second half quad (if field.ne.0)       
                 write(iou,'(A8)') KLE(6)
                 qb2=0.1*qb2*apt/2.
                 write(iou,*) ql2,qb2,apt/2.
               endif
cc*case2*end EoTL.eq.0
             endif
           endif           
           if(ktw.eq.8) then
c TRACEWIN FREQ
             READ(inline(5:128),*) param(1)
c frequency in MHz
             freq=param(1)
             write(iou,'(A4)') KLE(18)
             fh=param(1)*1000000.
             write(iou,*) fh
           endif           
           if(ktw.eq.9) then
c TRACEWIN GAP
             READ(inline(4:128),*) (param(i),i=1,10)
             write(iou,'(A7)') KLE(16)
             etl=param(1)/1000000.
             phi=param(2)
             apt=param(3)/10.
             write(iou,*) etl,phi,' 1 ',apt
           endif           
           if(ktw.eq.10) then
c TRACEWIN END
             write(iou,'(A4)') KLE(59)
           endif           
           if(ktw.eq.11) then
c TRACEWIN FIELD_MAP
             READ(inline(10:128),*) (param(i),i=1,8),fieldn
             ncavf=ncavf+1
             ncav=ncav+1
             if(ncavf.eq.1) then
               ilen=LEN_TRIM(fieldn)
               fieldn(ilen+1:ilen+4)='.txt'
               write(iou,'(A5)') KLE(14)
               write(iou,'(A)') fieldn(1:ilen+4)
c assume field file is in MV/m; convert to V/m               
               write(iou,'(A)') ' 1000000.'               
             endif  
c             
             write(iou,'(A6)') KLE(57)
             write(iou,*) ncavf
             geom=param(1)
             flen=param(2)/10.
             phin=param(3)
             apt=param(4)/10.
             bamp=param(5)
             eamp=(param(6)-1.)*100.
             sscf=param(7)
             aptf=param(8)
             fdum=0.
c use phase from file 
             phad=pha(ncav)
             write(iou,*) fdum,phad,eamp,' 8 1'
           endif           
           if(ktw.eq.12) then
c TRACEWIN LATTICE
c initialize some parameters
c counter for number of DTL cells
             ndtlcl=0
c counter for number of electric (cavity) fields       
             ncavf=0
           endif           
         ENDIF
       ENDDO       
999    WRITE(6,*) 'All read'
       GOTO 1000
777    WRITE(6,*) 'I/O error # ', N, ', on ',in2
       STOP     
789    WRITE(6,*) 'I/O error # ', N, ', on ',in4
       STOP
799    WRITE(6,*) 'I/O error # ', N, ', on ',in3
       STOP
888    WRITE(6,*) 'I/O error # ', N, ', on ',in
       STOP
1000   CONTINUE
       END         
c#define DRIFT 1
c#define QUAD 2
c#define EDGE 3
c#define BEND 4
c#define DTL_CEL 5
c#define NCELLS 6
c#define MULTIPOLE 7
c#define FREQ 8
c#define GAP 9
c#define END 100
