       BLOCK DATA
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       COMMON/RADIA/TRT,RMOY,XINTF,CRAE
       data CRAE,XINTF /2.81793910e-13,.86967/
       DATA VL,XMAT,RPEL,QST/2.99792458E10,938.27231,28.17938E-14,1./
       DATA H /.040484004,.092121499,.138873510,.178145981,
     1  .207816048,.226283180,.232551553,.226283180,
     2  .207816048,.178145981,.138873510,.092121499,
     3  .040484004/
       DATA T /-.984183055,-.917598399,-.801578091,
     1  -.642349339,-.448492751,-.230458316,0.,
     2  .230458316,.448492751,.642349339,.801578091,
     3  .917598399,.984183055/
       DATA H1 /-.990575473,-.950675522,-.880239154,
     1          -.781514004,-.657671159,-.512690537,
     2          -.351231763,-.178484181,  0.,
     3           .178484181, .351231763,.512690537,
     4           .657671159, .781514004,.880239154,
     5           .950675522, .990575473/
       DATA T1 /.024148303,.055459529,.085036148,
     1          .111883847,.135136368,.154045761,
     2          .168004102,.176562705,.179446470,
     3          .176562705,.168004102,.154045761,
     4          .135136368,.111883847,.085036148,
     5          .055459529,.024148303/
       common/randu/ck(15),kmax
c       randu(j):contains the Chebitcheff coefficients Ck in table 20.
c       kmax is the total number of these coefficients(from 1)
       data (ck(j),j=1,7)/.98933556,-.68838689,.28191718,
     *       -.66389307E-01,.87406854E-02,-.59534602E-03,
     *       .16300617E-04/
cxx        data (ck(j),j=1,10)/.99927015,-.78398394,.44577741,
cxx     *        -.18252873,.51928922e-01,-.99821428e-02,
cxx     *        .1259644e-02,-.99460265e-04,.44427251e-05,
cxx     *        -.85528336e-07/
cxx        data (ck(j),j=1,12)/.99989097,-.79469098,.48376369,
cxx     *        -.23183581,.85333907e-01,-.23354050e-01,
cxx     *        .46284273e-02,-.64840101e-03,.62211306e-04,
cxx     *        -.38752189e-05,.14087289e-06,-.22650639e-08/
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
c   pi2=pi*pi, sqpi=(pi/2)**3/2 , pwtpi=(2*pi)**3/2 ,sqpi2=sqrt(pi/2.)
       data pi2,sqpi,pwtpi,sqpi2/9.86960440,1.96870124,
     *                         15.74960995,1.25331414/
       data kmax/7/
       END
c *******************************************************************
C      PROGRAM DYNAC
C      THIS SOFTWARE WAS ORIGINALLY PRODUCED BY CERN/PS, CEN/SACLAY
C
C      @version 6.0R15
C      @date 13-Dec-2015
C
C      @author P. LAPOSTOLLE  CONSULTANT (Paris, France)
C      @author E. TANKE       ESS (Lund, Sweden)
C      @author S. VALERO      CONSULTANT (Paris, France)
C
C      Modified and maintained by :
C       - TANKE Eugene
C       - VALERO Saby
c *******************************************************************
       PROGRAM dynac
       implicit real*8 (a-h,o-z)
       parameter(ncards=64,iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
       COMMON/DYNI/VREFI,TREFI,FHINIT
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/CARAC/CARA(10)
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/compt1/ndtl,ncavmc,ncavnm
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/faisc/f(10,iptsz),imax,ngood
       common/tapes/in,ifile,meta
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/TILT/TIPHA,TIX,TIY,SHIFW,SHIFP
       COMMON /BLOC21/ BE, APB(2), LAYL, LAYX, RABT
       real*8 LAYL, LAYX
       COMMON/PORO/IROT1,IROT2
       LOGICAL IROT1,IROT2
       common/etchas/fractx,fracty,fractl
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/davprt/shortl
       common/shortl/davprt
       COMMON/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/etcha3/ichxyz(iptsz)
       COMMON/DCSPA/IESP
       COMMON/TSTNT/intgct
       COMMON/TABSCH/IDEG,NN,PCHOIX
       COMMON/CMPTE/IELL
       common/cptemit/xltot(maxcell1),nbemit
       common/rfield/ifield
       common/posc/xpsc
       common/rander/ialin
       common/apel/iapel
       common/qskew/qtwist,iqrand,itwist,iaqu
       common/femt/iemgrw,iemqesg
       common/mode/eflvl,rflvl
       common/aerp/vphase,vfield,ierpf
       common/blvl/bflvl
       common/rec/irec
       COMMON/HISTO/CENTRE(6)
       COMMON/GRPARM/GLIM(4,2),GLIM1(4,2),GLIM2(4,2),PATITL,
     *        ngraphs(100),idwdp,igrprm,ngrafs
       common/zones/frms(6),nzone
       LOGICAL CHASIT,SHIFT,ITVOL,IMAMIN
       LOGICAL ICHAES,IESP,ifield,ialin,itwist,iemgrw
       character*80 cara,cmnt,text,patitl,ofeldf,ofelds
       character*80 davprt(maxcell1)
       character*8 KLE(ncards),kley
       character*12 fdrft
       character*160 titre
       character iitime*30
       COMMON/CESPCH/NCHGE,ICHSP,NPPI
       logical ichsp
       common/drfq/prfq(9)
       COMMON/ALIN/XL,YL,XPL,YPL
       COMMON/TRCMP/PCHOIXA
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/mcs/imcs,ncstat,cstat(20)
       COMMON/SECDR/ISEOR
       LOGICAL ISEOR
       COMMON/RAYSHY/IRAYSH
       common/trfq/icour,ncell
       common/newref/dephas,dewref,iref,irefw
       common/tofev/ttvols
       common/conti/irfqp
       logical irfqp
       common/rf1ptq/tvolt,avolt,fph,mlc,nceltot
       common/strip/atm,qs,atms,ths,qop,sqst(6),anp,nqst
       common/grot/rzot,izrot
c allow plotting the beam after sector nsprint ******
       common/isector/nsector,nsprint
       common/mingw/mg
       common/xposi/xpost(10),xlce(2),xpax(2),iscx(2)
       common/fcont/ifcont
       logical ifcont,mg,ffound
c********************************************
c reservation for TRACE3D related stuff
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       common/t3dfld/fldctr(15),zend(15),T3D
       logical T3D,xiset
       character*128 trace3h,trace3t,tif,tifa,tifb
c********************************************
       logical izrot
       LOGICAL IRAYSH
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
     *          'QUAFK','CAVNUM','EDFLEC','EMITL','RFKICK',
     *          'FIRORD','DCBEAM','T3D','STOP'/
       integer narg,i
       character*80 inarg,myarg(10),myfile,wfile,infiln,shortl
       character*10 plane
       logical g77,gfortran
c if mg=.true., use MINGW on windows, which has a different result for
c ctime function than standard gfortran
c default is mg=.false. This can be set by giving mingw as argument on
c the command line
        mg=.false.
c get arguments from the command line
c format for dynac:
c dynac file1 [-h]
c where: -file1 is the input file, describing the beamline
c        -h  print help info
        gfortran=.true.
c g77 in case of g77, set g77=.true. ; in case of gfortran,
c set g77=.false.
        g77=.false.
        if(g77) gfortran=.false.
        narg=0
        if (gfortran) then
          DO
            call GET_COMMAND_ARGUMENT(narg,inarg)
            if(LEN_TRIM(inarg).eq.0) exit
            narg=narg+1
c g77 in case of g77, comment out the next line; in case of
c gfortran, leave it in
             myarg(narg)=TRIM(inarg)
          ENDDO
        else
          write(6,*) 'Compatible with g77'
        endif
c       INPUT ARGUMENTS:
c *******************************************************************
        if (gfortran) then
          ffound=.false.
          do i=2,narg
            text=myarg(i)
            if(text(1:1).ne.'-') then
c the input argument is the name of the input file
              write(6,2917)myarg(i)
2917          format('Input file: ',A)
              infiln=myarg(i)
C 7 or 'in' is the unit corresponding to the dynac input file
c describing the accelerator or transport line
              open(7,file=myarg(i),status='unknown')
              ffound=.true.
            else
#ifdef MINGW
c if mg=.true., use MINGW on windows, which has a different result
c for ctime function than standard gfortran
                mg=.true.
                write(6,*) 'Using MINGW gfortran format on MSWindows'
#endif
              if(myarg(i).eq.'-h') then
c     print out of help message, starting with DYNAC version
                WRITE(6,3101)
                write(6,*) 'Command format:'
                write(6,*) 'dynacv6_0 [-h] [-mingw] [file1]'
                write(6,*) 'where file1 is the input file, describing ',
     *                     'the beamline'
                write(6,*) 'Optional arguments:'
                write(6,*) '-h will list the argument options (this ',
     *                     'list)'
              stop
              endif
            endif
          enddo
          if(.not.ffound) then
            write(6,*) 'Error: Input file name required'
            write(6,*) 'Type'
            write(6,*) 'dynacv6_0 -h'
            write(6,*) 'for syntax'
            stop
          endif
        endif
c *******************************************************************
c  OUTPUT FILES:
c
c 'dynac.long'  : print extensive computations information
c
c 'dynac.short' : print essentials of beam dynamics information
c
c 'dynac.print' : printout of the envelope and the emittance of the
c                 beam in the directions  x, y and z at each position
c                 of the optical and accelerating elements, as well as
c                 the number of good particles
c
c 'dynac.dmp'   : printout of cavity cell related data:
c                 cell number,synchronous phase (deg),relativistic
c                 beta (output), output energy (MeV),
c                 hor. emittance (mm.mrad, normalized),
c                 ver. emittance (mm.mrad, normalized),
c                 longitudinal emittance (ns.kev)
c
c 'beam_core.dst': print the coordinates of the particles kept by
c                  CHASE
c
c 'beam_remove.dst': print the coordinates of the particles
c                    removed by CHASE
c
c 'dynac_in_pr.dst': print of the coordinates of the particles at the
c                    input of the machine
c
c 'emit.plot'      : datafile used for the plots
c
c *************************************************************************************************
       open(16,file='dynac.long',status='unknown')
       open(12,file='dynac.short',status='unknown')
       open(71,file='dynac.print',status='unknown')
       open(50,file='dynac.dmp',status='unknown')
       open(61,file='beam_core.dst',status='unknown')
       open(60,file='beam_remove.dst',status='unknown')
       open(11,file='dynac_in_pr.dst',status='unknown')
       open(66,file='emit.plot',status='unknown')
c ****************************************************************************************************
c*******files deactivated in the code
c
c     'emlg.data': print at each position of space charge computation:
c                  the length(m), the longitudinal emittance growth (ns.keV),
c                  the kinetic energy (MeV), the number of good particles
c
c     'emtr.data': print at each position of space charge computation:
c                  the length(m), the hor. and ver. normalized emittances (mm.mrad)
c
c     'chemlg.data': as 'emlg.data' but with CHASE
c
c     'chemtr.data': as 'emtr.data' but with CHASE
c
c     'ntxyz.data' : print at each position of space charge computation:
c             the coordinates of the particles in the 3D plane (x(cm), y(cm) and z(cm))
c             This file is rewound each time, only the last result is kept
c
c     'rms.size' :print at each position of space charge computation the RMS(m) of the bunch
c                  in x, y and z-direction
c
c     'egun_prtcl.data'   : print of the coordinates of the particle icont (see RDBEAM) along the DC EGUN
c
c     'champ_sol.data'    : plot  field in the solenoid
c   ****************************************************************************************
c  SV 24/08/2015
       open(13,file='emlg.data',status='unknown')
cold       open(14,file='emtr.data',status='unknown')
comment       open(15,file='chemlg.data',status='unknown')
comment       open(17,file='rms.size',status='unknown')
comment       open(18,file='chemtr.data',status='unknown')
comment       open(19,file='nxyz.data',status='unknown')
comment       open(21,file='ntxyz.data',status='unknown')
comment       open(49,file='egun_prtcl.data',status='unknown')
c   ****************************************************
comment       write(14,7787)
comment       write(18,7787)
comment7787   format(4x,'N',6x,'length',5x,'emitx',5x,'emity')
comment       write(13,7788)
comment       write(15,7788)
comment7788   format(4x,'N',6x,'length',5x,'emitz',6x,'energy',1x,'part.good')
c********************************************************
c  SV 24/10/2015
cc       write(13,99)
cc99     format(3x,'ncav',2x,'ncell',2x,'pos I(cm)',4x,'pos S(cm)',6x,
cc     *'L cell(cm)',5x,'W(Mev)',7x,'dw(Mev)',6x,'TOF(dg)',7x,'Ph RF(dg)'
cc     *,6x,'pos M(cm)',5x,'WM(Mev)',7x,'dwM(Mev)',6x,
cc     *'TOFM(dg)',5x,'PhM RF(dg)')
c**************************************************************************************
c    initialize constants
C    in corresponds to dynac input file
       in=7
C    itout corresponds to TRACE3D input file
       itout=48
c
       PI=4.*ATAN(1.)
c initialize the number of charge states to 1
       ncstat=1
c initialize the # of zones (see ZONES card) to zero
       izrot=.false.
       ofeldf=' '
       ofelds=' '
       icour=0
c initialize trace3d related stuff
c frequency indicator
       fid=1.
       T3D=.false.
       xiset=.false.
       kt3h=1
       trace3h(kt3h)=' $DATA'
       kt3t=0
c DONE:
c      DATA KLE/'GEBEAM','INPUT',        ,'ETAC','DRIFT',
c     1*        'QUADRUPO',         ,         ,'SOLENO',        ,
c     2*        'BMAGNET',          ,'CAVSC'  ,        ,        ,
c     3*        'BUNCHER',          ,'NEWF'   ,        ,        ,
c TO LOOK AT:
c      DATA KLE/                ,'RDBEAM',
c     1*                 ,'SEXTUPO','QUADSXT',        ,'SOQUAD',
c     2*        '       ',  'CAVMC' ,         ,'FIELD' , 'HARM' ,
c     3*                 ,'RFQCL',    ,'NREF','SCDYNAC',
c     4*        'SCDYNEL','SCPOS','TILT','TILZ','CHANGREF',
c     5*        'TOF','REJECT','ZROT','ALINER',
c     7*        'RWFIELD','RANDALI','TWQA','EMIPRT',
c     8*        'MMODE','RFQPTQ','STRIPPER','STEER',
c     9*        'PROFGR','SECORD','RASYN','FDRIFT','FSOLE',
c    10*        'EGUN','COMPRES','REFCOG','FPART','QUAELEC',
c    11*        'QUAFK','CAVNUM','EDFLEC','EMITL','RFKICK'
C
c IGNORE: 'ZONES','STOP','T3D','DCBEAM','FIRORD','COMMENT'
C         'EMIT','EMITGR','WRBEAM','ENVEL','REJECT','CHASE'
C         'ACCEPT',
       irfqp=.false.
       ttvols=0.
       iprf=1
       ncell=0
       nzone=0
       eflvl=0.
       rflvl=0.
       bflvl=0.
       nbemit=0
       iell=0
       imcs=0
       IFW=1
       WDISP=1000.
       WPHAS=400.*pi
       WX=100.
       WY=100.
       RLIM=140.
       NRTRE=0
       do i=1,6
         centre(i)=0.
       enddo
       do i=1,maxcell1
         davprt(i)=""
       enddo
       DO i=1,iptsz
         ichas(i)=1
         ichxyz(i)=1
       ENDDO
       SHIFT=.FALSE.
       ICHAES=.FALSE.
       ISEOR=.FALSE.
       IRAYSH=.FALSE.
       ICONT=iptsz+5
       IPRIN=1
       NRRES=0
       NDTL=0
       NCAVNM=0
       NCAVMC=0
       NRBUNC=0
       DTIPH=0.
       icont=150000
       intgct=0
       itye=1000
       ierpf=0
       ialin=.false.
       itwist=.false.
       ichsp=.FALSE.
       iesp=.FALSE.
       itvol=.FALSE.
       imamin=.false.
       chasit=.false.
       ifield=.false.
       iemgrw=.false.
       idav=0
       davtot=0.
       iapel=1
       iaqu=1
       xpsc=.5
       FRACTX=1.
       FRACTY=1.
       FRACTL=1.
       ISCSP=0
       IFILE=11
       IMAX=0
c default for beams is to be bunched (not continuous)
       ifcont=.FALSE.
c   *******************************************************
       WRITE(6,3101)
       WRITE(16,3101)
       WRITE(12,3101)
3101   FORMAT('****** DYNAC V6.0R15 (Beta test), 13-Dec-2015 *******')
       if (gfortran) then
         write(16,*) 'Input file: ',infiln
         write(12,*) 'Input file: ',infiln
       endif
       text=' '
       call mytime(iitime)
       if (mg) then
c using MINGW style gfortran format: 03/30/10 20:51:06 (10 is 2010)
         text(1:4)='    '
         if(iitime(1:2).eq.'01')text(5:7)='Jan'
         if(iitime(1:2).eq.'02')text(5:7)='Feb'
         if(iitime(1:2).eq.'03')text(5:7)='Mar'
         if(iitime(1:2).eq.'04')text(5:7)='Apr'
         if(iitime(1:2).eq.'05')text(5:7)='May'
         if(iitime(1:2).eq.'06')text(5:7)='Jun'
         if(iitime(1:2).eq.'07')text(5:7)='Jul'
         if(iitime(1:2).eq.'08')text(5:7)='Aug'
         if(iitime(1:2).eq.'09')text(5:7)='Sep'
         if(iitime(1:2).eq.'10')text(5:7)='Oct'
         if(iitime(1:2).eq.'11')text(5:7)='Nov'
         if(iitime(1:2).eq.'12')text(5:7)='Dec'
         text(8:8)=' '
         text(9:10)=iitime(4:5)
         text(11:13)=' 20'
         text(14:15)=iitime(7:8)
         text(16:19)=' at '
         text(20:27)=iitime(10:17)
       else
c standard gfortran format: Tue Mar 30 20:51:06 2010
         text(1:11)=iitime(1:11)
         text(12:15)=iitime(21:24)
         text(16:19)=' at '
         text(20:27)=iitime(12:19)
       endif
       write(6,789) text(1:27)
789    format('Started on ',A27)
       write(16,*) 'Started on ',text(1:27)
       write(12,*) 'Started on ',text(1:27)
       call cpu_time(exstrt)
C  Read title
       READ (in,3333) titre(1:80)
       WRITE(16,3334) titre(1:80)
       WRITE(12,3334) titre(1:80)
3333   FORMAT(A80)
3334   FORMAT(1x,A80)
200    CONTINUE
       READ (in,3333) cmnt(1:80)
       if (cmnt(1:1).eq.';') then
         write(16,3334) cmnt(1:80)
         goto 200
       else
         kley=cmnt(1:8)
       endif
       DO I=1,ncards
         IF(KLEY.EQ.KLE(I))  GO TO 220
       ENDDO
       indic=999
       WRITE(66,*) indic
       WRITE(16,111)KLEY
       WRITE(6,111) KLEY
c key not found, go to STOP
       go to 1010
220    CONTINUE
       GO TO (1,2,3,4,5,6,7,8,9,10,
     *        11,12,13,14,15,16,17,18,19,20,
     *        21,22,23,24,25,26,27,28,29,30,
     *        31,32,33,34,35,36,37,38,39,40,
     *        41,42,43,44,45,46,47,48,49,50,
     *        51,52,54,55,56,57,58,59,60,61,
     *        62,63,64,53),I
c*********************************************************************
c  GEBEAM     1  INPUT     2  RDBEAM    3   ETAC    4   DRIFT     5
c  QUADRUPO   6  SEXTUPO   7  QUADSXT   8   SOLENO  9   SOQUAD    10
c  BMAGNET    11 CAVMC     12 CAVSC     13  FIELD   14  HARM      15
c  BUNCHER    16 RFQCL     17 NEWF      18  NREF    19  SCDYNAC   20
c  SCDYNEL    21 SCPOS     22 TILT      23  TILZ    24  CHANGREF  25
c  TOF        26 REJECT    27 ZROT      28  ALINER  29  ACCEPT    30
c  EMIT       31 EMITGR    32 COMMENT   33  WRBEAM  34  ENVEL     35
c  CHASE      36 RWFIELD   37 RANDALI   38  TWQA    39  EMIPRT    40
c  MMODE      41 RFQPTQ    42 STRIPPER  43  STEER   44  ZONES     45
c  PROFGR     46 SECORD    47 RASYN     48  FDRIFT  49  FSOLE     50
c  EGUN       51 COMPRES   52 REFCOG    54  FPART   55  QUAELEC   56
c  QUAFK      57 CAVNUM    58 EDFLEC    59  EMITL   60  RFKICK    61
c  FIRORD     62 DCBEAM    63 T3D       64  STOP    53
c*********************************************************************
1      CONTINUE
C AFTER GEBEAM:Generates randomly the input beam
       write(16,*) ' TYPE CODE:GEBEAM********* '
       CALL MONTE
       write(16,*) '********************************'
       GO TO 200
2      CONTINUE
C AFTER INPUT: define the dynamics at input
c --- INPUT must be preceded by GEBEAM
c --- the reference (synchronous) particle is the c.o.g of the bunch
c --- ENTRY:
c ---- 1) uem, atm, qst
c ---    uem : unit of Rest mass in MeV
c  Examples:
c          proton:938.27231  MeV
c          H-    :939.3145   MeV
c          mesons:33.9093    MeV
c          pions :139.5685   MeV
c          kaons :493.667    MeV
c          electrons : 0.511 MeV
c --- atm : Atomic number
c --- qst : charge (unit of charge) of the reference
c
c ----  2) enedep,tof
c ---      enedep: Kinetic energy of the reference
c ---      tofini: Time of flight of the reference at input (sec)
       write(16,*) ' TYPE CODE:INPUT ********'
       CALL ENTRE
       write(16,*) '********************************'
       GO TO 200
3      CONTINUE
C   AFTER RDBEAM: Read the input beam in the disk
       write(16,*) 'TYPE CODE:RDBEAM **********'
c --- ENTRY
c ---- 1) filen : filename of file containing the particle distribution
c ---- 2) iflag : flag type of distribution file
c ---     iflag = 0 particle coordinates x,xp,y,yp,z,zp with z in rad
c ---     iflag = 1 particle coordinates x,xp,y,yp,z,zp,q,m0 with z in rad
c ---     iflag = 2 particle coordinates x,xp,y,yp,z,zp,q with z in rad
c ---     iflag = 10 particle coordinates x,xp,y,yp,z,zp with z in ns
c ---     iflag = 11 particle coordinates x,xp,y,yp,z,zp,q,m0 with z in ns
c ---     iflag = 12 particle coordinates x,xp,y,yp,z,zp,q with z in ns
c ---     add 100 to the above iflg values to read a file that has Wref as
c ---     4th parameter on line 1
c ---- 3) freq, tof
c ---     freq: RF frquency (MHz)
c ---     tof: phase offset to be applied  to the beam (deg)
c
c ---- 4) uem, atm, qst
c         like ENTRY in type code INPUT
c
       READ(in,3333) myfile(1:80)
       write(16,*) 'Distribution file: ',myfile(1:80)
       open(55,file=myfile,status='unknown')
       CALL ADJRFQ
       write(16,*) '********************************'
       GO TO 200
C   AFTER ETAC :possibility of a multiple charge state beam
4      CONTINUE
       write(16,*) 'TYPE CODE:ETAC **********'
       kt3h=kt3h+1
       trace3h(kt3h)='ERROR: ETAC multi-charge states not supported'
       CALL ETAC
       write(16,*) '********************************'
       GO TO 200
C  AFTER DRIFT : dl:drift length (cm),if dl is negatif no space charge effect
5      CONTINUE
       write(16,*) 'TYPE CODE:DRIFT **********'
       READ(IN,*)DL
       kt3t=kt3t+1
       write(tif,6001) kt3t,kt3t,10.*DL
6001   format('  nt(',i4,')= 1, a(1,',i4,')=',f12.6)
       trace3t(kt3t)=tif
       CALL DRIFT(DL)
       write(16,*) '********************************'
       GO TO 200
6      CONTINUE
C    AFTER QUADRUPO: quadrupole
c   BQUAD: field at pole tip (kG)
C   XLQUA: EFFECTIVE LENGHT (cm )
C   RG:    APERTURE RADIUS (cm)
       write(16,*) 'TYPE CODE:QUADRUPO**********'
       READ(IN,*)XLQUA,BQUAD,RG
       kt3t=kt3t+1
       write(tif,6002) kt3t,kt3t,10.*BQUAD/RG,XLQUA*10.
6002   format('  nt(',i4,')= 3, a(1,',i4,')= ',f9.5,' , ',f9.5)
       trace3t(kt3t)=tif
       call qalva(bquad,xlqua,rg)
       write(16,*) '********************************'
       GO TO 200
7      CONTINUE
C    AFTER SEXTUPO: sextupole
c     IMK2: IFLAG (see arg)
c     ARG = KS2 (cm-3) if IMK2 = 0, otherwise ARG = BSEX (kG)
C     XLSEX : EFFECTIVE LENGHT (CM )
C     RG :APERTURE RADIUS (CM)
       write(16,*) 'TYPE CODE:SEXTUPO**********'
       READ(IN,*)IMK2,ARG,XLSEX,RG
       call sextu(imk2,arg,xlsex,rg)
       write(16,*) '********************************'
       GO TO 200
8      CONTINUE
C    AFTER QUADSXT:
C       quadrupole associated with sextupole
c    IKSQ: IFLAG (see ARGS and ARGQ)
c    ARGS: strength of SEXTUPOLE
c      IKSQ = 0, then ARGS = KS2 (cm-3), otherwise ARGS = FIELD BS(kG)
c    ARGQ: strength of QUADRUPOLE
c     If IKSQ = 0, then ARGQ = K2 (cm-2), otherwise ARGQ = FIELD BQ(kG)
C    XLQUA : EFFECTIVE LENGHT OF THE LENS(cm)
C    RG : APERTURE RADIUS OF THE LENS (cm)
       write(16,*) 'TYPE CODE:QUADSXT**********'
       READ(IN,*)iksq,args,argq,xlqua,rg
       call qasex(iksq,args,argq,xlqua,rg)
       write(16,*) '********************************'
       GO TO 200
C   AFTER SOLENO: solenoid
9      CONTINUE
C     IMKS: IFLAG (see ARG)
C     ARG: IMKS = 0 then ARG is the strength K (cm-1), otherwise ARG is the field B (kG)
C     XLSOL : EFFECTIVE LENGHT (CM )
       write(16,*) 'TYPE CODE SOLENO**********'
       READ(IN,*)IMKS,XLSOL,ARG
       call solnoid(imks,arg,xlsol)
       write(16,*) '********************************'
       GO TO 200
10     CONTINUE
C   AFTER SOQUAD: solenoid associated with quadrupole
c --IKSQ: IFLAG
c --ARGS: STRENGTH or FIELD OF SOLENOID
c     If IKSQ = 0 then ARGS = K (cm-1), otherwise ARGS = B(kG)
c --ARGQ: STRENGTH or FIELD of QUADRUPOLE
c     If IKSQ = 0 then ARGQ = K2 (cm-2), otherwise ARGQ = B(kG)
c    SIGN CONVENTIONS:
c     SOLENOID: K positive => rotate the transverse coordinates about the
c                             z-axis in the clockwise direction.
c     QUADRUPOLE: K2 positive => focusing in the H plane (x,z)
C --XLSOL : EFFECTIVE LENGHT OF THE LENS(cm)
C --RG : APERTURE RADIUS OF THE LENS (cm)
        write(16,*) 'TYPE CODE:SOQUAD**********'
        READ(IN,*) IKSQ,ARGS,ARGQ,XLSOL,RG
        call solquad(iksq,args,argq,xlsol,rg)
        write(16,*) '********************************'
        GO TO 200
11      CONTINUE
c AFTER BMAGNET (sign conventions as in the code TRANSPORT)
c NSECTOR: number of sectors dividing the bending magnet
C WEDGE BENDING MAGNET
C   ANGL : DEG  bend angle of the central trajectory
C   RMO  : CM   radius of curvature of the central trajectory
C   BAIM : KG   field of the bending magnet
c     BAIM = 0  the field is computed from the momentum of the reference
c               otherwise the momentum is computed from the field
C   XN   :      FIELD GRADIENT (dimensionless,TRANSPORT: n)
C   XB   :      NORMALIZED SECOND DERIVATIVE OF B (TRANSPORT : beta)
C   AP(1) = AP(2) CM vertical half aperture (only if IPOLE = 0)
C ENTRANCE FACE
C  PENT1 EK1 EK2 RAB1
c   PENT1: DEG   angle of pole face rotation  (deg)
c   RAB1 : CM    radius of curvature
c   EK1  :       integral related to the extent of the fringing field (TRANSPORT K1)
c   EK2  :       integral related to the extent of the fringing field (TRANSPORT K2)
C   AP(1) : CM   vertical half aperture
C EXIT FACE
c  PENT2 SK1 SK2  RAB2
c   PENT2: DEG   angle of pole face rotation
c   RAB2 : CM  radius of curvature
c   SK1  :     integral related to the extent of the fringing field
c   SK2  :     integral related to the extent of the fringing field
C   AP(2) : CM   vertical half aperture
c
c   SPACE CHARGE COMPUTATION is automatically provided in the routine
c   CAUTION: with several states charges in the beam, ONLY THE SCHEFF METHOD MUST BE USED in the bending magnet
c --------------------------------------------------------------------------------------------------------------
c   nsprint: allow plotting the beam after the sector number nsprint in file 13 ('emlg.data')
ccc        read(in,*) nsector,nsprint
        read(in,*) nsector
        READ(IN,*)ANGL,RMO,BAIM,XN,XB
        READ(IN,*)PENT1,RAB1,EK1,EK2,APB(1)
c if apb(1) ne 0 and if ek1 = 0 the program inserts a default value of ek1 = 0.5
        if(apb(1).ne.0..and.ek1.lt.0.) ek1=0.5
        READ(IN,*)PENT2,RAB2,SK1,SK2,APB(2)
c if apb(2) ne 0 and  sk1 = 0 the program inserts a default value of sk1 = 0.5
        if(apb(2).ne.0..and.sk1.lt.0.) sk1=0.5
c trace3d
c edge
       kt3t=kt3t+1
       write(tif,6005) kt3t,kt3t,PENT1,10.*abs(RMO),
     *  20.*APB(1),EK1,EK2
6005    format('  nt(',i4,')= 9, a(1,',i4,')=',f9.5,' , ',f9.2,
     *  ' , ',f9.2,' , ',f9.2,' , ',f9.2)
       trace3t(kt3t)=tif
c bend
       kt3t=kt3t+1
       write(tif,6006) kt3t,kt3t,ANGL,10.*abs(RMO)
6006    format('  nt(',i4,')= 8, a(1,',i4,')=',f9.5,' , ',f9.2,
     *  ' , 0 , 0 ')
       trace3t(kt3t)=tif
c edge
       kt3t=kt3t+1
       write(tif,6005) kt3t,kt3t,PENT2,10.*abs(RMO),
     *  20.*APB(2),SK1,SK2
       trace3t(kt3t)=tif
c

        CALL aimalv (ANGL,RMO,BAIM,XN,XB,EK1,EK2,PENT1,RAB1,
     *                   SK1,SK2,PENT2,RAB2)
        write(16,*) '********************************'
        GO TO 200
12      CONTINUE
C   AFTER CAVMC: Multicell accelerating element
c           the electromagnetic field can be read:
c           in the file 'field.txt' in the form (z,E(z))
c           or from the type code: HARM in the form of a Fourier series expansion
       write(16,*) 'TYPE CODE:CAVMC  **********'
       CALL RESTAY
       write(16,*) '********************************'
       GO TO 200
C   AFTER CAVSC: single cell accelerating element
c                The transit time factors describe the eletromagnetic field
13    CONTINUE
      write(16,*) 'TYPE CODE:CAVSC  **********'
      CALL ETGAP
      write(16,*) '********************************'
      GO TO 200
14    CONTINUE
c  AFTER FIELD
c --- the axial field of the cavity is read from file in the form (z,E(z)) with a step size h in z
c     in the SUPERFISH format: z(m)  E(z) (V/m)
c  PART: the step size h may be divided in 'part' elements (with PART (= >) 1)
c  ATT:  the field E(z) is multiplied by ATT and converted in MV/cm
c
       write(16,*) 'TYPE CODE: FIELD  **********'
       write(16,*) 'ELECTRIC FIELD (z, E(z) ) '
c      get filename of input file of the electromagnetic field in the form (z,E(z))
       READ(in,3333) myfile(1:80)
c       write(16,*) 'Electromagnetic field file: ',myfile(1:80)
c      check if file is already open
       iostats = int(FTELL(20))
       if(iostats.eq.-1) then
c file not yet open
         open(20,file=myfile,status='unknown')
         write(16,*) 'Opening field file: ',myfile(1:80)
         ofeldf=myfile
       else
         if(ofeldf.ne.myfile) then
           write(16,*) 'Closing field file: ',ofeldf
           write(16,*) 'Opening field file: ',myfile(1:80)
           close(20)
           open(20,file=myfile,status='unknown')
           ofeldf=myfile
         endif
       endif
       read(in,*) att
comment       write(16,*) ' partition of a step h: ',part,' field factor: ', att
c     conversion V/m (SUPERFISH) --> MV/cm
       att=1.e-08*att
       write(16,*) ' ** Read the cavity field from ',myfile(1:80)
       call fieldcav(att)
       write(16,*) '********************************'
       go to 200
15    CONTINUE
c AFTER HARM: the field is read on the form of a Fourier series expansion
       write(16,*) ' TYPE CODE: HARM  **********'
       write(16,*) 'ELECTRIC FIELD (Fourier series expansion) ****'
       call rharm
       write(16,*) '********************************'
       go to 200
16     CONTINUE
C  After BUNCHER: buncher as a single element
c      PV: Voltage
c      PDP: PHASE OF RF (deg)
c      PRLIM: aperture radius (cm)
c      PHARM: harmonic factor (bucher fq.)/(DTL freq.)
       write(16,*) 'TYPE CODE:BUNCHER ************************'
       READ(IN,*) pv,pdp,pharm,prlim
       write(16,7777) pv,pdp,prlim
7777   format(' BUNCHER CAVITY ',/,' Voltage ',E12.5,' MV',/,
     X ' RF Phase ',e12.5,' deg',' Aperture Radius ',e12.5,' cm')
       kt3t=kt3t+1
       write(tif,6007) kt3t,kt3t,pv,pdp,pharm
6007   format('  nt(',i4,')=10, a(1,',i4,')=',f9.5,' , ',f9.2,
     &        ' , 1 , 1 , ',f5.1)
       trace3t(kt3t)=tif
       pdp=pdp*pi/180.
       call bunparm(pv,pdp,pharm,prlim)
       write(16,*) '********************************'
       GO TO 200
17     CONTINUE
C   AFTER RFQCL (single cell of a RFQ)
      write(16,*) 'TYPE CODE:RFQCL ************************'
c --- The parameter definitions of the cell are identical to the ones of TRACE3-D
c       VR02: Maximum intervane potentiel difference divided by square of average vanne displacement kV/(mm2)
c       AV  : product of accelerating efficience and maximum intervane voltage (kV)
c       XLRFQ: cell length (mm)
c       XPHRFQ: phase of RF (deg)
c       TYPE: TYPE is one of the following
c         0 = standard cell, no acceleration
c         1 = standard cell, acceleration
c         2 = fringing field, no acceleration
c         3 = fringing field, acceleration
c
       READ(IN,*) VR02,AV,XLRFQ,XPHRFQ,TYPE
c -----     convert in MV , m  in the array prfq(9)
c ---- prfq(1) : VR02 ( MV/(m*m) )
c ---- prfq(2) : AV (MV)
c ---- prfq(3) : cell length (m)
c ---- prfq(4) : phase of RF (deg)
c ---- prfq(5) : TYPE
      prfq(1)=vr02*1.e03
      prfq(2)=av*1.e-03
      prfq(3)=xlrfq*1.e-03
      prfq(4)=xphrfq
      prfq(5)=type
c --- The space charge routine SCHEFF only is available
       if(ichaes) then
        if(iscsp.lt.3) then
        write(6,*)'*** HERSC and SCHERM cannot be used in the RFQ'
        write(16,*)'*** HERSC and SCHERM cannot be used in the RFQ'
         go to 53
        endif
        write(16,*)'***** beam current: ',beamc,' mA'
       endif
      call rfq_o3
      write(16,*) '********************************'
      GO TO 200
C  AFTER NEWF:define a new frequency (Hertz)
18    CONTINUE
      IF (IMAX.EQ.0) THEN
        WRITE(16,*)' INIT TOF HAS TO BE PRECEEDED BY GBEAM OR RDBEAM'
      GO TO 200
      ENDIF
      write(16,*) 'TYPE CODE:NEWF ************************'
      READ (IN,*) FH1
c record frequency change in view of trace3d file
      FH1=2.*PI*FH1
      fid=fh1/fh
c adjust beam current
      beamc=beamc*fh1/fh
      rflvl=rflvl*fh1/fh
      fh=fh1
comment       do i=1,ngood
comment        f(6,i)=f(6,i)-tref
comment       enddo
comment       TREF=0.
comment       WRITE(16,*)'******** NEWF: TOF HAS BEEN RESET TO ZERO ********'
      WRITE(16,*)' NEW FREQUENCY : ',FH/(2.*PI),' Hertz'
      GO TO 200
C   AFTER NREF: define a new syncronous particle
19     CONTINUE
C ---- DEPHAS: the change of phase (DEG)
C ---- DEW   : the change of kinetic energy (see IREFW)
C ---- follow two falgs, IREF and IREFW
C ----  IF IREF=0: RELATIVE TO the previous synchronous particle
C ----  IF IREF=1: RELATIVE TO the previous COG
C ----  IF IREFW=0: DEWREF is in dW/W
C ----  IF IREFW=1: DEWREF is in dW (MeV)
       write(16,*) ' TYPE CODE:NREF  ***********'
       read (in,*) dephas,dewref,iref,irefw
       CALL REFER
       write(16,*) '********************************'
       GO TO 200
20     CONTINUE
c AFTER SCDYNAC
C   ISCSP: METHOD FOR SPACE CHARGE COMPUTATIONS
c     ISCSP=1  HERSC METHOD
C     ISCSP=2  SCHERM METHOD
C     ISCSP=3  SCHEFF METHOD
       write(16,*) ' TYPE CODE:SCDYNAC  ***********'
       read(in,*) iscsp
c      sce10 =1 : call in quads,solenoids,accelarating elements
c      sce10 =2 : call in drifts,accelarating elements
c      sce10 =3 : call in quads,solenoids,drifts,accelarating elements
C      BEAM CURRENT IN ma
       READ(IN,*) BEAMC,sce10
       if(.not.xiset) then
         kt3h=kt3h+1
         write(tif,7001) BEAMC
7001     format(' XI= ',f9.4)
         trace3h(kt3h)=tif
         xiset=.true.
       endif
       if(iscsp.le.1) write(16,*) 'HERSC method  '
       if(iscsp.eq.2) write(16,*) 'SCHERM method  '
       if(iscsp.eq.3) write(16,*) 'SCHEFF method  '
       if(iscsp.gt.3) then
        write(16,*) 'Error in SCDYNAC iscsp: ',iscsp
        stop
       endif
       WRITE(16,*)' Beam current : ',BEAMC,' mA'
       ECT=4.
       ICHAES=.TRUE.
       if(iscsp.le.1) then
c        initialise the routine HERSC
         if(iscsp.eq.1)ini=0
         if(iscsp.lt.1)ini=-1
         call hersc(ini)
         iscsp=1
       endif
c---  SCHERM, read third line as dummy
       if(iscsp.eq.2) read(in,*) idum
c --- SCHEFF
c      initialise the mesh of routine SCHEFF
       if(iscsp.eq.3) then
        call schfdyn
       endif
       if (beamc.eq.0.) ichaes=.FALSE.
       write(16,*) '****************************'
       GO TO 200
21     CONTINUE
C   AFTER SCDYNEL: space charge computations in the current position (i.e.bending magnet)
C     XTRANS is the acting lenght of the beam self-fields (cm)
       write(16,*) ' TYPE CODE:SCDYNEL *******'
       READ(IN,*) XTRANS
       CALL CESP(XTRANS)
       write(16,*) '****************************'
       GO TO 200
22    CONTINUE
c AFTER SCPOS : change the position of space charge computation in gaps or cavities
       write(16,*) ' TYPE CODE:SCPOS(space charge position)**********'
       read(in,*) xpsc
       if(xpsc.ge.1.)xpsc=.5
       write(16,*) '****************************'
       go to 200
23     CONTINUE
C   AFTER TILT:rotation and shift of beam ellipsoid
c    ICG    : = 1 => REFERENCE PARTICLE IS THE C.O.G. OF THE BEAM
c    ICG    : = 0 => IT IS DISTINCT FROM THE C.O.G.
       write(16,*) ' TYPE CODE:TILT *********************'
       READ(IN,*)ICG
       READ(IN,*)TIPHA,TIX,TIY,SHIFW,SHIFP
       CALL TILTBM(ICG)
       write(16,*) '****************************'
       GO TO 200
24     CONTINUE
C   AFTER TILZ: tilt in the plane (x,z) around the c.o.g. of the upright ellipse
       write(16,*) ' TYPE CODE:TILZ *********************'
       READ(IN,*) TILTA
       call tiltz(TILTA)
       write(16,*) '****************************'
       GO TO 200
25     CONTINUE
C  AFTER CHANGREF: change of reference frame
       write(16,*) ' TYPE CODE:CHANGREF *********************'
       CALL CHREFE
       GO TO 200
26     CONTINUE
C After TOF:
c ---  the T.O.F may be activated in the dynamics of bunchers, cavities and acc. gaps
c ---  Entry: indic and icor
c ---  indic = 0 : the T.O.F is activated, otherwise it is passive
c ---  icor = 0  : no adjustement on the phase offset, otherwise adjustments are automatically made on the phase offset
       write(16,*) ' TYPE CODE: TOF ************'
       CALL rmami
       write(16,*) '****************************'
       GOTO 200
C   AFTER REJECT: defining limits in X,X',Y,Y',Z,Z'
c ---- aperture of the Beam
c ----  WDISP: in half dispersion
C ------  IF IFW=0  ==> WDISP in (+-) dW/W
C ------  IF IFW=1  ==> WDISP in (+-) dW (MeV)
c ------  IF IFW=10 ==> wdisp = dW/W relative to REF
c ------  IF IFW=11 ==> wdisp = dW relative to REF
c ----  WPHAS: in half phase (+-) deg
c ----  WX   : in x-direction (+-) cm
c ----  WY   : in y-direction (+-) cm
c ----  RLIM : in radius (cm)
27     CONTINUE
       write(16,*) ' TYPE CODE:REJECT*********************'
       READ(IN,*) IFW,WDISP,WPHAS,WX,WY,RLIM
       if(ifw.eq.0)WRITE(16,1050)WDISP,WPHAS,WX,WY,RLIM
       if(ifw.ne.0)WRITE(16,1051)WDISP,WPHAS,WX,WY,RLIM
1050   FORMAT(5X,' *** BEAM SIZE LIMITS ',/,
     X 4X,' 1/2 dW/W :',E12.5,'  1/2 PHASE(DEG) :'
     X ,E12.5,/,4X,' 1/2 x (cm) :',E12.5,
     X ' 1/2 y(cm) :',E12.5,
     X ' RADIUS (cm) :',E12.5)
1051   FORMAT(5X,' *** BEAM SIZE LIMITS ',/,
     X 4X,'1/2 dW (MeV) :',E12.5,'  1/2 PHASE(DEG) :'
     X ,E12.5,/,4X,' 1/2 x (cm) :',E12.5,
     X ' 1/2 y(cm) :',E12.5,
     X ' RADIUS (cm) :',E12.5)
c ---- convert WPHAS in rad
       WPHAS=WPHAS*PI/180.
       write(16,*) '****************************'
       GO TO 200
28     CONTINUE
c AFTER ZROT : beam rotation
c  The transverse coordinates x and y may be rotated through an
c  angle ZROTA(deg).The positive sense of ratation is clockwise
c  about the positive z axis
      write(16,*) ' TYPE CODE:ZROT*********************'
      READ(IN,*)ZROTA
      CALL ZROTAT(ZROTA)
      write(16,*) '****************************'
      GO TO 200
29    CONTINUE
C   after ALINER: ALIGNMENT errors IN X,X',Y,Y'
C     XL,YL (cm)    XPL,YPL (mrad)
      write(16,*) ' TYPE CODE:ALINER*********************'
      read(in,*) XL,YL,XPL,YPL
      CALL ALINER
      write(16,*) '****************************'
      GO TO 200
30    CONTINUE
C    AFTER ACCEPT: Determination of the input acceptance for the structure
       write(16,*) ' TYPE CODE:ACCEPT*********************'
       CALL ACCEPT
       write(16,*) '****************************'
       GO TO 200
31     CONTINUE
C    AFTER EMIT  Print emittance data in the file 'dynac.short'
       CALL emiprt(0)
       GO TO 200
32     CONTINUE
C    AFTER EMITGR: emittance plots
C    PLOTS IN XX', YY', XY AND ZZ'
      write(16,*) ' TYPE CODE:EMITGR*********************'
      igrprm=0
      call ytzp
      write(16,*) '****************************'
      GO TO 200
33    CONTINUE
C    AFTER COMMENT
C allows for comments in the input data file
       write(16,'(a8)') kley
       READ (in,'(A)') cmnt(1:80)
       WRITE(16,'(A)')  cmnt(1:80)
       GO TO 200
34     CONTINUE
C   AFTER WRBEAM: prints coordinates of particles on output files
       write(16,*)'WRBEAM output coordinates of particles'
c      irec:flag   irec=0 the phase is recentered with regard to the c.o.g
c                  irec<>0 the phase is not recentered around the c.o.g
       read(in,'(A)') wfile
       write(16,'(A)') 'Distribution will be written to ',wfile
       read(in,*) irec,iflg
c       iflag=0
       call prbeam(iflg,wfile)
       write(16,*) '****************************'
       GO TO 200
35     CONTINUE
C   AFTER ENVEL: plot the longitudinal and the tranverse envelope of the beam
       write(16,*) ' ENVEL *********************'
       CALL PROFIL
       write(16,*) '****************************'
       GO TO 200
36     CONTINUE
C   AFTER CHASE:Temporary elimination of most distant particles for statistical purposes
       write(16,*) ' TYPE CODE:CHASE ********************'
       CALL CHASE
       write(16,*) '****************************'
       GO TO 200
37     CONTINUE
c  AFTER RWFIELD : rewinds the file 'field.txt'
       rewind(20)
       go to 200
38     CONTINUE
c  AFTER RANDALI : generates random errors in alignments
C     XL,YL (cm)    XPL,YPL (mrad)
c     ilier = 0 stop the effects of the random misalignment
      write(16,*) ' TYPE CODE:RANDALI*********************'
      read(in,*) ilier
      if(ilier.eq.0) then
       ialin=.false.
       write(16,*) '****************************'
       go to 200
      else
       ialin=.true.
       read(in,*) XL,YL,XPL,YPL
      endif
      write(16,*) '****************************'
      go to 200
39    CONTINUE
c  AFTER TWQA generates systematic or random twist of quadrupoles
c  QTWIST: rotation about Y axis (deg)
c  IQRAND: =0 systematic twist, otherwise random twist
       write(16,*) ' TYPE CODE TWQA*********************'
       read(in,*) iqrand,qtwist
       itwist=.true.
       if(abs(qtwist).le.1.e-20) itwist=.false.
       write(16,*) '****************************'
       go to 200
40    CONTINUE
c  AFTER EMIPRT : print the beam characteristics in the disk(tape12, file='dynac.short')
c --- the beam characteristics are systematically print after:
c       cavities, accelerating gaps, bunchers, electrons gun, rfq
c
c  IEMQESG :  =0  stop the prints for all optical lenses
c             =1  after all optical lenses apart from positive drifts
c             =2  after all optical lenses and positive drifts
c             =3  after quads, solen.,positive drifts and accel. elements
       iemgrw=.true.
       read(in,*) iemqesg
       if(iemqesg.eq.0) iemgrw=.false.
       go to 200
41    CONTINUE
c AFTER MMODE: systematic or random error on the phase offset and on the level of the field
c ---          MMODE is only acting on  particles in the bunch  (no change of the reference)
c --- ENTRY: ierpf , vphase , vfield
c --- ierpf: flag
c ---   IF:  ierpf = 0  ===> stop the type code effects
c ---        ierpf = 1  ===> systematic error
c ---        ierpf > 1  ===> random error
c ---  vphase (deg): error added to the nominal phase offset
c ---  vfield (%)  : error added to the level of the electric field
c ---  (new phase offset) = (previous phase offset) + vphase
c ---  (new level of field) = (previous level of field) * (1.+ vfield/100)
       write(16,*) ' TYPE CODE MMODE*********************'
       read(in,*) ierpf,vphase,vfield
       if(ierpf.eq.0) then
        vphase=0.
        vfield=0.
       endif
       if(ierpf.eq.1) write(16,4279) vphase,vfield
4279   format(2x,'systematic error on phase offset: ',e12.2,' deg',/,
     *        2x,'systematic error on level of field: ',e12.5,' %')
       if(ierpf.gt.1) write(16,4290) vphase,vfield
4290   format(2x,'maximun randon error on phase offset: ',
     * e12.2,' deg',/, 2x,'maximum random error on level of field: ',
     * e12.5,' %')
       vfield=vfield/100.
       write(16,*) '****************************'
       go to 200
42    CONTINUE
c --- AFTER RFQPTQ
      write(16,*) 'TYPE CODE:RFQPTQ ************************'
      write(6,*)
c -- ENTRIES:
c ---  ENTRY 1: input file 'myfile' contains the geometry of the RFQ (unit 27)
c ---  ENTRY 2: nceltot
c        nceltot: number of cells (may be less than the total number of cells)
c ---  ENTRY 3: tvolt avolt fph
c ---- tvolt: factor applied to intervane-voltage Vref of the synchronous particle (in %)
c ---- avolt: factor applied to intervane-voltage Vpart for particles (in %)
c        Vref = Vref(1 + tvolt/100)
c        Vpart = Vpart(1 + avolt/100)
c --- fph: factor applied to the phase at entrance of cells in the file myfile  in %)
c        (phase at entrance of cells)= (1 + fph/100) X (phase in myfile)
c  NOTE: fph is only available for cells of type = 0, type = 2 (type = E) or type = 5 (type R)
c --- pib > 0:  shift the particles at the entrance of the RFQ inside (+/-) pi w.r.t.the synchronous particle
c --- pib = 0: no action
c
c
c related files:
c
c     'rfq_list.data'  : list of parameters of the RFQ
c
c     'rfq_listmid.data' : list of phases at input & middle of each RFQ cell
c
c     'rfq_lost.data'  : lists where particles are lost in the RFQ
c
c      'rfq_coeflist.data' list of coefficients
      READ(in,3333) myfile(1:80)
      write(16,*) 'RFQ input data file: ',myfile(1:80)
      open(27,file=myfile,status='unknown')
      open(70,file='rfq_list.data',status='unknown')
      open(75,file='rfq_coef.data',status='unknown')
      open(49,file='rfq_lost.data',status='unknown')
      open(89,file='rfq_listmid.data',status='unknown')
      read(in,*) nceltot
      read(in,*) tvolt,avolt,fph,pib
      write(16,5279)nceltot,tvolt,avolt
5279  format(' RFQ number of cells: ',i5,/,
     *  '  factor on intervane voltage (reference):',e12.5,' %',/,
     *  '  factor on intervane voltage (bunch):',e12.5,' % ')
       tvolt=tvolt/100.
       avolt=avolt/100.
       fph=(1.+fph/100.)
       if(ichaes) then
        write(16,*)'***** beam current: ',beamc,' mA'
        if(iscsp.lt.3) then
        write(6,*)'*** HERSC and SCHERM are not available in the RFQ'
        write(16,*)'*** HERSC and SCHERM are not available in the RFQ'
         go to 53
        endif
       endif
       call cpardyn(pib)
       write(6,*)
       write(16,*) '****************************'
       go to 200
43    CONTINUE
c --- AFTER STRIPPER: stripper foils
c ---  ( based on the works of D.A. Eastham, ref.   )
c --- available for 'slow' hadron particles
c --- STRIPPER FOILS PARAMETERS:
c --- qs : charge (unit of charge)
c --- atms: atomic mass
c --- ths : thickness of the foils (g/cm**2)
c --- PARTICLES
c*2010-12-02 Use Baron formula for charge state distribution in case of carbon foils
c --- anp : atomic number of the projectile
c --- nqst: number of charge states after the stripper
c --- sqst: array holding the nsqt charge states
c --- qop : charge of particles after crossing the stripper foils (unit of charge)
c --- the atomic mass of particles (atm) is the one given in INPUT or in RDBEAM
       write(16,*) 'TYPE CODE:STRIPPER ************************'
       read(in,*)qs,atms,ths,anp
       call stripp
       write(16,*) '****************************'
       go to 200
44     CONTINUE
c   AFTER STEER: thin steering element
C---- PARAMETERS ARE INTEGRATED FIELD fld, nvf
c     fld in units of (Tm) for magnetic steerer
c     fld in units of (kV*m/m) for electrostatic steerer
c     (Voltage * length / plate separation)
c     if nvf=0,  horizontal magnetic steerer
c     if nvf=1,  vertical magnetic steerer
C     if nvf=2,  horzontal electrostatic steerer
C     if nvf=3,  vertical electrostatic steerer
       write(16,*) ' TYPE CODE:STEER*********************'
       read(in,*) fld,nvf
       call steer(fld,nvf)
       write(16,*) '****************************'
c  start of writes in file '.short' for steerer
       idav=idav+1
       iitem(idav)=9
       dav1(idav,1)=fld
       dav1(idav,2)=float(nvf)
       dav1(idav,3)=davtot*10.
c  end daves
       go to 200
45    CONTINUE
C  AFTER ZONES: specify zones of different colours in the bunch
       write(16,*) ' TYPE CODE:ZONES*********************'
       init=1
       call area(init)
       write(16,*) '****************************'
       go to 200
46    CONTINUE
C  AFTER PROFGR: X-Z and Y-Z scatter plots,(X,Y,Z) and (Xp,Yp,Zp) profiles
       write(16,*) ' TYPE CODE:PROFGR*********************'
       igrprm=0
       if (igrprm.eq.0) then
C        READ GRAPH TITLE
         READ(IN,6620) TEXT
6620     FORMAT(A)
c idwdp=0 cog=ref in XZ,YZ plots (for instance for Alvarez structure)
c idwdp=1 cog<>ref in XZ,YZ plots (for instance for IH structrure)
c iskale=0 vertical scale on profile plots is NOT a log scale
c iskale=1 vertical scale on profile plots IS a log scale
         READ(IN,*) idwdp,iskale
C        READ GRAPH LIMITS INTO GLIM(J,K), J=GRAPH NUMBER
C        K=1 HOR. LIMIT , K=2 VERT. LIMIT
         READ(IN,*) GLIM(3,1),GLIM(3,2),GLIM(4,1),GLIM(4,2)
       endif
       call grcomp(text,iskale)
       write(16,*) '****************************'
       go to 200
47    CONTINUE
C  AFTER SECORD : second order matrix for optical lenses
       write(16,*) '*****************************************'
       write(16,*) ' SECOND ORDER IN BEAM TRANSPORT**********'
       ISEOR=.TRUE.
       write(16,*) '*****************************************'
       go to 200
48    CONTINUE
C  AFTER RASYN : Radiation exitation in bending magnets (only for electrons)
       write(16,*) ' SYNCHRTRON RADIATION IN BENDING MAGNET**********'
       IRAYSH=.TRUE.
       write(16,*) '*************************************************'
       go to 200
49    CONTINUE
C  AFTER FDRIFT : Divide a drift length in partial drifts (for space charge computations)
       write(16,*) ' TYPE CODE:FDRIFT*********************'
C  XL    : Total drift length (cm)
C  NPART : number of partial drifts
C  IMIT  : if IMIT not equal zero,print emittance datas in the file 'short.txt'
       read(in,*) xl,npart,imit
       dl=xl/float(npart)
       write(16,*) ' total drift length : ',xl,' cm divided in : '
     *   ,npart,' drifts of : ',dl,' cm'
       call fdrift(xl,npart,imit)
       write(16,*) '*************************************************'
       go to 200
CAFTER FSOLE: solenoid the magnetic field is read from disk in the form (z,B(z))
C              unities: z (m), B(z) kG
C   BCRET : coefficient multiplier, the magnetic is: B(z) * BCRET
C   sign field convention like in the code TRANSPORT
50     CONTINUE
       write(16,*) 'TYPE CODE FSOLE*****************'
       READ(in,3333) myfile(1:80)
c     unit 25 corresponds to an input file giving solenoid magnetic field(s) in the form (z,B(z))
c      check if file is already open
       iostats = int(FTELL(25))
       if(iostats.eq.-1) then
c file not yet open
         open(25,file=myfile,status='unknown')
         write(16,*) 'Opening solenoid field file: ',myfile(1:80)
         ofelds=myfile
       else
         if(ofeldf.ne.myfile) then
           write(16,*) 'Closing solenoid field file: ',ofelds
           write(16,*) 'Opening solenoid field file: ',myfile(1:80)
           close(25)
           open(25,file=myfile,status='unknown')
           ofelds=myfile
         endif
       endif
       read(in,*)bcret,intgr
       write(16,990) intgr
990    format('*** SOLENOID WITH ARBITRARY MAGNETIC FIELD ',/,5x,
     *        'PARTITION IN: ',i4,' ELEMENTARY SOLENOIDS')
       call solfield(bcret,intgr)
       write(16,*) '********************************'
       GO TO 200
51     CONTINUE
C AFTER EGUN:
c     Only the SCHEFF routine is avalaible in the routine EGUN
comment       ichaes=.false.
       write(16,*) 'TYPE CODE EGUN*****************'
c     unit 22 corresponds to an input file giving egun field (z,E(z)) E(z) normalised
       READ(in,3333) myfile(1:80)
       write(16,*) 'Egun field file: ',myfile(1:80)
       open(22,file=myfile,status='unknown')
c  follows the particle ifpt (not active)
c not actived       read(in,*)ifpt
c    fmult: field factor
       read(in,*)fmult,indp
       if(.not.ichaes) indp=1
c  indp: number of space charge computations to be made in the EGUN
c   indp = 1 : 8 space charge computations. The EGUN field is divided in 16 elements
c   indp = 2 : 16 space charge computations.                          in 32 elements
c   indp = 3 : 32 space charge computations.                          in 64 elements
       if(ichaes) then
        write(16,*)'***** beam current: ',beamc,' mA'
        if(iscsp.lt.3) then
        write(6,*)'*** HERSC and SCHERM are not available in the RFQ'
        write(16,*)'*** HERSC and SCHERM are not available in the RFQ'
         go to 53
        endif
       endif
       call egun(fmult,indp)
       write(16,*) '*************************************************'
       go to 200
52     CONTINUE
C AFTER COMPRES
c  Do so by shifting particles belonging to the same bunch outside the (+/-) pib/2 (deg) window
c   w.r.t.the COG to inside the (+/-) pib/2 window w.r.t.the COG
       write(16,*) ' TYPE CODE: COMPRES ********'
c --- pib (deg)
       read(in,*) pib
       pib=pib/2.
       write(16,1890) pib
       call compress(pib)
1890   format('*** shift particles inside +/- ',e12.5, 'deg')
       write(16,*) '*************************************************'
      go to 200
54    continue
C   AFTER REFCOG:
c    ISHIFT  = 0: the synchronous particle and the cog are coinciding (shift = false)
c    ISHIFT  = 1: the synchronous particle and the cog are separated   (shift = true)
c    ISHIFT  = 2: at the start the TOF of synchronous particle  is the TOF of the cog
c                   after synchronous particle and cog are separated (shift = true)
c    (at the begining shift = false)
       write(16,*) ' TYPE CODE:REFCOG *********************'
       read(in,*) ishift
       if(ishift.eq.0) then
        shift=.false.
        write(16,*)' Synchronous particle is the COG of the bunch'
        obref=vref/vl
        otref=tref
        tcog=0.
        wcog=0.
        do ijp=1,ngood
         wcog=wcog+f(7,ijp)
         tcog=tcog+f(6,ijp)
        enddo
        wcog=wcog/float(ngood)
        tcog=tcog/float(ngood)
        gcog=wcog/xmat
        bcog=sqrt(gcog*gcog-1.)/gcog
        vref=bcog*vl
        tref=tcog
        fcpi=fh*180./pi
        write(16,5433)obref,otref*fcpi,vref/vl,tref*fcpi,bcog,tcog*fcpi
5433    format(' old ref. beta: ',e12.5,' TOF: ',e12.5,' deg',/,
     *         ' new ref. beta: ',e12.5,' TOF: ',e12.5,' deg',/,
     *         ' COG      beta: ',e12.5,' TOF: ',e12.5,' deg')
c        write(16,*) 'old TTVOL (deg): ',ottvols*fcpi
c        write(16,*) 'new TTVOL (deg): ',ttvols*fcpi
       endif
       if(ishift.eq.1) then
        shift=.true.
        write(16,*)
     * ' Synchronous particle and COG of the bucnch are independent'
       endif
       if(ishift.gt.1) then
        write(16,*)
     * ' Synchronous particle and COG of the bucnch are independent,',
     * ' but initially TOF and energy of the synchronous particle ',
     * ' are the ones of the bunch'
        shift=.true.
c --- the reference is the cog
        obref=vref/vl
        otref=tref
        wcog=0.
        tcog=0.
        do ijp=1,ngood
         wcog=wcog+f(7,ijp)
         tcog=tcog+f(6,ijp)
        enddo
        wcog=wcog/float(ngood)
        tcog=tcog/float(ngood)
        gcog=wcog/xmat
        bcog=sqrt(gcog*gcog-1.)/gcog
        vref=bcog*vl
        tref=tcog
        ottvols=ttvols
        if(itvol) ttvols=tref
        fcpi=fh*180./pi
c        write(16,5420)obref,otref*fcpi,obref,tref*fcpi
        write(16,5420)obref,otref*fcpi,vref,tref*fcpi
5420    format(' old ref. beta: ',e12.5,' tref: ',e12.5,' deg',/,
     *         ' new ref. beta: ',e12.5,' tref: ',e12.5,' deg')
        write(16,*) 'old TTVOL (deg): ',ottvols*fcpi
        write(16,*) 'new TTVOL (deg): ',ttvols*fcpi
       endif
       write(16,*) '****************************'
       GO TO 200
55     continue
C   AFTER FPART:
c --- a given particle is followed in the accelerating elements
c --- ICONT : the number of the particle followed in accelerating elements
       write(16,*) ' TYPE CODE:FPART *********************'
       read(in,*) icont
       write(16,*) ' the particle:',icont,' is followed'
       write(16,*) '****************************'
       GO TO 200
56     CONTINUE
C    AFTER QUAELEC: electric quadrupole
c       VOLT: voltage at pole tip (kV)
C       XLQUA: effective length  (cm)
C       RS: radial distance of pole tip from axis (cm)
       write(16,*) 'TYPE CODE:QUALEC**********'
       READ(IN,*)XLQUA,VOLT,RS
       call qelec(volt,xlqua,rs)
       write(16,*) '********************************'
       GO TO 200
57     continue
C    AFTER QUAFK:  quadrupole (magnetic or electric)
c       ITYQU: ITYQU = 0 electric quadrupole , otherwise magnetic quadrupole
c       ARG: K2 (cm-2)
C       XLQUA: effective length  (cm)
C       RS: radial distance of pole tip from the axis (cm)
       write(16,*) 'TYPE CODE:QUAFK**********'
       READ(IN,*)ITYQU,ARG,XLQUA,RS
       call qfk (ityqu,arg,xlqua,rs)
       write(16,*) '********************************'
       GO TO 200
58     continue
c    AFTER CAVNUM: numerical computation of multicell cavities
c           the electromagnetic field can be read:on the disk in the form (z,E(z))
c           or from the type code: HARM in the form of a Fourier series expansion
       write(16,*) 'TYPE CODE:CAVNUM  **********'
       call cavnum
       write(16,*) '********************************'
       GO TO 200
59    continue
c     AFTER EDFLEC:  ELECTROSTATIC DEFLECTOR
c --- Input parameters
c   radial radius (cm)
c   bend angle (deg)
c   radii: vertical (radial) radii of curvature (cm)
       call e_deflec
       write(16,*) '********************************'
       GO TO 200
60     CONTINUE
C    AFTER EMITL  Print emittance datas in the file 'dynac.short' but
C                 also read label to written to dynac.short
       read(in,'(A)') shortl
       call emiprt(1)
       GO TO 200
61     continue
C  After RFKICK: Electric RF kicker
c      PV: Voltage Factor (Voltage (kV) * electrode length (m) / gap (m))
c      PDP: PHASE OF RF (deg)
c      PHARM: harmonic factor (bucher fq.)/(DTL freq.)
C      NVF: 0 = horizontal, 1 = vertical
       write(16,*) 'TYPE CODE:RFKICK ************************'
       READ(IN,*) pv,pdp,pharm,nvf
       plane='horizontal'
       if(nvf.eq.1) plane='vertical  '
       write(16,7779) pv,pdp,plane
7779   format(' RF Kicker ',/,' Voltage Factor',E12.5,' kV*m/m',/,
     X ' RF Phase ',e12.5,' deg',' Type: ',a10)
       pdp=pdp*pi/180.
       call rfkick(pv,pdp,pharm,nvf)
       write(16,*) '********************************'
       GO TO 200
62     CONTINUE
C  AFTER FIRORD : first order matrix for optical lenses
       write(16,*) '*****************************************'
       write(16,*) '  FIRST ORDER IN BEAM TRANSPORT**********'
       ISEOR=.FALSE.
       write(16,*) '*****************************************'
       go to 200
63     CONTINUE
C  AFTER DCBEAM : select DC or bunched beam
       READ(IN,*) iscont
       write(16,*) '****************************'
       if(iscont.eq.1) then
         ifcont=.true.
         write(16,*) ' DC BEAM SELECTED    **********'
       elseif(iscont.eq.0) then
         ifcont=.false.
         write(16,*) ' BUNCHED BEAM SELECTED ********'
       endif
       write(16,*) '****************************'
       go to 200
64     CONTINUE
C  AFTER T3D : request TRACE3D input file
       T3D=.true.
       write(16,*) '***********************************'
       write(16,*) ' TRACE3D INPUT FILE WILL BE WRITTEN'
       write(16,*) '***********************************'
       go to 200
53     CONTINUE
c AFTER STOP: mandatory, stop the computations
       write(16,100) kle(ncards)
c write trace file (if requested)
       if(T3D) then
c write trace file (header part)
         open(itout,file='for_trace3d.t3d',status='unknown')
         write(itout,'(A)') trace3h(1)
         write(itout,1990) kt3t,kt3t,kt3t
1990     format(' N1= 1, N2= ',I4,', NEL1= 1, NEL2= ',I4,
     *        ', NP1= 1, NP2= ',I4)
         do i=2,kt3h
           write(itout,'(A)') trace3h(i)
         enddo
         write(itout,1992)
1992     format(' PQEXT= 2.5, ICHROM= 0, IBS= 0, SMAX= 2.0')
         write(itout,1993)
1993     format(' XM=  15.00, XPM=  50.00, YM=  15.00')
         write(itout,1994)
1994     format(' XMI= 15.00, XPMI= 50.00, XMF= 15.00, XPMF= 25.00')
         write(itout,1995)
1995     format(' DPM=  90.0, DWM=  50.00, DPP= 90.00')
         write(itout,1996)
1996     format(' DPMI= 90.0, DPMF= 35.00, DWMI= 50.0, DWMF= 200.0')
c write trace file (beam line part)
c
c deal with negative drifts (typical within DTLs)
c if the current drift is negative, then see if there is a drift before or
c after it. If so, add the 2 drifts, devide by 2 and put this result in both
c drifts
         do k=2,kt3t-1
           trdrifa=0.
           trdrifb=0.
           tif=trace3t(k)
           tifb=trace3t(k-1)
           tifa=trace3t(k+1)
           if(tif(12:13).eq.' 1' .and. tifb(12:13).eq.' 1') then
             read(tif(28:50),*) trdrift
             if(trdrift.lt.0.) then
               read(tifb(28:50),*) trdrifb
               if(trdrifb.gt.0.) then
                 aver=(trdrifb+trdrift)/2.
                 trdrifb=aver
                 trdrift=aver
c                 write(6,*) k,'beforeb',tif(12:13),trdrift,trdrifb
                 write(fdrft,555) trdrift
555              format(f12.6)
                 tif(28:40)=fdrft
                 write(fdrft,555) trdrifb
                 tifb(28:40)=fdrft
                 trace3t(k)=tif
                 trace3t(k-1)=tifb
               endif
             endif
           elseif(tif(12:13).eq.' 1' .and. tifa(12:13).eq.' 1') then
             read(tif(28:50),*) trdrift
             if(trdrift.lt.0.) then
               read(tifa(28:50),*) trdrifa
               if(trdrifa.gt.0.) then
                 aver=(trdrifa+trdrift)/2.
                 trdrifa=aver
                 trdrift=aver
c                 write(6,*) k,'aftera',tif(12:13),trdrift,trdrifa
                 write(fdrft,555) trdrift
                 tif(28:40)=fdrft
                 write(fdrft,555) trdrifa
                 tifa(28:40)=fdrft
                 trace3t(k)=tif
                 trace3t(k+1)=tifa
               endif
             endif
           endif
         enddo
         do k=1,kt3t
           write(itout,'(A)') trace3t(k)
         enddo
         write(itout,'(A)') ' $END'
       endif
1010   CONTINUE
       call daves
       call eugwrt
       call cpu_time(exfin)
       exfin=exfin-exstrt
       call mytime(iitime)
       write(6,*)
       if (mg) then
c using MINGW style gfortran format: 03/30/10 20:51:06 (10 is 2010)
         text(1:4)='    '
         if(iitime(1:2).eq.'01')text(5:7)='Jan'
         if(iitime(1:2).eq.'02')text(5:7)='Feb'
         if(iitime(1:2).eq.'03')text(5:7)='Mar'
         if(iitime(1:2).eq.'04')text(5:7)='Apr'
         if(iitime(1:2).eq.'05')text(5:7)='May'
         if(iitime(1:2).eq.'06')text(5:7)='Jun'
         if(iitime(1:2).eq.'07')text(5:7)='Jul'
         if(iitime(1:2).eq.'08')text(5:7)='Aug'
         if(iitime(1:2).eq.'09')text(5:7)='Sep'
         if(iitime(1:2).eq.'10')text(5:7)='Oct'
         if(iitime(1:2).eq.'11')text(5:7)='Nov'
         if(iitime(1:2).eq.'12')text(5:7)='Dec'
         text(8:8)=' '
         text(9:10)=iitime(4:5)
         text(11:13)=' 20'
         text(14:15)=iitime(7:8)
         text(16:19)=' at '
         text(20:27)=iitime(10:17)
       else
cet2010s
         text(1:11)=iitime(1:11)
         text(12:15)=iitime(21:24)
         text(16:19)=' at '
         text(20:27)=iitime(12:19)
       endif
       write(12,*) 'Stopped on ',text(1:27)
       write(16,*) 'Stopped on ',text(1:27)
       write(6,'(A11,A27)') 'Stopped on ',text(1:27)
       write(6,'(A12,F14.6,A4)') 'Executed in ',exfin,' sec'
cet2010e
       write(16,101)
c      close input files
       close(10)
       close(20)
       close(22)
       close(25)
       close(27)
       close(55)
c  close output files
c   files actived in the code
       close(16)
       close(12)
       close(71)
       close(50)
       close(61)
       close(60)
       close(11)
       close(66)
       close(75)
       close(49)
c  files not actived in the code
c  SV 24/08/2015
       close(13)
cold       close(14)
comment       close(15)
comment       close(17)
comment       close(18)
comment       close(19)
comment       close(21)
comment       close(49)
       close(70)
       stop
111    FORMAT('STOP ON KEY: ',A8,' (invalid key)')
100    FORMAT(/,40X,' STOP on key  : ',A8,//)
101    FORMAT('*******************************************************',
     1 '*****************')
       END
       SUBROUTINE mytime(iitime)
       implicit real*8 (a-h,o-z)
       character iitime*30
       integer*8 inttim
       inttim=time8()
       iitime=ctime(inttim)
       return
       end
c *******************************************************************
c       SUBROUTINE rmami
c activate the time of flight for bunchers, cavities and accelerating
c gaps
c - indic and icor:integer flags
c - itvol and imamin: logical flags
c - itvol = true =>   the time of flight is activated, otherwise
c - itvol = false =>  the time of flight is passive
c - imamin = true =>  adjustments are automatically made on the phase
c                     of bunchers, cavities and acc. gaps
c - imamin = false => no adjustments on the phase of accelerating
c                     elements
c - indic = 0 => itvol = true, in this case the time of flight is
c                activited for accelerating elements
c - indic (<>) 0 => itvol = false and imamin = false
c - icor = 0 => imamin = false
c - icor (<>) 0 => imamin = true
c *******************************************************************
       SUBROUTINE rmami
       implicit real*8 (a-h,o-z)
       parameter (maxcell1=3000)
       common/tapes/in,ifile,meta
       common/dyn/tref,vref
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/itvole/itvol,imamin
       common/tofev/ttvols
       logical itvol,imamin
       fcpi=fh*180./pi
       ttvols=0.
       read(in,*) indic,icor
       if(indic.eq.0) then
        itvol=.true.
        ttvols=tref
        write(16,10) ttvols*fcpi,davtot,tref*fcpi
       else
        itvol=.false.
        write(16,*) 'time of flight passive '
       endif
10     FORMAT(' ** time of flight activated at: ',e12.5,
     * ' deg at position: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg')
       imamin = .false.
       if(itvol.and.icor.ne.0) imamin = .true.
       if(imamin)
     * write(16,*)' Adjustments on phase offset of acc. elements'
       if(.not.imamin)
     *  write(16,*) ' No adjustments on phase offset'
       RETURN
       END
       SUBROUTINE kick
       implicit real*8 (a-h,o-z)
C     ..................................
C     ALIGNMENT DEFAULTS :
C       HORIZONTAL : XL(cm)
C       VERTICAL   : YL(cm)
C     ..................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DESA/XL,YL,IFLAG
       common/faisc/f(10,iptsz),imax,ngood
       WRITE(16,100) XL,YL
100    FORMAT(/,5X,' KICK x(cm) y(cm): ',2(e12.5,2x),/)
C         CM-MRD
       DO II=1,ngood
         if(iflag.ne.1) then
           F(2,II)=F(2,II) + XL
           F(4,II)=F(4,II) + YL
         else
           F(2,II)=F(2,II) - XL
           F(4,II)=F(4,II) - YL
         endif
       enddo
       RETURN
       END
       SUBROUTINE shuffle
       implicit real*8 (a-h,o-z)
c  ...............................................................................
c  Reshuffles f(i,j) array so that the "good" particles
c  are on top of the stack.  The number of "good" particles
c  (ngood) is passed back.
c  ...............................................................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/mcs/imcs,ncstat,cstat(20)
       common/macro/ratei
       k=ngood
       ngood1=0
  1    if(ngood1.ge.k)  goto 5
       if(f(8,ngood1+1).eq.1.) goto 4
       do j=1,10
         c=f(j,ngood1+1)
         f(j,ngood1+1)=f(j,k)
         f(j,k)=c
       enddo
       k=k-1
       goto 1
  4    ngood1=ngood1+1
       goto 1
  5    continue
       ngood=ngood1
       if(ngood.eq.0) then
         ratei=0.
       else
         ratei=float(imax)/float(ngood)
       endif
c*ets*2015/05/20
       ncstat=1
       cstat(1)=f(9,1)
       do j=2,ngood
         mcstat=0
         do k=1,ncstat
           if(f(9,j).eq.cstat(k)) then
             mcstat=1
           endif
         enddo
         if(mcstat.eq.0) then
           ncstat=ncstat+1
           cstat(ncstat)=f(9,j)
         endif
       enddo
       write(16,*) 'Number of charge states after shuffle: ',ncstat
       write(16,*) 'Charge states: ',(cstat(j),j=1,ncstat)
       imcs=0
       if(ncstat.gt.1) imcs=1
c*ete*2015/05/20
comment       write(16,*) '***reshuffle: ',ngood,' good particles'
       if(ngood.lt.10) then
         write(16,*)'Less than 10 particles left, statistics too low'
         write(6,*)
         write(6,*)'Less than 10 particles left, statistics too low'
         stop
       endif
       return
       end
       FUNCTION xitl0(GAMI,GAMS,BETR,SAPHI,QQC)
c   ....................................................................
c       called by RESTAY and ETGAP  dynamics computations
c  ......................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *              AA,BB,CC,DD,EE,PCREST,SQCTTF
       FH0=FH/VL
       CGI=QQC/XMAT
       BETI=SQRT(1.-1./(GAMI*GAMI))
       BETS=SQRT(1.-1./(GAMS*GAMS))
       XK1=FH0/BETI
       XK2=FH0/BETS
       XKM=FH0/BETR
       TILTA2=PHSLIP/(2.*EQVL)
       PAVPH=1./10.*(XK1- XK2)*EQVL
     X +(XK1-XKM)*ASDL
       XKC1=-FH0/(BETI**3 * GAMI**3)
       XKC2=-FH0/(BETS**3 * GAMS**3)
       DO I=1,2
         PHIT10=SAPHI-PHSLIP/2.+PAVPH
         DAZ0=COS(PHIT10)*TILTA2
         DBZ0=SIN(PHIT10)*TILTA2
         DGZ0=CGI*(TK*DAZ0-SK*DBZ0)
         DGZ0=DGZ0/SIN(PHSLIP/2.)
         XKP1=XKC1*DGZ0
         PHIT11=SAPHI+PHSLIP/2.+PAVPH
         DAZ1=COS(PHIT11)*TILTA2
         DBZ1=SIN(PHIT11)*TILTA2
         DGZ1=CGI*(TK*DAZ1-SK*DBZ1)
         DGZ1=DGZ1/SIN(PHSLIP/2.)
         XKP2=XKC2*DGZ1
         PAVPH=1./10.*(XK1- XK2)*EQVL
     X         +(XK1-XKM)*ASDL
         PAVPH=PAVPH+(XKP1+XKP2)*EQVL**2 / 120.
       ENDDO
       XK11=XK1-XKM
       XK22=XK2-XKM
       AA=XK11
       BB=XKP1/2.
       CCL1=-(4.*XK22+6.*XK11)/(EQVL**2)
       CCL2=-(3./2.*XKP1-XKP2/2.)/EQVL
       CC=CCL1+CCL2
       DDL1=(7.*XK22+8.*XK11)/(EQVL**3)
       DDL2=(3./2.*XKP1-XKP2)/(EQVL**2)
       DD=DDL1+DDL2
       EEL1=-(3.*XK22+3.*XK11)/(EQVL**4)
       EEL2=-(XKP1/2.-XKP2/2.)/(EQVL**3)
       EE=EEL1+EEL2
       PHIT0=SAPHI+PAVPH
       GIT = CGI * (
     X TK*COS(PHIT0)-SK*SIN(PHIT0))
       XITL0=GAMI+GIT
       RETURN
       END
       FUNCTION xitl2(GAMI,GAMS,BETR,SAPHI,QQC)
c    ..................................................................
c     called by RESTAY and ETGAP
C       INTEGRAL of S ( EZG * Z/(BETA*GAMA)**3  *DZ)
C       PHASE JUMP
c    ....................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       DIMENSION H(8),T(8)
       DATA H /.101228536,.222381034,.313706646,.362683783,
     1  .362683783,.313706646,.222381034,.101228536/
       DATA T /-.960289856,-.796666477,-.525532409,
     1  -.183434642,.183434642,.525532409,
     2  .796666477,.960289856/
       FH0=FH/VL
       CGI=QQC/XMAT
       XITL2=0.
       BETI=SQRT(1.-1./(GAMI*GAMI))
       BETS=SQRT(1.-1./(GAMS*GAMS))
       XK1=FH0/BETI
       XK2=FH0/BETS
       XKM=FH0/BETR
       TILTA2=PHSLIP/(2.*EQVL)
       DO I=1,8
         XCC= EQVL*(1.+T(I))/2.
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
         GIT = CGI * (
     X   TK*COS(PHIT0)-SK*SIN(PHIT0)  )
         GI=GAMI+GIT *
     X   SIN(XCC*TILTA2)/SIN(PHSLIP/2.)
         BI=SQRT(1.-1./(GI*GI))
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         DAZ=COS(PHIT1)*TILTA2
         DBZ=SIN(PHIT1)*TILTA2
         DGZ=CGI*(TK*DAZ-SK*DBZ)
         DGZ=DGZ/SIN(PHSLIP/2.)
         EZ = XMAT/QQC * DGZ
         XINT=1./(BI*GI)
         XCC1=XCC+ASDL
         XITL2=XITL2+H(I)*XINT**3 *XCC1*EZ
       ENDDO
       XITL2=XITL2/2. *EQVL
       RETURN
       END
       FUNCTION xitl3(GAMI,GAMS,BETR,NIT,SAPHI,QQC)
c  .....................................................................
c     called by RESTAY and ETGAP
C       INTEGRAL  S ( EZG * Z/(BETA*GAMA)**3  *DZ)
C       PHASE AND ENERGY AT THE MIDDLE OF THE GAP
c  .....................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/gaus13/H(13),T(13)
       common/sgcos/xkpc
       FH0=FH/VL
       CGI=QQC/XMAT
       XITL3=0.
       BETI=SQRT(1.-1./(GAMI*GAMI))
       BETS=SQRT(1.-1./(GAMS*GAMS))
       XK1=FH0/BETI
       XK2=FH0/BETS
       XKM=FH0/BETR
       XK11=XK1-XKM
       TILTA2=PHSLIP/(2.*EQVL)
       DO I=1,13
         XCC= EQVL*(1.+T(I))/2.
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
         GIT = CGI * (TK*COS(PHIT0)-SK*SIN(PHIT0))
         GI=GAMI+GIT * SIN(XCC*TILTA2)/SIN(PHSLIP/2.)
         BI=SQRT(1.-1./(GI*GI))
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         DAZ=COS(PHIT1)*TILTA2
         DBZ=SIN(PHIT1)*TILTA2
         DGZ=CGI*(TK*DAZ-SK*DBZ)
         DGZ=DGZ/SIN(PHSLIP/2.)
         EZ = XMAT/QQC * DGZ
         XINT=1./(BI*GI)
         XCC1=XCC+ASDL
         XITL3=XITL3+H(I)*XINT**3 *XCC1*EZ
C     ENERGY AND PHASE AT THE MIDDLE OF THE GAP
         IF(NIT.EQ.3.AND.I.EQ.7) THEN
c          ENMIL=XMAT*(GI-GAMI)
           ENMIL=XMAT*(GI-1.)
           VAPMI=(XK11*ASDL+SAPHI+XKM*XCC1+AA*XCC+
     *     BB*XCC*XCC+CC*XCC**3+DD*XCC**4+EE*XCC**5)*180./PI
         ENDIF
       ENDDO
       XITL3=XITL3/2. *EQVL
       RETURN
       END
       SUBROUTINE xtypl2(GAMI,SAPHI,QSC,DCG)
c  .....................................................................
c     called by RESTAY and ETGAP
c    integrals of the second derivative of the functions HA0(z) and HB0(z)
C       ZB = Z+ASDL
c  .....................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/THAD2/H0AKI,H0AKIM,H0AKM,H0BKI,H0BKIM,H0BKM,
     X              H1AKI,H1AKIM,H1AKM,H1BKI,H1BKIM,H1BKM
       common/gaus17/H1(17),T1(17)
       FH0=FH/VL
       CGI=QSC/XMAT
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
C     circular functions in cos
       H0AKI=0.
       H0AKIM=0.
       H0AKM=0.
       H1AKI=0.
       H1AKIM=0.
       H1AKM=0.
C     circular functions in sin
       H0BKI=0.
       H0BKIM=0.
       H0BKM=0.
       H1BKI=0.
       H1BKIM=0.
       H1BKM=0.
c
       DTILK=EQVL
       TILTA2=PHSLIP/(2.*EQVL)
       CGAM10=((GAMI*GAMI-1.)**3)/(FH0*FH0)
       DGAM10=GAMI*((GAMI*GAMI-1.)**2)/(FH0*FH0)
       PHCRTK=(T1K*SK-S1K*TK)/(TK*TK+SK*SK)
       DPHC1=(T2K*SK-S2K*TK)/(TK*TK+SK*SK)
       DPHC2=(T1K*TK+S1K*SK)*(T1K*SK-S1K*TK)/
     X ((TK*TK+SK*SK)**2)
       DPHCRTK=DPHC1-2.*DPHC2
       BIM1=BETI
       GAKM1=0.
       GAIT=GAMI
       do i=1,17
         XCC= EQVL*(1.+H1(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C        Function GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS =  SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GITC=GIT * GIS
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C        FONCTION DERIVE G0(Z) RELATIF A K MOYEN
         PHIT0K=-DTILK*(1.-XCC/EQVL)/2.
C        DERIVE PREMIERE FONGTION GAMMA
         IF(PHSLIP.NE.0.) THEN
           GIC =  COS(XCC*TILTA2)
           GAK1=DTILK*COS(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAK2=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAK3=DTILK*COS(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAK =CGI * SQCTTF* (-GAK1-GAK2+GAK3)
         ELSE
           GAK1=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAK =-CGI * SQCTTF* GAK1
         ENDIF
C      Second derivative of GAMMA(z)
         if(gi.ne.gait) then
           DGAK=(GAK-GAKM1)/(GI-GAIT) *GAK
         else
           dgak=0.
         endif
         GAKM1=GAK
         GAIT=GI
         XCC1=XCC+ASDL
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         PHTZ0=(XCC/EQVL-.5)*DTILK
         DPHTZ0=DTILK/EQVL
         PHCRZ0=(PHTZ0-PHCRTK)
C      INTEGRALES FONCTIONS HAKI(Z) (a multipler par (k1i-k10)**2 )
         HAKI1=SQCTTF*COS(PHIT1-PCREST)/((GI*GI-1.)**2.5)
         HAKI2=SQCTTF*COS(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
         HAKI3=SQCTTF*COS(PHIT1-PCREST)*GI*GI/((GI*GI-1.)**3.5)
C        n=0
         H0AKI=H0AKI+T1(I)*CGAM10*(-3.*HAKI1+15*HAKI3)-
     X         T1(I)*DGAM10*9.*HAKI2
C        n=1
         H1AKI=H1AKI+T1(I)*CGAM10*XCC1*(-3.*HAKI1+15.*HAKI3)-
     X         T1(I)*DGAM10*9.*HAKI2*XCC1
C  Integral of HAKIM(Z) (to be multiplied by (k1i-k10)*(kmi-km0) )
         HAKIM1=SQCTTF*COS(PHIT1-PCREST)*GAK/((GI*GI-1.)**2.5)
         HAKIM2=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GAK/((GI*GI-1.)**3.5)
         HAKIM2=HAKIM2
         HAKIM3=SQCTTF*SIN(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
C        n=0
         H0AKIM=H0AKIM+T1(I)*SQRT(CGAM10)*(6.*HAKIM1-30.*HAKIM2 -
     X          3.*PHCRZ0*HAKIM3)
C        n=1
         H1AKIM=H1AKIM+T1(I)*SQRT(CGAM10)*(6.*HAKIM1-30.*HAKIM2 -
     X          3.*PHCRZ0*HAKIM3)*XCC1
C  Integral of HAKM(Z) (to be multiplied by (kmi-km0)**2  )
         HAKM1=SQCTTF*COS(PHIT1-PCREST)*GAK*GAK/((GI*GI-1.)**2.5)
         HAKM2=SQCTTF*COS(PHIT1-PCREST)*DGAK*GI/((GI*GI-1.)**2.5)
         HAKM3=SQCTTF*COS(PHIT1-PCREST)*GAK*GAK*GI*GI/
     X         ((GI*GI-1.)**3.5)
         HAKM4=SQCTTF*SIN(PHIT1-PCREST)*GAK*GI/
     X         ((GI*GI-1.)**2.5)
         HAKM5=SQCTTF*COS(PHIT1-PCREST)/
     X         ((GI*GI-1.)**1.5)
         HAKM6=SQCTTF*SIN(PHIT1-PCREST)/
     X         ((GI*GI-1.)**1.5)
C        n=0
         H0AKM=H0AKM+T1(I)*(-3.*HAKM1-3.*HAKM2+15.*HAKM3 +
     X         3.*PHCRZ0*HAKM4 -
     X         PHCRZ0*PHCRZ0*HAKM5+DPHCRTK*HAKM6)
C        n=1
         H1AKM=H1AKM+T1(I)*XCC1*(-3.*HAKM1-3.*HAKM2+15.*HAKM3 +
     X         3.*PHCRZ0*HAKM4 -
     X         PHCRZ0*PHCRZ0*HAKM5+DPHCRTK*HAKM6)
C
C      INTEGRALES FONCTIONS HBKI(Z) (a multipler par (k1i-k10)**2 )
         HBKI1=SQCTTF*SIN(PHIT1-PCREST)/((GI*GI-1.)**2.5)
         HBKI2=SQCTTF*SIN(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
         HBKI3=SQCTTF*SIN(PHIT1-PCREST)*GI*GI/((GI*GI-1.)**3.5)
C        n=0
         H0BKI=H0BKI+T1(I)*CGAM10*(-3.*HBKI1+15.*HBKI3) -
     X               T1(I)*DGAM10*9.*HBKI2
C        n=1
         H1BKI=H1BKI+T1(I)*CGAM10*XCC1*(-3.*HBKI1+15.*HBKI3) -
     X              T1(I)*DGAM10*9.*HBKI2*XCC1
C  Integral of HBKIM(Z) (to be multiplied by (k1i-k10)*(kmi-km0) )
         HBKIM1=SQCTTF*SIN(PHIT1-PCREST)*GAK/((GI*GI-1.)**2.5)
         HBKIM2=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GAK/((GI*GI-1.)**3.5)
         HBKIM2=HBKIM2
         HBKIM3=SQCTTF*COS(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
C        n=0
         H0BKIM=H0BKIM+T1(I)*SQRT(CGAM10)*(6.*HAKIM1-30.*HAKIM2 +
     X          3.*PHCRZ0*HAKIM3)
C        n=1
         H1BKIM=H1BKIM+T1(I)*SQRT(CGAM10)*(6.*HAKIM1-30.*HAKIM2 +
     X          3.*PHCRZ0*HAKIM3)*XCC1
C  Integral of HBKM(Z) (to be multiplied by (kmi-km0)**2  )
         HBKM1=SQCTTF*SIN(PHIT1-PCREST)*GAK*GAK/((GI*GI-1.)**2.5)
         HBKM2=SQCTTF*SIN(PHIT1-PCREST)*DGAK*GI/((GI*GI-1.)**2.5)
         HBKM3=SQCTTF*SIN(PHIT1-PCREST)*GAK*GAK*GI*GI/
     X         ((GI*GI-1.)**3.5)
         HBKM4=SQCTTF*COS(PHIT1-PCREST)*GAK*GI/
     X         ((GI*GI-1.)**2.5)
         HBKM5=SQCTTF*SIN(PHIT1-PCREST)/
     X         ((GI*GI-1.)**1.5)
         HBKM6=SQCTTF*COS(PHIT1-PCREST)/
     X         ((GI*GI-1.)**1.5)
C       n=0
         H0BKM=H0BKM+T1(I)*(-3.*HBKM1-3.*HBKM2+15.*HBKM3 -
     X         3.*PHCRZ0*HBKM4 -
     X         PHCRZ0*PHCRZ0*HBKM5-DPHCRTK*HBKM6)
C       n=1
         H1BKM=H1BKM+T1(I)*XCC1*(-3.*HBKM1-3.*HBKM2+15.*HBKM3-
     X         3.*PHCRZ0*HBKM4 -
     X         PHCRZ0*PHCRZ0*HBKM5-DPHCRTK*HBKM6)
       enddo
200    CONTINUE
C    integrals in cos
       H0AKI=H0AKI/2. *EQVL
       H0AKIM=H0AKIM/2. *EQVL
       H0AKM=H0AKM/2. *EQVL
       H1AKI=H1AKI/2. *EQVL
       H1AKIM=H1AKIM/2. *EQVL
       H1AKM=H1AKM/2. *EQVL
C    integrals in sin
       H0BKI=H0BKI/2. *EQVL
       H0BKIM=H0BKIM/2. *EQVL
       H0BKM=H0BKM/2. *EQVL
       H1BKI=H1BKI/2. *EQVL
       H1BKIM=H1BKIM/2. *EQVL
       H1BKM=H1BKM/2. *EQVL
       RETURN
       END
       SUBROUTINE xtyplp1(GAMI,SAPHI,QSC,DCG)
c   ...............................................................
c     called by RESTAY and ETGAP
C   Integrals of functions  HA0(Z) and HB0(Z)
c   Integrals of the first and second derivative of HA0(Z) et HB0(Z)(with PH01)
c   Integrals of the third derivative of HA0(Z) et HB0(Z)(with PH01)
C       zb = z+ASDL
c   ...............................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/JACOB/GAKS,GAPS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/TYPLP1/YH1P1,YH2P1,HAPI,HBPI
       COMMON/TYPLP2/HAPPI,HBPPI
       common/gaus17/H1(17),T1(17)
       FH0=FH/VL
       CGI=QSC/XMAT
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
       DTILK=EQVL
C     circular functions in cos
       YH1P1=0.
       HAPI=0.
       HAPPI=0.
C     circular functions in sin
       YH2P1=0.
       HBPI=0.
       HBPPI=0.
       TILTA2=PHSLIP/(2.*EQVL)
       IF(PHSLIP.NE.0.) DESY=PHSLIP/SIN(PHSLIP/2.)
       IF(PHSLIP.EQ.0.) DESY=2.
       do i=1,17
         XCC= EQVL*(1.+H1(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C        FONCTION GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS =  SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C        FONCTION DERIVE G0(Z) RELATIF A PH01
         IF(PHSLIP.NE.0.) THEN
           GAP=-CGI*SQCTTF*SIN(PHIT0-PCREST)*GIS/PHSLIP
C          DERIVE TROISIEME
           DDGAP=-GAP
         ELSE
           GAP = SIN(PHIT0-PCREST)*GIS
           GAP = -CGI * SQCTTF* GAP
           DDGAP=-GAP
         ENDIF
         IF(I.EQ.17) GAPS=GAP
C        Second derivative of G0(Z) in relation with  PH01
         IF(PHSLIP.NE.0.) THEN
           DGAP=-CGI*SQCTTF*COS(PHIT0-PCREST)*GIS/PHSLIP
         ELSE
           DGAP = COS(PHIT0-PCREST)*GIS
           DGAP = -CGI * SQCTTF*DGAP
         ENDIF
C        INTEGRALES FONCTIONS HA0(Z) et HB0(Z)
         XINT=1./(BI*BI*BI*GI*GI*GI)
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         PHTZ0=(XCC/EQVL-.5)*DTILK
C        INTEGRALES DERIVES HA0(Z)
         DHA01=SQCTTF*COS(PHIT1-PCREST)*GI*GAP/((GI*GI-1.)**2.5)
         DHA02=0.
C        n=1
         YH1P1=YH1P1+T1(I)*XCC1*(-6.*DHA01-2.*DHA02)
C        INTEGRALES DERIVES  HB0(Z)
         DHB01=SQCTTF*SIN(PHIT1-PCREST)*GI*GAP/((GI*GI-1.)**2.5)
         DHB02=0.
C        n=1
         YH2P1=YH2P1+T1(I)*XCC1*(-6.*DHB01+2.*DHB02)
C        INTEGRALES HAPI(Z) (multiplies par (ph1i-ph10)**2 )
         HAPI1=SQCTTF*COS(PHIT1-PCREST)*GAP*GAP/((GI*GI-1.)**2.5)
         HAPI2=SQCTTF*COS(PHIT1-PCREST)*GI*DGAP/((GI*GI-1.)**2.5)
         HAPI3=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GAP*GAP/
     X         ((GI*GI-1.)**3.5)
C        n=1
         HAPI=HAPI+T1(I)*XCC1*(-3.*HAPI1-3.*HAPI2+15.*HAPI3)
C        INTEGRALES HAPPI(Z) (multiplies par (ph1i-ph10)**3 )/3
         HAPPI1=SQCTTF*COS(PHIT1-PCREST)*GAP*DGAP/
     X         ((GI*GI-1.)**2.5)
         HAPPI2=SQCTTF*COS(PHIT1-PCREST)*GAP*GAP*GAP*GI/
     X         ((GI*GI-1.)**3.5)
         HAPPI3=SQCTTF*COS(PHIT1-PCREST)*GAP*DGAP/
     X         ((GI*GI-1.)**2.5)
         HAPPI4=SQCTTF*COS(PHIT1-PCREST)*GI*DDGAP/
     X         ((GI*GI-1.)**2.5)
         HAPPI5=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GAP*DGAP/
     X         ((GI*GI-1.)**3.5)
         HAPPI6=SQCTTF*COS(PHIT1-PCREST)*GI*GAP**3/
     X         ((GI*GI-1.)**3.5)
         HAPPI7=SQCTTF*COS(PHIT1-PCREST)*GI*GI*DGAP*GAP/
     X         ((GI*GI-1.)**3.5)
         HAPPI8=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GI*GAP*GAP*GAP/
     X         ((GI*GI-1.)**4.5)
C        n=1
         HAPPI=HAPPI+T1(I)*XCC1*(-6.*HAPPI1+15.*HAPPI2-3.*HAPPI3 -
     X   3.*HAPPI4 +15.*HAPPI5 +30.*HAPPI6 +30.*HAPPI7 -105.*HAPPI8)
C
C        INTEGRAL OF HBPI(Z) (to be multiplied by (ph1i-ph10)**2 )
         HBPI1=SQCTTF*SIN(PHIT1-PCREST)*GAP*GAP/((GI*GI-1.)**2.5)
         HBPI2=SQCTTF*SIN(PHIT1-PCREST)*GI*DGAP/((GI*GI-1.)**2.5)
         HBPI3=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GAP*GAP/
     X         ((GI*GI-1.)**3.5)
C        n=1
         HBPI=HBPI+T1(I)*XCC1*(-3.*HBPI1-3.*HBPI2+15.*HBPI3 )
C        INTEGRALS of HBPPI(Z) (to be multiplied by(ph1i-ph10)**3 )/3
         HBPPI1=SQCTTF*SIN(PHIT1-PCREST)*GAP*DGAP/
     X         ((GI*GI-1.)**2.5)
         HBPPI2=SQCTTF*SIN(PHIT1-PCREST)*GAP*GAP*GAP*GI/
     X         ((GI*GI-1.)**3.5)
         HBPPI3=SQCTTF*SIN(PHIT1-PCREST)*GAP*DGAP/
     X         ((GI*GI-1.)**2.5)
         HBPPI4=SQCTTF*SIN(PHIT1-PCREST)*GI*DDGAP/
     X         ((GI*GI-1.)**2.5)
         HBPPI5=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GAP*DGAP/
     X         ((GI*GI-1.)**3.5)
         HBPPI6=SQCTTF*SIN(PHIT1-PCREST)*GI*GAP**3/
     X         ((GI*GI-1.)**3.5)
         HBPPI7=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*DGAP*GAP/
     X         ((GI*GI-1.)**3.5)
         HBPPI8=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GI*GAP*GAP*GAP/
     X         ((GI*GI-1.)**4.5)
C        n=1
         HBPPI=HAPPI+T1(I)*XCC1*(-6.*HBPPI1+15.*HBPPI2-3.*HBPPI3 -
     X   3.*HBPPI4 +15.*HBPPI5 +30.*HBPPI6 +30.*HBPPI7 -105.*HBPPI8)
       enddo
200    CONTINUE
C    INTEGRALS to be multiplied by cos
       YH1P1=YH1P1/2. *EQVL
       HAPI=HAPI/2. *EQVL
       HAPPI=HAPPI/2. *EQVL
C    INTEGRALS to be multiplied by sin
       YH2P1=YH2P1/2. *EQVL
       HBPI=HBPI/2. *EQVL
       HBPPI=HBPPI/2. *EQVL
       RETURN
       END
       SUBROUTINE xtylpk(GAMI,SAPHI,QSC,DCG)
c  .......................................................................
c   called by RESTAY and ETGAP
C   INTEGRALS of HA0(Z) and HB0(Z)
C   INTEGRALS of the second derivative of HA0(Z) and HB0(Z)
C       zb = z+ASDL
c  ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/TYPLPK/YH10PK,YH11PK,YH20PK,YH21PK
       common/gaus17/H1(17),T1(17)
       FH0=FH/VL
       CGI=QSC/XMAT
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
       DTILK=EQVL
C     CIRCULAIRES EN COS
       YH10PK=0.
       YH11PK=0.
C     CIRCULAIRES EN SIN
       YH20PK=0.
       YH21PK=0.
       TILTA2=PHSLIP/(2.*EQVL)
       IF(PHSLIP.NE.0.) DESY=PHSLIP/SIN(PHSLIP/2.)
       IF(PHSLIP.EQ.0.) DESY=2.
       CGAM10=((GAMI*GAMI-1.)**1.5)/FH0
       PHCRTK=(T1K*SK-S1K*TK)/(TK*TK+SK*SK)
       do i=1,17
         XCC= EQVL*(1.+H1(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C        FONCTION GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS =  SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C        Derivative of G0(Z) in relation with: PH01
         IF(PHSLIP.NE.0.) THEN
           GAP=-CGI*SQCTTF*SIN(PHIT0-PCREST)*GIS/PHSLIP
         ELSE
           GAP = SIN(PHIT0-PCREST)*GIS
           GAP = -CGI * SQCTTF* GAP
         ENDIF
C        Derivative of G0(Z) in relation with the equivalent K
         PHIT0K=-DTILK*(1.-XCC/EQVL)/2.
         IF(PHSLIP.NE.0.) THEN
           GIC =COS(XCC*TILTA2)
           GAK1=DTILK*COS(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAK2=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAK3=DTILK*COS(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAK =CGI * SQCTTF* (-GAK1-GAK2+GAK3)
         ELSE
           GAK1=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAK =-CGI * SQCTTF* GAK1
         ENDIF
C        derivative of G0(Z) in relation with the equivalent k and the phase PHI
         IF(PHSLIP.NE.0.) THEN
           GIC = COS(XCC*TILTA2)
           GAKP1=DTILK*SIN(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAKP2=COS(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAKP3=DTILK*SIN(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAKP =CGI * SQCTTF* (GAKP1-GAKP2-GAKP3)
         ELSE
           GAKP1 = COS(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAKP = -CGI * SQCTTF* GAKP1
         ENDIF
C        INTEGRALS of HA0(Z) and HB0(Z)
         XINT=1./(BI*BI*BI*GI*GI*GI)
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
C        INTEGRALES DERIVE HA0(Z) HB0(Z)  k0 ET PHI
         DHA01=SQCTTF*COS(PHIT1-PCREST)*GAP/((GI*GI-1.)**2.5)
         DHA02=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GAP/
     X         ((GI*GI-1.)**3.5)
         DHB01=SQCTTF*SIN(PHIT1-PCREST)*GAP/((GI*GI-1.)**2.5)
         DHB02=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GAP/
     X         ((GI*GI-1.)**3.5)
C        n=1
         YH10PK=YH10PK+T1(I)*XCC1*CGAM10*(3.*DHA01-15.*DHA02)
         YH20PK=YH20PK+T1(I)*XCC1*CGAM10*(3.*DHB01-15.*DHB02)
C        Integrals of HA(Z) (multiplied by (ph1i-ph10)*(kmi-km0)
         HAPI1=SQCTTF*COS(PHIT1-PCREST)*GAP*GAK/((GI*GI-1.)**2.5)
         HAPI2=SQCTTF*COS(PHIT1-PCREST)*GI*GAKP/((GI*GI-1.)**2.5)
         HAPI3=SQCTTF*COS(PHIT1-PCREST)*GI*GI*GAP*GAK/
     X         ((GI*GI-1.)**3.5)
C        n=1
         YH11PK=YH11PK+T1(I)*XCC1*(-3.*HAPI1-3.*HAPI2+15.*HAPI3)
C        Integrals of HB(Z) (multiplied by (ph1i-ph10)*(kmi-km0)
         HBPI1=SQCTTF*SIN(PHIT1-PCREST)*GAP*GAK/((GI*GI-1.)**2.5)
         HBPI2=SQCTTF*SIN(PHIT1-PCREST)*GI*GAKP/((GI*GI-1.)**2.5)
         HBPI3=SQCTTF*SIN(PHIT1-PCREST)*GI*GI*GAP*GAK/
     X         ((GI*GI-1.)**3.5)
C        n=1
         YH21PK=YH21PK+T1(I)*XCC1*(-3.*HBPI1-3.*HBPI2+15.*HBPI3)
       enddo
200    CONTINUE
C    Cosinus integrals
       YH10PK=YH10PK/2. *EQVL
       YH11PK=YH11PK/2. *EQVL
C    Sinus integrals
       YH20PK=YH20PK/2. *EQVL
       YH21PK=YH21PK/2. *EQVL
       RETURN
       END
       SUBROUTINE xtypj(GAMI,SAPHI,QSC,DCG)
c  ...............................................................
c    called by RESTAY and ETGAP
C INTEGRALS of g(G)*Ez(**2) *z**n   n=0,1,2
c  ...............................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/TYPJ/YFSK0,YFSK1,YFSK2,YFSP0,YFSP1,YFSP2,
     X               YFSKC0,YFSKC1,YFSKC2,
     X               YFSCK0,YFSCK1,YFSCK2,
     X               YFSCP0,YFSCP1,YFSCP2,
     X               YFS0,YFS1,YFS2
       common/gaus17/H1(17),T1(17)
       FH0=FH/VL
       CGI=QSC/XMAT
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
       DTILK=EQVL
       YFS0=0.
       YFS1=0.
       YFS2=0.
       YFSKC0=0.
       YFSKC1=0.
       YFSKC2=0.
       YFSK0=0.
       YFSK1=0.
       YFSK2=0.
       YFSCK0=0.
       YFSCK1=0.
       YFSCK2=0.
       YFSCP0=0.
       YFSCP1=0.
       YFSCP2=0.
       YFSP0=0.
       YFSP1=0.
       YFSP2=0.
       TILTA2=PHSLIP/(2.*EQVL)
       IF(PHSLIP.NE.0.) DESY=PHSLIP/SIN(PHSLIP/2.)
       IF(PHSLIP.EQ.0.) DESY=2.
       CGAM10=((GAMI*GAMI-1.)**1.5)/FH0
       PHCRTK=(T1K*SK-S1K*TK)/(TK*TK+SK*SK)
       do i=1,17
         XCC= EQVL*(1.+H1(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C        Function GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS =  SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C        Derivative of G0(Z) with regard to the average k
         PHIT0K=-DTILK*(1.-XCC/EQVL)/2.
         IF(PHSLIP.NE.0.) THEN
           GIC =COS(XCC*TILTA2)
           GAK1=DTILK*COS(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAK2=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAK3=DTILK*COS(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAK =CGI * SQCTTF* (-GAK1-GAK2+GAK3)
         ELSE
           GAK1=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAK =-CGI * SQCTTF* GAK1
         ENDIF
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         PHTZ0=(XCC/EQVL-.5)*DTILK
         XINT=(GI*GI+2.)/((GI*GI-1.)**2)
         XFK1=2.*GI*(1.-2.*(GI*GI+2.)/(GI*GI-1.))/((GI*GI-1.)**2)
         HA0=COS(PHIT1-PCREST)
         HB0=SIN(PHIT1-PCREST)
C    n=0
         YFS0=YFS0+T1(I)*HA0*HA0*XINT
         YFSKC0=YFSKC0-T1(I)*HA0*HA0*XFK1*CGAM10
         YFSK0=YFSK0+T1(I)*HA0*HA0*XFK1*GAK
         YFSCK0=YFSCK0-2.*T1(I)*HA0*HB0*XINT*(PHTZ0-PHCRTK)
C    n=1
         YFS1=YFS1+T1(I)*HA0*HA0*XINT*XCC
         YFSKC1=YFSKC1-T1(I)*HA0*HA0*XFK1*CGAM10*XCC
         YFSK1=YFSK1+T1(I)*HA0*HA0*XFK1*GAK*XCC
         YFSCK1=YFSCK1-2.*T1(I)*HA0*HB0*XINT*(PHTZ0-PHCRTK)*XCC
C    n=2
         YFS2=YFS2+T1(I)*HA0*HA0*XINT*XCC*XCC
         YFSKC2=YFSKC2-T1(I)*HA0*HA0*XFK1*CGAM10*XCC*XCC
         YFSK2=YFSK2+T1(I)*HA0*HA0*XFK1*GAK*XCC*XCC
         YFSCK2=YFSCK2-2.*T1(I)*HA0*HB0*XINT*(PHTZ0-PHCRTK)*XCC*XCC
C      Derivative of G0(Z) with regard to PH01
         IF(PHSLIP.NE.0.) THEN
           GAP = -CGI*SQCTTF*SIN(PHIT0-PCREST)*GIS/PHSLIP
         ELSE
           GAP=SIN(PHIT0-PCREST)*GIS
           GAP=-CGI*SQCTTF*GAP
         ENDIF
C    n=0
         YFSP0=YFSP0+T1(I)*XFK1*GAP*HA0*HA0
         YFSCP0=YFSCP0-2.*T1(I)*XINT*HA0*HB0
C    n=1
         YFSP1=YFSP1+T1(I)*XFK1*GAP*HA0*HA0*XCC
         YFSCP1=YFSCP1-2.*T1(I)*XINT*HA0*HB0*XCC
C    n=2
         YFSP2=YFSP2+T1(I)*XFK1*GAP*HA0*HA0*XCC*XCC
         YFSCP2=YFSCP2-2.*T1(I)*XINT*HA0*HB0*XCC*XCC
       enddo
200    CONTINUE
c       IST=IAST
       YFS0=YFS0/2. *EQVL
       YFS1=YFS1/2. *EQVL
       YFS2=YFS2/2. *EQVL
       YFSKC0=YFSKC0/2. *EQVL
       YFSKC1=YFSKC1/2. *EQVL
       YFSKC2=YFSKC2/2. *EQVL
       YFSK0=YFSK0/2. *EQVL
       YFSK1=YFSK1/2. *EQVL
       YFSK2=YFSK2/2. *EQVL
       YFSCK0=YFSCK0/2. *EQVL
       YFSCK1=YFSCK1/2. *EQVL
       YFSCK2=YFSCK2/2. *EQVL
       YFSCP0=YFSCP0/2. *EQVL
       YFSCP1=YFSCP1/2. *EQVL
       YFSCP2=YFSCP2/2. *EQVL
       YFSP0=YFSP0/2. *EQVL
       YFSP1=YFSP1/2. *EQVL
       YFSP2=YFSP2/2. *EQVL
       RETURN
       END
       SUBROUTINE xtypm(GAMI,SAPHI,QSC,DCG)
c   ....................................................................
c    called by RESTAY and ETGAP
C   INTEGRALS of g(G)*z**n  with n=0,1,2
c   ....................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/TYPM/YNSK0,YNSK1,YNSK2,YNSP0,YNSP1,YNSP2,
     X               YNSK0C,YNSK1C,YNSK2C,
     X               YNS0,YNS1,YNS2
       common/gaus17/H1(17),T1(17)
       FH0=FH/VL
       CGI=QSC/XMAT
       DTILK=EQVL
       YNS0=0.
       YNS1=0.
       YNS2=0.
       YNSK0C=0.
       YNSK1C=0.
       YNSK2C=0.
       YNSK0=0.
       YNSK1=0.
       YNSK2=0.
       YNSP0=0.
       YNSP1=0.
       YNSP2=0.
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
       TILTA2=PHSLIP/(2.*EQVL)
       IF(PHSLIP.NE.0.) DESY=PHSLIP/SIN(PHSLIP/2.)
       IF(PHSLIP.EQ.0.) DESY=2.
       CGAM10=((GAMI*GAMI-1.)**1.5)/FH0
       PHCRTK=(T1K*SK-S1K*TK)/(TK*TK+SK*SK)
       do i=1,17
         XCC= EQVL*(1.+H1(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C        Function GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT=CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS=SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C        Derivative of G0(Z) with regard to the average k
         PHIT0K=-DTILK*(1.-XCC/EQVL)/2.
         IF(PHSLIP.NE.0.) THEN
           GIC =COS(XCC*TILTA2)
           GAK1=DTILK*COS(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAK2=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAK3=DTILK*COS(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAK=CGI * SQCTTF* (-GAK1-GAK2+GAK3)
         ELSE
           GAK1=SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAK =-CGI * SQCTTF* GAK1
         ENDIF
         XINT=(GI*GI+2.)/((GI*GI-1.)**2)
         XNK1=2.*GI*(1.-2.*(GI*GI+2.)/(GI*GI-1.))/((GI*GI-1.)**2)
C    n=0
         YNS0=YNS0+T1(I)*XINT
         YNSK0C=YNSK0C-T1(I)*XNK1*CGAM10
         YNSK0=YNSK0+T1(I)*XNK1*GAK
C    n=1
         YNS1=YNS1+T1(I)*XINT*XCC
         YNSK1C=YNSK1C-T1(I)*XNK1*CGAM10*XCC
         YNSK1=YNSK1+T1(I)*XNK1*GAK*XCC
C    n=2
         YNS2=YNS2+T1(I)*XINT*XCC*XCC
         YNSK2C=YNSK2C-T1(I)*XNK1*CGAM10*XCC*XCC
         YNSK2=YNSK2+T1(I)*XNK1*GAK*XCC*XCC
C      DERIVE G0(Z) RELATIF A PH01
         IF(PHSLIP.NE.0.) THEN
           GAP = -CGI*SQCTTF*SIN(PHIT0-PCREST)*GIS/PHSLIP
         ELSE
           GAP=SIN(PHIT0-PCREST)*GIS
           GAP=-CGI*SQCTTF*GAP
         ENDIF
         YNSP0=YNSP0+T1(I)*XNK1*GAP
         YNSP1=YNSP1+T1(I)*XNK1*GAP*XCC
         YNSP2=YNSP2+T1(I)*XNK1*GAP*XCC*XCC
       enddo
200    CONTINUE
       YNS0=YNS0/2. *EQVL
       YNS1=YNS1/2. *EQVL
       YNS2=YNS2/2. *EQVL
       YNSK0C=YNSK0C/2. *EQVL
       YNSK1C=YNSK1C/2. *EQVL
       YNSK2C=YNSK2C/2. *EQVL
       YNSK0=YNSK0/2. *EQVL
       YNSK1=YNSK1/2. *EQVL
       YNSK2=YNSK2/2. *EQVL
       YNSP0=YNSP0/2. *EQVL
       YNSP1=YNSP1/2. *EQVL
       YNSP2=YNSP2/2. *EQVL
       RETURN
       END
       FUNCTION gamci(PHI,PCRESI,GAMI,IST,QSC)
c  ........................................................................
c   called by RESTAY and ETGAP
C       CURRENT GAMMA VALUE (the POSITION IS GIVEN BY IST)
c  ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *               AA,BB,CC,DD,EE,PCREST,SQCTTF
       common/gaus17/H1(17),T1(17)
       GAMCI=0.
       IF(IST.GT.17) RETURN
       CGI=QSC/XMAT
       TILTA2=PHSLIP/(2.*EQVL)
       I=IST
       XCC= EQVL*(1.+H1(I))/2.
       PHIT0=PHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)
C      Function GAMMA (Z)
       IF(PHSLIP.NE.0.) THEN
         GIT = CGI * SQCTTF*COS(PHIT0-PCRESI)/PHSLIP
         GIS =  SIN(XCC*TILTA2)
       ELSE
         GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
         GIS =  XCC/(2.*EQVL)
       ENDIF
       GAMCI=GAMI+GIT * GIS
       RETURN
       END
       SUBROUTINE intga(npt,ireca)
c   ........................................................................
c Calculate the bean self fields acting on each particle (called by SCHERM)
c        Gauss quadrature
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/HERMRR/AFXRR(20),AFYRR(20),AFZRR(20)
       COMMON/SIZR/XRMS3,YRMS3,ZRMS3,ZCGR3
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/ELCG/XCGD,YCGD,ZCGD,XCGR,YCGR,ZCGR
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/INTGRT/ex,ey,ez
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/twcst/epsilon
       COMMON/ECOR/const
       COMMON/NPART/IMAXR
       common/macro/ratei
       LOGICAL ichaes
       DIMENSION UI(6),WI(6)
       DATA (UI(J),J=1,6)/.033765,.169395,.380690,.619310,.830605,
     *      .966234/
       DATA (WI(J),J=1,6)/.085662,.180381,.233957,.233957,.180381,
     *      .085662/
c
c Initialize some constants and variables
c freq. in MHz
       freq=fh*0.5e-06/pi
c qmpart=coul , xrms=meters
c    ireca=0 : first ellipsoid  over nmaxy terms
c    ireca=4 : one ellipsoid over the first term
       xrmsc=0.
       yrmsc=0.
       zrmsc=0.
       xgc=0.
       ygc=0.
       zgc=0.
       if(ireca.eq.0.or.ireca.eq.4) then
         xrmsc=xrms1
         yrmsc=yrms1
         zrmsc=zrms1
         xgc=xcgd
         ygc=ycgd
         zgc=zcgd
       endif
c    ireca=1 : second ellipsoid
       if(ireca.eq.1) then
         xrmsc=xrms2
         yrmsc=yrms2
         zrmsc=zrms2
         xgc=xcgr
         ygc=ycgr
         zgc=zcgr
       endif
c    ireca=2 : third elllipsoid
       if(ireca.eq.2) then
         xrmsc=xrms3
         yrmsc=yrms3
         zrmsc=zrms3
         xgc=0.
         ygc=0.
         zgc=zcgr3
       endif
       qmpart=1.0e-09*beamc/(float(imax)*freq)
comment        const=qmpart*xrmsc*yrmsc*zrmsc/(2.*epsilon)
       qmpart=qmpart*ratei
       const=qmpart/(2.*epsilon)
       dnorm=(xrmsc*yrmsc*zrmsc)**.333333333
       dsq=dnorm*dnorm
       xsq=(xc(npt)-xgc)*(xc(npt)-xgc)
       ysq=(yc(npt)-ygc)*(yc(npt)-ygc)
       zsq=(zc(npt)-zgc)*(zc(npt)-zgc)
       if(ireca.eq.0) zc1=zc(npt)-zgc
       if(ireca.eq.1) zc2=zc(npt)-zgc
       ex=0.
       if(ireca.eq.2) zc2=zc(npt)-zgc
       ey=0.
       ez=0.
c initialize integrals to 0.
c integrate all 3 components (x,y,z)
       DO J=1,6
         a1=xrmsc*xrmsc-dsq+dsq/ui(j)
         a2=yrmsc*yrmsc-dsq+dsq/ui(j)
         a3=zrmsc*zrmsc-dsq+dsq/ui(j)
         t1=xsq/a1
         t2=ysq/a2
         t3=zsq/a3
         txyz=sqrt(t1+t2+t3)
         ff1=drxyz(nmaxy,txyz,ireca)/(ui(j)*ui(j)*sqrt(a1*a2*a3))
         fxn=ff1/a1
         fyn=ff1/a2
         fzn=ff1/a3
         ex=ex+wi(j)*fxn*dsq
         ey=ey+wi(j)*fyn*dsq
         ez=ez+wi(j)*fzn*dsq
       ENDDO
c Field components are in Newton/Coulomb
       ex=ex*const*xc(npt)
       ey=ey*const*yc(npt)
       ez=ez*const*(zc(npt)-zgc)
       RETURN
       END
       SUBROUTINE sizcor(ect,xrms,yrms,zrms,imaxd)
c   ........................................................................
c   Computes the R.M.S. of the bunch at positions of space charge computation
c   .......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/CGRMS/xsum,ysum,zsum
       common/faisc/f(10,iptsz),imax,ngood
       imaxf=0
       imaxx=ngood
       xsum=0.
       ysum=0.
       zsum=0.
       xrmsp=xrms
       yrmsp=yrms
       zrmsp=zrms
       xsqsum=0.
       ysqsum=0.
       zsqsum=0.
       if(imaxd.gt.0) imaxx=imaxd
       do  i=1,imaxx
         xcoup=abs(xc(i)/xrmsp)
         ycoup=abs(yc(i)/yrmsp)
         zcoup=abs(zc(i)/zrmsp)
         if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
           xsum=xsum+xc(i)
           ysum=ysum+yc(i)
           zsum=zsum+zc(i)
           xsqsum=xsqsum+xc(i)*xc(i)
           ysqsum=ysqsum+yc(i)*yc(i)
           zsqsum=zsqsum+zc(i)*zc(i)
           imaxf=imaxf+1
         endif
       enddo
       xsum=xsum/imaxf
       ysum=ysum/imaxf
       zsum=zsum/imaxf
       xsqsum=xsqsum/imaxf
       ysqsum=ysqsum/imaxf
       zsqsum=zsqsum/imaxf
       xrms=SQRT(xsqsum-xsum*xsum)
       yrms=SQRT(ysqsum-ysum*ysum)
       zrms=SQRT(zsqsum-zsum*zsum)
       RETURN
       END
       SUBROUTINE sizrms(imaxd,xrms,yrms,zrms,zmin)
c   ........................................................................
c   partial R.M.S. (called by SCHERM)
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/CGRMS/xsum,ysum,zsum
       common/faisc/f(10,iptsz),imax,ngood
       imaxx=ngood
       xsum=0.
       ysum=0.
       zsum=0.
       if(imaxd.gt.0) then
         do i=1,imaxd
           zc(i)=zc(i)-zmin
         enddo
         imaxx=imaxd
       endif
       xsqsum=0.
       ysqsum=0.
       zsqsum=0.
       do i=1,imaxx
         xsum=xsum+xc(i)
         ysum=ysum+yc(i)
         zsum=zsum+zc(i)
         xsqsum=xsqsum+xc(i)*xc(i)
         ysqsum=ysqsum+yc(i)*yc(i)
         zsqsum=zsqsum+zc(i)*zc(i)
       enddo
       xsum=xsum/imaxx
       ysum=ysum/imaxx
       zsum=zsum/imaxx
       xsqsum=xsqsum/float(imaxx)
       ysqsum=ysqsum/float(imaxx)
       zsqsum=zsqsum/float(imaxx)
       xrms=SQRT(xsqsum-xsum*xsum)
       yrms=SQRT(ysqsum-ysum*ysum)
       zrms=SQRT(zsqsum-zsum*zsum)
       RETURN
       END
       FUNCTION snzt(CC,DD)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       common/gaus17/H1(17),T1(17)
       SNZT=0.
       DO I=1,17
         Z= (CC+DD)/2.+(DD-CC)*H1(I)/2.
         DENZ=DENSZ(NMAZ,Z,0)
         IF(DENZ.LT.0.) DENZ=0.
         SNZT=SNZT+T1(I)*DENZ
       ENDDO
       SNZT=SNZT*(DD-CC)/2.
       RETURN
       END
       FUNCTION snzd(cc,dd)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       common/gaus17/H1(17),T1(17)
       SNZD=0.
       DO I=1,17
         Z= (CC+DD)/2.+(DD-CC)*H1(I)/2.
         DENZ=DENSZ(NMAZ,Z,0)
         IF(DENZ.LT.0.) GO TO 13
         ZZ=(Z-CC)
         SNZD=SNZD+T1(I)*DENZ
       ENDDO
13     CONTINUE
       SNZD=SNZD*(DD-CC)
       RETURN
       END
       FUNCTION vaprz(CC,DD)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       common/gaus17/H1(17),T1(17)
       VAPRZ=0.
       VAR1=0.
       VAR2=0.
       DO I=1,17
         Z= (CC+DD)/2.+(DD-CC)*H1(I)/2.
         DENZ=DENSZ(NMAZ,Z,0)
         IF(DENZ.LT.0.) GO TO 13
         ZZ=(Z-CC)
         VAR1=VAR1+T1(I)*ZZ*ZZ*DENZ
         VAR2=VAR2+T1(I)*DENZ
       ENDDO
13     CONTINUE
       VAPRZ=VAR1/VAR2
       RETURN
       END
       FUNCTION prinz(CC,DD,KAP,ZRMSS1)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       common/gaus17/H1(17),T1(17)
       PRINZ=0.
       K=KAP-1
       DO I=1,17
         DENZ=0.
         Z= (CC+DD)/2.+(DD-CC)*H1(I)/2.
         DENZ=DENSZ(NMAZ,Z,0)
         IF(DENZ.LT.0.) GO TO 13
         ZZ=(Z-CC)/ZRMSS1
         PRINZ=PRINZ+T1(I)*HERM(2*K,ZZ)*DENZ
       ENDDO
13     CONTINUE
       PRINZ=PRINZ*(DD-CC)/2.
       RETURN
       END
       SUBROUTINE rchsom(ZI,ZF,NMAZ)
c   ........................................................................
c    specific function called by SCHERM
c     Look for the top of the partial distributions
c   ........................................................................
       implicit real*8 (a-h,o-z)
       XPZ=ABS((ZF-ZI))/400.
1      CONTINUE
       ZTEST=ABS(ZF-ZI)
       IF(ZTEST.LE.XPZ)RETURN
       Z1=(ZF+ZI)/2.
       Z2=(Z1+ZI)/2.
       T1=DENSZ(NMAZ,Z1,0)
       T2=DENSZ(NMAZ,Z2,0)
C    new interval: Z2,ZF
       IF(T1.GT.T2) ZI=Z2
C    new interval  : ZI,Z1
       IF(T1.LT.T2) ZF=Z1
       GO TO 1
       END
       SUBROUTINE rchsor(AA,BB,CC,DD,ZS)
c   ........................................................................
c    specific function called by SCHERM
c     Look for the top of the partial distributions
c   ........................................................................
       implicit real*8 (a-h,o-z)
       ZI=AA
       ZF=CC
       XPZ=ABS((ZF-ZI))/400.
1      CONTINUE
       ZTEST=ABS(ZF-ZI)
       IF(ZTEST.LE.XPZ) THEN
         ZS=(ZI+ZF)/2.
         RETURN
       ENDIF
       Z1=(ZF+ZI)/2.
       Z2=(Z1+ZI)/2.
       T1=DENDIF(Z1,AA,BB,CC,DD)
       T2=DENDIF(Z2,AA,BB,CC,DD)
C    NEW INTERVAL : Z2,ZF
       IF(T1.GT.T2) ZI=Z2
C    NEW INTERVAL : ZI,Z1
       IF(T1.LT.T2) ZF=Z1
       GO TO 1
       END
       FUNCTION herm(M,X)
c   ........................................................................
c     Hermite polynomials
c   ........................................................................
       implicit real*8 (a-h,o-z)
       DIMENSION HE(30)
       IF(M.EQ.0) THEN
         HERM=1.
         RETURN
       ENDIF
       IF(M.EQ.1) THEN
         HERM=X
         RETURN
       ENDIF
       HE(1)=1.
       HE(2)=X
       M1=M-1
       DO K=1,M1
         HE(K+2)=X*HE(K+1)-FLOAT(K)*HE(K)
       ENDDO
       HERM=HE(M+1)
       RETURN
       END
       FUNCTION densz(M,Z,IRECA)
c   ........................................................................
c     called by SCHERM
c     Look for the distribution in z position
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       DENSZ=0.
       DO K=1,M
         KAP=K-1
         IF(IRECA.EQ.0)
     X    DENSZ=DENSZ+EXP(-Z*Z/2.)*AFZT(K)*HERM(KAP,Z)
         IF(IRECA.EQ.2)
     X    DENSZ=DENSZ+EXP(-Z*Z/2.)*AFZM(K)*HERM(2*KAP,Z)
         IF(IRECA.EQ.3)
     X    DENSZ=DENSZ+EXP(-Z*Z/2.)*AFZR(K)*HERM(2*KAP,Z)
       ENDDO
       RETURN
       END
       FUNCTION codsy(BB,CC,DD,EE,KAP)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/gaus13/H(13),T(13)
       CODSY=0.
       K=KAP-1
       FF=2.*EE-CC
       DO I=1,13
         DENDIFR=0.
         Z= (CC+FF)/2.+(CC-FF)*T(I)/2.
         IF(Z.GE.FF.AND.Z.LT.EE) THEN
           ZS=2.*EE-Z
           Z1=DD+BB-ZS
           IF(ZS.LT.BB) DENDIFR=DENSZ(NMAZ,ZS,0)
           IF(ZS.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,ZS,0)-DENSZ(NMAZ,Z1,0)
         ENDIF
         IF(Z.GE.EE) THEN
           IF(Z.LT.BB) DENDIFR=DENSZ(NMAZ,Z,0)
           Z1=DD+BB-Z
           IF(Z.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,Z,0)-DENSZ(NMAZ,Z1,0)
         ENDIF
         IF(DENDIFR.LT.0.) DENDIFR=0.
         ZZ=Z-EE
C        ZRMS2(EFFECTIF)=ZRMS2(CALCULE)*ZRMS
         ZZ=ZZ/ZRMS2
         CODSY=CODSY+H(I)*HERM(2*K,ZZ)*DENDIFR
       ENDDO
       CODSY=CODSY*(CC-FF)/2.
       RETURN
       END
       FUNCTION codif(BB,CC,DD,EE,EE1,KAP)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/SIZR/XRMS3,YRMS3,ZRMS3,ZCGR3
       COMMON/gaus13/H(13),T(13)
       CODIF=0.
       K=KAP-1
       FF=2.*EE1-EE
       FF1=2.*EE-CC
       DENDIFR=0.
       DO I=1,13
         Z=(EE+FF)/2.+(EE-FF)*T(I)/2.
         IF(Z.LE.FF1) DENDIFR=DENSZ(NMAZ,Z,0)
         IF(Z .GE.FF1)  THEN
           ZS=2.*EE-Z
           Z1=DD+BB-ZS
           IF(ZS.LT.BB) DENDIFR=DENSZ(NMAZ,ZS,0)
           IF(ZS.GE.BB)
     X       DENDIFR=DENSZ(NMAZ,ZS,0)-DENSZ(NMAZ,Z1,0)
           DENDIFR=DENSZ(NMAZ,Z,0)-DENDIFR
         endif
         IF(Z.Gt.0.) dendifr=0.
         IF(DENDIFR.LT.0.) DENDIFR=0.
         ZZ=Z-EE1
C        ZRMS2(EFFECTIF)=ZRMS2(CALCULE)*ZRMS
         ZZ=ZZ/ZRMS3
         CODIF=CODIF+H(I)*HERM(2*K,ZZ)*DENDIFR
       ENDDO
       CODIF=-CODIF*FF/2.
       RETURN
       END
       FUNCTION varia(BB,CC,DD,EE)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/gaus13/H(13),T(13)
       VARIA=0.
       CODI1=0.
       CODI2=0.
       FF=2.*EE-CC
       DO I=1,13
         DENDIFR=0.
         Z= (CC+FF)/2.+(CC-FF)*T(I)/2.
         IF(Z.GE.FF.AND.Z.LT.EE) THEN
           ZS=2.*EE-Z
           Z1=DD+BB-ZS
           IF(ZS.LT.BB) DENDIFR=DENSZ(NMAZ,ZS,0)
           IF(ZS.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,ZS,0)-DENSZ(NMAZ,Z1,0)
         ENDIF
         IF(Z.GE.EE) THEN
           IF(Z.LT.BB) DENDIFR=DENSZ(NMAZ,Z,0)
           Z1=DD+BB-Z
           IF(Z.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,Z,0)-DENSZ(NMAZ,Z1,0)
         ENDIF
         ZZ=Z-EE
         IF(DENDIFR.LT.0.) DENDIFR=0.
         CODI1=CODI1+H(I)*ZZ*ZZ*DENDIFR
         CODI2=CODI2+H(I)*DENDIFR
       ENDDO
       VARIA=CODI1/CODI2
       RETURN
       END
       FUNCTION variz(BB,CC,DD,EE,EE1)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/gaus13/H(13),T(13)
       VARIZ=0.
       CODI1=0.
       CODI2=0.
       FF=2.*EE1-EE
       FF1=2.*EE-CC
       DO I=1,13
         DENDIFR=0.
         Z= (EE+FF)/2.+(EE-FF)*T(I)/2.
         IF(Z.LE.FF1) DENDIFR=DENSZ(NMAZ,Z,0)
         IF(Z .GE.FF1)  THEN
           ZS=2.*EE-Z
           Z1=DD+BB-ZS
           IF(ZS.LT.BB) DENDIFR=DENSZ(NMAZ,ZS,0)
           IF(ZS.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,ZS,0)-DENSZ(NMAZ,Z1,0)
           DENDIFR=DENSZ(NMAZ,Z,0)-DENDIFR
         endif
         IF(Z.Gt.EE) dendifr=0.
         IF(DENDIFR.LT.0.) DENDIFR=0.
         ZZ=Z-EE1
         CODI1=CODI1+H(I)*ZZ*ZZ*DENDIFR
         CODI2=CODI2+H(I)*DENDIFR
       ENDDO
       VARIZ=CODI1/CODI2
       RETURN
       END
       FUNCTION grz(AA,BB,CC,DD,EE)
c   ........................................................................
c    specific function called by SCHERM
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/gaus13/H(13),T(13)
       GRZ=0.
       GZ=0.
       GS=0.
       DO I=1,13
         Z= (EE+AA)/2.+(EE-AA)*T(I)/2.
         dif=dendir(z,aa,bb,cc,dd,ee)
         if(dif.lt.0.)dif=0.
         GZ=GZ+H(I)*Z*DIF
         GS=GS+H(I)*DIF
       ENDDO
       IF(GS.LE.0.) GRZ=0.
       IF(GS.GT.0.) GRZ=GZ/GS
       RETURN
       END
       FUNCTION dendir(Z,AA,BB,CC,DD,EE)
c   ........................................................................
c    Specific function called by SCHERM
c     Calculate the value (nt(z)-nm(z))
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       DENDIR=0.
       FF=2.*EE-CC
       IF(Z.LT.AA)DENDIR=0.
       IF(Z.GE.AA.AND.Z.LT.FF) DENDIR=DENSZ(NMAZ,Z,0)
       IF(Z.GE.FF.AND.Z.LT.EE) THEN
         ZS=2.*EE-Z
         Z1=DD+BB-ZS
         DENDIFR=0.
         IF(ZS.LT.BB) DENDIFR=DENSZ(NMAZ,ZS,0)
         IF(ZS.GE.BB)
     X     DENDIFR=DENSZ(NMAZ,ZS,0)-DENSZ(NMAZ,Z1,0)
         DENDIR=DENSZ(NMAZ,Z,0)-DENDIFR
       ENDIF
       IF(Z.GE.EE) DENDIR=0.
       RETURN
       END
       FUNCTION dendif(Z,AA,BB,CC,DD)
c   ........................................................................
c    Specific function called by SCHERM
c     Calculate the value (nt(z)-nm(z))
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       DENDIF=0.
       IF(Z.LT.AA)DENDIF=0.
       IF(Z.GE.AA.AND.Z.LT.BB) DENDIF=DENSZ(NMAZ,Z,0)
       IF(Z.GE.BB.AND.Z.LT.CC) THEN
         Z1=DD+BB-Z
         DENDIF=DENSZ(NMAZ,Z,0)-DENSZ(NMAZ,Z1,0)
       ENDIF
       IF(Z.GE.CC) DENDIF=0.
       RETURN
       END
       FUNCTION denpd(XYZ,NMAXY,NMAZ)
c   ........................................................................
c    Specific function called by SCHERM
c     Represents the distribution:(n(x)+n(y)+n(z))/3
c   ........................................................................
       implicit real*8 (a-h,o-z)
       DENPD=(DENSX(NMAXY,XYZ,1)+DENSY(NMAXY,XYZ,1)+
     X        DENSZ(NMAZ,XYZ,2))/3.
       RETURN
       END
       FUNCTION drxyz(M,XYZ,IRECA)
c   ........................................................................
c    Specific function called by SCHERM
c     Calculate the derivatives of:(n(x)+n(y)+n(z))/3
C      IRECA=0 : for the first ellipse
C      IRECA=1 : for the second ellipse
C      IRECA=2 : for the third ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/HERMRR/AFXRR(20),AFYRR(20),AFZRR(20)
       COMMON/SIZR/XRMS3,YRMS3,ZRMS3,ZCGR3
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       DERIV1=0.
       DERIV2=0.
       IF(ABS(XYZ).GT.13.) XYZ=13.
       fe=EXP(-XYZ*XYZ/2.)
       axyz=ABS(xyz)
       MSTO=M
       DO K=2,M
         KAP=K-1
         KAH=2*KAP-1
         IF(IRECA.EQ.0) THEN
           DXXYZ=fe*FLOAT(KAP)*AFXM(K)*HERS(KAH,axyz)
           DYXYZ=fe*FLOAT(KAP)*AFYM(K)*HERS(KAH,axyz)
           DZXYZ=fe*FLOAT(KAP)*AFZM(K)*HERS(KAH,axyz)
           DERIV1=DERIV1+(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
         IF(IRECA.EQ.1) THEN
           DXXYZ=fe*FLOAT(KAP)*AFXR(K)*HERS(KAH,axyz)
           DYXYZ=fe*FLOAT(KAP)*AFYR(K)*HERS(KAH,axyz)
           DZXYZ=fe*FLOAT(KAP)*AFZR(K)*HERS(KAH,axyz)
           DERIV1=DERIV1+(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
         IF(IRECA.EQ.2) THEN
           DXXYZ=fe*FLOAT(KAP)*AFXRR(K)*HERS(KAH,axyz)
           DYXYZ=fe*FLOAT(KAP)*AFYRR(K)*HERS(KAH,axyz)
           DZXYZ=fe*FLOAT(KAP)*AFZRR(K)*HERS(KAH,axyz)
           DERIV1=DERIV1+(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
       ENDDO
       DO K=1,M
         KAP=2*(K-1)
         IF(IRECA.EQ.0) THEN
           DXXYZ=-fe*AFXM(K)*HERM(KAP,XYZ)
           DYXYZ=-fe*AFYM(K)*HERM(KAP,XYZ)
           DZXYZ=-fe*AFZM(K)*HERM(KAP,XYZ)
           DERIV2=DERIV2+.5*(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
         IF(IRECA.EQ.1) THEN
           DXXYZ=-fe*AFXR(K)*HERM(KAP,XYZ)
           DYXYZ=-fe*AFYR(K)*HERM(KAP,XYZ)
           DZXYZ=-fe*AFZR(K)*HERM(KAP,XYZ)
           DERIV2=DERIV2+.5*(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
         IF(IRECA.EQ.2) THEN
           DXXYZ=-fe*AFXRR(K)*HERM(KAP,XYZ)
           DYXYZ=-fe*AFYRR(K)*HERM(KAP,XYZ)
           DZXYZ=-fe*AFZRR(K)*HERM(KAP,XYZ)
           DERIV2=DERIV2+.5*(DXXYZ+DYXYZ+DZXYZ)
         ENDIF
       ENDDO
       DRXYZ=-(DERIV1+DERIV2)/(3.*PI)
       M=MSTO
       RETURN
       END
       FUNCTION hers(M,X)
c   ........................................................................
c    Specific Hermite polynomials called by the function DRXYZ
c   ........................................................................
       implicit real*8 (a-h,o-z)
       DIMENSION HE(30)
       IF(M.EQ.1) THEN
         HERS=1.
         RETURN
       ENDIF
       IF(M.EQ.3) THEN
         HERS=X*X-3.
         RETURN
       ENDIF
       HE(1)=1.
       HE(2)=X*X-3.
       XM1=FLOAT((M+1)/2)+.01
       M1=INT(XM1)-2
       DO K=1,M1
         HE(K+2)=HERM(2*(K+2)-2,ABS(X))-FLOAT(2*(K+2)-2)*HE(K+1)
       ENDDO
       HERS=HE(M1+2)
       RETURN
       END
       FUNCTION copdr(XI,XF,KAP)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the Hermite coefficients for :(n(x)+n(y)+n(z))/3.
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/gaus13/H(13),T(13)
       COPDR=0.
       K=KAP-1
       DO I=1,13
         Z= (XI+XF)/2.+(XF-XI)*T(I)/2.
         DEND=DENRS(Z)
         COPDR=COPDR+H(I)*HERM(2*K,Z)*DEND
       ENDDO
       COPDR=COPDR*(XF-XI)/2.
       RETURN
       END
       FUNCTION densx(M,X,IRECA)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the distribution :(n(x)+n(y)+n(z))/3. at the x position
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       DENSX=0.
       DO K=1,M
         KAP=K-1
         IF(IRECA.EQ.0)
     X   DENSX=DENSX+EXP(-X*X/2.)*AFXT(K)*HERM(2*KAP,ABS(X))
         IF(IRECA.EQ.1)
     X   DENSX=DENSX+EXP(-X*X/2.)*AFXM(K)*HERM(2*KAP,ABS(X))
       ENDDO
       RETURN
       END
       FUNCTION densy(M,Y,IRECA)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the distribution :(n(x)+n(y)+n(z))/3. at the y position
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       DENSY=0.
       DO K=1,M
         KAP=K-1
         IF(IRECA.EQ.0)
     X   DENSY=DENSY+EXP(-Y*Y/2.)*AFYT(K)*HERM(2*KAP,ABS(Y))
         IF(IRECA.EQ.1)
     X   DENSY=DENSY+EXP(-Y*Y/2.)*AFYM(K)*HERM(2*KAP,ABS(Y))
       ENDDO
       RETURN
       END
       FUNCTION denrs(XYZ)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the distribution :(n(x)+n(y)+n(z))/3. for the third ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       SP1=(DENSX(NMAXY,XYZ,0)+DENSY(NMAXY,XYZ,0))/3.
       SP2=(DENSX(NMAXY,XYZ,1)+DENSY(NMAXY,XYZ,1))/3.
       DENRS=SP1-SP2+DENSZ(NMAZR,XYZ,3)/3.
       RETURN
       END
       FUNCTION scgx(XI,XF)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the c.o.g in the x-direction for the third ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/gaus13/H(13),T(13)
       CGX=0.
       CGXX=0.
       DO I=1,13
         X= (XI+XF)/2.+(XF-XI)*T(I)/2.
         DEND=(DENSX(NMAXY,X,0)-DENSX(NMAXY,X,1))
         CGX=CGX+H(I)*DEND
         CGXX=CGXX+H(I)*DEND*X
       ENDDO
       CGX=CGX*(XF-XI)/2.
       CGXX=CGXX*(XF-XI)/2.
       SCGX=CGXX/CGX
       RETURN
       END
       FUNCTION scgy(XI,XF)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the c.o.g in the y-direction for the third ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/gaus13/H(13),T(13)
       CGY=0.
       CGYY=0.
       DO I=1,13
         Y=(XI+XF)/2.+(XF-XI)*T(I)/2.
         DEND=(DENSY(NMAXY,Y,0)-DENSY(NMAXY,Y,1))
         CGY=CGY+H(I)*DEND
         CGYY=CGYY+H(I)*DEND*Y
       ENDDO
       CGY=CGY*(XF-XI)/2.
       CGYY=CGYY*(XF-XI)/2.
       SCGY=CGYY/CGY
       RETURN
       END
       FUNCTION corxy(XI,XF,KAP,IK,XYRMS)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the Hermite coefficients for the third ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/gaus13/H(13),T(13)
        CORXY=0.
        K=KAP-1
        DO I=1,13
          Z= (XI+XF)/2.+(XF-XI)*T(I)/2.
          IF(IK.EQ.0)
     X    DEND=DENSX(NMAXY,Z,0)-DENSX(NMAXY,Z,1)
          IF(IK.EQ.1)
     X    DEND=DENSY(NMAXY,Z,0)-DENSY(NMAXY,Z,1)
          IF(DEND.LT.0.) DEND=0.
          Z=Z/XYRMS
          CORXY=CORXY+H(I)*HERM(2*K,Z)*DEND
       ENDDO
       CORXY=CORXY*(XF-XI)/2.
       RETURN
       END
       FUNCTION varxy(XI,XF,IK)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the Hermite coefficients for the residual part of the bunch
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/gaus13/H(13),T(13)
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
        VARXY=0.
        CORXY1=0.
        CORXY2=0.
        DO I=1,13
         Z= (XI+XF)/2.+(XF-XI)*T(I)/2.
         IF(IK.EQ.0)
     X   DEND=DENSX(NMAXY,Z,0)-DENSX(NMAXY,Z,1)
         IF(IK.EQ.1)
     X   DEND=DENSY(NMAXY,Z,0)-DENSY(NMAXY,Z,1)
         CORXY1=CORXY1+H(I)*Z*Z*DEND
         CORXY2=CORXY2+H(I)*DEND
       ENDDO
       VARXY=CORXY1/CORXY2
       RETURN
       END
       FUNCTION varzr(EE,CC,NMAZR)
       implicit real*8 (a-h,o-z)
c   ........................................................................
c    Specific function called by SCHERM
c    Calculate the rms sizes of th residual part of the bunch
c   ........................................................................
       VARZR=0.
       SMT=DENSZ(NMAZR,0.0d0,3)/2.
       if(smt.gt.0.) then
         SMTE=SMT/1000.
         V1=0.
         V2=CC-EE
20       Z1=(V1+V2)/2.
         SPL=DENSZ(NMAZR,Z1,3)
         IF(ABS(SPL-SMT).LE.SMTE) GO TO 10
         IF(SPL.GT.SMT) THEN
           V1=Z1
           GO TO 20
         ENDIF
         IF(SPL.LT.SMT) THEN
           V2=Z1
           GO TO 20
         ENDIF
10       CONTINUE
         VARZR=2.*Z1/2.36
       else
         VARZR=0.
       endif
       RETURN
       END
       SUBROUTINE cdg(IDCH)
c   ........................................................................
c    Calculate the c.o.g. of the bunch
C      IDCH = 1:WITH  CHASE
C      IDCH NE 1:OTHERWISE
c       cog(1) : Energy(MeV)
c       cog(3) : t.o.f. (sec)
c       cog(4) : x-direction (cm)
c       cog(5) : xp(mrd)
c       cog(6) : y-direction (cm)
c       cog(7) : yp(mrd)
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/DYN/TREF,VREF
       LOGICAL chasit
       cog(1)=0.
       cog(3)=0.
       cog(4)=0.
       cog(5)=0.
       cog(6)=0.
       cog(7)=0.
       imaxf=0
       IF(IDCH.EQ.1.) THEN
         DO i=1,ngood
           IF(ICHAS(I).EQ.1) THEN
             cog(1)=cog(1)+f(7,i)
             cog(3)=cog(3)+f(6,i)
             cog(4)=cog(4)+f(2,i)
             cog(5)=cog(5)+f(3,i)
             cog(6)=cog(6)+f(4,i)
             cog(7)=cog(7)+f(5,i)
             imaxf=imaxf+1
           ENDIF
         ENDDO
       ELSE
         DO i=1,ngood
           cog(1)=cog(1)+f(7,i)
           cog(3)=cog(3)+f(6,i)
           cog(4)=cog(4)+f(2,i)
           cog(5)=cog(5)+f(3,i)
           cog(6)=cog(6)+f(4,i)
           cog(7)=cog(7)+f(5,i)
           imaxf=imaxf+1
         ENDDO
       ENDIF
       cog(1)=cog(1)/imaxf
       cog(3)=cog(3)/imaxf
       cog(4)=cog(4)/imaxf
       cog(5)=cog(5)/imaxf
       cog(6)=cog(6)/imaxf
       cog(7)=cog(7)/imaxf
       return
       end
       SUBROUTINE ext2d(IDCH)
c   ........................................................................
c EXT2D looks for average extensions squared and returns them in array exten
c used in the routines: EMIPRT ETGAP RESTAY STATIS
c
C      IDCH = 1:WITH  CHASE TEST
C      IDCH NE 1:OTHERWISE
c
c       cog(1) : Energy(MeV)
c       cog(3) : t.o.f. (sec)
c       cog(4) : x-direction (cm)
c       cog(5) : xp(mrd)
c       cog(6) : y-direction (cm)
c       cog(7) : yp(mrd)
c
c  exten(1) : Sum( dE*dE )  MeV*MeV
c  exten(2) : Sum( dE*dPHase ) MeV*rad
c  exten(3) : Sum( dPHase*dPHase ) rad*rad
c  exten(4) : Sum( x*x )   cm*cm
c  exten(5) : Sum( xp*xp )   mrad*mrad
c  exten(6) : Sum( y*y )   cm*cm
c  exten(7) : Sum( yp*yp )  mrad*mrad
c  exten(8) : Sum( x*xp )   cm*mrad
c  exten(9) : Sum( y*yp )   cm*mrad
c  exten(12): Sum( (dp/p)*x(i) )  (cm), with p of REF
c  exten(13): Sum( (dp/p)*y(i) )  (cm), with p of REF
c  exten(14): Sum( (dp/p)*(dp/p) )  compared to p of REF
c  exten(15): Sum( (dp/p)*x(i) )  (cm), with p of COG
c  exten(16): Sum( (dp/p)*y(i) )  (cm), with p of COG
c  exten(17): Sum( (dp/p)*(dp/p) )  compared to p of COG
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/DYN/TREF,VREF
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       dimension ffd(iptsz)
       LOGICAL chasit
       exten(1)=0.
       exten(2)=0.
       exten(3)=0.
       exten(4)=0.
       exten(5)=0.
       exten(6)=0.
       exten(7)=0.
       exten(8)=0.
       exten(9)=0.
       exten(12)=0.
       exten(13)=0.
       exten(14)=0.
       exten(15)=0.
       exten(16)=0.
       exten(17)=0.
       qmoy=0.
c --- imaxf: particles keep in the bunch. With CHASE imaxf may be different from ngood(i.e. imaxf < ngood)
       imaxf=0
       gcog=cog(1)/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       bref=vref/vl
       gref=1./sqrt(1.-bref*bref)
       do i=1,ngood
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         fd(i)=(gpai*bpai)/(gref*bref)
         ffd(i)=(gpai*bpai)/(gcog*bcog)
       enddo
       do i=1,ngood
         FDPP=fd(i)-1.
         FFDPP=ffd(i)-1.
         fdp=f(7,i)-cog(1)
         TRPH1=FH*(f(6,i)-cog(3))
         TRXF=F(2,I)-COG(4)
         TRTF=F(3,I)-cog(5)
         TRYF=F(4,I)-COG(6)
         TRPF=F(5,I)-cog(7)
c --- with CHASE
         IF(IDCH.EQ.1.AND.ICHAS(I).EQ.1) THEN
           exten(1)=exten(1)+FDP**2
           exten(2)=exten(2)+TRPH1*FDP
           exten(3)=exten(3)+TRPH1**2
           EXTEN(4)=EXTEN(4)+TRXF**2
           exten(5)=exten(5)+TRTF**2
           EXTEN(6)=EXTEN(6)+TRYF**2
           exten(7)=exten(7)+TRPF**2
           EXTEN(8)=EXTEN(8)+TRXF*TRTF
           EXTEN(9)=EXTEN(9)+TRYF*TRPF
           exten(12)=exten(12)+FDPP*TRXF
           exten(13)=exten(13)+FDPP*TRYF
           exten(14)=exten(14)+FDPP**2
           exten(15)=exten(15)+FFDPP*TRXF
           exten(16)=exten(16)+FFDPP*TRYF
           exten(17)=exten(17)+FFDPP**2
           qmoy=qmoy+f(9,i)
           imaxf=imaxf+1
         ENDIF
c --- without CHASE
         IF(IDCH.NE.1) THEN
           exten(1)=exten(1)+FDP**2
           exten(2)=exten(2)+TRPH1*FDP
           exten(3)=exten(3)+TRPH1**2
           EXTEN(4)=EXTEN(4)+TRXF**2
           exten(5)=exten(5)+TRTF**2
           EXTEN(6)=EXTEN(6)+TRYF**2
           exten(7)=exten(7)+TRPF**2
           EXTEN(8)=EXTEN(8)+TRXF*TRTF
           EXTEN(9)=EXTEN(9)+TRYF*TRPF
           exten(12)=exten(12)+FDPP*TRXF
           exten(13)=exten(13)+FDPP*TRYF
           exten(14)=exten(14)+FDPP**2
           exten(15)=exten(15)+FFDPP*TRXF
           exten(16)=exten(16)+FFDPP*TRYF
           exten(17)=exten(17)+FFDPP**2
           qmoy=qmoy+f(9,i)
           imaxf=imaxf+1
         ENDIF
       enddo
       do i=1,9
         exten(i)=exten(i)/float(imaxf)
       enddo
       do i=12,17
         exten(i)=exten(i)/float(imaxf)
       enddo
       qmoy=qmoy/float(imaxf)
       return
       end
       SUBROUTINE ext2(IDCH)
c   ........................................................................
c EXT2 looks for average extensions squared and returns them in array exten
c used in the routines: stapl  tiltbm  ytzp
c
C      IDCH = 1:WITH  CHASE TEST
C      IDCH NE 1:OTHERWISE
c
c       cog(1) : Energy(MeV)
c       cog(3) : t.o.f. (sec)
c       cog(4) : x-direction (cm)
c       cog(5) : xp(mrd)
c       cog(6) : y-direction (cm)
c       cog(7) : yp(mrd)
c
c  exten(1) : Sum( (dp/p)*(dp/p) ) compared to p of COG
c  exten(2) : Sum( (dp/p)*phase )  (rad)
c  exten(3) : Sum( phase*phase) ) (rad*rad)
c  exten(4) : Sum( x(i)*x(i) )   cm*cm
c  exten(5) : Sum( xp(i)*xp(i) )   mrad*mrad
c  exten(6) : Sum( y(i)*y(i) )   cm*cm
c  exten(7) : Sum( yp(i)*yp(i) )  mrad*mrad
c  exten(8) : Sum( x(i)*xp(i) )   cm*mrad
c  exten(9) : Sum( y(i)*yp(i) )   cm*mrad
c  exten(10): Sum( dE*dE )        (MeV*MeV)
c  exten(11): Sum( dE*phase )     (MeV*rad)
c  exten(12): Sum( (dp/p)*x(i) )  (cm), with p of REF
c  exten(13): Sum( (dp/p)*y(i) )  (cm), with p of REF
c  exten(14): Sum( (dp/p)*(dp/p) ) compared to p of REF
c  exten(15): Sum( (dp/p)*x(i) )  (cm), with p of COG
c  exten(16): Sum( (dp/p)*y(i) )  (cm), with p of COG
c  exten(17): Sum( (dp/p)*(dp/p) )  compared to p of COG
c   ........................................................................

c    Calculate the emittances of the bunch
C      IDCH = 1:WITH  CHASE TEST
C      IDCH NE 1:OTHERWISE
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/DYN/TREF,VREF
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       dimension fdd(iptsz)
       LOGICAL chasit
c EXT2 looks for average extensions squared and returns them in exten(14)
       exten(1)=0.
       exten(2)=0.
       exten(3)=0.
       exten(4)=0.
       exten(5)=0.
       exten(6)=0.
       exten(7)=0.
       exten(8)=0.
       exten(9)=0.
       exten(10)=0.
       exten(11)=0.
       exten(12)=0.
       exten(13)=0.
       exten(14)=0.
       exten(15)=0.
       exten(16)=0.
       exten(17)=0.
       qmoy=0.
       imaxf=0
       gcog=cog(1)/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       bref=vref/vl
       gref=1./sqrt(1.-bref*bref)
       do i=1,ngood
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         fd(i)=(gpai*bpai)/(gcog*bcog)
         fdd(i)=(gpai*bpai)/(gref*bref)
       enddo
       do i=1,ngood
         fdp=fd(i)-1.
         FDDP=fdd(i)-1.
         dener=f(7,i)-cog(1)
         TRPH1=FH*(f(6,i)-cog(3))
         TRXF=F(2,I)-COG(4)
         TRTF=F(3,I)-cog(5)
         TRYF=F(4,I)-COG(6)
         TRPF=F(5,I)-cog(7)
         IF(IDCH.EQ.1.AND.ICHAS(I).EQ.1) THEN
           exten(1)=exten(1)+FDP**2
           exten(2)=exten(2)+TRPH1*FDP
           exten(3)=exten(3)+TRPH1**2
           EXTEN(4)=EXTEN(4)+TRXF**2
           exten(5)=exten(5)+TRTF**2
           EXTEN(6)=EXTEN(6)+TRYF**2
           exten(7)=exten(7)+TRPF**2
           EXTEN(8)=EXTEN(8)+TRXF*TRTF
           EXTEN(9)=EXTEN(9)+TRYF*TRPF
           exten(10)=exten(10)+dener*dener
           exten(11)=exten(11)+dener*trph1
           exten(12)=exten(12)+FDDP*TRXF
           exten(13)=exten(13)+FDDP*TRYF
           exten(14)=exten(14)+FDDP**2
           exten(15)=exten(15)+FDP*TRXF
           exten(16)=exten(16)+FDP*TRYF
           exten(17)=exten(17)+FDP**2
           qmoy=qmoy+f(9,i)
           imaxf=imaxf+1
         ENDIF
         IF(IDCH.NE.1) THEN
           exten(1)=exten(1)+FDP**2
           exten(2)=exten(2)+TRPH1*FDP
           exten(3)=exten(3)+TRPH1**2
           EXTEN(4)=EXTEN(4)+TRXF**2
           exten(5)=exten(5)+TRTF**2
           EXTEN(6)=EXTEN(6)+TRYF**2
           exten(7)=exten(7)+TRPF**2
           EXTEN(8)=EXTEN(8)+TRXF*TRTF
           EXTEN(9)=EXTEN(9)+TRYF*TRPF
           exten(10)=exten(10)+dener*dener
           exten(11)=exten(11)+dener*trph1
           exten(12)=exten(12)+FDDP*TRXF
           exten(13)=exten(13)+FDDP*TRYF
           exten(14)=exten(14)+FDDP**2
           exten(15)=exten(15)+FDP*TRXF
           exten(16)=exten(16)+FDP*TRYF
           exten(17)=exten(17)+FDP**2
           qmoy=qmoy+f(9,i)
           imaxf=imaxf+1
         ENDIF
       enddo
       do i=1,17
         exten(i)=exten(i)/float(imaxf)
       enddo
       qmoy=qmoy/float(imaxf)
       return
       end
c       SUBROUTINE mfordre(rc,ra,rb)
c .................................................................
c  Calculates RC=RA*RB
c .................................................................
c       implicit real*8 (a-h,o-z)
c       DIMENSION RA(6,6), RB(6,6), RC(6,6)
c       DO I1 = 1, 6
c         DO I2 = 1, 6
c           GHOST = 0.0
c           DO I3 = 1, 6
c             GHOST = GHOST + RA(I1,I3)*RB(I3,I2)
c           ENDDO
c           RC(I1,I2) = GHOST
c         ENDDO
c       ENDDO
c       RETURN
c       END
       SUBROUTINE chrefe
c .................................................................
c    change of reference frame
c    ENTRY :
c     XC YC A
c    XC   : DISPLACEMENT IN THE HORIZONTAL DIRECTION (CM)
c    YC   : DISPLACEMENT IN THE VERTICAL DIRECTION   (CM)
c    A    : ROTATION ABOUT THE vertical AXIS         (DEG)
c .................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/TAPES/IN,IFILE,META
       READ(IN,*)XC,YC,A
       a=a*pi/180.
       WRITE(16,100)XC,YC,A
       do i=1,ngood
         x=f(2,i)
         xp=f(3,i)*0.001
         yp=f(5,i)*0.001
         y=f(4,i)
         x0=x
         x=((x-yc)*cos(xp)+xc*sin(xp))/cos(xp-a)
         xp=xp-a
         xl=xc-x*sin(a)
         yl=yc-x0+x*cos(a)
         dl=sqrt(xl*xl+yl*yl)
         dl=sign(dl,xl)
         y=y+dl*tan(yp)
         f(2,i)=x
         f(3,i)=xp*1000.
         f(4,i)=y
         f(5,i)=yp*1000.
       enddo
  100  FORMAT(' New reference frame  XC =',F6.2,' CM , YC =',
     *         F6.2,' CM , A =',F6.4,' RADIAN',///)
       return
       end
       SUBROUTINE etac
c   ..............................................................................
c     Several charge state in the  bunch
c                   generated randomly
c      ENTRY :
c        N       : Number of charge states (maximum 6 different charge states)
c        CHARGE(I)  PCENT(I)  EOFF (I = 1 to N )
c        CHARGE(I)    : charge state
c        PCENT(I)     : percentage of charge state
c        EOFF(I)      : absolute energy offset of charge state w.r.t. COG (MeV)
c        ix           : RANDOM NUMBER GENERATOR
c   ..............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       common/tapes/in,ifile,meta
       common/mcs/imcs,ncstat,cstat(20)
       dimension charge(20),pcent(20),charm(20),pc(20),eoff(20),vecx(1)
       dimension foo(20,9),NDP(20)
       character*80 myfile
       READ(IN,*) NCSTAT
       N=NCSTAT
       IF (N.GT.20) THEN
         WRITE (16,140) N
         STOP
       ENDIF
       imcs=1
       IF(N.EQ.0) then
c read charge state distribution from file
         READ(in,3333) myfile(1:80)
3333     FORMAT(A80)
         write(16,*) 'Charge state distribution file: ',myfile(1:80)
         open(56,file=myfile,status='unknown')
         read(56,*) ntot
         write(16,*) 'Maximum number of particles:',imax
         write(16,*) 'Number  of   good particles:',ngood
         write(16,*) 'Number of particles in charge state file:',ntot
         IF(ntot.lt.ngood) then
           write(16,*) 'Not enough particles in charge state file'
           STOP
         ENDIF
         DO J=1,ntot
           read(56,*) chstate
           f(9,i)=chstate
         enddo
         close(56)
       ELSE
c generate charge state distribution
         POURC=0.
         DO I=1,N
c for each charge state read charge, percentage and energy offset
           READ(IN,*) charm(I),PC(I),eoff(i)
           CSTAT(i)=charm(i)
           POURC=POURC+PC(I)
           IF (POURC.GT.100.) THEN
             WRITE (16,100) I,POURC
             STOP
           ENDIF
         enddo
         POURC=0.
         DO I=1,N
c           WRITE(16,110) charm(I),PC(I),eoff(i)
           POURC=POURC+PC(I)
         enddo
         IF (POURC.NE.100.) THEN
           WRITE (16,120) POURC
           STOP
         ENDIF
         j=1
25       continue
         ts=500.
         is=1
         do i=1,n
           if(ts.ge.charm(i)) then
             ts=charm(i)
             is=i
           endif
         enddo
         charge(j)=ts
         pcent(j)=pc(is)
         charm(is)=1000.
         j=j+1
         if(j.le.n) goto 25
         write(16,*) '**********************'
c for each charge state write charge, percentage and energy offset
         jjj=0
         do i=1,n
           WRITE(16,110) CHARGE(I),PCENT(I),eoff(i)
           if(charge(i).eq.qst) jjj=1
         enddo
         len=1
C FIRST TRAJECTORY HAS CHARGE STATE AS DEFINED BY INPUT
         PCENT(1)=PCENT(1)/100.
         DO I=2,N
           PCENT(I)=PCENT(I-1)+PCENT(I)/100.
         enddo
         DO I=1,IMAX
           call rlux(vecx,len)
           XARPHA=VECX(1)
           IF (XARPHA.LE.PCENT(1)) THEN
             f(9,i)=CHARGE(1)
             f(7,i)=f(7,i)+eoff(1)
           ELSE
             DO J=1,N-1
               IF (XARPHA.GT.PCENT(J) .AND. XARPHA.LE.PCENT(J+1))then
                 f(7,i)=f(7,i)+eoff(J+1)
                 f(9,i)=CHARGE(J+1)
               ENDIF
             enddo
           ENDIF
         ENDDO
       ENDIF
c print energy, boro for each charge state
       DO k=1,ncstat
         NDP(k)=0
         do j=2,7
           foo(k,j)=0.
         enddo
       ENDDO
       DO i=1,imax
         DO k=1,ncstat
           if(f(9,i).eq.cstat(k)) then
             NDP(k)=NDP(k)+1
             do j=2,7
               foo(k,j)=foo(k,j)+f(j,i)
             enddo
           endif
         ENDDO
       ENDDO
       DO k=1,ncstat
         do j=2,7
           foo(k,j)=foo(k,j)/float(NDP(k))
         enddo
       ENDDO
       DO k=1,ncstat
         gref=foo(k,7)/xmat
         bref=sqrt(1.-1./(gref*gref))
         xe=(gref-1.)*xmat
c   magnetic rigidity
         bor=3.3356*xmat*bref*gref/cstat(k)
         write(16,*) ' Q: ',cstat(k),' COG : energy ',xe,
     *     ' MeV  momentum ',bor,' kG.cm'
       ENDDO
100    FORMAT(3X,' WRONG PERCENTAGE IN CHARGE STATE DISTRIBUTION',/,
     * 4X,' CHARGE STATE ',I3,' PERCENTAGE ',E12.5)
110    FORMAT(3X,' CHARGE STATE ',F6.1,'  PERCENTAGE ',E12.5,' %',
     * 4X,' ENERGY OFFSET ',E12.5,' MeV')
120    FORMAT(3X,' TOTAL PERCENTAGE OF ALL CHARGE STATES < 100 %',/,
     * 4X,' PERCENTAGE ',E12.5)
140    FORMAT(3X,' NUMBER OF CHARGE STATES : ',I3,' GREATER THAN 20')
       RETURN
       END
       SUBROUTINE crest(betr,eqvl,xpos,bkcr,ffield)
c  .....................................................................
c    called by RESTAY
c    Look for the beta giving the maximun energy gain
c     iterative method
c  .....................................................................
       implicit real*8 (a-h,o-z)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       XLHE=YLG
       ITR=0
       BE1=BETR
       BE2=BETR-8.333e-03*BETR
       XLEQ=XLHE
456    CONTINUE
       ITR=ITR+1
       XK1=FH/(VL*BE1)
       XK2=FH/(VL*BE2)
       T1=TTA0(BE1)*ffield
       S1=TSB0(BE1)*ffield
       TP1=TTA1(BE1)*ffield
       SP1=TSB1(BE1)*ffield
       T2=TTA0(BE2)*ffield
       S2=TSB0(BE2)*ffield
       TP2=TTA1(BE2)*ffield
       SP2=TSB1(BE2)*ffield
       DTS=(T1*TP1+S1*SP1)/(T1*T1+S1*S1)
       A1K12=(T1*TP1+S1*SP1)/(T2*TP2+S2*SP2)
       A2K12=(T2*T2+S2*S2)/(T1*T1+S1*S1)
       AK12=A1K12*A2K12
       AK12=1./AK12
       BK12=(XK2-XK1)/(AK12-1.)
       BK12=BK12*DTS
       DESY=-4.*ATAN(DTS*3.2/XLEQ)
C     improve  DESYNCHRONISATION
       EPSRD=1.e-04
       IF(ABS(DESY).GE.EPSRD) THEN
comment       write(6,*) ' desy ',itr,abs(desy)
comment       pause
         TIL2=DESY/2.
         DO III=1,3
           FTIL=TIL2/TAN(TIL2)-1.-BK12
           DFTIL=-TIL2/(SIN(TIL2)*SIN(TIL2))+1./TAN(TIL2)
           TIL2=TIL2-FTIL/DFTIL
         ENDDO
         DESY=TIL2*2.
       ENDIF
       IF(ABS(DESY).LT.EPSRD) THEN
         XPOS=(T1*SP1-S1*TP1)/(T1*T1+S1*S1)
         BKCR=SQRT(T1*T1+S1*S1)
         EQVL=XLEQ
         RETURN
       ENDIF
C    calculates the equivalent length
       XLEQ=DESY*(AK12-1.)/(XK2-XK1)
       DELTK=DESY/XLEQ
       XKCRT=XK1-DELTK
       BECRT=FH/(VL*XKCRT)
       BE1=BECRT
       BE2=BECRT-BECRT/120.
       GO TO 456
       END
       SUBROUTINE gaus(r1,r2,z1,z2,opt,er,ez)
c  .....................................................................
c    called by SCHEFF
c---calculate er and ez at the r and z location given in fldcom
c---by gauss quadrature integration over the double
c---interval from r1 to r2 and from z1 to z2.
c---if opt.gt.0, determine number of integration points as follows.
c---let rat = max(cr/cz,cz/cr), where cr=r2-r1, and cz=z2-z1
c---     if rat.le.2, use 2 x 2 point array
c---     if rat.gt.2, use 2 x 4 point array
c---     if rat.gt.4, use 2 x 6 point array
c---if opt.lt.0, use 2 x 2 point array regardless of rat
c  .....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       dimension r(6),z(6),wr(6),wz(6),xx(3,3),wx(3,3)
       data ((xx(i,j),i=1,3),j=1,3)/.2113248654,0.0,0.0,
     1   .06943184420,
     1  .33000947820,0.0,.03376524290,.16939530680,.3806904070/
       data ((wx(i,j),i=1,3),j=1,3)/.50,0.0,0.0,.17392742260,
     1   .32607257740,
     1  0.0, .085662246190,.1803807865000,.2339569673000/
       cr=r2-r1
       cz=z2-z1
       ir=1
       jz=1
       m=1
       if (opt.lt.0.) go to 20
c---determine number of integration points
       rat=abs(cz/cr)
       l=0
       if (rat.ge.1.)go to 10
       rat=1./rat
       l=1
   10  if (rat.gt.2.)m=2
       if (rat.gt.4.)m=3
       if (l.eq.0) jz=m
       if (l.eq.1) ir=m
   20  do i=1,ir
         k=2*i-1
         r(k)=r1+cr*xx(i,ir)
         r(k+1)=r2-cr*xx(i,ir)
         wr(k)=wx(i,ir)
         wr(k+1)=wx(i,ir)
       enddo
       do j=1,jz
         k=2*j-1
         z(k)=z1+cz*xx(j,jz)
         z(k+1)=z2-cz*xx(j,jz)
         wz(k)=wx(j,jz)
         wz(k+1)=wx(j,jz)
       enddo
       ser=0.
       sez=0.
       kr=2*ir
       kz=2*jz
       do i=1,kr
         do j=1,kz
           call flds(r(i),z(j),er1,ez1)
           ser=ser+wr(i)*wz(j)*er1*r(i)
           sez=sez+wr(i)*wz(j)*ez1*r(i)
         enddo
       enddo
       er=cr*cz*ser
       ez=cr*cz*sez
       return
       end
       SUBROUTINE flds(r,z,er,ez)
c   ..................................................................
c    called by SCHEFF
c     evaluate fields at r1,z1 due to ring of charge at r,z.
c     er=(pi/2)*r/d**3.     ez=(pi/2)*z/d**3.
c   ..................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/fldcom/ rp, zp,pl,opt,nip
       r1=rp
       z1=zp
c **********************************
       d=z1-z
       c=(r1-r)**2
       b=(r1+r)**2
       a=4.0*r*r1/(b+d**2)
       call eint (a,ee,ek)
       er1=0.0
       a=sqrt(b+d**2)
       if (r1.eq.0.) go to 10
       er1=(ek-(r**2-r1**2+d**2)*ee/(c+d**2))/(2.0*r1*a)
   10  ez1=d*ee/(a*(c+d**2))
       if (nip.eq.0) go to 50
       do i=1,nip
         xi=i
         do j=1,2
           d=z1-(z+xi*pl)
           a=4.0*r*r1/(b+d**2)
           call eint (a,ee,ek)
           a=sqrt(b+d**2)
           if (r1.eq.0.) go to 20
           er1=er1+(ek-(r**2-r1**2+d**2)*ee/(c+d**2))/(2.0*r1*a)
20         ez1=ez1+d*ee/(a*(c+d**2))
           xi=-xi
           enddo
       enddo
   50  er=er1
       ez=ez1
       return
       end
       SUBROUTINE eint(a,ee,ek)
c   .............................................................................
c      evaluate elliptic integrals  ( called by SCHEFF)
c   .............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       b=1.0-a
       c=log(b)
       ee=1.0+b*(.4630106-0.2452740*c+b*(0.1077857-0.04125321*c))
       ek=1.38629436-.5*c+b*(0.1119697-0.1213486*c+
     *    b*(.07253230-.028874721*c))
       return
       end
       SUBROUTINE tiltz(tilta)
c  ......................................................................
c      Skew the right ellipse generated by GEBEAM in the phase plane (x,z)
c  ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/DYN/TREF,VREF
       common/tapes/in,ifile,meta
       common/etcom/cog(8),exten(17),fd(iptsz)
       WRITE(16,1)TILTA
1      FORMAT(' tilt in the plane (x,z) around the c.o.g',
     X '    ANGLE :',E12.5,' DEG',/)
        TILTA=TILTA*PI/180.
       tref=0.
       vref=0.
       DO I=1,imax
         gpai=f(7,i)/xmat
         vref=vref+vl*sqrt(1.-1./(gpai*gpai))
         tref=tref+f(6,i)
       enddo
       vref=vref/float(imax)
       tref=tref/float(imax)
       vref1=0.
       tref1=0.
       DO I=1,IMAX
         gpai=f(7,i)/xmat
         vpai=sqrt(1.-1./(gpai*gpai))*vl
         trot=(f(6,i)-tref)*cos(tilta)-sin(tilta)*f(2,i)/vpai
         xrot=(f(6,i)-tref)*sin(tilta)*vpai+cos(tilta)*f(2,i)
         f(6,i)=trot
         f(2,i)=xrot
         tref1=tref1+f(6,i)
         vref1=vref1+vpai
       ENDDO
       tref=tref1/float(imax)
       vref=vref1/float(imax)
       RETURN
       END
      SUBROUTINE rfq_o3
c*****************************************************************
C---- Dynamics through a single cell of a RFQ with multipolar expansion
C---  units: MeV , m , sec
C----    p(1): V/(r0*r0)  (MV/m2)
C----    p(2): AV (MV)
C-----   p(3): cell length CL (m)
C-----   p(4): phase of RF  at the entrance of the cell (deg)
c-----   p(5): TYPE
C----      TYPE IS 0, 1., 2., OR 3., INDICATING AS FOLLOWS
C----      0, STANDARD CELL, NO ACCELERATION
C----      1, STANDARD CELL, ACCELERATION
C----      2, FRINGE-FIELD, NO ACCELERATION
C----      3, FRINGE FIELD, ACCELERATION
c****************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/RIGID/BORO
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/tapes/in,ifile,meta
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DCSPA/IESP
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SHIF/DTIPH,SHIFT
       common/femt/iemgrw,iemqesg
       common/posc/xpsc
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       common/drfq/p(9)
       character*1 cr
       COMMON/AZLIST/ICONT,IPRIN
       common/trfq/icour,ncell
       logical iesp,ichaes,shift,iemgrw,iflag
       common/itvole/itvol,imamin
       common/tofev/ttvols
       logical itvol,imamin
c allow for print out on terminal of gap# on one and the same line
       cr=char(13)
comment       WRITE(6,8254) NRTRE,NRRES,cr
c       write(6,8254) nrtre,ncell,cr
c8254   format('Transport element:',i5,
c     *        '      RFQ cell            :',i5,a1,$)
       iflag=.false.
       radian=pi/180.
       ilost=0
       twopi=2.*pi
c  convert vl in m
       vlm=vl/100.
C      STATISTIQUES FOR PLOTS
       if(iprf.eq.1) call stapl(davtot*10.)
       VORSQ=P(1)
       AV=P(2)
       CL=P(3)
       TYPE=P(5)
       write(16,*)'*** cell :',ncell+1,' length (cm): ',cl*100.
       write(16,*)'*** V/r0**2  (kV/mm**2): ',vorsq/1000.
       write(16,*)'*** AV (kV): ',av*1000.
       if(type.eq.0.)write(16,*)'*** no acceleration, standard cell '
       if(type.eq.1.)write(16,*)'*** acceleration, standard cell '
       if(type.eq.2.)write(16,*)
     *            '*** no acceleration, fringing field region '
       if(type.eq.3.)write(16,*)
     *            '*** acceleration, fringing field region '
       wavel=2.*pi*vlm/fh
       er=xmat
       CAY=PI/CL
       NS=18
       xl=cl/float(ns)
       hl=.5*xl
c----  scl: space charge length (in cm SCHEFF unit)
       scl=cl*100.
c----  c.o.g of the bunch at the entrance of the cell
        tcog=0.
        ecog=0.
cold        qcog=0.
        do i =1,ngood
         tcog=tcog+f(6,i)
         ecog=ecog+f(7,i)
cold         qcog=qcog+f(9,i)
        enddo
        tcog=tcog/float(ngood)
        ecog=ecog/float(ngood)
cold        qcog=qcog/float(ngood)
        gcog=ecog/er
        bcog=sqrt(1.-1./(gcog*gcog))
        wcog=ecog-er
       if(ncell.eq.0)then
c ---- shift = .false. ==> synchronous particle coincide with cog in the cell ncell = 0
        if(.not.shift) then
         write(16,*) '*** ref. part. and cog coincide in ncell = 0'
         tref=tcog
         bref=bcog
         vref=bref*vl
         gref=gcog
         wref=wcog
         wrefi=wref
        else
c ---- shift = .true. ==> the reference particle and the cog are different at the entrance of the RFQ
         write(16,*) '*** ref. part. and cog separated  in ncell = 0'
         bref=vref/vl
         gref=1./sqrt(1.-bref*bref)
         wref=er*(gref-1.)
         wrefi=wref
        endif
       endif
       IF (TYPE.GT.1.) then
        CAY=.5*CAY
        NS=int(36.*CL/(BREF*WAVEL))
       endif
        write(16,178)
178    format(/,' Dynamics at input',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       write(16,1788) bcog,gcog,wcog,tcog*fh*180./pi,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       WRITE(16,165) bref,gref,wref,tref*fh*180./pi,tref
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
c  start prints in file 'short.data'
       idav=idav+1
       iitem(idav)=14
       dav1(idav,1)=cl*1000.
       dav1(idav,2)=vorsq/1000.
       dav1(idav,3)=av*1000
       davtot=davtot+cl*100.
       dav1(idav,4)=davtot*10.
       dav1(idav,5)=type
c----  phini: phase of the reference at input of the cell
       phini=-tref*fh+p(4)*radian
       ph0=phini*180./pi
       write(16,3945)ph0,p(4)
3945   format('phase offset at input : ',
     *        e12.5,' deg phase shift: ',e12.5,' deg')
c --- the cell is divided in ns elements of length xl(m)
       z=0.
       do n=1,ns
        z=z+hl
c---- change of synchronous particle over the half step hl = xl/2
        tref=tref+hl/(bref*vlm)
        if(itvol) ttvols=tref
        phref=tref*fh+phini
        skz=sin(cay*z)
        ckz=cos(cay*z)
c---    cell with fringing field (TYPE > 1)
        if(type.gt.1.) then
          c3kz=cos(3.*cay*z)
          skz=.75*(skz+sin(3.*cay*z))
        endif
C------  CHANGE OF ENERGY OVER STEP XL
        if (type.ne.0..and.type.ne.2.) then
         sp=sin(phref)
         dwref=.5*qst*cay*av*skz*sp*xl
         wrefm=wref+0.5*dwref
         grefm=wrefm/er+1.
         brefm=sqrt(1.-1./(grefm*grefm))
         wref=wref+dwref
         gref=wref/er+1.
         bref=sqrt(1.-1./(gref*gref))
cc         dez=.5*qst*cay*av*skz*sp
cc         dref=.5*(dez/er) * xl*xl/(brefm**3*grefm**3*vlm)
        endif
c   start computations of particles
       do ip=1,ngood
c     convert in m and rad
         xi=f(2,ip)*1.e-02
         xpi=f(3,ip)*1.e-03
         yi=f(4,ip)*1.e-02
         ypi=f(5,ip)*1.e-03
         ww=f(7,ip)-er
         gi=ww/er+1.
         bi=sqrt(1.-1./(gi*gi))
         bgi=bi*gi
         tim=f(6,ip)+hl/(bi*vlm)
         phi=phini+fh*tim
         qq=abs(f(9,ip))
         sp=sin(phi)
         cp=cos(phi)
         bav=bi
         gav=gi
         bgav=bgi
         bg=bgi
         beta=bi
         delt=0.
         amort=1.
         xm=xi+xpi*hl
         ym=yi+ypi*hl
         rm=sqrt(xm*xm+ym*ym)
c ---- the particle is lost if rm>rlim (rlim is the limit of REJECT)
c        rlim is converted in m
c         rlimm=rlim*1.e-02
          if(rm.gt.rlimm) then
           f(8,ip)=0.
           iflag=.true.
           ilost=ilost+1
           write(16,5556) ip,rm,rlimm
5556       format(' particle lost: ',i5,' radius (m): ',e12.5,
     *           ' barrier (m):',e12.5)
           go to 6525
          endif
         theta=0.
         xml=xm
         yml=ym
         if(abs(xm).gt.1.e-10) theta=atan(ym/xm)
         if(abs(xm).le.1.e-10) then
          if(abs(ym).gt.1.e-10) then
           if(xm.ge.0..and.ym.gt.0.) theta=pi/2
           if(xm.ge.0..and.ym.lt.0.) theta=-pi/2
           if(xm.lt.0..and.ym.lt.0.) theta=pi/2
           if(xm.lt.0..and.ym.gt.0.) theta=-pi/2
          endif
          if(abs(ym).le.1.e-10) theta=0.
         endif
         zrm=cay*rm
c ----- Bessel functions I0 and I1
comment          bi0=1.+zrm*zrm/4.+zrm**4/64.+zrm**6/2304.+zrm**8/1.47456e05
          bi0=1.+zrm*zrm/4.+zrm**4/64.+zrm**6/2304.
          bi1=zrm/2.+zrm**3/16.
c         transverse fields ex and ey
          erf=vorsq*cos(2.*theta)*2.*rm+
     *        cay*(av*bi1)*ckz
          erf=-erf/2.
          etf=vorsq*sin(theta)*2.*rm
          etf=etf/2.
          ex=erf*cos(theta)-etf*sin(theta)
          ey=erf*sin(theta)+etf*cos(theta)
C----  change the energy over step xl
         IF (TYPE.NE.0..AND.TYPE.NE.2.) then
c ----  energy at the middle of the element
c------ standard cell
          if(type.lt.2.) then
           ez=0.5*(av*bi0)*skz*cay
           dw=qq*ez*sp*xl
c------ fringe-field cell
          else
           dw=.5*qq*cay*avb*skz*sp*xl
          endif
          WAV=WW+.5*DW
          GA=WAV/ER
          BGAV=SQRT(GA*(2.+GA))
          GAV=1.+GA
          BAV=BGAV/GAV
c ---- energy over the step xl
          WW=WW+DW
          GA=WW/ER
          GAM=1.+GA
          BG=SQRT(GA*(2.+GA))
          beta=sqrt(1.-1/(gam*gam))
c       calculate the jump of phase (in sec)
          if(type.lt.2.) dez=qq*ez*sp
          if(type.gt.2.) dez=.5*qq*cay*avb*skz*sp
          delt=.5*(dez/er) * xl*xl/(bav**3*gav**3*vlm)
          amort=bgi/bg
         endif
         BGFAC=GAV*BAV**2
c ----- IF TYPE LT 2: standard cell
         if(type.lt.2.) then
          cc=qq*xl*sp/(bgfac*er)
          signx=1.
          signy=1.
          if(theta.gt.0.) then
           if(xm.lt.0.) signx=-1.
           if(ym.lt.0.) signy=-1.
          endif
          if(theta.lt.0.) then
           if(xm.lt.0.) signx=-1.
           if(ym.gt.0.) signy=-1.
          endif
          if(theta.eq.0.) then
           signx=0.
           signy=0.
          endif
          rr1=cc*ex*signx
          rr2=cc*ey*signy
          xpm=xpi*amort+rr1
          ypm=ypi*amort+rr2
         else
c ----- fringe-field cell (first order in transverse directions)
c*********************************************************
c         C1 = (1/m**2) * (m) = (1/m)
c         C2 = (1/m**2) * (m) = (1/m)
c         RF1 = (MeV/(MeV*m**2)) = (1/m**2)
c         RF2 =  (MeV/MeV) * (1/m**2) = (1/m**2)
c*******************************************************
           RF1=QQ*VORB/ER
           RF2=.25*QQ*CAY*CAY*AVB/ER
           C1=RF1*SP*XL/BGFAC
           C2=RF2*CKZ*SP*XL/BGFAC
           C1=C1*.75*(CKZ+C3KZ/3.)
           C2=C2*.75*(CKZ+3.*C3KZ)
           RR1=-(C1+C2)
           RR2=(C1-C2)
           xpm=xpi*amort+rr1*xm
           ypm=ypi*amort+rr2*ym
         endif
         xf=xm+xpm*hl
         yf=ym+ypm*hl
c      convert m->cm , rad->mrad
         f(2,ip)=xf*100.
         f(4,ip)=yf*100.
         f(3,ip)=xpm*1000.
         f(5,ip)=ypm*1000.
c   tof over the length xl
         f(6,ip)=f(6,ip)+hl/(bi*vlm)+hl/(beta*vlm)+delt
         f(7,ip)=ww+er
c----- print in file 49 coordinates of particle icont (not active)
comment         if(ip.eq.icont) then
comment	       write(49,8888)n,icour,z,signx,signy,ex,ey,theta,xml,yml,rr1,rr2
comment          icour=icour+1
comment         endif
comment8888      format(2(2x,i5),10(2x,e12.5))
c---- ********************************************************************
6525     continue
        enddo
c ----- reshuffle good particles
         if(iflag) then
          call shuffle
          iflag=.false.
         endif
c    space charge at the middle of the cell  (revoir pour fringe field cell)
        if(n.eq.9) then
         if(ichaes) then
C      Charge space
          iesp=.true.
          call cesp(scl)
          iesp=.false.
         endif
        endif
        tref=tref+hl/(bref*vlm)+dref
        if(itvol) ttvols=tref
        vref=bref*vl
        z=z+hl
c     Change  dp/p over the cell
        call disp
       enddo
c---  kept in the bunch particles such that phase RF > 180 deg or phase RF < -180 deg
        do i=1,ngood
         dtvl=(f(6,i)-tref)*fh
         if(dtvl.gt.pi) f(6,i)=f(6,i)-2.*pi/fh
         if(dtvl.lt.-pi) f(6,i)=f(6,i)+2.*pi/fh
        enddo
c----  c.o.g of the bunch at the output of the cell
        tcog=0.
        ecog=0.
        do i =1,ngood
         tcog=tcog+f(6,i)
         ecog=ecog+f(7,i)
        enddo
        tcog=tcog/float(ngood)
        ecog=ecog/float(ngood)
        gcog=ecog/er
        bcog=sqrt(1.-1./(gcog*gcog))
        wcog=ecog-er
c---  window control relative to the energy of the c.o.g of the bunch
c ---- ifw = 0 ===> wdisp = dW/W
c ---- ifw = 1 ===> wdisp = dW (MeV)
c ----- convert wdisp in dp/p
       if(ifw.eq.0) dispr=gcog*gcog*wdisp/(gcog*(gcog+1.))
       if(ifw.eq.1) dispr=gcog*gcog*wdisp/(gcog*(gcog+1.)*wcog)
       iflag=.false.
       do i=1,ngood
        dese=abs(fd(i)-1.)
        if(dese.gt.dispr) then
         ilost=ilost+1
         f(8,i)=0.
         write(16,*) ' particle lost: ',i,' dp/p: ',dese,
     *             ' in window :',dispr
         iflag=.true.
        endif
       enddo
       if(iflag) then
        call shuffle
c----  c.o.g of the bunch after shuffle
        tcog=0.
        ecog=0.
        do i =1,ngood
         tcog=tcog+f(6,i)
         ecog=ecog+f(7,i)
        enddo
        tcog=tcog/float(ngood)
        ecog=ecog/float(ngood)
        gcog=ecog/er
        bcog=sqrt(1.-1./(gcog*gcog))
        wcog=ecog-er
       endif
       write(16,179)
179    format(/,' Dynamics at the output',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       write(16,1788) bcog,gcog,wcog,tcog*fh*180./pi,tcog
       write(16,165) bref,gref,wref,tref*fh*180./pi,tref
       write(16,*) '  time of flight: ',ttvols*fh*180./pi,' deg'

cet2010s
c dphete,dav1(idav,16),dav1(idav,21) and dav1(idav,12) still to be assigned correct value
       dphete=0.
       emns=0.
       tofprt=fh*tcog*180./pi
       n2kp=int(tofprt/360.)
       tofprt=tofprt-float(n2kp)*360.
       if(tofprt.gt.180.) tofprt=tofprt-360.
c cavity number, transmission (%), synchronous phase (deg), time of flight (deg) (within –180 deg and 180 deg),
c COG relativistic beta (@ output), COG output energy (MeV), REF relativistic beta (@ output), REF output energy (MeV),
c horizontal emittance (mm.mrad, RMS normalized), vertical emittance (mm.mrad, RMS normalized),
c longitudinal emittance (RMS, ns.keV) <- still to be implemented (=emns)
       trnsms=100.*float(ngood)/float(imax)
       if(ncell.eq.1) write(50,*) '# rfq_o3.dmp'
       write(50,7023) ncell,trnsms,p(4),tofprt,bcog,wccog,
     *  bets,wref,0.25*dav1(idav,16),0.25*dav1(idav,21),
     *  0.25*emns
7023   format(1x,i4,1x,f6.2,2(1x,f8.3),2(1x,f7.5,1x,f11.4),
     *        3(1x,f11.4))
cet2010e

C----  new magnetic rigidity of the reference
        xmor=xmat*bref*gref
        boro=33.356*xmor*1.e-01/qst
        dav1(idav,6)=(gref-1.)*er-wrefi
        dav1(idav,36)=ngood
C   plots
        CALL STAPL(davtot*10.)
        call emiprt(0)
        ncell=ncell+1
      RETURN
      END
      SUBROUTINE tdens(m,ireca,iacc)
c   ........................................................................
c     called by SCHERM
c      Look for the shape of the distribution n(t)
c       ireca=0 : for the first ellipse
c       ireca=1 : for the second ellipse
c   ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       iacc=0
       T=1.e-05
       dt=5.e-02
       fc=drxyz(m,t,ireca)
       f1=fc
       DO 10 I=1,100
         t=t+dt
         fu=drxyz(m,t,ireca)
         if(fu.ge.fc) then
           fc=fu
           go to 10
         endif
         go to 11
10     continue
11     continue
       IF(ABS(FC/F1).GT.1.2) IACC=1
       RETURN
       END
       SUBROUTINE schfdyn
c   ........................................................................
c  Input datas for SCHEFF
c      SCE(2)=r extension in rms multiples
c      SCE(3)=z half extension in rms multiples
c      SCE(4)=no. of radial mesh intervals (le 20)
c      SCE(5)=no. of longitudinal mesh intervals (le 40)
c      sce(6)=no. of adjacent bunches, applicable for buncher studies
c                and should be 0 for linac dynamics
c      sce(7)=pulse length, if not beta lambda.(transport studies)
c              distance bewteen beam pulses
c            input zero to get default "beta lambda"
c             units are cm
c      sce(8)=dummy
c      sce(9)=option to integrate space charge forces over box
c                  if.eq.0. no integration.  see sub gaus for further
c                  explanation.
c   ........................................................................
       implicit real*8 (a-h,o-z)
       common/tapes/in,ifile,meta
       common/rcshef/sce(20)
       read(in,*) iread
       if(iread.eq.0) then
c  standard SCHEFF parameters
         sce(2)=4.
         sce(3)=4.
         sce(4)=20
         sce(5)=40
         sce(6)=0
         sce(9)=0.
         sce(7)=0.
       else
c    read  SCHEFF parameters
         READ(IN,*) sce(2),sce(3),sce(4),sce(5),sce(6),sce(7),sce(9)
       endif
       sce(3)=sce(3)*2.
       sce(8)=0.
       call schefini
       RETURN
       END
       SUBROUTINE intg3(npt)
c  .....................................................................
c   called by SCHERM
c Calculate the electric field components acting on each particle
c        Gauss method
c  .....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/HERMRR/AFXRR(20),AFYRR(20),AFZRR(20)
       COMMON/SIZR/XRMS3,YRMS3,ZRMS3,ZCGR3
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZZT/XRMSZ,YRMSZ,ZRMSZ
       COMMON/ELCG/XCGD,YCGD,ZCGD,XCGR,YCGR,ZCGR
       COMMON/INTGRT/ex,ey,ez
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ichaes
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/twcst/epsilon
       COMMON/ECOR/const
       COMMON/CGRMS/xsum,ysum,zsum
       COMMON/NPART/IMAXR
       DIMENSION UI(6),WI(6)
       DATA (UI(J),J=1,6)/.033765,.169395,.380690,.619310,.830605,
     *      .966234/
       DATA (WI(J),J=1,6)/.085662,.180381,.233957,.233957,.180381,
     *      .085662/
c
c Initialize some constants and variables
c freq. in MHz
       freq=fh*0.5e-06/pi
       xrmsc=xrmsz
       yrmsc=yrmsz
       zrmsc=zrmsz
       xgc=xsum
       ygc=ysum
       zgc=zsum
       qmpart=1.0e-09*beamc/(float(ngood)*freq)
comment        const=qmpart*xrmsc*yrmsc*zrmsc/(2.*epsilon)
       const=qmpart/(2.*epsilon)
       dnorm=(xrmsc*yrmsc*zrmsc)**.333333333
       dsq=dnorm*dnorm
       xsq=(xc(npt)-xgc)*(xc(npt)-xgc)
       ysq=(yc(npt)-ygc)*(yc(npt)-ygc)
       zsq=(zc(npt)-zgc)*(zc(npt)-zgc)
       ex=0.
       ey=0.
       ez=0.
c initialize integrals to 0.
c integrate all 3 components (x,y,z)
       DO J=1,6
         a1=xrmsc*xrmsc-dsq+dsq/ui(j)
         a2=yrmsc*yrmsc-dsq+dsq/ui(j)
         a3=zrmsc*zrmsc-dsq+dsq/ui(j)
         t1=xsq/a1
         t2=ysq/a2
         t3=zsq/a3
         txyz=sqrt(t1+t2+t3)
         if(abs(txyz).gt.13.)txyz=13.
         ff1=exp(-txyz*txyz/2.)*afzt(1)
         ff1=ff1/(2.*pi)
         fxn=ff1/(ui(j)*ui(j)*sqrt(a1)*a1*sqrt(a2)*sqrt(a3))
         fyn=ff1/(ui(j)*ui(j)*sqrt(a1)*a2*sqrt(a2)*sqrt(a3))
         fzn=ff1/(ui(j)*ui(j)*sqrt(a1)*a3*sqrt(a2)*sqrt(a3))
         ex=ex+wi(j)*fxn*dsq
         ey=ey+wi(j)*fyn*dsq
         ez=ez+wi(j)*fzn*dsq
       ENDDO
c Field components are in Newton/Coulomb
       ex=ex*const*(xc(npt)-xgc)
       ey=ey*const*(yc(npt)-ygc)
       ez=ez*const*(zc(npt)-zgc)
       RETURN
       END
       SUBROUTINE prbeam1
c   .....................................................................
C      PRINT OF PARTICLE COORDINATES
c   .....................................................................
       implicit real*8 (a-h,o-z)
       parameter(iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/sc3/beamc,scdist,sce10,cplm,ectt,apl,ichaes,iscsp
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       common/mcs/imcs,ncstat,cstat(20)
       common/tapes/in,ifile,meta
       common/etcha3/ichxyz(iptsz)
       common/rec/irec
       common/etcom/cog(8),exten(17),fd(iptsz)
       logical chasit
       character*80 wfile
       ENTRY compress(pib)
c  Do so by shifting particles belonging to the same bunch outside the (+/-) pib window to
c   inside the +/- pib window.
       pib=pib*pi/180
       do ite=1,3
c Find cog of bunch in time
         tcog=0.
         do i=1,ngood
           tcog=tcog+f(6,i)
         enddo
         tcog=tcog/float(ngood)
         do i=1,ngood
           drad=(f(6,i)-tcog)*fh
           if(drad.gt.pib) then
             f(6,i)=(f(6,i)-2.*pi/fh)
           endif
           if(drad.lt.-pib) then
             f(6,i)=(f(6,i)+2.*pi/fh)
           endif
         enddo
       enddo
c Find cog of bunch in time after shifting particles
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       write(16,59) tcog*fh*180./pi
59     format(' cog in time after shifting particles: ',e13.7,' deg')
       return
c    ****************************************
       ENTRY prbeam(iflg,wfile)
C      PRINT OF PARTICLE COORDINATES
       ecog=0.
       tcog=0.
       xav=0.
       xpav=0.
       yav=0.
       ypav=0.
c **********************************************
c   cog of the bunch
       do i=1,ngood
         ecog=ecog+f(7,i)
         tcog=tcog+f(6,i)
         xav=xav+f(2,i)
         xpav=xpav+f(3,i)
         yav=yav+f(4,i)
         ypav=ypav+f(5,i)
       enddo
       ecog=ecog/float(ngood)
       tcog=tcog/float(ngood)
       xav=xav/float(ngood)
       xpav=xpav/float(ngood)
       yav=yav/float(ngood)
       ypav=ypav/float(ngood)
C write particle coordinates to files,  phase in radian
       BEREF=VREF/VL
       GAMREF=1./SQRT(1.-(BEREF*BEREF))
       ENREF=XMAT*GAMREF
       open(58,file=wfile,status='unknown')
       if(iflg.lt.100) then
         write(58,*) ngood,beamc, fh/(2000000.*pi)
         if (chasit) then
           write(60,*) ngood,beamc, fh/(2000000.*pi)
           write(61,*) ngood,beamc, fh/(2000000.*pi)
         endif
       else
         write(58,*) ngood,beamc, fh/(2000000.*pi),enref-xmat
         if (chasit) then
           write(60,*) ngood,beamc, fh/(2000000.*pi),enref-xmat
           write(61,*) ngood,beamc, fh/(2000000.*pi),enref-xmat
         endif
       endif
       f2=0.
       f3=0.
       f4=0.
       f5=0.
       do i=1,ngood
         if(irec.eq.2) then
c coordinates relative to the reference
           eprt=f(7,i)-enref
           tprt=fh*(f(6,i)-tref)
           f2=f(2,i)
           f3=f(3,i)
           f4=f(4,i)
           f5=f(5,i)
         endif
         if(irec.eq.1) then
c absolute values for phase and energy
           eprt=f(7,i)-xmat
           tprt=fh*(f(6,i)-tcog)
c **********************************
           f2=f(2,i)
           f3=f(3,i)
           f4=f(4,i)
           f5=f(5,i)
         endif
         if(irec.eq.0) then
c values for phase and energy relative to the COG
           tprt=fh*(f(6,i)-tcog)
           eprt=f(7,i)-ecog
           f2=f(2,i)-xav
           f3=f(3,i)-xpav
           f4=f(4,i)-yav
           f5=f(5,i)-ypav
         endif
         if(iflg.eq.0 .or. iflg.eq.100)
     *    write(58,100)f2,f3/1000.,f4,f5/1000.,tprt,eprt
         if(iflg.eq.1 .or. iflg.eq.101)
     *    write(58,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(1,i)
         if(iflg.eq.2 .or. iflg.eq.102)
     *    write(58,102)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,i)
         if(iflg.eq.3 .or. iflg.eq.103)
     *    write(58,103)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,i),f(1,i)
         if(iflg.eq.10 .or. iflg.eq.110)
     *   write(58,100)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt
         if(iflg.eq.11 .or. iflg.eq.111)
     *   write(58,101)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(1,i)
         if(iflg.eq.12 .or. iflg.eq.112)
     *   write(58,102)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(9,i)
         if(iflg.eq.13 .or. iflg.eq.113)
     *   write(58,103)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,
     *   f(9,i),f(1,i)
c  ************** only with CHASE, write to file particles removed
         if (chasit) then
          if(ichxyz(i).eq.0) then
           if(iflg.eq.0 .or. iflg.eq.100)
     *   write(60,100)f2,f3/1000.,f4,f5/1000.,tprt,eprt
           if(iflg.eq.10 .or. iflg.eq.110)
     *   write(60,100)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt
           if(iflg.eq.1 .or. iflg.eq.101)
     *   write(60,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(1,i)
           if(iflg.eq.11 .or. iflg.eq.111)
     *   write(60,102)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(1,i)
           if(iflg.eq.2 .or. iflg.eq.102)
     *   write(60,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,i)
           if(iflg.eq.12 .or. iflg.eq.112)
     *   write(60,102)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(9,i)
          endif
c  ************** only with CHASE, write to file particles kept
          if(ichxyz(i).eq.1) then
           if(iflg.eq.0 .or. iflg.eq.100)
     *   write(61,100)f2,f3/1000.,f4,f5/1000.,tprt,eprt
           if(iflg.eq.10 .or. iflg.eq.110)
     *   write(61,100)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt
           if(iflg.eq.1 .or. iflg.eq.101)
     *   write(61,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(1,i)
           if(iflg.eq.11 .or. iflg.eq.111)
     *   write(61,102)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(1,i)
           if(iflg.eq.2 .or. iflg.eq.102)
     *   write(61,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,i)
           if(iflg.eq.12 .or. iflg.eq.112)
     *   write(61,102)f2,f3/1000.,f4,f5/1000.,1.e09*tprt/fh,eprt,f(9,i)
          endif
         endif
       enddo
       close(58)
100    format(6(e13.6,1x))
101    format(7(e13.6,1x))
102    format(7(e13.6,1x))
103    format(8(e13.6,1x))
       return
       end
       SUBROUTINE gausse
c  ........................................................................
c    called by HERSC
c    GAUSS method
c  ........................................................................
       implicit real*8 (a-h,o-z)
       common/gauss1/absg(40),wg(40),igaus
       common/randu/ck(15),kmax
       dimension ui(40),u9(9),u10(10),u12(12)
       dimension w9(9),w10(10),w12(12)
c   GAUSS n=9 de -1. a 1
       DATA (U9(J),J=1,9)/-.9681602,-.8360311,-.6133714,-.3242534,0.,
     1                     .3242534,.6133714,.8360311,.9681602/
       data (w9(j),j=1,9)/.0812744,.1806482,.2606107,.3123471,.3302394,
     1                     .3123471,.2606107,.1806482,.0812744/
c   GAUSS n=10 de -1. a 1
       DATA (U10(J),J=1,10)/-.9739065,-.8650634,-.6794096,-.4333954,
     1   -.1488743,.1488743,.4333954,.6794096,.8650634,.9739065/
       data (w10(j),j=1,10)/.0666713,.1494513,.2190864,.2692667,
     1   .2955242,.2955242,.2692667,.2190864,.1494513,.0666713/
c   GAUSS n=12 de -1. a 1
       DATA (U12(J),J=1,12)/-.9815606,-.9041173,-.7699027,-.5873180,
     1   -.3678315,-.1252334,.1252334,.3678315,.5873180,.7699027,
     2    .9041173,.9815606/
       data (w12(j),j=1,12)/.0471753,.1069393,.1600783,.2031674,
     1   .2334925,.2491470,.2491470,.2334925,.2031674,.1600783,
     2   .1069393,.0471753/
c    built the abscissas from (-1,1) to (1,0)
       if(igaus.eq.9) then
         do i=1,igaus
           ui(i)=u9(i)
           wg(i)=w9(i)
         enddo
       endif
c   following options in case one wishes to use 10 or 12 steps
       if(igaus.eq.10) then
         do i=1,igaus
           ui(i)=u10(i)
           wg(i)=w10(i)
         enddo
       endif
       if(igaus.eq.12) then
         do i=1,igaus
           ui(i)=u12(i)
           wg(i)=w12(i)
         enddo
       endif
       do i=1,igaus
         absg(i)=(1.+ui(i))/2.
         wg(i)=wg(i)/2.
       enddo
       return
       end
       SUBROUTINE table(lbmax,mbmax,nbmax)
c  .........................................................................
c   called by HERSC
c  arrays of variables independent of the coordinates of particles
c   lbmax,mbmax,nbmax are the maximum degrees of the coefficients A(l,m,n)
c  ..........................................................................
       implicit real*8 (a-h,o-z)
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/randu/ck(15),kmax
       common/hass/carg(100),sarg(100),argip(100)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/consta/vl,pi,xmat,rpel,qst
       common/factor/fpir(40,40),fect(30)
c       maximum of subscripts : ideg
c igaus: order of the Gauss integration ******
c  maximum of lbmax,mmax,nbmax -->ideg
       ideg=nbmax
       if(lbmax.ge.ideg) ideg=lbmax
       if(mbmax.ge.ideg) ideg=mbmax
       ideg=ideg+1
       idegp2=ideg+2
c  arrays of the power of the circular functions fo the Gauss positions absg(i)
c   the maximum power is ideg+2
       do i=1,igaus
         co(i,1)=1.
         sn(i,1)=1.
       enddo
       do i=1,igaus
         aco=pi*absg(i)/2.
         cod=cos(aco)
         snd=sin(aco)
         do j=2,idegp2
           co(i,j)=co(i,j-1)*cod
           sn(i,j)=sn(i,j-1)*snd
         enddo
       enddo
c    arrays used in the integrals of tables 75 and 76
c    absg(i)**2(t1+t2+t3+1)+1-->ragp(i,j) i: Gauss positions, j: is the power
       idgp=3*ideg+kmax
       kemax=kmax/2
       do i=1,igaus
         ragp(i,1)=1.
         ragm1(i,1)=1.
         absm2=(1.-absg(i))*(1.-absg(i))
         do j=2,idgp
           ragp(i,j)=ragp(i,j-1)*absg(i)
         enddo
c   storage of (absg(i)-1)**2*kemax -->ragm1
         do j=2,kemax+1
           ragm1(i,j)=ragm1(i,j-1)*absm2
         enddo
       enddo
       do i=1,idgp+1
         carg(i)=sqrt((4.*float(i-1)+1.)/2.)
         sarg(i)=sqrt((4.*float(i-1)+3.)/2.)
         argip(i)=sqrt(2.*float(i-1)+1.)
       enddo
c   store the factorials
       do i=1,40
         do j=1,i
           fpir(i,j)=fper(i-1,j-1)
         enddo
       enddo
       fj=1.
       fect(1)=1.
       do i=1,23
         fi=float(i)
         fj=fj*fi
         fect(i+1)=fj
       enddo
       return
       end
       FUNCTION fpar(i,j)
c  ....................................................................
c    Factorial function called by HERSC
c  ....................................................................
       implicit real*8 (a-h,o-z)
       common/factor/fpir(40,40),fect(30)
       ii=i+1
       jj=j+1
       fpar=fpir(ii,jj)
       return
       end
       FUNCTION fper(i,j)
c  ....................................................................
c      called by HERSC
c       i*(i-1)*(i-2)*...(i-j+1)/(1.2....j)
c  ....................................................................
       implicit real*8 (a-h,o-z)
       faci=1.
       facj=1.
       fper=1.
       if(i.eq.0) return
       if(j.eq.0) return
       do  k=1,j
         facj=facj*float(k)
       enddo
       ii=i
       do k=1,j
         faci=faci*float(ii)
         ii=ii-1
       enddo
       fper=faci/facj
       return
       end
       FUNCTION hermint(s,ihd)
c  .....................................................................
c    called by HERSC
c  .....................................................................
       implicit real*8 (a-h,o-z)
       dimension he(100)
       hermint=0.
       if(ihd.eq.0) then
         hermint=1.
         return
       endif
       if(ihd.eq.1) then
         hermint=s
         return
       endif
       he(1)=1.
       he(2)=s
       m1=ihd-1
       do k=1,m1
         he(k+2)=s*he(k+1)-float(k)*he(k)
       enddo
       hermint=he(ihd+1)*exp(-s*s/2.)
       return
       end
       FUNCTION fact(m)
c  ....................................................................
C      FACTORIAL of M
c  ....................................................................
       implicit real*8 (a-h,o-z)
       FACT=1.
       IF(M.EQ.0) RETURN
       DO  K=1,M
         FACT=FACT*FLOAT(K)
       ENDDO
       RETURN
       END
       FUNCTION factd(m)
c  ....................................................................
C      Calculate of (-1)**m *(2m-1)!!
c  ....................................................................
       implicit real*8 (a-h,o-z)
       dimension he(100)
       factd=1.
       if(m.eq.0) factd=1.
       if(m.eq.1) then
         factd=-1.
         return
       endif
       he(1)=1.
       do k=1,2*m-1,2
         he(k+2)=-float(k)*he(k)
       enddo
       factd=he(2*m+1)
       continue
       return
       end
       SUBROUTINE bhdist
c  ....................................................................
c    called by HERSC
c       computes: the coefficients A(l,m,n)
c                 the rms sizes in x, y and z-direction
c       selects the significants terms in the Hermite series expansion
c       lmax, mmax and nmax are the maximum values of l, m, n for these coefficients
c  ....................................................................
       implicit real*8 (a-h,o-z)
       common/coef/a(30,30,30),xa,xb,xc
       common/ind/lmax,mmax,nmax
       common/indin/lmaxi,mmaxi,nmaxxi
       common/hcgrms/xcdg,ycdg,zcdg,ect,eps
c       character iitime*30
c      hermite degrees
       lmax=lmaxi
       mmax=mmaxi
       nmax=nmaxxi
       ect=4.
c       call mytime(iitime)
c       write(16,*) 'PINTFAST started at ',iitime
       call pintfast
c       call mytime(iitime)
c       write(16,*) 'PINTFAST end, HCOEF started at ',iitime
       call hcoef
c       call mytime(iitime)
c       write(16,*) 'HCOEF ended at ',iitime
       return
       end
       SUBROUTINE trms(isucc)
c  ....................................................................
c     called by HERSC
c   storage  of the variables depending only of the rms sizes
c   lmax,mmax,nmax : maximum of the subscripts: l,m,n
c                    for the significants Almn
c  isucc:order of succession of the integrals in table 1-b
c     isucc=1: order of succesion x-->y-->z
c     isucc=2: order of succesion y-->z-->x
c     isucc=3: order of succesion z-->x-->y
c  ....................................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/coef/a(30,30,30),xa,xb,xc
       common/ind/lmax,mmax,nmax
c  lmax,mmax,nmax from zero
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/randu/ck(15),kmax
c  minimum of the rms sizes xa,xb,xc
       rmin=xa
       if(rmin.gt.xb) rmin=xb
       if(rmin.gt.xc) rmin=xc
c  order of the succession of the integrals
       if(rmin.eq.xc) isucc=3
       if(rmin.eq.xb) isucc=2
       if(rmin.eq.xa) isucc=1
c   maximum of lmax,mmax,nmax
       maxi=lmax
       if(maxi.lt.mmax) maxi=mmax
       if(maxi.lt.nmax) maxi=nmax
       maxt=2*(maxi+1)+3
       x3=0.
       x2=0.
       x1=0.
c  order of succession x->y->z
       if(isucc.eq.1) then
         x3=xa
         x2=xb
         x1=xc
       endif
c  order of succession y->z->x
       if(isucc.eq.2) then
        x3=xb
        x2=xc
        x1=xa
       endif
c  order of succession z->x->y
       if(isucc.eq.3) then
         x3=xc
         x2=xa
         x1=xb
       endif
c   array rms(j,i) Powers of the rms sizes  j=1,2,3, the value i is the power
       rms(3,1)=1.
       rms(2,1)=1.
       rms(1,1)=1.
       do i=2,maxt
         rms(3,i)=rms(3,i-1)*x3
         rms(2,i)=rms(2,i-1)*x2
         rms(1,i)=rms(1,i-1)*x1
       enddo
       j1=kmax/2
       i1m=(lmax+1)/2
       i2m=(mmax+1)/2
       i3m=(nmax+1)/2
       im=i1m+i2m+i3m+j1+4
       if(im.ge.40) then
         write(16,*) ' overlap the array blam with im= ',im
         stop
       endif
       do i=1,igaus
         blam(i,1)=(rms(1,3)*co(i,3)+rms(2,3)*sn(i,3))/rms(3,3)
         do ii=2,im
           blam(i,ii)=blam(i,ii-1)*blam(i,1)
         enddo
       enddo
       return
       end
       SUBROUTINE uvrms
c  ..........................................................................
c   called by HERSC
c     storage of the variables depending only of the coordinate s3
c  ..........................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *              akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/ind/lmax,mmax,nmax
       common/indttal/lmnt
       common/randu/ck(15),kmax
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       s32=s3*s3
       s22=s2*s2
       s12=s1*s1
       exs3=exp(-s32/2.)
       exs2=exp(-s22/2.)
       exs1=exp(-s12/2.)
       as3=abs(s3)
       sgns3=1.
       if(s3.lt.0.) sgns3=-1.
       as31=as3*s1
       as32=as3*s2
       s1rms=s1*rms(1,2)
       s2rms=s2*rms(2,2)
c  power of as3, used in the functions spii,......
       s3pw(1)=1.
       do i=2,kmax
         s3pw(i)=s3pw(i-1)*as3
       enddo
c    memory array hsint(ig,arg,indice) for functions sipp,...
       do ig=1,igaus
         wgpi=wg(ig)*pi/2.
         sqblam(ig)=sqrt(blam(ig,1))
         xlblam=sqblam(ig)*rms(3,2)
         do iarg=1,2
           arghm1=s1rms*co(ig,2)/xlblam
           arghm2=s2rms*sn(ig,2)/xlblam
           if(iarg.eq.1) then
             arg=arghm1+arghm2
             earg=exp(-arg*arg/2.)
           endif
           if(iarg.eq.2) then
           arg=arghm1-arghm2
             earg=exp(-arg*arg/2.)
           endif
           hsint(ig,iarg,1)=1.*earg*wgpi
           hsint(ig,iarg,2)=arg*earg*wgpi
           do inhs=3,lmnt
             hsint(ig,iarg,inhs)=arg*hsint(ig,iarg,inhs-1)
     *                     -float(inhs-2)*hsint(ig,iarg,inhs-2)
           enddo
         enddo
       enddo
       r13=rms(1,2)/rms(3,2)
       r23=rms(2,2)/rms(3,2)
       do j=1,igaus
         do i=1,igaus
           aeps1=ragp(i,3)*(blam(j,1)-1.)/2.
           aeps1=s32*(aeps1+ragp(i,2))
           aeps2=ragp(i,3)*blam(j,1)/2.
           aeps2=aeps2*s32
           akc1=ragp(i,2)*r13*co(j,2)*as31
           akc2=ragp(i,2)*r23*sn(j,2)*as32
           aks1=ragp(i,2)*r13*co(j,2)*as31
           aks2=ragp(i,2)*r23*sn(j,2)*as32
           epsi1(i,j)=exp(-aeps1)*wg(i)
           epsi2(i,j)=exp(-aeps2)*wg(i)
           akpcc(i,j)=cos(akc1)*cos(akc2)
           akpcs(i,j)=cos(akc1)*sin(akc2)
           akpsc(i,j)=sin(akc1)*cos(akc2)
           akpss(i,j)=sin(akc1)*sin(akc2)
         enddo
       enddo
c    Hermite functions
       hs1(1)=exs1
       hs2(1)=exs2
       hs3(1)=exs3
       hs1(2)=s1*exs1
       hs2(2)=s2*exs2
       hs3(2)=s3*exs3
       do ihe=3,lmnt
         hs1(ihe)=s1*hs1(ihe-1)-float(ihe-2)*hs1(ihe-2)
         hs2(ihe)=s2*hs2(ihe-1)-float(ihe-2)*hs2(ihe-2)
         hs3(ihe)=s3*hs3(ihe-1)-float(ihe-2)*hs3(ihe-2)
       enddo
       return
       end
       SUBROUTINE fielde(lc,mc,nc,isucc)
c********************************************************************
c     beam self-fields computation called by HERSC
c   look for the parity of the currents lc,mc and nc
c   compute the corresponding field components
c      isucc=1: a<b,c
c      isucc=2: b<a,c
c      isucc=3: c<a,b
c      x,y,z are the scaling coordinates: x/a,y/b,z/c
c the corresponding analytical equations are in the number of the tables
c   the table number in 'comments' refers to the corresponding analytical equations
c********************************************************************
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/field/ex,ey,ez
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),sgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1),pw3as3=as3**(2*it3))
       common/sgpth/mksgi,mksgp
       common/ftth/makti,maktp
       common/fsth/maksi,maksp
       logical maksi,maksp
       logical makti,maktp
       logical mksgi,mksgp
       logical lpl,lpm,lpn
c         look for the parity of the currents lc,mc,nc
       lpl=.false.
       lpm=.false.
       lpn=.false.
       makti=.false.
       maktp=.false.
       maksi=.false.
       maksp=.false.
       mksgi=.false.
       mksgp=.false.
       xlc=lc
       tl=xlc-2.*int(xlc/2.+0.0001)
       if(tl.eq.0.)lpl=.true.
       xmc=mc
       tm=xmc-2.*int(xmc/2.+0.0001)
       if(tm.eq.0.)lpm=.true.
       xnc=nc
       tn=xnc-2.*int(xnc/2.+0.0001)
       if(tn.eq.0.)lpn=.true.
       if(lpl.and.lpm.and.lpn) itpar=1
       if(lpl.and.lpm.and..not.lpn) itpar=2
       if(lpl.and..not.lpm.and.lpn) itpar=3
       if(lpl.and..not.lpm.and..not.lpn) itpar=4
       if(.not.lpl.and.lpm.and.lpn) itpar=5
       if(.not.lpl.and.lpm.and..not.lpn) itpar=6
       if(.not.lpl.and..not.lpm.and.lpn) itpar=7
       if(.not.lpl.and..not.lpm.and..not.lpn) itpar=8
c
c  values of the variables it1,it2,it3 in relation with lc,mc,nc
c  the initial values of these lc,mc,nc is zero
       if (isucc.eq.1) then
c  a<b,c
         if(itpar.eq.1) then
c       E E E and  a<b,c:
c                  table 67-a: Ex=E(2t3+1,2t2,2t1)
c                              Ey=E(2t3,2t2+1,2t1)
c                              Ez=E(2t3,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-a
           it3=lc/2
           it2=mc/2
           it1=nc/2
           if(s3.ne.0.) then
c            it=2*(it1+it2+it3)+1
             it=lc+mc+nc+1
             pwas3=as3**it
           endif
           ex=eipp(it1,it2,it3)
           ey=epip(it1,it2,it3)
           ez=eppi(it1,it2,it3)
           return
         endif
         if(itpar.eq.2) then
c       E E I and  a<b,c
c                  table 67-g: Ex=E(2t3+1,2t2,2t1+1)
c                              Ey=E(2t3,2t2+1,2t1+1)
c                              Ez=E(2t3,2t2,2t1)
c        values of it1,it2,it3 in table 3-g
c        x-direction and y-direction
           it3=lc/2
           it2=mc/2
           it1=(nc-1)/2
           if(s3.ne.0.) then
c	  it=2*(it1+it2+it3)+1
             it=lc+mc+nc
             pwas3=as3**it
           endif
           ex=eipi(it1,it2,it3)
           ey=epii(it1,it2,it3)
c        z-direction
           it1=(nc+1)/2
           if(s3.ne.0.)pwas3=pwas3*s32
c     it=2*(it1+it2+it3)+1
           ez=eppp(it1,it2,it3)
           return
         endif
         if(itpar.eq.3) then
c       E I E and  a<b,c
c               in table 67-h: Ex=E(2t3+1,2t2+1,2t1)
c                              Ey=E(2t3,2t2,2t1)
c                              Ez=E(2t3,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-h
c        x-direction
           it3=lc/2
           it2=(mc-1)/2
           it1=nc/2
           if(s3.ne.0.) then
c	 it=2*(it1+it2+it3)+1
             it=lc+mc+nc
             apwas3=as3**it
             pwas3=apwas3
           endif
           ex=eiip(it1,it2,it3)
c        y-direction
           it3=lc/2
           it2=(mc+1)/2
           it1=nc/2
c       it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=pwas3*s32
           ey=eppp(it1,it2,it3)
c        z-direction
           it3=lc/2
           it2=(mc-1)/2
           it1=nc/2
           if(s3.ne.0.)pwas3=apwas3
c       it=2*(it1+it2+it3)+1
           ez=epii(it1,it2,it3)
           return
         endif
         if(itpar.eq.4) then
c       E I I and  a<b,c
c               in table 67-f: Ex=E(2t3+1,2t2+1,2t1+1)
c                              Ey=E(2t3,2t2,2t1+1)
c                              Ez=E(2t3,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-f
c        x-direction
           it3=lc/2
           it2=(mc-1)/2
           it1=(nc-1)/2
           if(s3.ne.0.) then
c	it=2*(it1+it2+it3)+1
             it=lc+mc+nc-1
             pwas3=as3**it
           endif
           ex=eiii(it1,it2,it3)
c        y-direction
           it3=lc/2
           it2=(mc+1)/2
           it1=(nc-1)/2
c         it=2*(it1+it2+it3)+1
           if(s3.ne.0.)pwas3=pwas3*s32
           ey=eppi(it1,it2,it3)
c        z-direction
           it3=lc/2
           it2=(mc-1)/2
           it1=(nc+1)/2
           ez=epip(it1,it2,it3)
           return
         endif
         if(itpar.eq.5) then
c       I E E and  a<b,c
c               in table 67-c: Ex=E(2t3,2t2,2t1)
c                              Ey=E(2t3+1,2t2+1,2t1)
c                              Ez=E(2t3+1,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-c
c        x-direction
           it3=(lc+1)/2
           it2=mc/2
           it1=nc/2
           if(s3.ne.0.) then
c	it=2*(it1+it2+it3)+1
             it=lc+mc+nc
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eppp(it1,it2,it3)
c        y-direction
           it3=(lc-1)/2
           it2=mc/2
           it1=nc/2
           if(s3.ne.0.) then
c       it=2*(it1+it2+it3)+1
             pwas3=apwas3
             pw3as3=apwas3
           endif
           ey=eiip(it1,it2,it3)
c        z-direction
           it3=(lc-1)/2
           it2=mc/2
           it1=nc/2
           ez=eipi(it1,it2,it3)
           return
         endif
         if(itpar.eq.6) then
c       I E I and  a<b,c
c               in table 67-e: Ex=E(2t3,2t2,2t1+1)
c                              Ey=E(2t3+1,2t2+1,2t1+1)
c                              Ez=E(2t3+1,2t2,2t1)
c        values of it1,it2,it3 in table 3-e
c        x-direction
           it3=(lc+1)/2
           it2=mc/2
           it1=(nc-1)/2
           if(s3.ne.0.) then
c	it=2*(it1+it2+it3)+1
             it=lc+mc+nc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eppi(it1,it2,it3)
c        y-direction
           it3=(lc-1)/2
           it2=mc/2
           it1=(nc-1)/2
           if(s3.ne.0.) then
c             it=2*(it1+it2+it3)+1
             pwas3=apwas3
           endif
           ey=eiii(it1,it2,it3)
c        z-direction
           it3=(lc-1)/2
           it2=mc/2
           it1=(nc+1)/2
c      it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=pwas3*s32
           ez=eipp(it1,it2,it3)
           return
         endif
         if(itpar.eq.7) then
c       I I E and  a<b,c
c               in table 67-d: Ex=E(2t3,2t2+1,2t1)
c                              Ey=E(2t3+1,2t2,2t1)
c                              Ez=E(2t3+1,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-d
c        x-direction
           it3=(lc+1)/2
           it2=(mc-1)/2
           it1=nc/2
           if(s3.ne.0.) then
c  	it=2*(it1+it2+it3)+1
             it=lc+mc+nc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=epip(it1,it2,it3)
c        y-direction
           it3=(lc-1)/2
           it2=(mc+1)/2
           it1=nc/2
c             it=2*(it1+it2+it3)+1
           ey=eipp(it1,it2,it3)
c        z-direction
           it3=(lc-1)/2
           it2=(mc-1)/2
           it1=nc/2
c         it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3
           ez=eiii(it1,it2,it3)
           return
         endif
         if(itpar.eq.8) then
c       I I I and  a<b,c
c               in table 67-b: Ex=E(2t3,2t2+1,2t1+1)
c                              Ey=E(2t3+1,2t2,2t1+1)
c                              Ez=E(2t3+1,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-b
c        x-direction
           it3=(lc+1)/2
           it2=(mc-1)/2
           it1=(nc-1)/2
           if(s3.ne.0.) then
c        it=2*(it1+it2+it3)+1
             it=lc+mc+nc
             pwas3=as3**it
           endif
           ex=epii(it1,it2,it3)
c        y-direction
           it3=(lc-1)/2
           it2=(mc+1)/2
           it1=(nc-1)/2
c     it=2*(it1+it2+it3)+1
           ey=eipi(it1,it2,it3)
c        z-direction
           it3=(lc-1)/2
           it2=(mc-1)/2
           it1=(nc+1)/2
           ez=eiip(it1,it2,it3)
           return
         endif
c    endif of isucc=1----> a<b,c
       endif
       if (isucc.eq.2) then
c       b<a,c
         if(itpar.eq.1) then
c       E E E and  b<a,c
c                  table 67-a: Ex=E(2t3,2t2,2t1+1)
c                              Ey=E(2t3+1,2t2,2t1)
c                              Ez=E(2t3,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-a
           it3=mc/2
           it2=nc/2
           it1=lc/2
           if(s3.ne.0.) then
c 	it=2*(it1+it2+it3)+1
             it=mc+nc+lc+1
             pwas3=as3**it
           endif
           ex=eppi(it1,it2,it3)
           ey=eipp(it1,it2,it3)
           ez=epip(it1,it2,it3)
           return
         endif
         if(itpar.eq.2) then
c       E E I and  b<a,c
c                  table 67-h: Ex=E(2t3,2t2+1,2t1+1)
c                              Ey=E(2t3+1,2t2+1,2t1)
c                              Ez=E(2t3,2t2,2t1)
c        values of it1,it2,it3 in table 3-h
c        x-direction and y-direction
           it3=mc/2
           it2=(nc-1)/2
           it1=lc/2
           if(s3.ne.0) then
c	it=2*(it1+it2+it3)+1
             it=mc+nc+lc
             pwas3=as3**it
           endif
           ex=epii(it1,it2,it3)
           ey=eiip(it1,it2,it3)
c        z-direction
           it2=(nc+1)/2
c     it=2*(it1+it2+it3)+1
           if(s3.ne.0.)pwas3=pwas3*s32
           ez=eppp(it1,it2,it3)
           return
         endif
         if(itpar.eq.3) then
c       E I E and  b<a,c
c               in table 67-c: Ex=E(2t3+1,2t2,2t1+1)
c                              Ey=E(2t3,2t2,2t1)
c                              Ez=E(2t3+1,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-c
c        x-direction
           it3=(mc-1)/2
           it2=nc/2
           it1=lc/2
           if(s3.ne.0.) then
c     it=2*(it1+it2+it3)+1
             it=mc+nc+lc
             apwas3=as3**it
             pwas3= apwas3
           endif
           ex=eipi(it1,it2,it3)
c        y-direction
           it3=(mc+1)/2
           it2=nc/2
           it1=lc/2
           if(s3.ne.0.) then
c      it=2*(it1+it2+it3)+1
             pwas3=pwas3*s32
           endif
           ey=eppp(it1,it2,it3)
c        z-direction
           it3=(mc-1)/2
           it2=nc/2
           it1=lc/2
           if(s3.ne.0.) pwas3= apwas3
           ez=eiip(it1,it2,it3)
           return
         endif
         if(itpar.eq.4) then
c       E I I and  b<a,c
c               in table 67-d: Ex=E(2t3+1,2t2+1,2t1+1)
c                              Ey=E(2t3,2t2+1,2t1)
c                              Ez=E(2t3+1,2t2,2t1)
c        values of it1,it2,it3 in table 3-d
c        x-direction
           it3=(mc-1)/2
           it2=(nc-1)/2
           it1=lc/2
           if(s3.ne.0.) then
c      it=2*(it1+it2+it3)+1
             it=mc+nc+lc-1
             pwas3=as3**it
           endif
           ex=eiii(it1,it2,it3)
c        y-direction
           it3=(mc+1)/2
           it2=(nc-1)/2
           it1=lc/2
c    it=2*(it1+it2+it3)+1
           if(s3.ne.0.) then
             pwas3=pwas3*s32
           endif
           ey=epip(it1,it2,it3)
c        z-direction
           it3=(mc-1)/2
           it2=(nc+1)/2
           it1=lc/2
c             it=2*(it1+it2+it3)+1
           ez=eipp(it1,it2,it3)
           return
         endif
         if(itpar.eq.5) then
c       I E E and  b<a,c
c               in table 67-g: Ex=E(2t3,2t2,2t1)
c                              Ey=E(2t3+1,2t2,2t1+1)
c                              Ez=E(2t3,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-g
c        x-direction
           it3=mc/2
           it2=nc/2
           it1=(lc+1)/2
           if(s3.ne.0.) then
c       it=2*(it1+it2+it3)+1
             it=mc+nc+lc
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eppp(it1,it2,it3)
c        y-direction
           it3=mc/2
           it2=nc/2
           it1=(lc-1)/2
c    it=2*(it1+it2+it3)+1
           if(s3.ne.0.)pwas3=apwas3
           ey=eipi(it1,it2,it3)
c        z-direction
           it3=mc/2
           it2=nc/2
           it1=(lc-1)/2
           ez=epii(it1,it2,it3)
           return
         endif
         if(itpar.eq.6) then
c       I E I and  b<a,c
c               in table 67-f: Ex=E(2t3,2t2+1,2t1)
c                              Ey=E(2t3+1,2t2+1,2t1+1)
c                              Ez=E(2t3,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-f
c        x-direction
           it3=mc/2
           it2=(nc-1)/2
           it1=(lc+1)/2
           if(s3.ne.0.) then
c    it=2*(it1+it2+it3)+1
             it=mc+nc+lc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=epip(it1,it2,it3)
c        y-direction
           it3=mc/2
           it2=(nc-1)/2
           it1=(lc-1)/2
c      it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3
           ey=eiii(it1,it2,it3)
c        z-direction
           it3=mc/2
           it2=(nc+1)/2
           it1=(lc-1)/2
c      it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3*s32
           ez=eppi(it1,it2,it3)
           return
         endif
         if(itpar.eq.7) then
c       I I E and  b<a,c
c               in table 66-e: Ex=E(2t3+1,2t2,2t1)
c                              Ey=E(2t3,2t2,2t1+1)
c                              Ez=E(2t3+1,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-e
c        x-direction
           it3=(mc-1)/2
           it2=nc/2
           it1=(lc+1)/2
           if(s3.ne.0.) then
c    it=2*(it1+it2+it3)+1
             it=mc+nc+lc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eipp(it1,it2,it3)
c        y-direction
           it3=(mc+1)/2
           it2=nc/2
           it1=(lc-1)/2
c             it=2*(it1+it2+it3)+1
           ey=eppi(it1,it2,it3)
c        z-direction
           it3=(mc-1)/2
           it2=nc/2
           it1=(lc-1)/2
           if(s3.ne.0.) pwas3=apwas3
c       it=2*(it1+it2+it3)+1
           ez=eiii(it1,it2,it3)
          return
         endif
         if(itpar.eq.8) then
c       I I I and  b<a,c
c               in table 67-b: Ex=E(2t3+1,2t2+1,2t1)
c                              Ey=E(2t3,2t2+1,2t1+1)
c                              Ez=E(2t3+1,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-b
c        x-direction
           it3=(mc-1)/2
           it2=(nc-1)/2
           it1=(lc+1)/2
           if(s3.ne.0.) then
c    it=2*(it1+it2+it3)+1
             it=mc+nc+lc
             pwas3=as3**it
           endif
           ex=eiip(it1,it2,it3)
c        y-direction
           it3=(mc+1)/2
           it2=(nc-1)/2
           it1=(lc-1)/2
c     it=2*(it1+it2+it3)+1
           ey=epii(it1,it2,it3)
c        z-direction
           it3=(mc-1)/2
           it2=(nc+1)/2
           it1=(lc-1)/2
c             it=2*(it1+it2+it3)+1
           ez=eipi(it1,it2,it3)
           return
         endif
c    endif of isucc=2----> b<a,c
       endif
       if (isucc.eq.3) then
c       c<a,b
         if(itpar.eq.1) then
c       E E E and  c<a,b
c                  table 67-a: Ex=E(2t3,2t2+1,2t1)
c                              Ey=E(2t3,2t2,2t1+1)
c                              Ez=E(2t3+1,2t2,2t1)
c        values of it1,it2,it3 in table 3-a
           it3=nc/2
           it2=lc/2
           it1=mc/2
           if(s3.ne.0.) then
c    it=2*(it1+it2+it3)+1
             it=nc+lc+mc+1
             pwas3=as3**it
           endif
           ex=epip(it1,it2,it3)
           ey=eppi(it1,it2,it3)
           ez=eipp(it1,it2,it3)
           return
         endif
         if(itpar.eq.2) then
c       E E I and  c<a,b
c                  table 67-c: Ex=E(2t3+1,2t2+1,2t1)
c                              Ey=E(2t3+1,2t2,2t1+1)
c                              Ez=E(2t3,2t2,2t1)
c        values of it1,it2,it3 in table 3-c
c        x-direction and y-direction
           it3=(nc-1)/2
           it2=lc/2
           it1=mc/2
           if(s3.ne.0.) then
c     it=2*(it1+it2+it3)+1
             it=nc+lc+mc
             pwas3=as3**it
           endif
           ex=eiip(it1,it2,it3)
           ey=eipi(it1,it2,it3)
c        z-direction
           it3=(nc+1)/2
           if(s3.ne.0.) pwas3=pwas3*s32
c       it=2*(it1+it2+it3)+1
           ez=eppp(it1,it2,it3)
           return
         endif
         if(itpar.eq.3) then
c       E I E and  c<a,b
c               in table 67-g: Ex=E(2t3,2t2+1,2t1+1)
c                              Ey=E(2t3,2t2,2t1)
c                              Ez=E(2t3+1,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-g
c        x-direction
           it3=nc/2
           it2=lc/2
           it1=(mc-1)/2
           if(s3.ne.0.) then
c     it=2*(it1+it2+it3)+1
             it=nc+lc+mc
             apwas3=as3**it
             pwas3=apwas3
           endif
           ex=epii(it1,it2,it3)
c        y-direction
           it3=nc/2
           it2=lc/2
           it1=(mc+1)/2
c     it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=pwas3*s32
           ey=eppp(it1,it2,it3)
c        z-direction
           it3=nc/2
           it2=lc/2
           it1=(mc-1)/2
           if(s3.ne.0.) pwas3=apwas3
           ez=eipi(it1,it2,it3)
           return
         endif
         if(itpar.eq.4) then
c       E I I and  c<a,b
c               in table 67-e: Ex=E(2t3+1,2t2+1,2t1+1)
c                              Ey=E(2t3+1,2t2,2t1)
c                              Ez=E(2t3,2t2,2t1+1)
c        values of it1,it2,it3 in table 3-e
c        x-direction
           it3=(nc-1)/2
           it2=lc/2
           it1=(mc-1)/2
           if(s3.ne.0.) then
c   it=2*(it1+it2+it3)+1
             it=nc+lc+mc-1
             pwas3=as3**it
           endif
           ex=eiii(it1,it2,it3)
c        y-direction
           it3=(nc-1)/2
           it2=lc/2
           it1=(mc+1)/2
c      it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=pwas3*s32
           ey=eipp(it1,it2,it3)
c        z-direction
           it3=(nc+1)/2
           it2=lc/2
           it1=(mc-1)/2
c        it=2*(it1+it2+it3)+1
           ez=eppi(it1,it2,it3)
           return
         endif
         if(itpar.eq.5) then
c       I E E and  c<a,b
c               in table 67-h: Ex=E(2t3,2t2,2t1)
c                              Ey=E(2t3,2t2+1,2t1+1)
c                              Ez=E(2t3+1,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-h
c        x-direction
           it3=nc/2
           it2=(lc+1)/2
           it1=mc/2
           if(s3.ne.0.) then
c        it=2*(it1+it2+it3)+1
             it=nc+lc+mc
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eppp(it1,it2,it3)
c        y-direction
           it3=nc/2
           it2=(lc-1)/2
           it1=mc/2
c         it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3
           ey=epii(it1,it2,it3)
c        z-direction
           it3=nc/2
           it2=(lc-1)/2
           it1=mc/2
           ez=eiip(it1,it2,it3)
           return
         endif
         if(itpar.eq.6) then
c       I E I and  c<a,b
c               in table 67-d: Ex=E(2t3+1,2t2,2t1)
c                              Ey=E(2t3+1,2t2+1,2t1+1)
c                              Ez=E(2t3,2t2+1,2t1)
c        values of it1,it2,it3 in table 3-d
c        x-direction
           it3=(nc-1)/2
           it2=(lc+1)/2
           it1=mc/2
           if(s3.ne.0.) then
c        it=2*(it1+it2+it3)+1
             it=nc+lc+mc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eipp(it1,it2,it3)
c        y-direction
           it3=(nc-1)/2
           it2=(lc-1)/2
           it1=mc/2
c       it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3
           ey=eiii(it1,it2,it3)
c        z-direction
           it3=(nc+1)/2
           it2=(lc-1)/2
           it1=mc/2
c       it=2*(it1+it2+it3)+1
           if(s3.ne.0.)  pwas3=pwas3*s32
           ez=epip(it1,it2,it3)
           return
         endif
         if(itpar.eq.7) then
c       I I E and  c<a,b
c               in table 67-f: Ex=E(2t3,2t2,2t1+1)
c                              Ey=E(2t3,2t2+1,2t1)
c                              Ez=E(2t3+1,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-f
c        x-direction
           it3=nc/2
           it2=(lc+1)/2
           it1=(mc-1)/2
           if(s3.ne.0.) then
c          it=2*(it1+it2+it3)+1
             it=nc+lc+mc-1
             apwas3=as3**it
             pwas3=apwas3*s32
           endif
           ex=eppi(it1,it2,it3)
c        y-direction
           it3=nc/2
           it2=(lc-1)/2
           it1=(mc+1)/2
c       it=2*(it1+it2+it3)+1
           ey=epip(it1,it2,it3)
c        z-direction
           it3=nc/2
           it2=(lc-1)/2
           it1=(mc-1)/2
c      it=2*(it1+it2+it3)+1
           if(s3.ne.0.) pwas3=apwas3
           ez=eiii(it1,it2,it3)
           return
         endif
         if(itpar.eq.8) then
c       I I I and  c<a,b
c               in table 67-b: Ex=E(2t3+1,2t2,2t1+1)
c                              Ey=E(2t3+1,2t2+1,2t1)
c                              Ez=E(2t3,2t2+1,2t1+1)
c        values of it1,it2,it3 in table 3-b
c        x-direction
         it3=(nc-1)/2
           it2=(lc+1)/2
           it1=(mc-1)/2
c        it=2*(it1+it2+it3)+1
           if(s3.ne.0.) then
             it=nc+lc+mc
             pwas3=as3**it
           endif
           ex=eipi(it1,it2,it3)
c        y-direction
           it3=(nc-1)/2
           it2=(lc-1)/2
           it1=(mc+1)/2
c        it=2*(it1+it2+it3)+1
           ey=eiip(it1,it2,it3)
c        z-direction
           it3=(nc+1)/2
           it2=(lc-1)/2
           it1=(mc-1)/2
c        it=2*(it1+it2+it3)+1
           ez=eppi(it1,it2,it3)
           return
         endif
c    endif of isucc=3----> c<a,b
       endif
       end
       FUNCTION eppp(it1,it2,it3)
c   .................................................................
c      E(2it3,2it2,2it1)   table 77-a-1
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/randu/ck(15),kmax
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),isgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1)
       eppp=0.
       isgnw=4*(it2+it1)+5*it3
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eppp=wsng*pi2*exs3*rms(3,3)
       eppp=eppp*(tppp(it1,it2,it3)+sppp(it1,it2,it3))
       eppp=eppp-8.*wsng*rms(3,3)*sqpi*sgppp(it1,it2,it3)
       if(s3.ne.0.) then
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+1
         ipt22=2*it2+1
         ipt212=2*(it1+it2)
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas3*rint*uppp(it1,it2,it3)
         scum=0.
         kj=1
         pcas3=pwas3
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vppp(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3
         eppp=eppp+e1+e2
       endif
       return
       end
       FUNCTION epip(it1,it2,it3)
c   .................................................................
c    E(2it3,2it2+1,2it1)  table 77-a-2
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),isgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1)
       epip=0.
       isgnw=4*(it2+it1)+5*it3+2
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       epip=wsng*pi2*exs3*rms(3,3)
       epip=epip*(tpip(it1,it2,it3)+spip(it1,it2,it3))
       epip=epip-8.*wsng*rms(3,3)*sqpi*sgpip(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*as3
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+2
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1.
         endif
         ipt12=2*it1+1
         ipt22=2*it2+2
         ipt212=2*(it1+it2)+1
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*upip(it1,it2,it3)
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vpip(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3
         epip=epip+(e1+e2)
       endif
       return
       end
       FUNCTION eppi(it1,it2,it3)
c   .................................................................
c    E(2it3,2it2,2it1+1) table 77-a-2
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),isgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1),pw3as3=as3**(2*it3)
       eppi=0.
       isgnw=4*(it2+it1)+5*it3+2
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eppi=wsng*pi2*exs3*rms(3,3)
       eppi=eppi*(tppi(it1,it2,it3)+sppi(it1,it2,it3))
       eppi=eppi-8.*wsng*rms(3,3)*sqpi*sgppi(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*as3
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+2
         xsgn1=isgn1
         pari=xsgn1-2.*int(xsgn1/2.+0.0001)
         if(pari.eq.0.) then
           sgn1=1.
           sgn2=-1.
         endif
         ipt12=2*it1+2
         ipt22=2*it2+1
         ipt212=2*(it1+it2)+1
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*uppi(it1,it2,it3)
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vppi(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3
         eppi=eppi+(e1+e2)
       endif
       return
       end
       FUNCTION epii(it1,it2,it3)
c   .................................................................
c     E(2*it3,2*it2+1,2*it1+1) table 77-a-1
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),isgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1),pw3as3=as3**(2*it3)
       epii=0.
       isgnw=4*(it2+it1)+5*it3+4
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       epii=wsng*pi2*exs3*rms(3,3)
       epii=epii*(tpii(it1,it2,it3)+spii(it1,it2,it3))
       epii=epii-8.*wsng*rms(3,3)*sqpi*sgpii(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*s32
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+4
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+2
         ipt22=2*it2+2
         ipt212=2*(it1+it2+1)
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*upii(it1,it2,it3)
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vpii(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3
         epii=epii+e1+e2
       endif
       return
       end
       FUNCTION eipp(it1,it2,it3)
c   .................................................................
c      E(2it3+1,2it2,2it1) table 77-b-2
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),sgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1),pw3as3=as3**(2*it3)
       eipp=0.
       isgnw=4*(it2+it1)+5*it3+3
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eipp=wsng*sgns3*pi2*exs3*rms(3,3)
       eipp=eipp*(tipp(it1,it2,it3)+sipp(it1,it2,it3))
       eipp=eipp+8.*wsng*rms(3,3)*sqpi*sgipp(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*as3
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+3
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+1
         ipt22=2*it2+1
         ipt212=2*(it1+it2)
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*uipp(it1,it2,it3)*sgns3
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vipp(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3*sgns3
         eipp=eipp-(e1+e2)
       endif
       return
       end
       FUNCTION eiip(it1,it2,it3)
c   .................................................................
c    E(2it3+1,2it2+1,2it1)  table 77-b-1
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),sgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1)
       eiip=0.
       isgnw=4*(it2+it1)+5*it3+5
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eiip=wsng*sgns3*pi2*exs3*rms(3,3)
       eiip=eiip*(tiip(it1,it2,it3)+siip(it1,it2,it3))
       eiip=eiip+8.*wsng*rms(3,3)*sqpi*sgiip(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*s32
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+5
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+1
         ipt22=2*it2+2
         ipt212=2*(it1+it2)+1
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*uiip(it1,it2,it3)*sgns3
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*viip(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3*sgns3
         eiip=eiip-(e1+e2)
       endif
       return
       end
       FUNCTION eipi(it1,it2,it3)
c   .................................................................
c     E(2*it3+1,2*it2,2*it1+1)  table 77-b-1
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),sgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1)
       eipi=0.
       isgnw=4*(it2+it1)+5*it3+5
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eipi=wsng*sgns3*pi2*exs3*rms(3,3)
       eipi=eipi*(tipi(it1,it2,it3)+sipi(it1,it2,it3))
       eipi=eipi+8.*wsng*rms(3,3)*sqpi*sgipi(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*s32
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+5
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+2
         ipt22=2*it2+1
         ipt212=2*(it1+it2)+1
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*uipi(it1,it2,it3)*sgns3
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*vipi(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3*sgns3
         eipi=eipi-(e1+e2)
       endif
       return
       end
       FUNCTION eiii(it1,it2,it3)
c   .................................................................
c    E(2*it3+1,2*it2+1,2*it1+1) table 77-b-2
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/partcl/x,y,z
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/randu/ck(15),kmax
c           pi2=pi*pi, sqpi=(pi/2)**3/2
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
c         exs3=exp(-s3*s3/2),abs3=abs(s3),sgns3=sign s3,s32=s3*s3
c         pwas3=as3**(2(it1+it2+it3)+1),pw3as3=as3**(2*it3)
       eiii=0.
       isgnw=4*(it2+it1)+5*it3+5
       ipar=isgnw-2*(isgnw/2)
       wsng=-1.
       if(ipar.eq.0)wsng=1.
       eiii=wsng*sgns3*pi2*exs3*rms(3,3)
       eiii=eiii*(tiii(it1,it2,it3)+siii(it1,it2,it3))
       eiii=eiii+8.*wsng*rms(3,3)*sqpi*sgiii(it1,it2,it3)
       if(s3.ne.0.) then
         pwas31=pwas3*s32*as3
         sgn1=-1.
         sgn2=1.
         isgn1=3*(it1+it2)+4*it3+5
         ipar=isgn1-2*(isgn1/2)
         if(ipar.eq.0) then
           sgn1=1.
           sgn2=-1
         endif
         ipt12=2*it1+2
         ipt22=2*it2+2
         ipt212=2*(it1+it2)+2
         rint=rms(1,ipt12+1)*rms(2,ipt22+1)/rms(3,ipt212+1)
         e1=2.*sgn1*pi2*pwas31*rint*uiii(it1,it2,it3)*sgns3
         scum=0.
         kj=1
         pcas3=pwas31
         do k=1,kmax,2
           scum=scum+ck(k)*pcas3*viii(kj)
           kj=kj+1
           pcas3=pcas3*s32
         enddo
         e2=scum*sgn2*2.*pi2*rint*exs3*sgns3
         eiii=eiii-(e1+e2)
       endif
       return
       end
       FUNCTION tipp(it1,it2,it3)
c   .................................................................
c   Part of W**(2t3+1,2t2,2t1)  table 61-b
c    with T(2t3+1,2t2,2t1,j1)   table 41
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tipp=0.
       if(.not.makti) then
         makti=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           xsj1=1.
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
         do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1i(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+1)
               mm=mm+1
             enddo
c      s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2i(kk,jj1)=2.*fpar(km1,2*jm1)*stoc
             tt1=tt1+tt*stc2i(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tipp=tt1*ck(k)+tipp
           kk=kk+1
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+1)
               mm1=mm1+1
             enddo
             tt1=tt1+tt*stc2i(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tipp=tt1*ck(k)+tipp
           kk=kk+1
         enddo
         return
       endif
       end
       FUNCTION tiip(it1,it2,it3)
c   .................................................................
c    part of W**(2t3+1,2t2+1,2t1)  table 61-b
c      with  T(2t3+1,2t2+1,2t1,j1) table 41
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tiip=0.
       if(.not.makti) then
         makti=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           xsj1=1.
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1i(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+1)
               mm=mm+1
             enddo
c      s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2i(kk,jj1)=2.*fpar(km1,2*jm1)*stoc
             tt1=tt1+tt*stc2i(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tiip=tt1*ck(k)+tiip
           kk=kk+1
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+1)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2i(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tiip=tt1*ck(k)+tiip
           kk=kk+1
         enddo
         return
       endif
       end
       FUNCTION tipi(it1,it2,it3)
c   .................................................................
c   part of W**(2t3+1,2t2,2t1+1) in table 61-b
c     with T(2t3+1,2t2,2t1+1,j1) in table 41
c   .................................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tipi=0.
       if(.not.makti) then
         makti=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           xsj1=1.
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2+1,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1i(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+2)
               mm=mm+1
             enddo
c      s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2i(kk,jj1)=2.*fpar(km1,2*jm1)*stoc
             tt1=tt1+tt*stc2i(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tipi=tt1*ck(k)+tipi
           kk=kk+1
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+2)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2i(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tipi=tt1*ck(k)+tipi
           kk=kk+1
         enddo
         return
       endif
       end
       FUNCTION tiii(it1,it2,it3)
c ...............................................................
c   part of W**(2t3+1,2t2+1,2t1+1) in table 61-b
c     with T(2t3+1,2t2+1,2t1+1,j1) in table 41
c ...............................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tiii=0.
       if(.not.makti) then
         makti=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           xsj1=1.
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2+1,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1i(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+2)
               mm=mm+1
             enddo
c      s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2i(kk,jj1)=2.*fpar(km1,2*jm1)*stoc
             tt1=tt1+tt*stc2i(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tiii=tt1*ck(k)+tiii
           kk=kk+1
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1
c   term T(2t3+1,2t2+1,2t1+1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1i(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+2)
               mm1=mm1+1
             enddo
             tt1=tt1+tt*stc2i(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tiii=tt1*ck(k)+tiii
           kk=kk+1
         enddo
         return
       endif
       end
       FUNCTION tppp(it1,it2,it3)
c  ...............................................................
c   part of W**(2t3,2t2,2t1) in table 61-a
c        with  T(2t3,2t2,2t1,j1) in table 41
c  ...............................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tppp=0.
       if(.not.maktp) then
         maktp=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           xsj1=1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 100
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1p(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+1)
               mm=mm+1
             enddo
c    s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2p(kk,jj1)=2.*fpar(km1,2*jm1+1)*stoc
             tt1=tt1+tt*stc2p(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tppp=tt1*ck(k)+tppp
           kk=kk+1
100        continue
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 200
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+1)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2p(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tppp=tt1*ck(k)+tppp
           kk=kk+1
200        continue
         enddo
         return
       endif
       end
       FUNCTION tpip(it1,it2,it3)
c  ...............................................................
c   the part of W**(2t3,2t2,2t1+1) in table 61-a
c        with T(2t3,2t2,2t1+1,j1)  in table 41
c  ...............................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tpip=0.
       if(.not.maktp) then
         maktp=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           xsj1=1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 100
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1p(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+1)
               mm=mm+1
             enddo
c    s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2p(kk,jj1)=2.*fpar(km1,2*jm1+1)*stoc
             tt1=tt1+tt*stc2p(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tpip=tt1*ck(k)+tpip
           kk=kk+1
100        continue
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 200
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+1)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2p(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tpip=tt1*ck(k)+tpip
           kk=kk+1
200        continue
         enddo
         return
       endif
       end
       FUNCTION tppi(it1,it2,it3)
c   .............................................................
c    part of W**(2t3,2t2,2t1+1) in table 61-a
c        with T(2t3,2t2,2t1+1,j1) in table 41
c   .............................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tppi=0.
       if(.not.maktp) then
         maktp=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           xsj1=1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 100
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1p(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+2)
               mm=mm+1
             enddo
c    s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2p(kk,jj1)=2.*fpar(km1,2*jm1+1)*stoc
             tt1=tt1+tt*stc2p(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tppi=tt1*ck(k)+tppi
           kk=kk+1
100        continue
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 200
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+2)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2p(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tppi=tt1*ck(k)+tppi
           kk=kk+1
200        continue
         enddo
         return
       endif
       end
       FUNCTION tpii(it1,it2,it3)
c  ............................................................
c       part of W**(2t3,2t2+1,2t1+1) in table 61-a
c         with T(2t3,2t2+1,2t1+1,j1) in table 41
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/ftsk/stc1i(8,8,40),stc1p(8,8,40),stc2i(8,8),stc2p(8,8)
       common/ftth/makti,maktp
       logical makti,maktp
       tpii=0.
       if(.not.maktp) then
         maktp=.true.
         kk=1
         do k=2,kmax,2
           km1=k-1
           xsj1=1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 100
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               stoc=rms(1,2*mm1+1)*rms(2,2*(it3j1-mm1)+1)
               stc1p(kk,jj1,mm)=fpar(it3j1,mm1)/stoc
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+2)*hs1(it1m+2)
               mm=mm+1
             enddo
c    s3pw(j1km+1)=as3**j1km
             stoc=rms(3,2*it3j1+1)*s3pw(j1km+1)*xsj1
             stc2p(kk,jj1)=2.*fpar(km1,2*jm1+1)*stoc
             tt1=tt1+tt*stc2p(kk,jj1)
             xsj1=-xsj1
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tpii=tt1*ck(k)+tpii
           kk=kk+1
100        continue
         enddo
         return
       else
         kk=1
         do k=2,kmax,2
           km1=k-1
           tt1=0.
           jm1=0
           jj1=1
           do j1=1,km1,2
             j1km=km1-2*jm1-1
c     j1km must be always greather or equal to zero
             if(j1km.lt.0) go to 200
c    T(2t3,2t2,2t1,j1) in table 41
             it3j1=it3+jm1
             tt=0.
             mm=1
             do m=1,it3j1+1
               mm1=m-1
               it23jm=2*(it2+it3+jm1-mm1)
               it1m=2*(it1+mm1)
               tt=tt+stc1p(kk,jj1,mm)*hs2(it23jm+1)*hs1(it1m+2)
               mm=mm+1
             enddo
             tt1=tt1+tt*stc2p(kk,jj1)
             jm1=jm1+1
             jj1=jj1+1
           enddo
           tpii=tt1*ck(k)+tpii
           kk=kk+1
200        continue
         enddo
         return
       endif
       end
       FUNCTION sipp(it1,it2,it3)
c  ............................................................
c   second part of W**(2t3+1,2t2,2t1) in table 61-b
c     with S(2t3+1,2t2,2t1,j1)  in table 70-b
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/hass/carg(100),sarg(100),argip(100)
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       sipp=0.
       it1p=2*it1+1
       it2p=2*it2+1
       it12p=2*(it1+it2+1)+1
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         xsj1=-1.
         km1=k-1
         tt1=0.
         jm1=0
         jj2=1
         do j1=1,km1,2
           jj1=2*jm1+1
           j1km=km1-jj1
           if(j1km.lt.0) go to 100
c    S(2t3+1,2t2,2t1,j1) computation
c      Gauss quadrature in table 70-b
           i1123j1=it1+it2+it3+jm1+1
           i2123j1=2*i1123j1
           tt=0.
           do ig=1,igaus
c Hermite functions are in the table hsint(ig,,ind )
c   caution!!! in the table hsint the indice ind is starting from 1
c              i2123j1 is starting from zero
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base/sqblam(ig)*(htm0+htp0)
           enddo
           if(.not.maksi) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstci(kk,jj2)=fpar(km1,jj1)*stock
             xsj1=-xsj1
           endif
           tt1=sstci(kk,jj2)*tt+tt1
           jj2=jj2+1
           jm1=jm1+1
         enddo
         sipp=sipp+tt1*ck(k)
         kk=kk+1
c  enddo from k (k=1,kmax+1)
100      continue
       enddo
       sipp=sipp*bsp
       maksi=.true.
       return
       end
       FUNCTION siip(it1,it2,it3)
c  ............................................................
c     W**(2t3+1,2t2+1,2t1) in table 61-b
c     with S(2t3+1,2t2+1,2t1,j1)in table 70-b
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       siip=0.
       it1p=2*it1+1
       it2p=2*it2+2
       it12p=2*(it1+it2+2)
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         xsj1=-1.
         km1=k-1
         tt1=0.
         jm1=0
         jj2=1
         do j1=1,km1,2
           jj1=2*jm1+1
           j1km=km1-jj1
           if(j1km.lt.0) go to 100
c    S(2t3+1,2t2,2t1,j1) computation
c      Gauss quadrature in table 70-b
           i1123j1=it1+it2+it3+jm1+2
           i2123j1=2*i1123j1-1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htp0-htm0)
           enddo
           if(.not.maksi) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstci(kk,jj2)=fpar(km1,jj1)*stock
             xsj1=-xsj1
           endif
           tt1=sstci(kk,jj2)*tt+tt1
           jj2=jj2+1
           jm1=jm1+1
         enddo
         siip=siip+tt1*ck(k)
         kk=kk+1
c  enddo from k (k=1,kmax+1)
100      continue
       enddo
       siip=siip*bsp
       maksi=.true.
       return
       end
       FUNCTION sipi(it1,it2,it3)
c  ............................................................
c  part of W**(2t3+1,2t2,2t1+1) in table 61-b
c   with S(2t3+1,2t2,2t1+1,j1) given in table 70-b
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       sipi=0.
       it1p=2*it1+2
       it2p=2*it2+1
       it12p=2*(it1+it2+2)
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         xsj1=-1.
         km1=k-1
         tt1=0.
         jm1=0
         jj2=1
         do j1=1,km1,2
           jj1=2*jm1+1
           j1km=km1-jj1
           if(j1km.lt.0) go to 100
c    S(2t3+1,2t2,2t1,j1) computation
c      Gauss quadrature in table 70-b
           i1123j1=it1+it2+it3+jm1+2
           i2123j1=2*i1123j1-1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htp0+htm0)
           enddo
           if(.not.maksi) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstci(kk,jj2)=fpar(km1,jj1)*stock
             xsj1=-xsj1
           endif
           tt1=sstci(kk,jj2)*tt+tt1
           jj2=jj2+1
           jm1=jm1+1
         enddo
         sipi=sipi+tt1*ck(k)
         kk=kk+1
c  enddo from k (k=1,kmax+1)
100      continue
       enddo
       sipi=sipi*bsp
       maksi=.true.
       return
       end
       FUNCTION siii(it1,it2,it3)
c  ............................................................
c   part of W**(2t3+1,2t2+1,2t1+1) in table 61-b
c    with S(2t3+1,2t2+1,2t1+1,j1)in table 70-b
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       siii=0.
       it1p=2*it1+2
       it2p=2*it2+2
       it12p=2*(it1+it2+2)+1
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         xsj1=-1.
         km1=k-1
         tt1=0.
         jm1=0
         jj2=1
         do j1=1,km1,2
           jj1=2*jm1+1
           j1km=km1-jj1
           if(j1km.lt.0) go to 100
c    S(2t3+1,2t2,2t1,j1) computation
c      Gauss quadrature in table 70-b
           i1123j1=it1+it2+it3+jm1+2
           i2123j1=2*i1123j1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htm0-htp0)/sqblam(ig)
           enddo
           if(.not.maksi) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstci(kk,jj2)=fpar(km1,jj1)*stock
             xsj1=-xsj1
           endif
           tt1=sstci(kk,jj2)*tt+tt1
           jj2=jj2+1
           jm1=jm1+1
         enddo
         siii=siii+tt1*ck(k)
c  enddo from k (k=1,kmax+1)
         kk=kk+1
100      continue
       enddo
       siii=-siii*bsp
       maksi=.true.
       return
       end
       FUNCTION sppp(it1,it2,it3)
c  ............................................................
c      part of W**(2t3,2t2,2t1) in table 61-a
c        with S(2t3,2t2,2t1,j1) given in table 70-a
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       sppp=0.
       it1p=2*it1+1
       it2p=2*it2+1
       it12p=2*(it1+it2+1)+1
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         km1=k-1
         xsj1=1
         tt1=0.
         jj1=1
         do j1=1,k
           jm1=j1-1
           j1km=km1-2*jm1
           if(j1km.lt.0) go to 100
c     j1km must be always greather or equal to zero
c      Gauss quadrature in table 70-a
           i1123j1=it1+it2+it3+jm1
           i2123j1=2*i1123j1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htm0+htp0)/sqblam(ig)
           enddo
           if(.not.maksp) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstcp(kk,jj1)=fpar(km1,2*jm1)*stock
             xsj1=-xsj1
           endif
           tt1=sstcp(kk,jj1)*tt+tt1
c  enddo from  j1
           jj1=jj1+1
         enddo
100      continue
         sppp=sppp+tt1*ck(k)
c  enddo from k (k=1,kmax+1)
          kk=kk+1
       enddo
       sppp=sppp*bsp
       maksp=.true.
       return
       end
       FUNCTION spip(it1,it2,it3)
c  ............................................................
c      part of W**(2t3,2t2,2t1+1) in table 61-a
c     with S(2t3,2t2,2t1+1,j1) given in table 70-a
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       spip=0.
       it1p=2*it1+1
       it2p=2*it2+2
       it12p=2*(it1+it2+1)+2
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         km1=k-1
         xsj1=1
         tt1=0.
         jj1=1
         do j1=1,k
           jm1=j1-1
           j1km=km1-2*jm1
           if(j1km.lt.0) go to 100
c     j1km must be always greather or equal to zero
c      Gauss quadrature in table 70-a
           i1123j1=it1+it2+it3+jm1
           i2123j1=2*i1123j1+1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1+1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htp0-htm0)
           enddo
           if(.not.maksp) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstcp(kk,jj1)=fpar(km1,2*jm1)*stock
             xsj1=-xsj1
           endif
           tt1=sstcp(kk,jj1)*tt+tt1
c  enddo from  j1
           jj1=jj1+1
         enddo
100      continue
         spip=spip+tt1*ck(k)
         kk=kk+1
c  enddo from k (k=1,kmax+1)
       enddo
       spip=spip*bsp
       maksp=.true.
       return
       end
       FUNCTION sppi(it1,it2,it3)
c  ............................................................
c      part of W**(2t3,2t2,2t1+1) in table 61-a
c     with S(2t3,2t2,2t1+1,j1) given in table 70-a
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       sppi=0.
       it1p=2*it1+2
       it2p=2*it2+1
       it12p=2*(it1+it2+1)+2
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         km1=k-1
         xsj1=1.
         tt1=0.
         jj1=1
         do j1=1,k
           jm1=j1-1
           j1km=km1-2*jm1
c     j1km must be always greather or equal to zero
           if(j1km.lt.0) go to 100
c      Gauss quadrature in table 70-a
           i1123j1=it1+it2+it3+jm1
           i2123j1=2*i1123j1+1
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1+1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htm0+htp0)
           enddo
           if(.not.maksp) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstcp(kk,jj1)=fpar(km1,2*jm1)*stock
             xsj1=-xsj1
           endif
           tt1=sstcp(kk,jj1)*tt+tt1
c  enddo from  j1
           jj1=jj1+1
         enddo
100      continue
         sppi=sppi+tt1*ck(k)
c  enddo from k (k=1,kmax+1)
         kk=kk+1
       enddo
       sppi=sppi*bsp
       maksp=.true.
       return
       end
       FUNCTION spii(it1,it2,it3)
c  ............................................................
c       part of W**(2t3,2t2,2t1+1) in table 61-a
c         with S(2t3,2t2+1,2t1+1,j1) given in table 70-a
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/randu/ck(15),kmax
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/rms/rms(3,50),s1,s2,s3
       common/comtab/hsint(40,2,60),sqblam(40),s3pw(15)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/fssk/sstci(8,8),sstcp(8,8)
       common/fsth/maksi,maksp
       logical maksi,maksp
       spii=0.
       it1p=2*it1+2
       it2p=2*it2+2
       it12p=2*(it1+it2+2)+1
       bsp=rms(1,it1p+1)*rms(2,it2p+1)/(sq2pi*rms(3,it12p))
       kk=1
       do k=1,kmax,2
         km1=k-1
         xsj1=1.
         tt1=0.
         jj1=1
         do j1=1,k
           jm1=j1-1
           j1km=km1-2*jm1
           if(j1km.lt.0) go to 100
c     j1km must be always greather or equal to zero
c      Gauss quadrature in table 70-a
           i1123j1=it1+it2+it3+jm1
           i2123j1=2*i1123j1+2
           tt=0.
           do ig=1,igaus
c Hermite functions
             base=co(ig,it1p)*sn(ig,it2p)/blam(ig,i1123j1+1)
             htm0=hsint(ig,2,i2123j1+1)
             htp0=hsint(ig,1,i2123j1+1)
             tt=tt+base*(htm0-htp0)/sqblam(ig)
           enddo
           if(.not.maksp) then
c   s3pw(j1km+1)=as3**j1km
             stock=2.*xsj1*s3pw(j1km+1)
             sstcp(kk,jj1)=fpar(km1,2*jm1)*stock
             xsj1=-xsj1
           endif
           tt1=sstcp(kk,jj1)*tt+tt1
c  enddo from  j1
           jj1=jj1+1
         enddo
100      continue
         spii=spii+tt1*ck(k)
         kk=kk+1
c  enddo from k (k=1,kmax+1)
       enddo
       maksp=.true.
       spii=-spii*bsp
       return
       end
c****************************************************************
c  funtions sigma in table 14
c   these functions are used in tables 77-a-1 to 77-b-2
c   sgppp -->l,m and n even
c   sgpip --> l even, m odd, n even
c   .................................
c*****************************************************************
        FUNCTION sgppp(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3,2t2,2t1) in table 14
c      in E(2t3,2t2,2t1) in table 77-a-1
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgppp=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgp) then
         mksgp=.true.
         do ik=1,it3
           sg1=0.
           km2=2*ik-2
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgp3(ik)=sgn*hs3(km2+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgrp(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+1)
           enddo
           sgppp=sgppp+sgp3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+1)
           enddo
           sgppp=sgppp+sgp3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgpip(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3,2t+1,2t1) in table 14
c     is in E(2t3,2t2+1,2t1) given in table 77-a-2
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgpip=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgp) then
         mksgp=.true.
         do ik=1,it3
           sg1=0.
           km2=2*ik-2
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
comment             hs30=hs3(km2+1)*rms(3,it3k+1)
           sgp3(ik)=sgn*hs3(km2+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgrp(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+1)
           enddo
           sgpip=sgpip+sgp3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+1)
           enddo
           sgpip=sgpip+sgp3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgppi(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3,2t2,2t1+1) in table 14
c     is in E(2t3,2t2,2t1+1) in table 77-a-2
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgppi=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgp) then
         mksgp=.true.
         do ik=1,it3
           sg1=0.
           km2=2*ik-2
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgp3(ik)=sgn*hs3(km2+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgrp(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+2)
           enddo
           sgppi=sgppi+sgp3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+2)
           enddo
           sgppi=sgppi+sgp3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgpii(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3,2t2+1,2t1+1) in table 14
c     correspond to E(2t3,2t2+1,2t1+1) in table 77-a-1
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgpii=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgp) then
         mksgp=.true.
         do ik=1,it3
           sg1=0.
           km2=2*ik-2
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgp3(ik)=sgn*hs3(km2+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgrp(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+2)
           enddo
           sgpii=sgpii+sgp3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgrp(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+2)
           enddo
           sgpii=sgpii+sgp3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgipp(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3+1,2t2,2t1) in table 14
c     E(2t3+1,2t2,2t1) is given in table 77-b-2
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgipp=0.
       if(it3.le.0) return
       sgn=-1.
       if(.not.mksgi) then
         mksgi=.true.
         do ik=1,it3
           sg1=0.
           km1=2*ik-1
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgi3(ik)=sgn*hs3(km1+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgri(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+1)
           enddo
           sgipp=sgipp+sgi3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+1)
           enddo
           sgipp=sgipp+sgi3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgiip(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3+1,2t2+1,2t1) in table 14
c     inside E(2t3+1,2t2+1,2t1) given in table 77-b-1
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgiip=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgi) then
         mksgi=.true.
         do ik=1,it3
           sg1=0.
           km1=2*ik-1
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgi3(ik)=sgn*hs3(km1+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgri(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+1)
           enddo
           sgiip=sgiip+sgi3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+1)
           enddo
           sgiip=sgiip+sgi3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgipi(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3+1,2t2,2t1+1)  in table 14
c    in E(2t3+1,2t2,2t1+1) in table 77-b-1
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgipi=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgi) then
         mksgi=.true.
         do ik=1,it3
           sg1=0.
           km1=2*ik-1
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgi3(ik)=sgn*hs3(km1+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgri(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+2)
           enddo
           sgipi=sgipi+sgi3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+1)*hs1(jt1j+2)
           enddo
           sgipi=sgipi+sgi3(ik)*sg1
         enddo
         return
       endif
       end
       FUNCTION sgiii(it1,it2,it3)
c  ............................................................
c   FUNCTION sigma(2t3+1,2t2+1,2t1+1) in table 14
c     in E(2t3+1,2t2+1,2t1+1) in table 77-b-2
c  ............................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/herfun/hs1(60),hs2(60),hs3(60)
       common/cars3/exs3,exs2,exs1,as3,sgns3,s32,pwas3,pw3as3
       common/sgsk/sgp3(30),sgi3(30),sgrp(30,30),sgri(30,30)
       common/sgpth/mksgi,mksgp
       logical mksgi,mksgp
       sgiii=0.
       if(it3.eq.0) return
       sgn=-1.
       if(.not.mksgi) then
         mksgi=.true.
         do ik=1,it3
           sg1=0.
           km1=2*ik-1
           it3kp1=it3-ik+1
           it3k=2*(it3-ik)
           sgi3(ik)=sgn*hs3(km1+1)*rms(3,it3k+1)
           do jk=1,it3kp1
             jkm1=jk-1
             jt3kj=2*(it3-ik-jkm1)
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             stoc=rms(1,2*jkm1+1)*rms(2,jt3kj+1)
             sgri(ik,jk)=fpar(it3-ik,jkm1)/stoc
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+2)
           enddo
           sgiii=sgiii+sgi3(ik)*sg1
           sgn=-sgn
         enddo
         return
       else
         do ik=1,it3
           sg1=0.
           it3kp1=it3-ik+1
           do jk=1,it3kp1
             jkm1=jk-1
             jt23kj=2*(it3+it2-ik-jkm1)
             jt1j=2*(it1+jkm1)
             sg1=sg1+sgri(ik,jk)*hs2(jt23kj+2)*hs1(jt1j+2)
           enddo
           sgiii=sgiii+sgi3(ik)*sg1
         enddo
         return
       endif
       end
c*****************************************************************************
c   Functions U and V given in tables 75 and 76
c   uppp ---> l,m,n  even
c   upip ---> l even, m odd, n even
c   .................................
c*****************************************************************************
       FUNCTION uppp(it1,it2,it3)
c  ..............................................................................
c      the integral U(2t3,2t2,2t1,thet) is given in table 75
c     the summation with cos(thet)**2*t1,cos(thet)**2*t2 is found in table 77-a-1
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       uppp=0.
       ind=2*(it1+it2+it3)+1
       idt1=2*it1+1
       idt2=2*it2+1
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arcc(i,j)=ragp(i,ind)*akpcc(i,j)
           u=u+arcc(i,j)*epsi1(i,j)
         enddo
         uppp=uppp+cs(j)*u
       enddo
       return
       end
       FUNCTION upip(it1,it2,it3)
c  ..............................................................................
c   the integral U(2t3,2t2+1,2t1,thet) is given in table 75
c   the summation with cos(thet)**2t1,cos(thet)**(2t2+1) is in table 77-a-2
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       upip=0.
       ind=2*(it1+it2+it3)+2
       idt1=2*it1+1
       idt2=2*it2+2
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arcs(i,j)=ragp(i,ind)*akpcs(i,j)
           u=u+arcs(i,j)*epsi1(i,j)
         enddo
         upip=upip+cs(j)*u
       enddo
       return
       end
       FUNCTION uppi(it1,it2,it3)
c  ..............................................................................
c   the integral U(2t3,2t2,2t1+1,thet)  in table 75
c   the summation with cos(thet)**(2*t1+1),cos(thet)**2*t2  in table 77-a-2
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       uppi=0.
       ind=2*(it1+it2+it3)+2
       idt1=2*it1+2
       idt2=2*it2+1
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arsc(i,j)=ragp(i,ind)*akpsc(i,j)
           u=u+arsc(i,j)*epsi1(i,j)
         enddo
         uppi=uppi+cs(j)*u
       enddo
       return
       end
       FUNCTION upii(it1,it2,it3)
c  ..............................................................................
c   the integral U(2*t3,2*t2+1,2*t1+1,thet)  table 75
c   the summation with cos(thet)**(2*t1+1),cos(thet)**(2*t2+1)  table 77-a-1
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       upii=0.
       ind=2*(it1+it2+it3+1)+1
       idt1=2*it1+2
       idt2=2*it2+2
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arss(i,j)=ragp(i,ind)*akpss(i,j)
           u=u+arss(i,j)*epsi1(i,j)
         enddo
         upii=upii+cs(j)*u
       enddo
       return
       end
       FUNCTION uipp(it1,it2,it3)
c  ..............................................................................
c   the integral U(2*t3+1,2*t2,2*t1,thet)  table 75
c   for the summation with cos(thet)**2*t1,cos(thet)**2*t2 see table 77-b-2
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       uipp=0.
       ind=2*(it1+it2+it3+1)
       idt1=2*it1+1
       idt2=2*it2+1
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arcc(i,j)=ragp(i,ind)*akpcc(i,j)
           u=u+arcc(i,j)*epsi1(i,j)
         enddo
         uipp=uipp+cs(j)*u
       enddo
       return
       end
       FUNCTION uiip(it1,it2,it3)
c  ..............................................................................
c   the integral U(2*t3+1,2*t2+1,2*t1,thet) table 75
c   summation over cos(thet)**2*t1,cos(thet)**(2*t2+1) table 77-b-1
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       uiip=0.
       ind=2*(it1+it2+it3+1)+1
       idt1=2*it1+1
       idt2=2*it2+2
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arcs(i,j)=ragp(i,ind)*akpcs(i,j)
           u=u+arcs(i,j)*epsi1(i,j)
         enddo
         uiip=uiip+cs(j)*u
       enddo
       return
       end
       FUNCTION uipi(it1,it2,it3)
c  ..............................................................................
c   the integral U(2t3+1,2t2,2t1+1,thet)  table 75
c   summation over cos(thet)**(2*t1+1),cos(thet)**2*t2  table 77-b-1
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *             arss(40,40),cs(40)
       uipi=0.
       ind=2*(it1+it2+it3+1)+1
       idt1=2*it1+2
       idt2=2*it2+1
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arsc(i,j)=ragp(i,ind)*akpsc(i,j)
           u=u+arsc(i,j)*epsi1(i,j)
         enddo
         uipi=uipi+cs(j)*u
       enddo
       return
       end
       FUNCTION uiii(it1,it2,it3)
c  .............................................................................
c   the integral U(2*t3+1,2*t2+1,2*t1+1,thet)  table 75
c   summation with cos(thet)**(2*t1+1),cos(thet)**(2*t2+1) (see table 77-b-2)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       uiii=0.
       ind=2*(it1+it2+it3+1)+2
       idt1=2*it1+2
       idt2=2*it2+2
       do j=1,igaus
         cs(j)=co(j,idt1)*sn(j,idt2)*wg(j)
         u=0.
         do i=1,igaus
           arss(i,j)=ragp(i,ind)*akpss(i,j)
           u=u+arss(i,j)*epsi1(i,j)
         enddo
         uiii=uiii+cs(j)*u
       enddo
       return
       end
       FUNCTION vppp(k)
c  ..............................................................................
c    the integral V(2t3,2t2,2t1,k,thet)  table 76
c    sum with cos(thet)**2*t1,cos(thet)**2*t2 (see table 77-a-1)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *              akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *               arss(40,40),cs(40)
       vppp=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arcc(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vppp=vppp+cs(j)*u
       enddo
       return
       end
       FUNCTION vpip(k)
c  ..............................................................................
c   the integral V(2t3,2t2+1,2t1,k,thet)  table 76
c   sum with cos(thet)**2*t1,cos(thet)**(2*t2+1) (table 77-a-2)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *              akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       vpip=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arcs(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vpip=vpip+cs(j)*u
       enddo
       return
       end
       FUNCTION vppi(k)
c  ..............................................................................
c  the integral V(2t3,2t2,2t1+1,k,thet)  table 76
c  sum wih cos(thet)**2*t1,cos(thet)**(2*t2+1) (table 77-a-2)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *              akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       vppi=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arsc(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vppi=vppi+cs(j)*u
       enddo
       return
       end
       FUNCTION vpii(k)
c  ..............................................................................
c   the integral V(2*t3,2*t2+1,2*t1+1,k,thet)   table 76
c   sum with cos(thet)**(2*t1+1),cos(thet)**(2*t2+1) (table 77-a-1)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       vpii=0.
       do j=1,igaus
       u=0.
         do i=1,igaus
           u=u+arss(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vpii=vpii+cs(j)*u
       enddo
       return
       end
       FUNCTION vipp(k)
c  ..............................................................................
c   the integral V(2*t3+1,2*t2,2*t1,thet)  table 76
c   sum with cos(thet)**2*t1,cos(thet)**2*t2 ( table 77-b-2)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *              akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       vipp=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arcc(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vipp=vipp+cs(j)*u
       enddo
       return
       end
       FUNCTION viip(k)
c  ..............................................................................
c   the integral V(2*t3+1,2*t2+1,2*t1,thet)   table 76
c   sum with cos(thet)**2*t1,cos(thet)**(2*t2+1) (table 77-b-1)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       viip=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arcs(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         viip=viip+cs(j)*u
       enddo
       return
       end
       FUNCTION vipi(k)
c  ..............................................................................
c  the integral V(2*t3+1,2*t2,2*t1+1,k,thet)  table 76
c  sum with cos(thet)**(2*t1+1),cos(thet)**2*t2 (table 77-b-1)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       vipi=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arsc(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         vipi=vipi+cs(j)*u
       enddo
       return
       end
       FUNCTION viii(k)
c  ..............................................................................
c   the integral V(2*t3+1,2*t2+1,2*t1+1,thet) in table 76
c   sum with cos(thet)**(2*t1+1),cos(thet)**(2*t2+1) (table 77-b-2)
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       common/rms/rms(3,50),s1,s2,s3
       common/gauss1/absg(40),wg(40),igaus
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/uvtab/epsi1(40,40),epsi2(40,40),akpcc(40,40),
     *             akpcs(40,40),akpsc(40,40),akpss(40,40)
       common/expmod/ragp(40,100),ragm1(40,40)
       common/uvint/arcc(40,40),arcs(40,40),arsc(40,40),
     *              arss(40,40),cs(40)
       viii=0.
       do j=1,igaus
         u=0.
         do i=1,igaus
           u=u+arss(i,j)*ragm1(i,k)*epsi2(i,j)
         enddo
         viii=viii+cs(j)*u
       enddo
       return
       end
       SUBROUTINE pintfast
       implicit real*8 (a-h,o-z)
c  ..............................................................................
C      Particles too far from the C. of G. of the bunch are eliminated
c      for the Almn computations with HERSC
c  ..............................................................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/vl,pi,xmat,rpel,qst
       common/dyn/tref,vref
       common/sc3/beamc,scdist,sce10,cplm,ectt,apl,ichaes,iscsp
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/part/xc(iptsz),yc(iptsz),zc(iptsz)
       common/dimens/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       common/nume/nuelm
       common/vpintim/gcg,bcg
       common/coef/a(30,30,30),xrmsz,yrmsz,zrmsz
       common/macro/ratei
       nuelm=nuelm+1
       write(16,*)' space charge with HERSC call number: ',nuelm
       trmoy=0.
       wcg=0.
       xcg=0.
       xcg=0.
       ycg=0.
       do i=1,ngood
         trmoy=trmoy+f(6,i)
         wcg=wcg+f(7,i)
         xcg=xcg+f(2,i)
         ycg=ycg+f(4,i)
       enddo
       trmoy=trmoy/float(ngood)
       wcg=wcg/float(ngood)
       gcg=wcg/xmat
       bcg=sqrt(1.-1./(gcg*gcg))
       xcg=xcg/float(ngood)
       ycg=ycg/float(ngood)
       zcg=trmoy*fh
c  Isochronism correction in bending magnet
c  apl is the angle due to the lack of isochronicity in the plane (x,z)
c  (see : A modified space charge routine for high intensity bunched
c         beam,P.Lapostolle and 5 authors,NIM A Vol.379,pp21-40,
c         September 1996)
c    does not work with SCHEFF method (iscsp=3)
       apl=0.
       if(iscsp.le.2) then
         xb2x=0.
         xb2z=0.
         xbxz=0.
         imaxx=0
         do np=1,ngood
           gnp=f(7,np)/xmat
           bnp=sqrt(1.-1./(gnp*gnp))
           zc(np)=(trmoy-f(6,np))*bnp*vl/100.
c   *   correct. relativiste
comment           zc(np)=zc(np)*gcg
c   *
           xc(np)=(f(2,np)-xcg)/100.
           xb2z=xb2z+zc(np)*zc(np)
           xb2x=xb2x+xc(np)*xc(np)
           xbxz=xbxz+zc(np)*xc(np)
           imaxx=imaxx+1
         enddo
         xb2z=xb2z/float(imaxx)
         xb2x=xb2x/float(imaxx)
         xbxz=xbxz/float(imaxx)
         apl=atan(-2.*xbxz/(xb2x-xb2z))/2.
       endif
       write(16,*)'*slope of the bunch in plane(Oz,Ox):',apl,' radian'
c bunch at the space charge position
       xbar=0.
       ybar=0.
       zbar=0.
       imaxx=0
c  Divide by 100. to convert from centimeters to meters
       do np=1,ngood
         gnp=f(7,np)/xmat
         bnp=sqrt(1.-1./(gnp*gnp))
         znp=(trmoy-f(6,np))*bnp*vl
c   *   correct. relativiste valero
comment         znp=znp*gcg
c   *
         xnp=f(2,np)
         zc(np)=znp*cos(apl)+xnp*sin(apl)
         xnp=xnp*cos(apl)-znp*sin(apl)
C  convert from mrad to rad
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
C  convert from cm   to m
         xc(np)=(xnp+zc(np)*f3)/100.
         yc(np)=(f(4,np)+zc(np)*f5)/100.
         zc(np)=zc(np)/100.
c evaluate xbar , ybar , zbar
         xbar=xbar+xc(np)
         ybar=ybar+yc(np)
         zbar=zbar+zc(np)
       enddo
       eng=float(ngood)
       xbar=xbar/eng
       ybar=ybar/eng
       zbar=zbar/eng
       do np=1,ngood
         xc(np)=xc(np)-xbar
         yc(np)=yc(np)-ybar
         zc(np)=zc(np)-zbar
       enddo
       return
       end
       SUBROUTINE hcoef
c  ....................................................................
c    the significant Hermite coefficients (called by HERSC)
c  ....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/vl,pi,xmat,rpel,qst
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/hcgrms/xcdg,ycdg,zcdg,ect,eps
       common/ind/lmax,mmax,nmax
       common/indttal/lmnt
       common/randu/ck(15),kmax
       common/coef/a(30,30,30),xrmsz,yrmsz,zrmsz
       common/part/xc(iptsz),yc(iptsz),zc(iptsz)
       common/isoch/apl
       common/dimens/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       common/nume/nuelm
       common/macro/ratei
       common/facrms/fxrms,fyrms,fzrms
       common/rcoef/rdcf
       common/factor/fpir(40,40),fect(30)
       dimension hx(30),hy(30),hz(30),tran(1)
c calculation of Hermite coefficients on the principal axis in x,y,z
       do i=1,ngood
         zcp(i)=zc(i)
         xcp(i)=xc(i)
         ycp(i)=yc(i)
       enddo
       do kn=1,nmax
         do km=1,mmax
           do kl=1,lmax
             a(kl,km,kn)=0.
           enddo
         enddo
       enddo
       zcg=0.
       xcg=0.
       ycg=0.
       do i=1,ngood
         zcg=zcg+zc(i)
         xcg=xcg+xc(i)
         ycg=ycg+yc(i)
       enddo
       xcg=xcg/float(ngood)
       ycg=ycg/float(ngood)
       zcg=zcg/float(ngood)
       xsqsum=0.
       ysqsum=0.
       zsqsum=0.
c    evaluate the rms sizes
       do j=1,ngood
         xcj=xc(j)-xcg
         ycj=yc(j)-ycg
         zcj=zc(j)-zcg
         xsqsum=xsqsum+xcj*xcj
         ysqsum=ysqsum+ycj*ycj
         zsqsum=zsqsum+zcj*zcj
       enddo
       xrmsz=xsqsum/float(ngood)
       yrmsz=ysqsum/float(ngood)
       zrmsz=zsqsum/float(ngood)
       xrmsz=sqrt(xrmsz)
       yrmsz=sqrt(yrmsz)
       zrmsz=sqrt(zrmsz)
       write(16,*) '*RMS of the bunch (m): ',xrmsz,yrmsz,zrmsz
comment       write(17,25) nuelm,xrmsz,yrmsz,zrmsz
comment25     format(2x,i5,3(2x,e12.5))
c     ***TEST control s.c.
comment        rewind(19)
comment        irep=0
comment        rewind(21)
comment        write(21,*) '      x(cm)         y(cm)         z(cm)'
c   coordinates of the bunch  x, y, z
comment        do j=1,ngood
comment          xcoup=(xcp(j)-xcg)/xrmsz
comment          ycoup=(ycp(j)-ycg)/yrmsz
comment          zcoup=(zcp(j)-zcg)/zrmsz
comment          write(21,1188) xcoup,ycoup,zcoup
comment          if((abs(xcoup).lt.fxrms).and.(abs(ycoup).lt.fyrms).and.
comment     *      (abs(zcoup).lt.fzrms)) then
comment           write(19,1188) xcoup,ycoup,zcoup
comment           irep=irep+1
comment          endif
comment        enddo
comment1188     format(3(2x,e12.5))
C    ***END TEST control s.c.
       if(rdcf.ge.1.) then
         hz(1)=1.
         hy(1)=1.
         hx(1)=1.
         irct=0
         do j=1,ngood
           xc(j)=(xcp(j)-xcg)/xrmsz
           yc(j)=(ycp(j)-ycg)/yrmsz
           zc(j)=(zcp(j)-zcg)/zrmsz
           if((abs(xc(j)).lt.fxrms).and.(abs(yc(j)).lt.fyrms).and.
     *        (abs(zc(j)).lt.fzrms)) then
             irct=irct+1
             ach=abs(f(9,j))
             hz(2)=zc(j)
             hy(2)=yc(j)
             hx(2)=xc(j)
             do kn=1,nmax
               if(kn.gt.2) hz(kn)=zc(j)*hz(kn-1)-float(kn-2)*hz(kn-2)
               do km=1,mmax
                 if(km.gt.2) hy(km)=yc(j)*hy(km-1)-float(km-2)*hy(km-2)
                 do kl=1,lmax
                   if(kl.gt.2) hx(kl)=xc(j)*hx(kl-1)-float(kl-2)*
     *                                hx(kl-2)
                   xherm=hx(kl)*hy(km)*hz(kn)/(fect(kl)*fect(km)*
     *                   fect(kn))
                   a(kl,km,kn)=a(kl,km,kn)+xherm/pwtpi*ach
                 enddo
               enddo
             enddo
           endif
         enddo
         rate=float(ngood)/float(irct)
c         write(16,*) ' particles kept in Almn computation: ',irct
         do kn=1,nmax
           do km=1,mmax
             do kl=1,lmax
               a(kl,km,kn)=a(kl,km,kn)*rate
             enddo
           enddo
         enddo
       endif
c*********************************
       if(rdcf.lt.1.) then
         rdcfc=rdcf*float(imax)/float(ngood)
         if(rdcfc.gt.1.) rdcfc=1.
         len=1
         irct=0
         hz(1)=1.
         hy(1)=1.
         hx(1)=1.
         do j=1,ngood
           call rlux(tran,len)
           if(tran(1).le.rdcfc) then
             xc(j)=(xcp(j)-xcg)/xrmsz
             yc(j)=(ycp(j)-ycg)/yrmsz
             zc(j)=(zcp(j)-zcg)/zrmsz
             if((abs(xc(j)).lt.fxrms).and.(abs(yc(j)).le.fyrms).and.
     *         (abs(zc(j)).lt.fzrms)) then
               ach=abs(f(9,j))
               irct=irct+1
               hz(2)=zc(j)
               hy(2)=yc(j)
               hx(2)=xc(j)
               do kn=1,nmax
                 if(kn.gt.2) hz(kn)=zc(j)*hz(kn-1)-float(kn-2)*hz(kn-2)
                 do km=1,mmax
                   if(km.gt.2) hy(km)=yc(j)*hy(km-1)-float(km-2)*
     *                                hy(km-2)
                   do kl=1,lmax
                     if(kl.gt.2) hx(kl)=xc(j)*hx(kl-1)-float(kl-2)*
     *                                  hx(kl-2)
                     xherm=hx(kl)*hy(km)*hz(kn)/(fect(kl)*fect(km)*
     *                     fect(kn))
                     a(kl,km,kn)=a(kl,km,kn)+xherm/pwtpi*ach
                   enddo
                 enddo
               enddo
             endif
           endif
         enddo
         rate=float(ngood)/float(irct)
         write(16,*) ' particles kept in Almn: ',irct
         do kn=1,nmax
           do km=1,mmax
             do kl=1,lmax
                a(kl,km,kn)=a(kl,km,kn)*rate
             enddo
           enddo
         enddo
       endif
c118    format(3(2x,e12.5))
c   *  cesaro  transformation
       lsup=lmax
       msup=mmax
       nsup=nmax
         do kn=1,nsup
           do km=1,msup
             do kl=1,lsup
               cesl=(1.-float(kl-1)/float(lsup))
               cesm=(1.-float(km-1)/float(msup))
               cesn=(1.-float(kn-1)/float(nsup))
               ces=cesl*cesm*cesn
cesaro force             a(kl,km,kn)=a(kl,km,kn)*ces*ces
               a(kl,km,kn)=a(kl,km,kn)*ces
             enddo
           enddo
         enddo
c    *
c   select the significant coefficients
c999    continue
       fond=abs(a(1,1,1))
       itot=0
       iret=0
       do kn=1,nmax
         n=kn-1
         ipar=n-2*int(n/2)
         if(ipar.eq.0)zz=0.
         if(ipar.ne.0) then
         if(n.eq.1) zz=1.
           if(n.eq.3) zz=.75
           if(n.eq.5) zz=.625
         if(n.gt.5.and.n.le.11) zz=0.50
           if(n.gt.11) zz=0.375
         endif
         do km=1,mmax
           m=km-1
         ipar=m-2*int(m/2)
         if(ipar.eq.0)yy=0.
           if(ipar.ne.0) then
           if(m.eq.1) yy=1.
             if(m.eq.3) yy=.75
             if(m.eq.5) yy=.625
           if(m.gt.5.and.m.le.11) yy=0.50
             if(m.gt.11) yy=0.375
           endif
           do kl=1,lmax
             l=kl-1
             ipar=l-2*int(l/2)
             if(ipar.eq.0)xx=0.
             if(ipar.ne.0) then
             if(l.eq.1) xx=1.
             if(l.eq.3) xx=.75
             if(l.eq.5) xx=.625
             if(l.gt.5.and.l.le.11) xx=0.50
             if(l.gt.11) xx=0.375
             endif
             itot=itot+1
             xherm=hermint(xx,l)*hermint(yy,m)*hermint(zz,n)
             ab=abs(a(kl,km,kn)*xherm)/fond
             if(ab.ge.eps) then
             iret=iret+1
           else
             a(kl,km,kn)=0.
           endif
           enddo
         enddo
       enddo
       rpeps=float(iret)/float(itot)
       write(16,*)'*significant terms in Hermite series expansion: ',
     *             iret,' total of terms :',itot
       if(rpeps.ge..3) then
         write(16,*) ' problem in space charge : rpeps gt .3 ',rpeps
         stop
       endif
c   select the maximum values l, m and n for the significant coefficients
       lsup=0
       msup=0
       nsup=0
       lmnt=0
       do kn=1,nmax
         do km=1,mmax
           do kl=1,lmax
             if(a(kl,km,kn).ne.0.) then
               itm=kl+km+kn-3
               if(itm.ge.lmnt) lmnt=itm
               if(lsup.le.kl) lsup=kl
               if(msup.le.km) msup=km
               if(nsup.le.kn) nsup=kn
             endif
            enddo
         enddo
       enddo
       lmnt=lmnt+kmax+4+3
       write(16,*) ' maximum of n m l for the significants terms ',
     *            nsup-1,msup-1,lsup-1
       write(16,*) ' maximun of (t) for the significants terms ',lmnt
c   *  cesaro  transformation
ccc         do kn=1,nsup
ccc           do km=1,msup
ccc             do kl=1,lsup
ccc               cesl=(1.-float(kl-1)/float(lsup))
ccc               cesm=(1.-float(km-1)/float(msup))
ccc               cesn=(1.-float(kn-1)/float(nsup))
ccc               ces=cesl*cesm*cesn
ccc              a(kl,km,kn)=a(kl,km,kn)*ces
ccc             enddo
ccc           enddo
ccc         enddo
c    *
       return
       end
       SUBROUTINE hersc(ini)
c  ........................................................................
c       space charge method: HERSC
c  ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/vl,pi,xmat,rpel,qst
       common/sc3/beamc,scdist,sce10,cplm,ectt,apl,ichaes,iscsp
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/cdek/dwp(iptsz)
       common/part/xc(iptsz),yc(iptsz),zc(iptsz)
       common/coef/a(30,30,30),xrmsz,yrmsz,zrmsz
       common/hcgrms/xcdg,ycdg,zcdg,ect,eps
       common/ind/lmax,mmax,nmax
       common/indin/lmaxi,mmaxi,nmaxxi
       common/indttal/lmnt
       common/rms/rms(3,50),s1,s2,s3
       common/randu/ck(15),kmax
       common/circu/co(40,50),sn(40,50),blam(40,100)
       common/gauss1/absg(40),wg(40),igaus
       common/field/ex,ey,ez
       common/expmod/ragp(40,100),ragm1(40,40)
       common/const/pi2,sqpi,pwtpi,sqpi2,sq2pi
       common/facrms/fxrms,fyrms,fzrms
       common/dimens/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       common/beamsa/fs(7,iptsz)
       common/dcspa/iesp
       common/cmpte/iell
       common/npart/imaxr
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rcoef/rdcf
       common/macro/ratei
       common/tapes/in,ifile,meta
       common/posc/xpsc
       common/locher/rx,ry,rz,nx,ny,nz
       logical ichaes,iesp
       dimension exk(15,15,30),eyk(15,15,30),ezk(15,15,30)
       data lbmax,mbmax,nbmax/23,23,23/
c ini = 0 :initialisation of the routine
       if(ini.le.0) then
         sq2pi=sqrt(2.*pi)
c   input parameters
         igaus=9
c   eps allows to select the significant terms in the HERMITE series expansion
         eps=8.e-03
c   maximun of the upper limits in the HERMITE series expansion
         lmaxi=11
         mmaxi=11
         nmaxxi=11
c   The HERMITE coefficients Almn are computed with the particles include in a cube
c   with regard to the RMS of the bunch. The sizes of the cube are:RMSx*fxrms, ........
         fxrms=2.5
         fyrms=2.5
         fzrms=2.5
         if(ini.lt.0) then
           read(in,*)lmaxi,mmaxi,nmaxxi
           write(16,*)'upper limits ',lmaxi,mmaxi,nmaxxi
           read(in,*)fxrms,fyrms,fzrms
           write(16,*) ' rms factors',fxrms,fyrms,fzrms
           read(in,*)eps
           write(16,*) ' select the Hermite cefficients with: ',eps
         endif
c    MESH**********
c      limits in (+-) RMS
         rx=4.
         ry=4.
         rz=4.
c      number of nodes in each direction
         nx=8
         ny=8
         nz=8
c  rdcf   : in %/100 , allows to select a reduced number of particles for the Almn.
         call shuffle
         call gausse
         read(in,*) rdcf
         if(rdcf.gt.1.) rdcf=1.
         if(rdcf.le.0.) then
           if(ngood.gt.15000) then
             rdcf=15000./float(ngood)
           else
             rdcf=1.
           endif
         endif
         call table(lbmax,mbmax,nbmax)
         return
       endif
c **********************************************************************
       if(ini.eq.1) then
c    .................................................................................
c       seeks: the significant coefficients  in the Hermite series expansions
c              ---> storage array a(l,m,n)
c              maximum values of l,m and n of the significants a(l,m,n) --->lmax,mmax,nmax
c              RMS of the bunch ---> xa,xb and xc
c       order of the succesives integrals
c       isucc=1: xa<xb,xc , rms(3,2)=xa, rms(2,2)=xb,  rms(1,2)=xc
c                           rms(3,k)=xa**(k-1), ..........
c       isucc=2: xb<xa,xc , rms(3,2)=xb, rms(2,2)=xc,  rms(1,2)=xa
c       isucc=1: xc<xa,xb , rms(3,2)=xc, rms(2,2)=xa,  rms(1,2)=xb
c       compute the array blam(i,j)
c        blam(i,j): FUNCTION beta**2 in the fifth relation in table 31
c                   blam(i,j)=beta**2j,j=1,...,im=(lmax+mmax+nmax+kmax)/2+2
c                   i is given from the i-Gauss position in absg(i)
c   ..................................................................................
         call bhdist
         call trms(isucc)
         sta=a(1,1,1)
         icoa=0
         do ican=1,nmax
           do icam=1,mmax
             do ical=1,lmax
               if(a(ican,icam,ical).ne.0.) icoa=icoa+1
             enddo
           enddo
         enddo
       return
       endif
c **********************************************************************
c   beam self-fields
       if(ini.eq.2) then
c   limits of the mesch : rx ry rz
c   length of the steps
         delx=2.*rx/float(nx)
         dely=2.*ry/float(ny)
         delz=2.*rz/float(nz)
c    Initialize constants
c   xmass-Kg,wavel=cm,charge=coul.
c   clight=cm/sec MHz/Hz ,freq=Mhz ,xmat=MeV
c freq. in MHz
         freq=fh*0.5e-06/pi
         epsilon=8.854189586e-12
         const3=1.E-06
         const2=1.E-06/xmat
         dxp=0.
         dyp=0.
         dw=0.
         dz=scdist/100.
         dz1=dz
         zsot1=0.
c   qmpart=macro-charge
         qmpart=1.0e-9*beamc/(float(imax)*freq)
         if(ratei.le.0.) then
           write(16,*) ' all the particles are lost '
           stop
         endif
         qmpart=qmpart*ratei
         vrms=xrmsz*yrmsz*zrmsz
         cmacro=qmpart/(epsilon*vrms)
         wcg=0.
         do i=1,ngood
          wcg=wcg+f(7,i)
         enddo
         wcg=wcg/float(ngood)
         gcg=wcg/xmat
         gmoy=gcg
         bcg=sqrt(1.-1./(gcg*gcg))
         bmoy=bcg
c  relativistic correction
         cmacro=cmacro/gmoy
c
         cmacrxy=cmacro/(bmoy*bmoy*gmoy*gmoy)
c   beam self-fields at the nodes of the mesh
         imail=0
         rcy=-ry
         do j=1,ny
           rcx=-rx
           do i=1,nx
             rcz=-rz
             do k=1,nz
               if(isucc.eq.1) then
                 s3=rcx
                 s2=rcy
                 s1=rcz
                 ax=pwtpi*rms(3,2)
                 ay=pwtpi*rms(2,2)
                 az=pwtpi*rms(1,2)
                 rrx=rms(3,2)
                 rry=rms(2,2)
                 rrz=rms(1,2)
               endif
               if(isucc.eq.2) then
                 s3=rcy
                 s2=rcz
                 s1=rcx
                 ax=pwtpi*rms(1,2)
                 ay=pwtpi*rms(3,2)
                 az=pwtpi*rms(2,2)
                 rrx=rms(1,2)
                 rry=rms(3,2)
                 rrz=rms(2,2)
               endif
               if(isucc.eq.3) then
                 s3=rcz
                 s2=rcx
                 s1=rcy
                 ax=pwtpi*rms(2,2)
                 ay=pwtpi*rms(1,2)
                 az=pwtpi*rms(3,2)
                 rrx=rms(2,2)
                 rry=rms(1,2)
                 rrz=rms(3,2)
               endif
             call uvrms
c     fields
c       loop over the l,m and n
             exk(i,j,k)=0.
             eyk(i,j,k)=0.
             ezk(i,j,k)=0.
             do jn=1,nmax
                 jn1=jn-1
                 do jm=1,mmax
                   jm1=jm-1
                   do jl=1,lmax
                     jl1=jl-1
                     if(a(jl,jm,jn).ne.0.) then
                       call fielde(jl1,jm1,jn1,isucc)
c     the beam self-fields are  in tables 67-a to 67-h
                       exk(i,j,k)=a(jl,jm,jn)/ax*ex+exk(i,j,k)
                       eyk(i,j,k)=a(jl,jm,jn)/ay*ey+eyk(i,j,k)
                       ezk(i,j,k)=a(jl,jm,jn)/az*ez+ezk(i,j,k)
                     endif
                   enddo
                 enddo
               enddo
               rcz=rcz+delz
               imail=imail+1
             enddo
             rcx=rcx+delx
           enddo
           rcy=rcy+dely
         enddo
c xi in Amps, ibeam in mA
         dxp=0.
         dyp=0.
         dw=0.
c   save the  particles coordinates
         do i=1,ngood
           xc(i)=xcp(i)
           yc(i)=ycp(i)
           zc(i)=zcp(i)
         enddo
c Do integration to determine Ex,Ey,Ez for each macro particle
         nprint=1
         insd=0
         iout=0
         ickl=0
         do ic=1,ngood
c     position in the mesh
           ickl=ickl+1
c    *    valero
           u=xc(ic)/xrmsz
           v=yc(ic)/yrmsz
           w=zc(ic)/zrmsz
comment           u=xc(ic)/rrx
comment           v=yc(ic)/rry
comment           w=zc(ic)/rrz
           i=int((u+rx)/delx)+1
           j=int((v+ry)/dely)+1
           k=int((w+rz)/delz)+1
c   the particle is in the mesh
           if(i.gt.0.and.i.le.nx.and.j.gt.0.and.j.le.ny.
     *                        and.k.gt.0.and.k.le.nz) then
             xnd1=-rx+float(i-1)*delx
             ynd1=-ry+float(j-1)*dely
             znd1=-rz+float(k-1)*delz
             delu=u-xnd1
             delv=v-ynd1
             delw=w-znd1
             delux=delu/delx
             delvy=delv/dely
             delwz=delw/delz
c *******interpollation for particle inside a cube
c %%%%%%     in plane 1:
c node 1: (i,j,k)  node 2:(i+1,j,k) node 3:(i+1,j,k+1) node 4:(i,j,k+1)
c  axis (nd1,nd2)
             ex12=(exk(i+1,j,k)-exk(i,j,k))*delux+exk(i,j,k)
c  axis (nd4,nd3)
             ex43=(exk(i+1,j,k+1)-exk(i,j,k+1))*delux+exk(i,j,k+1)
c    plane 1
             exp1=(ex43-ex12)*delwz+ex12
c %%%%%%%     in plane 3:
c  node 5:(i,j+1,k) node 6:(i+1,j+1,k) node 7:(i+1,j+1,k+1) node 8:(i,j+1,k+1)
c  axis(nd5,nd6)
             ex56=(exk(i+1,j+1,k)-exk(i,j+1,k))*delux+exk(i,j+1,k)
c  axis(nd8,nd7)
             ex87=(exk(i+1,j+1,k+1)-exk(i,j+1,k+1))*delux+exk(i,j+1,
     *             k+1)
c   plane 3
             exp3=(ex87-ex56)*delwz+ex56
c  $$$$$ plane1+plane3
             exp13=(exp3-exp1)*delvy+exp1
c %%%%%    in plane 2:
c node 1: (i,j,k)  node 2:(i+1,j,k) node 5:(i,j+1,k) node 6:(i+1,j+1,k)
c axis (nd1,nd2) : ex12
c axis (nd5,nd6) : ex56
c  plane 2
             exp2=(ex56-ex12)*delvy+ex12
c  %%%%%%    in plane 5
c  node 4:(i,j,k+1) node 3:(i+1,j,k+1) node 8:(i,j+1,k+1) node 7:(i+1,j+1,k+1)
c   axis(nd4,nd3)
             ex43=(exk(i+1,j,k+1)-exk(i,j,k+1))*delux+exk(i,j,k+1)
c   axis(nd8,nd7)
             ex87=(exk(i+1,j+1,k+1)-exk(i,j+1,k+1))*delux+exk(i,j+1,
     *             k+1)
c  plane 5
             exp5=(ex87-ex43)*delvy+ex43
c  $$$$ plane2+plane5
             exp25=(exp5-exp2)*delwz+exp2
c   &&&& (plane1+plane3)+(plane2+plane5) --> field Ex(u,v,w)
             ext=(exp13+exp25)/2.
c ******* field component Ey
c %%%%%%     in plane 2:
c node 1: (i,j,k)  node 2:(i+1,j,k) node 5:(i,j+1,k) node 6:(i+1,j+1,k)
c  axis(nd1,nd5)
             ey15=(eyk(i,j+1,k)-eyk(i,j,k))*delvy+eyk(i,j,k)
c  axis (nd2,nd6)
             ey26=(eyk(i+1,j+1,k)-eyk(i+1,j,k))*delvy+eyk(i+1,j,k)
c   plane 2
             eyp2=(ey26-ey15)*delu/delx+ey15
c %%%%%%     in plane 5:
c  node 4:(i,j,k+1) node 3:(i+1,j,k+1) node 8:(i,j+1,k+1) node 7:(i+1,j+1,k+1)
c  axis(nd4,nd8)
             ey48=(eyk(i,j+1,k+1)-eyk(i,j,k+1))*delvy+eyk(i,j,k+1)
c  axis(nd3,nd7)
             ey37=(eyk(i+1,j+1,k+1)-eyk(i+1,j,k+1))*delvy+eyk(i+1,j,
     *             k+1)
c   plane 5
             eyp5=(ey37-ey48)*delux+ey48
c $$$$$$ plane2+plane5
             eyp25=(eyp5-eyp2)*delwz+eyp2
c %%%%%%     in plane 6:
c  node 1:(i,j,k) node 5:(i,j+1,k) node 4:(i,j,k+1) node 8:(i,j+1,k+1)
c  axis(nd1,nd5): ey15
c  axis(nd4,nd8): ey48
c  plane 6
             eyp6=(ey48-ey15)*delwz+ey15
c %%%%%%     in plane 4:
c  node 2:(i+1,j,k) node 3:(i+1,j,k+1) node 6:(i+1,j+1,k) node 7:(i+1,j+1,k+1)
c   axis(nd2,nd6): ey26
c   axis(nd3,nd7): ey37
c    plane 4
             eyp4=(ey37-ey26)*delwz+ey26
c $$$$plane6+plane4
             eyp46=(eyp4-eyp6)*delux+eyp6
c  £££ (plane2+plane5)+(plane4+plane6)--> field Ez(u,v,w)
             eyt=(eyp25+eyp46)/2.
c *******interpollation of the field component Ez
c %%%%%%     in plane 1:
c  node 1:(i,j,k) node 2:(i+1,j,k) node 3:(i+1,j,k+1) node 4:(i,j,k+1)
c axis(nd1,nd4)
             ez14=(ezk(i,j,k+1)-ezk(i,j,k))*delwz+ezk(i,j,k)
c axis(nd2,nd3)
             ez23=(ezk(i+1,j,k+1)-ezk(i+1,j,k))*delwz+ezk(i+1,j,k)
c   plane 1
             ezp1=(ez23-ez14)*delux+ez14
c %%%%%%     in plane 3:
c  node 5:(i,j+1,k) node 6:(i+1,j+1,k) node 7:(i+1,j+1,k+1) node 8:(i,j+1,k+1)
c  axis (nd5,nd8)
             ez58=(ezk(i,j+1,k+1)-ezk(i,j+1,k))*delwz+ezk(i,j+1,k)
c  axis (nd6,nd7)
             ez67=(ezk(i+1,j+1,k+1)-ezk(i+1,j+1,k))*delwz+ezk(i+1,
     *             j+1,k)
c   plane 3
             ezp3=(ez67-ez58)*delux+ez58
c  $$$$$$ plane1+plane3
             ezp13=(ezp3-ezp1)*delvy+ezp1
c %%%%%%     in plane 6:
c  node 1:(i,j,k) node 5:(i,j+1,k) node 4:(i,j,k+1) node 8:(i,j+1,k+1)
c  axis (nd1,nd4)
             ez14=(ezk(i,j,k+1)-ezk(i,j,k))*delwz+ezk(i,j,k)
c  axis (nd5,nd8)
             ez58=(ezk(i,j+1,k+1)-ezk(i,j+1,k))*delwz+ezk(i,j+1,k)
c   plane 6
             ezp6=(ez58-ez14)*delvy+ez14
c %%%%%%     in plane 4:
c  node 2:(i+1,j,k) node 3:(i+1,j,k+1) node 6:(i+1,j+1,k) node 7:(i+1,j+1,k+1)
c  axis(nd2,nd3)
             ez23=(ezk(i+1,j,k+1)-ezk(i+1,j,k))*delwz+ezk(i+1,j,k)
c  axis(nd6,nd7)
             ez67=(ezk(i+1,j+1,k+1)-ezk(i+1,j+1,k))*delwz+ezk(i+1,
     *             j+1,k)
c     plane 4
             ezp4=(ez67-ez23)*delvy+ez23
c $$$$$$ plane6+plane4
             ezp64=(ezp4-ezp6)*delux+ezp6
c £££££££ (P1+P3)+(P6+P4)----> Ez(u,v,w)
             ezt=(ezp13+ezp64)/2.
             insd=insd+1
c calculate kicks in x',y' and z' (energy)
             eztp=ezt*cos(apl)-ext*sin(apl)
             extp=ezt*sin(apl)+ext*cos(apl)
             ezt=eztp
             ext=extp
             gsc=f(7,ic)/xmat
             bsc=sqrt(1.-1./(gsc*gsc))
comment             cmacrxy=cmacro/(bmoy*bmoy*gmoy*gmoy*gmoy)
comment             cmacrxy=cmacro/(bsc*bsc*gsc*gsc)
             dxp=const2*ext*dz*cmacrxy*abs(f(9,ic))
             dyp=const2*eyt*dz*cmacrxy*abs(f(9,ic))
             dw=const3*ezt*dz*cmacro*abs(f(9,ic))
c ****test************************************
comment              dw=dw/gmoy
c  ****************************************
             if(.not.iesp) then
c     load the entrance beam parameters for cavities or gaps
               do js=1,7
                 f(js,ic)=fs(js,ic)
               enddo
               f(3,ic)=f(3,ic)+dxp*1000.
               f(5,ic)=f(5,ic)+dyp*1000.
               f(2,ic)=f(2,ic)-dz1*dxp*100.*xpsc
               f(4,ic)=f(4,ic)-dz1*dyp*100.*xpsc
comment               f(2,ic)=f(2,ic)-dz1*dxp*100.
comment               f(4,ic)=f(4,ic)-dz1*dyp*100.
               dwp(ic)=dw
             else
               f(3,ic)=f(3,ic)+dxp*1000.
               f(5,ic)=f(5,ic)+dyp*1000.
               f(7,ic)=f(7,ic)+dw
             endif
           else
c  the particle is not in the mesh
c   computed the beam self-fields as made at the nodes
             ax=0.
             ay=0.
             az=0.
             if(isucc.eq.1) then
               s3=xc(ic)/rms(3,2)
               s2=yc(ic)/rms(2,2)
               s1=zc(ic)/rms(1,2)
               ax=pwtpi*rms(3,2)
               ay=pwtpi*rms(2,2)
               az=pwtpi*rms(1,2)
               rrx=rms(3,2)
               rry=rms(2,2)
               rrz=rms(1,2)
             endif
             if(isucc.eq.2) then
               s3=yc(ic)/rms(3,2)
               s2=zc(ic)/rms(2,2)
               s1=xc(ic)/rms(1,2)
               ax=pwtpi*rms(1,2)
               ay=pwtpi*rms(3,2)
               az=pwtpi*rms(2,2)
               rrx=rms(1,2)
               rry=rms(3,2)
               rrz=rms(2,2)
             endif
             if(isucc.eq.3) then
               s3=zc(ic)/rms(3,2)
               s2=xc(ic)/rms(2,2)
               s1=yc(ic)/rms(1,2)
               ax=pwtpi*rms(2,2)
               ay=pwtpi*rms(1,2)
               az=pwtpi*rms(3,2)
               rrx= rms(2,2)
               rry=rms(1,2)
               rrz=rms(3,2)
             endif
c       storage arrays of the functions in the integrals in tables 75 and 76
c       these values given in table 73 are independent of l, m and n
c       They are stored in: epsi1(i,j),epsi2(i,j), akpc1(i,j), akpc2(i,j),
c          akps1(i,j), akps2(i,j) where i,j are the Gauss positions
             call uvrms
c     field computation
c     loop over the l,m and n
             ext=0.
             eyt=0.
             ezt=0.
             do jn=1,nmax
               jn1=jn-1
               do jm=1,mmax
                 jm1=jm-1
                 do jl=1,lmax
                   jl1=jl-1
                   if(a(jl,jm,jn).ne.0.) then
c     in output the SUBROUTINE field returns the values of the FUNCTION E*(l,m,n)
c     shown in tables 77-a-1 to 77-b-2 in x,y and z-directions-->ex,ey,ez
                     call fielde(jl1,jm1,jn1,isucc)
c     the corresonding field components are obtained from tables 67-a to 67-h
                     ext=a(jl,jm,jn)/ax*ex+ext
                     eyt=a(jl,jm,jn)/ay*ey+eyt
                     ezt=a(jl,jm,jn)/az*ez+ezt
                   endif
                 enddo
               enddo
             enddo
             iout=iout+1
c    kicks
c    kicks computation
c calculate kick in x',y' and z' (energy)
c   isochronism correction
             eztp=ezt*cos(apl)-ext*sin(apl)
             extp=ezt*sin(apl)+ext*cos(apl)
             ezt=eztp
             ext=extp
             gsc=f(7,ic)/xmat
             bsc=sqrt(1.-1./(gsc*gsc))
comment             cmacrxy=cmacro/(bmoy*bmoy*gmoy*gmoy)
comment             cmacrxy=cmacro/(bsc*bsc*gsc*gsc)
             dxp=const2*ext*dz*cmacrxy*abs(f(9,ic))
             dyp=const2*eyt*dz*cmacrxy*abs(f(9,ic))
             dw=const3*ezt*dz*cmacro*abs(f(9,ic))
             if(.not.iesp) then
c     load the beam at the input at cavities or gaps
               do js=1,7
                 f(js,ic)=fs(js,ic)
               enddo
               f(3,ic)=f(3,ic)+dxp*1000.
               f(5,ic)=f(5,ic)+dyp*1000.
               f(2,ic)=f(2,ic)-dz1*dxp*100.*xpsc
               f(4,ic)=f(4,ic)-dz1*dyp*100.*xpsc
comment               f(2,ic)=f(2,ic)-dz1*dxp*100.
comment               f(4,ic)=f(4,ic)-dz1*dyp*100.
               dwp(ic)=dw
             else
               f(3,ic)=f(3,ic)+dxp*1000.
               f(5,ic)=f(5,ic)+dyp*1000.
               f(7,ic)=f(7,ic)+dw
             endif
           endif
c end of if particle not in mesh
           nprint=nprint+1
c end of the loop ic
         enddo
         write(16,*) ' particles in the mesh:',insd,' outside: ',iout
c end of if when ini=2
       endif
c1000   format(4(2x,e12.5))
       return
       end
       SUBROUTINE xtypl1(GAMI,SAPHI,QSC,DCG)
c  ...........................................................................
c  called by ETGAP and RESTAY
C INTEGRALS  E(z)*(BG)**-3 *z**n   n=0,1
C INTEGRALS  dE(z)/dT*(BG)**-3 *z**n   n=0,1,2
C INTEGRALS FONCTIONS HA0(Z) et HB0(Z)
C INTEGRALS OF THE FIRST DERIVATIVES OF HA0(Z) et HB0(Z)
c  ...........................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/POSI/IST
       COMMON/JACOB/GAKS,GAPS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/TYPL1/YH1K0,YH1K1,YP1K1,YP1K2,
     X              YH1K00,YH1K01,YP1K01,YP1K02,
     X              YH10,YH11,YP11,YP12
       COMMON/TYPL2/YH2K0,YH2K1,
     X              YP2K1,YP2K2,
     X              YH2K00,YH2K01,YP2K01,YP2K02,
     X              YH20,YH21,YP21,YP22
C      TRANSVERSAL dE(z)/dt
       COMMON/TYPI1/YE1K0,YE1K1,YE1K2,
     X              YE1KC0,YE1KC1,YE1KC2,
     X              YE10,YE11,YE12
       COMMON/TYPI2/YE2K0,YE2K1,YE2K2,
     X              YE2KC0,YE2KC1,YE2KC2,
     X              YE20,YE21,YE22
       DIMENSION H(17),T(17)
       DATA T /-.990575473,-.950675522,-.880239154,
     1          -.781514004,-.657671159,-.512690537,
     2          -.351231763,-.178484181,  0.,
     3           .178484181, .351231763,.512690537,
     4           .657671159, .781514004,.880239154,
     5           .950675522, .990575473/
       DATA H /.024148303,.055459529,.085036148,
     1         .111883847,.135136368,.154045761,
     2         .168004102,.176562705,.179446470,
     3         .176562705,.168004102,.154045761,
     4         .135136368,.111883847,.085036148,
     5         .055459529,.024148303/
       FH0=FH/VL
c  valero 08/08/07
comment       CGI=ABS(QSC/XMAT)
       CGI=QSC/XMAT
c   *
C      Circular cosinus functions
C     In longitudinal direction
       YH1K0=0.
       YH1K1=0.
       YP1K1=0.
       YP1K2=0.
       YH1K00=0.
       YH1K01=0.
       YP1K01=0.
       YP1K02=0.
       YH10=0.
       YH11=0.
       YP11=0.
       YP12=0.
C      In transverse direction ( dE(z)/dt )
       YE1K0=0.
       YE1K1=0.
       YE1K2=0.
       YE1KC0=0.
       YE1KC1=0.
       YE1KC2=0.
       YE10=0.
       YE11=0.
       YE12=0.
C      Circular sinus functions
C     In longitudinal direction
       YH2K0=0.
       YH2K1=0.
       YP2K1=0.
       YP2K2=0.
       YH2K00=0.
       YH2K01=0.
       YP2K01=0.
       YP2K02=0.
       YH20=0.
       YH21=0.
       YP21=0.
       YP22=0.
C      In transverse direction ( dE(z)/dt )
       YE2K0=0.
       YE2K1=0.
       YE2K2=0.
       YE2KC0=0.
       YE2KC1=0.
       YE2KC2=0.
       YE20=0.
       YE21=0.
       YE22=0.
C    Calculates the integrals
       DTILK=EQVL
       GAM2=GAMI**2
       BETI=SQRT(1. - 1./GAM2 )
       XK1=FH0/BETI
       TILTA2=PHSLIP/(2.*EQVL)
       CGAM10=((GAMI*GAMI-1.)**1.5)/FH0
       PHCRTK=(T1K*SK-S1K*TK)/(TK*TK+SK*SK)
       ist=0
       do i=1,17
         ist=ist+1
         XCC= EQVL*(1.+T(I))/2.
         XCC1=XCC+ASDL
         IF(XCC1.GT.DCG) GO TO 200
         PHIT0=SAPHI-PHSLIP*(EQVL-XCC)/(2.*EQVL)+PAVPH
C       FONCTION GAMMA (Z)
         IF(PHSLIP.NE.0.) THEN
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)/PHSLIP
           GIS =  SIN(XCC*TILTA2)
         ELSE
           GIT = CGI * SQCTTF*COS(PHIT0-PCREST)
           GIS =  XCC/(2.*EQVL)
         ENDIF
         GI=GAMI+GIT * GIS
         BI=SQRT(1.-1./(GI*GI))
C       Derivative of  G0(Z) relative to the equivalent k
         PHIT0K=-DTILK*(1.-XCC/EQVL)/2.
         IF(PHSLIP.NE.0.) THEN
           GIC =  COS(XCC*TILTA2)
           GAK1 = DTILK*COS(PHIT0-PCREST)*GIS/(PHSLIP*PHSLIP)
           GAK2 = SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)/PHSLIP
           GAK3 = DTILK*COS(PHIT0-PCREST)*XCC*GIC/(2.*PHSLIP*EQVL)
           GAK = CGI * SQCTTF* (-GAK1-GAK2+GAK3)
         ELSE
           GAK1 = SIN(PHIT0-PCREST)*GIS*(PHIT0K-PHCRTK)
           GAK = -CGI * SQCTTF* GAK1
         ENDIF
         IF(I.EQ.17) GAKS=GAK
C      Calculate the integrals of  HA0(Z) and HB0(Z)
         XINT=1./(BI*BI*BI*GI*GI*GI)
         PHIT1=PHIT0+XCC*PHSLIP/(2.*EQVL)
         PHTZ0=(XCC/EQVL-.5)*DTILK
         HA0=2.*SQCTTF*COS(PHIT1-PCREST)*XINT
         HB0=2.*SQCTTF*SIN(PHIT1-PCREST)*XINT
C      n=0
C       Longitudinal direction
         YH10=YH10+H(I)*HA0
         YH20=YH20+H(I)*HB0
C        Transverse direction
         YE10=YH20
         YE20=YH10
C      n=1
C        Longitudinal direction
         YH11=YH11+H(I)*HA0*XCC1
         YH21=YH21+H(I)*HB0*XCC1
C        Transverse direction
         YE11=YE11+H(I)*HB0*XCC
         YE21=YE21+H(I)*HA0*XCC
C      n=2
         YE12=YE12+H(I)*HB0*XCC*XCC
         YE22=YE22+H(I)*HA0*XCC*XCC
C      Calculate the integrals of the derivative of HA0(Z)
         DHA01=SQCTTF*COS(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
         DHA02=SQCTTF*COS(PHIT1-PCREST)*GI*GAK/((GI*GI-1.)**2.5)
         DHA03=SQCTTF*SIN(PHIT1-PCREST)/((GI*GI-1.)**1.5)
C      n=0
C      Longitudinal direction
         YH1K00=YH1K00+H(I)*DHA01*6.*CGAM10
         YH1K0=YH1K0+H(I)*(-6.*DHA02+
     X         2.*(PHCRTK-PHTZ0)*DHA03)
C      Transverse direction
         YE2KC0=YH1K00
         YE2K0=YH1K0
C      n=1
C      Longitudinal direction
         YH1K01=YH1K01+H(I)*DHA01*XCC1*6.*CGAM10
         YH1K1=YH1K1+H(I)*XCC1*(-6.*DHA02+
     X        2.*(PHCRTK-PHTZ0)*DHA03)
C      Transverse direction
         YE2KC1=YE2KC1+H(I)*DHA01*XCC*6.*CGAM10
         YE2K1=YE2K1+H(I)*XCC*(-6.*DHA02-
     X       2.*(PHCRTK-PHTZ0)*DHA03)
C      n=2
         YE2KC2=YE2KC2+H(I)*DHA01*XCC*XCC*6.*CGAM10
         YE2K2=YE2K2+H(I)*XCC*XCC*(-6.*DHA02-
     X       2.*(PHCRTK-PHTZ0)*DHA03)
C      INTEGRALS Of  HB0(Z)
         DHB01=SQCTTF*SIN(PHIT1-PCREST)*GI/((GI*GI-1.)**2.5)
         DHB02=SQCTTF*SIN(PHIT1-PCREST)*GI*GAK/((GI*GI-1.)**2.5)
         DHB03=SQCTTF*COS(PHIT1-PCREST)/((GI*GI-1.)**1.5)
C      n=0
C      LONGITUDINAL
         YH2K00=YH2K00+H(I)*DHB01*6.*CGAM10
         YH2K0=YH2K0+H(I)*(-6.*DHB02-
     X          2.*(PHCRTK-PHTZ0)*DHB03)
C      TRANSVERSE
         YE1KC0=YH2K00
         YE1K0=YH2K0
C      n=1
C      LONGITUDINAL
         YH2K01=YH2K01+H(I)*DHB01*XCC1*6.*CGAM10
         YH2K1=YH2K1+H(I)*XCC1*(-6.*DHB02-
     X       2.*(PHCRTK-PHTZ0)*DHB03)
C      TRANSVERSE
         YE1KC1=YE1KC1+H(I)*DHB01*XCC*6.*CGAM10
         YE1K1=YE1K1+H(I)*XCC*(-6.*DHB02-
     X       2.*(PHCRTK-PHTZ0)*DHB03)
C      n=2
         YE1KC2=YE1KC2+H(I)*DHB01*XCC*XCC*6.*CGAM10
         YE1K2=YE1K2+H(I)*XCC*XCC*(-6.*DHB02-
     X       2.*(PHCRTK-PHTZ0)*DHB03)
C      Calculate the integrals of PA0(Z) et PB0(Z)
         PA0=2.*SQCTTF*COS(PHIT1-PCREST)*XINT*XINT
         PB0=2.*SQCTTF*SIN(PHIT1-PCREST)*XINT*XINT
C      n=1
         YP11=YP11+H(I)*PA0*XCC1
         YP21=YP21+H(I)*PB0*XCC1
C      n=2
         YP12=YP12+H(I)*PA0*XCC1*XCC1
         YP22=YP22+H(I)*PB0*XCC1*XCC1
C      Calculate the integrals of the derivatives of PA0(Z)
         DPA01=SQCTTF*COS(PHIT1-PCREST)*GI/((GI*GI-1.)**4)
         DPA02=SQCTTF*COS(PHIT1-PCREST)*GI*GAK/((GI*GI-1.)**4)
         DPA03=SQCTTF*SIN(PHIT1-PCREST)/((GI*GI-1.)**3)
C      n=1
         YP1K01=YP1K01+H(I)*DPA01*12.*CGAM10*XCC1
         YP1K1=YP1K1+H(I)*XCC1*(-12.*DPA02+
     X         2.*(PHCRTK-PHTZ0)*DPA03)
C      n=2
         YP1K02=YP1K02+H(I)*DPA01*12.*CGAM10*XCC1*XCC1
         YP1K2=YP1K2+H(I)*XCC1*XCC1*(-12.*DPA02 +
     X         2.*(PHCRTK-PHTZ0)*DPA03)
C      INTEGRALES DERIVES PB0(Z)
         DPB01=SQCTTF*SIN(PHIT1-PCREST)*GI/((GI*GI-1.)**4)
         DPB02=SQCTTF*SIN(PHIT1-PCREST)*GI*GAK/((GI*GI-1.)**4)
         DPB02=DPB02
         DPB03=SQCTTF*COS(PHIT1-PCREST)/((GI*GI-1.)**3)
C      n=1
         YP2K01=YP2K01+H(I)*DPB01*12.*CGAM10*XCC1
         YP2K1=YP2K1+H(I)*XCC1*(-12.*DPB02-
     X         2.*(PHCRTK-PHTZ0)*DPB03)
C      n=2
         YP2K02=YP2K02+H(I)*DPB01*12.*CGAM10*XCC1*XCC1
         YP2K2=YP2K2+H(I)*XCC1*XCC1*(-12.*DPB02 -
     X         2.*(PHCRTK-PHTZ0)*DPB03)
       enddo
200    CONTINUE
C    in COS
C     LONGITUDINAL INTEGRALS
       YH1K00=YH1K00/2. *EQVL
       YH1K01=YH1K01/2. *EQVL
       YH1K0=YH1K0/2. *EQVL
       YH1K1=YH1K1/2. *EQVL
       YP1K1=YP1K1/2. *EQVL
       YP1K2=YP1K2/2. *EQVL
       YP1K01=YP1K01/2. *EQVL
       YP1K02=YP1K02/2. *EQVL
       YH10=YH10*EQVL/2.
       YH11=YH11*EQVL/2.
       YP11=YP11*EQVL/2.
       YP12=YP12*EQVL/2.
C     Transverse integrals
       YE1K0=YE1K0/2. *EQVL
       YE1KC0=YE1KC0/2. *EQVL
       YE1K1=YE1K1/2. *EQVL
       YE1KC1=YE1KC1/2. *EQVL
       YE1K2=YE1K2/2. *EQVL
       YE1KC2=YE1KC2/2. *EQVL
       YE10=YE10*EQVL/2.
       YE11=YE11*EQVL/2.
       YE12=YE12*EQVL/2.
C    in SIN
C     LONGITUDINAL INTEGRALS
       YH2K00=YH2K00/2. *EQVL
       YH2K01=YH2K01/2. *EQVL
       YH2K0=YH2K0/2. *EQVL
       YH2K1=YH2K1/2. *EQVL
       YP2K1=YP2K1/2. *EQVL
       YP2K2=YP2K2/2. *EQVL
       YP2K01=YP2K01/2. *EQVL
       YP2K02=YP2K02/2. *EQVL
       YH20=YH20*EQVL/2.
       YH21=YH21*EQVL/2.
       YP21=YP21*EQVL/2.
       YP22=YP22*EQVL/2.
C     TRANSVERSE INTEGRALS
       YE2K0=YE2K0/2. *EQVL
       YE2KC0=YE2KC0/2. *EQVL
       YE2K1=YE2K1/2. *EQVL
       YE2KC1=YE2KC1/2. *EQVL
       YE2K2=YE2K2/2. *EQVL
       YE2KC2=YE2KC2/2. *EQVL
       YE20=YE20*EQVL/2.
       YE21=YE21*EQVL/2.
       YE22=YE22*EQVL/2.
c100    CONTINUE
       RETURN
       END
       SUBROUTINE chasel
c  ...................................................................................
c  analysis of emittance by elimination of remote particles in  z-direction
c  ...................................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TAPES/IN,IFILE,META
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       LOGICAL CHASIT
       common/etchas/fractx,fracty,fractl
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/pool/zl(iptsz),ipin(iptsz)
       if(fractl.ge.1.) return
       do i=1,ngood
         ichas(i)=1
         zl(i)=0.
         ipin(i)=1
       enddo
       nl=int(float(ngood)*fractl)
       inz=0
       imaxf=ngood
       inx=0
c  ellipsoid of concentration (phase-dispersion)
       tx2=0.
       txp2=0.
       txxp=0.
       do i=1,ngood
         tx2=tx2+f(6,i)*f(6,i)
         txp2=txp2+fd(i)*fd(i)
         txxp=txxp+f(6,i)*fd(i)
       enddo
       tx2=tx2/float(ngood)
       txp2=txp2/float(ngood)
       txxp=txxp/float(ngood)
       delxxp=tx2*txp2-txxp*txxp
       ikept=0
       do i=1,ngood
         theta=pi/2.
         if(fd(i).ne.0.)theta=atan(f(6,i)/fd(i))
         rpart=f(6,i)*f(6,i)+fd(i)*fd(i)
         cos2=cos(theta)*cos(theta)
         sin2=sin(theta)*sin(theta)
         denom=tx2*cos2+txp2*sin2-2.*txxp*cos(theta)*sin(theta)
         relpse=2.5*delxxp/denom
         if(fractl.ge..97) then
           relpse=3.5*delxxp/denom
         else
           if(fractl.ge..95) relpse=3.*delxxp/denom
         endif
         if(rpart.le.relpse) then
           ipin(i)=0
           ikept=ikept+1
         endif
       enddo
       write(16,*) ' CHASEL:',fractl,'  % over: ',ngood-ikept,
     *             ' particles'
       do j=1,ngood
         if(ipin(j).eq.1) then
           inz=inz+1
           if(imaxf.le.nl) go to 9990
           imaxx=0
           tx=0.
           txp=0.
           txxp=0.
           tx2=0.
           txp2=0.
           do i=1,ngood
             if(ichas(i).eq.1) then
               tx=f(6,i)+tx
               txp=fd(i)+txp
               tx2=tx2+f(6,i)*f(6,i)
               txp2=txp2+fd(i)*fd(i)
               txxp=txxp+f(6,i)*fd(i)
               imaxx=imaxx+1
             endif
           enddo
           tx=tx/float(imaxx)
           txp=txp/float(imaxx)
           tx2=tx2/float(imaxx)
           txp2=txp2/float(imaxx)
           txxp=txxp/float(imaxx)
           xcg=tx
           zcg=txp
c       betatron parameters
           xxl=tx2-tx*tx
           zzl=txp2-txp*txp
           xzl=txxp-tx*txp
           if(inz.eq.1) go to 8880
           flcrit=2.*fl2rms*log(2.*imaxx)
           if(zlma.lt.flcrit) goto 7770
8880       emil=sqrt(xxl*zzl-xzl*xzl)
           if (emil.eq.0.) then
             al=0.
             bl=0.
             cl=0.
           else
             bl=sqrt(xxl/emil)
             cl=1./bl
             al=-xzl/emil
           endif
           tlx0=tx
           tlz0=txp
7770       fl2rms=(1.+al*al)*xxl/bl+2.*al*xzl+bl*zzl
           zlma=0.
           do i=1,ngood
             zl(i)=0.
             if(ichas(i).eq.1) then
               psx=f(6,i)-tlx0
               psz=fd(i)-tlz0
               zl(i)=psx*psx*cl*cl+(psx*al/bl+psz*bl)**2
               if (zlma.lt.zl(i)) then
                 zlma=zl(i)
                 izlma=i
               endif
             endif
           enddo
c    particle is eliminated
           imaxf=0
           if (zlma.eq.0.) zlma=1.e10
           do i=1,ngood
             if(ichas(i).eq.1 .and. zl(i).lt.zlma) then
               imaxf=imaxf+1
             else
               ichas(i)=0
             endif
           enddo
         endif
       enddo
9990   continue
       return
       end
       SUBROUTINE chasex
c  ...................................................................................
c  analysis of emittance by elimination of remote particles in x-direction
c  ...................................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/tapes/in,ifile,meta
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       LOGICAL CHASIT
       common/etchas/fractx,fracty,fractl
       common/pool/zl(iptsz),ipin(iptsz)
       if(fractx.ge.1.) return
       do i=1,ngood
         ichas(i)=1
         zl(i)=0.
         ipin(i)=1
       enddo
       nl=int(float(ngood)*fractx)
       inz=0
       imaxf=ngood
       inx=0
c  ellipsoid of concentration (x-xp)
       tx2=0.
       txp2=0.
       txxp=0.
       do i=1,ngood
         f2=f(2,i)*1.e-02
         f3=f(3,i)*1.e-03
         tx2=tx2+f2*f2
         txp2=txp2+f3*f3
         txxp=txxp+f2*f3
       enddo
       tx2=tx2/float(ngood)
       txp2=txp2/float(ngood)
       txxp=txxp/float(ngood)
       delxxp=tx2*txp2-txxp*txxp
       ikept=0
       do i=1,ngood
         theta=pi/2.
         f2=f(2,i)*1.e-02
         f3=f(3,i)*1.e-03
         if(f3.ne.0.)theta=atan(f2/f3)
         rpart=f2*f2+f3*f3
         cos2=cos(theta)*cos(theta)
         sin2=sin(theta)*sin(theta)
         denom=tx2*cos2+txp2*sin2-2.*txxp*cos(theta)*sin(theta)
         relpse=2.5*delxxp/denom
         if(fractx.ge..97) then
           relpse=3.5*delxxp/denom
         else
           if(fractx.ge..95) relpse=3.*delxxp/denom
         endif
         if(rpart.le.relpse) then
           ipin(i)=0
           ikept=ikept+1
         endif
       enddo
       write(16,*) ' CHASEX:',fractx,'  % over: ',ngood-ikept,
     *             ' particles'
       do j=1,ngood
         if(ipin(j).eq.1) then
           inz=inz+1
           if(imaxf.le.nl) go to 9990
           imaxx=0
           tx=0.
           txp=0.
           txxp=0.
           tx2=0.
           txp2=0.
           do i=1,ngood
             if(ichas(i).eq.1) then
               tx=f(2,i)+tx
               txp=f(3,i)+txp
               tx2=tx2+f(2,i)*f(2,i)
               txp2=txp2+f(3,i)*f(3,i)
               txxp=txxp+f(2,i)*f(3,i)
               imaxx=imaxx+1
             endif
           enddo
           tx=tx/float(imaxx)
           txp=txp/float(imaxx)
           tx2=tx2/float(imaxx)
           txp2=txp2/float(imaxx)
           txxp=txxp/float(imaxx)
           xcg=tx
           zcg=txp
c       betatron parameters
           xxl=tx2-tx*tx
           zzl=txp2-txp*txp
           xzl=txxp-tx*txp
           if(inz.eq.1) go to 8880
           flcrit=2.*fl2rms*log(2.*imaxx)
           if(zlma.lt.flcrit) goto 7770
8880       emil=sqrt(xxl*zzl-xzl*xzl)
           if (emil.eq.0.) then
             al=0.
             bl=0.
             cl=0.
           else
             bl=sqrt(xxl/emil)
             cl=1./bl
             al=-xzl/emil
           endif
           tlx0=tx
           tlz0=txp
7770       fl2rms=(1.+al*al)*xxl/bl+2.*al*xzl+bl*zzl
           zlma=0.
           do i=1,ngood
             zl(i)=0.
             if(ichas(i).eq.1) then
               psx=f(2,i)-tlx0
               psz=f(3,i)-tlz0
               zl(i)=psx*psx*cl*cl+(psx*al/bl+psz*bl)**2
               if (zlma.lt.zl(i)) then
                 zlma=zl(i)
                 izlma=i
               endif
             endif
           enddo
c    particle is eliminated
           imaxf=0
           if (zlma.eq.0.) zlma=1.e10
           do i=1,ngood
             if(ichas(i).eq.1 .and. zl(i).lt.zlma) then
               imaxf=imaxf+1
             else
               ichas(i)=0
             endif
           enddo
         endif
       enddo
9990   continue
       return
       end
       SUBROUTINE chase
c  ...................................................................................
c  analysis of emittance by elimination of remote particles (read the parameters)
c  ...................................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/tapes/in,ifile,meta
       common/etchas/fractx,fracty,fractl
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       LOGICAL CHASIT
       READ(in,*)FRACTX,FRACTY,FRACTL
       WRITE(16,1) FRACTX*100.,FRACTY*100.,FRACTL*100.
1      FORMAT('    ARE KEPT IN THE BUNCH ',/,
     X         '       (x,xp)   : ',F7.3,' %',/,
     X         '       (y,yp)   : ',F7.3,' % ',/,
     X         '       (w,phase):',F7.3,' %')
C
       CHASIT=.FALSE.
       if(fractx.lt.1.) CHASIT=.TRUE.
       if(fracty.lt.1.) CHASIT=.TRUE.
       if(fractl.lt.1.) CHASIT=.TRUE.
       RETURN
       END
       SUBROUTINE chasey
c  ...................................................................................
c  analysis of emittance by elimination of remote particles in y-direction
c  ...................................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/tapes/in,ifile,meta
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       LOGICAL CHASIT
       common/etchas/fractx,fracty,fractl
       common/pool/zl(iptsz),ipin(iptsz)
       if(fracty.ge.1.) return
       do i=1,ngood
         ichas(i)=1
         zl(i)=0.
         ipin(i)=1
       enddo
c  ellipsoid of concentration (y-yp)
       ty2=0.
       typ2=0.
       tyyp=0.
       do i=1,ngood
         f4=f(4,i)*1.e-02
         f5=f(5,i)*1.e-03
         ty2=ty2+f4*f4
         typ2=typ2+f5*f5
         tyyp=tyyp+f4*f5
       enddo
       ty2=ty2/float(ngood)
       typ2=typ2/float(ngood)
       tyyp=tyyp/float(ngood)
       delyyp=ty2*typ2-tyyp*tyyp
       ikept=0
       do i=1,ngood
         f4=f(4,i)*1.e-02
         f5=f(5,i)*1.e-03
         theta=pi/2.
         if(f5.ne.0.)theta=atan(f4/f5)
         rpart=f4*f4+f5*f5
         cos2=cos(theta)*cos(theta)
         sin2=sin(theta)*sin(theta)
         denom=ty2*cos2+typ2*sin2-2.*tyyp*cos(theta)*sin(theta)
         relpse=2.5*delyyp/denom
         if(fracty.ge..97) then
           relpse=3.5*delyyp/denom
         else
           if(fracty.ge..95) relpse=3.*delyyp/denom
         endif
         if(rpart.le.relpse) then
           ipin(i)=0
           ikept=ikept+1
         endif
       enddo
       write(16,*) ' CHASEY:',fracty,'  % over: ',ngood-ikept,
     *             ' particles'
       imaxf=ngood
       nl=int(float(ngood)*fracty)
       inz=0
       do j=1,ngood
         if(ipin(j).eq.1) then
           inz=inz+1
           if(imaxf.le.nl) go to 9990
           imaxx=0
           ty=0.
           typ=0.
           tyyp=0.
           ty2=0.
           typ2=0.
           do i=1,ngood
             if(ichas(i).eq.1) then
               ty=f(4,i)+ty
               typ=f(5,i)+typ
               ty2=ty2+f(4,i)*f(4,i)
               typ2=typ2+f(5,i)*f(5,i)
               tyyp=tyyp+f(4,i)*f(5,i)
               imaxx=imaxx+1
             endif
           enddo
           ty=ty/float(imaxx)
           typ=typ/float(imaxx)
           ty2=ty2/float(imaxx)
           typ2=typ2/float(imaxx)
           tyyp=tyyp/float(imaxx)
           xcg=ty
           zcg=typ
c       betatron parameters
           xxl=ty2-ty*ty
           zzl=typ2-typ*typ
           xzl=tyyp-ty*typ
           if(inz.eq.1) go to 8880
           flcrit=2.*fl2rms*log(2.*imaxx)
           if(zlma.lt.flcrit) goto 7770
8880       emil=sqrt(xxl*zzl-xzl*xzl)
           if (emil.eq.0.) then
             al=0.
             bl=0.
             cl=0.
           else
             bl=sqrt(xxl/emil)
             cl=1./bl
             al=-xzl/emil
           endif
           tlx0=ty
           tlz0=typ
7770       fl2rms=(1.+al*al)*xxl/bl+2.*al*xzl+bl*zzl
           zlma=0.
           do i=1,ngood
             zl(i)=0.
             if(ichas(i).eq.1) then
               psx=f(4,i)-tlx0
               psz=f(5,i)-tlz0
               zl(i)=psx*psx*cl*cl+(psx*al/bl+psz*bl)**2
               if (zlma.lt.zl(i)) then
                 zlma=zl(i)
                 izlma=i
               endif
             endif
           enddo
c    particle is eliminated
           imaxf=0
           if (zlma.eq.0.) zlma=1.e10
           do i=1,ngood
             if(ichas(i).eq.1 .and. zl(i).lt.zlma) then
               imaxf=imaxf+1
             else
               ichas(i)=0
             endif
           enddo
         endif
       enddo
9990   continue
       return
       end
       SUBROUTINE corre(n,nall)
c  ...........................................................................
c     correction over the beam generated by MONTE
c  ...........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /COM4/CORD(IPTSZ,6)
       DIMENSION A(6,6),B(6,6),C(6,6),D(6),E(6),F(6),G(6)
       DO I =1,6
         D(I)=0.0
         E(I)=0.0
         F(I)=0.0
         G(I)=0.0
         DO J =1,6
           A(I,J)=0.0
           B(I,J)=0.0
           C(I,J)=0.0
         ENDDO
       ENDDO
       DO J =1,6
         DO I =1,N
           D(J)=D(J)+CORD(I,J)
         ENDDO
         D(J)=D(J)/N
       ENDDO
       DO J =1,6
         DO I=1,NALL
           CORD(I,J)=CORD(I,J)-D(J)
         ENDDO
       ENDDO
       DO J=1,6
         DO K=1,J
           DO I=1,N
             A(J,K)=A(J,K)+ CORD(I,J)*CORD(I,K)
           ENDDO
           A(J,K)=A(J,K)/N
         ENDDO
       ENDDO
C WE MAKE MATRIX A=B*BT WHERE B,BT ARE TRIANGULAR
       DO I=1,6
         DO J=1,I
           H=A(I,J)
           J1=J-1
           IF(J.EQ.1)GO TO 8
           DO K=1,J1
             H=H-A(I,K)*A(J,K)
           ENDDO
8          IF(I.EQ.J)GO TO 10
           A(I,J)=H/A(J,J)
           GO TO 7
10         A(I,J)=SQRT(H)
7          B(I,J)=A(I,J)
         ENDDO
       ENDDO
C  WE INVERT B TO GIVE AN UPDATED B
       E(1)=1.0/B(1,1)
       DO I=2,6
         E(I)= 1.0/B(I,I)
         J1=I-1
         DO JJ=1,J1
           J=I-JJ
           J3=J+1
           S=0.0
           IF(JJ.EQ.1)GO TO 13
           DO K=J3,J1
             S=S-B(K,I)*B(K,J)
           ENDDO
13         B(J,I)=(S-B(I,J)*E(I))/B(J,J)
         ENDDO
       ENDDO
       DO I=1,6
         B(I,I)=E(I)
       ENDDO
       DO I=1,6
         DO J=1,I
           DUM=B(I,J)
           B(I,J)=B(J,I)
           B(J,I)=DUM
         ENDDO
       ENDDO
C  WE CONVERT CORD(**) SO THAT ITS SIGMA MATRIX IS UNITY
       DO I=1,NALL
         DO K=1,6
           F(K)=0.0
           DO J=1,K
             F(K)=F(K)+ B(K,J)*CORD(I,J)
           ENDDO
         ENDDO
         DO K=1,6
           CORD(I,K)=F(K)
         ENDDO
       ENDDO
C WE TEST THE MEANS AND SIGMA MATRIX
       DO J=1,6
         DO I=1,N
           G(J)=G(J)+CORD(I,J)
         ENDDO
         G(J)=G(J)/N
       ENDDO
       DO J=1,6
         DO I=1,N
           CORD(I,J)=CORD(I,J)-G(J)
         ENDDO
       ENDDO
       DO J=1,6
         DO K=1,6
           DO I=1,N
             C(J,K)=C(J,K)+ CORD(I,J)*CORD(I,K)
           ENDDO
           C(J,K)=C(J,K)/N
         ENDDO
       ENDDO
       RETURN
       END
       SUBROUTINE pintim
c   ..................................................................
C   Shifts particle coordinates to a single point in time. Uses
C   a linear shift
C   Divide by 100. to convert from cm to meters
c    called by SCHEFF or SCHERM
c   ..................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/PART/XC(iptsz),YC(iptsz),ZC(iptsz)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/AZLIST/ICONT,IPRIN
       grmoy=0.
       trmoy=0.
       xbax=0.
       do i=1,ngood
         grmoy=grmoy+f(7,i)/xmat
         trmoy=trmoy+f(6,i)
         xbax=xbax+f(2,i)
       enddo
       trmoy=trmoy/float(ngood)
       grmoy=grmoy/float(ngood)
       brmoy=sqrt(1.-1./(grmoy*grmoy))
       xbax=xbax/float(ngood)
       apl=0.
c  Isochronism correction  (bending magnet) only with SCHERM
c    does not work with  with SCHEFF  (iscsp=3)
       if(iscsp.eq.2) then
         xb2x=0.
         xb2z=0.
         xbxz=0.
         do np=1,ngood
           gpai=f(7,np)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           zc(np)=(trmoy-f(6,np))*bpai*vl/100.
           xc(np)=(f(2,np)-xbax)/100.
           xb2z=xb2z+zc(np)*zc(np)
           xb2x=xb2x+xc(np)*xc(np)
           xbxz=xbxz+zc(np)*xc(np)
         enddo
         xb2z=xb2z/float(ngood)
         xb2x=xb2x/float(ngood)
         xbxz=xbxz/float(ngood)
         apl=atan(-2.*xbxz/(xb2x-xb2z))/2.
         write(16,*) 'slope of the bunch in plane(Oz,Ox):',apl,
     *               ' radian'
       endif
       do np=1,ngood
         gpai=f(7,np)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
c      iscsp = 3 Lorentz transformation (only with scheff)
comment         if(iscsp.eq.3) znp=(trmoy-f(6,np))*bpai*vl*grmoy
comment         if(iscsp.eq.2) znp=(trmoy-f(6,np))*bpai*vl
         znp=(trmoy-f(6,np))*bpai*vl
         xnp=f(2,np)
         zc(np)=znp*cos(apl)+xnp*sin(apl)
         xnp=xnp*cos(apl)-znp*sin(apl)
C        convert from mrad to rad
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
C        convert from cm   to m
         xc(np)=(xnp+zc(np)*f3)/100.
         yc(np)=(f(4,np)+zc(np)*f5)/100.
         zc(np)=zc(np)/100.
       enddo
       xbar=0.
       ybar=0.
       zbar=0.
       do np=1,ngood
c      evaluate xbar , ybar , zbar
         xbar=xbar+xc(np)
         ybar=ybar+yc(np)
         zbar=zbar+zc(np)
       enddo
       xbar=xbar/float(ngood)
       ybar=ybar/float(ngood)
       zbar=zbar/float(ngood)
c  Translate distribution by center of mass coordinates to shift
c  coordinate origin to (0,0,0)
       do np=1,ngood
         xc(np)=xc(np)-xbar
         yc(np)=yc(np)-ybar
         zc(np)=zc(np)-zbar
       enddo
       return
       end
       SUBROUTINE schermi1
      implicit real*8 (a-h,o-z)
c   ..................................................................
c     Called by SCHERMI when the bunch can be represented
c      by a simple ellipse in the longitudinal direction
c   ..................................................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/ELCG/XCGD,YCGD,ZCGD,XCGR,YCGR,ZCGR
       COMMON/INTGRT/ex,ey,ez
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/CDEK/DWP(iptsz)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/NPART/IMAXR
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/twcst/epsilon
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/DCSPA/IESP
       COMMON/CMPTE/IELL
       COMMON/CGRMS/xsum,ysum,zsum
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/DIMENS/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       common/posc/xpsc
       LOGICAL ichaes,iesp
       dimension afx(20),afy(20)
       IF (beamc.eq.0. .OR. scdist.eq.0.) return
       IELL=IELL+1
       imaxf=ngood
       wavel=2.*pi*vl/fh
       xmass=xmat*1.78267581E-30
       NMAXY=5
       dxp=0.
       dyp=0.
       dw=0.
       dz=scdist/100.
       dz1=dz
       if(ngood.eq.0) then
         write(16,*) ' all the particles are lost '
         stop
       endif
c xi in Amps, beamc in mA
c epsilon =(coul*coul)/nt*(m*m)
       epsilon=8.854189586e-12
       c1=1./(3.*pi*sqrt(5.))
       cl=vl/100.
c charge per macro particle
comment       const1=c1*xi*xi*wavel*wavel/(10000.*imaxf*xmass*
c             *         epsilon*cl*cl*cl*cl)
       const3=1.E-06
       const2=1.E-06/xmat
c calculate rms beam size for beam in one point in time
       CALL sizrms(0,xrms,yrms,zrms,zz)
c6875   format(2x,e12.5,2x,e12.5,2x,e12.5)
       xrmsp=xrms
       yrmsp=yrms
       zrmsp=zrms
       CALL sizcor(ect,xrms,yrms,zrms,0)
       xrmsc=xrms
       yrmsc=yrms
       zrmsc=zrms
       xrms1=xrmsp
       yrms1=yrmsp
       zrms1=zrmsp
       xrms=xrmsp
       yrms=yrmsp
       zrms=zrmsp
       xcgd=xsum
       ycgd=ysum
       zcgd=zsum
       do i=1,ngood
         ZCP(I)=ZC(I)
         XCP(I)=XC(I)
         YCP(I)=YC(I)
       enddo
c      limits in z-direction
       zmat=0.
       zmit=1000.
       do i=1,ngood
         IF(ZCP(I).GE.ZMAT) ZMAT=ZCP(I)
         IF(ZCP(I).LT.ZMIT) ZMIT=ZCP(I)
       enddo
       zmat=zmat/zrms
       zmit=zmit/zrms
c    extends zmat
       zmat=zmat+zmat*.50
       zmit=zmit+zmit*.50
       if(zmat.gt.ect) zmat=ect
       if(abs(zmit).gt.ect) zmit=-ect
c    Hermite coefficients in x and y-direction
       nmaz=0
c6876   format(2x,i3,2x,i3)
       do k=1,20
         afzt(k)=0.
         afxt(k)=0.
         afyt(k)=0.
         afzm(k)=0.
         afxm(k)=0.
         afym(k)=0.
       enddo
       do k=1,nmaxy
         kap=k-1
         do j=1,ngood
           xcoup=abs(xcp(j)/xrms)
           ycoup=abs(ycp(j)/yrms)
           zcoup=abs(zcp(j)/zrms)
           if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
             xc(j)=xcp(j)/xrmsc
             yc(j)=ycp(j)/yrmsc
             zc(j)=zcp(j)/zrmsc
             afxm(k)=afxm(k)+herm(2*kap,xc(j))
             afym(k)=afym(k)+herm(2*kap,yc(j))
             AFZM(K)=AFZM(K)+HERM(2*KAP,ZC(J))
           endif
         enddo
         afxm(k)=afxm(k)/(fact(2*kap)*sqrt(2.*pi))
         afym(k)=afym(k)/(fact(2*kap)*sqrt(2.*pi))
         AFZm(K)=AFZm(K)/(FACT(2*KAP)*SQRT(2.*PI))
       enddo
       NMAZ=10
c   Hermite coefficients in z-direction
       zcdg=0.
       imaxx=0
       do j=1,ngood
         xcoup=abs(xcp(j)/xrms)
         ycoup=abs(ycp(j)/yrms)
         zcoup=abs(zcp(j)/zrms)
         if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
           zcdg=zcdg+zcp(j)
           imaxx=imaxx+1
         endif
       enddo
       zcdg=zcdg/float(imaxx)
       zsqsum=0.
       zcub=0.
       zcub1=0.
       do j=1,ngood
         xcoup=abs(xcp(j)/xrms)
         ycoup=abs(ycp(j)/yrms)
         zcoup=abs(zcp(j)/zrms)
         if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
           zc(j)=zcp(j)-zcdg
           zsqsum=zsqsum+zc(j)*zc(j)
           zcub=zcub+zc(j)*zc(j)*zc(j)
           zcub1=zcub1+zc(j)
         endif
       enddo
       zrmsz=zsqsum/float(imaxx)
       zrmsz=sqrt(zrmsz)
       zcub3=zcub
       zcub=zcub/(zrmsz*zrmsz*zrmsz)-3.*zcub1/zrmsz
       zcub=zcub/(6.*sqrt(2.*pi))
       do k=1,nmaz
         afzt(k)=0.
         afx(k)=0.
         afy(k)=0.
         kap=k-1
         do j=1,ngood
           xcoup=abs(xcp(j)/xrms)
           ycoup=abs(ycp(j)/yrms)
           zcoup=abs(zcp(j)/zrms)
           if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
             zc(j)=zcp(j)/zrmsc
             xc(j)=xcp(j)/xrmsc
             yc(j)=ycp(j)/yrmsc
             afzt(k)=afzt(k)+herm(kap,zc(j))
             afx(k)=afx(k)+herm(kap,xc(j))
             afy(k)=afy(k)+herm(kap,yc(j))
           endif
         enddo
         afzt(k)=afzt(k)/(fact(kap)*sqrt(2.*pi))
         afx(k)=afx(k)/(fact(kap)*sqrt(2.*pi))
         afy(k)=afy(k)/(fact(kap)*sqrt(2.*pi))
       enddo
c Do Hermite integration to determine Ex,Ey,Ez for each macro particle
c and apply space charge kick. Field components are passed through the
c common INTGRT in units Newton/Coulomb
       DO i=1,ngood
c      reprise coordonnees spaciales pour calcul des champs
         xc(i)=xcp(i)
         yc(i)=ycp(i)
         zc(i)=zcp(i)
         CALL INTGA(i,0)
         ext=ex
         eyt=ey
         ezt=ez
c calculate kick in x',y' and z' (energy)
c calculate kick in x',y' and z' (energy)
c   isochronism correction
         eztp=ezt*cos(apl)-ext*sin(apl)
         extp=ezt*sin(apl)+ext*cos(apl)
         ezt=eztp
         ext=extp
         gsc=f(7,i)/xmat
         bsc=sqrt(1.-1./(gsc*gsc))
c        *  valero
         dxp=const2*ext*dz/(bsc*bsc*gsc*gsc*gsc)*abs(f(9,i))
         dyp=const2*eyt*dz/(bsc*bsc*gsc*gsc*gsc)*abs(f(9,i))
         dw=const3*ezt*dz*abs(f(9,i))/gsc
c        *
         if(.not.iesp) then
c     load the entrance beam in cavities or gaps
           do js=1,7
             f(js,i)=fs(js,i)
           enddo
           f(3,i)=f(3,i)+dxp*1000.
           f(5,i)=f(5,i)+dyp*1000.
           f(2,i)=f(2,i)-dz1*dxp*100.*xpsc
           f(4,i)=f(4,i)-dz1*dyp*100.*xpsc
           dwp(i)=dw
         else
           f(3,i)=f(3,i)+dxp*1000.
           f(5,i)=f(5,i)+dyp*1000.
           f(7,i)=f(7,i)+dw
         endif
       enddo
       RETURN
       END
       SUBROUTINE schermi
c   ..................................................................
c    SCHERM space charge method
c     See NIM A 309(1996) 21-40
c   ..................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/HERMD/AFXM(20),AFYM(20),AFZM(20)
       COMMON/HERMR/AFXR(20),AFYR(20),AFZR(20)
       COMMON/HERMRR/AFXRR(20),AFYRR(20),AFZRR(20)
       COMMON/SIZR/XRMS3,YRMS3,ZRMS3,ZCGR3
       COMMON/SIZT/XRMS,YRMS,ZRMS
       COMMON/SIZP/XRMS1,YRMS1,ZRMS1,XRMS2,YRMS2,ZRMS2,IMAXD
       COMMON/ELCG/XCGD,YCGD,ZCGD,XCGR,YCGR,ZCGR
       COMMON/INTGRT/ex,ey,ez
       COMMON/DEGHERM/NMAZ,NMAZR,NMAXY
       COMMON/CDEK/DWP(iptsz)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/NPART/IMAXR
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/twcst/epsilon
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/DCSPA/IESP
       COMMON/CMPTE/IELL
       COMMON/CHAMP/FXRMS(10,15),FYRMS(10,15),FZRMS(10,15),NCHAMP(10),
     X NCCHAM(10),NCHPAS,JCHAM,ITYE
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/posc/xpsc
       LOGICAL ichaes,iesp
C *************************************************
c      nquad : number of quad in the transport line
       COMMON /CHQUA/ICQD,nquad
       LOGICAL ICQD
C *************************************************
       COMMON/DIMENS/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       IF (beamc.eq.0. .OR. scdist.eq.0.) return
c    dummy: in order to have the same entry that HERSC and SCHEFF
c Initialize some constants and variables c wavel in cm
       ect=4.
       call shuffle
       call pintim
       IELL=IELL+1
       imaxf=ngood
       wavel=2.*pi*vl/fh
       xmass=xmat*1.78267581E-30
       NMAZ=14
       NMAZR=8
       NMAXY=4
       dxp=0.
       dyp=0.
       dw=0.
       dz=scdist/100.
       dz1=dz
       zsot1=0.
       if(ngood.eq.0) then
         write(16,*) ' all the particles are lost '
         stop
       endif
       WRITE(16,*)' call SCHERM N: ',iell
c xi in Amps, beamc in mA
c epsilon =(coul*coul)/nt*(m*m)
       epsilon=8.854189586e-12
       c1=1./(3.*pi*sqrt(5.))
       cl=vl/100.
c charge per macro particle
       const3=1.E-06
       const2=1.E-06/xmat
c  normalized emittances(normalized) in m.radian
c calculate rms beam size for beam in one point in time
       CALL sizrms(0,xrms,yrms,zrms,zz)
       xrmsp=xrms
       yrmsp=yrms
       zrmsp=zrms
       CALL sizcor(ect,xrms,yrms,zrms,0)
       write(16,*) ' bunch RMS(m): ',xrms,yrms,zrms
c6875   format(2x,e12.5,2x,e12.5,2x,e12.5)
       xrmsc=xrms
       yrmsc=yrms
       zrmsc=zrms
       xrms=xrmsp
       yrms=yrmsp
       zrms=zrmsp
c   Total bunch charge densities in x,y,z
c      calculation of Hermite coefficients
       do i=1,ngood
         zcp(i)=zc(i)
         xcp(i)=xc(i)
         ycp(i)=yc(i)
       enddo
c      limits in z-direction
       zmat=0.
       zmit=1000.
       do i=1,ngood
         if(zcp(i).ge.zmat) zmat=zcp(i)
         if(zcp(i).lt.zmit) zmit=zcp(i)
       enddo
       zmat=zmat/zrms
       zmit=zmit/zrms
c    extend zmat
       zmat=zmat+zmat*.50
       zmit=zmit+zmit*.50
c    (xmat,ymat,zmat) >0
       if(zmat.gt.ect) zmat=ect
       if(abs(zmit).gt.ect) zmat=-ect
c   Hermite coefficients
       do k=1,nmaxy
         afxt(k)=0.
         afyt(k)=0.
         kap=k-1
         do j=1,ngood
           xcoup=abs(xcp(j)/xrms)
           ycoup=abs(ycp(j)/yrms)
           zcoup=abs(zcp(j)/zrms)
           if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
             xc(j)=xcp(j)/xrmsc
             yc(j)=ycp(j)/yrmsc
             afxt(k)=afxt(k)+herm(2*kap,xc(j))
             afyt(k)=afyt(k)+herm(2*kap,yc(j))
           endif
         enddo
         afxt(k)=afxt(k)/(fact(2*kap)*sqrt(2.*pi))
         afyt(k)=afyt(k)/(fact(2*kap)*sqrt(2.*pi))
       enddo
c6876   format(2x,i5,3x,e12.5,3x,e12.5)
       do k=1,nmaz
         afzt(k)=0.
         kap=k-1
         do j=1,ngood
           xcoup=abs(xcp(j)/xrms)
           ycoup=abs(ycp(j)/yrms)
           zcoup=abs(zcp(j)/zrms)
           if(xcoup.le.ect.and.ycoup.le.ect.and.zcoup.le.ect) then
             zc(j)=zcp(j)/zrms
             afzt(k)=afzt(k)+herm(kap,zc(j))
           endif
         enddo
         afzt(k)=afzt(k)/(fact(kap)*sqrt(2.*pi))
       enddo
c6877   format(3x,e12.5)
C     Seek for surface of the bunch
       szbt=snzt(zmit,zmat)
c      Seek for the vertex of the distribution n(z)
c        between  0 and zmat/2
       zf=zmat/2.
       zi=0.
       call rchsom(zi,zf,nmaz)
       zsot=(zi+zf)/2.
c9999   CONTINUE
       zcgd=zsot*zrms
C   Main ellipsoid (in the right of the vertex)
       imaxd=1
       do i=1,ngood
         if(zcp(i).ge.zcgd) then
           xc(imaxd)=xcp(i)
           yc(imaxd)=ycp(i)
           zc(imaxd)=zcp(i)
           imaxd=imaxd+1
         endif
       enddo
       imaxd=imaxd-1
       imaxr=ngood-2*imaxd
       if(imaxr.lt.0) then
c if imaxr<0 => use one ellips only, rather than 2
         pcent1=float(2*imaxd)/float(ngood)
         pcent2=float(imaxr)/float(ngood)
         if(icqd)NQUAD=NQUAD-1
         iell=iell-1
         write(16,*) ' one ellipsoid in z-direction '
         call pintim
         call schermi1
         return
       endif
       IF(25*IMAXR.LT.IMAXF.OR.IMAXR.LE.30) THEN
         pcent1=float(2*imaxd)/float(ngood)
         pcent2=float(imaxr)/float(ngood)
         iell=iell-1
         if(icqd)NQUAD=NQUAD-1
c if imaxr<.04*imaxf => use one ellips only, rather than 2
          write(16,*) 'one ellipsoid in z-direction '
         call pintim
         call schermi1
         IF (12*IMAXR.LT.ngood)
     x     write(16,*) 'one ellipsoid in z-direction '
         IF (IMAXR.LE.40)
     x     write(16,*) ' one ellipsoid in z-direction '
         return
       endif
       i1elli=0
c    C.O.G. of the principal ellipsoid (in cm)
c      longitudinal RMS of the pricipal ellipsoid
       zrmss1=sqrt(vaprz(zsot,zmat))
       zrms1=zrmss1*zrms
c   RMS in the transverse directions
       xrms1=xrms
       yrms1=yrms
c      surface n(z)
       szbd=snzd(zsot,zmat)
       rsnz=szbd/szbt
c   coefficients Hermite en x,y,z sur partie droite en z
c    principal ellipsoid
       do k=1,nmaxy
         kap=k-1
         afxm(k)=afxt(k)*rsnz
         afym(k)=afyt(k)*rsnz
         afzm(k)=prinz(zsot,zmat,k,zrmss1)
         afzm(k)=2.*afzm(k)/(fact(2*kap)*sqrt(2.*pi))
       enddo
c    limits in x,y,z of the principal ellipsoid
       do i=1,imaxd
         xc(i)=xc(i)/xrms
         yc(i)=yc(i)/yrms
         zc(i)=zc(i)/zrms
       enddo
       xmam=xc(1)
       ymam=yc(1)
       zmam=zc(1)
       do i=1,imaxd
         if(xc(i).ge.xmam) xmam=xc(i)
         if(yc(i).ge.ymam) ymam=yc(i)
         if(zc(i).ge.zmam) zmam=zc(i)
       enddo
       if(abs(xmam).ge.ect) xmam=ect
       if(abs(ymam).ge.ect) ymam=ect
       if(abs(zmam).ge.ect) zmam=ect
       xmim=xmam
       ymim=ymam
       zmim=zmam
       do i=1,imaxd
         if(xc(i).lt.xmim) xmim=xc(i)
         if(yc(i).lt.ymim) ymim=yc(i)
         if(zc(i).lt.zmim) zmim=zc(i)
       enddo
       if(abs(xmim).ge.ect) xmim=-ect
       if(abs(ymim).ge.ect) ymim=-ect
       if(abs(zmim).ge.ect) zmim=-ect
c      partition commune @ x,y,z
       xymam=ymam
       xymim=ymim
       if(xmam.ge.xymam) xymam=xmam
       if(zmam.ge.xymam) xymam=zmam
       if(xmim.lt.xymim) xymim=xmim
       if(zmim.lt.xymim) xymim=zmim
c      symmetrisation de l'intervalle
       if(abs(xymim).ge.xymam) then
         xymam=abs(xymim)
       else
         xymim=-xymam
       endif
c    vertex of the second ellipsoid
       aa=zmit
       bb=2.*zsot-zmat
       cc=zsot
       dd=zmat
       call rchsor(aa,bb,cc,dd,ee)
c second ellipsoid  around ee
c    sz2e : surface
       zrms2=varia(bb,cc,dd,ee)
       zrms2=sqrt(zrms2)
       zcgr=ee*zrms
       sz2e=codsy(bb,cc,dd,ee,1)
       afzr(1)=sz2e/sqrt(2.*pi)
       stm12=abs(afzt(1)-afzm(1)-afzr(1))
       rs2e=sz2e/szbt
       xrms2=xrms
       yrms2=yrms
       do k=1,nmaxy
         kap=k-1
         afxr(k)=afxt(k)*rs2e
         afyr(k)=afyt(k)*rs2e
         afzr(k)=codsy(bb,cc,dd,ee,k)
         afzr(k)=afzr(k)/(fact(2*kap)*sqrt(2.*pi))
       enddo
c   two ellipsoids:
       if(stm12*10.le.afzr(1)) inint=2
       if(stm12*10.gt.afzr(1)) then
c   3 ellipsoids;the principal ellipsoid,the second ellipsoid is symmetrized
c    the third ellipsoid is defined around the c. of g. of the residu
         ee1=grz(aa,bb,cc,dd,ee)
         if(ee1.gt.ee) then
           inint=2
           go to 1968
         endif
         zrms3=variz(bb,cc,dd,ee,ee1)
         zrms3=sqrt(zrms3)
         xrms3=xrms
         yrms3=yrms
         zcgr3=ee1*zrms
C      C.O.G. of the second ellipsoid
c        on suppose xcgr et ycgr nuls pour le residu
         xcgr=0.
         ycgr=0.
         sz3e=codif(bb,cc,dd,ee,ee1,1)
         rs3e=sz3e/szbt
c        correction complementaire
         srtot=sz3e/sqrt(2.*pi)+afzr(1)+afzm(1)
         srtot=srtot-afzt(1)
C       HErMITE coefficients over the residual  ellips
         do k=1,nmaxy
           kap=k-1
           afxrr(k)=afxt(k)*rs3e
           afyrr(k)=afyt(k)*rs3e
           afzrr(k)=codif(bb,cc,dd,ee,ee1,k)
           afzrr(k)=afzrr(k)/(fact(2*kap)*sqrt(2.*pi))
           if(k.eq.1) then
             tzrr=afzrr(k)-srtot
             if(tzrr.ge.0..and.abs(tzrr).ge.(afzrr(k)/10.)) then
               afxrr(k)=afxrr(k)-srtot
               afyrr(k)=afyrr(k)-srtot
               afzrr(k)=afzrr(k)-srtot
             endif
           endif
         enddo
c end of calculus of 2nd ellips
         inint=3
       endif
       zrms3=zrms3*zrms
1968   continue
       zrms2=zrms2*zrms
c    surface of the ellipsoids in pourcent
       write(16,*) ' surface of the ellipsoids in % of the bunch:'
       if(inint.eq.2) then
         pcent1=afzm(1)/afzt(1)
         pcent2=afzr(1)/afzt(1)
         write(16,7777) iell,pcent1,pcent2
       endif
       if(inint.eq.3) then
         pcent1=afzm(1)/afzt(1)
         pcent2=afzr(1)/afzt(1)
         pcent3=afzrr(1)/afzt(1)
         write(16,7778) iell,pcent1,pcent2,pcent3
       endif
7777   format(2x,i4,2x,f7.4,2x,f7.4)
7778   format(2x,i4,2x,f7.4,2x,f7.4,2x,f7.4)
c Do Hermite integration to determine Ex,Ey,Ez for each macro particle
c and apply space charge kick. Field components are passed through the
c common INTGRT in units Newton/Coulomb
c       gsct introduit pour calcul specifique impulsion
       gsct=0.
       igsct=0
       DO i=1,ngood
c      reprise coordonnees spaciales pour calcul des champs
         xc(i)=xcp(i)
         yc(i)=ycp(i)
         zc(i)=zcp(i)
c   principal ellipsoid
         CALL INTGA(i,0)
         ext=ex
         eyt=ey
         ezt=ez
c  second ellipsoid
         CALL INTGA(i,1)
         ext=ext+ex
         eyt=eyt+ey
         ezt=ezt+ez
c  third ellipsoid
         if(inint.eq.3) then
           CALL INTGA(i,2)
           ext=ext+ex
           eyt=eyt+ey
           ezt=ezt+ez
         endif
c calculate kick in x',y' and z' (energy)
c   isochronism correction
         eztp=ezt*cos(apl)-ext*sin(apl)
         extp=ezt*sin(apl)+ext*cos(apl)
         ezt=eztp
         ext=extp
         gsc=f(7,i)/xmat
         gsct=gsct+gsc
         igsct=igsct+1
         bsc=sqrt(1.-1./(gsc*gsc))
         dxp=const2*ext*dz/(bsc*bsc*gsc*gsc*gsc)*abs(f(9,i))
         dyp=const2*eyt*dz/(bsc*bsc*gsc*gsc*gsc)*abs(f(9,i))
         dw=const3*ezt*dz/gsc
         if(.not.iesp) then
c     load the entrance beam parameters for cavities or gaps
           do js=1,7
             f(js,i)=fs(js,i)
           enddo
           f(3,i)=f(3,i)+dxp*1000.
           f(5,i)=f(5,i)+dyp*1000.
           f(2,i)=f(2,i)-dz1*dxp*100.*xpsc
           f(4,i)=f(4,i)-dz1*dyp*100.*xpsc
comment           f(2,i)=f(2,i)-dz1*dxp*100.
comment           f(4,i)=f(4,i)-dz1*dyp*100.
           dwp(i)=dw
         else
           f(3,i)=f(3,i)+dxp*1000.
           f(5,i)=f(5,i)+dyp*1000.
           f(7,i)=f(7,i)+dw
         endif
       enddo
       RETURN
       END
       SUBROUTINE entre
c   ............................................................................
c    ******define the dynamics at input
c   uem : Rest mass in MeV
c
c          proton:938.27231  MeV
c          H_    :939.3145   MeV
c          mesons:33.9093    MeV
c          pions :139.5685   MeV
c          kaons :493.667    MeV
c          electrons : 0.511 MeV
c
c   atm : Atomic number
c   qst : charge
c
c   enedep: Kinetic energy
c   tofini: Time of flight
c
c    REMARK : After INPUT the reference coincides with the c.o.g.
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/speda/dave,idave
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON/DYNI/VREFI,TREFI,FHINIT
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/HISTO/CENTRE(6)
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/MASTRP/XMA(2,2),XMB(2,2),XMC(2,2)
       COMMON/STIS/suryth,surzph,enedep,ecogde,TESTCA
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/tapes/in,ifile,meta
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/TILT/TIPHA,TIX,TIY,SHIFW,SHIFP
       COMMON/newtlt/twissa(3),itwiss
       common/tof/tofini
       common/strip/atm,qs,atms,ths,qop,sqst(6),anp,nqst
       common/mcs/imcs,ncstat,cstat(20)
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       common/trace3e/tracebi(6),traceei(3)
       character*128 trace3h,trace3t,tif
       logical dave
c       dimension foo(9)
        DAVE=.FALSE.
       read(in,*)uem,atm,qst
       if(ncstat.eq.1) cstat(1)=qst
c    input energy(MeV) and initial time of flight(deg)
       read(in,*) enedep,tofini
c 19/03/2009        read(in,*) idave
c 19/03/2009      read(in,*)idum
c 19/03/2009   if icog=0 reference (vref,tref) is not cog
c 19/03/2009   if icog=1 reference (vref,tref) is cog
c --- xmat= rest mass
       xmat=uem*atm
       write(16,101) uem,atm,xmat,qst
101    format('  **** unit mass: ',e12.5,' MeV mass units: ',
     *     f5.1,' rest mass: ',e12.5,' MeV charge ',f4.1)
       write(16,102)enedep,tofini
102    format('  **** energy: ',e12.5,' MeV initial tof: ',e12.5,
     *        ' deg')
       tofini=tofini*pi/(180.*fh)
       gdep=enedep/xmat+1.
       bdep=sqrt(1.-1./(gdep*gdep))
c
       tracebi(6)=tracebi(6)*sqrt(atm)
       write(tif,901) (tracebi(i),i=1,6)
901    format(' BEAMI(1)= ',5(f12.6,1x),f12.4)
       kt3h=kt3h+1
       trace3h(kt3h)=tif
       traceei(3)=traceei(3)/sqrt(atm)
       write(tif,902) (traceei(i),i=1,3)
902    format(' EMITI(1)= ',f12.6,1x,f12.6,1x,f12.4)
       kt3h=kt3h+1
       trace3h(kt3h)=tif
c
c       write(tif,903) uem,int(qst),fhinit/2./pi/1.E6
c903    format(' ER= ',E16.7,', Q= ',I2,', FREQ= ',E14.7)
       write(tif,903) uem*atm/qst,fhinit/2./pi/1.E6
903    format(' ER= ',E16.7,', Q=1 FREQ= ',E14.7)
       kt3h=kt3h+1
       trace3h(kt3h)=tif
       write(tif,7001) enedep/qst
7001   format(' W= ',f12.4)
       kt3h=kt3h+1
       trace3h(kt3h)=tif
c fo(index,i), index=1 initial particle # ,
c              index=2 x
c              index=3 xp
c              index=4 y
c              index=5 yp
c              index=6 time of flight (tof)
c              index=7 energy
c              index=8 if = 0 , then particle is lost
c              index=9 charge
c              i= particle #
       do i=1,ngood
c if itwiss=1,  apply tofini in tiltbm routine
         if(itwiss.ne.1) fo(6,i)=fo(6,i)+tofini
         fo(9,i)=qst
         fo(7,i)=enedep+fo(1,i)+xmat
c  17/08/2009
         if(fo(7,i).lt.xmat) fo(7,i)=xmat
         fo(1,i)=float(i)
       enddo
c --- the reference particle and the cog coincide
       bref=0.
       tref=0.
       encog=0.
       do i=1,ngood
         encog=encog+fo(7,i)
         gai=fo(7,i)/xmat
         bref=bref+sqrt(1.-1/(gai*gai))
         tref=tref+fo(6,i)
       enddo
       encog=encog/float(ngood)
       bref=bref/float(ngood)
       vref=bref*vl
       tref=tref/float(ngood)
       vrefi=vref
       trefi=tref
c       pack the table f(,) of the current beam
       do i=1,ngood
         do j=1,9
           f(j,i)=fo(j,i)
         enddo
       enddo
c   momentum of the  reference (i.e. the c.o.g.)
       gcog=1./sqrt(1.-bref*bref)
       boro=3.3356*xmat*bref*gcog/abs(qst)
       write(16,3450) boro
3450   format('  **** momentum of c.o.g. (kG.cm): ',e12.5)
comment       write(16,*)'before TILT fh, TREF=',fh,tref,tref*fh*180./pi
       if (itwiss.eq.1) then
c Beam was defined in MCOBJET with Twiss parameters
c call TILTBM to apply Twiss alpha
         TIPHA=twissa(3)
         TIX=twissa(1)
         TIY=twissa(2)
         SHIFW=.000000
         SHIFP=.000000
         icg=1
         call tiltbm(icg)
c  write the input beam in file 'input.beam'
         dum=0.
         write(11,*) ngood,dum,fh/(2000000.*pi)
         do i=1,ngood
           f(2,i)=f(2,i)+centre(2)
           f(3,i)=f(3,i)+centre(3)
           f(4,i)=f(4,i)+centre(4)
           f(5,i)=f(5,i)+centre(5)
           f(6,i)=f(6,i)+centre(6)
           f(7,i)=f(7,i)+centre(1)
           etphas=fh*(f(6,i)-tref)
           etener=f(7,i)-xmat
           write(11,777) f(2,i),f(3,i)/1000.,f(4,i),f(5,i)/1000.,
     *                   etphas,etener
         enddo
777      format(6(F13.8,1x))
       else
c  in case of itwiss.ne.1 apply CENTRE and write input beam in TILTBM
         TIPHA=0.
         TIX=0.
         TIY=0.
         SHIFW=.000000
         SHIFP=.000000
         icg=1
         call tiltbm(icg)
       endif
       call emiprt(0)
       return
       end
       SUBROUTINE monte
c  ****************************************************************************************
c  This routine is using the CERN random # routine RLUX (named ranlux in the CERN library)
c
c  ****************************************************************************************
c  ..................................................................................
c   ***********generates randomly the 6D coordinates of the clouds of particles
C   LOI : 1 IMAX particles are generated randomly in three phase plane ellipse
c           with uniform distribution in real space (x,y,z), then xp,yp, and zp
c           from within each phase-plane ellipse. z,zp are converted to phi,w
C   LOI : 2 IMAX particles are generated randomly in a six dimensional ellipsoid
C   LOI : 3 IMAX particles are generated randomly in three phase plane ellipse
c           in real space (x,y,z),with distribution corrresponding to equilibrium
c           stationnary sphere at the limit of the of the current acceptable
c          (see help_DYNAC).Then then xp,yp, and zp from within each phase-plane
c           each phase-plane ellipse. z,zp are converted to phi,w
c   LOI : 4 IMAX particles are generated randomly in a six dimensional ellipsoid from
c           a distribution corrresponding to equilibrium stationnary sphere
c   LOI : 5 IMAX particles are generated randomly in a six dimensional cylinder (axis in z-direction)
c           with uniform distribution in transverse directions
c   LOI : 6 IMAX particles are generated randomly in a six dimensional cylinder (axis in z-direction)
c           with gaussian distribution in transverse directions
c   ITWISS=1   read Twiss parameters for emittance definition
c   ITWISS<>1  reading emittance bounderies for upright ellips
c  .....................................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DCSPA/IESP
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /COM4/CORD(IPTSZ,6)
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/QMOYEN/QMOY
       COMMON/HISTO/CENTRE(6)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/tapes/in,ifile,meta
       common/ranec1/dummy(6)
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/newtlt/twissa(3),itwiss
       COMMON/DYNI/VREFI,TREFI,FHINIT
       common/speda/dave,idave
       common/ragau/ntir
c********************************************
c v28/04/2015
       common/fcont/ifcont
       logical ifcont
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       common/trace3e/tracebi(6),traceei(3)
       character*128 trace3h,trace3t,tif
       logical iesp,chasit,dave
c       dimension fn(9,iptsz),foo(9),so(6),vecx(1)
       dimension vecx(1)
       read(in,*)LOI,itwiss
       WRITE(16,*) ' Generate particles based on law ',loi
       if(loi.gt.6) then
         write(16,*) ' ERROR in GEBEAM: law .gt. 6 invalid ! '
         stop
       endif
       read(in,*)FH,IMAX
       WRITE(16,102) FH
 102   FORMAT(//,30x,' FREQUENCY : ',E12.5,' Hz',//)
       FH=FH*2.*pi
       FHINIT=FH
C   Center of the beam ellipsoid
C    X X' Y Y' Z' Z
       read(in,*)(CENTRE(J),J=2,5),CENTRE(1),CENTRE(6)
       if (itwiss.ne.1) then
         read(in,*)YMAX,TMAX,ZMAX,PMAX,DMAX,TTMAX
         kt3h=kt3h+1
         trace3h(kt3h)='ERROR: GEBEAM ITWISS.NE.1 not yet implemented'
       else
         read(in,*) alphax,betax,emitx
         read(in,*) alphay,betay,emity
         read(in,*) alphaz,betaz,emitz
         gammax=(1.+alphax*alphax)/betax
         gammay=(1.+alphay*alphay)/betay
         ymax=0.1*sqrt(emitx/gammax)
         tmax=sqrt(emitx*gammax)
         zmax=0.1*sqrt(emity/gammay)
         pmax=sqrt(emity*gammay)
         twissa(1)=-alphax*ymax
         twissa(2)=-alphay*zmax
         if(loi.ne.5) then
           gammaz=(1.+alphaz*alphaz)/betaz
           dmax=0.001*sqrt(emitz*gammaz)
           ttmax=pi*sqrt(emitz/gammaz)/(fh*180.)
           twissa(3)=alphaz*sqrt(emitz/gammaz)
         endif
         if(loi.ge.5) then
           dmax=alphaz
           twissa(3)=0.
           ttmax=pi/fh
         endif
         if(loi.eq.6) then
           sig=betaz
         endif
c store parameters for trace3d, write them in ENTRE
         tracebi(1)=alphax
         tracebi(2)=betax
         tracebi(3)=alphay
         tracebi(4)=betay
         tracebi(5)=alphaz
         tracebi(6)=betaz
         traceei(1)=emitx
         traceei(2)=emity
         traceei(3)=emitz
       endif
c --- beam centre is activated in SUBROUTINE entre
       cph=centre(6)*fh*180/pi
       WRITE(16,123) (CENTRE(J),J=2,5),CENTRE(1),CENTRE(6),cph
123    FORMAT(' *** Beam centre defined as:',/,3X,
     X  4X,' Transverse direction :',/,6X,
     X' HORZ PLANE X(CM) = ',E12.5,'       XP(MRD) = ',E12.5,/,6X,
     X' VERT PLANE Y(CM) = ',E12.5,'       YP(MRD) = ',E12.5,/,4X,
     X' LONGITUDINAL :',/,5X,
     X' DELTA ENERGY(MeV) = ',E12.5,'     TIME(SEC) = ',E12.5,/,
     X  41x,' PHASE(DEG) = ',E12.5,//)
       PTMAX=ttmax*fh*180./pi
      if(itwiss.ne.1) then
        WRITE(16,99)YMAX,TMAX,ZMAX,PMAX,DMAX,TTMAX,PTMAX
99      FORMAT(3X,' *** Limits of the random distribution ',/,
     X  4X,' Transverse direction :',/,6X,
     X' HORZ PLANE X(CM) = ',E12.5,'       XP(MRD) = ',E12.5,/,6X,
     X' VERT PLANE Y(CM) = ',E12.5,'       YP(MRD) = ',E12.5,/,4X,
     X' LONGITUDINAL :',/,5X,
     X' DELTA ENERGY(MeV) = ',E12.5,'     TIME(SEC) = ',E12.5,/,
     X  41x,' PHASE(DEG) = ',E12.5,//)
      else
c        if(loi.eq.5) then
        if(loi.ge.5) then
          write(16,*) ' Beam distribution based on Twiss parameters'
          write(16,199) alphax,betax,emitx,alphay,betay,emity
          write(16,399) ptmax,dmax
399       format(4x,' Continuous beam in the longitudinal direction :',
     *    /,6x,' half phase length (deg): ',e12.5,/,6x,
     *    ' half energy width (MeV): ',e12.5)
c continuous beam
          ifcont=.true.
        else
          write(16,*) ' Beam distribution based on Twiss parameters'
          write(16,199) alphax,betax,emitx,alphay,betay,emity
199       format(4x,' Transverse direction :',/,6x,
     *    ' Horz plane: alpha: ',e12.5,' beta(mm/mrad): ',e12.5,
     *    ' emit(pi*mm*mrad): ',e12.5,/,6x,
     *    ' Vert plane: alpha: ',e12.5,' beta(mm/mrad): ',e12.5,
     *    ' emit(pi*mm*mrad): ',e12.5)
          write(16,299) alphaz,betaz,emitz
299       format(4x,' Longitudinal direction :',/,6x,
     *    ' alpha: ',e12.5,' beta(deg/keV): ',e12.5,
     *    ' emit(pi*deg*keV): ',e12.5)
        endif
      endif
C     len : starting value of the random vector for the routine rlux
       len=1
       do j=1,6
         cord(1,j)=0.0
       enddo
       if(loi.eq.1) then
         do i=2,imax
150        continue
           call rlux(vecx,len)
           r1=2.*vecx(1)-1.
           call rlux(vecx,len)
           r3=2.*vecx(1)-1.
           call rlux(vecx,len)
           r6=2.*vecx(1)-1.
c  make round in x,y,z plane
           rho=r1**2+r3**2+r6**2
           if (rho.gt.1) go to 150
152        call rlux(vecx,len)
           r2=2.*vecx(1)-1.
c make beam round in x,x' plane
           if((r1*r1+r2*r2).gt.1.) go to 152
153        call rlux(vecx,len)
           r4=2.*vecx(1)-1.
c  make beam round in y,yp plane
           if((r4*r4+r3*r3).gt.1.) go to 153
180        call rlux(vecx,len)
           r5=2.*vecx(1)-1.
c   make beam round in plane z, zp
           if((r6*r6+r5*r5).gt.1.) go to 180
c        store random numbers in preparation for rms correction
           cord(i,1)=r1
           cord(i,2)=r2
           cord(i,3)=r3
           cord(i,4)=r4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
       endif
       if(loi.eq.2) then
         do i=2,imax
1500       call rlux(vecx,len)
           r1=2.*vecx(1)-1.
           call rlux(vecx,len)
           r2=2.*vecx(1)-1.
c   make round x,xp
           rho=r1*r1+r2*r2
           if(rho.gt.1) go to 1500
1530       call rlux(vecx,len)
           r3=2.*vecx(1)-1.
           call rlux(vecx,len)
           r4=2.*vecx(1)-1.
c  make beam round in y,yp plane
           if((r4*r4+r3*r3).gt.1.) go to 1530
1800       call rlux(vecx,len)
           r5=2.*vecx(1)-1.
           call rlux(vecx,len)
           r6=2.*vecx(1)-1.
c   make beam round in plane z, zp
           if((r6*r6+r5*r5).gt.1.) go to 1800
c   make round x,xp,y,yp,z,zp
           rho=r1**2+r2**2+r3**2+r4**2+r5**2+r6**2
           if (rho.gt.1.0) go to 1500
c        store random numbers in preparation for rms correction
           cord(i,1)=r1
           cord(i,2)=r2
           cord(i,3)=r3
           cord(i,4)=r4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
       endif
       if(loi.eq.3) then
         ntir=24
         s=.2493
         am=0.
         do i=2,imax
14         call gcern (len,s,am,vec)
           r1=vec
           call gcern (len,s,am,vec)
           r3=vec
           call gcern (len,s,am,vec)
           r6=vec
c  make round in x,y,z plane
           rho=r1**2+r3**2+r6**2
           if (rho.gt.1) go to 14
16         call gcern (len,s,am,vec)
           r2=vec
c make beam round in x,x' plane
           if((r1*r1+r2*r2).gt.1.) go to 16
17         call gcern (len,s,am,vec)
           r4=vec
c  make beam round in y,yp plane
           if((r4*r4+r3*r3).gt.1.) go to 17
22         call gcern (len,s,am,vec)
           r5=vec
c   make beam round in plane z, zp
           if((r6*r6+r5*r5).gt.1.) go to 22
c        store random numbers in preparation for rms correction
           cord(i,1)=r1
           cord(i,2)=r2
           cord(i,3)=r3
           cord(i,4)=r4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
       endif
       if(loi.eq.4) then
         ntir=12
         s=.2493
         am=0.
         do i=2,imax
101        call gcern (len,s,am,vec)
           r1=vec
           call gcern (len,s,am,vec)
           r2=vec
c make beam round in x,xp plane
           rho=r1*r1+r2*r2
           if (rho.gt.1) go to 101
112        call gcern (len,s,am,vec)
           r3=vec
           call gcern (len,s,am,vec)
           r4=vec
c make beam round in y,yp plane
           if((r3*r3+r4*r4).gt.1.) go to 112
113        call gcern (len,s,am,vec)
           r5=vec
           call gcern (len,s,am,vec)
           r6=vec
c   make beam round in plane z, zp
           if((r6*r6+r5*r5).gt.1.) go to 113
c        store random numbers in preparation for rms correction
c   make round x,xp,y,yp,z,zp
           rho=r1**2+r2**2+r3**2+r4**2+r5**2+r6**2
           if (rho.gt.1.0) go to 101
           cord(i,1)=r1
           cord(i,2)=r2
           cord(i,3)=r3
           cord(i,4)=r4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
       endif
       if(loi.eq.5) then
         do i=2,imax
1566       continue
           call rlux(vecx,len)
           r1=2.*vecx(1)-1.
           call rlux(vecx,len)
           r3=2.*vecx(1)-1.
c  make round in x,y plane
           rho=r1**2+r3**2
           if (rho.gt.1) go to 1566
1525        call rlux(vecx,len)
           r2=2.*vecx(1)-1.
c make beam round in x,x' plane
           if((r1*r1+r2*r2).gt.1.) go to 1525
1535        call rlux(vecx,len)
           r4=2.*vecx(1)-1.
c  make beam round in y,yp plane
           if((r4*r4+r3*r3).gt.1.) go to 1535
c   make round x,xp,y,yp
           rho=r1**2+r2**2+r3**2+r4**2
           if (rho.gt.1.0) go to 1566
c*et*2010-11-23  do NOT make beam round in z,zp plane
cet1835       call rlux(vecx,len)
           call rlux(vecx,len)
           r5=2.*vecx(1)-1.
           call rlux(vecx,len)
           r6=2.*vecx(1)-1.
cet           if((r6*r6+r5*r5).gt.1.) go to 1835
c        store random numbers in preparation for rms correction
           cord(i,1)=r1
           cord(i,2)=r2
           cord(i,3)=r3
           cord(i,4)=r4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
       endif
       if(loi.eq.6) then
c Gaussian
c*et*2014-Mar-09
         do i=2,imax
           CALL RGAUS2(SIG,Y1,Y2,Y3,Y4)
c do NOT make beam gaussien in z,zp plane
           call rlux(vecx,len)
           r5=2.*vecx(1)-1.
           call rlux(vecx,len)
           r6=2.*vecx(1)-1.
c
c        store random numbers in preparation for rms correction
           cord(i,1)=y1
           cord(i,2)=y2
           cord(i,3)=y3
           cord(i,4)=y4
           cord(i,5)=r5
           cord(i,6)=r6
         enddo
         y1x=abs(cord(2,1))
         y2x=abs(cord(2,2))
         y3x=abs(cord(2,3))
         y4x=abs(cord(2,4))
         do i=3,imax
           if(abs(cord(i,1)).gt.y1x)y1x=abs(cord(i,1))
           if(abs(cord(i,2)).gt.y2x)y2x=abs(cord(i,2))
           if(abs(cord(i,3)).gt.y3x)y3x=abs(cord(i,3))
           if(abs(cord(i,4)).gt.y4x)y4x=abs(cord(i,4))
         enddo
         do i=2,imax
           cord(i,1)=cord(i,1)/y1x
           cord(i,2)=cord(i,2)/y2x
           cord(i,3)=cord(i,3)/y3x
           cord(i,4)=cord(i,4)/y4x
         enddo
       endif
       call corre(imax,imax)
c       fimax=ttmax*fh*180./pi
c maximum extent in case of continous beam is +/-pi (i.e. +/-180 deg)
       tcorct=abs(.5*cord(2,6)*ttmax)
c   in fo(1,) is stored the energy extent (MeV)
       do i=2,imax
         fo(1,i)=.5*cord(i,5)*dmax
         fo(2,i)=.5*cord(i,1)*ymax
         fo(3,i)=.5*cord(i,2)*tmax
         fo(4,i)=.5*cord(i,3)*zmax
         fo(5,i)=.5*cord(i,4)*pmax
         fo(6,i)=.5*cord(i,6)*ttmax
         if(abs(fo(6,i)).gt.tcorct) tcorct=abs(fo(6,i))
       enddo
c    first particle
cccc        imax=imax+1
        ichas(1)=1
        fo(8,1)=1.
        do j=1,6
          fo(j,1)=0.
        enddo
c        if(loi.ne.5) then
        if(loi.lt.5) then
          do i=2,imax
            fo(8,i)=1
            ICHAS(I)=1
          enddo
        else
c correct the phase length for a continuous beam (force +/- 180 deg)
          tcorct=ttmax/tcorct
          do i=2,imax
            fo(8,i)=1
            ICHAS(I)=1
            fo(6,i)=fo(6,i)*tcorct
          enddo
        endif
        ngood=imax
        imaxo=imax
        WRITE(16,8)LOI,IMAX
8       FORMAT(8X ,'  ****law ',i2,'  with ',i6,' particles',/)
        RETURN
        END
        subroutine rgaus2(sigma,y1,y2,y3,y4)
        implicit real*8 (a-h,o-z)
        dimension vecx(1)
        len=1
        do while ((w1.ge.1.0).or.(w1.eq.0.))
          call rlux(vecx,len)
          x1 = 2.0 * vecx(1) - 1.0
          call rlux(vecx,len)
          x3 = 2.0 * vecx(1) - 1.0
          w1 = x1 * x1 + x3 * x3
        enddo
        do while ((w2.ge.1.0).or.(w2.eq.0.))
          call rlux(vecx,len)
          x2 = 2.0 * vecx(1) - 1.0
          call rlux(vecx,len)
          x4 = 2.0 * vecx(1) - 1.0
          w2 = x2 * x2 + x4 * x4
        enddo
        w1 = sigma*sqrt( (-2.0 * log( w1 ) ) / w1 )
        w2 = sigma*sqrt( (-2.0 * log( w2 ) ) / w2 )
        y1 = x1 * w1
        y2 = x2 * w2
        y3 = x3 * w1
        y4 = x4 * w2
        return
        end
       SUBROUTINE rlux(RVEC,LENV)
       implicit real*8 (a-h,o-z)
C* ranlux.f Rev 1.2  1997/09/22 13:45:47  mclareni
C* Correct error in initializing RANLUX by using RLUXIN with the output of
C* RLUXUT from a previous run.
C* CERN Mathlib gen
C*
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C         Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       DIMENSION RVEC(LENV)
       DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
       PARAMETER (MAXLEV=4, LXDFLT=3)
       DIMENSION NDSKIP(0:MAXLEV)
       DIMENSION NEXT(24)
       PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
       PARAMETER (ITWO24=2**24, ICONS=2147483563)
       SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
       SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
       INTEGER LUXLEV
       LOGICAL NOTYET
       DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
       DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
       DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
       IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT
         INSEED = JSEED
         WRITE(16,'(A,I12)')' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(16,'(A,I2,A,I4)')' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
         TWOM24 = 1.
         DO I= 1, 24
           TWOM24 = TWOM24 * 0.5
           K = JSEED/53668
           JSEED = 40014*(JSEED-K*53668) -K*12211
           IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
           ISEEDS(I) = MOD(JSEED,ITWO24)
         ENDDO
         TWOM12 = TWOM24 * 4096.
         DO I= 1,24
           SEEDS(I) = REAL(ISEEDS(I))*TWOM24
           NEXT(I) = I-1
         ENDDO
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
       ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
       DO IVEC= 1, LENV
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
           UNI = UNI + 1.0
           CARRY = TWOM24
         ELSE
           CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
         RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
         IF (UNI .LT. TWOM12)  THEN
           RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C          and zero is forbidden in case someone takes a logarithm
           IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
         ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
         IN24 = IN24 + 1
         IF (IN24 .EQ. 24)  THEN
           IN24 = 0
           KOUNT = KOUNT + NSKIP
           DO ISK= 1, NSKIP
             UNI = SEEDS(J24) - SEEDS(I24) - CARRY
             IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
             ELSE
               CARRY = 0.
             ENDIF
             SEEDS(I24) = UNI
             I24 = NEXT(I24)
             J24 = NEXT(J24)
           ENDDO
         ENDIF
       ENDDO
       KOUNT = KOUNT + LENV
       IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
       ENDIF
       RETURN
C
C           Entry to input and float integer seeds from previous run
       ENTRY RLUXIN(ISDEXT)
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO I= 1, 24
           NEXT(I) = I-1
           TWOM24 = TWOM24 * 0.5
         ENDDO
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
         WRITE(16,'(A)')'FULL INITIALIZATION OF RANLUX WITH 25',
     *                  ' INTEGERS'
         WRITE(16,'(5X,5I12)') ISDEXT
         DO I= 1, 24
           SEEDS(I) = REAL(ISDEXT(I))*TWOM24
         ENDDO
         CARRY = 0.
         IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
         ISD = IABS(ISDEXT(25))
         I24 = MOD(ISD,100)
         ISD = ISD/100
         J24 = MOD(ISD,100)
         ISD = ISD/100
         IN24 = MOD(ISD,100)
         ISD = ISD/100
         LUXLEV = ISD
         IF (LUXLEV .LE. MAXLEV) THEN
           NSKIP = NDSKIP(LUXLEV)
           WRITE (6,'(A,I2)')'RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
         ELSE  IF (LUXLEV .GE. 24) THEN
           NSKIP = LUXLEV - 24
           WRITE (6,'(A,I5)')'RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
         ELSE
           NSKIP = NDSKIP(MAXLEV)
           WRITE (6,'(A,I5)')'RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
           LUXLEV = MAXLEV
         ENDIF
         INSEED = -1
       RETURN
C
C                    Entry to output seeds as integers
       ENTRY RLUXUT(ISDEXT)
         DO I= 1, 24
           ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
         ENDDO
         ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
         IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
       RETURN
C
C                    Entry to output the "convenient" restart point
       ENTRY RLUXAT(LOUT,INOUT,K1,K2)
         LOUT = LUXLEV
         INOUT = INSEED
         K1 = KOUNT
         K2 = MKOUNT
       RETURN
C
C                    Entry to initialize from one or three integers
       ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
           LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
           LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
           LUXLEV = MAXLEV
           WRITE (6,'(A,I7)')'RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
           LUXLEV = LUX
           DO ILX= 0, MAXLEV
             IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
           ENDDO
         ENDIF
         IF (LUXLEV .LE. MAXLEV)  THEN
           NSKIP = NDSKIP(LUXLEV)
           WRITE(16,'(A,I2,A,I4)')'RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
         ELSE
           NSKIP = LUXLEV - 24
           WRITE(16,'(A,I5)')'RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
         ENDIF
         IN24 = 0
         IF (INS .LT. 0)  WRITE (6,'(A)')
     +   ' Illegal initialization by RLUXGO, negative input seed'
         IF (INS .GT. 0)  THEN
           JSEED = INS
           WRITE(16,'(A,3I12)')'RANLUX INITIALIZED BY ',
     +     'RLUXGO FROM SEEDS',JSEED, K1,K2
         ELSE
           JSEED = JSDFLT
           WRITE(16,'(A)')'RANLUX INITIALIZED BY RLUXGO FROM DEFAULT',
     +                    ' SEED'
         ENDIF
         INSEED = JSEED
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO I= 1, 24
           TWOM24 = TWOM24 * 0.5
           K = JSEED/53668
           JSEED = 40014*(JSEED-K*53668) -K*12211
           IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
           ISEEDS(I) = MOD(JSEED,ITWO24)
         ENDDO
         TWOM12 = TWOM24 * 4096.
         DO I= 1,24
           SEEDS(I) = REAL(ISEEDS(I))*TWOM24
           NEXT(I) = I-1
         ENDDO
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
         KOUNT = K1
         MKOUNT = K2
         IF (K1+K2 .NE. 0)  THEN
           DO IOUTER= 1, K2+1
             INNER = IGIGA
             IF (IOUTER .EQ. K2+1)  INNER = K1
             DO ISK= 1, INNER
               UNI = SEEDS(J24) - SEEDS(I24) - CARRY
               IF (UNI .LT. 0.)  THEN
                 UNI = UNI + 1.0
                 CARRY = TWOM24
               ELSE
                 CARRY = 0.
               ENDIF
               SEEDS(I24) = UNI
               I24 = NEXT(I24)
               J24 = NEXT(J24)
             ENDDO
           ENDDO
C         Get the right value of IN24 by direct calculation
           IN24 = MOD(KOUNT, NSKIP+24)
           IF (MKOUNT .GT. 0)  THEN
             IZIP = MOD(IGIGA, NSKIP+24)
             IZIP2 = MKOUNT*IZIP + IN24
             IN24 = MOD(IZIP2, NSKIP+24)
           ENDIF
C       Now IN24 had better be between zero and 23 inclusive
           IF (IN24 .GT. 23) THEN
             WRITE (6,'(A/A,3I11,A,I5)')
     +    '  Error in RESTARTING with RLUXGO:',' The values',INS,
     +       K1, K2, ' cannot occur at luxury level', LUXLEV
             IN24 = 0
           ENDIF
         ENDIF
       RETURN
       END
       SUBROUTINE gcern(len,s,am,v)
c generateur aleatoire selon une loi normale
c        s : ecart-type de la distribution
c        am: moyenne de ladistribution
c        v : nombre aléatoire selon la loi normale
       implicit real*8 (a-h,o-z)
       common/ragau/ntir
       dimension vecx(1)
       a=0.
       do i=1,ntir
         call rlux(vecx,len)
         y=vecx(1)
         a=a+y
       enddo
       v=(a-float(ntir)/2.)*s+am
       return
       end
       SUBROUTINE adjrfq
c   ..............................................................................
c    read the coordinates of particles from file
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/speda/dave,idave
       COMMON/DYN/TREF,VREF
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/tapes/in,ifile,meta
       common/isxpyp/iflag
       common/mcs/imcs,ncstat,cstat(20)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYNI/VREFI,TREFI,FHINIT
       common/strip/atm,qs,atms,ths,qop,sqst(6),anp,nqst
       logical dave
       read(in,*) iflag
       if(iflag.eq.0 .or. iflag.eq.100)
     *   write(16,*)'Standard file, phase in rad'
       if(iflag.eq.1 .or. iflag.eq.101)
     *   write(16,*)'File with charge state and rest ',
     *               'mass, phase in rad'
       if(iflag.eq.2 .or. iflag.eq.102)
     *   write(16,*)'File with several charge states',
     *               ', phase in rad'
       if(iflag.eq.10 .or. iflag.eq.110)
     *   write(16,*)'Standard file, phase in ns'
       if(iflag.eq.11 .or. iflag.eq.111)
     *   write(16,*)'File with charge state and rest',
     *               ' mass, phase in ns'
       if(iflag.eq.12 .or. iflag.eq.112)
     *   write(16,*)'File with several charge states',
     *               ', phase in ns'
c ---  input: freq.(MHertz) tofini: phase offset (deg) to be applied both to the reference and the beam
       read(in,*) freq,tofini
       write(16,*) 'Frequency [MHz]:',freq
       fh=2.*pi*freq*1.e06
       FHINIT=FH
       read(in,*)uem,atm
       xmat=uem*atm
c --- reference:
c   enedep:energy(MeV), qst: charge of the reference
c   when iflag =0 or iflag = 1 qst is the charge of the beam
       read(in,*) enedep,qst
       tofini=tofini*pi/(180.*fh)
       gdep=enedep/xmat+1.
       bdep=sqrt(1.-1./(gdep*gdep))
       vref=bdep*vl
       tref=tofini
       ncstat=1
c*et*2012*March*02 use reference to define boro
       boro=3.3356*xmat*bdep*gdep/qst
c       write(16,111) uem,atm,enedep,qst,tofini
c111    format(' **** mass unit: ',e12.5,'MeV, mass unit: ',f6.1,/,
c     *      ' **** reference : energy ',e12.5,' Mev  charge ',f4.1,
c     *      ' time of flight ',e12.5,' sec')
       write(16,101) uem,atm,qst,tofini,enedep,boro
101    format(' **** unit mass: ',e12.5,' MeV, mass units: ',f6.1,/,
     *      ' **** reference charge ',f4.1,' time of flight ',e12.5,
     *      ' sec',/,' **** reference : energy ',e12.5,
     *      ' MeV  momentum ',e12.5,' kG.cm')

       read(55,*) imax,dum,dum
       if(imax+2.gt.iptsz) then
         WRITE(16,*) 'too many particles '
         stop
       endif
c  iflag = 0 standard file:f(1,)=x, f(2,)=xp, f(3,)=y, f(4,)=yp, f(5,)=phase, f(6,)=kinetic energy
        if(iflag.eq.0 .or. iflag.eq.10) read(55,*)
     *   ((f(i,j),i=1,6),j=1,imax)
c  iflag = 1  File with rest mass: f(1,)=x, f(2,)=xp, f(3,)=y, f(4,)=yp, f(5,)=phase, f(6,)=kinetic energy, dum1, dum2
        if(iflag.eq.1 .or. iflag.eq.11) read(55,*)
     *   ((f(i,j),i=1,6),dum1,dum2,j=1,imax)
c  iflag = 2 beam with different charges: f(7,) = charge; figure out how many different ones there are and store them
c*ets*13/12/2010
       if(ncstat.eq.1) cstat(1)=qst
       if(iflag.eq.2 .or. iflag.eq.12) then
         ncstat=1
         read(55,*) (f(i,1),i=1,7)
         cstat(1)=f(7,1)
         do j=2,imax
           read(55,*) (f(i,j),i=1,7)
           mcstat=0
           do k=1,ncstat
c*ets*23/12/2010
             if(f(7,j).eq.cstat(k)) then
               mcstat=1
             endif
           enddo
           if(mcstat.eq.0) then
c*ete*23/12/2010
             ncstat=ncstat+1
             cstat(ncstat)=f(7,j)
           endif
         enddo
         write(16,*) 'Number of charge states: ',ncstat
         write(16,*) 'Charge states: ',(cstat(j),j=1,ncstat)
         if(ncstat.gt.1) imcs=1
       endif
       call intfac(tofini)
       RETURN
       END
       SUBROUTINE intfac(tofini)
c   ............................................................................
c  convert the particles  in DYNAC units
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/tapes/in,ifile,meta
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/speda/dave,idave
       common/isxpyp/iflag
       common/mcs/imcs,ncstat,cstat(20)
       dimension foo(20,9),NDP(20)
       logical shift,chasit,dave
C      FH :INITIAL FREQENCY (hertz)
c Use frequency from part. dist. file to calculate f(6,i)
c 19/03/2009        dave=.true.
c    iflag = 0 : standard file, phase in rad
c    iflag = 1 : special file, phase in rad
c    iflag = 2 : several charge states, phase in rad
c    iflag = 10 : standard file, phase in ns
c    iflag = 11 : special file, phase in ns
c    iflag = 12 : several charge states, phase in ns
       imcs=0
       if(iflag.eq.2 .or. iflag.eq.12)  imcs=1
       QMOY=QST
       j=1
c   the table fo(,) is built from the input beam stored in the table f(,)
       if(iflag.le.2) then
         do i=1,imax
           fo(7,i)=f(6,j)+xmat
           if(iflag.le.1) fo(9,i)=qst
           if(iflag.eq.2) fo(9,i)=f(7,i)
           ichas(i)=1
           fo(8,i)=1.
           fo(1,i)=float(j+1)
           fo(2,i)=f(1,j)
             fo(3,i)=f(2,j)*1000.
           fo(4,i)=f(3,j)
           fo(5,i)=f(4,j)*1000.
           fo(6,i)=tofini+f(5,j)/fh
           j=j+1
         enddo
       else
         do i=1,imax
           fo(7,i)=f(6,j)+xmat
           if(iflag.le.11) fo(9,i)=qst
           if(iflag.eq.12) fo(9,i)=f(7,i)
           ichas(i)=1
           fo(8,i)=1.
           fo(1,i)=float(j+1)
           fo(2,i)=f(1,j)
           fo(3,i)=f(2,j)*1000.
           fo(4,i)=f(3,j)
           fo(5,i)=f(4,j)*1000.
           fo(6,i)=tofini+f(5,j)*1.e-09
           j=j+1
         enddo
       endif
c   cog of the beam
c*et*2012*March*02 changed foo(9) to foo(2,9)
       if(iflag.eq.0.or.iflag.eq.1.or.iflag.eq.10.or.iflag.eq.11)then
c single charge state
         NDP(1)=imax
         do j=2,7
           foo(1,j)=0.
         enddo
         do i=1,imax
           do j=2,7
             foo(1,j)=foo(1,j)+fo(j,i)
           enddo
         enddo
       else
c multi charge state
         DO k=1,ncstat
           NDP(k)=0
           do j=2,7
             foo(k,j)=0.
           enddo
         ENDDO
         DO i=1,imax
           DO k=1,ncstat
             if(fo(9,i).eq.cstat(k)) then
               NDP(k)=NDP(k)+1
               do j=2,7
                 foo(k,j)=foo(k,j)+fo(j,i)
               enddo
             endif
           ENDDO
         ENDDO
       endif
       DO k=1,ncstat
         do j=2,7
           foo(k,j)=foo(k,j)/float(NDP(k))
         enddo
       ENDDO
       ngood=imax
       if(iflag.eq.0.or.iflag.eq.1.or.iflag.eq.10.or.iflag.eq.11)then
c --- COG
         gref=foo(1,7)/xmat
         bref=sqrt(1.-1./(gref*gref))
         xe=(gref-1.)*xmat
c   magnetic rigidity
         bor=3.3356*xmat*bref*gref/qst
         write(16,*) '**** COG : energy ',xe,' MeV  momentum '
     *  ,boro,' kG.cm'
       else
         DO k=1,ncstat
           gref=foo(k,7)/xmat
           bref=sqrt(1.-1./(gref*gref))
           xe=(gref-1.)*xmat
c   magnetic rigidity
           bor=3.3356*xmat*bref*gref/cstat(k)
           write(16,*) ' Q: ',cstat(k),' COG : energy ',xe,
     *     ' MeV  momentum ',bor,' kG.cm'
         ENDDO
       endif
c  add the last particle equal to the center of gravity
c  cor 12/1/09
c*et*2012*March*02 suppress the imax+1 particle
c       ngood=imax+1
c       fo(1,ngood)=float(ngood)
c       fo(8,ngood)=1
c       fo(9,ngood)=qst
c       ichas(ngood)=1
c       do j=2,7
c         fo(j,ngood)=foo(j)/float(imax)
c       enddo
c       do j=2,7
c         fo(j,ngood)=foo(j)/float(imax)
c       enddo
c       imax=imax+1
       vrefi=vref
       trefi=tref
c  TEST sv
c   first particle forced to x=xp=y=yp=0,tref,vref
ccc       fo(1,1)=1.
ccc       ichas(1)=1
ccc       do j=2,5
ccc         fo(j,1)=0.
ccc       enddo
ccc       fo(6,1)=tref
ccc       fo(7,1)=enedep+xmat
c **************************************************
c now save data back to f
       do i=1,ngood
         do k=1,9
           f(k,i)=fo(k,i)
         enddo
       enddo
c9999   format(6(2x,e12.5))
       imaxo=ngood
       call emiprt(0)
       RETURN
       END
       SUBROUTINE stapl(zpos)
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/pltprf/sprfx(3000),sprfy(3000),sprfl(3000),sprfw(3000),
     *               sprfp(3000),SPRNG(3000),iprf
       common/pltprf1/sprww(3000),eprfw(3000),eprnx(3000),eprny(3000),
     *                sprfz(3000)
       common/pltprf2/sxmn(3000),sxmx(3000),symn(3000),symx(3000),
     *                stmn(3000),stmx(3000),spmn(3000),spmx(3000),
     *                swmn(3000),swmx(3000),disprx(3000),dispry(3000),
     *                dispcx(3000),dispcy(3000)
       common/etcom/cog(8),exten(17),fd(iptsz)
       common /consta/ vl, pi, xmat, rpel, qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/faisc/f(10,iptsz),imax,ngood
       common/grot/rzot,izrot
       logical izrot
c       *************************************************************************************
c       the statistics in EXT2
c
c       sprfx(cm): 2.*sqrt( sum(x*x) )
c       sprfy(cm): 2.*sqrt( sum(y*y) )
c       sprfw:     2.*sqrt( sum(dp/p * dp/p) )*beta*beta = (energy spread)/(energy of c.o.g)
c       sprfp(deg):2.*sqrt( sum(dphi * dphi) )
c       sprfz(cm) :   sqrt( sum(dt * dt) )
c       sprfl(m) : t.o.f of the reference
c       sprww(MeV):kinetic energy of the c.o.g
c       eprfw(ns.keV) : longitudinal emittance
c       eprnx(mm.mrad): normalized emittance in x-direction
c       eprny(mm.mrad): normalized emittance in y-direction
c       *************************************************************************************
c       iprf: pointer
c cor 08/01/2009
comment       iprf=iprf+1
c  izrot: logical flag set as .true. in the routine ZROTA
       if(izrot) call zrotap(-rzot)
       iarg=0
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       CALL EXT2(iarg)
       qdisp=2.*sqrt(exten(1))
       qmd=exten(1)*exten(3)-exten(2)*exten(2)
       delw=encog*encog*bcog**4
       qmdv=qmd*delw
       sqmdv=4.*pi*sqrt(qmdv)
       eprfw(iprf)=sqmdv*1.e12/(pi*fh)
       sprfx(iprf)=2.*sqrt(exten(4))
       sprfy(iprf)=2.*sqrt(exten(6))
       trqtx=exten(4)*exten(5)-exten(8)*exten(8)
       trqpy=exten(6)*exten(7)-exten(9)*exten(9)
       surxth=4.*pi*sqrt(trqtx)
       suryph=4.*pi*sqrt(trqpy)
       eprnx(iprf)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
       eprny(iprf)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
c horizontal and vertical dispersion in meter (w.r.t. reference)
c       disprx(iprf)=0.01*exten(12)/exten(14)
c       dispry(iprf)=0.01*exten(13)/exten(14)
c horizontal and vertical dispersion in meter (w.r.t. COG)
       dispcx(iprf)=0.01*exten(15)/exten(17)
       dispcy(iprf)=0.01*exten(16)/exten(17)
       trqfi=0.
       tof=0.
       chxmx=f(2,1)
       chxmn=f(2,1)
       chymx=f(4,1)
       chymn=f(4,1)
       chpmx=f(6,1)-tcog
       chpmn=f(6,1)-tcog
       chwmx=f(7,1)-encog
       chwmn=f(7,1)-encog
       do i=1,ngood
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         fdp=(gpai*bpai)/(gcog*bcog)-1.
         trqfi=trqfi+fdp*fdp
         tf=(tcog-f(6,i))*bpai*vl
         tof=tof+tf*tf
         if(f(2,i).gt.chxmx)chxmx=f(2,i)
         if(f(2,i).lt.chxmn)chxmn=f(2,i)
         if(f(4,i).gt.chymx)chymx=f(4,i)
         if(f(4,i).lt.chymn)chymn=f(4,i)
         if(f(6,i)-tcog.gt.chpmx) then
           chpmx=f(6,i)-tcog
           chdmx=chpmx*fh*180./pi
         endif
         if(f(6,i)-tcog.lt.chpmn) then
           chpmn=f(6,i)-tcog
           chdmn=chpmn*fh*180./pi
         endif
         if(f(7,i)-encog.gt.chwmx)chwmx=f(7,i)-encog
         if(f(7,i)-encog.lt.chwmn)chwmn=f(7,i)-encog
       enddo
       trqfi=trqfi/float(ngood)
       tof=tof/float(ngood)
       sprfz(iprf)=sqrt(tof)
       cmult=(gcog+1.)/gcog
       sprfw(iprf)=2.*sqrt(trqfi)*cmult
c v29/06/2010       sprfw(iprf)=2.*sqrt(trqfi)*bcog*bcog
       sprfp(iprf)=2.*sqrt(exten(3))*180./pi
       sprww(iprf)=cog(1)-xmat
c   t.o.f of the reference in m
       sprfl(iprf)=zpos/1000.
c*et*2013-Aug-21 set up ngood particle counter for dynac.print file
       sprng(iprf)=ngood
c*et*2014-Sep-25 Add recording of min. and max. positions in x,y and z,z'
       sxmx(iprf)=chxmx*10.
       sxmn(iprf)=chxmn*10.
       symx(iprf)=chymx*10.
       symn(iprf)=chymn*10.
       stmx(iprf)=chpmx
       stmn(iprf)=chpmn
       spmx(iprf)=chdmx
       spmn(iprf)=chdmn
       swmx(iprf)=chwmx
       swmn(iprf)=chwmn
c evolve pointer iprf (cor. 08/01/2009)
       iprf=iprf+1
       if(izrot) call zrotap(-rzot)
       return
       end
       SUBROUTINE eugwrt
       implicit real*8 (a-h,o-z)
       common/pltprf/sprfx(3000),sprfy(3000),sprfl(3000),sprfw(3000),
     *               sprfp(3000),SPRNG(3000),iprf
       common/pltprf1/sprww(3000),eprfw(3000),eprnx(3000),eprny(3000),
     *                sprfz(3000)
       common/pltprf2/sxmn(3000),sxmx(3000),symn(3000),symx(3000),
     *                stmn(3000),stmx(3000),spmn(3000),spmx(3000),
     *                swmn(3000),swmx(3000),disprx(3000),dispry(3000),
     *                dispcx(3000),dispcy(3000)
c l(m) x(mm) y(mm) z(deg) z(mm) Ex(mm.mrd) Ey(mm.mrd) Ez(KeV.ns) Wcog(MeV) xmin(mm) xmax(mm) ymin(mm) ymax(mm)
c  tmin(ns) tmax(ns) pmin(deg) pmax(deg) wmin(MeV) wmax(MeV) Dxref(m) Dyref(m) Dxcog(m) Dycog(m)
       write(71,99)
99     format(1x,'#    l(m)  ',1x,'     x(mm)  ',1x,'     y(mm)  ',
     *        3x,'    z(deg) ',1x,'     z(mm)  ',
     *        2x,'emx(mm.mrd)',2x,'emy(mm.mrd) ',
     *        1x,'emz(KeV.ns)',2x,'energy(MeV) ',1x,'#particles',
     *        2x,'xmin(mm)',5x,'xmax(mm)',5x,'ymin(mm)',5x,'ymax(mm)',
     *        5x,'tmin(s)',6x,'tmax(s)',5x,'phmin(deg)',
     *        4x,'phmax(deg)',3x,'Wmin(MeV)',4x,'Wmax(MeV)',
     *        4x,'Dx(m)',8x,'Dy(m)')
       iprf1=iprf-1
       do i=1,iprf1
c     x and y: cm-->mm
        sprx=sprfx(i)*10./2.
        spry=sprfy(i)*10./2.
        sprp=sprfp(i)/2.
        sprz=sprfz(i)*10.
        write(71,100) sprfl(i),sprx,spry,sprp,sprz,eprnx(i)/4.,
     *               eprny(i)/4.,eprfw(i),sprww(i),int(sprng(i)),
     *               sxmn(i),sxmx(i),symn(i),symx(i),stmn(i),stmx(i),
     *               spmn(i),spmx(i),swmn(i),swmx(i),
     *               dispcx(i),dispcy(i)
c     *               disprx(i),dispry(i),dispcx(i),dispcy(i)
       enddo
100    format(9(1x,e12.5),1x,i7,3x,12(1x,e12.5))
c100    format(9(1x,e12.5),1x,i7,3x,14(1x,e12.5))
       return
       end
      SUBROUTINE etgap
c   ............................................................................
c   single cell of a DTL
c
c      etcell(1)=cell#; etcell(2)=energy (MeV),,etcell(3)=beta
c      etcell(4)=cell length (cm), etcell(5)=T, etcell(6)=TP,
c      etcell (7)=S,etcell(8)=SP
c      etcell(9)=quad length (cm), etcell (10)=quad. strength (kG/cm);
c      etcell(11)=Eo (MV/m) , etcell(12)= phase of RF at middle(deg),
c      etcell(13)= actual length(cm), etcell(14)=TPP
c      etcell(15)=frequency (MHz), etcell(16)=field factor
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFS/DYNT(MAXCELL),DYNTP(MAXCELL),DYNTPP(MAXCELL),
     *   DYNE0(MAXCELL),DYNPH(MAXCELL),DYNLG(MAXCELL),FHPAR,NC
       COMMON/POSI/IST
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/AZMTCH/DLG,XMCPH,XMCE
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
C      TRANSIT TIME COEFFICIENTS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/TTFC1/T3K,T4K,S3K,S4K
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/JACOB/GAKS,GAPS
       common/iter1/DXDKI,DPHII,PHI,DKMSKE,DKMSPHI,RETPH,XKMI,XKM,
     *              DXK00,TKE,T1KE,SKE,S1KE,PHIWC,XK1I,XK1II,XK2II
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/CDEK/DWP(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/compt1/ndtl,ncavmc,ncavnm
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TAPES/IN,IFILE,META
       COMMON/RANEC1/DUMMY(6)
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/speda/dave,idave
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/DCSPA/IESP
       COMMON/fct/FAKT
       common/mode/eflvl,rflvl
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/appel/irstay,ilost,iavp,ispcel
       common/posc/xpsc
       common/femt/iemgrw,iemqesg
       common/aerp/vphase,vfield,ierpf
       common/tofev/ttvols
C*bugfix*2010_07_03 : next line was missing; caused z=0 in envelope plot
       common/pstpla/tstp
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       character*128 trace3h,trace3t,tif
c --- synchronous particle:
c *****     DWRFS(MeV): energy gain
c *****     SPHRFS(rad): phase jump
c *****     PHRFS(rad):phase RF
c *****       common/parmrf/DWRFS,SPHRFS,PHRFS,ngdrf
       logical iesp,irstay,iavp,ispcel,iemgrw
       CHARACTER*1 cr
C         SHIFT =TRUE: cog and synchronous particle are independent
C         SHIFT =FALSE: cog and synchronous particle are coinciding
       LOGICAL SHIFT,CHASIT,ITVOL,IMAMIN,DAVE,ICHAES
       dimension etcell(maxcell)
c --- iesp, irstay, ispcel logical flags used in the routines of space sharge computations
c       iesp=false: accelerating element
c       ispcel = .true.: space charge computation
       iesp=.false.
       irstay=.false.
       ispcel=.true.
       ilost=0
       FACT=1.
       NDTL=NDTL+1
C*bugfix*2010_07_03 : next line was missing; caused # of accelerating elements to be zero
       nrres=nrres+1
c allow for print out on terminal of gap# on one and the same line
       cr=char(13)
#if F77
       WRITE(6,8254) NRTRE,NDTL,cr
8254   format('Transport element:',i5,
     *        '      Accelerating gap    :',i5,a1, $)
#else
       WRITE(6,8254,advance='no') NRTRE,NDTL,cr
8254   format('Transport element:',i5,
     *        '      Accelerating gap    :',i5,a1)
#endif
       WRITE(16,*)'ACCELERATING GAP N :',NDTL
       read(IN,*) (etcell(iet),iet=1,16)
c trace3d
c split cell into drif-gap-drift
       kt3t=kt3t+1
       write(tif,6001) kt3t,kt3t,10.*etcell(4)/2.
6001   format('  nt(',i4,')= 1, a(1,',i4,')=',f12.6)
       trace3t(kt3t)=tif
c
       tre0tl=etcell(11)*etcell(5)*etcell(4)/100.
       trphase=etcell(12)
       kt3t=kt3t+1
       write(tif,6005) kt3t,kt3t,tre0tl,trphase,fid
6005   format('  nt(',i4,')=10, a(1,',i4,')=',f9.5,' , ',f9.2,
     &        ', 0., 1.,',f5.3,',')
       trace3t(kt3t)=tif
c
       kt3t=kt3t+1
       write(tif,6001) kt3t,kt3t,10.*etcell(4)/2.
       trace3t(kt3t)=tif
c end trace3d drift-gap-drift
       FH=etcell(15)
       FAKT=etcell(16)
       if (FAKT.eq.0.) FAKT=1.e-12
       FH=FH*2.*PI*1000000.
C --- YLG : CELL LENGTH (cm)
C --- YE0 : ELECTRIC FIELD (MV/cm)
       YLG=etcell(4)
       scdist=ylg
       YE0=etcell(11)/100.
C      TRANSIT TIME COEFFICIENTS (at the middle OF THE CELL)
       T0=etcell(5)*YLG*YE0
       TP0=-etcell(6)*YLG*YLG*YE0
       TPP0=-etcell(14)*YLG*YLG*YLG*YE0
C      TRANSIT TIME factors at the ENTRANCE OF THE CELL
C       with kg=2*PI/ylg:
C         TK0=T0*COS(kg*ylg/2) = -T0
C         SK0=T0*SIN(kg*ylg/2) = 0
C         TPK0=d(TK0)/dk
C         SPK0=d(SK0)/dk
C         TPPK0=d(TPK0)/dk
C         SPPK0=d(SPK0)/dk
       TK0=-T0
       SK0=0.
       TPK0=-TP0
       SPK0=-YLG*T0/2.
       TPPK0=YLG*YLG*T0/4.-TPP0
       SPPK0=-YLG*TP0/2.
       TK=TK0
       T1K=TPK0
       T2K=TPPK0
       SK=SK0
       S1K=SPK0
       S2K=SPPK0
       TP3K0=0.
       TP4K0=0.
       SP3K0=0.
       SP4K0=0.
C  MULTIPLY TRANSIT TIME factors WITH FAKT
       T0=T0*FAKT
       TP0=TP0*FAKT
       TPP0=TPP0*FAKT
       TK0=-T0
       TPK0=-TP0
       SPK0=-YLG*T0/2.
       TPPK0=YLG*YLG*T0/4.-TPP0
       SPPK0=-YLG*TP0/2.
       TK=TK0
       T1K=TPK0
       T2K=TPPK0
       SK=SK0
       S1K=SPK0
       S2K=SPPK0
c       IPOINR=IPOINR+1
c   print in file: 'short.data'
C --- ylg : CELL LENGTH (cm) ==> (mm)
C --- ye0 : ELECTRIC FIELD (MV/cm) ==> (Kv/mm)
c --- davtot (mm)
       idav=idav+1
       iitem(idav)=17
       dav1(idav,1)=ylg*10.
       dav1(idav,2)=ye0*100.
C*bugfix*2010_07_03 : next line was missing; caused z=0 in envelope plot
       tstp=(davtot+ylg*xpsc)*10.
       davtot=davtot+ylg
       dav1(idav,24)=davtot*10.
c 21.11.09       dav1(idav,40)=fh
       FH0=FH/VL
C      STATISTICS FOR PLOT
       IF(IPRF.EQ.1) CALL STAPL(dav1(idav,24))
C      reference particle
       iarg=1
       call cdg(iarg)
       ecog=cog(1)
       enold=ecog
       gcog=ecog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       IF(SHIFT) THEN
C        reference particle and cog are independent
         BEREF=VREF/VL
         GAMREF=1./SQRT(1.-(BEREF*BEREF))
         ENREF=XMAT*GAMREF
         TREFDG=TREF*FH*180./PI
         tcogdg=tcog*fh*180./pi
         dav1(idav,3)=1.
       ELSE
C        reference and c.o.g. are coinciding
         vref=bcog*vl
         tref=tcog
         gamref=gcog
         beref=bcog
         enref=cog(1)
         TREFDG=TREF*FH*180./PI
         tcogdg=tcog*fh*180./pi
         dav1(idav,3)=0.
       ENDIF
c --- ttvol: time of flight at entrance (sec)
       ttvol=0.
       if(itvol)ttvol=ttvols*fh
c ***** reference is placed in the position ngdrf=ngood+1 in array f(10,i)
c *****         ngdrf=ngood+1
c *****         BEREF=VREF/VL
c *****         GAMREF=1./SQRT(1.-(BEREF*BEREF))
c *****         ENREF=XMAT*GAMREF
c ****         f(1,ngdrf)=ngdrf
c ****         f(2,ngdrf)=0.
c ****         f(3,ngdrf)=0.
c ****         f(4,ngdrf)=0.
c ****         f(5,ngdrf)=0.
c ****         f(6,ngdrf)=tref
c ****         f(7,ngdrf)=enref
c ****         f(8,ngdrf)=1.
c ****         f(9,ngdrf)=qst
c ****         f(10,ngdrf)=0.
         if(dav1(idav,3).eq.1.) write(16,*)
     *   ' ****reference and cog are different'
         if(dav1(idav,3).eq.0.) write(16,*)
     *   ' **** reference and cog coincide '
       WRITE(16,178)
178    FORMAT(/,' DYNAMICS AT THE INPUT ',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       WRITE(16,1788) bcog,gcog,ecog-xmat,tcogdg,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       WRITE(16,165) beref,gamref,enref-xmat,TREFDG,TREF
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       tofprt=tref
c      iprint=1: full print ( in case of pbmes)
       iprint=0
       call statis
       XK1=FH/VREF
c ---  predictor of the average beta:BEMY.
       BEREF=VREF/VL
       GAMREF = 1./SQRT(1.-(BEREF*BEREF))
C      DPHASE : PHASE  AT THE MIDDLE OF THE CELL (rad)
       dphete=etcell(12)
       dphase=etcell(12)*pi/180.
c 15/12/09       dpherd=dphase
c 15/12/09       ttvol=0.
c --- The option TOF is passive in the DTL
c   t.o.f. at the middle of the cell
c 15/12/09       if(itvol) then
c 15/12/09        tvolm=ylg/(2.*vref)
c 15/12/09        ttvol=(ttvols+tvolm)
c 15/12/09       endif
c 15/12/09 adjust the phase of RF w.r.t. TOF
c 15/12/09       if(itvol.and.imamin) then
c 15/12/09        odphase=dphase
c 15/12/09        ottvol=fh*ttvol*180./pi
c 15/12/09        attvol=ottvol
c 15/12/09        xkpi=ottvol/360.
c 15/12/09        ixkpi=int(xkpi)
c 15/12/09        xkpi=(xkpi-float(ixkpi))*360.
c 15/12/09        dphase=dphase-xkpi*pi/180.
c 15/12/09       endif
       aqst=abs(qst)
       ddw=aqst*t0*cos(dphase)
       ENREFS=ENREF+DDW
       GAMS= ENREFS/XMAT
       BETS=SQRT(1.-1./(GAMS*GAMS))
       XK2=FH/(BETS*VL)
       BEMY=(GAMS*BETS+GAMREF*BEREF)/(GAMS+GAMREF)
       XKM=FH/(BEMY*VL)
       XKG=2.*PI/YLG
c --- average phase of RF at entrance
       saphi=dphase
       SAPHO=SAPHI
       EQVL=YLG
       DKG=(XKM-XKG)
       DTS=TP0/T0
       FK1=2.*DTS
       FPK0=(TP0*TP0+T0*TPP0)/(T0*T0)
       FPK1=2.*TP0**2/(T0*T0)
       FPK=2.*(FPK0-FPK1)
       PCREST=ATAN(-SK/TK)
       DDW=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))/2.
       IF(DDW.LT.0.) PCREST=PCREST+PI
       PEQVL=YLG/2.
       SCLEN=YLG
c 15/12/09       if(.not.imamin) then
        WRITE(16,1555)FH/(2.*pi),YLG,DPHASE*180./PI
1555    FORMAT(4X,'FREQENCY :',E12.5,' Hertz',/,4x,
     x        'GAP LENGTH :',e12.5,' cm',/,4x,
     x        'PHASE of RF (middle of the gap)  :',e12.5,'deg',/)
c 15/12/09       else
c 15/12/09        WRITE(16,1556)FH/(2.*pi),YLG,DPHASE*180./PI,odphase
c 15/12/09 1556    FORMAT(4X,'FREQENCY :',E12.5,' Hertz',/,4x,
c 15/12/09     x        'GAP LENGTH :',e12.5,' cm',/,4x,
c 15/12/09     x        'phase of RF after adjustement:',e12.5,'deg',/,4x,
c 15/12/09     *        'phase of RF before adjustment',e12.5,'deg')
c 15/12/09       endif
C ---  Follow ITERATIONS giving:
c            The equivalent field length   (cm)
c            The asociated drift length    (cm)
c            The slip of phase             (rd)
c            The energy gain               (MeV)
c            The phase jump                (rd)
c            The average k                 (cm-1)
c            The transit time coefficients (MeV,cm)
c            The phase crest               (rad)
c            The phase offset at entrance  (rad)
       EQVL=YLG
       DKG=(XKM-XKG)
       FPK=2.*(FPK0-FPK1)
       TIL2=0.
       DO IT=1,3
C        slip of phase and equivalent field length
         IF(IT.EQ.1) PHSLIP=-4.*ATAN(3.2*DTS/EQVL)
         IF(PHSLIP.NE.0.) THEN
           TIL2=PHSLIP/2.
           DO IIII=1,4
             GX=1./TAN(TIL2)-1./TIL2
             GPX=-1./(SIN(TIL2)*SIN(TIL2)) + 1./(TIL2*TIL2)
             GPPX=2.*COS(TIL2)/(SIN(TIL2)**3) - 2./(TIL2*TIL2*TIL2)
             HX=GPX/(GX*GX) -2.*FPK/(FK1**2)
             DHX=-(2.*GX*GPX*GPX-GX*GX*GPPX)/(GX**4)
             TIL2=TIL2-HX/DHX
             EQVL=SQRT(ABS(2.*FPK/GPX))
             IF(ABS(HX).LE.1.E-05) GO TO 556
           ENDDO
556        CONTINUE
           PHSLIP=TIL2*2.
         ENDIF
         PEQVL=YLG/2.
         ASDL=PEQVL-EQVL/2.
C      ENERGY GAIN AND PHASE JUMP (i.e. DELPHR)
c         saphi=sapho-pcrest+ttvol*fh
         saphi=sapho-pcrest+ttvol
c TESTsv05.05.2011         saphi=sapho-pcrest
         F0=XITL0(GAMREF,GAMS,BEMY,SAPHI,AQST)
         DELWRM=(F0-GAMREF)*XMAT
         ENRS=ENREF+DELWRM
         GAMS= ENRS/XMAT
         bets=sqrt(1.-1./(gams*gams))
         xk2=fh0/bets
         coeph=fh*aqst/(vl*xmat)
         F2=XITL2(GAMREF,GAMS,BEMY,SAPHI,AQST)
         DELPHR= COEPH * F2
         XKM=DELPHR/EQVL + XK2*(1.+ASDL/EQVL) - XK1*ASDL/EQVL
         BEMY=FH0/XKM
C    TRANSIT TIME FACTORS
         DKG=(XKM-XKG)
         TK=TK0+DKG*TPK0+DKG*DKG*TPPK0/2.+DKG**3*TP3K0/6.+
     X     DKG**4*TP4K0/24.
         T1K=TPK0+DKG*TPPK0+DKG*DKG*TP3K0/2.+DKG**3*TP4K0/6.
         T2K=TPPK0+DKG*TP3K0+DKG*DKG*TP4K0/2.
         T3K=TP3K0+DKG*TP4K0
         T4K=TP4K0
         SK=SK0+DKG*SPK0+DKG*DKG*SPPK0/2.+DKG**3*SP3K0/6.+
     X     DKG**4*SP4K0/24.
         S1K=SPK0+DKG*SPPK0+DKG*DKG*SP3K0/2.+DKG**3*SP4K0/6.
         S2K=SPPK0+DKG*SP3K0+DKG*DKG*SP4K0/2.
         S3K=SP3K0+DKG*SP4K0
         S4K=SP4K0
C      PHASE CREST
         PCREST=ATAN(-SK/TK)
         DDW=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))/2.
        IF(DDW.LT.0.) PCREST=PCREST+PI
       ENDDO
C      CREST VALUE = SQCTTF*(PHSLIP/2)/SIN(PHSLIP/2)
       SQCTTF=TIL2*SQRT(TK*TK+SK*SK)/SIN(TIL2) *2.
C      THE synchronous particle IS BASED ON CHARGE STATE : QMOY =AQST
       CFH=FH/(VL*2.*XMAT)
       CKH=QMOY*QMOY/(4.*XMAT*XMAT)
c ******---  save the energy and the T.O.F of the reference before the routine:GAP
c *****       enri=f(7,ngdrf)
c *****       trefi=f(6,ngdrf)
       call gap(gamref,saphi,gams,delphr)
c *****--- trefs and enrs: T.O.F. and  energy of the reference at output of the acc. element
c ****       trefs=f(6,ngdrf)
c ******       enrs=f(7,ngdrf)
c *****       grefs=f(7,ngdrf)/xmat
c *****       bets=sqrt(grefs*grefs-1.)/grefs
c ******       delwrm=enrs-enri
C      output of the element
c  new PHARES,TREFS to be in accordance with GENAC
       PHARES=SAPHI+FH*YLG/VREF+DELPHR
       trefs=tref+ylg/(bets*vl)+delphr/fh
       PHARED=(PHARES-SAPHI)*180./PI
       TREDG=fh*TREFS *180./PI
       write(16,*) ' PARAMETERS RELATING TO THE REFERENCE PARTICLE '
       write(16,*) '************************************************'
       write(16,*) ' ENERGY GAIN(MeV): ',delwrm,' TOF(DEG) ',tredg
c       write(16,*) ' PHASE JUMP(DEG): ',sphrfs*180./PI
       write(16,*) ' CREST PHASE OF RF (DG): ',
     *               PCREST*180./PI
       write(16,*) ' PHASE OF RF AT THE MIDDLE (DG): ',sapho*180./pi
       write(16,*) ' PHASE OF RF AT THE ENTRANCE (DG): ',saphi*180./pi
       write(16,*) ' AVERAGE k (cm-1) (freq./velocity): ',XKM
       write(16,*) ' TRANSIT TIME FACTORS (MeV-cm):'
       write(16,*) ' T dT/dk d2T/dk2 ',TK,T1K,T2K
       WRITE(16,*) ' S dS/dk d2S/dk2 ',SK,S1K,S2K
       write(16,*) ' PHASE SLIP(DEG) ',PHSLIP*180./PI
       write(16,*)
       write(16,*) ' PARAMETERS RELATING TO THE EQUIVALENT FIELD '
       write(16,*) '************************************************'
       WRITE(16,171)EQVL
171    FORMAT(' length :',e12.5,' cm ')
       write(16,*) ' Associated drift length: ',asdl,' cm'
       WRITE(16,*)' magnitude: ',SQCTTF,' MV/cm'
        iarg=1
        call cdg(iarg)
        encog=cog(1)
        gcog=encog/xmat
        bcog=sqrt(1.-1./(gcog*gcog))
        tcog=cog(3)
        CALL EXT2D(1)
c sup        phnew=-(int(tcog*fh/pi+0.5)-tcog*fh/pi)*180.
c sup        dav1(idav,7)=phnew
        dav1(idav,37)=saphi*180./pi
c 15/12/09        if(itvol) then
c 15/12/09        dav1(idav,38)=dphete
c 15/12/09        dav1(idav,39)=dphase*180./pi
c 15/12/09        else
        dav1(idav,38)=dphete
c 15/12/09        endif
c  end print in the file: 'short.data'
       WRITE(16,3777)
3777   FORMAT(/,3X,3(1H*),' DYNAMICS AT THE OUTPUT: ',/,
     2 5X,'   BETA     dW(MeV)    ENERGY(MeV) ',
     4 '   TOF(deg)     TOF(sec)')
       write(16,3473) bets,delwrm,enrs-xmat,fh*trefs*180./pi,trefs
3473   FORMAT(' REF ',f7.5,3x,f10.6,3x,f8.3,3x,e12.5,3x,e12.5)
       WRITE(16,1789) bcog,encog-enold,encog-xmat,tcog*fh*180./pi,tcog
1789   FORMAT(' COG ',f7.5,3x,f10.6,3x,f8.3,3x,e12.5,3x,e12.5)
       TESTCA=exten(1)*exten(2)*exten(3)
       epsil=1.E-40
       IF(abs(TESTCA).gt.epsil) THEN
         qdisp=2.*sqrt(exten(1))
         qmd=exten(1)*exten(3)-exten(2)**2
         SQMDV=4.*PI*SQRT(QMD)
         SURM=4.*PI*SQRT(QMD)*180./PI
         qdp=2.*sqrt(exten(3))
         cor12=exten(2)/sqrt(exten(1)*exten(3))
         QDPDE=QDP*180./PI
       ELSE
         QDISP=0.
         QMD=0.
         SQMDV=0.
         SURM=0.
         QDP=0.
         COR12=0.
         PENT12=0.
         PENT21=0.
         QDPDE=0.
       ENDIF
       TRQTX=exten(4)*exten(5)-exten(8)**2
       TRQPY=exten(6)*exten(7)-exten(9)**2
       QDITAX=2.*SQRT(exten(4))
       QDIANT=2.*SQRT(exten(5))
       QDITAY=2.*SQRT(exten(6))
       QDIANP=2.*SQRT(exten(7))
       SURXTH=4.*PI*SQRT(TRQTX)
       SURYPH=4.*PI*SQRT(TRQPY)
       IF(SHIFT) THEN
         vref=bets*vl
         tref=trefs
       ELSE
         vref=bcog*vl
         tref=tcog
       ENDIF
       if(itvol) then
        ttvols=tref
c 15/12/09        attvol=fh*ttvols*180./pi
c 15/12/09        write(16,7456) ottvol,attvol
       endif
c 15/12/09 7456   format(2x,'***tof at input: ',e12.5,' deg',/,
c 15/12/09     *           2x,'***tof at output: ',e12.5,' deg')
       call statis
C      PROFIL (plot)
       CALL STAPL(dav1(idav,24))
       dltaw=qdisp*xmat*bcog*bcog/sqrt(1.-bcog*bcog)
c sup       WRITE(16,9998) SQMDV
c sup9998   FORMAT(2X,'   EMITTANCE (norm): ',
c sup     *        E12.5,' PI*MEV*RAD')
c   print in the file: 'dynac.dmp':
c   gap number, phase offset(deg), relativistic beta, energy(MeV), horz. emit.(mm*mrd,norm), vert. emit.(mm*mrd,norm),long. emit(keV*sec)
c  dav1(idav,16): Emittance(norm)  x-xp (mm*mrad)
       dav1(idav,16)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
c  dav1(idav,21): Emittance(norm)  y-yp (mm*mrad)
       dav1(idav,21)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
       dav1(idav,25)=ndtl
       emns=1.e12*sqmdv/(pi*fh)
cet2010s
       trfprt=fh*tref*180./pi
       tcgprt=fh*tcog*180./pi
c       n2kp=int(tofprt/360.)
c       tofprt=tofprt-float(n2kp)*360.
c       if(tofprt.gt.180.) tofprt=tofprt-360.
c cavity number, z(m), transmission (%), synchronous phase (deg), time of flight (deg) (cog), COG relativistic beta (@ output)
c COG output energy (MeV), time of flight (deg) (REF), REF relativistic beta (@ output), REF output energy (MeV),
c horizontal emittance (mm.mrad, RMS normalized), vertical emittance (mm.mrad, RMS normalized),
c longitudinal emittance (RMS, ns.keV)
       trnsms=100.*float(ngood)/float(imax)
       if(ndtl.eq.1) write(50,*) '# gap.dmp'
       if(ndtl.eq.1) write(50,*) '# gap     Z       trans   ',
     *   'PHIs     TOF(COG)    COG        Wcog          TOF(REF)   ',
     *   '    REF        Wref       Ex,RMS,n     Ey,RMS,n     El,RMS'
       if(ndtl.eq.1) write(50,*) '#  #     (m)       (%)    ',
     *  '(deg)     (deg)      beta       (MeV)          (deg)      ',
     *  '   beta       (MeV)      (mm.mrad)    (mm.mrad)    (ns.keV)'
       write(50,7023) ndtl,0.01*davtot,trnsms,dphete,tcgprt,
     *  bcog,encog-xmat,trfprt,bets,enrs-xmat,
     *  0.25*dav1(idav,16),0.25*dav1(idav,21),0.25*emns
7023     format(1x,i4,1x,e12.5,1x,f6.2,1x,f7.2,1x,
     *   2(e14.7,1x,f7.5,1x,e14.7,1x),3(e12.5,1x))
cet2010e
C      RIGIDITY of the ref. prtcle
       gref=1./sqrt(1.-bets*bets)
       XMOR=XMAT*bets*gref
       BORO=33.356*XMOR*1.E-01/QST
       WRITE(16,*) ilost,' particles are lost in element ',ndtl
       write(16,*)
       call emiprt(0)
       return
       end
       SUBROUTINE gap(GAMREF,SAPHI,GAMS,DELPHR)
c   ............................................................................
c      ETGAP or RESTAY ==> GAP
c        dynamics in the accelerating element
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFS/DYNT(MAXCELL),DYNTP(MAXCELL),DYNTPP(MAXCELL),
     *   DYNE0(MAXCELL),DYNPH(MAXCELL),DYNLG(MAXCELL),FHPAR,NC
       COMMON/POSI/IST
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/AZMTCH/DLG,XMCPH,XMCE
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
C      TRANSIT TIME COEFFICIENTS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
C ***************************************************************
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/JACOB/GAKS,GAPS
       common/iter1/DXDKI,DPHII,PHI,DKMSKE,DKMSPHI,RETPH,XKMI,XKM,
     *              DXK00,TKE,T1KE,SKE,S1KE,PHIWC,XK1I,XK1II,XK2II
C         functions in COS(DXDPT)
       common/iterco/YH11T,YH1K1T,YH1K01T,YH10PKT,YH11PKT,YH1P1T,
     *                H1AKIT,H1AKIMT,H1AKMT,HAPIT,HAPPIT
C         functions in SIN(DXDPT)
       common/itersi/YH21T,YH2K1T,YH2K01T,YH2P1T,YH20PKT,YH21PKT,
     *                H1BKIT,H1BKIMT,H1BKMT,HBPIT,HBPPIT
       common/tranrs/SA11,SA12,SA21,SA22,SACT11,SACT12,SACT21,SACT22
C --- routine XTYPL1
       COMMON/TYPL1/YH1K0,YH1K1,YP1K1,YP1K2,YH1K00,YH1K01,YP1K01,
     X              YP1K02,YH10,YH11,YP11,YP12
       COMMON/TYPL2/YH2K0,YH2K1,YP2K1,YP2K2,YH2K00,YH2K01,YP2K01,
     X              YP2K02,YH20,YH21,YP21,YP22
       COMMON/TYPI1/YE1K0,YE1K1,YE1K2,YE1KC0,YE1KC1,YE1KC2,
     X              YE10,YE11,YE12
       COMMON/TYPI2/YE2K0,YE2K1,YE2K2,YE2KC0,YE2KC1,YE2KC2,
     X              YE20,YE21,YE22
       COMMON/THAD2/H0AKI,H0AKIM,H0AKM,H0BKI,H0BKIM,H0BKM,
     X              H1AKI,H1AKIM,H1AKM,H1BKI,H1BKIM,H1BKM
C --- routine XTYPLP1
       COMMON/TYPLP1/YH1P1,YH2P1,HAPI,HBPI
C --- routine XTYPL2
       COMMON/TYPLP2/HAPPI,HBPPI
C --- routine XTYLPK
       COMMON/TYPLPK/YH10PK,YH11PK,YH20PK,YH21PK
C     Integrals of E(z)**2
C --- routine XTYPJ
       COMMON/TYPJ/YFSK0,YFSK1,YFSK2,YFSP0,YFSP1,YFSP2,
     X              YFSKC0,YFSKC1,YFSKC2,YFSCK0,YFSCK1,YFSCK2,
     X              YFSCP0,YFSCP1,YFSCP2,YFS0,YFS1,YFS2
C --- routine XTYPM
       COMMON/TYPM/YNSK0,YNSK1,YNSK2,YNSP0,YNSP1,YNSP2,
     X              YNSK0C,YNSK1C,YNSK2C,YNS0,YNS1,YNS2
C ********************************************************************
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/CDEK/DWP(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON /CONSTA/ VL, PI, XMAT, RPEL,QST
       COMMON/DYN/TREF,VREF
       COMMON/TAPES/IN,IFILE,META
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/speda/dave,idave
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/DCSPA/IESP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/appel/irstay,ilost,iavp,ispcel
       common/posc/xpsc
       common/pstpla/tstp
       common/rander/ialin
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       logical iesp,ichaes,irstay,iavp,ispcel,ialin
       LOGICAL SHIFT,CHASIT,ITVOL,IMAMIN,DAVE
c *****      DWRFS(MeV):  gain of energy of the fictitious reference
c *****     SPHRFS(rad): phase jump
c *****      PHRFS(rad):  phase
c *****     NGDRF:  position of the reference in the array f(10,iptsz)
c ****       common/parmrf/DWRFS,SPHRFS,PHRFS,ngdrf
c      CHARACTER*1 CHOPT
       FH0=FH/VL
c   iesp is used in s.c. routines : iesp=.false. == > accel. element
       iesp=.false.
       iavp=.true.
       ipas=2
       DCUM=YLG
c       random errors in alignment
       if(ialin) call randali
C       Random variation on the phase for each particle
       VARPHA=0.
C   INTEGRALS REQUIRED BY THE EXPANSIONS IN LONGITUDINAL MOTION (ALL THE ELEMENT)
c   ( computations of the average k and slip of phase from which are computed deltk delphi,..)
c    *   valero 08/08/07
       aqmoy=abs(qst)
comment       aqmoy=qst
c    *
       CALL XTYPL1(GAMREF,SAPHI,AQMOY,DCUM)
C  see here equations 88 and 106 of Part. Acc. 1994 vol44 pp215-255
       CXLG=aqmoy/(4.*XMAT*EQVL)
       GAMI=GAMREF
       DKMP= (GAMI*GAMI-1.)**(1.5)*(GAMS*GAMS-1.)**(-1.5)
       DKMS= DKMP*(1.+ASDL/EQVL)
     X       +YH1K01*FH0*CXLG/EQVL -ASDL/EQVL
       DKM1= -GAKS*(GAMS*GAMS-1.)**(-1.5)*FH0*(1.+ASDL/EQVL)
C    DKMSKE :(eq.106, see article of Part.Acc.)
       DKMSKE=DKMS/(1.-YH1K1*CXLG*FH0/EQVL-DKM1)
       CALL XTYPL2(GAMREF,SAPHI,AQMOY,DCUM)
       CALL XTYPLP1(GAMREF,SAPHI,AQMOY,DCUM)
       CALL XTYLPK(GAMREF,SAPHI,AQMOY,DCUM)
       dphsph1=(yh1p1-yh21)*cxlg*fh0
       DKMSPHI=-FH0*(GAMS*GAMS-1.)**(-1.5)*GAPS*(1.+ASDL/EQVL)+
     X          DPHSPH1/EQVL
c   These integrals are saved
       gakst=gaks
       gapst=gaps
C         functions in COS(DXDPT)
       YH11T=YH11
       YH1K1T=YH1K1
       YH1K01T=YH1K01
       YH10PKT=YH10PK
       YH11PKT=YH11PK
       YH1P1T=YH1P1
       H1AKIT=H1AKI
       H1AKIMT=H1AKIM
       H1AKMT=H1AKM
       HAPIT=HAPI
       HAPPIT=HAPPI
C         functions in SIN(DXDPT)
       YH21T=YH21
       YH2K1T=YH2K1
       YH2K01T=YH2K01
       YH2P1T=YH2P1
       YH20PKT=YH20PK
       YH21PKT=YH21PK
       H1BKIT=H1BKI
       H1BKIMT=H1BKIM
       H1BKMT=H1BKM
       HBPIT=HBPI
       HBPPIT=HBPPI
       if(ichaes.and.ispcel) then
         ipas=1
         write(16,*)'  SPACE CHARGE ACTING ON LENGTH: ',scdist,' CM'
         dcum=ylg*xpsc
         write(16,*)'  POSITION OF S.C. COMPUTATION:  ',dcum,' CM'
c     computation of the integrals in the middle of the cavity
C INTEGRALS  E(z)*(BG)**-3 *z**n   n=0,1
C INTEGRALS  dE(z)/dt*(BG)**-3 *z**n   n=0,1,2
         call xtypl1(gamref,saphi,aqmoy,dcum)
         istm=ist-1
C     INTEGRALS ON SECOND DERIVATIVES k HA0(Z) & HB0(Z)
         call xtypl2(gamref,saphi,aqmoy,dcum)
C      1st,2nd,3rd DERIVATIVES ON PHASE,HA0(Z) & HB0(Z)
         call xtyplp1(gamref,saphi,aqmoy,dcum)
C    DERIVATIVES  COUPLED ON PHASE, K FUNCTIONS HA0(Z) & HB0(Z)
         call xtylpk(gamref,saphi,aqmoy,dcum)
C     TRANSVERSE INTEGRALS  TYPE J & M
         call xtypj(gamref,saphi,aqmoy,dcum)
         call xtypm(gamref,saphi,aqmoy,dcum)
       endif
1026   continue
       if(ipas.eq.2) then
         dcum=ylg
         gaks=gakst
         gaps=gapst
c --- COS(DXDPT)
         yh11=yh11t
         yh1k1=yh1k1t
         yh1k01=yh1k01t
         yh10pk=yh10pkt
         yh11pk=yh11pkt
         yh1p1=yh1p1t
         h1aki=h1akit
         h1akim=h1akimt
         h1akm=h1akmt
         hapi=hapit
         happi=happit
c ---  SIN(DXDPT)
         yh21=yh21t
         yh2k1=yh2k1t
         yh2k01=yh2k01t
         yh2p1=yh2p1t
         yh20pk=yh20pkt
         yh21pk=yh21pkt
         h1bki=h1bkit
         h1bkim=h1bkimt
         h1bkm=h1bkmt
         hbpi=hbpit
         hbppi=hbppit
         call xtypl1(gamref,saphi,aqmoy,dcum)
C     INTEGRALS ON SECOND DERIVATIVES k HA0(Z) & HB0(Z)
         call xtypl2(gamref,saphi,aqmoy,dcum)
C      1st,2nd,3rd DERIVATIVES ON PHASE,HA0(Z) & HB0(Z)
         call xtyplp1(gamref,saphi,aqmoy,dcum)
C    DERIVATIVES  COUPLED ON PHASE, K FUNCTIONS HA0(Z) & HB0(Z)
         call xtylpk(gamref,saphi,aqmoy,dcum)
C       TRANSVERSE INTEGRALS  TYPE J & M
         call xtypj(gamref,saphi,aqmoy,dcum)
         call xtypm(gamref,saphi,aqmoy,dcum)
       endif
       SA11=1.
       SA12=0.
       SA21=0.
       SA22=1.
       SACT11=1.
       SACT12=0.
       SACT21=0.
       SACT22=1.
c    SAVE BEAM
       do is=1,ngood
         do js=1,7
           fs(js,is)=f(js,is)
         enddo
       enddo
c1055   CONTINUE
       call boucle(ipas,gamref,saphi,dcum,delphr)
c  Reshuffles f(i,j) array after boucle
       call shuffle
       IF(IPAS.EQ.1) THEN
C   compute the space charge on the beam (except the reference)
c  call stapl at the position of space charge computation
c*et*26-Jul-2014         call stapl(tstp)
         if(iscsp.eq.1) then
           ini=1
           call hersc(ini)
           ini=2
           call hersc(ini)
         endif
         if(iscsp.eq.2) CALL SCHERMI
         if(iscsp.eq.3) CALL SCHEFF1(1)
C -----  window control
         write(16,*) 'Checking for lost particles'
         call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c    second step of the gap
         ipas=2
         iavp=.false.
         go to 1026
C      end of space charge computation
       ENDIF
C       charateristics of the beam
c 18/03/2009       if (dave) then
c 18/03/2009         gimax=f(7,ngood)/xmat
c 18/03/2009         bimax=sqrt(1.-1./(gimax*gimax))
c 18/03/2009         dav1(idav,3)=bimax
c 18/03/2009         dav1(idav,4)=f(7,ngood)-xmat
c 18/03/2009         dav1(idav,5)=-(INT(f(6,ngood)*FH/PI+0.5)-
c 18/03/2009     *                f(6,ngood)*FH/PI)*180.
c 18/03/2009         dav1(idav,38)=dphase*180./pi
c 18/03/2009       endif
C ----       WINDOW CONTROL
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       RETURN
       END
       SUBROUTINE cogetc
c      COG of TOF with regard the various charge states
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       itot=0
       n=1
       do i=1,20
         nbch(i)=0
       enddo
       charm(n)=f(9,1)
100    continue
       do i=1,ngood
         if(f(9,i).eq.charm(n)) nbch(n)=nbch(n)+1
       enddo
c   the following charge states
       itot=itot+nbch(n)
       if(itot.ge.ngood) go to 110
       do i=1,ngood
         do j=1,n
           if(f(9,i).eq.charm(j)) go to 120
         enddo
         if(f(9,i).ne.charm(n)) then
           n=n+1
           charm(n)=f(9,i)
           go to 100
         endif
120      continue
       enddo
110    continue
c   compute the cog of TOF for each charge state
       do i=1,n
         cgtdv(i)=0.
         do j=1,ngood
           if(f(9,j).eq.charm(i)) cgtdv(i)=cgtdv(i)+f(6,j)
         enddo
         cgtdv(i)=cgtdv(i)/float(nbch(i))
c        write(16,*)'i,nch,cgtdv=',i,nbch(i),cgtdv(i),charm(i)
       enddo
       netac=n
       return
       end
       SUBROUTINE boucle(ipas,gamref,saphi,dcum,delphr)
c   ............................................................................
c      RESTAY or ETGAP ==> GAP ==> BOUCLE
c      compute the dynamics of the fictitiuos reference and of particles
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/POSI/IST
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/AZMTCH/DLG,XMCPH,XMCE
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/ITVOLE/ITVOL,IMAMIN
       common/iter1/DXDKI,DPHII,PHI,DKMSKE,DKMSPHI,RETPH,XKMI,XKM,
     *              DXK00,TKE,T1KE,SKE,S1KE,PHIWC,XK1I,XK1II,XK2II
       common/tranrs/SA11,SA12,SA21,SA22,SACT11,SACT12,SACT21,SACT22
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       common/blvl/bflvl
C      TRANSIT TIME COEFFICIENTS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/TTFCB/T3K,T4K,S3K,S4K
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/JACOB/GAKS,GAPS
       COMMON/TYPL1/YH1K0,YH1K1,YP1K1,YP1K2,YH1K00,YH1K01,YP1K01,
     X              YP1K02,YH10,YH11,YP11,YP12
       COMMON/TYPL2/YH2K0,YH2K1,YP2K1,YP2K2,YH2K00,YH2K01,YP2K01,
     X              YP2K02,YH20,YH21,YP21,YP22
       COMMON/TYPI1/YE1K0,YE1K1,YE1K2,YE1KC0,YE1KC1,YE1KC2,
     X              YE10,YE11,YE12
       COMMON/TYPI2/YE2K0,YE2K1,YE2K2,YE2KC0,YE2KC1,YE2KC2,
     X              YE20,YE21,YE22
       COMMON/THAD2/H0AKI,H0AKIM,H0AKM,H0BKI,H0BKIM,H0BKM,
     X              H1AKI,H1AKIM,H1AKM,H1BKI,H1BKIM,H1BKM
C --- routine XTYPLP1
       COMMON/TYPLP1/YH1P1,YH2P1,HAPI,HBPI
C --- routine XTYPL2
       COMMON/TYPLP2/HAPPI,HBPPI
C --- routine XTYLPK
       COMMON/TYPLPK/YH10PK,YH11PK,YH20PK,YH21PK
C     integrals of E(z)**2
C --- routine XTYPJ
       COMMON/TYPJ/YFSK0,YFSK1,YFSK2,YFSP0,YFSP1,YFSP2,
     X              YFSKC0,YFSKC1,YFSKC2,YFSCK0,YFSCK1,YFSCK2,
     X              YFSCP0,YFSCP1,YFSCP2,YFS0,YFS1,YFS2
C     E(z)**2 INTEGRALS (COMPLEMENTARY ELECTRIC FIELD)
C --- routine XTYPM
       COMMON/TYPM/YNSK0,YNSK1,YNSK2,YNSP0,YNSP1,YNSP2,
     X              YNSK0C,YNSK1C,YNSK2C,YNS0,YNS1,YNS2
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/CDEK/DWP(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/CONSTA/ VL, PI, XMAT, RPEL,QST
       COMMON/DYN/TREF,VREF
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TAPES/IN,IFILE,META
       COMMON/DCSPA/IESP
       common/tcav/SV1P(iptsz),SV2P(iptsz),SXV1P(iptsz),SXV2P(iptsz),
     *       DWCIS(iptsz),BEINI1(iptsz),PHIP(iptsz),
     *       TEGLP(iptsz),DXDPIP(iptsz),DXDKIP(iptsz),
     *       DXDPTP(iptsz),DXK00P(iptsz),DPHIIP(iptsz),sauphcs(iptsz)
       common/iterco/YH11T,YH1K1T,YH1K01T,YH10PKT,YH11PKT,YH1P1T,
     *                H1AKIT,H1AKIMT,H1AKMT,HAPIT,HAPPIT
       common/itersi/YH21T,YH2K1T,YH2K01T,YH2P1T,YH20PKT,YH21PKT,
     *                H1BKIT,H1BKIMT,H1BKMT,HBPIT,HBPPIT
       common/appel/irstay,ilost,iavp,ispcel
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       common/tofev/ttvols
       common/aerp/vphase,vfield,ierpf
c *****     DWRFS(MeV): gain of energy of the fictitious reference
c *****     SPHRFS(rad):phase jump
c *****     PHRFS(rad):phase
c *****     NGDRF:  position of the reference in the array f(10,iptsz)
c *****       common/parmrf/DWRFS,SPHRFS,PHRFS,ngdrf
       logical itvol,imamin,ispcel
       logical iesp,iavp,ichaes,irstay
       dimension vecx(1)
c      CHARACTER*1 CHOPT
       ttvol=ttvols*fh
       FH0=FH/VL
       beref=sqrt(1.-1./(gamref*gamref))
       call cogetc
       tcog=0.
       gcog=0.
       do i=1,ngood
         gcog=gcog+f(7,i)/xmat
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       gcog=gcog/float(ngood)
       bcog=sqrt(1.-1./(gcog*gcog))
       wcg=(gcog-1.)*xmat
c ----- convert wdisp in dp/p (window control)
c ---- ifw = 0 ===> wdisp = dW/W
c ---- ifw = 1 ===> wdisp = dW (MeV)
c ----- convert wdisp in dp/p
       if(ifw.eq.0) dispr=gcog*gcog*wdisp/(gcog*(gcog+1.))
       if(ifw.eq.1) dispr=gcog*gcog*wdisp/(gcog*(gcog+1.)*wcg)
       DO 19 I=1,ngood
         IF(I.EQ.ICONT) THEN
           WRITE(16,*)'******************************************'
           WRITE(16,*)'*** FOLLOWED PARTICLE NUMBER: ',I
         IF(IPAS.EQ.1)
     *         WRITE(16,*) ' AT SPACE CHARGE POSITION'
         IF(IPAS.EQ.2) WRITE(16,*) ' AT OUTPUT '
           WRITE(16,*)'******************************************'
         ENDIF
         if(iavp) then
           gini=f(7,i)/xmat
           beini=sqrt(1.-1./(gini*gini))
           fd(i)=(gini*beini)/(gcog*bcog)
           f6i=0.
           do istc=1,netac
             if(f(9,i).eq.charm(istc))f6i=f(6,i)-cgtdv(istc)
           enddo
C       LONGITUDINAL WINDOW CONTROL
           if(fh*abs(f6i).ge.wphas) f(8,i)=0.
           IF(ABS(FD(I)-1.).GE.DISPR) f(8,i)=0.
C       TRANSVERSE WINDOW CONTROL
           RADIU=SQRT(F(2,I)*F(2,I)+F(4,I)*F(4,I))
           IF(RADIU.GE.RLIM) f(8,i)=0.
           IF(abs(f(2,i)).gt.wx) f(8,i)=0.
           IF(abs(f(4,i)).gt.wy) f(8,i)=0.
           if(f(8,i).eq.0) then
             write(16,3928) i,int(f(1,i)),f(2,i),f(3,i),f(4,i),f(5,i),
     *                      f6i*fh*180./pi,f(7,i)-xmat,int(f(9,i))
3928         FORMAT(' # ',i5,1x,i5,1x,6(f10.2,1x),1x,i2)
             ilost=ilost+1
             if(ilost.ge.ngood) stop
             go to 19
           endif
C   COMPUTATION OF AVERAGE K AND JUMP OF PHASE FOR THE CELL
C  START IF BLOCK on ICONT
           IF(I.EQ.ICONT) THEN
             WRITE(16,558) F(2,I),F(3,I),F(4,I),F(5,I)
558          FORMAT(1x,'* INPUT OF THE ELEMENT: ',/,1x,
     X       '* X :',E12.5,' CM  XP :',E12.5,' MRD',/,1x,
     X       '* Y :',E12.5,' CM  YP: ',E12.5,' MRD',/,1x,'*')
              f6dg=fh*f(6,i)*180./pi
             WRITE(16,*)' Tof(deg): ',f6dg,' ENER(MeV) ',f(7,i)-xmat
           ENDIF
           RADIU=SQRT(F(2,I)*F(2,I)+F(4,I)*F(4,I))
           IF(RADIU.LT.1.e-06) THEN
             DRADIU=.001*SQRT(F(3,I)*F(3,I)+F(5,I)*F(5,I))
           ELSE
             DRADIU=F(3,I)*.001*F(2,I)/RADIU+F(4,I)*F(5,I)*.001/RADIU
           ENDIF
c335        CONTINUE
C --- retph: phase delay between the actual particle  and the reference
           RETPH=FH*(f(6,i)-TREF)
c --- systematic or random defaults on the phase RF (not for the reference)
           if(ierpf.ne.0) then
             vphasi=vphase*pi/180.
c --- systematic default on the phase offset
             if(ierpf.eq.1) retph=retph+vphasi
c --- random error on  phase offset
             if(ierpf.gt.1) then
               len=1
               call rlux(vecx,len)
               r1=(2.*vecx(1)-1.)*vphasi
               retph=retph+r1
             endif
           endif
           PHI=SAPHI+RETPH
           IF(I.EQ.ICONT) write(16,*)
     *    '* PHASE DELAY RELATIVE TO REFERENCE ',retph*180./pi,' DEG'
c --- the TTF are the ones of the reference particle based on the charge state qst (input)
c ---  They are corrected in order to take into account the charge state of the current particle
c ****       cort=-(1.-(qst/f(9,i)))
           cort=0.
           TKC=TK*cort+TK
           T1KC=T1K*cort+T1K
           T2KC=T2K*cort+T2k
           T3KC=T3K*cort+T3k
           T4KC=T4K*cort+T4k
           SKC=SK*cort+SK
           S1KC=S1K*cort+S1K
           S2KC=S2K*cort+S2K
           S3KC=S3K*cort+S3K
           S4KC=S4K*cort+S4K
c       cort1=qst/f(9,i)
           cort1=1.
c ************************************
c   Predictor
           PHASE=PHI
C  DXKI0 : difference of energy betwen the reference and the actual particle
           DXKI0=FH0*(1./beini-1./BEREF)
           DXDTE=DXKI0
           TKE=TKC+DXDTE*T1KC+DXDTE*DXDTE*T2KC/2.+DXDTE**3*T3KC/6.+
     X     DXDTE**4*T4KC/24.
           SKE=SKC+DXDTE*S1KC+DXDTE*DXDTE*S2KC/2.+DXDTE**3*S3KC/6.+
     X     DXDTE**4*S4KC/24.
c  TEST*******
ccc           TKE=TKC
ccc           SKE=SKC
c  ******************************************
c --- systematic or random defaults on the field level (not for the reference)
           if(ierpf.ne.0) then
             if(ierpf.eq.1) then
               TKE=TKE*(1.+vfield)
               SKE=SKE*(1.+vfield)
             endif
             if(ierpf.gt.1) then
               len=1
               call rlux(vecx,len)
               r1=(2.*vecx(1)-1.)*vfield
               TKE=TKE*(1.+r1)
               SKE=SKE*(1.+r1)
             endif
           endif
           PHIWC=PHI+PAVPH
c  TEST***********
cccc           PHIWC=PHI
c *****************************************************
           DDWP=abs(f(9,i))*(TKE*COS(PHIWC)-SKE*SIN(PHIWC))
           ENPMT=f(7,i)+DDWP
           GAMPS=ENPMT/XMAT
           if(gamps.le.1.) f(8,i)=0.
           if(f(8,i).eq.0.) then
             ilost=ilost+1
             if(ilost.ge.ngood) stop
             go to 19
           endif
           BETPS=SQRT(1.-1./(GAMPS*GAMPS))
           XK1II=FH0/BEINI
           XK2II=FH0/BETPS
           XKMI=XK2II+(XK2II-XK1II)*ASDL/EQVL + DELPHR/EQVL
           XK1I=XK1II-XKMI
           XK2I=XK2II-XKMI
           BEMPY=FH0/XKMI
           SAUPHC=DELPHR
           DXDKI=XKMI-XKM
           DPHII=(XK1II-XK2II)*EQVL/10.+(XKP1+XKP2)/120.*EQVL**2
     X          + XK1I*ASDL
           DO IJK=1,3
c --- boucle IJK to improve  TTF
             TKE=TKC+DXDKI*T1KC+DXDKI*DXDKI*T2KC/2.+DXDKI**3*T3KC/6.+
     X           DXDKI**4*T4KC/24.
             T1KE=T1KC+DXDKI*T2KC+DXDKI*DXDKI*T3KC/2.+DXDKI**3*T4KC/6.
             SKE=SKC+DXDKI*S1KC+DXDKI*DXDKI*S2KC/2.+DXDKI**3*S3KC/6.+
     X           DXDKI**4*S4KC/24.
             S1KE=S1KC+DXDKI*S2KC+DXDKI*DXDKI*S3KC/2.+DXDKI**3*S4KC/6.
c  TEST*********
ccc             TKE=TKC
ccc             T1KE=T1KC
ccc             SKE=SKC
ccc             S1KE=S1KC
c **********************************************************
c --- systematic or random defaults on the field level (not for the reference)
             if(ierpf.eq.1) then
               TKE=TKE*(1.+vfield)
               SKE=SKE*(1.+vfield)
               T1KE=T1KE*(1.+vfield)
               S1KE=S1KE*(1.+vfield)
             endif
             if(ierpf.gt.1) then
               len=1
               call rlux(vecx,len)
               r1=(2.*vecx(1)-1.)*vfield
               TKE=TKE*(1.+r1)
               SKE=SKE*(1.+r1)
               T1KE=T1KE*(1.+r1)
               S1KE=S1KE*(1.+r1)
             endif
             PCRESI=0.
             DPHCI0=0.
c             if(f(9,i).ne.qst) then
c --- new crest phase is PCRESI
ccc              PCRESI=ATAN(-SKE/TKE)
ccc              DDWC=abs(f(9,i))*(TKE*COS(PCRESI)-SKE*SIN(PCRESI))
ccc              IF(DDWC.LT.0.) PCRESI=PCRESI+PI
C --- DPHCI0 is the difference between the previous crest phase and the present one
ccc              DPHCI0=PCREST-PCRESI
c             endif
c **************************************
ccc             PHIWC=PHI+DPHII-DPHCI0
c TEST********
             PHIWC=PHI+DPHII
c ***********************************
c ****             PHIWC=PHI+DPHII
             DWCI=abs(f(9,i))*(TKE*COS(PHIWC)-SKE*SIN(PHIWC))
             ENRC=f(7,i)+DWCI
             GACR=ENRC/XMAT
             BECR=SQRT(1.-1./(GACR*GACR))
             XK2II=FH0/BECR
             CXLG=abs(f(9,i))/(4.*XMAT*EQVL)
             DXDPI=RETPH-DPHCI0
c --- see Part. Acc., 1994, vol 44., pp. 215-255
C ----     relation 86
             DXDPT=DXDPI+DXK00*(1.-DKMSKE)*ASDL-DXDPI*DKMSPHI*ASDL
             DXK00=FH0*(1./BEINI-1./BEREF)
             XLH11=(YH11T+DXK00*(DKMSKE*YH1K1T+YH1K01T))*COS(DXDPT)
             XLH11=XLH11+DXDPI*YH1P1T*COS(DXDPT)
C        second derivative of Ha0(Z) (division by 2 has been made)
             XLH112=DXK00*DXK00*(H1AKIT+H1AKIMT*DKMSKE+
     X              H1AKMT*DKMSKE*DKMSKE)*COS(DXDPT)
             XLH112=XLH112+DXDPI*DXDPI*HAPIT*COS(DXDPT)
             XLH112=XLH112+DXDPI*DXK00*(YH10PKT+
     X            DKMSKE*YH11PKT)*COS(DXDPT)
C        third derivative of Ha0(Z)
             XLH113=(DXDPI**3)/3.*HAPPIT*COS(DXDPT)
             XLH11=XLH11+XLH112+XLH113
C        first derivative of Hb0(Z)
             XLH21=(YH21T+DXK00*(DKMSKE*YH2K1T+YH2K01T))*SIN(DXDPT)
             XLH21=XLH21+DXDPI*YH2P1T*SIN(DXDPT)
C        second derivative of Hb0(Z) (division by 2 has been made)
             XLH212=DXK00*DXK00*(H1BKIT+H1BKIMT*DKMSKE +
     X              H1BKMT*DKMSKE**2)*SIN(DXDPT)
             XLH212=XLH212+DXDPI*DXDPI*HBPIT*SIN(DXDPT)
             XLH212=XLH212+DXDPI*DXK00*(YH20PKT+
     X              DKMSKE*YH21PKT)*SIN(DXDPT)
C      third derivative of of Hb0(Z)
             XLH213=(DXDPI**3)/3.*HBPPIT*SIN(DXDPT)
             XLH21=XLH21+XLH212+XLH213
             XLH1I=CXLG*(XLH11-XLH21)
c --- SAUPHC is the jump of phase
             SAUPHC=FH0*XLH1I
C ---     XKMI:  AVERAGE factor k (k = frequency/velovity)
             XKMI=XK2II+SAUPHC/EQVL+(XK2II-XK1II)*ASDL/EQVL
             BEMPY=FH0/XKMI
             XK1I=XK1II-XKMI
             DXDKI=XKMI-XKM
             DPHII=(XK1II-XK2II)*EQVL/10.+(XKP1+XKP2)/120.*EQVL**2
     X            + XK1I*ASDL
           ENDDO
c TEST****       ENDDO IJK=1,3
C    Compute the shift of phase: PHSLIL
           DTS=(TKE*T1KE+SKE*S1KE)/(TKE*TKE+SKE*SKE)
           TILTAL=-4.*ATAN(DTS*3.2/EQVL)
           IF(TILTAL.NE.0.) THEN
             TIL2=TILTAL/2.
             XLREI=EQVL
             DO IIII=1,4
               FTIL=1./TAN(TIL2) - 1./TIL2 -DTS*2./XLREI
               DFTIL=-1./(SIN(TIL2)*SIN(TIL2)) + 1./(TIL2*TIL2)
               IF(DFTIL.NE.6.*0.) then
                 TIL2=TIL2 - FTIL/DFTIL
                 GX=1./TAN(TIL2) - 1./TIL2
                 XLREI=2.*DTS/GX
               endif
             ENDDO
             TILTAL=TIL2*2.
           ENDIF
           PHSLIL=TILTAL
           TEGL1=PHSLIL*PHSLIL/(SIN(PHSLIL/2.)*SIN(PHSLIL/2.))
           TEGL2=(TKE*TKE+SKE*SKE)
           TEGL=TEGL1*TEGL2/(EQVL*32.)
           TEGL=TEGL/EQVL
C ---   Transverse motion  (field E(z)**2)
           PHITTI=PHIWC+XK1I*ASDL-PHSLIP/2.
           PHITSI=PHIWC+XK2I*ASDL+PHSLIP/2.
           CETF1=f(9,i)*f(9,i)/(16.*XMAT*XMAT*EQVL*EQVL)
*  TESTsv ***************************
c OLD****           CETF1=CETF1*cort1*cort1
c **************************
           CFV1=(GINI*GINI+2.)/((GINI*GINI-1.)**2)
           CFV2=(GACR*GACR+2.)/((GACR*GACR-1.)**2)
           SV1=CFV1*(TKE*COS(PHITTI)-SKE*SIN(PHITTI))**2
           SV1=CETF1*SV1*(PHSLIL/SIN(PHSLIL/2.))**2
           SV2=CFV2*(TKE*COS(PHITSI)-SKE*SIN(PHITSI))**2
           SV2=CETF1*SV2*(PHSLIL/SIN(PHSLIL/2.))**2
C       Transverse motion  (field dE(z)/dt)
           CETI=FH0*abs(f(9,i))/(4.*XMAT*EQVL)
*  TESTsv ***************************
c OLD****           CETI=CETI*cort1
c **************************
           CXV1=(GINI*GINI-1.)**1.5
           CXV2=(GACR*GACR-1.)**1.5
           SXV1=(TKE*SIN(PHITTI)+SKE*COS(PHITTI))*PHSLIL/SIN(PHSLIL/2.)
           SXV1=-CETI*SXV1/CXV1
           SXV2=(TKE*SIN(PHITSI)+SKE*COS(PHITSI))*PHSLIL/SIN(PHSLIL/2.)
           SXV2=-CETI*SXV2/CXV2
C --- save all parameters
           SV1P(I)=SV1
           SV2P(I)=SV2
           SXV1P(I)=SXV1
           SXV2P(I)=SXV2
           DWCIS(I)=DWCI
           sauphcs(i)=SAUPHC
           BEINI1(I)=BEINI
           PHIP(I)=PHI
           TEGLP(I)=TEGL
           DXDPIP(I)=DXDPI
           DXDKIP(I)=DXDKI
           DXDPTP(I)=DXDPT
           DXK00P(I)=DXK00
           DPHIIP(I)=DPHII
C  END IF BLOCK  on APV
         ENDIF
c --- ipas = 2: the dynamics is computed over all the element
         if(ipas.eq.2) then
C    Recover all previous parameters
           SAUPHC=sauphcs(i)
           XLH1I=SAUPHC/FH0
           DWCI=DWCIS(I)
           SV1=SV1P(I)
           SV2=SV2P(I)
           SXV1=SXV1P(I)
           SXV2=SXV2P(I)
           BEINI=BEINI1(I)
           GINI = 1./SQRT(1.-BEINI*BEINI)
           PHI=PHIP(I)
           TEGL=TEGLP(I)
           DXDPI=DXDPIP(I)
           DXDKI=DXDKIP(I)
           DXDPT=DXDPTP(I)
           DXK00=DXK00P(I)
           DPHII=DPHIIP(I)
           GO TO 5678
         endif
c  compute the jump of phase at position of space charge computation
         CXLG=abs(f(9,i))/(4.*XMAT*EQVL)
C     FUNCTION Ha0(Z)
         XLH11=(YH11+DXK00*(DKMSKE*YH1K1+YH1K01))*COS(DXDPT)
         XLH11=XLH11+DXDPI*YH1P1*COS(DXDPT)
C     SECOND DERIVATIVES OF Ha0(Z) (division by 2 has been made)
         XLH112=DXK00*DXK00*(H1AKI+H1AKIM*DKMSKE+
     X          H1AKM*DKMSKE*DKMSKE)*COS(DXDPT)
         XLH112=XLH112+DXDPI*DXDPI*HAPI*COS(DXDPT)
         XLH112=XLH112+DXDPI*DXK00*(YH10PK+
     X          DKMSKE*YH11PK)*COS(DXDPT)
C      THIRD DERIVATIVE of Ha0(Z)
         XLH113=(DXDPI**3)/3.*HAPPI*COS(DXDPT)
         XLH11=XLH11+XLH112+XLH113
C     FUNCTION Hb0(Z)
         XLH21=(YH21+DXK00*(DKMSKE*YH2K1+YH2K01))*SIN(DXDPT)
         XLH21=XLH21+DXDPI*YH2P1*SIN(DXDPT)
C     SECOND DERIVATIVE OF Hb0(Z) (division by 2 has been made)
         XLH212=DXK00*DXK00*(H1BKI +H1BKIM*DKMSKE+
     X          H1BKM*DKMSKE**2)*SIN(DXDPT)
         XLH212=XLH212+DXDPI*DXDPI*HBPI*SIN(DXDPT)
         XLH212=XLH212+DXDPI*DXK00*(YH20PK+
     X          DKMSKE*YH21PK)*SIN(DXDPT)
C      THIRD DERIVATIVE OF Hb0(Z)
         XLH213=(DXDPI**3)/3.*HBPPI*SIN(DXDPT)
         XLH21=XLH21+XLH212+XLH213
         XLH1I=CXLG*(XLH11-XLH21)
         SAUPHC=FH0*XLH1I
5678     CONTINUE
         CXLG=abs(f(9,i))/(4.*XMAT*EQVL)
         XLH01=(YH10+DXDKI*YH1K0+DXK00*YH1K00)*COS(DXDPT)
         XLH02=(YH20+DXDKI*YH2K0+DXK00*YH2K00)*SIN(DXDPT)
         XLH0I=CXLG*(XLH01-XLH02)
C -- COUPLING INTEGRALS
         XLP11=(YP11+ DXDKI*YP1K1 + DXK00*YP1K01)*COS(DXDPT)
         XLP21=(YP21+ DXDKI*YP2K1 + DXK00*YP2K01)*SIN(DXDPT)
         XLP1I=CXLG*(XLP11-XLP21)
         XLP12=(YP12+ DXDKI*YP1K2 + DXK00*YP1K02)*COS(DXDPT)
         XLP22=(YP22+ DXDKI*YP2K2 + DXK00*YP2K02)*SIN(DXDPT)
         XLP2I=CXLG*(XLP12-XLP22)
         RADIU=SQRT(F(2,I)*F(2,I)+F(4,I)*F(4,I))
         IF(RADIU.LT.1.e-06) THEN
           DRADIU=.001*SQRT(F(3,I)*F(3,I)+F(5,I)*F(5,I))
         ELSE
           DRADIU=F(3,I)*.001*F(2,I)/RADIU+F(4,I)*F(5,I)*.001/RADIU
         ENDIF
         RP=RADIU
         RPP=DRADIU
C        Picht coordinates
         RRP=RP*SQRT(BEINI*GINI)
         RRPP=RPP*SQRT(BEINI*GINI)
c --- ipas = 1: the dynamics is computed until the position of space charge computation
C       gain of energy at the position of space charge computation (not for the reference)
         if(ipas.eq.1) then
c *****           PHIWC=PHI+DPHII-DPHCI0
           PHIWC=PHI+DPHII
c       istm: flag ==> the energy is compute at the position of space charge computation
           istm=ist-1
           GACR=GAMCI(PHIWC,PCRESI,GINI,ISTM,abs(f(9,i)))
           DWCI=(GACR-GINI)*XMAT
         endif
         DWPI=DWCI+XMAT*FH0*FH0*RRP*RRP/4.*XLH0I +
     X        XMAT*RRP*RRPP*FH0*FH0/2. *XLH1I
C       PHASE JUMP
         DELPHI=SAUPHC+ FH0**3*RRP*RRP/4. *XLP1I +
     X          FH0**3*RRP*RRPP/2. *XLP2I
         IF(IPAS.EQ.2) THEN
c    kicks of energy from space charge effects are in dwp(i)
           IF(ICHAES) THEN
             F(7,i)=F(7,I)+DWPI+DWP(I)
             GAMSOR = F(7,i)/XMAT
             IF(GAMSOR.LE.1.) f(8,i)=0.
C     particle is lost
             IF(F(8,i).EQ.0.) go to 19
             BESOR=SQRT(1. -1./(GAMSOR*GAMSOR))
C       phase jump resulting from space charge
             DELGAM=DWP(I)/XMAT
             GAMKK0=F(7,I)/XMAT
             BEKK0=SQRT(1. -1./(GAMKK0*GAMKK0))
             DBEK21=DELGAM/(BEKK0**3 * GAMKK0**3)
             DELSC=FH0*SCDIST*DBEK21/2.
             DELPHI=DELPHI+DELSC
             DITEMP=YLG/(BESOR*VL)+DELPHI/FH
             f(6,i)=fs(6,i)+ditemp
           ELSE
             f(7,i)=f(7,i)+DWPI
             GAMSOR = f(7,i)/XMAT
             IF(GAMSOR.LE.1.) f(8,i)=0.
C     particle is lost
             if(f(8,i).eq.0.) then
               f6i=f(6,i)-tcog
               write(16,3928) i,int(f(1,i)),f(2,i),f(3,i),f(4,i),
     *               f(5,i),f6i*fh*180./pi,f(7,i)-xmat,int(f(9,i))
               if(ilost.ge.ngood) stop
               ilost=ilost+1
               go to 19
             ENDIF
             BESOR=SQRT(1. -1./(GAMSOR*GAMSOR))
             DITEMP=YLG/(besor*vl)+DELPHI/FH
             f(6,i)=f(6,i)+DITEMP
           ENDIF
         ENDIF
         IF(ipas.eq.1) THEN
           f(7,i)=f(7,i)+DWPI
           GAMSOR = f(7,i)/XMAT
           if(GAMSOR.LE.1.) f(8,i)=0.
           if(f(8,i).eq.0.) then
             write(16,3928) i,int(f(1,i)),f(2,i),f(3,i),f(4,i),
     *            f(5,i),f(6,i)*fh*180./pi,f(7,i)-xmat,int(f(9,i))
             if(ilost.ge.ngood) stop
             ilost=ilost+1
             go to 19
           endif
           BESOR=SQRT(1. -1./(GAMSOR*GAMSOR))
c   TEST******
c           DITEMP=YLG/(BESOR*VL)+DELPHI/FH
           DITEMP=YLG/(2.*BESOR*VL)+DELPHI/FH
c2014-Aug-08           DITEMP=ylg/(BESOR*VL)+DELPHI/FH
ccc             DITEMP=(2.*XK1II*ASDL+XKMI*(EQVL-ASDL)+
ccc     x              XK2II*(YLG-(EQVL+ASDL)))/FH+DELPHI/FH
c ****************************
           f(6,i)=f(6,i)+DITEMP
         ENDIF
         IF(I.EQ.ICONT) THEN
           IF(IPAS.EQ.2) THEN
             WRITE(16,*) '* DYNAMICS AT THE OUTPUT :'
             WRITE(16,994) DELPHI*180./PI,DELSC*180./PI,DWP(I)
994          FORMAT(1x,'* PHASE JUMP ',E12.5,' DEG  CORRECTED BY :',
     X       E12.5,' DEG ',' SC KICK(MEV) ',E12.5)
             ENRPRIN=f(7,i)-XMAT
           WRITE(16,88) DWPI,ENRPRIN,BESOR,f(9,i),DITEMP
88           FORMAT(1x,'* ENERGY GAIN : ',E14.7,' MEV','   ENERGY :',
     X       E14.7,' MEV',/,1X,'* BETA :',E12.5,/,
     X       1X,'* CHARGE :',F5.0,' TRANSIT TIME :',E12.5,' SEC',/,
     X       1X,'*')
            write(16,*) 'TK TKE ',TK,TKE
            write(16,*) 'SK SKE ',SK,SKE
            write(16,*) 'PCREST  PCRESI',PCREST,PCRESI
           ENDIF
           if(ipas.eq.1) then
             WRITE(16,*) '* DYNAMICS AT THE POSITION OF SPACE CHARGE :'
             ENRPRIN=f(7,i)-XMAT
             WRITE(16,9944) DELPHI*180./PI,DWPI,ENRPRIN,DITEMP
9944         FORMAT(1x,'* PHASE JUMP ',E12.5,' DEG ',' ENERGY GAIN : ',
     X       E14.7,' MeV',' ENERGY :',E14.7,' MEV',/,
     X       ' TRANSIT TIME :',E12.5,' SEC',/,
     X       1X,'*')
           ENDIF
         ENDIF
         AMORT=SQRT(BEINI*GINI/(BESOR*GAMSOR))
C     INTEGRALS OF THE TRANSVERSE MOTION (not computed for the reference)
C     ******************
c ******       if(i.ne.ngdrf) then
C ----  integrals of E(Z)**2
C ----     (G**2+2/(G**2-1)**2)*E(z)**2
         CETF=f(9,i)*f(9,i)*SQCTTF*SQCTTF/(16.*XMAT*XMAT*EQVL*EQVL)
c TESTsv ********
c OLD****         CETF=CETF*cort1*cort1
c *************************************
         XJF0I=YFS0+DXDKI*(YFSK0+YFSCK0) + DXDPI*(YFSP0 +YFSCP0) +
     X         DXK00*YFSKC0
         XJF0I=XJF0I*CETF
         XJF1I=YFS1+DXK00*DKMSKE*(YFSK1+YFSCK1)+DXDPI*(YFSP1+YFSCP1)+
     X         DXK00*YFSKC1
         XJF1I=XJF1I*CETF
         XJF2I=YFS2+DXK00*DKMSKE*(YFSK2+YFSCK2)+DXDPI*(YFSP2+YFSCP2)+
     X         DXK00*YFSKC2
         XJF2I=XJF2I*CETF
         V1=SV1
         V2=SV2
C ----     (G**2+2/(G**2-1)**2)
         CETM=f(9,i)*f(9,i)*TEGL/(XMAT*XMAT)
c TESTsv ********
c OLD****         CETM=CETM*cort1*cort1
c *************************************
         XMN0I=YNS0+DXDKI*YNSK0+DXDPI*YNSP0+DXK00*YNSK0C
         XMN0I=XMN0I*CETM
         XMN1I=YNS1+DXK00*DKMSKE*YNSK1+DXDPI*YNSP1+DXK00*YNSK1C
         XMN1I=XMN1I*CETM
         XMN2I=YNS2+DXK00*DKMSKE*YNSK2 + DXDPI*YNSP2 +DXK00*YNSK2C
         XMN2I=XMN2I*CETM
         F0=XJF0I + XMN0I
         F1=XJF1I + XMN1I
         F2=XJF2I + XMN2I
C ---     1/(B*G)**3 *dE(z)/dt
         CETI=FH0*abs(f(9,i))/(8.*XMAT*EQVL)
*  TESTsv ***************************
c OLD****         CETI=CETI*cort1
c **************************
         XIE01= (YE10+DXDKI*YE1K0+DXK00*YE1KC0)*COS(DXDPT)
         XIE02= (YE20+DXDKI*YE2K0+DXK00*YE2KC0)*SIN(DXDPT)
         XIE0I=-CETI*(XIE01+XIE02)
         XIE11= (YE11+DXDKI*YE1K1+DXK00*YE1KC1)*COS(DXDPT)
         XIE12= (YE21+DXDKI*YE2K1+DXK00*YE2KC1)*SIN(DXDPT)
         XIE1I=-CETI*(XIE11+XIE12)
         XIE21= (YE12+DXDKI*YE1K2+DXK00*YE1KC2)*COS(DXDPT)
         XIE22= (YE22+DXDKI*YE2K2+DXK00*YE2KC2)*SIN(DXDPT)
         XIE2I=-CETI*(XIE21+XIE22)
         XV1=SXV1
         XV2=SXV2
         XI0=XIE0I
         XI1=XIE1I
         XI2=XIE2I
C      Transport matrix in 'PICHT' coordonates
         XQ0=XI0-F0
         XQ1=XI1-F1
         XQ2=XI2-F2
         XQ01=(XQ1+ASDL*XQ0)
         XQ12=(XQ2+ASDL*XQ1)
         V1=XV1-V1
         V2=XV2-V2
         A11=-XQ01*(1.+ (V1+V2)*EQVL*EQVL/120.)
         A12=-(XQ2+2.*ASDL*XQ1+ASDL*ASDL*XQ0 +
     X       EQVL*EQVL*((V1+V2)/120. + EQVL*V2/120.)*XQ01)
         ZA=-(XQ12/EQVL +V2*EQVL*EQVL*XQ01/120.)
         ZB=-((EQVL+ASDL)*XQ12/EQVL -
     X       EQVL*XQ01/10. +V2*(EQVL+ASDL)*EQVL*EQVL*XQ01/120.)
         A21=XQ0*(1.+ (V1+V2)*EQVL*EQVL/120.)
         A22= XQ1 + XQ0*(ASDL+ASDL*EQVL*EQVL*(V1+V2)/120. +
     X       EQVL**3*V2/120.)
         ZC=XQ1/EQVL + V2*EQVL*EQVL*XQ0/120.
         ZD= (ASDL+EQVL)*XQ1/EQVL -
     X      ( EQVL/10. + V2*(EQVL+ASDL)*EQVL*EQVL/120.)*XQ0
C
         TMA=1./(1.-ZA-ZC*ZB/(1.-ZD))
         T11= (A11+ZB*A21/(1.-ZD)) *TMA
         T12=( A12+ZB*A22/(1.-ZD) )*TMA
         T21=(A21+ZC*T11)/(1.-ZD)
         T22=(A22+ZC*T12)/(1.-ZD)
         VR11=(1.+T11+DCUM*T21)
         VR12=(T12+DCUM*(1.+T22))
         VR21=T21
         VR22=1.+T22
         IF(I.EQ.ICONT) THEN
           DETRE=VR11*VR22-VR12*VR21
           WRITE(16,8921) VR11,VR12,VR21,VR22,DETRE
8921       FORMAT(2X,' TRANSVERSE CANONICAL MATRIX:(cm,radian) ',/,
     X     2X,' VR11:',E12.5,'   VR12:',E12.5,/,
     X     2X,' VR21:',E12.5,'   VR22:',E12.5,/,
     X     2X,' DETERMINANT       :',E12.5,//)
         ENDIF
C       REAL MATRIX
         A11=VR11*AMORT
         A12=VR12*AMORT
         A21=VR21*AMORT
         A22=VR22*AMORT
Comment       ** cumulative matrix (particle 1)
Comment          tables R(,) ET T(,,,)
comment         IF (I .EQ. 1.AND.IPAS.EQ.2) THEN
comment           STA11=VR11*SA11+VR12*SA21
comment           STA12=VR11*SA12+VR12*SA22
comment           STA21=VR21*SA11+VR22*SA21
comment           STA22=VR21*SA12+VR22*SA22
comment           SA11=STA11
comment           SA12=STA12
comment           SA21=STA21
comment           SA22=STA22
comment           SAA11=SA11*AMORT
comment           SAA12=SA12*AMORT
comment           SAA21=SA21*AMORT
comment           SAA22=SA22*AMORT
comment           RS(1,1) =SAA11
comment           RS(1,2) =SAA12
comment           RS(2,1) =SAA21
comment           RS(2,2) =SAA22
comment           RS(3,3) =SAA11
comment           RS(3,4) =SAA12
comment           RS(4,3) =SAA21
comment           RS(4,4) =SAA22
comment           DO IA=1,6
comment             DO IB=1,6
comment               R(IA,IB)=RCUL(IA,IB)
comment             ENDDO
comment           ENDDO
comment           CALL MFORDRE(RCUL,RS,R)
comment         ENDIF
C      BEAM COORDINATES
         FXT1=A11*F(2,I)+A12*F(3,I)*1.E-03
         FXT2=A21*F(2,I)+A22*F(3,I)*1.E-03
         FXT3=A11*F(4,I)+A12*F(5,I)*1.E-03
         FXT4=A21*F(4,I)+A22*F(5,I)*1.E-03
         F(2,I)=FXT1
         F(3,I)=FXT2*1.E03
         F(4,I)=FXT3
         F(5,I)=FXT4*1.E03
         IF(I.EQ.ICONT) THEN
C        ** CURRENT MATRIX
           STTA11=VR11*SACT11+VR12*SACT21
           STTA12=VR11*SACT12+VR12*SACT22
           STTA21=VR21*SACT11+VR22*SACT21
           STTA22=VR21*SACT12+VR22*SACT22
           SACT11=STTA11
           SACT12=STTA12
           SACT21=STTA21
           SACT22=STTA22
           SAA11=SACT11*AMORT
           SAA12=SACT12*AMORT
           SAA21=SACT21*AMORT
           SAA22=SACT22*AMORT
           DET=A11*A22-A12*A21
           WRITE(16,992) A11,A12*1.E-3,A21*1.E3,A22,DET,AMORT
992        FORMAT(1x,'*  TRANSVERSE MATRIX (cm,mrd)',/,1x,
     X           '*',E12.5,3X,E12.5,/,1x,
     X           '*',E12.5,3X,E12.5,/,1x,
     X           '* DETERMINANT :',E12.5,' DUMPING OF ENERGY :',E12.5)
           write(16,*)'*'
           WRITE(16,559) F(2,I),F(3,I),F(4,I),F(5,I)
559        FORMAT(
     X     ' * TRANVERSE COORDINATES AT OUTPUT  ',/,1x,
     X     '* X :',E12.5,' CM  XP :',E12.5,' MRD ',/,1x,
     X     '* Y :',E12.5,' CM  YP :',E12.5,' MRD')
           IF(IPAS.EQ.2)
     X       WRITE(16,*) '********** END OF FOLLOWED PARTICLE ********'
         ENDIF
c *****        endif
19     CONTINUE
       return
       end
       SUBROUTINE bunparm(v,dp,harm,prlim)
c   ............................................................................
C          BUNCHER  (NO SPACE CHARGE EFFECT)
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/RIGID/BORO
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/corec/tref1
       COMMON/QMOYEN/QMOY
       common/aerp/vphase,vfield,ierpf
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SHIF/DTIPH,SHIFT
       common/tofev/ttvols
       common/fcont/ifcont
       logical ifcont
       character*1 cr
       dimension vecx(1)
       LOGICAL chasit,itvol,imamin,shift
C      ENVELOPE
       call stapl(davtot*10.)
cxx       ilost=0
       twopi=2.*pi
       freq=fh/twopi
       wavel=vl/freq
       fcpi=fh*180./pi
c print out on terminal of transport element # on one and the same line
       nrbunc=nrbunc+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) nrtre,nrbunc,cr
8254   format('Transport element:',i5,
     *        '      Buncher             :',i5,a1,$)
#else
       WRITE(6,8254,advance='no') nrtre,nrbunc,cr
8254   format('Transport element:',i5,
     *        '      Buncher             :',i5,a1)
#endif
       if (harm.le.0.) harm=1.
       rhs=prlim*prlim
c   test window
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       tcog=0.
       bcog=0.
       do np=1,ngood
         tcog=tcog+f(6,np)
         gpa=f(7,np)/xmat
         bcog=sqrt(1.-1./(gpa*gpa))+bcog
       enddo
       tcog=tcog/float(ngood)
       bcog=bcog/float(ngood)
       gcog=1./sqrt(1.-bcog*bcog)
       encog=xmat*gcog-xmat
c adjustement of the phase of RF w.r.t. the T.O.F.
      xkpi=0.
      if(imamin) then
       ttvpi=harm*ttvols*fcpi
       xkpi=ttvpi/360.
       ixkpi=int(xkpi)
       xkpi=(xkpi-float(ixkpi))*360.
       write(16,*) ' *** TOF correction:',-xkpi,' deg'
       dp=dp-xkpi*pi/180.
       write(16,*)' ***phase of RF adjusted : ',dp*180./pi,' deg'
      endif
c 20/08/2009    delay of phase of the reference at input w.r.t. the synchronous phase
c 20/08/2009       ttvpi=harm*ttvols*fcpi
c 20/08/2009       xkpi=ttvpi/180.
c 20/08/2009       ixkpi=xkpi+0.5
c 20/08/2009       ixkpi=ixkpi*180
c 20/08/2009       xkpc=cos(ixkpi*pi/180.)
c 20/08/2009       xkpi=ttvpi-float(ixkpi)
c 20/08/2009       if(itvol.and.imamin)
c 20/08/2009     *  write(16,8975) dp*180./pi,xkpi
c 20/08/2009 8975   format('***previous phase offset: ',e12.5,' deg',/,
c 20/08/2009     * '***new phase offset: ',e12.5,' deg')
c  start of write to file '.short' for buncher
       idav=idav+1
       iitem(idav)=8
       dav1(idav,1)=v
       dav1(idav,2)=dp*180./pi
       dav1(idav,3)=prlim
       dav1(idav,4)=davtot*10.
       if(itvol) dav1(idav,5)=-xkpi
c  end
       WRITE(16,178)
178    FORMAT(/,' Dynamics at the input',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       WRITE(16,1788) bcog,gcog,encog,tcog*fcpi,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       e0t=harm*v/(bcog*wavel)
c    random or systematic error on  RF level
       if(ierpf.gt.0) then
c  systematic error on RF level
         if(ierpf.eq.1) e0t=e0t*(1.+vfield)
c  random error on RF level
         if(ierpf.ge.2) then
           len=1
           call rlux(vecx,len)
           r1=(2.*vecx(1)-1.)*vfield
           e0t=e0t*(1.+r1)
         endif
       endif
       cay=harm*twopi/(bcog*gcog*wavel)
       caysq=cay**2
       con=twopi*e0t*qmoy/xmat
       rad=pi/180.
c  systematic error on  phase
       if(ierpf.eq.1) dp=dp+vphase*rad
c  random error on phase
       if(ierpf.ge.2) then
         len=1
         call rlux(vecx,len)
         r1=(2.*vecx(1)-1.)*vphase*rad
         dp=dp+r1
       endif
c shift=true => reference and COG seperated, otherwise reference=COG
c --- save the reference
       ovref=vref
       otref=tref
c --- shift = false: the reference particle is the cog
       if(shift) then
         ovref=vref
         beref=vref/vl
         gamref=1./sqrt(1.-beref*beref)
         older=xmat*gamref
       else
         tref=tcog
         vref=bcog*vl
         ovref=vref
         beref=bcog
         gamref=1./sqrt(1.-beref*beref)
         older=xmat*gamref
       endif
c ---  if imamin = false: phase setting has been forced equal to dp, otherwise phase setting has been adjusted
       dgr=v*cos(harm*ttvols*fh+dp)*qmoy
       ewer=older+dgr
       gor=ewer/xmat
       vref=vl*sqrt(1.-1./(gor*gor))
       ENRPRIN=older-xmat
       WRITE(16,165) beref,gamref,ENRPRIN,tref*fh*180./pi,tref
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       wsync=0.
       bcour=0.
       tcog=0.
       do np=1,ngood
         rs=f(2,np)**2+f(4,np)**2
         a=harm*(f(6,np)-tref+ttvols)*fh+dp
         rs=rs*1.e-04
c 20/08/2009         s=sin(a)*xkpc
         s=sin(a)
         w=f(7,np)-xmat
         bg=sqrt(w/xmat*(2.+w/xmat))
c    conversion cm-->m and mrd--->rd
         f2=f(2,np)
         f4=f(4,np)
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
         bgx=f3*bg
         bgy=f5*bg
         cayrsq=caysq*rs
         arg=cayrsq/4.
         xi0=1.+arg*(1.+arg*(.25+arg/36.))
c           dw=v*cos(a)*xi0*qmoy
         dw=v*cos(a)*xi0*f(9,np)
c 20/08/2009         dw=dw*xkpc
         wb=w+.5*dw
         bgav=sqrt(wb/xmat*(2.+wb/xmat))
         gav=1.+wb/xmat
         bav=bgav/gav
         bcour=bcour+bav
         wf=w+dw
         bgf=sqrt(wf/xmat*(2.+wf/xmat))
         xi1okr=.5+.25*arg+arg**2/24.
         del=-con*s*(1.-bav*bcog)*xi1okr/bav
         tcog=tcog+f(6,np)
         f3=(bgx+del*f2)/bgf
         f(3,np)=f3*1.e03
         f5=(bgy+del*f4)/bgf
         f(5,np)=f5*1.e03
         f(7,np)=wf+xmat
         wsync=wsync+wf
       enddo
       wsync=wsync/float(ngood)
       bcour=bcour/float(ngood)
       tcog=tcog/float(ngood)
c Test window
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c new output
c shift=true => reference and COG seperated, otherwise reference=COG
       engain=wsync-encog
       if(shift) then
         beref=vref/vl
         gamref=1./sqrt(1.-beref*beref)
         enref=ewer-xmat
       else
         tref=tcog
         vref=bcour*vl
         beref=bcour
         enref=wsync
         dgr=engain
       endif
       if(itvol) ttvols=tref
       WRITE(16,3777)
3777   FORMAT(/,' Dynamics at the output',/,
     1 5X,'   BETA     dW(MeV)    ENERGY(MeV) ',
     2 '   TOF(deg)     TOF(sec)')
       WRITE(16,3473) beref,dgr,enref,fh*TREF*180./PI,TREF
3473   FORMAT(' REF ',f7.5,3x,f10.6,3x,f8.3,3x,e12.5,3x,e12.5)
       WRITE(16,1789) bcour,engain,wsync,tcog*fh*180./pi,tcog
1789   FORMAT(' COG ',f7.5,3x,f10.6,3x,f8.3,3x,e12.5,3x,e12.5)
c dave start for buncher
       dav1(idav,36)=ngood
c as beam passed through buncher, assume bunched beam
       ifcont=.false.
       WRITE(16,*) 'After buncher, bunched beam assumed'
       call emiprt(0)
c dave end
cxx3928   FORMAT(' # ',i5,1x,7(f10.2,1x))
       return
       end
       SUBROUTINE refer
c   ............................................................................
c  change the longitudinal position of the reference
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/newref/dephas,dewref,iref,irefw
       common/tapes/in,ifile,meta
       common/itvole/itvol,imamin
       common/tofev/ttvols
       logical itvol,imamin
c*************************************************************************************************
c ---- IREF =0 and IREFW = 0:  dewref is dW/W where W is the kinetic energy of the old reference
c ---- IREF =0 and IREFW = 1:  dewref is dW (MeV) relative to the kinetic energy of the old reference
c ---- IREF =0 and IREFW = 2:  dewref is dW (MeV) is new reference energy and dbref new reference phase (deg) in abs. units
c ---- IREF =1 and IREFW = 0:  dewref is dW/W where W is the kinetic energy of the c.o.g
c ---- IREF =1 and IREFW = 1:  dewref is dW (MeV) relative to the kinetic energy of the c.o.g.
c*************************************************************************************************
c ---- save  reference
      avref=vref
      atref=tref
      attvols=ttvols
      if(irefw.eq.2) then
        wnref=dewref
        gnref=wnref/xmat+1.
        bref=sqrt(gnref*gnref-1.)/gnref
c ---- set tref in sec based on dephas (deg)
        delt=dephas*pi/(fh*180.)
c ------ the new reference is:
        tref=delt
        vref=bref*vl
        if(itvol)ttvols=tref
      endif
      if(iref.eq.0) then
       if(irefw.eq.0) then
        bref=vref/vl
        gref=1./sqrt(1-bref*bref)
        dbref=bref*dewref/(gref*(gref+1.))
c ---- change dephas (deg) in delt (sec)
        delt=dephas*pi/(fh*180.)
c ------ the new reference is:
         tref=tref+delt
         vref=vref+dbref*vl
         if(itvol)ttvols=tref
        endif
        if(irefw.eq.1) then
         bref=vref/vl
         gref=1./sqrt(1.-bref*bref)
         wref=(gref-1.)*xmat
         wnref=wref+dewref
         gnref=wnref/xmat+1.
         bref=sqrt(gnref*gnref-1.)/gnref
c ---- change dephas (deg) in delt (sec)
         delt=dephas*pi/(fh*180.)
c ------ the new reference is:
         tref=tref+delt
         vref=bref*vl
         if(itvol)ttvols=tref
        endif
       endif
      if(iref.eq.1) then
c  ----- c.o.g. of the bunch is the new reference
        tcog=0.
        bcog=0.
        do i=1,ngood
         tcog=tcog+f(6,i)
         gpai=f(7,i)/xmat
c-bugfix-R10-2013-11-12
c         bcog=bbref+sqrt(1.-1./(gpai*gpai))
         bcog=bcog+sqrt(1.-1./(gpai*gpai))
        enddo
        tcog=tcog/float(ngood)
        bcog=bcog/float(ngood)
        gcog=1./sqrt(1-bcog*bcog)
        wcog=(gcog-1.)*xmat
c ---- change dephas (deg) in delt (sec)
        delt=dephas*pi/(fh*180.)
        if(irefw.eq.0) then
c-bugfix-R10-2013-11-12
c          dbcog=bcog*dewref/(gcog*(gcog+1.))
          wrefn=wcog+wcog*dewref/100.
c ------ the new reference is:
          tref=tcog+delt
c-bugfix-R10-2013-11-12
c          vref=(bcog+dbcog)*vl
           grefn=wrefn/xmat+1.
           vref=vl*sqrt(grefn*grefn-1.)/grefn
          if(itvol)ttvols=tref
        endif
        if(irefw.eq.1) then
          wncog=wcog+dewref
          gncog=wncog/xmat+1.
          bcog=sqrt(gncog*gncog-1.)/gncog
c ------ the new reference is:
          tref=tcog+delt
          vref=bcog*vl
          if(itvol)ttvols=tref
        endif
      endif
      baref=avref/vl
      garef=1./sqrt(1.-baref*baref)
      waref=(garef-1.)*xmat
      bnref=vref/vl
      gnref=1./sqrt(1.-bnref*bnref)
      wnref=(gnref-1.)*xmat
      fcpi=fh*180/pi
      write(16,20)atref*fcpi,attvols*fcpi,waref
20    format(3x,'**before NREF',/,5x,'tof of the reference: ',e12.5,
     *' deg tof for adjustments: ',e12.5,' deg energy of reference: '
     * ,e12.5,' MeV')
      write(16,21)tref*fcpi,ttvols*fcpi,wnref
21    format(3x,'**after NREF',/,5x,'tof of the reference: ',e12.5,
     *' deg tof for adjustments: ',e12.5,' deg energy of reference: '
     * ,e12.5,' MeV')
       RETURN
       END
       SUBROUTINE steer(fld,nvf)
       IMPLICIT REAL*8 (A-H,O-Z)
C---- TRANSFORM BEAM THRU THIN STEERER
C     4/15/14 - Daniel Alt: Added electrostatic steerers.
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       character*1 cr
       logical iesp,ichaes
C---- MAGNETIC STEERER
C---- PARAMETERS: fld (Tm), nvf
c     where fld is the integrated field
c     if nvf=0,  horizontal magnetic steerer
c     if nvf=1,  vertical magnetic steerer
C
C---- ELECTROSTATIC STEERER
c     This is a zero length element: the length shown below is
c     ONLY used for kick calculation.
C---- PARAMETERS: fld (kV*m/m), nvf
c     where fld is Plate Voltage * plate length/plate separation.
C     if nvf=2, horziontal electrostatic steerer
C     if nvf=3, vertical electrostatic steerer
C
C-----ANGULAR DISPLACEMENTS DUE TO STEERING
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#if F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *       '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,advance='no') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *       '      Accelerating element:',i5,a1)
#endif
       if (nvf.eq.0) then
c horizontal magnetic steerer
         write(16,*) 'Horizontal magnetic steerer: ',fld,' Tm'
         do i=1,ngood
           const=xmat*1.e8/(f(9,i)*vl)
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           DISPX=fld/(const*gpai*bpai)*1000.
           f(3,i)=f(3,i)+dispx
         enddo
       else if (nvf.eq.1) then
c vertical magnetic steerer
         write(16,*) 'Vertical magnetic steerer: ',fld,' Tm'
         do i=1,ngood
           const=xmat*1.e8/(f(9,i)*vl)
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           DISPY=fld/(const*gpai*bpai)*1000.
           f(5,i)=f(5,i)+dispy
         enddo
       else if (nvf.eq.2) then
C horizontal electrostatic steerer
         write (16,*) 'Horizontal electrostatic steerer: ',fld,' kV*m/m'
         do i=1,ngood
           gpai=f(7,i)/xmat
           const=(gpai/(gpai*gpai-1.))*f(9,i)
           DISPX=const*fld/xmat
           f(3,i)=f(3,i)+dispx
         enddo
       else if (nvf.eq.3) then
C vertical electrostatic steerer
         write (16,*) 'Vertical electrostatic steerer: ',fld,' kV*m/m'
         do i=1,ngood
           gpai=f(7,i)/xmat
           const=(gpai*gpai/(gpai*gpai-1.))*f(9,i)
           DISPY=const*fld/xmat
           f(5,i)=f(5,i)+dispy
         enddo
       else
c error on input
         write(16,*) 'Wrong value for NVF in STEERER'
         STOP
       endif
       RETURN
       END
       SUBROUTINE emiprt(L)
c  ......................................................................
c  following EMIT or EMITL card, store data in arrays to be printed by
c  subroutine "daves"
c  L=0 corresponds to EMIT (do not read and print a label)
c  L=1 corresponds to EMITL (do read and print a label in dynac.short)
c  looks for the statistics with EXT2D and returns them in array dav1 and dav2
C      IDCH = 1:WITH  CHASE TEST
C      IDCH NE 1:OTHERWISE
c
c       cog(1) : Energy(MeV)
c       cog(3) : t.o.f. (sec)
c       cog(4) : x-direction (cm)
c       cog(5) : xp(mrd)
c       cog(6) : y-direction (cm)
c       cog(7) : yp(mrd)
c
c
c  exten(1) : Sum( dE*dE )  MeV*MeV
c  exten(2) : Sum( dE*dPHase ) MeV*rad
c  exten(3) : Sum( dPHase*dPHase ) rad*rad
c  exten(4) : Sum( x*x )   cm*cm
c  exten(5) : Sum( xp*xp )   mrad*mrad
c  exten(6) : Sum( y*y )   cm*cm
c  exten(7) : Sum( yp*yp )  mrad*mrad
c  exten(8) : Sum( x*xp )   cm*mrad
c  exten(9) : Sum( y*yp )   cm*mrad
c
c  ......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/etchas/fractx,fracty,fractl
       common/dyn/tref,vref
       common/faisc/f(10,iptsz),imax,ngood
       common/davcom/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/davprt/shortl
       common/qmoyen/qmoy
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/etcha3/ichxyz(iptsz)
       common/speda/dave,idave
       common/cptemit/xltot(maxcell1),nbemit
       common /consta/ vl, pi, xmat, rpel, qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/mcs/imcs,ncstat,cstat(20)
       common/strip/atm,qs,atms,ths,qop,sqst(6),anp,nqst
       common/shortl/davprt
       logical dave,chasit
       dimension foo(20,9),NDP(20)
       character*80 davprt(maxcell1),shortl
       nbemit=nbemit+1
       do i=1,ngood
         ichxyz(i)=1
       enddo
       xltot(nbemit)=davtot
       idav=idav+1
       iitem(idav)=3
       dav1(idav,40)=fh
       do i=1,ngood
         ichas(i)=1
       enddo
       iarg=1
       if (L.eq.1) davprt(idav)=shortl
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       CALL EXT2D(1)
c --- qdisp : average dispersion (MeV)
c --- sqmdv:  emittance (MeV*rad)
       qdisp=2.*sqrt(exten(1))
       qmdv=exten(1)*exten(3)-exten(2)*exten(2)
       sqmdv=4.*pi*sqrt(qmdv)
c --- qdp : average extension in phase (rad)
       qdp=2.*sqrt(exten(3))
c --- cor12: coefficient of correlation in (dE, dPHI)
       cor12=exten(2)/sqrt(exten(1)*exten(3))
c sup       pent12=sqrt(exten(1)/exten(3))/cor12
c sup       pent21=sqrt(exten(1)/exten(3))*cor12
c sup       qdpde=qdp*180./pi
c ---  particle reference
c      dav1(idav,3): relativistic beta
c      dav1(idav,4): Kinetic energy (MeV)
c      dav1(idav,5): phase (in deg. w.r.t. k*pi)
        beref=vref/vl
        gref=1./sqrt(1.-beref*beref)
        dav1(idav,3)=beref
        dav1(idav,4)=xmat*(gref-1.)
        dav1(idav,5)=-(int(tref*fh/pi+0.5)-tref*fh/pi)*180.
c*et*2013*change to print out absolute TOF
        dav1(idav,5)=(tref*fh/pi)*180.
c ---  c.o.g of the bunch
c      dav1(idav,6): Kinetic energy (MeV)
c      dav1(idav,7): phase (in deg. w.r.t. k*pi)
       phnw=-(int(tcog*fh/pi+0.5)-tcog*fh/pi)*180.
c*et*2013*change to print out absolute TOF
       phnw=(tcog*fh/pi)*180.
       dav1(idav,6)=encog-xmat
       dav1(idav,7)=phnw
c ---- deviation between the fictious reference and the c.o.g. of the bunch
c     dav1(idav,8) : deviation in energy (MeV)
c     dav1(idav,9) : deviation of phase (deg)
       dav1(idav,8)=encog-xmat-dav1(idav,4)
       dav1(idav,9)=phnw-dav1(idav,5)
c ---- statistics in z-zp
c --- dav1(idav,10) :  extension dPHI (deg)
c --- dav1(idav,11) :  dispersion dE (MeV)
c --- dav1(idav,12) :  emittance (MeV*rad)
c --- dav1(idav,23) :  correlation in  between dE an dPHI
       dav1(idav,10)=qdp*180./pi
       dav1(idav,11)=qdisp
       dav1(idav,12)=sqmdv/pi
       dav1(idav,23)=cor12
c sup       dav1(idav,39)=sqmdv*180./(pi*pi)
c sup       Ez(ns.keV)
c sup       dav1(idav,12)=sqmdv*(1.E09)/(pi*fh)
c ---- statistics in x-xp and y-yp
       trqtx=exten(4)*exten(5)-exten(8)*exten(8)
       trqpy=exten(6)*exten(7)-exten(9)*exten(9)
       surxth=4.*pi*sqrt(trqtx)
       suryph=4.*pi*sqrt(trqpy)
       qditax=2.*sqrt(exten(4))
       qdiant=2.*sqrt(exten(5))
       qditay=2.*sqrt(exten(6))
       qdianp=2.*sqrt(exten(7))
c  dav1(idav,13): extension in x (mm)
c  dav1(idav,14): extension in xp (mrad)
c --- dav1(idav,15): correlation between x and xp
       dav1(idav,13)=qditax*10.
       dav1(idav,14)=qdiant
       dav1(idav,15)=0.
       if (exten(4).ne.0. .and. exten(5).ne.0.)
     *   dav1(idav,15)=exten(8)/sqrt(exten(4)*exten(5))
c    Emittance(norm)  x-xp (mm*mrad)
       dav1(idav,16)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
c    Emittance(non norm) x-xp (mm*mrad)
       dav1(idav,17)=surxth*10./pi
c --- dav1(idav,18): y-extension  (mm)
c --- dav1(idav,19): yp-extension (mrad)
c --- dav1(idav,20): correlation between y and yp
       dav1(idav,20)=0.
       if (exten(6).ne.0. .and. exten(7).ne.0.)
     *   dav1(idav,20)=exten(9)/sqrt(exten(6)*exten(7))
       dav1(idav,18)=qditay*10.
       dav1(idav,19)=qdianp
c    dav1(idav,21) : Emittance(norm) y-yp (pi*mm*mrad)
       dav1(idav,21)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
c    dav1(idav,22) : Emittance(non norm) y-yp (mm*mrad)
       dav1(idav,22)=suryph*10./pi
       dav1(idav,30)=float(ngood)
       dav1(idav,31)=cog(4)*10.
       dav1(idav,32)=cog(5)
       dav1(idav,33)=cog(6)*10.
       dav1(idav,34)=cog(7)
c next card to indicate no chase
       dav1(idav,26)=0.
c --- statistics with chase
       if (chasit) then
         dav2(idav,31)=fractx
         dav2(idav,32)=fracty
         dav2(idav,33)=fractl
c   longitudinal direction
         call chasel
         do i=1,ngood
           ichxyz(i)=ichas(i)
         enddo
         iarg=1
         call cdg(iarg)
         encog=cog(1)
         gcog=encog/xmat
         bcog=sqrt(1.-1./(gcog*gcog))
         enprt=encog-xmat
         CALL EXT2D(1)
c     qdisp : average dispersion dE (MeV)
c     sqmdv: longitudinal emittance (MeV*rad)
       qdisp=2.*sqrt(exten(1))
       qmdv=exten(1)*exten(3)-exten(2)*exten(2)
       sqmdv=4.*pi*sqrt(qmdv)
c     qdp : average extension dPHI (rad)
       qdp=2.*sqrt(exten(3))
c     cor12: coefficient of correlation in (dE, dPHI)
       cor12=exten(2)/sqrt(exten(1)*exten(3))
c sup       pent12=sqrt(exten(1)/exten(3))/cor12
c sup       pent21=sqrt(exten(1)/exten(3))*cor12
c sup       qdpde=qdp*180./pi
c ---  fictitious reference
c      dav2(idav,3): relativistic beta
c      dav2(idav,4): Kinetic energy (MeV)
c      dav2(idav,5): phase (in deg. w.r.t. k*pi)
         beref=vref/vl
         gref=1./sqrt(1.-beref*beref)
         dav2(idav,3)=beref
         dav2(idav,4)=xmat*(gref-1.)
         dav2(idav,5)=(int(tref*fh/pi+0.5)-tref*fh/pi)*180.
c ---  c.o.g of the bunch
c      dav2(idav,6): Kinetic energy (MeV)
c      dav2(idav,7): phase (in deg. w.r.t. k*pi)
       phnw=-(int(tcog*fh/pi+0.5)-tcog*fh/pi)*180.
       dav2(idav,6)=encog-xmat
       dav2(idav,7)=phnw
c ---- deviation between fictitious reference and c.o.g.
c     dav2(idav,8) : deviation in energy (MeV)
c     dav2(idav,9) : deviation of phase (deg)
       dav2(idav,8)=encog-xmat-dav2(idav,4)
       dav2(idav,9)=phnw-dav2(idav,5)
c ---- caracteristics of the bunch in longitudinal plane (dE, dPHI)
c       dav2(idav,10) :  dPHI extension (deg)
c       dav2(idav,11) :  dispersion dE (MeV)
c       dav2(idav,12) :  emittance (MeV*rad)
c       dav2(idav,23) :  coefficient of correlation in  (dE, dPHI)
       dav2(idav,10)=qdp*180./pi
       dav2(idav,11)=qdisp
       dav2(idav,12)=sqmdv/pi
       dav2(idav,23)=cor12
c sup       dav2(idav,39)=sqmdv*180./(pi*pi)
c sup       Ez(ns.keV)
c sup       dav2(idav,12)=sqmdv*(1.E09)/(pi*fh)
c ---  chase in x-xp
         call chasex
         do i=1,ngood
           ichxyz(i)=ichas(i)*ichxyz(i)
         enddo
         iarg=1
         call cdg(iarg)
         encog=cog(1)
         gcog=encog/xmat
         bcog=sqrt(1.-1/(gcog*gcog))
         dav2(idav,26)=cog(4)*10.
         dav2(idav,27)=cog(5)
         CALL EXT2D(1)
         trqty=exten(4)*exten(5)-exten(8)*exten(8)
         surxth=4.*pi*sqrt(trqty)
         qditax=2.*sqrt(exten(4))
         qdiant=2.*sqrt(exten(5))
ccc below for fourth line of EMIT output
c ---- caracteristics of the bunch in x-direction
c        dav2(idav,13): x-extension (mm)
c        dav2(idav,14): xp-extension (mrad)
c        dav2(idav,15): coefficient of correlation in plane (x, xp)
       dav2(idav,13)=qditax*10.
       dav2(idav,14)=qdiant
       dav2(idav,15)=0.
       if (exten(4).ne.0. .and. exten(5).ne.0.)
     *   dav1(idav,15)=exten(8)/sqrt(exten(4)*exten(5))
c    Emittance(norm) x-xp (pi*mm*mrad)
       dav2(idav,16)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
c    Emittance(non norm) in x-xp (mm*mrad)
       dav2(idav,17)=surxth*10./pi
ccc below for fourth line of EMIT output
c ---  chase in x-xp
         call chasey
         do i=1,ngood
           ichxyz(i)=ichas(i)*ichxyz(i)
         enddo
         iarg=1
         call cdg(iarg)
         encog=cog(1)
         gcog=encog/xmat
         bcog=sqrt(1.-1./(gcog*gcog))
         dav2(idav,28)=cog(6)*10.
         dav2(idav,29)=cog(7)
         CALL EXT2D(1)
         trqpz=exten(6)*exten(7)-exten(9)*exten(9)
         suryph=4.*pi*sqrt(trqpz)
         qditay=2.*sqrt(exten(6))
         qdianp=2.*sqrt(exten(7))
c     dav2(idav,18): y-extension  (mm)
c     dav2(idav,19): yp-extension (mrad)
c     dav2(idav,20): coefficient of correlation in the plane (y, yp)
       dav2(idav,20)=0.
       if (exten(6).ne.0. .and. exten(7).ne.0.)
     *   dav2(idav,20)=exten(9)/sqrt(exten(6)*exten(7))
       dav2(idav,18)=qditay*10.
       dav2(idav,19)=qdianp
c    dav2(idav,21) : Emittance(norm) in y-yp (mm*mrad)
       dav2(idav,21)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
c    dav1(idav,22) : Emittance(non norm) in y-yp (mm*mrad)
       dav2(idav,22)=suryph*10./pi
c next card to indicate chase
         dav1(idav,26)=1.
         do i=1,ngood
           ichas(i)=1
         enddo
       endif
c*et*2012*March*02 print energy, boro etc for each charge state in long
       if(ncstat.gt.1) then
         write(16,'(1x,a4)') 'EMIT'
         DO k=1,ncstat
           NDP(k)=0
           do j=2,7
             foo(k,j)=0.
           enddo
         ENDDO
         DO i=1,imax
           if(f(8,i).eq.1) then
             DO k=1,ncstat
               if(f(9,i).eq.cstat(k)) then
                 NDP(k)=NDP(k)+1
                 do j=2,7
                   foo(k,j)=foo(k,j)+f(j,i)
                 enddo
               endif
             ENDDO
           endif
         ENDDO
         DO k=1,ncstat
           if(NDP(k).ne.0) then
             do j=2,7
               foo(k,j)=foo(k,j)/float(NDP(k))
             enddo
           endif
         ENDDO
         write(16,*) ' Q     Particles  beta       Wcog(MeV)',
     *     '    Wcog(MeV/u)  Pcog(kG.cm)     TOF(deg)   TOF(sec)   ',
     *     '     X_avg(cm)     Xp_avg(mrad)  ',
     *     'Y_avg(cm)     Yp_avg(mrad)'
         DO k=1,ncstat
           if(NDP(k).ne.0) then
             gref=foo(k,7)/xmat
             bref=sqrt(1.-1./(gref*gref))
             xe=(gref-1.)*xmat
c   magnetic rigidity
             bor=3.3356*xmat*bref*gref/cstat(k)
             write(16,'(2x,f5.2,3x,I5,5x,F9.7,
     *       5(1x,E12.5),1x,4(F12.5,2x))')
     *       cstat(k),NDP(k),bref,xe,xe/atm,bor,foo(k,6)*180.*fh/pi,
     *       foo(k,6),foo(k,2),foo(k,3),foo(k,4),foo(k,5)
           else
             write(16,'(2x,f5.2,3x,I5)') cstat(k),NDP(k)
           endif
         ENDDO
         write(16,*)
       endif
       return
       end
       SUBROUTINE statis
c  ......................................................................
c       statitics of the 6-d ellipsoid (for print)
c   EXT2D:
c  exten(1) : Sum( dE(i)*dE(i) )  MeV*MeV
c  exten(2) : Sum( dE(i)*dPHase(i) ) MeV*rad
c  exten(3) : Sum( dPHase(i)*dPHase(i) ) rad*rad
c  exten(4) : Sum( x(i)*x(i) )   cm*cm
c  exten(5) : Sum( xp(i)*xp(i) )   mrad*mrad
c  exten(6) : Sum( y(i)*y(i) )   cm*cm
c  exten(7) : Sum( yp(i)*yp(i) )  mrad*mrad
c  exten(8) : Sum( x(i)*xp(i) )   cm*mrad
c  exten(9) : Sum( y(i)*yp(i) )   cm*mrad
c  ......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/STIS/surxth,suryph,enedep,ecogde,TESTCA
       common/etcom/cog(8),exten(17),fd(iptsz)
c EXT2D looks for average extensions squared and returns them in EXTEN
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       BREF=VREF/VL
       GREF=1./SQRT(1.-BREF*BREF)
       ENREF=XMAT*GREF
       CCGP=(tcog-tref)*fh*180./pi
       CCGD=ENCOG-ENREF
       DO I=1,ngood
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         fd(i)=bpai/bcog * gpai/gcog
       enddo
       CALL EXT2D(1)
       TESTCA=exten(1)*exten(2)*exten(3)
       IF(abs(TESTCA).gt.1.e-40) THEN
         QDISP=2.*SQRT(exten(1))
         QMD=exten(1)*exten(3)-exten(2)*exten(2)
         SURM=4.*PI*SQRT(QMD)*180./PI
         QDP=2.*SQRT(exten(3))
         COR12=exten(2)/sqrt(exten(1)*exten(3))
c  sup         PENT12=SQRT(exten(1)/exten(3))/COR12
c  sup         PENT21=SQRT(exten(1)/exten(3))*COR12
         QDPDE=QDP*180./PI
       ELSE
         QDISP=0.
         QMD=0.
         SURM=0.
         QDP=0.
         COR12=0.
         PENT12=0.
         PENT21=0.
         QDPDE=0.
       ENDIF
       TRQTX=EXTEN(4)*EXTEN(5)-EXTEN(8)*EXTEN(8)
       TRQPY=EXTEN(6)*EXTEN(7)-EXTEN(9)*EXTEN(9)
       QDITAX=2.*SQRT(EXTEN(4))
       QDIANT=2.*SQRT(exten(5))
       QDITAY=2.*SQRT(EXTEN(6))
       QDIANP=2.*SQRT(exten(7))
       SURXTH=4.*PI*SQRT(TRQTX)
       SURYPH=4.*PI*SQRT(TRQPY)
       SQMDV=4.*PI*SQRT(QMD)
       WRITE(16,52)IMAX,ngood
52     FORMAT(4X,' TOTAL NUMBER OF PARTICLES :',
     X    I5,'   NUMBER OF PARTICLES CONSIDERED :',I6,/)
       WRITE(16,1557)
1557   FORMAT(5X,' *** TRANSVERSE AND LONGITUDINAL STATISTICS')
       WRITE(16,1553) cog(4),cog(5)
1553   FORMAT(4X,' COG COORD   X : ',E12.5,' CM  XP :',E12.5,' MRD')
       WRITE(16,1556) cog(6),cog(7)
1556   FORMAT(4X,' COG COORD   Y : ',E12.5,' CM  YP :',E12.5,' MRD')
       WRITE(16,14)CCGP,CCGD
14     FORMAT(4X,' COG COORD dPHI: ',E12.5,' deg dW :',E12.5,' MeV')
       WRITE(16,1552) QDITAX,QDIANT,SURXTH/pi
1552   FORMAT(4X,' X :',E12.5,' CM  XP :',E12.5,' MRD  EMITTANCE :',
     * E15.8,' CM.MRD')
       WRITE(16,1555) QDITAY,QDIANP,SURYPH/pi
1555   FORMAT(4X,' Y :',E12.5,' CM  YP :',E12.5,' MRD  EMITTANCE :',
     * E15.8,' CM.MRD')
       WRITE(16,154) QDPDE,qdisp,SQMDV*180./(pi*pi)
154    FORMAT(4X,' dPHI :  ',F7.3,' deg dW :',3x,f7.3,3x,
     * 'MeV  EMITTANCE :',E15.8,' MeV.deg',/)
c sup       write(16,*) ' *** ellips (dE,dPHI) '
c sup       write(16,19) surm,qdp,qdisp,cor12
c sup19     format(' area: ',e12.5,' MeV*deg  half phase extent: ',e12.5,
c sup     *        ' deg half energy extend: ',e12.5,' MeV',/,
c sup     *        ' correlation coef: ',e12.5)
c sup     *
       RETURN
       END
       SUBROUTINE tiltbm(icg)
c  ......................................................................
c   tilt and shift of the beam with regard to the cog
c  ......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/TILT/TIPHA,TIX,TIY,SHIFW,SHIFP
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/DYN/TREF,VREF
       COMMON/STIS/suryth,surzph,enedep,ecogde,TESTCA
       COMMON/newtlt/twissa(3),itwiss
       COMMON/HISTO/CENTRE(6)
       COMMON/SHIF/DTIPH,SHIFT
       common/tof/tofini
       LOGICAL SHIFT
       common/tapes/in,ifile,meta
       common/etcom/cog(8),exten(17),fd(iptsz)
C      TIPHA : Shift on the phase axis(DEG)
C      TIX   : Shift in the x-direction (CM)
C      TIY   : Shift in the y-direction   (CM)
C      SHIFW : Change the energy position of the c.o.g.(MeV)
C      SHIFP : Change the phase position of the c.o.g.(deg)
C      DTIPH : Change the position of the phase (radian)
C      ICG  : = 0 VREF and TREF are the ones of sync. particle
C      ICG  :. NE. 0 VREF ET TREF are the ones of the c. of g.
       WRITE(16,1)TIPHA,TIX,TIY,SHIFW,SHIFP
1      FORMAT('  Shift the position of the bunch',/,
     X '    in the z-direction :',E12.5,' DEG',/,
     X '    in the x-direction :',E12.5,' CM ',/,
     X '    in the y-direction :',E12.5,' CM ',/,
     X '    shift the energy of the cog with  :',E12.5,' MeV',/,
     X '    shift the phase  of the cog with  :',E12.5,' deg',/)
       TIPHA=TIPHA*PI/180.
       SHIFP=SHIFP*PI/180.
       deg=fh*180./pi
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       tcog=cog(3)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
c --- this routine is call by the routine ENTRE with icg = 1
c       the reference particle is the cog of the bunch
        if(icg.ne.0) then
         tref=tcog
         vref=vl*bcog
         vcour=vref
         gcour=gcog
         bcour=bcog
         WRITE(16,*)'*** Before tilt and shift '
         WRITE(16,24)bcog,tref,tref*deg,encog-xmat
24       FORMAT(2X,'*** the reference particle is the cog:',/,
     X   ' REF AND COG: BETA :',E12.5,2X,'TOF :',E12.5,
     *   ' SEC OR: ',e12.5,' deg',
     x   2x,' ENERGY :',E12.5,' MeV',/)
       endif
       if(icg.eq.0) then
         vcour=vl*bcog
         bcour=bcog
         gcour=gcog
         bvref=vref/vl
         gamref=1./sqrt(1.-bvref*bvref)
         wvref=(gamref-1.)*xmat
         WRITE(16,*)'*** Before tilt and shift '
         WRITE(16,*)' the reference particle and cog are distinct'
         WRITE(16,16)bvref,tref,tref*deg,wvref,bcog,tcog,
     *      tcog*deg,encog-xmat
16     FORMAT(2X,'  REF: BETA ',E12.5,
     X ' T.O.F. ',E12.5,' SEC  OR ',E12.5,' DG',
     X ' ENERGY ',E12.5,' MeV',//,2x,
     X '  COG: BETA ',E12.5,' T.O.F. ',
     X  E12.5,' SEC  OR ',E12.5,' DG',
     X ' ENERGY ',E12.5,' MeV',/)
       endif
       DO I=1,ngood
         gapi=f(7,i)/xmat
         bepi=sqrt(1.-1./(gapi*gapi))
         fd(i)=bepi/bcour*gapi/gcour
       ENDDO
       CALL EXT2(1)
c        shift of the ellipse in the longitudinal plane with regard to the cog (TILT)
c        the slip of phase is : TIPHA
c     the ellipsoid generated by GEBEAM and INPUT is upright
c       qdispw=dispersion in energy
       QDISPW=2.*SQRT(exten(10))
       encrt=encog+qdispw
       gamcrt=encrt/xmat
       bcrt=sqrt(1.-1./(gamcrt*gamcrt))
       vcrt=bcrt*vl
       delv=vcrt-vcour
C    shift of the ellipse in the transverse planes with regard to the cog
C  DELTH ET DELPH: slips in xp and yp (mrd)
c     horizontal plane
       delxp=2.*sqrt(exten(5))
C    vertical plane
       delyp=2.*sqrt(exten(7))
       write(16,22) qdispw,delxp,delyp
22     format(' half size in energy ',e12.5,' MeV',/,
     *        ' half size in xp ',e12.5,' mrd',/,
     *        ' half size in yp ',e12.5,' mrd')
C    Shift of the ellipses in the phase spaces
       tcrt=0.
       do i=1,ngood
         gapi=f(7,i)/xmat
         bpai=sqrt(1.-1./(gapi*gapi))
         vpai=vl*bpai
         dv=(vpai-vcour)*tipha/(delv*fh)
         f(6,i)=f(6,i)-dv
         tcrt=tcrt+f(6,i)
         dlx=0.
         dly=0.
         if(delxp.ne.0.) dlx=f(3,i)*tix/delxp
         if(delyp.ne.0.) dly=f(5,i)*tiy/delyp
         f(2,i)=f(2,i)+dlx
         f(4,i)=f(4,i)+dly
       enddo
       tcrt=tcrt/float(ngood)
       DTIPH=0.
       if(shifw.ne.0. .or. shifp.ne.0.) then
         shtref=shifp/fh
         enshift=encog+shifw
         gshift=enshift/xmat
         bshift=sqrt(1.-1./(gshift*gshift))
         eshift=enshift-xmat
         vshift=vl*bshift
         tshift=tcrt+shtref
         deltav=vshift-vcour
         deltat=tshift-tcog
       endif
       if(shifw.eq.0. .and. shifp.eq.0.) then
         enshift=encog
         eshift=enshift-xmat
         gshift=gcour
         vshift=vcour
         bshift=vshift/vl
         tshift=tcrt
         shtref=0.
         deltav=0.
         deltat=0.
       endif
       do i=1,ngood
         f(7,i)=f(7,i)+shifw
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         f(6,i)=f(6,i)+shtref
         vapi=vl*bpai
       enddo
       if(itwiss.eq.1) tref=tref+tofini
       WRITE(16,*)' ***After tilt and shift '
       bvref=vref/vl
       gamref=1./sqrt(1.-bvref*bvref)
       wvref=(gamref-1.)*xmat
       write(16,16) bvref,tref,tref*deg,wvref,bshift,tshift,
     *              tshift*deg,eshift
       do i=1,ngood
         fo(7,i)=f(7,i)
         fo(2,i)=f(2,i)
         fo(3,i)=f(3,i)
         fo(4,i)=f(4,i)
         fo(5,i)=f(5,i)
         fo(6,i)=f(6,i)
       enddo
       if(itwiss.ne.1) then
         dum=0.
         write(11,*) ngood,dum,fh/(2000000.*pi)
         do i=1,ngood
           f(2,i)=f(2,i)+centre(2)
           f(3,i)=f(3,i)+centre(3)
           f(4,i)=f(4,i)+centre(4)
           f(5,i)=f(5,i)+centre(5)
           f(6,i)=f(6,i)+centre(6)
           f(7,i)=f(7,i)+centre(1)
           etphas=FH*(f(6,i)-tcog)
coption          etener=f(7,i)-encog
           etener=f(7,i)-xmat
           write(11,777) f(2,i),f(3,i)/1000.,f(4,i),f(5,i)/1000.,
     *                   etphas,etener
         enddo
777      format(6(F13.8,1x))
       endif
C      ENVELOPE
       CALL STAPL(davtot*10.)
       RETURN
       END
       SUBROUTINE tiltbm_bis(icg)
c  ......................................................................
c   Change the positions of the beam in the phase planes
c  ......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/TILT/TIPHA,TIX,TIY,SHIFW,SHIFP
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/objet/fo(9,iptsz),imaxo
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/DYN/TREF,VREF
       COMMON/STIS/suryth,surzph,enedep,ecogde,TESTCA
       COMMON/newtlt/twissa(3),itwiss
       COMMON/HISTO/CENTRE(6)
       COMMON/SHIF/DTIPH,SHIFT
comment      common/tof/tofini
       LOGICAL SHIFT
       common/tapes/in,ifile,meta
       common/etcom/cog(8),exten(17),fd(iptsz)
C      TIPHA : Shift with regard to the phase axis(DEG)
C      TIX   : Shift in the x-direction (CM)
C      TIY   : Shift in the y-direction   (CM)
C      SHIFW : Change the energy position of the c.o.g.(MeV)
C      SHIFP : Change the phase position of the c.o.g.(deg)
C      DTIPH : Change the position of the phase (radian)
C      ICG  : = 0 VREF and TREF are the ones of sync. particle
C      ICG  :. NE. 0 VREF ET TREF are the ones of the c. of g.
       WRITE(16,1)TIPHA,TIX,TIY,SHIFW,SHIFP
1      FORMAT(' shift the position around the c.o.g. of the bunch',/,
     X '    with regard to the phase axis :',E12.5,' DEG',/,
     X '    in the x-direction            :',E12.5,' CM ',/,
     X '    in the y-direction            :',E12.5,' CM ',/,
     X '   Change of energy position of the c.o.g.  :',E12.5,' MEV',/,
     X '   Change of phase  position of the c.o.g.  :',E12.5,' DEG',/)
       shtref=0.
       TIPHA=TIPHA*PI/180.
       SHIFP=SHIFP*PI/180.
       deg=fh*180./pi
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       tcog=cog(3)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       IF(ICG.NE.0 .AND. (SHIFW.EQ.0. .OR. SHIFP.EQ.0.)) THEN
C      the reference and the cog coincide
         tref=tcog
         vref=vl*bcog
         vcour=vref
         gcour=gcog
         bcour=bcog
         WRITE(16,*)' Before shift '
         WRITE(16,24)VREF,TREF,encog-xmat
24       FORMAT(2X,' Note : reference will coincide with the cog',/,
     X   2X,' velocity :',E12.5,' CM/SEC',2X,'tof :',E12.5,' SEC',/,
     x   2x,' ENERGY :',E12.5,' MeV',/)
       ELSE
C         the reference is distinct from the cog
         vcour=vl*bcog
         bcour=bcog
         gcour=gcog
         bvref=vref/vl
         gamref=1./sqrt(1.-bvref*bvref)
         wvref=(gamref-1.)*xmat
         WRITE(16,*)' Before shift '
         WRITE(16,16)VREF,TREF,tref*deg,wvref,vcour,tcog,
     *      tcog*deg,encog-xmat
       ENDIF
       DO I=1,ngood
         gapi=f(7,i)/xmat
         bepi=sqrt(1.-1./(gapi*gapi))
         fd(i)=bepi/bcour*gapi/gcour
       ENDDO
       CALL EXT2(1)
c        shift of the ellipse in the longitudinal plane with regard to the cog (TILT)
c        the slip of phase is : TIPHA
c     the elipsoid generated from MONTE end ENTRE is upright
c       qdispw=dispersion in energy
       QDISPW=2.*SQRT(exten(10))
       encrt=encog+qdispw
       gamcrt=encrt/xmat
       bcrt=sqrt(1.-1./(gamcrt*gamcrt))
       vcrt=bcrt*vl
       delv=vcrt-vcour
C    shift of the ellipse in the transverse planes with regard to the cog
C  DELTH ET DELPH: slips in xp and yp (MRD)
c     horizontal plane
       delxp=2.*sqrt(exten(5))
C    vertical plane
       delyp=2.*sqrt(exten(7))
       write(16,22) qdispw,delxp,delyp
22     format(' half size in energy ',e12.5,' MeV',/,
     *        ' half size in xp ',e12.5,' mrd',/,
     *        ' half size in yp ',e12.5,' mrd')
C    Shift of the ellipses in the phase spaces
       tcrt=0.
       do i=1,ngood
         gapi=f(7,i)/xmat
         bpai=sqrt(1.-1./(gapi*gapi))
         vpai=vl*bpai
         dv=(vpai-vcour)*tipha/(delv*fh)
         f(6,i)=f(6,i)-dv
         tcrt=tcrt+f(6,i)
         dlx=0.
         dly=0.
         if(delxp.ne.0.) dlx=f(3,i)*tix/delxp
         if(delyp.ne.0.) dly=f(5,i)*tiy/delyp
         f(2,i)=f(2,i)+dlx
         f(4,i)=f(4,i)+dly
       enddo
       tcrt=tcrt/float(ngood)
       DTIPH=0.
c     shift=.true.: The reference and COG are independent
c     shift=.false. The reference and COG are the same
       IF(SHIFW.NE.0. .OR. SHIFP.NE.0.) THEN
         SHIFT=.TRUE.
         SHTREF=SHIFP/FH
         enshift=encog+shifw
         gshift=enshift/xmat
         bshift=sqrt(1.-1./(gshift*gshift))
         eshift=enshift-xmat
         vshift=vl*bshift
         tshift=tcrt+shtref
       ELSE
         SHIFT=.FALSE.
         enshift=encog
         eshift=enshift-xmat
         gshift=gcour
         vshift=vcour
         tshift=tcrt
       ENDIF
       if(icg.eq.0) then
         deltav=vshift-vcour
         deltat=tshift-tcog
       else
         deltav=0.
         deltat=0.
         shifw=0.
         shtref=0.
       endif
       write(16,*) 'tofini,tcrt=',tofini,tcrt,tcrt*deg
       if(shift) then
         write(16,*) 'Reference and COG are independent',ICG
       else
         write(16,*) 'Reference and COG coincide',ICG
       endif
       do i=1,ngood
         f(7,i)=f(7,i)+shifw
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         f(6,i)=f(6,i)+shtref
         vapi=vl*bpai
       enddo
       if(itwiss.eq.1) tref=tref+tofini
       WRITE(16,*)' After shift '
       bvref=vref/vl
       gamref=1./sqrt(1.-bvref*bvref)
       wvref=(gamref-1.)*xmat
       write(16,16)vref,tref,tref*deg,wvref,vshift,tshift,
     *              tshift*deg,eshift
16     FORMAT(2X,'  REFERENCE: VELOCITY :',E12.5,
     X ' CM/SEC, T.O.F. :',E12.5,' SEC  OR ',E12.5,' DG',/,
     X 3x,' ENERGY : ',E12.5,' MeV',/,2x,
     X '  C.O.G.:    VELOCITY :',E12.5,' CM/SEC, T.O.F. :',
     X E12.5,' SEC OR ',E12.5,' DG',/,
     X 3x,' ENERGY : ',E12.5,' MeV',//)
       bpai=0.
       do i=1,ngood
         fo(7,i)=f(7,i)
         fo(2,i)=f(2,i)
         fo(3,i)=f(3,i)
         fo(4,i)=f(4,i)
         fo(5,i)=f(5,i)
         fo(6,i)=f(6,i)
       enddo
       if(itwiss.ne.1) then
         write(11,*) ngood
         do i=1,ngood
           f(2,i)=f(2,i)+centre(2)
           f(3,i)=f(3,i)+centre(3)
           f(4,i)=f(4,i)+centre(4)
           f(5,i)=f(5,i)+centre(5)
           f(6,i)=f(6,i)+centre(6)
           f(7,i)=f(7,i)+centre(1)
           etphas=FH*(f(6,i)-tcog)
coption          etener=f(7,i)-encog
           etener=f(7,i)-xmat
           write(11,777) f(2,i),f(3,i),f(4,i),f(5,i),etphas,etener
         enddo
777      format(6(F13.8,1x))
       endif
C      ENVELOPE
       CALL STAPL(davtot*10.)
       RETURN
       END
       SUBROUTINE accept
c ..........................................................
c   computes a physical beam acceptance at the entrance of the machine
c ..........................................................
       implicit real*8 (a-h,o-z)
       parameter(iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/TAPES/IN,IFILE,META
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/DYN/TREF,VREF
       COMMON/DYNI/VREFI,TREFI,FHINIT
       COMMON/QMOYEN/QMOY
       common/objet/fo(9,iptsz),imaxo
       COMMON/CONSTA/VL,PI,XMAT,RPEL,QST
       COMMON/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       DIMENSION BACK(9,iptsz)
       logical chasit
       common/mcs/imcs,ncstat,cstat(20)
c make statistics on initial beam using IEX at the point where
c the ACCEPT card has been placed
       write(16,*)
       write(16,*) 'Physical acceptance has:'
       write(16,*) imaxo,' particles at origin'
       iprint=1
       savfh=fh
       satref=tref
       savref=vref
c recall initial frequency
       fh=fhinit
       tref=trefi
       vref=vrefi
c save current beam
       do i=1,9
         do j=1,ngood
           back(i,j)=f(i,j)
         enddo
       enddo
c original particle number is in f(1,j)
c F (I,8): =1 the particle is good,   =0 the particle is lost
C recuperate original particle numbers and their coordinates
C next loop for graphics and file of good particles
       do j=1,imax
         tprt=fh*fo(6,j)
         eprt=fo(7,j)-xmat
         f2=fo(2,j)
         f3=fo(3,j)/1000.
         f4=fo(4,j)
         f5=fo(5,j)/1000.
         nold=int(f(1,j))
       enddo
       open(23,file='input_kept.dst',status='unknown')
       dummy=0.
       write(23,*) ngood,dummy, fh/(2000000.*pi)
       do j=1,ngood
         nold=int(f(1,j))
         do jj=1,9
           f(jj,j)=fo(jj,nold)
         enddo
         tprt=fh*f(6,j)
         eprt=f(7,j)-xmat
         f2=f(2,j)
         f3=f(3,j)
         f4=f(4,j)
         f5=f(5,j)
         if(ncstat.gt.1) then
           write(23,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,j)
         else
           write(23,100)f2,f3/1000.,f4,f5/1000.,tprt,eprt
         endif
       enddo
       close(23)
c now make graphics
       write(16,*) 'Starting good particles graphics for ACCEPT card'
       igrprm=0
       call ytzp
c print out beam statistics
       call emiprt(0)
C next loop for graphics and file of lost particles
       open(23,file='input_lost.dst',status='unknown')
       dummy=0.
       write(23,*) imax-ngood,dummy, fh/(2000000.*pi)
       j=1
       ngd=ngood
       do k=ngood+1,imax
         nold=int(f(1,k))
         do jj=1,9
           f(jj,j)=fo(jj,nold)
         enddo
         tprt=fh*f(6,j)
         eprt=f(7,j)-xmat
         f2=f(2,j)
         f3=f(3,j)
         f4=f(4,j)
         f5=f(5,j)
         if(ncstat.gt.1) then
           write(23,101)f2,f3/1000.,f4,f5/1000.,tprt,eprt,f(9,j)
         else
           write(23,100)f2,f3/1000.,f4,f5/1000.,tprt,eprt
         endif
         j=j+1
       enddo
       ngood=imax-ngood
100    format(6(F15.8,1x))
101    format(7(F15.8,1x))
       close(23)
c now make graphics
       write(16,*) 'Starting lost particles graphics for ACCEPT card'
       igrprm=0
       call ytzp
c print out beam statistics
       call emiprt(0)
c recall original coordinates
       ngood=ngd
       do i=1,9
         do j=1,ngood
           f(i,j)=back(i,j)
         enddo
       enddo
       tref=satref
       vref=savref
       fh=savfh
       return
       end
       SUBROUTINE ytzp
c ..........................................................
C      STORAGE OF PARTICLE COORDINATES FOR PLOTS
c ..........................................................
       implicit real*8 (a-h,o-z)
       parameter(iptsz=100002,maxcell=3000,maxcell1=3000)
       CHARACTER*80 TEXT,PATITL
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/STIS/suryth,surzph,enedep,ecogde,TESTCA
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
       common/tapes/in,ifile,meta
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/GRPARM/GLIM(4,2),GLIM1(4,2),GLIM2(4,2),PATITL,
     *        ngraphs(100),idwdp,igrprm,ngrafs
       common/mcs/imcs,ncstat,cstat(20)
       common/zones/frms(6),nzone
       LOGICAL chasit
       DIMENSION oldcog(7),slim(4,2)
c       real*4 xx(iptsz),yy(iptsz),cs(iptsz),axx,ayy
       dimension xx(iptsz),yy(iptsz),cs(iptsz)
       satref=tref
       savref=vref
       init=0
       call area(init)
       if (igrprm.eq.0) then
C        READ GRAPH TITLE
         READ(IN,20) TEXT
20       FORMAT(A)
c idwdp=0 cog=ref in ZZ' plot (for instance for Alvarez structure)
c idwdp=1 cog<>ref in ZZ' plot (for instance for IH structrure)
         READ(IN,*) idwdp,rmsmtp
C       READ GRAPH LIMITS INTO GLIM(J,K), J=GRAPH NUMBER
C       K=1 HOR. LIMIT , K=2 VERT. LIMIT
         READ(IN,*) ((GLIM(J,K), K=1,2), J=1,4)
       endif
       if (igrprm.eq.1) then
         text=patitl
       endif
       if (igrprm.eq.2) then
         text=patitl
c  save limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             slim(j,k)=glim(j,k)
             glim(j,k)=glim1(j,k)
           enddo
         enddo
       endif
       if (igrprm.eq.3) then
         text=patitl
c  save limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             slim(j,k)=glim(j,k)
             glim(j,k)=glim2(j,k)
           enddo
         enddo
       endif
       do i=1,ngood
         ichas(i)=1
       enddo
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       oldcog(1)=cog(4)
       oldcog(2)=cog(5)
       oldcog(3)=cog(6)
       oldcog(4)=cog(7)
       oldcog(5)=cog(3)
       oldcog(6)=cog(8)
       oldcog(7)=cog(1)
       enprt=encog-xmat
       teps=1.e-08
       rmssz=sqrt(rmsmtp)
       iarg=1
       call ext2(iarg)
       QDISP=rmssz*SQRT(exten(1))
       QMD=exten(1)*exten(3)-exten(2)*exten(2)
       QMDW=exten(10)*exten(3)-exten(11)*exten(11)
       SURM=rmsmtp*180.*SQRT(QMD)
       QDP=rmssz*SQRT(exten(3))
       COR12=exten(2)/sqrt(exten(1)*exten(3))
       PENT12=SQRT(exten(1)/exten(3))/COR12
       PENT21=SQRT(exten(1)/exten(3))*COR12
       QDPDE=QDP*180./pi
       TRQTY=exten(4)*exten(5)-exten(8)*exten(8)
       TRQPZ=exten(6)*exten(7)-exten(9)*exten(9)
       SURYTH=rmsmtp*PI*SQRT(TRQTY)
       SURZPH=rmsmtp*PI*SQRT(TRQPZ)
       QDITAY=rmssz*SQRT(exten(4))
       QDIANT=rmssz*SQRT(exten(5))
       QDITA=rmssz*SQRT(exten(6))
       QDIANP=rmssz*SQRT(exten(7))
       write(16,*) ' *** PLOT Ellips for ',rmsmtp,' RMS'
       write(16,'(a)') text
       WRITE(16,1557) IMAX,ngood
1557   FORMAT(1x,' *** GRAPH, TOTAL NUMBER OF PARTICLES : ',I6,
     1 ' PARTICLES KEPT : ',I6,//,' ***  HORIZONTAL phase plane ',/)
       WRITE(16,1553) cog(4),cog(5)
1553   FORMAT(4X,' C.O.G. :',
     1 5X,' X : ',E12.5,' CM  XP :',E12.5,' MRD',/)
       if(rmsmtp.gt.teps) then
         WRITE(16,1552) QDITAY,QDIANT,SURYTH
1552     FORMAT(4X,' 1/2 EXTENSION X : ',E12.5,' CM',/,4X,
     1   ' 1/2 EXTENSION XP : ',E12.5,' MRD',4X,
     2   ' SURFACE : ',E15.8,' CM.MRD',/)
       endif
       WRITE(16,1554)
1554   FORMAT(' ***  VERTICAL phase plane ',/)
       WRITE(16,1556) cog(6),cog(7)
1556   FORMAT(4X,' C.O.G :',
     1   5X,' Y : ',E12.5,' CM  YP :',E12.5,' MRD',/)
       if(rmsmtp.gt.teps) then
         WRITE(16,1555) QDITA,QDIANP,SURZPH
1555     FORMAT(4X,' 1/2 EXTENSION Y : ',E12.5,' CM',/,4X,
     1   ' 1/2 EXTENSION YP : ',E12.5,' MRD',4X,
     2   ' SURFACE : ',E15.8,' CM.MRD',/)
       endif
C Store header and particle coordinates in binary file for
C graphics post-processor
C
C igrtyp is type of graph
C        igrtyp=1  for ytzp emittance plots
C        igrtyp=6  for ytzp emittance plots for multi-charge state beam
C        igrtyp=11 for ytzp emittance plots with ZONES card
C
       igrtyp=1
       if(nzone.ne.0)igrtyp=11
       if(imcs.eq.1) igrtyp=6
       WRITE(66,*) igrtyp
       if (igrtyp.eq.6)then
         write(66,*) ncstat
         write(66,*) (cstat(j),j=1,ncstat)
       endif
       if (igrtyp.eq.11)then
         write(66,*) nzone
         write(66,*) (frms(j),j=2,nzone),' 0.'
       endif
       WRITE(66,*) text
       xx(1)=GLIM(1,1)
       yy(1)=GLIM(1,2)
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
C     first store ellips coordinates XX'
       step=qdiant/50.
       TTA=exten(5)
       do 345 i=1,101
         xii=-qdiant+step*float(i-1)
         yy(i)=xii+cog(5)
         yy(202-i)=yy(i)
         TTB=exten(8)*XII
         TTC=exten(4)*XII**2-TRQTY*rmsmtp
         TTCB=TTB**2-TTC*TTA
         IF(TTCB.LE.0.) TTCB=0.
         IF(TTA.eq.0.) then
           YI=0.
           YII=0.
         ELSE
           QUOT=TTCB/TTA**2
           YI=TTB/TTA-SQRT(quot)
           YII=TTB/TTA+SQRT(quot)
         ENDIF
         xx(i)=yi+cog(4)
         xx(202-i)=yii+cog(4)
345    continue
       do i=1,201
         WRITE(66,*) xx(i),yy(i)
       enddo
       do i=1,ngood
         xx(i)=F(2,I)
         yy(i)=F(3,I)
       enddo
       if (imcs.eq.1) then
         do i=1,ngood
           cs(i)=f(9,i)
         enddo
       endif
C write particle coordinates to graphics file
       WRITE(66,*) ngood
       if (imcs.eq.0) then
         if(nzone.eq.0) then
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i)
           enddo
         else
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i),f(10,i)
           enddo
         endif
       else
         do i=1,ngood
          WRITE(66,*) xx(i),yy(i),cs(i)
         enddo
       endif
C YY' next
       xx(1)=GLIM(2,1)
       yy(1)=GLIM(2,2)
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
C      first write ellips coordinates
       step=qdianp/50.
       TTA=exten(7)
       do 346 i=1,101
         xii=-qdianp+step*float(i-1)
         yy(i)=xii+cog(7)
         yy(202-i)=yy(i)
         TTB=exten(9)*XII
c        TTC=exten(6)*XII**2-TRQPZ*4.
         TTC=exten(6)*XII**2-TRQPZ*rmsmtp
         TTCB=TTB**2-TTC*TTA
         IF(TTCB.LE.0.) TTCB=0.
         IF(TTA.eq.0.) then
           YI=0.
           YII=0.
         ELSE
           QUOT=TTCB/TTA**2
           YI=TTB/TTA-SQRT(quot)
           YII=TTB/TTA+SQRT(quot)
         ENDIF
         xx(i)=yi+cog(6)
         xx(202-i)=yii+cog(6)
346    continue
       do i=1,201
        WRITE(66,*) xx(i),yy(i)
       enddo
       do i=1,ngood
         xx(i)=F(4,I)
         yy(i)=F(5,I)
       enddo
C write particle coordinates to graphics file
       WRITE(66,*) ngood
       if (imcs.eq.0) then
         if(nzone.eq.0) then
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i)
           enddo
         else
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i),f(10,i)
           enddo
         endif
       else
         do i=1,ngood
          WRITE(66,*) xx(i),yy(i),cs(i)
         enddo
       endif
C
c1400   CONTINUE
C     TRACE GRAPHE ZZ'
       xx(1)=GLIM(3,1)
       yy(1)=GLIM(3,2)
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
c
       xx(1)=GLIM(4,1)
       yy(1)=GLIM(4,2)
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
C       TRACE DE L ELLIPSE  DE CONCENTRATION
       bref=vref/vl
       gref=1./sqrt(1.-bref*bref)
       WREF=XMAT*(gref-1.)
       gcog=1./sqrt(1.-bcog*bcog)
       WCOG=XMAT*(gcog-1.)
       WRITE(16,22) WREF,TREF,WCOG,tcog
22     FORMAT(' *** LONGITUDINAL phase plane ',/,
     1 6X,' REFERENCE : ',
     2 ' ENERGY: ',E15.8,' (MeV), TOF: ',E15.8,' (SEC)',
     3 /,   6X,' COG       : ',
     4 ' ENERGY: ',E15.8,' (MeV), TOF: ',E15.8,' (SEC)',/)
       WRITE(16,167) QMD,SURM,QDP,QDPDE,QDISP,COR12,PENT12,PENT21
167    FORMAT(3X,' ***',
     2 ' 2nd ORDER MOMENTS :',E12.5,' (RD DP/P)**2',
     3 '    SURFACE : ',E12.5,' (DEG DP/P)',/,
     4 6X,' 1/2 EXTENSION PHASE : ',E12.5,' RD ',
     5 ' OR ',E12.5,' DEG',/,6X,
     6 ' 1/2 EXTENSION DISPERSION : ',E15.8,' IN DP/P ',/,6X,
     7 ' CORRELATION COEF : ',E15.8,/,6X,
     8 ' DISPERSION SLOPE: ',E15.8,' (DP/P)/RD ',/,6X,
     9 ' PHASE SLOPE : ',E15.8,' (DP/P)/RD ')
C      GRAPHE DE L ELLIPSE
       step=qdpde/50.
       TTA=exten(3)*180.*180./(pi*pi)
       do 347 i=1,101
         xii=-qdpde+step*float(i-1)
         xx(i)=xii
         xx(202-i)=xx(i)
         TTB=exten(11)*XII*180./pi
         TTC=exten(10)*XII**2-QMDW*rmsmtp*180.*180./(pi*pi)
         TTCB=TTB**2-TTC*TTA
         IF(TTCB.LE.0.) TTCB=0.
         IF(TTA.eq.0.) then
           YI=0.
           YII=0.
         ELSE
           QUOT=TTCB/TTA**2
           YI=TTB/TTA-SQRT(quot)
           YII=TTB/TTA+SQRT(quot)
         ENDIF
         yy(i)=yi
         yy(202-i)=yii
347    continue
       enihrf=0.
       if(idwdp.eq.0) then
c Alvarez
         do i=1,201
         WRITE(66,*) xx(i),yy(i)
         enddo
       else
c IH
         bref=vref/vl
         gref=1./sqrt(1.-bref*bref)
         ENIHRF=XMAT*(gref-1.)
         PHIHRF=FH*(tcog-tref)*180./pi
         do i=1,201
           axx=phihrf
           ayy=encog-enihrf-xmat
       WRITE(66,*) xx(i)+axx,yy(i)+ayy
         enddo
       endif
C  write particle coordinates to graphics file
       do i=1,ngood
         if (idwdp.eq.0) then
c Alvarez like
           xx(i)=fh*(f(6,i)-tcog)*180./pi
           yy(i)=f(7,i)-encog
           if (imcs.eq.1) cs(i)=f(9,i)
         endif
         IF(idwdp.eq.1) THEN
c IH like
           xx(i)=FH*(f(6,i)-tref)*180./pi
           yy(i)=f(7,i)-enihrf-xmat
           if (imcs.eq.1) cs(i)=f(9,i)
         endif
       enddo
C write particle coordinates to graphics file
       WRITE(66,*) ngood
       if (imcs.eq.0) then
         if(nzone.eq.0) then
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i)
           enddo
         else
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i),f(10,i)
           enddo
         endif
       else
         do i=1,ngood
          WRITE(66,*) xx(i),yy(i),cs(i)
         enddo
       endif
C plot emittance  with chase (if applicable)
c      if (chasit) call ytzprc(glim,oldcog,idwdp)
       if (igrprm.eq.2.or.igrprm.eq.3) then
c  restore limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             glim(j,k)=slim(j,k)
           enddo
         enddo
       endif
       tref=satref
       vref=savref
       RETURN
       END
       SUBROUTINE grcomp(text,iskale)
c ..........................................................
c      routine for plots
c ..........................................................
       implicit real*8 (a-h,o-z)
       parameter(iptsz=100002,maxcell=3000,maxcell1=3000)
       CHARACTER*80 TEXT,PATITL
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/STIS/suryth,surzph,enedep,ecogde,TESTCA
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       common/tapes/in,ifile,meta
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/GRPARM/GLIM(4,2),GLIM1(4,2),GLIM2(4,2),PATITL,
     *       ngraphs(100),idwdp,igrprm,ngrafs
       common/mcs/imcs,ncstat,cstat(20)
       common/zones/frms(6),nzone
       common/hist/xpos(200),xn(200),ypos(200),yn(200),
     *             zpos(200),zn(200),ixt,iyt,izt
       common/hist1/xps(200),xpn(200),yps(200),ypn(200),
     *             zps(200),zpn(200),ixpt,iypt,izpt
       LOGICAL chasit
c       real*4 xx(iptsz),yy(iptsz),cs(iptsz)
       dimension xx(iptsz),yy(iptsz),cs(iptsz)
       DIMENSION SLIM(4,2)
       tcog=0.
       do i=1,ngood
         tcog=f(6,i)+tcog
       enddo
       tcog=tcog/float(ngood)
       if (igrprm.eq.1) then
         text=patitl
       endif
       if (igrprm.eq.2) then
         text=patitl
c  save limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             slim(j,k)=glim(j,k)
             glim(j,k)=glim1(j,k)
           enddo
         enddo
       endif
       if (igrprm.eq.3) then
         text=patitl
c  save limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             slim(j,k)=glim(j,k)
             glim(j,k)=glim2(j,k)
           enddo
         enddo
       endif
       write(16,*)
       write(16,*) 'LIMITS',((GLIM(J,K), K=1,2), J=1,4)
c X-Z start
       XYLIZ=GLIM(4,1)*VREF*PI/(180.*FH)
       ZLIX=XYLIZ
       XLIX=GLIM(3,1)
       DSTRLY=GLIM(4,2)
C
C Store header and particle coordinates in binary file for
C graphics post-processor
C
C igrtyp is type of graph
c        igrtyp=2  for grcomp extra plots
c        igrtyp=7  for grcomp extra plots with multi charge state beam
C        igrtyp=12 for grcomp extra plots with ZONES card
       igrtyp=2
       if(nzone.ne.0) igrtyp=12
       if (imcs.eq.1)igrtyp=7
c igrtyp=2,7,12   --> normal scale in envelopes
c igrtyp=17,22,27 --> log    scale in envelopes
       if(iskale.eq.1) then
         WRITE(66,*) igrtyp+15
         write(66,*) dstrly
       else
         WRITE(66,*) igrtyp
       endif
       if (igrtyp.eq.7)then
         write(66,*) ncstat
         write(66,*) (cstat(j),j=1,ncstat)
       endif
       if (igrtyp.eq.12)then
         write(66,*) nzone
         write(66,*) (frms(j),j=2,nzone),' 0.'
       endif
       WRITE(66,*) text
cnew
       xx(1)=zlix
       yy(1)=xlix
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
c120    CONTINUE
       if (imcs.eq.1) then
         do i=1,ngood
           cs(i)=f(9,i)
         enddo
       endif
       if(idwdp.eq.0) then
c Alvarez type (plot w.r.t. reference)
         DO I=1,ngood
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           xx(i)=(tcog-f(6,i))*vl*bpai
           yy(i)=F(2,I)
         ENDDO
       else
c IH type (plot w.r.t. centre of gravity)
         DO I=1,ngood
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           xx(i)=(tref-f(6,i))*vl*bpai
           yy(i)=F(2,I)
         ENDDO
       endif
C write particle coordinates to graphics file
       WRITE(66,*) ngood
       if (imcs.eq.0) then
         if(nzone.eq.0) then
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i)
           enddo
         else
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i),f(10,i)
           enddo
         endif
       else
         do i=1,ngood
          WRITE(66,*) xx(i),yy(i),cs(i)
         enddo
       endif
C Y-Z start
       ZLIY=XYLIZ
       ZLIY=GLIM(3,2)
       xx(1)=zlix
       yy(1)=zliy
       WRITE(66,*) -xx(1),xx(1),-yy(1),yy(1)
c1200   CONTINUE
       if(idwdp.eq.0) then
c Alvarez type (plot w.r.t. reference)
         do i=1,ngood
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           xx(i)=(tcog-f(6,i))*vl*bpai
           yy(i)=F(4,I)
         enddo
       else
c IH type (plot w.r.t. centre of gravity)
         do i=1,ngood
           gpai=f(7,i)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           xx(i)=(tref-f(6,i))*vl*bpai
           yy(i)=F(4,I)
         enddo
       endif
C write particle coordinates to graphics file
       WRITE(66,*) ngood
       if (imcs.eq.0) then
         if(nzone.eq.0) then
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i)
           enddo
         else
           do i=1,ngood
             WRITE(66,*) xx(i),yy(i),f(10,i)
           enddo
         endif
       else
         do i=1,ngood
          WRITE(66,*) xx(i),yy(i),cs(i)
         enddo
       endif
c1400   CONTINUE
c beam profile plots for X,Y & Z
       call histgrm
       write(66,*) ixt
       do i=1,ixt
         write(66,*) xpos(i),xn(i)
       enddo
       write(66,*) iyt
       do i=1,iyt
         write(66,*) ypos(i),yn(i)
       enddo
       write(66,*) izt
       do i=1,izt
         write(66,*) zpos(i),zn(i)
       enddo
c beam profile plots for Xp,Yp & Zp
       write(66,*) ixpt
       do i=1,ixpt
         write(66,*) xps(i),xpn(i)
       enddo
       write(66,*) iypt
       do i=1,iypt
         write(66,*) yps(i),ypn(i)
       enddo
       write(66,*) izpt
       do i=1,izpt
         write(66,*) zps(i),zpn(i)
       enddo
       if (igrprm.eq.2.or.igrprm.eq.3) then
c  restore limits in GLIM(j,k)
         do j=1,4
           do k=1,2
             glim(j,k)=slim(j,k)
           enddo
         enddo
       endif
       RETURN
       END
       SUBROUTINE restay
c   ............................................................................
c   motion of particles in a cavity
c    the field can read on the disk on the form: (z,E(z)
c    or it read in the command list on the form of a Fourier series expansion
c   ............................................................................
       implicit real*8 (a-h,o-z)
C      ****************************************************
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFS/DYNT(MAXCELL),DYNTP(MAXCELL),DYNTPP(MAXCELL),
     *   DYNE0(MAXCELL),DYNPH(MAXCELL),DYNLG(MAXCELL),FHPAR,NC
       COMMON/POSI/IST
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/AZMTCH/DLG,XMCPH,XMCE
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
C      TRANSIT TIME COEFFICIENTS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/TTFCB/T3K,T4K,S3K,S4K
       COMMON/BEDYCP/PHSLIP,EQVL,ASDL,PEQVL,PAVPH,XKP1,XKP2,
     *                AA,BB,CC,DD,EE,PCREST,SQCTTF
       COMMON/JACOB/GAKS,GAPS
       common/iter1/DXDKI,DPHII,PHI,DKMSKE,DKMSPHI,RETPH,XKMI,XKM,
     *              DXK00,TKE,T1KE,SKE,S1KE,PHIWC,XK1I,XK1II,XK2II
       common/faisc/f(10,iptsz),imax,ngood
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/rfield/ifield
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/CDEK/DWP(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/compt1/ndtl,ncavmc,ncavnm
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TAPES/IN,IFILE,META
       COMMON/RANEC1/DUMMY(6)
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/speda/dave,idave
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/DCSPA/IESP
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/appel/irstay,ilost,iavp,ispcel
       common/posc/xpsc
       common/pstpla/tstp
       common/femt/iemgrw,iemqesg
       common/mode/eflvl,rflvl
       common/aerp/vphase,vfield,ierpf
       common/tofev/ttvols
c ****    reference
c ****     DWRFS(MeV): energy gain
c ****     SPHRFS(rad):phase jump
c ****     PHRFS(rad):phase
c ****       common/parmrf/DWRFS,SPHRFS,PHRFS,ngdrf
       logical iesp,ichaes,irstay,iavp,ispcel,ifield,iemgrw
       LOGICAL SHIFT,CHASIT,ITVOL,IMAMIN,DAVE
       CHARACTER*1 cr
C************************************************************
C    XESLN : NEGATIVE LENGHT OF THE DRIFT FOLLOWING THE GAP
C    IF XESLN N.E.0 THEN THE CHARGE SPACE EFFECT IMPLIES THE
C    LENGTH (YLG-XESLN)
       NRRES=NRRES+1
       ncavmc=ncavmc+1
c allow for print out on terminal of gap# on one and the same line
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,advance='no') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       WRITE(16,*)' CAVITY N :',NRRES
       ilost=0
       aqst=abs(qst)
       qmoy=aqst
c --- the frequency fh may be changed with delfh
       oldfh=fh
C        IDUM : dummy variable (indicate in input file the number of the cavity)
         READ (IN,*) IDUM
c
c --- XESLN (cm): difference between the length of the field YLG and the physical length of the cavity
c      (The space charge is acting on the length (YLG-XESLN) )
c      dphase (deg): phase offset from the phase crest (giving the maximum of energy gain)
c      FFIELD : in percent;
c                    (electric field)=(initial electric field)*(1.+FFIELD/100)
c         isec: flag
c           isec = 0 : The crest phase (or the valley) is adjusted at the entrance of the accelerating element
c           isec = 1 : The crest phase (or the valley)is adjusted at the middle of the accelerating element
c
c           idum is for compatability with CAVNUM
         read (in,*) xesln,dphase,ffield,isec,idum
         ffield=1.+ffield/100.
       if(ifield) then
c --- The field is read on the disk in file 'field.txt' in the form:
c          z   (cm)   E(z) MV/cm
c     fhc: frequency of the cavity (Hertz) (read in the file 'field.txt' )
         fh=fhc*2.*pi
         ncel=ncell
         atte=att
         ye0=atte
c    flength : length of the field (cm)
         ylg=flength
         SCDIST=YLG-ABS(XESLN)
       else
c --- The field is read in the input list in the form of a Fourier series
         oldfh=fh
c     atte: factor acting on the amplitude of the field (read in the input list)
         ye0=atte
         SCDIST=YLG-ABS(XESLN)
       endif
       dphete=dphase
       if(itvol.and.imamin) then
c ---  adjustment of the phase offset w.r.t. the t.o.f.  (deg)
        ottvol=fh*ttvols*180./pi
        attvol=ottvol
        xkpi=ottvol/360.
        ixkpi=int(xkpi)
        xkpi=(xkpi-float(ixkpi))*360.
        dphase=dphase-xkpi
       endif
c --- iesp, irstay and ispcel: logical flags for space charge computations
       iesp=.false.
       irstay=.true.
       ispcel=.true.
c ---  dwp(*): array reserved to space charge computations
       do i=1,iptsz
         dwp(i)=0.
       enddo
       WRITE(16,150)FH/(2.*pi),YLG,ATTE,ffield,NCEL
150    FORMAT(4X,'FREQUENCY :',E12.5,' Hertz',/,4x,
     x        'FIELD LENGTH :',e12.5,' cm',/,4x,
     x        'FIELD FACTOR (UNITS CONVERSION) :',e12.5,/,4x,
     x        'FIELD FACTOR (ATTENUATION)      :',f12.6,/,4x,
     X        'FIELD DIVIDED IN: ',I4,' SECTIONS ')
       if(.not.imamin) write(16,*) '    PHASE OFFSET: ',dphete,' DEG'
       if(imamin) write(16,1501)dphete,DPHASE,xkpi
1501   format(4x,
     x        'PHASE OFFSET (before adjustment): ',e12.5,' deg',/,4x,
     x        'PHASE OFFSET (after adjustment): ',e12.5,' deg',/,4x,
     *        'ADJUSTMENT ON THE PHASE OFFSET: ',e12.5,' deg')
       BEREF=VREF/VL
       fh0=fh/vl
C --- prediction of transit time factors TK and SK  based on the velocity at the entrance
       TK=TTA0(BEREF)/2. * FFIELD
       SK=TSB0(BEREF)/2. * FFIELD
C --- prediction of PCREST (phase of RF giving the maximum of energy gain in the cavity)
       PCREST=ATAN(-SK/TK)
       DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       IF(DDWC.LT.0.) PCREST=PCREST+PI
c --- ttvol: time of flight at entrance (sec)
       ttvol=0.
       if(itvol)ttvol=ttvols*fh
c  start file 'short.data'
c --- dav1(idav,3)=0: the particle reference and the cog coincide at the input
c --- dav1(idav,3)=1: the particle reference and the cog are independent
       dav1(idav,3)=0.
       idav=idav+1
       iitem(idav)=1
       dav1(idav,1)=ylg*10.
       dav1(idav,2)=ye0*100.
       tstp=(davtot+ylg*xpsc)*10.
       davtot=davtot+ylg
       dav1(idav,24)=davtot*10.
       dav1(idav,40)=fh
       IF(IPRF.EQ.1) CALL STAPL(dav1(idav,24))
       iarg=1
       call cdg(iarg)
       enold=cog(1)
       encog=enold
       gcog=enold/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       IF(SHIFT) THEN
c --- the reference particle and the cog are independent
         BEREF=VREF/VL
         GAMREF=1./SQRT(1.-(BEREF*BEREF))
         ENREF=XMAT*GAMREF
         TREFDG=TREF*FH*180./PI
         dav1(idav,3)=1.
       ELSE
c --- the reference particle and the cog are coinciding
         beref=bcog
         vref=bcog*vl
         tref=tcog
         gamref=gcog
         enref=cog(1)
         dav1(idav,3)=0.
       ENDIF
c ---   the reference particle is put in the array f(10,iptsz) at the position ngdrf = ngood + 1
c ****         ngdrf=ngood+1
c ****         BEREF=VREF/VL
c ****         GAMREF=1./SQRT(1.-(BEREF*BEREF))
c ****         ENREF=XMAT*GAMREF
c ****         f(1,ngdrf)=ngdrf
c ****         f(2,ngdrf)=0.
c ****         f(3,ngdrf)=0.
c ****         f(4,ngdrf)=0.
c ****         f(5,ngdrf)=0.
c ****         f(6,ngdrf)=tref
c ****         f(7,ngdrf)=enref
c ****         f(8,ngdrf)=1.
c ****         f(9,ngdrf)=qst
c ****         f(10,ngdrf)=0.
         if(dav1(idav,3).eq.1.) write(16,*)
     *   ' ****reference and cog evolve independently'
         if(dav1(idav,3).eq.0.) write(16,*)
     *   ' **** the reference is the cog '
       WRITE(16,178)
178    FORMAT(/,' Dynamics at the input',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       WRITE(16,1788) bcog,gcog,encog-xmat,tcog*fh*180./pi,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       ENRPRIN=ENREF-XMAT
       WRITE(16,165) beref,gamref,ENRPRIN,tref*fh*180./pi,tref
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       tofprt=tref
       iprint=0
       call statis
       XK1=FH/VREF
c --- prediction of transit time factors based on an average value of velocity
       DDW=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       ENREFS=ENREF+DDW
       GAMS= ENREFS/XMAT
       BETS=SQRT(1.-1./(GAMS*GAMS))
       BEMY=(GAMS+GAMREF)/(GAMS*BETS+GAMREF*BEREF)
       BEMY=1./BEMY
       TK0=TTA0(BEMY)/2.  * FFIELD
       TPK0=TTA1(BEMY)/2. * FFIELD
       TPPK0=TTA2(BEMY)/2.* FFIELD
       TP3K0=TTA3(BEMY)/2.* FFIELD
       TP4K0=TTA4(BEMY)/2.* FFIELD
       SK0=TSB0(BEMY)/2.  * FFIELD
       SPK0=TSB1(BEMY)/2. * FFIELD
       SPPK0=TSB2(BEMY)/2.* FFIELD
       SP3K0=TSB3(BEMY)/2.* FFIELD
       SP4K0=TSB4(BEMY)/2.* FFIELD
       TK=TK0
       T1K=TPK0
       T2K=TPPK0
       T3K=TP3K0
       T4K=TP4K0
       SK=SK0
       S1K=SPK0
       S2K=SPPK0
       S3K=SP3K0
       S4K=SP4K0
c ----  prediction of PCREST (crest phase) based on the actual coefficients factors T and S
       PCREST=ATAN(-SK0/TK0)
       DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       IF(DDWC.LT.0.) PCREST=PCREST+PI
c --- the routine crest is computing the equivalent field lenght: EQVL
       call crest(bemy,eqvl,xpos,sqcttf,ffield)
C ---  follow computations of the parameters:
c            equivalent field length   (cm)
c            asociated drift length    (cm)
c            slip of phase             (rad)
c            energy gain               (MeV)
c            phase jump                (rad)
c            average k (=frequency/velocity) (cm-1)
c            transit time coefficients (MeV,cm)
c            crest phase               (rd)
c            phase of RF at entrance   (rd)
       SAPHI=PCREST
c --- start iterations: improve the average velocity and the transit time factors
       DO 770 IT=1,3
         DTS=(TK*T1K+SK*S1K)/(TK*TK+SK*SK)
         FK1=2.*DTS
         EQVLP=EQVL
c --- computation of the phase slip: PHSLIP
         PHSLIP=-4.*ATAN(3.2*DTS/EQVL)
         IF(PHSLIP.NE.0.) THEN
           TIL2=PHSLIP/2.
           do iiii=1,4
             GX=1./TAN(TIL2)-1./TIL2-FK1/EQVLP
             GPX=-1./(SIN(TIL2)*SIN(TIL2)) + 1./(TIL2*TIL2)
             TIL2=TIL2-GX/GPX
             HX=1./TAN(TIL2)-1./TIL2
             EQVLP=FK1/HX
           enddo
           PHSLIP=TIL2*2.
         ENDIF
         PEQVL=XPOS
         ASDL=PEQVL-EQVL/2.
         F0=XITL0(GAMREF,GAMS,BEMY,SAPHI,AQST)
         DELWRM=(F0-GAMREF)*XMAT
         enrs=enref+delwrm
         gams=enrs/xmat
         bets=sqrt(1.-1./(gams*gams))
C ---  computation of the jump of phase:DELPHR
         COEPH =FH*AQST/(VL*XMAT)
         F3=XITL3(GAMREF,GAMS,BEMY,IT,SAPHI,AQST)
         DELPHR= COEPH * F3
         XK2=FH0/BETS
         XKM=DELPHR/EQVL + XK2*(1.+ASDL/EQVL) - XK1*ASDL/EQVL
         BEMY=FH0/XKM
c --- computations of transit time factors
         TK=TTA0(BEMY)/2.  * FFIELD
         T1K=TTA1(BEMY)/2. * FFIELD
         SK=TSB0(BEMY)/2.  * FFIELD
         S1K=TSB1(BEMY)/2. * FFIELD
770    CONTINUE
c --- crest phase PCREST (after iterations)
       PCREST=ATAN(-SK/TK)
       DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       IF(DDWC.LT.0.) PCREST=PCREST+PI
       DCEMD=0.
c --- isec = 0: the phase offset is given relative to the entrance of the cavity
c --- isec = 1: the phase offset is given relative to the middle of the cavity
c --- (vapmi has been computed in the function XITL3 (in deg) )
       if(isec.ne.0) then
        IMEDI=INT(VAPMI/360.+.4)
        DCEMD=VAPMI-360.*IMEDI
        WRITE(16,773) NRRES,ENMIL,VAPMI,DCEMD
773     FORMAT(2X,' AT THE MIDDLE OF THE CAVITY:',i4,/,2x,
     * ' *ENERGY :',E12.5,' MEV  *PHASE :',E12.5,' DEG',2X,
     * ' *SLIP OF PHASE  :',e12.5,' deg',/)
        DCEMD=DCEMD*PI/180.
      endif
C      SAPHI: phase of RF w.r.t. the phase offset
       DPHASE=DPHASE*PI/180.
       SAPHI=PCREST+DPHASE+TTVOL-DCEMD
       DDW=AQST*(TK0*COS(SAPHI)-SK0*SIN(SAPHI))
       enrs=enref+ddw
       gams=enrs/xmat
C --- start new iterations to improve:
c         transit time factors
c         phase slip
c         jump of phase
c         crest phase
c         phase of RF
c
       SAPHO=SAPHI
       DO 78 IT=1,3
         DTS=(TK*T1K+SK*S1K)/(TK*TK+SK*SK)
         FK1=2.*DTS
         EQVLP=EQVL
         PHSLIP=-4.*ATAN(3.2*DTS/EQVL)
         IF(PHSLIP.NE.0.) THEN
           TIL2=PHSLIP/2.
           do iiii=1,3
             GX=1./TAN(TIL2)-1./TIL2-FK1/EQVLP
             GPX=-1./(SIN(TIL2)*SIN(TIL2)) + 1./(TIL2*TIL2)
             TIL2=TIL2-GX/GPX
             HX=1./TAN(TIL2)-1./TIL2
             EQVLP=FK1/HX
           enddo
           PHSLIP=TIL2*2.
         ENDIF
         PEQVL=XPOS
         ASDL=PEQVL-EQVL/2.
C --- ENERGY GAIN AND PHASE JUMP
         F0=XITL0(GAMREF,GAMS,BEMY,SAPHI,AQST)
         DELWRM=(F0-GAMREF)*XMAT
         ENRS=ENREF+DELWRM
         gams=enrs/xmat
         bets=sqrt(1.-1./(gams*gams))
         xk2=fh0/bets
C ---    DELPHR: jump of phase
         COEPH =FH*AQST/(VL*XMAT)
         F2=XITL3(GAMREF,GAMS,BEMY,IT,SAPHI,AQST)
         DELPHR= COEPH * F2
C ---   XKM:  average k =  frequency/velocity)
         XKM=DELPHR/EQVL + XK2*(1.+ASDL/EQVL) - XK1*ASDL/EQVL
         BEMY=FH0/XKM
C ---  TRANSIT TIME FACTORS (based on the velocity BEMY)
         TK=TTA0(BEMY)/2.  * FFIELD
         T1K=TTA1(BEMY)/2. * FFIELD
         T2K=TTA2(BEMY)/2. * FFIELD
         T3K=TTA3(BEMY)/2. * FFIELD
         T4K=TTA4(BEMY)/2. * FFIELD
         SK=TSB0(BEMY)/2.  * FFIELD
         S1K=TSB1(BEMY)/2. * FFIELD
         S2K=TSB2(BEMY)/2. * FFIELD
         S3K=TSB3(BEMY)/2. * FFIELD
         S4K=TSB4(BEMY)/2. * FFIELD
         PCREST=ATAN(-SK/TK)
         DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
         IF(DDWC.LT.0.) PCREST=PCREST+PI
         DPHII=(XK1-XK2)*EQVL/10.+(XKP1+XKP2)/120.*EQVL*EQVL
     X         + (XK1-XKM)*ASDL
c --- phase of RF
c*2011*July*16 DPHII sign
c*2011 SV        SAPHI=PCREST+DPHASE+TTVOL-DCEMD+DPHII
c*2011 ET        SAPHI=PCREST+DPHASE+TTVOL-DCEMD-DPHII
       SAPHI=PCREST+DPHASE+TTVOL-DCEMD+DPHII
78     CONTINUE
       savph=saphi*180./pi
C      REFERENCE PARTICLE BASED ON CHARGE STATE: QMOY = ABS(QST)
       CFH=FH/(VL*2.*XMAT)
       CKH=QMOY*QMOY/(4.*XMAT*XMAT)
       DPHII=(XK1-XK2)*EQVL/10.+(XKP1+XKP2)/120.*EQVL*EQVL
     X         + (XK1-XKM)*ASDL
       PHARES=SAPHI+XK2*YLG+DELPHR
c   ***** TEST
c       DITEMP=(2.*XK1*ASDL+XKM*(EQVL-ASDL)+
c     x        XK2*(YLG-(EQVL+ASDL)))/FH+DELPHR/FH
ccc       DITEMP=(XK1*ASDL+XKM*(EQVL)+
ccc     x        XK2*(YLG-(EQVL+ASDL)))/FH+DELPHR/FH
       TREFS=TREF+(XK2*YLG+DELPHR)/FH
ccc       TREFS1=TREF+DITEMP
ccc       write(6,*) 'trefs trefs1 ',trefs,trefs1
ccc       pause
c ****************************
       PHARED=(PHARES-SAPHI)*180./PI
       TREDG=fh*TREFS *180./PI
C ******  REFERENCE BASED ON THE AVERAGE CHARGE STATE QMOY (if several charges state)
c *******       CFH=FH/(VL*2.*XMAT)
c ******       CKH=QMOY*QMOY/(4.*XMAT*XMAT)
c ******  save the energy and the T.O.F of the particle reference at the input of the cavity
c ******       enri=f(7,ngdrf)
c *****       trefi=f(6,ngdrf)
c *****       call gap(gamref,saphi,gams,delphr)
c ***** trefs and enrs: time of flight and energy of the reference and the output of the cavity
c *****       trefs=f(6,ngdrf)
c *****       enrs=f(7,ngdrf)
c *****       grefs=f(7,ngdrf)/xmat
c *****       bets=sqrt(grefs*grefs-1.)/grefs
c *****       dwrfs=enrs-enri
c *****       PHARES=SAPHI+XK2*YLG+DELPHR
c *****       TREFS=TREF+(XK2*YLG+DELPHR)/FH
c *****       TREDG=fh*TREFS *180./PI
       write(16,*) ' PARAMETERS RELATING TO THE REFERENCE PARTICLE '
       write(16,*) '***********************************************'
       write(16,*) ' ENERGY GAIN(MeV) ',DELWRM,' TOF ',tredg,' DEG'
       write(16,*) ' PHASE JUMP(DG) ',DELPHR*180./PI
       write(16,*) ' SLIP OF PHASE AT THE INPUT(DG) ',SAPHO*180./PI
       write(16,*) ' PHASE OF RF AT ENTRANCE(DG) ',savph
       write(16,*) ' AVERAGE k (cm-1) (freq./velocity): ',XKM
       write(16,*) ' Associated drift length ',asdl,' (cm)'
       write(16,*) ' Equivalent field length ',eqvl,' cm center at ',
     *               xpos,' cm'
       write(16,*) ' TRANSIT TIME FACTORS AND DERIVATIVES (MeV,cm):'
       write(16,*) ' T ',TK,T1K,T2K,T3K,T4K
       WRITE(16,*) ' S ',SK,S1K,S2K,S3K,S4K
       write(16,*) ' PHASE SLIP(DEG) ',PHSLIP*180./PI
       write(16,*) ' CREST PHASE OF RF (DEG) ',PCREST*180./PI
       WRITE(16,*) ' MAGNITUDE ',SQCTTF,' MV/cm'
       t0s=sqrt(tk*tk+sk*sk)
       WRITE(16,*) ' T0 ',T0S
c *************************************************************************
       call gap(gamref,saphi,gams,delphr)
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       CALL EXT2D(1)
c print in file 'short.data'
c 3.12.09       phnew=-(int(tcog*fh/pi+0.5)-tcog*fh/pi)*180.
c 3.12.09       dav1(idav,6)=encog-xmat
c 3.12.09        dav1(idav,7)=phnew
       if(itvol) then
        dav1(idav,38)=dphete
        dav1(idav,39)=dphase*180./pi
       else
        dav1(idav,38)=dphete
       endif
       WRITE(16,3777)
3777   FORMAT(/,' Dynamics at the output',/,
     1 5X,'   BETA     dW(MeV)     ENERGY(MeV) ',
     2 '   TOF(deg)     TOF(sec)')
       engain=encog-enold
       WRITE(16,3473) bets,DELWRM,enrs-xmat,fh*TREFS *180./PI,TREFS
3473   FORMAT(' REF ',f7.5,3x,f10.6,3x,e12.5,3x,e12.5,3x,e12.5)
       WRITE(16,1789) bcog,engain,encog-xmat,tcog*fh*180./pi,tcog
1789   FORMAT(' COG ',f7.5,3x,f10.6,3x,e12.5,3x,e12.5,3x,e12.5)
       TESTCA=exten(1)*exten(2)*exten(3)
c       epsil=1.E-20
       epsil=1.E-40
       IF(abs(TESTCA).gt.epsil) THEN
         QDISP=2.*SQRT(exten(1))
         QMD=exten(1)*exten(3)-exten(2)**2
         SQMDV=4.*PI*SQRT(QMD)
         SURM=4.*PI*SQRT(QMD)*180./PI
         QDP=2.*SQRT(exten(3))
         COR12=exten(2)/sqrt(exten(1)*exten(3))
comment         PENT12=SQRT(exten(1)/exten(3))/COR12
comment         PENT21=SQRT(exten(1)/exten(3))*COR12
         QDPDE=QDP*180./PI
       ELSE
         QDISP=0.
         QMD=0.
         SQMDV=0.
         SURM=0.
         QDP=0.
         COR12=0.
         PENT12=0.
         PENT21=0.
         QDPDE=0.
       ENDIF
       TRQTX=exten(4)*exten(5)-exten(8)**2
       TRQPY=exten(6)*exten(7)-exten(9)**2
       QDITAX=2.*SQRT(exten(4))
       QDIANT=2.*SQRT(exten(5))
       QDITAY=2.*SQRT(exten(6))
       QDIANP=2.*SQRT(exten(7))
       SURXTH=4.*PI*SQRT(TRQTX)
       SURYPH=4.*PI*SQRT(TRQPY)
       IF(SHIFT) THEN
         vref=bets*vl
         tref=trefs
       ELSE
         vref=bcog*vl
         tref=tcog
       ENDIF
       if(itvol) then
        ttvols=tref
comment        dphete=dgphr
comment        attvol=fh*ttvols*180./pi
comment        write(16,7456) ottvol,attvol
       endif
comment 7456   format(2x,'***tof at the input: ',e12.5,' deg',/,
comment     *           2x,'***tof at the output: ',e12.5,' deg')
       call statis
C      ENVEL
       CALL STAPL(dav1(idav,24))
comment       WRITE(16,9998) SQMDV
comment9998   FORMAT(2X,'   EMITTANCE (norm): ',
comment     *        E12.5,' PI*MEV*RAD')
       dav1(idav,16)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
c 3.12.09       dav1(idav,17)=surxth*10./pi
c
       dav1(idav,21)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
       dav1(idav,25)=nrres
       dav1(idav,30)=ngood
c
c   print in the file: 'dynac.dmp':
c   gap number, phase offset(deg), relativistic beta, energy(MeV), horz. emit.(mm*mrd,norm), vert. emit.(mm*mrd,norm),long. emit(keV*sec)
c
c --- dav1(idav,16): Emittance(norm)  x-xp (mm*mrad)
c --- dav1(idav,21): Emittance(norm)  y-yp (mm*mrad)
       emns=1.e12*sqmdv/(pi*fh)
cet2010s
       tcgprt=fh*tcog*180./pi
       trfprt=fh*tref*180./pi
c cavity number, z (m), transmission (%), synchronous phase (deg), time of flight (deg) (reference),
c COG relativistic beta (@ output), COG output energy (MeV), REF relativistic beta (@ output), REF output energy (MeV),
c horizontal emittance (mm.mrad, RMS normalized), vertical emittance (mm.mrad, RMS normalized),
c longitudinal emittance (RMS, ns.keV)
       trnsms=100.*float(ngood)/float(imax)
       if(ncavmc.eq.1) write(50,*) '# cavmc.dmp'
       if(ncavmc.eq.1) write(50,*) '# cav     Z       trans   ',
     *   'PHIs     TOF(COG)    COG        Wcog          TOF(REF)   ',
     *   '    REF        Wref       Ex,RMS,n     Ey,RMS,n     El,RMS'
       if(ncavmc.eq.1) write(50,*) '#  #     (m)       (%)    ',
     *  '(deg)     (deg)      beta       (MeV)          (deg)      ',
     *  '   beta       (MeV)      (mm.mrad)    (mm.mrad)    (ns.keV)'
       write(50,7023) nrres,0.001*dav1(idav,24),trnsms,dphete,
     *  tcgprt,bcog,encog-xmat,trfprt,bets,enrs-xmat,
     *  0.25*dav1(idav,16),0.25*dav1(idav,21),0.25*emns
7023     format(1x,i4,1x,e12.5,1x,f6.2,1x,f7.2,1x,
     *   2(e14.7,1x,f7.5,1x,e14.7,1x),3(e12.5,1x))
cet2010e
       fh=oldfh
C       new magnetic rigidity of the reference
       gref=1./sqrt(1.-bets*bets)
       xmor=xmat*bets*gref
       BORO=33.356*XMOR*1.E-01/AQST
       WRITE(16,*) ilost,' particles lost in cavity ',nrres
       call emiprt(0)
       return
       end
       SUBROUTINE fieldcav(atte)
c.....................................................................
c  read from disk the electromagnetic field in the form (z,E(z))
c  SUPERFISH units: z(m) E(z) (Volt/m)
c  converted to: z(cm)  E(z) (MVolt/cm)
c.....................................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/rfield/ifield
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/mode/eflvl,rflvl
       common/cavnum1/xnh,xpas,ffield,npt
       common/kcell/avrg(15)
       logical ifield
       part=1.
       ifield=.true.
       read(20,*)fhc
cnov02       att=att*(1.+eflvl)
       att=atte
       read(20,*) xspl(1),yspl(1)
       yspl(1)=yspl(1)*att
       xspl(1)=xspl(1)*100.
       npt=2
       do i=2,4000
         read(20,*) xspl(i),yspl(i)
         if (xspl(i).eq.0.) go to 10
         xspl(i)=xspl(i)*100.
         yspl(i)=yspl(i)*att
         npt=npt+1
       enddo
10     continue
       npt=npt-1
c **************************************
        tdep=xspl(1)
        do it=1,npt
         xspl(it)=xspl(it)-tdep
        enddo
c **************************************
       call deriv2(npt)
       xpas=(xspl(3)-xspl(2))/part
       xcour=xspl(1)
       nfpt=1
       do i=1,10000
         if(xcour.gt.xspl(npt) ) go to 20
         yf(i)=spline(npt,xcour)
         xf(i)=xcour
         xcour=xcour+xpas
         nfpt=nfpt+1
       enddo
20     continue
c   *  valero  mars 2006
       xlimf=xspl(npt)
       xf(nfpt)=0.
comment       yf(nfpt)=0.
c look for the number of cells,the limits of the cells, the number of coordinates in each cell
       do i=1,15
         npoint(i)=0
       enddo
       ncell=1
       xlim(ncell)=xf(1)
       do i=2,nfpt
         if (xf(i).eq.0.) then
           ncell=ncell+1
c  * valero mars 2006
comment      xlim(ncell)=xf(i-1)
           xlim(ncell)=xlimf
           go to 30
         endif
         if(yf(i)*yf(i-1).lt.0.) then
           ncell=ncell+1
           xlim(ncell)=xf(i)
         else
           npoint(ncell)=npoint(ncell)+1
         endif
       enddo
30     continue
       flength=xlim(ncell)
       ncell=ncell-1
c       write(16,*) ' ******Read the field of the cavity************'
       write(16,100) ncell,flength,att,fhc
100    format(' Number of cells: ',i3,' field length: ',e12.5,'cm',
     *        ' field factor: ',e12.5,' frequency: ',e12.5,' Hz')
c ******************************************
c  sv 28/10/2015
       call celint
c *****************************************
       do i=1,ncell
         write(16,200) i,xlim(i),xlim(i+1),avrg(i)
       enddo
200    format(' Cell number ',i3,' lower limit ',e12.5,' cm ',
     *        ' upper limit ',e12.5,' cm ','average ',e12.5,' cm')
       return
       end
       SUBROUTINE celint
c ********************************************
c   average position in the cell
c ********************************************
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/consta/vl,pi,xmat,rpel,qst
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/kcell/avrg(15)
       xlpos=0.
       isce=20
       xlcum=0.
      do inc=1,ncell
       xint1=0.
       xint2=0.
        xnh=0.
        xlcel=xlim(inc+1)-xlim(inc)
        xlpos=xlpos+xlcel
        xpas=xlcel/float(isce)
        xnhc=0
       estop=xpas/10.
       istop=0
20     continue
        xpat=xnhc*xpas
        if(xpat.lt.(xlcel-estop)) then
         z0=xnh*xpas
         z1=(xnh+0.20)*xpas
         z2=(xnh+0.40)*xpas
         z3=(xnh+0.60)*xpas
         z4=(xnh+0.80)*xpas
         z5=(xnh+1.00)*xpas
         fpos0=xnh*xpas+xlcum
         fpos1=(xnh+0.20)*xpas+xlcum
         fpos2=(xnh+0.40)*xpas+xlcum
         fpos3=(xnh+0.60)*xpas+xlcum
         fpos4=(xnh+0.80)*xpas+xlcum
         fpos5=(xnh+1.0)*xpas+xlcum
         tspl0=spline(npt,fpos0)
         tspl1=spline(npt,fpos1)
         tspl2=spline(npt,fpos2)
         tspl3=spline(npt,fpos3)
         tspl4=spline(npt,fpos4)
         tspl5=spline(npt,fpos5)
         xspl0=z0*tspl0
         xspl1=z1*tspl1
         xspl2=z2*tspl2
         xspl3=z3*tspl3
         xspl4=z4*tspl4
         xspl5=z5*tspl5
c     integral E(z)*z over the cell
         tspl11=19.*xspl0+75.*xspl1+50.*xspl2+50.*xspl3+75.*xspl4
     *     +19.*xspl5
         xi1=xpas/288. * tspl11
         xint1=xint1+xi1
         xspl0=tspl0
         xspl1=tspl1
         xspl2=tspl2
         xspl3=tspl3
         xspl4=tspl4
         xspl5=tspl5
c     integral E(z) over the cell
         tspl2=19.*xspl0+75.*xspl1+50.*xspl2+50.*xspl3+75.*xspl4
     *     +19.*xspl5
         xi2=xpas/288. * tspl2
         xint2=xint2+xi2
         xnh=xnh+1.
         xnhc=xnhc+1
         go to 20
        endif
        avrg(inc)=xint1/xint2+xlcum
        xlcum=xlcum+xlcel
       enddo
       return
       end
       FUNCTION fcav(xc,nrc)
c  .......................................................................
c       electromagnetic field  at the position xc
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/fcv/jsav
       fcav=0.
       if(nrc.eq.1) then
         j=1
       else
         j=jsav
       endif
       do k=1,npoint(nrc)
         tk=xc-xf(j)
         if(tk.eq.0.)then
           fcav=yf(j)
           go to 10
         endif
         if(tk.lt.0.) then
           a=(xc-xf(j-1))/(xf(j)-xf(j-1))
           b=(xf(j)-xc)/(xf(j)-xf(j-1))
           fcav=b*yf(j-1)+a*yf(j)
           go to 10
         endif
         j=j+1
       enddo
10     continue
       jsav=j
       return
       end
       FUNCTION ta0(betr,nrc)
c  .......................................................................
c    Transit time factor t(k) (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/consta/vl,pi,xmat,rpel,qst
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       common/fcv/jsav
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR+T1(I)*RFONC*COS(XK*XC)
       enddo
       ta0=ar*(xc2-xc1)
       return
       END
       FUNCTION tta0(betr)
c  .......................................................................
c    Transit time factor t(k) (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta0=0.
         do i=1,ncell
           tta0=tta0+ta0(betr,i)
         enddo
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta0=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR+H(I)*RFONC*COS(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA0=TTA0+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       return
       end
       FUNCTION ta1(betr,nrc)
c  .......................................................................
c    Transit time factor dT(k)/dk (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR-T1(I)*XC*RFONC*SIN(XK*XC)
       enddo
       ta1=ar*(xc2-xc1)
       return
       end
       FUNCTION tta1(BETR)
c  .......................................................................
c    Transit time factor dT(k)/dk (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta1=0.
         do i=1,ncell
           tta1=tta1+ta1(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta1=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR-H(I)*XC*RFONC*SIN(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA1=TTA1+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION ta2(betr,nrc)
c  .......................................................................
c    Transit time factor d2T(k)/dk2 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR-T1(I)*XC*XC*RFONC*COS(XK*XC)
       enddo
       TA2=AR*(XC2-XC1)
       return
       END
       FUNCTION tta2(BETR)
c  .......................................................................
c    Transit time factor d2T(k)/dk2 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta2=0.
         do i=1,ncell
           tta2=tta2+ta2(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta2=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR-H(I)*XC*XC*RFONC*COS(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA2=TTA2+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION ta3(betr,nrc)
c  .......................................................................
c    Transit time factor d3T(k)/dk3 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR+T1(I)*XC*XC*XC*RFONC*SIN(XK*XC)
       enddo
       TA3=AR*(XC2-XC1)
       return
       END
       FUNCTION tta3(BETR)
c  .......................................................................
c    Transit time factor d3T(k)/dk2 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta3=0.
         do i=1,ncell
           tta3=tta3+ta3(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta3=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR+H(I)*XC*XC*XC*RFONC*SIN(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA3=TTA3+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION ta4(betr,nrc)
c  .......................................................................
c    Transit time factor d4T(k)/dk4 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR+T1(I)*XC**4*RFONC*COS(XK*XC)
       enddo
       TA4=AR*(XC2-XC1)
       return
       END
       FUNCTION tta4(BETR)
c  .......................................................................
c    Transit time factor d4T(k)/dk4 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta4=0.
         do i=1,ncell
           tta4=tta4+ta4(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta4=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR+H(I)*XC**4*RFONC*COS(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA4=TTA4+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION ta5(betr,nrc)
c  .......................................................................
c    Transit time factor d5T(k)/dk5 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       ar=0.
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         AR=AR-T1(I)*XC**5*RFONC*SIN(XK*XC)
       enddo
       TA5=AR*(XC2-XC1)
       return
       END
       FUNCTION tta5(BETR)
c  .......................................................................
c    Transit time factor d5T(k)/dk5 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tta5=0.
         do i=1,ncell
           tta5=tta5+ta5(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tta5=0.
         ar=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           AR=AR-H(I)*XC**5*RFONC*SIN(XK*XC)
         enddo
         AR=AR*(XC2-XC1)
         TTA5=TTA5+AR
         AR=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb0(BETR,nrc)
c  .......................................................................
c    Transit time factor s(k) (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR+T1(I)*RFONC*SIN(XK*XC)
       enddo
       SB0=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb0(BETR)
c  .......................................................................
c    Transit time factor s(k) (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb0=0.
         do i=1,ncell
           tsb0=tsb0+sb0(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb0=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR+H(I)*RFONC*SIN(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb0=tsb0+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb1(BETR,nrc)
c  .......................................................................
c    Transit time factor dS(k)/dk (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR+T1(I)*XC*RFONC*COS(XK*XC)
       enddo
       SB1=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb1(BETR)
c  .......................................................................
c    Transit time factor ds(k)/dk (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb1=0.
         do i=1,ncell
           tsb1=tsb1+sb1(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb1=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR+H(I)*XC*RFONC*COS(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb1=tsb1+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb2(BETR,nrc)
c  .......................................................................
c    Transit time factor d2S(k)/dk2 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR-T1(I)*XC*XC*RFONC*SIN(XK*XC)
       enddo
       SB2=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb2(BETR)
c  .......................................................................
c    Transit time factor d2S(k)/dk2 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb2=0.
         do i=1,ncell
           tsb2=tsb2+sb2(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb2=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR-H(I)*XC*XC*RFONC*SIN(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb2=tsb2+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb3(BETR,nrc)
c  .......................................................................
c    Transit time factor d3S(k)/dk3 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR-T1(I)*XC*XC*XC*RFONC*COS(XK*XC)
       enddo
       SB3=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb3(BETR)
c  .......................................................................
c    Transit time factor d3S(k)/dk3 (multi-cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb3=0.
         do i=1,ncell
           tsb3=tsb3+sb3(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb3=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR-H(I)*XC*XC*XC*RFONC*COS(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb3=tsb3+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb4(BETR,nrc)
c  .......................................................................
c    Transit time factor d4S(k)/dk5 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       do i=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR+T1(I)*XC**4*RFONC*SIN(XK*XC)
       enddo
       SB4=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb4(BETR)
c  .......................................................................
c    Transit time factor d4S(k)/dk5 (multi cells)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb4=0.
         do i=1,ncell
           tsb4=tsb4+sb4(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb4=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR+H(I)*XC**4*RFONC*SIN(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb4=tsb4+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION sb5(BETR,nrc)
c  .......................................................................
c    Transit time factor d5S(k)/dk5 (single cell)
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/gaus13/H(13),T(13)
       common/gaus17/H1(17),T1(17)
       xk=fhc*2.*pi/(betr*vl)
       xc1=xlim(nrc)
       xc2=xlim(nrc+1)
       br=0.
       DO I=1,17
         XC= (XC2+XC1)/2.+(XC2-XC1)*H1(I)/2.
         rfonc=fcav(xc,nrc)
         BR=BR+T1(I)*XC**5*RFONC*COS(XK*XC)
       enddo
       SB5=BR*(XC2-XC1)
       return
       END
       FUNCTION tsb5(BETR)
c  .......................................................................
c    Transit time factor d5S(k)/dk5
c  .......................................................................
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/rfield/ifield
       common/gaus13/H(13),T(13)
       dimension xmin(15)
       logical ifield
       if(ifield) then
         tsb5=0.
         do i=1,ncell
           tsb5=tsb5+sb5(betr,i)
         enddo
         return
       else
         xk=fh/(betr*vl)
         xmin(1)=0.
         do i=2,ncel+2
           xmin(i)=0.
         enddo
         do i=2,ncel+1
           xmin(i)=ylg*(i-1)/ncel
         enddo
         ipas=1
         tsb5=0.
         br=0.
101      continue
         if(xmin(ipas+1).eq.0.or.ipas.gt.(ncel+2)) return
         do i=1,13
           xc1=xmin(ipas)
           xc2=xmin(ipas+1)
           XC= (XC2+XC1)/2.+(XC2-XC1)*T(I)/2.
           RFONC= FONE(XC)
           BR=BR+H(I)*XC**5*RFONC*COS(XK*XC)
         enddo
         BR=BR*(XC2-XC1)
         tsb5=tsb5+br
         br=0.
         IPAS=IPAS+1
         GO TO 101
       endif
       end
       FUNCTION fone(Z)
c  .......................................................................
c   Electromagnetic field at the longitudinal point (z,0)
c   The field harmonics are  stored in  A(200)
c ........................................................................
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       FONE=0.
       DO J=1,NHARM
         XL=PI*(J-1)/YLG
         FONE=FONE+A(J)*cos(XL*Z)
       enddo
       RETURN
       end
       SUBROUTINE rharm
c   ..................................................
c    the field is in the form of a Fourier series expansion
c   ..................................................
       implicit real*8 (a-h,o-z)
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/TAPES/IN,IFILE,META
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/rfield/ifield
       common/mode/eflvl,rflvl
       logical ifield
       ifield=.false.
c    ylg: length of the field (cm)
c    fh: frequency (hertz)
c    atte: field factor
c    ncel: number of cells in the cavity
       read(in,*) ylg,fh,atte,ncel
       read(in,*) nharm
       read(in,*) (a(i),i=1,nharm)
       write(16,100) ncel,ylg,fh,atte
100    format(' number of cells in the cavity: ',i3,/,
     *  ' field length: ',e12.5, 'cm',/,
     *  ' freq. ',e12.5,' Hertz',' field factor ',e12.5)
       write(16,*) ' number of harmonics: ',nharm
       write(16,200) (a(i),i=1,nharm)
200    format(3(2x,e12.5))
       do i=1,nharm
         a(i)=a(i)*atte
       enddo
       fh=fh*2.*pi
c **********************************************
c  sv 24/10/2015
       xp=ylg/100.
       open(18,file='chemtr.txt',status='unknown')
       zx=0.
       yh=0.
98     continue
       do inh=1,nharm
cc        xl=pi*float(inh-1)/ylg
cc        yh=yh+a(inh)*cos(xl*zx)
        ff=fone(zx)
       enddo
       write(18,99) zx/100.,ff
99     format(3(2x,e12.5))
       zx=zx+xp
       if(zx.lt.ylg) go to 98
       close(18)
c **********************************************
       return
       end
       SUBROUTINE profil
c...................................................................
C Store header and envelopes in binary file for
C graphics post-processor
c....................................................................
       implicit real*8 (a-h,o-z)
cet2010May*start
       COMMON/DPLT/ZDEB,ZFIN,YWMAX,YPMAX,RMSN
       COMMON/TAPES/IN,IFILE,META
       CALL PLPRF1
       CALL PLPRF2
       RETURN
       END
       SUBROUTINE plprf1
c.........................................................................
C     IPRF  : POINTEUR
C     RMSN  :envelope size in multiples of RMS size
C     SPRFX :half horizontal extent (cm)
C     SPRFY :half vertical extent (cm)
C     SPRFW :half energy extent (MeV)
C     SPRFP :half phase extent    (deg)
C     SPRFL :position along Z(m)
c.........................................................................
       implicit real*8 (a-h,o-z)
       dimension xx(3000),yy(3000)
C     RMSN  :envelope size in multiples of RMS size
C     ZDEB : Starting position of the plot
C     ZFIN : End of the plot
       CHARACTER*80 CAR,text
       COMMON/DPLT/ZDEB,ZFIN,YWMAX,YPMAX,RMSN
       COMMON/PLTPRF/SPRFY(3000),SPRFZ(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/TAPES/IN,IFILE,META
       common/prof/car
       READ(IN,69)CAR(1:80)
69     FORMAT(A80)
       READ(IN,*) RMSN
       READ(IN,*) ZDEB,ZFIN
       READ(IN,*) XXMAX,XYMAX,YWMAX,YPMAX
       IPRF=IPRF-1
       IDEB=1
       IFIN=IPRF
       IF(ZFIN.GT.SPRFL(IPRF))ZFIN=SPRFL(IPRF)
       write(16,*) ' ******* PROFIL ***************** '
       write(16,*) ' IPRF ZFIN ',iprf,zfin
       DO I=2,IPRF
         IF((ZDEB.GT.SPRFL(I-1)).AND.(ZDEB.LE.SPRFL(I)))IDEB=I
         IF((ZFIN.GE.SPRFL(I-1)).AND.(ZFIN.LT.SPRFL(I)))IFIN=I-1
       ENDDO
       if (xxmax.le.0.) then
         XXMAX=0.
         DO I=IDEB,IFIN
           IF(0.5*SPRFY(I)*RMSN .GT. XXMAX) XXMAX=0.5*SPRFY(I)*RMSN
         ENDDO
       endif
       if (xymax.le.0.) then
         XYMAX=0.
         DO I=IDEB,IFIN
           IF(0.5*SPRFZ(I)*RMSN.GT.XYMAX) XYMAX=0.5*SPRFZ(I)*RMSN
         ENDDO
       endif
C Store header and envelopes in binary file for
C graphics post-processor
C
C igrtyp is type of graph (igrtyp=3 for x,y envelope plots)
       igrtyp=3
       text='X and Y envelopes   '
       text(21:80)=car(1:60)
       WRITE(66,*) igrtyp
       WRITE(66,*) text
       xx(1)=zdeb
       xx(2)=zfin
       yy(1)=-xymax
       yy(2)=xxmax
       write(16,*) ' XMAX YMAX ',xxmax,xymax
       WRITE(66,*) xx(1),xx(2),yy(1),yy(2)
c       GRADX=XXMAX/3.
c       YFAXE=XXMAX
       icnt=0
       DO I=IDEB,IFIN
         icnt=icnt+1
         xx(icnt)=sprfl(i)
         yy(icnt)=0.5*sprfy(i)*RMSN
       ENDDO
C write envelope coordinates to graphics file
       WRITE(66,*) icnt
       DO I=1,icnt
         WRITE(66,*) xx(i),yy(i)
       ENDDO
       icnt=0
       DO I=IDEB,IFIN
         icnt=icnt+1
         xx(icnt)=sprfl(i)
         yy(icnt)=-0.5*sprfz(i)*RMSN
       ENDDO
C write envelope coordinates to graphics file
       WRITE(66,*) icnt
       DO I=1,icnt
         WRITE(66,*) xx(i),yy(i)
       ENDDO
       RETURN
       END
       SUBROUTINE plprf2
       implicit real*8 (a-h,o-z)
       dimension xx(3000),yy(3000)
       CHARACTER*80 CAR,text
       common/prof/car
       COMMON/DPLT/ZDEB,ZFIN,YWMAX,YPMAX,RMSN
       COMMON/PLTPRF/SPRFY(3000),SPRFZ(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       IDEB=1
       IFIN=IPRF
       DO I=2,IPRF
         IF((ZDEB.GT.SPRFL(I-1)).AND.(ZDEB.LE.SPRFL(I)))IDEB=I
         IF((ZFIN.GE.SPRFL(I-1)).AND.(ZFIN.LT.SPRFL(I)))IFIN=I-1
       ENDDO
       xxmax=ywmax/1000.
       xymax=ypmax
       if (xxmax.le.0.) then
         XXMAX=0.
         DO I=IDEB,IFIN
           IF(0.5*SPRFW(I)*RMSN.GT.XXMAX) XXMAX=0.5*SPRFW(I)*RMSN
         ENDDO
       endif
       if (xymax.le.0.) then
         XYMAX=0.
         DO I=IDEB,IFIN
           IF(0.5*SPRFP(I)*RMSN.GT.XYMAX) XYMAX=0.5*SPRFP(I)*RMSN
         ENDDO
       endif
cet2010May*end
c Store header and envelopes in binary file for
C graphics post-processor
C
C igrtyp is type of graph (igrtyp=4 for dW/W envelope plots)
       igrtyp=4
       text='dW/W envelope       '
       text(21:80)=car(1:60)
       WRITE(66,*) igrtyp
       WRITE(66,*) text
       xx(1)=zdeb
       xx(2)=zfin
       yy(1)=0.
       yy(2)=xxmax*1000.
       WRITE(66,*) xx(1),xx(2),yy(1),yy(2)
       icnt=0
       DO I=IDEB,IFIN
         icnt=icnt+1
         xx(icnt)=sprfl(i)
         yy(icnt)=0.5*sprfw(i)*1000.*RMSN
       ENDDO
C write envelope coordinates to graphics file
       WRITE(66,*) icnt
       DO I=1,icnt
         WRITE(66,*) xx(i),yy(i)
       ENDDO
C Store header and envelopes in binary file for
C graphics post-processor
C
C igrtyp is type of graph (igrtyp=5 for dPHI envelope plots)
       igrtyp=5
       text='dPHI envelope       '
       text(21:80)=car(1:60)
       WRITE(66,*) igrtyp
       WRITE(66,*) text
       xx(1)=zdeb
       xx(2)=zfin
c 14/08/2009       yy(1)=-xymax
       yy(1)=0.
       yy(2)=xymax
       write(16,*) ' dW/WMAX dPhiMAX ',xxmax,xymax
       WRITE(66,*) xx(1),xx(2),yy(1),yy(2)
       icnt=0
       DO I=IDEB,IFIN
         icnt=icnt+1
         xx(icnt)=sprfl(i)
         yy(icnt)=0.5*sprfp(i)*RMSN
       ENDDO
C write envelope coordinates to graphics file
       WRITE(66,*) icnt
       DO I=1,icnt
         WRITE(66,*) xx(i),yy(i)
       ENDDO
       IPRF=IPRF+1
       RETURN
       END
       FUNCTION slope(N,XV)
c   ..........................................
c    first derivative of the spline function
c   ..........................................
       implicit real*8 (a-h,o-z)
       common/spl/x(4000),y(4000),s(5000),p(5000),q(5000)
       do i=2,n
         xtvi=xv-x(i)
         if(xtvi.gt.0.) go to 4
         if(xtvi.lt.0.) go to 2
         if(xtvi.eq.0.00) go to 3
4        continue
       enddo
3       I=I-1
        AVX=X(I+1)-X(I)
        SLOPE=S(I+1)*AVX/3.+S(I)*AVX/6.+(Y(I+1)-Y(I))/AVX
        return
2       I=I-1
        DGX=XV-X(I)
        DDX=X(I+1)-XV
        AVX=X(I+1)-X(I)
        SLOPE=-(S(I)*DDX*DDX)/(2.*AVX)+(S(I+1)*DGX*DGX)/(2.*AVX)
     *          +((Y(I+1)-Y(I))/AVX)-(AVX*(S(I+1)-S(I))/6.)
        return
        end
       FUNCTION spline (N,XV)
C    ..................................................................
C      SPLINE FUNCTION
C    ..................................................................
       implicit real*8 (a-h,o-z)
       COMMON /SPL/X(4000),Y(4000),S(5000),P(5000),Q(5000)
       spline=y(1)
       xtv1=xv-x(1)
       if(xtv1.lt.0.) then
         SPLINE=Y(1)+((Y(2)-Y(1))/(X(2)-X(1))-S(2)*(X(2)-X(1))/6.)
     *          *(XV-X(1))
         return
       endif
       if(xtv1.eq.0.00) then
         SPLINE=Y(1)
         return
       endif
       if(xtv1.gt.0.) then
        xtvn=xv-x(n)
        if(xtvn.eq.0.00) then
          SPLINE=Y(N)
          return
        endif
        if(xtvn.gt.0.) then
         SPLINE=Y(N)+((Y(N)-Y(N-1))/(X(N)-X(N-1))+S(N-1)*(X(N)-X(N-1))
     *           /6.)*(XV-X(N))
         return
        endif
        if(xtvn.lt.0.) then
         do i=2,n
          xtvi=xv-x(i)
          if(xtvi.gt.0.) go to 11
          if(xtvi.lt.0.) go to 2
          if(xtvi.eq.0.) go to 3
11        continue
         enddo
3        spline=y(i)
         return
2        I=I-1
         DGX=XV-X(I)
         DDX=X(I+1)-XV
         AVX=X(I+1)-X(I)
         SPLINE=S(I)*DDX**3/(6.*AVX)+S(I+1)*DGX**3/(6.*AVX)
     *   +(Y(I+1)/AVX-S(I+1)*AVX/6.)*DGX+(Y(I)/AVX-S(I)*AVX/6.)*DDX
         return
        endif
       endif
       end
       SUBROUTINE DERIV2(N)
C   ...................................................................
C     second derivative of spline functions at position (x,y)
C   ...................................................................
       implicit real*8 (a-h,o-z)
       common/spl/x(4000),y(4000),s(5000),p(5000),q(5000)
       AVXN=X(N)-X(N-1)
       AVVXN=X(N-1)-X(N-2)
       AVYN=Y(N)-Y(N-1)
       AVVYN=Y(N-1)-Y(N-2)
       F=AVXN-(AVVXN**2)/AVXN
       P(N-1)=1.
       Q(N-1)=0.
       if(f.ne.0.) then
         P(N-1)=(-2.*AVXN-3.*AVVXN-AVVXN**2/AVXN)/F
         Q(N-1)=6.*(AVYN/AVXN-AVVYN/AVVXN)/F
       endif
       NM1=N-1
       DO J=2,NM1
         I=N-J
         AVX=X(I+1)-X(I)
         AVVX=X(I+2)-X(I+1)
         AVY=Y(I+1)-Y(I)
         AVVY=Y(I+2)-Y(I+1)
         D=2.*(AVX+AVVX)+AVVX*P(I+1)
         P(I)=-AVX/D
         Q(I)=(6.*(AVVY/AVVX-AVY/AVX)-AVVX*Q(I+1))/D
       ENDDO
       AVX1=X(2)-X(1)
       AVVX1=X(3)-X(2)
       G1=(AVVX1/AVX1)+1.-P(2)-(Q(2)/Q(1))
       G2=(AVVX1/(AVX1*P(1)))-(AVVX1/AVX1)-1.+P(2)
       S(1)=(Q(1)*G1)/(P(1)*G2)
       DO I=1,NM1
         S(I+1)=P(I)*S(I)+Q(I)
       ENDDO
       RETURN
       END
       SUBROUTINE DERIF2(N)
C   ...................................................................
C     second derivative of spline functions at position (x,y)
C   ...................................................................
       implicit real*8 (a-h,o-z)
       common/spff/x(400),y(400),s(500),p(500),q(500)
       AVXN=X(N)-X(N-1)
       AVVXN=X(N-1)-X(N-2)
       AVYN=Y(N)-Y(N-1)
       AVVYN=Y(N-1)-Y(N-2)
       F=AVXN-(AVVXN**2)/AVXN
       P(N-1)=1.
       Q(N-1)=0.
       if(f.ne.0.) then
         P(N-1)=(-2.*AVXN-3.*AVVXN-AVVXN**2/AVXN)/F
         Q(N-1)=6.*(AVYN/AVXN-AVVYN/AVVXN)/F
       endif
       NM1=N-1
       DO J=2,NM1
         I=N-J
         AVX=X(I+1)-X(I)
         AVVX=X(I+2)-X(I+1)
         AVY=Y(I+1)-Y(I)
         AVVY=Y(I+2)-Y(I+1)
         D=2.*(AVX+AVVX)+AVVX*P(I+1)
         P(I)=-AVX/D
         Q(I)=(6.*(AVVY/AVVX-AVY/AVX)-AVVX*Q(I+1))/D
       ENDDO
       AVX1=X(2)-X(1)
       AVVX1=X(3)-X(2)
       G1=(AVVX1/AVX1)+1.-P(2)-(Q(2)/Q(1))
       G2=(AVVX1/(AVX1*P(1)))-(AVVX1/AVX1)-1.+P(2)
       S(1)=(Q(1)*G1)/(P(1)*G2)
       DO I=1,NM1
         S(I+1)=P(I)*S(I)+Q(I)
       ENDDO
       RETURN
       END
       FUNCTION splinf (N,XV)
C    ..................................................................
C      SPLINE FUNCTION
C    ..................................................................
       implicit real*8 (a-h,o-z)
       COMMON /SPFF/X(400),Y(400),S(500),P(500),Q(500)
       splinf=y(1)
       xtv1=xv-x(1)
       if(xtv1.lt.0.) then
         SPLINF=Y(1)+((Y(2)-Y(1))/(X(2)-X(1))-S(2)*(X(2)-X(1))/6.)
     *          *(XV-X(1))
         return
       endif
       if(xtv1.eq.0.00) then
         SPLINF=Y(1)
         return
       endif
       if(xtv1.gt.0.) then
        xtvn=xv-x(n)
        if(xtvn.eq.0.00) then
          SPLINF=Y(N)
          return
        endif
        if(xtvn.gt.0.) then
         SPLINF=Y(N)+((Y(N)-Y(N-1))/(X(N)-X(N-1))+S(N-1)*(X(N)-X(N-1))
     *           /6.)*(XV-X(N))
         return
        endif
        if(xtvn.lt.0.) then
         do i=2,n
          xtvi=xv-x(i)
          if(xtvi.gt.0.) go to 11
          if(xtvi.lt.0.) go to 2
          if(xtvi.eq.0.) go to 3
11        continue
         enddo
3        splinf=y(i)
         return
2        I=I-1
         DGX=XV-X(I)
         DDX=X(I+1)-XV
         AVX=X(I+1)-X(I)
         SPLINF=S(I)*DDX**3/(6.*AVX)+S(I+1)*DGX**3/(6.*AVX)
     *   +(Y(I+1)/AVX-S(I+1)*AVX/6.)*DGX+(Y(I)/AVX-S(I)*AVX/6.)*DDX
         return
        endif
       endif
       end
       FUNCTION slopf(N,XV)
c   ..........................................
c    first derivative of the spline function
c   ..........................................
       implicit real*8 (a-h,o-z)
       COMMON /SPFF/X(400),Y(400),S(500),P(500),Q(500)
       do i=2,n
         xtvi=xv-x(i)
         if(xtvi.gt.0.) go to 4
         if(xtvi.lt.0.) go to 2
         if(xtvi.eq.0.00) go to 3
4        continue
       enddo
3       I=I-1
        AVX=X(I+1)-X(I)
        SLOPF=S(I+1)*AVX/3.+S(I)*AVX/6.+(Y(I+1)-Y(I))/AVX
        return
2       I=I-1
        DGX=XV-X(I)
        DDX=X(I+1)-XV
        AVX=X(I+1)-X(I)
        SLOPF=-(S(I)*DDX*DDX)/(2.*AVX)+(S(I+1)*DGX*DGX)/(2.*AVX)
     *          +((Y(I+1)-Y(I))/AVX)-(AVX*(S(I+1)-S(I))/6.)
        return
        end
       SUBROUTINE area(init)
       implicit real*8 (a-h,o-z)
c  ..............................................................................
C    selection of the regions in the space x/a, y/b and z/c. The values a,b and c
c    are the RMS of the bunch in the space x, y, and z
c    The particles lying in a region are affected with the same color.
c  ..............................................................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/part/xc(iptsz),yc(iptsz),zc(iptsz)
       common/tapes/in,ifile,meta
       common/zones/frms(6),nzone
       dimension inzonn(6)
       if(init.eq.1) then
c  set to -1. the last column of the table f(10,iptsz)
         do i=1,ngood
           f(10,i)=-1.
         enddo
c read: ityp and nzone; the number of regions
c       ityp=0, the zones are defined in the 3-d space (x/a,y/b,x/c)
c       ityp=1, the zones are defined in the plane (x/a,y/b)
         read(in,*) ityp,nzone
         write(16,*)'Number of zones selected:',nzone
         if(nzone.gt.5) then
           write(16,*) 'Number of zones is greater than 5 ',nzone
           stop
         endif
c The regions are selected in the space x/xrms,y/yrms and z/zrms
c limits of the zones ; read the upper limits of the regions
         frms(1)=0.
         if(nzone.ne.0) then
           read(in,*) (frms(i), i=2,nzone)
           frms(nzone+1)=100.
           do i=1,nzone
             write(16,*) 'Zone: ',i,' lower limit: ',frms(i),
     *                 ' upper limit:',frms(i+1)
           enddo
         else
           read(in,*) frms(2)
         endif
       endif
       if(nzones.ge.2) then
         trmoy=0.
         do i=1,ngood
           trmoy=trmoy+f(6,i)
         enddo
         trmoy=trmoy/float(ngood)
         xbar=0.
         ybar=0.
         zbar=0.
         imaxx=0
c   Divide by 100. to convert from centimeters to meters
         do np=1,ngood
           gnp=f(7,np)/xmat
           vnp=vl*sqrt(1.-1./(gnp*gnp))
           zc(np)=(trmoy-f(6,np))*vnp
C   convert from mrad to rad
           f3=f(3,np)*1.e-03
           f5=f(5,np)*1.e-03
C   convert from cm   to m
           xc(np)=(f(2,np)+zc(np)*f3)/100.
           yc(np)=(f(4,np)+zc(np)*f5)/100.
           zc(np)=zc(np)/100.
c evaluate xbar , ybar , zbar
           xbar=xbar+xc(np)
           ybar=ybar+yc(np)
           zbar=zbar+zc(np)
         enddo
         eng=float(ngood)
         xbar=xbar/eng
         ybar=ybar/eng
         zbar=zbar/eng
         do np=1,ngood
           xc(np)=xc(np)-xbar
           yc(np)=yc(np)-ybar
           zc(np)=zc(np)-zbar
         enddo
c    evaluate the rms sizes
         xsqsum=0.
         ysqsum=0.
         zsqsum=0.
         do j=1,ngood
           xcj=xc(j)
           ycj=yc(j)
           zcj=zc(j)
           xsqsum=xsqsum+xcj*xcj
           ysqsum=ysqsum+ycj*ycj
           zsqsum=zsqsum+zcj*zcj
         enddo
         xrmsz=xsqsum/float(ngood)
         yrmsz=ysqsum/float(ngood)
         zrmsz=zsqsum/float(ngood)
         xrmsz=sqrt(xrmsz)
         yrmsz=sqrt(yrmsz)
         zrmsz=sqrt(zrmsz)
c  select the particles in the regions and count them
         if(ityp.eq.0) then
c       ityp=0, the zones are defined in the 3-d space (x/a,y/b,x/c)
           do i=1,nzone
             inzonn(i)=0
             do j=1,ngood
               xcp=xc(j)/xrmsz
               ycp=yc(j)/yrmsz
               zcp=zc(j)/zrmsz
               rxyz=sqrt((xcp*xcp+ycp*ycp+zcp*zcp)/3.)
               if(rxyz.lt.frms(i+1).and.rxyz.ge.frms(i)) then
                 inzonn(i)=inzonn(i)+1
                 if(init.eq.1) f(10,j)=frms(i+1)
               endif
               if(f(10,j).eq.100. .and. init.eq.1) f(10,j)=0.
             enddo
             if(init.eq.1) then
               write(16,*) inzonn(i),' particles initially in zone ',i
             else
               write(16,*) inzonn(i),' particles in zone ',i
             endif
           enddo
         else
c       ityp=1, the zones are defined in the plane (x/a,y/b)
           do i=1,nzone
             inzonn(i)=0
             do j=1,ngood
               xcp=xc(j)/xrmsz
               ycp=yc(j)/yrmsz
               zcp=zc(j)/zrmsz
               rxyz=sqrt((xcp*xcp+ycp*ycp)/2.)
               if(rxyz.lt.frms(i+1).and.rxyz.ge.frms(i)) then
                 inzonn(i)=inzonn(i)+1
                 if(init.eq.1) f(10,j)=frms(i+1)
               endif
               if(f(10,j).eq.100. .and. init.eq.1) f(10,j)=0.
             enddo
             if(init.eq.1) then
               write(16,*) inzonn(i),' particles initially in zone ',i
             else
               write(16,*) inzonn(i),' particles in zone ',i
             endif
           enddo
         endif
       endif
       return
       end
       subroutine histgrm
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002)
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/vl,pi,xmat,rpel,qst
       common/tapes/in,ifile,meta
       common/dyn/tref,vref
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/part/xc(iptsz),yc(iptsz),zc(iptsz)
       common/hist/xpos(200),xn(200),ypos(200),yn(200),
     *             zpos(200),zn(200),ixt,iyt,izt
       common/hist1/xps(200),xpn(200),yps(200),ypn(200),
     *             zps(200),zpn(200),ixpt,iypt,izpt
       trmoy=0.
       wcg=0.
       xcg=0.
       xcg=0.
       ycg=0.
       fnstp=100.
       do i=1,ngood
         trmoy=trmoy+f(6,i)
         wcg=wcg+f(7,i)
         xcg=xcg+f(2,i)
         ycg=ycg+f(4,i)
       enddo
       trmoy=trmoy/float(ngood)
       wcg=wcg/float(ngood)
       xcg=xcg/float(ngood)
       ycg=ycg/float(ngood)
       zcg=trmoy*fh
c  Isochronism correction (in case the bunch is not erect)
       xb2x=0.
       xb2z=0.
       xbxz=0.
       imaxx=0
       do np=1,ngood
         gnp=f(7,np)/xmat
         vnp=vl*sqrt(1.-1./(gnp*gnp))
         zc(np)=(trmoy-f(6,np))*vnp/100.
         xc(np)=(f(2,np)-xcg)/100.
         xb2z=xb2z+zc(np)*zc(np)
         xb2x=xb2x+xc(np)*xc(np)
         xbxz=xbxz+zc(np)*xc(np)
         imaxx=imaxx+1
       enddo
       xb2z=xb2z/float(imaxx)
       xb2x=xb2x/float(imaxx)
       xbxz=xbxz/float(imaxx)
       apl=atan(-2.*xbxz/(xb2x-xb2z))/2.
c coordinates of the particles at the point of time position
       xbar=0.
       ybar=0.
       zbar=0.
       imaxx=0
c  Divide by 100. to convert from centimeters to meters
       do np=1,ngood
         gnp=f(7,np)/xmat
         vnp=vl*sqrt(1.-1./(gnp*gnp))
         znp=(trmoy-f(6,np))*vnp
         xnp=f(2,np)
         zc(np)=znp*cos(apl)+xnp*sin(apl)
         xnp=xnp*cos(apl)-znp*sin(apl)
C  convert from mrad to rad
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
C  convert from cm   to m
         xc(np)=(xnp+zc(np)*f3)/100.
         yc(np)=(f(4,np)+zc(np)*f5)/100.
         zc(np)=zc(np)/100.
c evaluate xbar , ybar , zbar
         xbar=xbar+xc(np)
         ybar=ybar+yc(np)
         zbar=zbar+zc(np)
       enddo
       eng=float(ngood)
       xbar=xbar/eng
       ybar=ybar/eng
       zbar=zbar/eng
       do np=1,ngood
         xc(np)=xc(np)-xbar
         yc(np)=yc(np)-ybar
         zc(np)=zc(np)-zbar
       enddo
c  evaluate the rms sizes
       xsqsum=0.
       ysqsum=0.
       zsqsum=0.
       do j=1,ngood
         xsqsum=xsqsum+xc(j)*xc(j)
         ysqsum=ysqsum+yc(j)*yc(j)
         zsqsum=zsqsum+zc(j)*zc(j)
       enddo
       xrmsz=xsqsum/float(ngood)
       yrmsz=ysqsum/float(ngood)
       zrmsz=zsqsum/float(ngood)
       xrmsz=sqrt(xrmsz)
       yrmsz=sqrt(yrmsz)
       zrmsz=sqrt(zrmsz)
c normalize the coordinates x, y and z
       do j=1,ngood
         xc(j)=xc(j)/xrmsz
         yc(j)=yc(j)/yrmsz
         zc(j)=zc(j)/zrmsz
       enddo
c   look for the limits of the cloud of particles in the plane (x/a, y/b, z/c)
c   these limits are assumed included within at most +/- 5 RMS
       fract=5.
       xinf=0.
       yinf=0.
       zinf=0.
c    lower limits
       do i=1,ngood
         if((abs(xc(i)).le.fract).and.(abs(yc(i)).le.fract).and.
     *       (abs(zc(i)).le.fract)) then
           if(xinf.gt.xc(i)) xinf=xc(i)
           if(yinf.gt.yc(i)) yinf=yc(i)
           if(zinf.gt.zc(i)) zinf=zc(i)
         endif
       enddo
c     upper limits
       xsup=xinf
       ysup=yinf
       zsup=zinf
       do i=1,ngood
         if((abs(xc(i)).le.fract).and.(abs(yc(i)).le.fract).and.
     *       (abs(zc(i)).le.fract)) then
           if(xsup.lt.xc(i)) xsup=xc(i)
           if(ysup.lt.yc(i)) ysup=yc(i)
           if(zsup.lt.zc(i)) zsup=zc(i)
         endif
       enddo
c maximal sizes in x, y, and z-directions
       pax=(xsup-xinf)
       pay=(ysup-yinf)
       paz=(zsup-zinf)
c  histogram in x-direction,the step (stepx) is: pax/50
       stepx=pax/fnstp
       do i=1,200
         xn(i)=0.
       enddo
       xtot=0.
       x0=xinf-stepx
       x1=x0+stepx
       j=1
150    continue
c   xpos(j): position of the elementary cylinder j
c   xn(j)  : number of particles in the elementary cylinder j
c   xtot   : total number of particles in x-direction
       if(x1.gt.xsup+stepx) go to 160
       do i=1,ngood
         if(xc(i).gt.x0.and.xc(i).le.x1) xn(j)=xn(j)+1.
       enddo
       xtot=xtot+xn(j)
       xpos(j)=x0+stepx/2.
       j=j+1
       sta=x1
       x1=x1+stepx
       x0=sta
       go to 150
160    continue
c  normalize the number of particles in each step with regard to max.(xn(j))
       j=j-1
       ixt=j
       xnor=0.
       do i=1,j
         if(xnor.lt.xn(i)) xnor=xn(i)
       enddo
       do i=1,j
         xn(i)=xn(i)/xnor
       enddo
c   histogram in y-direction,the step (stepy) is: pay/50
       stepy=pay/fnstp
       do i=1,200
         yn(i)=0.
       enddo
       ytot=0.
       y0=yinf-stepy
       y1=y0+stepy
       j=1
c   ypos(j): position of the step j
c   yn(j)  : number of particles lying in the  j
c   ytot   : total number of particles in y-direction
151    continue
       if(y1.gt.ysup+stepy) go to 161
       do i=1,ngood
         if(yc(i).gt.y0.and.yc(i).le.y1) yn(j)=yn(j)+1.
       enddo
       ytot=ytot+yn(j)
       ypos(j)=y0+stepy/2.
       j=j+1
       sta=y1
       y1=y1+stepy
       y0=sta
       go to 151
161    continue
       j=j-1
       iyt=j
       y0=yinf
c  normalize the number of particles in step with regard to max(yn(j))
       ynor=0.
       do i=1,j
         if(ynor.lt.yn(i)) ynor=yn(i)
       enddo
       do i=1,j
         yn(i)=yn(i)/ynor
       enddo
c  histogram in z-direction,the step (stepz) is: paz/50
       do i=1,200
         zn(i)=0.
       enddo
       stepz=paz/fnstp
       ztot=0.
       z0=zinf-stepz
       z1=z0+stepz
       j=1
c  of length :stepz and radius:ray
c   zpos(j): position of the step j
c   zn(j)  : number of particles in the step j
c   ztot   : total number of particles in z-direction
152    continue
       if(z1.gt.zsup+2.*stepz) go to 162
       do i=1,ngood
         if(zc(i).gt.z0.and.zc(i).le.z1)zn(j)=zn(j)+1.
       enddo
       ztot=ztot+zn(j)
       zpos(j)=z0+stepz/2.
       j=j+1
       sta=z1
       z1=z1+stepz
       z0=sta
       go to 152
162    continue
       j=j-1
       izt=j
c*et*2011 added next line
       z0=zinf
c  normalize the number of particles in each step with regard to max(zn(j))
       znor=0.
       do i=1,j
         if(znor.lt.zn(i)) znor=zn(i)
       enddo
       do i=1,j
         zn(i)=zn(i)/znor
       enddo
c   look for the limits of xp, yp, zp
       xpinf=f(3,1)
       ypinf=f(5,1)
       zpinf=f(6,1)-trmoy
c    lower limits
       do i=1,ngood
         f3=f(3,i)
         f5=f(5,i)
         f6=f(6,i)-trmoy
         if(xpinf.gt.f3) xpinf=f3
         if(ypinf.gt.f5) ypinf=f5
         if(zpinf.gt.f6) zpinf=f6
       enddo
c     upper limits
       xpsup=xpinf
       ypsup=ypinf
       zpsup=zpinf
       do i=1,ngood
         f3=f(3,i)
         f5=f(5,i)
         f6=f(6,i)-trmoy
         if(xpsup.lt.f3) xpsup=f3
         if(ypsup.lt.f5) ypsup=f5
         if(zpsup.lt.f6) zpsup=f6
       enddo
c maximal sizes in xp, yp, and zp-directions
       paxp=(xpsup-xpinf)
       payp=(ypsup-ypinf)
       pazp=(zpsup-zpinf)
c  evaluate the rms sizes in xp, yp, zp
       xpsum=0.
       ypsum=0.
       zpsum=0.
       do i=1,ngood
         f3=f(3,i)
         f5=f(5,i)
         f6=f(6,i)-trmoy
         xpsum=xpsum+f3*f3
         ypsum=ypsum+f5*f5
         zpsum=zpsum+f6*f6
       enddo
       xpsum=xpsum/float(ngood)
       ypsum=ypsum/float(ngood)
       zpsum=zpsum/float(ngood)
       xpsum=sqrt(xpsum)
       ypsum=sqrt(ypsum)
       zpsum=sqrt(zpsum)
c  histogram in xp-direction,the step (stxp) is: paxp/50
       stpx=paxp/fnstp
       do i=1,200
         xpn(i)=0.
       enddo
       xp0=xpinf-stpx
       xp1=xp0+stpx
       j=1
250    continue
c   xps(j): position of the step j
c   xpn(j): number of particles lying in the step j
       if(xp1.gt.xpsup+stpx) go to 260
       do i=1,ngood
         if(f(3,i).gt.xp0.and.f(3,i).le.xp1) xpn(j)=xpn(j)+1.
       enddo
       xps(j)=xp0+stpx/2.
       j=j+1
       sta=xp1
       xp1=xp1+stpx
       xp0=sta
       go to 250
260    continue
c   normalize the number of particles in each step with regard to max.(xpn(j))
       j=j-1
       ixpt=j
       xnor=0.
       do i=1,j
         if(xnor.lt.xpn(i)) xnor=xpn(i)
       enddo
       do i=1,j
         xpn(i)=xpn(i)/xnor
         xps(i)=xps(i)/xpsum
       enddo
c  histogram in yp-direction,the step (stpy) is: payp/50
       stpy=payp/fnstp
       do i=1,200
         ypn(i)=0.
       enddo
       yp0=ypinf-stpy
       yp1=yp0+stpy
       j=1
c   yps(j): position of the step j in mrd
c   ypn(j)  : number of particles lying in the  j
251    continue
       if(yp1.gt.ypsup+stpy) go to 261
       do i=1,ngood
         if(f(5,i).gt.yp0.and.f(5,i).le.yp1) ypn(j)=ypn(j)+1.
       enddo
       yps(j)=yp0+stpy/2.
       j=j+1
       sta=yp1
       yp1=yp1+stpy
       yp0=sta
       go to 251
261    continue
       j=j-1
       iypt=j
       yp0=ypinf
c  normalize the number of particles in step with regard to max(ypn(j))
       ynor=0.
       do i=1,j
         if(ynor.lt.ypn(i)) ynor=ypn(i)
       enddo
       do i=1,j
         ypn(i)=ypn(i)/ynor
         yps(i)=yps(i)/ypsum
       enddo
c  histogram in zp-direction,the step (stpz) is: pazp/50
       do i=1,200
         zpn(i)=0.
       enddo
       stpz=pazp/fnstp
       zp0=zpinf-stpz
       zp1=zp0+stpz
       j=1
c  of length :stepz and radius:ray
c   zps(j): position of the step j (deg)
c   zpn(j)  : number of particles in the step j
252    continue
       if(zp1.gt.zpsup+2.*stpz) go to 262
       do i=1,ngood
         f6=f(6,i)-trmoy
         if(f6.gt.zp0.and.f6.le.zp1) zpn(j)=zpn(j)+1.
       enddo
       zps(j)=zp0+stpz/2.
       j=j+1
       sta=zp1
       zp1=zp1+stpz
       zp0=sta
       go to 252
262    continue
c  normalize the number of particles in each step with regard to max(zn(j))
       j=j-1
       izpt=j
       zpnor=0.
       do i=1,j
         if(zpnor.lt.zpn(i)) zpnor=zpn(i)
       enddo
       do i=1,j
         zpn(i)=zpn(i)/zpnor
         zps(i)=zps(i)/zpsum
       enddo
       return
       end
       SUBROUTINE daves
c
c  ......................................................................
c    Write the characteristics of the beam to the disk
c
c ---- statistics (dE-dPHI)
c    dav1(i,10) :  extension of phase  dPHI (deg)
c    dav1(i,11) :  dispersion of energy dE (MeV)
c old dav1(i,12)=(1.e09)*sqmdv/(fh*pi)
c    dav1(i,12) :  emittance (MeV*rad)
c    dav1(i,23) :  correlation  between dE an dPHI
c
c ---- statistics (x-xp)
c
c    dav1(i,13): x-extension (mm)
c    dav1(i,14): xp-extension (mrad)
c    dav1(i,15): correlation between x and xp
c    dav1(i,16): Emittance(norm)  x-xp (mm*mrad)
c    dav1(i,17): Emittance(non norm) x-xp (mm*mrad)
c
c ---- statistics (y-yp)
c     dav1(i,18): y-extension  (mm)
c     dav1(i,19): yp-extension (mrad)
c     dav1(i,20): correlation between y and yp
c     dav1(i,21) : Emittance(norm) y-yp (mm*mrad)
c     dav1(i,22) : Emittance(non norm) y-yp (mm*mrad)
c
c  with CHASE These parameters are in the array dav2 like above
c
c  ......................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/SPEDA/DAVE,idave
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/cptemit/xltot(maxcell1),nbemit
       COMMON/PORO/IROT1,IROT2
       LOGICAL IROT1,IROT2
       COMMON/SECDR/ISEOR
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/DYNI/VREFI,TREFI,FHINIT
       common/qfkd/ityq
       common/shortl/davprt
       LOGICAL ISEOR
       LOGICAL DAVE,CHASIT,ICHAES,ITVOL,IMAMIN,ityq
       character*80 davprt(maxcell1)
c dave end
       WRITE(12,3334)
3334   FORMAT('*****************************************************')
c34     FORMAT(10A8)
c ********************************************************************
       write(12,*) ' Energies are in [MeV], phases in [deg]',
     * ' lengths in [mm] ,tof in [deg]'
       write(12,*) ' ** For lenses followed by :',
     * ' Cummulative length, element type, length '
c
       write(12,*) ' ** For emit followed by'
       write(12,*) ' * Line 1:',
     * ' Particle reference:  beta, energy, tof ',
     * '  COG: energy, tof, energy offset, tof offset'
       write(12,*) ' * Line 2:',
     * ' COG coordinates for x  xp  y  yp (mm and mrad)'
       write(12,*) ' * Line 3:',
     * ' alpha-x beta-x(mm/mrad) alpha-y beta-y(mm/mrad)',
     * ' alpha-z beta-z(ns/keV)'
       write(12,*) ' * Line 4:',
     * ' alpha-z beta-z(deg/keV) emit-z(non norm.,keV.deg) f(MHz)'
       write(12,*) ' * Line 5:',
     * ' dPHI(deg)  dW(keV)   r12   long. emittance',
     * '  (keV.ns) particles left'
       write(12,*) ' * Line 6:',
     * '   x(mm)    xp(mrad)  r12   hor.  emittance',
     * '  (norm & non norm, mm.mrad)'
       write(12,*) ' * Line 7:',
     * '   y(mm)    yp(mrad)  r12   vert. emittance',
     * '  (norm & non norm, mm.mrad)'
c ********************************************************************
       WRITE(12,3334)
       write(12,*) '            Simulation with ',imax,' particles'
       IF(ISEOR) THEN
        write(12,*) '            Second order transport matrix    '
       ELSE
        write(12,*) '            First order transport matrix    '
       ENDIF
       IF(ICHAES) THEN
         write(12,*) '            Beam intensity  ',beamc,' mA'
         if(iscsp.eq.1) write(12,*)
     *   '            Space charge calculations with HERSC   '
         if(iscsp.eq.2) write(12,*)
     *   '            Space charge calculations with SCHERM  '
         if(iscsp.eq.3.or.iscsp.eq.4) write(12,*)
     *   '            Space charge calculations with SCHEFF  '
         if(sce10.eq.1) write(12,*)
     *   'Space charge calculated for all relevant elements,
     *    but not at drifts'
         if(sce10.eq.2) write(12,*)
     *   'Space charge calculated for accelerating elements only'
         if(sce10.eq.3) write(12,*)
     *   'Space charge calculated for all relevant elements'
       ENDIF
       if(itvol) write(12,*)
     *  ' TOF is operational in accelerating elements  '
       if(imamin) write(12,*)
     *  'Phase adjustments for accelerating elements active'
       WRITE(12,3334)
       ifirst=1
       iit6=0
       nemit=0
       do 10 i=1,idav
         if(davprt(i).ne."") write(12,'(A)') davprt(i)
         if(iitem(i) .eq. 1) then
c  cavity
         n=int(dav1(i,25)+.5)
         if(itvol.and.imamin) then
          write(12,1000) dav1(i,24),n,dav1(i,1),dav1(i,38),dav1(i,39)
1000     format(f9.2,' mm Cavity ',i3,' length ',f7.2,' mm',/,
     *   ' phase offset: before adjustement ',e12.5,' deg',
     *   ' after adjustement ',e14.7,' deg')
         else
          write(12,2789) dav1(i,24),n,dav1(i,1),dav1(i,38)
         endif
2789     format(f9.2,' mm Cavity ',i3,' length ',f7.2,' mm',
     *   ' phase offset: ',e12.5,' deg')
      endif
c
       if(iitem(i) .eq. 2) then
c quadrupole (magnetic)
         write(12,1010) dav1(i,4),dav1(i,1),dav1(i,7),dav1(i,2),
     *   dav1(i,3),dav1(i,5),dav1(i,6),dav1(i,36)
1010     format(f9.2,' mm  Quadrupole:  length = ',e12.5,' mm',
     *   '  aperture radius = ',e12.5,' mm',/,
     *   '   field = ',e12.5,' kG  K2 = ',e12.5,' cm-2 gradient = ',
     *   e12.5,' kG/cm',/,
     *   '   momentum = ',e12.5,' kG.cm  particles left  ',f7.0)
         write(12,*)
       endif
c
       if(iitem(i) .eq. 3) then
c emiprt
         write(12,305)
305      FORMAT('**********          beam (emit card)         ',
     *   '**********')
         nemit=nemit+1
comment         write(13,556) nemit,xltot(nemit)/100.,dav1(i,12)*1000.,
comment     *                 dav1(i,6),dav1(i,30)
comment         write(14,557) nemit,xltot(nemit)/100.,dav1(i,16),dav1(i,21)
comment556      format(2x,i3,3(3x,f8.3),3x,f6.0)
comment557      format(2x,i3,3(3x,f8.3))
c   ****************************************
         write(12,1001) (dav1(i,j), j=3,9)
1001     format(2x,f7.5,4(1x,e14.7),2(2x,e12.5),' MeV-deg')
         write(12,2003) (dav1(i,j), j=31,34)
2003     format(4(2x,f7.3),' mm and mrad ')
c --- following lines describe Courant-Snyder parameters
         fh=dav1(i,40)
c     1)  alpz betz
c
c    1-a) emz: emittance (keV*deg) betz(deg/keV)  gamz(keV/deg)
         emz=dav1(i,12)*1000.*(180./pi)
         betz=0.
         if(emz.gt.1.e-10) betz=dav1(i,10)*dav1(i,10)/emz
         dez=dav1(i,11)*1000.
         gamz=0.
         if(emz.gt.1.e-10) gamz=dez*dez/emz
         alpz=0.
         if(betz*gamz.ge.1.) alpz=sqrt(betz*gamz-1.)
         if(dav1(i,23).gt.0.)alpz=-alpz
c
c   1-b) emzz: emittance dE-dphi (keV*ns) betzz(ns/keV) gamzz(keV/ns)
         emzz=1.e12*dav1(i,12)/fh
         dphizz=1.e09 * dav1(i,10)/fh *(pi/180.)
         betzz=0.
         if(emzz.gt.1.e-10) betzz=dphizz*dphizz/emzz
         gamzz=0.
         if(emzz.gt.1.e-10) gamzz=dez*dez/emzz
         alpzz=0.
         if(betzz*gamzz.ge.1.) alpzz=sqrt(betzz*gamzz-1.)
         if(dav1(i,23).gt.0.)alpzz=-alpzz
c
c    2)  alpx btx
c       betx(mm/mrad)  gamx (mrad/mm)
         betx=0.
         emx=dav1(i,17)
         if(emx.gt.1e-10) betx=dav1(i,13)*dav1(i,13)/emx
         gamx=0.
         if(emx.gt.1e-10) gamx=dav1(i,14)*dav1(i,14)/emx
         alpx=0.
         if(betx*gamx.ge.1.) alpx=sqrt(betx*gamx-1.)
         if(dav1(i,15).gt.0.)alpx=-alpx
c
c    3)  alpy bety
c       bety(mm/mrad)  gamy (mrad/mm)
         bety=0.
         emy=dav1(i,22)
         if(emy.gt.1e-10) bety=dav1(i,18)*dav1(i,18)/emy
         gamy=0.
         if(emy.gt.1e-10) gamy=dav1(i,19)*dav1(i,19)/emy
         alpy=0.
         if(bety*gamy.ge.1.) alpy=sqrt(bety*gamy-1.)
         if(dav1(i,20).gt.0.)alpy=-alpy
c
c    betzz: ns/keV
         write(12,3213) alpx,betx,alpy,bety,alpzz,betzz
3213     format(3(2x,e12.5,2x,e12.5),2x,e12.5)
c    emittance (keV*deg) betz(deg/keV)
         write(12,597)alpz,betz,emz,fh/(2.*pi*1.e6)
597      format(2x,f8.4,2x,e13.6,2x,e13.6,' keV.deg',2x,f8.3,' MHz')
         if (emzz.gt.1000.) then
c ns.MeV
          write(12,6332) dav1(i,10),dez,dav1(i,23),
     *     emzz/1000.,DAV1(I,30)
6332     format(2x,f7.3,1x,f10.2,2x,f8.4,3x,f8.3,' ns.MeV   ',
     *         f7.0,' particles left')
         else
c ns.keV
           write(12,1002) dav1(i,10),dez,dav1(i,23),
     *     emzz,DAV1(I,30)
1002     format(2x,e12.5,2x,f7.2,2x,f8.4,2x,f7.3,' ns.keV   ',
     *         f7.0,' particles left')
         endif
         write(12,1003) (dav1(i,j), j=13,22)
1003     format(2(2x,f7.3,3x,f8.3,2x,f8.4,2x,e12.5,' mm.mrad (norm)',
     *   2x,f7.3,' (non norm)',/))
         write(12,*)
         if (dav1(i,26).eq.1.) then
           write(12,8333) (dav2(i,j), j=31,33)
8333       format('********** With chase',3(1x,f6.4),' **********')
           write(12,1001) (dav1(i,j), j=3,9)
           write(12,2003) (dav2(i,j), j=26,29)
c following lines describe Courant-Snyder parameters
c     1)  alpz betz
c
c    1-a) emz: emittance (keV*deg) betz(deg/keV)  gamz(keV/deg)
         emz=dav2(i,12)*1000.*(180./pi)
         betz=0.
         if(emz.gt.1.e-10) betz=dav2(i,10)*dav2(i,10)/emz
         dez=dav2(i,11)*1000.
         gamz=0.
         if(emz.gt.1.e-10) gamz=dez*dez/emz
         alpz=0.
         if(betz*gamz.ge.1.) alpz=sqrt(betz*gamz-1.)
         if(dav2(i,23).gt.0.) alpz=-alpz
c
c   1-b) emzz: emittance dE-dphi (keV*ns) betzz(ns/keV) gamzz(keV/ns)
         emzz=1.e12*dav2(i,12)/fh
         dphizz=1.e09 * dav2(i,10)/fh *(pi/180.)
         betzz=0.
         if(emzz.gt.1.e-10) betzz=dphizz*dphizz/emzz
         gamzz=0.
         if(emzz.gt.1.e-10) gamzz=dez*dez/emzz
         alpzz=0.
         if(betzz*gamzz.ge.1.) alpzz=sqrt(betzz*gamzz-1.)
         if(dav2(i,23).gt.0.) alpzz=-alpzz
c
c    2)  alpx btx
c       betx(mm/mrad)  gamx (mrad/mm)
         betx=0.
         emx=dav1(i,17)
         if(emx.gt.1e-10) betx=dav2(i,13)*dav2(i,13)/emx
         gamx=0.
         if(emx.gt.1e-10) gamx=dav2(i,14)*dav2(i,14)/emx
         alpx=0.
         if(betx*gamx.ge.1.) alpx=sqrt(betx*gamx-1.)
         if(dav2(i,15).gt.0.) alpx=-alpx
c
c    3)  alpy bety
c       bety(mm/mrad)  gamy (mrad/mm)
         bety=0.
         emy=dav2(i,22)
         if(emy.gt.1e-10) bety=dav2(i,18)*dav2(i,18)/emy
         gamy=0.
         if(emy.gt.1e-10) gamy=dav2(i,19)*dav2(i,19)/emy
         alpy=0.
         if(bety*gamy.ge.1.) alpy=sqrt(bety*gamy-1.)
         if(dav2(i,20).gt.0.) alpy=-alpy
c
c    betzz: ns/keV
         write(12,3213) alpx,betx,alpy,bety,alpzz,betzz
         if (emzz.gt.1000.) then
c ns.MeV
           write(12,6332) dav2(i,10),dez,dav2(i,23),
     *     emzz/1000.,DAV2(I,30)
         else
c ns.keV
           write(12,1002) dav2(i,10),dez,dav2(i,23),
     *     emzz,DAV2(I,30)
         endif
           ifirst=0
comment           write(15,1556) nemit,xltot(nemit)/100.,dav2(i,12)*1000.
comment1556       format(2x,i3,2(3x,f8.3))
comment           write(18,1557) nemit,xltot(nemit)/100.,dav2(i,16),dav2(i,21)
comment1557       format(2x,i3,3(3x,f8.3))
comment656        format(2x,i3,3x,f7.3)
           write(12,1003) (dav2(i,j), j=13,22)
           write(12,*)
         endif
       endif
       if(iitem(i) .eq. 4) then
c bending magnet
         write(12,1025) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,3)
1025     format(f9.2,' mm bending magnet: central trajectory: '
     * ,f8.2,' mm',/,
     *  '     bend angle: ',f7.3,' deg bending radius: ',e12.5,' mm')
         write(12,1029) dav1(i,16),dav1(i,14),dav1(i,15)
1029     format('     field: ',f7.3,' T  n: ',f8.3,
     *  ' beta: ',f8.3)
         write(12,1026) dav1(i,6),dav1(i,9),dav1(i,7),dav1(i,8),
     *                  dav1(i,5)
1026     format( '  *Entrance ',/,'     pole-face rotation:',
     *   f8.3,'  deg curvature: ',f8.3,' mm',/,
     *   '     fringe field corrections: K1 ',f8.3,' K2 ',f8.3,/,
     *   '     vertical half-aperture: ',f8.3,' mm')
         write(12,1027) dav1(i,10),dav1(i,13),dav1(i,11),dav1(i,12),
     *                  dav1(i,17)
1027     format( '   *Exit ',/,'     pole-face rotation  :',
     *   f8.3,' deg  curvature: ',f8.3,' mm',/,
     *   '      fringe field correction: K1 ',f7.3,' K2 ',f7.3,/,
     *   '      vertical half-aperture: ',f8.3,' mm')
         write(12,1028)dav1(i,37)
1028   format('  particles left  ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 5) then
c solenoid
         write(12,5010) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,3),
     *     dav1(i,5),dav1(i,36)
5010     format(f9.2,' mm  Solenoid:  length = ',f8.3,' mm',
     *   '   field = ',e12.5,' kG  K = ',e12.5,' cm-1',/,
     *'      momentum = ',e12.5,' kG.cm   particles left ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 6) then
C From adjust or from entre(not yet introduced) input beam reference
         write(12,303)
303      FORMAT('**********     INITIAL BEAM          **********')
         write(12,2003) (dav1(i,j), j=31,34)
c     1)  alpz betz
c
c    1-a) emz: emittance (keV*deg) betz(deg/keV)  gamz(keV/deg)
         emz=dav1(i,12)*1000.*(180./pi)
         betz=0.
         if(emz.gt.1.e-10) betz=dav1(i,10)*dav1(i,10)/emz
         dez=dav1(i,11)*1000.
         gamz=0.
         if(emz.gt.1.e-10) gamz=dez*dez/emz
         alpz=0.
         if(betz*gamz.ge.1.) alpz=sqrt(betz*gamz-1.)
c
c   1-b) emzz: emittance dE-dphi (keV*ns) betzz(ns/keV) gamzz(keV/ns)
         emzz=1.e12*dav1(i,12)/fhinit
         dphizz=1.e09 * dav1(i,10)/fhinit *(pi/180.)
         betzz=0.
         if(emzz.gt.1.e-10) betzz=dphizz*dphizz/emzz
         gamzz=0.
         if(emzz.gt.1.e-10) gamzz=dez*dez/emzz
         alpzz=0.
         if(betzz*gamzz.ge.1.) alpzz=sqrt(betzz*gamzz-1.)
c
c    2)  alpx btx
c       betx(mm/mrad)  gamx (mrad/mm)
         betx=0.
         emx=dav1(i,17)
         if(emx.gt.1e-10) betx=dav1(i,13)*dav1(i,13)/emx
         gamx=0.
         if(emx.gt.1e-10) gamx=dav1(i,14)*dav1(i,14)/emx
         alpx=0.
         if(betx*gamx.ge.1.) alpx=sqrt(betx*gamx-1.)
c
c    3)  alpy bety
c       bety(mm/mrad)  gamy (mrad/mm)
         bety=0.
         emy=dav1(i,22)
         if(emy.gt.1e-10) bety=dav1(i,18)*dav1(i,18)/emy
         gamy=0.
         if(emy.gt.1e-10) gamy=dav1(i,19)*dav1(i,19)/emy
         alpy=0.
         if(bety*gamy.ge.1.) alpy=sqrt(bety*gamy-1.)
c
c    betzz: ns/keV
         write(12,3213) alpx,betx,alpy,bety,alpzz,betzz
         fh=dav1(i,40)
c    emittance (keV*deg) betz(deg/keV)
         write(12,597)alpz,betz,emz,fh/(2.*pi*1.e6)
         if (emzz.gt.1000.) then
c ns.MeV
          write(12,6332) dav1(i,10),dez,dav1(i,23),
     *     emzz/1000.,DAV1(I,30)
         else
c ns.keV
           write(12,1002) dav1(i,10),dez,dav1(i,23),
     *     emzz,DAV1(I,30)
         endif
         write(12,1003) (dav1(i,j), j=13,22)
         iit6=1
         write(12,304)
304      FORMAT('**********            With chase             ',
     *          '**********')
         write(12,2003) (dav2(i,j), j=26,29)
c     1)  alpz betz
c
c    1-a) emz: emittance (keV*deg) betz(deg/keV)  gamz(keV/deg)
         emz=dav2(i,12)*1000.*(180./pi)
         betz=0.
         if(emz.gt.1.e-10) betz=dav2(i,10)*dav2(i,10)/emz
         dez=dav2(i,11)*1000.
         gamz=0.
         if(emz.gt.1.e-10) gamz=dez*dez/emz
         alpz=0.
         if(betz*gamz.ge.1.) alpz=sqrt(betz*gamz-1.)
c
c   1-b) emzz: emittance dE-dphi (keV*ns) betzz(ns/keV) gamzz(keV/ns)
         emzz=1.e12*dav2(i,12)/fh
         dphizz=1.e09 * dav2(i,10)/fh *(pi/180.)
         betzz=0.
         if(emzz.gt.1.e-10) betzz=dphizz*dphizz/emzz
         gamzz=0.
         if(emzz.gt.1.e-10) gamzz=dez*dez/emzz
         alpzz=0.
         if(betzz*gamzz.ge.1.) alpzz=sqrt(betzz*gamzz-1.)
c
c    2)  alpx btx
c       betx(mm/mrad)  gamx (mrad/mm)
         betx=0.
         emx=dav1(i,17)
         if(emx.gt.1e-10) betx=dav2(i,13)*dav2(i,13)/emx
         gamx=0.
         if(emx.gt.1e-10) gamx=dav2(i,14)*dav2(i,14)/emx
         alpx=0.
         if(betx*gamx.ge.1.) alpx=sqrt(betx*gamx-1.)
c
c    3)  alpy bety
c
c       bety(mm/mrad)  gamy (mrad/mm)
         bety=0.
         emy=dav2(i,22)
         if(emy.gt.1e-10) bety=dav2(i,18)*dav2(i,18)/emy
         gamy=0.
         if(emy.gt.1e-10) gamy=dav2(i,19)*dav2(i,19)/emy
         alpy=0.
         if(bety*gamy.ge.1.) alpy=sqrt(bety*gamy-1.)
c
c    betzz: ns/keV
         write(12,3213) alpx,betx,alpy,bety,alpzz,betzz
         if (emzz.gt.1000.) then
c ns.MeV
           write(12,6332) dav2(i,10),dez,dav2(i,23),
     *     emzz/1000.,DAV2(I,30)
         else
c ns.keV
           write(12,1002) dav2(i,10),dez,dav2(i,23),
     *     emzz,DAV2(I,30)
         endif
           ifirst=0
           write(12,1003) (dav2(i,j), j=13,22)
           write(12,*)
         endif
       if(iitem(i) .eq. 7) then
c drift
         write(12,7010) dav1(i,4),dav1(i,1),dav1(i,36)
7010     format(f9.2,' mm  Drift:  length ',f10.3,' mm ',/,
     *   '   particles left  ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 8) then
c buncher
       if(.not.imamin)
     * write(12,8010) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,3),
     *                dav1(i,36)
8010     format(f9.2,' mm  Buncher ',f9.3,' MV ',
     *   '  RF Phase ',f9.3,' deg  Aperture radius',f5.1,' cm',/,
     *   '   particles left  ',f7.0)
       if(imamin)
     * write(12,8110) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,5),
     *                dav1(i,3),dav1(i,36)
8110     format(f9.2,' mm  Buncher ',f9.3,' MV ',
     *   '  RF Phase ',f9.3,' deg correction ',f9.3,' deg',
     *   '  Aperture radius',f5.1,' cm',/,
     *   '   particles left  ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 9) then
c steerer
         if (dav1(i,2).eq.0.) then
           write(12,8020) dav1(i,3),dav1(i,1)
8020       format(f9.2,' mm  Hor. Mag. Steerer ',E12.5,' Tm ')
         else if (dav1(i,2).eq.1.) then
           write(12,8021) dav1(i,3),dav1(i,1)
8021       format(f9.2,' mm  Ver. Mag. Steerer ',f12.5,' Tm ')
         else if (dav1(i,2).eq.2.) then
           write(12,8022) dav1(i,3),dav1(i,1)
8022       format(f9.2,' mm  Hor. El. Steerer ',E12.5,' kVm/m ')
         else if (dav1(i,2).eq.3.) then
           write(12,8023) dav1(i,3),dav1(i,1)
8023       format(f9.2,' mm  Ver. El. Steerer ',E12.5,' kVm/m ')
         endif
         write(12,*)
       endif
       if(iitem(i) .eq. 10) then
c  sextupole
         write(12,1011) dav1(i,4),dav1(i,1),dav1(i,6),dav1(i,3),
     *                  dav1(i,5),dav1(i,2),
     *                  dav1(i,7),dav1(i,36)
1011     format(f9.2,' mm  Sextupole:  length = ',f7.3,' mm',
     *       '  aperture radius = ',e12.5,' cm',/,
     *       '   field = ',e12.5,' kG  KS2 = ',e12.5,' cm-3',
     *       '   gradient = ',e12.5,' kG/cm2',/,
     *       '   momentum = ',e12.5,' kG.cm   particles left  ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 11) then
c solenoid+quadrupole
         write(12,5011) dav1(i,4),dav1(i,1),dav1(i,7),dav1(i,2),
     *   dav1(i,6),dav1(i,3),dav1(i,5),dav1(i,8),dav1(i,36)
5011     format(f9.2,' mm Sol+Quad: length = ',f7.3,
     *   ' mm aperture radius= ',e12.5,' mm',/,
     *   ' Solenoid: field = ',e12.5,' kG  K = ',e12.5,' cm-1',/,
     *   ' Quadrupole: field ',e12.5,' kG  K2 = ',e12.5,' cm-2',/,
     *   ' momentum = ',e12.5,' kG.cm   particles left ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 12) then
c quadrupole + sextupole
         write(12,5021) dav1(i,4),dav1(i,1),dav1(i,6),dav1(i,2),
     *   dav1(i,3),dav1(i,7),dav1(i,8),dav1(i,10),dav1(i,36)
5021     format(f9.2,' mm  Quad+Sext: length = ',e12.5,' mm ',
     *   ' aperture radius = ',e12.5,' mm',/,
     *   '  Quadrupole: B = ',e12.5,' kG  K2 = ',e12.5,' cm-2',/,
     *   '  Sextupole:  B = ',e12.5,' kG  K2 = ',e12.5, ' cm-3',/,
     *   ' momentum = ',e12.5,' kG.cm   particles left ',f7.0)
         write(12,*)
       endif
       if(iitem(i) .eq. 13) then
c electron gun
         write(12,5031) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,3),
     *   dav1(i,6),dav1(i,36)
5031     format(f9.2,' mm  DC egun  length:',f7.3,' mm',/,
     *   ' Crest field:',f8.3,' MV/m',3x,
     *   ' field stength:',f8.3,' kV',/,
     *   ' beta ( output):',e12.5,' particles left ',f7.0)
c 14/01/2010         write(12,5032) dav1(i,5),dav1(i,7)
c 14/10/2010 5032     format(' beam (emit card) at ',e12.5,
c 14/10/2010     *          ' mm from the cathode, field (MV/m) ',e12.5)
         write(12,*)
       endif
c rfqcl (RFQ cell)
       if(iitem(i) .eq. 14) then
         ncell=int(dav1(i,7))
         write(12,5041) dav1(i,4),ncell,dav1(i,1),dav1(i,2),
     *   dav1(i,3),dav1(i,5),dav1(i,6),dav1(i,36)
5041     format(f9.2,' mm rfq cell:',i5,' length: ',f7.3,' mm',/,
     *   '  V/r02: ',e12.5,' kV/mm**2  AV:',e12.5,' kV ',
     *   '  type: ',f3.0,/,
     *   '  energy(output): ',e12.5,' MeV ','  particles left ',f7.0)
         write(12,*)
       endif
c rfqptq (RFQ)
       if(iitem(i) .eq. 15) then
         ncell=int(dav1(i,7))
         write(12,5042) dav1(i,4),ncell,dav1(i,5),dav1(i,9),
     *   dav1(i,8),dav1(i,6),dav1(i,36)
5042     format(f9.2,' mm rfq: number of cells:',i5,
     *    ' total length: ',e12.5,' mm',/,
     *   '  intervane voltage (reference): ',e12.5,' kV',/,
     *   '  intervane voltage (bunch): ',e12.5,' kV',/,
     *   '  energy (output): ',e12.5,' MeV',/,
     *   '  particles left ',f7.0)
         write(12,*)
      endif
c stipper foils
       if(iitem(i) .eq. 16) then
        write(12,5043) dav1(i,4),dav1(i,1),dav1(i,2),dav1(i,3),
     *                 dav1(i,5),dav1(i,6),dav1(i,36)
5043    format(f9.2,' mm stripper: atomic number: ',f4.0,
     *     '   atomic mass : ',f4.0,'  thickness :',e12.5,
     *     ' g/cm**2',/,4x,'particles charge : ',
     *     f4.0,2x,'energy loss: ',e12.5,' MeV',/,
     *        '    particles left ',f7.0)
         write(12,*)
       endif
c accelerating gap
       if(iitem(i) .eq. 17) then
         n=int(dav1(i,25)+.5)
c 16/12/09         if(itvol.and.imamin) then
c 16/12/09          write(12,1009) dav1(i,24),n,dav1(i,1),dav1(i,2),
c 16/12/09     *                   dav1(i,38),dav1(i,39)
c 16/12/09 1009     format(f9.2,' mm Ac. gap ',i3,' length ',f7.2,
c 16/12/09     *                 'mm field ',e12.5,' kV/mm',/,
c 16/12/09     *   '    phase offset adjusted with TOF',/,
c 16/12/09     *   '    before adjustement ',e12.5,' deg',
c 16/12/09     *   ' after adjustement ',e14.7,' deg')
c 16/12/09         else
          write(12,1008) dav1(i,24),n,dav1(i,1),dav1(i,2),dav1(i,38)
c 16/12/09         endif
1008     format(f9.2,' mm Ac. gap ',i3,' length ',f7.2,' mm field ',
     *   e12.5,' kV/mm  phase of RF (middle): ',e12.5,' deg')
         write(12,*)
       endif
c QUAELEC electric quadrupole
       if(iitem(i) .eq. 18) then
         write(12,2010) dav1(i,4),dav1(i,1),dav1(i,7), dav1(i,2),
     *   dav1(i,6),dav1(i,5),dav1(i,3),dav1(i,36)
2010     format(f9.2,' mm  Quadrupole (electric): length = '
     *   ,e12.5,' mm aperture radius = ',e12.5,' mm',/,
     *   '  voltage = ',f8.3,' kV  K2 = ',e12.5,' cm-2  gradient = ',
     *   e12.5,' kV/(cm*cm) ',/,
     *   '  rigidity = ',e12.5,' kV  particles left  ',f7.0)
         write(12,*)
       endif
c  QUAFK  quadrupole(magnetic or electric)
       if(iitem(i) .eq. 19) then
c  electric quadrupole
        if(ityq) then
         write(12,2110) dav1(i,4),dav1(i,1),dav1(i,7),
     *   dav1(i,2),dav1(i,6),dav1(i,5),dav1(i,3),dav1(i,36)
2110    format(f9.2,' mm  Quadrupole (electric): length = '
     *  ,e12.5,' mm   aperture radius = ',e12.5,' mm',/,' voltage = '
     *  ,e12.5,' kV   K2 = ',e12.5,' cm-2  gradient = ',e12.5,
     * ' kV/cm2',/,' rigidity = ',e12.5,' kV   particles left  ',f7.0)
        else
c  magnetic quadrupole
         write(12,2111) dav1(i,4),dav1(i,1),dav1(i,7),dav1(i,2),
     *   dav1(i,6),dav1(i,5),dav1(i,3),dav1(i,36)
2111     format(f9.2,' mm  Quadrupole (magnetic): length = '
     *  ,e12.5,' mm  aperture radius= ',e12.5,' mm',/,' field = ',e12.5,
     *  ' kG   K2 = ',e12.5,' cm-2  gradient = ',e12.5,' kG/cm',/,
     *  ' momentum = ',e12.5,' kG.cm   particles left  ',f7.0)
       endif
       write(12,*)
      endif
       if(iitem(i) .eq. 20) then
c  rotating the transverse coordinates about the z axis
        write(12,2112) dav1(i,1)
2112    format(9x,'rotating the transverse coordinates',
     *       ' about the z-axis by an angle: ',e12.5,' deg')
       write(12,*)
      endif
c EDFLEC electric deflector
       if(iitem(i) .eq. 21) then
         write(12,3010) dav1(i,4),dav1(i,1),dav1(i,2), dav1(i,3),
     *   dav1(i,5),dav1(i,6),dav1(i,7),dav1(i,8),dav1(i,36)
3010     format(f9.2,' mm  Deflector (electric): length = '
     *   ,e12.5,' mm bend angle = ',e12.5,' deg',/,
     *   '  bend radius = ',e12.5,' mm  radii = ',e12.5,
     *   '  field index = ',e12.5,/,
     *   '  rigidity = ',e12.5,' kV field = ',e12.5,' kV/mm',
     *   '  particles left  ',f7.0)
         write(12,*)
       endif
c   end (big boucle)
10    continue
      return
      end
      SUBROUTINE secord
C     FIRST AND SECOND order MATRIX R AND T
C **** WARNING : IN THIS ROUTINE ALL VARIABLES ARE REAL UNLESS
C                OTHERWIZE DECLARED *********
C
      IMPLICIT REAL*8(A-Z)
C      -- POLE FACE ROTATION
      COMMON /BLOC21/ BE, APB(2), LAYL, LAYX, RABT
      COMMON /BLOC23/ H, DEVI, NB, BDB,L
      COMMON /BLOC11/ R(6,6), T(6,6,6)
      COMMON/SECDR/ISEOR
      LOGICAL ISEOR
C      H      : INVERSE OF THE MAGNET RADIUS (1./CM)
C      NB     : DIMENSIONLESS TRANSPORT n coefficient
C      BDB    : TRANSPORT beta coefficient
C      APB( ) : HALF GAP (CM)
C      L      : LENGTH OF THE MAGNET (CM)
C      BE     : ANGLE OF INLINATION (RAD)
C      LAYL   : K1 TRANSPORT COEFFICIENT
C      LAYX   : K2 TRANSPORT COEFFICIENT
c      RABT   : INVERSE OF THE RADIUS OF CURVATURE (1/CM)
C
C  ENTRANCE OR EXIT OF THE BENDING MAGNET, MATRIX R AND T ARE IN (M,RD)
      ENTRY POFAR1(GAP)
C     SAVE H, GAP AND RABT BEFORE CONVERSION IN M UNITS
      SAH=H
      SARABT=RABT
      SAGAP=GAP
      H=H*100.
      GAP=GAP*1.E-02
      RABT=RABT*100.
      TB = TAN(BE)
      CB = COS(BE)
      SB = (1.0 + SIN(BE)**2)/CB
      TCOR = 2.0*H*GAP*LAYL
      BE1 = BE -TCOR* SB* (  1. -LAYX*TCOR * TB )
      TB1  =  TAN(BE1)
      TB2 = TB**2
      R(2,1) = H*TB
      R(4,3) = - H*TB1
      IF(.NOT.ISEOR) GO TO 3333
      SB=1./CB
      SB2 = SB**2
      SB3 = SB2*  SB
      T(1,1,1)=- 0.5*H*TB2
      T(1,3,3)=0.5*H*SB2
      T(2,1,1)=0.5*H*RABT*SB3- TB*NB*H**2
      T(2,1,2)=H*TB2
      T(2,1,6)=-H*TB
      T(2,3,3)=H**2*(NB+0.5+TB2)*TB-0.5*H*RABT*SB3
      T(2,3,4)=-H*TB2
      T(3,1,3)=H*TB2
      T(4,1,3)=-H*RABT*SB3+2.*H**2*NB*TB
      T(4,1,4)=-H*TB2
      T(4,2,3)=-H*SB2
      SEC2=COS(BE1)*COS(BE1)
      SEC2=1./SEC2
      T(4,3,6)=H*TB-H*TCOR*SEC2
3333  CONTINUE
C     RESTORE H , RABT, APB(2)
      H=SAH
      RABT=SARABT
      GAP=SAGAP
      RETURN
c***************************
      ENTRY POFAR2(GAP)
C     SAVE H, GAP AND RABT BEFORE CONVERSION IN M UNITY
      SAH=H
      SARABT=RABT
      SAGAP=GAP
      H=H*100.
      GAP=GAP*1.E-02
      RABT=RABT*100.
      TB = TAN(BE)
      CB = COS(BE)
      SB = (1.0 + SIN(BE)**2)/CB
      TCOR = 2.0*H*GAP*LAYL
      BE1 = BE -TCOR* SB* (  1. -LAYX*TCOR * TB )
      TB1  =  TAN(BE1)
      TB2 = TB**2
      R(2,1) = H*TB
      R(4,3) = - H*TB1
      IF(.NOT.ISEOR) GO TO 4444
      SB=1./CB
      SB2 = SB**2
      SB3 = SB2*  SB
      T(1,1,1)=0.5*H*TB2
      T(1,3,3)=-0.5*H*SB2
      T(2,1,1)=0.5*H*RABT*SB3- TB*(NB+0.5*TB2)*H**2
      T(2,1,2)=-H*TB2
      T(2,1,6)=-H*TB
      T(2,3,3)=H**2*(NB-0.5*TB2)*TB-0.5*H*RABT*SB3
      T(2,3,4)=H*TB2
      T(3,1,3)=-H*TB2
      T(4,1,3)=-H*RABT*SB3+H**2*(2.*NB+SB2)*TB
      T(4,1,4)=H*TB2
      T(4,2,3)=H*SB2
      SEC2=COS(BE1)*COS(BE1)
      SEC2=1./SEC2
      T(4,3,6)=H*TB-H*TCOR*SEC2
4444  CONTINUE
C     RESTORE H , RABT, APB(2)
      H=SAH
      RABT=SARABT
      GAP=SAGAP
      RETURN
C*******************************
C     IDEAL MAGNET, MATRIX R AND T ARE IN (M,RD)
c*et*2013*Apr*1
      ENTRY BENMAG(sbet,fdtot)
c sbet is relativistic beta of the charge state of interest
      sgam2=1./(1.-sbet*sbet)
C     SAVE H BEFORE CONVERSION IN m
      SAH=H
      H=H*100.
      FIELDN=NB
      BETA=BDB
      AL=L*1.E-02
      RAD=1.0/H
      H2=H*H
      H3=H2*H
      H4=H3*H
      H5=H4*H
      H6=H5*H
      KX2=(1.0-FIELDN)*H2
      KY2=FIELDN*H2
      KX=SQRT(ABS(KX2))
      KY=SQRT(ABS(KY2))
      KX3=KX2*KX
      ARGX=KX*AL
      ARGY=KY*AL
      SINX=SIN(ARGX)
      SHX=SINH(ARGX)
      AL2=AL*AL
      AL3=AL2*AL
      CX=0.
      SX=0.
      DX=0.
      J1XL=0.
      IF (kx2.lt.6.*0) then
       CX=COSH(ARGX)
       SX=SINH(ARGX)/KX
       DX=H*(1.0-CX)/KX2
       J1XL=(ARGX-SHX)/KX3
      ENDIF
      IF (KX2.EQ.6.*0) THEN
       CX=1
       SX=AL
       DX=H*AL2/2.0
       J1XL=AL3/6.0
      ENDIF
      IF (KX2.GT.6.*0) THEN
       CX=COS(ARGX)
       SX=SIN(ARGX)/KX
       DX=H*(1.0-CX)/KX2
       J1XL=(ARGX-SINX)/KX3
      ENDIF
      CPX=-KX2*SX
      SPX=CX
      DPX=H*SX
      CY=0.
      SY=0.
      IF(KY2.LT.6.*0) THEN
       CY=COSH(ARGY)
       SY=SINH(ARGY)/KY
      ENDIF
      IF(KY2.EQ.6.*0) THEN
       CY=1.0
       SY=AL
      ENDIF
      IF(KY2.GT.6.*0) THEN
       CY=COS(ARGY)
       SY=SIN(ARGY)/KY
      ENDIF
      CPY=-KY2*SY
      SPY=CY
C    First order Matrix R
      R(1,1)=CX
      R(1,2)=SX
      R(1,6)=DX
      R(2,1)=CPX
      R(2,2)=SPX
      R(2,6)=DPX
      R(3,3)=CY
      R(3,4)=SY
      R(4,3)=CPY
      R(4,4)=SPY
      R(5,1)=H*SX
      R(5,2)=DX
      R(5,6)=H2*J1XL
c  Carey eq.2.41 p 34
c      R(5,6)=R(5,6)-AL/sgam2
      R(5,6)=R(5,6)-AL*fdtot/sgam2
c*et*2013*Apr*1
c*et*2013*Apr*1 New term definition based on TRACE3D matrix
c original     R(5,6)=H2*J1XL
c      R(5,6)=-(H2/KX2)*(AL*sbet*sbet-SX)+(AL/(sgam2))*(1-H2/KX2)
c      R(5,6)=H2*(ARGX*sbet*sbet-KX*SX)/KX3+(AL/(sgam2))*(1-H2/KX2)
c  ****test
c*et*2013*May*19 Sign definition changed for next 3 lines
c originally next 3 lines commented out;
c      R(5,1)=-R(5,1)
c      R(5,2)=-R(5,2)
c      R(5,6)=-R(5,6)
c ************************
      IF(.NOT.ISEOR) GO TO 3334
C   Second order Matrix T
      COSX=COS(ARGX)
      CHX=COSH(ARGX)
      KX4=KX2*KX2
      KX6=KX4*KX2
      ARGY=KY*AL
      SINY=SIN(ARGY)
      COSY=COS(ARGY)
      SHY=SINH(ARGY)
      CHY=COSH(ARGY)
      KY3=KY2*KY
      KY4=KY2*KY2
      TARGX=ARGX+ARGX
      SIN2X=SIN(TARGX)
      COS2X=COS(TARGX)
      SH2X=SINH(TARGX)
      CH2X=COSH(TARGX)
      KX3=KX*KX2
      KX5=KX3*KX2
      KX7=KX5*KX2
      TARGY=ARGY+ARGY
      SIN2Y=SIN(TARGY)
      COS2Y=COS(TARGY)
      SH2Y=SINH(TARGY)
      CH2Y=COSH(TARGY)
      AL4=AL3*AL
      AL5=AL4*AL
      AL6=AL5*AL
      AL7=AL6*AL
      C=1.0/(KX2-4.0*KY2)
      j1l=0.
      j2l=0.
      j3l=0.
      j2xl=0.
      j3xl=0.
      j4xl=0.
      j5xl=0.
      j7xl=0.
      j9xl=0.
      j10xl=0.
      j11xl=0.
      j12xl=0.
      j13xl=0.
      j14xl=0.
      j15xl=0.
      j16xl=0.
      j17xl=0.
      IF(kx2.gt.6.*0) THEN
       CX=DCOS(KX*AL)
       SX=DSIN(KX*AL)/KX
       DX=H*(1.0D0-CX)/KX2
       J1XL=(ARGX-SINX)/KX3
       J2XL=(1.0-COSX-.5*ARGX*SINX)/KX4
       J3XL=.5*(SINX-ARGX*COSX)/KX3
       J4XL=(.5*ARGX-2.0*SINX/3.0+SIN2X/12.)/KX5
       J5XL=(.25D0-COSX/3.0+COS2X/12.0)/KX4
       J10XL=(ARGX-1.5*SINX+.5*ARGX*COSX)/KX5
       J11XL=(-2.0*ARGX+3.0*SINX-ARGX*COSX)/KX5
       J12XL=(4.0*ARGX-5.5*SINX+1.5*ARGX*COSX)/KX3
       J13XL=(.75-2.0*COSX/3.0-COS2X/12.0-.5*ARGX*SINX)
     <  /KX6
       J14XL=(1.5-4.0*COSX/3.0-COS2X/6.0-ARGX*SINX)/KX6
       J15XL=(-1.75+4.0*COSX/3.0+5.0*COS2X/12.0+1.5*ARGX
     <  *SINX)/KX4
       J16XL=(1.5*ARGX-7.0*SINX/3.0-SIN2X/12.0+ARGX*COSX)/
     <  KX7
       J17XL=(-1.75*ARGX+17.0*SINX/6.0+5.0*SIN2X/24.0-1.5
     <  *ARGX*COSX)/KX5
       J1L=(.5*ARGX-.25*SIN2X)/KX3
       J2L=(.5*ARGX+.25*SIN2X)/KX
       J3L=.25*(1.0-COS2X)/KX2
       if(ky2.lt.6.*0) then
        J7XL=.5*(ARGX-SINX+KX2*C*(SINX-.5*KX*SH2Y/KY))/(KX3*KY2)
        J9XL=C*((COSX-1.0)/KX2+(1.0-CH2Y)/(4.0*KY2))
       endif
       if(ky2.eq.6.*0) then
        J7XL=(2.0*(SINX-ARGX)/KX3+AL3/3.0)/KX2
        J9XL=AL2/(2.0*KX2)+(COSX-1.0)/KX4
       endif
       if(ky2.gt.6.*0) then
        J7XL=.5*(ARGX-SINX+KX2*C*(SINX-.5*KX*SIN2Y/KY))/(KX3*KY2)
        J9XL=C*((COSX-1.0)/KX2+(1.0-COS2Y)/(4.0*KY2))
       endif
c   endif kx2.gt.0.
      endif
      IF (kx2.lt.6.*0) then
       CX=COSH(KX*AL)
       SX=SINH(KX*AL)/KX
       DX=H*(1.0-CX)/KX2
       J1XL=(ARGX-SHX)/KX3
       J2XL=(1.0-CHX+.5*ARGX*SHX)/KX4
       J3XL=.5*(SHX-ARGX*CHX)/KX3
       J4XL=(.5*ARGX-2.0*SHX/3.0+SH2X/12.0)/KX5
       J5XL=(.25-CHX/3.0+CH2X/12.0)/KX4
       J10XL=(ARGX-1.5*SHX+.5*ARGX*CHX)/KX5
       J11XL=(-2.0*ARGX+3.0*SHX-ARGX*CHX)/KX5
       J12XL=(4.0*ARGX-5.5*SHX+1.5*ARGX*CHX)/KX3
       J13XL=(.75-2.0*CHX/3.0-CH2X/12.0+.5*ARGX*SHX)
     <  /KX6
       J14XL=(1.5-4.0*CHX/3.0-CH2X/6.0+ARGX*SHX)/KX6
       J15XL=(-1.75+4.0*CHX/3.0+5.0*CH2X/12.0-1.5*ARGX*SHX)
     <  /KX4
       J16XL=(1.5*ARGX-7.0*SHX/3.0-SH2X/12.0+ARGX*CHX)/
     <  KX7
       J17XL=(-1.75*ARGX+17.0*SHX/6.0+5.0*SH2X/24.0-1.5
     <  *ARGX*CHX)/KX5
       J1L=(.5*ARGX-.25*SH2X)/KX3
       J2L=(.5*ARGX+.25*SH2X)/KX
       J3L=.25*(1.0-CH2X)/KX2
       IF(KY2.LT.6.*0) THEN
         J7XL=.5*(ARGX-SHX+KX2*C*(SHX-.5*KX*SH2Y/KY))/(KX3*KY2)
         J9XL=C*((CHX-1)/KX2+(1.0D0-CH2Y)/(4.0D0*KY2))
       ENDIF
       IF(KY2.EQ.6.*0) THEN
         J7XL=(2.0*(SINX-ARGX)/KX3+AL3/3.0D0)/KX2
         J9XL=AL2/(2.0D0*KX2)+(CHX-1.0D0)/KX4
       ENDIF
       IF(KY2.GT.6.*0) THEN
        J7XL=.5*(ARGX-SHX+KX2*C*(SHX-.5D0*KX*SIN2Y/KY))/(KX3*KY2)
        J9XL=C*((CHX-1.0)/KX2+(1.0-COS2Y)/(4.0*KY2))
       ENDIF
c   end kx2.lt.0.
      endif
      IF (KX2.EQ.6.*0) THEN
       CX=1
       SX=AL
       DX=H*AL*AL/2.0
       J1XL=AL3/6.0
       J2XL=AL4/24.0
       J3XL=AL3/6.0
       J4XL=AL5/60.0
       J5XL=AL4/24.0
       J7XL=(AL3/12.0-AL/(8.0*KY2)-SIN2Y/(16.0*KY3))/KY2
       J10XL=-AL5/24.0
       J11XL=-AL5/60.0
       J12XL=AL3/6.0
       J9XL=AL2/(8.0*KY2)-(1.0-COS2Y)/(16.0*KY4)
       J13XL=AL6/240.0
       J14XL=AL6/1080.0
       J15XL=AL4/12.0
       J16XL=AL7/840.0
       J17XL=AL5/60.0
       J1L=AL3/3.0
       J2L=AL
       J3L=AL2/2.0
c   end kx2.eq.0.
      endif
      CPX=-KX2*SX
      SPX=CX
      DPX=H*SX
      j4l=0.
      j5l=0.
      j6l=0.
      IF (KY2.LT.6.*0) THEN
       CY=COSH(KY*AL)
       SY=SINH(KY*AL)/KY
       J4L=(.5*ARGY-.25*SH2Y)/KY3
       J5L=(.5*ARGY+.25*SH2Y)/KY
       J6L=.25*(1.0-CH2Y)/KY2
      ENDIF
      IF (KY2.EQ.6.*0) THEN
       CY=1.0
       SY=AL
       J4L=AL3/3.0
       J5L=AL
       J6L=AL2/2.0
      ENDIF
      IF (KY2.GT.6.*0) THEN
       CY=COS(KY*AL)
       SY=SIN(KY*AL)/KY
       J4L=(.5*ARGY-.25*SIN2Y)/KY3
       J5L=(.5*ARGY+.25*SIN2Y)/KY
       J6L=.25*(1.0-COS2Y)/KY2
      ENDIF
      CPY=-KY2*SY
      SPY=CY
      SY2=SY*SY
      A=2.0*FIELDN-1.0-BETA
      B=(2.0-FIELDN)
      BN1=2.0*FIELDN-1.0-BETA
      BN2=2.5*FIELDN-BETA-1.5
      BN3=2.0*BETA-FIELDN
c  tabulation of the integrals (notations TRANSPORT SLAC R-75 table VIb)
      I10=DX/H
      I11=0.5*AL*SX
      I111=1.0*(SX**2+DX*RAD)/3.0
      I112=SX*DX*RAD/3.0
      I133=DX/H-(KY2/(KX2-4.0*KY2))*(SY2-2.0*DX*RAD)
      I134=C*(SY*CY-SX)
      I144=(SY2-2.0*DX*RAD)*C
      I20=SX
      I21=(SX+AL*CX)/2.0
      I22=I11
      I211=SX*(1.0+2.0*CX)/3.0
      I212=(2.0*SX**2-DX/H)/3.0
      I222=2.0*SX*DX*RAD/3.0
      I233=SX-2.0*KY2*(SY*CY-SX)*C
      I234=(KX2*DX*RAD-2.0*KY2*SY2)*C
      I244=2.0*C*(SY*CY-SX)
      I33=0.5*AL*SY
      IF(KY2.EQ.6.*0) I34=AL3/6.0
      IF(KY2.NE.6.*0) I34=0.5*(SY-AL*CY)/KY2
      I314=(2.0*SX*CY-SY*(1.0+CX))*C
      I324=C*(2.0*CY*DX*RAD-SX*SY)
      I43=0.5*(SY+AL*CY)
      I44=I33
      I413=C*((KX2-2.0*KY2)*SX*CY-KY2*SY*(1.0+CX))
      I414=C*((KX2-2.0*KY2)*SX*SY-CY*(1.0-CX))
      I424=C*(CY*SX-CX*SY-2.0*KY2*SY*DX*RAD)
      I12=(SX-AL*CX)*0.5/KX2
      I27=(DX*RAD-.5*AL*SX)/KX2
      I313=C*(KX2*CY*DX*RAD-2.0*SX*SY*KY2)
      I314=(2.0*SX*CY-SY*(1.0+CX))*C
      I324=C*(2.0*CY*DX*RAD-SX*SY)
      I43=0.5*(SY+AL*CY)
      I44=I33
      I413=C*((KX2-2.0*KY2)*SX*CY-KY2*SY*(1.0+CX))
      I414=C*((KX2-2.0*KY2)*SX*SY-CY*(1.0-CX))
      I424=C*(CY*SX-CX*SY-2.0*KY2*SY*DX*RAD)
      IF (KX.NE.6.*0) THEN
       I12=(SX-AL*CX)*0.5/KX2
       I27=(DX*RAD-.5*AL*SX)/KX2
       I116=(0.5*AL*SX-(SX**2+DX/H)/3.0)*H/KX2
       I122=(2.0*DX/H-SX**2)/3.0/KX2
       I126=H*(SX+2.0*SX*CX-3.0*AL*CX)/6.0/KX2**2
       I166=H2*(4.0*DX*RAD/3.0+SX**2/3.0-AL*SX)/KX2**2
       I216=H*(AL*CX/2.0+SX/6.0-2.0*SX*CX/3.0)/KX2
       I226=H*(0.5*AL*SX-2.0*SX**2/3.0+DX*RAD/3.0)/KX2
       I266=H2*(SX/3.0+2.0*SX*CX/3.0-AL*CX)/KX2**2
       I323=C*(2.0*KY2*SY*(1.0+CX)/KX2-SX*CY)+SY/KX2
       I336=H*(0.5*AL*SY-C*(CY*(1.0-CX)-2.0*KY2*SX*SY))/KX2
       I346=H*(I34-C*(2.0*SX*CY-SY*(1.0+CX)))/KX2
       I423=C*(2.0*KY2*CY*(1.0+CX)/KX2-CX*CY-KY2*SX*SY)+CY/KX2
       I436=H*(0.5*AL*CY+0.5*SY+C*(KY2*SY*(1.0+CX)
     >-(KX2-2.0*KY2)*SX*CY))/KX2
       I446=H*(AL*SY*0.5-C*((KX2-2.0*KY2)*SX*SY-CY*(1.0-CX)))/KX2
       I26=I12*H
      ELSE
       I12=AL3/6.0
       I27=AL4/12.0
       I116=H*AL4/24.0
       I122=AL4/12.0
       I126=H*AL5/40.0
       I166=H2*AL2/120.0
       I216=H*AL3/6.0
       I226=H*AL4/8.0
       I266=H2*AL5/20.0
       I323=AL2*SY/4.0
       I336=H*AL*(AL2*SY/12.0+(AL*CY-SY)/(KY2*8.0))
       I346=H*AL2*(SY/(KY2*8.0)-AL*CY/(KY2*12.0))
       I423=(AL2*CY+AL*SY)/4.0
       I436=H*AL2*(SY/8.0+CY*AL/12.0)
       I446=H*AL*(AL2*SY/12.0+(SY-AL*CY)/(KY2*8.0))
       I26=I12*H
      ENDIF
C     MATRIX T(i,j,k), TRANSPORT SLAC R-75 table VIa
C
      T(1,1,1)=A*H3*I111+0.5*KX2**2*I122*H
      T(1,1,2)=2.0*A*H3*I112-KX2*H*I112+H*SX
      T(1,1,6)=B*H2*I11+2.0*A*H3*I116-KX2*H2*I122
      T(1,2,2)=A*H3*I122+0.5*H*I111
      T(1,2,6)=B*H2*I12+2.0*A*H3*I126+H2*I112
      T(1,3,3)=BETA*H3*I133-0.5*KY2*H*I10
      T(1,3,4)=2.0*BETA*H3*I134
      T(1,4,4)=BETA*H3*I144-0.5*H*I10
      T(1,6,6)=B*H2*H*I27+A*H3*I166+0.5D0*H3*I122-
     1           H*I10
C
      T(2,1,1)=A*H3*I211+0.5*KX2**2*H*I222-H*CX*CPX
      T(2,1,2)=H*SPX+2.0D0*A*H3*I212-KX2*H*I212
     >       -H*(CX*SPX+CPX*SX)
      T(2,1,6)=B*H2*I21+2.0*A*H3*I216-KX2*H2*I222-H*
     >(CX*DPX+CPX*DX)
      T(2,2,2)=A*H3*I222+0.5*H*I211-H*SX*SPX
      T(2,2,6)=B*H2*I22+2.0*A*H3*I226+H2*I212
     >        -H*(SX*DPX+SPX*DX)
      T(2,3,3)=BETA*H3*I233-0.5*KY2*H*I20
      T(2,3,4)=2.0*BETA*H3*I234
      T(2,4,4)=BETA*H3*I244-0.5*H*I20
      T(2,6,6)=B*H2*I26+A*H3*I266+0.5*H3*I222-H*DX*DPX-H*I20
C
C     VALUE OF "B" IS CHANGED.
C
      B=BETA-FIELDN
C
      T(3,1,3)=2.0*B*H3*I313+KX2*KY2*H*I324
      T(3,1,4)=H*SY+2.0D0*B*H3*I314-KX2*H*I323
      T(3,2,3)=2.0*B*H3*I323-KY2*H*I314
      T(3,2,4)=2.0*B*H3*I324+H*I313
      T(3,3,6)=KY2*I33+2.0*B*H3*I336-KY2*H2*I324
      T(3,4,6)=KY2*I34+2.0*B*H3*I346+H2*I323
C
      T(4,1,3)=2.0*B*H3*I413+KX2*KY2*H*I424-H*CX*CPY
      T(4,1,4)=H*SPY+2.0*H3*B*I414-KX2*H*I423-H*CX*SPY
      T(4,2,3)=2.0*B*H3*I423-KY2*H*I414-H*SX*CPY
      T(4,2,4)=2.0*B*H3*I424+H*I413-H*SX*SPY
      T(4,3,6)=KY2*I43+2.0*B*H3*I436-KY2*H2*I424-H*DX*CPY
      T(4,4,6)=KY2*I44+2.0*B*H3*I446+H2*I423-H*DX*SPY
C
      T(5,1,1)=H4*(BN1*J1XL-BN2*KX2*J4XL)+.5*KX4*J1L
      T(5,1,2)=H4*2.0*BN2*J5XL-KX2*J3L+H*DX
      T(5,1,6)=H5*J11XL+H3*J12XL+H*KX2*J3XL
     <   +H5*2.0*BN2*J4XL+2.0*BETA*H5*J10XL-H*KX2*J1L
c Carey eq. 7.45  p.143
      T(5,1,6)=T(5,1,6)-R(5,1)/sgam2
c
      T(5,2,2)=.5*(H2*J1XL+H4*2.0*BN2*J4XL+J2L)
      T(5,2,6)=-2.0*BETA*H5*J13XL+H5*J14XL+H3*J15XL+H*KX2*
     <  J2XL+H*J3L
c Carey eq. 7.45  p.143
      T(5,2,6)=T(5,2,6)-R(5,2)/sgam2
c
      T(5,3,3)=.5*(H4*(BN3*J1XL-2.0*BETA*KY2*J7XL)+KY4*J4L)
      T(5,3,4)=2.0*BETA*H4*J9XL-KY2*J6L
      T(5,4,4)=BETA*H4*J7XL-.5*(H2*J1XL-J5L)
      T(5,6,6)=(1.0-BETA)*H6*J16XL+H4*J17XL-H2*J3XL
     >        +.5*H2*J1L
c Carey eq. 7.45  p.143
      T(5,6,6)=T(5,6,6)-R(5,6)/sgam2
      T(5,6,6)=T(5,6,6)+AL*((1./sgam2)**2+1.5*sbet*sbet/sgam2)
c Carey eq.7.46 p.143
      T(5,5,6)=-1./sgam2
3334  CONTINUE
C     RESTORE H
      H=SAH
      RETURN
      END
       SUBROUTINE syrout(ii)
c***************************************************************
c  RADIATION  EXITATION FOR PARTICLE ii
C   The method is described in 'synchrotron radiation in DYNAC'
c **************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       COMMON/RADIA/TRT,RMOY,XINTF,CRAE
C*************************************************************************************
C                    Synchrotron radiation
C  COMMON/RADIA/TRT,RMOY,XINTF,CRAE
C   TRT   : t.o.f in sec for a particle crossing the bending magnet
C   RMOY  : average radius of the bending magnet in the magnetic plane
C   XINTF : integral of (u/uc) * S(u/uc) where S(u/uc) is the spectral function
C            ( see Synchrotron radiation in DYNAC)
C   CRAE    : classical electron radius (cm)
C***************************************************************************************
C   PGAM  : Instantaneous radiation power (MeV/sec)
C   ETA   : eta=u/uc
C   uc    :  critical quanta energy (eV)
C   u     :  quanta energy (eV) considered
C           the most significant quanta are assumed emitted between: 0.1<eta<3.
       e4ii=f(7,ii)**4
       gpaii=f(7,ii)/xmat
       cgam=(4.*pi/3.)*crae/(xmat**3)
       pgam=vl*cgam*e4ii/(2.*pi*rmoy*rmoy)
C   ELOST  : Total energy lost due to radiation (MeV)
       elost=xintf*pgam*trt
C      variation of the momentum (only available for relativistic electrons)
       dmo=-elost/f(7,ii)
C      change the total energy (MeV)
       f(7,ii)=f(7,ii)-elost
C      change f(2,ii) and f(3,ii)
C       dmo*r(1,6) is given in m and  dmo*r(2,6) in rad, convert in cm and mrad
       f(2,ii)=f(2,ii)+dmo*r(1,6)*100.
       f(3,ii)=f(3,ii)+dmo*r(2,6)*1000.
       return
       end
       SUBROUTINE syref
c***************************************************************
c  RADIATION  EXITATION FROM PARTICLE REFERENCE
C   The method is described in 'synchrotron radiation in DYNAC'
c **************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RADIA/TRT,RMOY,XINTF,CRAE
       COMMON/DYN/TREF,VREF
       COMMON /BLOC23/ H, DEVI, NB, BDB,L
       common/itvole/itvol,imamin
       common/tofev/ttvols
       logical itvol,imamin
       real*8 L,NB
       beref=vref/vl
       gamref=1./sqrt(1.-beref*beref)
       ener=xmat*gamref*beref
       trt=l/vref
       e4=ener**4
       cgam=(4.*pi/3.)*crae/(xmat**3)
       pgam=vl*cgam*e4/(2.*pi*rmoy*rmoy)
C   ELOST  : Total energy lost due to radiation (MeV)
C   change vref and tref of the reference
       elost=xintf*pgam*trt
       fener=ener-elost
       fgam=fener/xmat
       fberef=sqrt(1.-1./(fgam*fgam))
       vref=fberef*vl
       tref=tref+l/vref
       if(itvol) ttvols=tref
       write(16,250) elost,ener,fener
250    format(//,' REFERENCE AFTER RADIATION EXITATION*****',/,
     *        '    ENERGY LOST (MeV):  ',e12.5,/,
     *        '    OLD ENERGY  (MeV):  ',e12.5,/,
     *        '    NEW ENERGY  (MeV):  ',e12.5)
       return
       end
       SUBROUTINE sextu(imk2,arg,xlsex,rg)
c   ............................................................................
C       magnetic sextupole
c       space charge at the middle of the lens
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       COMMON/FAISC/F(10,iptsz),IMAX,NGOOD
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       LOGICAL IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       LOGICAL SHIFT
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rander/ialin
       LOGICAL IALIN
       common/qskew/qtwist,iqrand,itwist,iaqu
       LOGICAL ITWIST
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/QSEX/L,KQ2,KS2
       REAL*8 L,KQ2,KS2
       character*1 cr
       dimension trans(2)
       ilost=0
       write(16,*) ' ******SEXTUPOLE*********'
C      statistics
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
       NRTRE=NRTRE+1
       cr=char(13)
c print out on terminal of transport element # on one and the same line
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,advance='no') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       len=1
c  if itwist=.true. skews the sextupole before misalignments
       sqtwist=0.
       if(itwist.and.arg.ne.0.) then
         if(iqrand.eq.0) then
           qtwrad=qtwist
           sqtwist=qtwrad
           call zrotat(qtwrad)
         else
           rdcf=.5
           call rlux(trans,len)
           if(trans(1).le.rdcf) sign=-1.
           if(trans(1).gt.rdcf) sign=1.
           call rlux(trans,len)
           qtwrad=qtwist*sign*trans(1)
           sqtwist=qtwrad
           call zrotat(qtwrad)
         endif
       endif
       if(ialin) call randali
       FH0=FH/VL
       write(16,*)'TOF at input:',tref*fh*180./pi,' deg'
c
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       l=xlsex
       bpa=sqrt(1.-1./(gpa*gpa))
       xmco=xmat*bpa*gpa
       ri=33.356*xmco*1.e-01/qst
       if(imk2.ne.0) then
c   arg: field B (kG)
        fb=arg
        b=arg/(rg*rg)
        ks2=b/ri
       else
c  arg: KS2 (cm-3)
        ks2=arg
        b=ks2*ri
        fb=b*(rg*rg)
       endif
       write(16,3300)xlsex,rg,fb,ks2,b,ri
3300   format(' LENGTH = ',e12.3,' cm  APERTURE RADIUS= ',e12.5,' cm',
     *      /,' FIELD = ',e12.5,' kG  KS2 = ',e12.5,' cm-3',
     *        '  GRADIENT = ',e12.5,' kG/(cm*cm)',/,
     *        ' MOMENTUM = ',e12.5,' kG.cm',/)
       call clear
       call elsex
c    print transport matrix (c.o.g.)
       call matrix
c     start daves
       idav=idav+1
       iitem(idav)=10
       dav1(idav,1)=xlsex*10.
       dav1(idav,2)=b
       dav1(idav,3)=fb
       davtot=davtot+xlsex
       dav1(idav,4)=davtot*10.
       dav1(idav,5)=ks2
       dav1(idav,6)=rg
       dav1(idav,7)=ri
C      first half sextupole
       l=xlsex/2.
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KS2  =  B/RI
        CALL ELSEX
        CALL cobeam(ii,l)
C       the evolution of the t.o.f is made in routine cobeam)
comment      f(6,ii)=f(6,ii)+l/(bpai*vl)
       enddo
C      Charge space effect (if dl >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlsex)
           iesp=.false.
         endif
       endif
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- convert wdisp in dp/p (window control)
c ---- ifw = 0 ===> wdisp = dW/W
c ---- ifw = 1 ===> wdisp = dW (MeV)
       if(ifw.eq.0) dispr=gcour*gcour*wdisp/(gcour*(gcour+1.))
       if(ifw.eq.1) dispr=gcour*gcour*wdisp/(gcour*(gcour+1.)*wcg)
c ---- Change the dispersion with the new energy
      do i=1,ngood
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         fd(i)=bpai/bcour * gpai/gcour
       enddo
c Test window after the first half sextupole
       call cogetc
c  Change the t.o.f of the reference
       tref=tref+xlsex/(2.*vref)
       call reject(nlost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       call shuffle
c    second half sextupole
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KS2  =  B/RI
        CALL ELSEX
        CALL COBEAM(II,L)
       enddo
c Test window after the second half sextupole (only in transverse directions and phase)
       call cogetc
c  Change the t.o.f of the reference
       tref=tref+xlsex/(2.*vref)
       call reject(ilost)
       ilost=ilost+nlost
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       dav1(idav,36)=ngood
       write(16,*)'TOF at output:',tref*fh*180./pi,' deg'
       write(16,*)' particles lost in sextupole :',ilost
c  returns coordinates to the initial orientation
       if(itwist.and.b.ne.0.) then
         qtwrad=-sqtwist
         call zrotat(qtwrad)
       endif
       if(iemgrw) call emiprt(0)
       CALL STAPL(davtot*10.)
       RETURN
       END
       SUBROUTINE qalva(bquad,xlqua,rg)
c   ............................................................................
C       QUADRUPOLE
c        space charge computations at the middle of the lens
c   ............................................................................
c   BQUAD: field at pole tip (kG)
c    If BQUAD positive = focalisation in the H-plane (x,z)
C   XLQUA: EFFECTIVE LENGHT (cm )
C   RG:    APERTURE RADIUS (cm)
c
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       COMMON/FAISC/F(10,iptsz),IMAX,NGOOD
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       LOGICAL IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       LOGICAL SHIFT
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rander/ialin
       LOGICAL IALIN
       common/qskew/qtwist,iqrand,itwist,iaqu
       LOGICAL ITWIST
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/QSEX/L,KQ2,KS2
       common/tofev/ttvols
       COMMON/ITVOLE/ITVOL,IMAMIN
       logical itvol,imamin
       REAL*8 L,KQ2,KS2
       dimension trans(1)
       character*1 cr
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       ilost=0
C   B : GRADIENT IN kG/cm
       b=bquad/rg
C      statistics
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
       write(16,*) ' ***QUADRUPOLE (magnetic) ***'
       fcpi=fh*180./pi
       if(itvol) write(16,10) ttvols*fcpi,davtot
10     FORMAT(' ** TOF (input of the lens): ',e12.5,
     *        ' deg at: ',e12.5,' cm in the lattice')
c   if itwist=.true: skews the quadrupole before misalignments (if b ne 0)
c     len = 1 level in rlux routine
       len=1
       if(itwist.and.b.ne.0.) then
         if(iqrand.eq.0) then
           qtwrad=qtwist
           sqtwist=qtwrad
           call zrotat(qtwrad)
         else
           rdcf=.5
           call rlux(trans,len)
           if(trans(1).le.rdcf) sign=-1.
           if(trans(1).gt.rdcf) sign=1.
           call rlux(trans,len)
           qtwrad=qtwist*sign*trans(1)
           sqtwist=qtwrad
           call zrotat(qtwrad)
         endif
       else
         sqtwist=0.
       endif
       if(ialin) call randali
       fh0=fh/vl
c    print transport matrix (cog)
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       l=xlqua
       bpa=sqrt(1.-1./(gpa*gpa))
       xmco=xmat*bpa*gpa
       ri=33.356*xmco*1.e-01/qst
       kq2=b/ri
       call clear
       call elqua
       write(16,3300) xlqua,rg,bquad,kq2,b,ri
3300   format(' LENGTH = ',e12.5,' cm  APERTURE RADIUS= ',e12.5,' cm',
     *      /,' FIELD = ',e12.5,' kG  K2 = ',e12.5,' cm-2 ',
     *        ' GRADIENT = ',e12.5,' kG/cm',/,
     *        ' MOMENTUM = ',e12.5,' kG.cm',/)
       call matrix
c  Daves start
       idav=idav+1
       iitem(idav)=2
       dav1(idav,1)=xlqua*10.
       dav1(idav,2)=bquad
       davtot=davtot+xlqua
       dav1(idav,4)=davtot*10.
       dav1(idav,3)=kq2
       dav1(idav,5)=b
       dav1(idav,6)=ri
       dav1(idav,7)=rg*10.
C      first half quadrupole
       l=xlqua/2.
       do ii=1,ngood
        call clear
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        ri=33.356*xmco*1.e-01/f(9,ii)
        kq2=b/ri
        call elqua
        call cobeam(ii,l)
       enddo
C      Charge space computations (if dl >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlqua)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
c Test window after the first half quadrupole (after s.c. computations)
       call cogetc
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- window control
c  Change the t.o.f of the reference
       tref=tref+xlqua/(2.*vref)
       call reject(nlost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c    second half quadrupole
       do ii=1,ngood
        call clear
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        ri=33.356*xmco*1.e-01/f(9,ii)
        kq2=b/ri
        call elqua
C       ********CONTROL :print matrix for particle number 1 (only for tests)
comment        if(ii.eq.1) then
comment        write(16,*) ' *****second half quadrupole'
comment         xkq2=sqrt(abs(kq2))
comment         xkql=xkq2*l*57.29578
comment         write(16,3300) ri,xkq2,xkql
comment		    call matrix
comment        endif
        call cobeam(ii,l)
       enddo
c Test window after the second half quadrupole
       call cogetc
c  Change the t.o.f of the reference
       tref=tref+xlqua/(2.*vref)
       call reject(ilost)
       ilost=ilost+nlost
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       if(itvol) ttvols=tref
       tcog=0.
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       if(itvol) then
        write(16,11) ttvols*fcpi,davtot,tref*fcpi,tcog*fcpi
11      FORMAT(' ** tof: ',e12.5,
     *        ' deg at: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg tof of the cog: ',
     * e12.5,' deg')
       else
        write(16,12) tref*fcpi,tcog*fcpi
12      FORMAT(' ** tof of the reference: ',e12.5,
     *        ' deg tof of the cog: ',e12.5,' deg')
       endif
       dav1(idav,36)=ngood
       write(16,*)' particles lost in quadrupole :',ilost
c  returns coordinates to the initial orientation
       if(itwist.and.b.ne.0.) then
         qtwrad=-sqtwist
         call zrotat(qtwrad)
       endif
       if(iemgrw) call emiprt(0)
C  envelope
       call stapl(davtot*10.)
       return
       end
       SUBROUTINE qasex(iksq,args,argq,xlqua,rg)
c   ............................................................................
C       quadrupole associated sextupole field
c       space charge computation at the middle of the lens
c    IKSQ: IFLAG (see ARGS and ARGQ)
c    ARGS: strength of SEXTUPOLE
c      IKSQ = 0, then ARGS = KS2 (cm-3), otherwise ARGS = FIELD FS(kG)
c    ARGQ: strength of QUADRUPOLE
c     If IKSQ = 0, then ARGQ = KQ2 (cm-2), otherwise ARGQ = FIELD FQ(kG)
C    XLQUA : EFFECTIVE LENGHT OF THE LENS(cm)
C    RG : APERTURE RADIUS OF THE LENS (cm)
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       COMMON/FAISC/F(10,iptsz),IMAX,NGOOD
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       LOGICAL IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       LOGICAL SHIFT
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rander/ialin
       LOGICAL IALIN
       common/qskew/qtwist,iqrand,itwist,iaqu
       LOGICAL ITWIST
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/QSEX/L,KQ2,KS2
       REAL*8 L,KQ2,KS2
       character*1 cr
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       rgorge=rg
       ilost=0
C      statistics
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
c      if itwist=.true. skews the quasex before misalignments
comment       len=1
comment       if(itwist.and.b.ne.0.) then
comment         if(iqrand.eq.0) then
comment           qtwrad=qtwist
comment           sqtwist=qtwrad
comment           call zrotat(qtwrad)
comment         else
comment           rdcf=.5
comment           call rlux(trans,len)
comment           if(trans.le.rdcf) sign=-1.
comment           if(trans.gt.rdcf) sign=1.
comment           call rlux(trans,len)
comment           qtwrad=qtwist*sign*trans
comment           sqtwist=qtwrad
comment           call zrotat(qtwrad)
comment         endif
comment       endif
c   misalignments
       if(ialin) call randali
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       l=xlqua
       bpa=sqrt(1.-1./(gpa*gpa))
       xmco=xmat*bpa*gpa
       ri=33.356*xmco*1.e-01/qst
       if(iksq.ne.0) then
c  --- quadrupole: argq = field fq (kG)
        fq=argq
        bq=fq/rg
c kq2: strength (cm-2)
        kq2=bq/ri
c --- sextupole: args = field fs (kG)
        fs=args
        bs=fs/(rg*rg)
c ks2 strength (cm-3)
        ks2=bs/ri
       else
c --- quadrupole: argq is kq2 (cm-2)
        kq2=argq
        fq=kq2*rg*ri
        bq=kq2*ri
c --- sextupole: args is ks2 (cm-3)
        ks2=args
        bs=ks2*ri
        fs=bs*(rg*rg)
       endif
       call clear
       call elqsex
       write(16,*) ' *****LENS QUADRUPOLE+SEXTUPOLE ********'
       write(16,3300)xlqua,rg,fq,kq2,bq,fs,ks2,bs,ri
3300   format(' LENS: LENGTH = ',e12.5,' cm APERTURE RADIUS = ',e12.5,
     *  ' cm',/,' QUADRUPOLE: FIELD = ',e12.5,' kG   KQ2 = ',e12.5,
     * ' cm-2  GRADIENT = ',e12.5,' kG/cm',/,
     * ' SEXTUPOLE: FIELD = ',e12.5,' kG   KS2 = ',e12.5,' cm-3',
     * ' GRADIENT = ',e12.5,' kG/(cm*cm)',/,
     * ' MOMENTUM = ',e12.5,' kG.cm',/)
       FH0=FH/VL
       write(16,*)'TOF at input:',tref*fh*180./pi,' deg'
c    print transport matrix  (c.o.g.)
       call matrix
c  start Daves
       idav=idav+1
       iitem(idav)=12
       dav1(idav,1)=xlqua*10.
       dav1(idav,6)=rg*10.
       dav1(idav,2)=fq
       dav1(idav,3)=kq2
c gradient (quad)       dav1(idav,5)=bq
       dav1(idav,7)=fs
       dav1(idav,8)=ks2
c gradient (sol)       dav1(idav,9)=bs
       dav1(idav,10)=ri
       davtot=davtot+xlqua
       dav1(idav,4)=davtot*10.
C      first half quasex
       l=xlqua/2.
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KQ2  =  BQ/RI
        KS2  =  BS/RI
        CALL ELQSEX
C       ********CONTROL :print matrix for particle number 1 (only for tests)
comment        if(ii.eq.5) then
comment         write(16,*) ' *****first half quasex'
comment         xkq2=sqrt(abs(kq2))
comment         xkql=xkq2*l*57.29578
comment         xks=sqrt(abs(ks2))
comment         xksl=xks*l**(1.5)*57.29578
comment         write(16,3300) ri,xkq2,xkql,xks,xksl
comment		 call matrix
comment        endif
C       ********END CONTROL
        CALL cobeam(ii,l)
       enddo
C      Charge space effect (if dl >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlqua)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
c Test window after the first half quasex (after s.c. computations)
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- window control
       call cogetc
c  Change the t.o.f of the reference
       tref=tref+xlqua/(2.*vref)
       call reject(nlost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c    second half quasex
       do ii=1,ngood
        call clear
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        ri=33.356*xmco*1.e-01/f(9,ii)
        kq2  =  bq/ri
        ks2  =  bs/ri
        call elqsex
C       ********CONTROL :print matrix for particle number 1 (only for tests)
comment        if(ii.eq.5) then
comment         write(16,*) ' *****second half quasex'
comment         xkq2=sqrt(abs(kq2))
comment         xkql=xkq2*l*57.29578
comment         xks=sqrt(abs(ks2))
comment         xksl=xks*l**(1.5)*57.29578
comment         write(16,3300) ri,xkq2,xkql,xks,xksl
comment		    call matrix
comment        endif
C       ********END CONTROL
        call cobeam(ii,l)
       enddo
c Test window after the second half quasex (only in transverse directions and phase)
       call cogetc
c  Change the t.o.f of the reference
       tref=tref+xlqua/(2.*vref)
       call reject(ilost)
       ilost=ilost+nlost
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       dav1(idav,36)=ngood
       write(16,*)'TOF at output:',tref*fh*180./pi,' deg'
       write(16,*)' particles lost:',ilost
c  returns coordinates to the initial orientation
comment       if(itwist.and.b.ne.0.) then
comment         qtwrad=-sqtwist
comment         call zrotat(qtwrad)
comment       endif
       if(iemgrw) call emiprt(0)
       CALL STAPL(davtot*10.)
       RETURN
       END
      SUBROUTINE quasex
C...............................................................................
C     FIRST AND SECOND order TRANSPORT MATRIX R AND T
C     quadrupole field combined with a sextupole field
C...............................................................................
C **** WARNING : IN THIS ROUTINE ALL VARIABLES ARE REAL UNLESS
C                OTHERWIZE DECLARED *********
C
      IMPLICIT REAL*8(A-Z)
      COMMON/QSEX/L,KQ2,KS2
      COMMON /BLOC11/ R(6,6), T(6,6,6)
      COMMON/SECDR/ISEOR
      LOGICAL ISEOR
      ENTRY ELQUA
C      PURE QUADRUPOLE FIELD, MATRIX R AND T  ARE IN  (M,RD)
      AKQ2=KQ2*1.E04
      AL=L*1.e-02
      AKQ=SQRT(ABS(AKQ2))
      IF(AKQ.LT.1.E-13) AKQ=1.E-13
      ARG=AKQ*AL
      CSIN=SIN(ARG)
      CCOS=COS(ARG)
      HSIN=SINH(ARG)
      HCOS=COSH(ARG)
      IF(AKQ2.GT.0.0) THEN
       SSOK=CSIN/AKQ
       SC=CCOS
       BSOK=HSIN/AKQ
       BC=HCOS
      ELSE
       SSOK=HSIN/AKQ
       SC=HCOS
       BSOK=CSIN/AKQ
       BC=CCOS
      ENDIF
      R(1,1)=SC
      R(1,2)=SSOK
      R(2,1)=-AKQ2*SSOK
      R(2,2)=SC
      R(3,3)=BC
      R(3,4)=BSOK
      R(4,3)=AKQ2*BSOK
      R(4,4)=BC
      IF(.NOT.ISEOR) RETURN
c      T116=AMAT(12,1)
      T(1,1,6)=AL*AKQ2*SSOK/2.0
c      T126=AMAT(17,1)
      T(1,2,6)=(SSOK-AL*SC)/2.0
c      T216=AMAT(12,2)
      T(2,1,6)=AKQ2*(SSOK+AL*SC)/2.0
c      T226=AMAT(17,2)
      T(2,2,6)=T(1,1,6)
c      T336= AMAT(21,3) T346=AMAT(24,3)
      T(3,3,6)=-AL*AKQ2*BSOK/2.0
      T(3,4,6)=(BSOK-AL*BC)/2.0
c      T436=AMAT(21,4)  T446=AMAT(24,4)
      T(4,3,6)=-AKQ2*(BSOK+AL*BC)/2.0
      T(4,4,6)=T(3,3,6)
c     T511=AMAT(7,5)  T512=AMAT(8,5) T522=AMAT(13,5) T533=AMAT(18,5)
      T(5,1,1)=AKQ2*(AL-SC*SSOK)/4.0
      T(5,1,2)=-AKQ2*SSOK*SSOK/2.0
      T(5,2,2)=(AL+SC*SSOK)/4.0
      T(5,3,3)=-AKQ2*(AL-BC*BSOK)/4.0
c      T534=AMAT(19,5)  T544=AMAT(22,5)
      T(5,3,4)=AKQ2*(BSOK*BSOK)/2.0
      T(5,4,4)=(AL+BC*BSOK)/4.0
      RETURN
      ENTRY ELSEX
C    PURE SEXTUPOLE FIELD, MATRIX R AND T  ARE IN  (M,RD)
      AKS2=KS2*1.E06
      AL=L*1.e-02
1000  CONTINUE
      AKL1=AKS2*AL
      AKL2=AKL1*AL
      AKL3=AKL2*AL
      AKL4=AKL3*AL
      R(1,2)=AL
      R(3,4)=AL
      IF(.NOT.ISEOR) RETURN
C      T522=AMAT(13,5) T544=AMAT(22,5)
      T(5,2,2)=AL/2.0
      T(5,4,4)=AL/2.0
      IF(AKS2.EQ.0.0)RETURN
C     T111=AMAT( 7,1) T112=AMAT( 8,1) T122=AMAT(13,1)
      T(1,1,1)=-AKL2/2.0
      T(1,1,2)=-AKL3/3.0
      T(1,2,2)=-AKL4/12.0
C     T133=AMAT(18,1)  T134=AMAT(19,1)  T144=AMAT(22,1)
      T(1,3,3)= AKL2/2.0
      T(1,3,4)= AKL3/3.0
      T(1,4,4)= AKL4/12.0
C     T211=AMAT( 7,2) T212=AMAT( 8,2) T222=AMAT(13,2)
      T(2,1,1)=-AKL1
      T(2,1,2)=-AKL2
      T(2,2,2)=-AKL3/3.0
C     T233=AMAT(18,2) T234=AMAT(19,2) T244=AMAT(22,2)
      T(2,3,3)= AKL1
      T(2,3,4)= AKL2
      T(2,4,4)= AKL3/3.0
C     T313=AMAT( 9,3) T314=AMAT( 10,3)  T323=AMAT(14,3) T324=AMAT(15,3)
      T(3,1,3)= AKL2
      T(3,1,4)= AKL3/3.0
      T(3,2,3)= AKL3/3.0
      T(3,2,4)= AKL4/6.0
C     T413=AMAT(9,4) T414=AMAT( 10,4) T423=AMAT(14,4) T424=AMAT(15,4)
      T(4,1,3)= AKL1*2.0
      T(4,1,4)= AKL2
      T(4,2,3)= AKL2
      T(4,2,4)= AKL3*2.0/3.0
      RETURN
      ENTRY ELQSEX
C      QUADRUPOLE FIELD	+ SEXTUPOLE FIELD
      AKS2=KS2*1.E06
      AL=L*1.e-02
      AKQ2=KQ2*1.E04
      AKQ=SQRT(ABS(AKQ2))
      ARG=AKQ*AL
      CSIN=SIN(ARG)
      CCOS=COS(ARG)
      HSIN=SINH(ARG)
      HCOS=COSH(ARG)
      IF(AKQ2.GT.6.*0) THEN
       SSOK=CSIN/AKQ
       SC=CCOS
       BSOK=HSIN/AKQ
       BC=HCOS
      ELSE
       SSOK=HSIN/AKQ
       SC=HCOS
       BSOK=CSIN/AKQ
       BC=CCOS
      ENDIF
      IF(AKQ2.EQ.6.*0) THEN
       R(1,2)=AL
       R(3,4)=AL
      ELSE
       R(1,1)=SC
       R(1,2)=SSOK
       R(2,1)=-AKQ2*SSOK
       R(2,2)=SC
       R(3,3)=BC
       R(3,4)=BSOK
       R(4,3)=AKQ2*BSOK
       R(4,4)=BC
      ENDIF
      IF(.NOT.ISEOR) RETURN
comment       IF(KS2.EQ.0.0)RETURN
      IF(AKQ2.EQ.6.*0) GO TO 1000
c      T111=AMAT( 7,1)  T112=AMAT( 8,1) T116=AMAT(12,1)
      T(1,1,1)=-AKS2*(SSOK*SSOK+(1.0-SC)/AKQ2)/3.0
      T(1,1,2)=-2.0*AKS2*(SSOK*(1.0-SC)/AKQ2)/3.0
      T(1,1,6)=AL*AKQ2*SSOK/2.0
c       T122=AMAT(13,1) T126=AMAT(17,1)
      T(1,2,2)=-AKS2*(2.0*(1.0-SC)/AKQ2-SSOK*SSOK)/
     >(3.0*AKQ2)
      T(1,2,6)=(SSOK-AL*SC)/2.0
c      T133=AMAT(18,1) T134=AMAT(19,1) T144=AMAT(22,1)
      T(1,3,3)=AKS2*(BSOK*BSOK+3.0*(1.0-SC)/AKQ2)/5.0
      T(1,3,4)=2.0*AKS2*(BSOK*BC-SSOK)/(5.0*AKQ2)
      T(1,4,4)=AKS2*(BSOK*BSOK-2.0*(1.0-SC)/AKQ2)/
     >(5.0*AKQ2)
c      T211=AMAT( 7,2) T212= AMAT( 8,2) T216=AMAT(12,2)
      T(2,1,1)=-AKS2*(2.0*SSOK*SC+SSOK)/3.0
      T(2,1,2)=-2.0*AKS2*(SC*(1.0-SC)/AKQ2+SSOK*SSOK)/3.0
      T(2,1,6)=AKQ2*(SSOK+AL*SC)/2.0
c      T222=AMAT(13,2)  T226=AMAT(17,2)
      T(2,2,2)=-AKS2*(2.0*SSOK-2.0*SSOK*SC)/(3.0*AKQ2)
      T(2,2,6)=T(1,1,6)
c       T233=AMAT(18,2) T234=AMAT(19,2) T244=AMAT(22,2)
      T(2,3,3)=AKS2*(2.0*BSOK*BC+3.0*SSOK)/5.0
      T(2,3,4)=2.0*AKS2*(BC*BC+BSOK*BSOK*AKQ2-SC)/
     >       (5.0*AKQ2)
      T(2,4,4)=2.0*AKS2*(BSOK*BC-SSOK)/(5.0*AKQ2)
c       T313=AMAT(9,3)   T314=AMAT(10,3)  T323=AMAT(14,3)
      T(3,1,3)=2.0*AKS2*(BC*(1.0-SC)/AKQ2+
     >                2.0*SSOK*BSOK)/5.0
      T(3,1,4)=2.0*AKS2*(2.0*SSOK*BC-BSOK*(1.0+SC))
     >  /(5.0*AKQ2)
      T(3,2,3)=2.0*AKS2*(3.0*BSOK-2.0*BSOK*SC-
     >    SSOK*BC)/(5.0*AKQ2)
c      T324=AMAT(15,3) T336= AMAT(21,3) T346=AMAT(24,3)
      T(3,2,4)=2.0*AKS2*(2.0*BC*(1.0-SC)/AKQ2-
     >     SSOK*BSOK)/(5.0*AKQ2)
      T(3,3,6)=-AL*AKQ2*BSOK/2.0
      T(3,4,6)=(BSOK-AL*BC)/2.0
c      T413=AMAT( 9,4)  T414=AMAT(10,4)  T423=AMAT(14,4)    T424=AMAT(15,4)
      T(4,1,3)=2.0*AKS2*(BSOK*(1.0-SC)+BC*SSOK+
     >  2.0*SC*BSOK+2.0*SSOK*BC)/5.0
      T(4,1,4)=2.0*AKS2*(2.0*SC*BC+2.0*SSOK*BSOK*AKQ2-
     >   BC*(1.0+SC)+BSOK*SSOK*AKQ2)/(5.0*AKQ2)
      T(4,2,3)=2.0*AKS2*(3.0*BC-2.0*BC*SC
     > +2.0*BSOK*SSOK*AKQ2-SC*BC-SSOK*BSOK*AKQ2)/(5.0*AKQ2)
      T(4,2,4)=2.0*AKS2*(2.0*BSOK*(1.0-SC)+
     >  2.0*BC*SSOK-SC*BSOK-SSOK*BC)/(5.0*AKQ2)
c      T436=AMAT(21,4)  T446=AMAT(24,4)
      T(4,3,6)=-AKQ2*(BSOK+AL*BC)/2.0
      T(4,4,6)=T(3,3,6)
c     T511=AMAT(7,5)  T512=AMAT(8,5) T522=AMAT(13,5) T533=AMAT(18,5)
      T(5,1,1)=AKQ2*(AL-SC*SSOK)/4.0
      T(5,1,2)=-AKQ2*SSOK*SSOK/2.0
      T(5,2,2)=(AL+SC*SSOK)/4.0
      T(5,3,3)=-AKQ2*(AL-BC*BSOK)/4.0
c      T534=AMAT(19,5)  T544=AMAT(22,5)
      T(5,3,4)=AKQ2*(BSOK*BSOK)/2.0
      T(5,4,4)=(AL+BC*BSOK)/4.0
      RETURN
      END
       SUBROUTINE solfield(bcret,intgr)
c   ....................................................................
C           Solenoid with an arbitrary magnetic field
c   ....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
       COMMON/RIGID/BORO
       COMMON/DCSPA/IESP
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       common/rander/ialin
       common/femt/iemgrw,iemqesg
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       dimension rs(6,6),rcul(6,6)
       logical iesp,ichaes,shift,ialin,iemgrw
       character*1 cr
c    print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
c       read the magnetic field on the disk
       read(25,*) ncord
       if(ncord.eq.0) then
        rewind(25)
        read(25,*)ncord
       endif
c*et*2013*may*24 fix printout of bcrest(=attenuation) and peak field
       crest=0.
       do i=1,ncord
        read(25,*) xspl(i),yspl(i)
        if(yspl(i).gt.crest) crest=yspl(i)
       enddo
       zinf=xspl(1)
       zsup=xspl(ncord)
       call deriv2(ncord)
       npas=intgr
       dsol=(zsup-zinf)/float(npas)
       xlsol=(zsup-zinf)*100.
       write(16,101) xlsol,crest,bcret
101    FORMAT(5X,'Field length =',F7.3,' cm ',/,
     X  5X,'Crest of the field =',F10.4,' kG' ,/
     X  5X,'Attenuation factor =',F10.4,/)
       FH0=FH/VL
c  random errors in alignment
       if(ialin) call randali
C test window
       ilost=0
C      PLOT
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
c  start prints in file 'short.data'
       idav=idav+1
       iitem(idav)=5
       dav1(idav,1)=xlsol*10.
       dav1(idav,2)=crest
       davti=davtot
       davtot=davtot+xlsol
       dav1(idav,4)=davtot*10.
       zc=zinf
       iesp=.true.
       do ia=1,6
         do ib=1,6
          if(ia.eq.ib) rcul(ia,ib)=1.
          if(ia.ne.ib) rcul(ia,ib)=0.
         enddo
       enddo
       bisol=0.
       bisol2=0.
       fnpas=float(npas)
       do i=1,npas
         zcf=zc+dsol
         zcf2=zc+dsol/2.
         bsol=bcret*spline(ncord,zcf2)
comment        write(70,110) i,zcf2,bsol
         dscm=dsol*100.
         tref=tref+xlsol/(vref*fnpas)
         call fldsol(bsol,dscm)
         write(16,520) i,zc,zcf,bsol
520      format(2x,'**STEP: ',i2,' LIMITS: inf(m)= ',f7.5,
     *        ' sup(m)= ',f7.5,' AVERAGE FIELD(kG): ',e12.5)
         bisol=bisol+bsol*dsol
         bisol2=bisol2+bsol*bsol*dsol
c     save equivalent transport matrix
         do ia=1,6
           do ib=1,6
             rs(ia,ib)=rcul(ia,ib)
           enddo
         enddo
         call mfordre(rcul,r,rs)
c   ************TEST******************
c   print first order matrix
comment         write(16,*) ' EQUIVALENT FIRST order MATRIX TRANSFORM (m-rad)'
comment       skl=0.5*acos(2.*r(1,1)-1.)*57.29578
comment       write(16,*) '**** K*LENGTH: ',skl,' degrees'
comment       DO IA=1,6
comment        write(16,100) (r(ia,ib),ib=1,6)
comment       ENDDO
comment        write(16,*) '*******cumul***************'
comment         DO IA=1,6
comment         write(16,100) (rcul(ia,ib),ib=1,6)
comment         ENDDO
c      **********END TEST********************
c  SPACE CHARGE
         if(.not.iesp) then
           iesp=.true.
           go to 720
         endif
         if(ichaes.and.xlsol.gt.0.) then
           if(sce10.eq.1.or.sce10.eq.3.) then
             dscm2=dscm*2.
             call cesp(dscm2)
             iesp=.false.
c     dispersion dE/E with respect to the C.O.G of the bunch
             call disp
cet*2015-04-25*test*s
cet             call compress(360.)
cet*2015-04-25*test*e
           endif
         endif
720      continue
C  envelope sol
         write(16,*) 'distance',davti,dsol
         davti=davti+dsol*100.
         call stapl((davti)*10.)
         zc=zcf
       enddo
       write(16,922) zc,bisol,bisol2,bisol2/zc
922    format(/,'Field length                 (m): ',e12.5,/,
     *  'Field integral            (kG.m): ',e12.5,/,
     *  'Field squared integral (kG**2.m): ',e12.5,/,
     *  'Field squared integral/L (kG**2): ',e12.5,/)
c110    format(2x,i5,2(2x,e12.5))
c   print first order matrix
       write(16,*) ' EQUIVALENT FIRST order MATRIX TRANSFORM (m-rad)'
       skl=0.5*acos(2.*rcul(1,1)-1.)*57.29578
       write(16,*) ' ******* K*LENGTH: ',skl,' degrees'
       DO IA=1,6
        write(16,100) (rcul(ia,ib),ib=1,6)
       ENDDO
100    format(6(3x,e12.5))
c **************************************************
c   evolution of the t.o.f of the reference (moved inside loop - Alt)
c       TREF=TREF+XLSOL/VREF
       dav1(idav,36)=ngood
C   plots
       CALL STAPL(davtot*10.)
       if(iemgrw) call emiprt(0)
       WRITE(16,*) 'particles lost in solenoid',ilost
       RETURN
       END
       SUBROUTINE fldsol(dbs,step)
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/femt/iemgrw,iemqesg
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/SOLE/L,KL,KO
       REAL*8 L,KL,KO
       l=step
       b=dbs
       gpai=0.
       do ii=1,ngood
        gpai=gpai+f(7,ii)/xmat
       enddo
       gpai=gpai/float(ngood)
       bpai=sqrt(1.-1./(gpai*gpai))
       xmco=xmat*bpai*gpai
       ri=33.356*xmco*1.e-01/qst
       ko = b /ri
       IF(KO .EQ. 0.) KO=1.E-16
       kl = ko*l
c half step
comment       kl=kl/2.
       call clear
       call elsol
       do ii=1,ngood
c   **********
comment        call clear
comment        gpai=f(7,ii)/xmat
comment        bpai=sqrt(1.-1./(gpai*gpai))
comment        xmco=xmat*bpai*gpai
comment        RI=33.356*XMCO*1.E-01/f(9,ii)
comment        KO = B /RI
comment        IF(KO .EQ. 0.) KO=1.E-16
comment        kl=ko*l
comment        call elsol
c   ***************
        call cobeam(ii,l)
C     evolution of the t.o.f has been made in the routine cobeam
comment         f(6,ii)=f(6,ii)+step/(bpai*vl)
       enddo
C test window after the step (only in the transverse directions)
       call cogetc
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       return
       end
       SUBROUTINE solnoid(imks,arg,xlsol)
c   ....................................................................
C             SOLENOID (fringe-fields are included in the transport matrix)
c             space charge computations at the middle of the lens
c   ....................................................................
C     IMKS: IFLAG (see ARG)
C     ARG:  IMKS = 0 then ARG is K (cm-1), otherwise ARG is the field BSOL (kG)
C     XLSOL : EFFECTIVE LENGHT (CM )
C
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
       COMMON/RIGID/BORO
       COMMON/DCSPA/IESP
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       common/rander/ialin
       common/femt/iemgrw,iemqesg
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       common/itvole/itvol,imamin
       common/tofev/ttvols
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       character*128 trace3h,trace3t,tif
       logical itvol,imamin
c       dimension rs(6,6),rcul(6,6)
       logical iesp,ichaes,shift,ialin,iemgrw
       character*1 cr
       COMMON/SOLE/L,KL,KO
       REAL*8 L,KL,KO
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
c print out on terminal of transport element # on one and the same line
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       write(16,*) ' ****** SOLENOID *********'
C      PLOT
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
       ilost=0
       FH0=FH/VL
       fcpi=fh*180./pi
c  random errors in alignment
       if(ialin) call randali
c    print out transport matrix (cog)
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       l=xlsol
       bpa=sqrt(1.-1./(gpa*gpa))
       xmco=xmat*bpa*gpa
       ri=33.356*xmco*1.e-01/qst
       if(imks.ne.0) then
c   ARG is the field B (kG)
        b=arg
c   KO = 2 * K (see ELSOL)
        KO = B /RI
        IF(KO .EQ. 0.) KO=1.E-16
       else
c   ARG is the strength K (cm-1)
c   KO = 2 * K (see ELSOL)
        ko = 2.*arg
        IF(KO .EQ. 0.) KO=1.E-16
        b=ko*ri
       endif
       KL = KO*L
       call clear
       call elsol
       xkql=(kl/2.)*57.29578
c trace3d
       kt3t=kt3t+1
       write(tif,6002) kt3t,kt3t,B*1000.,10.*XLSOL
6002   format('  nt(',i4,')= 5, a(1,',i4,')= ',f12.5,' , ',f9.5)
       trace3t(kt3t)=tif
c
       WRITE(16,101) XLSOL,B,ko/2.,ri,xkql
101    FORMAT('  LENGTH = ',F7.3,' CM ',/,
     X        '  FIELD = ',F10.4,' KG' ,/,
     *        '  K = ',e12.5,' cm-1',/,
     *        '  MOMENTUM = ',e12.5,' kG.cm',/,
     *        '  TRANSVERSE COORDINATES ROTATION = ',e12.5,' deg',/)
       call matrix
       write(16,10) ttvols*fcpi,davtot
10     FORMAT(' ** time of flight (input): ',e12.5,
     *        ' deg  position: ',e12.5,' cm')
c  start prints in file 'short.data'
       idav=idav+1
       iitem(idav)=5
       dav1(idav,1)=xlsol*10.
       dav1(idav,2)=b
       dav1(idav,3)=ko/2.
       davtot=davtot+xlsol
       dav1(idav,4)=davtot*10.
       dav1(idav,5)=ri
C     first half solenoid
       dg=xlsol
       l=xlsol/2.
       do ii=1,ngood
        call clear
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KO = B /RI
        IF(KO .EQ. 0.) KO=1.E-16
        KL = KO*L
        CALL ELSOL
        CALL cobeam(ii,l)
C     evolution of the t.o.f is made in routine cbeam
comment         f(6,ii)=f(6,ii)+l/(bpai*vl)
       enddo
C      space charge computations (only if l >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlsol)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
C test window after the first half solenoid
       call cogetc
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- window control
       tref=tref+xlsol/(vref*2.)
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c  second half solenoid
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KO = B /RI
        IF(KO .EQ. 0.) KO=1.E-16
        KL = KO*L
        CALL ELSOL
        CALL cobeam(ii,l)
       enddo
c   t.o.f
       TREF=TREF+XLSOL/(VREF*2.)
       call reject(ilost)
       if(itvol) ttvols=tref
       tcog=0.
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       if(itvol) then
        write(16,11) ttvols*fcpi,davtot,tref*fcpi,tcog*fcpi
11      FORMAT(' ** tof for adjustments: ',e12.5,
     *        ' deg at position: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg tof of the cog: ',
     * e12.5,' deg')
       else
        write(16,12) tref*fcpi,tcog*fcpi
12      FORMAT(' ** tof of the reference: ',e12.5,
     *        ' deg tof of the cog: ',e12.5,' deg')
       endif
       dav1(idav,36)=ngood
C   plots
       CALL STAPL(davtot*10.)
       if(iemgrw) call emiprt(0)
       WRITE(16,*) 'particles lost in solenoid',ilost
       RETURN
       END
       SUBROUTINE elsol
c  ..........................................................................
c    first order and second order R and T solenoid matrix
c  ..........................................................................
       implicit real*8 (a-h,o-z)
       COMMON/SOLE/L,KL,KO
       COMMON/RIGID/BORO
       REAL*8 L,KL,KO
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       COMMON/SECDR/ISEOR
       LOGICAL ISEOR
C      PURE SOLENOID FIELD, MATRIX R AND T  ARE IN  (M,RD)
C      SAVE KO AND L
       SAKO=KO
       SAL=L
C      convert (cm,rd) ==> (m,rd)
       KO=KO*100.
       L=L*1.E-02
       SN  = SIN(KL)
       CS  = COS(KL)
       R(4,4)  =  0.5 + 0.5*CS
       R(3,3)   = R(4,4)
       R(2,2)   = R(4,4)
       R(1,1)   = R(4,4)
       R(1,4)    =   (1.- CS)/KO
       R(3,2)    =  -R(1,4)
       R(4,1)    =  0.25*KO*(1.- CS)
       R(2,3)    =  -R(4,1)
       R(4,2)    =  -0.5*SN
       R(3,1)    =  R(4,2)
       R(2,4)    =  -R(3,1)
       R(1,3)    =  R(2,4)
       R(3,4)    =  R(1,3)*2./KO
       R(1,2)    =  R(3,4)
       R(4,3)    =  -0.25*KO*SN
       R(2,1)    =  R(4,3)
c     see element 16;  3.  in TRANSPORT USER MANUAL for SM (here: SM=0)
comment       R(5,6) = R(5,6)  +  L * SM**2/(RI**2 + SM**2)
      IF(.NOT.ISEOR) GO TO 200
      TEMP = 0.5*KO*L*SN
      T(1,1,6) = TEMP
      T(2,2,6) = TEMP
      T(3,3,6) = TEMP
      T(4,4,6) = TEMP
      T(1,2,6) = SN/KO - L*CS
      T(3,4,6) = T(1,2,6)
      TEMP = - 0.5*KO*L*CS
      T(1,3,6) = TEMP
      T(2,4,6) = TEMP
      T(4,2,6) = - TEMP
      T(3,1,6) = - TEMP
      T(1,4,6) = (1.0 - CS)/KO - L*SN
      T(3,2,6) = -T(1,4,6)
      T(2,1,6) = 0.25*KO*(KO*L*CS + SN)
      T(4,3,6) =T(2,1,6)
      T(2,3,6) = 0.25*KO*(1.0 - CS + KO*L*SN)
      T(4,1,6) = - T(2,3,6)
      T(5,2,2) =  0.5*L
      T(5,4,4) =  0.5*L
200   CONTINUE
C     RESTORE KO AND L
      KO=SAKO
      L=SAL
C
C     PATH LENGTH TERMS
C see element 16;  3.  in TRANSPORT USER MANUAL for SM here one takes: SM=0
comment  T(5,6,6) = T(5,6,6) - L*(SM**2 + 1.5*RI**2)*SM**2/
comment     1 (SM**2 + RI**2)**2
comment  T(5,5,6) = SM**2/(SM**2 + RI**2)
      RETURN
      END
      SUBROUTINE mfordre(rc,ra,rb)
c .................................................................
c  Calculates RC=RA*RB
c .................................................................
       implicit real*8 (a-h,o-z)
       DIMENSION RA(6,6), RB(6,6), RC(6,6)
       DO I1 = 1, 6
         DO I2 = 1, 6
           GHOST = 0.0
           DO I3 = 1, 6
             GHOST = GHOST + RA(I1,I3)*RB(I3,I2)
           ENDDO
           RC(I1,I2) = GHOST
         ENDDO
       ENDDO
       RETURN
       END
       SUBROUTINE solquad(iksq,args,argq,xlsol,rg)
c   ....................................................................
C        SOLENOID FIELD ASSOCIATED WITH QUADRUPOLE FIELD
c        space charge computations at the middle of the lens
c   ....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
       COMMON/RIGID/BORO
       COMMON/DCSPA/IESP
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       common/rander/ialin
       common/femt/iemgrw,iemqesg
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       logical iesp,shift,ialin,iemgrw
       character*1 cr
       COMMON/SLQ/L,KSO,KQO
       REAL*8 L,KSO,KQO
c --IKSQ: IFLAG
c --ARGS: STRENGTH or FIELD OF SOLENOID
c     If IKSQ = 0 then ARGS is the STRENGTH (cm-1), otherwise ARGS is the FIELD (kG)
c --ARGQ: STRENGTH or FIELD of QUADRUPOLE
c     If IKSQ = 0 then ARGQ is the STRENGTH (cm-2), otherwise ARGQ is the FIELD (kG)
c    SIGN CONVENTIONS:
c     SOLENOID: ARGS positive = rotate the transverse coordinates about the z-axis in the clockwise direction.
c     QUADRUPOLE: ARGQ positive = focusing in the plane (x,z)
C --XLSOL : EFFECTIVE LENGHT OF THE LENS(cm)
C --RG : APERTURE RADIUS OF THE LENS (cm)
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       FH0=FH/VL
C      PLOT
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
c  random errors in alignment
       if(ialin) call randali
       ilost=0
c   magnetic rigidity (cog)
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       l=xlsol
       bpa=sqrt(1.-1./(gpa*gpa))
       xmco=xmat*bpa*gpa
       ri=33.356*xmco*1.e-01/qst
       if(iksq.ne.0) then
c  argq: quadrupole field (kG), args:solenoid field:(kG)
c   quadrupole
        fq=argq
        bq=fq/rg
        kqo=bq/ri
        strq=kqo
c   solenoid (KSO = 2 * K, see ELSQ)
        bs=args
        kso=bs/ri
        strs=kso/2
       else
c  argq: quadupole strength (cm-2), args: solenoid strength (cm-1)
c    quadrupole
        kqo=argq
        strq=argq
        bq=argq*ri
        fq=bq*rg
c   solenoid (KSO = 2 * K, see ELSQ)
        strs=args
        bs=2.*strs*ri
        kso=bs/ri
       endif
       xks=kso/2.
       xksl=xks*l*57.29578
       write(16,*) '****SOLENOID+QUADRUPOLE*******'
       write(16,101) xlsol,rg,bs,strs,xksl,fq,strq,ri
101    FORMAT('  LENGTH =',F7.3,' cm  APERTURE RADIUS=',e12.5,' cm',/,
     X        '  SOLENOID: FIELD = ',F10.4,' kG  K = ',e12.5,
     *        ' cm-1  ROTATING ANGLE = ',e12.5,' deg',/,
     X        '  QUADRUPOLE: FIELD  =',F10.4,' kG   K2 =',e12.5,
     *        ' cm-2 ',/,
     *        '  RIGIDITY = ',e12.5,' kG.cm',/)
c    print transport matrix (c.o.g.)
C     For convenience, the matrix R and T are computed for a positive KQO = (B/r)*(1/BRO)
C     if KQO is negative one has seting up 90 degree rotation on the beam
       CALL elsq
       write (16,*)
     *' The matrix R and T are shown for a positive strength'
       write (16,*)
     *' For a negative strength set up 90 deg rotation on the beam'
       call matrix
c  start prints in file 'short.data'
       idav=idav+1
       iitem(idav)=11
       dav1(idav,1)=xlsol*10.
       dav1(idav,2)=bs
       dav1(idav,3)=fq
       dav1(idav,5)=strq
       dav1(idav,6)=kso/2.
       dav1(idav,7)=rg*10.
       dav1(idav,8)=ri
       davtot=davtot+xlsol
       dav1(idav,4)=davtot*10.
C     first half lens
       dg=xlsol
       l=xlsol/2.
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
c     KSO =K/2, see ELSQ
        KSO = BS / RI
        IF(KSO .EQ. 6.*0) KSO=1.E-16
        KQO = BQ / RI
        CALL ELSQ
C     the matrix R and T are computed for a positive KQO = (B/rg)*(1/BRO)
C     if KQO is negative set up 90 degree rotation on the beam
        IF(KQO.LT.6.*0) THEN
         irot=ii
         CALL rotat(irot)
         CALL cobeam(ii,l)
C   return the coordinates to their initial orientation
         irot=-ii
         CALL rotat(irot)
        ELSE
         CALL cobeam(ii,l)
       ENDIF
       enddo
C      Charge space effect (if l >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlsol)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
C test window after the first half solquad
       call cogetc
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- convert wdisp in dp/p (for window control)
       tref=tref+xlsol/(2.*vref)
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c   beam after the second half lens
       do ii=1,ngood
        CALL CLEAR
        gpai=f(7,ii)/xmat
        bpai=sqrt(1.-1./(gpai*gpai))
        xmco=xmat*bpai*gpai
        RI=33.356*XMCO*1.E-01/f(9,ii)
        KSO = BS / RI
c       KSO =K/2, see ELSQ
        IF(KSO .EQ. 0.) KSO=1.E-16
        KQO = BQ / RI
        CALL ELSQ
C     the matrix R and T are computed for a positive KQO = (B/r)*(1/BRO)
C     if KQO is negative set up 90 degree rotation on the beam
        IF(KQO.LT.6.*0) THEN
         IROT=II
         CALL rotat(irot)
         CALL cobeam(ii,l)
         IROT=-II
         CALL rotat(irot)
        ELSE
         CALL cobeam(ii,l)
        ENDIF
       enddo
c    t.o.f of reference
       TREF=TREF+XLSOL/(VREF*2.)
       call reject(ilost)
       dav1(idav,36)=ngood
C   plots
       CALL STAPL(davtot*10.)
       if(iemgrw) call emiprt(0)
       WRITE(16,*) 'particles lost in solenoid',ilost
       RETURN
       END
       SUBROUTINE rotat(ii)
c    .........................................................................
C     (+-) 90 DEG BEAM ROTATION
c    .........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       DIMENSION  FF(6)
       IF(II.GE.0) THEN
C   90 degree rotation
         FF(1)=F(2,II)
         FF(2)=F(3,II)
         FF(3)=F(4,II)
         FF(4)=F(5,II)
         F(2,II)=FF(3)
         F(3,II)=FF(4)
         F(4,II)=-FF(1)
         F(5,II)=-FF(2)
       ELSE
C   -90 degree rotation
         II=-II
         FF(1)=F(2,II)
         FF(2)=F(3,II)
         FF(3)=F(4,II)
         FF(4)=F(5,II)
         F(2,II)=-FF(3)
         F(3,II)=-FF(4)
         F(4,II)=FF(1)
         F(5,II)=FF(2)
      ENDIF
       RETURN
       END
       SUBROUTINE elsq
C...............................................................................
C            FIRST AND SECOND order MATRIX R AND T  (m,rd)
C     quadrupole field combined with a solenoid field
C...............................................................................
C **** WARNING : IN THIS ROUTINE ALL VARIABLES ARE REAL UNLESS
C                OTHERWIZE DECLARED *********
      IMPLICIT REAL*8(A-Z)
      COMMON/RIGID/BORO
      COMMON/SLQ/L,KSO,KQO
      COMMON /BLOC11/ R(6,6), T(6,6,6)
      COMMON/SECDR/ISEOR
      LOGICAL ISEOR
      AL=L*1.e-02
      AKQ= KQO*1.E04
      AKS= KSO*1.e02
      AKS=AKS/2.
      AKQ4=AKQ*AKQ
      AKS2=AKS*AKS
      AKS3=AKS2*AKS
      AKS4=AKS2*AKS2
      Q2=DSQRT(AKQ4+4.0*AKS4)
      Q2I=1.0/Q2
      SK1=SQRT(2.0*AKS2+Q2)
      SK3=SQRT(ABS(Q2-2.0*AKS2))
      SS=SIN(SK1*AL)
      SC=COS(SK1*AL)
      BS=SINH(SK3*AL)
      BC=COSH(SK3*AL)
      SKP=0.5*(SK1+SK3)
      SKM=0.5*(SK1-SK3)
      SKP2=SKP*SKP
      SKP3=SKP2*SKP
      SKM2=SKM*SKM
      SKM3=SKM2*SKM
      R(1,1)=Q2I*(SKP2*SC+SKM2*BC)
      R(1,2)=Q2I*(SKP*SS-SKM*BS)
      R(1,3)=Q2I*AKS*(SKM*SS+SKP*BS)
      R(1,4)=Q2I*AKS*(BC-SC)
      R(2,1)=Q2I*(-SKP3*SS-SKM3*BS)
      R(2,2)=Q2I*(SKP2*SC+SKM2*BC)
      R(2,3)=Q2I*AKS3*(SC-BC)
      R(2,4)=Q2I*AKS*(SKP*SS-SKM*BS)
      R(3,1)=Q2I*AKS*(SKM*BS-SKP*SS)
      R(3,2)=Q2I*AKS*(SC-BC)
      R(3,3)=Q2I*(SKM2*SC+SKP2*BC)
      R(3,4)=Q2I*(SKM*SS+SKP*BS)
      R(4,1)=Q2I*AKS3*(BC-SC)
      R(4,2)=Q2I*AKS*(-SKM*SS-SKP*BS)
      R(4,3)=Q2I*(SKP3*BS-SKM3*SS)
      R(4,4)=Q2I*(SKM2*SC+SKP2*BC)
      IF(.NOT.ISEOR) RETURN
      DAKS=-AKS
      DAKQ=-AKQ
      DAKQ4=2.0*AKQ*DAKQ
      DAKS2=2.0*AKS*DAKS
      DAKS3=3.0*AKS2*DAKS
      DAKS4=4.0*AKS3*DAKS
      DQ2=0.5*Q2I*(DAKQ4+4.0*DAKS4)
      Q2I2=Q2I*Q2I
      DQ2I=-Q2I2*DQ2
      DSK1=(2.0*DAKS2+DQ2)/(2.0*SK1)
      DSK3=0.0
      IF(SK3.NE.6.*0)DSK3=(DQ2-2.0*DAKS2)/(2.0*SK3)
      DSS=SC*DSK1*AL
      DSC=-SS*DSK1*AL
      DBS=BC*DSK3*AL
      DBC=BS*DSK3*AL
      DSKP=0.5*(DSK1+DSK3)
      DSKM=0.5*(DSK1-DSK3)
      DSKP2=2.0*SKP*DSKP
      DSKP3=3.0*SKP2*DSKP
      DSKM2=2.0*SKM*DSKM
      DSKM3=3.0*SKM2*DSKM
      T(1,1,6)=DQ2I*(SKP2*SC+SKM2*BC)+
     > Q2I*(DSKP2*SC+SKP2*DSC+DSKM2*BC+SKM2*DBC)
      T(1,2,6)=DQ2I*(SKP*SS-SKM*BS)+
     > Q2I*(DSKP*SS+SKP*DSS-DSKM*BS-SKM*DBS)
      T(1,3,6)=(DQ2I*AKS+Q2I*DAKS)*(SKM*SS+SKP*BS)
     > +Q2I*AKS*(DSKM*SS+SKM*DSS+DSKP*BS+SKP*DBS)
      T(1,4,6)=(DQ2I*AKS+Q2I*DAKS)*(BC-SC)
     > +Q2I*AKS*(DBC-DSC)
      T(2,1,6)=DQ2I*(-SKP3*SS-SKM3*BS)
     > +Q2I*(-DSKP3*SS-SKP3*DSS-DSKM3*BS-SKM3*DBS)
      T(2,2,6)=T(1,1,6)
      T(2,3,6)=(DQ2I*AKS3+Q2I*DAKS3)*(SC-BC)
     > +Q2I*AKS3*(DSC-DBC)
      T(2,4,6)=(DQ2I*AKS+Q2I*DAKS)*(SKP*SS-SKM*BS)
     > +Q2I*AKS*(DSKP*SS+SKP*DSS-DSKM*BS-SKM*DBS)
      T(3,1,6)=-T(2,4,6)
      T(3,2,6)=-T(1,4,6)
      T(3,3,6)=DQ2I*(SKM2*SC+SKP2*BC)
     > +Q2I*(DSKM2*SC+SKM2*DSC+DSKP2*BC+SKP2*DBC)
      T(3,4,6)=DQ2I*(SKM*SS+SKP*BS)
     > +Q2I*(DSKM*SS+SKM*DSS+DSKP*BS+SKP*DBS)
      T(4,1,6)=-T(2,3,6)
comment      AMAT(17,4,MATADR)=-AMAT(21,1,MATADR)
      T(4,2,6)=-T(1,3,6)
      T(4,3,6)=DQ2I*(SKP3*BS-SKM3*SS)
     > +Q2I*(DSKP3*BS+SKP3*DBS-DSKM3*SS-SKM3*DSS)
      T(4,4,6)=T(3,3,6)
      AISSSC=0.0
      IF(SK1.NE.0.0)AISSSC=SS*SS/(2.0*SK1)
      AIBSBC=0.0
      IF(SK3.NE.0.0)AIBSBC=BS*BS/(2.0*SK3)
      AISSBC=Q2I*(SK3*SS*BS-SK1*SC*BC+SK1)
      AISCBS=Q2I*(SK3*SC*BC+SK1*SS*BS-SK3)
      AISS2=0.0
      AISC2=AL
      AIBS2=0.0
      AIBC2=AL
      IF(SK1.NE.0.0)AISS2=0.5*(AL-SS*SC/SK1)
      IF(SK1.NE.0.0)AISC2=0.5*(AL+SS*SC/SK1)
      IF(SK3.NE.0.0)AIBS2=0.5*(BS*BC/SK3 - AL)
      IF(SK3.NE.0.0)AIBC2=0.5*(AL+BS*BC/SK3)
      AISSBS=Q2I*(SK3*SS*BC-SK1*SC*BS)
      AISCBC=Q2I*(SK1*SS*BC+SK3*SC*BS)
      AKS5=AKS4*AKS
      AKS6=AKS5*AKS
      SKP4=SKP3*SKP
      SKP5=SKP4*SKP
      SKP6=SKP5*SKP
      SKM4=SKM3*SKM
      SKM5=SKM4*SKM
      SKM6=SKM5*SKM
      Q4I=Q2I*Q2I
      T(5,1,1)=Q4I*0.5*(SKP6*AISS2+SKM6*AIBS2-2.0*SKP3*SKM3
     >   *AISSBS+AKS6*(AISC2+AIBC2-2.0*AISCBC))
      T(5,1,2)=Q4I*(-SKP5*AISSSC-SKP3*SKM2*AISSBC
     >      -SKM3*SKP2*AISCBS
     >      -SKM5*AIBSBC-AKS4*(-SKM*AISSSC-SKP*AISCBS+
     >      SKM*AISSBC+SKP*AIBSBC))
      T(5,1,3)=Q4I*AKS3*((SKP3-SKM3)*(AISSBC-AISSSC)
     >  +(SKP3+SKM3)*(AIBSBC-AISCBS))
      T(5,1,4)=Q4I*(AKS*(-SKP4*AISS2-(SKP3*SKM+SKM3*SKP)*AISSBS
     >  +SKM4*AIBS2)+AKS3*(-SKM2*AISC2-(SKP2-SKM2)*AISCBC+SKP2*AIBC2))
      T(5,2,2)=Q4I*0.5*(SKP4*AISC2+2.0*SKP2*SKM2*AISCBC
     >  +SKM4*AIBC2+AKS2*(SKM2*AISS2+2.0*SKP*SKM*AISSBS+SKP2*AIBS2))
      T(5,2,3)=Q4I*(AKS3*(SKP2*AISC2-(SKP2-SKM2)*AISCBC
     >  -SKM2*AIBC2)-AKS*(-SKM4*AISS2+(SKM*SKP3-SKP*SKM3)*AISSBS+
     >  SKP4*AIBS2))
      T(5,2,4)=Q4I*AKS*((SKP3-SKM3)*AISSSC-(SKP2*SKM+SKP*SKM2)
     >  *AISCBS+(SKM2*SKP-SKM*SKP2)*AISSBC-(SKM3+SKP3)*AIBSBC)
      T(5,3,3)=Q4I*0.5*(AKS6*(AISC2-2.0*AISCBC+AIBC2)
     >  +SKM6*AISS2-2.0*SKM3*SKP3*AISSBS+SKP6*AIBS2)
      T(5,3,4)=Q4I*(AKS4*(SKP*AISSSC-SKM*AISCBS-SKP*AISSBC
     >  +SKM*AIBSBC)-SKM5*AISSSC-SKM3*SKP2*AISSBC+SKP3*SKM2*AISCBS
     >  +SKP5*AIBSBC)
      T(5,4,4)=Q4I*0.5*(AKS2*(SKP2*AISS2-2.0*SKP*SKM*AISSBS
     > +SKM2*AIBS2)+SKM4*AISC2+2.0*SKM2*SKP2*AISCBC+SKP4*AIBC2)
      RETURN
      END
      SUBROUTINE fdrift(xl,npart,imit)
c............................................................................
C          DIVIDE A DRIFT LENGHT OF :XL CM  IN :NPART PARTIAL DRIFTS
c          This will allow several space charge computations in the drift
c............................................................................
       implicit real*8 (a-h,o-z)
       dl=xl/float(npart)
       do i=1,npart
         call drift(dl)
         if(imit.ne.0.) call emiprt(0)
       enddo
       return
       end
              SUBROUTINE drift(dl)
c   ............................................................................
C          DRIFT LENGHT
c         space charge computation at the middle of the drift
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/RIGID/BORO
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/tapes/in,ifile,meta
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DCSPA/IESP
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/SHIF/DTIPH,SHIFT
       common/femt/iemgrw,iemqesg
       common/posc/xpsc
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       common/tofev/ttvols
       COMMON/ITVOLE/ITVOL,IMAMIN
       logical iesp,ichaes,shift,iemgrw,itvol,imamin
       ilost=0
       FH0=FH/VL
       fcpi=fh*180./pi
       write(16,*)'*** DRIFT of length ',dl,' cm'
       if(itvol)write(16,10) ttvols*fcpi,davtot
10     FORMAT(' ** tof for adjustments at input: ',e12.5,
     *        ' deg at position: ',e12.5,' cm in the lattice')
C      STATISTICS
       if(iprf.eq.1) call stapl(davtot*10.)
c       random error in position only if dl is positif
C        Beam at the half drift
       dg2=dl/2.
       do i=1,ngood
c   conversion mrad ===> rad
         F2=F(3,I)*.001
         F4=F(5,I)*.001
         F(2,I)=F(2,I)+DG2*TAN(F2)
         F(4,I)=F(4,I)+DG2*TAN(F4) / COS(F2)
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         vpai=vl*bpai
cold         f(6,i)=f(6,i)+dg2/vpai
         f(6,i)=f(6,i)+dg2/(vpai*cos(f2)*cos(f4))
       enddo
C      ENVEL
c  start prints in file 'short.data'
       davtot=davtot+dl
       if(dl.gt.0.) then
c         nlength=nlength+1
         idav=idav+1
         iitem(idav)=7
         dav1(idav,1)=dl*10.
         dav1(idav,4)=davtot*10.
       endif
C      Charge space  (if dl >1.e-04)
       if(ichaes.and.dl.ge.1.e-04) then
c 11.03.2010         if(sce10.eq.2..or.sce10.eq.3.) then
         if(sce10.eq.3.) then
           iesp=.true.
           write(16,*) ' space charge at the middle'
           call cesp(dl)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
c Test window after the first half drift
       tref=tref+dl/(vref*2.)
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       il=ilost
C        Beam for the second half drift
       do i=1,ngood
         f2=f(3,i)*.001
         f4=f(5,i)*.001
         f(2,i)=f(2,i)+dg2*tan(f2)
         f(4,i)=f(4,i)+dg2*tan(f4) / cos(f2)
         gpai=f(7,i)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
         vpai=vl*bpai
cold         f(6,i)=f(6,i)+dg2/vpai
         f(6,i)=f(6,i)+dg2/(vpai*cos(f2)*cos(f4))
       enddo
c Test window after the second half drift (only in transverse directions and phase)
c  TEST *****
       tref=tref+dl/(vref*2.)
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
cold         call shuffle
       il=il+ilost
C     change the reference and the TOF
       if(itvol) ttvols=tref
       if(dl.gt.0.)then
         dav1(idav,36)=ngood
C  envelope
         call stapl(davtot*10.)
       endif
       tcog=0.
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       if(itvol) then
        write(16,11) ttvols*fcpi,davtot,tref*fcpi,tcog*fcpi
11      FORMAT(' ** tof for adjustments : ',e12.5,
     *        ' deg at position: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg tof of the cog: ',
     * e12.5,' deg')
       else
        write(16,12) tref*fcpi,tcog*fcpi
12      FORMAT(' ** tof of the reference: ',e12.5,
     *        ' deg tof of the cog: ',e12.5,' deg')
       endif
       write(16,*) 'particles lost in drift: ',il
       if(iemgrw.and.dl.gt.0.) then
         if(iemqesg.eq.2) call emiprt(0)
       endif
       RETURN
       END
       SUBROUTINE clear
C    CLEAR MATRIX R AND T
       implicit real*8 (a-h,o-z)
      COMMON /BLOC11/ R(6,6), T(6,6,6)
      COMMON/SECDR/ISEOR
      LOGICAL ISEOR
C   CLEAR R
       DO IA=1,6
        DO IB=1,6
         R(IA,IB)=0.
         IF(IA.EQ.IB) R(IA,IB)=1.
        ENDDO
       ENDDO
C   CLEAR T
      IF(ISEOR) THEN
       DO IA=1,6
        DO IB=1,6
         DO IC=1,6
          T(IA,IB,IC)=0.
         ENDDO
        ENDDO
       ENDDO
      ENDIF
       RETURN
       END
       SUBROUTINE cobeam(ii,xl)
C       BEAM COMPUTATION
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RADIA/TRT,RMOY,XINTF,CRAE
       COMMON/SECDR/ISEOR
       LOGICAL ISEOR
       DIMENSION SF(6),SSF(6)
C     MATRIX R AND T ARE IN  M-RD
C   beam conversion:  CM ===> M      MRD ==> RD
       SF(1)=F(2,II)*1.E-02
       SF(2)=F(3,II)*1.E-03
       SF(3)=F(4,II)*1.E-02
       SF(4)=F(5,II)*1.E-03
       SF(5)=XL*1.E-02
       SF(5)=SF(5)/(cos(sf(2))*cos(sf(4)))
       SF(6)=(FD(II)-1.)
C       FIRST order COMPUTATION (M-RD)
       DO IA=1,6
        SSF(IA)=0.
        DO IB=1,6
         SSF(IA)=R(IA,IB)*SF(IB)+SSF(IA)
        ENDDO
       ENDDO
       IF(.NOT.ISEOR) GO TO 10
C    SECOND order COMPUTATION  (M-RD)
       DO IA=1,6
        DO IB=1,6
         DO IC=1,6
          SSF(IA)=SSF(IA)+T(IA,IB,IC)*SF(IB)*SF(IC)
         ENDDO
        ENDDO
       ENDDO
10     CONTINUE
       DO IA=1,4
         F(IA+1,II)=SSF(IA)
       ENDDO
C   CONVERT RD ==> MRD    M ==> CM
       F(3,II)=F(3,II)*1000.
       F(5,II)=F(5,II)*1000.
       F(2,II)=F(2,II)*100.
       F(4,II)=F(4,II)*100.
       GPAI=F(7,II)/XMAT
       BPAI=SQRT(1.-1./(GPAI*GPAI))
       VPAI=VL*BPAI
       TRT=100.*SSF(5)/VPAI
       F(6,II)=F(6,II)+TRT
      RETURN
      END
      SUBROUTINE disp
c     Compute the dispersion dE/E with respect to the center of gravity of the bunch
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/mcs/imcs,ncstat,cstat(20)
       common /consta/ vl,pi,xmat,rpel,qst
       gcog=0.
       do i=1,ngood
         gcog=f(7,i)/xmat+gcog
       enddo
       gcog=gcog/float(ngood)
       bcog=sqrt(1.-1./(gcog*gcog))
       do i=1,ngood
         gpai=f(7,i)/xmat
         if(gpai.ge.1.) then
           bpai=sqrt(1.-1./(gpai*gpai))
         else
           bpai=0.
         endif
         fd(i)=bpai/bcog * gpai/gcog
       enddo
       return
       end
cold    do ist=1,ncstat
cold         gcog=0.
cold         nii=0
cold         do ii=1,ngood
cold          if(f(9,ii).eq.cstat(ist)) then
cold            gcog=gcog+f(7,ii)/xmat
cold            nii=nii+1
cold           endif
cold          enddo
cold           gcog=gcog/float(nii)
cold           bcog=sqrt(1.-1./(gcog*gcog))
cold          DO II=1,NGOOD
cold           if(f(9,ii).eq.cstat(ist)) then
cold            gpai=f(7,ii)/xmat
cold            bpai=sqrt(1.-1./(gpai*gpai))
cold            fd(ii)=(gpai*bpai)/(gcog*bcog)
cold           endif
cold          enddo
cold         enddo
cold        return
cold        end
      SUBROUTINE zrotat(zrot)
c    .........................................................................
C     BEAM ROTATION
c     the transverse coordinates X and Y may be rotated
c     through an angle about the axis tangent to the
c     central trajectory at the point in question.
c    .........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/grot/rzot,izrot
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       logical izrot
       DIMENSION  RS(6,6),FF(6),FC(6)
       WRITE(16,100) ZROT
 100   FORMAT(/,20X,'BEAM ROTATION ',F10.4
     1 ,'  degrees  ABOUT THE POSITIVE Z-AXIS',/)
       if(.not.izrot) then
        izrot=.true.
        rzot=zrot
        go to 500
       endif
       if(izrot) izrot=.false.
500    continue
c  Daves start
       idav=idav+1
       iitem(idav)=20
       dav1(idav,1)=zrot
c*et*2015-07-07*s
c initialize
       DO IA=1,4
         DO IB=1,4
           RS(IA,IB)=0.
         ENDDO
       ENDDO
c*et*2015-07-07*e
C      Conversion DEG ==> RAD
       ZROT=  ZROT*pi/180.
       RS44  =  COS(ZROT)
       RS( 4,4 )  =    RS44
       RS( 3,3 )  =    RS44
       RS( 2,2 )  =    RS44
       RS( 1,1 )  =    RS44
       RS24     =  SIN(ZROT)
       RS( 2,4 )  =    RS24
       RS( 1,3 )  =    RS24
       RS( 4,2 )  =  - RS24
       RS( 3,1 )  =  - RS24
       DO II=1,ngood
         FF(1)=F(2,II)
         FF(2)=F(3,II)
         FF(3)=F(4,II)
         FF(4)=F(5,II)
         DO IA=1,4
           FC(IA)=0.
           DO IB=1,4
             FC(IA)=FC(IA)+FF(IB)*RS(IA,IB)
           ENDDO
         ENDDO
         F(2,II)=FC(1)
         F(3,II)=FC(2)
         F(4,II)=FC(3)
         F(5,II)=FC(4)
       ENDDO
c110    continue
       RETURN
       END
      SUBROUTINE zrotap(zrot)
c    .........................................................................
C     BEAM ROTATION
c     the transverse coordinates X and Y may be rotated
c     through an angle about the axis tangent to the
c     central trajectory at the point in question.
c    .........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       common/grot/rzot,izrot
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       logical izrot
       DIMENSION  RS(6,6),FF(6),FC(6)
       WRITE(16,100) ZROT
 100   FORMAT(/,20X,'BEAM ROTATION ',F10.4
     1 ,'  degrees  ABOUT THE POSITIVE Z-AXIS',/)
       rzot=zrot
C      Conversion DEG ==> RAD
       ZROT=  ZROT*pi/180.
       RS44  =  COS(ZROT)
       RS( 4,4 )  =    RS44
       RS( 3,3 )  =    RS44
       RS( 2,2 )  =    RS44
       RS( 1,1 )  =    RS44
       RS24     =  SIN(ZROT)
       RS( 2,4 )  =    RS24
       RS( 1,3 )  =    RS24
       RS( 4,2 )  =  - RS24
       RS( 3,1 )  =  - RS24
       DO II=1,ngood
         FF(1)=F(2,II)
         FF(2)=F(3,II)
         FF(3)=F(4,II)
         FF(4)=F(5,II)
         DO IA=1,4
           FC(IA)=0.
           DO IB=1,4
             FC(IA)=FC(IA)+FF(IB)*RS(IA,IB)
           ENDDO
         ENDDO
         F(2,II)=FC(1)
         F(3,II)=FC(2)
         F(4,II)=FC(3)
         F(5,II)=FC(4)
       ENDDO
       RETURN
       END
       SUBROUTINE aliner
       implicit real*8 (a-h,o-z)
C     ..................................
C     ALIGNMENT  errors:
C       HORIZONTAL : XL(cm) XLP(mrad)
C       VERTICAL   : YL(cm) YLP(mrad)
C     ..................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/ALIN/XL,YL,XPL,YPL
       common/faisc/f(10,iptsz),imax,ngood
       WRITE(16,100) XL,YL,XPL,YPL
100    FORMAT(/,5X,' KICK x(cm)    y(cm)   : ',2(e12.5,2x),/,
     *         5X,' KICK xp(mrad) yp(mrad): ',2(e12.5,2x),//)
C    CALCUL DES COORDONNEES DES TRAJECTOIRES
C       UNITES  =  CM-MRD
       DO II=1,ngood
         F(2,II)=F(2,II) + XL
         F(4,II)=F(4,II) + YL
         F(3,II)=F(3,II) + XPL
         F(5,II)=F(5,II) + YPL
       ENDDO
       RETURN
       END
       SUBROUTINE randali
       implicit real*8 (a-h,o-z)
C     ..................................
C      RANDOM ALIGNMENT  errors :
C       HORIZONTAL PLANE : XL(cm) XLP(mrad)
C       VERTICAL PLANE   : YL(cm) YLP(mrad)
C     ..................................
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/ALIN/XL,YL,XPL,YPL
       common/faisc/f(10,iptsz),imax,ngood
       common/apel/iapel
       dimension trans(1)
       write(16,*) ' random error in alignment with:'
       WRITE(16,100) XL,YL,XPL,YPL
100    FORMAT(/,5X,' KICK x(cm)    y(cm)   : ',2(e12.5,2x),/,
     *         5X,' KICK xp(mrad) yp(mrad): ',2(e12.5,2x),//)
       len=1
       rdcf=.5
       call rlux(trans,len)
       if(trans(1).le.rdcf) sign=-1.
       if(trans(1).gt.rdcf) sign=1.
       call rlux(trans,len)
       xla=xl*sign*trans(1)
       call rlux(trans,len)
       if(trans(1).le.rdcf) sign=-1.
       if(trans(1).gt.rdcf) sign=1.
       call rlux(trans,len)
       yla=yl*sign*trans(1)
       call rlux(trans,len)
       if(trans(1).le.rdcf) sign=-1.
       if(trans(1).gt.rdcf) sign=1.
       call rlux(trans,len)
       xpla=xpl*sign*trans(1)
       call rlux(trans,len)
       if(trans(1).le.rdcf) sign=-1.
       if(trans(1).gt.rdcf) sign=1.
       call rlux(trans,len)
       ypla=ypl*sign*trans(1)
       WRITE(16,100) XLA,YLA,XPLA,YPLA
       DO II=1,ngood
         f(2,ii)=f(2,ii)+xla
         f(4,ii)=f(4,ii)+yla
         f(3,ii)=f(3,ii)+xpla
         f(5,ii)=f(5,ii)+ypla
       enddo
       return
       end
       SUBROUTINE matrix
C    print first and second order matrix of a lens
       implicit real*8 (a-h,o-z)
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       COMMON/SECDR/ISEOR
       LOGICAL ISEOR
       write(16,*) ' TRANSPORT MATRIX (m-rd)'
       write(16,*) ' FIRST ORDER TRANSPORT********'
       DO IA=1,6
        write(16,100) (r(ia,ib),ib=1,6)
       ENDDO
100    format(6(3x,e12.5))
      write(16,*) ' *************************************************'
      IF(ISEOR) THEN
       write(16,*) ' SECOND ORDER TRANSPORT (m-rd)********'
       DO IA=1,6
        DO IB=1,6
         DO IC=1,6
          IF(T(IA,IB,IC).NE.6.*0) WRITE(16,101)IA,IB,IC,T(IA,IB,IC)
         ENDDO
        ENDDO
       ENDDO
      write(16,*) ' *************************************************'
      ENDIF
101   format(' T',3(i1),3x,e12.5)
      RETURN
      END
       SUBROUTINE egun(fmult,indp)
c **********************************************************************************
c    method:Bode's rule
c  read on the disk the axial field of the DC gun field
c     z (m)   E(z) is normalized
c **********************************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/ vl,pi, xmat,rpel,qst
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/DYN/TREF,VREF
       COMMON/RIGID/BORO
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/femt/iemgrw,iemqesg
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/dcspa/iesp
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       dimension gam(3000),xe(3000),xpe(3000),ye(3000),ype(3000)
       character*1 cr
       logical flgsc,iesp,ichaes,iemgrw
c    print out on terminal of transport element # on one and the same line
       nrres=nrres+1
       cr=char(13)
       WRITE(6,8254) NRTRE,NRRES
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5)
       write(6,*) 'EGUN calculation started'
c  energy at the entrance
       we=0.
       do i=1,ngood
        we=we+f(7,i)
       enddo
       we=we/float(ngood)-xmat
c  indp: define the number of space charge computations
c   indp = 1 : 8 space charge computations. The EGUN field is divided in 16 elements
c   indp = 2 : 16 space charge computation.                           in 32 elements
c   indp = 3 : 32 space charge computation.                           in 64 elements
       if(indp.gt.3) indp=3
       if(indp.eq.1) ipart=16
       if(indp.ge.2) ipart=32
       if(indp.eq.3) ipart=64
c research the number of steps w.r.t. the input energy
c    convert energy: MeV --> eV
       we=we*1.e06
       weinf=19.99
       if(we.le.weinf) then
        write(6,*) ' Energy at the cathode: ',we,
     *             ' eV is below the lower limit of 20 eV '
        stop
       endif
c    read egun field on the disk
       e0=xmat
       vlm=vl*1.e-02
       read(22,*)npt
       do i=1,npt
        read(22,*) xspl(i),yspl(i)
       enddo
       zinf=xspl(1)
       zsup=xspl(npt)
       egl=zsup-zinf
       call deriv2(npt)
       elgun=egl*100.
C      PLOT
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
c  start prints in file 'short.data'
       idav=idav+1
       idavs=idav
       iitem(idav)=13
       dav1(idav,1)=egl*1000.
       dav1(idav,2)=fmult
       davtot=davtot+elgun
       dav1(idav,4)=elgun*10.
c    divide the length in ipart partitions
c    the space charge computations are made at the middle of each partition
       eglp=egl/float(ipart)
c      eglsc :space charge acting lenght (cm)  (of each part)
       scl=float(ipart)/2.
       eglsc=egl*100./scl
c    compute the EGUN field (MV)
       npas=200
       xpas=eglp/float(npas)
       xnhf=0
       dcfld=0.
       do i=1,npas
         fpos=xnhf*xpas
         tspl=spline(npt,fpos)*fmult
         dcfld=dcfld+qst*tspl*xpas
         xnhf=xnhf+1.
       enddo
c4443   format(2x,i5,3(2x,e12.5))
       dav1(idav,3)=dcfld*1000.
       write(16,101) elgun,fmult,dcfld*1000.
101    FORMAT(5X,' FIELD LENGTH =',F7.3,' CM ',/,
     X  5X,'  FIELD CREST=',F10.4,' MV/m' ,/
     x  5x,'  FIELD STENGTH= ',E12.5,' kV',/)
comment       write(49,5557)
comment5557   format(6x,'z(m)',14x,'x(m)',14x,'xp(rad)',11x,'y(m)',14x,
comment     *        'yp(rad)',11x,'energy(MeV)')
c **** demarrage sur 0.5 mm (?) soit 0.5 e-03 m
c     thresold of energy
       swe1=19.99
       swe2=79.99
       npas=200
       if(we.lt.swe1) npas=400
       if(we.lt.swe2) npas=300
       xlstart=0.5e-03
       xpas=xlstart/float(npas)
       npas1=npas+1
       iflg=-1
       xnht=0.
c*et*2012*Jun*s
       flgsc=.false.
c*et*2012*Jun*e
501    continue
       if(iflg.eq.ipart) go to 500
       do j=1,ngood
         xnh=xnht
         qc=f(9,j)
         gam0=f(7,j)/e0
         gam(1)=gam0
         tof=f(6,j)
c       the tranverse coordinates are converted in (m,rad)
         x0=f(2,j)*1.e-02
         y0=f(4,j)*1.e-02
         t0=f(3,j)*1.e-03
         p0=f(5,j)*1.e-03
c        Pitch transformation
         xe0=x0*(gam0*gam0-1.)**0.25
         xpe0=t0*(gam0*gam0-1.)**0.25
         ye0=y0*(gam0*gam0-1.)**0.25
         ype0=p0*(gam0*gam0-1.)**0.25
         fpos=xnh*xpas
         tspl=spline(npt,fpos)*fmult
         dgam=(qc/e0)*tspl
         dgami=dgam
         xpe0=xpe0+.5*xe0*gam0*dgam/(gam0*gam0-1)
         ype0=ype0+.5*ye0*gam0*dgam/(gam0*gam0-1)
         a1=qc*qc/(4.*e0*e0)
         a2=qc/(e0*vlm)
         xe(1)=xe0
         xpe(1)=xpe0
         ye(1)=ye0
         ype(1)=ype0
         do i=2,npas1
          i1=i-1
          fpos1=xnh*xpas
          fpos2=(xnh+0.25)*xpas
          fpos3=(xnh+0.5)*xpas
          fpos4=(xnh+0.75)*xpas
          fpos5=(xnh+1.)*xpas
          tspl1=spline(npt,fpos1)*fmult
          tspl2=spline(npt,fpos2)*fmult
          tspl3=spline(npt,fpos3)*fmult
          tspl4=spline(npt,fpos4)*fmult
          tspl5=spline(npt,fpos5)*fmult
          cw=(qc/e0)*xpas/90.
          tspl=7.*tspl1+32.*tspl2+12.*tspl3+32.*tspl4+7.*tspl5
          gam(i)=cw*tspl+gam(i1)
          gam1=gam(i1)
          gam5=gam(i)
          dgam1=(qc/e0)*tspl1
          cof1=(gam5-gam1)/(xpas*xpas)
          cof2=dgam1/xpas
          cof=cof1-cof2
          gam2=gam1+dgam1*xpas/4.+cof*xpas*xpas/16.
          gam3=gam1+dgam1*xpas*0.5+cof*xpas*xpas/4.
          gam4=gam1+dgam1*xpas*0.75+cof*9.*xpas*xpas/16.
          gams1=gam1*gam1
          gams2=gam2*gam2
          gams3=gam3*gam3
          gams4=gam4*gam4
          gams5=gam5*gam5
comment          bgt1=(gams1-1.)**1.5
          bgt2=(gams2-1.)**1.5
          bgt3=(gams3-1.)**1.5
          bgt4=(gams4-1.)**1.5
          bgt5=(gams5-1.)**1.5
          tslpt=8.*tspl2/bgt2+6.*tspl3/bgt3+24.*tspl4/bgt4
     *          +7.*tspl5/bgt5
          dt=a2*xpas*xpas*tslpt/90.
          bet=sqrt(1.-1./gams1)
          tof=tof+xpas/(vlm*bet)+dt
          f(7,j)=gam(i)*e0
          f(6,j)=tof
          bg1=(gams1+2.)/((gams1-1.)*(gams1-1.))
          bg2=(gams2+2.)/((gams2-1.)*(gams2-1.))
          bg3=(gams3+2.)/((gams3-1.)*(gams3-1.))
          bg4=(gams4+2.)/((gams4-1.)*(gams4-1.))
          bg5=(gams5+2.)/((gams5-1.)*(gams5-1.))
          bgts1=bg1*tspl1*tspl1
          bgts2=bg2*tspl2*tspl2
          bgts3=bg3*tspl3*tspl3
          bgts4=bg4*tspl4*tspl4
          bgts5=bg5*tspl5*tspl5
c            gtpm=bgts1+3.*bgts2+3.*bgts3+bgts4
          gtpm=7.*bgts1+32.*bgts2+12.*bgts3+32.*bgts4+7.*bgts5
c                      gtm=bgts1+2.*bgts2+bgts3
          gtm=7.*bgts1+24.*bgts2+6.*bgts3+8.*bgts4
c                 gtpz=bgts2+2.*bgts3+bgts4
          gtpz=8.*bgts2+6.*bgts3+24.*bgts4+7.*bgts5
          gtm1=2.*bgts2+3.*bgts3+18.*bgts4+7.*bgts5
          de=-a1*xpas*xpas*gtm/90.
          de1=-a1*xpas*xpas*xpas*gtm1/90.
          dpe1=-a1*xpas*gtpm/90.
          dpe2=-a1*xpas*xpas*gtpz/90.
          dxpe1=dpe1*xe(i1)
          dype1=dpe1*ye(i1)
          dxpe2=dpe2*xpe(i1)
          dype2=dpe2*ype(i1)
          dxe=de*xe(i1)+de1*xpe(i1)
          dye=de*ye(i1)+de1*ype(i1)
          xpe(i)=xpe(i1)+dxpe1+dxpe2
          ype(i)=ype(i1)+dype1+dype2
          xe(i)=xe(i1)+dxe+xpe(i1)*xpas
          ye(i)=ye(i1)+dye+ype(i1)*xpas
c       back to the real variables and convert to (cm,mrad)
          gamm1=(gams5-1.)**0.25
          gamm2=(gams5-1.)**1.25
          dgam=(qc/e0)*tspl5
          xi=xe(i)/gamm1
          xpi=xpe(i)/gamm1-xe(i)*gam(i)*dgam/(gamm2*2.)
          yi=ye(i)/gamm1
          ypi=ype(i)/gamm1-ye(i)*gam(i)*dgam/(gamm2*2.)
c      convert in cm and mrd
          f(2,j)=xi*1.e02
          f(4,j)=yi*1.e02
          f(3,j)=xpi*1.e03
          f(5,j)=ypi*1.e03
c    *****  follow prtcl ifpt  not active ************
comment          if(j.eq.ifpt)
comment     *    write(49,4445) fpos5,xi,xpi,yi,ypi,e0*(gam5-1.)
comment4445      format(6(2x,e12.5))
          xnh=xnh+1.
         enddo
        enddo
        iflg=iflg+1
        xnht=xnh
        if(iflg.eq.0) then
         if(indp.eq.1) npas=96
         if(indp.eq.2) npas=48
         if(indp.eq.3) npas=24
         npas1=npas+1
         xlres=egl-xlstart
         xpas=xlres/(float(npas)*float(ipart))
         xnht=fpos5/xpas
         flgsc=.true.
         dav1(idavs,7)=tspl5
         dav1(idav,5)=xlstart*1000.
         call disp
c v29/06/2010         call emiprt(0)
         go to 501
        endif
        if(iflg.eq.1) then
         if(indp.eq.1) npas=48
         if(indp.eq.2) npas=24
         if(indp.eq.3) npas=12
         npas1=npas+1
         xlres=egl-xlstart
         xpas=xlres/(float(npas)*float(ipart))
         xnht=fpos5/xpas
        endif
        if(.not.flgsc) then
         flgsc=.true.
         call disp
         go to 501
        endif
        if(flgsc) then
         if(ichaes) then
C      Charge space
           iesp=.true.
           call cesp(eglsc)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
          endif
          call disp
          flgsc=.false.
c    reference ( c.o.g. of the bunch)
           tref=0.
           gref=0.
           do ij=1,ngood
            tref=tref+f(6,ij)
            gref=gref+f(7,ij)/e0
           enddo
           tref=tref/float(ngood)
           gref=gref/float(ngood)
           bets=sqrt(1.-1./(gref*gref))
           vref=bets*vl
comment    call emiprt(0)
           write(16,562) fpos5*1000,tref,bets,tspl5
562    format(' ref.(c.o.g.) at ',e12.5,' mm of the  cathode',/,
     * 5x,' tof: ',e12.5,' sec beta: ',e12.5,' field (MV/m) ',e12.5)
#ifdef F77
           write(6,1562) fpos5*1000,bets,cr
1562   format(' EGUN: at ',e12.5,' mm from the cathode; beta: ',e12.5,
     * a1,$)
#else
           write(6,1562,ADVANCE='NO') fpos5*1000,bets,cr
1562   format(' EGUN: at ',e12.5,' mm from the cathode; beta: ',e12.5,
     * a1)
#endif
         go to 501
        endif
500     continue
        write(6,*)
        write(6,*) 'EGUN calculation finished'
c    reference ( c.o.g. of the bunch)
        tref=0.
        gref=0.
        do i=1,ngood
         tref=tref+f(6,i)
         gref=gref+f(7,i)/e0
        enddo
        tref=tref/float(ngood)
        gref=gref/float(ngood)
        bets=sqrt(1.-1./(gref*gref))
        write(16,561) tref,bets
561     format(' ref. at output of the DC gun',/,
     *         5x,' tof: ',e12.5,' sec beta: ',e12.5)
        vref=bets*vl
C       new magnetic rigidity of the reference
        xmor=xmat*bets*gref
        boro=33.356*xmor*1.e-01/qst
        dav1(idavs,6)=bets
        dav1(idavs,36)=ngood
C   plots
        CALL STAPL(davtot*10.)
        call emiprt(0)
        return
        end
       SUBROUTINE scheff1_swesson(int)
C 21-Apr-2014: This routine is NOT used by DYNAC
comment       SUBROUTINE scheff1(int)
c  ..............................................................................
c    non relativistic SCHEFF space charge method (version swesson)
c         SCE(i) (for ISCSP=3 only)
c         int=0: initiate the arrays
c      SCE(1)=BEAMC (as above)
c      SCE(2)=r extension in rms multiples
c      SCE(3)=z extension in rms multiples
c      SCE(4)=no. of radial mesh intervals (le 20)
c      SCE(5)=no. of longitudinal mesh intervals (le 40)
c      SCE(6)=no. of adjacent bunches, applicable for buncher studies
c	     and should be 0 for linac dynamics
c      SCE(7)=pulse length, if not beta lambda.(transport studies)
c	     distance bewteen beam pulses input zero to get default
c	     "beta lambda"; units are cm
c      SCE(8)=determines frequency of calculating mesh size. see pard
c      SCE(9)=option to integrate space charge forces over box
c	     if.eq.0. no integration.  see sub gaus for further
c	     explanation.
c      SCE(10)=same meaning as SCE10 above
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/CMPTE/IELL
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/DIMENS/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/CDEK/DWP(iptsz)
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/DCSPA/IESP
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/posc/xpsc
       LOGICAL ichaes,iesp
       common/bg/bsc,gsc,phis,wsync
c   named common bg, contains beta and gamma as determined in pardyn.
c   they may be values at beginning of cell, end of cell, or mid-way.
       common/fldcom/ rp, zp,pl,opt,nip
       common/spacech1/ rm(21), zm(41), rs(20), ers(16800), ezs(16800),
     1 ez(861), aa(800), rssq(20), zzs(41)
     2 ,er(861),rss(20), ismax(40), iemax(41)
       common/rcshef/sce(20)
c          set up field tables
       if(int.eq.0) then
         beami=beamc/1000.0
         wavel=2.*pi*vl/fh
         freq=fh/(2.*pi)
         btazro=vref/vl
         frrms=sce(2)
         fzrms=sce(3)
         nr=idint(sce(4))
         nz=idint(sce(5))
         nip=idint(sce(6))
         opt=sce(9)
comment        if(ncell.gt.0) pl=cell(4,ncell)
         if(sce(7).gt. 0.) pl=sce(7)
         pl=btazro*wavel
         nr1=nr+1
         nz1=nz+1
         im1=nr*nz
         im2=nr1*nz1
         im3=nr1*nz
         na=1
         nb=ngood
         nq=nb-na+1
         return
       endif
       beami=beamc/1000.0
       IF (beami.eq.0. .OR. scdist.eq.0.) return
       IELL=IELL+1
       write(16,*) ' *call SCHEFF  ',iell
       call pintim
       CALL sizrms(0,xrms,yrms,zrms,zz)
       write(16,6875) xrms,yrms,zrms
c   write the size rms in the output file 17
comment       write(17,25) iell,xrms,yrms,zrms
comment25     format(2x,i5,3(2x,e12.5))
6875   format(' RMS size(m)',e12.5,2x,e12.5,2x,e12.5)
       rrms=sqrt(xrms*xrms+yrms*yrms)
c    change unit: m==>cm
       zrms1=zrms*100.
       rrms=rrms*100.
       dr=rrms*frrms/float(nr)
       dz=zrms1*fzrms/float(nz)
       rmax=float(nr)*dr
c          load rm, zm, rs, zs
       rm(1)=0.0
       do i=2,nr1
         rm(i)=float(i-1)*dr
         rssq(i-1)=.5*(rm(i-1)**2+rm(i)**2)
         rss(i-1)=0.5*(rm(i-1)+rm(i))
         rs(i-1)=sqrt(rssq(i-1))
       enddo
       zs=.5*dz
       do i=1,nz1
         zm(i)=float(i-1)*dz
         zzs(i)=zm(i)+zs
       enddo
       hl=float(nz)*zs
c          load ers and ezs
c     mesh dimensions are in cm. ers and ezs are in 1/cm.
c     c1, c2 and c3 are in cm., and c4 is in mev-cm.
c     q=coulombs/point.   (2/pi)*e/(4*pi*epsilon)=572167 cm mev/coul.
       q=beami/(freq*float(nq))
       c1=572167.*q/xmat
       l=0
       do k=1,nr
         rfac=(rm(k+1)**2-rm(k)**2)*dz/2.
         if(opt.eq.0.) rfac=1.
         do j=1,nz
           zp=zm(j+1)
           do i=1,nr1
             rp=rm(i)
             if(opt.eq.0.) call flds(rs(k),zs,er1,ez1)
             if(opt.eq.0.) go to 35
             call gaus (rm(k),rm(k+1),zm(1),zm(2),opt,er1,ez1)
35           l=l+1
             ers(l)=c1*er1/rfac
             ezs(l)=c1*ez1/rfac
           enddo
         enddo
       enddo
       IF (beamc.eq.0. .OR. scdist.eq.0.) return
       dz1=scdist/100.
       dist=scdist
       WRITE(16, *) ' fields acting length(cm): ',DIST
c          evaluate and apply space charge effects.
c     phimc=phi of mesh center.
C   Shifts particle coordinates to a single point in time. Uses
C   a linear shift
C      Beam c.g.
       xbar=0.
       ybar=0.
       zbar=0.
       brmoy=0.
       trmoy=0.
       do np=1,ngood
         gpai=f(7,np)/xmat
         brmoy=brmoy+sqrt(1.-1./(gpai*gpai))
         trmoy=trmoy+f(6,np)
       enddo
       trmoy=trmoy/float(ngood)
       phimc=trmoy*fh
       pbar=phimc
       beta=brmoy/float(ngood)
       gsc=1./sqrt(1.-beta*beta)
       bg=beta*gsc
       c2=beta*wavel/(2.*pi)
       gmsq=1.+bg**2
       c3=dist/(bg*beta*gmsq)
       c4=dist*xmat
       gam=sqrt(gmsq)
       c5=1./(gam*(gam+1.))
c          evaluate ng, xbar, ybar
       ng=0
       xbar=0.
       ybar=0.
       xsq=0.
       ysq=0.
       do np=1,ngood
         ng=ng+1
         x=f(2,np)
         y=f(4,np)
         xbar=xbar+x
         ybar=ybar+y
         xsq=xsq+x**2
         ysq=ysq+y**2
       enddo
       eng=float(ngood)
       xbar=xbar/eng
       ybar=ybar/eng
c  the mesh center is phi*syn
       xsq=xsq/eng
       ysq=ysq/eng
       epsq=sqrt((xsq-xbar*xbar)/(ysq-ybar*ybar))
       epsqi=1./epsq
       xfac=2./(epsq + 1.)
       yfac=epsq*xfac
c          clear and load bins
       ng=0
       do i=1,im1
         aa(i)=0.0
       enddo
       do 120 np=1,ngood
         rsq=(f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq
c     i=sqrt(rsq)/dr+1.
         r=sqrt(rsq)
         halfdr=dr*0.5
         i=idint(r/dr+1.0)
         if (i.gt.nr) go to 120
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if (abs(z).ge.hl) go to 120
c------distribute charge among adjacent bins.
         ng=ng+1
         zz=z+hl
         jm1=idint(zz/dz+1.)
         i1=i+1
c     if (rsq.lt.rssq(i)) i1=i-1
         if (rsq.lt.rss (i)) i1=i-1
         if (i1.lt.1) i1=1
         if (i1.gt.nr) i1=nr
         j1=jm1+1
         if (zz.lt.zzs(jm1)) j1=jm1-1
         if (j1.lt.1) j1=1
         if (j1.gt.nz) j1=nz
         a=1.
c     if (i1.ne.i) a=(rsq-rssq(i1))/(rssq(i)-rssq(i1))
         if (r.gt.halfdr)then
           rminsq=(r-halfdr)**2
           rmaxsq=(r+halfdr)**2
           if (i1.lt.i) then
             a=(rmaxsq-rm(i)**2)/(rmaxsq-rminsq)
           else
             a=(rm(i1)**2-rminsq)/(rmaxsq-rminsq)
           endif
         endif
         b=1.-a
         cc=1.
         if (j1.ne.jm1) cc=(zz-zzs(j1))/(zzs(jm1)-zzs(j1))
         d=1.-cc
         k=(jm1-1)*nr+i
         aa(k)=aa(k)+a*cc
         k=k+i1-i
         aa(k)=aa(k)+b*cc
         k=(j1-1)*nr+i
         aa(k)=aa(k)+a*d
         k=k+i1-i
         aa(k)=aa(k)+b*d
  120  continue
       eng=ng
c          find ismax for each j
c  v02/04/2010       do 140 j=1,nz
c  v02/04/2010         l=(j-1)*nr
c  v02/04/2010         k=nr
c  v02/04/2010         do 130 i=1,nr
c  v02/04/2010           m=l+k
c  v02/04/2010           if (aa(m)) 130,130,140
c  v02/04/2010  130    k=k-1
c  v02/04/2010  140  ismax(j)=k
      do j=1,nz
       l=(j-1)*nr
       k=nr
       do i=1,nr
        m=l+k
        if(aa(m).le.0.00) then
         k=k-1
         go to 130
        else
         go to 140
        endif
130     continue
       enddo
140       ismax(j)=k
      enddo
c        find iemax for each j
         iemax(1)=1+ismax(1)
         do j=2,nz
           iemax(j)=1+max0(ismax(j-1),ismax(j))
         enddo
         iemax(nz1)=1+ismax(nz)
c        set er and ez to zero
       do i=1,im2
         er(i)=0.0
         ez(i)=0.0
       enddo
c          sum up fields
       do 220 js=1,nz
       js1=js+1
       ism=ismax(js)
       if (ism.eq.0) go to 220
       do 210 is=1,ism
       l=(js-1)*nr+is
       a1=aa(l)
       if (a1.eq.0.) go to 210
       l=(is-1)*im3
       do 180 je=1,js
       k1=l+(js-je)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 180
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)-a1*ezs(k)
       enddo
  180  continue
       do 200 je=js1,nz1
       k1=l+(je-js1)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 200
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)+a1*ezs(k)
       enddo
  200 continue
  210 continue
  220 continue
c          evaluate and apply impulse
       rrmax=0.
       zzmax=0.
       zzmin=1000.
       npz=0
       npr=0
       do 270 np=1,ngood
         r=sqrt((f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq)
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if(z.ge.zzmax) zzmax=z
         if(z.lt.zzmin) zzmin=z
         if(r.ge.rrmax) rrmax=r
         if (r.eq.0.) r=.000001
         xor=(f(2,np)-xbar)*xfac/r
         yor=(f(4,np)-ybar)*yfac/r
         if (r.gt.rmax) then
           npr=npr+1
           go to 230
         endif
c         zph=f(6,np)*fh
c         z=-c2*(zph-phimc)
c         if(z.ge.zzmax) zzmax=z
c         if(z.lt.zzmin) zzmin=z
         if (abs(z).gt.hl) then
           npz=npz+1
           go to 230
         endif
c        interpolate impulse within mesh.
         rb=r/dr
         i=idint(1.0+rb)
         a=rb-float(i-1)
         b=1.0-a
         zb=(z+hl)/dz
         j=idint(1.0+zb)
         c=zb-float(j-1)
         d=1.0-c
         l=i+(j-1)*nr1
         m=l+nr1
         crp=c3*(d*(a*er(l+1)+b*er(l))+c*(a*er(m+1)+b*er(m)))
         cen=c4*(d*(a*ez(l+1)+b*ez(l))+c*(a*ez(m+1)+b*ez(m)))
         crp=crp*abs(f(9,np))
         cen=cen*abs(f(9,np))
         go to 260
c        estimate impulse based on point charge at xbar,ybar,pbar.
  230    continue
         d=sqrt(z**2+r**2)
         rod3=r/d**3
         zod3=z/d**3
         if (nip.eq.0) go to 250
c        include neighboring bunches.
         do i=1,nip
           xi=i
           do j=1,2
             s=z+xi*pl
             d=sqrt(s**2+r**2)
             rod3=rod3+r/d**3
             zod3=zod3+s/d**3
             xi=-xi
           enddo
         enddo
c        evaluate impulse.
  250    continue
         crp=eng*c1*c3*rod3*pi/2.
         cen=eng*c1*c4*zod3*pi/2.
         crp=crp*abs(f(9,np))
         cen=cen*abs(f(9,np))
c        apply impulse
  260    continue
C        convert from mrad to rad
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
         dwc=f(7,np)-xmat
         dxp=crp*xor - f3*cen*c5/dwc
         dyp=crp*yor - f5*cen*c5/dwc
         if(.not.iesp) then
c     load the entrance beam parameters for cavities or gaps
           do js=1,7
             f(js,np)=fs(js,np)
           enddo
           f(3,np)=f(3,np)+dxp*1000.
           f(5,np)=f(5,np)+dyp*1000.
           F(2,np)=F(2,np)-DZ1*DXP*100.*xpsc
           F(4,np)=F(4,np)-DZ1*DYP*100.*xpsc
           DWP(np)=cen
         else
           f(3,np)=f(3,np)+dxp*1000.
           f(5,np)=f(5,np)+dyp*1000.
           f(7,np)=f(7,np)+cen
         endif
  270  continue
       RETURN
       END
       SUBROUTINE cpardyn(pib)
c *******************************************************************************************************************************
c --- make up a list of cell by cell RFQ data based on the file generated by the external code
c
c     nc: cell number
c     ityp = 0: Standard accelerating cell
c     ityp = 1: Transition cell of type T
c     ityp = 2: Transition cell of type E
c     ityp = 3: Transition cell of type M
c     ityp = 4: Fringe-field region from the potential function (after type T, M or accelerating cell)
c     ityp = 5: Radial matching section based on coefficients
c     ityp = 6: Radial matching section based on RMS shape defined in an external file
c     ityp = 7: Fringe-field region; the profile is defined in an external file (after type T, M or accelerating cell)
c     a(1): cell length (cm)
c     a(2): coefficient A10 (no dimension)
c     a(3): smallest aperture  a  of vanes (cm)
c     a(4): modulation factor  m  (no dimension)
c     a(5): mean aperture ( r0 ) of vanes (middle of cells)   (cm)
c     a(6): transverse radius of curvature (rh0) of the surface of electrodes at the vane tip (cm)
c     a(7): phase (deg) (at the entrance of cells)
c     a(8): factor FVOLT to be applied at the intervane potentiel (only for particles)
c              VOLT =  (1 + FVOLT )* (intervane voltage)
c     a(9): intervane-voltage (KV)
c
c     ipar: cell parity (like for PARMTEQ cells)
c --- Note: the last line in file must be: 0  0  0. 0. 0. 0. 0. 0. 0. 0. 0. 0
c
c ***********************************************************************************************************************************
       implicit real*8 (a-h,o-z)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/rf1ptq/tvolt,avolt,fph,mlc,nceltot
       common/rf2ptq/rfq1(500),rfq2(500),rfq3(500),rfq4(500),rfq6(500)
     *        ,rfq7(500),rfq8(500),rfq9(500)
       common/rf5ptq/tdvolt,rfq10(500),rfq11(500)
       common/rfq3ptq/itype(500),ipari(500),evens,evenr
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       common/bonda/cbx(500),bbx(500),ablx(500),cby(500),
     *        bby(500),ably(500)
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/spff/xspf(400),yspf(400),sf(500),pf(500),qf(500)
       common/rms_prfl/npt,npf
       dimension a(16),vptq(16)
       character*80 vprof,fprof
       logical even,evens,evenr
c  Do so by shifting particles belonging to the same bunch outside the (+/-) pib/2 (rad) window
c   w.r.t.the REF to inside the (+/-) pib/2 window w.r.t.the COG
       if(abs(pib).gt.6.*0.) then
         pib=pi
         call accep_rfq(pib)
       endif
       erest=xmat
c --- for simulating PARMTEQ cells, odd cells have A01 positive and even cells have A01 negative
c ---  even: logical flag, even = true ==> A01 is positive, otherwise A01 is negative
       even=.true.
       netc=0
50     continue
c --- read the data of the RFQ from the input file
c ****************************************************************
c ---   read 11 parameters from unit 27), file 'myfile', in the form:
c      nc  ityp intervane-voltage (KV) cl(cm)  A10  a(cm)  m  r0(cm) rho(cm) phase(deg) fvolt ipar
       netc=netc+1
       read(27,*) nc,ityp,(vptq(j),j=1,9),ipar
       if(ityp.eq.6) then
c if the RMS is of type 6, get the filename of the RMS descriptor
         read(27,'(A80)') vprof
         write(16,*) 'Read RMS from file ',vprof
       endif
       if(ityp.eq.7) then
c if the Fringe Field is of type 7, get the filename of the FF descriptor
         read(27,'(A80)') fprof
         write(16,*) 'Read Fringe Field from file ',fprof
       endif
       a(1)=vptq(2)
       a(2)=vptq(3)
       a(3)=vptq(4)
       a(4)=vptq(5)
       a(5)=vptq(6)
       a(6)=vptq(7)
       a(7)=vptq(8)
       a(8)=vptq(9)
       a(9)=vptq(1)
c **************************************************
c convert intervane voltage in MV
       a(9)=a(9)*1.e-03
       if(nc.eq.1) tdvolt=a(9)
c ----  stop at the last line with nc = 0
       if(nc.eq.0) go to 60
       if(nc.gt.nceltot) go to 60
       itype(nc)=ityp
       ipari(nc)=ipar
       if(itype(nc).eq.5) then
c Radial matching section from potential function
c ---------------------------------------------------------------------------------------------
c Note: The total length L(rms) of the RMS must be the sum of the RMS-cell lengths in the
c       Parmteq file. rho is equal to the distance between the axis and the vane, r0, at the
c       exit of the RMS
c       we have used the approximation of the modified Bessel functions for small arguments
c       and we have negleted the terms of higher order in r greater than 2
c ----------------------------------------------------------------------------------------------
         write(16,*) 'RFQ RMS from potential function'
         cl=a(1)
         r0=a(5)
         rh0=a(6)
         phim=a(7)
         fact=1.+a(8)
         xk=pi/(2.*cl)
         xq0=1./6.*xk*xk*rh0*rh0
cold sv         aq=9./(8.*xq0)
         aq=1./xq0
c ---   limits x-vane y-vane
         cbx(nc)=1.5*r0*1.e-02
         cby(nc)=1.5*r0*1.e-02
         ablx(nc)=0.
         ably(nc)=0.
         bbx(nc)=0.
         bby(nc)=0.
c ---  convert all parameters in units (MeV, m)
c     rfq1(nc): Aq (no dimension)
c     rfq2(nc): not used
c     rfq3(nc): RMS length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor to be applied at the inter-vane potentiel (only for particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (kV)
c     rfq11(nc):intervane voltage applied to the particles (kV)
c
         rfq1(nc)=aq
         rfq2(nc)=0.
         rfq3(nc)=cl*1.e-02
         rfq4(nc)=phim
         rfq7(nc)=r0*1.e-02
         rfq9(nc)=fact
         rfq10(nc)=(1.+tvolt)*a(9)
         rfq11(nc)=(1.+avolt)*a(9)
         go to 50
       endif
       if(itype(nc).eq.6) then
c  read from disk the RMS profile in the form (z,R(z)), with z,R in meter
c  npt is the number of (z,R(z)) coordinates
         open(67,file=vprof,status='unknown')
c         read(67,*)npt
         npt=0
         do
           read(67,*,END=777) xspl(npt+1),yspl(npt+1)
           npt=npt+1
         enddo
777      write(16,*) 'RFQ RMS from file with ',npt,' data points'
         call deriv2(npt)
c cl; length of the RMS (in M)
         cl=xspl(npt)
         r0=yspl(npt)
         phim=a(7)
         fact=1.+a(8)
c ---   limits x-vane y-vane
         cbx(nc)=1.5*r0
         cby(nc)=1.5*r0
         ablx(nc)=0.
         ably(nc)=0.
         bbx(nc)=0.
         bby(nc)=0.
c     rfq1(nc) : not used
c     rfq2(nc): not used
c     rfq3(nc): RMS length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): distance axe-vanes at the end of the RMS (m)
c     rfq9(nc): factor to be applied at the inter-vane potentiel (only for particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (kV)
c     rfq11(nc):intervane voltage applied to the particles (kV)
c       rfq1(nc)=0.
        rfq2(nc)=0.
        rfq3(nc)=cl
        rfq4(nc)=phim
        rfq7(nc)=r0
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        close(67)
ccc        write(16,579)nc,npt,cl,phim,r0,rfq10(nc),rfq11(nc)
ccc579     format(2(2x,i4),5(2x,e12.5))
        go to 50
       endif
c  End Radial matching section
c --- Standard accelerationg cell
       if(itype(nc).eq.0) then
        cl=a(1)
        a10=a(2)
        xa=a(3)
        xm=a(4)
        r0=a(5)
        rh0=a(6)
c fph : factor affecting the phase at entrance of cells (only available for type = 0)
        a(7)=fph*a(7)
        phim=a(7)
        fact=1.+a(8)
c    coefficient A01 (1/(cm*cm))
        alpha=rh0/r0
        a01=3.*(1.+5.*alpha)/(2.*r0*r0*(1.+7.*alpha))
c    coefficient A03: (1/(cm**6) )
        a03=-(1.+alpha)/(2.*r0**6*(1.+7.*alpha))
c  coefficient A12: no dimensions
        xam=xa*xm
        a12=0.
        xk1=1.-a01*xa*xa-a03*(xa**6)
        yk2=-1.+a01*xam*xam+a03*(xam**6)
c   Bessel functions: I0(ka) I0(mka) I4(ka) I4(mka)
        xk=pi/cl
        za0=xk*xa
        zam=za0*xm
        no=0
        bi0=bint(no,za0)
cold        bi0=1.+za0*za0/4.+za0**4/64.+za0**6/2304.+za0**8/1.47456e05
cold        bi0=bi0+za0**10/1.47456e07
        bim=bint(no,zam)
cold        bim=1.+zam*zam/4.+zam**4/64.+zam**6/2304.+zam**8/1.47456e05
cold        bim=bim+zam**10/1.47456e07
        no=4
        bi4=bint(no,za0)
cold        bi4=(za0/2.)**4/24.
cold        zaa0=za0*za0/4.
cold        bi4=bi4+(za0/2.)**4*zaa0/120.
cold        bi4=bi4+(za0/2.)**4*zaa0*zaa0/1440.
        bim4=bint(no,zam)
cold        zamm=zam*zam/4.
cold        bim4=(zam/2.)**4/24.
cold        bim4=bim4+(zam/2.)**4*zamm/120.
cold        bim4=bim4+(zam/2.)**4*zamm*zamm/1440.
        den1=bim4*bi0-bim*bi4
        if(abs(den1).gt.1.e-09) a12=(yk2*bi0-xk1*bim)/den1
c ********************************************************************************************
c NOTE: the coefficient A10 is read in the file 'myfile' but it also can be computed in Dynac
c        from multipolar expansions (see the two following fortran lines):
ccc        a10=0.
ccc        if(abs(den1).gt.1.e-09) a10=(xk1*bim4-yk2*bi4)/den1
c ************************************
c        or from first order computations (see the two following lines):
ccc        dencc=xm*xm*bi0+bim
ccc        a10=(xm*xm-1.)/dencc
c ********************************************************************************************
c ---  convert all parameters in units (MeV, m)
c     rfq1(nc): A01 ( 1/(m*m) )
c     rfq2(nc): A10 (no dimension)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq6(nc): A12 (no dimension)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq8(nc): A03 (1/(m**6)
c     rfq9(nc): error factor F = 1 + a(8)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c ---  odds cells have a01 positive, even cells have a01 negative
c
        rfq1(nc)=a01*1.e04
c ---  limits x-vane y-vane
        ncf=(ipari(nc)/2)*2-ipari(nc)
        if(ncf.ne.0) then
c  odd number standard cell limits
         cbx(nc)=xa
         bbx(nc)=(4.*r0-xm*xa-3.*cbx(nc))/cl
         ablx(nc)=xa*xm/(cl*cl)-bbx(nc)/cl-cbx(nc)/(cl*cl)
         cby(nc)=xa*xm
         bby(nc)=(4.*r0-xa-3.*cby(nc))/cl
         ably(nc)=xa/(cl*cl)-bby(nc)/cl-cby(nc)/(cl*cl)
       else
c  even number standard cell limits
        cbx(nc)=xm*xa
        bbx(nc)=(4.*r0-xa-3.*cbx(nc))/cl
        ablx(nc)=xa/(cl*cl)-bbx(nc)/cl-cbx(nc)/(cl*cl)
        cby(nc)=xa
        bby(nc)=(4.*r0-xm*xa-3.*cby(nc))/cl
        ably(nc)=xm*xa/(cl*cl)-bby(nc)/cl-cby(nc)/(cl*cl)
      endif
c8526  continue
c   conversion: cm --> m
         cbx(nc)=cbx(nc)*1.e-02
         ablx(nc)=ablx(nc)*1.e02
         cby(nc)=cby(nc)*1.e-02
         ably(nc)=ably(nc)*1.e02
c  end limits ******************************
        rfq2(nc)=a10
        rfq3(nc)=cl*1.e-02
        rfq4(nc)=phim
        rfq6(nc)=a12
        rfq7(nc)=r0*1.e-02
        rfq8(nc)=a03*1.e12
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        go to 50
c  endif of itype(nc) = 0 (accelerating cell)
      endif
      if(itype(nc).eq.1) then
c ------ T-cell(ityp = 1) (must follows a standard cell and be followed by an M-cell or by fringe-field region)
        cl=a(1)
        xa=a(3)
        xm=a(4)
        r0=a(5)
        phim=a(7)
        fact=1.+a(8)
        xk=pi/(2.*cl)
        za0=xk*xa
        zam=za0*xm
c   Bessel functions
        no=0
        bi0=bint(no,za0)
        bim=bint(no,zam)
        bi3=bint(no,3.*za0)
        bim3=bint(no,3.*zam)
c   coefficients A10 and A30 (no dimensions)
        t10k=xm*xm*bi0+bim
        t30k=xm*xm*bi3+bim3
        zr0=xk*r0
        zr3=3.*zr0
        bir0=bint(no,zr0)
        bir3=bint(no,zr3)
        alpk=0.
        if(abs(bir3).ne.6.*0.) alpk=bir0/bir3
        dtk=t10k+alpk*t30k/3.
        a10=0.
        if(abs(dtk).ne.6.*0.) a10=(xm*xm-1.)/dtk
        a30=alpk*a10/3.
c -- limits x-vane  y-vane
        cbx(nc)=xm*xa
        bbx(nc)=(4.*r0-xa-3.*cbx(nc))/(2.*cl)
        ablx(nc)=xa/(4.*cl*cl)-bbx(nc)/(2.*cl)-cbx(nc)/(4.*cl*cl)
        cby(nc)=xa
        bby(nc)=(4.*r0-xm*xa-3.*cby(nc))/(2.*cl)
        ably(nc)=xm*xa/(4.*cl*cl)-bby(nc)/(2.*cl)-cby(nc)/(4.*cl*cl)
c   conversion: cm --> m
         cbx(nc)=cbx(nc)*1.e-02
         ablx(nc)=ablx(nc)*1.e02
         cby(nc)=cby(nc)*1.e-02
         ably(nc)=ably(nc)*1.e02
c  end limits ******************************
c ---  convert all parameters in units (MeV, m)
c     rfq1(nc): A30 (no dimension)
c     rfq2(nc): A10 (no dimension)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor to be applied at the inter-vane potentiel (only for particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
        rfq1(nc)=a30
        rfq2(nc)=a10
        rfq3(nc)=cl*1.e-02
        rfq4(nc)=phim
        rfq7(nc)=r0*1.e-02
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        go to 50
c  endif for T-cell
       endif
      if(itype(nc).eq.2) then
c --- E-cell(itype = 2) must follow a RMS and be folllowed by an even standard cell
        cl=a(1)
        xa=a(3)
        xm=a(4)
        r0=a(5)
        a(7)=fph*a(7)
        phim=a(7)
        fact=1.+a(8)
        xk=pi/(2.*cl)
        za0=xk*xa
        zam=za0*xm
        xam=xa*xm
c   Bessel functions
        no=0
        bi0=bint(no,za0)
        bim=bint(no,zam)
        bi3=bint(no,3.*za0)
        bim3=bint(no,3.*zam)
c   coefficients A10 and A30 (no dimensions)
        t10k=xm*xm*bi0+bim
        t30k=xm*xm*bi3+bim3
        zr0=xk*r0
        zr3=3.*zr0
        bir0=bint(no,zr0)
        bir3=bint(no,3.*zr3)
        alpk=0.
        if(abs(bir3).ne.6.*0.) alpk=bir0/bir3
        dtk=t10k+alpk*t30k/3.
        a10=0.
        if(abs(dtk).ne.6.*0.) a10=(xm*xm-1.)/dtk
cold   E-cell is a even cell number (followed by a odd number standard cell
cold        ncf=(nc/2)*2-nc
cold        if(ncf.eq.0) a10=-a10
        a30=-alpk*a10/3.
c  cell limits x-vane y-vane
        cbx(nc)=r0
        bbx(nc)=(xm*xa-r0)/cl
        ablx(nc)=0.
        cby(nc)=r0
        bby(nc)=(xa-r0)/cl
        ably(nc)=0.
c   conversion: cm --> m
         cbx(nc)=cbx(nc)*1.e-02
         cby(nc)=cby(nc)*1.e-02
c ---  convert in units (MeV, m)
c     rfq1(nc): A30 (no dimensions)
c     rfq2(nc): A10 (no dimensions)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor applied at the inter-vane potentiel
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
        rfq1(nc)=a30
        rfq2(nc)=a10
        rfq3(nc)=cl*1.e-02
        rfq4(nc)=phim
        rfq7(nc)=r0*1.e-02
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        go to 50
       endif
      if(itype(nc).eq.3) then
c   M-cell  (ityp = 3)
        cl=a(1)
        xa=a(3)
        xm=a(4)
        r0=a(5)
        phim=a(7)
        fact=1.+a(8)
c  cell limits x-vane = y-vane = average radius
        cbx(nc)=r0
        bbx(nc)=0.
        ablx(nc)=0.
        cby(nc)=r0
        bby(nc)=0.
        ably(nc)=0.
c   conversion: cm --> m
        cbx(nc)=cbx(nc)*1.e-02
        cby(nc)=cby(nc)*1.e-02
c  end vanes limits ******************************
c ---  convert the parameters in units (MeV, m)
c     rfq1(nc): A30 = 0
c     rfq2(nc): A10 = 0
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor applied at the inter-vane potentiel(F = 1 + a(15))
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)c
        rfq1(nc)=0.
        rfq2(nc)=0.
        rfq3(nc)=cl*1.e-02
        rfq4(nc)=phim
        rfq7(nc)=r0*1.e-02
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        go to 50
c endif M-cell
      endif
c   Fringe-field region after T-cell or M-cell  (Type = 4)
       if(itype(nc).eq.4) then
        cl=a(1)
        a10=a(2)
        xa=a(3)
        xm=a(4)
        r0=a(5)
        rh0=a(6)
        phim=a(7)
        fact=1.+a(8)
c  fringe field region limits
        cbx(nc)=r0*1.5
        bbx(nc)=0.
        ablx(nc)=0.
        cby(nc)=r0*1.5
        bby(nc)=0.
        ably(nc)=0.
c   conversion: cm --> m
        cbx(nc)=cbx(nc)*1.e-02
        cby(nc)=cby(nc)*1.e-02
c  end limits ******************************
c ---  convert the parameters in units (MeV, m)
c     rfq1(nc): A01 ( 1/(m*m) )
c     rfq2(nc): A10 (no dimensions)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor applied at the inter-vane potentiel (F = 1 + a(15))
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
cold        rfq1(nc)=a01*1.e04
        rfq2(nc)=a10
        rfq3(nc)=cl*1.e-02
        rfq4(nc)=phim
        rfq7(nc)=r0*1.e-02
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
       endif
c********************************************
c  fringe field region the profil is read from disk
       if(itype(nc).eq.7) then
c  read the ff profile in the form (z,R(z)) in M
c    npt is the number of coordinates (z,R(z)
         open(67,file=fprof,status='unknown')
c         read(67,*)npf
         npf=0
         do
           read(67,*,END=888) xspf(npf+1),yspf(npf+1)
           npf=npf+1
         enddo
888      write(16,*) 'RFQ FF from file with ',npf,' data points'
         reca=xspf(1)
         do i=1,npf
           xspf(i)=xspf(i)-reca
         enddo
c length of the fringe field  m
         cl=xspf(npf)
         call derif2(npf)
         phim=a(7)
         fact=1.+a(8)
         write(16,*) 'Fringe Field length=',cl,' m'
         write(16,*) 'Fringe Field Phase=',phim
         write(16,*) 'Fringe Field Field Factor=',fact
c ---   limits x-vane y-vane
         cbx(nc)=100.
         cby(nc)=100.
         ablx(nc)=0.
         ably(nc)=0.
         bbx(nc)=0.
         bby(nc)=0.
c     rfq1(nc) : not used
c     rfq2(nc): not used
c     rfq3(nc): Fringe field length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): distance axe-vanes at the start of the FF (m)
c     rfq9(nc): factor to be applied at the inter-vane potentiel (only for particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
        rfq1(nc)=0.
        rfq2(nc)=0.
        rfq3(nc)=cl
        rfq4(nc)=phim
        rfq7(nc)=r0
        rfq9(nc)=fact
        rfq10(nc)=(1.+tvolt)*a(9)
        rfq11(nc)=(1.+avolt)*a(9)
        close(67)
ccc        write(16,579)nc,npt,cl,phim,r0,rfq10(nc),rfq11(nc)
ccc579     format(2(2x,i4),5(2x,e12.5))
       endif
c********************************************
60     continue
       if(netc.lt.nceltot) then
         write(6,*) 'Error: Parameter NCELTOT after RFQPTQ entry in',
     *   ' DYNAC input file is ',NCELTOT
         write(6,*) 'This is larger than the',
     *   ' number of cells in the RFQ datafile, which is ',netc
         write(16,*) 'Error: Parameter NCELTOT after RFQPTQ entry in',
     *   ' DYNAC input file is ',NCELTOT
         write(16,*) 'This is larger than the',
     *   ' number of cells in the RFQ datafile, which is ',netc
         stop
       endif
       call rfq_parm
       return
       end
       function bint(n,z)
c ---  integral representation of modified Bessel functions
c      n integer order    z argument
        implicit real*8 (a-h,o-z)
        COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
        dimension ui(16),wi(16)
c   GAUSS n=16 de -1. a 1
       DATA (UI(J),J=1,16)/-.9894009,-.9445750,-.8656312,-.7554044,
     1   -.6178762,-.4580168,-.2816036,-.0950125,
     2    .0950125,.2816036,.4580168,.6178762,.7554044,.8656312,
     3    .9445750,.9894009/
        DATA (WI(J),J=1,16)/.0271524,.0622535,.0951585,.1246288,
     1   .1495960,.1691565,.1826034,.1894506,.1894506,.1826034,
     2   .1691565,.1495960,.1246288,.0951585,.0622535,.0271524/
        bint=0.
        do i=1,16
         thet=pi/2.*(1.+ui(i))
         fln=float(n)
         cthet=cos(thet)
         fonc=exp(cthet*z)*cos(fln*thet)
         bint=bint+fonc*wi(i)
        enddo
        bint=bint/2.
        return
        end
       SUBROUTINE accep_rfq(pib)
c   .....................................................................
c  Do so by shifting particles belonging to the same bunch from outside
c  to inside (+/-) pi w.r.t.the synchronous particle
c   .....................................................................
       implicit real*8 (a-h,o-z)
       parameter(iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common/faisc/f(10,iptsz),imax,ngood
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       common/sc3/beamc,scdist,sce10,cplm,ectt,apl,ichaes,iscsp
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       common/mcs/imcs,ncstat,cstat(20)
       common/tapes/in,ifile,meta
       common/etcha3/ichxyz(iptsz)
       common/rec/irec
       common/etcom/cog(8),exten(17),fd(iptsz)
       logical chasit
c       character*80 wfile
       do ite=1,3
         do i=1,ngood
           drad=(f(6,i)-tref)*fh
           if(drad.gt.pib) then
cold             f(6,i)=(f(6,i)-pi/fh)
             f(6,i)=(f(6,i)-2.*pi/fh)
           endif
           if(drad.lt.-pib) then
cold             f(6,i)=(f(6,i)+pi/fh)
             f(6,i)=(f(6,i)+2.*pi/fh)
           endif
         enddo
       enddo
       tcog=0.
c --- cog in time of bunch after reinjecting particles into the bunch
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
c       tref=tcog
       write(16,59) tcog*fh*180./pi
59     format(' Phase of COG after reinjecting particles into the ',
     *        'bunch: ',e13.7,' deg')
       return
       end
      SUBROUTINE rfq_parm
c***********************************************************************************************
c  NOTE:
c      The reference particle and particles may evolve separately or may be connected
c      Only the SCHEFF space charge method is available
c      Space charge computations are automaticaly made at the middle of each cell
c ------------------------------------------------------------------------------------------------
c  Radial matching section
c     rfq1(nc): Aq (no dimension)
c     rfq2(nc): not used
c     rfq3(nc): RMS length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor to be applied at the inter-vane potentiel (only for particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
c ---- standard accelerating cells (ityp = 0)
c     rfq1(nc): A01 ( 1/(m*m) )
c     rfq2(nc): A10 (no dimensions)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq6(nc): A12 (no dimension)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq8(nc): A03 (1/(m**6)
c     rfq9(nc): factor F = 1 + a() (only for inter-vane potentiel of particles)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c ---  odds cells have a01 positive, even cells have a01 negative
c
c ----  transition T-cell and E-cell
c     rfq1(nc): A30 (no dimensions)
c     rfq2(nc): A10 (no dimensions)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor F = 1 + a(8)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
c ----  M-cell
c     rfq1(nc): A30 = 0
c     rfq2(nc): A10 = 0
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor F = 1 + a(8)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
c ---- Fringe-field region from the potential function (F-cell)
c     rfq1(nc): A01 ( 1/(m*m) )
c     rfq2(nc): A10 (no dimensions)
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m)
c     rfq9(nc): factor F = 1 + a(8)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
c ---- Fringe-field region the profil is read on the disk (F-cell)
c     rfq1(nc): not used
c     rfq2(nc): not used
c     rfq3(nc): cell length (m)
c     rfq4(nc): phase RF (deg)
c     rfq7(nc): mean aperture of the vane r0 (m) (not used)
c     rfq9(nc): factor F = 1 + a(8)
c     rfq10(nc):intervane voltage applied to the synchronous particle (KV)
c     rfq11(nc):intervane voltage applied to the particles (KV)
c
c**********************************************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/RIGID/BORO
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/tapes/in,ifile,meta
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DCSPA/IESP
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SHIF/DTIPH,SHIFT
       common/femt/iemgrw,iemqesg
       common/posc/xpsc
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
c       common/drfq/p(9)
       COMMON/AZLIST/ICONT,IPRIN
       common/trfq/icour,ncell
       logical iesp,ichaes,shift,iemgrw,iflag
       common/itvole/itvol,imamin
       common/tofev/ttvols
       common/rf1ptq/tvolt,avolt,fph,mlc,nceltot
       common/rf2ptq/rfq1(500),rfq2(500),rfq3(500),rfq4(500),rfq6(500)
     *        ,rfq7(500),rfq8(500),rfq9(500)
       common/rf5ptq/tdvolt,rfq10(500),rfq11(500)
cold       common/rfq3ptq/itype(500),evens,evenr
       common/rfq3ptq/itype(500),ipari(500),evens,evenr
       common/bonda/cbx(500),bbx(500),ablx(500),cby(500),
     *        bby(500),ably(500)
       logical itvol,imamin
       common/conti/irfqp
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/spff/xspf(400),yspf(400),sf(500),pf(500),qf(500)
       common/rms_prfl/npt,npf
c********************************************
c v28/04/2015
       common/fcont/ifcont
       logical ifcont
c********************************************
c       logical irfqp,evens,evenr,evenf
       logical irfqp,evens,evenr
       dimension rfqdmp(500,10)
       character*1 cr
c allow for print out on terminal of cell# on one and the same line
       cr=char(13)
       irfqp=.true.
       iflag=.false.
       radian=pi/180.
       ilost=0
       twopi=2.*pi
c  convert vl in m
       vlm=vl/100.
       wavel=2.*pi*vlm/fh
       er=xmat
C      STATISTIQUES FOR PLOTS
       if(iprf.eq.1) call stapl(davtot*10.)
c  start prints in file 'short.data'
       tlgth=0.
       idav=idav+1
       iitem(idav)=15
       dav1(idav,9)=tdvolt*1000.
       dav1(idav,7)=float(nceltot)
       dav1(idav,8)=tdvolt*1000.
c --- end daves
c  ns = 18: number of steps in cells (apart for the fringe field region and RMS)
       xlrfq=0.
c   the integer nn is only used if the RMS profile is required
       nn=0
       do ncell=1,nceltot
#ifdef F77
         write(6,8254) nrtre,ncell,cr
8254     format('Transport element:',i5,'      RFQ cell :',
     *     i5,'            ',a1,$)
#else
         write(6,8254,ADVANCE='NO') nrtre,ncell,cr
8254     format('Transport element:',i5,'      RFQ cell :',
     *     i5,'            ',a1)
#endif
        ns=18
        if(itype(ncell).eq.5) ns=126
        if(itype(ncell).eq.6) ns=126
        nsm=ns/2
        r0=rfq7(ncell)
        cl=rfq3(ncell)
        davtot=davtot+cl*100.
        xlrfq=xlrfq+cl*100.
         if(itype(ncell).eq.0) then
c ------------- accelerating cells(ityp = 0)
c --- synchronous particle
           rtvolt=rfq10(ncell)
           vorsq=rfq1(ncell)*rtvolt
           av=rfq2(ncell)*rtvolt
           a12v=rfq6(ncell)*rtvolt
c --- particles
           pavolt=rfq11(ncell)
           vorb=rfq1(ncell)*pavolt*rfq9(ncell)
           avb=rfq2(ncell)*pavolt*rfq9(ncell)
           a12vb=rfq6(ncell)*pavolt*rfq9(ncell)
           a03vb=rfq8(ncell)*pavolt*rfq9(ncell)
         endif
         if(itype(ncell).eq.1) then
c ---------------- T-cell(Type = 1)
c ---- synchronous particle
           rtvolt=rfq10(ncell)
           a31v=rfq1(ncell)*rtvolt
           a10v=rfq2(ncell)*rtvolt
c ---- particles (field)
           pavolt=rfq11(ncell)
           a31vb=rfq1(ncell)*pavolt*rfq9(ncell)
           av10b=rfq2(ncell)*pavolt*rfq9(ncell)
         endif
         if(itype(ncell).eq.2) then
c --------   E-cell(Type = 2)
c ---- synchronous particle
           rtvolt=rfq10(ncell)
           a31v=rfq1(ncell)*rtvolt
           a10v=rfq2(ncell)*rtvolt
c ---- particles
           pavolt=rfq11(ncell)
           a31vb=rfq1(ncell)*pavolt*rfq9(ncell)
           a10vb=rfq2(ncell)*pavolt*rfq9(ncell)
         endif
         if(itype(ncell).eq.4) then
c -------, F-cell (ityp = 4)
c --- synchronous particle
           rtvolt=rfq10(ncell)
           av=rfq2(ncell)*rtvolt
c --- particles
           pavolt=rfq11(ncell)
           avb=rfq2(ncell)*pavolt*rfq9(ncell)
         endif
c********************************************
c  fringe field type 7
         if(itype(ncell).eq.7) then
c -------, F-cell (ityp = 7)
c --- synchronous particle
           rtvolt=rfq10(ncell)
c --- particles
           pavolt=rfq11(ncell)
           a31vb=pavolt*rfq9(ncell)
         endif
c********************************************
c -----Radial matching section from potential function(Type = 5)
         if(itype(ncell).eq.5) then
c ---- synchronous particle
           rtvolt=rfq10(ncell)
           a31v=rfq1(ncell)*rtvolt
           a10v=rfq2(ncell)*rtvolt
c ---- particles (field)
           pavolt=rfq11(ncell)
           a31vb=rfq1(ncell)*pavolt*rfq9(ncell)
           av10b=rfq2(ncell)*pavolt*rfq9(ncell)
         endif
c -----Radial matching section with vanes profile read from disk (Type = 6)
         if(itype(ncell).eq.6) then
c ---- potential applied to the synchronous particle
           rtvolt=rfq10(ncell)
c ---- potenatial applied to particles
           pavolt=rfq11(ncell)
           a31vb=pavolt*rfq9(ncell)
         endif
c   c.o.g of the bunch
         tcog=0.
         ecog=0.
         do i =1,ngood
           tcog=tcog+f(6,i)
           ecog=ecog+f(7,i)
         enddo
         tcog=tcog/float(ngood)
         ecog=ecog/float(ngood)
         gcog=ecog/er
         bcog=sqrt(1.-1./(gcog*gcog))
         wcog=ecog-er
         if(ncell.eq.1) then
c ---- shift = .false. ==> the synchronous particle is the center of gravity
           if(.not.shift) then
             tref=tcog
             bref=bcog
             vref=bref*vl
             gref=gcog
             wref=wcog
             wrefi=wref
           else
c ---- shift = .true. ==> the synchronous particle and the c.o.g are separated
             bref=vref/vl
             gref=1./sqrt(1.-bref*bref)
             wref=er*(gref-1.)
             wrefi=wref
           endif
         endif
c --- standard accelerating cell (itype = 0)
         if(itype(ncell).eq.0) cay=pi/cl
c ---  T-cell(itype = 1 , E-cell(itype = 2, fringe field region(itype = 4) RMS (itype = 5)
         if(itype(ncell).eq.1) cay=pi/(2.*cl)
         if(itype(ncell).eq.2) cay=pi/(2.*cl)
         if(itype(ncell).eq.5) cay=pi/(2.*cl)
         if(itype(ncell).eq.6) cay=pi/(2.*cl)
c --- fringe field region type 4
         if(itype(ncell).eq.4) then
           cay=pi/(2.*cl)
           NS=int(36.*CL/(BREF*WAVEL))
           if(ns.le.5) ns=6
           nsm=ns/2
         endif
c********************************************
c   Fringe field region type 7
         if(itype(ncell).eq.7) then
           cay=pi/(2.*cl)
           NS=int(36.*CL/(BREF*WAVEL))
           ns=3*ns
cold           if(ns.le.5) ns=6
           nsm=ns/2
         endif
c********************************************
c ---  M-cell (itype = 3)
cold sv         if(itype(ncell).eq.3) cay=pi/cl
         xl=cl/float(ns)
         hl=.5*xl
c----  scl: space charge length in SCHEFF unit (cm)
         scl=cl*100.
c----  phini: phase of the synchronous at input of the cell
         phini=-tref*fh+rfq4(ncell)*radian
         if(ncell.eq.1) then
           write(16,178)
178        format(/,' Dynamics at the input',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
           write(16,1788) bcog,gcog,wcog,tcog*fh*180./pi,tcog
1788       FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
           WRITE(16,165) bref,gref,wref,tref*fh*180./pi,tref
165        FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
           write(16,*)
           write(75,9627)
9627    format(5x,'ncell',4x,'A01(m-2)',8x,'A10',
     *         12x,'A12',11x,'r0(m)',11x,'A03(m-6)')
          write(70,9977)
9977      format(5x,'ncell',4x,'Zcell(m)',7x,'Z(m)',
     *    9x,'Phi(deg)',8x,'Pho(deg)',7x,'Wsyn',11x,'Wcog',9x,'ngood')
          write(89,9888)
9888      format(4x,'ncell',1x,'z(m) middle',4x,'phini(deg)',
     *           5x,'phmid(deg)')
         endif
         z=0.
c ---- iterations over steps xl
         nsp1=ns+1
         do n=1,nsp1
           z=z+hl
           if(z.gt.cl) then
            zl=z-hl
            tlgth=tlgth+zl
            phfin=tref*fh+phini
            phfin=phfin*180./pi
ccc            phfin=phref*180./pi
            if(ncell.eq.1) then
              rfqdmp(ncell,1)=zl
            else
              rfqdmp(ncell,1)=zl+rfqdmp(ncell-1,1)
            endif
            rfqdmp(ncell,2)=phdep
         write(70,9999)ncell,zl,tlgth,phdep,phfin,wref,wcog,ngood
9999         format(2x,i5,6(3x,e12.5),3x,i6)
             go to 9527
           endif
c --- change of reference over the half step hl
           tref=tref+hl/(bref*vlm)
           if(itvol) ttvols=tref
           phref=tref*fh+phini
           skz=sin(cay*z)
           ckz=cos(cay*z)
c ---  synchronous particle
c        change of energy over the step xl
           sp=sin(phref)
c standard accelerating cell
           if(itype(ncell).eq.0) then
             dwref=.5*qst*cay*av*skz*sp*xl
c print in file rfq_list1.data
            if(n.eq.1) phdep=rfq4(ncell)
           endif
c T-cell (type = 1)
           if(itype(ncell).eq.1) then
            skz3=sin(3.*cay*z)
            ckz3=cos(3.*cay*z)
            dwref=0.5*qst*cay*(a10v*skz+3.*a31v*skz3)*sp*xl
c print in file rfq_list1.data
            if(n.eq.1) phdep=rfq4(ncell)
           endif
c E-cell (type = 2)
           if(itype(ncell).eq.2) then
            skz3=sin(3.*cay*z)
            ckz3=cos(3.*cay*z)
            dwref=0.5*qst*cay*(a10v*ckz+3.*a31v*ckz3)*sp*xl
c print in file rfq_list1.data
            if(n.eq.1) phdep=rfq4(ncell)
           endif
c M-cell (type = 3)
           if(itype(ncell).eq.3) then
             dwref=0.
c print in file rfq_list1.data
             if(n.eq.1) phdep=rfq4(ncell)
           endif
c Fringe-field region (itype = 4)
           if(itype(ncell).eq.4) then
             rtvolt=rfq10(ncell)
             av=rfq2(ncell)*rtvolt
             c3kz=cos(3.*cay*z)
             skpz=.75*(skz+sin(3.*cay*z))
             dwref=.5*qst*cay*av*skpz*sp*xl
c print in file rfq_list1.data
             if(n.eq.1) phdep=rfq4(ncell)
           endif
c********************************************
c Fringe-field region (itype = 7)
            if(itype(ncell).eq.7) then
              dwref=0.
c   rprof: current fringe field profile at the azimutal position z
              c3kz=cos(3.*cay*z)
              rproff=splinf(npf,z)
cold              write(16,*) 'FF z,radius ',z,rproff
c print in file rfq_list1.data
            if(n.eq.1) phdep=rfq4(ncell)
           endif
c********************************************
c --- RMS (type = 5)
           if(itype(ncell).eq.5) then
             dwref=0.
c print in file rfq_list1.data
             if(n.eq.1) phdep=rfq4(ncell)
           endif
c --- RMS (type = 6)
           if(itype(ncell).eq.6) then
             dwref=0.
c   rprof: current RMS profile given at z
             rprof=spline(npt,z)
ccc            write(16,4527) n,z,rprof,a31vb
ccc4527        format(2x,i4,3(2x,e12.5))
c print in file rfq_list1.data
             if(n.eq.1) phdep=rfq4(ncell)
           endif
c --- gain of energy (synchronous particle)
           wrefm=wref+0.5*dwref
           grefm=wrefm/er+1.
           brefm=sqrt(1.-1./(grefm*grefm))
           wref=wref+dwref
           gref=wref/er+1.
           bref=sqrt(1.-1./(gref*gref))
           dez=0.
           dref=0.
c ---- Beam
c         coordinates x, xp,y, yp convert in m and rad
           do ip=1,ngood
             xi=f(2,ip)*1.e-02
             xpi=f(3,ip)*1.e-03
             yi=f(4,ip)*1.e-02
             ypi=f(5,ip)*1.e-03
             ww=f(7,ip)-er
c remove from the bunch particles that have energy < 0
             if(ww.lt.0.)then
               f(8,ip)=0.
               ilost=ilost+1
               write(49,5558) ip,ncell,ww
               iflag=.true.
               go to 6525
             endif
5558         format(' particle: ',i5,' cell: ',i5,
     *           ' energy: ',e12.5)
             gi=ww/er+1.
             if(gi.lt.1.) then
                f(8,ip)=0.
                ilost=ilost+1
c needs to be on energy
                write(49,*) ' particle lost: ',i,' W: ',ww,
     *             ' MeV'
                iflag=.true.
                go to 6525
             endif
             bi=sqrt(1.-1./(gi*gi))
             bgi=bi*gi
c move particles back in the bunch if they have abs(phi) > pi
             tim=f(6,ip)+hl/(bi*vlm)
             phi=fh*(tim-tref)
             if(phi.ge.pi) then
               f(6,ip)=f(6,ip)-2.*pi/fh
               tim=f(6,ip)+hl/(bi*vlm)
               phi=fh*(tim-tref)
               if(abs(phi).gt.pi) then
                f(8,ip)=0.
                ilost=ilost+1
                write(49,5559) ip,ncell,phi,f(9,ip)
                iflag=.true.
                go to 6525
               endif
             endif
             if(phi.le.-pi) then
               f(6,ip)=f(6,ip)+2.*pi/fh
               tim=f(6,ip)+hl/(bi*vlm)
               phi=fh*(tim-tref)
               if(abs(phi).gt.pi) then
                f(8,ip)=0.
                ilost=ilost+1
                write(49,5559) ip,ncell,phi,f(9,ip)
                iflag=.true.
                go to 6525
               endif
             endif
5559           format(' particle: ',i5,' cell: ',i5,
     *               ' phi(rad): ',e12.5,' charge state ',f6.2)
             tim=f(6,ip)+hl/(bi*vlm)
             phi=phini+fh*tim
             qq=abs(f(9,ip))
             sp=sin(phi)
             cp=cos(phi)
             bav=bi
             gav=gi
             bgav=bgi
             bg=bgi
             beta=bi
             delt=0.
             amort=1.
             xm=xi+xpi*hl
             ym=yi+ypi*hl
             rm=sqrt(xm*xm+ym*ym)
c ---- particle is lost if abs(xm) > vanx or abs(ym) > vany
            if(n.ne.ns/2) then
              vanx=ablx(ncell)*z*z+bbx(ncell)*z+cbx(ncell)
              vany=ably(ncell)*z*z+bby(ncell)*z+cby(ncell)
              vanx=vanx*0.7523
              vany=vany*0.7523
              if(abs(xm).ge.vanx) f(8,ip)=0.
              if(abs(ym).ge.vany) f(8,ip)=0.
              if(f(8,ip).eq.0.) then
               ilost=ilost+1
               write(49,5556) ip,ncell,abs(xm),vanx,abs(ym),vany,rm
5556           format(' particle: ',i5,' cell: ',i5,
     *       ' abs(x) (m): ',e12.5,' x-v(m): ',e12.5,
     *      ' abs(y) (m):',e12.5,' y-v(m): ',e12.5,
     *      ' radius (m) ',e12.5)
               iflag=.true.
               go to 6525
             endif
            endif
            if(n.eq.ns/2) then
             if(rm.gt.r0) then
               f(8,ip)=0.
               ilost=ilost+1
               write(49,5557) ip,ncell,rm,r0
5557           format(' particle: ',i5,' cell: ',i5,
     *           ' radius (ptcl) (m): ',e12.5,
     *           ' radius (cell) (m):',e12.5)
               iflag=.true.
               go to 6525
              endif
             endif
             theta=0.
             signx=0.
             signy=0.
             xml=xm
             yml=ym
             if(abs(xm).gt.1.e-10) then
              theta=atan(ym/xm)
              signx=1.
              signy=1.
              if(theta.gt.0.) then
               if(xm.lt.0.) signx=-1.
               if(ym.lt.0.) signy=-1.
              endif
              if(theta.lt.0.) then
              if(xm.lt.0.) signx=-1.
              if(ym.gt.0.) signy=-1.
             endif
            endif
            if(abs(xm).le.1.e-10) then
              if(abs(ym).gt.1.e-10) then
                if(xm.ge.0..and.ym.gt.0.) theta=pi/2
                if(xm.ge.0..and.ym.lt.0.) theta=-pi/2
                if(xm.lt.0..and.ym.lt.0.) theta=pi/2
                if(xm.lt.0..and.ym.gt.0.) theta=-pi/2
              endif
            endif
            if(theta.eq.0.) then
             signx=0.
             signy=0.
            endif
c ----  standard accelerating cell
           if(itype(ncell).eq.0) then
             zrm=cay*rm
c ----- Bessel functions: I0 to I4
             bi0=1.+zrm*zrm/4.+zrm**4/64.
cold     *     +zrm**6/2304 +zrm**8/1.47456e05
cold      bi0=bi0+zrm**10/1.47456e07
cold             bi1=zrm/2.+zrm**3/16.+zrm**5/128.
             bi1=zrm/2.+zrm**3/16.
             bi1p=0.
             if(rm.ne.6.*0.) bi1p=bi1/rm
             bi3=zrm**3/48.+zrm**5/768
             bi4=zrm**4/384.
             bi4r=0.
             if(rm.gt.1.e-06) bi4r=bi4/rm
             bi5=zrm**5/3840.
             c2t=cos(2.*theta)
             s2t=sin(2.*theta)
             c1t=cos(theta)
             s1t=sin(theta)
c ---   transverse fields (cylindrical coordinates)
             erf=vorb*c2t*2.*rm+
     *         cay*(avb*bi1+a12vb*(bi3+bi5)*cos(4.*theta)/2.)*ckz
             erf=-erf/2.
             etf=vorb*s2t*2.*rm+
     *         4.*a12vb*bi4r*sin(4.*theta)*ckz
             etf=etf/2.
cold sv            evens=.false.
cold            ncf=(ncell/2)*2-ncell
             ncf=(ipari(ncell)/2)*2-ipari(ncell)
             if(ncf.eq.0) then
               erf=-erf
               etf=-etf
cold sv            evens=.true.
             endif
c  cartesian fields Ex and Ey
             ex=erf*c1t-etf*s1t
             ey=erf*s1t+etf*c1t
             ex=signx*ex
             ey=signy*ey
           endif
c ---- T-cell (Type = 1)
           if(itype(ncell).eq.1) then
c   transverse fields (cylindrical coordinates)
             zrm=cay*rm
             zrm3=zrm*3.
             bi0=1.+zrm*zrm/4.+zrm**4/64.+zrm**6/2304.
cold     *              +zrm**8/1.47456e05
             bi03=1.+zrm3*zrm3/4.+zrm3**4/64.+zrm3**6/2304.
cold     *              +zrm**8/1.47456e05
             bi1=zrm/2.+zrm**3/16.+zrm**5/384.
cold             bi1p=0.
cold             if(rm.ne.6.*0.) bi1p=bi1/rm
             bi13=zrm3/2.+zrm3**3/16.+zrm3**5/384.
cold             bi13p=0.
cold             if(rm.ne.6.*0.) bi13p=bi13/rm
             c2t=cos(2.*theta)
             s2t=sin(2.*theta)
             c1t=cos(theta)
             s1t=sin(theta)
             pavolt=rfq11(ncell)
             rpv=pavolt*rfq9(ncell)
             ncf=(ipari(ncell)/2)*2-ipari(ncell)
             if(ncf.eq.0) then
               rpv=-rpv
               erf=-rpv/(r0*r0)*c2t*rm
               erf=erf+cay/2.*(av10b*bi1*ckz+3.*a31vb*bi13*ckz3)
               etf=rpv/(r0*r0)*s2t*rm
             else
               erf=-rpv/(r0*r0)*c2t*rm
               erf=erf-cay/2.*(av10b*bi1*ckz+3.*a31vb*bi13*ckz3)
               etf=rpv/(r0*r0)*s2t*rm
             endif
             ex=erf*c1t-etf*s1t
             ey=erf*s1t+etf*c1t
             ex=signx*ex
             ey=signy*ey
           endif
c   tof over the length xl
c ---- E-cell (Type = 2)
           if(itype(ncell).eq.2) then
c ---   transverse fields (cylindrical coordinates)
             zrm=cay*rm
             zrm3=zrm*3.
             bi0=1.+zrm*zrm/4.+zrm**4/64.+zrm**6/2304.
cold     *              +zrm**8/1.47456e05
             bi03=1.+zrm3*zrm3/4.+zrm3**4/64.+zrm3**6/2304.
cold     *              +zrm**8/1.47456e05
             bi1=zrm/2.+zrm**3/16.+zrm**5/384.
cold             bi1p=0.
cold             if(rm.ne.6.*0.) bi1p=bi1/rm
             bi13=zrm3/2.+zrm3**3/16.+zrm3**5/384.
cold             bi13p=0.
cold             if(rm.ne.6.*0.) bi13p=bi13/rm
             c2t=cos(2.*theta)
             s2t=sin(2.*theta)
             c1t=cos(theta)
             s1t=sin(theta)
             pavolt=rfq11(ncell)
             rpv=pavolt*rfq9(ncell)
             erf=-rpv/(r0*r0)*c2t*rm
c*et*2015-07-07*s
c            erf=erf+0.5*cay*(av10b*bi1*skz+3.*a31vb*bi13*skz3)
             erf=erf+0.5*cay*(a10vb*bi1*skz+3.*a31vb*bi13*skz3)
c*et*2015-07-07*e
             etf=rpv/(r0*r0)*s2t*rm
c  control polarity of the cell
cold            ncf=(ncell/2)*2-ncell
            ncf=(ipari(ncell)/2)*2-ipari(ncell)
            if(ncf.eq.0) then
             erf=-erf
             etf=-etf
            endif
             ex=erf*c1t-etf*s1t
             ey=erf*s1t+etf*c1t
             ex=signx*ex
             ey=signy*ey
            endif
c ---- M-cell (Type = 3)
            if(itype(ncell).eq.3) then
             c2t=cos(2.*theta)
             s2t=sin(2.*theta)
             c1t=cos(theta)
             s1t=sin(theta)
             pavolt=rfq11(ncell)
             rpv=pavolt*rfq9(ncell)
             erf=-rpv/(r0*r0)*c2t*rm
             etf=rpv/(r0*r0)*s2t*rm
c ---- M-cell (Type = 3)
cold            if(itype(ncell).eq.3) then
cold             c2t=cos(2.*theta)
cold             s2t=sin(2.*theta)
cold             c1t=cos(theta)
cold             s1t=sin(theta)
cold             pavolt=rfq11(ncell)
cold             rpv=pavolt*rfq9(ncell)
cold             erf=-rpv/(r0*r0)*c2t*rm
cold             etf=rpv/(r0*r0)*s2t*rm
cold             c3kz=cos(3.*cay*z)
cold            erf=erf*(ckz+c3kz/3.)*0.75
cold            etf=etf*(ckz+c3kz/3.)*0.75
c *************************************************
c  control polarity of the cell
cold            ncf=(ncell/2)*2-ncell
            ncf=(ipari(ncell)/2)*2-ipari(ncell)
            if(ncf.eq.0) then
             erf=-erf
             etf=-etf
            endif
             ex=erf*c1t-etf*s1t
             ey=erf*s1t+etf*c1t
             ex=signx*ex
             ey=signy*ey
            endif
c ---- Radial matching section (type = 5)
            if(itype(ncell).eq.5) then
             c2t=cos(2.*theta)
             s2t=sin(2.*theta)
             c1t=cos(theta)
             s1t=sin(theta)
             skz3=sin(3.*cay*z)
             ckz3=cos(3.*cay*z)
             zrm=cay*rm
             b2kr=zrm*zrm/8.
             b2kr3=(9./8.)*zrm*zrm
             qzrm=0.
             erf=0.
             etf=0.
             if(rm.gt.6*0.) then
              qzrm=(b2kr*skz-(1./27.)*b2kr3*skz3)
              erf=-1./8.*a31vb*cay*cay*(skz-1/3.*skz3)*c2t
              erf=erf*rm
              etf=a31vb*qzrm*s2t/rm
             endif
             ex=erf*c1t-etf*s1t
             ey=erf*s1t+etf*c1t
             ex=signx*ex
             ey=signy*ey
            endif
c ****************************************************
            if(itype(ncell).eq.6) then
              BGFAC=GAV*BAV**2
              cc=a31vb*qq*xl*sp/(bgfac*er)
              rr1=-cc/(rprof*rprof)
              rr2=cc/(rprof*rprof)
              rr1=rr1*xm
              rr2=rr2*ym
            endif
c********************************************
c Fringe field profil from the disk (type 7)
           if(itype(ncell).eq.7) then
             BGFAC=GAV*BAV**2
             c1=.75*(ckz+c3kz/3.)
             cc=a31vb*qq*xl*sp/(bgfac*er)
             rr1=-cc/(rproff*rproff)
             rr2=cc/(rproff*rproff)
c  test the cell parity
             ncf=(ipari(ncell)/2)*2-ipari(ncell)
             if(ncf.eq.0) then
               rr1=-rr1
               rr2=-rr2
             endif
             rr1=rr1*xm
             rr2=rr2*ym
cold       if(ip.eq.25) write(16,5287)n,z,rproff,rr1/xm,rr2/ym
c5287       format(2x,i4,4(2x,e12.5))
           endif
c********************************************
C----  CHANGE ENERGY OVER the STEP XL
c ------  standard accelerating cell or R-cell (Type = 0)
           if(itype(ncell).eq.0) then
             ez=0.5*(avb*bi0+a12vb*bi4*cos(4.*theta))*skz*cay
c 2014 a12test line above replaced by line below
c             ez=0.5*(avb*bi0)*skz*cay
             dw=qq*ez*sp*xl
           endif
c -------  T-cell (Type = 1)
           if(itype(ncell).eq.1) then
             ez=0.5*cay*(av10b*skz*bi0+3.*a31vb*skz3*bi03)
             dw=qq*ez*sp*xl
           endif
c -------  E-cell (Type = 2)
           if(itype(ncell).eq.2) then
             ez=0.5*cay*(a10vb*ckz*bi0+3.*a31vb*ckz3*bi03)
             dw=qq*ez*sp*xl
           endif
c ------- M-cell (itype = 3)
           if(itype(ncell).eq.3) dw=0.
c ----    fringe field region (itype = 4)
c********************************************
           if(itype(ncell).eq.4) dw=.5*qq*cay*avb*skz*sp*xl
           if(itype(ncell).eq.7) dw=0.
c********************************************
c --------RMS (itype = 5)
           if(itype(ncell).eq.5) then
             ez=-1./16.*a31vb*cay**3*rm*rm*(ckz-ckz3)*c2t
             dw=qq*ez*sp*xl
            endif
c --------RMS (itype = 6)
           if(itype(ncell).eq.6) then
             dw=0.
           endif
c ----  WAV: energy at the middle of the element
           WAV=WW+.5*DW
           GA=WAV/ER
           if(ga.le.0.) then
             f(8,ip)=0.
             ilost=ilost+1
c needs to be on energy
             write(49,*) ' particle lost: ',i,' W: ',wav,
     *           ' MeV'
             iflag=.true.
             go to 6525
           endif
           BGAV=SQRT(GA*(2.+GA))
           GAV=1.+GA
           BAV=BGAV/GAV
c ---- energy over step xl
           WW=WW+DW
           GA=WW/ER
           GAM=1.+GA
           if(gam.le.1.) then
             f(8,ip)=0.
             ilost=ilost+1
c needs to be on energy
             write(49,*) ' particle lost: ',i,' W: ',ww,
     *           ' MeV'
             iflag=.true.
             go to 6525
           endif
           BG=SQRT(GA*(2.+GA))
           beta=sqrt(1.-1/(gam*gam))
cold jump of phase (sec) of particles (only for standard accelerating cells)
           delt=0.
cold           if(itype(ncell).eq.0) then
cold             dez=.5*qq*cay*avb*skz*sp
cold             delt=.5*(dez/er) * xl*xl/(bav**3*gav**3*vlm)
cold           endif
           amort=bgi/bg
           BGFAC=GAV*BAV**2
           cc=qq*xl*sp/(bgfac*er)
c********************************************
           if(itype(ncell).ne.4) then
              if(itype(ncell).eq.6) go to 9766
              if(itype(ncell).eq.7) go to 9766
              rr1=cc*ex
              rr2=cc*ey
9766         continue
              xpm=xpi*amort+rr1
              ypm=ypi*amort+rr2
              xf=xm+xpm*hl
              yf=ym+ypm*hl
            endif
c********************************************
c ------- Fringe-field region (itype = 4)
           if(itype(ncell).eq.4) then
c*********************************************************
c         C1 = (1/m**2) * (m) = (1/m)
c         C2 = (1/m**2) * (m) = (1/m)
c         RF1 = (MeV/(MeV*m**2)) = (1/m**2)
c         RF2 =  (MeV/MeV) * (1/m**2) = (1/m**2)
c*******************************************************
cold              RF1=QQ*VORB/ER
cold              RF2=.25*QQ*CAY*CAY*AVB/ER
cold              C1=RF1*SP*XL/BGFAC
cold              C2=RF2*CKZ*SP*XL/BGFAC
cold              C1=C1*.75*(CKZ+C3KZ/3.)
cold              C2=C2*.75*(CKZ+3.*C3KZ)
cold              RR1=-(C1+C2)
cold              RR2=(C1-C2)
c rfq11(ncell) is (1.+avolt)*VV , where VV is intervane voltage seen by particles
              pavolt=rfq11(ncell)
c rfq9(ncell) is 1.+fvolt
              rpv=pavolt*rfq9(ncell)
              rf1=qq*rpv/(r0*r0*er)
c avb is a10*rpv, but a10 is usually zero in fringe field
              rf2=.25*qq*cay*cay*avb/er
c  test cell parity
              ncf=(ipari(ncell)/2)*2-ipari(ncell)
              if(ncf.eq.0) then
               rf1=-rf1
               rf2=-rf2
              endif
              c1=rf1*sp*xl/bgfac
              c2=rf2*ckz*sp*xl/bgfac
              c1=c1*.75*(ckz+c3kz/3.)
              c2=c2*.75*(ckz+3.*c3kz)
              rr1=-(c1+c2)
              rr2=(c1-c2)
              xpm=xpi*amort+rr1*xm
              ypm=ypi*amort+rr2*ym
              xf=xm+xpm*hl
              yf=ym+ypm*hl
            endif
c      restore coordinates x, xp, y, yp in cm and mrad
            f(2,ip)=xf*100.
            f(4,ip)=yf*100.
            f(3,ip)=xpm*1000.
            f(5,ip)=ypm*1000.
c   tof over the length xl
            if(beta.eq.0.) then
              f(6,ip)=f(6,ip)+hl/(bi*vlm)
            else
              f(6,ip)=f(6,ip)+hl/(bi*vlm)+hl/(beta*vlm)
            endif
            f(7,ip)=ww+er
6525        continue
c --- end do ip (particle loop)
           enddo
c ----- reshuffle the good particles at the end of each element
           if(iflag) then
             call shuffle
             iflag=.false.
           endif
c    space charge at the middle of the cell
           if(n.eq.nsm) then
             if(ichaes) then
C      Space Charge
               iesp=.true.
               call cesp(scl)
               iesp=.false.
             endif
           endif
c --- change the reference over the half step hl
           tref=tref+hl/(bref*vlm)
           if(itvol) ttvols=tref
           vref=bref*vl
           z=z+hl
c print in file rfq_list1.phase RF in the middle of the cell
           if(n.eq.ns/2) then
            phmil=tref*fh+phini
            phmil=phmil*180./pi
            write(89,9735) ncell,z,phdep,phmil
9735        format(2x,i5,3(2x,e12.5))
            endif
c     Change  dp/p over the cell
           call disp
c      end of do n=1,nsp1 (steps in cell)
         enddo
9527     continue
c----  c.o.g of the bunch at the output of the cell
         tcog=0.
         ecog=0.
         do i =1,ngood
           tcog=tcog+f(6,i)
           ecog=ecog+f(7,i)
         enddo
         tcog=tcog/float(ngood)
         ecog=ecog/float(ngood)
         gcog=ecog/er
         bcog=sqrt(1.-1./(gcog*gcog))
         wcog=ecog-er
c---  window control relative to the energy of the c.o.g of the bunch
c ---- ifw = 0 ===> wdisp = dW/W
c ---- ifw = 1 ===> wdisp = dW (MeV)
c ----- convert wdisp in dp/p
         if(ifw.eq.0 .or. ifw.eq.10) then
           dispr=gcog*gcog*wdisp/(gcog*(gcog+1.))
         endif
         if(ifw.eq.1 .or. ifw.eq.11) then
           dispr=gcog*gcog*wdisp/(gcog*(gcog+1.)*wcog)
         endif
         iflag=.false.
         do i=1,ngood
           dese=abs(fd(i)-1.)
           if(dese.gt.dispr) then
             ilost=ilost+1
             f(8,i)=0.
             write(49,*) ' particle lost: ',i,' dp/p: ',dese,
     *             ' window :',dispr
             iflag=.true.
           endif
         enddo
         if(iflag) then
           call shuffle
c----  c.o.g of the bunch after shuffle
           tcog=0.
           ecog=0.
           do i =1,ngood
             tcog=tcog+f(6,i)
             ecog=ecog+f(7,i)
           enddo
           tcog=tcog/float(ngood)
           ecog=ecog/float(ngood)
           gcog=ecog/er
           bcog=sqrt(1.-1./(gcog*gcog))
           wcog=ecog-er
         endif
         br0=rfq7(ncell)*rfq7(ncell)
         bff=(1./er)*wavel*wavel*tdvolt/br0
         write(75,5555)ncell,rfq1(ncell),rfq2(ncell),
     *   rfq6(ncell),rfq7(ncell),rfq8(ncell)
5555     format(3x,i5,5(3x,e12.5))
cet2010s
c dphete,dav1(idav,16),dav1(idav,21) and dav1(idav,12) still to be assigned correct value
         dphete=0.
         trfprt=fh*tref*180./pi
         tcgprt=fh*tcog*180./pi
c         n2kp=int(tofprt/360.)
c         tofprt=tofprt-float(n2kp)*360.
c         if(tofprt.gt.180.) tofprt=tofprt-360.
c cavity number, transmission (%), synchronous phase (deg), time of flight (deg) (reference),
c COG relativistic beta (@ output), COG output energy (MeV), REF relativistic beta (@ output), REF output energy (MeV),
c horizontal emittance (mm.mrad, RMS normalized), vertical emittance (mm.mrad, RMS normalized),
c longitudinal emittance (RMS, ns.keV)
         trnsms=100.*float(ngood)/float(imax)
c*et*2012*Jun*s
         call cdg(1)
         encog=cog(1)
         gcog=encog/xmat
         bcog=sqrt(1.-1./(gcog*gcog))
         tcog=cog(3)
         CALL EXT2D(1)
         SURXTH=SQRT(exten(4)*exten(5)-exten(8)**2)
         SURYPH=SQRT(exten(6)*exten(7)-exten(9)**2)
         sqmdv=sqrt(exten(1)*exten(3)-exten(2)*exten(2))
         exns=bcog*surxth*10./sqrt(1.-bcog*bcog)
         eyns=bcog*suryph*10./sqrt(1.-bcog*bcog)
         emns=sqmdv*1.e12/fh
c*et*2012*Jun*e
         if(ncell.eq.1) then
           write(50,*) '# rfqparm.dmp'
           write(50,*) '# cell    Z       trans   ',
     *   'PHIs     TOF(COG)    COG        Wcog          TOF(REF)   ',
     *   '    REF        Wref       Ex,RMS,n     Ey,RMS,n     El,RMS'
           write(50,*) '#  #     (m)       (%)    ',
     *  '(deg)     (deg)      beta       (MeV)          (deg)      ',
     *  '   beta       (MeV)      (mm.mrad)    (mm.mrad)    (ns.keV)'
         endif
         write(50,7023) ncell,rfqdmp(ncell,1),trnsms,
     *   rfqdmp(ncell,2),tcgprt,bcog,wcog,trfprt,bref,wref,
     *   exns,eyns,emns
7023     format(1x,i4,1x,e12.5,1x,f6.2,1x,f7.2,1x,
     *   2(e14.7,1x,f7.5,1x,e14.7,1x),3(e12.5,1x))
cet2010e
         if(ncell.eq.nceltot) then
           write(16,179)
179        format(/,' Dynamics at the output',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
           write(16,1788) bcog,gcog,wcog,tcog*fh*180./pi,tcog
           write(16,165) bref,gref,wref,tref*fh*180./pi,tref
           if(itvol) write(16,*) '  time of flight: ',
     *               ttvols*fh*180./pi,' deg'
C----  new magnetic rigidity of the reference
           xmor=xmat*bref*gref
           boro=33.356*xmor*1.e-01/qst
           dav1(idav,4)=davtot*10.
           dav1(idav,5)=xlrfq*10.
           dav1(idav,6)=(gref-1.)*er
           dav1(idav,36)=ngood
           irfqp=.false.
         endif
C   plots
         CALL STAPL(davtot*10.)
c --- end do ncell
       enddo
       call emiprt(0)
c********************************************
c v28/04/2015
       ifcont=.false.
       WRITE(16,*) 'After RFQ, bunched beam assumed'
c********************************************
       RETURN
       END
        SUBROUTINE stripp
c************************************************************************************************
c --- solid stripper foils for 'slow' hadron particles
c************************************************************************************************
        implicit real*8 (a-h,o-z)
        parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
        COMMON /SPL/X(4000),Y(4000),S(5000),P(5000),Q(5000)
        common/strip/atm,qs,atms,ths,qop,sqst(6),anp,nqst
        common/mcs/imcs,ncstat,cstat(20)
        COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
        COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X                SPRFP(3000),SPRNG(3000),IPRF
        COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
        COMMON/QMOYEN/QMOY
        COMMON/RIGID/BORO
        common/faisc/f(10,iptsz),imax,ngood
        dimension pc(20),npcent(20)
        common/cgtof/charm(20),cgtdv(20),nbch(20),netac
c                       dimension charge(6),pcent(6),charm(6),pc(6),eoff(6)
        character*1 cr
        dimension xeast(21),yeast(21),vecx(1)
c --- Eastham curve (reduced half angle over reduced thickness)
c ---- xeast: reduced thickness
c ---- yeast: reduced half angle
        data xeast/0.0,2.0,4.0,6.0,8.0,10.0,
     1             12.0,14.0,16.0,18.0,20.0,
     2             22.0,24.0,26.0,28.0,30.0,
     3             32.0,34.0,36.0,38.0,40.0/
        data yeast/0.00,0.60,1.04,1.56,1.80,2.08,2.32,
     1             2.72,2.96,3.20,3.44,3.624,3.808,
     2             3.992,4.176,4.360,4.544,4.728,4.912,
     3             5.096,5.28/
c ---- characteristics of the stripper:
c ---- qs : atomic number
c ---- atms   : Atomic mass (amu)
c ---- ths:thickness (g/cm**2)
c allow for print out on terminal of stripper number on one and the same line
        cr=char(13)
        nstrp=1
        imcs=1
#ifdef F77
        write(6,8254) nrtre,nstrp,cr
8254    format('Transport element:',i5,
     *        '      Charge Stripper     :',i5,a1,$)
#else
        write(6,8254,ADVANCE='NO') nrtre,nstrp,cr
8254    format('Transport element:',i5,
     *        '      Charge Stripper     :',i5,a1)
#endif
        write(16,101)anp,atm
101     format('***** Projectile ',/,4x,
     *  'atomic number: ',f4.0,' atomic mass : ',f4.0)
        write(16,100)qs,atms,ths
100     format('***** Charge Stripper ',/,4x,
     *  'atomic number: ',f4.0,' atomic mass : ',f4.0,
     1     ' thickness: ',e12.5,' g/cm**2')
c ---- characteristics of particles
c --------- qp=f(9,i)
c --------- uem: atomic mass (amu)
c --------- wp=f(7,i)-xmat
c ----------------------------
c ----      ENVELOPE
        if(iprf.eq.1) call stapl(davtot*10.)
c ---- start prints in file 'short.data'
        idav=idav+1
        iitem(idav)=16
        dav1(idav,1)=qs
        dav1(idav,2)=atms
        dav1(idav,3)=ths
        dav1(idav,4)=davtot*10.
c ---- kinetic energy at the input
        wicg=0.
        do i=1,ngood
          wicg=f(7,i)+wicg
        enddo
        wicg=wicg/float(ngood)
        gcog=wicg/xmat
        bcog=sqrt(gcog*gcog-1.)/gcog
        wicg=wicg-xmat
          np=21
        denes=0.
          do i=1,np
            x(i)=xeast(i)
            y(i)=yeast(i)
          enddo
c ----  deriv2: develop the matrix for the b-splines functions
          call deriv2(np)
          len=1
          do ip=1,ngood
            wp=f(7,ip)-xmat
c --- wpatm: MeV/atm
            wpatm=wp/atm
            qp=f(9,ip)
c ---- xa: screening distance (cm)
            qsp=qs**(2./3)+qp**(2./3)
            xa=4.68165e-9/sqrt(qsp)
c ---- reduced thickness: thck
            xn=6.022e23/atms
            thck=xn*pi*xa*xa*ths
c --- reduce half angle: rtheta(rad) from Eastham curve
c ----  b-splines interpollation if thck<40, otherwise: rthet = 9.2e-02*thck + 1.6 (linear interpolation)
            if (thck.lt.40.) then
              rthet=spline(np,thck)
            else
              rthet=9.2e-02*thck+1.6
            endif
c ---- half angle of diffusion: dthet (mrad)
            zps=qs*qp/(xa*wp)
            dthet=2.88e-10*zps*rthet
c ---  angle of diffusion (mrad) scattered from M.C. separatly in xp and yp directions from a uniform law
cunif         call rlux(vec,len)
cunif         rx=(2.*vec-1.)*dthet*2.
cunif         f(3,ip)=f(3,ip)+rx
cunif         call rlux(vec,len)
cunif         ry=(2.*vec-1.)*dthet*2.
cunif         f(5,ip)=f(5,ip)+ry
c ---  angle of diffusion (mrad) scattered from M.C. separatly in xp and yp directions from a Gaussian law
            ax=f(3,ip)
            sm=dthet*2.
            call randga(len,sm,ax,vx)
            f(3,ip)=vx
            ay=f(5,ip)
            call randga(len,sm,ay,vy)
            f(5,ip)=vy
c ---  closest distance of approach: xb(cm)
            aps=(atms+atm)/(atms*atm)
            xb=1.44e-13*aps*qs*qp/sqrt(wpatm)
            if(ip.eq.1) then
              alpha=1.576e-02*qp*qs/sqrt(wpatm)
              write(16,5830)xa,thck,rthet,xb,dthet,alpha
5830          format(4x,'screening distance: ',e12.5,' cm',/,
     *        4x,'reduced thickness: ',e12.5,' reduced half angle: ',
     *        e12.5,' rad ',/,4x,'closest distance of approach: ',
     *        e12.5,' cm',/,4x,
     *        'half angle of diffusion: ',e12.5,' mrad',/,4x,
     *        'Bohr parameter: ',e12.5)
            endif
c ----  loss of energy per scatter in the stripper (eq.16)
            wapc=4.*atm*atms/((atm+atms)**2)
            dene=wapc*xb*xb*rthet*rthet/(xa*xa)*wp
            denes=denes+dene
            f(7,ip)=f(7,ip)-dene
c ---- change the electric charge state of the particle over the foil
            f(9,ip)=qop
          enddo
          denes=denes/float(ngood)
        write(16,*)'dE(MeV) (Eastham): ',dene,denes
        if(qs.eq.6. .or. qs.eq.3.) then
c ---- change the electric charge state of the particles (carbon foil case)
c ---- Based on E.Baron et al, NIM A328 (1993) p.177-182
c calculate dX', dY', dZ'
          fksi=0.1535375*(qs/atms)*anp*ths/(bcog*bcog)
          des=0.5*0.001*(1.866+1.57*log(wicg/atm))*(anp/atm)
          des=des*sqrt(1000000.*ths*qs/atms)
          write(16,*)'dE(MeV) stripping: ',des,des*atm,wicg,wicg/atm
          write(16,*)'dE(MeV) ksi: ',fksi
c calculate the charge state distribution
          qbar=anp*(1.-exp(-83.275*bcog/(anp**0.447)))
          qavg=qbar*(1.-exp(-12.905+0.2124*anp-0.00122*anp*anp))
          yy=qbar/anp
          stdv=sqrt(qbar*(0.07535+0.19*yy-0.2654*yy*yy))
          con=1./(stdv*sqrt(2.*pi))
          fact=-1./(2.*stdv*stdv)
          pcsum=0.
          numchs=0
          qst=int(qavg)
c ---- only take charge states that have more than thresh % of the particles
          thresh=100./float(ngood)
          write(16,7830) thresh
7830      format(4x,'Carbon foil stripper. Charge state distribution',
     *    ' based on E.Baron et al, NIM A328 (1993) p.177-182',/,4x,
     *    'Threshhold for cutoff of the distribution: ',f12.7,' %')
          do i=1,100
c            pc(i)=100.*con*exp(fact*(sqst(i)-qavg)*(sqst(i)-qavg))
            pcent=100.*con*exp(fact*(float(i)-qavg)*(float(i)-qavg))
            if(pcent.gt.thresh) then
              numchs=numchs+1
              sqst(numchs)=float(i)
              pc(numchs)=pcent
              pcsum=pcsum+pcent
            endif
          enddo
          nqst=numchs
          f(9,1)=float(int(qavg))
          qavg=f(9,1)
          write(16,111) nqst,int(qavg)
          netac=nqst
111       format(4x,'Number of charge states after the foil ',I2,/,
     *           4x,'Average charge state: ',i3)
C FIRST TRAJECTORY HAS AVERAGE CHARGE STATE
          NTOT=0
          DO I=1,numchs
            NPCENT(I)=int(pc(i)*float(ngood)/100.)
            ntot=ntot+NPCENT(I)
          enddo
c*temp*2012 : before correction
          DO I=1,numchs
            write(16,122) sqst(i),npcent(i),pc(i)
          enddo
C add missing number of particles to central charge state
          ncstat=numchs
          write(16,*) '   ntot,ngood=',ntot,ngood,' particles'
          DO I=1,numchs
            if(int(sqst(i)) .eq. int(f(9,1))) then
              NPCENT(I)=NPCENT(I)+ngood-ntot
            endif
c*et*2012-May-02*s
            cstat(i)=float(int(sqst(i)))
            charm(i)=cstat(i)
c*et*2012-May-02*e
            write(16,122) sqst(i),npcent(i),pc(i)
          enddo
122       format(4x,'Charge=',f3.0,' with ',i5,' particles',
     *       ' or ',f12.7,' %')
          len=1
          i=2
          DO
            if(i.gt.ngood) exit
            DO
              call rlux(vecx,len)
              XARPHA=VECX(1)
              NCOUNT=int(XARPHA*(float(numchs)+0.5))
              IF(NCOUNT.gt.0) exit
            ENDDO
            IF(NPCENT(NCOUNT).gt.0) then
              NPCENT(NCOUNT)=NPCENT(NCOUNT)-1
c*et*2011-Nov-2 bug fixed in next line
c              f(9,i)=sqst(NCOUNT)
              f(9,i)=float(int(sqst(NCOUNT)))
              i=i+1
            ENDIF
          ENDDO
        endif
c*et*2012-May-02*s
        call cogetc
c*et*2012-May-02*e
c ----  Change  dp/p over the stripper
        call disp
c
c ----  the new reference is the cog
        qcg=0.
        wcg=0.
        do i=1,ngood
          wcg=f(7,i)+wcg
          qcg=f(9,i)+qcg
        enddo
        wcg=wcg/float(ngood)
        gcg=wcg/xmat
        qcg=qcg/float(ngood)
        qmoy=qcg
        bref=sqrt(gcg*gcg-1.)/gcg
        vref=bref*vl
C----  new magnetic rigidity
        xmor=xmat*bref*gcg
        boro=33.356*xmor*1.e-01/qcg
        diff=(wcg-xmat)-wicg
        dav1(idav,5)=qavg
        dav1(idav,6)=diff
        dav1(idav,36)=ngood
        write(16,5420) wicg,-diff
5420    format(4x,'energy of cog: at entrance: ',e12.5,' MeV',/,
     *  4x,'energy loss of cog: ',e12.5,' MeV')
C   plots
        CALL STAPL(davtot*10.)
        call emiprt(0)
        return
        end
        SUBROUTINE randga (len,s,am,v)
c generateur aleatoire selon une loi normale
c         s : ecart-type de la distribution
c         am: moyenne de la distribution
c         v : nombre aléatoire selon la loi normale
        implicit real*8 (a-h,o-z)
        dimension vecx(1)
        a=0.
        do i=1,24
          call rlux(vecx,len)
          y=vecx(1)
          a=a+y
        enddo
        v=(a-12.)*s+am
        return
        end
       SUBROUTINE qelec(volt,xlqua,rs)
c   ............................................................................
C         electrostatic quadrupole
c         space charge computation at the middle of the lens
c   ............................................................................
c       VOLT: electric voltage at pole tip (kV)
C       XLQUA: effective length  (cm)
C       RS: radial distance of pole tip from axis (cm)
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       COMMON/FAISC/F(10,iptsz),IMAX,NGOOD
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       LOGICAL IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       LOGICAL SHIFT
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rander/ialin
       LOGICAL IALIN
       common/qskew/qtwist,iqrand,itwist,iaqu
       LOGICAL ITWIST
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/QSEX/L,KQ2,KS2
       common/tofev/ttvols
       COMMON/ITVOLE/ITVOL,IMAMIN
       logical itvol,imamin
       REAL*8 L,KQ2,KS2
       dimension trans(1)
       character*1 cr
       ilost=0
C      statistics
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
       write(16,*) ' ***QUADRUPOLE (electrostatic)***'
       write(16,*)
       fcpi=fh*180./pi
       if(itvol) write(16,10) ttvols*fcpi,davtot
10     FORMAT(' ** tof for adjustments at input: ',e12.5,
     *        ' deg at position: ',e12.5,' cm in lattice')
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,advance='no') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
c  if itwist=.true. skews the quadrupole before misalignments (only if abs(volt) gt 1.e-13)
c     len = 1 level 1 in rlux routine
       len=1
       sqtwist=0.
       if(abs(volt).gt.1.e-13) then
         if(itwist) then
           if(iqrand.eq.0) then
             qtwrad=qtwist
             sqtwist=qtwrad
             call zrotat(qtwrad)
           else
             rdcf=.5
             call rlux(trans,len)
             if(trans(1).le.rdcf) sign=-1.
             if(trans(1).gt.rdcf) sign=1.
             call rlux(trans,len)
             qtwrad=qtwist*sign*trans(1)
             sqtwist=qtwrad
             call zrotat(qtwrad)
           endif
         endif
c   misalignments
         if(ialin) call randali
       endif
       idav=idav+1
       iitem(idav)=18
       dav1(idav,1)=xlqua*10.
       davtot=davtot+xlqua
       dav1(idav,4)=davtot*10.
       fh0=fh/vl
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       bet=sqrt(gpa*gpa-1.)/gpa
       l=xlqua
       eni=gpa*xmat
       qi=qst
c      electric rigidity (kV)
       rigid=eni*bet*bet/qi *1.e03
c     electric field gradient: xgrad (kV/cm-2)
        xgrad=2.*volt/rs**2
        dav1(idav,2)=volt
        dav1(idav,3)=rigid
        dav1(idav,5)=xgrad
c   qk2 (cm-2)
        kq2=xgrad/rigid
        dav1(idav,6)=kq2
        dav1(idav,7)=rs*10.
       write(16,100) xlqua,rs,volt,kq2,xgrad,rigid
100   format(' LENGTH = ',e12.5,' cm   APERTURE RADIUS= ',e12.5,' cm',
     *  /,' VOLTAGE = ',e12.5,' kV  K2 = ',e12.5,' cm-2  GRADIENT = ',
     *  e12.5,' kV/cm2',/,' RIGIDITY = ',e12.5,' kV',/)
        call clear
        call elqua
c    print out the transport matrix (cog)
        call matrix
C      first half quadrupole
       l=xlqua/2.
       do ii=1,ngood
        call clear
        gpa=f(7,ii)/xmat
        bet=sqrt(gpa*gpa-1.)/gpa
        qi=f(9,ii)
c      electric rigidity (kV)
        rigi=f(7,ii)*bet*bet/qi *1.e03
c   qk2 (cm-2)
        kq2=xgrad/rigi
        call elqua
        call cobeam(ii,l)
       enddo
C      space charge computations (if l >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle of the lens'
           call cesp(xlqua)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
c Test window after the first half quadrupole (after s.c. computations)
       call cogetc
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- window control
       tref=tref+xlqua/(2.*vref)
       call reject(nlost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c    second half quadrupole
       do ii=1,ngood
         call clear
         gpa=f(7,ii)/xmat
         bet=sqrt(gpa*gpa-1.)/gpa
         qi=f(9,ii)
         rigi=f(7,ii)*bet*bet/qi *1.e03
         kq2=xgrad/rigi
         call elqua
         call cobeam(ii,l)
       enddo
c Test window after the second half quadrupole
       call cogetc
       tref=tref+xlqua/(2.*vref)
       call reject(ilost)
       ilost=ilost+nlost
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c  Change the t.o.f
       if(itvol) ttvols=tref
C  envelope
       call stapl(davtot*10.)
       tcog=0.
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       if(itvol) then
        write(16,11) ttvols*fcpi,davtot,tref*fcpi,tcog*fcpi
11      FORMAT(' ** tof for adjustments: ',e12.5,
     *        ' deg at position: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg tof of the cog: ',
     * e12.5,' deg')
       else
        write(16,12) tref*fcpi,tcog*fcpi
12      FORMAT(' ** tof of the reference: ',e12.5,
     *        ' deg tof of the cog: ',e12.5,' deg')
       endif
       dav1(idav,36)=ngood
       write(16,*)' particles lost :',ilost
c  returns coordinates to the initial orientation
       if(itwist) then
        if(abs(volt).gt.1.e-13) then
         qtwrad=-sqtwist
         call zrotat(qtwrad)
        endif
       endif
       if(iemgrw) call emiprt(0)
c  envelope
       call stapl(davtot*10.)
       return
       end
       SUBROUTINE qfk (ityqu,arg,xlqua,rs)
c   ............................................................................
C         electrostatic or magnetic quadrupole based on the strength K2
c         space charge computation at the middle of the lens
c   ............................................................................
c       ITYQU: ITYQU = 0 electric quadrupole, otherwise magnetic quadrupole
c       ARG: strength K (cm-2)
C       XLQUA: effective length  (cm)
C       RS: radial distance of pole tip from the axis (cm)
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/DYN/TREF,VREF
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/RIGID/BORO
       COMMON/FAISC/F(10,iptsz),IMAX,NGOOD
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON/DCSPA/IESP
       LOGICAL IESP
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       LOGICAL ICHAES
       COMMON/TAPES/IN,IFILE,META
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/SHIF/DTIPH,SHIFT
       LOGICAL SHIFT
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/rander/ialin
       LOGICAL IALIN
       common/qskew/qtwist,iqrand,itwist,iaqu
       LOGICAL ITWIST
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/QSEX/L,KQ2,KS2
       common/tofev/ttvols
       COMMON/ITVOLE/ITVOL,IMAMIN
       common/qfkd/ityq
       logical itvol,imamin,ityq
       REAL*8 L,KQ2,KS2
       dimension trans(1)
       character*1 cr
       ilost=0
       if(ityqu.eq.0) ityq=.true.
       if(ityqu.ne.0) ityq=.false.
       if(ityq) then
         write(16,*) ' ***QUADRUPOLE (electrostatic)***'
       else
         write(16,*) ' ***QUADRUPOLE (magnetic)***'
       endif
C      statistics
       IF(IPRF.EQ.1) CALL STAPL(davtot*10.)
       xfqu=arg
       fcpi=fh*180./pi
       if(itvol) write(16,10) ttvols*fcpi,davtot
10     FORMAT(' ** tof at input: ',e12.5,
     *        ' deg  position in the lattice: ',e12.5,' cm ')
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
c  if itwist=.true. skews the quadrupole before misalignments (only if abs(arg) gt 1.e-13)
c     len = 1 level 1 in rlux routine
       len=1
       sqtwist=0.
       if(abs(xfqu).gt.1.e-13) then
         if(itwist) then
           if(iqrand.eq.0) then
             qtwrad=qtwist
             sqtwist=qtwrad
             call zrotat(qtwrad)
           else
             rdcf=.5
             call rlux(trans,len)
             if(trans(1).le.rdcf) sign=-1.
             if(trans(1).gt.rdcf) sign=1.
             call rlux(trans,len)
             qtwrad=qtwist*sign*trans(1)
             sqtwist=qtwrad
             call zrotat(qtwrad)
           endif
         endif
c     misalignments
         if(ialin) call randali
       endif
       idav=idav+1
       iitem(idav)=19
       dav1(idav,1)=xlqua*10.
       davtot=davtot+xlqua
       dav1(idav,4)=davtot*10.
       fh0=fh/vl
       gpa=0.
       do ii=1,ngood
        gpa=gpa+f(7,ii)/xmat
       enddo
       gpa=gpa/float(ngood)
       bet=sqrt(gpa*gpa-1.)/gpa
       l=xlqua
       eni=gpa*xmat
       qi=qst
c   electric quadrupole
       if(ityq) then
c  rigid: electric momentum (kV)
        rigid=eni*bet*bet/qi * 1.e03
c    compute the electric voltage VOLT (kV) from K(cm-2)
        volt=xfqu*rs*rs*rigid
        volt=volt/2.
c   kq2 (cm-2)
        kq2=xfqu
c     electric field gradient: xgrad (kV/cm2)
        xgrad=kq2*rigid
        dav1(idav,2)=volt
        dav1(idav,3)=rigid
        dav1(idav,5)=xgrad
        dav1(idav,6)=xfqu
        dav1(idav,7)=rs*10.
        call clear
        call elqua
       write(16,100) xlqua,rs,volt,xfqu,xgrad,rigid
100   format(' LENGTH = ',e12.5,' cm   APERTURE RADIUS= ',e12.5,' cm',
     * /,' VOLTAGE = ',e12.5,' kV  K2 = ',e12.5,' cm-2  GRADIENT = ',
     *  e12.5,' kV/cm2',/,' MOMENTUM = ',e12.5,' kV',/)
c    print the transport matrix (of the cog)
        call matrix
       endif
c     magnetic quadrupole
       if(.not.ityq) then
        xmco=xmat*bet*gpa
c       rigid: magnetic rigidity (kG.cm)
        rigid=33.356*xmco*1.e-01/qst
        kq2=xfqu
c   bgrad: gradient (kG/cm)
        bgrad=kq2*rigid
c       bgaus: field (kG)
        bgaus=bgrad*rs
        dav1(idav,2)=bgaus
        dav1(idav,3)=rigid
        dav1(idav,5)=bgrad
        dav1(idav,6)=xfqu
        dav1(idav,7)=rs*10.
        call clear
        call elqua
       write(16,3300) xlqua,rs,bgaus,xfqu,bgrad,rigid
3300   format(' LENGTH = ',e12.5,' cm   APERTURE RADIUS= ',e12.5,' cm',
     * /,' FIELD = ',e12.5,' kG  K2 = ',e12.5,' cm-2  GRADIENT = ',
     *  e12.5,' kG/cm',/,' MOMENTUM = ',e12.5,' kG.cm',/)
        call matrix
       endif
C      first half quadrupole
       l=xlqua/2.
       do ii=1,ngood
        call clear
        gpa=f(7,ii)/xmat
        bet=sqrt(gpa*gpa-1.)/gpa
        qi=f(9,ii)
c electric quadrupole
        if(ityq) then
c   rigi: momentum (kV)
         rigi=f(7,ii)*bet*bet/qi * 1.e03
c   qk2 (cm-2)
         kq2=xgrad/rigi
         call elqua
         call cobeam(ii,l)
        endif
c magnetic quadrupole
        if(.not.ityq) then
         xmco=xmat*bet*gpa
c    rigi: momentum (kG.cm)
         rigi=3.3356*xmco/f(9,ii)
c    kq2 (cm-2)
         kq2=bgrad/rigi
         call elqua
         call cobeam(ii,l)
        endif
       enddo
C      space charge computations (if l >0)
       if(ichaes.and.l.gt.0.) then
         if(sce10.eq.1.or.sce10.eq.3.) then
           iesp=.true.
           write(16,*) 'space charge at the middle  '
           call cesp(xlqua)
           iesp=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
           call disp
         endif
       endif
c Test window after the first half quadrupole (after s.c. computations)
       call cogetc
       bcour=0.
       do i=1,ngood
         gpai=f(7,i)/xmat
         bcour=sqrt(1.-1./(gpai*gpai))+bcour
       enddo
       bcour=bcour/float(ngood)
       gcour=1./sqrt(1.-bcour*bcour)
       wcg=(gcour-1.)*xmat
c ----- convert window control
       tref=tref+xlqua/(2.*vref)
       call reject(nlost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c    second half quadrupole
       do ii=1,ngood
        call clear
        gpa=f(7,ii)/xmat
        bet=sqrt(gpa*gpa-1.)/gpa
        qi=f(9,ii)
c electric quadrupole
        if(ityq) then
c     rigi (kV)
         rigi=f(7,ii)*bet*bet/qi * 1.e03
c   qk2 (cm-2)
         kq2=xgrad/rigi
         call elqua
         call cobeam(ii,l)
        endif
c magnetic quadrupole
        if(.not.ityq) then
         xmco=xmat*bet*gpa
c    rigi: momentum (kG.cm)
         rigi=3.3356*xmco/f(9,ii)
c    kq2 (cm-2)
         kq2=bgrad/rigi
         call elqua
         call cobeam(ii,l)
        endif
       enddo
c Test window after the second half quadrupole
       call cogetc
       tref=tref+xlqua/(2.*vref)
       call reject(ilost)
       ilost=ilost+nlost
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
c  Change the t.o.f
       if(itvol) ttvols=tref
C  envelope
       call stapl(davtot*10.)
       tcog=0.
       do i=1,ngood
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       if(itvol) then
        write(16,11) ttvols*fcpi,davtot,tref*fcpi,tcog*fcpi
11      FORMAT(' ** tof for adjustments: ',e12.5,
     *        ' deg at position: ',e12.5,' cm in the lattice',
     * /,3x,'tof of the reference: ',e12.5,' deg tof of the cog: ',
     * e12.5,' deg')
       else
        write(16,12) tref*fcpi,tcog*fcpi
12      FORMAT(' ** tof of the reference: ',e12.5,
     *        ' deg tof of the cog: ',e12.5,' deg')
       endif
       dav1(idav,36)=ngood
       write(16,*)' particles lost :',ilost
c  returns coordinates to the initial orientation
       if(itwist) then
        if(abs(xfqu).gt.1.e-13) then
         qtwrad=-sqtwist
         call zrotat(qtwrad)
        endif
       endif
       if(iemgrw) call emiprt(0)
c  envelope
       call stapl(davtot*10.)
       return
       end
       SUBROUTINE cavnum
c   ............................................................................
c   numerical computations of the dynamic in cavities or gap
c    the field can be read on the disk on the form: (z,E(z)
c    or it can be read in the command list on the form of a Fourier series expansion
c   ............................................................................
       implicit real*8 (a-h,o-z)
C      ****************************************************
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFS/DYNT(MAXCELL),DYNTP(MAXCELL),DYNTPP(MAXCELL),
     *   DYNE0(MAXCELL),DYNPH(MAXCELL),DYNLG(MAXCELL),FHPAR,NC
       COMMON/POSI/IST
       COMMON/MIDGAP/ENMIL,VAPMI
       COMMON/AZMTCH/DLG,XMCPH,XMCE
       COMMON/AZLIST/ICONT,IPRIN
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/FUNC/A(200),YLG,ATTE,NCEL,NHARM
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
C      TRANSIT TIME COEFFICIENTS
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/cavnum1/xnh,xpas,ffield,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/rfield/ifield
       COMMON/QMOYEN/QMOY
       COMMON/RIGID/BORO
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/compt1/ndtl,ncavmc,ncavnm
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/TAPES/IN,IFILE,META
       common/etcom/cog(8),exten(17),fd(iptsz)
       common/speda/dave,idave
       COMMON/SHIF/DTIPH,SHIFT
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/DCSPA/IESP
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/appel/irstay,ilost,iavp,ispcel
       common/femt/iemgrw,iemqesg
       common/mode/eflvl,rflvl
       common/aerp/vphase,vfield,ierpf
       common/tofev/ttvols
       common/elec/jelec
       common/step/istep
c  ****************************************************************
c   SV 09/09/2015
       common/xposi/xpost(10),xlce(2),xpax(2),iscx(2)
       common/kcell/avrg(15)
c ****************************************************************
C*ET*2015*Jul*S
c********************************************
c reservation for TRACE3D related stuff
       common/trace3d/trace3h(100),trace3t(maxcell1),tif,kt3h,kt3t,fid
       common/t3dfld/fldctr(15),zend(15),T3D
c       logical T3D,xiset
       logical T3D
c       character*128 trace3h,trace3t,tif,tifa,tifb
       character*128 trace3h,trace3t,tif
c********************************************
       common/xitrd3/bcour1(15),bcour2(15),tr3dw(15),tr3ph(15),
     *  t0tr3d(15)
C*ET*2015*Jul*E
       logical iesp,ichaes,irstay,iavp,ispcel,ifield,iemgrw
       LOGICAL SHIFT,CHASIT,ITVOL,IMAMIN,DAVE,JELEC
       CHARACTER*1 cr
C************************************************************
C    XESLN : NEGATIVE LENGHT OF THE DRIFT FOLLOWING THE GAP
C    IF XESLN N.E.0 THEN THE CHARGE SPACE EFFECT IMPLIES THE
C    LENGTH (YLG-XESLN)
       NRRES=NRRES+1
       ncavnm=ncavnm+1
c allow for print out on terminal of gap# on one and the same line
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
       WRITE(16,*)' CAVITY N :',NRRES
       ilost=0
       aqst=abs(qst)
       qmoy=aqst
c --- the frequency fh may be changed with delfh
       oldfh=fh
C        IDUM : dummy variable (indicate in input file the number of the cavity)
         READ (IN,*) IDUM
c
c      ielec = 0 acceleration for non relativistic particles with erest < 1 MeV (e.g. non-relativistic electrons)
c                otherwise acceleration for hadrons
c      dphase (deg): phase offset from the phase crest (giving the maximum of energy gain)
c      FFIELD : in percent;
c                    (electric field)=(initial electric field)*(1.+FFIELD/100)
        jelec=.false.
        read (in,*) dielec,dphase,ffield,istep,ielec
c *********************************************************
c  SV  28/10/2015
        if(ncavnm.eq.1) write(13,990)
990     format(3x,'ncav',2x,'ncell',2x,'pos I(cm)',4x,'pos S(cm)',6x,
     *'L cell(cm)',5x,'W(Mev)',7x,'dw(Mev)',6x,'TOF(dg)',7x,
     * 'Ph RF(dg)',5x,'avrg pos(cm)')
        if(ncavnm.eq.1) write(13,995)
995     format(3x,'ncell',3x,'E0TL(Mev)',4x,'T(k) Mev/q',4x,
     * 'S(k) Mev/q',5x,'dW(Mev)',6x,'PHASE(dg)',5x,
     * 'Ph crest(dg)')
        write(13,*)
c ***************************************************************
        if(ielec.eq.0) jelec=.true.
        ffield=1.+ffield/100.
        if(ffield.eq.0.) ffield=1.E-12
       if(ifield) then
c --- The field is read on the disk in the form:
c          z   (cm)   E(z) MV/cm
c     fhc: frequency of the cavity (Hertz) (read in the file 'field.txt' )
         fh=fhc*2.*pi
         ncel=ncell
         atte=att
         ye0=atte
c    flength : length of the field (cm)
         flength=xspl(npt)-xspl(1)
         ylg=flength
       else
c --- The field is read on the input list in the form of a Fourier series
         ncell=ncel
         oldfh=fh
c     atte: factor acting on the amplitude of the field (read in the input list)
         ye0=atte
       endif
       dphete=dphase
       if(itvol.and.imamin) then
c ---  adjustment of the phase offset w.r.t. the t.o.f.  (deg)
        ottvol=fh*ttvols*180./pi
        attvol=ottvol
        xkpi=ottvol/360.
        ixkpi=int(xkpi)
        xkpi=(xkpi-float(ixkpi))*360.
        dphase=dphase-xkpi
       endif
       WRITE(16,150)FH/(2.*pi),YLG,ATTE,ffield,NCEL,istep
150    FORMAT(4X,'FREQUENCY :',E12.5,' Hertz',/,4x,
     x        'FIELD LENGTH :',e12.5,' cm',/,4x,
     x        'FIELD FACTOR (UNITS CONVERSION) :',e12.5,/,4x,
     x        'FIELD FACTOR (ATTENUATION)      :',E12.5,/,4x,
     X        'FIELD DIVIDED IN: ',I4,' CELLS  STEPS BY CELL ',i5)
       if(.not.imamin) write(16,*) '   PHASE OFFSET: ',dphete,' DEG'
       if(imamin) write(16,1501)dphete,DPHASE,xkpi
1501   format(4x,
     x        'PHASE OFFSET (before adjustment): ',e12.5,' deg',/,4x,
     x        'PHASE OFFSET (after adjustment): ',e12.5,' deg',/,4x,
     *        'ADJUSTMENT ON THE PHASE OFFSET: ',e12.5,' deg')
       fh0=fh/vl
       BEREF=VREF/VL
c --- ttvol: time of flight at entrance (sec)
       ttvol=0.
       if(itvol)ttvol=ttvols*fh
c  start file 'short.data'
c --- dav1(idav,3)=0: the particle reference and the cog coincide at the input
c --- dav1(idav,3)=1: the particle reference and the cog are independent
       dav1(idav,3)=0.
       idav=idav+1
       iitem(idav)=1
       dav1(idav,1)=ylg*10.
       dav1(idav,2)=ye0*100.
       davtot=davtot+ylg
       dav1(idav,24)=davtot*10.
       dav1(idav,40)=fh
       IF(IPRF.EQ.1) CALL STAPL(dav1(idav,24))
       iarg=1
       call cdg(iarg)
       enold=cog(1)
       encog=enold
       gcog=enold/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       IF(SHIFT) THEN
c --- the reference particle and the cog are independent
         BEREF=VREF/VL
         GAMREF=1./SQRT(1.-(BEREF*BEREF))
         ENREF=XMAT*GAMREF
         TREFDG=TREF*FH*180./PI
         dav1(idav,3)=1.
       ELSE
c --- the reference particle and the cog are coinciding
         beref=bcog
         vref=bcog*vl
         tref=tcog
         gamref=gcog
         enref=cog(1)
         dav1(idav,3)=0.
       ENDIF
       if(dav1(idav,3).eq.1.) write(16,*)
     *   ' ****reference and cog evolve independently'
       if(dav1(idav,3).eq.0.) write(16,*)
     *   ' **** the reference is the cog '
       WRITE(16,178)
178    FORMAT(/,' Dynamics at the input',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       WRITE(16,1788) bcog,gcog,encog-xmat,tcog*fh*180./pi,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       ENRPRIN=ENREF-XMAT
       WRITE(16,165) beref,gamref,ENRPRIN,tref*fh*180./pi,tref
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       iprint=0
       call statis
       XK1=FH/VREF
C ---  transit time factors TK and SK  based on the velocity at the entrance
       TK=TTA0(BEREF)/2. * FFIELD
       SK=TSB0(BEREF)/2. * FFIELD
C --- prediction of PCREST (phase of RF giving the maximum of energy gain in the cavity)
       PCREST=ATAN(-SK/TK)
       DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       IF(DDWC.LT.0.) PCREST=PCREST+PI
c --- first correction of pcrest based on an average beta
       call phcrest(pcrest,ylg,ncell,zcrest)
c    avbet: average value of beta
       avbet=fh/(zcrest*vl)
       tk=tta0(avbet)/2. * ffield
       sk=tsb0(avbet)/2. * ffield
       PCREST=ATAN(-SK/TK)
       DDWC=AQST*(TK*COS(PCREST)-SK*SIN(PCREST))
       IF(DDWC.LT.0.) PCREST=PCREST+PI
c  final correction of the phase crest PCREST
       call phcrest1(pcrest,ylg,ncell)
c   phase at the entrance of the cavity
       DPHASE=DPHASE*PI/180.
       phi0=pcrest+dphase+ttvol
c --- compute energy and TOF of reference
c ******************************************
       call dwref(phi0,gams,ts)
c ********************************************
       enrs=gams*xmat
       ddw=enrs-enref
       if(t3d) then
         do klm=1,ncel

           trphase=tr3ph(klm)
           tre0tl=t0tr3d(klm)/aqst
cc split cell into drif-gap-drift
cc debugging to be done: include q/m !!! protons assumed for now
           kt3t=kt3t+1
c drift
           write(tif,6001) kt3t,kt3t,10.*fldctr(klm)
6001       format('  nt(',i4,')= 1, a(1,',i4,')=',f12.6)
           trace3t(kt3t)=tif
c
           kt3t=kt3t+1
c gap
           write(tif,6005) kt3t,kt3t,tre0tl,trphase,fid
6005       format('  nt(',i4,')=10, a(1,',i4,')=',f9.5,' , ',f9.2,
     &        ', 0., 1.,',f5.3,',')
           trace3t(kt3t)=tif
c
           kt3t=kt3t+1
c drift
           write(tif,6001) kt3t,kt3t,10.*zend(klm)
           trace3t(kt3t)=tif
         enddo
       endif
c end trace3d drift-gap-drift
c
       trefs=ts+tref
       bets=sqrt(gams*gams-1.)/gams
cold       bets=b5
       TREDG=fh*TREFS *180./PI
c  dynamic of the bunch
       call bcnum(phi0,ylg,ncell)
C -----  window control
         gcg=0.
         do i=1,ngood
           gcg=gcg+f(7,i)/xmat
         enddo
         gcg=gcg/float(ngood)
         bcg=sqrt(1.-1./(gcg*gcg))
         wcg=(gcg-1.)*xmat
         call cogetc
comment  twind=0.
         do i=1,ngood
           gpai=f(7,i)/xmat
           if(gpai.lt.1.) gpai=1.
           bpai=sqrt(1.-1./(gpai*gpai))
           fd(i)=bpai/bcg * gpai/gcg
comment    twind=twind+f(6,i)
         enddo
comment  twind=twind/float(ngood)
c ---- window control
       tref=trefs
       vref=bets*vl
       call reject(ilost)
c  Reshuffles f(i,j) array after window (now done in 'reject')
c         call shuffle
       write(16,*) ' PARAMETERS RELATING TO THE REFERENCE PARTICLE '
       write(16,*) '***********************************************'
       write(16,*) ' ENERGY GAIN(MeV) ',DDW,' TOF ',tredg,' DEG'
       write(16,*) ' PHASE OF RF AT ENTRANCE(DG) ',phi0*180./pi
       write(16,*) ' CREST PHASE OF RF (DEG) ',PCREST*180./PI
       iarg=1
       call cdg(iarg)
       encog=cog(1)
       gcog=encog/xmat
       bcog=sqrt(1.-1./(gcog*gcog))
       tcog=cog(3)
       CALL EXT2D(1)
c print in file 'short.data'
c 3.12.09       phnew=-(int(tcog*fh/pi+0.5)-tcog*fh/pi)*180.
c 3.12.09       dav1(idav,6)=encog-xmat
c 3.12.09        dav1(idav,7)=phnew
       if(itvol) then
        dav1(idav,38)=dphete
        dav1(idav,39)=dphase*180./pi
       else
        dav1(idav,38)=dphete
       endif
       WRITE(16,3777)
3777   FORMAT(/,' Dynamics at the output',/,
     1 5X,'   BETA     dW(MeV)     ENERGY(MeV) ',
     2 '   TOF(deg)     TOF(sec)')
       engain=encog-enold
       WRITE(16,3473) bets,ddw,enrs-xmat,fh*TREFS *180./PI,TREFS
3473   FORMAT(' REF ',f7.5,3x,f10.6,3x,e12.5,3x,e12.5,3x,e12.5)
       WRITE(16,1789) bcog,engain,encog-xmat,tcog*fh*180./pi,tcog
1789   FORMAT(' COG ',f7.5,3x,f10.6,3x,e12.5,3x,e12.5,3x,e12.5)
       TESTCA=exten(1)*exten(2)*exten(3)
c       epsil=1.E-20
       epsil=1.E-40
       IF(abs(TESTCA).gt.epsil) THEN
         QDISP=2.*SQRT(exten(1))
         QMD=exten(1)*exten(3)-exten(2)**2
         SQMDV=4.*PI*SQRT(QMD)
         SURM=4.*PI*SQRT(QMD)*180./PI
         QDP=2.*SQRT(exten(3))
         COR12=exten(2)/sqrt(exten(1)*exten(3))
comment         PENT12=SQRT(exten(1)/exten(3))/COR12
comment         PENT21=SQRT(exten(1)/exten(3))*COR12
         QDPDE=QDP*180./PI
       ELSE
         QDISP=0.
         QMD=0.
         SQMDV=0.
         SURM=0.
         QDP=0.
         COR12=0.
         PENT12=0.
         PENT21=0.
         QDPDE=0.
       ENDIF
       TRQTX=exten(4)*exten(5)-exten(8)**2
       TRQPY=exten(6)*exten(7)-exten(9)**2
       QDITAX=2.*SQRT(exten(4))
       QDIANT=2.*SQRT(exten(5))
       QDITAY=2.*SQRT(exten(6))
       QDIANP=2.*SQRT(exten(7))
       SURXTH=4.*PI*SQRT(TRQTX)
       SURYPH=4.*PI*SQRT(TRQPY)
       IF(SHIFT) THEN
         vref=bets*vl
         tref=trefs
       ELSE
         vref=bcog*vl
         tref=tcog
       ENDIF
       if(itvol) ttvols=tref
       call statis
C      ENVEL
       CALL STAPL(dav1(idav,24))
cold       WRITE(16,9998) SQMDV
cold 9998   FORMAT(2X,'   EMITTANCE (norm): ',
cold     *        E12.5,' PI*MEV*RAD')
       dav1(idav,16)=bcog*surxth*10./(pi*sqrt(1.-bcog*bcog))
c 3.12.09       dav1(idav,17)=surxth*10./pi
c
       dav1(idav,21)=bcog*suryph*10./(pi*sqrt(1.-bcog*bcog))
       dav1(idav,25)=nrres
       dav1(idav,30)=ngood
c
c   print in the file: 'dynac.dmp':
c   gap number, phase offset(deg), relativistic beta, energy(MeV), horz. emit.(mm*mrd,norm), vert. emit.(mm*mrd,norm),long. emit(keV*sec)
c
c --- dav1(idav,16): Emittance(norm)  x-xp (mm*mrad)
c --- dav1(idav,21): Emittance(norm)  y-yp (mm*mrad)
       emns=1.e12*sqmdv/(pi*fh)
cet2010s
       trfprt=fh*tref*180./pi
       tcgprt=fh*tcog*180./pi
c TEST       n2kp=int(tofprt/360.)
c TEST       tofprt=tofprt-float(n2kp)*360.
c TEST       if(tofprt.gt.180.) tofprt=tofprt-360.
c cavity number, z (m), transmission (%), synchronous phase (deg), time of flight (deg) (within –180 deg and 180 deg),
c COG relativistic beta (@ output), COG output energy (MeV), REF relativistic beta (@ output), REF output energy (MeV),
c horizontal emittance (mm.mrad, RMS normalized), vertical emittance (mm.mrad, RMS normalized),
c longitudinal emittance (RMS, ns.keV)
       trnsms=100.*float(ngood)/float(imax)
       if(ncavnm.eq.1) write(50,*) '# cavnum.dmp'
       if(ncavnm.eq.1) write(50,*) '# cav     Z       trans   ',
     *   'PHIs     TOF(COG)    COG        Wcog          TOF(REF)   ',
     *   '    REF        Wref       Ex,RMS,n     Ey,RMS,n     El,RMS'
       if(ncavnm.eq.1) write(50,*) '#  #     (m)       (%)    ',
     *  '(deg)     (deg)      beta       (MeV)          (deg)      ',
     *  '   beta       (MeV)      (mm.mrad)    (mm.mrad)    (ns.keV)'
       write(50,7023) nrres,0.001*dav1(idav,24),trnsms,dphete,
     *  tcgprt,bcog,encog-xmat,trfprt,bets,enrs-xmat,
     *  0.25*dav1(idav,16),0.25*dav1(idav,21),0.25*emns
7023     format(1x,i4,1x,e12.5,1x,f6.2,1x,f7.2,1x,
     *   2(e14.7,1x,f7.5,1x,e14.7,1x),3(e12.5,1x))
cet2010e
       fh=oldfh
C       new magnetic rigidity of the reference
       gref=1./sqrt(1.-bets*bets)
       xmor=xmat*bets*gref
       BORO=33.356*XMOR*1.E-01/AQST
       WRITE(16,*) ilost,' particles lost in cavity ',nrres
       call emiprt(0)
       return
       end
       SUBROUTINE phcrest(phi0,ylg,ncell,zcrest)
c ..........................................................
c     REFERENCE:
c        average k (cm-1)
c .................................................
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum3/bgt0,bgt1,bgt2,bgt3,bgt4,bgt5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       COMMON/DYN/TREF,VREF
c   ylg(cm) length of cavity, xpas(cm): step length
c   8 steps over a cell
       npas=ncell*10
       xpas=ylg/float(npas)
       e0=xmat
       b0=vref/vl
       b1=b0
       b2=b0
       b3=b0
       b4=b0
       b5=b0
       t0=0.
       gam0=1./sqrt(1.-b0*b0)
       xnh=0.
       estop=xpas/10.
20      continue
ccc        if(xnh*xpas.ge.ylg) go to 10
        xpat=xnh*xpas
       if(xpat.lt.(ylg-estop)) then
        call fposb
c ---  enegy gain over the step
c   predictor of energy gain
       b1=b0
       b2=b0
       b3=b0
       b4=b0
       b5=b0
       dgam=xi1(phi0,t0,t5)*qst/e0
       gam5=gam0+dgam
c   corrector of energy gain
        xpas2=xpas*xpas
        dgdz=qst/e0 * tspl0
        d2gdz2=dgam/xpas2-dgdz/xpas
        d2gdz2=2.*d2gdz2
        gam1=gam0+dgdz*xpas/5.+d2gdz2*xpas2/50.
        gam2=gam0+dgdz*xpas*2./5.+d2gdz2*xpas2*4./50.
        gam3=gam0+dgdz*xpas*3./5.+d2gdz2*xpas2*9./50.
        gam4=gam0+dgdz*xpas*4./5.+d2gdz2*xpas2*16./50.
        b1=sqrt(gam1*gam1-1.)/gam1
        b2=sqrt(gam2*gam2-1.)/gam2
        b3=sqrt(gam3*gam3-1.)/gam3
        b4=sqrt(gam4*gam4-1.)/gam4
        b5=sqrt(gam5*gam5-1.)/gam5
        dgam=xi1(phi0,t0,t5)*qst/e0
        gam5=gam0+dgam
        b5=sqrt(gam5*gam5-1.)/gam5
        b0=b5
        t0=t5
        gam0=gam5
        xnh=xnh+1.
        go to 20
       endif
c*et*2013 10     continue
c   compute an average k: zcrest
       zcrest=fh*t5/ylg
       return
       end
       SUBROUTINE phcrest1(phi0,ylg,ncell)
c ..........................................................
c     REFERENCE:
c       computation of the phase giving the maximum energy gain)
c .................................................
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum3/bgt0,bgt1,bgt2,bgt3,bgt4,bgt5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       COMMON/DYN/TREF,VREF
       dimension stx(400),sty(400),sts(400),stp(400),stq(400)
       dimension phc(100),wph(100)
c   ylg(cm) length of cavity, xpas(cm): step length
       npas=ncell*10
       xpas=ylg/float(npas)
       e0=xmat
       rad=pi/180.
       dph=rad
       dph1=rad
       dplim=rad*0.01
       pmax=phi0+10.*rad
       pmin=phi0-10.*rad
       ibcl=1
       phi=pmin
30     continue
       if(phi.ge.pmax) go to 50
       phc(ibcl)=phi
       b0=vref/vl
       b1=b0
       b2=b0
       b3=b0
       b4=b0
       b5=b0
       t0=0.
       gam0=1./sqrt(1.-b0*b0)
       wwref=(gam0-1.)*e0
       xnh=0.
20      continue
        if(xnh*xpas.ge.ylg) go to 10
        call fposb
c ---  enegy gain over the step
c   predictor of energy gain
       b1=b0
       b2=b0
       b3=b0
       b4=b0
       b5=b0
       dgam=xi1(phi,t0,t5)*qst/e0
       gam5=gam0+dgam
c   corrector
        xpas2=xpas*xpas
        dgdz=qst/e0 * tspl0
        d2gdz2=dgam/(xpas2)-dgdz/xpas
        d2gdz2=2.*d2gdz2
        gam1=gam0+dgdz*xpas/5.+d2gdz2*xpas2/50.
        gam2=gam0+dgdz*xpas*2./5.+d2gdz2*xpas2*4./50.
        gam3=gam0+dgdz*xpas*3./5.+d2gdz2*xpas2*9./50.
        gam4=gam0+dgdz*xpas*4./5.+d2gdz2*xpas2*16./50.
        b1=sqrt(gam1*gam1-1.)/gam1
        b2=sqrt(gam2*gam2-1.)/gam2
        b3=sqrt(gam3*gam3-1.)/gam3
        b4=sqrt(gam4*gam4-1.)/gam4
        b5=sqrt(gam5*gam5-1.)/gam5
        dgam=xi1(phi,t0,t5)*qst/e0
        gam5=gam0+dgam
        b5=sqrt(gam5*gam5-1.)/gam5
        b0=b5
        t0=t5
        gam0=gam5
        xnh=xnh+1.
        go to 20
10     continue
       wwpcr=(gam0-1.)*e0
       dwcpr=wwpcr-wwref
       wph(ibcl)=dwcpr
       phi=phi+dph
       ibcl=ibcl+1
       go to 30
50     continue
       ibcl=ibcl-1
c     save spline areas (partial save)
      do i=1,400
       stx(i)=xspl(i)
       sty(i)=yspl(i)
       sts(i)=s(i)
       stp(i)=p(i)
       stq(i)=q(i)
      enddo
c     padding spline areas
      do i=1,ibcl
       xspl(i)=phc(i)
       yspl(i)=wph(i)
      enddo
      call deriv2(ibcl)
      i=1
      phi=xspl(1)
      yfb=slope(ibcl,phi)/100.
70    continue
      if(phi.ge.xspl(ibcl))go to 60
      if(dph1.le.dplim) go to 60
      yf=slope(ibcl,phi)
ccccc      if(abs(yf).lt.yfb) go to 60
      if(yf.gt.0.) then
       phi=phi+dph1
       go to 70
      else
       phi=phi-dph1
       dph1=dph1/2.
       phi=phi+dph1
       go to 70
      endif
60    continue
      phi0=phi
c     restore spline areas (partial save)
      do i=1,400
       xspl(i)=stx(i)
       yspl(i)=sty(i)
       s(i)=sts(i)
       p(i)=stp(i)
       q(i)=stq(i)
      enddo
      return
      end
       SUBROUTINE dwref(phi0,gam5,t5)
c ..........................................................
c   compute the energy gain and the time of flight of the reference over the cavity(gap)
c    at the exit of the cavity(gap): gam4 = relativistic gamma,  t4 = tof (cavity)
c ....................................................................................
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/cavnum7/sspl0,sspl1,sspl2,sspl3,sspl4,sspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       COMMON/DYN/TREF,VREF
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/step/istep
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/kcell/avrg(15)
       common/xitrd3/bcour1(15),bcour2(15),tr3dw(15),tr3ph(15),
     *  t0tr3d(15)
       common/sphi/tcour1(15),tcour2(15)
       common/t3dfld/fldctr(15),zend(15),T3D
c  save previous xpas
       xpcel=xpas
c
       xlpos=0.
       isce=20
       xlcum=0.
       b0=vref/vl
       b1=b0
       b2=b0
       b3=b0
       b4=b0
       b5=b0
       t0=0.
c   xlcel (cm): length of current cell
c   xlim(inc) cm: limits of current cell
c   isce: step number in cell (forced isce = 20)
c   wdgam: energy gain over the cell
c **************************************
       do inc=1,ncell
        bcour1(inc)=b0
        tcour1(inc)=t0
        wdgams=0.
        xnh=0.
        xlcel=xlim(inc+1)-xlim(inc)
        xlpos=xlpos+xlcel
        xpas=xlcel/float(isce)
        e0=xmat
        xnhc=0
        gam0=1./sqrt(1.-b0*b0)
        estop=xpas/10.
        istop=0
20     continue
       xpat=xnhc*xpas
        if(xpat.lt.(xlcel-estop)) then
         call fposbb(xlcum)
c ---  energy gain over the step
         b1=b0
         b2=b0
         b3=b0
         b4=b0
         b5=b0
c predictor
         dgam=xi1(phi0,t0,t5)*qst/e0
         gam5=gam0+dgam
c   corrector
         xpas2=xpas*xpas
         dgdz=qst/e0 * tspl0
         d2gdz2=dgam/xpas2-dgdz/xpas
         d2gdz2=2.*d2gdz2
         gam1=gam0+dgdz*xpas/5.+d2gdz2*xpas2/50.
         gam2=gam0+dgdz*xpas*2./5.+d2gdz2*xpas2*4./50.
         gam3=gam0+dgdz*xpas*3./5.+d2gdz2*xpas2*9./50.
         gam4=gam0+dgdz*xpas*4./5.+d2gdz2*xpas2*16./50.
         b1=sqrt(gam1*gam1-1.)/gam1
         b2=sqrt(gam2*gam2-1.)/gam2
         b3=sqrt(gam3*gam3-1.)/gam3
         b4=sqrt(gam4*gam4-1.)/gam4
         b5=sqrt(gam5*gam5-1.)/gam5
         dgam=xi1(phi0,t0,t5)*qst/e0
         gam5=gam0+dgam
         eww=xmat*(gam5-1.)
         wdgam=xmat*dgam
         wdgams=wdgams+wdgam
cccc        write(13,*) ' inc wdgam eww t5 ',inc,wdgam,eww,t5
         b5=sqrt(gam5*gam5-1.)/gam5
         b0=b5
         t0=t5
         gam0=gam5
         xnh=xnh+1.
         xnhc=xnhc+1.
         go to 20
       endif
       tcour2(inc)=t5
       xlcum=xlcum+xlcel
       ttvol=t5*fh*180./pi
       ttphi=ttvol+tref*fh*180/pi
       bcour2(inc)=b5
       tr3dw(inc)=wdgams
c ************************************************
       write(13,101) NRRES,inc,xlim(inc),xlim(inc+1),
     *  xlcel,eww,wdgams,ttvol,ttphi,avrg(inc)
101    format(2(2x,i4),8(2x,e12.5))
         fldctr(inc)=avrg(inc)-xlim(inc)
         zend(inc)=xlcel-fldctr(inc)
       enddo
       call itrd3
c restore  xpas
       xpas=xpcel
       return
       end
c **************************************************************
       SUBROUTINE itrd3
       implicit real*8 (a-h,o-z)
       common/cell/yf(10000),xf(10000),xlim(15),flength,fhc,att,
     *             npoint(15),ncell
       common/consta/vl,pi,xmat,rpel,qst
       common/xitrd3/bcour1(15),bcour2(15),tr3dw(15),tr3ph(15),
     *  t0tr3d(15)
       common/cavnum1/xnh,xpas,ffield,npt
       common/sphi/tcour1(15),tcour2(15)
       common/kcell/avrg(15)
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       dimension ttr3d(15),str3d(15)
       do nrc=1,ncell
         bcour=(bcour1(nrc)+bcour2(nrc))/2.
         ttr3d(nrc)=ta0(bcour,nrc)
       enddo
       do nrc=1,ncell
         bcour=(bcour1(nrc)+bcour2(nrc))/2.
         str3d(nrc)=sb0(bcour,nrc)
       enddo
       do nrc=1,ncell
         t0tr3d(nrc)=qst*sqrt(ttr3d(nrc)**2+str3d(nrc)**2)*ffield/2.
         atr3d=0.
         tt=tr3dw(nrc)/t0tr3d(nrc)
         tph=atan(-str3d(nrc)/ttr3d(nrc))
         if(tt.lt.1.) atr3d=acos(tt)
         ddw=ttr3d(nrc)*cos(tph)-str3d(nrc)*sin(tph)
         if(ddw.lt.0.) then
           tph=tph+pi
         endif
c         if(tph.ge.pcrest) atr3d=-atr3d
         drift1=avrg(nrc)-xlim(nrc)
         drift2=xlim(nrc+1)-avrg(nrc)
         xsign1=tcour2(nrc)-tcour1(nrc)
         xsing=xsign1-(drift1/bcour1(nrc)+drift2/bcour2(nrc))/vl
         if(xsing.gt.0.) atr3d=-atr3d
         tr3ph(nrc)=atr3d*180./pi
         write(13,100)nrc,t0tr3d(nrc),ttr3d(nrc),str3d(nrc),tr3dw(nrc)
     *              ,atr3d*180./pi,tph*180./pi
100      format(2x,i4,6(2x,e12.5))
       enddo
       return
       end
c **************************************************************
       SUBROUTINE fposbbb(xlcum,fposs,jx)
c  SV 24/08/2015
c -------------------------------------------------------
c  electric field at the 6 Bode's positions in the step (see routine dwref)
c -------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/cavnum6/fpos0,fpos1,fpos2,fpos3,fpos4,fpos5
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/rfield/ifield
       common/xposi/xpost(10),xlce(2),xpax(2),iscx(2)
       logical ifield
c *******************************************
       if(fposs.eq.0.) then
        fpos0=xnh*xpas+xlcum+fposs
        fpos1=(xnh+0.20)*xpas+xlcum+fposs
        fpos2=(xnh+0.40)*xpas+xlcum+fposs
        fpos3=(xnh+0.60)*xpas+xlcum+fposs
        fpos4=(xnh+0.80)*xpas+xlcum+fposs
        fpos5=(xnh+1.0)*xpas+xlcum+fposs
       else
        fpos0=xnh*xpas+fposs
        fpos1=(xnh+0.20)*xpas+fposs
        fpos2=(xnh+0.40)*xpas+fposs
        fpos3=(xnh+0.60)*xpas+fposs
        fpos4=(xnh+0.80)*xpas+fposs
        fpos5=(xnh+1.0)*xpas+fposs
       endif
       if(jx.eq.1)then
        xnh1=xnh+1
        if(xnh1.eq.float(iscx(jx))) fposs=fpos5
       endif
       if(ifield) then
c --- the field is read from disk in MV/cm
         tspl0=spline(npt,fpos0)*fmult
         tspl1=spline(npt,fpos1)*fmult
         tspl2=spline(npt,fpos2)*fmult
         tspl3=spline(npt,fpos3)*fmult
         tspl4=spline(npt,fpos4)*fmult
         tspl5=spline(npt,fpos5)*fmult
c ****TEST
ccc         write(13,*) 'xnh fpos tspl '
c          write(13,700) xnh,fpos0,tspl0,xpas,fposs
c         write(13,700) xnh,fpos1,tspl1
c         write(13,700) xnh,fpos2,tspl2
c         write(13,700) xnh,fpos3,tspl3
c         write(13,700) xnh,fpos4,tspl4
c         write(13,700) xnh,fpos5,tspl5,xpas,fposs
c700      format(5(2x,e12.5))
c ***********************************
         go to 10
      else
c ++++ THIS OPTION IS NOT AVAILABLE IN THE CODE +++++++++++++
c the field (MV/cm) is given in the form of a Fourier series
         tspl0=fone(fpos0)*fmult
         tspl1=fone(fpos1)*fmult
         tspl2=fone(fpos2)*fmult
         tspl3=fone(fpos3)*fmult
         tspl4=fone(fpos4)*fmult
         tspl5=fone(fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
cc         write(70,700) fpos0,tspl0
cc         write(70,700) fpos1,tspl1
cc         write(70,700) fpos2,tspl2
cc         write(70,700) fpos3,tspl3
cc         write(70,700) fpos4,tspl4
cold         write(70,700) fpos5,tspl5
c ***********************************
       endif
10     continue
       return
       end
       SUBROUTINE fposbb(xlcum)
c  SV 24/08/2015
c -------------------------------------------------------
c  electric field at the 6 Bode's positions in the step  (see routine dwref1)
c -------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/cavnum6/fpos0,fpos1,fpos2,fpos3,fpos4,fpos5
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/rfield/ifield
       logical ifield
c *******************************************
       fpos0=xnh*xpas+xlcum
       fpos1=(xnh+0.20)*xpas+xlcum
       fpos2=(xnh+0.40)*xpas+xlcum
       fpos3=(xnh+0.60)*xpas+xlcum
       fpos4=(xnh+0.80)*xpas+xlcum
       fpos5=(xnh+1.0)*xpas+xlcum
       if(ifield) then
c --- the field is read from disk in MV/cm
         tspl0=spline(npt,fpos0)*fmult
         tspl1=spline(npt,fpos1)*fmult
         tspl2=spline(npt,fpos2)*fmult
         tspl3=spline(npt,fpos3)*fmult
         tspl4=spline(npt,fpos4)*fmult
         tspl5=spline(npt,fpos5)*fmult
c ****TEST
ccc         write(13,*) 'xnh fpos tspl '
cc         write(13,700) xnh,fpos0,tspl0,xpas
c         write(13,700) xnh,fpos1,tspl1
c         write(13,700) xnh,fpos2,tspl2
c         write(13,700) xnh,fpos3,tspl3
c         write(13,700) xnh,fpos4,tspl4
c         write(13,700) xnh,fpos5,tspl5,xpas
c700      format(4(2x,e12.5))
c ***********************************
         go to 10
      else
c the field (MV/cm) is given in the form of a Fourier series
         tspl0=fone(fpos0)*fmult
         tspl1=fone(fpos1)*fmult
         tspl2=fone(fpos2)*fmult
         tspl3=fone(fpos3)*fmult
         tspl4=fone(fpos4)*fmult
         tspl5=fone(fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
cc         write(70,700) fpos0,tspl0
cc         write(70,700) fpos1,tspl1
cc         write(70,700) fpos2,tspl2
cc         write(70,700) fpos3,tspl3
cc         write(70,700) fpos4,tspl4
cold         write(70,700) fpos5,tspl5
c ***********************************
       endif
10     continue
       return
       end
       SUBROUTINE fposb
c -------------------------------------------------------
c  electric field at the 6 Bode's positions in the step
c -------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/cavnum6/fpos0,fpos1,fpos2,fpos3,fpos4,fpos5
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/rfield/ifield
       logical ifield
       fpos0=xnh*xpas
       fpos1=(xnh+0.20)*xpas
       fpos2=(xnh+0.40)*xpas
       fpos3=(xnh+0.60)*xpas
       fpos4=(xnh+0.80)*xpas
       fpos5=(xnh+1.0)*xpas
       if(ifield) then
c --- the field is read from disk in MV/cm
         tspl0=spline(npt,fpos0)*fmult
         tspl1=spline(npt,fpos1)*fmult
         tspl2=spline(npt,fpos2)*fmult
         tspl3=spline(npt,fpos3)*fmult
         tspl4=spline(npt,fpos4)*fmult
         tspl5=spline(npt,fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
cc         write(13,700) xnh,fpos0,tspl0
cc         write(13,700) xnh,fpos1,tspl1
cc         write(13,700) xnh,fpos2,tspl2
cc         write(13,700) xnh,fpos3,tspl3
cc         write(13,700) xnh,fpos4,tspl4
cc         write(13,700) xnh,fpos5,tspl5
cc700      format(3(2x,e12.5))
c ***********************************
         go to 10
      else
c the field (MV/cm) is given in the form of a Fourier series
         tspl0=fone(fpos0)*fmult
         tspl1=fone(fpos1)*fmult
         tspl2=fone(fpos2)*fmult
         tspl3=fone(fpos3)*fmult
         tspl4=fone(fpos4)*fmult
         tspl5=fone(fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
cc         write(70,700) fpos0,tspl0
cc         write(70,700) fpos1,tspl1
cc         write(70,700) fpos2,tspl2
cc         write(70,700) fpos3,tspl3
cc         write(70,700) fpos4,tspl4
cold         write(70,700) fpos5,tspl5
c ***********************************
       endif
10     continue
       return
       end
       SUBROUTINE sposb
c -------------------------------------------------------
c  derivative of electric field at the 6 Bode's positions in the step
c    not used in the code
c -------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/cavnum6/fpos0,fpos1,fpos2,fpos3,fpos4,fpos5
       common/cavnum7/sspl0,sspl1,sspl2,sspl3,sspl4,sspl5
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       common/rfield/ifield
       logical ifield
       fpos0=xnh*xpas
       fpos1=(xnh+0.20)*xpas
       fpos2=(xnh+0.40)*xpas
       fpos3=(xnh+0.60)*xpas
       fpos4=(xnh+0.80)*xpas
       fpos5=(xnh+1.0)*xpas
       if(ifield) then
c --- the field is read from disk in MV/cm
         sspl0=slope(npt,fpos0)*fmult
         sspl1=slope(npt,fpos1)*fmult
         sspl2=slope(npt,fpos2)*fmult
         sspl3=slope(npt,fpos3)*fmult
         sspl4=slope(npt,fpos4)*fmult
         sspl5=slope(npt,fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
ccc         write(13,700) xnh,fpos0,sspl0
ccc         write(13,700) xnh,fpos1,sspl1
ccc         write(13,700) xnh,fpos2,sspl2
ccc         write(13,700) xnh,fpos3,sspl3
ccc         write(13,700) xnh,fpos4,sspl4
ccc         write(13,700) xnh,fpos5,sspl5
ccc700      format(3(2x,e12.5))
c ***********************************
         go to 10
      else
c the field (MV/cm) is given in the form of a Fourier series
         tspl0=fone(fpos0)*fmult
         tspl1=fone(fpos1)*fmult
         tspl2=fone(fpos2)*fmult
         tspl3=fone(fpos3)*fmult
         tspl4=fone(fpos4)*fmult
         tspl5=fone(fpos5)*fmult
c ****TEST
cold         write(6,*) 'xnh fmult ',xnh,fmult
cc         write(70,700) fpos0,tspl0
cc         write(70,700) fpos1,tspl1
cc         write(70,700) fpos2,tspl2
cc         write(70,700) fpos3,tspl3
cc         write(70,700) fpos4,tspl4
cold         write(70,700) fpos5,tspl5
c ***********************************
       endif
10     continue
       return
       end
       FUNCTION xi1(phi0,t0,t5)
c ----------------------------------------------
c     energy gain over the step
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       xi1=0.
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2
       b34=(b3+b4)/2
       b45=(b4+b5)/2
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
       xspl0=cos(fh*t0+phi0)*tspl0
       xspl1=cos(fh*t1+phi0)*tspl1
       xspl2=cos(fh*t2+phi0)*tspl2
       xspl3=cos(fh*t3+phi0)*tspl3
       xspl4=cos(fh*t4+phi0)*tspl4
       xspl5=cos(fh*t5+phi0)*tspl5
       tspl=19.*xspl0+75.*xspl1+50.*xspl2+50.*xspl3+75.*xspl4
     *     +19.*xspl5
       xi1=xpas/288. * tspl
       return
       end
       FUNCTION xi2(phi0,t0)
c ----------------------------------------------
c     coupling terms in R and R' (energy gain)
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2.
       b34=(b3+b4)/2.
       b45=(b4+b5)/2.
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
cold       xspl0=cos(fh*t0+phi0)*tspl0
       xspl1=cos(fh*t1+phi0)*tspl1
       xspl2=cos(fh*t2+phi0)*tspl2
       xspl3=cos(fh*t3+phi0)*tspl3
       xspl4=cos(fh*t4+phi0)*tspl4
       xspl5=cos(fh*t5+phi0)*tspl5
       tspl=15.*xspl1+20.*xspl2+30.*xspl3+60*xspl4+19.*xspl5
       xi2=xpas*xpas/288. * tspl
       return
       end
       FUNCTION xj1(phi0,t0)
c ----------------------------------------------
c     transverse motion field dE/dt
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum3/bgt0,bgt1,bgt2,bgt3,bgt4,bgt5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
c --- the field is read on the disk
       xj1=0.
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2.
       b34=(b3+b4)/2.
       b45=(b4+b5)/2.
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
       xspl0=-fh*sin(fh*t0+phi0)*tspl0/bgt0
       xspl1=-fh*sin(fh*t1+phi0)*tspl1/bgt1
       xspl2=-fh*sin(fh*t2+phi0)*tspl2/bgt2
       xspl3=-fh*sin(fh*t3+phi0)*tspl3/bgt3
       xspl4=-fh*sin(fh*t4+phi0)*tspl4/bgt4
       xspl5=-fh*sin(fh*t5+phi0)*tspl5/bgt5
       tspl=19.*xspl0+75.*xspl1+50.*xspl2+50.*xspl3+75.*xspl4
     *     +19.*xspl5
       xj1=xpas/288. * tspl
       return
       end
       FUNCTION xj2(phi0,t0)
c ----------------------------------------------
c     transverse motion field dE/dt
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum3/bgt0,bgt1,bgt2,bgt3,bgt4,bgt5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       xj2=0.
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2.
       b34=(b3+b4)/2.
       b45=(b4+b5)/2.
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
       xspl1=-fh*sin(fh*t1+phi0)*tspl1/bgt1
       xspl2=-fh*sin(fh*t2+phi0)*tspl2/bgt2
       xspl3=-fh*sin(fh*t3+phi0)*tspl3/bgt3
       xspl4=-fh*sin(fh*t4+phi0)*tspl4/bgt4
       xspl5=-fh*sin(fh*t5+phi0)*tspl5/bgt5
       tspl=15.*xspl1+20.*xspl2+30.*xspl3
     *          +60.*xspl4+19.*xspl5
       xj2=xpas*xpas/288. * tspl
       return
       end
       FUNCTION xe21(phi0,t0)
c ----------------------------------------------
c     transverse motion field E*E
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum4/bge0,bge1,bge2,bge3,bge4,bge5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       xe21=0.
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2.
       b34=(b3+b4)/2.
       b45=(b4+b5)/2.
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
       xspl0=cos(fh*t0+phi0)*tspl0
       xspl1=cos(fh*t1+phi0)*tspl1
       xspl2=cos(fh*t2+phi0)*tspl2
       xspl3=cos(fh*t3+phi0)*tspl3
       xspl4=cos(fh*t4+phi0)*tspl4
       xspl5=cos(fh*t5+phi0)*tspl5
       xspl0=xspl0*xspl0*bge0
       xspl1=xspl1*xspl1*bge1
       xspl2=xspl2*xspl2*bge2
       xspl3=xspl3*xspl3*bge3
       xspl4=xspl4*xspl4*bge4
       xspl5=xspl5*xspl5*bge5
       tspl=19.*xspl0+75.*xspl1+50.*xspl2+50.*xspl3+75.*xspl4
     *     +19.*xspl5
       xe21=xpas/288. * tspl
       return
       end
       FUNCTION xe22(phi0,t0)
c ----------------------------------------------
c     transverse motion for field E*E
c ----------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum4/bge0,bge1,bge2,bge3,bge4,bge5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/consta/vl,pi,xmat,rpel,qst
       common/spl/xspl(4000),yspl(4000),s(5000),p(5000),q(5000)
       COMMON/DYN/TREF,VREF
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/rfield/ifield
       logical ifield
       xe22=0.
       b01=(b0+b1)/2.
       b12=(b1+b2)/2.
       b23=(b2+b3)/2.
       b34=(b3+b4)/2.
       b45=(b4+b5)/2.
       t1=t0+xpas/(5.*b01*vl)
       t2=t1+xpas/(5.*b12*vl)
       t3=t2+xpas/(5.*b23*vl)
       t4=t3+xpas/(5.*b34*vl)
       t5=t4+xpas/(5.*b45*vl)
       xspl1=cos(fh*t1+phi0)*tspl1
       xspl2=cos(fh*t2+phi0)*tspl2
       xspl3=cos(fh*t3+phi0)*tspl3
       xspl4=cos(fh*t4+phi0)*tspl4
       xspl5=cos(fh*t4+phi0)*tspl5
       xspl1=xspl1*xspl1*bge1
       xspl2=xspl2*xspl2*bge2
       xspl3=xspl3*xspl3*bge3
       xspl4=xspl4*xspl4*bge4
       xspl5=xspl5*xspl5*bge5
       tspl=15.*xspl1+20.*xspl2+30.*xspl3
     *          +60.*xspl4+19.*xspl5
       xe22=xpas*xpas/288. * tspl
       return
       end
      SUBROUTINE bcnum(phref,ylg,ncell)
c .............................................
c  dynamics of the bunch
c ....................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/dcspa/iesp
       common/ttfc/tk,t1k,t2k,sk,s1k,s2k,fh
       common/faisc/f(10,iptsz),imax,ngood
       common/consta/ vl,pi, xmat,rpel,qst
       common/cavnum1/xnh,xpas,fmult,npt
       common/cavnum2/b0,b1,b2,b3,b4,b5
       common/cavnum3/bgt0,bgt1,bgt2,bgt3,bgt4,bgt5
       common/cavnum4/bge0,bge1,bge2,bge3,bge4,bge5
       common/cavnum5/tspl0,tspl1,tspl2,tspl3,tspl4,tspl5
       common/rfield/ifield
       common/TESTREF/trefs,ddw
       common/elec/jelec
       common/step/istep
       dimension gam(500),xe(500),xpe(500),ye(500),ype(500)
       dimension tcour(iptsz),phi(iptsz)
       logical ifield,flgsc,ichaes,iesp,jelec
c   ylg(cm) length of cavity, xpas(cm): step length
c   8 steps over a cell (a voir!!!)
       BEREF=VREF/VL
       GAMREF=1./SQRT(1.-(BEREF*BEREF))
       ENREF=XMAT*GAMREF
       e0=xmat
       npas=ncell*istep
       xpas=ylg/float(npas)
       npas1=npas+1
       xnh=0.
c   flgsc = true ---> s.c. computation
c   eglsc = 2*xpas : acting length of s.c. computation
       eglsc=2.*xpas
       flgsc=.false.
       do i=2,npas1
        i1=i-1
c --- seek the field E(z) values in the 6 positions in the step length xh
        call fposb
        do j=1,ngood
         qc=f(9,j)
         gam0=f(7,j)/e0
         gam(i1)=gam0
         if(i1.eq.1) then
           tcour(j)=0.
           tof=f(6,j)
C --- rphas: phase delay between the actual particle  and the reference (entrance of the cavity)
           rphas=fh*(tof-tref)
           phi(j)=phref+rphas
         endif
         t0=tcour(j)
         ddt=t0
c  predictor (energy gain)
          b0=sqrt(gam0*gam0-1.)/gam0
          b1=b0
          b2=b0
          b3=b0
          b4=b0
          b5=b0
          dgam=xi1(phi(j),t0,t5)*qc/e0
          gam5=gam(i1)+dgam
c   corrector (energy gain)
c   tspl0 = dE/dz  (MV/(cm*cm)
        xpas2=xpas*xpas
        dgdz=qst/e0 * tspl0
        d2gdz2=dgam/xpas2-dgdz/xpas
        d2gdz2=2.*d2gdz2
        gam1=gam0+dgdz*xpas/5.+d2gdz2*xpas2/50.
        gam2=gam0+dgdz*xpas*2./5.+d2gdz2*xpas2*4./50.
        gam3=gam0+dgdz*xpas*3./5.+d2gdz2*xpas2*9./50.
        gam4=gam0+dgdz*xpas*4./5.+d2gdz2*xpas2*16./50.
        b1=sqrt(gam1*gam1-1.)/gam1
        b2=sqrt(gam2*gam2-1.)/gam2
        b3=sqrt(gam3*gam3-1.)/gam3
        b4=sqrt(gam4*gam4-1.)/gam4
        b5=sqrt(gam5*gam5-1.)/gam5
        dgam=xi1(phi(j),t0,t5)*qc/e0
        gam5=gam(i1)+dgam
        b5=sqrt(gam5*gam5-1.)/gam5
c      tranverse coordinates in (cm,rad)
         x0=f(2,j)
         y0=f(4,j)
         xt0=f(3,j)*1.e-03
         yp0=f(5,j)*1.e-03
c        Picht transformation: xe0 and ye0 (cm) xpe0 and ype0 (rad)
         gamm0=(gam0*gam0-1.)**0.25
         xe0=x0*gamm0
         xpe0=xt0*gamm0
         ye0=y0*gamm0
         ype0=yp0*gamm0
         xpe0=xpe0+.5*xe0*gam0*dgdz/(gam0*gam0-1.)
         ype0=ype0+.5*ye0*gam0*dgdz/(gam0*gam0-1.)
         xe(i1)=xe0
         xpe(i1)=xpe0
         ye(i1)=ye0
         ype(i1)=ype0
c       transverse coupling terms
         gam00=gam0*gam0
         gam11=gam1*gam1
         gam22=gam2*gam2
         gam33=gam3*gam3
         gam44=gam4*gam4
         bgt0=(gam00-1.)**1.5
         xk1=fh*fh/(4.*vl*vl*bgt0)
         red=sqrt(xe(i1)*xe(i1)+ye(i1)*ye(i1))
         red2=red*red
         dred=0.
         if(red.gt.1.e-08) then
          dred=xe(i1)*xpe(i1)+ye(i1)*ype(i1)
          dred=dred/red
         endif
         rk1=xk1*red2*xi1(phi(j),t0,t5)*qc/e0
         rk2=red*dred*xk1*xi2(phi(j),t0)*qc/e0
         gam(i)=gam5+rk1+rk2
         gam55=gam(i)*gam(i)
         bgt1=(gam11-1.)**1.5
         bgt2=(gam22-1.)**1.5
         bgt3=(gam33-1.)**1.5
         bgt4=(gam44-1.)**1.5
         bgt5=(gam55-1.)**1.5
         bge0=(gam00+2.)/((gam00-1.)*(gam00-1.))
         bge1=(gam11+2.)/((gam11-1.)*(gam11-1.))
         bge2=(gam22+2.)/((gam22-1.)*(gam22-1.))
         bge3=(gam33+2.)/((gam33-1.)*(gam33-1.))
         bge4=(gam44+2.)/((gam44-1.)*(gam44-1.))
         bge5=(gam55+2.)/((gam55-1.)*(gam55-1.))
c ******************************************************
c   compute the jump of phase over the step
cold         dlt=(1.+red2*xk1)*xi3(phi(j),t0)
cold         dlt=dlt+red*dred*xk1*xi4(phi(j),t0)
cold         dlt=dlt*qc/(e0*vl)
cold         tof=tof+xpas/(vl*b0)+dlt
c *******************************************************
         f(7,j)=gam(i)*e0
         tcour(j)=t5
         ddt1=t5-ddt
         f(6,j)=ddt1+f(6,j)
c -- angular deviation
c    1) terms in dE/dt
         a1=qc/(2.*e0*vl)
         ttt1=xj1(phi(j),t0)
         ttt2=xj2(phi(j),t0)
cold         dxpe=xe(i1)*xj1(phi(j),t0)+xpe(i1)*xj2(phi(j),t0)
cold         dype=ye(i1)*xj1(phi(j),t0)+ype(i1)*xj2(phi(j),t0)
         dxpe1=xe(i1)*ttt1+xpe(i1)*ttt2
         dype1=ye(i1)*ttt1+ype(i1)*ttt2
         xpe(i)=xpe(i1)+a1*dxpe1
         ype(i)=ype(i1)+a1*dype1
c --- 2) terms in E*E (only for no-relativistic electrons)
       if(jelec) then
         ae2=qc/(2.*e0)
         ae2=ae2*ae2
         stt1=xe21(phi(j),t0)
         stt2=xe22(phi(j),t0)
         dxpe2=xe(i1)*stt1+xpe(i1)*stt2
         dype2=ye(i1)*stt1+ype(i1)*stt2
         xpe(i)=xpe(i1)+a1*dxpe1-ae2*dxpe2
         ype(i)=ype(i1)+a1*dype1-ae2*dype2
        endif
c  extension
cold         dxe=xe(i1)*xj2(phi,t0)+xpe(i1)*xj3(phi,t0)
cold         dye=ye(i1)*xj2(phi,t0)+ype(i1)*xj3(phi,t0)
cold         xe(i)=xe(i1)+a1*dxe+xpas*xpe(i1)
cold         ye(i)=ye(i1)+a1*dye+xpas*ype(i1)
         xe(i)=xe(i1)+xpas*(xpe(i1)+xpe(i))/2.
         ye(i)=ye(i1)+xpas*(ype(i1)+ype(i))/2.
c       back to the real variables and convert to (cm,mrad)
        dgdzr=qc/e0 * tspl5
        gamm1=(gam(i)*gam(i)-1.)**0.25
        gamm2=(gam(i)*gam(i)-1.)**1.25
        xi=xe(i)/gamm1
        xpi=xpe(i)/gamm1-xe(i)*gam(i)*dgdzr/(gamm2*2.)
        yi=ye(i)/gamm1
        ypi=ype(i)/gamm1-ye(i)*gam(i)*dgdzr/(gamm2*2.)
c      convert in cm and mrd
        f(2,j)=xi
        f(4,j)=yi
        f(3,j)=xpi*1.e03
        f(5,j)=ypi*1.e03
       enddo
c  space charge computation (only odd step numbers)
        if(.not.flgsc) then
         flgsc=.true.
         call disp
cet        endif
        else
cet        if(flgsc) then
           if(ichaes) then
C      Charge space (only SCHEFF is available)
           iesp=.true.
           call cesp(eglsc)
           iesp=.false.
           flgsc=.false.
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
          endif
          call disp
        endif
       xnh=xnh+1.
      enddo
      return
      end
C**********************************************************************************
       SUBROUTINE reject(ilost)
c   ............................................................................
C          reject particles outside window set by REJECT card
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       common/mcs/imcs,ncstat,cstat(20)
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DYN/TREF,VREF
c       common/tapes/in,ifile,meta
c       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
c       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
c       COMMON/DCSPA/IESP
c       COMMON/SHIF/DTIPH,SHIFT
c       common/femt/iemgrw,iemqesg
c       common/posc/xpsc
c       common/tofev/ttvols
c       COMMON/ITVOLE/ITVOL,IMAMIN
c       logical iesp,ichaes,shift,iemgrw,itvol,imamin
c ---- ifw = 0  ===> wdisp = dW/W relative to cog
c ---- ifw = 1  ===> wdisp = dW (MeV) relative to cog
c ---- ifw = 10 ===> wdisp = dW/W relative to REF
c ---- ifw = 11 ===> wdisp = dW (MeV) relative to REF
       ilost=0
       FH0=FH/VL
       fcpi=fh*180./pi
c Test window
       write(16,*)'Check if the ',ngood,' particles are within window'
       write(16,*) 'Number of charge states: ',ncstat
       if(ifw.le.1) then
C Relative to COG
         write(16,3900)
c ---- ifw = 0 ===> wdisp = dW/W
c ---- ifw = 1 ===> wdisp = dW (MeV)
         if (ncstat.gt.1) call cogetc
         bcour=0.
         cgtv=0.
         do i=1,ngood
           gpai=f(7,i)/xmat
           bcour=sqrt(1.-1./(gpai*gpai))+bcour
           cgtv=cgtv+f(6,i)
         enddo
         cgtv=cgtv/float(ngood)
         bcour=bcour/float(ngood)
         bcg=bcour
         gcour=1./sqrt(1.-bcour*bcour)
         gcg=gcour
         wcg=(gcour-1.)*xmat
         do i=1,ngood
           gpai=f(7,i)/xmat
           if(gpai.lt.1.) gpai=1.
           bcour=sqrt(1.-1./(gpai*gpai))
           fd(i)=bcour/bcg * gpai/gcg
         enddo
         if(ifw.eq.0) then
           dispr=gcour*gcour*wdisp/(gcour*(gcour+1.))
         else
           dispr=gcour*gcour*wdisp/(gcour*(gcour+1.)*wcg)
         endif
       else
C Relative to REF
         write(16,3901)
         bref=VREF/VL
         gref=1./sqrt(1.-bref*bref)
         wref=(gref-1.)*xmat
         do i=1,ngood
            gpai=f(7,i)/xmat
            bcour=sqrt(1.-1./(gpai*gpai))
            fd(i)=bcour/bref * gpai/gref
         enddo
         if(ifw.eq.10) then
           dispr=gref*gref*wdisp/(gref*(gref+1.))
         else
           dispr=gref*gref*wdisp/(gref*(gref+1.)*wref)
         endif
       endif
       write(16,3927)   rlim,wx,wy,wphas,wdisp
c     *         ,dispr,f6i*fh*180./pi,f(7,i)-xmat,int(f(9,i))
       f6i=0.
       do i=1,ngood
         ray=f(2,i)*f(2,i)+f(4,i)*f(4,i)
         ray=sqrt(ray)
         if(ray.gt.rlim) f(8,i)=0.
         if(abs(f(2,i)).gt.wx) f(8,i)=0.
         if(abs(f(4,i)).gt.wy) f(8,i)=0.
         if(ifw.le.1) then
C Relative to cog
           if (ncstat.gt.1) then
c---- each charge state has its own COG in phase
             do istc=1,ncstat
               if(f(9,i).eq.charm(istc))then
                 f6i=f(6,i)-cgtdv(istc)
c             tmp=cgtdv(istc)
               endif
             enddo
           else
             f6i=f(6,i)-cgtv
           endif
         else
C Relative to REF
           f6i=f(6,i)-TREF
         endif
         if(fh*abs(f6i).ge.wphas) then
           f(8,i)=0.
c           write(16,*) "Test1:",fh*abs(f6i),wphas,tmp,f(6,i)
         endif
         if(abs(fd(i)-1.).ge.dispr) then
           f(8,i)=0.
c           write(16,*) "Test2:",abs(fd(i)-1.)
         endif
         if(f(8,i).eq.0.) then
           write(16,3928) i,int(f(1,i)),f(2,i),f(3,i),f(4,i),f(5,i),
     *                  f6i*fh*180./pi,f(7,i)-xmat,f(9,i)
           ilost=ilost+1
         endif
       enddo
c  Reshuffles f(i,j) array after window
       call shuffle
ccccc       cstat=1
       ncstat=1
       cstat(1)=f(9,1)
       do j=2,ngood
         mcstat=0
         do k=1,ncstat
           if(f(9,j).eq.cstat(k)) then
             mcstat=1
           endif
         enddo
         if(mcstat.eq.0) then
           ncstat=ncstat+1
           cstat(ncstat)=f(9,j)
         endif
       enddo
       netac=ncstat
       write(16,*) 'Number of good particles left: ',ngood
       write(16,*) 'Number of charge states left : ',ncstat
       write(16,4030) (cstat(j),j=1,ncstat)
       imcs=0
       if(ncstat.gt.1) imcs=1
3900   FORMAT(' Window w.r.t. COG')
3901   FORMAT(' Window w.r.t. reference particle')
3927   FORMAT(' LIM R,X,Y ',3(f10.2,9x),'P,W ',e12.5,9x,e12.5)
3928   FORMAT(' # ',i5,1x,i5,1x,6(f10.2,1x),1x,f5.2)
4030   FORMAT('Charge state(s): ',20(f5.1,1x))
       RETURN
       END
C************************************************************************
C           FIRST AND SECOND order Bending Magnet                       *
C************************************************************************
       SUBROUTINE aimalv (ANGL,RMO,BAIM,XN,XB,EK1,EK2,PENT1,RAB1,
     *                   SK1,SK2,PENT2,RAB2)
c   ....................................................................
c      Bending magnet
c   ....................................................................
C WEDGE BENDING MAGNET
C   ANGL : DEG  bend angle of the central trajectory
C   RMO  : CM   radius of curvature of the central trajectory
C   BAIM : KG   field of the bending magnet
c     BAIM = 0  the field is computed from the momentum of the reference
c               otherwise the momentum is computed from the field
C   XN   :      FIELD GRADIENT (dimensionless,TRANSPORT: n)
C   XB   :      NORMALIZED SECOND DERIVATIVE OF B (TRANSPORT : beta)
C   AP(1) = AP(2) CM vertical half aperture (only if IPOLE = 0)
C ENTRANCE FACE
C  PENT1 EK1 EK2 RAB1
c   PENT1: DEG   angle of pole face rotation  (deg)
c   RAB1 : CM    radius of curvature
c   EK1  :       integral related to the extent of the fringing field (TRANSPORT K1)
c   EK2  :       integral related to the extent of the fringing field (TRANSPORT K2)
C   AP(1) : CM   vertical half aperture
C EXIT FACE
c  PENT2 SK1 SK2  RAB2
c   PENT2: DEG   angle of pole face rotation
c   RAB2 : CM  radius of curvature
c   SK1  :     integral related to the extent of the fringing field
c   SK2  :     integral related to the extent of the fringing field
C   AP(2) : CM   vertical half aperture
c ****************************************************************************
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/TAPES/IN,IFILE,META
       COMMON/DYN/TREF,VREF
       COMMON/RIGID/BORO
       common/faisc/f(10,iptsz),imax,ngood
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON /BLOC23/ H, DEVI, NB, BDB,L
       real*8 L,NB
       COMMON/PORO/IROT1,IROT2
       LOGICAL IROT1,IROT2
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       COMMON /BLOC21/ BE, APB(2), LAYL, LAYX, RABT
       real*8 LAYL, LAYX
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/rander/ialin
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/RADIA/TRT,RSYNC,XINTF,CRAE
       COMMON/RAYSHY/IRAYSH
       LOGICAL IRAYSH
       common/itvole/itvol,imamin
       common/tofev/ttvols
       common/mcs/imcs,ncstat,cstat(20)
c nsprint: control
       common/isector/nsector,nsprint
c --------------------------------------
       logical itvol,imamin,ichaes
       character*1 cr
       logical ialin
       dimension xmoy(20),ymoy(20),rmoy(20),rig(20),ncs(20)
       dimension xpmoy(20),ypmoy(20)
       dimension xcl1(20),xcl2(20),alp(20),sxeb1(20),charge(20)
       dimension xsa1b1(20),baims(20)
c*et*2013*Apr*1
       dimension sbeta(20)
       WRITE(16,100)
100    FORMAT(' ****** BENDING MAGNET: input list ****** ')
       if(baim.eq.0.0) then
c     BORO: momentum of reference (kG.cm)
         ri=boro
c --- BAIM: bend field (KG)
         baim=abs(ri/rmo)
         write(16,'(A,F4.1,A,F12.5,A)')' Based on reference charge ',
     *         qst,' momentum ',boro,' (kG.cm)'
       else
c --- RI momentum kG.cm
         ri=baim*rmo
       endif
       DEVI=ANGL
       NB=XN
       BDB=XB
       RSYNC=rmo
       gap=wy
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
cold       IF(IROT1) WRITE(16,1010) PENT1,RAB1,EK1,EK2,APB(1)
       if(ichaes) then
        write(16,*)'***** beam current: ',beamc,' mA'
       endif
       WRITE(16,1010) PENT1,RAB1,EK1,EK2,APB(1)
1010    FORMAT('    ENTRANCE FACE ******',/,
     1 '  ANGLE OF POLE FACE ROTATION  ',E12.5,' DEG',/,
     1 '  RADIUS OF CURVATURE          ',E12.5,' CM',/,
     1 '  FRINGE FIELD CORECTIONS K1  K2',2(2X,E12.5),/,
     1 '  VERTICAL HALF-APERTURE       ',E12.5,' CM')
       WRITE(16,1020)RMO,DEVI,BAIM,NB,BDB
1020    FORMAT('    WEDGE MAGNET************',/,
     2 '  BENDING RADIUS:  ',E12.5,' CM ',/,
     2 '  BEND ANGLE:      ',E12.5,' DEG',/,
     2 '  FIELD:           ',E12.5,' KG',/,
     2 '  FIELD GRADIENTS: N ',E12.5,' BETA:',E12.5)
cold       IF(IROT2) WRITE(16,1030) PENT2,RAB2,SK1,SK2,APB(2)
        WRITE(16,1030) PENT2,RAB2,SK1,SK2,APB(2)
1030    FORMAT('    EXIT FACE******',/,
     3 '  ANGLE OF POLE FACE ROTATION  ',E12.5,' DEG',/,
     3 '  RADIUS OF CURVATURE          ',E12.5,'CM',/,
     3 '  FRINGE FIELD CORRECTIONS K1  K2',2(2X,E12.5),/,
     3 '  VERTICAL HALF-APERTURE       ',E12.5,' CM')
C  start prints in file 'short.data'
       idav=idav+1
       iitem(idav)=4
       dav1(idav,2)=devi
       dav1(idav,3)=rmo*10.
       dav1(idav,5)=APB(1)*10.
       dav1(idav,6)=PENT1
       dav1(idav,7)=EK1
       dav1(idav,8)=EK2
       dav1(idav,9)=RAB1*10.
       dav1(idav,10)=PENT2
       dav1(idav,11)=SK1
       dav1(idav,12)=SK2
       dav1(idav,13)=RAB2*10.
       dav1(idav,14)=nb
       dav1(idav,15)=bdb
       dav1(idav,16)=baim*.1
       dav1(idav,17)=APB(2)*10.
       FH0=FH/VL
c  Conversion deg--->rad
       radia=pi/180.
       pent1=pent1*radia
       devi=devi*radia
       devr=devi
       l=devr*rmo
       dav1(idav,1)=l*10.
       sdavtot=davtot
       davtot=davtot+l
       dav1(idav,4)=davtot*10.
       pent2=pent2*radia
       devtot=devi
c ----------------------------
c  define the bending angle of the synchronous particle over each sector
       devit=devi/float(nsector)
c --- devit must be different from pent2
       if(devit.eq.pent2) nsector=nsector+1
c --- space charge computation: nsector must be g.t. 1
       if(ichaes.and.(nsector.eq.1) ) nsector = 2
       devi=devi/float(nsector)
       devr=devr/float(nsector)
c ----------------------------------------------------------
c -- scl effective length for space charge computation
c     scl is the path length of the C.T. over two following sectors
       scl=2.*devi*rmo
c   save pent1 , pent2 , devi, ek1, ek2 ,sk1 , sk2 ,
       pent1s=pent1
       pent2s=pent2
       devis=devi
       devrs=devr
       ek1s=ek1
       ek2s=ek2
       sk1s=sk1
       sk2s=sk2
c --------------------
       do ist=1,ncstat
        xcl2(ist)=0.
        rmoy(ist)=0.
        rig(ist)=0.
        baims(ist)=baim
       enddo
c -------------------
c ---- nsector: number of sectors in the bending magnet
      do nsec=1,nsector
c       write(6,*) '********'
c       write(6,*) ' BENDING MAGNET sector ',nsec
       devi=devis
       devr=devrs
       xlsy=devi*rmo
       sdavtot=sdavtot+xlsy
       if(nsector.gt.1) then
         if(nsec.eq.1) then
          pent1=pent1s
          pent2=0.
          ek1=ek1s
          ek2=ek2s
          sk1=0.
          sk2=0.
         endif
         if(nsec.eq.nsector) then
          pent1=0.
          pent2=pent2s
          ek1=0.
          ek2=0.
          sk1=sk1s
          sk2=sk2s
         endif
         if((nsec.gt.1).and.(nsec.lt.nsector)) then
          pent1=0.
          pent2=0.
          ek1=0.0001
          ek2=0.
          sk1=0.0001
          sk2=0.
         endif
        endif
c ------------------------------------------------------------
c --- nsctat: number of charges in the beam
        do ist=1,ncstat
          charge(ist)=cstat(ist)
          devi=devis
          devr=devrs
          xcl1(ist)=pent1-xcl2(ist)
          xmoy(ist)=0.
          xpmoy(ist)=0.
          ymoy(ist)=0.
          ypmoy(ist)=0.
          ncs(ist)=0
          rig(ist)=0.
c*et*2013*Apr*1
          sbeta(ist)=0.
          do i=1,ngood
            if(f(9,i).eq.charge(ist)) then
              xmoy(ist)=xmoy(ist)+f(2,i)
              ymoy(ist)=ymoy(ist)+f(4,i)
              xpmoy(ist)=xpmoy(ist)+f(3,i)
              ypmoy(ist)=ypmoy(ist)+f(5,i)
              gpai=f(7,i)/xmat
              bpai=sqrt(1.-1./(gpai*gpai))
c*et*2013*Apr*1
              sbeta(ist)=sbeta(ist)+bpai
              xmco=xmat*bpai*gpai
              rip=33.356*xmco*1.e-01/f(9,i)
              rig(ist)=rip+rig(ist)
cold            rmoy(ist)=rmoy(ist)+rip
              ncs(ist)=ncs(ist)+1
            endif
          enddo
c*et*2013*Apr*1
          sbeta(ist)=sbeta(ist)/float(ncs(ist))
          xmoy(ist)=xmoy(ist)/float(ncs(ist))
          ymoy(ist)=ymoy(ist)/float(ncs(ist))
          xpmoy(ist)=xpmoy(ist)/float(ncs(ist))
          ypmoy(ist)=ypmoy(ist)/float(ncs(ist))
          rig(ist)=rig(ist)/float(ncs(ist))
cold         rmoy(ist)=rig(ist)/baim
          rmoy(ist)=rig(ist)/baims(ist)
c eq.14
          ctan=cos(devi-pent2)/sin(devi-pent2)
          xep=rmo*(sin(devi)*ctan-cos(devi))
          xepc=xep+rmo-rmoy(ist)+xmoy(ist)
c eq.15 and eq.16
          argu=-xmoy(ist)*tan(pent1)/xepc
c eq.16
          omga=atan(argu)
c eq.15
          thet=omga+devi-pent2
c eq.18
          eo1=xepc/cos(omga)
c eq.17
          arg1=eo1*sin(thet)/rmoy(ist)
          eta=asin(arg1)
c eq.13
          xeb1=xepc*cos(thet)/cos(omga)+rmoy(ist)*cos(eta)
          sxeb1(ist)=xeb1
c eq.12
          xk2b1=-xmoy(ist)*tan(xcl1(ist))+xeb1*sin(devi-pent2)
c eq.11  (bend angle)
          alp(ist)=asin(xk2b1/rmoy(ist))
c eq.18
          xeo1=xepc/cos(omga)
c eq.19 (angle of inclination exit)
          argu=eo1/rmoy(ist) * sin(thet)
          xcl2(ist)=asin(argu)
          sa1b1=-rmo*sin(devr)
          sa1b1=sa1b1/sin(devr-pent2)
          xsa1b1(ist)=sa1b1+sxeb1(ist)
c --- field
c      first order
          baims(ist)=baim*(1.-nb*xsa1b1(ist)/rmo)
c       second order
          rih=1./(rmo*rmo)
          baims(ist)=baims(ist)+xb*rih*xsa1b1(ist)*xsa1b1(ist)
c -----------------------------------------------------
c ----  Transport matrix
          sbet=sbeta(ist)
          devi=alp(ist)
          AILONG=devi*rmoy(ist)
          WRITE(16,101) charge(ist),nsec,nsector,baims(ist),
     *         xsa1b1(ist),rmoy(ist),devi*180./pi,AILONG,rig(ist)
101       FORMAT(/,'  **************************************',/,
     *     '  *CENTRAL TRAJECTORY for charge: ',f4.1,' *',/,
     *     '  **************************************',/,
     *     '  SECTOR: ',i4,' SECTORS NUMBER: ',i5, /,
     *     '  BENDING FIELD:   ',e12.5,' kG  at: ',E12.5,' cm',/,
     *     '  BENDING RADIUS:  ',E12.5,' CM ',/,
     *     '  BENDING ANGLE:   ',E12.5,' DEG',/,
     *     '  length: ',e12.5,' cm  rigidity: ',E12.5,' kG.cm')
          L=AILONG
          H=1./rmoy(ist)
C   ENTRANCE FACE OF THE BENDING MAGNET
C   CLEAR R AND T
          CALL CLEAR
          GAP=APB(1)
          BE= xcl1(ist)
          LAYL =EK1
          LAYX =EK2
          RABT=0.
          IF(ABS(RAB1).GT.6.*0)RABT=1./RAB1
c -----------------------------------------
          gcog=0.
          nii=0
          do ii=1,ngood
            if(f(9,ii).eq.charge(ist)) then
              gcog=gcog+f(7,ii)/xmat
              nii=nii+1
            endif
          enddo
          gcog=gcog/float(nii)
          bcog=sqrt(1.-1./(gcog*gcog))
          fdtot=0.
          DO II=1,NGOOD
            if(f(9,ii).eq.charge(ist)) then
              tbe=tan(be)
              f(2,ii)=f(2,ii)-xmoy(ist)
              gpai=f(7,ii)/xmat
              bpai=sqrt(1.-1./(gpai*gpai))
              f(6,ii)=f(6,ii)+xmoy(ist)*tbe/(bpai*vl)
              fd(ii)=(gpai*bpai)/(gcog*bcog)
              fdtot=fdtot+fd(ii)
            endif
          enddo
c  TEST
          fdtot=fdtot/float(nii)-1.
c --------------------------------------------------
          CALL POFAR1(GAP)
          write(16,4502) be*180./pi,charge(ist)
4502      format('  ****INPUT FACE*** SLOPE: ',e12.5,' deg ',
     *           'CHARGE: ',f4.1)
          call matrix
          XLL=0.
          DO II=1,NGOOD
            if(f(9,ii).eq.charge(ist)) then
              CALL COBEAM(II,XLL)
            endif
          enddo
c *******************************************************
C  WEDGE BENDING MAGNET
C   CLEAR R AND T
          CALL CLEAR
          CALL BENMAG(sbet,fdtot)
C ---  :print the transport matrix
          write(16,4101) charge(ist)
4101      format('  ****BENDING MAGNET for charge ',f4.1)
          call matrix
          r51=r(5,1)
C  ----  transport of particles
          DO II=1,NGOOD
            if(f(9,ii).eq.charge(ist)) then
              CALL COBEAM(II,L)
            endif
C      synchrotron radiation (only for electrons, i.e. erest = 0.511 Mev)
            IF(IRAYSH.and.xmat.eq.0.511) CALL SYROUT(II)
          ENDDO
c *******************************************************
C --- EXIT FACE OF THE BENDING MAGNET
C   CLEAR R AND T
          CALL CLEAR
          xll=0.
          GAP=APB(2)
          BE= xcl2(ist)
          LAYL =SK1
          LAYX =SK2
          RABT=0.
          IF(ABS(RAB2).GT.1.E-10)RABT=1./RAB2
          CALL POFAR2(GAP)
C ---- :print the transport matrix
          write(16,4501) be*180./pi,charge(ist)
4501      format('  ****EXIT FACE*** SLOPE: ',e12.5,' deg CHARGE: ',
     *          f4.1)
          call matrix
C     transport of particles
          XLL=0.
          DO II=1,NGOOD
            if(f(9,ii).eq.charge(ist)) then
              CALL COBEAM(II,XLL)
            endif
          ENDDO
c -------------------------------------------------------------------------
c --- get back the particles coordinates in the coordinates system of the synchronous particle
          sa1b1=-rmo*sin(devr)
          sa1b1=sa1b1/sin(devr-pent2)
          ttt=xcl2(ist)-pent2
          ttt=ttt*1.e03
          do ii=1,ngood
            if(f(9,ii).eq.charge(ist)) then
c EQ.29
              a1b1=sa1b1+sxeb1(ist)
c EQ.28
              f(2,ii)=(a1b1+f(2,ii)/cos(xcl2(ist))) * cos(pent2)
c EQ.30
              f(3,ii)=f(3,ii)-ttt
c EQ.31
              gpai=f(7,ii)/xmat
              bpai=sqrt(1.-1./(gpai*gpai))
              f(6,ii)=f(6,ii)+r51*xmoy(ist)/(bpai*vl)
ccc              f(6,ii)=f(6,ii)-r51*xmoy(ist)/(bpai*vl)
c **** allow plotting the beam after the sector number nsprint in file 13 (see the MAIN)
ccc            if(nsec.eq.nsprint) then
ccc       write(13,2587) nsec,ist,charge(ist),f(2,ii),f(3,ii),f(4,ii),
ccc     *                   f(5,ii),xmoy(ist),ymoy(ist),rmoy(ist)
ccc2587   format(2x,i3,2x,i3,8(2x,e12.5))
ccc            endif
c ***************************************************************************
            endif
          enddo
c   enddo for ist (number of charges in the beam)
        enddo
c --------------------------------------------------------------
c  Space charge computation
        if(ichaes) then
c --- check the parity of nsec
          pnsec=float(nsec)/2.-nsec/2
c    nsec is odd:---> space charge computation
          if((pnsec.ne.0.).and.(nsec.lt.nsector)) then
            call cesp(scl)
            write(16,*) ' space charge after sector: ',nsec
          endif
        endif
c ------------------------------------------------------------
C    synchronous radiation (only for electrons i.e. erest = 0.511 Mev)
        IF(IRAYSH.and.xmat.eq.0.511) CALL SYREF
c     The routine SYREF changes vref and tref (reference)
c  envelope
        CALL STAPL(sdavtot*10.)
c  enddo for nsec (sectors numbers)
       enddo
c ------------------------------------------------------
c  random error in alignment
       if(ialin) call randali
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
       call disp
       call cogetc
c --- Test window after the bending magnet
       tcog=0.
       gcog=0.
       do i=1,ngood
         gcog=gcog+f(7,i)/xmat
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       gcog=gcog/float(ngood)
       bcog=sqrt(1.-1./(gcog*gcog))
       wcg=(gcog-1.)*xmat
c devtot: total angle of the bending magnet(computation of tref)
       devi=devtot
       ailong=devi*rmo
       IF(IRAYSH.and.xmat.eq.0.511) go to 2561
       TREF=TREF+AILONG/VREF
c    ----- window control
       call reject(ilost)
c ---------------------------------------------------
       dav1(idav,37)=ngood
       if(itvol) ttvols=tref
2561   continue
       BEREF=VREF/VL
       GAMREF=1./SQRT(1.-BEREF*BEREF)
       XMCO=XMAT*BEREF*GAMREF
       BORO=33.356*XMCO*1.E-01/QST
       TLONG = TREF*VREF
       WRITE(16,256) BEREF,GAMREF,TREF,TLONG,BORO,NGOOD
256    FORMAT(//,3X,' *** REFERENCE AT THE EXIT :',/,
     X '   BETA :',E12.5,' GAMMA :',E12.5,/,'   T.O.F (SEC): ',E12.5,
     X '   T.O.F (CM): ',E12.5,/,
     X '   RIGIDITY(KGAUSS.CM) :',E12.5,/,
     X '   NUMBER OF PARTICLES :',I6,/)
       if(itvol) write(16,*) ' tof for adjustments: ',ttvols,' sec'
       if(iemgrw) call emiprt(0)
cold       CALL STAPL(davtot*10.)
       RETURN
       END
      SUBROUTINE deflect(fdtot)
c ----------------------------------------------------------------
c Electrostatic deflector: Transport matrix
C     first order transport matrix (M,RD)
c     drad: horizontal deflector radius (cm)
c     l: length of the central trajectory (cm)
c     kx2, ky2: matrix terms arguments (cm-2)
c     avb: average relativistic beta
c        in this routine drad, l kx2, ky2 are converted to m
c -------------------------------------------------------------------
c*et*2013*Mar
c changed matrix elements from F.Hinterberger to Triumf note TRI-DN-05-7
c
      IMPLICIT REAL*8(A-Z)
      COMMON /BLOC11/ R(6,6), T(6,6,6)
      common /edef/avb,drad,kx2,ky2,l
c --- convert drad, kx2 , ky2 in m
      h=100./drad
      h2=h*h
      AL=L*1.E-02
      KX2=KX2*1.e04
      KY2=KY2*1.e04
      KX=SQRT(ABS(KX2))
      KY=SQRT(ABS(KY2))
      ARGX=KX*AL
      ARGY=KY*AL
      avb2=avb*avb
      avg=sqrt(1.-avb2)
      avg=1./avg
      avg2=avg*avg
      dx=(2.-avb2)*h/kx
      dxp=(2.-avb2)*h/kx2
c ***** variante:(???)
c      dx=h/kx
c      dxp=h/kx2
c ************************
c  kx2 < 0
      IF (kx2.lt.6.*0) then
       CX=COSH(ARGX)
       SX=SINH(ARGX)/KX
       SXP=sinh(argx)*kx
c ---------------------------------
C    First order Matrix R (plane (X,XP)
        R(1,1)=CX
        R(1,2)=SX
        R(1,6)=DXP*(1.-cx)
        R(2,1)=SXP
        R(2,2)=CX
c*et*2013*Mar
c        R(2,6)=DX*SX
        R(2,6)=DX*SX*KX
c*et*2013*Mar
c        R(5,1)=-H*SX
        R(5,1)=-DX*SX*KX
c*et*2013*Mar
c        R(5,2)=(CX-1.)*H/KX2
        R(5,2)=-DXP*(1.-CX)
        R(5,5)=1.
c*et*2013*Mar
c        R(5,6)=AL/avg2-DXP*H*(AL-SX)
C        R(5,6)=AL/avg2-(2.-avb2)*DXP*H*(AL-SX)
C ALT 5/28
        R(5,6)=fdtot*AL/avg2-(2.-avb2)*DXP*H*(AL-SX)
        R(6,6)=1.
      ENDIF
c  kx2 > 0
      IF (kx2.gt.6.*0) then
       CX=COS(ARGX)
       SX=SIN(ARGX)/KX
       SXP=sin(argx)*kx
C    First order Matrix R (plane (X,XP)
        R(1,1)=CX
        R(1,2)=SX
        R(1,6)=DXP*(1.-cx)
        R(2,1)=-SXP
        R(2,2)=CX
c*et*2013*Mar
c        R(2,6)=DX*SX
        R(2,6)=DX*SX*KX
c*et*2013*Mar
c        R(5,1)=-H*SX
C        R(5,1)=-DX*SX*KX
C ALT 5/28
        R(5,1)=DX*SX*KX
c*et*2013*Mar
c        R(5,2)=(CX-1.)*H/KX2
C        R(5,2)=-DXP*(1.-CX)
C ALT 5/28
        R(5,2)=DXP*(1.-CX)
        R(5,5)=1.
c*et*2013*Mar
c        R(5,6)=AL/avg2-DXP*H*(AL-SX)
C        R(5,6)=AL/avg2-(2.-avb2)*DXP*H*(AL-SX)
C ALT 5/28
        R(5,6)=fdtot*AL/avg2-(2.-avb2)*DXP*H*(AL-SX)
        R(6,6)=1.
      ENDIF
c  kx2 = 0
      IF (kx2.eq.6.*0) then
        R(1,1)=1.
        R(1,2)=AL
        R(1,6)=0.
        R(2,1)=0.
        R(2,2)=1.
c*et*2013*Mar
c        R(2,6)=L*H
        R(2,6)=L*H*(2.-avb2)
c*et*2013*Mar
c        R(5,1)=-L*H
        R(5,1)=-L*H*(2.-avb2)
        R(5,2)=0.
        R(5,5)=1.
        R(5,6)=AL/avg2
      ENDIF
c ky2 < 0
      IF (ky2.lt.6.*0) then
       CY=COSH(ARGY)
       SY=SINH(ARGY)/KY
       SYP=sinh(argy)*ky
       R(3,3)=CY
       R(3,4)=SY
       R(4,3)=SYP
       R(4,4)=CY
      ENDIF
c ky2 > 0
      IF (ky2.gt.6.*0) then
       CY=COS(ARGY)
       SY=SIN(ARGY)/KY
       SYP=sin(argy)*ky
       R(3,3)=CY
       R(3,4)=SY
       R(4,3)=-SYP
       R(4,4)=CY
      ENDIF
c  ky2 = 0
      IF (ky2.eq.6.*0) then
       CY=1.
       SY=AL
       SYP=0.
       R(3,3)=CY
       R(3,4)=SY
       R(4,3)=SYP
       R(4,4)=CY
      ENDIF
      return
      end
       SUBROUTINE e_deflec
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       COMMON/TAPES/IN,IFILE,META
       COMMON/DYN/TREF,VREF
       COMMON/ERIGID/edr0
       common/faisc/f(10,iptsz),imax,ngood
       common/femt/iemgrw,iemqesg
       LOGICAL IEMGRW
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/QMOYEN/QMOY
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       common /edef/avb,drad,kx2,ky2,L
       real*8 L,kx2,ky2
       COMMON /BLOC11/ R(6,6), T(6,6,6)
       COMMON/PLTPRF/SPRFX(3000),SPRFY(3000),SPRFL(3000),SPRFW(3000),
     X               SPRFP(3000),SPRNG(3000),IPRF
       common/rander/ialin
       common/cgtof/charm(20),cgtdv(20),nbch(20),netac
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/itvole/itvol,imamin
       common/tofev/ttvols
       COMMON/SECDR/ISEOR
       LOGICAL ISEOR,sseor
       common/mcs/imcs,ncstat,cstat(20)
       logical itvol,imamin,ichaes
       character*1 cr
       logical ialin
       dimension xmoy(20),ymoy(20),rmoy(20),rig(20),ncs(20)
       dimension xpmoy(20),ypmoy(20),avbt(20),charge(20),alp(20)
c -----------------------------------------------------------------
c print out on terminal of transport element # on one and the same line
       NRTRE=NRTRE+1
       cr=char(13)
#ifdef F77
       WRITE(6,8254) NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') NRTRE,NRRES,cr
8254   format('Transport element:',i5,
     *        '      Accelerating element:',i5,a1)
#endif
c ------------------------------------------------------------------
c   characteristics of the deflector (central trajectory)
c   input parameters
c   nsector: nombre of sectors in the deflector
c   rm0: radial radius (cm)
c   devtot: bend angle (deg)
c   radii: vertical (radial) radii of curvature (cm)
c   elecf: nominal electric field
c   characteristics of the synchronous particle
c   wt0 :total energy (MeV)
c   qst: electric charge
c ----------------------------------------------------------------------
c  ---- save iseor in sseor
        sseor=iseor
c   iseor = false ---> second ordre transport matrix not available for the deflector
       iseor=.false.
c
       read(in,*) nsector
       read(in,*) rm0,devtot,radii,elecf
       bd0=vref/vl
       gd0=sqrt(1.-bd0*bd0)
       gd0=1./gd0
       wt0=gd0*xmat
       wt0c=(wt0-xmat)
c   edr0: electric rigidity (kV)
       edr0=wt0*bd0*bd0*1.e03/qst
c     efd0: radial electric field (in kV/cm)
       efd0=edr0/rm0
       edfnom=efd0
       if(elecf.ge.0.) efd0=elecf
       WRITE(16,1020)RM0,DEVTOT,RADII,EDR0,EDFNOM,EFD0,WT0C
1020    FORMAT('    ELECTROSTATIC DEFLECTOR************',/,
     2 '  BENDING RADIUS:  ',E12.5,' cm ',/,
     2 '  BEND ANGLE:      ',E12.5,' deg',/,
     2 '  VERTICAL RADII OF CURVATURE: ',E12.5,' cm',/,
     2 '  RIGIDITY: ',e12.5,' kV ',/,
     2 '  RADIAL ELECTRIC FIELD (nominal): ',e12.5,' kV/cm',/,
     2 '  RADIAL ELECTRIC FIELD (applied): ',e12.5,' kV/cm',/,
     2 '  INPUT ENERGY: ',e12.5,' MeV',/)
       if(ichaes) then
        write(16,*)'***** beam current: ',beamc,' mA'
        if((iscsp.lt.3).and.(ncstat.gt.1)) then
         write(6,*) '****************************'
         write(6,2748)
2748     format(' CAUTION: In the case of multiple charge states',/,
     *  ' HERSC and SCHERM can not be used for electrostatic bends')
        endif
       endif
C  start prints in file 'short.data'
       findex=1.+rm0/radii
       idav=idav+1
       iitem(idav)=21
       dav1(idav,2)=devtot
       dav1(idav,3)=rm0*10.
       dav1(idav,5)=radii*10.
       dav1(idav,6)=findex
       dav1(idav,7)=edr0
       dav1(idav,8)=efd0*0.1
c   convert bend angle in rad
       devtot=devtot*pi/180.
       l=devtot*rm0
       dav1(idav,1)=l*10.
       sdavtot=davtot
       davtot=davtot+l
       dav1(idav,4)=davtot*10.
c --- space charge computation: nsector must be g.t. 1
       if(ichaes.and.(nsector.eq.1) ) nsector = 2
       devi=devtot/float(nsector)
       devr=devtot/float(nsector)
c ----------------------------------------------------------
c -- scl effective length for space charge computation
c     scl is the path length of the C.T. over two following sectors
       scl=2.*devi*rm0
       do ist=1,ncstat
        rmoy(ist)=0.
        rig(ist)=0.
c*et*2015-May-13-s
        alp(ist)=devi
c*et*2015-May-13-e
       enddo
c -------------------
c ---- nsector: number of sectors in the deflector
      do nsec=1,nsector
        write(6,*) '********'
        write(6,*) ' deflector sector ',nsec
        xlsy=devi*rm0
        sdavtot=sdavtot+xlsy
c ------------------------------------------------------------
c --- nsctat: number of charge states in the beam
        do ist=1,ncstat
          charge(ist)=cstat(ist)
          xmoy(ist)=0.
          xpmoy(ist)=0.
          ymoy(ist)=0.
          ypmoy(ist)=0.
          ncs(ist)=0
          rig(ist)=0.
          avbt(ist)=0.
          do i=1,ngood
            if(f(9,i).eq.charge(ist)) then
              xmoy(ist)=xmoy(ist)+f(2,i)
              ymoy(ist)=ymoy(ist)+f(4,i)
              xpmoy(ist)=xpmoy(ist)+f(3,i)
              ypmoy(ist)=ypmoy(ist)+f(5,i)
              gpai=f(7,i)/xmat
              bpai=sqrt(1.-1./(gpai*gpai))
              avbt(ist)=avbt(ist)+bpai
c      electric rigidity (kV)
              rip=f(7,i)*bpai*bpai/f(9,i) *1.e03
              rig(ist)=rip+rig(ist)
              ncs(ist)=ncs(ist)+1
            endif
          enddo
          xmoy(ist)=xmoy(ist)/float(ncs(ist))
          ymoy(ist)=ymoy(ist)/float(ncs(ist))
          xpmoy(ist)=xpmoy(ist)/float(ncs(ist))
          ypmoy(ist)=ypmoy(ist)/float(ncs(ist))
          rig(ist)=rig(ist)/float(ncs(ist))
          rmoy(ist)=rig(ist)/efd0
          avbt(ist)=avbt(ist)/float(ncs(ist))
c   dispersion in dp/p relative to the cog of the bunch
          gcog=sqrt(1.-avbt(ist)*avbt(ist))
C ALT 5/28
          fdtot=0.
          nii=0
C ALT --
          do i=1,ngood
            if(f(9,i).eq.charge(ist)) then
              gpai=f(7,i)/xmat
              bpai=sqrt(1.-1./(gpai*gpai))
              fd(i)=(gpai*bpai)/(gcog*avbt(ist))
C ALT 5/28
              fdtot=fdtot+fd(i)
              nii=nii+1
C ALT --
            endif
          enddo
C ALT 5/28
          fdtot=fdtot/float(nii)
c --- local deflector
c       alp(ist): angle of the local deflector
c        parametres kx2 and ky2
c  eq.11
         oo1=rm0-rmoy(ist)+xmoy(ist)
c  eq.10
c*et*2015-May-13-s
c         abet=oo1*sin(devi)/rmoy(ist)
         abet=oo1*sin(alp(ist))/rmoy(ist)
c*et*2015-May-13-e
         abet=asin(abet)
c  eq.12 (angle of the local central trajectory)
c*et*2015-May-13-s
c         alp(ist)=devi+abet
         alp(ist)=alp(ist)+abet
c*et*2015-May-13-e
c  eq.14 (field index)
         findex=1.+rmoy(ist)/radii
c  eq.15 (parameters kx, ky)
         kx2=3.-findex-avbt(ist)*avbt(ist)
         rmoy2=rmoy(ist)*rmoy(ist)
         kx2=kx2/rmoy2
         ky2=(findex-1.)/rmoy2
c -----------------------------------------------------
c ----  Transport matrix
         devi=alp(ist)
         AILONG=devi*rmoy(ist)
         L=AILONG
c     drad: horizontal deflector radius (cm)
c     l: length of the central trajectory (cm)
c     kx2, ky2: DIMENSIONLESS coefficients depending on the field indice
c     avb: average relativistic beta
       avb=avbt(ist)
       drad=rmoy(ist)
C --- deflector matrix
         WRITE(16,4101) charge(ist),nsec,nsector,efd0,rig(ist),
     *         findex,kx2,ky2,rmoy(ist),devi*180./pi,AILONG
4101      FORMAT(/,'  **************************************',/,
     *     '  *CENTRAL TRAJECTORY for charge: ',f4.1,' *',/,
     *     '  **************************************',/,
     *     '  SECTOR: ',i4,' SECTORS NUMBER: ',i5, /,
     *     '  RADIAL FIELD:   ',e12.5,' kV*cm-1: ',/,
     *     '  RIGIDITY:   ',e12.5,' kV ',/,
     *     '  FIELD INDEX: ',e12.5,' PARAMETER Kx: ',e12.5,
     *     ' cm-2  PARAMETER Ky: ',e12.5,' cm-2',/,
     *     '  BENDING RADIUS:  ',E12.5,' cm ',
     *     '  BENDING ANGLE:   ',E12.5,' deg',/,
     *     '  LENGTH: ',e12.5,' cm',/)
C   CLEAR R AND T
       CALL CLEAR
C ALT 5/28
       CALL deflect(fdtot)
C CALL deflect
C ---  :print the transport matrix
        call matrix
        r51=r(5,1)
C  ----  transport of particles
        DO II=1,NGOOD
          if(f(9,ii).eq.charge(ist)) then
            CALL COBEAM(II,L)
          endif
         enddo
c -------------------------------------------------------------------------
c --- get back the particles coordinates in the coordinates system of the synchronous particle
c   eq.33
       ec=-rmoy(ist)*cos(abet)-oo1*sin(devi)+rm0
       do ii=1,ngood
         if(f(9,ii).eq.charge(ist)) then
c   eq.36
           f(2,ii)=f(2,ii)*cos(abet)-ec
c   eq.41
           f(3,ii)=f(3,ii)-abet*1.e03
c   eq.45
           gpai=f(7,ii)/xmat
           bpai=sqrt(1.-1./(gpai*gpai))
           f(6,ii)=f(6,ii)+r51*xmoy(ist)/(bpai*vl)
         endif
        enddo
c   enddo for ist (number of charges in the beam)
       enddo
c --------------------------------------------------------------
c  Space charge computation
         if(ichaes) then
c --- check the parity of nsec
          pnsec=float(nsec)/2.-nsec/2
c    nsec is odd:---> space charge computation
          if((pnsec.ne.0.).and.(nsec.lt.nsector)) then
            call cesp(scl)
            write(6,*) ' space charge after sector: ',nsec
          endif
        endif
c ------------------------------------------------------------
c  enveloppe
        CALL STAPL(sdavtot*10.)
c  enddo for nsec (sectors numbers)
       enddo
c ------------------------------------------------------
c  random error in alignment
       if(ialin) call randali
c     Change the dispersion dE/E with respect to the C.O.G of the bunch
       call disp
       call cogetc
c --- Test window after the bending magnet
       tcog=0.
       gcog=0.
       do i=1,ngood
         gcog=gcog+f(7,i)/xmat
         tcog=tcog+f(6,i)
       enddo
       tcog=tcog/float(ngood)
       gcog=gcog/float(ngood)
       bcog=sqrt(1.-1./(gcog*gcog))
       wcg=(gcog-1.)*xmat
c devtot: total angle of the deflector (computation of tref)
       devi=devtot
       ailong=devi*rm0
       TREF=TREF+AILONG/VREF
c    ----- window control
       call reject(ilost)
c ---------------------------------------------------
       dav1(idav,36)=ngood
       if(itvol) ttvols=tref
       if(itvol) write(16,*) ' tof for adjustments: ',ttvols,' sec'
       if(iemgrw) call emiprt(0)
cold       CALL STAPL(davtot*10.)
c  ---- restore iseor (from sseor)
        iseor=sseor
       RETURN
       END
c  SPECIAL SHEFF (bunches are separated)*******************************
       SUBROUTINE cesp(xlqua)
c   ....................................................................
C       select the space charge method (optical lenses)
c   ....................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/CDEK/DWP(iptsz)
       COMMON/DCSPA/IESP
       common/faisc/f(10,iptsz),imax,ngood
       common/mcs/imcs,ncstat,cstat(20)
       logical ichaes,iesp,isepa
       DO I=1,ngood
         DWP(I)=0.
       enddo
C      Space charge
       SCDIST=0.
       IF(.NOT.ICHAES) RETURN
C       XLQUA: length of space charge effect
        if((iscsp.ne.3).and.(ncstat.gt.1)) then
         write(6,*) '****************************'
         write(6,2748)
2748     format(' ERROR: Wrong space charge model chosen',/,
     *          ' With multiple charge states in the beam',/,
     *          ' only the SCHEFF routine should be used')
         write(16,2748)
         STOP
        endif
       SCDIST=XLQUA
       write(16,*) 'space charge length(cm): ',scdist
       iesp=.true.
       if(iscsp.le.1) then
         ini=1
         call hersc(ini)
         ini=2
         call hersc(ini)
       endif
       if(iscsp.eq.2)  call schermi
       if(iscsp.eq.3) then
        if(ncstat.eq.1) call scheff1(1)
c --- otherwise: ncstat > 1 check if the bunches are separated or not
        if(ncstat.gt.1) then
          isepa=.false.
          call b_sep(isepa)
c isepa = true  call special scheff --->scheff_sep
c isepa = false call usual scheff ----> scheff1(1)
          if(isepa) call scheff_sep
          if(.not.isepa) call scheff1(1)
        endif
       endif
       RETURN
       END
       SUBROUTINE sizer(ist,xrms,yrms,zrms)
c   ........................................................................
c   partial R.M.S. (called by SCHEFF_sep)
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/CGRMS/xsum,ysum,zsum
       common/faisc/f(10,iptsz),imax,ngood
       common/mcs/imcs,ncstat,cstat(20)
       xsum=0.
       ysum=0.
       zsum=0.
       xsqsum=0.
       ysqsum=0.
       zsqsum=0.
       ngist=0
       do i=1,ngood
        if(f(9,i).eq.cstat(ist)) then
         ngist=ngist+1
         xsum=xsum+xc(i)
         ysum=ysum+yc(i)
         zsum=zsum+zc(i)
         xsqsum=xsqsum+xc(i)*xc(i)
         ysqsum=ysqsum+yc(i)*yc(i)
         zsqsum=zsqsum+zc(i)*zc(i)
        endif
       enddo
       xsum=xsum/float(ngist)
       ysum=ysum/float(ngist)
       zsum=zsum/float(ngist)
       xsqsum=xsqsum/float(ngist)
       ysqsum=ysqsum/float(ngist)
       zsqsum=zsqsum/float(ngist)
       xrms=SQRT(xsqsum-xsum*xsum)
       yrms=SQRT(ysqsum-ysum*ysum)
       zrms=SQRT(zsqsum-zsum*zsum)
       RETURN
       END
       SUBROUTINE pintim1(ist)
c   ..................................................................
C   Shifts particle coordinates to a single point in time. Uses
C   a linear shift
C   Divide by 100. to convert from cm to meters
c    called by SCHEFF or SCHERM
c   ..................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/QMOYEN/QMOY
       COMMON/PART/XC(iptsz),YC(iptsz),ZC(iptsz)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       COMMON/AZLIST/ICONT,IPRIN
       common/mcs/imcs,ncstat,cstat(20)
       grmoy=0.
       trmoy=0.
       xbax=0.
       ngist=0
       do i=1,ngood
        if(f(9,i).eq.cstat(ist)) then
         ngist=ngist+1
         grmoy=grmoy+f(7,i)/xmat
         trmoy=trmoy+f(6,i)
         xbax=xbax+f(2,i)
        endif
       enddo
       trmoy=trmoy/float(ngist)
       grmoy=grmoy/float(ngist)
       brmoy=sqrt(1.-1./(grmoy*grmoy))
       xbax=xbax/float(ngist)
cccc       apl=0.
c  Isochronism correction  (bending magnet) only with SCHERM
c    does not work with  with SCHEFF  (iscsp=3)
cccc       if(iscsp.eq.2) then
cccc         xb2x=0.
cccc         xb2z=0.
cccc         xbxz=0.
cccc         do np=1,ngood
cccc           gpai=f(7,np)/xmat
cccc           bpai=sqrt(1.-1./(gpai*gpai))
cccc           zc(np)=(trmoy-f(6,np))*bpai*vl/100.
cccc           xc(np)=(f(2,np)-xbax)/100.
cccc           xb2z=xb2z+zc(np)*zc(np)
cccc           xb2x=xb2x+xc(np)*xc(np)
cccc           xbxz=xbxz+zc(np)*xc(np)
cccc         enddo
cccc         xb2z=xb2z/float(ngood)
cccc         xb2x=xb2x/float(ngood)
cccc         xbxz=xbxz/float(ngood)
cccc         apl=atan(-2.*xbxz/(xb2x-xb2z))/2.
cccc         write(16,*) 'slope of the bunch in plane(Oz,Ox):',apl,' radian'
cccc       endif
       do np=1,ngood
        if(f(9,np).eq.cstat(ist)) then
         gpai=f(7,np)/xmat
         bpai=sqrt(1.-1./(gpai*gpai))
c      iscsp = 3 Lorentz transformation (only with scheff)
comment         if(iscsp.eq.3) znp=(trmoy-f(6,np))*bpai*vl*grmoy
comment         if(iscsp.eq.2) znp=(trmoy-f(6,np))*bpai*vl
         znp=(trmoy-f(6,np))*bpai*vl
         xnp=f(2,np)
         zc(np)=znp*cos(apl)+xnp*sin(apl)
         xnp=xnp*cos(apl)-znp*sin(apl)
C        convert from mrad to rad
         f3=f(3,np)*1.e-03
         f5=f(5,np)*1.e-03
C        convert from cm   to m
         xc(np)=(xnp+zc(np)*f3)/100.
         yc(np)=(f(4,np)+zc(np)*f5)/100.
         zc(np)=zc(np)/100.
        endif
       enddo
       xbar=0.
       ybar=0.
       zbar=0.
       do np=1,ngood
c      evaluate xbar , ybar , zbar
        if(f(9,np).eq.cstat(ist)) then
         xbar=xbar+xc(np)
         ybar=ybar+yc(np)
         zbar=zbar+zc(np)
        endif
       enddo
       xbar=xbar/float(ngist)
       ybar=ybar/float(ngist)
       zbar=zbar/float(ngist)
c  Translate distribution by center of mass coordinates to shift
c  coordinate origin to (0,0,0)
       do np=1,ngood
        if(f(9,np).eq.cstat(ist)) then
         xc(np)=xc(np)-xbar
         yc(np)=yc(np)-ybar
         zc(np)=zc(np)-zbar
        endif
       enddo
       return
       end
       SUBROUTINE b_sep(isepa)
c   ........................................................................
c     seek the bunches in the beam separated or not
c   ........................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/CGRMS/xsum,ysum,zsum
       common/faisc/f(10,iptsz),imax,ngood
       common/mcs/imcs,ncstat,cstat(20)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       dimension d11(2),d22(2),d12(2),rp12(2),xpint(2),xint(2)
       logical isepa
       csmax1=0.
       icms1=0
       do i=1,ncstat
        if(cstat(i).gt.csmax1) then
         icsm1=i
         csmax1=cstat(i)
        endif
       enddo
       csmax2=0.
       do i=1,ncstat
        if(i.ne.icsm1) then
         if(cstat(i).gt.csmax2) then
         csmax2=cstat(i)
         endif
        endif
       enddo
       xg1=0.
       xpg1=0.
       xg2=0.
       xpg2=0.
       imax1=0
       imax2=0
c   COG over state charges csmax1 and csmax2
       do i=1,ngood
        if(f(9,i).eq.csmax1) then
         xg1=xg1+f(2,i)
         xpg1=xpg1+f(3,i)
         imax1=imax1+1
        endif
        if(f(9,i).eq.csmax2) then
         xg2=xg2+f(2,i)
         xpg2=xpg2+f(3,i)
         imax2=imax2+1
        endif
       enddo
       xg1=xg1/float(imax1)
       xpg1=xpg1/float(imax1)
       xg2=xg2/float(imax2)
       xpg2=xpg2/float(imax2)
c ------------------------------------
       d11(1)=0.
       d22(1)=0.
       d12(1)=0.
       d11(2)=0.
       d22(2)=0.
       d12(2)=0.
       do i=1,ngood
        if(f(9,i).eq.csmax1) then
         d11(1)=d11(1)+(f(3,i)-xpg1)**2
         d22(1)=d22(1)+(f(2,i)-xg1)**2
         d12(1)=d12(1)+(f(3,i)-xpg1)*(f(2,i)-xg1)
        endif
        if(f(9,i).eq.csmax2) then
         d11(2)=d11(2)+(f(3,i)-xpg2)**2
         d22(2)=d22(2)+(f(2,i)-xg2)**2
         d12(2)=d12(2)+(f(3,i)-xpg2)*(f(2,i)-xg2)
        endif
       enddo
       d11(1)=d11(1)/float(imax1)
       d22(1)=d22(1)/float(imax1)
       d12(1)=d12(1)/float(imax1)
       d11(2)=d11(2)/float(imax2)
       d22(2)=d22(2)/float(imax2)
       d12(2)=d12(2)/float(imax2)
       rp12(1)=d12(1)/sqrt(d11(1)*d22(1))
       rp12(2)=d12(2)/sqrt(d11(2)*d22(2))
       xpint(1)=sqrt(d11(1)*(1.-rp12(1)))
       xint(1)=sqrt(d22(1)*(1.-rp12(1)))
       xpint(2)=sqrt(d11(2)*(1.-rp12(2)))
       xint(2)=sqrt(d22(2)*(1.-rp12(2)))
c ------------------------------------------------------------
       elip1=xpg1+xpint(1)
       elip2=xpg2-xpint(2)
       if(elip1.lt.elip2)isepa=.true.
c TEST*********
cold       write(6,*)'xpg1 xpint(1) elip1 ',xpg1,xpint(1),elip1
cold       write(6,*)'xpg2 xpint(2) elip2 ',xpg2,xpint(2),elip2
cold       write(6,*)'isepa ',isepa
cold       write(6,*) ' ************************************'
c ***************************************************************
       RETURN
       END
      SUBROUTINE schefini
c  ..............................................................................
c    SCHEFINI set up field tables for SCHEFF1 and SCHEFF_sep
c ...............................................................................
c     input data
c         sce(2)=radial extension in rms multiples
c         sce(3)=longitudinal extension in rms multiples
c         sce(4)=no. of radial mesh intervals (le 20)
c         sce(5)=no. of longitudinal mesh intervals (le 40)
c         sce(6)=no. of adjacent bunches, applicable for buncher studies
c                and should be 0 for linac dynamics
c         sce(7)=distance between adjacent beam pulses in cm (transport studies)
c                input zero to get (beta*lambda) default
c         sce(8)=desactived
c         sce(9)=option to integrate space charge forces over box
c                  if.eq.0. no integration  see sub gaus for further
c                  explanation.
c         sce(10) =1 : call in quads,solenoids,accelarating elements
c         sce(10) =2 : call in drifts,accelarating elements
c         sce(10) =3 : call at both
c  standard SCHEFF parameters (see user guide)
c         sce(2)=4
c         sce(3)=4
c         sce(4)=20
c         sce(5)=40
c         sce(6)=0
c         sce(7)=0
c         sce(9)=0
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/faisc/f(10,iptsz),imax,ngood
       LOGICAL ichaes
c  common modifiés *********************
       common/stsc/beami,wavel,freq,btazero,frrms,fzrms,nr,nz
       common/stsc1/beams,im1,im2,im3,nr1,nz1,nq
       common/fldcom/ rp, zp,pl,opt,nip
c ***********************************
       common/rcshef/sce(20)
       common/conti/irfqp
       logical irfqp
         gmoy=0.
         do np=1,ngood
           gmoy=f(7,np)/xmat+gmoy
         enddo
         gmoy=gmoy/float(ngood)
         bgmoy=sqrt(gmoy*gmoy-1.)
         beams=beamc/1000.0
         wavel=2.*pi*vl/fh
         freq=fh/(2.*pi)
         frrms=sce(2)
         fzrms=sce(3)
         nr=idint(sce(4))
         nz=idint(sce(5))
         nip=idint(sce(6))
         opt=sce(9)
         pl=bgmoy*wavel
         if(irfqp) pl=pl/2.
c         sce(7)=pulse length, if not beta lambda.(transport studies), units are cm
         if(sce(7).gt. 0.) pl=sce(7)*gmoy
         nr1=nr+1
         nz1=nz+1
         im1=nr*nz
         im2=nr1*nz1
         im3=nr1*nz
         na=1
         nb=ngood
         nq=nb-na+1
         return
         end
       SUBROUTINE scheff1(idum)
c  ..............................................................................
c    SCHEFF space charge method
c    remark: In this version int is a dummy parameter
c     This version of SCHEFF, starting from Swesson version, has modifications made to include
c     corrections for relativistic beams. The dynamics
c     have been modified to transform to the beam rest frame,
c     calculate the space-charge kicks in this frame, and then
c     transform back to the lab frame.
c     input data
c         sce(1)=beam current in ma.
c         sce(2)=radial extension in rms multiples
c         sce(3)=longitudinal extension in rms multiples
c         sce(4)=no. of radial mesh intervals (le 20)
c         sce(5)=no. of longitudinal mesh intervals (le 40)
c         sce(6)=no. of adjacent bunches, applicable for buncher studies
c                and should be 0 for linac dynamics
c         sce(7)=distance between adjacent beam pulses in cm (transport studies)
c                input zero to get (beta*lambda) default
c         sce(8)=desactived
c         sce(9)=option to integrate space charge forces over box
c                  if.eq.0. no integration  see sub gaus for further
c                  explanation.
c         sce(10) =1 : call in quads,solenoids,accelarating elements
c         sce(10) =2 : call in drifts,accelarating elements
c         sce(10) =3 : call at both
c  standard SCHEFF parameters (see sub schfdyn and user guide)
c         sce(2)=4
c         sce(3)=4
c         sce(4)=20
c         sce(5)=40
c         sce(6)=0
c         sce(7)=0
c         sce(9)=0
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/CMPTE/IELL
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/DIMENS/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/CDEK/DWP(iptsz)
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/DCSPA/IESP
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/posc/xpsc
       LOGICAL ichaes,iesp
       common/bg/bsc,gsc,phis,wsync
       common/stsc/beami,wavel,freq,btazero,frrms,fzrms,nr,nz
       common/stsc1/beams,im1,im2,im3,nr1,nz1,nq
       common/fldcom/ rp, zp,pl,opt,nip
       common/spacech1/ rm(21), zm(41), rs(20), ers(16800), ezs(16800),
     1 ez(861), aa(800), rssq(20), zzs(41)
     2 ,er(861),rss(20), ismax(40), iemax(41)
       common/rcshef/sce(20)
       common/conti/irfqp
c********************************************
c v28/04/2015
       common/fcont/ifcont
       logical ifcont
c********************************************
       logical irfqp
       idm=idum
c     set up field tables with int=0
comment       beami=beamc/1000.0
c********************************************
c v28/04/2015
       if(ifcont) then
         pib=180.
         call compress(pib)
       endif
c********************************************
       gmoy=0.
       do np=1,ngood
         gmoy=f(7,np)/xmat+gmoy
       enddo
       gmoy=gmoy/float(ngood)
       beami=beams/gmoy
       IF (beami.eq.0. .OR. scdist.eq.0.) return
       IELL=IELL+1
c       write(16,*) ' *call SCHEFF  ',iell
       call pintim
c       write(16,*) ' *after pintim '
       CALL sizrms(0,xrms,yrms,zrms,zz)
c       write(16,*) ' *after sizrms ',zz
       write(16,6875) iell,xrms,yrms,zrms
c   write the size rms in the output file 17
comment       write(17,25) iell,xrms,yrms,zrms
comment25     format(2x,i5,3(2x,e12.5))
6875   format(' Cell ',i4,' RMS size(m)',e12.5,2x,e12.5,2x,e12.5)
       rrms=sqrt(xrms*xrms+yrms*yrms)
c    change unit: m==>cm
       zrms1=zrms*100.
       rrms=rrms*100.
       dr=rrms*frrms/float(nr)
       dz=zrms1*fzrms/float(nz)
       rmax=float(nr)*dr
c          load rm, zm, rs, zs
       rm(1)=0.0
       do i=2,nr1
         rm(i)=float(i-1)*dr
         rssq(i-1)=.5*(rm(i-1)**2+rm(i)**2)
         rss(i-1)=0.5*(rm(i-1)+rm(i))
         rs(i-1)=sqrt(rssq(i-1))
       enddo
       zs=.5*dz
       do i=1,nz1
         zm(i)=float(i-1)*dz
         zzs(i)=zm(i)+zs
       enddo
       hl=float(nz)*zs
c          load ers and ezs
c     mesh dimensions are in cm. ers and ezs are in 1/cm.
c     c1, c2 and c3 are in cm., and c4 is in mev-cm.
c     q=coulombs/point.   (2/pi)*e/(4*pi*epsilon)=572167 cm mev/coul.
       q=beami/(freq*float(nq))
       c1=572167.*q/xmat
       l=0
       do 40 k=1,nr
       rfac=(rm(k+1)**2-rm(k)**2)*dz/2.
       if(opt.eq.0.) rfac=1.
       do 40 j=1,nz
       zp=zm(j+1)
       do 40 i=1,nr1
       rp=rm(i)
       if(opt.eq.0.) call flds(rs(k),zs,er1,ez1)
       if(opt.eq.0.) go to 35
       call gaus (rm(k),rm(k+1),zm(1),zm(2),opt,er1,ez1)
 35    l=l+1
       ers(l)=c1*er1/rfac
ccc       write(16,*) ' ersl ',ers(l),' er1 ',er1,' l ',l
ccc       write(16,*) ' rsk ',rs(k),' k ',k,' zs ',zs
ccc       pause
   40  ezs(l)=c1*ez1/rfac
       IF (beamc.eq.0. .OR. scdist.eq.0.) return
       dz1=scdist/100.
       dist=scdist
c sup       WRITE(16, *) ' fields acting length(cm): ',DIST
c          evaluate and apply space charge effects.
c     phimc=phi of mesh center.
C   Shifts particle coordinates to a single point in time. Uses
C   a linear shift
C      Beam c.g.
       xbar=0.
       ybar=0.
       zbar=0.
       brmoy=0.
       trmoy=0.
       do np=1,ngood
         gpai=f(7,np)/xmat
         brmoy=brmoy+sqrt(1.-1./(gpai*gpai))
         trmoy=trmoy+f(6,np)
       enddo
       trmoy=trmoy/float(ngood)
       phimc=trmoy*fh
       pbar=phimc
       beta=brmoy/float(ngood)
       gsc=1./sqrt(1.-beta*beta)
       bg=beta*gsc
       c3=dist/bg
       c4=dist*xmat
c          evaluate ng, xbar, ybar, and pbar.
       c2=beta*wavel/(2.*pi)
       gmsq=1.+bg**2
       gam=sqrt(gmsq)
c          evaluate ng, xbar, ybar
       ng=0
       xbar=0.
       ybar=0.
       xsq=0.
       ysq=0.
       do np=1,ngood
         ng=ng+1
         x=f(2,np)
         y=f(4,np)
         xf=f(3,np)
         yf=f(5,np)
         xbar=xbar+x
         ybar=ybar+y
         xsq=xsq+x**2
         ysq=ysq+y**2
       enddo
       eng=float(ngood)
       xbar=xbar/eng
       ybar=ybar/eng
c  the mesh center is phi*syn
       xsq=xsq/eng
       ysq=ysq/eng
       epsq=sqrt((xsq-xbar*xbar)/(ysq-ybar*ybar))
       epsqi=1./epsq
       xfac=2./(epsq + 1.)
       yfac=epsq*xfac
c          clear and load bins
       ng=0
       do i=1,im1
         aa(i)=0.0
       enddo
       do 120 np=1,ngood
         rsq=(f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq
c     i=sqrt(rsq)/dr+1.
         r=sqrt(rsq)
         halfdr=dr*0.5
         i=idint(r/dr+1.0)
         if (i.gt.nr) go to 120
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if (abs(z).ge.hl) go to 120
c------distribute charge among adjacent bins.
         ng=ng+1
         zz=z+hl
         jm1=idint(zz/dz+1.)
         i1=i+1
c     if (rsq.lt.rssq(i)) i1=i-1
         if (rsq.lt.rss (i)) i1=i-1
         if (i1.lt.1) i1=1
         if (i1.gt.nr) i1=nr
         j1=jm1+1
         if (zz.lt.zzs(jm1)) j1=jm1-1
         if (j1.lt.1) j1=1
         if (j1.gt.nz) j1=nz
         a=1.
c     if (i1.ne.i) a=(rsq-rssq(i1))/(rssq(i)-rssq(i1))
      if (i1.eq.i) then
        a=1.
        else
        rdr2=rsq/dr**2
        sqr=sqrt(4.*rdr2-1.)
        rminsq=(halfdr*(sqr-1.))**2
        rmaxsq=(halfdr*(sqr+1.))**2
        if (i1.lt.i) then
          a=(rmaxsq-rm(i)**2)/(rmaxsq-rminsq)
        else
          a=(rm(i1)**2-rminsq)/(rmaxsq-rminsq)
        endif
      endif
comment         if (r.gt.halfdr)then
comment           rminsq=(r-halfdr)**2
comment           rmaxsq=(r+halfdr)**2
comment           if (i1.lt.i) then
comment             a=(rmaxsq-rm(i)**2)/(rmaxsq-rminsq)
comment           else
comment             a=(rm(i1)**2-rminsq)/(rmaxsq-rminsq)
comment           endif
comment         endif
         b=1.-a
         cc=1.
         if (j1.ne.jm1) cc=(zz-zzs(j1))/(zzs(jm1)-zzs(j1))
         d=1.-cc
         k=(jm1-1)*nr+i
         aa(k)=aa(k)+a*cc
         k=k+i1-i
         aa(k)=aa(k)+b*cc
         k=(j1-1)*nr+i
         aa(k)=aa(k)+a*d
         k=k+i1-i
         aa(k)=aa(k)+b*d
  120  continue
       eng=float(ng)
c          find ismax for each j
c  v02/04/2010       do 140 j=1,nz
c  v02/04/2010         l=(j-1)*nr
c  v02/04/2010         k=nr
c  v02/04/2010         do 130 i=1,nr
c  v02/04/2010           m=l+k
c  v02/04/2010           if (aa(m)) 130,130,140
c  v02/04/2010  130    k=k-1
c  v02/04/2010  140  ismax(j)=k
      do j=1,nz
       l=(j-1)*nr
       k=nr
       do i=1,nr
        m=l+k
        if(aa(m).le.0.00) then
         k=k-1
         go to 130
        else
         go to 140
        endif
130     continue
       enddo
140       ismax(j)=k
       enddo
c        find iemax for each j
         iemax(1)=1+ismax(1)
         do j=2,nz
           iemax(j)=1+max0(ismax(j-1),ismax(j))
         enddo
         iemax(nz1)=1+ismax(nz)
c        set er and ez to zero
       do i=1,im2
         er(i)=0.0
         ez(i)=0.0
       enddo
c          sum up fields
       do 220 js=1,nz
       js1=js+1
       ism=ismax(js)
       if (ism.eq.0) go to 220
       do 210 is=1,ism
       l=(js-1)*nr+is
       a1=aa(l)
       if (a1.eq.0.) go to 210
       l=(is-1)*im3
       do 180 je=1,js
       k1=l+(js-je)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 180
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)-a1*ezs(k)
       enddo
  180  continue
       do 200 je=js1,nz1
       k1=l+(je-js1)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 200
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)+a1*ezs(k)
       enddo
  200 continue
  210 continue
  220 continue
c          evaluate and apply impulse
       rrmax=0.
       zzmax=0.
       zzmin=1000.
       npz=0
       npr=0
       do np=1,ngood
c
c  Transforming to the bunch reference frame
c
         dwc=f(7,np)-xmat
         gm1=dwc/xmat
c    convert xp an yp from mrad to rad
         f3np=f(3,np)*1.e-03
         f5np=f(5,np)*1.e-03
comment   gm1*(2.+gm1)=(gam-1)*(gam+1)=gam*gam-1=beta*beta*gam*gam
         bgz=sqrt(gm1*(2.+gm1))
         bgx=bgz*f3np
         bgy=bgz*f5np
         gamma=1.+gm1
c  Particle momentum in the bunch frame
c
         bgzstar=gam*(bgz-beta*gamma)
c
c  Particle energy in bunch frame
c
         gstar=gam*(gamma-beta*bgz)
c
         r=sqrt((f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq)
         if(r.ge.rrmax) rrmax=r
         if (r.eq.0.) r=.000001
         xor=(f(2,np)-xbar)*xfac/r
         yor=(f(4,np)-ybar)*yfac/r
         if (r.gt.rmax) then
           npr=npr+1
           go to 230
         endif
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if(z.ge.zzmax) zzmax=z
         if(z.lt.zzmin) zzmin=z
         if (abs(z).gt.hl) then
           npz=npz+1
           go to 230
         endif
c          interpolate impulse within mesh.
         rb=r/dr
         i=idint(1.+rb)
         a=rb-float(i-1)
         b=1.-a
         zb=(z+hl)/dz
         j=idint(1.+zb)
         c=zb-float(j-1)
         d=1.-c
         l=i+(j-1)*nr1
         m=l+nr1
         cbgr=c3*(d*(a*er(l+1)+b*er(l))+c*(a*er(m+1)+b*er(m)))
         cbgzs=c3*(d*(a*ez(l+1)+b*ez(l))+c*(a*ez(m+1)+b*ez(m)))
c *******************
ccc         write(14,5755)np,cbgr,cbgzs,c3,a,b,c,d
ccc         write(14,5755)np,l,m,er(l+1),er(l),er(m+1),er(m)
ccc5755   format(2x,i5,7(2x,e12.5))
ccc5755   format(3(2x,i5),4(2x,e12.5))
c *******************
c     different space charge in the bunch (valero)
         cbgr=cbgr*abs(f(9,np))
         cbgzs=cbgzs*abs(f(9,np))
         go to 260
c        estimate impulse based on point charge at xbar,ybar,pbar.
c          estimate impulse based on point charge at xbar,ybar,pbar.
230      continue
         d=sqrt(z**2+r**2)
         rod3=r/d**3
         zod3=z/d**3
         if (nip.eq.0) go to 250
c          include neighboring bunches.
         do i=1,nip
           xi=i
           do j=1,2
             s=z+xi*pl
             d=sqrt(s**2+r**2)
             rod3=rod3+r/d**3
             zod3=zod3+s/d**3
             xi=-xi
           enddo
         enddo
c  Evaluate impulse.
c
250      cbgr=eng*c1*c3*rod3*pi/2.
         cbgzs=eng*c1*c3*zod3*pi/2.
c     different charges in the bunch (valero)
         cbgr=cbgr*abs(f(9,np))
         cbgzs=cbgzs*abs(f(9,np))
c
c  Apply impulse and transform back to lab frame.
c
260      bgx=bgx+cbgr*xor
         bgy=bgy+cbgr*yor
         pbgzstar=bgzstar
         bgzstar=bgzstar+cbgzs
         gstar=1.+0.5*bgzstar**2
         bgzf=gam*(bgzstar+beta*gstar)
         f3=bgx/bgzf
         f5=bgy/bgzf
         dww=f(7,np)-xmat
         dws=dww*((gamma+1.)/gamma)*(bgzf-bgz)/bgz
ccc         write(14,5755)np,pbgzstar,cbgzs,bgzstar
ccc5755   format(2x,i5,3(2x,e12.5))
c  ********************
         if(.not.iesp) then
c     load the entrance beam parameters for cavities or gaps
           do js=1,7
             f(js,np)=fs(js,np)
           enddo
c     dxp and dyp are the jumps of xp and yp (in rad) at the position dz1*xpsc (in m)
           dxp=f3-f3np
           dyp=f5-f5np
c   correction of xp and yp ( in rad)
           f(3,np)=f(3,np)+dxp*1000.
           f(5,np)=f(5,np)+dyp*1000.
           f(2,np)=f(2,np)-dz1*100.*dxp*xpsc
           f(4,np)=f(4,np)-dz1*100.*dyp*xpsc
           dwp(np)=dws
         else
           f(3,np)=f3*1000.
           f(5,np)=f5*1000.
           f(7,np)=f(7,np)+dws
         endif
       enddo
       return
       end
       SUBROUTINE scheff_sep
c  ..............................................................................
c    SCHEFF_sep special space charge method
c     This version of SCHEFF, starting from Swesson version, has modifications made to include
c     corrections for relativistic beams. The dynamics
c     have been modified to transform to the beam rest frame,
c     calculate the space-charge kicks in this frame, and then
c     transform back to the lab frame.
c     input data
c         sce(1)=beam current in ma.
c         sce(2)=radial extension in rms multiples
c         sce(3)=longitudinal extension in rms multiples
c         sce(4)=no. of radial mesh intervals (le 20)
c         sce(5)=no. of longitudinal mesh intervals (le 40)
c         sce(6)=no. of adjacent bunches, applicable for buncher studies
c                and should be 0 for linac dynamics
c         sce(7)=distance between adjacent beam pulses in cm (transport studies)
c                input zero to get (beta*lambda) default
c         sce(8)=desactived
c         sce(9)=option to integrate space charge forces over box
c                  if.eq.0. no integration  see sub gaus for further
c                  explanation.
c         sce(10) =1 : call in quads,solenoids,accelarating elements
c         sce(10) =2 : call in drifts,accelarating elements
c         sce(10) =3 : call at both
c  standard SCHEFF parameters (see sub schfdyn and user guide)
c         sce(2)=4
c         sce(3)=4
c         sce(4)=20
c         sce(5)=40
c         sce(6)=0
c         sce(7)=0
c         sce(9)=0
c  ..............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/DYN/TREF,VREF
       COMMON/CMPTE/IELL
       COMMON/part/xc(iptsz),yc(iptsz),zc(iptsz)
       COMMON/DIMENS/zcp(iptsz),xcp(iptsz),ycp(iptsz)
       COMMON/HERMT/AFXT(22),AFYT(22),AFZT(22)
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/SC3/BEAMC,SCDIST,SCE10,CPLM,ECT,APL,ICHAES,ISCSP
       common/faisc/f(10,iptsz),imax,ngood
       COMMON/CDEK/DWP(iptsz)
       COMMON/BEAMSA/FS(7,iptsz)
       COMMON/DCSPA/IESP
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       common/posc/xpsc
       LOGICAL ichaes,iesp
       common/bg/bsc,gsc,phis,wsync
       common/stsc/beami,wavel,freq,btazero,frrms,fzrms,nr,nz
       common/stsc1/beams,im1,im2,im3,nr1,nz1,nq
       common/fldcom/ rp, zp,pl,opt,nip
       common/spacech1/ rm(21), zm(41), rs(20), ers(16800), ezs(16800),
     1 ez(861), aa(800), rssq(20), zzs(41)
     2 ,er(861),rss(20), ismax(40), iemax(41)
       common/rcshef/sce(20)
       common/mcs/imcs,ncstat,cstat(20)
       common/conti/irfqp
       logical irfqp
c     set up field tables with int=0
       IF (beams.eq.0. .OR. scdist.eq.0.) return
c  isepa true ---> the ncstat bunches are separated
        IELL=IELL+1
        write(16,*) ' ****SCHEFF  ',iell
        write(16,*) ' states charges ',ncstat
        do isp=1,ncstat
          ngisp=0
          gmoy=0.
          do np=1,ngood
            if (f(9,np).eq.cstat(isp)) then
              ngisp=ngisp+1
              gmoy=f(7,np)/xmat+gmoy
            endif
          enddo
          gmoy=gmoy/float(ngisp)
cold          beamc=beams/ncstat
          beamc=beams*ngisp/ngood
          beami=beamc/gmoy
          call pintim1(isp)
          CALL sizer(isp,xrms,yrms,zrms)
       write(16,6875) cstat(isp),beamc,xrms,yrms,zrms
c   write the size rms in the output file 17
cold       write(17,25) iell,xrms,yrms,zrms
cold25     format(2x,i5,3(2x,e12.5))
6875   format(' charge: ',f8.0,' bunch intensity: ',e12.5,' amp',/,
     *        ' with RMS size(m)',e12.5,2x,e12.5,2x,e12.5)
         rrms=sqrt(xrms*xrms+yrms*yrms)
c    change unit: m==>cm
         zrms1=zrms*100.
         rrms=rrms*100.
         dr=rrms*frrms/float(nr)
         dz=zrms1*fzrms/float(nz)
         rmax=float(nr)*dr
c          load rm, zm, rs, zs
         rm(1)=0.0
         do i=2,nr1
           rm(i)=float(i-1)*dr
           rssq(i-1)=.5*(rm(i-1)**2+rm(i)**2)
           rss(i-1)=0.5*(rm(i-1)+rm(i))
           rs(i-1)=sqrt(rssq(i-1))
         enddo
         zs=.5*dz
         do i=1,nz1
           zm(i)=float(i-1)*dz
           zzs(i)=zm(i)+zs
         enddo
         hl=float(nz)*zs
c          load ers and ezs
c     mesh dimensions are in cm. ers and ezs are in 1/cm.
c     c1, c2 and c3 are in cm., and c4 is in mev-cm.
c     q=coulombs/point.   (2/pi)*e/(4*pi*epsilon)=572167 cm mev/coul.
         q=beami/(freq*float(nq))
         c1=572167.*q/xmat
         l=0
         do 40 k=1,nr
           rfac=(rm(k+1)**2-rm(k)**2)*dz/2.
           if(opt.eq.0.) rfac=1.
           do 40 j=1,nz
             zp=zm(j+1)
             do 40 i=1,nr1
               rp=rm(i)
               if(opt.eq.0.) call flds(rs(k),zs,er1,ez1)
               if(opt.eq.0.) go to 35
               call gaus(rm(k),rm(k+1),zm(1),zm(2),opt,er1,ez1)
 35            l=l+1
               ers(l)=c1*er1/rfac
   40    ezs(l)=c1*ez1/rfac
         dz1=scdist/100.
         dist=scdist
c sup       WRITE(16, *) ' fields acting length(cm): ',DIST
c          evaluate and apply space charge effects.
c     phimc=phi of mesh center.
C   Shifts particle coordinates to a single point in time. Uses
C   a linear shift
C      Beam c.g.
       xbar=0.
       ybar=0.
       zbar=0.
       brmoy=0.
       trmoy=0.
       do np=1,ngood
         if(f(9,np).eq.cstat(isp)) then
           gpai=f(7,np)/xmat
           brmoy=brmoy+sqrt(1.-1./(gpai*gpai))
           trmoy=trmoy+f(6,np)
         endif
       enddo
       trmoy=trmoy/float(ngisp)
       phimc=trmoy*fh
       pbar=phimc
       beta=brmoy/float(ngisp)
       gsc=1./sqrt(1.-beta*beta)
       bg=beta*gsc
       c3=dist/bg
       c4=dist*xmat
c          evaluate ng, xbar, ybar, and pbar.
       c2=beta*wavel/(2.*pi)
       gmsq=1.+bg**2
       gam=sqrt(gmsq)
c          evaluate ng, xbar, ybar
       ng=0
       xbar=0.
       ybar=0.
       xsq=0.
       ysq=0.
       do np=1,ngood
         if(f(9,np).eq.cstat(isp)) then
           ng=ng+1
           x=f(2,np)
           y=f(4,np)
           xf=f(3,np)
           yf=f(5,np)
           xbar=xbar+x
           ybar=ybar+y
           xsq=xsq+x**2
           ysq=ysq+y**2
         endif
       enddo
       eng=float(ngisp)
       xbar=xbar/eng
       ybar=ybar/eng
c  the mesh center is phi*syn
       xsq=xsq/eng
       ysq=ysq/eng
       epsq=sqrt((xsq-xbar*xbar)/(ysq-ybar*ybar))
       epsqi=1./epsq
       xfac=2./(epsq + 1.)
       yfac=epsq*xfac
c          clear and load bins
       ng=0
       do i=1,im1
         aa(i)=0.0
       enddo
       do 120 np=1,ngood
         if(f(9,np).eq.cstat(isp)) then
         rsq=(f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq
c     i=sqrt(rsq)/dr+1.
         r=sqrt(rsq)
         halfdr=dr*0.5
         i=idint(r/dr+1.0)
         if (i.gt.nr) go to 120
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if (abs(z).ge.hl) go to 120
c------distribute charge among adjacent bins.
         ng=ng+1
         zz=z+hl
         jm1=idint(zz/dz+1.)
         i1=i+1
c     if (rsq.lt.rssq(i)) i1=i-1
         if (rsq.lt.rss (i)) i1=i-1
         if (i1.lt.1) i1=1
         if (i1.gt.nr) i1=nr
         j1=jm1+1
         if (zz.lt.zzs(jm1)) j1=jm1-1
         if (j1.lt.1) j1=1
         if (j1.gt.nz) j1=nz
         a=1.
c     if (i1.ne.i) a=(rsq-rssq(i1))/(rssq(i)-rssq(i1))
       if (i1.eq.i) then
        a=1.
       else
        rdr2=rsq/dr**2
        sqr=sqrt(4.*rdr2-1.)
        rminsq=(halfdr*(sqr-1.))**2
        rmaxsq=(halfdr*(sqr+1.))**2
        if (i1.lt.i) then
          a=(rmaxsq-rm(i)**2)/(rmaxsq-rminsq)
        else
          a=(rm(i1)**2-rminsq)/(rmaxsq-rminsq)
        endif
       endif
       b=1.-a
       cc=1.
       if (j1.ne.jm1) cc=(zz-zzs(j1))/(zzs(jm1)-zzs(j1))
         d=1.-cc
         k=(jm1-1)*nr+i
         aa(k)=aa(k)+a*cc
         k=k+i1-i
         aa(k)=aa(k)+b*cc
         k=(j1-1)*nr+i
         aa(k)=aa(k)+a*d
         k=k+i1-i
         aa(k)=aa(k)+b*d
       endif
  120  continue
       eng=float(ng)
c          find ismax for each j
      do j=1,nz
       l=(j-1)*nr
       k=nr
       do i=1,nr
        m=l+k
        if(aa(m).le.0.00) then
         k=k-1
         go to 130
        else
         go to 140
        endif
130     continue
       enddo
140       ismax(j)=k
      enddo
c        find iemax for each j
         iemax(1)=1+ismax(1)
         do j=2,nz
           iemax(j)=1+max0(ismax(j-1),ismax(j))
         enddo
         iemax(nz1)=1+ismax(nz)
c        set er and ez to zero
       do i=1,im2
         er(i)=0.0
         ez(i)=0.0
       enddo
c          sum up fields
       do 220 js=1,nz
       js1=js+1
       ism=ismax(js)
       if (ism.eq.0) go to 220
       do 210 is=1,ism
       l=(js-1)*nr+is
       a1=aa(l)
       if (a1.eq.0.) go to 210
       l=(is-1)*im3
       do 180 je=1,js
       k1=l+(js-je)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 180
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)-a1*ezs(k)
       enddo
  180  continue
       do 200 je=js1,nz1
       k1=l+(je-js1)*nr1
       n1=(je-1)*nr1
       iem=iemax(je)
       if (iem.le.1) go to 200
       do ie=1,iem
         n=n1+ie
         k=k1+ie
         er(n)=er(n)+a1*ers(k)
         ez(n)=ez(n)+a1*ezs(k)
       enddo
  200 continue
  210 continue
  220 continue
c          evaluate and apply impulse
       rrmax=0.
       zzmax=0.
       zzmin=1000.
       npz=0
       npr=0
       do np=1,ngood
c  Transforming to the bunch reference frame
c
         if(f(9,np).eq.cstat(isp)) then
         dwc=f(7,np)-xmat
         gm1=dwc/xmat
c    convert xp an yp from mrad to rad
         f3np=f(3,np)*1.e-03
         f5np=f(5,np)*1.e-03
comment   gm1*(2.+gm1)=(gam-1)*(gam+1)=gam*gam-1=beta*beta*gam*gam
         bgz=sqrt(gm1*(2.+gm1))
         bgx=bgz*f3np
         bgy=bgz*f5np
         gamma=1.+gm1
c  Particle momentum in the bunch frame
c
         bgzstar=gam*(bgz-beta*gamma)
c
c  Particle energy in bunch frame
c
         gstar=gam*(gamma-beta*bgz)
c
         r=sqrt((f(2,np)-xbar)**2*epsqi+(f(4,np)-ybar)**2*epsq)
         if(r.ge.rrmax) rrmax=r
         if (r.eq.0.) r=.000001
         xor=(f(2,np)-xbar)*xfac/r
         yor=(f(4,np)-ybar)*yfac/r
         if (r.gt.rmax) then
           npr=npr+1
           go to 230
         endif
         zph=f(6,np)*fh
         z=-c2*(zph-phimc)
         if(z.ge.zzmax) zzmax=z
         if(z.lt.zzmin) zzmin=z
         if (abs(z).gt.hl) then
           npz=npz+1
           go to 230
         endif
c          interpolate impulse within mesh.
         rb=r/dr
         i=idint(1.+rb)
         a=rb-float(i-1)
         b=1.-a
         zb=(z+hl)/dz
         j=idint(1.+zb)
         c=zb-float(j-1)
         d=1.-c
         l=i+(j-1)*nr1
         m=l+nr1
         cbgr=c3*(d*(a*er(l+1)+b*er(l))+c*(a*er(m+1)+b*er(m)))
         cbgzs=c3*(d*(a*ez(l+1)+b*ez(l))+c*(a*ez(m+1)+b*ez(m)))
         cbgr=cbgr*abs(f(9,np))
         cbgzs=cbgzs*abs(f(9,np))
         go to 260
c        estimate impulse based on point charge at xbar,ybar,pbar.
c          estimate impulse based on point charge at xbar,ybar,pbar.
230      continue
         d=sqrt(z**2+r**2)
         rod3=r/d**3
         zod3=z/d**3
         if (nip.eq.0) go to 250
c          include neighboring bunches.
         do i=1,nip
           xi=i
           do j=1,2
             s=z+xi*pl
             d=sqrt(s**2+r**2)
             rod3=rod3+r/d**3
             zod3=zod3+s/d**3
             xi=-xi
           enddo
         enddo
c  Evaluate impulse.
c
250      cbgr=eng*c1*c3*rod3*pi/2.
         cbgzs=eng*c1*c3*zod3*pi/2.
         cbgr=cbgr*abs(f(9,np))
         cbgzs=cbgzs*abs(f(9,np))
c
c  Apply impulse and transform back to lab frame.
c
260      bgx=bgx+cbgr*xor
         bgy=bgy+cbgr*yor
         bgzstar=bgzstar+cbgzs
         gstar=1.+0.5*bgzstar**2
         bgzf=gam*(bgzstar+beta*gstar)
         f3=bgx/bgzf
         f5=bgy/bgzf
         dww=f(7,np)-xmat
         dws=dww*((gamma+1.)/gamma)*(bgzf-bgz)/bgz
c  ********************
         if(.not.iesp) then
c     load the entrance beam parameters for cavities or gaps
           do js=1,7
             f(js,np)=fs(js,np)
           enddo
c     dxp and dyp are the jumps of xp and yp (in rad) at the position dz1*xpsc (in m)
           dxp=f3-f3np
           dyp=f5-f5np
c   correction of xp and yp ( in rad)
           f(3,np)=f(3,np)+dxp*1000.
           f(5,np)=f(5,np)+dyp*1000.
           f(2,np)=f(2,np)-dz1*100.*dxp*xpsc
           f(4,np)=f(4,np)-dz1*100.*dyp*xpsc
           dwp(np)=dws
         else
           f(3,np)=f3*1000.
           f(5,np)=f5*1000.
           f(7,np)=f(7,np)+dws
         endif
        endif
       enddo
c big enddo of isp
      enddo
      return
      end
       SUBROUTINE rfkick(v,dp,harm,nvf)
c   ............................................................................
c          RFKICK  (NO SPACE CHARGE EFFECT)
c
c	Contributing Author: Daniel Alt, NSCL/MSU, East Lansing,  MI, USA
c   Date: 23-May-2014
c
c	Electric RF kicker.  Simulates a sine wave chopper consisting
c	of two plates deflecting the beam in the transverse direction.
c	This is a zero length element.
c
c	V: Voltage Factor.  Consists of plate voltage in kV, times
c		the electrode length divided by the gap between them.
c	DP: RF phase offset in radians
c	Harm: Kicker harmonic number relative to current reference frequency
c	NVF: Kicker axis: 0 = horziontal, 1= vertical (use negative
c		voltage for negative deflection)
c   ............................................................................
       implicit real*8 (a-h,o-z)
       parameter (iptsz=100002,maxcell=3000,maxcell1=3000)
       COMMON/RIGID/BORO
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       COMMON/DYN/TREF,VREF
       COMMON/DAVCOM/dav1(maxcell1,40),davtot,iitem(maxcell1),idav
       common/etcha1/dav2(maxcell1,33),ichas(iptsz),chasit
       COMMON/TTFC/TK,T1K,T2K,SK,S1K,S2K,FH
       common/faisc/f(10,iptsz),imax,ngood
       common/etcom/cog(8),exten(17),fd(iptsz)
       COMMON/FENE/WDISP,WPHAS,WX,WY,RLIM,IFW
       common/corec/tref1
       COMMON/QMOYEN/QMOY
       common/aerp/vphase,vfield,ierpf
       COMMON/ITVOLE/ITVOL,IMAMIN
       COMMON/COMPT/NRRES,NRTRE,NRBUNC,NRDBUN
       COMMON/SHIF/DTIPH,SHIFT
       common/tofev/ttvols
       character*1 cr
c       dimension vecx(1)
       LOGICAL chasit,itvol,imamin,shift
C      ENVELOPE
       call stapl(davtot*10.)
cxx       ilost=0
       twopi=2.*pi
       freq=fh/twopi
       wavel=vl/freq
       fcpi=fh*180./pi
c print out on terminal of transport element # on one and the same line
       cr=char(13)
#if F77
       WRITE(6,8254) nrtre,cr
8254   format('Transport element:',i5,a1,$)
#else
       WRITE(6,8254,ADVANCE='NO') nrtre,cr
8254   format('Transport element:',i5,a1)
#endif
       if (harm.le.0.) harm=1.
c   test window - unneccesary at START of element - Alt
c       call reject(ilost)
C   Calculate tof, beta, gamma, k.e. for c.o.g.
       tcog=0.
       bcog=0.
       do np=1,ngood
         tcog=tcog+f(6,np)
         gpa=f(7,np)/xmat
         bcog=sqrt(1.-1./(gpa*gpa))+bcog
       enddo
       tcog=tcog/float(ngood)
       bcog=bcog/float(ngood)
       gcog=1./sqrt(1.-bcog*bcog)
       encog=xmat*gcog-xmat
c adjustement of the phase of RF w.r.t. the T.O.F.
       xkpi=0.
       if(imamin) then
         ttvpi=harm*ttvols*fcpi
         xkpi=ttvpi/360.
         ixkpi=int(xkpi)
         xkpi=(xkpi-float(ixkpi))*360.
         write(16,*) ' *** TOF correction:',-xkpi,' deg'
         dp=dp-xkpi*pi/180.
         write(16,*)' ***phase of RF adjusted : ',dp*180./pi,' deg'
       endif
c  start of write to file '.short' for kicker
       idav=idav+1
       iitem(idav)=22
       dav1(idav,1)=v
       dav1(idav,2)=dp*180./pi
       dav1(idav,3)=nvf
       dav1(idav,4)=davtot*10.
C       if(itvol) dav1(idav,5)=-xkpi
       dav1(idav,5)=harm
c  end
       WRITE(16,178)
178    FORMAT(/,' Longitudinal parameters',/,
     2 5X,'   BETA     GAMMA      ENERGY(MeV) ',
     3 '       TOF(deg)     TOF(sec)')
       WRITE(16,1788) bcog,gcog,encog,tcog*fcpi,tcog
1788   FORMAT(' COG ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       e0t=harm*v/(bcog*wavel)
       cay=harm*twopi/(bcog*gcog*wavel)
       caysq=cay**2
       con=twopi*e0t*qmoy/xmat
       rad=pi/180.
c shift=true => reference and COG seperated, otherwise reference=COG
c --- save the reference
       ovref=vref
       otref=tref
c --- shift = false: the reference particle is the cog
       if(shift) then
         ovref=vref
         beref=vref/vl
         gamref=1./sqrt(1.-beref*beref)
         older=xmat*gamref
       else
         tref=tcog
         vref=bcog*vl
         ovref=vref
         beref=bcog
         gamref=1./sqrt(1.-beref*beref)
         older=xmat*gamref
       endif
c ---  if imamin = false: phase setting has been forced equal to dp, otherwise phase setting has been adjusted
       ENRPRIN=older-xmat
       WRITE(16,165) beref,gamref,ENRPRIN,tref*fh*180./pi,tref
165    FORMAT(' REF ',F7.5,3X,E12.5,1X,E12.5,3X,E12.5,3x,E12.5)
       wsync=0.
       bcour=0.
       tcog=0.
       do np=1,ngood
         a=harm*(f(6,np)-tref+ttvols)*fh+dp
         s=sin(a)
         w=f(7,np)-xmat
         gpai=f(7,np)/xmat
         bg=sqrt(w/xmat*(2.+w/xmat))
         bpai=bg/gpai
         const=(gpai/(gpai*gpai-1.))*f(9,np)
         disp=const*v/xmat*s
         bcour=bcour+bpai
         tcog=tcog+f(6,np)
         if (nvf.eq.0) then
           f(3,np)=f(3,np)+disp
         else if (nvf.eq.1) then
           f(5,np)=f(5,np)+disp
         else
           write(6,*) 'Invalid parameter NVF in RFKICK'
         endif
         wsync=wsync+w
       enddo
       wsync=wsync/float(ngood)
       bcour=bcour/float(ngood)
       tcog=tcog/float(ngood)
c Test window
       call reject(ilost)
c dave start for kicker
       dav1(idav,36)=ngood
       call emiprt(0)
       return
       end
