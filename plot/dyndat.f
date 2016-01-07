        PROGRAM dyndat
        implicit real*8 (a-h,o-z)
        CHARACTER*11 pfnm(20),strng
c        CHARACTER*11 pfnm(20),prnm(20),strng
        character*5 cccst(20)
        common/prtcnt/imax
        common/wfil2/title,labels
        common/fichier/filenm,pfnm
        common/opsys/opsy
        common/chstat1/cst(250002),cstat(20),fcstat(20),ncstat,mcstat
        common/chstat2/cccst
        common/pscl/yminsk,iskale
        DIMENSION dx(250002),dxp(250002),data(10)
        CHARACTER*80 title,command,junk
        CHARACTER*40 labels(20)
        character filenm*23
        character*1 fsave
        character*6 fmt
        character*132 ras
        CHARACTER backslash
        parameter (backslash="\\")
        dimension icstat(20),zstat(20)
        dimension cx(250002),cxp(250002),ctrx1(250002),ctry1(250002)
        dimension cy(250002),cyp(250002),ctrx2(250002),ctry2(250002)
        dimension cz(250002),czp(250002),ctrx3(250002),ctry3(250002)
        common/wfil20/xmin(10),xmax(10),ymin(10),ymax(10)
        common/mingw/mg
        logical mg
c if mg=.true., use MINGW on windows, which has a different result for ctime function than standard gfortran
c default is mg=.false.        
        mg=.false.
        if(iargc() .ge. 1) then
          call getArg(1, command)
        else
          write(6,'(a)') "For LINUX   GNUPLOT enter L or l"
          write(6,'(a)') "For MAC     GNUPLOT enter M or m"
          write(6,'(a)') "For WINDOWS GNUPLOT enter W or w"
          write(6,'(a$)') "Enter GNUPLOT version "
          read(5,*) command
        endif
        if(command(1:1).eq.'L' .or. command(1:1).eq.'l') opsy=1
        if(command(1:1).eq.'W' .or. command(1:1).eq.'w') then
          opsy=2
          if(command(2:3).eq.'MG' .or. command(2:3).eq.'mg') mg=.true.
        endif  
        if(command(1:1).eq.'M' .or. command(1:1).eq.'m') opsy=3      
c opsy=1 --> LINUX   GNUPLOT version
c opsy=2 --> WINDOWS GNUPLOT version
c opsy=3 --> MAC     GNUPLOT version
        if(opsy.eq.1) then
          write(6,'(A)') 'PLOTIT V2.6 13-Nov-2013 LINUX'
        elseif(opsy.eq.2) then
          write(6,'(A)') 'PLOTIT V2.6 13-Nov-2013 WINDOWS'
        elseif(opsy.eq.3) then
          write(6,'(A)') 'PLOTIT V2.6 13-Nov-2013 MAC'
        else
          write(6,*) 'Error in operating system type entry'
          stop
        endif
c V1.2 Original released version, compatible with g77
c V2.0 Formatting changed to be compatible with gfortran and with WGNUPLOT gp440win32
c V2.1 Replaced 'dots' with 'points' syntax
c V2.2 Fixed compiler warnings related to array sizing of ctrx1,2,3 and ctry1,2,3
c V2.3 Mods to add MAC as a valid operating system for gnuplot
c V2.4 Allow for non-integer charge states
c V2.5 Allow for longer file names under the save option
c      Made dyndat compatible with MINGW gfortran in view of different result for ctime
c      function than standard gfortran
c      Fixed bug related to MAC: on a MAC 02 Apr is shown as 2 Apr (on windows as
c      02 Apr); fill the blank with the character '0' to avoid error message on file 
c      copy
c      Fix color of points used in particle plots for MAC (changed from grey to black)
c V2.6 Change from 100k to 250k macro particles
         write(6,*)
C igrtyp is type of graph (there is no igrtyp=8,9,10,13,14,15 or 16)
C single charge state, no zones:
C        igrtyp=1  for xx'-yy'-xy-zz' plots
C        igrtyp=2  for z-x, z-y plots & profiles
C        igrtyp=3  for x-z, y-z envelopes
C        igrtyp=4  for dW envelope
C        igrtyp=5  for dPHI envelope    
C multi charge state, no zones:
C        igrtyp=6  for xx'-yy'-xy-zz' plots for multi-charge state beam
C        igrtyp=7  for z-x, z-y plots & profiles for multi-charge state beam
C with zones, single charge state:
C        igrtyp=11 for xx'-yy'-xy-zz' plots with ZONES card
C        igrtyp=12 for z-x, z-y plots & profiles with ZONES card
C
C        igrtyp=17 or 22 or 27 log scale in bunch profiles
        ncstat=1
        data(1)=0.
        pfnm(1)='dynac01.plt'
        pfnm(2)='dynac02.plt'
        pfnm(3)='dynac03.plt'
        pfnm(4)='dynac04.plt'
        pfnm(5)='dynac05.plt'
        pfnm(6)='dynac06.plt'
        pfnm(7)='dynac07.plt'
        pfnm(8)='dynac08.plt'
        pfnm(9)='dynac09.plt'
        pfnm(10)='dynac10.plt'
        pfnm(11)='dynac11.plt'
        pfnm(12)='dynac12.plt'
        pfnm(13)='dynac13.plt'
        pfnm(14)='dynac14.plt'
        pfnm(15)='dynac15.plt'
        pfnm(16)='dynac16.plt'
        pfnm(17)='dynac17.plt'
        pfnm(18)='dynac18.plt'
        pfnm(19)='dynac19.plt'
        pfnm(20)='dynac20.plt'
        if(opsy.eq.1 .or. opsy.eq.3) then
c LINUX or MAC
          command(1:40)="test ! -e savedplots && mkdir savedplots"
          CALL System(COMMAND(1:40))
        else
C WINDOWS
          command(1:29)="if not exist savedplots\\*.* "
          command(30:58)="echo creating plots directory"
          CALL System(COMMAND(1:58))
          command(30:45)="mkdir savedplots"
          CALL System(COMMAND(1:45))
          command(30:57)="copy ..\\bin\\tst savedplots"
          CALL System(COMMAND(1:57))
        endif
        command=""
        IF (data(1).eq.0.) THEN
           OPEN(unit=66,file='emit.plot')
           data(1)=1.
           data(2)=0.
        ENDIF
10      READ(66,*,end=20) igrtyp
        iskale=0
        IF(igrtyp.eq.17 .or. igrtyp.eq.22 .or. igrtyp.eq.27)THEN
c log scale in bunch profiles
          igrtyp=igrtyp-15
          iskale=1
          read(66,*) yminsk
        ENDIF
        IF (igrtyp.eq.1 .or. igrtyp.eq.6 .or. igrtyp.eq.11) THEN
c x-x', y-y', x-y, z-z' plots
          data(2)=data(2)+1.
          data(5)=igrtyp
          nplot=int(data(2))
          if (igrtyp.eq.6) then
            read(66,*) ncstat
            read(66,*) (cstat(j),j=1,ncstat)
            do j=1,ncstat
              icstat(j)=int(cstat(j))
              write(command(2:5),'(I2)') icstat(j)
              command(1:1)=' '
              cccst(j)=command(1:5)
            enddo
          endif
          if (igrtyp.eq.11) then
            read(66,*) ncstat
            read(66,*) (zstat(j),j=1,ncstat)
            do j=1,ncstat
              cstat(j)=zstat(j)
              write(command(1:3),'(F3.1)') zstat(j)
c              command(1:1)=' '
              cccst(j)=command(1:3)
            enddo
          endif
          READ(66,1000) title(1:80)
1000      format(a80)
          id=nplot
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(1)=uxmin
          xmax(1)=uxmax
          ymin(1)=uymin
          ymax(1)=uymax
          DO i=1,201
            READ(66,*) dx(i),dxp(i)
            ctrx1(i)=dx(i)
            ctry1(i)=dxp(i)
          ENDDO
          READ(66,*) imax
          if (igrtyp.eq.1) then
            DO i=1,imax
              READ(66,*) dx(i),dxp(i)
              cx(i)=dx(i)
              cxp(i)=dxp(i)
            ENDDO
          else
            DO i=1,imax
              READ(66,*) dx(i),dxp(i),cst(i)
              cx(i)=dx(i)
              cxp(i)=dxp(i)
            ENDDO
          endif
          labels(1)='x (cm)'
          labels(2)='xp (mrad)'
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(2)=uxmin
          xmax(2)=uxmax
          ymin(2)=uymin
          ymax(2)=uymax
          DO i=1,201
            READ(66,*) dx(i),dxp(i)
            ctrx2(i)=dx(i)
            ctry2(i)=dxp(i)
          ENDDO
          READ(66,*) imax
          if (igrtyp.eq.1) then
            DO i=1,imax
              READ(66,*) dx(i),dxp(i)
              cy(i)=dx(i)
              cyp(i)=dxp(i)
            ENDDO
          else
            DO i=1,imax
              READ(66,*) dx(i),dxp(i),cst(i)
              cy(i)=dx(i)
              cyp(i)=dxp(i)
            ENDDO
          endif
          labels(3)='y (cm)'
          labels(4)='yp (mrad)'
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(3)=uxmin
          xmax(3)=uxmax
          ymin(3)=uymin
          ymax(3)=uymax
          labels(5)='x (cm)'
          labels(6)='y (cm)'
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(4)=uxmin
          xmax(4)=uxmax
          ymin(4)=uymin
          ymax(4)=uymax
          DO i=1,201
            READ(66,*) dx(i),dxp(i)
            ctrx3(i)=dx(i)
            ctry3(i)=dxp(i)
          ENDDO
          READ(66,*) imax
          if (igrtyp.eq.1) then
            DO i=1,imax
              READ(66,*) dx(i),dxp(i)
              cz(i)=dx(i)
              czp(i)=dxp(i)
            ENDDO
          else
            DO i=1,imax
              READ(66,*) dx(i),dxp(i),cst(i)
              cz(i)=dx(i)
              czp(i)=dxp(i)
            ENDDO
          endif
          labels(7)='z (deg)'
          labels(8)='zp (MeV)'
          call wfile10(0,imax,cx,cxp,cy,cyp,cz,czp,igrtyp)
          call wfile10(1,imax,ctrx1,ctry1,ctrx2,ctry2,ctrx3,ctry3,
     *                 igrtyp)
          isave=0
          call wfile20(isave,igrtyp)
          if(opsy.eq.1) then
C LINUX
            command(1:27)="gnuplot -noraise -geometry "
	      command(28:51)="500x515-250+25 dynac.gnu"
            CALL System(COMMAND(1:51))
          elseif(opsy.eq.3) then
C MAC
            command(1:27)="gnuplot -noraise dynac.gnu "
            CALL System(COMMAND(1:27))
          else
C WINDOWS          
            command(1:18)="wgnuplot dynac.gnu"
            CALL System(COMMAND(1:18))
          endif
        ENDIF
        IF (igrtyp.eq.2 .or. igrtyp.eq.7 .or. igrtyp.eq.12) THEN
c z-x, z-y plots & profiles
          data(5)=igrtyp
          data(2)=data(2)+1.
          nplot=int(data(2))
          if (igrtyp.eq.7) then
            read(66,*) ncstat
            read(66,*) (cstat(j),j=1,ncstat)
            do j=1,ncstat
              icstat(j)=int(cstat(j))
              write(command(2:5),'(I2)') icstat(j)
              command(1:1)=' '
              cccst(j)=command(1:5)
            enddo
          endif
          if (igrtyp.eq.12) then
            read(66,*) ncstat
            read(66,*) (zstat(j),j=1,ncstat)
            do j=1,ncstat
              cstat(j)=zstat(j)
              write(command(1:3),'(f3.1)') zstat(j)
c              command(1:1)=' '
              cccst(j)=command(1:3)
            enddo
          endif
          READ(66,1000) title(1:80)
          id=nplot
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(1)=uxmin
          xmax(1)=uxmax
          ymin(1)=uymin
          ymax(1)=uymax
          READ(66,*) imax
          if (igrtyp.eq.2) then
            DO i=1,imax
              READ(66,*) dx(i),dxp(i)
              cx(i)=dx(i)
              cxp(i)=dxp(i)
            ENDDO
          else
            DO i=1,imax
              READ(66,*) dx(i),dxp(i),cst(i)
              cx(i)=dx(i)
              cxp(i)=dxp(i)
            ENDDO
          endif
          labels(1)='z (cm)'
          labels(2)='x (cm)'
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(2)=uxmin
          xmax(2)=uxmax
          ymin(2)=uymin
          ymax(2)=uymax
          READ(66,*) imax
          if (igrtyp.eq.2) then
            DO i=1,imax
              READ(66,*) dx(i),dxp(i)
              cy(i)=dx(i)
              cyp(i)=dxp(i)
            ENDDO
          else
            DO i=1,imax
              READ(66,*) dx(i),dxp(i),cst(i)
              cy(i)=dx(i)
              cyp(i)=dxp(i)
            ENDDO
          endif
          labels(3)='z (cm)'
          labels(4)='y (cm)'
          call wfile11(imax,cx,cxp,cy,cyp,igrtyp)
          isave=0
          call wfile21(isave,igrtyp)
c profiles
c          prnm(1)='dynac01.pro'
c          prnm(2)='dynac02.pro'
c          prnm(3)='dynac03.pro'
c          prnm(4)='dynac04.pro'
c          prnm(5)='dynac05.pro'
c          prnm(6)='dynac06.pro'
c          prnm(7)='dynac07.pro'
c          prnm(8)='dynac08.pro'
c          prnm(9)='dynac09.pro'
c          prnm(10)='dynac10.pro'
c          prnm(11)='dynac11.pro'
c          prnm(12)='dynac12.pro'
c          prnm(13)='dynac13.pro'
c          prnm(14)='dynac14.pro'
c          prnm(15)='dynac15.pro'
c          prnm(16)='dynac16.pro'
c          prnm(17)='dynac17.pro'
c          prnm(18)='dynac18.pro'
c          prnm(19)='dynac19.pro'
c          prnm(20)='dynac20.pro'
          OPEN(unit=70,file='dynac01.pro')
          OPEN(unit=71,file='dynac02.pro')
          OPEN(unit=72,file='dynac03.pro')
          OPEN(unit=73,file='dynac04.pro')
          OPEN(unit=74,file='dynac05.pro')
          OPEN(unit=75,file='dynac06.pro')
          rewind(unit=70)
          rewind(unit=71)
          rewind(unit=72)
          rewind(unit=73)
          rewind(unit=74)
          rewind(unit=75)
c          READ(66,'(a)') junk
c          write(6,*) junk
c          write(6,*) junk
c          write(6,*) junk
          READ(66,*) imx
          DO i=1,imx
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(70,*)dx(i),dxp(i)
          ENDDO
          READ(66,*) imy
          DO i=1,imy
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(71,*)dx(i),dxp(i)
          ENDDO
          READ(66,*) imz
          DO i=1,imz
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(72,*)dx(i),dxp(i)
          ENDDO
          READ(66,*) imxp
          DO i=1,imxp
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(73,*)dx(i),dxp(i)
          ENDDO
          READ(66,*) imyp
          DO i=1,imyp
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(74,*)dx(i),dxp(i)
          ENDDO
          READ(66,*) imzp
          DO i=1,imzp
             READ(66,*) dx(i),dxp(i)
             if (dxp(i).eq.0.) dxp(i)=1.e-8
             write(75,*)dx(i),dxp(i)
          ENDDO
          close(70)
          close(71)
          close(72)
          close(73)
          close(74)
          close(75)
          if(opsy.eq.1) then
C LINUX
            command(1:27)="gnuplot -noraise -geometry "
	      command(28:51)="500x515-250+25 dynac.gnu"
            CALL System(COMMAND(1:51))
          elseif(opsy.eq.3) then
C MAC
            command(1:27)="gnuplot -noraise dynac.gnu "
            CALL System(COMMAND(1:27))
          else
C WINDOWS          
            command(1:18)="wgnuplot dynac.gnu"
            CALL System(COMMAND(1:18))
          endif
        ENDIF
        IF (igrtyp.eq.3) THEN
c x,y envelopes as f(z)
          data(5)=igrtyp
          data(2)=data(2)+1.
          nplot=int(data(2))
          READ(66,1000) title(1:80)
          id=nplot
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(1)=uxmin
          xmax(1)=uxmax
          ymin(1)=uymin
          ymax(1)=uymax
          READ(66,*) imax
          DO i=1,imax
            READ(66,*) dx(i),dxp(i)
            cx(i)=dx(i)
            cxp(i)=dxp(i)
          ENDDO
          READ(66,*) imax
          DO i=1,imax
            READ(66,*) dx(i),dxp(i)
            cy(i)=dx(i)
            cyp(i)=dxp(i)
          ENDDO
          labels(1)='z (m)'
          labels(2)='x,y (cm)'
          call wfile3(imax,cx,cxp,cy,cyp)
          isave=0
          call wfile2(isave,2,uxmin,uxmax,uymin,uymax)
          if(opsy.eq.1) then
C LINUX
            command(1:27)="gnuplot -noraise -geometry "
	      command(28:51)="500x515-250+25 dynac.gnu"
            CALL System(COMMAND(1:51))
          elseif(opsy.eq.3) then
C MAC
            command(1:27)="gnuplot -noraise dynac.gnu "
            CALL System(COMMAND(1:27))
          else
C WINDOWS          
            command(1:18)="wgnuplot dynac.gnu"
            CALL System(COMMAND(1:18))
          endif
        ENDIF
        IF (igrtyp.eq.4) THEN
c dW/W envelope as f(z)
          data(5)=igrtyp
          data(2)=data(2)+1.
          nplot=int(data(2))
          READ(66,1000) title(1:80)
          id=nplot
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(1)=uxmin
          xmax(1)=uxmax
          ymin(1)=uymin
          ymax(1)=uymax
          READ(66,*) imax
          DO i=1,imax
            READ(66,*) dx(i),dxp(i)
            cx(i)=dx(i)
            cxp(i)=dxp(i)
          ENDDO
          labels(1)='z (m)'
          labels(2)='dW/W (per mille)'
          call wfile1(imax,cx,cxp,cy,cyp)
          isave=0
          call wfile2(isave,3,uxmin,uxmax,uymin,uymax)
          if(opsy.eq.1) then
C LINUX
            command(1:27)="gnuplot -noraise -geometry "
	      command(28:51)="500x515-250+25 dynac.gnu"
            CALL System(COMMAND(1:51))
          elseif(opsy.eq.3) then
C MAC
            command(1:27)="gnuplot -noraise dynac.gnu "
            CALL System(COMMAND(1:27))
          else
C WINDOWS          
            command(1:18)="wgnuplot dynac.gnu"
            CALL System(COMMAND(1:18))
          endif
        ENDIF
        IF (igrtyp.eq.5) THEN
c dPHI envelope as f(z)
          data(5)=igrtyp
          data(2)=data(2)+1.
          nplot=int(data(2))
          READ(66,1000) title(1:80)
          id=nplot
          READ(66,*) uxmin,uxmax,uymin,uymax
          xmin(1)=uxmin
          xmax(1)=uxmax
          ymin(1)=uymin
          ymax(1)=uymax
          READ(66,*) imax
          DO i=1,imax
            READ(66,*) dx(i),dxp(i)
            cx(i)=dx(i)
            cxp(i)=dxp(i)
          ENDDO
          labels(1)='z (m)'
          labels(2)='dPHI (deg)'
          call wfile1(imax,cx,cxp,cy,cyp)
          isave=0
          call wfile2(isave,4,uxmin,uxmax,uymin,uymax)
          if(opsy.eq.1) then
C LINUX
            command(1:27)="gnuplot -noraise -geometry "
	      command(28:51)="500x515-250+25 dynac.gnu"
            CALL System(COMMAND(1:51))
          elseif(opsy.eq.3) then
C MAC
            command(1:27)="gnuplot -noraise dynac.gnu "
            CALL System(COMMAND(1:27))
          else
C WINDOWS          
            command(1:18)="wgnuplot dynac.gnu"
            CALL System(COMMAND(1:18))
          endif
        ENDIF
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
         write(6,'(a$)')
     *   'Save plot file (y/n/p(rint)/e(xit)/q(uit)/<cr>=n)? '
        else
C WINDOWS          
         write(6,'(a$)')'Save plot file (y/n/e(xit)/q(uit)/<cr>=n)? '
        endif
        read(5,'(A)') fsave
        fsave=fsave(1:1)
        if (fsave.eq.' ' .or. fsave.eq.'n' .or. fsave.eq.'N')then
            write(6,*) ' No files saved'
        elseif (fsave.eq.'e' .or. fsave.eq.'E' .or.
     *          fsave.eq.'q' .or. fsave.eq.'Q') then
          stop  
        elseif (fsave.eq.'p' .or. fsave.eq.'P') then
C currently for linux and MAC only
          if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC 
            OPEN(unit=50,file='dynac.gnu')
            OPEN(unit=51,file='dynacp.gnu')
            write(51,7334)
            write(51,7335)
7333        read(50,'(a)',END=3456) ras
            if (ras(1:5).ne.'pause') then
              k=len_trim(ras)
              call mkfrmt(k,fmt)
              write(51,fmt) ras(1:k)
            endif
7334        FORMAT('set output "dynac.ps"')
7335        FORMAT('set term postscript color')
            goto 7333
3456        continue
            close(50)
            close(51)
            command(1:18)="gnuplot dynacp.gnu"
            CALL System(COMMAND(1:18))
            command(1:13)="lpr dynac.ps"
            CALL System(COMMAND(1:13))
          endif
        elseif (fsave.eq.'y' .or. fsave.eq.'Y') then
c filename format: sXMmmDDYYYYHHMMSSaa.eee (see user guide)        
          isave=1
          IF (igrtyp.eq.1 .or. igrtyp.eq.6 .or. igrtyp.eq.11) THEN
c x-xp', y-xp', x-y, z-zp' plots
            call wfile20(isave,igrtyp)
            if(igrtyp.eq.1) then
              filenm(18:21)='.plt'
              if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                command(1:24)="cp dynac.plt savedplots/"
                command(25:45)=filenm(1:21)
                CALL System(COMMAND(1:45))
                filenm(18:21)='.cnt' 
                command(1:24)="cp dynac.cnt savedplots/"
                command(25:45)=filenm(1:21)         
                CALL System(COMMAND(1:45))
              else
C WINDOWS          
                command(1:25)="copy dynac.plt savedplots"
                command(26:26)=backslash
                command(27:47)=filenm(1:21)
                CALL System(COMMAND(1:47))
                filenm(18:21)='.cnt' 
                command(1:25)="copy dynac.cnt savedplots"
                command(27:47)=filenm(1:21)         
                CALL System(COMMAND(1:47))
              endif
            else
              if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                do j=1,mcstat
                  strng=pfnm(j)
                  filenm(18:23)=strng(6:11)
                  command(1:3)="cp "
                  command(4:14)=pfnm(j)
                  command(15:26)=' savedplots/'
                  command(27:49)=filenm(1:23)
                  CALL System(COMMAND(1:49))
                enddo
                filenm(18:21)='.cnt' 
                command(1:24)="cp dynac.cnt savedplots/"
                command(25:45)=filenm(1:21)         
                CALL System(COMMAND(1:45))
              else
C WINDOWS           
                do j=1,mcstat
                  strng=pfnm(j)
                  filenm(18:23)=strng(6:11)
                  command(1:5)="copy "
                  command(6:16)=pfnm(j)
                  command(17:27)=' savedplots'
                  command(28:28)=backslash
                  command(29:51)=filenm(1:23)
                  CALL System(COMMAND(1:51))
                enddo
                filenm(18:21)='.cnt' 
                command(1:25)="copy dynac.cnt savedplots"
                command(26:26)=backslash
                command(27:47)=filenm(1:21)         
                CALL System(COMMAND(1:47))
              endif
            endif
          ENDIF
          IF (igrtyp.eq.2 .or. igrtyp.eq.7 .or. igrtyp.eq.12) THEN
c z-x, z-y plots
            call wfile21(isave,igrtyp)
            if(igrtyp.eq.2) then
              filenm(18:21)='.plt'
              if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                command(1:24)="cp dynac.plt savedplots/"
                command(25:45)=filenm(1:21)
                CALL System(COMMAND(1:45))
              else
C WINDOWS          
                command(1:25)="copy dynac.plt savedplots"
                command(26:26)=backslash
                command(27:47)=filenm(1:21)
                CALL System(COMMAND(1:47))
              endif
            else
              if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                do j=1,mcstat
                  strng=pfnm(j)
                  filenm(18:23)=strng(6:11)
                  command(1:3)="cp "
                  command(4:14)=pfnm(j)
                  command(15:26)=' savedplots/'
                  command(27:49)=filenm(1:23)
                  CALL System(COMMAND(1:49))
                enddo
              else
C WINDOWS           
                do j=1,mcstat
                  strng=pfnm(j)
                  filenm(18:23)=strng(6:11)
                  command(1:5)="copy "
                  command(6:16)=pfnm(j)
                  command(17:27)=' savedplots'
                  command(28:28)=backslash
                  command(29:51)=filenm(1:23)
                  CALL System(COMMAND(1:51))
                enddo
              endif
            endif
c profiles
            filenm(20:23)='.pro'
            if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
              command(1:26)="cp dynac01.pro savedplots/"
              filenm(18:19)='01'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
              command(1:26)="cp dynac02.pro savedplots/"
              filenm(18:19)='02'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
              command(1:26)="cp dynac03.pro savedplots/"
              filenm(18:19)='03'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
              command(1:26)="cp dynac04.pro savedplots/"
              filenm(18:19)='04'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
              command(1:26)="cp dynac05.pro savedplots/"
              filenm(18:19)='05'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
              command(1:26)="cp dynac06.pro savedplots/"
              filenm(18:19)='06'
              command(27:49)=filenm(1:23)
              CALL System(COMMAND(1:49))
            else
C WINDOWS    
c s2Jan14200420535101.pro
              command(1:27)="copy dynac01.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='01'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
              command(1:27)="copy dynac02.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='02'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
              command(1:27)="copy dynac03.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='03'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
              command(1:27)="copy dynac04.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='04'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
              command(1:27)="copy dynac05.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='05'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
              command(1:27)="copy dynac06.pro savedplots"
              command(28:28)=backslash
              filenm(18:19)='06'
              command(29:51)=filenm(1:23)
              CALL System(COMMAND(1:51))
            endif
          ENDIF
          IF (igrtyp.eq.3) THEN
c x,y envelopes as f(z)
            call wfile2(isave,2,uxmin,uxmax,uymin,uymax)
            filenm(18:21)='.plt'
            if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
              command(1:24)="cp dynac.plt savedplots/"
              command(25:45)=filenm(1:21)
              CALL System(COMMAND(1:45))
            else
C WINDOWS          
              command(1:25)="copy dynac.plt savedplots"
              command(26:26)=backslash
              command(27:47)=filenm(1:21)
              CALL System(COMMAND(1:47))
            endif
          ENDIF
          IF (igrtyp.eq.4) THEN
c dW/W envelope as f(z)
            call wfile2(isave,3,uxmin,uxmax,uymin,uymax)
            filenm(18:21)='.plt'
            if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
              command(1:24)="cp dynac.plt savedplots/"
              command(25:45)=filenm(1:21)
              CALL System(COMMAND(1:45))
            else
C WINDOWS          
              command(1:25)="copy dynac.plt savedplots"
              command(26:26)=backslash
              command(27:47)=filenm(1:21)
              CALL System(COMMAND(1:47))
            endif
          ENDIF
          IF (igrtyp.eq.5) THEN
c dPHI envelope as f(z)
            call wfile2(isave,4,uxmin,uxmax,uymin,uymax)
            filenm(18:21)='.plt'
            if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
              command(1:24)="cp dynac.plt savedplots/"
              command(25:45)=filenm(1:21)
              CALL System(COMMAND(1:45))
            else
C WINDOWS          
              command(1:25)="copy dynac.plt savedplots"
              command(26:26)=backslash
              command(27:47)=filenm(1:21)
              CALL System(COMMAND(1:47))
            endif
          ENDIF
          write(6,*) ' Saved ',filenm(1:17),' in savedplots directory'
          isave=0
        ENDIF
        write(6,7733) int(data(2))
7733    format('Plot ',i3,' has been plotted')
        write(6,*) 
        GOTO 10
20      write(6,*) int(data(2)),' plots in total'
        data(4)=1.
        CLOSE(66)
30      data(3)=float(imax)
40      END
        SUBROUTINE wfile1(imax,x,xp,cx,cy)  
        implicit real*8 (a-h,o-z)
        common/opsys/opsy
c This routine writes the data points which will be plotted by GNU to a file
        dimension x(250002),xp(250002),cx(300),cy(300)
        CHARACTER*80 title,command
        command=""
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          command="rm -f dynac.cnt"
        else
C WINDOWS
          command="if exist dynac.cnt del dynac.cnt"
        endif
        CALL System(COMMAND)
        OPEN(unit=48,file='dynac.cnt')
        DO i=1,201
          write(48,*) cx(i),cy(i)
        ENDDO
        CLOSE(48)
        command=""
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          command="rm -f dynac.plt"
        else
C WINDOWS
          command="if exist dynac.plt del dynac.plt"
        endif
        CALL System(COMMAND)
        OPEN(unit=49,file='dynac.plt')
        DO i=1,imax
          write(49,*) x(i),xp(i)
        ENDDO
        CLOSE(49)
        RETURN
        END
        SUBROUTINE wfile2(isave,icontr,xmin,xmax,ymin,ymax)
        implicit real*8 (a-h,o-z)
c this routine will write the .GNU file containing the commands to be
c executed by GNUPLOT
        common/wfil2/title,labels
        common/fichier/filenm,pfnm
        common/opsys/opsy
        CHARACTER*80 command,title
        CHARACTER*40 labels(20)
        CHARACTER*11 pfnm(20)
        CHARACTER backslash,filenm*23,paf*33
        parameter (backslash="\\")
        command=""
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          command="rm -f dynac.gnu"
        else
C WINDOWS
          command="if exist dynac.gnu del dynac.gnu"
        endif
        CALL System(COMMAND)
        if (isave.eq.0) then
          OPEN(unit=50,file='dynac.gnu')
        else
          call fn
          if (icontr.eq.2) filenm(1:2)='s3' 
          if (icontr.eq.3) filenm(1:2)='s4' 
          if (icontr.eq.4) filenm(1:2)='s5' 
          filenm(18:21)='.gnu'
          paf(1:10)='savedplots'
          if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
            paf(11:11)="/"
          else
C WINDOWS
            paf(11:11)=backslash
          endif
          paf(12:32)=filenm(1:21)
          OPEN(unit=50,file=paf(1:32))
        endif
        write(50,1000)
        ytitle=0.99
        if(opsy.eq.1) then
C LINUX
        elseif(opsy.eq.3) then
C MAC        
          write(50,6000)
          ytitle=0.985
        else
C WINDOWS        
          write(50,7000)
        endif
        write(50,1202) 
        write(50,1001) title,ytitle
        write(50,1002) labels(1)
        write(50,1003) labels(2)
        write(50,1004) xmin,xmax
        write(50,1005) ymin,ymax
        write(50,1006) 
        write(50,1007) 
        if (icontr.eq.2) then
          if (isave.eq.0) then
            write(50,1012) backslash
            write(50,1013)
          else
            filenm(18:21)='.plt'
            write(50,2012) filenm,backslash
            write(50,2013) filenm
          endif
        endif
        if (icontr.eq.3) then
          if (isave.eq.0) then
            write(50,1014)
          else
            filenm(18:21)='.plt'
            write(50,3014) filenm
          endif
        endif
        if (icontr.eq.4) then
          if (isave.eq.0) then
            write(50,1014)
          else
            filenm(18:21)='.plt'
            write(50,4014) filenm
          endif
        endif
        write(50,1010) 
1000    format('set style data points')
1001    format('set label "',A80,'" at screen 0.1 ,',f5.3)
1002    format('set xlabel "',A40,'"')
1003    format('set ylabel "',A40,'"')
1004    format('set xrange [',f8.2,':',f8.2,']')
1005    format('set yrange [',f12.5,':',f12.5,']')
1006    format('set size 1., 1.')
1007    format('set samples 50')   
1008    format('plot "dynac.plt" using 1:2 with points pt 0 lc 0, ',
     *         A1)
1011    format('plot "dynac.plt" using 1:2 with points pt 0 lc 0')
1012    format('plot "dynac.plt" using 1:2 with lines, ',
     *         A1)
1013    format('"dynac.plt" using 3:4 with lines')
1014    format('plot "dynac.plt" using 1:2 with lines')
1009    format('"dynac.cnt" using 1:2 with lines')
1010    format('pause -1 "hit return to continue"')
1202    format('set nokey')
2012    format('plot "',a21,'" using 1:2 with lines, ',
     *         A1)
2013    format('"',a21,'" using 3:4 with lines')
3014    format('plot "',a21,'" using 1:2 with lines')
4014    format('plot "',a21,'" using 1:2 with lines')
6000    format('set terminal aqua title "DYNAC" size 900 500') 
7000    format('set term wxt size 900,500')     
        close (50)
        RETURN
        END
        SUBROUTINE wfile3(imax,x,xp,y,yp)  
        implicit real*8 (a-h,o-z)
        common/opsys/opsy
c This routine writes the data points which will be plotted by GNU to a file
        dimension x(250002),xp(250002),y(250002),yp(250002)
        CHARACTER*80 title,command
        command=""
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          command="rm -f dynac.plt"
        else
C WINDOWS
          command="if exist dynac.plt del dynac.plt"
        endif
        CALL System(COMMAND)
        OPEN(unit=48,file='dynac.plt')
        DO i=1,imax
          write(48,*) x(i),xp(i),y(i),yp(i)
        ENDDO
        CLOSE(48)
        RETURN
        END
        SUBROUTINE wfile10(icont,imax,x,xp,y,yp,z,zp,igrtyp)  
        implicit real*8 (a-h,o-z)
c This routine writes the data points which will be plotted by GNU to a file
        character*5 cccst(20)
        common/fichier/filenm,pfnm
        common/opsys/opsy
        common/chstat1/cst(250002),cstat(20),fcstat(20),ncstat,mcstat
        common/chstat2/cccst
        dimension x(250002),xp(250002),y(250002),yp(250002)
        dimension z(250002),zp(250002)
        CHARACTER*80 title,command
        CHARACTER*11 pfnm(20)
        character filenm*23
        if (icont.eq.1) then
c store ellips contour
          command=""
          if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
            command="rm -f dynac.cnt"
          else
C WINDOWS
            command="if exist dynac.cnt del dynac.cnt"
          endif
          CALL System(COMMAND)
          OPEN(unit=51,file='dynac.cnt')
          DO i=1,201
            write(51,*) x(i),xp(i),y(i),yp(i),z(i),zp(i)
          ENDDO
          CLOSE(51)
        endif
        if (icont.eq.0) then
c store particle coordinates 
          if(igrtyp.eq.1) then
            command=""
            if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
              command="rm -f dynac.plt"
            else
C WINDOWS
              command="if exist dynac.plt del dynac.plt"
            endif
            CALL System(COMMAND)
            OPEN(unit=52,file='dynac.plt')
            DO i=1,imax
              write(52,*) x(i),xp(i),y(i),yp(i),z(i),zp(i)
            ENDDO
            CLOSE(52)
          else
            mcstat=0
            DO j=1,ncstat
              klm=0              
              DO i=1,imax
                if(cst(i).eq.cstat(j)) then
                  klm=klm+1
                endif
              ENDDO
              if(igrtyp.eq.6) write(6,111) klm,cstat(j)
111           format(i6,' particles with charge state ',f4.1)
              if(igrtyp.eq.11) then
                if(j.eq.1) write(6,222) klm,cstat(j)
222           format(i6,' particles originally within zone ',
     *               'delimited by  0.0 and ',f4.1,'*RMS')
                if(j.gt.1 .and. j.lt.ncstat)
     *               write(6,223)klm,cstat(j-1),cstat(j)
223           format(i6,' particles originally within zone ',
     *               'delimited by ',f4.1,' and ',f4.1,'*RMS')
                if(j.eq.ncstat) write(6,224) klm,cstat(j-1)
224           format(i6,' particles beyond ',f4.1,'*RMS')
              endif
              if(klm.ne.0) then
                mcstat=mcstat+1
                iunit=20+j
                if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                  command(1:6)="rm -f "
                  command(7:17)=pfnm(j)
                  CALL System(COMMAND(1:17))
                else
C WINDOWS
                  command(1:9)="if exist "
                  command(10:20)=pfnm(j)
                  command(21:25)=" del "  
                  command(26:36)=pfnm(j)              
                  CALL System(COMMAND(1:36))
                endif
                fcstat(mcstat)=cstat(j)
                OPEN(unit=iunit,file=pfnm(j))
                DO i=1,imax
                  if(cst(i).eq.cstat(j)) then
                    write(iunit,*) x(i),xp(i),y(i),yp(i),z(i),zp(i)
                  endif
                ENDDO
                CLOSE(iunit)
              endif
            ENDDO
          endif
        endif
        RETURN
        END
        SUBROUTINE wfile11(imax,x,xp,y,yp,igrtyp)  
        implicit real*8 (a-h,o-z)
        character*5 cccst(20)
        common/fichier/filenm,pfnm
        common/opsys/opsy
        common/chstat1/cst(250002),cstat(20),fcstat(20),ncstat,mcstat
        common/chstat2/cccst
c This routine writes the data points which will be plotted by GNU to a file
        dimension x(250002),xp(250002),y(250002),yp(250002)
        CHARACTER*80 title,command
        CHARACTER*11 pfnm(20)
        character filenm*23
c store particle coordinates    
        IF(igrtyp.eq.2) then
          command=""
          if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
            command="rm -f dynac.plt"
          else
C WINDOWS
            command="if exist dynac.plt del dynac.plt"
          endif
          CALL System(COMMAND)
          OPEN(unit=52,file='dynac.plt')  
          DO i=1,imax
            write(52,*) x(i),xp(i),y(i),yp(i)
          ENDDO
          CLOSE(52)
        ELSE
          mcstat=0
          DO j=1,ncstat
            klm=0              
            DO i=1,imax
              if(cst(i).eq.cstat(j)) then
                klm=klm+1
              endif
            ENDDO
            if(igrtyp.eq.7) write(6,111) klm,cstat(j)
111         format(i6,' particles with charge state ',f3.0)
            if(igrtyp.eq.12) then
              if(j.eq.1) write(6,222) klm,cstat(j)
222           format(i6,' particles originally within zone ',
     *               'delimited by  0.0 and ',f4.1,'*RMS')
                if(j.gt.1 .and. j.lt.ncstat)
     *               write(6,223)klm,cstat(j-1),cstat(j)
223           format(i6,' particles originally within zone ',
     *               'delimited by ',f4.1,' and ',f4.1,'*RMS')
              if(j.eq.ncstat) write(6,224) klm,cstat(j-1)
224           format(i6,' particles beyond ',f4.1,'*RMS')
            endif
            if(klm.ne.0) then
              mcstat=mcstat+1
              iunit=20+j
              if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
                command(1:6)="rm -f "
                command(7:17)=pfnm(j)
                CALL System(COMMAND(1:17))
              else
C WINDOWS
                command(1:9)="if exist "
                command(10:20)=pfnm(j)
                command(21:25)=" del "  
                command(26:36)=pfnm(j)              
                CALL System(COMMAND(1:36))
              endif
              fcstat(mcstat)=cstat(j)
              OPEN(unit=iunit,file=pfnm(j))
              DO i=1,imax
                if(cst(i).eq.cstat(j)) then
                  write(iunit,*) x(i),xp(i),y(i),yp(i)
                endif
              ENDDO
              CLOSE(iunit)
            endif
          ENDDO
        ENDIF
        RETURN
        END
        SUBROUTINE wfile21(isave,igrtyp)
        implicit real*8 (a-h,o-z)
        character*5 cccst(20)
        common/prtcnt/imax
        common/wfil2/title,labels
        common/fichier/filenm,pfnm
        common/opsys/opsy
        common/chstat1/cst(250002),cstat(20),fcstat(20),ncstat,mcstat
        common/chstat2/cccst
        common/pscl/yminsk,iskale
        CHARACTER*80 command,title
        CHARACTER*40 labels(20),labels3,labels4
        CHARACTER*11 pfnm(20),strng,fnm
        CHARACTER*6 parcnt
        CHARACTER backslash,filenm*23,paf*33,indx*2,indxx*3
        common/wfil20/xmin(10),xmax(10),ymin(10),ymax(10)
        parameter (backslash="\\")
        command=""
        write(6,7) imax
7       format(i6,' particles total')
        write(parcnt,'(I6)') imax
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          command="rm -f dynac.gnu"
        else
C WINDOWS
          command="if exist dynac.gnu del dynac.gnu"
        endif
        CALL System(COMMAND)
        if(isave.eq.0) then
          OPEN(unit=50,file='dynac.gnu')
        else
          call fn
          filenm(1:2)='s2' 
          filenm(18:21)='.gnu'
          paf(1:10)='savedplots'
          if(opsy.eq.1) then
C LINUX
            paf(11:11)="/"
          elseif(opsy.eq.3) then
C MAC
            paf(11:11)="/"
          else
C WINDOWS
            paf(11:11)=backslash
          endif
          paf(12:32)=filenm(1:21)
          OPEN(unit=50,file=paf(1:32))
        endif
        write(50,1000)
        ytitle=0.99
        if(opsy.eq.1) then
C LINUX
        elseif(opsy.eq.3) then
C MAC        
          write(50,6000)
          ytitle=0.985
        else
C WINDOWS        
          write(50,7000)
        endif
        if(ncstat.eq.1) then
          write(50,1203)
        else
          if(opsy.eq.1 .or. opsy.eq.3) then
            if(igrtyp.eq.12) then
              write(50,1202)
            else
              write(50,1204)
            endif
          else
            if(igrtyp.eq.12) then
              write(50,1214)
            else
              write(50,1202)
            endif
          endif
        endif
        write(50,1001) title,ytitle
        if(opsy.eq.1 .or. opsy.eq.3) then
C LINUX or MAC
          write(50,1510) parcnt
        else
C WINDOWS
          write(50,1500) parcnt
        endif
        write(50,1200)   
c x-z
        write(50,1006)
        write(50,1101)   
        write(50,1002) labels(1)
        write(50,1003) labels(2)
        write(50,1004) xmin(1),xmax(1)
        write(50,1005) ymin(1),ymax(1)  
        if(igrtyp.eq.2) then
          if(isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8008)
            else
              write(50,1008)
            endif
          else
            filenm(18:21)='.plt'
            if(opsy.eq.3) then
              write(50,8808) filenm(1:21)
            else
              write(50,2008) filenm(1:21)
            endif
          endif
        else
          strng=pfnm(1)
          indxx=cccst(1)
          if(isave.eq.0) then
            write(50,3003) pfnm(1),indxx,backslash
            do j=2,mcstat-1
              strng=pfnm(j)
              indxx=cccst(j)
              write(50,3033) pfnm(j),indxx,j,backslash
            enddo
            strng=pfnm(mcstat)
            if(igrtyp.eq.12) then
              indxx=cccst(mcstat-1)
              write(50,3013) pfnm(j),indxx,mcstat
            else
              indxx=cccst(mcstat)
              write(50,3043) pfnm(j),indxx,mcstat
            endif
          else
            filenm(18:23)=strng(6:11)
            write(50,5003) filenm,indxx,backslash
            do j=2,mcstat-1
              strng=pfnm(j)
              indxx=cccst(j)
              filenm(18:23)=strng(6:11)
              write(50,5033) filenm,indxx,j,backslash
            enddo
            strng=pfnm(mcstat)
            filenm(18:23)=strng(6:11)
            if(igrtyp.eq.12) then
              indxx=cccst(mcstat-1)
              write(50,5013) filenm,indxx,mcstat
            else
              indxx=cccst(mcstat)
              write(50,5043) filenm,indxx,mcstat
            endif
          endif
        endif
        write(50,1203)
        write(50,1300)
c y-z
        write(50,1006)
        write(50,1103)   
        write(50,1002) labels(3)
        write(50,1003) labels(4)
        write(50,1004) xmin(2),xmax(2)
        write(50,1005) ymin(2),ymax(2)
        if(igrtyp.eq.2) then
          if(isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8010)
            else
              write(50,1010)
            endif
          else
            filenm(18:21)='.plt'
            if(opsy.eq.3) then
              write(50,8810) filenm(1:21)
            else
              write(50,2010) filenm(1:21)
            endif
          endif
        else
          if(isave.eq.0) then
            write(50,4003) pfnm(1),backslash
            do j=2,mcstat-1
              write(50,4033) pfnm(j),j,backslash
            enddo
            write(50,4013) pfnm(mcstat),mcstat
          else
            strng=pfnm(1)
            filenm(18:23)=strng(6:11)
            write(50,6003) filenm,backslash
            do j=2,mcstat-1
              strng=pfnm(j)
              filenm(18:23)=strng(6:11)
              write(50,6033) filenm,j,backslash
            enddo
            strng=pfnm(mcstat)
            filenm(18:23)=strng(6:11)
            write(50,6013) filenm,mcstat
          endif
        endif
c x,y,z profiles
c new start
        labels3="X, Y, Z (RMS multiples)"
        labels4="N (normalized)"
        xmin2=-5.
        xmax2=5.
        ymin2=0.
        if(iskale.eq.1) then
          if (yminsk.gt.0. .and. yminsk.lt.1.) then
            ymin2=yminsk
          else
            write(6,*) '***         Error in logscale          ***'
            write(6,*) '** log scale minimum defaults to 1.E-06 **'
            ymin2=1.e-6
          endif
        endif
        ymax2=1.
c new end
        write(50,1006)
        write(50,1105)
        write(50,1002) labels3
        write(50,1003) labels4
        if(iskale.eq.1) then
          write(50,9000)
          write(50,9012)
        endif
        write(50,1004) xmin2,xmax2
        write(50,1005) ymin2,ymax2
        if(opsy.eq.1) then
          write(50,1208)
        elseif(opsy.eq.3) then
          write(50,1207)
        else
          write(50,1206)
        endif
        if(isave.eq.0) then
          fnm=pfnm(1)
          fnm(9:11)='pro'
          indx=' X'
          write(50,4303) fnm,indx,backslash
          fnm=pfnm(2)
          fnm(9:11)='pro'
          indx=' Y'
          write(50,4305) fnm,indx,backslash
          fnm=pfnm(3)
          fnm(9:11)='pro'
          indx=' Z'
          write(50,4313) fnm,indx
        else
          strng=pfnm(1)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx=' X'
          write(50,6303) filenm,indx,backslash
          strng=pfnm(2)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx=' Y'
          write(50,6305) filenm,indx,backslash
          strng=pfnm(3)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx=' Z'
          write(50,6313) filenm,indx
        endif
        write(50,1203)
        if(iskale.eq.1) then
          write(50,9010)
          write(50,9014)
        endif
c
c xp,yp,zp profiles
c new start
        labels3="Xp, Yp, Zp (RMS multiples)"
        labels4="N (normalized)"
        xmin2=-5.
        xmax2=5.
        ymin2=0.
        if(iskale.eq.1) then
          if (yminsk.gt.0. .and. yminsk.lt.1.) then
            ymin2=yminsk
          else
            ymin2=1.e-6
          endif
        endif
        ymax2=1.
c new end
        write(50,1006)
        write(50,1107)
        write(50,1002) labels3
        write(50,1003) labels4
        if(iskale.eq.1) then
          write(50,9000)
          write(50,9012)
        endif
        write(50,1004) xmin2,xmax2
        write(50,1005) ymin2,ymax2
        if(opsy.eq.1) then
          write(50,1212)
        elseif(opsy.eq.3) then
          write(50,1211)
        else
          write(50,1210)
        endif
        if(isave.eq.0) then
          fnm=pfnm(4)
          fnm(9:11)='pro'
          indx='Xp'
          write(50,4303) fnm,indx,backslash
          fnm=pfnm(5)
          fnm(9:11)='pro'
          indx='Yp'
          write(50,4305) fnm,indx,backslash
          fnm=pfnm(6)
          fnm(9:11)='pro'
          indx='Zp'
          write(50,4313) fnm,indx
        else
          strng=pfnm(4)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx='Xp'
          write(50,6303) filenm,indx,backslash
          strng=pfnm(5)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx='Yp'
          write(50,6305) filenm,indx,backslash
          strng=pfnm(6)
          filenm(18:20)=strng(6:8)
          filenm(21:23)='pro'
          indx='Zp'
          write(50,6313) filenm,indx
        endif
        write(50,1203)
        if(iskale.eq.1) then
          write(50,9010)
          write(50,9014)
        endif
c
        write(50,1201)
        write(50,9020)        
1000    format('set style data points',/,'set pointsize 0.01')
1001    format('set label "',A80,'" at screen 0.13 ,',f5.3)
1500    format('set label "',A6,' particles" at screen 0.45,0.52')
1510    format('set label "',A6,' particles" at screen 0.45,0.49')
1300    format('set nolabel')
1002    format('set xlabel "',A40,'"')
1003    format('set ylabel "',A40,'"')
1004    format('set xrange [',f8.2,':',f8.2,']')
1005    format('set yrange [',f12.6,':',f12.6,']')
1006    format('set size 0.5,0.5')
1100    format('set size 1.0, 1.0')
1101    format('set origin 0.,0.5')
1103    format('set origin 0.5,0.5')
1105    format('set origin 0.,0.')
1107    format('set origin 0.5,0.')
1007    format('set samples 50')   
1008    format('plot "dynac.plt" using 1:2 with points pt 0 lc 0')
8008    format('plot "dynac.plt" using 1:2 with points pt 0 lc 9')
1010    format('plot "dynac.plt" using 3:4 with points pt 0 lc 0')
8010    format('plot "dynac.plt" using 3:4 with points pt 0 lc 9')
1200    format('set multiplot')
1201    format('set nomultiplot')
1202    format('set key at screen 0.54, 0.95 spacing 0.8', 
     *         ' samplen 1 textcolor rgb variable ')
1204    format('set key at screen 0.52, 0.90 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1206    format('set key at screen 0.535, 0.45 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1207    format('set key at screen 0.535, 0.475 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1208    format('set key at screen 0.5, 0.40 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1210    format('set key at screen 0.535, 0.35 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1211    format('set key at screen 0.535, 0.375 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1212    format('set key at screen 0.5, 0.30 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1214    format('set key at screen 0.56, 0.95 spacing 0.8', 
     *         ' samplen 1 textcolor rgb variable ')
1203    format('unset key')
2008    format('plot "',a21,'" using 1:2 with points pt 0 lc 0')
2010    format('plot "',a21,'" using 3:4 with points pt 0 lc 0')
8808    format('plot "',a21,'" using 1:2 with points pt 0 lc 9')
8810    format('plot "',a21,'" using 3:4 with points pt 0 lc 9')
3003    format('plot "',a11,'" using 1:2 title "',A3,
     *         '" with points pt 1, ',A1)
3033    format('     "',a11,'" using 1:2 title "',A3,
     *         '" with points pt ',I2,',',A1)
3013    format('     "',a11,'" using 1:2 title " >',A3,
     *         '" with points pt ',I2)
3043    format('     "',a11,'" using 1:2 title " ',A3,
     *         '" with points pt ',I2)
4003    format('plot "',a11,'" using 3:4',
     *         ' with points pt 1, ',A1)
4303    format('plot "',a11,'" using 1:2 title "',A2,
     *         '" with lines, ',A1)
4305    format('     "',a11,'" using 1:2 title "',A2,
     *         '" with lines ls 2, ',A1)
4313    format('     "',a11,'" using 1:2 title "',A2,
     *         '" with lines ls 3')
4033    format('     "',a11,'" using 3:4',
     *         ' with points pt ',I2,',',A1)
4013    format('     "',a11,'" using 3:4',
     *         ' with points pt ',I2)
5003    format('plot "',a23,'" using 1:2 title "',A3,
     *         '" with points pt 1, ',A1)
5033    format('     "',a23,'" using 1:2 title "',A3,
     *         '" with points pt ',I2,',',A1)
5013    format('     "',a23,'" using 1:2 title " >',A3,
     *         '" with points pt ',I2)
5043    format('     "',a23,'" using 1:2 title " ',A3,
     *         '" with points pt ',I2)
6003    format('plot "',a23,'" using 3:4',
     *         ' with points pt 1, ',A1)
6033    format('     "',a23,'" using 3:4',
     *         ' with points pt ',I2,',',A1)
6013    format('     "',a23,'" using 3:4',
     *         ' with points pt ',I2)
6303    format('plot "',a23,'" using 1:2 title "',A2,
     *         '" with lines, ',A1)
6305    format('     "',a23,'" using 1:2 title "',A2,
     *         '" with lines ls 2, ',A1)
6313    format('     "',a23,'" using 1:2 title "',A2,
     *         '" with lines ls 3')
6000    format('set terminal aqua title "DYNAC" size 750 675') 
7000    format('set term wxt size 750,675')     
9000    format('set logscale y')
9010    format('set nologscale y')
9012    format('set format y "%.0t.E%+02T"')
9014    format('set format y "%g"')
9020    format('pause -1 "hit return to continue"')
        close (50)
        RETURN
        END
        SUBROUTINE wfile20(isave,igrtyp)
        implicit real*8 (a-h,o-z)
        character*5 cccst(20)
        common/prtcnt/imax
        common/chstat1/cst(250002),cstat(20),fcstat(20),ncstat,mcstat
        common/chstat2/cccst
        common/wfil2/title,labels
        common/fichier/filenm,pfnm
        common/opsys/opsy
        CHARACTER*80 command,title
        CHARACTER*40 labels(20)
        CHARACTER*11 pfnm(20),strng
        CHARACTER*6 parcnt
        CHARACTER backslash,filenm*23,paf*33,indx*2,indxx*3
        common/wfil20/xmin(10),xmax(10),ymin(10),ymax(10)
        parameter (backslash="\\")
        command=""
        write(parcnt,'(I6)') imax
        write(6,7) imax
7       format(i6,' particles total')
        ytitle=0.993
        if(opsy.eq.1) then
C LINUX
          command="rm -f dynac.gnu"
        elseif(opsy.eq.3) then
C MAC
          command="rm -f dynac.gnu"
          ytitle=0.985          
        else
C WINDOWS
          command="if exist dynac.gnu del dynac.gnu"
        endif
        CALL System(COMMAND)
        IF (isave.eq.0) then
          OPEN(unit=50,file='dynac.gnu')
        ELSE
          call fn
          filenm(1:2)='s1' 
          filenm(18:21)='.gnu'
          paf(1:10)='savedplots'
          if(opsy.eq.1) then
C LINUX
            paf(11:11)="/"
          elseif(opsy.eq.3) then
C MAC
            paf(11:11)="/"
          else
C WINDOWS
            paf(11:11)=backslash
          endif
          paf(12:32)=filenm(1:21)
          OPEN(unit=50,file=paf(1:32))
        ENDIF
        write(50,1000)
        if(mcstat.eq.1) then
          write(50,1203)
        else
          if(opsy.eq.1 .or. opsy.eq.3) then
            if(igrtyp.eq.11) then
              write(50,1202)
            else
              write(50,1204)
            endif
          else
            if(igrtyp.eq.11) then
              write(50,1206)
            else
              write(50,1202)
            endif
          endif
        endif
        write(50,1100)   
        write(50,1001) title,ytitle
        if(opsy.eq.1) then
C LINUX
          write(50,1510) parcnt
        elseif(opsy.eq.3) then
C MAC        
          write(50,1510) parcnt
          write(50,6000)
        else
C WINDOWS
          write(50,1500) parcnt
          write(50,7000)
        endif
        write(50,1200)   
c x-xp
        write(50,1006)
        write(50,1101)   
        write(50,1002) labels(1)
        write(50,1003) labels(2)
        write(50,1004) xmin(1),xmax(1)
        write(50,1005) ymin(1),ymax(1)  
        IF(igrtyp.eq.1) then
          IF (isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8008) backslash
            else
              write(50,1008) backslash
            endif
            write(50,1009)
          ELSE
            filenm(18:21)='.plt'
            write(50,2022) filenm(1:21),backslash
            filenm(18:21)='.cnt'
	      write(50,2019) filenm(1:21)
          ENDIF
        ELSE
          IF (isave.eq.0) then
            strng=pfnm(1)
            indxx=cccst(1)
            write(50,3001) pfnm(1),indxx,backslash
            do j=2,mcstat
              strng=pfnm(j)
              if(j.eq.mcstat .and. igrtyp.eq.11) then
                indxx=cccst(j-1)
                write(50,3018) pfnm(j),indxx,j,backslash
              else
                indxx=cccst(j)
                write(50,3008) pfnm(j),indxx,j,backslash
              endif
            enddo
            write(50,1009)
          ELSE
            strng=pfnm(1)
            indxx=cccst(1)
            filenm(18:23)=strng(6:11)
            write(50,4001) filenm,indxx,backslash
            do j=2,mcstat
              strng=pfnm(j)
              filenm(18:23)=strng(6:11)
              if(j.eq.mcstat .and. igrtyp.eq.11) then
                indxx=cccst(j-1)
                write(50,4018) filenm,indxx,j,backslash
              else
                indxx=cccst(j)
                write(50,4008) filenm,indxx,j,backslash
              endif
            enddo
            filenm(18:21)='.cnt'
   	      write(50,2019) filenm(1:21)
          ENDIF
        ENDIF
        write(50,1203)
        write(50,1300)
c y-yp
        write(50,1006)
        write(50,1102)   
        write(50,1002) labels(3)
        write(50,1003) labels(4)
        write(50,1004) xmin(2),xmax(2)
        write(50,1005) ymin(2),ymax(2)
        IF(igrtyp.eq.1) then        
c        IF(ncstat.eq.1) then
          IF (isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8011) backslash
            else
              write(50,1011) backslash
            endif
            write(50,1012) 
          ELSE
            filenm(18:21)='.plt'
            write(50,2031) filenm(1:21),backslash
            filenm(18:21)='.cnt'
            write(50,2042) filenm(1:21)
          ENDIF
        ELSE
          IF (isave.eq.0) then
            write(50,3002) pfnm(1),backslash
            do j=2,mcstat
              write(50,3011) pfnm(j),j,backslash
            enddo
            write(50,1012) 
          ELSE
            strng=pfnm(1)
            filenm(18:23)=strng(6:11)
            write(50,4002) filenm,backslash
            do j=2,mcstat
              strng=pfnm(j)
              filenm(18:23)=strng(6:11)
              write(50,4011) filenm,j,backslash
            enddo
            filenm(18:21)='.cnt'
            write(50,2042) filenm(1:21)
          ENDIF
        ENDIF
c x-y
        write(50,1006)
        write(50,1103)   
        write(50,1002) labels(5)
        write(50,1003) labels(6)
        write(50,1004) xmin(3),xmax(3)
        write(50,1005) ymin(3),ymax(3)
        IF(igrtyp.eq.1) then
          IF (isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8013)
            else
              write(50,1013)
            endif
          ELSE
            filenm(18:21)='.plt'
            write(50,2033) filenm(1:21)
          ENDIF
        ELSE
          IF (isave.eq.0) then
            write(50,3003) pfnm(1),backslash
            do j=2,mcstat-1
              write(50,3033) pfnm(j),j,backslash
            enddo
            write(50,3013) pfnm(mcstat),mcstat
          ELSE
            strng=pfnm(1)
            filenm(18:23)=strng(6:11)
            write(50,4003) filenm,backslash
            do j=2,mcstat-1
              strng=pfnm(j)
              filenm(18:23)=strng(6:11)
              write(50,4033) filenm,j,backslash
            enddo
            strng=pfnm(mcstat)
            filenm(18:23)=strng(6:11)
            write(50,4013) filenm,mcstat
          ENDIF
        ENDIF
c dW-dPHI
        write(50,1006)
        write(50,1104)   
        write(50,1002) labels(7)
        write(50,1003) labels(8)
        write(50,1004) xmin(4),xmax(4)
        write(50,1005) ymin(4),ymax(4)
        IF(igrtyp.eq.1) then
          IF (isave.eq.0) then
            if(opsy.eq.3) then
              write(50,8014) backslash
            else
              write(50,1014) backslash
            endif
            write(50,1015) 
          ELSE
            filenm(18:21)='.plt'
            write(50,2024) filenm(1:21),backslash
            filenm(18:21)='.cnt'
            write(50,2025) filenm(1:21)
          ENDIF
        ELSE
          IF (isave.eq.0) then
            write(50,3004) pfnm(1),backslash
            do j=2,mcstat
              write(50,3014) pfnm(j),j,backslash
            enddo
            write(50,1015) 
          ELSE
            strng=pfnm(1)
            filenm(18:23)=strng(6:11)
            write(50,4004) filenm,backslash
            do j=2,mcstat
              strng=pfnm(j)
              filenm(18:23)=strng(6:11)
              write(50,4014) filenm,j,backslash
            enddo
            filenm(18:21)='.cnt'
	    write(50,2025) filenm(1:21)
          ENDIF
        ENDIF
        write(50,1201)
        write(50,1020) 
1000    format('set style data points',/,'set pointsize 0.01')
1001    format('set label "',A80,'" at screen 0.13 ,',F5.3)
1500    format('set label "',A6,' particles" at screen 0.45,0.52')
1510    format('set label "',A6,' particles" at screen 0.45,0.49')
1300    format('set nolabel')
1002    format('set xlabel "',A40,'"')
1003    format('set ylabel "',A40,'"')
1004    format('set xrange [',f8.2,':',f8.2,']')
1005    format('set yrange [',f12.5,':',f12.5,']')
1006    format('set size 0.5, 0.5')
1100    format('set size 1.0, 1.0')
1101    format('set origin 0.,0.5')
1102    format('set origin 0.5,0.5')
1103    format('set origin 0.,0.')
1104    format('set origin 0.5,0.')
1007    format('set samples 50')   
1008    format('plot "dynac.plt" using 1:2 title "" with ',
     *         'points pt 0 lc 0, ',A1)
8008    format('plot "dynac.plt" using 1:2 title "" with ',
     *         'points pt 0 lc 9, ',A1)
1009    format('     "dynac.cnt" using 1:2 title "" with lines')
1011    format('plot "dynac.plt" using 3:4 with',
     *         ' points pt 0 lc 0, ',A1)
8011    format('plot "dynac.plt" using 3:4 with',
     *         ' points pt 0 lc 9, ',A1)
1012    format('     "dynac.cnt" using 3:4 with lines')
1013    format('plot "dynac.plt" using 1:3',
     *         ' with points pt 0 lc 0')
8013    format('plot "dynac.plt" using 1:3',
     *         ' with points pt 0 lc 9')
1014    format('plot "dynac.plt" using 5:6 with',
     *         ' points pt 0 lc 0, ',A1)
8014    format('plot "dynac.plt" using 5:6 with',
     *         ' points pt 0 lc 9, ',A1)
1015    format('     "dynac.cnt" using 5:6 with lines')
1020    format('pause -1 "hit return to continue"')
1200    format('set multiplot')
1201    format('set nomultiplot')
1202    format('set key at screen 0.54, 0.95 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1204    format('set key at screen 0.52, 0.98 spacing 0.8 maxcols 1',
     *         ' samplen 1 horizontal textcolor rgb variable ')
1206    format('set key at screen 0.56, 0.95 spacing 0.8',
     *         ' samplen 1 textcolor rgb variable ')
1203    format('set nokey')
2002    format('plot "',a23,'" using 1:2 with points pt 0, ',
     *         A1)
2022    format('plot "',a21,'" using 1:2 title "" with ',
     *         'points pt 0, ',A1)
2008    format('     "',a23,'" using 1:2 with points pt 1, ',
     *         A1)
2009    format('     "',a23,'" using 1:2 title "" with lines')
2019    format('     "',a21,'" using 1:2 title "" with lines')
2011    format('plot "',a23,'" using 3:4 with',
     *         ' points pt 0, ',A1)
2031    format('plot "',a21,'" using 3:4 with',
     *         ' points pt 0, ',A1)
2012    format('     "',a23,'" using 3:4 with lines')
2042    format('     "',a21,'" using 3:4 with lines')
2013    format('plot "',a23,'" using 1:3',
     *         ' with points pt 0 lc 0')
2033    format('plot "',a21,'" using 1:3',
     *         ' with points pt 0 lc 0')
2014    format('plot "',a23,'" using 5:6 with',
     *         ' points pt 0, ',A1)
2015    format('     "',a23,'" using 5:6 with lines')
2024    format('plot "',a21,'" using 5:6 with',
     *         ' points pt 0, ',A1)
2025    format('     "',a21,'" using 5:6 with lines')
3001    format('plot "',a11,'" using 1:2 title "',A5,
     *         '" with points pt 1, ',A1)
3008    format('     "',a11,'" using 1:2 title "',A5,
     *         '" with points pt ',I2,',',A1)
3018    format('     "',a11,'" using 1:2 title " >',A3,
     *         '" with points pt ',I2,',',A1)
3002    format('plot "',a11,'" using 3:4 with',
     *         ' points pt 1, ',A1)
3011    format('     "',a11,'" using 3:4 with',
     *         ' points pt ',I2,',',A1)
3003    format('plot "',a11,'" using 1:3',
     *         ' with points pt 1, ',A1)
3033    format('     "',a11,'" using 1:3',
     *         ' with points pt ',I2,',',A1)
3013    format('     "',a11,'" using 1:3',
     *         ' with points pt ',I2)
3004    format('plot "',a11,'" using 5:6 with',
     *         ' points pt 1, ',A1)
3014    format('     "',a11,'" using 5:6 with',
     *         ' points pt ',I2,',',A1)
4001    format('plot "',a23,'" using 1:2 title "',A3,
     *         '" with points pt 1, ',A1)
4008    format('     "',a23,'" using 1:2 title "',A3,
     *         '" with points pt ',I2,',',A1)
4018    format('     "',a23,'" using 1:2 title " >',A3,
     *         '" with points pt ',I2,',',A1)
4002    format('plot "',a23,'" using 3:4 with',
     *         ' points pt 1, ',A1)
4011    format('     "',a23,'" using 3:4 with',
     *         ' points pt ',I2,',',A1)
4003    format('plot "',a23,'" using 1:3',
     *         ' with points pt 1, ',A1)
4033    format('     "',a23,'" using 1:3',
     *         ' with points pt ',I2,',',A1)
4013    format('     "',a23,'" using 1:3',
     *         ' with points pt ',I2)
4004    format('plot "',a23,'" using 5:6 with',
     *         ' points pt 1, ',A1)
4014    format('     "',a23,'" using 5:6 with',
     *         ' points pt ',I2,',',A1)
5000    format('set label "',A3,'" at screen ',f4.2,',',f4.2)
6000    format('set terminal aqua title "DYNAC" size 750 675') 
7000    format('set term wxt size 750,675')     
        close (50)
        RETURN
        END
        SUBROUTINE fn
        implicit real*8 (a-h,o-z)
        common/fichier/filenm,pfnm
        common/mingw/mg
        CHARACTER*11 pfnm(20)
        character iitime*30,filenm*23
        logical mg
        inttim=time8()
        iitime=ctime(inttim)
        if (mg) then
c using MINGW style gfortran 03/30/10 20:51:06 (10 is 2010)
          if(iitime(1:2).eq.'01')filenm(3:5)='Jan'
          if(iitime(1:2).eq.'02')filenm(3:5)='Feb'
          if(iitime(1:2).eq.'03')filenm(3:5)='Mar'
          if(iitime(1:2).eq.'04')filenm(3:5)='Apr'
          if(iitime(1:2).eq.'05')filenm(3:5)='May'
          if(iitime(1:2).eq.'06')filenm(3:5)='Jun'
          if(iitime(1:2).eq.'07')filenm(3:5)='Jul'
          if(iitime(1:2).eq.'08')filenm(3:5)='Aug'
          if(iitime(1:2).eq.'09')filenm(3:5)='Sep'
          if(iitime(1:2).eq.'10')filenm(3:5)='Oct'
          if(iitime(1:2).eq.'11')filenm(3:5)='Nov'
          if(iitime(1:2).eq.'12')filenm(3:5)='Dec'
          filenm(6:7)=iitime(4:5)
          filenm(8:9)='20'
          filenm(10:11)=iitime(7:8)          
          filenm(12:13)=iitime(10:11)
          filenm(14:15)=iitime(13:14)
          filenm(16:17)=iitime(16:17)
        else
c using standard gfortran Tue Mar 30 20:51:06 2010
          filenm(3:5)=iitime(5:7)
          filenm(6:7)=iitime(9:10)
          filenm(8:11)=iitime(21:24)
          filenm(12:13)=iitime(12:13)
          filenm(14:15)=iitime(15:16)
          filenm(16:17)=iitime(18:19)
        endif
c on a MAC 02 Apr is shown as 2 Apr (on windows as 02 Apr)
c fill the blank with the character '0' to avoid error message on file copy
        if(filenm(6:6).eq.' ') filenm(6:6)='0'
        return
        end
        subroutine mkfrmt(i,fmt)
        implicit real*8 (a-h,o-z)
c this routine makes a (variable) A format
c e.g. if the integer i=11, the character format fmt will
c be fmt=(A11)
        character fmt*6
        fmt(1:1)='('
        fmt(2:2)='A'
        fmt(3:6)=''
        if (i.lt.10) then
          fmt(3:3)=char(i+48)
          fmt(4:4)=')'
        elseif (i.lt.100) then
          j1=i/10
          j2=i-j1*10
          fmt(3:3)=char(j1+48)
          fmt(4:4)=char(j2+48)
          fmt(5:5)=')'
        elseif (i.lt.1000) then
          j1=i/100
          j=i-j1*100
          j2=j/10
          j3=j-j2*10
          fmt(3:3)=char(j1+48)
          fmt(4:4)=char(j2+48)
          fmt(5:5)=char(j3+48)
          fmt(6:6)=')'
        endif
        return
        end