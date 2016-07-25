       PROGRAM convert
       implicit real*8 (a-h,o-z)
c ******************************************************************
C Version 1.0  31-Oct-2013
C      THIS SOFTWARE WAS ORIGINALLY PRODUCED BY CERN/PS, CEN/SACLAY
C      AUTHORS :E. TANKE
c ******************************************************************
C Converts TRACEWIN binary particle distribution file to DYNAC ascii
C standard
C
C     Syntax:  convert input_file output_file
C
C     7 or 'in' is the unit corresponding to the TRACEWIN .dst file describing the
C     particle distribution in binary format
C     16 or 'iou' is the unit corresponding to the DYNAC .dst file describing the
C     particle distribution in ascii format
c ******************************************************************
       parameter(ncards=59,iptsz=1000000,maxcell=3000,maxcell1=3000)
       COMMON /CONSTA/ VL, PI, XMAT, RPEL, QST
       DATA VL,XMAT,RPEL,QST/2.99792458E10,938.27231,28.17938E-14,1./
       common/faisc/f(6,iptsz),imax,ngood
       common/tapes/in,ifile,meta
       character*80 text
c       character mychar
       logical ffound
c********************************************
       INTEGER buff(13)
       integer narg,i,status,imax
       integer*1 ishort
c       real*4 cur
       character*80 inarg,myarg(10),infiln,oufiln
       in=7
       iou=16
c get arguments from the command line
c format for convert:
c convert file1 file2 [-h]
c where: file1 is the TRACEWIN file (input to the convert program)
c where: file2 is the DYNAC file (output from the convert program)
c        -h  print help info
       DO
         call GET_COMMAND_ARGUMENT(narg,inarg)
         if(LEN_TRIM(inarg).eq.0) exit
         narg=narg+1
         myarg(narg)=TRIM(inarg)
       ENDDO
c       INPUT ARGUMENTS:
c ********************************************************************************
       ffound=.false.
c       write(6,*) myarg(1)
c       write(6,*) myarg(2)
c       write(6,*) '# of args ',narg
c       goto 20
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
             WRITE(6,3101)
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
3101   FORMAT('******** CONVERTC V1.0R1 (Beta), 31-Oct-2013 *********')
       if(.not.ffound) then
         write(6,*) 'Error: Input file name required'
         write(6,*) 'Type'
         write(6,*) 'convert -h'
         stop
       endif
c ****************************************************************************************************
c fromat in the TRACEWIN input file:
c 2xCHAR+INT(Np)+DOUBLE(Ib(mA))+DOUBLE(freq(MHz))+CHAR+
c Np x [6xDOUBLE(x(cm),x'(rad),y(cm),y'(rad),phi(rad),Energie(MeV))]+
c DOUBLE(mc2(MeV))
Comments:
c  CHAR is 1 byte long ,
c  INT is 4 bytes long,
c  DOUBLE is 8 bytes long.
c  Np is the number of particles,
c  Ib is the beam current,
c  freq is the bunch frequency,
c  mc2 is the particle rest mass.
       CALL STAT(myarg(2), buff, status)
       IF (status == 0) THEN
c         WRITE (6,*) 'Device ID:',buff(1)
c         WRITE (6,*) 'Inode number:',buff(2)
c         WRITE (6,*) 'File mode (octal):',buff(3)
c         WRITE (6,*) 'Number of links:',buff(4)
c         WRITE (6,*) 'Owner’s uid:',buff(5)
c         WRITE (6,*) 'Owner’s gid:',buff(6)
c         WRITE (6,*) 'Device where located:',buff(7)
         WRITE (6,*) 'Input File size:',buff(8)
c         WRITE (6,*)’Last access time:’, T30, A19)") CTIME(buff(9))
c         WRITE (6,*)’Last modification time’, T30, A19)") CTIME(buff(10))
c         WRITE (6,*)’Last status change time:’, T30, A19)") CTIME(buff(11))
c         WRITE (6,*) 'Preferred block size:',buff(12)
c         WRITE (6,*) 'No. of blocks allocated:',buff(13)
       END IF
       open(in,file=infiln,form='unformatted',recl=100,
     *      access='direct', status='unknown')
       READ(in, REC=1) ISHORT,ISHORT,imax,cur,freq,ISHORT
       close(in)
c       open(in,file=myarg(2),form='unformatted',recl=4800031,
       open(in,file=infiln,form='unformatted',recl=buff(8),
     *      access='direct', status='unknown')
       open(iou,file=oufiln,status='unknown')
       READ( 7, REC=1) ISHORT,ISHORT,imax,cur,freq,ISHORT,
     * (f(1,km),f(2,km),f(3,km),f(4,km),f(5,km),f(6,km),km=1,imax),
     *  energy
       write(6,*) imax,' sets of particle coordinates read'
       write(6,*) 'enter factor to multiply phase with'
       read(5,*) pf
       write(iou,*) imax,' 0. 0.'
       do k=1,imax
         write(iou,*)f(1,k),f(2,k),f(3,k),f(4,k),pf*f(5,k),f(6,k)
       enddo
       GOTO 10
c8      WRITE( *, * ) 'I/O error # ', N, ', on 7'
       goto 20
10     write(6,*) 'np,I,freq,Wrest=',imax,cur,freq,energy
20     close(in)
       close(iou)
       STOP
       END
