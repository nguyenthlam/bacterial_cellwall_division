MODULE MODS_GROW


CONTAINS

!===================================================================================

          SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
          END SUBROUTINE

!==================================================

   subroutine paras(nthreads,nsyn,n_edhold,ntrap,lenmax,newlengthave,oldtypmax,joldbond,k_g,l_g,ktheta,theta_0,k_p,l_p, &
              LSTART,LREP,LSWITCH,KSWITCH,pres,pi,delta,invdelta,beta,ksur,kgtase,ktpase,ltpase, &
              dreact,dreact2,kedase,ledase,l_edcap,kgttp,lgttp13,lgttp2,kgted,lgted,kside,lside,ktped,ltped, &
              kwall,klead,kpair,lpair,invmpbp,wthick,pthick,kaxis,sigma,pconstrict0,delx,invdelx,lcons, &
              lambda,prmbond,growth,p_acti,p_deact,pterm,p_glyin,p_glyout,p_transfast,nprint, &
              p_transslow,p_pepout,p_edhold,pcleave_spon,pedcap_spon,pedcap_inact,plyt,pgbelowbondmax)

   implicit none

   character chara*32

   integer n

   integer nthreads,nsyn,n_edhold,ntrap,lenmax,oldtypmax,joldbond,newlengthave,pgbelowbondmax,nprint
   real(kind=8)::k_g,l_g,ktheta,theta_0,k_p,l_p,LSTART,LREP,LSWITCH,KSWITCH,SLOPE,pres
   real(kind=8)::pi,delta,invdelta,beta
   real(kind=8)::ksur,kgtase,ktpase,ltpase,dreact,dreact2,kedase,ledase,kgttp,lgttp13,lgttp2
   real(kind=8)::kgted,lgted,kside,lside,kwall,klead,kpair,lpair,invmpbp,ktped,ltped,l_edcap
   real(kind=8)::wthick,pthick,kaxis,sigma,pconstrict0,delx,invdelx,lambda,lcons
   real(kind=8)::prmbond,growth
   real(kind=8)::p_acti,p_deact,pterm,p_glyin,p_glyout,p_transfast,p_transslow,p_pepout
   real(kind=8)::p_edhold,pcleave_spon,pedcap_spon,pedcap_inact,plyt


16 format(A32)

   open(1,file='paras.info')


   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='NTHREAD') exit
   end do

   read(1,*)nthreads

   do n=1,1000
      read(1,16)chara
      if(chara(1:2)=='PI') exit
   end do
   read(1,*)pi

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='DELTA') exit
   end do

   read(1,*)delta

   INVDELTA=1.0D0/DELTA

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='BETA') exit
   end do

   read(1,*)beta


!---------------------------------
   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='GLYSPRING') exit
   end do

   read(1,*)k_g

   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='GLYLENGTH') exit
   end do

   read(1,*)l_g

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='GLYBEND') exit
   end do

   read(1,*)ktheta

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='GLYTHETA') exit
   end do

   read(1,*)theta_0

   THETA_0=THETA_0*PI/180.0D0

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='PEPFORCE') exit
   end do

   read(1,*)k_p

   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='PEPLENGTH') exit
   end do

   read(1,*)l_p

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='PRESSURE') exit
   end do

   read(1,*)pres

   LSTART=1.0D0

   LREP=L_P-LSTART


   LSWITCH=0.9D0*LREP


   SLOPE=K_P*(LREP/4.0D0/(1.0D0-LSWITCH/LREP)**2-LREP/4.0D0+LSWITCH)

   KSWITCH=SLOPE/2.0D0/LSWITCH

!------------------------------------------


! Force constants for enzymes:

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='KSUR') exit
   end do

   read(1,*)ksur

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KGTASE') exit
   end do

   read(1,*)kgtase

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KTPASE') exit
   end do

   read(1,*)ktpase

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='LTPASE') exit
   end do

   read(1,*)ltpase


   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='DREACT') exit
   end do

   read(1,*)dreact

   DREACT2=DREACT*DREACT

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KEDASE') exit
   end do

   read(1,*)kedase

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='LEDASE') exit
   end do

   read(1,*)ledase

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='L_EDCAP') exit
   end do

   read(1,*)l_edcap


   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KGTTP') exit
   end do

   read(1,*)kgttp

   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='LGTTP_ONE') exit
   end do

   read(1,*)lgttp13

   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='LGTTP_TWO') exit
   end do

   read(1,*)lgttp2

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KGTED') exit
   end do

   read(1,*)kgted

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LGTED') exit
   end do

   read(1,*)lgted

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KTPED') exit
   end do

   read(1,*)ktped

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LTPED') exit
   end do

   read(1,*)ltped


   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KSIDE') exit
   end do

   read(1,*)kside

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LSIDE') exit
   end do

   read(1,*)lside

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KWALL') exit
   end do

   read(1,*)kwall

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KLEAD') exit
   end do

   read(1,*)klead

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KPAIR') exit
   end do

   read(1,*)kpair

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LPAIR') exit
   end do

   read(1,*)lpair

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='INVMPBP') exit
   end do

   read(1,*)invmpbp

!----------------------------------

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='WTHICK') exit
   end do

   read(1,*)wthick

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='PTHICK') exit
   end do

   read(1,*)pthick

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KAXIS') exit
   end do

   read(1,*)kaxis

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='SIGMA') exit
   end do

   read(1,*)sigma

   do n=1,1000
      read(1,16)chara
      if(chara(1:10)=='PCONSTRICT') exit
   end do

   read(1,*)pconstrict0

!   do n=1,1000
!      read(1,16)chara
!      if(chara(1:6)=='DPCONS') exit
!   end do

!   read(1,*)dpcons

!  constriction range

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LCONS') exit
   end do

   read(1,*)lcons

!  set grid to calculate constriction pressure

   delx=0.2d0

   invdelx=1.0d0/delx

!   do n=1,1000

!      pconstrictold(n)=pconstrict0*exp(-n*n*delx*delx/sigma/sigma)*invdelta

!   end do

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='LAMBDA') exit
   end do

   read(1,*)lambda

!-----------------------

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='JOLDBOND') exit
   end do

   read(1,*)joldbond

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='PRMBOND') exit
   end do

   read(1,*)prmbond

   do n=1,1000
      read(1,16)chara
      if(chara(1:9)=='OLDTYPMAX') exit
   end do

   read(1,*)oldtypmax

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='PGBELOW') exit
   end do

   read(1,*)pgbelowbondmax

!----------------------------------------

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='SYNCOMP') exit
   end do

   read(1,*)nsyn

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='GROWTH') exit
   end do

   read(1,*)growth

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='NPRINT') exit
   end do

   read(1,*)nprint


!----------------------------------------

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='P_ACTI') exit
   end do

   read(1,*)p_acti

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='P_DEACT') exit
   end do

   read(1,*)p_deact

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='PTERM') exit
   end do

   read(1,*)pterm

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='P_GLYIN') exit
   end do

   read(1,*)p_glyin

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='P_GLYOUT') exit
   end do

   read(1,*)p_glyout

   do n=1,1000
      read(1,16)chara
      if(chara(1:11)=='P_TRANSFAST') exit
   end do

   read(1,*)P_TRANSFAST

   do n=1,1000
      read(1,16)chara
      if(chara(1:11)=='P_TRANSSLOW') exit
   end do

   read(1,*)P_TRANSSLOW

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='P_PEPOUT') exit
   end do

   read(1,*)p_pepout

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='P_EDHOLD') exit
   end do

   read(1,*)p_edhold

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='N_EDHOLD') exit
   end do

   read(1,*)N_edhold

   do n=1,1000
      read(1,16)chara
      if(chara(1:12)=='PCLEAVE_SPON') exit
   end do

   read(1,*)PCLEAVE_SPON

   do n=1,1000
      read(1,16)chara
      if(chara(1:11)=='PEDCAP_SPON') exit
   end do

   read(1,*)PEDCAP_SPON

   do n=1,1000
      read(1,16)chara
      if(chara(1:12)=='PEDCAP_INACT') exit
   end do

   read(1,*)PEDCAP_INACT

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='PLYT') exit
   end do

   read(1,*)plyt

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='NTRAP') exit
   end do

   read(1,*)ntrap

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='LENMAX') exit
   end do

   read(1,*)lenmax

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='AVELEN') exit
   end do

   read(1,*)newlengthave




   end subroutine

!==================================================

   SUBROUTINE WHATTIME(TIMERUN)

   IMPLICIT NONE
   INTEGER TIMES(10),YEARS,MONS,TIMERUN,DAYS,HOURS,MINS,SECS,LASTTIME,TIMESTART,NM

   call date_and_time(values=times)

   mins=times(6); hours=times(5); days=times(3)-1; mons=times(2); years=times(1)
   timerun=mins+60*hours+60*24*days

   if(mons>1)then
      do nm=1,mons-1
         if((mod(nm,2)==1.and.nm<=7).or.(mod(nm,2)==0.and.nm>=8))then
            timerun=timerun+60*24*31
         elseif(nm==2.and.mod(years,4)==0)then
            timerun=timerun+60*24*29
         elseif(nm==2.and.mod(years,4)/=0)then
            timerun=timerun+60*24*28
         else
            timerun=timerun+60*24*30
         end if
      end do
   end if

   if(mod(years-1,4)==0)then
      timerun=timerun+60*24*366*(years-1)
   else
      timerun=timerun+60*24*365*(years-1)
   end if

   END SUBROUTINE


!========================================================================

   SUBROUTINE GETINFO(NSTART,NATOMOLD,NPGOLD,PGLENMAX,LOOPLENMAX,caplen)


   IMPLICIT NONE

   INTEGER N,NSTART,NATOMOLD,NPGOLD,PGLENMAX,LOOPLENMAX,caplen
   CHARACTER (LEN=64) FILECONFIG
   CHARACTER ZERO*1,CHARID1*1,CHARID2*2,CHARID3*3,CHARID4*4
   CHARACTER CHARA*512



   WRITE(ZERO,'(I1)')0

   IF(NSTART<10)THEN
      WRITE(CHARID1,'(I1)')NSTART
      FILECONFIG='config'//ZERO//ZERO//ZERO//CHARID1//'.dat'
   ELSE IF(NSTART<100)THEN
      WRITE(CHARID2,'(I2)')NSTART
      FILECONFIG='config'//ZERO//ZERO//CHARID2//'.dat'
   ELSE IF(NSTART<1000)THEN
      WRITE(CHARID3,'(I3)')NSTART
      FILECONFIG='config'//ZERO//CHARID3//'.dat'
   ELSE IF(NSTART<10000)THEN
      WRITE(CHARID4,'(I4)')NSTART
      FILECONFIG='config'//CHARID4//'.dat'
   ELSE
      PRINT*,'NSTART IS TOO BIG! STOP NOW.'
      STOP
   END IF

! TOPOLOGY:

   OPEN(1,FILE=FILECONFIG)

9  READ(1,*)CHARA
   IF(CHARA(1:4)/='ATOM')THEN
      GOTO 9
   END IF

   READ(1,*)NATOMOLD

3  READ(1,*)CHARA
   IF(CHARA(1:4)/='CAPS')THEN
      GOTO 3
   END IF

   READ(1,*)caplen

12 READ(1,*)CHARA
   IF(CHARA(1:5)/='OLDPG')THEN
      GOTO 12
   END IF

   READ(1,*)NPGOLD,PGLENMAX


7  READ(1,*)CHARA
   IF(CHARA(1:4)/='LOOP')THEN
      GOTO 7
   END IF

   READ(1,*)N,LOOPLENMAX

   CLOSE(1)

   END SUBROUTINE GETINFO

!========================================================================

   SUBROUTINE GETCAPS(NSTART,NCAP,CAPBOUND1,CAPBOUND2,CAPID,CAPLEN)

   IMPLICIT NONE

   INTEGER NSTART,NCAP,CAPBOUND1,CAPBOUND2,N,M,I
   INTEGER,DIMENSION(:),ALLOCATABLE::CAPLEN
   INTEGER,DIMENSION(:,:),ALLOCATABLE::CAPID

   CHARACTER CHARA*512
   CHARACTER (LEN=64) FILECONFIG
   CHARACTER ZERO*1,CHARID1*1,CHARID2*2,CHARID3*3,CHARID4*4

   WRITE(ZERO,'(I1)')0

   IF(NSTART<10)THEN
      WRITE(CHARID1,'(I1)')NSTART
      FILECONFIG='config'//ZERO//ZERO//ZERO//CHARID1//'.dat'
   ELSE IF(NSTART<100)THEN
      WRITE(CHARID2,'(I2)')NSTART
      FILECONFIG='config'//ZERO//ZERO//CHARID2//'.dat'
   ELSE IF(NSTART<1000)THEN
      WRITE(CHARID3,'(I3)')NSTART
      FILECONFIG='config'//ZERO//CHARID3//'.dat'
   ELSE IF(NSTART<10000)THEN
      WRITE(CHARID4,'(I4)')NSTART
      FILECONFIG='config'//CHARID4//'.dat'
   ELSE
      PRINT*,'NSTART IS TOO BIG! STOP NOW.'
      STOP
   END IF

   OPEN(1,FILE=FILECONFIG)

2  READ(1,*)CHARA
   IF(CHARA(1:4)/='CAPS')THEN
      GOTO 2
   END IF

   READ(1,*)M,N,CAPBOUND1,CAPBOUND2

   DO I=1,NCAP
      READ(1,*)CAPLEN(I),CAPID(1:CAPLEN(I),I)
   END DO

   CLOSE(1)

   END SUBROUTINE

!========================================================================


 SUBROUTINE PG_INPUT(NSTART,NATOMOLD,NATOM,NATOMDEL,OLDNATOMDEL,DNOR,ATOR,PEPDIR,X,Y,Z, &
              NPGOLD,PGID,PGLEN,NPG,NEWPGID,NEWPGLEN,NBONDGLY,NBONDGLYOLD,NBONDPEP,NBONDDEL,caplen, &
               noldbond,oldbond,oldtyp,PARTNER,BONDGLY,BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                jsyndir,na_red,NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDCAP,EDHOLD,CRLKAGE, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD, &
                  npress,pconstrict,delx,pconstrict0,sigma,invdelta,oldglytip,oldglysec,temlink,ntemlink)


   IMPLICIT NONE

   INTEGER NATOMOLD,NATOM,NSYN,NATOMDEL,OLDNATOMDEL,NPGOLD,NPG,length,n1,n2
   INTEGER NBONDGLY,NBONDGLYOLD,NBONDPEP,NBONDDEL,ntemlink
   INTEGER NLOOP,NLOOPDEL,noldbond
   INTEGER N,M,N0,I,J,K,IPG,LOOPLENMAX
   INTEGER NSTART,CAPlen,npress,jsyndir,na_red

   INTEGER,DIMENSION(:,:),ALLOCATABLE::PGID,NEWPGID,BONDGLY,BONDPEP,LOOP,SYNTHESIS,SYNPG,GLYTIP,GLYSEC
   INTEGER,DIMENSION(:),ALLOCATABLE::PGLEN,NEWPGLEN,BONTYP,DNOR,ATOR,PEPDIR,oldtyp,GTLOAD,GTATRANS
   INTEGER,DIMENSION(:,:),ALLOCATABLE::SIGCROSS,EDHOLD,CRLKAGE,EDPEP,PARTNER,oldbond,oldglytip,oldglysec,temlink
   INTEGER,DIMENSION(:),ALLOCATABLE::LOOPLEN,LOOPTYP,SYNDIR,SYNLOOP,JDEACT,SIGCLEAVE,EDCAP
   INTEGER,DIMENSION(:,:,:),ALLOCATABLE::TPPEP

   DOUBLE PRECISION,VALUE::delx,pconstrict0,sigma,invdelta

   DOUBLE PRECISION XCAP1,XCAP2,XCEN,X0

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: pconstrict

   CHARACTER CHARA*512
   CHARACTER (LEN=64) FILECOOR,FILECONFIG
   CHARACTER ZERO*1,CHARID1*1,CHARID2*2,CHARID3*3,CHARID4*4


   WRITE(ZERO,'(I1)')0

   IF(NSTART<10)THEN
      WRITE(CHARID1,'(I1)')NSTART
      FILECOOR='coor'//ZERO//ZERO//ZERO//CHARID1//'.dat'
      FILECONFIG='config'//ZERO//ZERO//ZERO//CHARID1//'.dat'
   ELSE IF(NSTART<100)THEN
      WRITE(CHARID2,'(I2)')NSTART
      FILECOOR='coor'//ZERO//ZERO//CHARID2//'.dat'
      FILECONFIG='config'//ZERO//ZERO//CHARID2//'.dat'
   ELSE IF(NSTART<1000)THEN
      WRITE(CHARID3,'(I3)')NSTART
      FILECOOR='coor'//ZERO//CHARID3//'.dat'
      FILECONFIG='config'//ZERO//CHARID3//'.dat'
   ELSE IF(NSTART<10000)THEN
      WRITE(CHARID4,'(I4)')NSTART
      FILECOOR='coor'//CHARID4//'.dat'
      FILECONFIG='config'//CHARID4//'.dat'
   ELSE
      PRINT*,'NSTART IS TOO BIG! STOP NOW.'
      STOP
   END IF

   PRINT*,'=================================================================='

   PRINT*,'STARTING UP THE SYSTEM AT NSTART =',NSTART
   WRITE(*,*)

   OPEN(1,FILE=FILECOOR,FORM='UNFORMATTED')


   READ(1)NATOM

   READ(1)X(1:NATOM)
   READ(1)Y(1:NATOM)
   READ(1)Z(1:NATOM)

   IF(NSTART>0)THEN
    READ(1)NSYN
    READ(1)XGTASE(1,1:NSYN),XGTASE(2,1:NSYN),XTPASE(1,1:NSYN),XTPASE(2,1:NSYN),XTPASE(3,1:NSYN),XEDASE(1:NSYN),XEDASEOLD(1:NSYN)
    READ(1)YGTASE(1,1:NSYN),YGTASE(2,1:NSYN),YTPASE(1,1:NSYN),YTPASE(2,1:NSYN),YTPASE(3,1:NSYN),YEDASE(1:NSYN),YEDASEOLD(1:NSYN)
    READ(1)ZGTASE(1,1:NSYN),ZGTASE(2,1:NSYN),ZTPASE(1,1:NSYN),ZTPASE(2,1:NSYN),ZTPASE(3,1:NSYN),ZEDASE(1:NSYN),ZEDASEOLD(1:NSYN)
   END IF


!  constriction pressure

!   if(nstart==0)then

      do n=1,npress-1

         pconstrict(n)=pconstrict0*exp(-n*n*delx*delx/sigma/sigma)*invdelta

      end do

      pconstrict(npress)=0.0d0

!   else

!      read(1)pconstrict(1:npress)

!   end if


   CLOSE(1)

   WRITE(*,*)'NUMBER OF ATOMS IN THE SYSTEM =',NATOM
   WRITE(*,*)

!-------------------------------------------------------
   OPEN(1,FILE=FILECONFIG)


42 READ(1,*)CHARA
   IF(CHARA(1:4)/='ATOM')THEN
      GOTO 42
   END IF

   READ(1,*)NATOMOLD,M,NATOMDEL,OLDNATOMDEL
   READ(1,*)DNOR(1:NATOM)
   READ(1,*)ATOR(1:NATOM)
   READ(1,*)PEPDIR(1:NATOM)

!-----------------------------------------------

!31 READ(1,*)CHARA
!   IF(CHARA(1:4)/='CAPS')THEN
!      GOTO 31
!   END IF

!   READ(1,*)caplen

!   read(1,*)capopen(1:caplen)
!   read(1,*)capopen(caplen+1:2*caplen)



!   if(nstart==0)then

      NBONDGLY=0

3     READ(1,*)CHARA
      IF(CHARA(1:5)/='OLDPG')THEN
         GOTO 3
      END IF

      READ(1,*)N,length

      DO I=1,NPGOLD
         READ(1,*)PGLEN(I),PGID(1:PGLEN(I),I)

         LENGTH=PGLEN(i)

         BONDGLY(1,NBONDGLY+1:NBONDGLY+LENGTH-1)=PGID(1:LENGTH-1,i)
         BONDGLY(2,NBONDGLY+1:NBONDGLY+LENGTH-1)=PGID(2:LENGTH,i)

         NBONDGLY=NBONDGLY+LENGTH-1

      END DO

!   else

!      NBONDGLY=NBONDGLYOLD

!   end if


33 READ(1,*)CHARA
   IF(CHARA(1:5)/='NEWPG')THEN
      GOTO 33
   END IF

   READ(1,*)NPG,length

   if(npg>0)then

      do i=1,npg
         read(1,*)newpglen(i),newpgid(1:newpglen(i),i)

         length=newpglen(i)

         IF(LENGTH>1)THEN

            BONDGLY(1,NBONDGLY+1:NBONDGLY+LENGTH-1)=newPGID(1:LENGTH-1,i)
            BONDGLY(2,NBONDGLY+1:NBONDGLY+LENGTH-1)=newPGID(2:LENGTH,i)

            NBONDGLY=NBONDGLY+LENGTH-1

         END IF

      end do

   end if



   WRITE(*,*)'NUMBER OF STRANDS ON THE ROD = ',NPGOLD+NPG
   WRITE(*,*)

   CALL FIXPEPDIR(NPG,newPGID,newPGLEN,PEPDIR)



!   NBONDGLY=0

!   do n=1,caplen-1

!      BONDGLY(1,n)=n
!      BONDGLY(2,n)=n+1

!   end do

!   NBONDGLY=caplen


!   BONDGLY(1,NBONDGLY)=caplen
!   BONDGLY(2,NBONDGLY)=1

!   do n=caplen+1,2*caplen-1

!      BONDGLY(1,N)=n
!      BONDGLY(2,N)=n+1

!   end do

!   NBONDGLY=NBONDGLY+caplen

!   BONDGLY(1,NBONDGLY)=caplen*2
!   BONDGLY(2,NBONDGLY)=1+caplen



!   DO N=1,NPG

!      LENGTH=PGLEN(N)

!      IF(LENGTH>1)THEN

!         BONDGLY(1,NBONDGLY+1:NBONDGLY+LENGTH-1)=PGID(1:LENGTH-1,N)
!         BONDGLY(2,NBONDGLY+1:NBONDGLY+LENGTH-1)=PGID(2:LENGTH,N)

!         NBONDGLY=NBONDGLY+LENGTH-1

!      END IF


!   END DO


!-----------------------------------------------

32 READ(1,*)CHARA
   IF(CHARA(1:7)/='BONDPEP')THEN
      GOTO 32
   END IF


   READ(1,*)NBONDPEP,NBONDDEL


11   FORMAT(I6,4X,I1,4X,I6,4X,I6)

   DO N=1,NBONDPEP
      READ(1,11)N0,BONTYP(N),BONDPEP(1,N),BONDPEP(2,N)
      IF(N/=N0)THEN
         PRINT*,'BONDS READING ERROR!','N=',N,'N0=',N0
         STOP
      END IF

      IF(BONTYP(N)>0)THEN

         N1=BONDPEP(1,N)
         N2=BONDPEP(2,N)

         IF(PARTNER(1,N1)==0)THEN
            PARTNER(1,N1)=N2
         ELSE
            PARTNER(2,N1)=N2
         END IF

         IF(PARTNER(1,N2)==0)THEN
            PARTNER(1,N2)=N1
         ELSE
            PARTNER(2,N2)=N1
         END IF


      END IF

   END DO

   if(nstart>0)then

      read(1,*)noldbond
      do n=1,noldbond
         read(1,*)n0,oldtyp(n),oldbond(1,n),oldbond(2,n)
      end do

   end if

!   WRITE(*,*)'NUMBER OF GLY-GLY BONDS =',NBONDGLY
   WRITE(*,*)'NUMBER OF GLY-PEP BONDS =',NBONDPEP

!------------------------------------------------
   WRITE(*,*)

7  READ(1,*)CHARA
   IF(CHARA(1:4)/='LOOP')THEN
      GOTO 7
   END IF

   READ(1,*)NLOOP,N,NLOOPDEL

   DO N=1,NLOOP

      READ(1,*)LOOPLEN(N),LOOPTYP(N),(LOOP(I,N),I=1,LOOPLEN(N))

   END DO

   WRITE(*,*)'NUMBER OF LOOPS =',NLOOP

!-------------------------------------------



!-------------------------------------------

12 FORMAT(2X,I2,2(2X,I1),10(2X,I7))


   IF(NSTART==0)THEN


      CALL ALIGNCELL(NATOM,X,Y,Z,CAPlen)



      PRINT*,'NUMBER OF COMPLEXES TO START =',NSYN

      CALL SYNSETUP(NATOM,X,Y,Z,NSYN,SYNDIR, &
                    XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE)


   ELSE

41    READ(1,*)CHARA
      IF(CHARA(1:7)/='SYNCOMP')THEN
         GOTO 41
      END IF

      READ(1,*)NSYN!,SYNRATIO   ! THIS ARE NUMBER OF SYNTHESIS COMPLEXES AND THE ATOM/SYNTHESIS RATIO

!      READ(1,*)
!      READ(1,*)COUPLING(1:NSYN)

      READ(1,*)
      READ(1,*)SYNDIR(1:NSYN)

      READ(1,*)
      READ(1,*)SYNTHESIS(1,1:NSYN),SYNTHESIS(2,1:NSYN)

      READ(1,*)
      READ(1,*)GTLOAD(1:NSYN)!,GTLOAD(2,1:NSYN)

      READ(1,*)
      READ(1,*)SYNPG(1,1:NSYN),SYNPG(2,1:NSYN)

      READ(1,*)
      READ(1,*)SYNLOOP(1:NSYN)

      DO N=1,NSYN

         IF(SYNLOOP(N)/=0)THEN
            LOOPTYP(SYNLOOP(N))=2
         END IF

      END DO

      READ(1,*)
      READ(1,*)GLYTIP(1,1:NSYN),GLYTIP(2,1:NSYN)

      READ(1,*)
      READ(1,*)GLYSEC(1,1:NSYN),GLYSEC(2,1:NSYN)

      READ(1,*)
      READ(1,*)TPPEP(1,1,1:NSYN),TPPEP(1,2,1:NSYN),TPPEP(1,3,1:NSYN)

      READ(1,*)
      READ(1,*)TPPEP(2,1,1:NSYN),TPPEP(2,2,1:NSYN),TPPEP(2,3,1:NSYN)

      READ(1,*)
      READ(1,*)EDPEP(1,1:NSYN),EDPEP(2,1:NSYN)

      READ(1,*)
      READ(1,*)GTATRANS(1:NSYN)!,GTATRANS(2,1:NSYN)

      READ(1,*)
      READ(1,*)JDEACT(1:NSYN)

      READ(1,*)
      READ(1,*)SIGCROSS(1,1:NSYN),SIGCROSS(2,1:NSYN),SIGCROSS(3,1:NSYN)

      READ(1,*)
      READ(1,*)SIGCLEAVE(1:NSYN)

!      READ(1,*)
!      READ(1,*)EDLOCKIN(1:NSYN)

      READ(1,*)
      READ(1,*)EDCAP(1:NSYN)

      READ(1,*)
      READ(1,*)EDHOLD(1,1:NSYN),EDHOLD(2,1:NSYN),EDHOLD(3,1:NSYN)

      READ(1,*)
      READ(1,*)CRLKAGE(1,1:NSYN),CRLKAGE(2,1:NSYN)


46    READ(1,*)CHARA
      IF(CHARA(1:7)/='TEMLINK')THEN
         GOTO 46
      END IF


      read(1,*)ntemlink
      do n=1,ntemlink
         read(1,*)temlink(1:4,n)
      end do

      read(1,*)oldglytip(1,1:nsyn),oldglytip(2,1:nsyn),oldglysec(1,1:nsyn),oldglysec(2,1:nsyn)

47    READ(1,*)CHARA
      IF(CHARA(1:5)/='REDIS')THEN
         GOTO 47
      END IF
  
      read(1,*)jsyndir,na_red

   END IF

   CLOSE(1)

!   DO N=1,NSYN

!      IF(SYNLOOP(N)/=0)THEN
!         LOOPTYP(SYNLOOP(N))=2
!      END IF

!   END DO

!-----------------------------------------------------------

print*,'done PG Input'

   END SUBROUTINE PG_INPUT

!==================================================================

   SUBROUTINE ALIGNCELL(NATOM,X,Y,Z,caplen)

   IMPLICIT NONE

   INTEGER NATOM,N,caplen
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z

!   INTEGER,DIMENSION(:,:),ALLOCATABLE::PGID
!   INTEGER,DIMENSION(:),ALLOCATABLE::PGLEN
   DOUBLE PRECISION XCAP1,XCAP2,YCAP1,YCAP2,ZCAP1,ZCAP2,DX,DY,DZ,SINP,COSP,XCEN,YCEN,ZCEN,XR,XA,YA,ZA




   XCAP1=SUM(X(1:caplen))/caplen
   YCAP1=SUM(Y(1:caplen))/caplen
   ZCAP1=SUM(Z(1:caplen))/caplen


   XCAP2=SUM(X(1+caplen:caplen*2))/caplen
   YCAP2=SUM(Y(1+caplen:caplen*2))/caplen
   ZCAP2=SUM(Z(1+caplen:caplen*2))/caplen



!----------------------------------------------
!  ROTATE ABOUT Y:

   DX=XCAP2-XCAP1
   DY=YCAP2-YCAP1
   DZ=ZCAP2-ZCAP1

   SINP=DZ/SQRT(DX*DX+DZ*DZ)
   COSP=DX/SQRT(DX*DX+DZ*DZ)

   DO N=1,NATOM

      XA=X(N)
      ZA=Z(N)

      X(N)=XA*COSP+ZA*SINP

      Z(N)=ZA*COSP-XA*SINP

   END DO

!----------------------------------------------
!  ROTATE ABOUT Z:

!   XCAP1=SUM(X(capid(1:caplen,1)))/caplen
!   ZCAP1=SUM(Z(capid(1:caplen,1)))/caplen

   XCAP1=SUM(X(1:caplen))/caplen
   ZCAP1=SUM(Z(1:caplen))/caplen


!   XCAP2=SUM(X(capid(1:caplen,2)))/caplen
!   ZCAP2=SUM(Z(capid(1:caplen,2)))/caplen

   XCAP2=SUM(X(1+caplen:caplen*2))/caplen
   ZCAP2=SUM(Z(1+caplen:caplen*2))/caplen


   DX=XCAP2-XCAP1
   DY=YCAP2-YCAP1

   SINP=-DY/SQRT(DX*DX+DY*DY)
   COSP=DX/SQRT(DX*DX+DY*DY)

   DO N=1,NATOM

      XA=X(N)
      YA=Y(N)

      X(N)=XA*COSP-YA*SINP

      Y(N)=YA*COSP+XA*SINP

   END DO

!----------------------------------------------
!  MOVE THE CENTER TO (0,0,0):

!   XCAP1=SUM(X(capid(1:caplen,1)))/caplen
!   YCAP1=SUM(Y(capid(1:caplen,1)))/caplen


!   XCAP2=SUM(X(capid(1:caplen,2)))/caplen
!   YCAP2=SUM(Y(capid(1:caplen,2)))/caplen

   XCAP1=SUM(X(1:caplen))/caplen
   YCAP1=SUM(Y(1:caplen))/caplen


   XCAP2=SUM(X(1+caplen:caplen*2))/caplen
   YCAP2=SUM(Y(1+caplen:caplen*2))/caplen


   XCEN=(XCAP1+XCAP2)/2
   YCEN=(YCAP1+YCAP2)/2
   ZCEN=(ZCAP1+ZCAP2)/2


   DO N=1,NATOM
      X(N)=X(N)-XCEN
      Y(N)=Y(N)-YCEN
      Z(N)=Z(N)-ZCEN
   END DO

!   XCAP1=SUM(X(capid(1:caplen,1)))/caplen
!   YCAP1=SUM(Y(capid(1:caplen,1)))/caplen
!   ZCAP1=SUM(Z(capid(1:caplen,1)))/caplen


!   XCAP2=SUM(X(capid(1:caplen,2)))/caplen
!   YCAP2=SUM(Y(capid(1:caplen,2)))/caplen
!   ZCAP2=SUM(Z(capid(1:caplen,2)))/caplen

!print*,xcap1,ycap1,zcap1
!print*,xcap2,ycap2,zcap2
!stop

   END SUBROUTINE


!=======================================================================

   SUBROUTINE RESETCAPS(NATOM,X,Y,Z,ORAD,PGID,PGLEN,CAPBOUND1,CAPBOUND2,XCAP1,XCAP2)

   IMPLICIT NONE

   INTEGER NATOM,CAPBOUND1,CAPBOUND2,N,J
   INTEGER COUNT1(10),COUNT2(10)
   INTEGER,DIMENSION(:),ALLOCATABLE::PGLEN
   INTEGER,DIMENSION(:,:),ALLOCATABLE::PGID

   DOUBLE PRECISION ORAD,XCAP1,XCAP2,DELT,X1,X2
   DOUBLE PRECISION RAD1(10),RAD2(10)
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z


   XCAP1=SUM(X(PGID(1:PGLEN(CAPBOUND1),CAPBOUND1)))/PGLEN(CAPBOUND1)
   XCAP2=SUM(X(PGID(1:PGLEN(CAPBOUND2),CAPBOUND2)))/PGLEN(CAPBOUND2)

   DELT=ORAD/10

   X1=XCAP1+ORAD

   X2=XCAP2-ORAD

   COUNT1=0

   COUNT2=0

   RAD1=0.0D0

   RAD2=0.0D0

   DO N=1,NATOM

      IF(X(N)<X1.AND.X(N)>XCAP1)THEN

         J=(X(N)-XCAP1)/DELT+1

         COUNT1(J)=COUNT1(J)+1

         RAD1(J)=RAD1(J)+SQRT(Y(N)*Y(N)+Z(N)*Z(N))

      END IF

      IF(X(N)>X2.AND.X(N)<XCAP2)THEN

         J=(XCAP2-X(N))/DELT+1

         COUNT2(J)=COUNT2(J)+1

         RAD2(J)=RAD2(J)+SQRT(Y(N)*Y(N)+Z(N)*Z(N))

      END IF

   END DO

   DO J=1,10

      XCAP1=XCAP1+DELT

      IF(RAD1(J)/COUNT1(J)>0.5D0*ORAD)THEN
         EXIT
      END IF

   END DO

   DO J=1,10

      XCAP2=XCAP2-DELT

      IF(RAD2(J)/COUNT2(J)>0.5D0*ORAD)THEN
         EXIT
      END IF

   END DO



   END SUBROUTINE

!=======================================================================

   SUBROUTINE SYNSETUP(NATOM,X,Y,Z,NSYN,SYNDIR, &
                    XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE)

   IMPLICIT NONE

   INTEGER,VALUE:: NATOM,NSYN
!   INTEGER SYNRATIO!,NS,NSEG
   INTEGER N,NA!,JS

   INTEGER,DIMENSION(:),ALLOCATABLE::SYNDIR

!   DOUBLE PRECISION,VALUE::XCAP1,XCAP2
   DOUBLE PRECISION CYLEN,R,DELT,XLIM1,XLIM2,R2(2),X0,y0,z0,radius,phi,dphi,pi

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE

!   SYNRATIO=NATOM/NSYN

!  SETUP LOCALIZATION PROBABILITY AS A FUNCTION OF RADIUS:

!   CYLEN=XCAP2-XCAP1

!   DELT=5.0D0!CYLEN/NSYN/2

!   NSEG=CYLEN/DELT

!   DELT=CYLEN/NSEG

   radius=sqrt(y(1)*y(1)+z(1)*z(1))

!  DIRECTION OF SYNTHESIS:

   DO N=1,NSYN

!      CALL RANDOM_NUMBER(R)


!      IF(R<0.5D0)THEN
         SYNDIR(N)=1
!      ELSE
!         SYNDIR(N)=-1
!      END IF

   END DO

   pi=4*atan(1.0d0)
   dphi=2*pi/nsyn


   DO N=1,NSYN

      phi=(n-1)*dphi

      y0=radius*cos(phi)

      z0=radius*sin(phi)

!      phi2=phi1+dphi

!      CALL RANDOM_NUMBER(R)

!      JS=NSEG*R+1



!            XLIM1=XCAP1+(JS-1)*DELT

!            XLIM2=XLIM1+DELT

      xlim1=-2.0d0

      xlim2=-xlim1



10    CALL RANDOM_NUMBER(R)


      NA=NATOM*R+1

      IF(X(NA)<XLIM1.OR.X(NA)>XLIM2)THEN
         GOTO 10
      END IF

      if(abs(y(na)-y0)>10.0) goto 10

      if(abs(z(na)-z0)>10.0) goto 10

      CALL RANDOM_NUMBER(R2)



      XGTASE(1,N)=X(NA)+0.5D0*SYNDIR(N)
      YGTASE(1,N)=Y(NA)+R2(1)-0.5D0
      ZGTASE(1,N)=Z(NA)+R2(2)-0.5D0

      CALL RANDOM_NUMBER(R2)



      XGTASE(2,N)=X(NA)-0.5D0*SYNDIR(N)
      YGTASE(2,N)=Y(NA)+R2(1)-0.5D0
      ZGTASE(2,N)=Z(NA)+R2(2)-0.5D0

      CALL RANDOM_NUMBER(R2)



      XTPASE(1,N)=X(NA)+1.0D0*SYNDIR(N)
      YTPASE(1,N)=Y(NA)+R2(1)-0.5D0
      ZTPASE(1,N)=Z(NA)+R2(2)-0.5D0

      CALL RANDOM_NUMBER(R2)



      XTPASE(2,N)=X(NA)
      YTPASE(2,N)=Y(NA)+R2(1)-0.5D0
      ZTPASE(2,N)=Z(NA)+R2(2)-0.5D0

      CALL RANDOM_NUMBER(R2)



      XTPASE(3,N)=X(NA)-1.0D0*SYNDIR(N)
      YTPASE(3,N)=Y(NA)+R2(1)-0.5D0
      ZTPASE(3,N)=Z(NA)+R2(2)-0.5D0

      CALL RANDOM_NUMBER(R2)



      XEDASE(N)=X(NA)
      YEDASE(N)=Y(NA)+R2(1)-0.5D0
      ZEDASE(N)=Z(NA)+R2(2)-0.5D0

   END DO

   PRINT*,'COMPLEXES DISTRIBUTED AT'

1  FORMAT(I2,3(2X,F7.1))

   DO N=1,NSYN
      WRITE(*,1)N,XEDASE(N),YEDASE(N),ZEDASE(N)
   END DO


   END SUBROUTINE

!==================================================================

   SUBROUTINE SETPARA(L_G,K_G,L_P,K_P,THETA_0,KTHETA,PRES,LSTART,LREP,LSWITCH,KSWITCH,MSWITCH, &
              KGTASE,KTPASE,LTPASE,DREACT,DREACT2,KSUR,KEDASE,LEDASE,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED, &
              KSIDE,LSIDE,KWALL,KLEAD,KPAIR,LPAIR,DELTA,INVDELTA,BETA,PI,INVMPBP, &
              wthick,pthick,kaxis,sigma,pconstrict,dpcons,invdelx,prmbond,oldtypmax,nsyn,growth)

   IMPLICIT NONE

   integer oldtypmax,nsyn,n
   DOUBLE PRECISION growth

   DOUBLE PRECISION L_G,K_G,L_P,K_P,THETA_0,KTHETA,PRES,LSTART,LREP,DREACT,DREACT2
   DOUBLE PRECISION KSWITCH,LSWITCH,MSWITCH,SLOPE,ESWITCH,KSUR,wthick,pthick,kaxis,sigma
   DOUBLE PRECISION delx,pconstrict0
   DOUBLE PRECISION pconstrict(1000)
   DOUBLE PRECISION dpcons,invdelx
   DOUBLE PRECISION KGTASE,KTPASE,LTPASE,KSIDE,LSIDE,PI,INVMPBP,prmbond
   DOUBLE PRECISION KEDASE,LEDASE,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED,KWALL,KLEAD,KPAIR,LPAIR,DELTA,INVDELTA,BETA
   CHARACTER CHARA*512

   PI=3.141592653589793239D0

!   twopi=2*pi

!-----------------------
!  FOR SURFACE CONSTRAINTS:
   KSUR=50.0D0

!  PARAMETERS FOR ESYN SUBROUTINE:

   KGTASE=2.0D0

   KTPASE=5.0D0
   LTPASE=2.0D0

   DREACT=2.0D0

   DREACT2=DREACT*DREACT

   KEDASE=20.0D0
   LEDASE=1.0D0

   KGTTP=10.0D0
   LGTTP13=1.0D0
   LGTTP2=0.5D0

   KGTED=2.0D0
   LGTED=2.0D0

   KSIDE=10.0D0
   LSIDE=0.5D0

   KWALL=20.0D0

   KLEAD=10.0D0

   KPAIR=5.0D0
   LPAIR=1.5D0

   DELTA=0.000001D0
   INVDELTA=1.0D0/DELTA

   BETA=0.00000000000001D0


   INVMPBP=0.25D0

!  wall thickness
   wthick=2.0d0

!  periplasm thickness
   pthick=0.5d0

!  force constant to constrain cell axis
   kaxis=0.1d0

!  set width of midcell in nm
   sigma=10.0d0

!---------------

   OPEN(1,FILE='paras.info')

   WRITE(*,*)
   PRINT*,'READING PARAMETERS'

15 READ(1,16)CHARA

   IF(CHARA(1:4)/='Bond')THEN
      GOTO 15
   END IF

   READ(1,*)CHARA,L_G,K_G

   READ(1,*)CHARA,L_P,K_P

17 READ(1,16)CHARA
   IF(CHARA(1:4)/='Angl')THEN
      GOTO 17
   END IF

   READ(1,*)CHARA,THETA_0,KTHETA

   THETA_0=THETA_0*PI/180.0D0

18 READ(1,16)CHARA
   IF(CHARA(1:4)/='Pres')THEN
      GOTO 18
   END IF

   READ(1,*)PRES

28 READ(1,16)CHARA
   IF(CHARA(1:4)/='Cons')THEN
      GOTO 28
   END IF

   read(1,*)pconstrict0
   read(1,*)dpcons


!  set grid to calculate constriction pressure

   delx=0.2d0

   invdelx=1.0d0/delx

   do n=1,1000

      pconstrict(n)=pconstrict0*exp(-n*n*delx*delx/sigma/sigma)*invdelta

   end do

38 READ(1,16)CHARA
   IF(CHARA(1:4)/='Remo')THEN
      GOTO 38
   END IF

   read(1,*)prmbond

   read(1,*)oldtypmax

48 READ(1,16)CHARA
   IF(CHARA(1:4)/='SYNC')THEN
      GOTO 48
   END IF

   read(1,*)nsyn
   read(1,*)growth

!   alpha0=alpha0*pi/180.0d0

!   sina0=sin(alpha0)

!   cosa0=cos(alpha0)

!   dalpha=10.0*pi/180.0


   CLOSE(1)


!----------------------------------------------------

16 format(A50)

20 FORMAT(5X,A5,1X,F7.4,1X,A2,19X,A5,E10.3,1X,A15)

   WRITE(*,20)'L_G =',L_G,'nm','K_G =',K_G,'* 1e-20 J/nm**2'
   WRITE(*,20)'L_P =',L_P,'nm','K_P =',K_P,'* 1e-20 J/nm**2'

21 FORMAT(1X,A9,1X,F7.4,1X,A3,14X,A9,E10.3,1X,A9)

   WRITE(*,21)'theta_0 =',THETA_0,'rad','Ktheta =',KTHETA,'* 1e-20 J'

22 FORMAT(A10,1X,F7.4,1X,A8)
   WRITE(*,22)'Pressure =',PRES,'* 10 MPa'  !,'Temperature =',TEM,'K'


   WRITE(*,*)


! -- ASSUME THERE IS NO FORCE FOR EXTENSION < CONTOUR LENGTH / 4. USE LREP INSTEAD OF LBOND:
   LSTART=1.0D0

   LREP=L_P-LSTART

! -- ALSO USE ALTERNATIVE BOND CONSTANT:
!      KREP=KBOND*LBOND/LREP

   LSWITCH=0.9D0*LREP

   ESWITCH=K_P*LSWITCH**2*(LREP/4.0D0/(LREP-LSWITCH)+0.5D0)

   SLOPE=K_P*(LREP/4.0D0/(1.0D0-LSWITCH/LREP)**2-LREP/4.0D0+LSWITCH)

   KSWITCH=SLOPE/2.0D0/LSWITCH
   MSWITCH=ESWITCH-KSWITCH*LSWITCH**2


! BEFORE SWITCHING, POTENTIAL OF PEPTIDE: E(XBOND) = KBOND*XBOND**2*(1/2+1/4(1-XBOND/LBOND))
! AFTER SWITCHING: E(XBOND)=KSWITCH*XBOND**2+MSWITCH



   END SUBROUTINE

!========================================================================


   SUBROUTINE VISUALPSF(npgold,NATOMOLD,NATOMNEW,NPEPMAX,nsyn,noldmax,pglen,pgid)

   IMPLICIT NONE

   INTEGER,VALUE:: NPGOLD,NATOMOLD,NATOMNEW,nsyn,NPEPMAX,noldmax
   INTEGER NTOTAL,length
   INTEGER IPG,N,J,I,JATOM,NLINE,NBOND
   INTEGER J1,J2,J3,J4,J5,J6,J7,J8

   INTEGER, ALLOCATABLE, DIMENSION(:,:),intent(in)::pgid
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::BOND
   INTEGER, ALLOCATABLE, DIMENSION(:),intent(in)::pglen

   CHARACTER(8)TEX,TYP,RES,SEGID,RESNO,ATOM
   DOUBLE PRECISION CHARG,MASS,W1,W2,GROWTH




   NTOTAL=natomold+2*NATOMNEW  +2*NPEPMAX+6*nsyn+ 4*nsyn + 12*nsyn+  2*nsyn    +4*nsyn+2*noldmax
!                         |                                 |          |          |            |          |            |
            !             |                                 |          |          |            |          |            |
!                    PG beads                          pep bonds     PBPs     GTASE-TIP  crlking    cleaving      EDASE-holding


   NBOND=    NATOMold+natomnew      +NPEPMAX+ 2*nsyn+ 6*nsyn  + nsyn   +2*nsyn  +noldmax
!               |                    |             |         |           |        |
!               |                    |             |         |           |        |
!         gly-bonds               pep-bonds   GTASE-TIP   crlking     cleaving  EDASE-holding

!
!  EACH PEPTIDE BOND IS VISUALIZED WITH 2 BEADS
!  EACH PBP IS VISUALIZED WITH A BEAD
!  A CROSSLINKING EVENT IS VISUALIZED WITH 2 BEADS
!  A CLEAVING EVENT IS VISUALIZED WITH 2 BEADS
!  EDASE HOLDING A PEPTIDE IS VISUALIZED WITH 2 BEADS


   OPEN(1,FILE='visual.psf')

20  FORMAT(I8,1X,A4,1X,A4,1X,A3,2X,A3,2X,A4,2X,F9.6,6X,F8.4,11X,I1)
21  FORMAT(I8,1X,A)

   ALLOCATE(BOND(2,NBOND))


   WRITE(1,21)NTOTAL,'!NATOM'

!-- FIRST, CAP UNITS:

   RES='GLY'

   ATOM='ATOM'
   CHARG=0.0
   MASS=1.0
   TYP='GLY'

   IPG=1
   WRITE(RESNO,'(I1)')IPG


   JATOM=0

   NBOND=0


!   segid='CAP'

!   do n=1,caplen

!      JATOM=JATOM+1

!      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

!      NBOND=NBOND+1

!      BOND(1,NBOND)=JATOM


!      if(n<caplen)then

!         BOND(2,NBOND)=JATOM+1

!      else

!         BOND(2,NBOND)=jatom-caplen+1

!      end if

!      if(capopen(bond(1,nbond))==1.and.capopen(bond(2,nbond))==1) nbond=nbond-1

!   end do

!   do n=1,caplen

!      JATOM=JATOM+1

!      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

!      NBOND=NBOND+1

!      BOND(1,NBOND)=JATOM


!      if(n<caplen)then

!         BOND(2,NBOND)=JATOM+1

!      else

!         BOND(2,NBOND)=jatom-caplen+1

!      end if

!      if(capopen(bond(1,nbond))==1.and.capopen(bond(2,nbond))==1) nbond=nbond-1

!   end do




   SEGID='OLD'

   DO N=1,Natomold

      WRITE(1,20)n,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   ENDDO

   DO N=1,NPGOLD

      length=pglen(n)

      do j=1,length-1

         nbond=nbond+1

         bond(1,nbond)=pgid(j,n)
         bond(2,nbond)=pgid(j+1,n)

       end do

   end do

   jatom=natomold



   SEGID='NEW'

   DO N=NATOMOLD+1,NATOMold+natomnew


         JATOM=JATOM+1

         WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

         JATOM=JATOM+1

         WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0



            NBOND=NBOND+1

            BOND(1,NBOND)=JATOM-1

            BOND(2,NBOND)=JATOM


   ENDDO




! - LIST OF PEPTIDE CROSS-LINKS:

   RES='PEP'
   SEGID='PEP'
   TYP='PEP'

   IPG=2
   WRITE(RESNO,'(I1)')IPG

   DO N=1,NPEPMAX

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO


!- LIST OF PBPS:

   RES='PBP'
   SEGID='GT1'
   TYP='PBP'

   IPG=3
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

   SEGID='GT2'

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

   SEGID='TP1'

   IPG=4
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

   SEGID='TP2'

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

   SEGID='TP3'

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

   SEGID='EDA'

   IPG=5
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

   END DO

!  CONNECTION BETWEEN GTASES AND STRAND TIPS:

   RES='TIP'
   SEGID='TIP'
   TYP='TIP'
   IPG=6
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO


!  LIST OF PBP CROSSLINKING:

   RES='CRL'
   SEGID='CRL'
   TYP='CRL'

   IPG=7
   WRITE(RESNO,'(I1)')IPG

   DO N=1,6*nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO

!  LIST OF DELETED PEP-BONDS:

   RES='DEL'
   SEGID='DEL'
   TYP='DEL'

   IPG=8
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO

!  EDASE HOLDING PEPTIDES:

   RES='HOD'
   SEGID='HOD'
   TYP='HOD'

   IPG=9
   WRITE(RESNO,'(I1)')IPG

   DO N=1,nsyn

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO

!  old bonds

   res='OBO'
   segid='OBO'
   typ='OBO'
   ipg=10
   WRITE(RESNO,'(I2)')IPG

   do n=1,noldmax

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      JATOM=JATOM+1

      WRITE(1,20)JATOM,SEGID,RESNO,RES,TYP,ATOM,CHARG,MASS,0

      NBOND=NBOND+1
      BOND(1,NBOND)=JATOM-1
      BOND(2,NBOND)=JATOM

   END DO



!-------------------------------------------
! -- WRITE LIST OF BONDS:

   NLINE=NBOND/4

   WRITE(1,*)
   WRITE(1,21)NBOND,'!NBOND: bonds'

   DO N=1,NLINE
      I=1+4*(N-1)
      J1=BOND(1,I); J2=BOND(2,I); J3=BOND(1,I+1); J4=BOND(2,I+1)
      J5=BOND(1,I+2); J6=BOND(2,I+2); J7=BOND(1,I+3); J8=BOND(2,I+3)

      WRITE(1,'(8I8)')J1,J2,J3,J4,J5,J6,J7,J8
   END DO

   IF(MOD(NBOND,4)==1)THEN
      J1=BOND(1,NBOND); J2=BOND(2,NBOND)
      WRITE(1,'(2I8)')J1,J2

   ELSE IF(MOD(NBOND,4)==2)THEN
      J1=BOND(1,NBOND-1); J2=BOND(2,NBOND-1)
      J3=BOND(1,NBOND); J4=BOND(2,NBOND)
      WRITE(1,'(4I8)')J1,J2,J3,J4

   ELSE IF(MOD(NBOND,4)==3)THEN
      J1=BOND(1,NBOND-2); J2=BOND(2,NBOND-2)
      J3=BOND(1,NBOND-1); J4=BOND(2,NBOND-1)
      J5=BOND(1,NBOND); J6=BOND(2,NBOND)
      WRITE(1,'(6I8)')J1,J2,J3,J4,J5,J6

   END IF

   CLOSE(1)

   END SUBROUTINE

!==============================================================================

   SUBROUTINE DCDHEADER(JUNIT,JFILE,NTOTAL)

   IMPLICIT NONE

   CHARACTER*4 COOR
   CHARACTER*80 STRING1,STRING2
   INTEGER IFIRST,NFRAME,NFREQ,ZEROS5(5),PEROFF,ZEROS7(7),TWO,TWENTYFOUR,NTOT,JUNIT,JFILE,NTOTAL
   INTEGER*8 JDELTA
   CHARACTER (LEN=64) FILEDCD
   CHARACTER ZERO*1,CHARID1*1,CHARID2*2,CHARID3*3

   WRITE(ZERO,'(I1)')0

   IF(JFILE<10)THEN
      WRITE(CHARID1,'(I1)')JFILE
      FILEDCD='visual'//ZERO//ZERO//CHARID1//'.dcd'
   ELSE IF(JFILE<100)THEN
      WRITE(CHARID2,'(I2)')JFILE
      FILEDCD='visual'//ZERO//CHARID2//'.dcd'
   ELSE IF(JFILE<1000)THEN
      WRITE(CHARID3,'(I3)')JFILE
      FILEDCD='visual'//CHARID3//'.dcd'
   ELSE
      PRINT*,'TOO MANY DCD FILES, STOP NOW'
      STOP
   END IF

   OPEN(JUNIT,FILE=FILEDCD,FORM='UNFORMATTED')

   COOR='CORD'; NFRAME=10000; IFIRST=0; NFREQ=1; NTOT=100
   ZEROS5=0; JDELTA=1; PEROFF=0; ZEROS7=0; TWENTYFOUR=24; TWO=2
   STRING1='HELLOOOOOOO'; STRING2='WHAT THE HELL!'

   WRITE(JUNIT)COOR,NFRAME,IFIRST,NFREQ,NTOT,ZEROS5,JDELTA,PEROFF,ZEROS7,TWENTYFOUR
   WRITE(JUNIT)TWO,STRING1,STRING2
   WRITE(JUNIT)NTOTAL

!   PRINT*,'OPEN DCD FILE #',JFILE

   END SUBROUTINE

!==============================================================================

   SUBROUTINE WRITEDCD(JUNIT,NTOTAL,NATOM,NATOMOLD,NATOMnew,NPG,DNOR,ATOR,PEPDIR,X,Y,Z,GLYTIP, &
              newPGID,newPGLEN,NBONDPEP,BONDPEP,BONTYP,NPEPMAX,NSYN,TPPEP,noldbond,oldbond,oldtyp, &
              EDHOLD,EDPEP,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,NFRAME)

   IMPLICIT NONE

   INTEGER,VALUE:: JUNIT,NTOTAL,NATOM,NATOMOLD,NATOMnew,NPG,NBONDPEP,NPEPMAX,NSYN,noldbond
   INTEGER NFRAME,JCYCLE,NA,JW,NR,JR,NS,NB,JSTART,JSTOP,CHECK,JS,JEXIT,NNEW,N0
   INTEGER, ALLOCATABLE, DIMENSION(:),intent(in)::newPGLEN,BONTYP,ATOR,DNOR,PEPDIR,oldtyp
   INTEGER, ALLOCATABLE, DIMENSION(:)::mark
   INTEGER, ALLOCATABLE, DIMENSION(:,:),intent(in)::newPGID,BONDPEP,EDHOLD,GLYTIP,EDPEP,oldbond
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:),intent(in)::TPPEP
   REAL, ALLOCATABLE, DIMENSION(:)::XW,YW,ZW
   REAL R(2),shift
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),intent(in)::X,Y,Z,XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),intent(in)::XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE

   NFRAME=NFRAME+1

   ALLOCATE(XW(NTOTAL),YW(NTOTAL),ZW(NTOTAL))

!---------------

!  old PG

   xw(1:natomold)=x(1:natomold)
   yw(1:natomold)=y(1:natomold)
   zw(1:natomold)=z(1:natomold)


   JW=natomold

   n0=jw



!  new PG

   IF(NPG==0)THEN

      N0=n0+2*NATOMnew

      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0

      JW=N0

!print*,'new',n0

      GOTO 10
   END IF

   DO NR=1,NPG


      DO JR=1,newPGLEN(NR)-1

         na=newPGID(JR,NR)

         IF(ATOR(na)==-1.AND.DNOR(na)==-1)THEN
            CYCLE
         END IF

         JW=JW+1
         XW(JW)=X(na)
         YW(JW)=Y(na)
         ZW(JW)=Z(na)

         na=newPGID(JR+1,NR)

         IF(ATOR(na)/=-1.OR.DNOR(na)/=-1)THEN

            JW=JW+1
            XW(JW)=X(na)
            YW(JW)=Y(na)
            ZW(JW)=Z(na)

         ELSE

            na=newPGID(JR,NR)

            JW=JW+1
            XW(JW)=X(na)
            YW(JW)=Y(na)
            ZW(JW)=Z(na)

            EXIT

         END IF

      END DO




   END DO

   N0=n0+2*natomnew


   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF

   JW=N0

!--------------------
!  PEPTIDE CROSSLINKS:

10 DO NB=1,NBONDPEP
      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      XW(JW+1:JW+2)=X(BONDPEP(1:2,NB))
      YW(JW+1:JW+2)=Y(BONDPEP(1:2,NB))
      ZW(JW+1:JW+2)=Z(BONDPEP(1:2,NB))
      JW=JW+2
   END DO


!stop
!--------------------

!  SHOW FREE HOOKS:

   allocate(mark(natom))

   mark=0

   if(noldbond>0)then

      do nb=1,noldbond

         if(oldtyp(nb)==0)cycle

         mark(oldbond(1:2,nb))=1

      end do

   end if

   CALL RANDOM_NUMBER(R)



   DO NA=1,NATOM

      IF(DNOR(NA)==0.OR.ATOR(NA)/=1.or.mark(na)==1)THEN
         CYCLE
      END IF

      JCYCLE=0

      JEXIT=0

      DO NS=1,NSYN

         IF(TPPEP(1,1,NS)==NA.OR.TPPEP(2,1,NS)==NA)THEN
            JCYCLE=1
            EXIT
         END IF

         IF(TPPEP(1,2,NS)==NA.OR.TPPEP(2,2,NS)==NA)THEN
            JCYCLE=1
            EXIT
         END IF

         IF(TPPEP(1,3,NS)==NA.OR.TPPEP(2,3,NS)==NA)THEN
            JCYCLE=1
            EXIT
         END IF

      END DO

      IF(JCYCLE==1)THEN
         CYCLE
      END IF

      DO NS=1,NSYN
         IF(EDHOLD(1,NS)==NA.OR.EDHOLD(2,NS)==NA)THEN
            JCYCLE=1
            EXIT
         END IF
      END DO

      IF(JCYCLE==1)THEN
         CYCLE
      END IF

      XW(JW+1)=X(NA)
      YW(JW+1)=Y(NA)
      ZW(JW+1)=Z(NA)


!      CALL RANDOM_NUMBER(R)

      if(pepdir(na)==1)then

         XW(JW+2)=X(NA)+0.5
         YW(JW+2)=Y(NA)+(R(1)-0.5)/10
         ZW(JW+2)=Z(NA)+(R(2)-0.5)/10

!      elseif(pepdir(na)==2)then

!         XW(JW+2)=X(NA)+(R(1)-0.5)/10

!         shift=0.5/sqrt(y(na)*y(na)+z(na)*z(na))

!         YW(JW+2)=Y(NA)+y(na)*shift
!         ZW(JW+2)=Z(NA)+z(na)*shift

!      elseif(pepdir(na)==3)then

      else

         XW(JW+2)=X(NA)-0.5
         YW(JW+2)=Y(NA)+(R(1)-0.5)/10
         ZW(JW+2)=Z(NA)+(R(2)-0.5)/10

!      else

!         XW(JW+2)=X(NA)+(R(1)-0.5)/10

!         shift=0.5/sqrt(y(na)*y(na)+z(na)*z(na))

!         YW(JW+2)=Y(NA)-y(na)*shift
!         ZW(JW+2)=Z(NA)-z(na)*shift

      end if

      JW=JW+2

   END DO


   N0=N0+2*NPEPMAX


!  FUTURE CROSSLINKS ARE KEPT AT THE CENTER:

   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF

   JW=N0

!---------------------
!  TRANSGLYCOSYLASES:

   XW(JW+1:JW+NSYN)=XGTASE(1,1:NSYN)
   YW(JW+1:JW+NSYN)=YGTASE(1,1:NSYN)
   ZW(JW+1:JW+NSYN)=ZGTASE(1,1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

   XW(JW+1:JW+NSYN)=XGTASE(2,1:NSYN)
   YW(JW+1:JW+NSYN)=YGTASE(2,1:NSYN)
   ZW(JW+1:JW+NSYN)=ZGTASE(2,1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

!---------------------
!  TRANSPEPTIDASES:

   XW(JW+1:JW+NSYN)=XTPASE(1,1:NSYN)
   YW(JW+1:JW+NSYN)=YTPASE(1,1:NSYN)
   ZW(JW+1:JW+NSYN)=ZTPASE(1,1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

   XW(JW+1:JW+NSYN)=XTPASE(2,1:NSYN)
   YW(JW+1:JW+NSYN)=YTPASE(2,1:NSYN)
   ZW(JW+1:JW+NSYN)=ZTPASE(2,1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

   XW(JW+1:JW+NSYN)=XTPASE(3,1:NSYN)
   YW(JW+1:JW+NSYN)=YTPASE(3,1:NSYN)
   ZW(JW+1:JW+NSYN)=ZTPASE(3,1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

!---------------------
!  ENDOPEPTIDASES:

   XW(JW+1:JW+NSYN)=XEDASE(1:NSYN)
   YW(JW+1:JW+NSYN)=YEDASE(1:NSYN)
   ZW(JW+1:JW+NSYN)=ZEDASE(1:NSYN)

!   XW(JW+NSYN:JW+nsyn)=XW(JW+NSYN)
!   YW(JW+NSYN:JW+nsyn)=YW(JW+NSYN)
!   ZW(JW+NSYN:JW+nsyn)=ZW(JW+NSYN)

   JW=JW+NSYN!MAX

   N0=N0+6*nsyn

!---------------------

!   XW(JW+1:NTOTAL)=0.0
!   YW(JW+1:NTOTAL)=0.0
!   ZW(JW+1:NTOTAL)=0.0

!---------------------

!  VISUALIZE THE CONNECTION BETWEEN GTASES AND STRAND TIPS:

   DO NS=1,NSYN

      DO JS=1,2

         IF(GLYTIP(JS,NS)>0)THEN

            NA=GLYTIP(JS,NS)
            JW=JW+1
            XW(JW)=X(NA)
            YW(JW)=Y(NA)
            ZW(JW)=Z(NA)

            JW=JW+1

            XW(JW)=XGTASE(JS,NS)
            YW(JW)=YGTASE(JS,NS)
            ZW(JW)=ZGTASE(JS,NS)

         END IF

      END DO

   END DO

   N0=N0+4*NSYN!MAX

   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF


   JW=N0

!-------------------------

!  CROSSLINKING DONORS:

   DO NS=1,NSYN

      IF(TPPEP(1,1,NS)>0)THEN

         XW(JW+1)=XTPASE(1,NS)
         YW(JW+1)=YTPASE(1,NS)
         ZW(JW+1)=ZTPASE(1,NS)

         XW(JW+2)=X(TPPEP(1,1,NS))
         YW(JW+2)=Y(TPPEP(1,1,NS))
         ZW(JW+2)=Z(TPPEP(1,1,NS))

         JW=JW+2

      END IF

      IF(TPPEP(1,2,NS)>0)THEN

         XW(JW+1)=XTPASE(2,NS)
         YW(JW+1)=YTPASE(2,NS)
         ZW(JW+1)=ZTPASE(2,NS)

         XW(JW+2)=X(TPPEP(1,2,NS))
         YW(JW+2)=Y(TPPEP(1,2,NS))
         ZW(JW+2)=Z(TPPEP(1,2,NS))

         JW=JW+2

      END IF

      IF(TPPEP(1,3,NS)>0)THEN

         XW(JW+1)=XTPASE(3,NS)
         YW(JW+1)=YTPASE(3,NS)
         ZW(JW+1)=ZTPASE(3,NS)

         XW(JW+2)=X(TPPEP(1,3,NS))
         YW(JW+2)=Y(TPPEP(1,3,NS))
         ZW(JW+2)=Z(TPPEP(1,3,NS))

         JW=JW+2

      END IF

   END DO

!---------------------
!  CROSSLINKING ACCEPTORS:

   DO NS=1,NSYN

      IF(TPPEP(2,1,NS)>0)THEN

         XW(JW+1)=XTPASE(1,NS)
         YW(JW+1)=YTPASE(1,NS)
         ZW(JW+1)=ZTPASE(1,NS)

         XW(JW+2)=X(TPPEP(2,1,NS))
         YW(JW+2)=Y(TPPEP(2,1,NS))
         ZW(JW+2)=Z(TPPEP(2,1,NS))

         JW=JW+2

      END IF

      IF(TPPEP(2,2,NS)>0)THEN

         XW(JW+1)=XTPASE(2,NS)
         YW(JW+1)=YTPASE(2,NS)
         ZW(JW+1)=ZTPASE(2,NS)

         XW(JW+2)=X(TPPEP(2,2,NS))
         YW(JW+2)=Y(TPPEP(2,2,NS))
         ZW(JW+2)=Z(TPPEP(2,2,NS))

         JW=JW+2

      END IF

      IF(TPPEP(2,3,NS)>0)THEN

         XW(JW+1)=XTPASE(3,NS)
         YW(JW+1)=YTPASE(3,NS)
         ZW(JW+1)=ZTPASE(3,NS)

         XW(JW+2)=X(TPPEP(2,3,NS))
         YW(JW+2)=Y(TPPEP(2,3,NS))
         ZW(JW+2)=Z(TPPEP(2,3,NS))

         JW=JW+2

      END IF

   END DO

   N0=N0+12*NSYN!MAX

   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF


   JW=N0

!   JW=NTOTAL-4*nsyn

!---------------------

!  CLEAVING PEPTIDE BONDS:

   DO NS=1,NSYN

      IF(EDHOLD(1,NS)>0)THEN
         XW(JW+1)=X(EDHOLD(1,NS))
         YW(JW+1)=Y(EDHOLD(1,NS))
         ZW(JW+1)=Z(EDHOLD(1,NS))

         XW(JW+2)=X(EDHOLD(2,NS))
         YW(JW+2)=Y(EDHOLD(2,NS))
         ZW(JW+2)=Z(EDHOLD(2,NS))

         JW=JW+2

      END IF
   END DO

   n0=n0+2*nsyn

   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF

   JW=N0

!   JW=NTOTAL-2*NSYN!MAX

!---------------------

!  EDASE IS HOLDING PEPTIDE:

   DO NS=1,NSYN

      IF(EDPEP(1,NS)>0.AND.EDPEP(1,NS)/=TPPEP(2,1,NS).AND.EDPEP(1,NS)/=TPPEP(2,3,NS))THEN

         XW(JW+1)=XEDASE(NS)
         YW(JW+1)=YEDASE(NS)
         ZW(JW+1)=ZEDASE(NS)

         XW(JW+2)=X(EDPEP(1,NS))
         YW(JW+2)=Y(EDPEP(1,NS))
         ZW(JW+2)=Z(EDPEP(1,NS))

         JW=JW+2

      END IF

      IF(EDPEP(2,NS)>0.AND.EDPEP(2,NS)/=TPPEP(2,1,NS).AND.EDPEP(2,NS)/=TPPEP(2,3,NS))THEN

         XW(JW+1)=XEDASE(NS)
         YW(JW+1)=YEDASE(NS)
         ZW(JW+1)=ZEDASE(NS)

         XW(JW+2)=X(EDPEP(2,NS))
         YW(JW+2)=Y(EDPEP(2,NS))
         ZW(JW+2)=Z(EDPEP(2,NS))

         JW=JW+2

      END IF

   END DO

   n0=n0+4*nsyn

   IF(JW<N0)THEN
      XW(JW+1:N0)=0.0
      YW(JW+1:N0)=0.0
      ZW(JW+1:N0)=0.0
   END IF

   JW=N0

!  old bonds

!print*,'oldbond',noldbond

   do nb=1,noldbond

      if(oldtyp(nb)==0) cycle

      xw(jw+1:jw+2)=x(oldbond(1:2,nb))
      yw(jw+1:jw+2)=y(oldbond(1:2,nb))
      zw(jw+1:jw+2)=z(oldbond(1:2,nb))

      jw=jw+2

   end do

   if(jw<ntotal)then
      xw(jw+1:ntotal)=0.0
      yw(jw+1:ntotal)=0.0
      zw(jw+1:ntotal)=0.0
   end if

!---------------------

   WRITE(JUNIT)XW(1:NTOTAL)
   WRITE(JUNIT)YW(1:NTOTAL)
   WRITE(JUNIT)ZW(1:NTOTAL)

   DEALLOCATE(XW,YW,ZW,mark)

   END SUBROUTINE

!==============================================================================

   SUBROUTINE SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XSYN,YSYN,ZSYN, &
                         NEIGLY,GLYNUM)

   IMPLICIT NONE

   INTEGER NBONDGLY,NBONDPEP,NS,N,N1,N2,LENGTH
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::BONDGLY,BONDPEP
   INTEGER, ALLOCATABLE, DIMENSION(:)::BONTYP,GLYNUM
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)::NEIGLY
   DOUBLE PRECISION DIST2,XBOND,YBOND,ZBOND
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z,XSYN,YSYN,ZSYN
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE

   LENGTH=0

   DO N=1,NBONDGLY

      IF(LENGTH==100)THEN
         PRINT*,'WARNING: COMPLEX',NS,'HAS MORE THAN 50 GLYCAN NEIGHBORS'
         EXIT
      END IF

      N1=BONDGLY(1,N)
      N2=BONDGLY(2,N)

      XBOND=(X(N1)+X(N2))/2
      YBOND=(Y(N1)+Y(N2))/2
      ZBOND=(Z(N1)+Z(N2))/2

      DIST2=(XSYN(NS)-XBOND)**2+(YSYN(NS)-YBOND)**2+(ZSYN(NS)-ZBOND)**2

      IF(DIST2<36.0D0)THEN
         LENGTH=LENGTH+1

!         NEIBOND(1,LENGTH,NS)=N1
!         NEIBOND(2,LENGTH,NS)=N2

         NEIGLY(1,LENGTH,NS)=N1
         NEIGLY(2,LENGTH,NS)=N2

      END IF
   END DO

   GLYNUM(NS)=LENGTH

!   DO N=1,NBONDPEP

!      IF(LENGTH==100)THEN
!         PRINT*,'WARNING: COMPLEX',NS,'HAS MORE THAN 100 NEIGHBORS'
!         EXIT
!      END IF

!      IF(BONTYP(N)==0)THEN
!         CYCLE
!      END IF

!      N1=BONDPEP(1,N)
!      N2=BONDPEP(2,N)

!      XBOND=(X(N1)+X(N2))/2
!      YBOND=(Y(N1)+Y(N2))/2
!      ZBOND=(Z(N1)+Z(N2))/2

!      DIST2=(XSYN(NS)-XBOND)**2+(YSYN(NS)-YBOND)**2+(ZSYN(NS)-ZBOND)**2

!      IF(DIST2<36.0D0)THEN
!         LENGTH=LENGTH+1

!         NEIBOND(1,LENGTH,NS)=N1
!         NEIBOND(2,LENGTH,NS)=N2

!      END IF
!   END DO

!   NEINUM(NS)=LENGTH

   END SUBROUTINE
!==============================================================================

   SUBROUTINE LOCALRAD(NATOM,X,Y,Z,XCEN,XSIZE,LORAD)

   IMPLICIT NONE
   INTEGER,VALUE::NATOM
   INTEGER NTEMP,N,CHECK,NCOUNT

   REAL(KIND=8),VALUE::XCEN,XSIZE

   REAL(KIND=8)::LORAD

   REAL(KIND=8)::X0,Y0,Z0,DY,DZ
   REAL(KIND=8)::F,FY,FZ,FYOLD,FZOLD,DEV2,SIGMA2,REP,BETA,INVBETA

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)::XTEMP,YTEMP,ZTEMP

   ALLOCATE(XTEMP(NATOM),YTEMP(NATOM),ZTEMP(NATOM))

   BETA=0.000000001D0
   INVBETA=1.0D0/BETA

   X0=XCEN


      NTEMP=0

      DO N=1,NATOM
         IF(X(N)>X0-XSIZE.AND.X(N)<X0+XSIZE)THEN
            NTEMP=NTEMP+1
            XTEMP(NTEMP)=X0
            YTEMP(NTEMP)=Y(N)
            ZTEMP(NTEMP)=Z(N)
         END IF
      END DO

      FY=0.0D0; FZ=0.0D0
      CHECK=0

      DO N=1,100000000
         FYOLD=FY; FZOLD=FZ

         CALL DEVIATION(X0,Y0,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
         DEV2=SIGMA2
!--------
         REP=Y0+BETA
         CALL DEVIATION(X0,REP,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FY=-(SIGMA2-DEV2)*INVBETA
!-------
         REP=Z0+BETA
         CALL DEVIATION(X0,Y0,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

         F=SQRT(FY**2+FZ**2)
         Y0=Y0+0.002D0*FY/F
         Z0=Z0+0.002D0*FZ/F

         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
            CHECK=CHECK+1
            IF(CHECK==10)THEN
               EXIT
            END IF
         ELSE
            CHECK=0
         END IF

      END DO


   DEALLOCATE(XTEMP,YTEMP,ZTEMP)

      LORAD=0.0D0

      NCOUNT=0


      DO N=1,NATOM

         IF(ABS(X(N)-XCEN)>XSIZE)THEN
            CYCLE
         END IF

         NCOUNT=NCOUNT+1

         DY=Y(N)-Y0
         DZ=Z(N)-Z0


         LORAD=LORAD+SQRT(DY*DY+DZ*DZ)

      END DO

      LORAD=LORAD/NCOUNT

   END SUBROUTINE

!==============================================================================

   SUBROUTINE FINDCENTER(NATOM,X,Y,Z,XCEN,YCEN,ZCEN,XSIZE)

   IMPLICIT NONE
   INTEGER,VALUE::NATOM
   INTEGER NTEMP,N,CHECK

   REAL(KIND=8),VALUE::XCEN,XSIZE

   REAL(KIND=8)::YCEN,ZCEN

   REAL(KIND=8)::X0,Y0,Z0
   REAL(KIND=8)::F,FY,FZ,FYOLD,FZOLD,DEV2,SIGMA2,REP,BETA,INVBETA

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)::XTEMP,YTEMP,ZTEMP

   ALLOCATE(XTEMP(NATOM),YTEMP(NATOM),ZTEMP(NATOM))


   BETA=0.000000001D0
   INVBETA=1.0D0/BETA

   X0=XCEN


      NTEMP=0

      DO N=1,NATOM
         IF(X(N)>X0-XSIZE.AND.X(N)<X0+XSIZE)THEN
            NTEMP=NTEMP+1
            XTEMP(NTEMP)=X0
            YTEMP(NTEMP)=Y(N)
            ZTEMP(NTEMP)=Z(N)
         END IF
      END DO

      FY=0.0D0; FZ=0.0D0
      CHECK=0

      DO N=1,100000000
         FYOLD=FY; FZOLD=FZ

         CALL DEVIATION(X0,Y0,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
         DEV2=SIGMA2
!--------
         REP=Y0+BETA
         CALL DEVIATION(X0,REP,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FY=-(SIGMA2-DEV2)*INVBETA
!-------
         REP=Z0+BETA
         CALL DEVIATION(X0,Y0,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

         F=SQRT(FY**2+FZ**2)
         Y0=Y0+0.002D0*FY/F
         Z0=Z0+0.002D0*FZ/F

         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
            CHECK=CHECK+1
            IF(CHECK==10)THEN
               EXIT
            END IF
         ELSE
            CHECK=0
         END IF

      END DO

      YCEN=Y0; ZCEN=Z0

   DEALLOCATE(XTEMP,YTEMP,ZTEMP)

   END SUBROUTINE

!==============================================================================

   SUBROUTINE FINDAXIS(NATOM,X,Y,Z,XCEN,YCEN,ZCEN,XSIZE,UX,UY,UZ)

   IMPLICIT NONE
   INTEGER,VALUE::NATOM
   INTEGER NTEMP,N,CHECK

   REAL(KIND=8),VALUE::XCEN,XSIZE

   REAL(KIND=8)::YCEN,ZCEN,UX,UY,UZ

   REAL(KIND=8)::XSIZE_BY_2,X0,Y0,Z0,X00,Y00,Z00
   REAL(KIND=8)::F,FY,FZ,FYOLD,FZOLD,DEV2,SIGMA2,REP,DIST,INVDIST,BETA,INVBETA

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)::XTEMP,YTEMP,ZTEMP

   ALLOCATE(XTEMP(NATOM),YTEMP(NATOM),ZTEMP(NATOM))

   BETA=0.000000001D0
   INVBETA=1.0D0/BETA

   X0=XCEN-XSIZE; Y0=0.0D0; Z0=0.0D0

   XSIZE_BY_2=0.5D0*XSIZE

      NTEMP=0

      DO N=1,NATOM
         IF(X(N)>X0-XSIZE_BY_2.AND.X(N)<X0+XSIZE_BY_2)THEN
            NTEMP=NTEMP+1
            XTEMP(NTEMP)=X0
            YTEMP(NTEMP)=Y(N)
            ZTEMP(NTEMP)=Z(N)
         END IF
      END DO

      FY=0.0D0; FZ=0.0D0
      CHECK=0

      DO N=1,100000000
         FYOLD=FY; FZOLD=FZ

         CALL DEVIATION(X0,Y0,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
         DEV2=SIGMA2
!--------
         REP=Y0+BETA
         CALL DEVIATION(X0,REP,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FY=-(SIGMA2-DEV2)*INVBETA
!-------
         REP=Z0+BETA
         CALL DEVIATION(X0,Y0,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

         F=SQRT(FY**2+FZ**2)
         Y0=Y0+0.002D0*FY/F
         Z0=Z0+0.002D0*FZ/F

         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
            CHECK=CHECK+1
            IF(CHECK==10)THEN
               EXIT
            END IF
         ELSE
            CHECK=0
         END IF

      END DO

      X00=XCEN+XSIZE; Y00=0.0D0; Z00=0.0D0

      NTEMP=0

      DO N=1,NATOM
         IF(X(N)>X00-XSIZE_BY_2.AND.X(N)<X00+XSIZE_BY_2)THEN
            NTEMP=NTEMP+1
            XTEMP(NTEMP)=X00
            YTEMP(NTEMP)=Y(N)
            ZTEMP(NTEMP)=Z(N)
         END IF
      END DO

      FY=0.0D0; FZ=0.0D0
      CHECK=0

      DO N=1,100000000
         FYOLD=FY; FZOLD=FZ

         CALL DEVIATION(X00,Y00,Z00,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
         DEV2=SIGMA2
!--------
         REP=Y00+BETA
         CALL DEVIATION(X00,REP,Z00,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FY=-(SIGMA2-DEV2)*INVBETA
!-------
         REP=Z00+BETA
         CALL DEVIATION(X00,Y00,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

         F=SQRT(FY**2+FZ**2)
         Y00=Y00+0.002D0*FY/F
         Z00=Z00+0.002D0*FZ/F

         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
            CHECK=CHECK+1
            IF(CHECK==10)THEN
               EXIT
            END IF
         ELSE
            CHECK=0
         END IF

      END DO

      UX=X00-X0; UY=Y00-Y0; UZ=Z00-Z0

      DIST=SQRT(UX**2+UY**2+UZ**2)

      INVDIST=1.0D0/DIST

      UX=UX*INVDIST
      UY=UY*INVDIST
      UZ=UZ*INVDIST

      YCEN=0.5D0*(Y0+Y00)
      ZCEN=0.5D0*(Z0+Z00)

   DEALLOCATE(XTEMP,YTEMP,ZTEMP)

   END SUBROUTINE

!==============================================================================

   SUBROUTINE RESETFRAME(NATOM,NSYN,SYNDIR,ORAD,X,Y,Z,XEDASE,YEDASE,ZEDASE, &
              XGTASE,XTPASE,P_CURV,XLEAD,YLEAD,ZLEAD,YCENBI1,ZCENBI1,YCENBI2,ZCENBI2, &
              YCENMON,ZCENMON,YCENEDO,ZCENEDO)

   IMPLICIT NONE

   INTEGER,VALUE::NATOM,NSYN
   INTEGER NS,NTEMP,N,NCOUNT
   INTEGER,ALLOCATABLE, DIMENSION(:),INTENT(IN)::SYNDIR
   DOUBLE PRECISION,VALUE::ORAD

   DOUBLE PRECISION X0,Y0,Z0,X00,Y00,Z00,DIST,INVDIST!,BETA,INVBETA
   DOUBLE PRECISION DX,DY,DZ,PROJ,LRAD,XCEN,YCEN,ZCEN,UX,UY,UZ
!   DOUBLE PRECISION XTEMP(2000),YTEMP(2000),ZTEMP(2000)
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z,XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::XGTASE,XTPASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::P_CURV,XLEAD,YLEAD,ZLEAD
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::YCENBI1,ZCENBI1,YCENBI2,ZCENBI2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::YCENMON,ZCENMON,YCENEDO,ZCENEDO



!   BETA=0.000000001D0
!   INVBETA=1.0D0/BETA

!   CYLEN=XCAP2-XCAP1
!   MIDCEL=XCAP1+CYLEN/2
!   L0=CYLEN/8

!OMP PARALLEL &

!OMP DEFAULT(NONE) &

!OMP PRIVATE(NS,NTEMP,N,CHECK,NCOUNT,XTEMP,YTEMP,ZTEMP,DX,DY,DZ,PROJ) &
!OMP PRIVATE(X0,Y0,Z0,X00,Y00,Z00,F,FY,FZ,FYOLD,FZOLD,DEV2,SIGMA2,REP,DIST,INVDIST) &
!OMP SHARED(NATOM,NSYN,BETA,INVBETA,X,Y,Z,XEDASE,YEDASE,ZEDASE,XCEN,YCEN,ZCEN,UX,UY,UZ,ORAD,P_CURV,L0,MIDCEL)

!OMP DO

   DO NS=1,NSYN

!     center for surface constraint for bifunctional 1:

      X0=XGTASE(1,NS)

      Y0=0.0D0; Z0=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN

            NTEMP=NTEMP+1

            Y0=Y0+Y(N)

            Z0=Z0+Z(N)

         END IF

      END DO


      YCENBI1(NS)=Y0/NTEMP

      ZCENBI1(NS)=Z0/NTEMP

!     center for surface constraint for bifunctional 2:

      X0=XGTASE(2,NS)

      Y0=0.0D0; Z0=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN

            NTEMP=NTEMP+1

            Y0=Y0+Y(N)

            Z0=Z0+Z(N)
      
         END IF

      END DO


      YCENBI2(NS)=Y0/NTEMP

      ZCENBI2(NS)=Z0/NTEMP

!     center for surface constraint for monofunctional:

      X0=XTPASE(3,NS)

      Y0=0.0D0; Z0=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN

            NTEMP=NTEMP+1

            Y0=Y0+Y(N)

            Z0=Z0+Z(N)

         END IF

      END DO


      YCENMON(NS)=Y0/NTEMP

      ZCENMON(NS)=Z0/NTEMP


!     center for surface constraint for endopeptidase:

      X0=XEDASE(NS)

      Y0=0.0D0; Z0=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN

            NTEMP=NTEMP+1

            Y0=Y0+Y(N)

            Z0=Z0+Z(N)

         END IF

      END DO


      YCENEDO(NS)=Y0/NTEMP

      ZCENEDO(NS)=Z0/NTEMP



!-------------------------------------------------

!      X0=XEDASE(NS)-5.0D0; Y0=0.0D0; Z0=0.0D0

!      NTEMP=0

!      DO N=1,NATOM
!         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN
!            NTEMP=NTEMP+1
!            XTEMP(NTEMP)=X0
!            YTEMP(NTEMP)=Y(N)
!            ZTEMP(NTEMP)=Z(N)
!         END IF
!      END DO

!      FY=0.0D0; FZ=0.0D0
!      CHECK=0

!      DO N=1,100000000
!         FYOLD=FY; FZOLD=FZ

!         CALL DEVIATION(X0,Y0,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
!         DEV2=SIGMA2
!--------
!         REP=Y0+BETA
!         CALL DEVIATION(X0,REP,Z0,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

!         FY=-(SIGMA2-DEV2)*INVBETA
!-------
!         REP=Z0+BETA
!         CALL DEVIATION(X0,Y0,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

!         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

!         F=SQRT(FY**2+FZ**2)
!         Y0=Y0+0.002D0*FY/F
!         Z0=Z0+0.002D0*FZ/F

!         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
!            CHECK=CHECK+1
!            IF(CHECK==10)THEN
!               EXIT
!            END IF
!         ELSE
!            CHECK=0
!         END IF

!      END DO

!-----------------------------------

!      X00=XEDASE(NS)+5.0D0; Y00=0.0D0; Z00=0.0D0

!      NTEMP=0

!      DO N=1,NATOM
!         IF(X(N)>X00-5.0D0.AND.X(N)<X00+5.0D0)THEN
!            NTEMP=NTEMP+1
!            XTEMP(NTEMP)=X00
!            YTEMP(NTEMP)=Y(N)
!            ZTEMP(NTEMP)=Z(N)
!         END IF
!      END DO

!      FY=0.0D0; FZ=0.0D0
!      CHECK=0

!      DO N=1,100000000
!         FYOLD=FY; FZOLD=FZ

!         CALL DEVIATION(X00,Y00,Z00,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)
!         DEV2=SIGMA2
!--------
!         REP=Y00+BETA
!         CALL DEVIATION(X00,REP,Z00,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

!         FY=-(SIGMA2-DEV2)*INVBETA
!-------
!         REP=Z00+BETA
!         CALL DEVIATION(X00,Y00,REP,XTEMP,YTEMP,ZTEMP,NTEMP,SIGMA2)

!         FZ=-(SIGMA2-DEV2)*INVBETA
!----------

!         F=SQRT(FY**2+FZ**2)
!         Y00=Y00+0.002D0*FY/F
!         Z00=Z00+0.002D0*FZ/F

!         IF(FY*FYOLD+FZ*FZOLD<0.0D0)THEN
!            CHECK=CHECK+1
!            IF(CHECK==10)THEN
!               EXIT
!            END IF
!         ELSE
!            CHECK=0
!         END IF

!      END DO



!---------------------------------

!  NEW METHOD TO CALCULATE LEADING DIRECTION:


      X0=XEDASE(NS)-10.0D0

      Y0=0.0D0; Z0=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X0-5.0D0.AND.X(N)<X0+5.0D0)THEN

            NTEMP=NTEMP+1

            Y0=Y0+Y(N)

            Z0=Z0+Z(N)

         END IF

      END DO


      Y0=Y0/NTEMP

      Z0=Z0/NTEMP


      X00=XEDASE(NS)+10.0D0

      Y00=0.0D0; Z00=0.0D0

      NTEMP=0

      DO N=1,NATOM

         IF(X(N)>X00-5.0D0.AND.X(N)<X00+5.0D0)THEN

            NTEMP=NTEMP+1

            Y00=Y00+Y(N)

            Z00=Z00+Z(N)

         END IF

      END DO


      Y00=Y00/NTEMP

      Z00=Z00/NTEMP


      UX=X00-X0; UY=Y00-Y0; UZ=Z00-Z0

      DIST=SQRT(UX**2+UY**2+UZ**2)

      INVDIST=1.0D0/DIST

      UX=UX*INVDIST
      UY=UY*INVDIST
      UZ=UZ*INVDIST

      XCEN=XEDASE(NS)

      YCEN=0.5D0*(Y0+Y00)
      ZCEN=0.5D0*(Z0+Z00)


      DX=0.0D0!XEDASE(NS)-XCEN(NS)
      DY=YEDASE(NS)-YCEN
      DZ=ZEDASE(NS)-ZCEN

      XLEAD(NS)=UY*DZ-UZ*DY
      YLEAD(NS)=UZ*DX-UX*DZ
      ZLEAD(NS)=UX*DY-UY*DX

      DIST=SQRT(XLEAD(NS)**2+YLEAD(NS)**2+ZLEAD(NS)**2)

      INVDIST=1.0D0/DIST

      XLEAD(NS)=XLEAD(NS)*SYNDIR(NS)*INVDIST
      YLEAD(NS)=YLEAD(NS)*SYNDIR(NS)*INVDIST
      ZLEAD(NS)=ZLEAD(NS)*SYNDIR(NS)*INVDIST


!---------------------------------

!      IF(A22==1)THEN
!         P_CURV(NS)=1.0D0
!         CYCLE
!      END IF

!     LOCAL RADIUS:

      LRAD=0.0D0
      NCOUNT=0

      DO N=1,NATOM

         IF(ABS(X(N)-XCEN)>5.0D0)THEN
            CYCLE
         END IF

         NCOUNT=NCOUNT+1

         DX=X(N)-XCEN
         DY=Y(N)-YCEN
         DZ=Z(N)-ZCEN

         PROJ=DX*UX+DY*UY+DZ*UZ

         LRAD=LRAD+SQRT(DX*DX+DY*DY+DZ*DZ-PROJ*PROJ)

      END DO

      LRAD=LRAD/NCOUNT

!     DEPENDENCE ON RADIUS:

      P_CURV(NS)=EXP(-100*(ORAD/LRAD-1.0D0)**2)

!      IF(ORAD<LRAD)THEN

!         P_CURV(NS)=EXP(-25*(ORAD/LRAD-1.0D0)**2)

         IF(P_CURV(NS)<0.001D0)THEN
            P_CURV(NS)=0.001D0
         END IF

!         MreBlen(NS)=MreBlen0*EXP(-50*(ORAD/LRAD-1.0D0)**2)

!         IF(MreBlen(NS)<0.05D0)then
!            MreBlen(NS)=0.05D0
!         END IF

!      ELSE

!         P_CURV(NS)=EXP(-100*(ORAD/LRAD-1.0D0)**2)

!         MreBlen(NS)=MreBlen0*EXP(-40*(ORAD/LRAD-1.0D0)**2)

!      END IF

!      IF(P_CURV(NS)<0.001D0)THEN
!         P_CURV(NS)=0.001D0
!      END IF

!      P_UCURV(NS)=1.0D0/P_CURV(NS)

!      P_UCURV(NS)=EXP(10*ABS(ORAD/LRAD-1.0D0))

!      P_UCURV(NS)=1.0

!     DEPENDENCE ON LONGITUDE:
!      P_CURV(NS)=P_CURV(NS)*EXP(-ABS(XEDASE(NS)-MIDCEL)/L0)

   END DO

!OMP END DO NOWAIT

!OMP END PARALLEL



   END SUBROUTINE

!========================================================================

   SUBROUTINE DEVIATION(CENX,CENY,CENZ,X,Y,Z,LENGTH,SIGMA2)

   IMPLICIT NONE

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::R
   DOUBLE PRECISION CENX,CENY,CENZ,SIGMA2,RMEAN
   DOUBLE PRECISION X(2000),Y(2000),Z(2000)
   INTEGER LENGTH,N

!======================  
   ALLOCATE(R(LENGTH))

   RMEAN=0.0D0

   DO N=1,LENGTH
      R(N)=SQRT((X(N)-CENX)**2+(Y(N)-CENY)**2+(Z(N)-CENZ)**2)
      RMEAN=RMEAN+R(N)
   END DO

   RMEAN=RMEAN/LENGTH

   SIGMA2=0.0D0
   DO N=1,LENGTH
      SIGMA2=SIGMA2+(R(N)-RMEAN)**2
   END DO

   SIGMA2=SIGMA2/LENGTH

   DEALLOCATE(R)

   END SUBROUTINE


!==============================================================================

!- CALCULATE LEADING DIRECTION:

   SUBROUTINE LEADIR(NSYN,SYNDIR,XEDASE,YEDASE,ZEDASE,XCEN,YCEN,ZCEN,UX,UY,UZ,XLEAD,YLEAD,ZLEAD)

   IMPLICIT NONE

   INTEGER NS,NSYN

   INTEGER, ALLOCATABLE, DIMENSION(:)::SYNDIR
!   INTEGER, ALLOCATABLE, DIMENSION(:,:)::PGID

   DOUBLE PRECISION DIST,INVDIST,DX,DY,DZ
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::XEDASE,YEDASE,ZEDASE,XCEN,YCEN,ZCEN,UX,UY,UZ,XLEAD,YLEAD,ZLEAD


!OMP PARALLEL &

!OMP DEFAULT(NONE) &

!OMP PRIVATE(NS,DX,DY,DZ,DIST,INVDIST) &
!OMP SHARED(NSYN,SYNDIR,XEDASE,YEDASE,ZEDASE,XCEN,YCEN,ZCEN,UX,UY,UZ,XLEAD,YLEAD,ZLEAD)

!OMP DO

   DO NS=1,NSYN

      DX=XEDASE(NS)-XCEN(NS)
      DY=YEDASE(NS)-YCEN(NS)
      DZ=ZEDASE(NS)-ZCEN(NS)

      XLEAD(NS)=UY(NS)*DZ-UZ(NS)*DY
      YLEAD(NS)=UZ(NS)*DX-UX(NS)*DZ
      ZLEAD(NS)=UX(NS)*DY-UY(NS)*DX

      DIST=SQRT(XLEAD(NS)**2+YLEAD(NS)**2+ZLEAD(NS)**2)

      INVDIST=1.0D0/DIST

      XLEAD(NS)=XLEAD(NS)*SYNDIR(NS)/DIST
      YLEAD(NS)=YLEAD(NS)*SYNDIR(NS)/DIST
      ZLEAD(NS)=ZLEAD(NS)*SYNDIR(NS)/DIST

   END DO

!OMP END DO NOWAIT

!OMP END PARALLEL
   END SUBROUTINE


!==============================================================================

   SUBROUTINE ACTIVATE(NS,SYNTHESIS,SYNLOOP,XGTASE,YGTASE,ZGTASE,X,Y,Z,NEIGLY,GLYNUM, &
                          NLOOP,LOOP,LOOPLEN,LOOPTYP,PSTART,jredist)

   IMPLICIT NONE

   INTEGER NLOOP,N,NS,CHECK,J,N1,N2,NL,ILOOP,JCHECK,GETLOOP,JEXIT,JG
   integer jredist
   INTEGER,DIMENSION(:,:,:),ALLOCATABLE::NEIGLY
   INTEGER,DIMENSION(:,:),ALLOCATABLE:: LOOP,SYNTHESIS
   INTEGER,DIMENSION(:),ALLOCATABLE:: LOOPLEN,LOOPTYP,SYNLOOP,GLYNUM
   DOUBLE PRECISION XLOOP,YLOOP,ZLOOP,DIST,PSTART,R
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::XGTASE,YGTASE,ZGTASE

   jredist=0

!  CHECK IF THE COMPLEX IS IN A POTENTIAL TRAP DUE TO STERIC HINDRANCE:

   XP1=1.5D0*XGTASE(1,NS)-0.5D0*XGTASE(2,NS)
   YP1=2*YGTASE(1,NS)
   ZP1=2*ZGTASE(1,NS)

   XP2=1.5D0*XGTASE(2,NS)-0.5D0*XGTASE(1,NS)
   YP2=2*YGTASE(2,NS)
   ZP2=2*ZGTASE(2,NS)

   XP3=0.5D0*(XP1+XP2)
   YP3=0.0D0; ZP3=0.0D0

   JEXIT=0

   DO JG=1,GLYNUM(NS)

      N1=NEIGLY(1,JG,NS)
      N2=NEIGLY(2,JG,NS)

      XL1=X(N1); YL1=Y(N1); ZL1=Z(N1)
      XL2=X(N2); YL2=Y(N2); ZL2=Z(N2)

      CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

      IF(CHECK==1)THEN
         JEXIT=1
         EXIT
      END IF

   END DO

   IF(JEXIT==1)THEN
      RETURN
   END IF

   IF(SYNTHESIS(1,NS)==1)THEN
      GOTO 10
   END IF

!  THE SECOND GTASE IS ALREADY ACTIVE:

   IF(SYNTHESIS(2,NS)==1)THEN

!     CHECK IF PROBABILITY OF ACTIVATION IS HIGH:

!      JRSTART(NS)=JRSTART(NS)+1

      CALL RANDOM_NUMBER(R)


      IF(PSTART<R)THEN
         RETURN
      END IF

      NL=SYNLOOP(NS)

!     CHECK IF THE FIRST GTASE IS IN THE LOOP:

      JCHECK=0

      XP3=X(LOOP(1,NL))
      YP3=Y(LOOP(1,NL))
      ZP3=Z(LOOP(1,NL))

      XL1=XGTASE(1,NS); YL1=2*YGTASE(1,NS); ZL1=2*ZGTASE(1,NS)

      XL2=XGTASE(1,NS); YL2=0.0D0; ZL2=0.0D0

      DO J=2,LOOPLEN(NL)-1

         N1=LOOP(J,NL)

         N2=LOOP(J+1,NL)

         IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
            CYCLE
         END IF

         XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

         XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==0)THEN
         RETURN
      END IF

      SYNTHESIS(1,NS)=1

      RETURN

   END IF

!---------------------------------------------------

!  CHECK IF PROBABILITY OF ACTIVATION FOR THE FIRST GTASE IS HIGH:

!   JRSTART(NS)=JRSTART(NS)+1

   CALL RANDOM_NUMBER(R)


   IF(PSTART<R)THEN
      GOTO 10
   END IF

! --  DETERMINE THE LOOP ENCIRCLING THE FIRST TRANSGLYCOSYLASE:

   DIST=5.0D0

   GETLOOP=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      XLOOP=SUM(X(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(XGTASE(1,NS)-XLOOP)>DIST)THEN
         CYCLE
      END IF

      YLOOP=SUM(Y(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(YGTASE(1,NS)-YLOOP)>DIST)THEN
         CYCLE
      END IF

      ZLOOP=SUM(Z(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(ZGTASE(1,NS)-ZLOOP)>DIST)THEN
         CYCLE
      END IF

!     NOW CHECK THIS LOOP:

      JCHECK=0

      XP3=X(LOOP(1,NL))
      YP3=Y(LOOP(1,NL))
      ZP3=Z(LOOP(1,NL))

      XL1=XGTASE(1,NS); YL1=2*YGTASE(1,NS); ZL1=2*ZGTASE(1,NS)

      XL2=XGTASE(1,NS); YL2=0.0D0; ZL2=0.0D0

      DO J=2,LOOPLEN(NL)-1

         N1=LOOP(J,NL)

         N2=LOOP(J+1,NL)

         IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
            CYCLE
         END IF

         XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

         XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==1)THEN
         ILOOP=NL
         GETLOOP=1
         EXIT
      ENDIF


   END DO

   IF(GETLOOP==0)THEN
      PRINT*,'WARNING: CANNOT IDENTIFY SYNLOOP 1'

      jredist=1
      RETURN
   END IF

!----

!  CHECK IF THERE IS ANOTHER COMPLEX IN THIS LOOP:

   IF(LOOPTYP(ILOOP)==2)THEN
      GOTO 10
   END IF

!  CHECK IF THE OTHER TRANSGLYCOSYLASE IS IN THE SAME LOOP:

   NL=ILOOP

   JCHECK=0

   XP3=X(LOOP(1,NL))
   YP3=Y(LOOP(1,NL))
   ZP3=Z(LOOP(1,NL))

   XL1=XGTASE(2,NS); YL1=2*YGTASE(2,NS); ZL1=2*ZGTASE(2,NS)

   XL2=XGTASE(2,NS); YL2=0.0D0; ZL2=0.0D0

   DO J=2,LOOPLEN(NL)-1

      N1=LOOP(J,NL)

      N2=LOOP(J+1,NL)

      IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
         CYCLE
      END IF

      XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

      XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

      CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

      IF(CHECK==1)THEN
         JCHECK=JCHECK+1
      END IF

   END DO

   IF(MOD(JCHECK,2)==1)THEN
      SYNTHESIS(1,NS)=1
      LOOPTYP(ILOOP)=2
      SYNLOOP(NS)=ILOOP
   END IF

!===========================

!  IF THE FIRST GTASE IS ALREADY ACTIVE:

10 IF(SYNTHESIS(1,NS)==1)THEN

!     CHECK IF PROBABILITY OF ACTIVATION IS HIGH:

!      JRSTART(NS)=JRSTART(NS)+1

      CALL RANDOM_NUMBER(R)


      IF(PSTART<R)THEN
         RETURN
      END IF

      NL=SYNLOOP(NS)

!     CHECK IF THE SECOND GTASE IS IN THE LOOP:

      JCHECK=0

      XP3=X(LOOP(1,NL))
      YP3=Y(LOOP(1,NL))
      ZP3=Z(LOOP(1,NL))

      XL1=XGTASE(2,NS); YL1=2*YGTASE(2,NS); ZL1=2*ZGTASE(2,NS)

      XL2=XGTASE(2,NS); YL2=0.0D0; ZL2=0.0D0

      DO J=2,LOOPLEN(NL)-1

         N1=LOOP(J,NL)

         N2=LOOP(J+1,NL)

         IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
            CYCLE
         END IF

         XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

         XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==0)THEN
         RETURN
      END IF

      SYNTHESIS(2,NS)=1

      RETURN

   END IF
!----------------------------------------------------

!  CHECK IF PROBABILITY OF ACTIVATION IS HIGH:

!   JRSTART(NS)=JRSTART(NS)+1

   CALL RANDOM_NUMBER(R)


   IF(PSTART<R)THEN
      RETURN
   END IF

! --  DETERMINE THE LOOP ENCIRCLING THE SECOND TRANSGLYCOSYLASE:

   DIST=5.0D0

   GETLOOP=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      XLOOP=SUM(X(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(XGTASE(2,NS)-XLOOP)>DIST)THEN
         CYCLE
      END IF

      YLOOP=SUM(Y(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(YGTASE(2,NS)-YLOOP)>DIST)THEN
         CYCLE
      END IF

      ZLOOP=SUM(Z(LOOP(1:LOOPLEN(NL),NL)))/LOOPLEN(NL)

      IF(ABS(ZGTASE(2,NS)-ZLOOP)>DIST)THEN
         CYCLE
      END IF

!     NOW CHECK THIS LOOP:

      JCHECK=0

      XP3=X(LOOP(1,NL))
      YP3=Y(LOOP(1,NL))
      ZP3=Z(LOOP(1,NL))

      XL1=XGTASE(2,NS); YL1=2*YGTASE(2,NS); ZL1=2*ZGTASE(2,NS)

      XL2=XGTASE(2,NS); YL2=0.0D0; ZL2=0.0D0

      DO J=2,LOOPLEN(NL)-1

         N1=LOOP(J,NL)

         N2=LOOP(J+1,NL)

         IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
            CYCLE
         END IF

         XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

         XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==1)THEN
         ILOOP=NL
         GETLOOP=1
         EXIT
      ENDIF

   END DO

   IF(GETLOOP==0)THEN
      PRINT*,'WARNING: CANNOT IDENTIFY SYNLOOP 2'
      jredist=1
      RETURN
   END IF

!-----------------------------------------------

!  CHECK IF THERE IS ANOTHER COMPLEX IN THIS LOOP:

   IF(LOOPTYP(ILOOP)==2)THEN
      RETURN
   END IF

!  CHECK IF THE OTHER TRANSGLYCOSYLASE IS IN THE SAME LOOP:

   NL=ILOOP

   JCHECK=0

   XP3=X(LOOP(1,NL))
   YP3=Y(LOOP(1,NL))
   ZP3=Z(LOOP(1,NL))

   XL1=XGTASE(1,NS); YL1=2*YGTASE(1,NS); ZL1=2*ZGTASE(1,NS)

   XL2=XGTASE(1,NS); YL2=0.0D0; ZL2=0.0D0

   DO J=2,LOOPLEN(NL)-1

      N1=LOOP(J,NL)

      N2=LOOP(J+1,NL)

      IF(LOOP(1,NL)==N1.OR.LOOP(1,NL)==N2)THEN
         CYCLE
      END IF

      XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

      XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

      CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

      IF(CHECK==1)THEN
         JCHECK=JCHECK+1
      END IF

   END DO

   IF(MOD(JCHECK,2)==1)THEN
      SYNTHESIS(2,NS)=1
      LOOPTYP(ILOOP)=2
      SYNLOOP(NS)=ILOOP

   END IF


   END SUBROUTINE

!====================================================================

   SUBROUTINE FIRSTPG(NS,NPG,newPGID,newPGLEN,NATOM,DNOR,ATOR,PEPDIR,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                         SYNPG,GTLOAD,SYNDIR,GLYTIP,GTATRANS)

   IMPLICIT NONE

   INTEGER NPG,NATOM,NS
   INTEGER, ALLOCATABLE, DIMENSION(:):: newPGLEN,DNOR,ATOR,PEPDIR,SYNDIR,GTLOAD,GTATRANS
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::newPGID,SYNPG,GLYTIP
   DOUBLE PRECISION R,dy,dz,invrad
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z!,XLEAD,YLEAD,ZLEAD
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE

!  first strand

   NATOM=NATOM+1

   DNOR(NATOM)=1
   ATOR(NATOM)=1


   CALL RANDOM_NUMBER(R)

   IF(R<0.5)THEN
      PEPDIR(NATOM)=-1
   ELSE
      PEPDIR(NATOM)=1
   END IF

!   pepdir(natom)=4*r+1


   NPG=NPG+1
   newPGID(1,NPG)=NATOM
   newPGLEN(NPG)=1
!   PGTYP(NPG)=3

   GLYTIP(1,NS)=NATOM

   SYNPG(1,NS)=NPG


   X(NATOM)=XGTASE(1,NS)

   invrad=1.0d0/sqrt(YGTASE(1,NS)*YGTASE(1,NS)+ZGTASE(1,NS)*ZGTASE(1,NS))
   dy=-syndir(ns)*ZGTASE(1,NS)*invrad
   dz=syndir(ns)*YGTASE(1,NS)*invrad

   Y(NATOM)=YGTASE(1,NS)-0.5D0*dy
   Z(NATOM)=ZGTASE(1,NS)-0.5D0*dz


!  second strand

   NATOM=NATOM+1

   DNOR(NATOM)=1
   ATOR(NATOM)=1

   pepdir(natom)=-pepdir(natom-1)

!   if(pepdir(natom-1)>2)then
!      pepdir(natom)=pepdir(natom-1)-2
!   else
!      pepdir(natom)=pepdir(natom-1)+2
!   end if

   NPG=NPG+1
   newPGID(1,NPG)=NATOM
   newPGLEN(NPG)=1
!   PGTYP(NPG)=3

   GLYTIP(2,NS)=NATOM

   SYNPG(2,NS)=NPG


   X(NATOM)=XGTASE(2,NS)

   Y(NATOM)=YGTASE(2,NS)-0.5D0*dy
   Z(NATOM)=ZGTASE(2,NS)-0.5D0*dz


   GTLOAD(NS)=0


   GTATRANS(NS)=0



   END SUBROUTINE

!====================================================================

   SUBROUTINE ELONGATE(NS,NPG,SYNPG,newPGID,newPGLEN,NBONDGLY,BONDGLY,X,Y,Z,XGTASE,YGTASE,ZGTASE,L_G,LGTASE, &
                          GLYTIP,GLYSEC,GTLOAD,NATOM,DNOR,ATOR,PEPDIR,GTATRANS,syndir, &
                            ntemlink,ntemlinkmax,oldglysec,oldglytip,temlink)

   IMPLICIT NONE

   INTEGER NS,NPG,NATOM,IPG1,IPG2,NA1,NA2,NBONDGLY,ntemlink,ntemlinkmax,n
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::newPGID,BONDGLY,SYNPG,GLYTIP,GLYSEC,oldglysec,oldglytip,temlink
   INTEGER, ALLOCATABLE, DIMENSION(:)::newPGLEN,DNOR,ATOR,PEPDIR,syndir,GTATRANS,GTLOAD
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE,LGTASE
   DOUBLE PRECISION DX1,DY1,DZ1,DIST1,DX2,DY2,DZ2,DIST2,L_G
   DOUBLE PRECISION R,dy,dz,invrad

   IPG1=SYNPG(1,NS)

   IPG2=SYNPG(2,NS)

   if(ipg1>0)then

      NA1=newPGID(newPGLEN(IPG1),IPG1)

      DX1=XGTASE(1,NS)-X(NA1)
      DY1=YGTASE(1,NS)-Y(NA1)
      DZ1=ZGTASE(1,NS)-Z(NA1)

      DIST1=SQRT(DX1**2+DY1**2+DZ1**2)

      IF(DIST1<=L_G+0.1D0) RETURN

   END IF

   if(ipg2>0)then

      NA2=newPGID(newPGLEN(IPG2),IPG2)

      DX2=XGTASE(2,NS)-X(NA2)
      DY2=YGTASE(2,NS)-Y(NA2)
      DZ2=ZGTASE(2,NS)-Z(NA2)

      DIST2=SQRT(DX2**2+DY2**2+DZ2**2)

      IF(DIST2<=L_G+0.1D0) RETURN

   END IF

!  first strand elongating


   if(ipg1>0)then

      NATOM=NATOM+1

      DNOR(NATOM)=1
      ATOR(NATOM)=1


      PEPDIR(NATOM)=-PEPDIR(NA1)!+syndir(ns)

!      if(PEPDIR(NATOM)>4) PEPDIR(NATOM)=PEPDIR(NATOM)-4

!      if(pepdir(natom)<1) pepdir(natom)=pepdir(natom)+4


      newPGLEN(IPG1)=newPGLEN(IPG1)+1

      newPGID(newPGLEN(IPG1),IPG1)=NATOM

      X(NATOM)=X(NA1)+DX1*L_G/DIST1
      Y(NATOM)=Y(NA1)+DY1*L_G/DIST1
      Z(NATOM)=Z(NA1)+DZ1*L_G/DIST1

      NBONDGLY=NBONDGLY+1

      BONDGLY(1,NBONDGLY)=NA1
      BONDGLY(2,NBONDGLY)=NATOM

      GLYTIP(1,NS)=NATOM
      GLYSEC(1,NS)=NA1


      if(newpglen(ipg1)==2.and.oldglytip(1,ns)>0)then

         if(ntemlink==ntemlinkmax)then

            do n=1,ntemlink-1

               temlink(1:4,n)=temlink(1:4,n+1)

            end do

            ntemlink=ntemlink-1

         end if

         ntemlink=ntemlink+1

         temlink(1,n)=oldglysec(1,ns)

         temlink(2,n)=oldglytip(1,ns)

         temlink(3,n)=newpgid(1,ipg1)

         temlink(4,n)=natom

      end if



   end if


!  second strand elongating


   if(ipg2>0)then

      NATOM=NATOM+1

      DNOR(NATOM)=1
      ATOR(NATOM)=1

      PEPDIR(NATOM)=-PEPDIR(NA2)!+syndir(ns)

!      if(PEPDIR(NATOM)>4) PEPDIR(NATOM)=PEPDIR(NATOM)-4

!      if(pepdir(natom)<1) pepdir(natom)=pepdir(natom)+4


      newPGLEN(IPG2)=newPGLEN(IPG2)+1

      newPGID(newPGLEN(IPG2),IPG2)=NATOM

      X(NATOM)=X(NA2)+DX2*L_G/DIST2
      Y(NATOM)=Y(NA2)+DY2*L_G/DIST2
      Z(NATOM)=Z(NA2)+DZ2*L_G/DIST2

      NBONDGLY=NBONDGLY+1

      BONDGLY(1,NBONDGLY)=NA2
      BONDGLY(2,NBONDGLY)=NATOM

      GLYTIP(2,NS)=NATOM
      GLYSEC(2,NS)=NA2

      if(newpglen(ipg2)==2.and.oldglytip(2,ns)>0)then

         if(ntemlink==ntemlinkmax)then

            do n=1,ntemlink-1

               temlink(1:4,n)=temlink(1:4,n+1)

            end do

            ntemlink=ntemlink-1

         end if

         ntemlink=ntemlink+1

         temlink(1,n)=oldglysec(2,ns)

         temlink(2,n)=oldglytip(2,ns)

         temlink(3,n)=newpgid(1,ipg2)

         temlink(4,n)=natom

      end if

   end if


!  first strand initiated

   if(ipg1==0)then

      NATOM=NATOM+1

      DNOR(NATOM)=1
      ATOR(NATOM)=1

      pepdir(natom)=-pepdir(natom-1)

!      if(pepdir(natom-1)>2)then
!         pepdir(natom)=pepdir(natom-1)-2
!      else
!         pepdir(natom)=pepdir(natom-1)+2
!      end if

      NPG=NPG+1
      newPGID(1,NPG)=NATOM
      newPGLEN(NPG)=1
!      PGTYP(NPG)=3

      GLYTIP(1,NS)=NATOM

      SYNPG(1,NS)=NPG

      X(NATOM)=XGTASE(1,NS)

      invrad=1.0d0/sqrt(YGTASE(1,NS)*YGTASE(1,NS)+ZGTASE(1,NS)*ZGTASE(1,NS))
      dy=-syndir(ns)*ZGTASE(1,NS)*invrad
      dz=syndir(ns)*YGTASE(1,NS)*invrad

      Y(NATOM)=YGTASE(1,NS)-0.5D0*dy
      Z(NATOM)=ZGTASE(1,NS)-0.5D0*dz

   end if

!  second strand initiated

   if(ipg2==0)then

      NATOM=NATOM+1

      DNOR(NATOM)=1
      ATOR(NATOM)=1

      pepdir(natom)=-pepdir(natom-1)

!      if(pepdir(natom-1)>2)then
!         pepdir(natom)=pepdir(natom-1)-2
!      else
!         pepdir(natom)=pepdir(natom-1)+2
!      end if

      NPG=NPG+1
      newPGID(1,NPG)=NATOM
      newPGLEN(NPG)=1
!      PGTYP(NPG)=3

      GLYTIP(2,NS)=NATOM

      SYNPG(2,NS)=NPG

      X(NATOM)=XGTASE(2,NS)

      invrad=1.0d0/sqrt(YGTASE(2,NS)*YGTASE(2,NS)+ZGTASE(2,NS)*ZGTASE(2,NS))
      dy=-syndir(ns)*ZGTASE(2,NS)*invrad
      dz=syndir(ns)*YGTASE(2,NS)*invrad

      Y(NATOM)=YGTASE(2,NS)-0.5D0*dy
      Z(NATOM)=ZGTASE(2,NS)-0.5D0*dz

   end if


   GTLOAD(NS)=0


   GTATRANS(NS)=0

   LGTASE(1:2,NS)=0.5D0

   END SUBROUTINE

!====================================================================


   SUBROUTINE TP1CRLK(NS,npgold,npg,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,newPGID,newPGLEN,BONDPEP,BONTYP,partner, &
               edpep,PEPDIR,syndir,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,JRANDS,RANDS,DREACT,JSKIP)

   IMPLICIT NONE

   INTEGER,value:: NS,npgold,npg!,nbondpep
   integer NA,IPG1,IPG2,ILOOP,JRANDS,JSKIP
   INTEGER JCYCLE,JR,N,J,NB,LOLEN,JL,NP,JCHECK,CHECK,N1,N2

   integer naback,ipg,jexit,npartner,jlen,jl1,jl2,j1

   integer npep

   INTEGER, ALLOCATABLE, DIMENSION(:,:,:):: TPPEP
   INTEGER, ALLOCATABLE, DIMENSION(:,:),intent(in):: pgid,newPGID,LOOP,BONDPEP,SYNPG,EDHOLD,partner,edpep
   INTEGER, ALLOCATABLE, DIMENSION(:):: pglen,newPGLEN,LOOPLEN,DNOR,ATOR,PEPDIR,syndir,SYNLOOP,BONTYP

   INTEGER, ALLOCATABLE, DIMENSION(:):: ncan

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: X,Y,Z,RANDS
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: XTPASE,YTPASE,ZTPASE

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::prob

   DOUBLE PRECISION DX,DY,DZ,DIST,psum,R,DREACT
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

!print*,'tp1'

   ILOOP=SYNLOOP(NS)

   LOLEN=LOOPLEN(ILOOP)

   IPG1=SYNPG(1,NS)

   IPG2=SYNPG(2,NS)

   NA=TPPEP(1,1,NS)


   npep=0

   allocate(prob(lolen),ncan(lolen))
   prob=0.0d0
   ncan=0

   DO J=1,LOLEN

      N=LOOP(J,ILOOP)

      IF(ATOR(N)/=1.or.pepdir(n)==pepdir(na)) cycle

      IF((X(NA)-X(N))*SYNDIR(ns)>1.0D0) cycle

      jcycle=0


      if(n==edpep(1,ns).or.n==edpep(2,ns)) goto 47


!-----------------------------------------------------------
!     MAKE SURE THE PEPTIDE POINTING INWARD:

      IF(DNOR(N)/=0)THEN

         XL1=X(N)-0.5D0*synDIR(ns)
         YL1=2*Y(N)
         ZL1=2*Z(N)

         XL2=XL1
         YL2=0.0D0; ZL2=0.0D0

         XP3=X(LOOP(1,ILOOP))
         YP3=Y(LOOP(1,ILOOP))
         ZP3=Z(LOOP(1,ILOOP))

         JCHECK=0

         DO JL=2,LOLEN-1

            N1=LOOP(JL,ILOOP)

            N2=LOOP(JL+1,ILOOP)

            IF(LOOP(1,ILOOP)==N1.OR.LOOP(1,ILOOP)==N2)THEN
               CYCLE
            END IF

            XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

            XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

            CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

            IF(CHECK==1)THEN
               JCHECK=JCHECK+1
            END IF

         END DO

         IF(MOD(JCHECK,2)==0)THEN
            CYCLE
         END IF

      ELSE ! SO THIS COULD BE TRIMERIC CROSSLINK FORMATION


         NB=EDHOLD(3,NS)

         if(nb==0) cycle

         IF(BONDPEP(1,NB)/=N)THEN
            CYCLE
         END IF

      END IF

!-----------------------------------------------------------
!     prevent crosslinking to the accompany strand

!      JCYCLE=0

      DO JR=1,newPGLEN(IPG1)
         IF(newPGID(JR,IPG1)==N)THEN
            JCYCLE=1
            EXIT
         END IF
      END DO

      IF(JCYCLE==1)THEN
         CYCLE
      END IF

      IF(IPG2>0)THEN

         DO JR=1,newPGLEN(IPG2)
            IF(newPGID(JR,IPG2)==N)THEN
               JCYCLE=1
               EXIT
            END IF
         END DO

         IF(JCYCLE==1)THEN
            CYCLE
         END IF

      END IF


!---------------------------------------------------

!     2017jan18A: to prevent crosslinking across the other strand if the current tip is on LOOP and the acceptor bead is crosslinked

      if(dnor(n)==0)then

         do jl=1,lolen

            if(na==loop(jl,iloop))then


               do jl1=1,lolen

                  if(n==loop(jl1,iloop))then

                     if(jl1==1)then
                        jl2=lolen
                     else
                        jl2=jl1-1
                     end if

                     npartner=loop(jl2,iloop)

                     if(npartner/=partner(1,n).and.npartner/=partner(2,n)) jcycle=1


!                     if(jget==1.and.js==1) jcycle=1

!                     if(jget==0.and.js==3) jcycle=1

                     exit

                  end if

               end do

               exit

            end if

         end do

      end if

      if(jcycle==1) cycle



!     to prevent forming two crosslinks in a row on the same two strands:

      do jr=1,newpglen(ipg1)

         if(jr>1.and.newpgid(jr,ipg1)==na)then

            naback=newpgid(jr-1,ipg1)


            npartner=partner(1,naback)

            if(npartner>0)then

               jexit=0

!              check old PG

               do ipg=1,npgold

                  do jlen=1,pglen(ipg)

                     if(pgid(jlen,ipg)==npartner)then

                        jexit=1

                        if(jlen<pglen(ipg))then
                           if(pgid(jlen+1,ipg)==n) jcycle=1 
                        end if

                        if(jlen>1)then
                           if(pgid(jlen-1,ipg)==n) jcycle=1  
                        end if


                        exit

                     end if

                  end do

                  if(jexit==1) exit

               end do

               if(jcycle==1) exit

!              check new PG

               do ipg=1,npg

                  do jlen=1,newpglen(ipg)

                     if(newpgid(jlen,ipg)==npartner)then

                        jexit=1

                        if(jlen<newpglen(ipg))then
                           if(newpgid(jlen+1,ipg)==n) jcycle=1
                        end if

                        if(jlen>1)then
                           if(newpgid(jlen-1,ipg)==n) jcycle=1
                        end if


                        exit

                     end if

                  end do

                  if(jexit==1) exit
                  
               end do


            end if

         end if

         if(jcycle==1) exit

      end do

      if(jcycle==1) cycle

!---------------------------------

47    DX=XTPASE(1,NS)-X(N)
      DY=YTPASE(1,NS)-Y(N)
      DZ=ZTPASE(1,NS)-Z(N)

      DIST=SQRT(DX**2+DY**2+DZ**2)

      IF(DIST>DREACT)THEN
         CYCLE
      END IF

      npep=npep+1

      PROB(npep)=(1.0-DIST/DREACT)**2*JSKIP

      ncan(npep)=n

!      JRANDS=JRANDS+1

!      R=RANDS(JRANDS)


!      IF(PROB>R)THEN

!         TPPEP(2,1,NS)=N

!         EXIT

!      END IF

   END DO

   if(npep==0) return

   psum=sum(prob(1:npep))

   JRANDS=JRANDS+1

   R=RANDS(JRANDS)

   if(psum<r) return

   prob(1:npep)=prob(1:npep)/psum

   JRANDS=JRANDS+1

   R=RANDS(JRANDS)

   psum=0.0d0

   do j=1,npep

      psum=psum+prob(j)

      if(psum>r)then

         TPPEP(2,1,NS)=ncan(j)

         deallocate(prob,ncan)

         exit

      end if

   end do

!print*,'tp1 end'


   END SUBROUTINE

!====================================================================

   SUBROUTINE TP3CRLK(NS,npgold,npg,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,newPGID,newPGLEN,BONDPEP,BONTYP,partner, &
               edpep,PEPDIR,syndir,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,JRANDS,RANDS,DREACT,JSKIP)

   IMPLICIT NONE

   INTEGER,value:: NS,npgold,npg!,nbondpep
   integer NA,IPG1,IPG2,ILOOP,JRANDS,JSKIP
   INTEGER JCYCLE,JR,N,J,NB,LOLEN,JL,NP,JCHECK,CHECK,N1,N2

   integer naback,ipg,jexit,npartner,jlen,jl1,jl2,j1

   integer npep

   INTEGER, ALLOCATABLE, DIMENSION(:,:,:):: TPPEP
   INTEGER, ALLOCATABLE, DIMENSION(:,:),intent(in):: pgid,newPGID,LOOP,BONDPEP,SYNPG,EDHOLD,partner,edpep
   INTEGER, ALLOCATABLE, DIMENSION(:):: pglen,newPGLEN,LOOPLEN,DNOR,ATOR,PEPDIR,syndir,SYNLOOP,BONTYP

   INTEGER, ALLOCATABLE, DIMENSION(:):: ncan

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: X,Y,Z,RANDS
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: XTPASE,YTPASE,ZTPASE

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::prob

   DOUBLE PRECISION DX,DY,DZ,DIST,psum,R,DREACT
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

!print*,'tp3'

   ILOOP=SYNLOOP(NS)

   LOLEN=LOOPLEN(ILOOP)

   IPG1=SYNPG(2,NS)

   IPG2=SYNPG(1,NS)

   NA=TPPEP(1,3,NS)

   npep=0

   allocate(prob(lolen),ncan(lolen))
   prob=0.0d0
   ncan=0

   DO J=1,LOLEN

      N=LOOP(J,ILOOP)

      IF(ATOR(N)/=1.or.pepdir(n)==pepdir(na)) cycle

      IF((X(N)-X(NA))*SYNDIR(ns)>1.0D0) cycle

      jcycle=0

      if(n==edpep(1,ns).or.n==edpep(2,ns)) goto 47

!-----------------------------------------------------------
!     MAKE SURE THE PEPTIDE POINTING INWARD:

      IF(DNOR(N)/=0)THEN

         XL1=X(N)+0.5D0*synDIR(ns)
         YL1=2*Y(N)
         ZL1=2*Z(N)

         XL2=XL1
         YL2=0.0D0; ZL2=0.0D0

         XP3=X(LOOP(1,ILOOP))
         YP3=Y(LOOP(1,ILOOP))
         ZP3=Z(LOOP(1,ILOOP))

         JCHECK=0

         DO JL=2,LOLEN-1

            N1=LOOP(JL,ILOOP)

            N2=LOOP(JL+1,ILOOP)

            IF(LOOP(1,ILOOP)==N1.OR.LOOP(1,ILOOP)==N2)THEN
               CYCLE
            END IF

            XP1=X(N1); YP1=Y(N1); ZP1=Z(N1)

            XP2=X(N2); YP2=Y(N2); ZP2=Z(N2)

            CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

            IF(CHECK==1)THEN
               JCHECK=JCHECK+1
            END IF

         END DO

         IF(MOD(JCHECK,2)==0)THEN
            CYCLE
         END IF

      ELSE ! SO THIS COULD BE TRIMERIC CROSSLINK FORMATION

         NB=EDHOLD(3,NS)

         if(nb==0)then
            cycle
         end if

         IF(BONDPEP(1,NB)/=N)THEN
            CYCLE
         END IF

      END IF

!-----------------------------------------------------------

!      JCYCLE=0

      DO JR=1,newPGLEN(IPG1)
         IF(newPGID(JR,IPG1)==N)THEN
            JCYCLE=1
            EXIT
         END IF
      END DO

      IF(JCYCLE==1)THEN
         CYCLE
      END IF

      IF(IPG2>0)THEN

         DO JR=1,newPGLEN(IPG2)
            IF(newPGID(JR,IPG2)==N)THEN
               JCYCLE=1
               EXIT
            END IF
         END DO

         IF(JCYCLE==1)THEN
            CYCLE
         END IF

      END IF

!---------------------------------------------------

!     2017jan18A: to prevent crosslinking across the other strand if the current tip is on LOOP and the acceptor bead is crosslinked

      if(dnor(n)==0)then

         do jl=1,lolen

            if(na==loop(jl,iloop))then


               do jl1=1,lolen

                  if(n==loop(jl1,iloop))then

                     if(jl1==1)then
                        jl2=lolen
                     else
                        jl2=jl1-1
                     end if

                     npartner=loop(jl2,iloop)

                     if(npartner==partner(1,n).or.npartner==partner(2,n)) jcycle=1


                     exit

                  end if

               end do

               exit

            end if

         end do

      end if

      if(jcycle==1) cycle


!     to prevent forming two crosslinks in a row on the same two strands:

      do jr=1,newpglen(ipg1)

         if(jr>1.and.newpgid(jr,ipg1)==na)then

            naback=newpgid(jr-1,ipg1)

            npartner=partner(1,naback)

            if(npartner>0)then

               jexit=0

!              check old PG

               do ipg=1,npgold

                  do jlen=1,pglen(ipg)

                     if(pgid(jlen,ipg)==npartner)then

                        jexit=1

                        if(jlen<pglen(ipg))then
                           if(pgid(jlen+1,ipg)==n) jcycle=1
                        end if

                        if(jlen>1)then
                           if(pgid(jlen-1,ipg)==n) jcycle=1
                        end if

                        exit

                     end if

                  end do

                  if(jexit==1) exit

               end do

               if(jcycle==1) exit

!              check new PG

               do ipg=1,npg

                  do jlen=1,newpglen(ipg)

                     if(newpgid(jlen,ipg)==npartner)then

                        jexit=1

                        if(jlen<newpglen(ipg))then
                           if(newpgid(jlen+1,ipg)==n) jcycle=1
                        end if

                        if(jlen>1)then
                           if(newpgid(jlen-1,ipg)==n) jcycle=1
                        end if

                        exit

                     end if

                  end do

                  if(jexit==1) exit

               end do

            end if

         end if

         if(jcycle==1) exit

      end do

      if(jcycle==1) cycle

!---------------------------------

47    DX=XTPASE(3,NS)-X(N)
      DY=YTPASE(3,NS)-Y(N)
      DZ=ZTPASE(3,NS)-Z(N)

      DIST=SQRT(DX**2+DY**2+DZ**2)

      IF(DIST>DREACT)THEN
         CYCLE
      END IF

      npep=npep+1

      PROB(npep)=(1.0-DIST/DREACT)**2*JSKIP

      ncan(npep)=n

!      PROB=(1.0-DIST/DREACT)**2*JSKIP

!      JRANDS=JRANDS+1

!      R=RANDS(JRANDS)


!      IF(PROB>R)THEN

!         TPPEP(2,3,NS)=N

!         EXIT

!      END IF

   END DO

   if(npep==0) return

   psum=sum(prob(1:npep))

   JRANDS=JRANDS+1

   R=RANDS(JRANDS)

   if(psum<r) return

   prob(1:npep)=prob(1:npep)/psum

   JRANDS=JRANDS+1

   R=RANDS(JRANDS)

   psum=0.0d0

   do j=1,npep

      psum=psum+prob(j)

      if(psum>r)then

         TPPEP(2,3,NS)=ncan(j)

         deallocate(prob,ncan)

         exit

      end if

   end do

!print*,'tp3 end'

   END SUBROUTINE



!====================================================================

   SUBROUTINE TP2CRLK(NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,glytip,syndir,PEPDIR,X,Y,Z,ATOR,JRANDS,RANDS,DREACT,JSKIP)

   IMPLICIT NONE
   INTEGER NS,IPG,NA,N,JRANDS,JSKIP
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)::TPPEP
   INTEGER, ALLOCATABLE, DIMENSION(:,:):: SYNPG,glytip
   INTEGER, ALLOCATABLE, DIMENSION(:):: PEPDIR,syndir
   INTEGER, ALLOCATABLE, DIMENSION(:):: ATOR
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: X,Y,Z,RANDS
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XTPASE,YTPASE,ZTPASE
   DOUBLE PRECISION DX,DY,DZ,DIST,PROB,R,DREACT

!print*,'tp2'

   IPG=SYNPG(2,NS)

   IF(IPG==0)THEN
      RETURN
   END IF

!print*,'tp2'


   NA=TPPEP(1,2,NS)

   N=glytip(2,ns)!PGID(PGLEN(IPG),IPG)


   IF(PEPDIR(N)==PEPDIR(NA).or.ATOR(N)==0)THEN
      RETURN
   END IF

!   IF(ATOR(N)==0)THEN
!      RETURN
!   END IF

   DX=XTPASE(2,NS)-X(N)
   DY=YTPASE(2,NS)-Y(N)
   DZ=ZTPASE(2,NS)-Z(N)

   DIST=SQRT(DX**2+DY**2+DZ**2)

   IF(DIST<DREACT)THEN

      PROB=(1.0D0-DIST/DREACT)**2*JSKIP

      JRANDS=JRANDS+1

      R=RANDS(JRANDS)


!      CALL RANDOM_NUMBER(R)

      IF(PROB>R)THEN
         TPPEP(2,2,NS)=N
      END IF

   END IF

!print*,'tp2 end'


   END SUBROUTINE

!====================================================================


   SUBROUTINE CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,xtpase,ytpase,ztpase,SIGCROSS,DELTA)

   IMPLICIT NONE

   INTEGER NS,JS,JS1,JS2,N1,N2,CHECK
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::SIGCROSS
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)::TPPEP

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),intent(in)::XGTASE,YGTASE,ZGTASE,xtpase,ytpase,ztpase
   DOUBLE PRECISION DX,DY,DZ,DIST2,DELTA
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   N1=TPPEP(1,JS,NS)
   N2=TPPEP(2,JS,NS)

   DX=X(N2)-X(N1)
   DY=Y(N2)-Y(N1)
   DZ=Z(N2)-Z(N1)

   DIST2=DX*DX+DY*DY+DZ*DZ

   IF(DIST2>16.0D0)THEN
      RETURN
   END IF


   IF(JS==1.OR.JS==3)THEN

!     check if the potential pep-bond is almost straight:

      dx=0.5d0*(x(n1)+x(n2))-xtpase(js,ns)
      dy=0.5d0*(y(n1)+y(n2))-ytpase(js,ns)
      dz=0.5d0*(z(n1)+z(n2))-ztpase(js,ns)

      dist2=dx*dx+dy*dy+dz*dz

      if(dist2>0.25d0) return

      JS1=(JS+1)/2
      JS2=(5-JS)/2

      XL1=XGTASE(JS1,NS)
      YL1=YGTASE(JS1,NS)
      ZL1=ZGTASE(JS1,NS)

      XL2=XGTASE(JS2,NS)
      YL2=YGTASE(JS2,NS)
      ZL2=ZGTASE(JS2,NS)

      XL1=XL1+(XL1-XL2)*DELTA
      YL1=YL1+(YL1-YL2)*DELTA
      ZL1=ZL1+(ZL1-ZL2)*DELTA


!      XP1=X(N1); YP1=2*Y(N1); ZP1=2*Z(N1)

!      XP2=X(N2); YP2=2*Y(N2); ZP2=2*Z(N2)

      xp1=x(n1); yp1=y(n1); zp1=z(n1)
      xp2=x(n2); yp2=y(n2); zp2=z(n2)

!      x0=0.5d0*(xp1+xp2)
!      y0=0.5d0*(yp1+yp2)
!      z0=0.5d0*(zp1+zp2)

!      xp1=xp1+xp1-x0
!      yp1=yp1+yp1-y0
!      zp1=zp1+zp1-z0

      xp2=xp2+xp2-xp1
      yp2=yp2+yp2-yp1
      zp2=zp2+zp2-zp1


      XP3=0.5d0*(XP1+XP2)
      YP3=0.0D0
      ZP3=0.0D0

      XP1=2*XP1-XP3
      yp1=yp1+yp1
      zp1=zp1+zp1

      XP2=2*XP2-XP3
      yp2=yp2+yp2
      zp2=zp2+zp2

      CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

      IF(CHECK==1)THEN
         RETURN
      END IF

   END IF

   SIGCROSS(JS,NS)=1

   END SUBROUTINE

!====================================================================

   SUBROUTINE SEARCHPATH(NATOM,ILOOP,LOOP,LOOPLEN,NBONDGLY,BONDGLY,NBONDPEP,BONDPEP,BONTYP,NSTART,NSTOP,IPATH,IPLEN)

   IMPLICIT NONE
   INTEGER NATOM,ILOOP,NBONDGLY,NBONDPEP,NSTART,NSTOP,IPLEN,IPATH(100)
   INTEGER ILEN,NB,N1,N2,NPATH,JSTOP,JSTEP,JMOVE,NP,NPATH0,LENGTH,NTIP,J,NA
   INTEGER, ALLOCATABLE, DIMENSION(:)::LOOPLEN,BONTYP,MARKLOOP,PARTLEN,CHECK,PATHLEN,PATHTYP
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::LOOP,BONDGLY,BONDPEP,PARTNER,PATH

   ALLOCATE(MARKLOOP(NATOM))

   MARKLOOP=0

   MARKLOOP(LOOP(1:LOOPLEN(ILOOP),ILOOP))=1

   IF(MARKLOOP(NSTART)==1)THEN
      IPATH(1)=NSTART
      IPLEN=1
      DEALLOCATE(MARKLOOP)
      RETURN
   END IF

   IPLEN=0

   ALLOCATE(PARTNER(4,NATOM),PARTLEN(NATOM))
   PARTNER=0
   PARTLEN=0

   DO NB=1,NBONDGLY

      N1=BONDGLY(1,NB)
      N2=BONDGLY(2,NB)

      PARTLEN(N1)=PARTLEN(N1)+1
      PARTNER(PARTLEN(N1),N1)=N2

      PARTLEN(N2)=PARTLEN(N2)+1
      PARTNER(PARTLEN(N2),N2)=N1

   END DO

   DO NB=1,NBONDPEP

      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      N1=BONDPEP(1,NB)
      N2=BONDPEP(2,NB)

      PARTLEN(N1)=PARTLEN(N1)+1
      PARTNER(PARTLEN(N1),N1)=N2

      PARTLEN(N2)=PARTLEN(N2)+1
      PARTNER(PARTLEN(N2),N2)=N1

   END DO

!---------

   ALLOCATE(PATH(100,100),PATHLEN(100),PATHTYP(100))

   ALLOCATE(CHECK(NATOM))
   CHECK=0

   NPATH=1

   PATH(1,1)=NSTART ! THIS IS THE START

   PATHLEN(1)=1

   PATHTYP(1)=1

   CHECK(NSTART)=1
   CHECK(NSTOP)=1

   JSTOP=0

   DO JSTEP=1,1000

      JMOVE=0

      NPATH0=NPATH

      DO NP=1,NPATH0

         IF(PATHTYP(NP)==0)THEN
            CYCLE
         END IF

         LENGTH=PATHLEN(NP)

         NTIP=PATH(PATHLEN(NP),NP)

         IF(PARTLEN(NTIP)>0)THEN

            DO J=1,PARTLEN(NTIP)

               NA=PARTNER(J,NTIP)

               IF(CHECK(NA)==1)THEN
                  CYCLE
               END IF

               CHECK(NA)=1

               IF(LENGTH==PATHLEN(NP))THEN

                  PATHLEN(NP)=LENGTH+1

                  PATH(LENGTH+1,NP)=NA

                  IF(MARKLOOP(NA)==1)THEN
                     JSTOP=1
                     PATHTYP(NP)=2
                     JMOVE=1
                     EXIT
                  END IF

               ELSE

                  NPATH=NPATH+1

                  PATH(1:LENGTH,NPATH)=PATH(1:LENGTH,NP)

                  PATH(LENGTH+1,NPATH)=NA

                  PATHTYP(NPATH)=1

                  PATHLEN(NPATH)=LENGTH+1

                  IF(MARKLOOP(NA)==1)THEN
                     JSTOP=1
                     PATHTYP(NPATH)=2
                     JMOVE=1
                     EXIT
                  END IF

               END IF

               IF(JSTOP==1)THEN
                  EXIT
               END IF

               JMOVE=1

            END DO

         END IF

         IF(LENGTH==PATHLEN(NP))THEN
            PATHTYP(NP)=0
         END IF

      END DO

      IF(JMOVE==0)THEN
!         PRINT*,'DEAD END AT JSTEP=',JSTEP
         EXIT
      END IF

      IF(JSTOP==1)THEN
!         PRINT*,'FOUND PATH AT STEP',JSTEP
!         PRINT*,'NPATH',NPATH
         EXIT
      END IF

   END DO

   DO NP=1,NPATH

      IF(PATHTYP(NP)==2)THEN
!         PRINT*,'PATH',PATH(1:PATHLEN(NP),NP)
         IPLEN=PATHLEN(NP)
         IPATH(1:IPLEN)=PATH(1:IPLEN,NP)
      END IF

   END DO

   DEALLOCATE(MARKLOOP,PARTNER,PARTLEN,CHECK,PATH,PATHLEN,PATHTYP)

   END SUBROUTINE

!====================================================================

   SUBROUTINE POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                       NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,EDPEP,JDEACT,PARTNER, &
                       X,Y,Z,XGTASE,YGTASE,ZGTASE,jrestart)


   IMPLICIT NONE

   integer jrestart

   INTEGER JS,NS,NATOM,NLOOP,NLOOPDEL,NBONDGLY,NBONDPEP,NBONDDEL,NB,ILOOP,NBSIG
   INTEGER NA1,NA2,NSTART1,NSTART2,NSTOP1,NSTOP2,LEN1,LEN2,J1,J2,JL,J,JDEL,LOLEN
   INTEGER PATH1(100),PATH2(100)

   INTEGER, ALLOCATABLE, DIMENSION(:)::SYNLOOP,DNOR,ATOR,LOOPLEN,LOOPTYP,BONTYP,JDEACT
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::CRLKAGE,SIGCROSS,LOOP,BONDGLY,BONDPEP,EDPEP,PARTNER
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)::TPPEP

   INTEGER JCHECK,CHECK,n1,n2,js1,jcount

   DOUBLE PRECISION X0,y0,z0,r(3),radius

   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   DOUBLE PRECISION,ALLOCATABLE,INTENT(IN), DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE



   NA1=TPPEP(1,JS,NS)
   NA2=TPPEP(2,JS,NS)

   N1=NA1
   N2=NA2

   TPPEP(1:2,JS,NS)=0

   NBONDPEP=NBONDPEP+1

!print*,'add new bond',nbondpep

   BONTYP(NBONDPEP)=1

   BONDPEP(1,NBONDPEP)=NA1

   BONDPEP(2,NBONDPEP)=NA2

   DNOR(NA1)=0

   ATOR(NA2)=0

   IF(PARTNER(1,NA1)==0)THEN
      PARTNER(1,NA1)=NA2
   ELSE
      PARTNER(2,NA1)=NA2
   END IF

   IF(PARTNER(1,NA2)==0)THEN
      PARTNER(1,NA2)=NA1
   ELSE
      PARTNER(2,NA2)=NA1
   END IF


   NBSIG=0

   IF(DNOR(NA2)==0)THEN

      DO NB=1,NBONDPEP-1

         IF(BONTYP(NB)==0)THEN
            CYCLE
         END IF

         IF(BONDPEP(1,NB)==NA2)THEN
            BONTYP(NB)=2
            NBSIG=NB
            EXIT
         END IF

      END DO

   END IF

   IF(JS==1)THEN

      CRLKAGE(1,NS)=CRLKAGE(1,NS)+1

   ELSEIF(JS==2)THEN

      CRLKAGE(1,NS)=CRLKAGE(1,NS)+1
      CRLKAGE(2,NS)=CRLKAGE(2,NS)+1


   ELSE

      CRLKAGE(2,NS)=CRLKAGE(2,NS)+1

   END IF

   SIGCROSS(JS,NS)=0

   IF(EDPEP(1,NS)==NA2)THEN
      EDPEP(1,NS)=0
   END IF

   IF(EDPEP(2,NS)==NA2)THEN
      EDPEP(2,NS)=0
   END IF

!---------------------------------------

   ILOOP=SYNLOOP(NS)

   IF(JS==1)THEN

      NSTART1=NA2

      NSTOP1=NA1

      NSTART2=NA1

      NSTOP2=NA2

   ELSE

      NSTART1=NA1

      NSTOP1=NA2

      NSTART2=NA2

      NSTOP2=NA1

   END IF

   CALL SEARCHPATH(NATOM,ILOOP,LOOP,LOOPLEN,NBONDGLY,BONDGLY,NBONDPEP,BONDPEP,BONTYP,NSTART1,NSTOP1,PATH1,LEN1)


   IF(LEN1==0)THEN
      RETURN
   END IF

   CALL SEARCHPATH(NATOM,ILOOP,LOOP,LOOPLEN,NBONDGLY,BONDGLY,NBONDPEP,BONDPEP,BONTYP,NSTART2,NSTOP2,PATH2,LEN2)

   IF(LEN2==0)THEN

      RETURN
   END IF

!  NA1 AND NA2 ARE IN DIFFERENT ROLES:

   NA1=PATH1(LEN1)

   NA2=PATH2(LEN2)

   J1=0
   J2=0

   DO JL=1,LOOPLEN(ILOOP)

      IF(LOOP(JL,ILOOP)==NA1)THEN
         J1=JL
      ELSEIF(LOOP(JL,ILOOP)==NA2)THEN
         J2=JL
      END IF

   END DO

   IF(J1==0.OR.J2==0)THEN
      PRINT*,'CHECK ERROR'
      PRINT*,'PATH 1',LEN1,PATH1(1:LEN1)
      PRINT*,'PATH 2',LEN2,PATH2(1:LEN2)
      PRINT*,'LOOP',LOOP(1:LOOPLEN(ILOOP),ILOOP)

      IF(BONTYP(NBONDPEP)==0)THEN
         PRINT*,'ERROR: BOND ALREADY DELETED',NBONDPEP,BONDPEP(1:2,NBONDPEP)
         STOP
      END IF

      BONTYP(NBONDPEP)=0
      NBONDDEL=NBONDDEL+1

      NA1=BONDPEP(1,NBONDPEP)
      NA2=BONDPEP(2,NBONDPEP)

      IF(PARTNER(1,NA1)==NA2)THEN
         PARTNER(1,NA1)=PARTNER(2,NA1)
      END IF

      PARTNER(2,NA1)=0


      IF(PARTNER(1,NA2)==NA1)THEN
         PARTNER(1,NA2)=PARTNER(2,NA2)
      END IF

      PARTNER(2,NA2)=0



!print*,'nbonddel=',nbonddel

      DNOR(NA1)=-1
      ATOR(NA2)=1

      JDEACT(NS)=1

      IF(NBSIG/=0)THEN
         BONTYP(NBSIG)=1
      END IF

      RETURN

   END IF


!  THE FIRST LOOP:

   NLOOP=NLOOP+1

   LOOPTYP(NLOOP)=1

   LOLEN=LEN1

   LOOP(1:LOLEN,NLOOP)=PATH1(1:LEN1)

   JDEL=J2-J1

   IF(JDEL<0)THEN
      JDEL=JDEL+LOOPLEN(ILOOP)
   END IF

   DO J=1,JDEL-1

      JL=J+J1

      IF(JL>LOOPLEN(ILOOP))THEN
         JL=JL-LOOPLEN(ILOOP)
      END IF

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=LOOP(JL,ILOOP)

   END DO

   DO J=LEN2,1,-1

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=PATH2(J)

   END DO

   LOOPLEN(NLOOP)=LOLEN


!  THE SECOND LOOP, ALSO THE NEW SYNLOOP:

   NLOOP=NLOOP+1

   LOOPTYP(NLOOP)=2

   SYNLOOP(NS)=NLOOP

   LOLEN=LEN2

   LOOP(1:LOLEN,NLOOP)=PATH2(1:LEN2)

   JDEL=J1-J2

   IF(JDEL<0)THEN
      JDEL=JDEL+LOOPLEN(ILOOP)
   END IF

   DO J=1,JDEL-1

      JL=J+J2

      IF(JL>LOOPLEN(ILOOP))THEN
         JL=JL-LOOPLEN(ILOOP)
      END IF

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=LOOP(JL,ILOOP)

   END DO

   DO J=LEN1,1,-1

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=PATH1(J)

   END DO

   LOOPLEN(NLOOP)=LOLEN

!  DELETE ILOOP:

   IF(LOOPTYP(ILOOP)==0)THEN
      PRINT*,'ERROR: LOOP ALREADY DELETED'
      PRINT*,LOOP(1:LOOPLEN(ILOOP),ILOOP)
      STOP
   END IF

   NLOOPDEL=NLOOPDEL+1

   LOOPTYP(ILOOP)=0


   if(js==2) return


   return

!  check if Gtases are out of SYNLOOP

      ILOOP=SYNLOOP(NS)

      XP3=X(LOOP(1,ILOOP))
      YP3=Y(LOOP(1,ILOOP))
      ZP3=Z(LOOP(1,ILOOP))

      if(js==1)then
         js1=1
      else
         js1=2
      end if

      x0=XGTASE(js1,NS)
      y0=yGTASE(js1,NS)
      z0=zGTASE(js1,NS)

      jcount=0

      radius=0.5d0

13    XL1=XGTASE(js1,NS)
      YL1=2*YGTASE(js1,NS)
      ZL1=2*ZGTASE(js1,NS)

      XL2=XL1; YL2=0.0D0; ZL2=0.0D0

      JCHECK=0

      DO J=2,LOOPLEN(ILOOP)-1

         NA1=LOOP(J,ILOOP)

         NA2=LOOP(J+1,ILOOP)

         IF(LOOP(1,ILOOP)==NA1.OR.LOOP(1,ILOOP)==NA2)THEN
            CYCLE
         END IF

         XP1=X(NA1); YP1=Y(NA1); ZP1=Z(NA1)

         XP2=X(NA2); YP2=Y(NA2); ZP2=Z(NA2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==0)THEN

         call random_number(r(3))


         xgtase(js1,ns)=x0+radius*(r(1)-0.5d0)
         ygtase(js1,ns)=y0+radius*(r(2)-0.5d0)
         zgtase(js1,ns)=z0+radius*(r(3)-0.5d0)

         jcount=jcount+1

         if(jcount>1000)then

!            print*,'increase radius to adjust position of gtase',js1,ns

            print*,'cannot move Gtase inside SYNLOOP',js1,ns

!            jprint=2

            xgtase(js1,ns)=x0
            ygtase(js1,ns)=y0
            zgtase(js1,ns)=z0

            jrestart=1

            return
!            jcount=0

!            radius=radius+0.1d0

         end if

         goto 13

      end if


   END SUBROUTINE


!====================================================================


   SUBROUTINE SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

   IMPLICIT NONE
   INTEGER NSIG,NBONDPEP,NB
   INTEGER, ALLOCATABLE, DIMENSION(:)::BONTYP
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::BONDPEP,SIGBOND

   NSIG=0

   DO NB=1,NBONDPEP

      IF(BONTYP(NB)==2)THEN
         NSIG=NSIG+1
         SIGBOND(1:2,NSIG)=BONDPEP(1:2,NB)
      END IF

   END DO

   END SUBROUTINE

!====================================================================

   SUBROUTINE ENDOPEP(NS,NSYN,SYNLOOP,LOOP,LOOPLEN,NBONDPEP,BONDPEP,BONTYP,SYNPG,newPGID,newPGLEN,PARTNER,EDHOLD, &
              pepdir,syndir,EDCAP,X,Y,Z,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD,JRANDS,RANDS,JSKIP,glytip,l_edcap)


   IMPLICIT NONE

   INTEGER,VALUE:: NS,NSYN,JSKIP,NBONDPEP
   INTEGER NS1,N1,N2,ILOOP,IPG1,IPG2,JR,JCYCLE,JL,NB,CHECK,JRANDS,NBGET,LENGTH

   INTEGER, ALLOCATABLE, DIMENSION(:),INTENT(IN)::SYNLOOP,LOOPLEN,BONTYP,newPGLEN,EDCAP,pepdir,syndir
   INTEGER, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::LOOP,BONDPEP,newPGID,SYNPG,PARTNER,glytip
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::EDHOLD

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z,XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::RANDS,XEDASEOLD,YEDASEOLD,ZEDASEOLD
   DOUBLE PRECISION DX,DY,DZ,R
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   DOUBLE PRECISION,value::l_edcap

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::pcan
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::cand

   DOUBLE PRECISION x0,y0,z0,xc,yc,zc,psum,d2,d2cap

   INTEGER ncan,ntip1,ntip2,j

!print*,'endopep'

   d2cap=l_edcap*l_edcap

   XL1=XEDASE(NS)
   YL1=YEDASE(NS)
   ZL1=ZEDASE(NS)

   XL2=XEDASEOLD(NS)
   YL2=YEDASEOLD(NS)
   ZL2=ZEDASEOLD(NS)

!  IF THE COMPLEX IS INACTIVE:

   IF(EDCAP(NS)==2)THEN

      DO NB=1,NBONDPEP

         IF(BONTYP(NB)==0)THEN
            CYCLE
         END IF

         JCYCLE=0

         DO NS1=1,NSYN
            IF(EDHOLD(3,NS1)==NB)THEN
               JCYCLE=1
               EXIT
            END IF
         END DO

         IF(JCYCLE==1)THEN
            CYCLE
         END IF

         N1=BONDPEP(1,NB)
         N2=BONDPEP(2,NB)

         DX=XL1-0.5D0*(X(N1)+X(N2))
         DY=YL1-0.5D0*(Y(N1)+Y(N2))
         DZ=ZL1-0.5D0*(Z(N1)+Z(N2))

         IF(DX*DX+DY*DY+DZ*DZ<0.01D0)THEN

!         XP3=0.5D0*(X(N1)+X(N2))
!         YP3=0.0D0
!         ZP3=0.0D0

!         XP1=2*X(N1)-XP3
!         YP1=2*Y(N1)
!         ZP1=2*Z(N1)

!         XP2=2*X(N2)-XP3
!         YP2=2*Y(N2)
!         ZP2=2*Z(N2)

!         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

!         IF(CHECK==1)THEN

            EDHOLD(1,NS)=N1
            EDHOLD(2,NS)=N2
            EDHOLD(3,NS)=NB

            EXIT

         END IF

      END DO

      RETURN
   END IF

!----------------------------------


   ILOOP=SYNLOOP(NS)


   IPG1=SYNPG(1,NS)
   IPG2=SYNPG(2,NS)

   LENGTH=LOOPLEN(ILOOP)

   N2=LOOP(LENGTH,ILOOP)

   allocate(pcan(length),cand(3,length))

   ncan=0

   ntip1=glytip(1,ns)

   ntip2=glytip(2,ns)

   x0=0.5*(x(ntip1)+x(ntip2))
   y0=0.5*(y(ntip1)+y(ntip2))
   z0=0.5*(z(ntip1)+z(ntip2))


   DO JL=1,LENGTH

      N1=N2

      N2=LOOP(JL,ILOOP)


!     CHECK IF THIS IS A PEPTIDE BOND:

      IF(N1/=PARTNER(1,N2).AND.N1/=PARTNER(2,N2))THEN
         CYCLE
      END IF

!     pepdir of the first bead must be the same as SYNDIR of the complex

      if(pepdir(n1)/=syndir(ns)) cycle


!     NOT CLEAVING NEW BONDS:

      JCYCLE=0

      IF(IPG1>0)THEN

         DO JR=1,newPGLEN(IPG1)

            IF(newPGID(JR,IPG1)==N1.OR.newPGID(JR,IPG1)==N2)THEN
               JCYCLE=1
               EXIT
            END IF

         END DO

         IF(JCYCLE==1)THEN
            CYCLE
         END IF

      END IF

      IF(IPG2>0)THEN

         DO JR=1,newPGLEN(IPG2)

            IF(newPGID(JR,IPG2)==N1.OR.newPGID(JR,IPG2)==N2)THEN
               JCYCLE=1
               EXIT
            END IF

         END DO

         IF(JCYCLE==1)THEN
            CYCLE
         END IF

      END IF



      JRANDS=JRANDS+1

      R=RANDS(JRANDS)



!     GET THE BOND:


      DO NB=1,NBONDPEP

         IF(BONTYP(NB)==0)THEN
            CYCLE
         END IF

         IF(BONDPEP(1,NB)==N1.AND.BONDPEP(2,NB)==N2)THEN
!            JGET=1
            NBGET=NB
            EXIT
         END IF

         IF(BONDPEP(1,NB)==N2.AND.BONDPEP(2,NB)==N1)THEN
!            JGET=1
            NBGET=NB
            EXIT
         END IF

      END DO

      ncan=ncan+1

      cand(1,ncan)=n1
      cand(2,ncan)=n2
      cand(3,ncan)=nbget


      xc=0.5*(x(n1)+x(n2))
      yc=0.5*(y(n1)+y(n2))
      zc=0.5*(z(n1)+z(n2))

      d2=(x0-xc)**2+(y0-yc)**2+(z0-zc)**2

      if(d2>4*d2cap) cycle

      pcan(ncan)=1.0/d2

      if(d2>d2cap) pcan(ncan)=pcan(ncan)*1e-10

!     to make deviation angle from the track a factor:

      pcan(ncan)=pcan(ncan)/max(1.0d0,(x0-xc)**2)


   end do


   if(ncan==0) return


   if(ncan==1)then

      edhold(1:3,ns)=cand(1:3,ncan)

   else

      psum=sum(pcan(1:ncan))

      pcan(1:ncan)=pcan(1:ncan)/psum

      psum=0.0

      do j=1,ncan

         psum=psum+pcan(j)

         if(psum>r)then

            edhold(1:3,ns)=cand(1:3,j)

            exit

         end if

      end do

   end if


!      NB=PARTNER(2,N1)



!      CALL RANDOM_NUMBER(R)

!      JRANDS=JRANDS+1

!      R=RANDS(JRANDS)

!      IF(0.1D0*JSKIP>R)THEN

!         EDHOLD(1,NS)=N1
!         EDHOLD(2,NS)=N2
!         EDHOLD(3,NS)=NBGET

!         EXIT
!      END IF

!   END DO

!print*,'endopep end'

   END SUBROUTINE

!====================================================================

   SUBROUTINE CLEAVESIGNAL(NS,EDHOLD,X,Y,Z,SIGCLEAVE,EDLOCKIN,BONTYP)

   IMPLICIT NONE

   INTEGER NS,N1,N2,NB
   INTEGER, ALLOCATABLE, DIMENSION(:)::EDLOCKIN,SIGCLEAVE,BONTYP
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::EDHOLD
   DOUBLE PRECISION DX,DY,DZ,DIST
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE
!   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3


   N1=EDHOLD(1,NS)
   N2=EDHOLD(2,NS)
   NB=EDHOLD(3,NS)

   IF(EDLOCKIN(NS)==1)THEN

!      IF(BONTYP(NB)==2)THEN
         SIGCLEAVE(NS)=1
!      END IF

      RETURN
   END IF

   DX=X(N1)-X(N2)
   DY=Y(N1)-Y(N2)
   DZ=Z(N1)-Z(N2)

   DIST=SQRT(DX*DX+DY*DY+DZ*DZ)


   IF(DIST>1.5D0)THEN
      EDLOCKIN(NS)=1
   END IF


   END SUBROUTINE

!====================================================================

   SUBROUTINE POSTCLEAVE(NS,NSYN,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOP,NLOOPDEL,EDHOLD,SIGCLEAVE, &
              JDEACT,X,Y,Z,NBONDDEL,NBONDPEP,BONDPEP,BONTYP,DNOR,ATOR,EDPEP,EDCAP,JSTEP,J_EDHOLD,PARTNER, &
              noldbond,oldbond,oldtyp)

   IMPLICIT NONE

   INTEGER(kind=8)::jstep
   INTEGER NS,NS1,NSYN,NLOOP,NLOOPDEL,NBONDDEL,NBONDPEP,NB,NB0,NB1,NB2
   INTEGER ILOOP,N1,N2,JL,NL,IPICK,J,J1,J2,JDEL,JSTART,JA1,JA2,NL1,NL2
   INTEGER NPICK,LOLEN,CHECK,JCHECK,LEN1,LEN2,JALTER,NA1,NA2

   INTEGER, ALLOCATABLE, DIMENSION(:,:):: LOOP,EDHOLD,BONDPEP,EDPEP,PARTNER
   INTEGER, ALLOCATABLE, DIMENSION(:)::SYNLOOP,LOOPLEN,LOOPTYP,PART1,PART2,BONTYP,JDEACT
   INTEGER, ALLOCATABLE, DIMENSION(:)::SIGCLEAVE,DNOR,ATOR,EDCAP
   INTEGER(kind=8),ALLOCATABLE, DIMENSION(:)::J_EDHOLD

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   integer noldbond
   INTEGER, ALLOCATABLE, DIMENSION(:,:)::oldbond
   INTEGER, ALLOCATABLE, DIMENSION(:)::oldtyp


   SIGCLEAVE(NS)=0

!   EDLOCKIN(NS)=0

!   EDCAP(NS)=0

!   ILOOP=SYNLOOP(NS)

   N1=EDHOLD(1,NS)
   N2=EDHOLD(2,NS)
   NB=EDHOLD(3,NS)

   IF(N1==0.OR.N2==0.OR.NB==0)THEN
      PRINT*,'NO HOLD FOR ENDOPEP',N1,N2,NB
      PRINT*,'COMPLEX',NS
      STOP
   END IF

!  IF THE COMPLEX IS INACTIVE:

   IF(EDCAP(NS)==2)THEN

      EDCAP(NS)=0

      EDHOLD(1:3,NS)=0

      N1=BONDPEP(1,NB)
      N2=BONDPEP(2,NB)

      IF(BONTYP(NB)==0)THEN
         PRINT*,'ERROR: BOND ALREADY DELETED',NB,BONDPEP(1:2,NB)
         STOP
      END IF

      BONTYP(NB)=0

!     make it oldbond

      noldbond=noldbond+1

      oldbond(1:2,noldbond)=BONDPEP(1:2,NB)

      oldtyp(noldbond)=1

      NBONDDEL=NBONDDEL+1

      IF(PARTNER(1,N1)==N2)THEN
         PARTNER(1,N1)=PARTNER(2,N1)
      END IF

      PARTNER(2,N1)=0


      IF(PARTNER(1,N2)==N1)THEN
         PARTNER(1,N2)=PARTNER(2,N2)
      END IF

      PARTNER(2,N2)=0


      print*,'cleaving pep bond while inactive'

      DNOR(N1)=-1
      ATOR(N2)=1

!     turn off signal bond:

      IF(DNOR(N2)==0)THEN

         DO NB0=1,NBONDPEP

            IF(BONTYP(NB0)/=2)THEN
               CYCLE
            END IF

            IF(BONDPEP(1,NB0)==N2)THEN
               BONTYP(NB0)=1
               EXIT
            END IF

         END DO

      END IF

      NL1=0; NL2=0

!     THEN MERGE TWO LOOPS:

      DO NL=1,NLOOP

         IF(LOOPTYP(NL)==0)THEN
            CYCLE
         END IF

         DO JL=1,LOOPLEN(NL)

            JA1=LOOP(JL,NL)

            IF(JL==LOOPLEN(NL))THEN
               JA2=LOOP(1,NL)
            ELSE
               JA2=LOOP(JL+1,NL)
            END IF

            IF(JA1==N1.AND.JA2==N2)THEN
               NL1=NL

               IF(JL==LOOPLEN(NL))THEN
                  J1=1
               ELSE
                  J1=JL+1
               END IF

               EXIT
            END IF

            IF(JA2==N1.AND.JA1==N2)THEN
               NL2=NL

               IF(JL==LOOPLEN(NL))THEN
                  J2=1
               ELSE
                  J2=JL+1
               END IF

               EXIT
            END IF

         END DO

         IF(NL1>0.AND.NL2>0)THEN
            EXIT
         END IF

      END DO

      IF(NL1==0.OR.NL2==0)THEN
         RETURN
      END IF

!     IF THE BOND IS NOT CONNECTING TWO LOOPS:

      IF(NL1==NL2)THEN

         LOLEN=LOOPLEN(NL1)

         J1=0; J2=0

         DO JL=1,LOLEN

            JA1=LOOP(JL,NL1)

            IF(JL==LOLEN)THEN
               JA2=LOOP(1,NL1)
            ELSE
               JA2=LOOP(JL+1,NL1)
            END IF

            IF(JA1==N2.AND.JA2==N1)THEN

               IF(JL==LOLEN)THEN
                  J1=1
               ELSE
                  J1=JL+1
               END IF

            END IF

            IF(JA1==N1.AND.JA2==N2)THEN

               IF(JL==LOLEN)THEN
                  J2=1
               ELSE
                  J2=JL+1
               END IF

            END IF

         END DO

         IF(J1==0.OR.J2==0)THEN
            PRINT*,'COULD NOT FIND LOOP INDICES with inactive complex'
            STOP
         END IF

         ALLOCATE(PART1(LOLEN),PART2(LOLEN))

         LEN1=0

         DO JL=1,LOLEN

            J=JL+J1-1

            IF(J>LOLEN)THEN
               J=J-LOLEN
            END IF

            LEN1=LEN1+1

            PART1(LEN1)=LOOP(J,NL1)

            IF(PART1(LEN1)==N2)THEN
               LEN1=LEN1-2
               EXIT
            END IF

         END DO

         LEN2=0

         DO JL=1,LOLEN

            J=JL+J2-1

            IF(J>LOLEN)THEN
               J=J-LOLEN
            END IF

            LEN2=LEN2+1

            PART2(LEN2)=LOOP(J,NL1)

            IF(PART2(LEN2)==N1)THEN
               LEN2=LEN2-2
               EXIT
            END IF

         END DO

         IF(LEN1<LEN2)THEN
            LOOP(1:LEN2,NL1)=PART2(1:LEN2)
            LOOPLEN(NL1)=LEN2
            DEALLOCATE(PART1,PART2)

         ELSE
            LOOP(1:LEN1,NL1)=PART1(1:LEN1)
            LOOPLEN(NL1)=LEN1
            DEALLOCATE(PART1,PART2)

         ENDIF

         RETURN

      END IF

!--------------------------

   IF(LOOPTYP(NL1)==0.OR.LOOPTYP(NL2)==0)THEN
      PRINT*,'ERROR: AT LEAST ONE LOOP ALREADY DELETED'
      PRINT*,LOOPTYP(NL1),LOOP(1:LOOPLEN(NL1),NL1)
      PRINT*,LOOPTYP(NL2),LOOP(1:LOOPLEN(NL2),NL2)
      STOP
   END IF

      LOOPTYP(NL1)=0
      LOOPTYP(NL2)=0
      NLOOPDEL=NLOOPDEL+2

      NLOOP=NLOOP+1

      LOOPTYP(NLOOP)=1

      LOLEN=0

      DO JL=1,LOOPLEN(NL1)

         J=JL+J1-1

         IF(J>LOOPLEN(NL1))THEN
            J=J-LOOPLEN(NL1)
         END IF

         LOLEN=LOLEN+1

         LOOP(LOLEN,NLOOP)=LOOP(J,NL1)

      END DO

      DO JL=2,LOOPLEN(NL2)-1

         J=JL+J2-1

         IF(J>LOOPLEN(NL2))THEN
            J=J-LOOPLEN(NL2)
         END IF

         LOLEN=LOLEN+1

         LOOP(LOLEN,NLOOP)=LOOP(J,NL2)

      END DO

      LOOPLEN(NLOOP)=LOLEN

      DO NS1=1,NSYN

         IF(SYNLOOP(NS1)==NL1.OR.SYNLOOP(NS1)==NL2)THEN
            SYNLOOP(NS1)=NLOOP

            LOOPTYP(NLOOP)=2
         END IF

      END DO

      RETURN

   END IF

!--------------------------------------------------------

!  CLEAVE BOND WHEN COMPLEX IS ACTIVE:


!   EDLOCKIN(NS)=0

   EDCAP(NS)=0

   ILOOP=SYNLOOP(NS)


   JSTART=0

   DO JL=1,LOOPLEN(ILOOP)

      NA2=LOOP(JL,ILOOP)

      IF(JL==1)THEN
         NA1=LOOP(LOOPLEN(ILOOP),ILOOP)
      ELSE
         NA1=LOOP(JL-1,ILOOP)
      END IF

      IF(NA1==N1.AND.NA2==N2)THEN
         JSTART=JL
         EXIT
      END IF

   END DO

   IF(JSTART==0)THEN
      PRINT*,'WARNING: EDHOLD DOES NOT MATCH SYNLOOP FOR SYNCOMP',NS
      PRINT*,'EDHOLD',EDHOLD(1:3,NS)
      PRINT*,'LOOP',LOOP(1:LOOPLEN(ILOOP),ILOOP)
      EDHOLD(1:3,NS)=0
      RETURN
   END IF

   IPICK=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      DO JL=1,LOOPLEN(NL)

         NA1=LOOP(JL,NL)

         IF(JL==LOOPLEN(NL))THEN
            NA2=LOOP(1,NL)
         ELSE
            NA2=LOOP(JL+1,NL)
         END IF

         IF(NA1==N2.AND.NA2==N1)THEN
            IPICK=NL
            EXIT
         END IF

      END DO

      IF(IPICK>0)THEN
         EXIT
      END IF

   END DO

   IF(IPICK==0)THEN
      PRINT*,'CANNOT FIND IPICK FOR COMPLEX',NS,N1,N2
      STOP
   END IF

!----------------------------

!  IF TWO COMPLEXES MEET:

!   IF(ILOOP/=IPICK.AND.LOOPTYP(IPICK)==2)THEN


!         DO NS1=1,NSYN

!            IF(SYNLOOP(NS1)==IPICK.AND.NS1/=NS)THEN
!               JDEACT(NS1)=1
!               EXIT
!            END IF

!         END DO


!   END IF

!-------------------------


      IF(BONTYP(NB)==0)THEN
         PRINT*,'ERROR: BOND ALREADY DELETED',NB,BONDPEP(1:2,NB)
         STOP
      END IF

      noldbond=noldbond+1

      oldbond(1:2,noldbond)=BONDPEP(1:2,NB)

      oldtyp(noldbond)=1


   NBONDDEL=NBONDDEL+1

!print*,'nbonddel=',nbonddel

   BONTYP(NB)=0

   NB1=BONDPEP(1,NB)
   NB2=BONDPEP(2,NB)

      IF(PARTNER(1,NB1)==NB2)THEN
         PARTNER(1,NB1)=PARTNER(2,NB1)
      END IF

      PARTNER(2,NB1)=0


      IF(PARTNER(1,NB2)==NB1)THEN
         PARTNER(1,NB2)=PARTNER(2,NB2)
      END IF

      PARTNER(2,NB2)=0


   DNOR(NB1)=-1
   ATOR(NB2)=1

   EDHOLD(1:3,NS)=0

   IF(ATOR(NB1)==1)THEN
      EDPEP(1,NS)=NB1
   END IF

   EDPEP(2,NS)=NB2

   J_EDHOLD(NS)=JSTEP

!   IF(ATOR(N1)==0)THEN
!      N0=N1
!   ELSE
!      N0=N2
!   END IF

   IF(DNOR(NB2)==0)THEN

      DO NB0=1,NBONDPEP

         IF(BONTYP(NB0)/=2)THEN
            CYCLE
         END IF

         IF(BONDPEP(1,NB0)==NB2)THEN
            BONTYP(NB0)=1
            EXIT
         END IF

      END DO

   END IF

!   IF(DNOR(N1)==0)THEN
!      DNOR(N1)=-1
!      ATOR(N2)=1
!   ELSE
!      ATOR(N1)=1
!      DNOR(N2)=-1
!   END IF

!---------------------------------------
!  IF THE BOND IS NOT CONNECTING TWO LOOPS:

   IF(IPICK==ILOOP)THEN

      PRINT*,'SYSTEM MAY BE A MESS AT SYNCOMP',NS
      LOLEN=LOOPLEN(ILOOP)

!----------------------------
!     New method to separate two parts of the loop:

      J1=0; J2=0

      DO JL=1,LOLEN

         JA1=LOOP(JL,ILOOP)

         IF(JL==LOLEN)THEN
            JA2=LOOP(1,ILOOP)
         ELSE
            JA2=LOOP(JL+1,ILOOP)
         END IF

         IF(JA1==N2.AND.JA2==N1)THEN

            IF(JL==LOLEN)THEN
               J1=1
            ELSE
               J1=JL+1
            END IF

         END IF

         IF(JA1==N1.AND.JA2==N2)THEN

            IF(JL==LOLEN)THEN
               J2=1
            ELSE
               J2=JL+1
            END IF

         END IF

      END DO

      IF(J1==0.OR.J2==0)THEN
         PRINT*,'COULD NOT FIND LOOP INDICES'
         PRINT*,'CROSSLINK AT',N1,N2
         PRINT*,'LOOP',LOOP(1:LOLEN,ILOOP)
         STOP
      END IF

      ALLOCATE(PART1(LOLEN),PART2(LOLEN))

      LEN1=0

      DO JL=1,LOLEN

         J=JL+J1-1

         IF(J>LOLEN)THEN
            J=J-LOLEN
         END IF

         LEN1=LEN1+1

         PART1(LEN1)=LOOP(J,ILOOP)

         IF(PART1(LEN1)==N2)THEN
            LEN1=LEN1-2
            EXIT
         END IF

      END DO

      LEN2=0

      DO JL=1,LOLEN

         J=JL+J2-1

         IF(J>LOLEN)THEN
            J=J-LOLEN
         END IF

         LEN2=LEN2+1

         PART2(LEN2)=LOOP(J,ILOOP)

         IF(PART2(LEN2)==N1)THEN
            LEN2=LEN2-2
            EXIT
         END IF

      END DO

!----------------------------

!      J1=0

!      J2=0

!      DO JL=1,LOLEN

!         IF(LOOP(JL,ILOOP)==N1.AND.J1==0)THEN
!            J1=JL

!         ELSEIF(LOOP(JL,ILOOP)==N1.AND.J1>0)THEN
!            J2=JL
!         END IF

!      END DO

!      IF(J2==0)THEN

!         LOOP(J1:LOLEN-2,ILOOP)=LOOP(J1+2:LOLEN,ILOOP)

!         LOOPLEN(ILOOP)=LOLEN-2

!         RETURN

!      END IF

!      ALLOCATE(PART1(LOLEN))

!      JALTER=0

!      JDEL=J2-J1

!      IF(JDEL<0)THEN
!         JDEL=JDEL+LOLEN
!      END IF

!      LEN1=0

!      DO JL=1,JDEL

!         J=JL+J1+1

!         IF(J>LOLEN)THEN
!            J=J-LOLEN
!         END IF

!         IF(LOOP(J,ILOOP)==N2)THEN
!            JALTER=1
!            EXIT
!         END IF

!         LEN1=LEN1+1

!         PART1(LEN1)=LOOP(J,ILOOP)

!      END DO

!      IF(JALTER==1)THEN

!         JDEL=J1-J2

!         IF(JDEL<0)THEN
!            JDEL=JDEL+LOLEN
!         END IF

!         LEN1=0

!         DO JL=1,JDEL

!            J=JL+J2-1

!            IF(J>LOLEN)THEN
!               J=J-LOLEN
!            END IF

!            LEN1=LEN1+1

!            PART1(LEN1)=LOOP(J,ILOOP)

!         END DO

!      END IF

!-----

!      J1=0

!      J2=0

!      DO JL=1,LOLEN

!         IF(LOOP(JL,ILOOP)==N2.AND.J1==0)THEN
!            J1=JL

!         ELSEIF(LOOP(JL,ILOOP)==N2.AND.J1>0)THEN
!            J2=JL
!         END IF

!      END DO
 
!      IF(J2==0)THEN

!         LOOP(J1:LOLEN-2,ILOOP)=LOOP(J1+2:LOLEN,ILOOP)

!         LOOPLEN(ILOOP)=LOLEN-2

!         DEALLOCATE(PART1)

!         RETURN

!      END IF

!      ALLOCATE(PART2(LOLEN))

!      JALTER=0

!      JDEL=J2-J1

!      IF(JDEL<0)THEN
!         JDEL=JDEL+LOLEN
!      END IF

!      LEN2=0

!      DO JL=1,JDEL

!         J=JL+J1-1

!         IF(J>LOLEN)THEN
!            J=J-LOLEN
!         END IF

!         IF(LOOP(J,ILOOP)==N1)THEN
!            JALTER=1
!            EXIT
!         END IF

!         LEN2=LEN2+1

!         PART2(LEN2)=LOOP(J,ILOOP)

!      END DO

!      IF(JALTER==1)THEN

!         JDEL=J1-J2

!         IF(JDEL<0)THEN
!            JDEL=JDEL+LOLEN
!         END IF

!         LEN2=0

!         DO JL=1,JDEL

!            J=JL+J2-1

!            IF(J>LOLEN)THEN
!               J=J-LOLEN
!            END IF

!            LEN2=LEN2+1

!            PART2(LEN2)=LOOP(J,ILOOP)

!         END DO

!      END IF

!-----

      IF(LEN1<4)THEN
         LOOP(1:LEN2,ILOOP)=PART2(1:LEN2)
         LOOPLEN(ILOOP)=LEN2
         DEALLOCATE(PART1,PART2)
         RETURN
      ELSEIF(LEN2<4)THEN
         LOOP(1:LEN1,ILOOP)=PART1(1:LEN1)
         LOOPLEN(ILOOP)=LEN1
         DEALLOCATE(PART1,PART2)
         RETURN
      ENDIF

!     CHECK IF PART1 IS A LOOP:

      XL1=0.5D0*(X(N1)+X(N2))
      YL1=Y(N1)+Y(N2)
      ZL1=Z(N1)+Z(N2)

      XL2=XL1; YL2=0.0D0; ZL2=0.0D0

      XP3=X(PART1(1))
      YP3=Y(PART1(1))
      ZP3=Z(PART1(1))

      JCHECK=0

      DO J=2,LEN1-1

         NA1=PART1(J)
         NA2=PART1(J+1)

         IF(PART1(1)==NA1.OR.PART1(1)==NA2)THEN
            CYCLE
         END IF

         XP1=X(NA1)
         YP1=Y(NA1)
         ZP1=Z(NA1)

         XP2=X(NA2)
         YP2=Y(NA2)
         ZP2=Z(NA2)

         CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

         IF(CHECK==1)THEN
            JCHECK=JCHECK+1
         END IF

      END DO

      IF(MOD(JCHECK,2)==1)THEN
         LOOP(1:LEN1,ILOOP)=PART1(1:LEN1)
         LOOPLEN(ILOOP)=LEN1
      ELSE
         LOOP(1:LEN2,ILOOP)=PART2(1:LEN2)
         LOOPLEN(ILOOP)=LEN2
      END IF

      DEALLOCATE(PART1,PART2)

      RETURN

   END IF

!----------------------------------------

!  LET ENDOPEPTIDASE HOLD NB2:


   IF(LOOPTYP(ILOOP)==0.OR.LOOPTYP(IPICK)==0)THEN
      PRINT*,'ERROR: LOOP ALREADY DELETED'
      PRINT*,LOOPTYP(ILOOP),LOOP(1:LOOPLEN(ILOOP),ILOOP)
      PRINT*,LOOPTYP(IPICK),LOOP(1:LOOPLEN(IPICK),IPICK)
      STOP
   END IF

   NLOOPDEL=NLOOPDEL+2
   LOOPTYP(ILOOP)=0
   LOOPTYP(IPICK)=0

   NLOOP=NLOOP+1
   LOOPTYP(NLOOP)=2
   SYNLOOP(NS)=NLOOP

   DO NS1=1,NSYN

      IF(SYNLOOP(NS1)==IPICK)THEN
         SYNLOOP(NS1)=NLOOP
         EXIT
      END IF

   END DO


!   DO JL=1,LOOPLEN(ILOOP)

!      NA2=LOOP(JL,ILOOP)

!      IF(JL==1)THEN
!         NA1=LOOP(LOOPLEN(ILOOP),ILOOP)
!      ELSE
!         NA1=LOOP(JL-1,ILOOP)
!      END IF

!      IF(NA1==N1.AND.NA2==N2)THEN
!         JSTART=JL
!         EXIT
!      END IF

!   END DO

   LOLEN=0

   DO JL=1,LOOPLEN(ILOOP)

      J=JL-1+JSTART

      IF(J>LOOPLEN(ILOOP))THEN
         J=J-LOOPLEN(ILOOP)
      END IF

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=LOOP(J,ILOOP)

   END DO

!--
   DO JL=1,LOOPLEN(IPICK)

      NA2=LOOP(JL,IPICK)

      IF(JL==1)THEN
         NA1=LOOP(LOOPLEN(IPICK),IPICK)
      ELSE
         NA1=LOOP(JL-1,IPICK)
      END IF

      IF(NA1==N2.AND.NA2==N1)THEN
         JSTART=JL
         EXIT
      END IF

   END DO

   DO JL=2,LOOPLEN(IPICK)-1

      J=JL-1+JSTART

      IF(J>LOOPLEN(IPICK))THEN
         J=J-LOOPLEN(IPICK)
      END IF

      LOLEN=LOLEN+1

      LOOP(LOLEN,NLOOP)=LOOP(J,IPICK)

   END DO

   LOOPLEN(NLOOP)=LOLEN


   END SUBROUTINE

!====================================================================

   SUBROUTINE PGCLEAN(IPG,newPGID,newPGLEN,ILOOP,LOOP,LOOPLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,DNOR,ATOR)

   IMPLICIT NONE

   INTEGER IPG,ILOOP,NBONDPEP,NBONDDEL
   INTEGER CHECK,JL,NA,J,NB
   INTEGER,DIMENSION(:,:),ALLOCATABLE::newPGID,BONDPEP,LOOP
   INTEGER,DIMENSION(:),ALLOCATABLE::newPGLEN,LOOPLEN,BONTYP,DNOR,ATOR

   CHECK=0

   DO JL=1,LOOPLEN(ILOOP)

      NA=LOOP(JL,ILOOP)

      DO J=1,newPGLEN(IPG)

         IF(newPGID(J,IPG)==NA)THEN
            CHECK=1
            EXIT
         END IF

      END DO

      IF(CHECK==1)THEN
         EXIT
      END IF

   END DO

   IF(CHECK==1)THEN
      RETURN
   END IF

   DO J=1,newPGLEN(IPG)

      NA=newPGID(J,IPG)

      IF(DNOR(NA)==0.OR.ATOR(NA)==0)THEN

         DO NB=1,NBONDPEP

            IF(BONTYP(NB)==0)THEN
               CYCLE
            END IF

            IF(BONDPEP(1,NB)==NA.OR.BONDPEP(2,NB)==NA)THEN

               BONTYP(NB)=0

               NBONDDEL=NBONDDEL+1

               DNOR(BONDPEP(1,NB))=-1

               ATOR(BONDPEP(2,NB))=1

               EXIT

            END IF

         END DO

      END IF

   END DO



   END SUBROUTINE

!====================================================================


   SUBROUTINE RMPEPBOND(NATOM,NPG,newPGID,newPGLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP, &
                        DNOR,ATOR,NSYN,SYNPG,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOP)

   IMPLICIT NONE

   INTEGER NATOM,NPG,NR,JR,NB,NB1,JCOUNT,NPICK,NBONDPEP,NBONDDEL,NSYN,NS,NL,JL,NLOOP,JGET,LENGTH,JCYCLE,NA
   INTEGER,DIMENSION(:,:),ALLOCATABLE::newPGID,BONDPEP,SYNPG,LOOP
   INTEGER,DIMENSION(:),ALLOCATABLE::newPGLEN,BONTYP,DNOR,ATOR,MARK,PGMARK,LOOPMARK,SYNLOOP,LOOPLEN,LOOPTYP

   if(npg==0) return

   ALLOCATE(MARK(NATOM),PGMARK(NPG),LOOPMARK(NATOM))

   MARK=0

   DO NB=1,NBONDPEP

      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      MARK(BONDPEP(1:2,NB))=1

   END DO

   PGMARK=0
   LOOPMARK=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      LOOPMARK(LOOP(1:LOOPLEN(NL),NL))=1

   END DO

   DO NS=1,NSYN

      IF(SYNPG(1,NS)>0)THEN
         PGMARK(SYNPG(1,NS))=1
      END IF

      IF(SYNPG(2,NS)>0)THEN
         PGMARK(SYNPG(2,NS))=1
      END IF

      IF(SYNLOOP(NS)>0)THEN
         NL=SYNLOOP(NS)
         LOOPMARK(LOOP(1:LOOPLEN(NL),NL))=2
      END IF

   END DO

   DO NR=1,NPG

      IF(PGMARK(NR)==1)THEN
         CYCLE
      END IF

      JCOUNT=0

      DO JR=1,newPGLEN(NR)

         NA=newPGID(JR,NR)

         IF(MARK(NA)==1)THEN
            NB1=NA
            JCOUNT=JCOUNT+1
         END IF

      END DO

      IF(JCOUNT/=1)THEN
         CYCLE
      END IF

      IF(LOOPMARK(NB1)==2)THEN
         CYCLE
      END IF

      NPICK=0

      DO NB=1,NBONDPEP

         IF(BONTYP(NB)==0)THEN
            CYCLE
         END IF

         IF(BONDPEP(1,NB)==NB1)THEN
            NPICK=NB
            EXIT
         END IF

         IF(BONDPEP(2,NB)==NB1)THEN
            NPICK=NB
            EXIT
         END IF

      END DO

      IF(BONTYP(NPICK)==0)THEN
         PRINT*,'ERROR: BOND ALREADY DELETED',NPICK,BONDPEP(1:2,NPICK)
         STOP
      END IF

      NBONDDEL=NBONDDEL+1

!print*,'nbonddel=',nbonddel

      BONTYP(NPICK)=0

      DNOR(BONDPEP(1,NPICK))=-1
      ATOR(BONDPEP(2,NPICK))=1

      MARK(BONDPEP(1,NPICK))=0
      MARK(BONDPEP(2,NPICK))=0

      IF(LOOPMARK(NB1)==0)THEN
         CYCLE
      END IF

!     OTHERWISE MODIFY LOOPS:

      DO NL=1,NLOOP

         IF(LOOPTYP(NL)==0)THEN
            CYCLE
         END IF

         JGET=0

         DO JL=1,LOOPLEN(NL)

            IF(LOOP(JL,NL)==NB1)THEN
               JGET=1
               EXIT
            END IF

         END DO

         IF(JGET==0)THEN
            CYCLE
         END IF

         LENGTH=0

         NA=0

         DO JL=1,LOOPLEN(NL)

            JCYCLE=0

            DO JR=1,newPGLEN(NR)

               IF(LOOP(JL,NL)==newPGID(JR,NR))THEN
                  JCYCLE=1
                  EXIT
               END IF

            END DO

            IF(JCYCLE==1)THEN
               CYCLE
            END IF

            IF(LOOP(JL,NL)==NA)THEN
               CYCLE
            END IF

            NA=LOOP(JL,NL)

            LENGTH=LENGTH+1

            LOOP(LENGTH,NL)=NA

         END DO

         IF(LOOP(LENGTH,NL)==LOOP(1,NL))THEN
            LENGTH=LENGTH-1
         END IF

         LOOPLEN(NL)=LENGTH

         EXIT

      END DO

   END DO

   DEALLOCATE(MARK,PGMARK,LOOPMARK)


   END SUBROUTINE

!====================================================================

   SUBROUTINE LYTGTASE(NATOM,NPG,newPGID,newPGLEN,DNOR,ATOR,PLYT,NATOMDEL,NSYN,SYNPG,SYNLOOP,LOOPLEN,LOOP, &
              noldbond,oldbond,oldtyp)

   IMPLICIT NONE

   INTEGER NATOM,NPG,NATOMDEL,NR,JR,JFREE1,JFREE2,JCRLK1,JCRLK2,JDEL,JCHECK,NSYN,NS,JCYCLE,ILOOP,IPG,JCOUNT,na,nb,j
   INTEGER,ALLOCATABLE,DIMENSION(:)::DNOR,ATOR,MARK,MARKLOOP,newPGLEN,LOOPLEN,SYNLOOP
   INTEGER,ALLOCATABLE,DIMENSION(:,:)::newPGID,LOOP,SYNPG
   DOUBLE PRECISION PLYT,R

   integer,value::noldbond
   INTEGER,ALLOCATABLE,DIMENSION(:,:),intent(in)::oldbond
   INTEGER,ALLOCATABLE,DIMENSION(:)::oldtyp

   if(npg==0) return

   ALLOCATE(MARK(NPG),MARKLOOP(NATOM))

   MARK=0
   MARKLOOP=0

   DO NS=1,NSYN

      ILOOP=SYNLOOP(NS)

      IF(ILOOP/=0)THEN

         MARKLOOP(LOOP(1:LOOPLEN(ILOOP),ILOOP))=1

      END IF

      IPG=SYNPG(1,NS)

      IF(IPG/=0)THEN

         JCHECK=SUM(MARKLOOP(newPGID(1:newPGLEN(IPG),IPG)))

         IF(JCHECK==0)THEN

            MARK(IPG)=1

         ELSE

            MARK(IPG)=2

         END IF

      END IF

      IPG=SYNPG(2,NS)

      IF(IPG/=0)THEN

         JCHECK=SUM(MARKLOOP(newPGID(1:newPGLEN(IPG),IPG)))

         IF(JCHECK==0)THEN

            MARK(IPG)=1

         ELSE

            MARK(IPG)=2

         END IF

      END IF


   END DO


   DO NR=1,NPG

      IF(MARK(NR)==1)THEN
         CYCLE
      END IF

!      JCYCLE=0

!      DO JR=1,PGLEN(NR)

!         IF(MARK(PGID(JR,NR))==1)THEN
!            JCYCLE=1
!            EXIT
!         END IF

!      END DO

!      IF(JCYCLE==1)THEN
!         CYCLE
!      END IF

!     FIRST DEAL WITH FREE STRANDS:

      JCHECK=0
      JDEL=0

      DO JR=1,newPGLEN(NR)

         IF(DNOR(newPGID(JR,NR))==0.OR.ATOR(newPGID(JR,NR))==0)THEN
            JCHECK=1
         ELSEIF(ATOR(newPGID(JR,NR))==-1)THEN
            JDEL=JDEL+1
         END IF

      END DO

      IF(JCHECK==0)THEN

         IF(JDEL==newPGLEN(NR))THEN
            CYCLE
         END IF

         CALL RANDOM_NUMBER(R)


         IF(PLYT<R)THEN
            CYCLE
         END IF

         IF(newPGLEN(NR)==1)THEN

            NATOMDEL=NATOMDEL+1

            na=newPGID(1,NR)

            DNOR(na)=-1
            ATOR(na)=-1


            do nb=1,noldbond

               if(na==oldbond(1,nb).or.na==oldbond(2,nb)) oldtyp(nb)=0

            end do

            CYCLE

         END IF

         DO JR=1,newPGLEN(NR)

            IF(ATOR(newPGID(JR,NR))/=-1)THEN
               JFREE1=JR
               EXIT
            END IF

         END DO

         DO JR=newPGLEN(NR),1,-1

            IF(ATOR(newPGID(JR,NR))/=-1)THEN
               JFREE2=JR
               EXIT
            END IF

         END DO

         CALL RANDOM_NUMBER(R)


         JDEL=(JFREE2-JFREE1+1)*R+1

         NATOMDEL=NATOMDEL+JDEL

         do j=JFREE1,JFREE1+JDEL-1

            na=newPGID(j,nr)

            DNOR(na)=-1
            ATOR(na)=-1


            do nb=1,noldbond

               if(na==oldbond(1,nb).or.na==oldbond(2,nb)) oldtyp(nb)=0
            
            end do

         end do

         CYCLE

      END IF

!---------------------------

!     CHECK FROM THE TAIL:

      JCRLK1=1

      DO JR=1,newPGLEN(NR)
         IF(DNOR(newPGID(JR,NR))==0.OR.ATOR(newPGID(JR,NR))==0)THEN
            JCRLK1=JR
            EXIT
         END IF
      END DO

      IF(JCRLK1==1)THEN
         GOTO 10
      END IF


      JFREE1=0

      DO JR=1,JCRLK1-1
         IF(ATOR(newPGID(JR,NR))==1)THEN
            JFREE1=JR
            EXIT
         END IF
      END DO

      IF(JFREE1==0)THEN
         GOTO 10
      END IF

      JCOUNT=SUM(MARKLOOP(newPGID(JFREE1:JCRLK1-1,NR)))
      IF(JCOUNT>0)THEN
         GOTO 10
      END IF

      CALL RANDOM_NUMBER(R)


      IF(PLYT<R)THEN
         GOTO 10
      END IF

      CALL RANDOM_NUMBER(R)


      JDEL=(JCRLK1-JFREE1)*R+1

      NATOMDEL=NATOMDEL+JDEL


      do j=JFREE1,JFREE1+JDEL-1

         na=newpgid(j,nr)

         DNOR(na)=-1
         ATOR(na)=-1


         do nb=1,noldbond
         
            if(na==oldbond(1,nb).or.na==oldbond(2,nb)) oldtyp(nb)=0

         end do

      end do
!--------------------------

!     CHECK FROM THE TIP:

10    IF(MARK(NR)==2)THEN
         CYCLE
      END IF

      JCRLK2=newPGLEN(NR)

      DO JR=newPGLEN(NR),1,-1

         IF(DNOR(newPGID(JR,NR))==0.OR.ATOR(newPGID(JR,NR))==0)THEN
            JCRLK2=JR
            EXIT
         END IF

      END DO

      IF(JCRLK2==newPGLEN(NR))THEN
         CYCLE
      END IF

      JFREE2=0

      DO JR=newPGLEN(NR),JCRLK2+1,-1

         IF(ATOR(newPGID(JR,NR))/=-1)THEN
            JFREE2=JR
            EXIT
         END IF

      END DO

      IF(JFREE2==0)THEN
         CYCLE
      END IF

      JCOUNT=SUM(MARKLOOP(newPGID(JCRLK2+1:JFREE2,NR)))
      IF(JCOUNT>0)THEN
         CYCLE
      END IF

      CALL RANDOM_NUMBER(R)


      IF(PLYT<R)THEN
         CYCLE
      END IF

      CALL RANDOM_NUMBER(R)


      JDEL=(JFREE2-JCRLK2)*R+1

      NATOMDEL=NATOMDEL+JDEL

      do j=JFREE2-JDEL+1,JFREE2

         na=newpgid(j,nr)

         DNOR(na)=-1
         ATOR(na)=-1


         do nb=1,noldbond

            if(na==oldbond(1,nb).or.na==oldbond(2,nb)) oldtyp(nb)=0

         end do

      end do

   END DO


   DEALLOCATE(MARK,MARKLOOP)

   END SUBROUTINE

!====================================================================

   subroutine updateoldbond(noldbond,oldbond,oldtyp,oldtypmax,joldbond,pgbelowbond, &
              pgbelowbondmax,nsyn,x,y,z,xedase,yedase,zedase,prmbond,ator)

   implicit none

   integer,value::nsyn,oldtypmax,joldbond,pgbelowbondmax

   integer noldbond,nnew,n,n1,n2,ns,jcleave

   integer, allocatable,dimension(:,:)::oldbond
   integer, allocatable,dimension(:)::oldtyp
   integer, allocatable,dimension(:),intent(in)::ator,pgbelowbond

   real(kind=8),allocatable,dimension(:),intent(in)::x,y,z,xedase,yedase,zedase

!   real(kind=8)::alpha0,sina0,cosa0,dalpha,delta
   real xb,yb,zb,dx,dy,dz!,alpha,arg,yr,zr

   real(kind=8),value::prmbond
   real(kind=8)::r


   if(noldbond==0) return


   jcleave=0

   do n=1,noldbond

      if(oldtyp(n)==0)then
         jcleave=1
         cycle
      end if

!     if the effect of old-bonds on new strands are not considered:

      if(joldbond==0)then

         n1=oldbond(1,n)
         n2=oldbond(2,n)

         if(ator(n1)==0.and.ator(n2)==0)then
            oldtyp(n)=0
            jcleave=1
            cycle

         end if
      end if



      oldtyp(n)=oldtyp(n)+1

      if(oldtyp(n)>oldtypmax) then


         call random_number(r)

         if(prmbond>r.or.prmbond*oldtyp(n)>1.0d0)then
            oldtyp(n)=0
            jcleave=1

         elseif(pgbelowbond(n)>pgbelowbondmax)then
            oldtyp(n)=0
            jcleave=1

         else

            n1=oldbond(1,n)
            n2=oldbond(2,n)

            xb=0.5*(x(n1)+x(n2))
            yb=0.5*(y(n1)+y(n2))
            zb=0.5*(z(n1)+z(n2))

            do ns=1,nsyn

               dx=xb-xedase(ns)
               dy=yb-yedase(ns)
               dz=zb-zedase(ns)

               if(dx*dx+dy*dy+dz*dz<4.0)then
                  oldtyp(n)=0
                  jcleave=1
                  exit
               end if

            end do

         end if

      end if

   end do




   if(jcleave==0) return


   nnew=0

   do n=1,noldbond

      if(oldtyp(n)==0) cycle

      nnew=nnew+1

      oldbond(1:2,nnew)=oldbond(1:2,n)

      oldtyp(nnew)=oldtyp(n)

   end do

   oldtyp(nnew+1:noldbond)=0

   noldbond=nnew

!   oldtyp(1:noldbond)=1

!   print*,'after removal, noldbond=',noldbond

!   deallocate(typ)

   end subroutine

!====================================================================

   SUBROUTINE UPDATESYS(natomold,NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,newPGID,newPGLEN, &
              NSYN,SYNPG,GLYTIP,GLYSEC,TPPEP,EDPEP,EDHOLD,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
              NBONDGLY,NBONDGLYOLD,NBONDCAP,NBONDPEP,NBONDDEL,BONDGLY,BONDPEP,BONTYP,SYNLOOP,PARTNER, &
              noldbond,oldbond,jpgnew)

   IMPLICIT NONE

   integer, value::natomold
   INTEGER NATOM,NDEL,NPG,NSYN,NLOOP,NBONDGLY,NBONDGLYOLD,NBONDCAP,NBONDPEP,NBONDDEL,noldbond,jpgnew
   INTEGER N,N1,N2,NTEMP,NR,JR,LENGTH,NS,NL,JL,NA,NPGTEMP,J,NB,NNEW,NLOOPDEL
   INTEGER,ALLOCATABLE,DIMENSION(:)::DNOR,ATOR,PEPDIR,newPGLEN,LOOPLEN,LOOPTYP,MAP
   INTEGER,ALLOCATABLE,DIMENSION(:,:)::SYNPG,GLYTIP,GLYSEC,EDPEP,PARTNER
   INTEGER,ALLOCATABLE,DIMENSION(:)::BONTYP,DNORTEMP,ATORTEMP,SYNLOOP
   INTEGER,ALLOCATABLE,DIMENSION(:,:)::newPGID,LOOP,BONDGLY,BONDPEP,EDHOLD,oldbond
   INTEGER,ALLOCATABLE,DIMENSION(:,:,:)::TPPEP
   DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::X,Y,Z


   ALLOCATE(MAP(NATOM),DNORTEMP(NATOM-NDEL),ATORTEMP(NATOM-NDEL))

   do n=1,natomold
      map(n)=n
   end do

   NTEMP=natomold

   DO N=natomold+1,NATOM
      IF(ATOR(N)/=-1)THEN
         NTEMP=NTEMP+1
         MAP(N)=NTEMP
         DNORTEMP(NTEMP)=DNOR(N)
         ATORTEMP(NTEMP)=ATOR(N)
         PEPDIR(NTEMP)=PEPDIR(N)
         X(NTEMP)=X(N)
         Y(NTEMP)=Y(N)
         Z(NTEMP)=Z(N)
      END IF
   END DO

   IF(NATOM-NDEL/=NTEMP)THEN
      PRINT*,'ERROR IN COUNTING DELETED ATOMS',NATOM,NDEL,NTEMP
      STOP
   END IF

   NATOM=NTEMP

!-----------

!  shorten the strands:

   DO NR=1,NPG

      LENGTH=0


      DO JR=1,newPGLEN(NR)
         IF(ATOR(newPGID(JR,NR))==-1)THEN
            CYCLE
         END IF

         LENGTH=LENGTH+1

         NA=newPGID(JR,NR)

         newPGID(LENGTH,NR)=MAP(NA)

      END DO

      newPGLEN(NR)=LENGTH

   END DO

!----------------------------------
!  CLEAN UP UNCROSSLINKED STRANDS:

   NPGTEMP=0

   DO NR=1,NPG

      IF(newPGLEN(NR)==0)THEN
         CYCLE
      END IF

      NPGTEMP=NPGTEMP+1

      newPGLEN(NPGTEMP)=newPGLEN(NR)

!      PGTYP(NPGTEMP)=PGTYP(NR)

!      PGDIR(NPGTEMP)=PGDIR(NR)

      newPGID(1:newPGLEN(NR),NPGTEMP)=newPGID(1:newPGLEN(NR),NR)


      DO NS=1,NSYN

         IF(SYNPG(1,NS)==NR)THEN
            SYNPG(1,NS)=NPGTEMP
         END IF

         IF(SYNPG(2,NS)==NR)THEN
            SYNPG(2,NS)=NPGTEMP
         END IF

      END DO

   END DO

!print*,'update',npgtemp,npg

   NPG=NPGTEMP

!   do nr=1,npg
!      if(pgtyp(nr)==3)then
!         jpgnew=nr
!         exit
!      end if
!   end do
!-----------------------------

!  updating glycan bonds:

   NBONDGLY=NBONDGLYOLD

   DO NR=1,NPG

      LENGTH=newPGLEN(NR)-1

      IF(LENGTH>0)THEN

         BONDGLY(1,NBONDGLY+1:NBONDGLY+LENGTH)=newPGID(1:LENGTH,NR)
         BONDGLY(2,NBONDGLY+1:NBONDGLY+LENGTH)=newPGID(2:LENGTH+1,NR)

         NBONDGLY=NBONDGLY+LENGTH

      END IF


   END DO

!-------------
   NNEW=0

   DO NB=1,NBONDPEP

      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1

      BONDPEP(1:2,NNEW)=MAP(BONDPEP(1:2,NB))

      BONTYP(NNEW)=BONTYP(NB)

      DO NS=1,NSYN

         IF(EDHOLD(3,NS)==NB)THEN
            EDHOLD(3,NS)=NNEW
            EXIT
         END IF

      END DO


   END DO


   IF(NNEW/=NBONDPEP-NBONDDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED BONDS',NNEW,NBONDPEP,NBONDDEL
      STOP
   END IF


!   PRINT*,'DELETED BONDS',NBONDDEL

   NBONDPEP=NNEW
   NBONDDEL=0
!print*,'nbonddel=',nbonddel

!   PRINT*,'UPDATED BONDS',NBONDPEP

   PARTNER=0

   DO NB=1,NBONDPEP

         N1=BONDPEP(1,NB)
         N2=BONDPEP(2,NB)

         IF(PARTNER(1,N1)==0)THEN
            PARTNER(1,N1)=N2
         ELSE
            PARTNER(2,N1)=N2
         END IF

         IF(PARTNER(1,N2)==0)THEN
            PARTNER(1,N2)=N1
         ELSE
            PARTNER(2,N2)=N1
         END IF

   END DO

!------------------------

!  old bonds update

   do nb=1,noldbond

      oldbond(1:2,nb)=map(oldbond(1:2,nb))

   end do

!------------


   DO NS=1,NSYN

      IF(GLYTIP(1,NS)/=0)THEN
         GLYTIP(1,NS)=MAP(GLYTIP(1,NS))
      END IF

      IF(GLYTIP(2,NS)/=0)THEN
         GLYTIP(2,NS)=MAP(GLYTIP(2,NS))
      END IF

!--

      IF(GLYSEC(1,NS)/=0)THEN
         GLYSEC(1,NS)=MAP(GLYSEC(1,NS))
      END IF

      IF(GLYSEC(2,NS)/=0)THEN
         GLYSEC(2,NS)=MAP(GLYSEC(2,NS))
      END IF

!--

      DO J=1,3

         IF(TPPEP(1,J,NS)/=0)THEN
            TPPEP(1,J,NS)=MAP(TPPEP(1,J,NS))
         END IF

         IF(TPPEP(2,J,NS)/=0)THEN
            TPPEP(2,J,NS)=MAP(TPPEP(2,J,NS))
         END IF

      END DO

!--
      IF(EDPEP(1,NS)/=0)THEN
         EDPEP(1,NS)=MAP(EDPEP(1,NS))
      END IF

      IF(EDPEP(2,NS)/=0)THEN
         EDPEP(2,NS)=MAP(EDPEP(2,NS))
      END IF

      IF(EDHOLD(1,NS)/=0)THEN
         EDHOLD(1,NS)=MAP(EDHOLD(1,NS))
         EDHOLD(2,NS)=MAP(EDHOLD(2,NS))
      END IF

   END DO


!-----------------

   NNEW=0

   DO NL=1,NLOOP
      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF


      LENGTH=0

      DO JL=1,LOOPLEN(NL)

         NA=LOOP(JL,NL)

         IF(ATOR(NA)==-1)THEN
            CYCLE
         END IF

         IF(LENGTH>0)THEN
            IF(NA==LOOP(LENGTH,NL))THEN
               CYCLE
            END IF
         END IF


         LENGTH=LENGTH+1

         LOOP(LENGTH,NL)=NA

      END DO

      IF(LENGTH<4)THEN
         PRINT*,'REMOVING UNSTRUCTURED LOOP, LENGTH=',LENGTH
         NLOOPDEL=NLOOPDEL+1
         CYCLE
      END IF

      IF(LOOP(1,NL)==LOOP(LENGTH,NL))THEN
         LENGTH=LENGTH-1
      END IF

      NNEW=NNEW+1

      LOOPLEN(NNEW)=LENGTH

      LOOP(1:LENGTH,NNEW)=MAP(LOOP(1:LENGTH,NL))

      LOOPTYP(NNEW)=LOOPTYP(NL)

      IF(LOOPTYP(NL)==2)THEN
         DO NS=1,NSYN
            IF(SYNLOOP(NS)==NL)THEN
               SYNLOOP(NS)=NNEW
            END IF
         END DO
      END IF

   END DO


   IF(NNEW/=NLOOP-NLOOPDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED LOOPS',NNEW,NLOOP,NLOOPDEL
      STOP
   END IF

   NLOOP=NNEW
   NLOOPDEL=0


!-------------

   DNOR(natomold+1:NATOM)=DNORTEMP(natomold+1:NATOM)
   ATOR(natomold+1:NATOM)=ATORTEMP(natomold+1:NATOM)

   DEALLOCATE(MAP,DNORTEMP,ATORTEMP)

   END SUBROUTINE

!====================================================================

   SUBROUTINE SYNCOOR(NSYN,XSYN,X0,XCAP1,XCAP2)

   INTEGER NSYN,NS,J

   INTEGER, ALLOCATABLE, DIMENSION(:)::COUNTS

   DOUBLE PRECISION XCAP1,XCAP2,X0,LENGTH,R,DELT,PSUM

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::PROB,XSYN


   LENGTH=XCAP2-XCAP1


   IF(NSYN==1)THEN

      CALL RANDOM_NUMBER(R)

      X0=XCAP1+LENGTH*R

      RETURN

   END IF

   ALLOCATE(COUNTS(NSYN),PROB(NSYN))

   COUNTS=0

   DELT=LENGTH/NSYN

   DO NS=1,NSYN-1

      J=(XSYN(NS)-XCAP1)/DELT+1

      IF(J<1) J=1

      IF(J>NSYN) J=NSYN

      COUNTS(J)=COUNTS(J)+1

   END DO

   DO NS=1,NSYN

      IF(COUNTS(NS)==0)THEN

         PROB(NS)=10.0D0

      ELSE

         PROB(NS)=1.0D0/COUNTS(NS)

      END IF

   END DO

   PSUM=SUM(PROB(1:NSYN))

   PROB(1:NSYN)=PROB(1:NSYN)/PSUM

   CALL RANDOM_NUMBER(R)

   PSUM=0.0D0

   DO NS=1,NSYN

      PSUM=PSUM+PROB(NS)

      IF(PSUM>R)THEN

         J=NS

         EXIT

      END IF

   END DO

   CALL RANDOM_NUMBER(R)

   X0=XCAP1+(J-1)*DELT+R*DELT

   DEALLOCATE(COUNTS,PROB)

   END SUBROUTINE

!====================================================================

   SUBROUTINE REDISTRIBUTE(jsyndir,na_red,NATOM,X,Y,Z,NS,XGTASE,YGTASE,ZGTASE, &
               XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,SYNDIR)

   IMPLICIT NONE

   INTEGER,VALUE::NATOM,NS
   INTEGER NA,jsyndir,na_red
   INTEGER,DIMENSION(:),ALLOCATABLE:: SYNDIR

   DOUBLE PRECISION R,R2(2),XLIM1,XLIM2

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN)::X,Y,Z
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: XEDASE,YEDASE,ZEDASE!,PROB



   if(jsyndir/=0)then

      syndir(ns)=jsyndir

      na=na_red

      jsyndir=0

      goto 20

   end if



   CALL RANDOM_NUMBER(R)

   IF(R>0.5D0)THEN
      SYNDIR(NS)=1
   ELSE
      SYNDIR(NS)=-1
   END IF

   jsyndir=-syndir(ns)


      xlim1=-2.0d0

      xlim2=-xlim1


10 CALL RANDOM_NUMBER(R)

   NA=NATOM*R+1

   IF(X(NA)<XLIM1.OR.X(NA)>XLIM2)THEN
      GOTO 10
   END IF

   na_red=na

20 CALL RANDOM_NUMBER(R2)

   XEDASE(NS)=X(NA)
   YEDASE(NS)=Y(NA)+R2(1)-0.5D0
   ZEDASE(NS)=Z(NA)+R2(2)-0.5D0




11 CALL RANDOM_NUMBER(R2)

   XGTASE(1,NS)=XEDASE(NS)+0.5D0*SYNDIR(NS)
   YGTASE(1,NS)=YEDASE(NS)+R2(1)-0.5D0
   ZGTASE(1,NS)=ZEDASE(NS)+R2(2)-0.5D0

   CALL RANDOM_NUMBER(R2)

   XGTASE(2,NS)=XEDASE(NS)-0.5D0*SYNDIR(NS)
   YGTASE(2,NS)=YEDASE(NS)+R2(1)-0.5D0
   ZGTASE(2,NS)=ZEDASE(NS)+R2(2)-0.5D0

   CALL RANDOM_NUMBER(R2)

   XTPASE(1,NS)=XEDASE(NS)+1.0D0*SYNDIR(NS)
   YTPASE(1,NS)=YEDASE(NS)+R2(1)-0.5D0
   ZTPASE(1,NS)=ZEDASE(NS)+R2(2)-0.5D0

   CALL RANDOM_NUMBER(R2)

   XTPASE(2,NS)=XEDASE(NS)
   YTPASE(2,NS)=YEDASE(NS)+R2(1)-0.5D0
   ZTPASE(2,NS)=ZEDASE(NS)+R2(2)-0.5D0

   CALL RANDOM_NUMBER(R2)

   XTPASE(3,NS)=XEDASE(NS)-1.0D0*SYNDIR(NS)
   YTPASE(3,NS)=YEDASE(NS)+R2(1)-0.5D0
   ZTPASE(3,NS)=ZEDASE(NS)+R2(2)-0.5D0


   PRINT*,'COMPLEX',NS, 'RE-DISTRIBUTED AT'

1  FORMAT(I2,3(2X,F7.1))

   WRITE(*,1)NS,XEDASE(NS),YEDASE(NS),ZEDASE(NS)


   END SUBROUTINE


!==========================================

   SUBROUTINE INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

   IMPLICIT NONE

   INTEGER CHECK

   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3
   DOUBLE PRECISION M11,M12,M13,M21,M22,M23,M31,M32,M33
   DOUBLE PRECISION I11,I12,I13,I21,I22,I23,I31,I32,I33
   DOUBLE PRECISION DET,INVDET,T,U,V,XB,YB,ZB

! -- FORM THE MATRIX FROM THE LINE THROUGH L1 L2 AND PLANE THROUGH P1 P2 P3:

   M11=XL1-XL2
   M21=YL1-YL2
   M31=ZL1-ZL2

   M12=XP2-XP1
   M22=YP2-YP1
   M32=ZP2-ZP1

   M13=XP3-XP1
   M23=YP3-YP1
   M33=ZP3-ZP1

   DET=M11*M22*M33+M12*M23*M31+M13*M21*M32-M11*M23*M32-M12*M21*M33-M13*M22*M31

   INVDET=1.0D0/DET

!--- MATRIX INVERSE:

   I11=(M22*M33-M23*M32)*INVDET
   I12=(M13*M32-M12*M33)*INVDET
   I13=(M12*M23-M13*M22)*INVDET

   I21=(M23*M31-M21*M33)*INVDET
   I22=(M11*M33-M13*M31)*INVDET
   I23=(M13*M21-M11*M23)*INVDET

   I31=(M21*M32-M22*M31)*INVDET
   I32=(M12*M31-M11*M32)*INVDET
   I33=(M11*M22-M12*M21)*INVDET

! --- BASE VECTOR:
   XB=XL1-XP1
   YB=YL1-YP1
   ZB=ZL1-ZP1

! -- INTERSECTION IS REPRESENTED BY (T U V):

   T=I11*XB+I12*YB+I13*ZB

   U=I21*XB+I22*YB+I23*ZB

   V=I31*XB+I32*YB+I33*ZB

! --- INTERSECTION OCCURS IF T =(0,1); U,V = (0,1); U+V = (0,1)
! --- FIRST ASSUME CHECK = 1:

   CHECK = 1

   IF(T<0.0D0.OR.T>1.0D0)THEN
      CHECK=0
   END IF

   IF(U<0.0D0.OR.U>1.0D0)THEN
      CHECK=0
   END IF

   IF(V<0.0D0.OR.V>1.0D0)THEN
      CHECK=0
   END IF

   IF(U+V<0.0D0.OR.U+V>1.0D0)THEN
      CHECK=0
   END IF


   END SUBROUTINE

!==================================================================

   SUBROUTINE CALSURF(NSYN,SYNPG,PGID,PGLEN,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE, &
                      XEDASE,YEDASE,ZEDASE,NATOM,X,Y,Z,EDRAD,BIRAD1,BIRAD2,MORAD)


   IMPLICIT NONE

   INTEGER,VALUE:: NSYN,NATOM
   INTEGER NS,N,NA,NCOUNT,J,IPG
   INTEGER, ALLOCATABLE, DIMENSION(:),INTENT(IN):: PGLEN
   INTEGER, ALLOCATABLE, DIMENSION(:)::MARK
   INTEGER, ALLOCATABLE, DIMENSION(:,:),INTENT(IN):: SYNPG,PGID
   DOUBLE PRECISION DIST,Y0,Z0
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::X,Y,Z,XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::YCENBI1,ZCENBI1,YCENBI2,ZCENBI2
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::YCENMON,ZCENMON,YCENEDO,ZCENEDO

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::EDRAD,BIRAD1,BIRAD2,MORAD

   ALLOCATE(MARK(NATOM))
   MARK=0

   DO NS=1,NSYN

      IF(SYNPG(1,NS)>0)THEN
         IPG=SYNPG(1,NS)

         MARK(PGID(1:PGLEN(IPG),IPG))=1
      END IF

      IF(SYNPG(2,NS)>0)THEN
         IPG=SYNPG(2,NS)

         MARK(PGID(1:PGLEN(IPG),IPG))=1
      END IF

   END DO



   DO NS=1,NSYN

      DIST=2.0D0

!      Y0=YCENEDO(NS)

!      Z0=ZCENEDO(NS)


      DO J=1,100

!        for endopeptidase:

         DIST=DIST+1.0D0

         if(dist>50.0d0) print*,'warning: complex in big hole > 50 nm in radius'

         NCOUNT=0

         EDRAD(NS)=0.0D0

         DO N=1,NATOM

            IF(ABS(X(N)-XEDASE(NS))>DIST.OR.ABS(Y(N)-YEDASE(NS))>DIST.OR.ABS(Z(N)-ZEDASE(NS))>DIST.OR.MARK(N)==1)THEN
               CYCLE
            END IF

            NCOUNT=NCOUNT+1

            EDRAD(NS)=EDRAD(NS)+SQRT(Y(N)**2+Z(N)**2)

         END DO

         IF(NCOUNT>=4)THEN
            EDRAD(NS)=EDRAD(NS)/NCOUNT
            EXIT
         END IF

      END DO

      IF(NCOUNT<4)THEN
         PRINT*,'COULD NOT CALCULATE EDRAD FOR SYNCOMP',NS
         PRINT*,XEDASE(NS),YEDASE(NS),ZEDASE(NS)
         STOP
      END IF

!     for bifunctional enzyme 1:

      DIST=2.0D0

!      Y0=YCENBI1(NS)

!      Z0=ZCENBI1(NS)


      DO J=1,100

         DIST=DIST+1.0D0

         if(dist>50.0d0) print*,'warning: complex in big hole > 50 nm in radius'

         NCOUNT=0

         BIRAD1(NS)=0.0D0


         DO N=1,NATOM

            IF(ABS(X(N)-XGTASE(1,NS))>DIST.OR.ABS(Y(N)-YGTASE(1,NS))>DIST.OR.ABS(Z(N)-ZGTASE(1,NS))>DIST.OR.MARK(N)==1) CYCLE

            NCOUNT=NCOUNT+1

            BIRAD1(NS)=BIRAD1(NS)+SQRT(Y(N)**2+Z(N)**2)

         END DO

         IF(NCOUNT>=4)THEN
            BIRAD1(NS)=BIRAD1(NS)/NCOUNT
            EXIT
         END IF

      END DO

      IF(NCOUNT<4)THEN
         PRINT*,'COULD NOT CALCULATE BIRAD1 FOR SYNCOMP',NS
         PRINT*,XGTASE(1,NS),YGTASE(1,NS),ZGTASE(1,NS)
         STOP
      END IF


!     for bifunctional enzyme 2:

      DIST=2.0D0

!      Y0=YCENBI2(NS)

!      Z0=ZCENBI2(NS)

      DO J=1,100

         DIST=DIST+1.0D0

         if(dist>50.0d0) print*,'warning: complex in big hole > 50 nm in radius'

         NCOUNT=0

         BIRAD2(NS)=0.0D0

         DO N=1,NATOM

            IF(ABS(X(N)-XGTASE(2,NS))>DIST.OR.ABS(Y(N)-YGTASE(2,NS))>DIST.OR.ABS(Z(N)-ZGTASE(2,NS))>DIST.OR.MARK(N)==1)THEN
               CYCLE
            END IF

            NCOUNT=NCOUNT+1

            BIRAD2(NS)=BIRAD2(NS)+SQRT(Y(N)**2+Z(N)**2)

         END DO

         IF(NCOUNT>=4)THEN
            BIRAD2(NS)=BIRAD2(NS)/NCOUNT
            EXIT
         END IF

      END DO


      IF(NCOUNT<4)THEN
         PRINT*,'COULD NOT CALCULATE BIRAD2 FOR SYNCOMP',NS
         PRINT*,XGTASE(2,NS),YGTASE(2,NS),ZGTASE(2,NS)
         STOP
      END IF


!     for monofuntional transpeptidase:

      DIST=2.0D0

!      Y0=YCENMON(NS)

!      Z0=ZCENMON(NS)

      DO J=1,20

         DIST=DIST+1.0D0

         if(dist>50.0d0) print*,'warning: complex in big hole > 50 nm in radius'

         NCOUNT=0

         MORAD(NS)=0.0D0

         DO N=1,NATOM

            IF(ABS(X(N)-XTPASE(3,NS))>DIST.OR.ABS(Y(N)-YTPASE(3,NS))>DIST.OR.ABS(Z(N)-ZTPASE(3,NS))>DIST.OR.MARK(N)==1)THEN
               CYCLE
            END IF

            NCOUNT=NCOUNT+1

            MORAD(NS)=MORAD(NS)+SQRT(Y(N)**2+Z(N)**2)

         END DO

         IF(NCOUNT>=4)THEN
            MORAD(NS)=MORAD(NS)/NCOUNT
            EXIT
         END IF

      END DO

      IF(NCOUNT<4)THEN
         PRINT*,'COULD NOT CALCULATE MORAD FOR SYNCOMP',NS
         PRINT*,XTPASE(3,NS),YTPASE(3,NS),ZTPASE(3,NS)
         STOP
      END IF

   END DO


   DEALLOCATE(MARK)


   END SUBROUTINE

!==========================================

   SUBROUTINE SURFLINK(NSYN,loop,looplen,synloop,y,z,YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE, &
              FYGTASE,FZGTASE,FYTPASE,FZTPASE,FYEDASE,FZEDASE,KSUR,pthick)

   IMPLICIT NONE

   INTEGER,VALUE:: NSYN
   INTEGER N,iloop,j,na

   INTEGER,ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::loop
   INTEGER,ALLOCATABLE, DIMENSION(:),INTENT(IN)::looplen,synloop

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::FYEDASE,FZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::YTPASE,ZTPASE,YGTASE,ZGTASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::y,z
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::YCENBI1,ZCENBI1,YCENBI2,ZCENBI2
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::YCENMON,ZCENMON,YCENEDO,ZCENEDO

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::FYTPASE,FZTPASE,FYGTASE,FZGTASE
   DOUBLE PRECISION,VALUE:: KSUR,pthick
   DOUBLE PRECISION RAD,F,radius!,DY,DZ


   DO N=1,NSYN

      iloop=synloop(n)

      if(iloop==0) cycle

      radius=0.0d0

      do j=1,looplen(iloop)

         na=loop(j,iloop)

         radius=radius+sqrt(y(na)*y(na)+z(na)*z(na))

      end do

      radius=radius/looplen(iloop)



!     CONSTRAINT ON TRANSGLYCOSYLASES:

!         DY=YGTASE(1,N)-YCENBI1(N)
!         DZ=ZGTASE(1,N)-ZCENBI1(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         RAD=SQRT(YGTASE(1,N)**2+ZGTASE(1,N)**2)

         if(rad>radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad<radius-pthick)then
            f=-KSUR*(RAD-radius+pthick)/RAD
         else
            f=0.0d0
         end if

         FYGTASE(1,N)=FYGTASE(1,N)+F*YGTASE(1,N)!GTASE(1,N)/RAD

         FZGTASE(1,N)=FZGTASE(1,N)+F*ZGTASE(1,N)!GTASE(1,N)/RAD


!         DY=YGTASE(2,N)-YCENBI2(N)
!         DZ=ZGTASE(2,N)-ZCENBI2(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         RAD=SQRT(YGTASE(2,N)**2+ZGTASE(2,N)**2)

         if(rad>radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad<radius-pthick)then
            f=-KSUR*(RAD-radius+pthick)/RAD
         else
            f=0.0d0
         end if

!         F=-KSUR*(RAD-radius)/RAD

         FYGTASE(2,N)=FYGTASE(2,N)+F*YGTASE(2,N)!GTASE(2,N)/RAD

         FZGTASE(2,N)=FZGTASE(2,N)+F*ZGTASE(2,N)!GTASE(2,N)/RAD

!     CONSTRAINT ON TRANSPEPTIDASES

!         DY=YTPASE(1,N)-YCENBI1(N)
!         DZ=ZTPASE(1,N)-ZCENBI1(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         RAD=SQRT(YTPASE(1,N)**2+ZTPASE(1,N)**2)

         if(rad>radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad<radius-pthick)then
            f=-KSUR*(RAD-radius+pthick)/RAD
         else
            f=0.0d0
         end if

!         F=-KSUR*(RAD-radius)/RAD

         FYTPASE(1,N)=FYTPASE(1,N)+F*YTPASE(1,N)!TPASE(1,N)/RAD

         FZTPASE(1,N)=FZTPASE(1,N)+F*ZTPASE(1,N)!TPASE(1,N)/RAD


         RAD=SQRT(YTPASE(2,N)**2+ZTPASE(2,N)**2)

!         DY=YTPASE(2,N)-YCENBI2(N)
!         DZ=ZTPASE(2,N)-ZCENBI2(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         if(rad>radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad<radius-pthick)then
            f=-KSUR*(RAD-radius+pthick)/RAD
         else
            f=0.0d0
         end if

!         F=-KSUR*(RAD-radius)/RAD

         FYTPASE(2,N)=FYTPASE(2,N)+F*YTPASE(2,N)!TPASE(2,N)/RAD

         FZTPASE(2,N)=FZTPASE(2,N)+F*ZTPASE(2,N)!TPASE(2,N)/RAD


         RAD=SQRT(YTPASE(3,N)**2+ZTPASE(3,N)**2)

!         DY=YTPASE(3,N)-YCENMON(N)
!         DZ=ZTPASE(3,N)-ZCENMON(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         if(rad>radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad<radius-pthick)then
            f=-KSUR*(RAD-radius+pthick)/RAD
         else
            f=0.0d0
         end if

!         F=-KSUR*(RAD-radius)/RAD

         FYTPASE(3,N)=FYTPASE(3,N)+F*YTPASE(3,N)!TPASE(3,N)/RAD

         FZTPASE(3,N)=FZTPASE(3,N)+F*ZTPASE(3,N)!TPASE(3,N)/RAD

!     CONSTRAINT ON ENDOPEPTIDASE:

         RAD=SQRT(YEDASE(N)**2+ZEDASE(N)**2)

!         DY=YEDASE(N)-YCENEDO(N)
!         DZ=ZEDASE(N)-ZCENEDO(N)

!         RAD=SQRT(DY*DY+DZ*DZ)

         if(rad<radius)then
            F=-KSUR*(RAD-radius)/RAD
         elseif(rad>radius+pthick)then
            f=-KSUR*(RAD-radius-pthick)/RAD
         else
            f=0.0d0
         end if

!         F=-KSUR*(RAD-radius)/RAD

         FYEDASE(N)=FYEDASE(N)+F*YEDASE(N)!EDASE(N)/RAD

         FZEDASE(N)=FZEDASE(N)+F*ZEDASE(N)!EDASE(N)/RAD

   END DO

   END SUBROUTINE

!==========================================
   subroutine axis(caplen,kaxis,y,z,fy,fz)

   IMPLICIT NONE

   integer,value::caplen

!   integer,allocatable,dimension(:,:),intent(in)::capid

   real(kind=8),value::kaxis

   real(kind=8)::f

   real(kind=8),allocatable,dimension(:),intent(in)::y,z

   real(kind=8),allocatable,dimension(:)::fy,fz


   f=-kaxis*SUM(Y(1:caplen))/caplen

   fy(1:caplen)=fy(1:caplen)+f

   f=-kaxis*SUM(Z(1:caplen))/caplen

   fz(1:caplen)=fz(1:caplen)+f


   f=-kaxis*SUM(Y(1+caplen:caplen*2))/caplen

   fy(1+caplen:caplen*2)=fy(1+caplen:caplen*2)+f

   f=-kaxis*SUM(Z(1+caplen:caplen*2))/caplen

   fz(1+caplen:caplen*2)=fz(1+caplen:caplen*2)+f


   END SUBROUTINE

!==========================================

   SUBROUTINE temlinkforce(ntemlink,temlink,KTHETA,THET0,BETA,INVDELTA,DELTA,l_g,k_g,X,Y,Z,FX,FY,FZ)

   implicit none

   INTEGER,VALUE::ntemlink
   integer n,n1,n2,n3,n4

   INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::temlink


   REAL(KIND=8),VALUE::KTHETA,THET0,BETA,INVDELTA,DELTA,l_g,k_g

   REAL(KIND=8)::DX1,DY1,DZ1,DX3,DY3,DZ3,INVDIST1,INVDIST3,COS_T0,COS_T,THET,F0
   REAL(KIND=8)::dx,dy,dz,DFX1,DFY1,DFZ1,DFX3,DFY3,DFZ3,DREP,dist,f

   REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(IN)::X,Y,Z
   REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::FX,FY,FZ


   do n=1,ntemlink

      n1=temlink(1,n)
      n2=temlink(2,n)
      n3=temlink(3,n)
      n4=temlink(4,n)

!     1st angle term

      dx1=x(n1)-x(n2)
      dy1=y(n1)-y(n2)
      dz1=z(n1)-z(n2)

      INVDIST1=1.0d0/sqrt(dx1**2+dy1**2+dz1**2)


      dx3=x(n3)-x(n2)
      dy3=y(n3)-y(n2)
      dz3=z(n3)-z(n2)

      dist=sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      invdist3=1.0d0/dist


      COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)*invDIST1*invDIST3

      THET=ACOS(COS_T0*(1.0D0-BETA))

      F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA

!---- FORCE ON N1 ALONG X:

      DX=DX1+DELTA

      DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

      COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP*invDIST3


      DFX1=F0*(COS_T-COS_T0)

!---- FORCE ON N1 ALONG Y:

      DY=DY1+DELTA

      DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

      COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP*invDIST3


      DFY1=F0*(COS_T-COS_T0)

!---- FORCE ON N1 ALONG Z:

      DZ=DZ1+DELTA

      DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP*invDIST3


      DFZ1=F0*(COS_T-COS_T0)

!---- FORCE ON N3 ALONG X:

      DX=DX3+DELTA

      DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

      COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)*invDIST1/DREP

      DFX3=F0*(COS_T-COS_T0)

!---- FORCE ON N3 ALONG Y:

      DY=DY3+DELTA

      DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

      COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)*invDIST1/DREP

      DFY3=F0*(COS_T-COS_T0)

!---- FORCE ON N3 ALONG Z:

      DZ=DZ3+DELTA

      DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)*invDIST1/DREP

      DFZ3=F0*(COS_T-COS_T0)


!     FORCES ON N1:

      FX(N1)=FX(N1)+DFX1

      FY(N1)=FY(N1)+DFY1

      FZ(N1)=FZ(N1)+DFZ1

!     FORCES ON N2:

      FX(N2)=FX(N2)-DFX1-DFX3

      FY(N2)=FY(N2)-DFY1-DFY3

      FZ(N2)=FZ(N2)-DFZ1-DFZ3

!     FORCES ON N3:

      FX(N3)=FX(N3)+DFX3

      FY(N3)=FY(N3)+DFY3

      FZ(N3)=FZ(N3)+DFZ3


!     2nd angle term

      dx1=-dx3
      dy1=-dy3
      dz1=-dz3

      invdist1=invdist3

      dx3=x(n4)-x(n3)
      dy3=y(n4)-y(n3)
      dz3=z(n4)-z(n3)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)


      COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)*invDIST1*invDIST3

      THET=ACOS(COS_T0*(1.0D0-BETA))

      F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA

!---- FORCE ON N2 ALONG X:

      DX=DX1+DELTA

      DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

      COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP*invDIST3


      DFX1=F0*(COS_T-COS_T0)

!---- FORCE ON N2 ALONG Y:

      DY=DY1+DELTA

      DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

      COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP*invDIST3


      DFY1=F0*(COS_T-COS_T0)


!---- FORCE ON N2 ALONG Z:

      DZ=DZ1+DELTA

      DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP*invDIST3


      DFZ1=F0*(COS_T-COS_T0)

!---- FORCE ON N4 ALONG X:

      DX=DX3+DELTA

      DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

      COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)*invDIST1/DREP

      DFX3=F0*(COS_T-COS_T0)


!---- FORCE ON N4 ALONG Y:

      DY=DY3+DELTA

      DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

      COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)*invDIST1/DREP

      DFY3=F0*(COS_T-COS_T0)


!---- FORCE ON N4 ALONG Z:

      DZ=DZ3+DELTA

      DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)*invDIST1/DREP

      DFZ3=F0*(COS_T-COS_T0)

!     FORCES ON N2:

      FX(N2)=FX(N2)+DFX1

      FY(N2)=FY(N2)+DFY1

      FZ(N2)=FZ(N2)+DFZ1

!     FORCES ON N3:

      FX(N3)=FX(N3)-DFX1-DFX3

      FY(N3)=FY(N3)-DFY1-DFY3

      FZ(N3)=FZ(N3)-DFZ1-DFZ3

!     FORCES ON N4:

      FX(N4)=FX(N4)+DFX3

      FY(N4)=FY(N4)+DFY3

      FZ(N4)=FZ(N4)+DFZ3


!     bond term

      if(dist<l_g)then

         F=K_G*(DIST-L_G)/DIST

         DFX1=F*DX1
         DFY1=F*DY1
         DFZ1=F*DZ1

         fx(n2)=fx(n2)-dfx1
         fy(n2)=fy(n2)-dfy1
         fz(n2)=fz(n2)-dfz1

         fx(n3)=fx(n3)+dfx1
         fy(n3)=fy(n3)+dfy1
         fz(n3)=fz(n3)+dfz1

      end if

   end do


   END SUBROUTINE

!==========================================

   SUBROUTINE constrict(natom,dnor,ator,fconstrict,sigma,x,y,z,fy,fz,ftotal)

   implicit none

   integer,value::natom
   integer n
   integer,allocatable,dimension(:),intent(in)::dnor,ator

   real(kind=8),value::fconstrict,sigma
   real(kind=8)::ftotal
   real(kind=8)::rad,f,sigma2
   real(kind=8),allocatable,dimension(:),intent(in)::x,y,z
   real(kind=8),allocatable,dimension(:)::fy,fz

   sigma2=0.25d0*sigma*sigma

   ftotal=0.0d0

!$omp parallel &
!$omp default(none) &
!$omp private(n,rad,f) &
!$omp shared(natom,dnor,ator,x,y,z,fy,fz,fconstrict,sigma2) &
!$omp reduction(+:ftotal)

!$omp do schedule(guided,32)

   do n=1,natom

      if(dnor(n)/=0.and.ator(n)/=0) cycle

      rad=sqrt(y(n)*y(n)+z(n)*z(n))

      f=fconstrict*exp(-x(n)*x(n)/sigma2)/rad

      fy(n)=fy(n)-f*y(n)

      fz(n)=fz(n)-f*z(n)

      ftotal=ftotal+f*rad

   end do

!$omp end do
!$omp end parallel

   END SUBROUTINE


!==========================================

   SUBROUTINE SURFCONSTR(natom,NPG,newPGLEN,DNOR,ATOR,newPGID,KSUR,X,Y,Z,FY,FZ,pthick)

   IMPLICIT NONE

   INTEGER,VALUE::NPG,natom

   INTEGER N,J,NA,N1,N2,NCOUNT,NS,JS,NB,jcycle

   INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::newPGLEN,DNOR,ATOR
   INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::newPGID

   REAL(KIND=8),VALUE::KSUR,pthick
   REAL(KIND=8)::RAD,DX,DY,DZ,D2,RAD_NA,F
   REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(IN)::X,Y,Z
   REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::FY,FZ

   DO N=1,NPG

!      IF(PGTYP(N)==0)THEN
!         CYCLE
!      END IF


      IF(newPGLEN(N)<3) cycle

      DO J=2,newPGLEN(N)-1

         NA=newPGID(J,N)

         N1=newPGID(J-1,N)

         N2=newPGID(J+1,N)

         IF(DNOR(NA)/=0.AND.ATOR(NA)/=0.AND.DNOR(N1)/=0.AND.ATOR(N1)/=0.AND.DNOR(N2)/=0.AND.ATOR(N2)/=0)THEN

            RAD=0.0D0

            NCOUNT=0

            do nb=1,natom

               jcycle=0

               do js=1,newPGLEN(N)

                  if(newPGID(JS,N)==nb)then
                     jcycle=1
                     exit
                  end if

               end do

               if(jcycle==1) cycle

!            DO NS=1,NPG

!               IF(NS==N)THEN
!                  CYCLE
!               END IF


!               DO JS=1,PGLEN(NS)

!                  NB=PGID(JS,NS)



               IF(DNOR(NB)/=0.AND.ATOR(NB)/=0)THEN
                  CYCLE
               END IF

               DX=X(NA)-X(NB)

               DY=Y(NA)-Y(NB)

               DZ=Z(NA)-Z(NB)

               D2=DX*DX+DY*DY+DZ*DZ

               IF(D2<36.0D0)THEN

                  NCOUNT=NCOUNT+1

                  RAD=RAD+SQRT(Y(NB)*Y(NB)+Z(NB)*Z(NB))

               END IF

            END DO

!            END DO

            IF(NCOUNT>0)THEN

               RAD=RAD/NCOUNT

               RAD_NA=SQRT(Y(NA)*Y(NA)+Z(NA)*Z(NA))

               if(rad_na>rad+pthick)then

                  F=-KSUR*(RAD_NA-RAD-pthick)/RAD_NA

               elseif(rad_na<rad-pthick)then

                  F=-KSUR*(RAD_NA-RAD+pthick)/RAD_NA

               else

                  cycle

               end if

               FY(NA)=FY(NA)+F*Y(NA)

               FZ(NA)=FZ(NA)+F*Z(NA)

            END IF

         END IF

      END DO

   END DO

   END SUBROUTINE


!==========================================

   SUBROUTINE SETRAND(RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRAND,NRAND)

   IMPLICIT NONE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::R1,R2,rx1,rx2,ry1,ry2,rz1,rz2
   INTEGER NRAND,JRAND,NS,K1,K2,K4,NHALF
   DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793239D0

!   INTEGER(KIND=8)::count_rate,cr,t0,t1,t2,t3,t4,t5
!   REAL(KIND=8)::rate

   integer nfull

!   real(kind=8)::arg,thet

!print*,'setrand'

!   CALL system_clock(count_rate=cr)
!   rate = REAL(cr)


   JRAND=0

   K1=70
   K2=60
   K4=50


   nfull=3*nrand/2
   nhalf=nrand/2

!   call system_clock(t2)

   allocate(r1(nfull),r2(nfull))

!  for Gtases and Edases

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RXGTASE(1:nrand)=K1*SQRT(-2*LOG(R1(1:nrand)))*SIN(2*PI*R2(1:nrand))
      RXGTASE(1+nrand:2*NRAND)=K1*SQRT(-2*LOG(R1(1:nrand)))*COS(2*PI*R2(1:nrand))

      RXEDASE(1:NHALF)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*SIN(2*PI*R2(1+nrand:nfull))
      RXEDASE(1+NHALF:NRAND)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*cos(2*PI*R2(1+nrand:nfull))

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RyGTASE(1:nrand)=K1*SQRT(-2*LOG(R1(1:nrand)))*SIN(2*PI*R2(1:nrand))
      RyGTASE(1+nrand:2*NRAND)=K1*SQRT(-2*LOG(R1(1:nrand)))*COS(2*PI*R2(1:nrand))

      RyEDASE(1:NHALF)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*SIN(2*PI*R2(1+nrand:nfull))
      RyEDASE(1+NHALF:NRAND)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*cos(2*PI*R2(1+nrand:nfull))

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RzGTASE(1:nrand)=K1*SQRT(-2*LOG(R1(1:nrand)))*SIN(2*PI*R2(1:nrand))
      RzGTASE(1+nrand:2*NRAND)=K1*SQRT(-2*LOG(R1(1:nrand)))*COS(2*PI*R2(1:nrand))

      RzEDASE(1:NHALF)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*SIN(2*PI*R2(1+nrand:nfull))
      RzEDASE(1+NHALF:NRAND)=K4*SQRT(-2*LOG(R1(1+nrand:nfull)))*cos(2*PI*R2(1+nrand:nfull))


!  for Tpases

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RXTPASE(1:nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*SIN(2*PI*R2(1:nfull))
      RXTPASE(1+nfull:2*nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*COS(2*PI*R2(1:nfull))

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RyTPASE(1:nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*SIN(2*PI*R2(1:nfull))
      RyTPASE(1+nfull:2*nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*COS(2*PI*R2(1:nfull))

      CALL RANDOM_NUMBER(R1)
      CALL RANDOM_NUMBER(R2)

      RzTPASE(1:nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*SIN(2*PI*R2(1:nfull))
      RzTPASE(1+nfull:2*nfull)=K2*SQRT(-2*LOG(R1(1:nfull)))*COS(2*PI*R2(1:nfull))

   DEALLOCATE(R1,R2)

!   call system_clock(t3)


!-----------------------------------------


!   call system_clock(t4)

!   allocate(rx1(nfull),rx2(nfull),ry1(nfull),ry2(nfull),rz1(nfull),rz2(nfull))


!  for Gtases and Edases

!   CALL RANDOM_NUMBER(Rx1)
!   CALL RANDOM_NUMBER(Rx2)
!   CALL RANDOM_NUMBER(Ry1)
!   CALL RANDOM_NUMBER(Ry2)
!   CALL RANDOM_NUMBER(Rz1)
!   CALL RANDOM_NUMBER(Rz2)



!omp parallel &
!omp default(none) &

!omp private(n,arg,thet,j) &
!omp shared(nrand,k1,rx1,rx2,ry1,ry2,rz1,rz2,rxgtase,rygtase,rzgtase,nfull,rxedase,ryedase,rzedase)

!omp do schedule(guided,64)

!   do n=1,nrand

!      arg=k1*sqrt(-2*log(rx1(n)))
!      thet=2*pi*rx2(n)

!      rxgtase(2*n-1)=arg*sin(thet); rxgtase(2*n)=arg*cos(thet)

!      arg=k1*sqrt(-2*log(ry1(n)))
!      thet=2*pi*ry2(n)

!      rygtase(2*n-1)=arg*sin(thet); rygtase(2*n)=arg*cos(thet)

!      arg=k1*sqrt(-2*log(rz1(n)))
!      thet=2*pi*rz2(n)

!      rzgtase(2*n-1)=arg*sin(thet); rzgtase(2*n)=arg*cos(thet)

!   end do

!omp end do nowait

!omp do schedule(guided,64)

!   do n=nrand+1,nfull

!      arg=k1*sqrt(-2*log(rx1(n)))
!      thet=2*pi*rx2(n)

!      j=n-nrand

!      rxedase(2*j-1)=arg*sin(thet); rxedase(2*j)=arg*cos(thet)

!      arg=k1*sqrt(-2*log(ry1(n)))
!      thet=2*pi*ry2(n)

!      ryedase(2*j-1)=arg*sin(thet); ryedase(2*j)=arg*cos(thet)


!      arg=k1*sqrt(-2*log(rz1(n)))
!      thet=2*pi*rz2(n)

!      rzedase(2*j-1)=arg*sin(thet); rzedase(2*j)=arg*cos(thet)


!   end do

!omp end do nowait
!omp end parallel




!  for Tpases

!   CALL RANDOM_NUMBER(Rx1)
!   CALL RANDOM_NUMBER(Rx2)
!   CALL RANDOM_NUMBER(Ry1)
!   CALL RANDOM_NUMBER(Ry2)
!   CALL RANDOM_NUMBER(Rz1)
!   CALL RANDOM_NUMBER(Rz2)


!omp parallel &
!omp default(none) &

!omp private(n,arg,thet) &
!omp shared(k1,rx1,rx2,ry1,ry2,rz1,rz2,nfull,rxtpase,rytpase,rztpase)

!omp do schedule(guided,64)

!   do n=1,nfull

!      arg=k1*sqrt(-2*log(rx1(n)))
!      thet=2*pi*rx2(n)

!      rxtpase(2*n-1)=arg*sin(thet); rxtpase(2*n)=arg*cos(thet)

!      arg=k1*sqrt(-2*log(ry1(n)))
!      thet=2*pi*ry2(n)

!      rytpase(2*n-1)=arg*sin(thet); rytpase(2*n)=arg*cos(thet)

!      arg=k1*sqrt(-2*log(rz1(n)))
!      thet=2*pi*rz2(n)

!      rztpase(2*n-1)=arg*sin(thet); rztpase(2*n)=arg*cos(thet)


!   end do


!omp end do nowait
!omp end parallel



!   call system_clock(t5)

!print*,(t1-t0)/rate,(t3-t2)/rate,(t5-t4)/rate


!stop
  END SUBROUTINE

!=========================================================================

   SUBROUTINE RANDFORCES(NSYN,FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE, &
                         RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRAND, &
                         YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE)


   IMPLICIT NONE

   INTEGER NSYN,JRAND,NS,J0

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FXEDASE,FYEDASE,FZEDASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::FXTPASE,FYTPASE,FZTPASE,FXGTASE,FYGTASE,FZGTASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::RZTPASE,RXEDASE,RYEDASE,RZEDASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::YGTASE,ZGTASE,YTPASE,ZTPASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN)::YEDASE,ZEDASE
   DOUBLE PRECISION ARG,DIST

!   JRAND=JRAND+1
   J0=2*JRAND


   FXGTASE(1,1:NSYN)=RXGTASE(J0+1:NSYN+J0)


   FYGTASE(1,1:NSYN)=RYGTASE(J0+1:NSYN+J0)


   FZGTASE(1,1:NSYN)=RZGTASE(J0+1:NSYN+J0)

!--------------------------------


   FXGTASE(2,1:NSYN)=RXGTASE(J0+NSYN+1:2*NSYN+J0)


   FYGTASE(2,1:NSYN)=RYGTASE(J0+NSYN+1:2*NSYN+J0)


   FZGTASE(2,1:NSYN)=RZGTASE(J0+NSYN+1:2*NSYN+J0)

   J0=3*JRAND


   FXTPASE(1,1:NSYN)=RXTPASE(J0+1:NSYN+J0)


   FYTPASE(1,1:NSYN)=RYTPASE(J0+1:NSYN+J0)


   FZTPASE(1,1:NSYN)=RZTPASE(J0+1:NSYN+J0)


   FXTPASE(2,1:NSYN)=RXTPASE(J0+NSYN+1:2*NSYN+J0)


   FYTPASE(2,1:NSYN)=RYTPASE(J0+NSYN+1:2*NSYN+J0)


   FZTPASE(2,1:NSYN)=RZTPASE(J0+NSYN+1:2*NSYN+J0)




   FXTPASE(3,1:NSYN)=RXTPASE(J0+2*NSYN+1:3*NSYN+J0)


   FYTPASE(3,1:NSYN)=RYTPASE(J0+2*NSYN+1:3*NSYN+J0)


   FZTPASE(3,1:NSYN)=RZTPASE(J0+2*NSYN+1:3*NSYN+J0)



   FXEDASE(1:NSYN)=RXEDASE(JRAND+1:NSYN+JRAND)


   FYEDASE(1:NSYN)=RYEDASE(JRAND+1:NSYN+JRAND)


   FZEDASE(1:NSYN)=RZEDASE(JRAND+1:NSYN+JRAND)



   JRAND=JRAND+NSYN

!---- REDUCE RANDOM FORCE ON THE PERPENDICULAR DIRECTION TO THE SURFACE 90%:

!OMP PARALLEL &
!OMP DEFAULT(NONE) &

!OMP PRIVATE(NS,ARG,DIST) &
!OMP SHARED(NSYN,FYGTASE,FZGTASE,YGTASE,ZGTASE,FYTPASE,FZTPASE,YTPASE,ZTPASE,FYEDASE,FZEDASE,YEDASE,ZEDASE)

!OMP DO 


   DO NS=1,NSYN

      DIST=YGTASE(1,NS)**2+ZGTASE(1,NS)**2


      ARG=(FYGTASE(1,NS)*YGTASE(1,NS)+FZGTASE(1,NS)*ZGTASE(1,NS))/DIST


      FYGTASE(1,NS)=FYGTASE(1,NS)-0.9D0*ARG*YGTASE(1,NS)!/DIST

      FZGTASE(1,NS)=FZGTASE(1,NS)-0.9D0*ARG*ZGTASE(1,NS)!/DIST

      DIST=YGTASE(2,NS)**2+ZGTASE(2,NS)**2

      ARG=(FYGTASE(2,NS)*YGTASE(2,NS)+FZGTASE(2,NS)*ZGTASE(2,NS))/DIST


      FYGTASE(2,NS)=FYGTASE(2,NS)-0.9D0*ARG*YGTASE(2,NS)!/DIST

      FZGTASE(2,NS)=FZGTASE(2,NS)-0.9D0*ARG*ZGTASE(2,NS)!/DIST

!--------

      DIST=YTPASE(1,NS)**2+ZTPASE(1,NS)**2

      ARG=(FYTPASE(1,NS)*YTPASE(1,NS)+FZTPASE(1,NS)*ZTPASE(1,NS))/DIST


      FYTPASE(1,NS)=FYTPASE(1,NS)-0.9D0*ARG*YTPASE(1,NS)!/DIST

      FZTPASE(1,NS)=FZTPASE(1,NS)-0.9D0*ARG*ZTPASE(1,NS)!/DIST


      DIST=YTPASE(2,NS)**2+ZTPASE(2,NS)**2

      ARG=(FYTPASE(2,NS)*YTPASE(2,NS)+FZTPASE(1,NS)*ZTPASE(2,NS))/DIST


      FYTPASE(2,NS)=FYTPASE(2,NS)-0.9D0*ARG*YTPASE(2,NS)!/DIST

      FZTPASE(2,NS)=FZTPASE(2,NS)-0.9D0*ARG*ZTPASE(2,NS)!/DIST


      DIST=YTPASE(3,NS)**2+ZTPASE(3,NS)**2

      ARG=(FYTPASE(3,NS)*YTPASE(3,NS)+FZTPASE(3,NS)*ZTPASE(3,NS))/DIST


      FYTPASE(3,NS)=FYTPASE(3,NS)-0.9D0*ARG*YTPASE(3,NS)!/DIST

      FZTPASE(3,NS)=FZTPASE(1,NS)-0.9D0*ARG*ZTPASE(3,NS)!/DIST

!---------

      DIST=YEDASE(NS)**2+ZEDASE(NS)**2

      ARG=(FYEDASE(NS)*YEDASE(NS)+FZEDASE(NS)*ZEDASE(NS))/DIST


      FYEDASE(NS)=FYEDASE(NS)-0.9D0*ARG*YEDASE(NS)!/DIST

      FZEDASE(NS)=FZEDASE(NS)-0.9D0*ARG*ZEDASE(NS)!/DIST

   END DO

!OMP END DO NOWAIT
!OMP END PARALLEL

   END SUBROUTINE

!=========================================================================

   SUBROUTINE SYNCOMP(NSYN,SYNDIR,edpep,tppep,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                      KPAIR,LPAIR,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED,KSIDE,LSIDE,DELTA,BETA,ktped,ltped, &
                      FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE)

   IMPLICIT NONE

   INTEGER NS
   INTEGER,VALUE::NSYN
   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::SYNDIR
   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::edpep
   INTEGER,DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN)::tppep

   DOUBLE PRECISION,VALUE:: KPAIR,LPAIR,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED,KSIDE,LSIDE,DELTA,BETA,ktped,ltped
   DOUBLE PRECISION DX,DY,DZ,DIST,FX,FY,FZ,F,invrad,ylead,zlead,proj,lgted_rep,ygtase0,zgtase0,l0

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN)::XEDASE,YEDASE,ZEDASE!,XLEAD,YLEAD,ZLEAD
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FXEDASE,FYEDASE,FZEDASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE


!OMP PARALLEL &
!OMP DEFAULT(NONE) &

!OMP PRIVATE(NS,DX,DY,DZ,FX,FY,FZ,F,DIST,INVDIST) &

!OMP SHARED(NSYN,SYNDIR,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE) &
!OMP SHARED(KPAIR,LPAIR,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED,KSIDE,LSIDE,DELTA,BETA,XLEAD,YLEAD,ZLEAD) &
!OMP SHARED(FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE)

!OMP DO 

   DO NS=1,NSYN

!---- constrain THE FIRST TPASE TO THE left of FIRST GTASE:

      DX=(XTPASE(1,NS)-XGTASE(1,NS))*syndir(ns)

      IF(DX<LSIDE)THEN

         F=SYNDIR(NS)*KSIDE*(DX-LSIDE)**2

         FXTPASE(1,NS)=FXTPASE(1,NS)+F

         FXGTASE(1,NS)=FXGTASE(1,NS)-F

      END IF

      if(dx>lgttp13)then

         F=SYNDIR(NS)*KSIDE*(DX-lgttp13)**2

         FXTPASE(1,NS)=FXTPASE(1,NS)-F

         FXGTASE(1,NS)=FXGTASE(1,NS)+F

      end if

!     pair THE FIRST TPASE TO THE FIRST GTASE on YZ:

!      DY=YTPASE(1,NS)-YGTASE(1,NS)

!      fy=-2*KPAIR*dy

!      FYGTASE(1,NS)=FYGTASE(1,NS)-FY
!      FYTPASE(1,NS)=FYTPASE(1,NS)+FY


!      DZ=ZTPASE(1,NS)-ZGTASE(1,NS)

!      fz=-2*KPAIR*dz

!      FZGTASE(1,NS)=FZGTASE(1,NS)-FZ
!      FZTPASE(1,NS)=FZTPASE(1,NS)+FZ




!     CONSTRAIN THE THIRD TPASE TO RIGHT SIDE OF THE SECOND GTASE:

      DX=(XGTASE(2,NS)-XTPASE(3,NS))*syndir(ns)

      IF(DX<LSIDE)THEN

         F=SYNDIR(NS)*KSIDE*(DX-LSIDE)**2

         FXGTASE(2,NS)=FXGTASE(2,NS)+F

         FXTPASE(3,NS)=FXTPASE(3,NS)-F

      END IF

      if(dx>lgttp13)then

         F=SYNDIR(NS)*KSIDE*(DX-lgttp13)**2

         FXGTASE(2,NS)=FXGTASE(2,NS)-F

         FXTPASE(3,NS)=FXTPASE(3,NS)+F

      end if

!     pair THE THIRD TPASE TO THE SECOND GTASE:

!      DY=YGTASE(2,NS)-YTPASE(3,NS)

!      fy=-2*KPAIR*dy

!      fygtase(2,ns)=fygtase(2,ns)+fy

!      fytpase(3,ns)=fytpase(3,ns)-fy

!      DZ=ZGTASE(2,NS)-ZTPASE(3,NS)

!      fz=-2*KPAIR*dz

!      fzgtase(2,ns)=fzgtase(2,ns)+fz

!      fztpase(3,ns)=fztpase(3,ns)-fz


!---- TETHER THE 1st TPASE TO THE 1st GTASE:


      DX=XGTASE(1,NS)-XTPASE(1,NS)
      DY=YGTASE(1,NS)-YTPASE(1,NS)
      DZ=ZGTASE(1,NS)-ZTPASE(1,NS)

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      if(tppep(1,1,ns)==0)then
         l0=0.5*LGTTP13
      else
         l0=LGTTP13
      end if

      IF(DIST>l0)THEN

         F=-KGTTP*(DIST-l0)/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXGTASE(1,NS)=FXGTASE(1,NS)+FX
         FYGTASE(1,NS)=FYGTASE(1,NS)+FY
         FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

         FXTPASE(1,NS)=FXTPASE(1,NS)-FX
         FYTPASE(1,NS)=FYTPASE(1,NS)-FY
         FZTPASE(1,NS)=FZTPASE(1,NS)-FZ

      END IF

!---- TETHER THE SECOND TPASE TO THE SECOND GTASE:

      DX=XGTASE(2,NS)-XTPASE(2,NS)
      DY=YGTASE(2,NS)-YTPASE(2,NS)
      DZ=ZGTASE(2,NS)-ZTPASE(2,NS)

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      IF(DIST>LGTTP2)THEN

         F=-KGTTP*(DIST-LGTTP2)/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXGTASE(2,NS)=FXGTASE(2,NS)+FX
         FYGTASE(2,NS)=FYGTASE(2,NS)+FY
         FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

         FXTPASE(2,NS)=FXTPASE(2,NS)-FX
         FYTPASE(2,NS)=FYTPASE(2,NS)-FY
         FZTPASE(2,NS)=FZTPASE(2,NS)-FZ

      END IF



!---- TETHER THE SECOND TPASE TO THE FIRST GTASE:

      DX=XGTASE(1,NS)-XTPASE(2,NS)
      DY=YGTASE(1,NS)-YTPASE(2,NS)
      DZ=ZGTASE(1,NS)-ZTPASE(2,NS)

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      IF(DIST>LGTTP2)THEN

         F=-KGTTP*(DIST-LGTTP2)/DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXGTASE(1,NS)=FXGTASE(1,NS)+FX
         FYGTASE(1,NS)=FYGTASE(1,NS)+FY
         FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

         FXTPASE(2,NS)=FXTPASE(2,NS)-FX
         FYTPASE(2,NS)=FYTPASE(2,NS)-FY
         FZTPASE(2,NS)=FZTPASE(2,NS)-FZ

      END IF


!---- TETHER THE THIRD TPASE TO THE SECOND GTASE:

      DX=XGTASE(2,NS)-XTPASE(3,NS)
      DY=YGTASE(2,NS)-YTPASE(3,NS)
      DZ=ZGTASE(2,NS)-ZTPASE(3,NS)

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      if(tppep(1,3,ns)==0)then
         l0=0.5*LGTTP13
      else
         l0=LGTTP13
      end if

      IF(DIST>l0)THEN

         F=-KGTTP*(DIST-l0)/DIST


         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXGTASE(2,NS)=FXGTASE(2,NS)+FX
         FYGTASE(2,NS)=FYGTASE(2,NS)+FY
         FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

         FXTPASE(3,NS)=FXTPASE(3,NS)-FX
         FYTPASE(3,NS)=FYTPASE(3,NS)-FY
         FZTPASE(3,NS)=FZTPASE(3,NS)-FZ

      END IF

!---- PAIR OF TRANSGLYCOSYLASES:

      DX=XGTASE(1,NS)-XGTASE(2,NS)
      DY=YGTASE(1,NS)-YGTASE(2,NS)
      DZ=ZGTASE(1,NS)-ZGTASE(2,NS)

      FX=-KPAIR*(DX-LPAIR*SYNDIR(NS))

      FXGTASE(1,NS)=FXGTASE(1,NS)+FX

      FXGTASE(2,NS)=FXGTASE(2,NS)-FX

!      DIST=SQRT(DY*DY+DZ*DZ)

!      INVDIST=1.0D0/(DIST+BETA)

      F=-2*KPAIR!*DIST

      FY=F*DY!*INVDIST
      FZ=F*DZ!*INVDIST

      FYGTASE(1,NS)=FYGTASE(1,NS)+FY
      FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

      FYGTASE(2,NS)=FYGTASE(2,NS)-FY
      FZGTASE(2,NS)=FZGTASE(2,NS)-FZ


!     CONSTRAIN THE FIRST TPASE TO LEFT SIDE OF THE FIRST GTASE:

!      DX=XTPASE(1,NS)-XGTASE(1,NS)

!      IF(DX*SYNDIR(NS)<LSIDE)THEN

!         F=SYNDIR(NS)*KSIDE*(DX-SYNDIR(NS)*LSIDE)**2

!         FXTPASE(1,NS)=FXTPASE(1,NS)+F

!         FXGTASE(1,NS)=FXGTASE(1,NS)-F

!      END IF

!---- CONSTRAIN THE SECOND TPASE IN BETWEEN 2 GTASES:

      DX=XTPASE(2,NS)-0.5D0*(XGTASE(1,NS)+XGTASE(2,NS))

      F=-KSIDE*DX

      FXTPASE(2,NS)=FXTPASE(2,NS)+F

      FXGTASE(1,NS)=FXGTASE(1,NS)-0.5D0*F

      FXGTASE(2,NS)=FXGTASE(2,NS)-0.5D0*F

!     CONSTRAIN THE THIRD TPASE TO RIGHT SIDE OF THE SECOND GTASE:

!      DX=XGTASE(2,NS)-XTPASE(3,NS)

!      IF(DX*SYNDIR(NS)<LSIDE)THEN

!         F=SYNDIR(NS)*KSIDE*(DX-SYNDIR(NS)*LSIDE)**2

!         FXGTASE(2,NS)=FXGTASE(2,NS)+F

!         FXTPASE(3,NS)=FXTPASE(3,NS)-F

!      END IF

!---- TETHER ENDOPEPTIDASE TO THE FIRST GTASE:

!      DX=XGTASE(1,NS)-XEDASE(NS)
!      DY=YGTASE(1,NS)-YEDASE(NS)
!      DZ=ZGTASE(1,NS)-ZEDASE(NS)

!      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

!      IF(DIST>LGTED)THEN

!         F=-KGTED*(DIST-LGTED)

!         INVDIST=1.0D0/DIST

!         FX=F*DX*INVDIST
!         FY=F*DY*INVDIST
!         FZ=F*DZ*INVDIST

!         FXGTASE(1,NS)=FXGTASE(1,NS)+FX
!         FYGTASE(1,NS)=FYGTASE(1,NS)+FY
!         FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

!         FXEDASE(NS)=FXEDASE(NS)-FX
!         FYEDASE(NS)=FYEDASE(NS)-FY
!         FZEDASE(NS)=FZEDASE(NS)-FZ

!      END IF

!---- TETHER ENDOPEPTIDASE TO THE SECOND GTASE:

!      DX=XGTASE(2,NS)-XEDASE(NS)
!      DY=YGTASE(2,NS)-YEDASE(NS)
!      DZ=ZGTASE(2,NS)-ZEDASE(NS)

!      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

!      IF(DIST>LGTED)THEN

!         F=-KGTED*(DIST-LGTED)

!         INVDIST=1.0D0/DIST

!         FX=F*DX*INVDIST
!         FY=F*DY*INVDIST
!         FZ=F*DZ*INVDIST

!         FXGTASE(2,NS)=FXGTASE(2,NS)+FX
!         FYGTASE(2,NS)=FYGTASE(2,NS)+FY
!         FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

!         FXEDASE(NS)=FXEDASE(NS)-FX
!         FYEDASE(NS)=FYEDASE(NS)-FY
!         FZEDASE(NS)=FZEDASE(NS)-FZ

!      END IF

!---- CONSTRAIN EDASE IN BETWEEN 2 GTASES:

!      DX=XEDASE(NS)-0.5D0*(XGTASE(1,NS)+XGTASE(2,NS))

!      F=-KSIDE*DX

!      FXEDASE(NS)=FXEDASE(NS)+F

!      FXGTASE(1,NS)=FXGTASE(1,NS)-0.5D0*F

!      FXGTASE(2,NS)=FXGTASE(2,NS)-0.5D0*F


!     CONSTRAIN EDASE ALWAYS IN between Gtase along X, and in FRONT OF GTASES along YZ:


      if(syndir(ns)==1)then

         if(xedase(ns)>xgtase(1,ns))then

            fx=kgted*(xedase(ns)-xgtase(1,ns))

            FXEDASE(NS)=FXEDASE(NS)-FX

            FXGTASE(1,NS)=FXGTASE(1,NS)+fx

         end if

         if(xedase(ns)<xgtase(2,ns))then

            fx=kgted*(xedase(ns)-xgtase(2,ns))

            FXEDASE(NS)=FXEDASE(NS)-FX

            FXGTASE(2,NS)=FXGTASE(2,NS)+fx

         end if

      else

         if(xedase(ns)<xgtase(1,ns))then

            fx=kgted*(xedase(ns)-xgtase(1,ns))

            FXEDASE(NS)=FXEDASE(NS)-FX

            FXGTASE(1,NS)=FXGTASE(1,NS)+fx

         end if

         if(xedase(ns)>xgtase(2,ns))then

            fx=kgted*(xedase(ns)-xgtase(2,ns))

            FXEDASE(NS)=FXEDASE(NS)-FX

            FXGTASE(2,NS)=FXGTASE(2,NS)+fx

         end if

      end if


      invrad=1.0d0/sqrt(YEDASE(NS)*YEDASE(NS)+ZEDASE(NS)*ZEDASE(NS))*syndir(ns)

      ylead=-ZEDASE(NS)*invrad
      zlead=YEDASE(NS)*invrad

!      DX=XEDASE(NS)-0.5D0*(XGTASE(1,NS)+XGTASE(2,NS))!-XLEAD(NS)*LGTED

      ygtase0=0.5D0*(YGTASE(1,NS)+YGTASE(2,NS))
      zgtase0=0.5D0*(ZGTASE(1,NS)+ZGTASE(2,NS))

      DY=YEDASE(NS)-ygtase0!0.5D0*(YGTASE(1,NS)+YGTASE(2,NS))!+ZEDASE(NS)*invrad*LGTED
      DZ=ZEDASE(NS)-zgtase0!0.5D0*(ZGTASE(1,NS)+ZGTASE(2,NS))!-YEDASE(NS)*invrad*LGTED

      proj=dy*ylead+dz*zlead

      if(edpep(1,ns)==0.and.edpep(2,ns)==0)then
         lgted_rep=lgted
      else
         lgted_rep=0.0d0
      end if

      if(proj<0.0d0)then

         f=-0.5d0*kgted*proj
         fy=f*ylead
         fz=f*zlead

         fyedase(ns)=fyedase(ns)+fy+fy
         fzedase(ns)=fzedase(ns)+fz+fz

         fygtase(1:2,ns)=fygtase(1:2,ns)-fy
         fzgtase(1:2,ns)=fzgtase(1:2,ns)-fz

      elseif(proj>lgted_rep)then

         f=-0.5d0*kgted*(proj-lgted_rep)
         fy=f*ylead
         fz=f*zlead

         fyedase(ns)=fyedase(ns)+fy+fy
         fzedase(ns)=fzedase(ns)+fz+fz

         fygtase(1:2,ns)=fygtase(1:2,ns)-fy
         fzgtase(1:2,ns)=fzgtase(1:2,ns)-fz

      end if


!     CONSTRAIN TPASE 1 and 3 ALWAYS in FRONT OF GTASES along YZ:

      dy=ytpase(1,ns)-ygtase0
      dz=ztpase(1,ns)-zgtase0

      proj=dy*ylead+dz*zlead

      if(proj<0.0d0)then

         f=-0.5d0*kgttp*proj
         fy=f*ylead
         fz=f*zlead

         fytpase(1,ns)=fytpase(1,ns)+fy+fy
         fztpase(1,ns)=fztpase(1,ns)+fz+fz

         fygtase(1:2,ns)=fygtase(1:2,ns)-fy
         fzgtase(1:2,ns)=fzgtase(1:2,ns)-fz

      end if

      dy=ytpase(3,ns)-ygtase0
      dz=ztpase(3,ns)-zgtase0

      proj=dy*ylead+dz*zlead

      if(proj<0.0d0)then

         f=-0.5d0*kgttp*proj
         fy=f*ylead
         fz=f*zlead

         fytpase(3,ns)=fytpase(3,ns)+fy+fy
         fztpase(3,ns)=fztpase(3,ns)+fz+fz

         fygtase(1:2,ns)=fygtase(1:2,ns)-fy
         fzgtase(1:2,ns)=fzgtase(1:2,ns)-fz

      end if



!     pull TPase to EDase for faster crosslinking

      if(tppep(1,1,ns)>0.and.tppep(2,1,ns)==0.and.edpep(1,ns)+edpep(2,ns)>0)then

         dx=xtpase(1,ns)-xedase(ns)
         dy=ytpase(1,ns)-yedase(ns)
         dz=ztpase(1,ns)-zedase(ns)

         DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

         f=ktped*(dist-ltped)/dist

         fx=f*dx
         fy=f*dy
         fz=f*dz

         fxtpase(1,ns)=fxtpase(1,ns)-fx
         fytpase(1,ns)=fytpase(1,ns)-fy
         fztpase(1,ns)=fztpase(1,ns)-fz

         fxedase(ns)=fxedase(ns)+fx
         fyedase(ns)=fyedase(ns)+fy
         fzedase(ns)=fzedase(ns)+fz

      end if

      if(tppep(1,3,ns)>0.and.tppep(2,3,ns)==0.and.edpep(1,ns)+edpep(2,ns)>0)then

         dx=xtpase(3,ns)-xedase(ns)
         dy=ytpase(3,ns)-yedase(ns)
         dz=ztpase(3,ns)-zedase(ns)

         DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

         f=ktped*(dist-ltped)/dist

         fx=f*dx
         fy=f*dy
         fz=f*dz

         fxtpase(3,ns)=fxtpase(3,ns)-fx
         fytpase(3,ns)=fytpase(3,ns)-fy
         fztpase(3,ns)=fztpase(3,ns)-fz

         fxedase(ns)=fxedase(ns)+fx
         fyedase(ns)=fyedase(ns)+fy
         fzedase(ns)=fzedase(ns)+fz

      end if



   END DO


!OMP END DO NOWAIT
!OMP END PARALLEL

   END SUBROUTINE

!=========================================================================

   SUBROUTINE GTAHOLD(NSYN,GLYTIP,GLYSEC,KGTASE,LGTASE,KLEAD,KTHETA,THETA_0,DELTA,INVDELTA,BETA,X,Y,Z, &
                      FXPG,FYPG,FZPG,XGTASE,YGTASE,ZGTASE,FXGTASE,FYGTASE,FZGTASE,syndir)

   IMPLICIT NONE

   INTEGER,VALUE:: NSYN
   INTEGER NS,N1,N2!,jsig

   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::GLYTIP,GLYSEC
   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::syndir

   DOUBLE PRECISION,VALUE:: KGTASE,KLEAD,KTHETA,THETA_0,DELTA,INVDELTA,BETA
   DOUBLE PRECISION DX,DY,DZ,F,FXGT,FYGT,FZGT,DIST,INVDIST,DX2,DY2,DZ2,DREP,REP
   DOUBLE PRECISION THET,COS_T,COS_T0,INVDIST2,F0,XL,YL,ZL,X1,Y1,Z1,invrad

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z!,XLEAD,YLEAD,ZLEAD
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FXPG,FYPG,FZPG
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::FXGTASE,FYGTASE,FZGTASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::XGTASE,YGTASE,ZGTASE,LGTASE




   DO NS=1,NSYN

!---  THE FIRST GTASE HOLDS GLYCAN TIP:

      IF(GLYTIP(1,NS)/=0)THEN

!jsig=jsig+1

!print*,'GTAhold'

         N1=GLYTIP(1,NS)

         x1=x(n1)
         y1=y(n1)
         z1=z(n1)

         DX=XGTASE(1,NS)-X1
         DY=YGTASE(1,NS)-Y1
         DZ=ZGTASE(1,NS)-Z1

         DIST=SQRT(DX**2+DY**2+DZ**2)

         INVDIST=1.0D0/DIST

         IF(LGTASE(1,NS)>1.0D0)THEN
            F=-KGTASE*(DIST-LGTASE(1,NS))
         ELSE
            F=-10*KGTASE*(DIST-LGTASE(1,NS))
         END IF

         F=F*INVDIST

!if(jsig>0) print*,'force 1',f

         FXgt=F*DX!*INVDIST
         FYgt=F*DY!*INVDIST
         FZgt=F*DZ!*INVDIST


!     CONSTRAIN STRAND TO LEADING DIRECTION AT THE TIP:

         xl=0.0d0!XLEAD(NS)

         invrad=1.0d0/sqrt(YGTASE(1,NS)*YGTASE(1,NS)+ZGTASE(1,NS)*ZGTASE(1,NS))*syndir(ns)

         yl=-ZGTASE(1,NS)*invrad!YLEAD(NS)
         zl=YGTASE(1,NS)*invrad!ZLEAD(NS)

         cos_t0=(DX*XL+DY*YL+DZ*ZL)*invdist

         THET=ACOS(COS_T0*(1.0D0-BETA))

         F0=Klead*THET/SIN(THET)*INVDELTA

!--------
         REP=DX+DELTA

         Drep=SQRT(REP**2+DY**2+DZ**2)

         cos_t=(REP*XL+DY*YL+DZ*ZL)/Drep

         f=F0*(COS_T-COS_T0)+fxgt

!if(jsig>0) print*,'force 11',f

         FXGTASE(1,NS)=FXGTASE(1,NS)+F

         FXPG(N1)=FXPG(N1)-F

!--------

         REP=DY+DELTA

         Drep=SQRT(DX**2+REP**2+DZ**2)

         cos_t=(DX*XL+REP*YL+DZ*ZL)/Drep


         f=F0*(COS_T-COS_T0)+fygt

!if(jsig>0) print*,'force 12',f

         FYGTASE(1,NS)=FYGTASE(1,NS)+F

         FYPG(N1)=FYPG(N1)-F

!--------

         REP=DZ+DELTA

         Drep=SQRT(DX**2+DY**2+REP**2)

         cos_t=(DX*XL+DY*YL+REP*ZL)/Drep


         f=F0*(COS_T-COS_T0)+fzgt

!if(jsig>0) print*,'force 13',f

         FZGTASE(1,NS)=FZGTASE(1,NS)+F

         FZPG(N1)=FZPG(N1)-F


      END IF




! --  TAKING INTO ACCOUNT THE ANGLE TERM AT THE NEW PG UNIT:

      IF(GLYSEC(1,NS)/=0)THEN

         N2=GLYSEC(1,NS)


         DX2=X(N2)-X1
         DY2=Y(N2)-Y1
         DZ2=Z(N2)-Z1

         invdist2=1.0d0/SQRT(DX2*dx2+DY2*dy2+DZ2*dz2)

         cos_t0=(DX*DX2+DY*DY2+DZ*DZ2)*invDIST*invdist2


         THET=ACOS(COS_T0*(1.0D0-BETA))

         F0=KTHETA*(THET-THETA_0)/SIN(THET)*INVDELTA

!-------

         REP=DX+DELTA

         DREP=SQRT(REP**2+DY**2+DZ**2)

         cos_t=(REP*DX2+DY*DY2+DZ*DZ2)/DREP*invdist2


         fxgt=F0*(COS_T-COS_T0)


         FXGTASE(1,NS)=FXGTASE(1,NS)+Fxgt


!--------

         REP=DY+DELTA

         DREP=SQRT(DX**2+REP**2+DZ**2)

         cos_t=(DX*DX2+REP*DY2+DZ*DZ2)/DREP*invdist2


         fygt=F0*(COS_T-COS_T0)


         FYGTASE(1,NS)=FYGTASE(1,NS)+Fygt


!-------

         REP=DZ+DELTA

         DREP=SQRT(DX**2+DY**2+REP**2)

         cos_t=(DX*DX2+DY*DY2+REP*DZ2)/DREP*invdist2

         fzgt=F0*(COS_T-COS_T0)

         FZGTASE(1,NS)=FZGTASE(1,NS)+Fzgt


!-------

         REP=DX2+DELTA

         DREP=SQRT(REP**2+DY2**2+DZ2**2)

         cos_t=(DX*REP+DY*DY2+DZ*DZ2)*invDIST/DREP

         f=F0*(COS_T-COS_T0)

         FXPG(N2)=FXPG(N2)+F

         FXPG(N1)=FXPG(N1)-F-fxgt

!--------

         REP=DY2+DELTA

         DREP=SQRT(DX2**2+REP**2+DZ2**2)

         cos_t=(DX*DX2+DY*REP+DZ*DZ2)*invDIST/DREP


         f=F0*(COS_T-COS_T0)


         FYPG(N2)=FYPG(N2)+F

         FYPG(N1)=FYPG(N1)-F-fygt

!--------

         REP=DZ2+DELTA

         DREP=SQRT(DX2**2+DY2**2+REP**2)

         cos_t=(DX*DX2+DY*DY2+DZ*REP)*invDIST/DREP


         f=F0*(COS_T-COS_T0)


         FZPG(N2)=FZPG(N2)+F

         FZPG(N1)=FZPG(N1)-F-fzgt

      END IF


!---- THE SECOND GTASE HOLDING THE SECOND STRAND TIP:

      IF(GLYTIP(2,NS)/=0)THEN

!jsig=jsig+1


         N1=GLYTIP(2,NS)

         x1=x(n1)
         y1=y(n1)
         z1=z(n1)

         DX=XGTASE(2,NS)-X1
         DY=YGTASE(2,NS)-Y1
         DZ=ZGTASE(2,NS)-Z1

         DIST=SQRT(DX**2+DY**2+DZ**2)

         INVDIST=1.0D0/DIST

         IF(LGTASE(2,NS)>1.0D0)THEN
            F=-KGTASE*(DIST-LGTASE(2,NS))
         ELSE
            F=-10*KGTASE*(DIST-LGTASE(2,NS))
         END IF

         F=F*INVDIST

!if(jsig>0) print*,'force 2',f


         FXgt=F*DX!*INVDIST
         FYgt=F*DY!*INVDIST
         FZgt=F*DZ!*INVDIST

!     CONSTRAIN STRAND TO LEADING DIRECTION AT THE TIP:

         xl=0.0d0!XLEAD(NS)

         invrad=1.0d0/sqrt(YGTASE(2,NS)*YGTASE(2,NS)+ZGTASE(2,NS)*ZGTASE(2,NS))*syndir(ns)

         yl=-ZGTASE(2,NS)*invrad!YLEAD(NS)
         zl=YGTASE(2,NS)*invrad!ZLEAD(NS)

         cos_t0=(DX*XL+DY*YL+DZ*ZL)*invdist

         THET=ACOS(COS_T0*(1.0D0-BETA))

         F0=Klead*THET/SIN(THET)*INVDELTA

!--------
         REP=DX+DELTA

         Drep=SQRT(REP**2+DY**2+DZ**2)

         cos_t=(REP*XL+DY*YL+DZ*ZL)/Drep


         f=F0*(COS_T-COS_T0)+fxgt

!if(jsig>0) print*,'force 21',f0,cos_t,cos_t0,fxgt



         FXGTASE(2,NS)=FXGTASE(2,NS)+F

         FXPG(N1)=FXPG(N1)-F

!--------

         REP=DY+DELTA

         Drep=SQRT(DX**2+REP**2+DZ**2)

         cos_t=(DX*XL+REP*YL+DZ*ZL)/Drep


         f=F0*(COS_T-COS_T0)+fygt
!if(jsig>0) print*,'force 22',f


         FYGTASE(2,NS)=FYGTASE(2,NS)+F

         FYPG(N1)=FYPG(N1)-F

!--------

         REP=DZ+DELTA

         Drep=SQRT(DX**2+DY**2+REP**2)

         cos_t=(DX*XL+DY*YL+REP*ZL)/Drep


         f=F0*(COS_T-COS_T0)+fzgt

!if(jsig>0) print*,'force 23',f

         FZGTASE(2,NS)=FZGTASE(2,NS)+F

         FZPG(N1)=FZPG(N1)-F

      END IF

! --  TAKING INTO ACCOUNT THE ANGLE TERM AT THE NEW PG UNIT:

      IF(GLYSEC(2,NS)/=0)THEN

         N2=GLYSEC(2,NS)

         DX2=X(N2)-X1
         DY2=Y(N2)-Y1
         DZ2=Z(N2)-Z1

         invdist2=1.0d0/SQRT(DX2*dx2+DY2*dy2+DZ2*dz2)

         cos_t0=(DX*DX2+DY*DY2+DZ*DZ2)*invDIST*invdist2

         THET=ACOS(COS_T0*(1.0D0-BETA))

         F0=KTHETA*(THET-THETA_0)/SIN(THET)*INVDELTA



!-------

         REP=DX+DELTA

         DREP=SQRT(REP**2+DY**2+DZ**2)

         cos_t=(REP*DX2+DY*DY2+DZ*DZ2)/DREP*invdist2


         fxgt=F0*(COS_T-COS_T0)


         FXGTASE(2,NS)=FXGTASE(2,NS)+Fxgt


!--------

         REP=DY+DELTA

         DREP=SQRT(DX**2+REP**2+DZ**2)

         cos_t=(DX*DX2+REP*DY2+DZ*DZ2)/DREP*invdist2


         fygt=F0*(COS_T-COS_T0)


         FYGTASE(2,NS)=FYGTASE(2,NS)+Fygt


!-------

         REP=DZ+DELTA

         DREP=SQRT(DX**2+DY**2+REP**2)

         cos_t=(DX*DX2+DY*DY2+REP*DZ2)/DREP*invdist2


         fzgt=F0*(COS_T-COS_T0)


         FZGTASE(2,NS)=FZGTASE(2,NS)+Fzgt


!-------

         REP=DX2+DELTA

         DREP=SQRT(REP**2+DY2**2+DZ2**2)

         cos_t=(DX*REP+DY*DY2+DZ*DZ2)*invDIST/DREP


         f=F0*(COS_T-COS_T0)



         FXPG(N2)=FXPG(N2)+F

         FXPG(N1)=FXPG(N1)-F-fxgt

!--------

         REP=DY2+DELTA

         DREP=SQRT(DX2**2+REP**2+DZ2**2)

         cos_t=(DX*DX2+DY*REP+DZ*DZ2)*invDIST/DREP


         f=F0*(COS_T-COS_T0)


         FYPG(N2)=FYPG(N2)+F

         FYPG(N1)=FYPG(N1)-F-fygt

!--------

         REP=DZ2+DELTA

         DREP=SQRT(DX2**2+DY2**2+REP**2)

         cos_t=(DX*DX2+DY*DY2+DZ*REP)*invDIST/DREP


         f=F0*(COS_T-COS_T0)


         FZPG(N2)=FZPG(N2)+F

         FZPG(N1)=FZPG(N1)-F-fzgt

      END IF


   END DO


   END SUBROUTINE

!===============================================================

   SUBROUTINE TPAHOLD(NSYN,X,Y,Z,FXPG,FYPG,FZPG,TPPEP,XTPASE,YTPASE,ZTPASE,FXTPASE,FYTPASE,FZTPASE,KTPASE,LTPASE)

   IMPLICIT NONE

   INTEGER,VALUE::NSYN
   INTEGER NS,N1,N2

   INTEGER,DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN)::TPPEP

   DOUBLE PRECISION,VALUE:: KTPASE,LTPASE
   DOUBLE PRECISION DX,DY,DZ,DIST,FX,FY,FZ,F

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FXPG,FYPG,FZPG
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::XTPASE,YTPASE,ZTPASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::FXTPASE,FYTPASE,FZTPASE




   DO NS=1,NSYN

!---- THE FIRST TPASE HOLDING A DONOR PEPTIDE:

      IF(TPPEP(1,1,NS)>0)THEN

         N1=TPPEP(1,1,NS)

         DX=XTPASE(1,NS)-X(N1)
         DY=YTPASE(1,NS)-Y(N1)
         DZ=ZTPASE(1,NS)-Z(N1)

         DIST=SQRT(DX**2+DY**2+DZ**2)

         IF(TPPEP(2,1,NS)==0)THEN
            F=-KTPASE*(DIST-LTPASE)
         ELSE
            F=-KTPASE*DIST
         END IF

         IF((DIST>LTPASE.AND.TPPEP(2,1,NS)==0).OR.TPPEP(2,1,NS)>0)THEN

            F=F/DIST

            FX=F*DX!*INVDIST
            FY=F*DY!*INVDIST
            FZ=F*DZ!*INVDIST

            FXTPASE(1,NS)=FXTPASE(1,NS)+FX
            FYTPASE(1,NS)=FYTPASE(1,NS)+FY
            FZTPASE(1,NS)=FZTPASE(1,NS)+FZ

            FXPG(N1)=FXPG(N1)-FX
            FYPG(N1)=FYPG(N1)-FY
            FZPG(N1)=FZPG(N1)-FZ
         END IF

      END IF

!---- THE FIRST TPASE HOLDING AN ACCEPTOR PEPTIDE:

      IF(TPPEP(2,1,NS)>0)THEN

         N1=TPPEP(2,1,NS)

         DX=XTPASE(1,NS)-X(N1)
         DY=YTPASE(1,NS)-Y(N1)
         DZ=ZTPASE(1,NS)-Z(N1)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXTPASE(1,NS)=FXTPASE(1,NS)+FX
         FYTPASE(1,NS)=FYTPASE(1,NS)+FY
         FZTPASE(1,NS)=FZTPASE(1,NS)+FZ

         FXPG(N1)=FXPG(N1)-FX
         FYPG(N1)=FYPG(N1)-FY
         FZPG(N1)=FZPG(N1)-FZ

      END IF

!     STRAIGHTEN UP THE POTENTIAL NEW PEPTIDE BOND:

      IF(TPPEP(1,1,NS)>0.AND.TPPEP(2,1,NS)>0)THEN

         N1=TPPEP(1,1,NS)

         N2=TPPEP(2,1,NS)

         DX=0.5D0*(X(N1)+X(N2))-XTPASE(1,NS)
         DY=0.5D0*(Y(N1)+Y(N2))-YTPASE(1,NS)
         DZ=0.5D0*(Z(N1)+Z(N2))-ZTPASE(1,NS)

 !        DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+0.5D0*FX
         FYPG(N1)=FYPG(N1)+0.5D0*FY
         FZPG(N1)=FZPG(N1)+0.5D0*FZ

         FXPG(N2)=FXPG(N2)+0.5D0*FX
         FYPG(N2)=FYPG(N2)+0.5D0*FY
         FZPG(N2)=FZPG(N2)+0.5D0*FZ

         FXTPASE(1,NS)=FXTPASE(1,NS)-FX
         FYTPASE(1,NS)=FYTPASE(1,NS)-FY
         FZTPASE(1,NS)=FZTPASE(1,NS)-FZ

      END IF

!---- THE SECOND TPASE HOLDING A DONOR PEPTIDE:

      IF(TPPEP(1,2,NS)>0)THEN

         N1=TPPEP(1,2,NS)

         DX=XTPASE(2,NS)-X(N1)
         DY=YTPASE(2,NS)-Y(N1)
         DZ=ZTPASE(2,NS)-Z(N1)

         DIST=SQRT(DX**2+DY**2+DZ**2)

         IF(TPPEP(2,1,NS)==0)THEN
            F=-KTPASE*(DIST-LTPASE)
         ELSE
            F=-KTPASE*DIST
         END IF

         IF((DIST>LTPASE.AND.TPPEP(2,2,NS)==0).OR.TPPEP(2,2,NS)>0)THEN

            F=F/DIST

            FX=F*DX!*INVDIST
            FY=F*DY!*INVDIST
            FZ=F*DZ!*INVDIST

            FXTPASE(2,NS)=FXTPASE(2,NS)+FX
            FYTPASE(2,NS)=FYTPASE(2,NS)+FY
            FZTPASE(2,NS)=FZTPASE(2,NS)+FZ

            FXPG(N1)=FXPG(N1)-FX
            FYPG(N1)=FYPG(N1)-FY
            FZPG(N1)=FZPG(N1)-FZ
         END IF

      END IF

!---- THE SECOND TPASE HOLDING AN ACCEPTOR PEPTIDE:

      IF(TPPEP(2,2,NS)>0)THEN

         N1=TPPEP(2,2,NS)

         DX=XTPASE(2,NS)-X(N1)
         DY=YTPASE(2,NS)-Y(N1)
         DZ=ZTPASE(2,NS)-Z(N1)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXTPASE(2,NS)=FXTPASE(2,NS)+FX
         FYTPASE(2,NS)=FYTPASE(2,NS)+FY
         FZTPASE(2,NS)=FZTPASE(2,NS)+FZ

         FXPG(N1)=FXPG(N1)-FX
         FYPG(N1)=FYPG(N1)-FY
         FZPG(N1)=FZPG(N1)-FZ

      END IF

!     STRAIGHTEN UP THE POTENTIAL NEW PEPTIDE BOND:

      IF(TPPEP(1,2,NS)>0.AND.TPPEP(2,2,NS)>0)THEN

         N1=TPPEP(1,2,NS)

         N2=TPPEP(2,2,NS)

         DX=0.5D0*(X(N1)+X(N2))-XTPASE(2,NS)
         DY=0.5D0*(Y(N1)+Y(N2))-YTPASE(2,NS)
         DZ=0.5D0*(Z(N1)+Z(N2))-ZTPASE(2,NS)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+0.5D0*FX
         FYPG(N1)=FYPG(N1)+0.5D0*FY
         FZPG(N1)=FZPG(N1)+0.5D0*FZ

         FXPG(N2)=FXPG(N2)+0.5D0*FX
         FYPG(N2)=FYPG(N2)+0.5D0*FY
         FZPG(N2)=FZPG(N2)+0.5D0*FZ

         FXTPASE(2,NS)=FXTPASE(2,NS)-FX
         FYTPASE(2,NS)=FYTPASE(2,NS)-FY
         FZTPASE(2,NS)=FZTPASE(2,NS)-FZ

      END IF

!---- THE THIRD TPASE HOLDING A DONOR PEPTIDE:

      IF(TPPEP(1,3,NS)>0)THEN

         N1=TPPEP(1,3,NS)

         DX=XTPASE(3,NS)-X(N1)
         DY=YTPASE(3,NS)-Y(N1)
         DZ=ZTPASE(3,NS)-Z(N1)

         DIST=SQRT(DX**2+DY**2+DZ**2)

         IF(TPPEP(2,3,NS)==0)THEN
            F=-KTPASE*(DIST-LTPASE)
         ELSE
            F=-KTPASE*DIST
         END IF

         IF((DIST>LTPASE.AND.TPPEP(2,3,NS)==0).OR.TPPEP(2,3,NS)>0)THEN

            F=F/DIST

            FX=F*DX!*INVDIST
            FY=F*DY!*INVDIST
            FZ=F*DZ!*INVDIST

            FXTPASE(3,NS)=FXTPASE(3,NS)+FX
            FYTPASE(3,NS)=FYTPASE(3,NS)+FY
            FZTPASE(3,NS)=FZTPASE(3,NS)+FZ

            FXPG(N1)=FXPG(N1)-FX
            FYPG(N1)=FYPG(N1)-FY
            FZPG(N1)=FZPG(N1)-FZ
         END IF

      END IF

!---- THE THIRD TPASE HOLDING AN ACCEPTOR PEPTIDE:

      IF(TPPEP(2,3,NS)>0)THEN

         N1=TPPEP(2,3,NS)

         DX=XTPASE(3,NS)-X(N1)
         DY=YTPASE(3,NS)-Y(N1)
         DZ=ZTPASE(3,NS)-Z(N1)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXTPASE(3,NS)=FXTPASE(3,NS)+FX
         FYTPASE(3,NS)=FYTPASE(3,NS)+FY
         FZTPASE(3,NS)=FZTPASE(3,NS)+FZ

         FXPG(N1)=FXPG(N1)-FX
         FYPG(N1)=FYPG(N1)-FY
         FZPG(N1)=FZPG(N1)-FZ

      END IF

!     STRAIGHTEN UP THE POTENTIAL NEW PEPTIDE BOND:

      IF(TPPEP(1,3,NS)>0.AND.TPPEP(2,3,NS)>0)THEN

         N1=TPPEP(1,3,NS)

         N2=TPPEP(2,3,NS)

         DX=0.5D0*(X(N1)+X(N2))-XTPASE(3,NS)
         DY=0.5D0*(Y(N1)+Y(N2))-YTPASE(3,NS)
         DZ=0.5D0*(Z(N1)+Z(N2))-ZTPASE(3,NS)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

         F=-KTPASE!*DIST

!         INVDIST=1.0D0/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+0.5D0*FX
         FYPG(N1)=FYPG(N1)+0.5D0*FY
         FZPG(N1)=FZPG(N1)+0.5D0*FZ

         FXPG(N2)=FXPG(N2)+0.5D0*FX
         FYPG(N2)=FYPG(N2)+0.5D0*FY
         FZPG(N2)=FZPG(N2)+0.5D0*FZ

         FXTPASE(3,NS)=FXTPASE(3,NS)-FX
         FYTPASE(3,NS)=FYTPASE(3,NS)-FY
         FZTPASE(3,NS)=FZTPASE(3,NS)-FZ

      END IF

   END DO


   END SUBROUTINE

!==========================================================

   SUBROUTINE EDAHOLD(NSYN,EDHOLD,EDPEP,KEDASE,LEDASE,X,Y,Z,FXPG,FYPG,FZPG, &
                      XEDASE,YEDASE,ZEDASE,FXEDASE,FYEDASE,FZEDASE)

   IMPLICIT NONE

   INTEGER,VALUE:: NSYN
   INTEGER NS,N1,N2
!   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::PEPDIR
   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::EDHOLD,EDPEP

   DOUBLE PRECISION,VALUE:: KEDASE,LEDASE
   DOUBLE PRECISION FX,FY,FZ,F,DX,DY,DZ,dist
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::FXPG,FYPG,FZPG
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FXEDASE,FYEDASE,FZEDASE
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN)::XEDASE,YEDASE,ZEDASE


   DO NS=1,NSYN

!---- ENDOPEPTIDASE HOLDING A BOND TO CLEAVE:

      IF(EDHOLD(1,NS)>0)THEN

         N1=EDHOLD(1,NS)
         N2=EDHOLD(2,NS)

!        FORCE HOLDING EDASE CLOSE TO THE CROSSLINK CENTER:

         DX=0.5D0*(X(N1)+X(N2))-XEDASE(NS)
         DY=0.5D0*(Y(N1)+Y(N2))-YEDASE(NS)
         DZ=0.5D0*(Z(N1)+Z(N2))-ZEDASE(NS)

!         DIST=SQRT(DX**2+DY**2+DZ**2)

!         INVDIST=1.0D0/DIST


         F=-KEDASE!*DIST


         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+0.5D0*FX
         FYPG(N1)=FYPG(N1)+0.5D0*FY
         FZPG(N1)=FZPG(N1)+0.5D0*FZ

         FXPG(N2)=FXPG(N2)+0.5D0*FX
         FYPG(N2)=FYPG(N2)+0.5D0*FY
         FZPG(N2)=FZPG(N2)+0.5D0*FZ

         FXEDASE(NS)=FXEDASE(NS)-FX
         FYEDASE(NS)=FYEDASE(NS)-FY
         FZEDASE(NS)=FZEDASE(NS)-FZ

!        FORCE STRETCHING THE CROSSLINK ON N1:

!         DX=X(N1)-XEDASE(NS)+PEPDIR(N1)*LEDASE
!         DY=Y(N1)-YEDASE(NS)
!         DZ=Z(N1)-ZEDASE(NS)


!         F=-KEDASE!*DIST

!         FX=F*DX!*INVDIST
!         FY=F*DY!*INVDIST
!         FZ=F*DZ!*INVDIST

!         FXPG(N1)=FXPG(N1)+FX
!         FYPG(N1)=FYPG(N1)+FY
!         FZPG(N1)=FZPG(N1)+FZ

!         FXEDASE(NS)=FXEDASE(NS)-FX
!         FYEDASE(NS)=FYEDASE(NS)-FY
!         FZEDASE(NS)=FZEDASE(NS)-FZ

!        FORCE STRETCHING THE CROSSLINK ON N2:

!         DX=X(N2)-XEDASE(NS)+PEPDIR(N2)*LEDASE
!         DY=Y(N2)-YEDASE(NS)
!         DZ=Z(N2)-ZEDASE(NS)

!         DIST=SQRT(DX**2+DY**2+DZ**2)
!         INVDIST=1.0D0/DIST

!         F=-KEDASE!*DIST

!         FX=F*DX!*INVDIST
!         FY=F*DY!*INVDIST
!         FZ=F*DZ!*INVDIST

!         FXPG(N2)=FXPG(N2)+FX
!         FYPG(N2)=FYPG(N2)+FY
!         FZPG(N2)=FZPG(N2)+FZ

!         FXEDASE(NS)=FXEDASE(NS)-FX
!         FYEDASE(NS)=FYEDASE(NS)-FY
!         FZEDASE(NS)=FZEDASE(NS)-FZ


      END IF

!---- ENDOPEPTIDASE TETHERED TO HOOKS FROM CLEAVED PEPTIDE BOND:

      IF(EDPEP(1,NS)>0)THEN

         N1=EDPEP(1,NS)

         DX=X(N1)-XEDASE(NS)
         DY=Y(N1)-YEDASE(NS)
         DZ=Z(N1)-ZEDASE(NS)

         DIST=SQRT(DX**2+DY**2+DZ**2)
!         INVDIST=1.0D0/DIST

         F=-KEDASE*(dist-ledase)/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+FX
         FYPG(N1)=FYPG(N1)+FY
         FZPG(N1)=FZPG(N1)+FZ

         FXEDASE(NS)=FXEDASE(NS)-FX
         FYEDASE(NS)=FYEDASE(NS)-FY
         FZEDASE(NS)=FZEDASE(NS)-FZ

      END IF

      IF(EDPEP(2,NS)>0)THEN

         N1=EDPEP(2,NS)

         DX=X(N1)-XEDASE(NS)!+PEPDIR(N1)*LEDASE
         DY=Y(N1)-YEDASE(NS)
         DZ=Z(N1)-ZEDASE(NS)

         DIST=SQRT(DX**2+DY**2+DZ**2)
!         INVDIST=1.0D0/DIST

         F=-KEDASE*(dist-ledase)/DIST

         FX=F*DX!*INVDIST
         FY=F*DY!*INVDIST
         FZ=F*DZ!*INVDIST

         FXPG(N1)=FXPG(N1)+FX
         FYPG(N1)=FYPG(N1)+FY
         FZPG(N1)=FZPG(N1)+FZ

         FXEDASE(NS)=FXEDASE(NS)-FX
         FYEDASE(NS)=FYEDASE(NS)-FY
         FZEDASE(NS)=FZEDASE(NS)-FZ

      END IF

   END DO




   END SUBROUTINE

!=========================================================================

   SUBROUTINE STERIC(NSYN,SYNTHESIS,SYNLOOP,LOOP,LOOPLEN,GLYNUM,NEIGLY,GLYTIP, &
              KWALL,X,Y,Z,FXPG,FYPG,FZPG,FXGTASE,FYGTASE,FZGTASE,XGTASE,YGTASE,ZGTASE)


   IMPLICIT NONE

   INTEGER,VALUE:: NSYN
   INTEGER NS,J,NA,NB,JL,ILOOP

   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::SYNTHESIS,GLYTIP,LOOP
   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::GLYNUM,SYNLOOP,LOOPLEN
   INTEGER,DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN)::NEIGLY

   DOUBLE PRECISION,value:: KWALL
   DOUBLE PRECISION F,PROJ,FX,FY,FZ,XPBP,YPBP,ZPBP,DX,DY,DZ,DIST
   DOUBLE PRECISION XA,YA,ZA,XB,YB,ZB,M1,M2,A,B,C,T,V,DET,DA,DB

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FXPG,FYPG,FZPG
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::FXGTASE,FYGTASE,FZGTASE
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::XGTASE,YGTASE,ZGTASE





   DO NS=1,NSYN


!---- INTERACTION WITH LIPOPROTEINS, CONSTRAINT FROM NEIGHBOR BONDS:

      IF(SYNTHESIS(1,NS)==0.AND.SYNTHESIS(2,NS)==0)THEN
         CYCLE
      END IF

      ILOOP=SYNLOOP(NS)

      DO JL=1,LOOPLEN(ILOOP)

         NA=LOOP(JL,ILOOP)

         IF(JL==LOOPLEN(ILOOP))THEN
            NB=LOOP(1,ILOOP)
         ELSE
            NB=LOOP(JL+1,ILOOP)
         END IF

!         IF(NA==GLYTIP(1,NS).OR.NA==GLYTIP(2,NS).OR.NB==GLYTIP(1,NS).OR.NB==GLYTIP(2,NS))THEN
!            CYCLE
!         END IF

         XA=X(NA); YA=Y(NA); ZA=Z(NA)

         XB=X(NB); YB=Y(NB); ZB=Z(NB)

!        FOR THE FIRST TRANSGLYCOSYLASE:

         IF(SYNTHESIS(1,NS)==1)THEN

            XPBP=XGTASE(1,NS); YPBP=YGTASE(1,NS); ZPBP=ZGTASE(1,NS)

            M1=(XPBP-XA)*(XB-XA)-YA*(YB-YA)-ZA*(ZB-ZA)

            M2=-YA*YPBP-ZA*ZPBP

            A=YPBP*(YB-YA)+ZPBP*(ZB-ZA)

            B=(XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2

            C=YPBP**2+ZPBP**2

            DET=A**2-B*C

            T=(M2*A-M1*C)/DET

!            IF(T<0.0D0.OR.T>1.0D0)THEN

!               PROJ=YA*YPBP+ZA*ZPBP

!               DX=XA-XPBP
!               DY=YA-PROJ*YPBP/(YPBP**2+ZPBP**2)
!               DZ=ZA-PROJ*ZPBP/(YPBP**2+ZPBP**2)

!               DIST=SQRT(DX**2+DY**2+DZ**2)

!               IF(DIST<0.5D0)THEN

!                  F=KWALL*(0.5D0-DIST)**2/DIST**3

!                  FX=F*DX!*INVDIST
!                  FY=F*DY!*INVDIST
!                  FZ=F*DZ!*INVDIST

!                  FXPG(NA)=FXPG(NA)+FX
!                  FYPG(NA)=FYPG(NA)+FY
!                  FZPG(NA)=FZPG(NA)+FZ

!                  FXGTASE(1,NS)=FXGTASE(1,NS)-FX
!                  FYGTASE(1,NS)=FYGTASE(1,NS)-FY
!                  FZGTASE(1,NS)=FZGTASE(1,NS)-FZ

!               END IF

!               PROJ=YB*YPBP+ZB*ZPBP

!               DX=XB-XPBP
!               DY=YB-PROJ*YPBP/(YPBP**2+ZPBP**2)
!               DZ=ZB-PROJ*ZPBP/(YPBP**2+ZPBP**2)

!               DIST=SQRT(DX**2+DY**2+DZ**2)

!               IF(DIST<0.5D0)THEN


!                  F=KWALL*(0.5D0-DIST)**2/DIST**3

!                  FX=F*DX!*INVDIST
!                  FY=F*DY!*INVDIST
!                  FZ=F*DZ!*INVDIST

!                  FXPG(NB)=FXPG(NB)+FX
!                  FYPG(NB)=FYPG(NB)+FY
!                  FZPG(NB)=FZPG(NB)+FZ

!                  FXGTASE(1,NS)=FXGTASE(1,NS)-FX
!                  FYGTASE(1,NS)=FYGTASE(1,NS)-FY
!                  FZGTASE(1,NS)=FZGTASE(1,NS)-FZ

!               END IF

!            ELSE

            if(t>0.0d0.and.t<1.0d0)then

               V=(M2*B-M1*A)/DET

               DX=XPBP-XA-T*(XB-XA)

               DY=-YA+V*YPBP-T*(YB-YA)

               DZ=-ZA+V*ZPBP-T*(ZB-ZA)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<0.5D0)THEN

!                  INVDIST=1.0D0/DIST

                  F=KWALL*(0.5D0-DIST)**2/DIST**3

                  FX=F*DX!*INVDIST
                  FY=F*DY!*INVDIST
                  FZ=F*DZ!*INVDIST

                  FXGTASE(1,NS)=FXGTASE(1,NS)+FX
                  FYGTASE(1,NS)=FYGTASE(1,NS)+FY
                  FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

                  FXPG(NA)=FXPG(NA)-0.5D0*FX
                  FYPG(NA)=FYPG(NA)-0.5D0*FY
                  FZPG(NA)=FZPG(NA)-0.5D0*FZ

                  FXPG(NB)=FXPG(NB)-0.5D0*FX
                  FYPG(NB)=FYPG(NB)-0.5D0*FY
                  FZPG(NB)=FZPG(NB)-0.5D0*FZ

               END IF

            END IF

         END IF

!        FOR THE SECOND TRANSGLYCOSYLASE:

         IF(SYNTHESIS(2,NS)==1)THEN

            XPBP=XGTASE(2,NS); YPBP=YGTASE(2,NS); ZPBP=ZGTASE(2,NS)

            M1=(XPBP-XA)*(XB-XA)-YA*(YB-YA)-ZA*(ZB-ZA)

            M2=-YA*YPBP-ZA*ZPBP

            A=YPBP*(YB-YA)+ZPBP*(ZB-ZA)

            B=(XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2

            C=YPBP**2+ZPBP**2

            DET=A**2-B*C

            T=(M2*A-M1*C)/DET

!            IF(T<0.0D0.OR.T>1.0D0)THEN

!               PROJ=YA*YPBP+ZA*ZPBP

!               DX=XA-XPBP
!               DY=YA-PROJ*YPBP/(YPBP**2+ZPBP**2)
!               DZ=ZA-PROJ*ZPBP/(YPBP**2+ZPBP**2)

!               DIST=SQRT(DX**2+DY**2+DZ**2)

!               IF(DIST<0.5D0)THEN

!                  F=KWALL*(0.5D0-DIST)**2/DIST**3

!                  FX=F*DX!*INVDIST
!                  FY=F*DY!*INVDIST
!                  FZ=F*DZ!*INVDIST

!                  FXPG(NA)=FXPG(NA)+FX
!                  FYPG(NA)=FYPG(NA)+FY
!                  FZPG(NA)=FZPG(NA)+FZ

!                  FXGTASE(2,NS)=FXGTASE(2,NS)-FX
!                  FYGTASE(2,NS)=FYGTASE(2,NS)-FY
!                  FZGTASE(2,NS)=FZGTASE(2,NS)-FZ

!               END IF

!               PROJ=YB*YPBP+ZB*ZPBP

!               DX=XB-XPBP
!               DY=YB-PROJ*YPBP/(YPBP**2+ZPBP**2)
!               DZ=ZB-PROJ*ZPBP/(YPBP**2+ZPBP**2)

!               DIST=SQRT(DX**2+DY**2+DZ**2)

!               IF(DIST<0.5D0)THEN


!                  F=KWALL*(0.5D0-DIST)**2/DIST**3

!                  FX=F*DX!*INVDIST
!                  FY=F*DY!*INVDIST
!                  FZ=F*DZ!*INVDIST

!                  FXPG(NB)=FXPG(NB)+FX
!                  FYPG(NB)=FYPG(NB)+FY
!                  FZPG(NB)=FZPG(NB)+FZ

!                  FXGTASE(2,NS)=FXGTASE(2,NS)-FX
!                  FYGTASE(2,NS)=FYGTASE(2,NS)-FY
!                  FZGTASE(2,NS)=FZGTASE(2,NS)-FZ

!               END IF

!            ELSE

            if(t>0.0d0.and.t<1.0d0)then

               V=(M2*B-M1*A)/DET

               DX=XPBP-XA-T*(XB-XA)

               DY=-YA+V*YPBP-T*(YB-YA)

               DZ=-ZA+V*ZPBP-T*(ZB-ZA)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<0.5D0)THEN

!                  INVDIST=1.0D0/DIST

                  F=KWALL*(0.5D0-DIST)**2/DIST**3

                  FX=F*DX!*INVDIST
                  FY=F*DY!*INVDIST
                  FZ=F*DZ!*INVDIST

                  FXGTASE(2,NS)=FXGTASE(2,NS)+FX
                  FYGTASE(2,NS)=FYGTASE(2,NS)+FY
                  FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

                  FXPG(NA)=FXPG(NA)-0.5D0*FX
                  FYPG(NA)=FYPG(NA)-0.5D0*FY
                  FZPG(NA)=FZPG(NA)-0.5D0*FZ

                  FXPG(NB)=FXPG(NB)-0.5D0*FX
                  FYPG(NB)=FYPG(NB)-0.5D0*FY
                  FZPG(NB)=FZPG(NB)-0.5D0*FZ

               END IF

            END IF

         END IF

      END DO

!----------------------------------------------------------------------------------

!     STERIC HINDRANCE PREVENTS OLD STRANDS GOING THROUGH BETWEEN TWO GTASES:

      IF(SYNTHESIS(1,NS)==0.OR.SYNTHESIS(2,NS)==0)THEN
         CYCLE
      END IF

!      ILOOP=SYNLOOP(NS)

!      DO JL=1,LOOPLEN(ILOOP)

      XA=XGTASE(1,NS); YA=YGTASE(1,NS); ZA=ZGTASE(1,NS)

      XB=XGTASE(2,NS); YB=YGTASE(2,NS); ZB=ZGTASE(2,NS)

      DO J=1,GLYNUM(NS)

!        THE FIRST BEAD ON THE BOND:

         NA=NEIGLY(1,J,NS)

         IF(GLYTIP(1,NS)==NA.OR.GLYTIP(2,NS)==NA)THEN
            GOTO 45
         END IF


!         XA=XGTASE(1,NS); YA=YGTASE(1,NS); ZA=ZGTASE(1,NS)

!         XB=XGTASE(2,NS); YB=YGTASE(2,NS); ZB=ZGTASE(2,NS)

         XPBP=X(NA); YPBP=Y(NA); ZPBP=Z(NA)

         DA=(XA-XPBP)**2+(YA-YPBP)**2+(ZA-ZPBP)**2
         DB=(XB-XPBP)**2+(YB-YPBP)**2+(ZB-ZPBP)**2

         IF(DA>4.0D0.AND.DB>4.0D0)THEN
            GOTO 45
         END IF

         M1=(XPBP-XA)*(XB-XA)-YA*(YB-YA)-ZA*(ZB-ZA)

         M2=-YA*YPBP-ZA*ZPBP

         A=YPBP*(YB-YA)+ZPBP*(ZB-ZA)

         B=(XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2

         C=YPBP**2+ZPBP**2

         DET=A**2-B*C

         T=(M2*A-M1*C)/DET

         IF(T<0.0D0.OR.T>1.0D0)THEN

            PROJ=YA*YPBP+ZA*ZPBP

            DX=XA-XPBP
            DY=YA-PROJ*YPBP/(YPBP**2+ZPBP**2)
            DZ=ZA-PROJ*ZPBP/(YPBP**2+ZPBP**2)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN
!               INVDIST=1.0D0/DIST
               F=KWALL*(0.5D0-DIST)**2/DIST**3

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)-FX
               FYPG(NA)=FYPG(NA)-FY
               FZPG(NA)=FZPG(NA)-FZ

               FXGTASE(1,NS)=FXGTASE(1,NS)+FX
               FYGTASE(1,NS)=FYGTASE(1,NS)+FY
               FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

            END IF

            PROJ=YB*YPBP+ZB*ZPBP

            DX=XB-XPBP
            DY=YB-PROJ*YPBP/(YPBP**2+ZPBP**2)
            DZ=ZB-PROJ*ZPBP/(YPBP**2+ZPBP**2)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN
!               INVDIST=1.0D0/DIST

               F=KWALL*(0.5D0-DIST)**2/DIST**3

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)-FX
               FYPG(NA)=FYPG(NA)-FY
               FZPG(NA)=FZPG(NA)-FZ

               FXGTASE(2,NS)=FXGTASE(2,NS)+FX
               FYGTASE(2,NS)=FYGTASE(2,NS)+FY
               FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

            END IF

         ELSE

            V=(M2*B-M1*A)/DET

            DX=XPBP-XA-T*(XB-XA)

            DY=-YA+V*YPBP-T*(YB-YA)

            DZ=-ZA+V*ZPBP-T*(ZB-ZA)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN

!               INVDIST=1.0D0/DIST

               F=KWALL*(0.5D0-DIST)**2/DIST**3

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)+FX
               FYPG(NA)=FYPG(NA)+FY
               FZPG(NA)=FZPG(NA)+FZ

               FXGTASE(1,NS)=FXGTASE(1,NS)-0.5D0*FX
               FYGTASE(1,NS)=FYGTASE(1,NS)-0.5D0*FY
               FZGTASE(1,NS)=FZGTASE(1,NS)-0.5D0*FZ

               FXGTASE(2,NS)=FXGTASE(2,NS)-0.5D0*FX
               FYGTASE(2,NS)=FYGTASE(2,NS)-0.5D0*FY
               FZGTASE(2,NS)=FZGTASE(2,NS)-0.5D0*FZ

            END IF

         END IF


!-------------
!        THE SECOND BEAD ON THE BOND:

45       NA=NEIGLY(2,J,NS)

         IF(GLYTIP(1,NS)==NA.OR.GLYTIP(2,NS)==NA)THEN
            CYCLE
         END IF


!         XA=XGTASE(1,NS); YA=YGTASE(1,NS); ZA=ZGTASE(1,NS)

!         XB=XGTASE(2,NS); YB=YGTASE(2,NS); ZB=ZGTASE(2,NS)

         XPBP=X(NA); YPBP=Y(NA); ZPBP=Z(NA)

         DA=(XA-XPBP)**2+(YA-YPBP)**2+(ZA-ZPBP)**2
         DB=(XB-XPBP)**2+(YB-YPBP)**2+(ZB-ZPBP)**2

         IF(DA>4.0D0.AND.DB>4.0D0)THEN
            CYCLE
         END IF

         M1=(XPBP-XA)*(XB-XA)-YA*(YB-YA)-ZA*(ZB-ZA)

         M2=-YA*YPBP-ZA*ZPBP

         A=YPBP*(YB-YA)+ZPBP*(ZB-ZA)

         B=(XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2

         C=YPBP**2+ZPBP**2

         DET=A**2-B*C

         T=(M2*A-M1*C)/DET

         IF(T<0.0D0.OR.T>1.0D0)THEN

            PROJ=YA*YPBP+ZA*ZPBP

            DX=XA-XPBP
            DY=YA-PROJ*YPBP/(YPBP**2+ZPBP**2)
            DZ=ZA-PROJ*ZPBP/(YPBP**2+ZPBP**2)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN
!               INVDIST=1.0D0/DIST
               F=KWALL*(0.5D0-DIST)**2/DIST**3

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)-FX
               FYPG(NA)=FYPG(NA)-FY
               FZPG(NA)=FZPG(NA)-FZ

               FXGTASE(1,NS)=FXGTASE(1,NS)+FX
               FYGTASE(1,NS)=FYGTASE(1,NS)+FY
               FZGTASE(1,NS)=FZGTASE(1,NS)+FZ

            END IF

            PROJ=YB*YPBP+ZB*ZPBP

            DX=XB-XPBP
            DY=YB-PROJ*YPBP/(YPBP**2+ZPBP**2)
            DZ=ZB-PROJ*ZPBP/(YPBP**2+ZPBP**2)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN
!               INVDIST=1.0D0/DIST

               F=KWALL*(0.5D0-DIST)**2/DIST**2

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)-FX
               FYPG(NA)=FYPG(NA)-FY
               FZPG(NA)=FZPG(NA)-FZ

               FXGTASE(2,NS)=FXGTASE(2,NS)+FX
               FYGTASE(2,NS)=FYGTASE(2,NS)+FY
               FZGTASE(2,NS)=FZGTASE(2,NS)+FZ

            END IF

         ELSE

            V=(M2*B-M1*A)/DET

            DX=XPBP-XA-T*(XB-XA)

            DY=-YA+V*YPBP-T*(YB-YA)

            DZ=-ZA+V*ZPBP-T*(ZB-ZA)

            DIST=SQRT(DX**2+DY**2+DZ**2)

            IF(DIST<0.5D0)THEN

!               INVDIST=1.0D0/DIST

               F=KWALL*(0.5D0-DIST)**2/DIST**3

               FX=F*DX!*INVDIST
               FY=F*DY!*INVDIST
               FZ=F*DZ!*INVDIST

               FXPG(NA)=FXPG(NA)+FX
               FYPG(NA)=FYPG(NA)+FY
               FZPG(NA)=FZPG(NA)+FZ

               FXGTASE(1,NS)=FXGTASE(1,NS)-0.5D0*FX
               FYGTASE(1,NS)=FYGTASE(1,NS)-0.5D0*FY
               FZGTASE(1,NS)=FZGTASE(1,NS)-0.5D0*FZ

               FXGTASE(2,NS)=FXGTASE(2,NS)-0.5D0*FX
               FYGTASE(2,NS)=FYGTASE(2,NS)-0.5D0*FY
               FZGTASE(2,NS)=FZGTASE(2,NS)-0.5D0*FZ

            END IF

         END IF

      END DO


   END DO


   END SUBROUTINE

!==========================================

   SUBROUTINE ESTRAND(NATOM,FXGLY,FYGLY,FZGLY,FXTHETA,FYTHETA,FZTHETA, &
               NCAP,NPG,PGID,PGLEN,PGTYP,X,Y,Z,LBOND,KBOND,THET0,KTHETA,JFORCE,DELTA,INVDELTA,BETA,PI)

   IMPLICIT NONE

   INTEGER,VALUE:: NCAP,NPG,JFORCE,NATOM
   INTEGER N,I,N1,N2,N3,L

   INTEGER,DIMENSION(:,:),ALLOCATABLE::PGID
   INTEGER,DIMENSION(:),ALLOCATABLE:: PGLEN,PGTYP

   DOUBLE PRECISION PI,LBOND,KBOND,THET0,KTHETA,DELTA,INVDELTA,BETA
   DOUBLE PRECISION DX,DY,DZ,DX1,DY1,DZ1,DX3,DY3,DZ3,DREP,DIST1,INVDIST,DIST3
   DOUBLE PRECISION COS_T,COS_T0,THET,F0,F,FX,FY,FZ
   DOUBLE PRECISION DFX1,DFY1,DFZ1,DFX3,DFY3,DFZ3
   DOUBLE PRECISION x1,y1,z1,x2,y2,z2,x3,y3,z3

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FXGLY,FYGLY,FZGLY,FXTHETA,FYTHETA,FZTHETA
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z


!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(N,I,N1,N2,x1,y1,z1,x2,y2,z2) &
!$OMP PRIVATE(DX1,DY1,DZ1,DIST1,INVDIST,f) &
!$OMP SHARED(NPG,LBOND,KBOND) &
!$OMP SHARED(FXGLY,FYGLY,FZGLY,X,Y,Z,PGID,PGLEN) 


!$OMP DO SCHEDULE(GUIDED,64)

   DO N=1,NPG


      N2=PGID(1,N)

      FXGLY(N2)=0.0D0
      FYGLY(N2)=0.0D0
      FZGLY(N2)=0.0D0

      IF(PGLEN(N)<2)THEN
         CYCLE
      END IF

      X2=X(N2)
      Y2=Y(N2)
      Z2=Z(N2)

      DO I=2,PGLEN(N)

         n1=n2

         x1=x2
         y1=y2
         z1=z2

         N2=PGID(I,N)

         x2=x(n2)
         y2=y(n2)
         z2=z(n2)

         dx1=x1-x2
         dy1=y1-y2
         dz1=z1-z2

         DIST1=SQRT(DX1**2+DY1**2+DZ1**2)

         INVDIST=1.0D0/DIST1


         F=KBOND*(DIST1-LBOND)

         FXGLY(N2)=F*DX1*INVDIST
         FYGLY(N2)=F*DY1*INVDIST
         FZGLY(N2)=F*DZ1*INVDIST

         FXGLY(N1)=FXGLY(N1)-FXGLY(N2)
         FYGLY(N1)=FYGLY(N1)-FYGLY(N2)
         FZGLY(N1)=FZGLY(N1)-FZGLY(N2)

      end do

   end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

   DO N=1,NCAP

      N1=PGID(1,N)

      N2=PGID(PGLEN(N),N)

      DX1=X(N1)-X(N2)
      DY1=Y(N1)-Y(N2)
      DZ1=Z(N1)-Z(N2)

      DIST1=SQRT(DX1**2+DY1**2+DZ1**2)

      INVDIST=1.0D0/DIST1


      F=KBOND*(DIST1-LBOND)

      FX=F*DX1*INVDIST
      FY=F*DY1*INVDIST
      FZ=F*DZ1*INVDIST

      FXGLY(N1)=FXGLY(N1)-FX
      FYGLY(N1)=FYGLY(N1)-FY
      FZGLY(N1)=FZGLY(N1)-FZ

      FXGLY(N2)=FXGLY(N2)+FX
      FYGLY(N2)=FYGLY(N2)+FY
      FZGLY(N2)=FZGLY(N2)+FZ

   END DO

!-------------------------------------

   IF(JFORCE==0)THEN
      RETURN
   END IF


!-------------------------------------

!  FORCES FROM BENDING

!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(N,I,N1,N2,N3,x1,y1,z1,x2,y2,z2,x3,y3,z3) &
!$OMP PRIVATE(DX,DY,DZ,DX1,DY1,DZ1,DX3,DY3,DZ3,DREP,DIST1,DIST3) &
!$OMP PRIVATE(COS_T,COS_T0,THET,F0,DFX1,DFY1,DFZ1) &
!$OMP SHARED(NPG,THET0,KTHETA,DELTA,INVDELTA,BETA) &
!$OMP SHARED(FXTHETA,FYTHETA,FZTHETA,X,Y,Z,PGID,PGLEN) 


!$OMP DO SCHEDULE(GUIDED,64)

   DO N=1,NPG

      n2=pgid(1,n)

      FXTHETA(N2)=0.0D0
      FYTHETA(N2)=0.0D0
      FZTHETA(N2)=0.0D0


      IF(PGLEN(N)<2)THEN
         CYCLE
      END IF

      n3=pgid(2,n)

      FXTHETA(N3)=0.0D0
      FYTHETA(N3)=0.0D0
      FZTHETA(N3)=0.0D0

      IF(PGLEN(N)<3)THEN
         CYCLE
      END IF

      x3=x(n3)
      y3=y(n3)
      z3=z(n3)

      dx3=x3-x(n2)
      dy3=y3-y(n2)
      dz3=z3-z(n2)

      dist3=sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      do i=2,pglen(n)-1

         dx1=-dx3
         dy1=-dy3
         dz1=-dz3

         dist1=dist3

         x2=x3
         y2=y3
         z2=z3

         n1=n2

         n2=n3

         N3=PGID(I+1,N)

         x3=x(n3)
         y3=y(n3)
         z3=z(n3)

         dx3=x3-x2
         dy3=y3-y2
         dz3=z3-z2

         DIST3=SQRT(DX3**2+DY3**2+DZ3**2)

         COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)/DIST1/DIST3

         THET=ACOS(COS_T0*(1.0D0-BETA))


         F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA


!------ FORCE ON N1 ALONG X:

         DX=DX1+DELTA

         DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

         COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP/DIST3


         DFX1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Y:

         DY=DY1+DELTA

         DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

         COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP/DIST3


         DFY1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Z:

         DZ=DZ1+DELTA

         DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

         COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP/DIST3


         DFZ1=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG X:

         DX=DX3+DELTA

         DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

         COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)/DIST1/DREP


         FXTHETA(N3)=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Y:

         DY=DY3+DELTA

         DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

         COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)/DIST1/DREP


         FYTHETA(N3)=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Z:

         DZ=DZ3+DELTA

         DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

         COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)/DIST1/DREP


         FZTHETA(N3)=F0*(COS_T-COS_T0)

!       FORCES ON N1:

         FXTHETA(N1)=FXTHETA(N1)+DFX1
         FYTHETA(N1)=FYTHETA(N1)+DFY1
         FZTHETA(N1)=FZTHETA(N1)+DFZ1

!       FORCES ON N2:

         FXTHETA(N2)=FXTHETA(N2)-DFX1-FXTHETA(N3)
         FYTHETA(N2)=FYTHETA(N2)-DFY1-FYTHETA(N3)
         FZTHETA(N2)=FZTHETA(N2)-DFZ1-FZTHETA(N3)


      END DO

   END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL


   DO N=1,NCAP

      L=PGLEN(N)

      DO I=1,2

         N1=PGID(L+1-I,N)

         N2=PGID(2-I+L*(I-1),N)

         N3=PGID(3-I,N)

         DX1=X(N1)-X(N2)
         DY1=Y(N1)-Y(N2)
         DZ1=Z(N1)-Z(N2)

         DIST1=SQRT(DX1**2+DY1**2+DZ1**2)


         DX3=X(N3)-X(N2)
         DY3=Y(N3)-Y(N2)
         DZ3=Z(N3)-Z(N2)


         DIST3=SQRT(DX3**2+DY3**2+DZ3**2)

         COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)/DIST1/DIST3

         THET=ACOS(COS_T0*(1.0D0-BETA))


         F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA

!------ FORCE ON N1 ALONG X:

         DX=DX1+DELTA

         DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

         COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP/DIST3


         DFX1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Y:

         DY=DY1+DELTA

         DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

         COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP/DIST3


         DFY1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Z:

         DZ=DZ1+DELTA

         DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

         COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP/DIST3


         DFZ1=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG X:

         DX=DX3+DELTA

         DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

         COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)/DIST1/DREP

         DFX3=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Y:

         DY=DY3+DELTA

         DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

         COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)/DIST1/DREP

         DFY3=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Z:

         DZ=DZ3+DELTA

         DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

         COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)/DIST1/DREP

         DFZ3=F0*(COS_T-COS_T0)

!       FORCES ON N1:

         FXTHETA(N1)=FXTHETA(N1)+DFX1
         FYTHETA(N1)=FYTHETA(N1)+DFY1
         FZTHETA(N1)=FZTHETA(N1)+DFZ1

!       FORCES ON N2:

         FXTHETA(N2)=FXTHETA(N2)-DFX1-DFX3
         FYTHETA(N2)=FYTHETA(N2)-DFY1-DFY3
         FZTHETA(N2)=FZTHETA(N2)-DFZ1-DFZ3

!       FORCES ON N3:

         FXTHETA(N3)=FXTHETA(N3)+DFX3
         FYTHETA(N3)=FYTHETA(N3)+DFY3
         FZTHETA(N3)=FZTHETA(N3)+DFZ3

      END DO

   END DO




   END SUBROUTINE
!========================================================================

   SUBROUTINE FGLYBONDS(NPG,PGID,PGLEN,L_G,K_G,FX,FY,FZ,X,Y,Z)

   IMPLICIT NONE

   INTEGER N,N1,N2,I
   INTEGER,VALUE::NPG

   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::PGID
   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::PGLEN

   DOUBLE PRECISION,VALUE:: L_G,K_G
   DOUBLE PRECISION DX,DY,DZ,DIST,DFX,DFY,DFZ,F,X1,Y1,Z1,X2,Y2,Z2

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FX,FY,FZ
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z



!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(N,I,N1,N2,X1,Y1,Z1,X2,Y2,Z2,DX,DY,DZ,DIST,DFX,DFY,DFZ,F) &
!$OMP SHARED(NPG,PGID,PGLEN,FX,FY,FZ,X,Y,Z,L_G,K_G)

!$OMP DO SCHEDULE(GUIDED,64)


DO N=1,NPG

   IF(PGLEN(N)<2)THEN
      CYCLE
   END IF

   N2=PGID(1,N)

   X2=X(N2)
   Y2=Y(N2)
   Z2=Z(N2)

   DO I=2,PGLEN(N)

      N1=N2

      X1=X2
      Y1=Y2
      Z1=Z2

      N2=PGID(I,N)

      X2=X(N2)
      Y2=Y(N2)
      Z2=Z(N2)

      DX=X1-X2
      DY=Y1-Y2
      DZ=Z1-Z2

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      F=K_G*(DIST-L_G)/DIST

      DFX=F*DX
      DFY=F*DY
      DFZ=F*DZ

      FX(N1)=FX(N1)-DFX

      FY(N1)=FY(N1)-DFY

      FZ(N1)=FZ(N1)-DFZ

      FX(N2)=FX(N2)+DFX

      FY(N2)=FY(N2)+DFY

      FZ(N2)=FZ(N2)+DFZ

   END DO

END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL



!   DO N=1,2

!      n2=caplen*n

!      x2=x(n2)
!      y2=y(n2)
!      z2=z(n2)

!      do i=1,caplen

!         n1=n2

!         x1=x2
!         y1=y2
!         z1=z2

!         n2=i+(n-1)*caplen

!         x2=x(n2)
!         y2=y(n2)
!         z2=z(n2)

!         if(capopen(n1)==1.and.capopen(n2)==1) cycle

!         dx=x1-x2
!         dy=y1-y2
!         dz=z1-z2



!         DIST=SQRT(DX**2+DY**2+DZ**2)

!         F=K_G*(DIST-L_G)/DIST

!         DFX=F*DX!*INVDIST
!         DFY=F*DY!*INVDIST
!         DFZ=F*DZ!*INVDIST

!         FX(N1)=FX(N1)-DFX
!         FY(N1)=FY(N1)-DFY
!         FZ(N1)=FZ(N1)-DFZ

!         FX(N2)=FX(N2)+DFX
!         FY(N2)=FY(N2)+DFY
!         FZ(N2)=FZ(N2)+DFZ

!      end do

!   END DO

   END SUBROUTINE



!==========================================


   SUBROUTINE FANGLES(NPG,PGID,PGLEN,KTHETA,THET0,BETA,DELTA,INVDELTA,X,Y,Z,FX,FY,FZ)

   INTEGER N,N1,N2,N3,I
   INTEGER,VALUE::NPG
   INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::PGID
   INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::PGLEN

   REAL(KIND=8),VALUE::KTHETA,THET0,BETA,INVDELTA,DELTA
   REAL(KIND=8)::DX,DY,DZ,DX1,DY1,DZ1,DX3,DY3,DZ3,INVDIST1,INVDIST3,COS_T0,COS_T,THET,F0
   REAL(KIND=8)::DFX1,DFY1,DFZ1,DFX3,DFY3,DFZ3,DREP,X2,Y2,Z2,X3,Y3,Z3
   REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(IN)::X,Y,Z
   REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::FX,FY,FZ

!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(N,I,N1,N2,N3,DX,DY,DZ,DX1,DY1,DZ1,DX3,DY3,DZ3,INVDIST1,INVDIST3) &
!$OMP PRIVATE(X2,Y2,Z2,X3,Y3,Z3) &
!$OMP PRIVATE(COS_T0,COS_T,THET,F0,DFX1,DFY1,DFZ1,DFX3,DFY3,DFZ3,DREP) &
!$OMP SHARED(NPG,PGID,PGLEN,KTHETA,THET0,BETA,INVDELTA,DELTA,X,Y,Z,FX,FY,FZ)

!$OMP DO SCHEDULE(GUIDED,64)

DO N=1,NPG

   IF(PGLEN(N)<3)THEN
      CYCLE
   END IF

   N2=PGID(1,N)


   N3=PGID(2,N)

   X3=X(N3)
   Y3=Y(N3)
   Z3=Z(N3)

   DX3=X3-X(N2)
   DY3=Y3-Y(N2)
   DZ3=Z3-Z(N2)

   INVDIST3=1.0D0/SQRT(DX3**2+DY3**2+DZ3**2)

   DO I=3,PGLEN(N)

      DX1=-DX3
      DY1=-DY3
      DZ1=-DZ3


      N1=N2
      N2=N3


      X2=X3
      Y2=Y3
      Z2=Z3

      INVDIST1=INVDIST3

      N3=PGID(I,N)

      X3=X(N3)
      Y3=Y(N3)
      Z3=Z(N3)

      DX3=X3-X2
      DY3=Y3-Y2
      DZ3=Z3-Z2

      INVDIST3=1.0D0/SQRT(DX3**2+DY3**2+DZ3**2)



      COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)*invDIST1*invDIST3

      THET=ACOS(COS_T0*(1.0D0-BETA))

      F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA


!---- FORCE ON N1 ALONG X:

      DX=DX1+DELTA

      DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

      COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP*invDIST3


      DFX1=F0*(COS_T-COS_T0)


!---- FORCE ON N1 ALONG Y:

      DY=DY1+DELTA

      DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

      COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP*invDIST3


      DFY1=F0*(COS_T-COS_T0)


!---- FORCE ON N1 ALONG Z:

      DZ=DZ1+DELTA

      DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP*invDIST3


      DFZ1=F0*(COS_T-COS_T0)


!---- FORCE ON N3 ALONG X:

      DX=DX3+DELTA

      DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

      COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)*invDIST1/DREP

      DFX3=F0*(COS_T-COS_T0)


!---- FORCE ON N3 ALONG Y:

      DY=DY3+DELTA

      DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

      COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)*invDIST1/DREP

      DFY3=F0*(COS_T-COS_T0)


!---- FORCE ON N3 ALONG Z:

      DZ=DZ3+DELTA

      DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

      COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)*invDIST1/DREP

      DFZ3=F0*(COS_T-COS_T0)

!     FORCES ON N1:

      FX(N1)=FX(N1)+DFX1

      FY(N1)=FY(N1)+DFY1

      FZ(N1)=FZ(N1)+DFZ1

!     FORCES ON N2:

      FX(N2)=FX(N2)-DFX1-DFX3

      FY(N2)=FY(N2)-DFY1-DFY3

      FZ(N2)=FZ(N2)-DFZ1-DFZ3

!     FORCES ON N3:

      FX(N3)=FX(N3)+DFX3

      FY(N3)=FY(N3)+DFY3

      FZ(N3)=FZ(N3)+DFZ3

   END DO

END DO



!$OMP END DO NOWAIT
!$OMP END PARALLEL


!------------------------------------------


!   DO N=1,2


!      n2=n*caplen-1

!      n3=n*caplen

!      X3=X(N3)
!      Y3=Y(N3)
!      Z3=Z(N3)

!      DX3=X3-X(N2)
!      DY3=Y3-Y(N2)
!      DZ3=Z3-Z(N2)

!      INVDIST3=1.0D0/SQRT(DX3**2+DY3**2+DZ3**2)

!      DO I=1,caplen

!         DX1=-DX3
!         DY1=-DY3
!         DZ1=-DZ3


!         N1=N2
!         N2=N3


!         X2=X3
!         Y2=Y3
!         Z2=Z3

!         INVDIST1=INVDIST3

!         n3=i+(n-1)*caplen

!         X3=X(N3)
!         Y3=Y(N3)
!         Z3=Z(N3)

!         DX3=X3-X2
!         DY3=Y3-Y2
!         DZ3=Z3-Z2


!         invDIST3=1.0D0/SQRT(DX3**2+DY3**2+DZ3**2)

!         if(capopen(n1)+capopen(n2)+capopen(n3)==2) cycle


!         COS_T0=(DX1*DX3+DY1*DY3+DZ1*DZ3)*invDIST1*invDIST3

!         THET=ACOS(COS_T0*(1.0D0-BETA))


!         F0=KTHETA*(THET-THET0)/SIN(THET)*INVDELTA

!------ FORCE ON N1 ALONG X:

!         DX=DX1+DELTA

!         DREP=SQRT(DX*DX+DY1*DY1+DZ1*DZ1)

!         COS_T=(DX*DX3+DY1*DY3+DZ1*DZ3)/DREP*invDIST3


!         DFX1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Y:

!         DY=DY1+DELTA

!         DREP=SQRT(DX1*DX1+DY*DY+DZ1*DZ1)

!         COS_T=(DX1*DX3+DY*DY3+DZ1*DZ3)/DREP*invDIST3


!         DFY1=F0*(COS_T-COS_T0)

!------ FORCE ON N1 ALONG Z:

!         DZ=DZ1+DELTA

!         DREP=SQRT(DX1*DX1+DY1*DY1+DZ*DZ)

!         COS_T=(DX1*DX3+DY1*DY3+DZ*DZ3)/DREP*invDIST3


!         DFZ1=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG X:

!         DX=DX3+DELTA

!         DREP=SQRT(DX*DX+DY3*DY3+DZ3*DZ3)

!         COS_T=(DX1*DX+DY1*DY3+DZ1*DZ3)*invDIST1/DREP

!         DFX3=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Y:

!         DY=DY3+DELTA

!         DREP=SQRT(DX3*DX3+DY*DY+DZ3*DZ3)

!         COS_T=(DX1*DX3+DY1*DY+DZ1*DZ3)*invDIST1/DREP

!         DFY3=F0*(COS_T-COS_T0)

!------ FORCE ON N3 ALONG Z:

!         DZ=DZ3+DELTA

!         DREP=SQRT(DX3*DX3+DY3*DY3+DZ*DZ)

!         COS_T=(DX1*DX3+DY1*DY3+DZ1*DZ)*invDIST1/DREP

!         DFZ3=F0*(COS_T-COS_T0)

!       FORCES ON N1:

!         FX(N1)=FX(N1)+DFX1
!         FY(N1)=FY(N1)+DFY1
!         FZ(N1)=FZ(N1)+DFZ1

!       FORCES ON N2:

!         FX(N2)=FX(N2)-DFX1-DFX3
!         FY(N2)=FY(N2)-DFY1-DFY3
!         FZ(N2)=FZ(N2)-DFZ1-DFZ3

!       FORCES ON N3:

!         FX(N3)=FX(N3)+DFX3
!         FY(N3)=FY(N3)+DFY3
!         FZ(N3)=FZ(N3)+DFZ3

!      END DO

!   END DO

   END SUBROUTINE



!========================================================================

   SUBROUTINE FPEPBONDS(FX,FY,FZ,X,Y,Z,NBOND,BOND,BONTYP,NSIG,SIGBOND, &
                        LBOND,KBOND,LSTART,LREP,LSWITCH,KSWITCH,noldbond,oldbond,oldtyp)


   IMPLICIT NONE

   INTEGER,VALUE:: NBOND,NSIG,noldbond
   INTEGER N,N1,N2,NS

   INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)::BOND,SIGBOND,oldbond
   INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::BONTYP,oldtyp!,SIGBOND

   DOUBLE PRECISION DFX,DFY,DFZ
   DOUBLE PRECISION,VALUE:: KBOND,LBOND,LSTART,LREP
   DOUBLE PRECISION DX,DY,DZ,DIST,F
   DOUBLE PRECISION,VALUE:: KSWITCH,LSWITCH

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FX,FY,FZ
   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z



!---------------------------------------------------------------------------------------
!   THIS HAS BEEN PRE-CALCULATED SO NO NEED TO DO IT AGAIN:
! -- ASSUME THERE IS NO FORCE FOR EXTENSION < CONTOUR LENGTH / 4. USE LREP INSTEAD OF LBOND:
!   LSTART=1.0D0

!   LREP=LBOND-LSTART

! -- ALSO USE ALTERNATIVE BOND CONSTANT:
!      KREP=KBOND*LBOND/LREP

!   LSWITCH=0.9D0*LREP

!   ESWITCH=KBOND*LSWITCH**2*(LREP/4.0D0/(LREP-LSWITCH)+0.5D0)

!   SLOPE=KBOND*(LREP/4.0D0/(1.0D0-LSWITCH/LREP)**2-LREP/4.0D0+LSWITCH)

!   KSWITCH=SLOPE/2.0D0/LSWITCH
!   MSWITCH=ESWITCH-KSWITCH*LSWITCH**2


! BEFORE SWITCHING, POTENTIAL OF PEPTIDE: E(XBOND) = KBOND*XBOND**2*(1/2+1/4(1-XBOND/LBOND))
! AFTER SWITCHING: E(XBOND)=KSWITCH*XBOND**2+MSWITCH
!------------------------------------------------------------------------------------------------


!   ALLOCATE(SIGBOND(2,100)

!  BECAUSE THERE ARE DOUBLE BONDS SO FIRST WE CALCULATE NORMAL BONDS:


!$OMP parallel  &
!$OMP DEFAULT(NONE) &
!$OMP private(N,N1,N2,DX,DY,DZ,DIST,F,DFX,DFY,DFZ) &
!$OMP shared(BOND,BONTYP,NBOND,X,Y,Z,FX,FY,FZ) &
!$OMP SHARED(KBOND,LREP,LSTART,KSWITCH,LSWITCH)

!$OMP DO SCHEDULE(GUIDED,64)

      DO N=1,NBOND

         IF(BONTYP(N)/=1)THEN
            CYCLE
         END IF


         N1=BOND(1,N)
         N2=BOND(2,N)


!       Calculate the distance between two ends of the bond:

         DX=X(N1)-X(N2)
         DY=Y(N1)-Y(N2)
         DZ=Z(N1)-Z(N2)

         DIST=SQRT(DX**2+DY**2+DZ**2)-LSTART


!       Energy of the peptide:

!         IF(DIST>0.0D0.AND.DIST<=LSWITCH)THEN
!            ENG=ENG+KBOND*DIST**2*(LREP/4.0D0/(LREP-DIST)+0.5D0)

!         ELSEIF(DIST>0.0D0.AND.DIST>LSWITCH)THEN
!            ENG=ENG+KSWITCH*DIST**2+MSWITCH
!         END IF

!         IF(JFORCE==0)THEN
!            CYCLE
!         END IF
!-------------------------------------------

         IF(DIST<0.0D0)THEN
            F=0.0D0

         ELSEIF(DIST<=LSWITCH)THEN

            F=KBOND*(LREP/4.0D0/(1.0D0-DIST/LREP)**2-LREP/4.0D0+DIST)

         ELSEIF(DIST>LSWITCH)THEN

            F=2.0D0*KSWITCH*DIST

         END IF

         F=F/(DIST+LSTART)

!         DIST=DIST+LSTART
!         INVDIST=1.0D0/DIST

         DFX=F*DX!*INVDIST
         DFY=F*DY!*INVDIST
         DFZ=F*DZ!*INVDIST

!         FX(N2)=F*DX*INVDIST
!         FY(N2)=F*DY*INVDIST
!         FZ(N2)=F*DZ*INVDIST

!         FX(N1)=-FX(N2)
!         FY(N1)=-FY(N2)
!         FZ(N1)=-FZ(N2)


         FX(N2)=FX(N2)+DFX
         FY(N2)=FY(N2)+DFY
         FZ(N2)=FZ(N2)+DFZ

         FX(N1)=FX(N1)-DFX
         FY(N1)=FY(N1)-DFY
         FZ(N1)=FZ(N1)-DFZ

      END DO

!$OMP END DO NOWAIT
!$OMP end parallel 

!---------------------------------

   if(noldbond==0) goto 10

   do n=1,noldbond

      if(oldtyp(n)==0) cycle

      n1=oldbond(1,n)
      n2=oldbond(2,n)

!       Calculate the distance between two ends of the bond:

         DX=X(N1)-X(N2)
         DY=Y(N1)-Y(N2)
         DZ=Z(N1)-Z(N2)

         DIST=SQRT(DX**2+DY**2+DZ**2)-LSTART

         IF(DIST<0.0D0)THEN
            F=0.0D0

         ELSEIF(DIST<=LSWITCH)THEN

            F=0.1d0*KBOND*(LREP/4.0D0/(1.0D0-DIST/LREP)**2-LREP/4.0D0+DIST)

         ELSEIF(DIST>LSWITCH)THEN

            F=0.2D0*KSWITCH*DIST

         END IF

         F=F/(DIST+LSTART)

         DFX=F*DX!*INVDIST
         DFY=F*DY!*INVDIST
         DFZ=F*DZ!*INVDIST

         FX(N2)=FX(N2)+DFX
         FY(N2)=FY(N2)+DFY
         FZ(N2)=FZ(N2)+DFZ

         FX(N1)=FX(N1)-DFX
         FY(N1)=FY(N1)-DFY
         FZ(N1)=FZ(N1)-DFZ


      END DO


!---------------------------------



10   IF(NSIG==0)THEN
        RETURN
     END IF

      DO N=1,NSIG

         N1=SIGBOND(1,N)
         N2=SIGBOND(2,N)

!       Calculate the distance between two ends of the bond:

         DX=X(N1)-X(N2)
         DY=Y(N1)-Y(N2)
         DZ=Z(N1)-Z(N2)

         DIST=SQRT(DX**2+DY**2+DZ**2)-LSTART


!-------------------------------------------

         IF(DIST<0.0D0)THEN
            F=0.0D0

         ELSEIF(DIST<=LSWITCH)THEN

            F=KBOND*(LREP/4.0D0/(1.0D0-DIST/LREP)**2-LREP/4.0D0+DIST)

         ELSEIF(DIST>LSWITCH)THEN

            F=2.0D0*KSWITCH*DIST

         END IF

         F=F/(DIST+LSTART)

!         DIST=DIST+LSTART
!         INVDIST=1.0D0/DIST

         DFX=F*DX!*INVDIST
         DFY=F*DY!*INVDIST
         DFZ=F*DZ!*INVDIST

         FX(N2)=FX(N2)+DFX
         FY(N2)=FY(N2)+DFY
         FZ(N2)=FZ(N2)+DFZ

         FX(N1)=FX(N1)-DFX
         FY(N1)=FY(N1)-DFY
         FZ(N1)=FZ(N1)-DFZ


      END DO



   END SUBROUTINE

!==================================================================

   subroutine oldbondpush(noldbond,npg,oldbond,oldtyp,pgbelowbond,newpgid,newpglen,wthick,kwall,x,y,z,fy,fz)

   implicit none

   integer,value::noldbond,npg
   integer nb,n1,n2,n,j,j1,j2,check

   integer,allocatable,dimension(:,:),intent(in)::oldbond,newpgid
   integer,allocatable,dimension(:),intent(in)::newpglen,oldtyp
   integer,allocatable,dimension(:)::pgbelowbond

   real(kind=8),value::wthick,kwall
   real(kind=8)::XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3
   real(kind=8)::y0,z0,radius,rad,f,dfy,dfz,xc,yc,zc

   real(kind=8),allocatable,dimension(:),intent(in)::x,y,z

   real(kind=8),allocatable,dimension(:)::fy,fz

   if(noldbond==0.or.npg==0) return


   YP3=0.0D0; ZP3=0.0D0

   do nb=1,noldbond

      if(oldtyp(nb)==0) cycle

      pgbelowbond(nb)=0

      n1=oldbond(1,nb)
      n2=oldbond(2,nb)

      xp1=x(n1)
      yp1=y(n1)
      zp1=z(n1)

      xp2=x(n2)
      yp2=y(n2)
      zp2=z(n2)

      xc=0.5d0*(xp1+xp2)
      yc=0.5d0*(yp1+yp2)
      zc=0.5d0*(zp1+zp2)

      radius=sqrt(yc*yc+zc*zc)-wthick

      xp3=0.5d0*(xp1+xp2)

      xp1=xp1+xp1-xp3
      yp1=yp1+yp1
      zp1=zp1+zp1

      xp2=xp2+xp2-xp3
      yp2=yp2+yp2
      zp2=zp2+zp2

      do n=1,npg

         if(newpglen(n)<2) cycle

         j2=newpgid(1,n)

         xl2=x(j2)
         yl2=y(j2)
         zl2=z(j2)

         do j=2,newpglen(n)

            j1=j2

            xl1=xl2
            yl1=yl2
            zl1=zl2

!            if(abs(yl1-y0)>5.0d0) cycle

!            if(abs(zl1-z0)>5.0d0) cycle

            j2=newpgid(j,n)

!            if(j1==n1.or.j1==n2.or.j2==n1.or.j2==n2) cycle

            xl2=x(j2)
            yl2=y(j2)
            zl2=z(j2)

            if(abs(xl1-xc)>5.0d0) cycle

            if(abs(yl1-yc)>5.0d0) cycle

            if(abs(zl1-zc)>5.0d0) cycle

            if(j1==n1.or.j1==n2.or.j2==n1.or.j2==n2) cycle

            CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

            if(check==0) cycle


            pgbelowbond(nb)=pgbelowbond(nb)+1

            y0=0.5d0*(yl1+yl2)
            z0=0.5d0*(zl1+zl2)

            rad=sqrt(y0*y0+z0*z0)

            if(rad>radius)then

               f=kwall*(rad-radius)/rad

               dfy=f*y0
               dfz=f*z0

               fy(n1)=fy(n1)+dfy
               fz(n1)=fz(n1)+dfz

               fy(n2)=fy(n2)+dfy
               fz(n2)=fz(n2)+dfz

               fy(j1)=fy(j1)-dfy
               fz(j1)=fz(j1)-dfz

               fy(j2)=fy(j2)-dfy
               fz(j2)=fz(j2)-dfz

            end if

         end do

      end do

   end do

   end subroutine


!==================================================================
   subroutine oldbondpush1(nthreads,natomold,natom,noldbond,npg,oldbond,oldtyp, &
                           pgbelowbond,newpgid,newpglen,wthick,kwall,x,y,z,fy,fz)

   implicit none

   integer,value::noldbond,npg,nthreads,natomold,natom
   integer nb,n1,n2,n,j,j1,j2,check,tid,omp_get_thread_num

   integer,allocatable,dimension(:,:),intent(in)::oldbond,newpgid
   integer,allocatable,dimension(:),intent(in)::newpglen,oldtyp
   integer,allocatable,dimension(:)::pgbelowbond

   real(kind=8),value::wthick,kwall
   real(kind=8)::XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3
   real(kind=8)::y0,z0,radius,rad,f,dfy,dfz,xc,yc,zc

   real(kind=8),allocatable,dimension(:),intent(in)::x,y,z

   real(kind=8),allocatable,dimension(:)::fy,fz
   real(kind=8),allocatable,dimension(:,:)::fytem,fztem

!   if(noldbond==0.or.npg==0) return

   allocate(fytem(nthreads,natom-natomold),fztem(nthreads,natom-natomold))

   fytem=0.0d0
   fztem=0.0d0

   YP3=0.0D0; ZP3=0.0D0


!$omp parallel &
!$omp default(none) &
!$omp private (nb,n1,n2,n,j,j1,j2,check,tid) &
!$omp private (XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3) &
!$omp private (y0,z0,radius,rad,f,dfy,dfz,xc,yc,zc) &
!$omp shared(noldbond,npg,oldbond,newpgid,newpglen,oldtyp,pgbelowbond,wthick,kwall,YP3,ZP3) &
!$omp shared(x,y,z,fy,fz,fytem,fztem,natomold)
!$omp do

   do nb=1,noldbond

      if(oldtyp(nb)==0) cycle

      pgbelowbond(nb)=0

      n1=oldbond(1,nb)
      n2=oldbond(2,nb)

      xp1=x(n1)
      yp1=y(n1)
      zp1=z(n1)

      xp2=x(n2)
      yp2=y(n2)
      zp2=z(n2)

      xc=0.5d0*(xp1+xp2)
      yc=0.5d0*(yp1+yp2)
      zc=0.5d0*(zp1+zp2)

      radius=sqrt(yc*yc+zc*zc)-wthick

      xp3=0.5d0*(xp1+xp2)

      xp1=xp1+xp1-xp3
      yp1=yp1+yp1
      zp1=zp1+zp1

      xp2=xp2+xp2-xp3
      yp2=yp2+yp2
      zp2=zp2+zp2

      do n=1,npg

         if(newpglen(n)<2) cycle

         j2=newpgid(1,n)

         xl2=x(j2)
         yl2=y(j2)
         zl2=z(j2)

         do j=2,newpglen(n)

            j1=j2

            xl1=xl2
            yl1=yl2
            zl1=zl2

!            if(abs(yl1-y0)>5.0d0) cycle

!            if(abs(zl1-z0)>5.0d0) cycle

            j2=newpgid(j,n)

!            if(j1==n1.or.j1==n2.or.j2==n1.or.j2==n2) cycle

            xl2=x(j2)
            yl2=y(j2)
            zl2=z(j2)

            if(abs(xl1-xc)>5.0d0) cycle

            if(abs(yl1-yc)>5.0d0) cycle

            if(abs(zl1-zc)>5.0d0) cycle

            if(j1==n1.or.j1==n2.or.j2==n1.or.j2==n2) cycle

            CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

            if(check==0) cycle


            pgbelowbond(nb)=pgbelowbond(nb)+1

            y0=0.5d0*(yl1+yl2)
            z0=0.5d0*(zl1+zl2)

            rad=sqrt(y0*y0+z0*z0)

            if(rad>radius)then

               f=kwall*(rad-radius)/rad

               dfy=f*y0
               dfz=f*z0

               fy(n1)=fy(n1)+dfy
               fz(n1)=fz(n1)+dfz

               fy(n2)=fy(n2)+dfy
               fz(n2)=fz(n2)+dfz

               tid=omp_get_thread_num()+1

               fytem(tid,j1-natomold)=fytem(tid,j1-natomold)-dfy
               fztem(tid,j1-natomold)=fztem(tid,j1-natomold)-dfz

               fytem(tid,j2-natomold)=fytem(tid,j2-natomold)-dfy
               fztem(tid,j2-natomold)=fztem(tid,j2-natomold)-dfz

            end if

         end do

      end do

   end do

!$omp end do
!$omp end parallel


!$omp parallel &
!$omp default(none) &
!$omp private (n,j) &
!$omp shared(natomold,natom,nthreads,fy,fz,fytem,fztem)
!$omp do

   do n=natomold+1,natom

      j=n-natomold

      fy(n)=fy(n)+sum(fytem(1:nthreads,j))
      fz(n)=fz(n)+sum(fztem(1:nthreads,j))

   end do

!$omp end do
!$omp end parallel

   deallocate(fytem,fztem)

   end subroutine




!==================================================================

 SUBROUTINE FPRES(NATOM,FX,FY,FZ,X,Y,Z,NLOOP,LOOP,LOOPLEN,LOOPTYP,PRES,DELTA,INVDELTA,caplen, &
                  npress,pconstrict,invdelx,fycon,fzcon,nsyn)

      IMPLICIT NONE

      INTEGER,VALUE:: NATOM,nsyn
      INTEGER,VALUE:: NLOOP,caplen,npress
      INTEGER I,J,N,NI,NJ,LENGTH,NA,psign(2),jdist,ns

      INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(IN)::LOOPLEN,LOOPTYP
      INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: LOOP

      DOUBLE PRECISION,VALUE:: PRES,DELTA,INVDELTA,invdelx
      DOUBLE PRECISION XCEN,YCEN,ZCEN
      DOUBLE PRECISION X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
      DOUBLE PRECISION V0,VLOOP,VX,VY,VZ,XNOR,YNOR,ZNOR,INVDIST,F,KBOND

      DOUBLE PRECISION dp,dv,turgor,xij

!      DOUBLE PRECISION pconstrict(1000)

      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: FX,FY,FZ,fycon,fzcon
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,INTENT(IN):: X,Y,Z,pconstrict

!      if(jprint==1) allocate(fxcon(natom),fycon(natom),fzcon(natom))

      XCEN=0.0D0!SUM(X(1:NATOM))/NATOM
      YCEN=0.0D0!SUM(Y(1:NATOM))/NATOM
      ZCEN=0.0D0!SUM(Z(1:NATOM))/NATOM

!--   PICK AN ARTIFICIAL BOND CONSTANT TO FLATTEN THE LOOPS:
      KBOND=0.15D0    ! -- THIS VALUE IS CHOSEN SO THE FREQUENCY = 1/50 COMPARED TO K_G

      FX(1:NATOM)=0.0D0; FY(1:NATOM)=0.0D0; FZ(1:NATOM)=0.0D0

      fycon(1:natom)=0.0d0; fzcon(1:natom)=0.0d0

!      sigma2=sigma*sigma

      turgor=pres*invdelta


!      fconstrict=0.0d0


!OMP PARALLEL &
!OMP DEFAULT(NONE) &
!OMP PRIVATE(N,NA,I,NI,NJ,X0,Y0,Z0,V0,V,VLOOP,VX,VY,VZ,XNOR,YNOR,ZNOR,INVDIST,F) &
!OMP PRIVATE(XI,YI,ZI,XJ,YJ,ZJ,LENGTH) &
!OMP SHARED(X,Y,Z,XCEN,YCEN,ZCEN,PRES,DELTA,INVDELTA,NLOOP,LOOP,LOOPLEN) &
!OMP SHARED(FX,FY,FZ,LOOPTYP,KBOND) &

!omp private(dp,da,dfy1,dfy2,dfz1,dfz2) &
!omp shared(pconstrict,sigma2,fycon,fzcon)
!omp reduction(+:fconstrict)

!OMP DO SCHEDULE(GUIDED,64)

   DO N=1,NLOOP

      IF(LOOPTYP(N)==0)THEN
         CYCLE
      END IF

      LENGTH=LOOPLEN(N)

! --- FIND THE LOOP CENTER:

      X0=0.0D0; Y0=0.0D0; Z0=0.0D0

      DO I=1,LENGTH

         NA=LOOP(I,N)

         X0=X0+X(NA)
         Y0=Y0+Y(NA)
         Z0=Z0+Z(NA)

      END DO

      X0=X0/LENGTH; Y0=Y0/LENGTH; Z0=Z0/LENGTH

!     pressure difference due to constriction

!      dp=pconstrict*exp(-x0*x0/sigma2)

! --- THESE ARE USED TO FIND THE LOOP NORMAL VECTOR:
      VLOOP=0.0D0; VX=0.0D0; VY=0.0D0; VZ=0.0D0

! --  NOW RUN THE LOOP:

      NJ=LOOP(LENGTH,N)

      XJ=X(NJ)
      YJ=Y(NJ)
      ZJ=Z(NJ)


!     constriction force

!      fycons=0.0d0
!      fzcons=0.0d0

      DO I=1,LENGTH

         NI=NJ

         XI=XJ
         YI=YJ
         ZI=ZJ


         NJ=LOOP(I,N)

         XJ=X(NJ);YJ=Y(NJ);ZJ=Z(NJ)

!        pressure difference due to constriction

         xij=(xi+xj)

         jdist=abs(xij)*invdelx+1

         if(jdist>npress) jdist=npress

         dp=pconstrict(jdist)

!         do ns=1,nsyn

!            jdist=abs(xij-xsyn(ns))*invdelx+1

!            dp=dp+pconstrict(jdist)

!         end do



!       Volume of the tetrahedron formed by I, J, tetragon
!       center and cell center:

         V0=VOL(XI,YI,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)

! -- TO CALCULATE FORCE ON NI, PLUS ADJUSTMENT DUE TO DIFFERENCE IN RADII FROM CENTER:

!        x component

         dV=VOL(XI+DELTA,YI,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0


         FX(NI)=FX(NI)+turgor*dv!PRES*(V-V0)*INVDELTA


!        y component

         dV=VOL(XI,YI+DELTA,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

!         da=(V-V0)*INVDELTA

!         dfy1=pres*da

!         dfy2=dp*da

         FY(NI)=FY(NI)+(turgor-dp)*dv!dfy1-dfy2!PRES*(V-V0)*INVDELTA

         fycon(ni)=fycon(ni)+dp*dv!dfy2



!        z component

         dV=VOL(XI,YI,ZI+DELTA,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

!         da=(V-V0)*INVDELTA

!         dfz1=pres*da

!         dfz2=dp*da

         FZ(NI)=FZ(NI)+(turgor-dp)*dv!dfz1-dfz2!PRES*(V-V0)*INVDELTA

         fzcon(ni)=fzcon(ni)+dp*dv!dfz2



! --  TO CALCULATE FORCE ON NJ:

!        x component

         dV=VOL(XI,YI,ZI,XJ+DELTA,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FX(NJ)=FX(NJ)+turgor*dv!PRES*(V-V0)*INVDELTA


!        y component

         dV=VOL(XI,YI,ZI,XJ,YJ+DELTA,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

!         da=(V-V0)*INVDELTA

!         dfy1=pres*da

!         dfy2=dp*da


         FY(NJ)=FY(NJ)+(turgor-dp)*dv!dfy1-dfy2!PRES*(V-V0)*INVDELTA

         fycon(nj)=fycon(nj)+dp*dv!dfy2


!        z component

         dV=VOL(XI,YI,ZI,XJ,YJ,ZJ+DELTA,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

!         da=(V-V0)*INVDELTA

!         dfz1=pres*da

!         dfz2=dp*da


         FZ(NJ)=FZ(NJ)+(turgor-dp)*dv!dfz1-dfz2!PRES*(V-V0)*INVDELTA

         fzcon(nj)=fzcon(nj)+dp*dv!dfz2


! --  TO FIND THE NORMAL VECTOR:
         VLOOP=VLOOP+V0

         VX=VX+VOL(XI,YI,ZI,XJ,YJ,ZJ,X0+DELTA,Y0,Z0,XCEN,YCEN,ZCEN)

         VY=VY+VOL(XI,YI,ZI,XJ,YJ,ZJ,X0,Y0+DELTA,Z0,XCEN,YCEN,ZCEN)

         VZ=VZ+VOL(XI,YI,ZI,XJ,YJ,ZJ,X0,Y0,Z0+DELTA,XCEN,YCEN,ZCEN)

      END DO

!      fconstrict=fconstrict+sqrt(fycons*fycons+fzcons*fzcons)


! --  THIS IS THE NORMAL VECTOR:

      XNOR=VX-VLOOP!)*INVDELTA
      YNOR=VY-VLOOP!)*INVDELTA
      ZNOR=VZ-VLOOP!)*INVDELTA

      INVDIST=1.0D0/SQRT(XNOR**2+YNOR**2+ZNOR**2)
!      INVDIST=1.0D0/DIST

      XNOR=XNOR*INVDIST
      YNOR=YNOR*INVDIST
      ZNOR=ZNOR*INVDIST

!---  NOW FLATTEN THE LOOP:

      DO I=1,LENGTH

         NI=LOOP(I,N)

         F=KBOND*(XNOR*(X(NI)-X0)+YNOR*(Y(NI)-Y0)+ZNOR*(Z(NI)-Z0))

!omp atomic
         FX(NI)=FX(NI)-F*XNOR

!omp atomic
         FY(NI)=FY(NI)-F*YNOR

!omp atomic
         FZ(NI)=FZ(NI)-F*ZNOR

      END DO

   END DO

!OMP END DO NOWAIT
!OMP END PARALLEL

!----------------------------------

!  presure on caps

   psign(1)=-1
   psign(2)=1


   do n=1,2

      x0=sum(x(1+(n-1)*caplen:caplen*n))/caplen
      y0=sum(y(1+(n-1)*caplen:caplen*n))/caplen
      z0=sum(z(1+(n-1)*caplen:caplen*n))/caplen

      nj=caplen*n

      XJ=X(NJ)
      YJ=Y(NJ)
      ZJ=Z(NJ)


      do i=1,caplen

         ni=nj

         XI=XJ
         YI=YJ
         ZI=ZJ

         nj=i+(n-1)*caplen

         XJ=X(NJ);YJ=Y(NJ);ZJ=Z(NJ)

!       Volume of the tetrahedron formed by I, J, tetragon

         V0=VOL(XI,YI,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)

! -- TO CALCULATE FORCE ON NI, PLUS ADJUSTMENT DUE TO DIFFERENCE IN RADII FROM CENTER:

         dV=VOL(XI+DELTA,YI,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FX(NI)=FX(NI)+turgor*dv*psign(n)

         dV=VOL(XI,YI+DELTA,ZI,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FY(NI)=FY(NI)+turgor*dv*psign(n)

         dV=VOL(XI,YI,ZI+DELTA,XJ,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FZ(NI)=FZ(NI)+turgor*dv*psign(n)

! --  TO CALCULATE FORCE ON NJ:

         dV=VOL(XI,YI,ZI,XJ+DELTA,YJ,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FX(NJ)=FX(NJ)+turgor*dv*psign(n)

         dV=VOL(XI,YI,ZI,XJ,YJ+DELTA,ZJ,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FY(NJ)=FY(NJ)+turgor*dv*psign(n)

         dV=VOL(XI,YI,ZI,XJ,YJ,ZJ+DELTA,X0,Y0,Z0,XCEN,YCEN,ZCEN)-v0

         FZ(NJ)=FZ(NJ)+turgor*dv*psign(n)

      END DO


   END DO



  END SUBROUTINE FPRES

!======================================

      DOUBLE PRECISION function vol(xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD)
      DOUBLE PRECISION,value:: xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD
      DOUBLE PRECISION x1,y1,z1,x2,y2,z2,x3,y3,z3,x23,y23,z23,a,b,c

      x1=xa-xd
      y1=ya-yd
      z1=za-zd

      x2=xb-xd
      y2=yb-yd
      z2=zb-zd

      x3=xc-xd
      y3=yc-yd
      z3=zc-zd

      x23=y2*z3-y3*z2
      y23=z2*x3-z3*x2
      z23=x2*y3-x3*y2

      vol=abs(x1*x23+y1*y23+z1*z23)/6

      a=(y1-y2)*(z2-z3)-(z1-z2)*(y2-y3)
      b=(z1-z2)*(x2-x3)-(x1-x2)*(z2-z3)
      c=(x1-x2)*(y2-y3)-(y1-y2)*(x2-x3)

      if(a*x1+b*y1+c*z1<0.0d0)then
         vol=-vol
      end if

      end function vol

!======================================

!       This function calculates the volume of a tetrahedron
!       formed of 4 points in space A, B, C and D.

      DOUBLE PRECISION function vol1(xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD)
      DOUBLE PRECISION xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD
      DOUBLE PRECISION Area,H,s
      DOUBLE PRECISION a,b,c,d,delta
      INTEGER ivol

! -- use a very small number to avoid SQRT of negative numbers:

      delta=0.0000000000001d0

!       Consider A, B and C form a base of the tetrahedron which has
!       area denoted as Area.
!       To calculate the area of triangle ABC, we use a,b,c to denote
!       the lengths of BC,AC,AB. 

      a=sqrt((xB-xC)**2+(yB-yC)**2+(zB-zC)**2)
      b=sqrt((xA-xC)**2+(yA-yC)**2+(zA-zC)**2)
      c=sqrt((xB-xA)**2+(yB-yA)**2+(zB-zA)**2)

!       The semiperimeter of the triangle is:

      s=(a+b+c)/2.0D0+delta

!       And here is the area:

      Area=sqrt(s*(s-a)*(s-b)*(s-c))

      if(Area<0.00001D0)then
         ivol=1
         goto 10
      else
         ivol=0
      end if

!       The distance from D to the ABC plane is the height of the
!       tetrahedron H.
!       The ABC plane can be discribed as ax+by+cz+d=0.
!       A normal vector to the plane can be written as n=(a,b,c).
!       n can be calculated as a cross product of vectors AB and BC
!       n = AB cross AC which give:

      a=(yA-yB)*(zB-zC)-(zA-zB)*(yB-yC)
      b=(zA-zB)*(xB-xC)-(xA-xB)*(zB-zC)
      c=(xA-xB)*(yB-yC)-(yA-yB)*(xB-xC)

!       We can find d by using a relation n dot (x-xA) = 0
!       Then here it is:

      d=-a*xA-b*yA-c*zA

!       So the height of the tetrahedron is calculated using a formula
!       in wolfram:

      H=abs((a*xD+b*yD+c*zD+d)/sqrt(a**2+b**2+c**2))

!       Volume of the tetrahedron is Area*Height/3

      vol1=Area*H/3.0D0

!       We want to give volume the sign of n cross DA:

      if(a*(xA-xD)+b*(yA-yD)+c*(zA-zD)<0.0d0)then
         vol1=-vol1
      end if

10    if(ivol==1)then
         vol1=0.0d0
      end if

      end function vol1

!==================================================================

   SUBROUTINE NEWCOOR(NATOM,X,Y,Z,FX,FY,FZ,NSYN,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                   FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE, &
                   FXEDASE,FYEDASE,FZEDASE,INVMPBP,GAM)

   IMPLICIT NONE

   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::XEDASE,YEDASE,ZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::XTPASE,YTPASE,ZTPASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)::FX,FY,FZ,FXEDASE,FYEDASE,FZEDASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::FXGTASE,FYGTASE,FZGTASE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(IN)::FXTPASE,FYTPASE,FZTPASE
   DOUBLE PRECISION,VALUE:: GAM,INVMPBP
   INTEGER NA,NS
   INTEGER,VALUE:: NATOM,NSYN!,JSCALE
!   INTEGER,ALLOCATABLE, DIMENSION(:),INTENT(IN)::SYNLOOP


!$omp parallel &
!$omp default(none) &
!$omp private(na,ns) &
!$omp shared(natom,x,y,z,gam,fx,fy,fz,nsyn,xgtase,ygtase,zgtase,xtpase,ytpase,ztpase) &
!$omp shared(XEDASE,YEDASE,ZEDASE,FXEDASE,FYEDASE,FZEDASE,FXGTASE,FYGTASE,FZGTASE) &
!$omp shared(FXTPASE,FYTPASE,FZTPASE,INVMPBP)
!$omp do schedule(guided,64)

   do na=1,natom
      X(NA)=X(NA)+GAM*FX(NA)
      Y(NA)=Y(NA)+GAM*FY(NA)
      Z(NA)=Z(NA)+GAM*FZ(NA)
   end do

!$omp enddo nowait

!$omp do

   DO NS=1,NSYN

!      IF(SYNLOOP(NS)==0)THEN
!         JG=JSCALE*5
!      ELSE
!         JG=1
!      END IF

      XGTASE(1:2,NS)=XGTASE(1:2,NS)+GAM*FXGTASE(1:2,NS)*INVMPBP!*JG
      YGTASE(1:2,NS)=YGTASE(1:2,NS)+GAM*FYGTASE(1:2,NS)*INVMPBP!*JG
      ZGTASE(1:2,NS)=ZGTASE(1:2,NS)+GAM*FZGTASE(1:2,NS)*INVMPBP!*JG

      XTPASE(1:3,NS)=XTPASE(1:3,NS)+GAM*FXTPASE(1:3,NS)*INVMPBP!*JG
      YTPASE(1:3,NS)=YTPASE(1:3,NS)+GAM*FYTPASE(1:3,NS)*INVMPBP!*JG
      ZTPASE(1:3,NS)=ZTPASE(1:3,NS)+GAM*FZTPASE(1:3,NS)*INVMPBP!*JG

!      XEDASEOLD(NS)=XEDASE(NS)
!      YEDASEOLD(NS)=YEDASE(NS)
!      ZEDASEOLD(NS)=ZEDASE(NS)

      XEDASE(NS)=XEDASE(NS)+GAM*FXEDASE(NS)*INVMPBP!*JG
      YEDASE(NS)=YEDASE(NS)+GAM*FYEDASE(NS)*INVMPBP!*JG
      ZEDASE(NS)=ZEDASE(NS)+GAM*FZEDASE(NS)*INVMPBP!*JG


   ENDDO

!$omp enddo nowait

!$omp end parallel


   END SUBROUTINE

!==================================================================

   SUBROUTINE deactivate(ns,nsyn,jstep,nbondpep,nbonddel,jdeact,newpglen,looplen,bontyp,dnor,ator,gtload,looptyp, &
              edcap,sigcleave,jtrap,synloop,synpg,newpgid,loop,bondpep,synthesis,glytip,glysec,crlkage, &
              sigcross,edhold,edpep,tppep)

   implicit none

   integer,value::ns,nsyn
   integer(kind=8),value::jstep
   integer nbondpep,nbonddel
   integer iloop,ipg,ns1
   integer, allocatable, dimension(:)::jdeact,newpglen,looplen,bontyp,dnor,ator,gtload,looptyp,edcap,sigcleave,synloop
   integer(kind=8),allocatable, dimension(:)::jtrap
   integer, allocatable, dimension(:,:)::synpg,newpgid,loop,bondpep,synthesis,glytip,glysec,crlkage
   integer, allocatable, dimension(:,:)::sigcross,edhold,edpep
   integer, allocatable, dimension(:,:,:)::tppep


            ILOOP=SYNLOOP(NS)

            IPG=SYNPG(1,NS)

            IF(IPG>0)THEN

               CALL PGCLEAN(IPG,newPGID,newPGLEN,ILOOP,LOOP,LOOPLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,DNOR,ATOR)

            END IF

            IPG=SYNPG(2,NS)

            IF(IPG>0)THEN

               CALL PGCLEAN(IPG,newPGID,newPGLEN,ILOOP,LOOP,LOOPLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,DNOR,ATOR)

            END IF

            JDEACT(NS)=0

            SYNTHESIS(1:2,NS)=0


            SYNPG(1:2,NS)=0

            GLYTIP(1:2,NS)=0

            GLYSEC(1:2,NS)=0

            GTLOAD(NS)=0

            CRLKAGE(1:2,NS)=0

            IF(ILOOP>0)THEN

               LOOPTYP(ILOOP)=1

               DO NS1=1,NSYN

                  IF(NS1/=NS.AND.SYNLOOP(NS1)==ILOOP)THEN
                     LOOPTYP(ILOOP)=2
                     EXIT
                  END IF

               END DO

            END IF

            SYNLOOP(NS)=0

            TPPEP(1:2,1,NS)=0

            TPPEP(1:2,2,NS)=0

            TPPEP(1:2,3,NS)=0

            SIGCROSS(1:3,NS)=0

            EDHOLD(1:3,NS)=0

            EDPEP(1:2,NS)=0

            EDCAP(NS)=0

            SIGCLEAVE(NS)=0


            JTRAP(NS)=JSTEP

 

   END SUBROUTINE

!==================================================================

   SUBROUTINE LOOPCHECK(NSYN,XGTASE,YGTASE,ZGTASE,X,Y,Z,SYNTHESIS,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,JDEACT)

   IMPLICIT NONE

   INTEGER NA1,NA2,NSYN,NS,JCHECK,CHECK,J,ILOOP,js

   INTEGER, ALLOCATABLE, DIMENSION(:)::SYNLOOP,LOOPLEN,LOOPTYP,JDEACT
   INTEGER, ALLOCATABLE,INTENT(IN), DIMENSION(:,:)::LOOP,SYNTHESIS

   DOUBLE PRECISION XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3

   DOUBLE PRECISION,ALLOCATABLE,INTENT(IN), DIMENSION(:)::X,Y,Z
   DOUBLE PRECISION,ALLOCATABLE,INTENT(IN), DIMENSION(:,:)::XGTASE,YGTASE,ZGTASE


   DO NS=1,NSYN

      IF(SYNTHESIS(1,NS)==0.OR.SYNTHESIS(2,NS)==0)THEN
         CYCLE
      END IF

      ILOOP=SYNLOOP(NS)

      XP3=X(LOOP(1,ILOOP))
      YP3=Y(LOOP(1,ILOOP))
      ZP3=Z(LOOP(1,ILOOP))

      do js=1,2

         XL1=XGTASE(js,NS)
         YL1=2*YGTASE(js,NS)
         ZL1=2*ZGTASE(js,NS)

         XL2=XL1; YL2=0.0D0; ZL2=0.0D0

         JCHECK=0

         DO J=2,LOOPLEN(ILOOP)-1

            NA1=LOOP(J,ILOOP)

            NA2=LOOP(J+1,ILOOP)

            IF(LOOP(1,ILOOP)==NA1.OR.LOOP(1,ILOOP)==NA2)THEN
               CYCLE
            END IF

            XP1=X(NA1); YP1=Y(NA1); ZP1=Z(NA1)

            XP2=X(NA2); YP2=Y(NA2); ZP2=Z(NA2)

            CALL INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

            IF(CHECK==1)THEN
               JCHECK=JCHECK+1
            END IF

         END DO

         IF(MOD(JCHECK,2)==0)THEN
            PRINT*,'out loop for complex',NS,js

            JDEACT(NS)=1


            exit
         END IF

      end do

   END DO

   END SUBROUTINE

!==================================================================

 SUBROUTINE PG_OUT(NSTART,NATOMOLD,NATOM,NATOMDEL,OLDNATOMDEL,DNOR,ATOR,PEPDIR,X,Y,Z, &
              NPG,NPGOLD,PGID,PGLEN,PGLENMAX,lenmax,newPGid,newpglen,NBONDPEP,NBONDDEL,noldbond,oldbond,oldtyp, &
               BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP,caplen, &
                jsyndir,na_red,NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDCAP,EDHOLD,CRLKAGE, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD, &
                  oldglytip,oldglysec,temlink,ntemlink)


   IMPLICIT NONE

   INTEGER,VALUE:: NATOMOLD,NATOM,NSYN,NATOMDEL,OLDNATOMDEL,NPG
   INTEGER,VALUE:: NBONDPEP,NBONDDEL,NPGOLD,noldbond,ntemlink
   INTEGER,VALUE:: NLOOP,NLOOPDEL,PGLENMAX,lenmax,caplen,jsyndir,na_red
   INTEGER N,N0,I,J,K,LOOPLENMAX,NSTART

   INTEGER,DIMENSION(:,:),ALLOCATABLE::PGID,newpgid,BONDPEP,LOOP,SYNTHESIS,SYNPG,GLYTIP,GLYSEC,oldbond
   INTEGER,DIMENSION(:),ALLOCATABLE::PGLEN,newpglen,BONTYP,DNOR,ATOR,PEPDIR,oldtyp,GTLOAD,GTATRANS
   INTEGER,DIMENSION(:,:),ALLOCATABLE::SIGCROSS,EDHOLD,CRLKAGE,EDPEP,oldglytip,oldglysec,temlink
   INTEGER,DIMENSION(:),ALLOCATABLE::LOOPLEN,LOOPTYP,SYNDIR,SYNLOOP,JDEACT,SIGCLEAVE,EDCAP
   INTEGER,DIMENSION(:,:,:),ALLOCATABLE::TPPEP

   DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: X,Y,Z,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: XTPASE,YTPASE,ZTPASE,XGTASE,YGTASE,ZGTASE
!   DOUBLE PRECISION ORAD

   CHARACTER CHARA*512
   CHARACTER (LEN=64) FILECOOR,FILECONFIG
   CHARACTER ZERO*1,CHARID1*1,CHARID2*2,CHARID3*3,CHARID4*4

   NSTART=NSTART+1

   WRITE(ZERO,'(I1)')0

   IF(NSTART<10)THEN
      WRITE(CHARID1,'(I1)')NSTART
      FILECOOR='coor'//ZERO//ZERO//ZERO//CHARID1//'.dat'
      FILECONFIG='config'//ZERO//ZERO//ZERO//CHARID1//'.dat'
   ELSE IF(NSTART<100)THEN
      WRITE(CHARID2,'(I2)')NSTART
      FILECOOR='coor'//ZERO//ZERO//CHARID2//'.dat'
      FILECONFIG='config'//ZERO//ZERO//CHARID2//'.dat'
   ELSE IF(NSTART<1000)THEN
      WRITE(CHARID3,'(I3)')NSTART
      FILECOOR='coor'//ZERO//CHARID3//'.dat'
      FILECONFIG='config'//ZERO//CHARID3//'.dat'
   ELSE IF(NSTART<10000)THEN
      WRITE(CHARID4,'(I4)')NSTART
      FILECOOR='coor'//CHARID4//'.dat'
      FILECONFIG='config'//CHARID4//'.dat'
   ELSE
      PRINT*,'NSTART IS TOO BIG! STOP NOW.'
      STOP
   END IF

!  ----------------------------------------------
   OPEN(1,FILE=FILECOOR,FORM='UNFORMATTED')

   WRITE(1)NATOM

   WRITE(1)X(1:NATOM)
   WRITE(1)Y(1:NATOM)
   WRITE(1)Z(1:NATOM)

   WRITE(1)NSYN

   WRITE(1)XGTASE(1,1:NSYN),XGTASE(2,1:NSYN),XTPASE(1,1:NSYN),XTPASE(2,1:NSYN),XTPASE(3,1:NSYN),XEDASE(1:NSYN),XEDASEOLD(1:NSYN)
   WRITE(1)YGTASE(1,1:NSYN),YGTASE(2,1:NSYN),YTPASE(1,1:NSYN),YTPASE(2,1:NSYN),YTPASE(3,1:NSYN),YEDASE(1:NSYN),YEDASEOLD(1:NSYN)
   WRITE(1)ZGTASE(1,1:NSYN),ZGTASE(2,1:NSYN),ZTPASE(1,1:NSYN),ZTPASE(2,1:NSYN),ZTPASE(3,1:NSYN),ZEDASE(1:NSYN),ZEDASEOLD(1:NSYN)


!  constriction pressure

!   write(1)pconstrict(1:npress)

   CLOSE(1)

!  ----------------------------------------------

   OPEN(1,FILE=FILECONFIG)

   WRITE(1,*)'ATOM INPUT'

   WRITE(1,*)NATOMOLD,NATOM,NATOMDEL,OLDNATOMDEL
   WRITE(1,*)DNOR(1:NATOM)
   WRITE(1,*)ATOR(1:NATOM)
   WRITE(1,*)PEPDIR(1:NATOM)
   WRITE(1,*)

!  ----------------------------------------------

   WRITE(1,*)'CAPS INPUT'

   WRITE(1,*)caplen

!   WRITE(1,*)capopen(1:caplen)
!   WRITE(1,*)capopen(1+caplen:caplen*2)

   WRITE(1,*)

   WRITE(1,*)'OLDPG INPUT'

!   PGLENMAX=MAXVAL(PGLEN(1:NPG))

   WRITE(1,*)NPGOLD,PGLENMAX

   DO I=1,NPGold
      WRITE(1,*)PGLEN(I),PGID(1:PGLEN(I),I)
   END DO
   WRITE(1,*)

   WRITE(1,*)'NEWPG INPUT'

   WRITE(1,*)NPG,LENMAX

   DO I=1,NPG
      WRITE(1,*)newPGLEN(I),newPGID(1:newPGLEN(I),I)
   END DO
   WRITE(1,*)

!  ----------------------------------------------


!  ----------------------------------------------

   WRITE(1,*)'BONDPEP INPUT'

   WRITE(1,*)NBONDPEP,NBONDDEL

11   FORMAT(I6,4X,I1,4X,I6,4X,I6)

   DO N=1,NBONDPEP
      WRITE(1,11)N,BONTYP(N),BONDPEP(1,N),BONDPEP(2,N)
   END DO

   write(1,*)noldbond
   do n=1,noldbond
      write(1,*)n,oldtyp(n),oldbond(1,n),oldbond(2,n)
   end do 

   WRITE(1,*)

!  ----------------------------------------------

   WRITE(1,*)'LOOP INPUT'

   LOOPLENMAX=MAXVAL(LOOPLEN(1:NLOOP))

   WRITE(1,*)NLOOP,LOOPLENMAX,NLOOPDEL

   DO N=1,NLOOP
      WRITE(1,*)LOOPLEN(N),LOOPTYP(N),LOOP(1:LOOPLEN(N),N)
   END DO
   WRITE(1,*)

!  ----------------------------------------------

!   WRITE(1,*)'ORIGINAL RADIUS'
!   WRITE(1,*)ORAD

!  ----------------------------------------------

12 FORMAT(2X,I2,2(2X,I1),10(2X,I7))

   WRITE(1,*)'SYNCOMP INPUT'

      WRITE(1,*)NSYN!,SYNRATIO   ! THIS ARE # OF SYNTHESIS COMPLEXES AND THE ATOM/SYNTHESIS RATIO

!      WRITE(1,*)'COUPLING(1:NSYN)'
!      WRITE(1,*)COUPLING(1:NSYN)

      WRITE(1,*)'SYNDIR(1:NSYN)'
      WRITE(1,*)SYNDIR(1:NSYN)

      WRITE(1,*)'SYNTHESIS(1,1:NSYN),SYNTHESIS(2,1:NSYN)'
      WRITE(1,*)SYNTHESIS(1,1:NSYN),SYNTHESIS(2,1:NSYN)

      WRITE(1,*)'GTLOAD(1,1:NSYN)'
      WRITE(1,*)GTLOAD(1:NSYN)!,GTLOAD(2,1:NSYN)

      WRITE(1,*)'SYNPG(1,1:NSYN),SYNPG(2,1:NSYN)'
      WRITE(1,*)SYNPG(1,1:NSYN),SYNPG(2,1:NSYN)

      WRITE(1,*)'SYNLOOP(1:NSYN)'
      WRITE(1,*)SYNLOOP(1:NSYN)

      WRITE(1,*)'GLYTIP(1,1:NSYN),GLYTIP(2,1:NSYN)'
      WRITE(1,*)GLYTIP(1,1:NSYN),GLYTIP(2,1:NSYN)

      WRITE(1,*)'GLYSEC(1,1:NSYN),GLYSEC(2,1:NSYN)'
      WRITE(1,*)GLYSEC(1,1:NSYN),GLYSEC(2,1:NSYN)

      WRITE(1,*)'TPPEP(1,1,1:NSYN),TPPEP(1,2,1:NSYN),TPPEP(1,3,1:NSYN)'
      WRITE(1,*)TPPEP(1,1,1:NSYN),TPPEP(1,2,1:NSYN),TPPEP(1,3,1:NSYN)

      WRITE(1,*)'TPPEP(2,1,1:NSYN),TPPEP(2,2,1:NSYN),TPPEP(2,3,1:NSYN)'
      WRITE(1,*)TPPEP(2,1,1:NSYN),TPPEP(2,2,1:NSYN),TPPEP(2,3,1:NSYN)


      WRITE(1,*)'EDPEP(1,1:NSYN),EDPEP(2,1:NSYN)'
      WRITE(1,*)EDPEP(1,1:NSYN),EDPEP(2,1:NSYN)

      WRITE(1,*)'GTATRANS(1:NSYN)'
      WRITE(1,*)GTATRANS(1:NSYN)!,GTATRANS(2,1:NSYN)

      WRITE(1,*)'JDEACT(1:NSYN)'
      WRITE(1,*)JDEACT(1:NSYN)

      WRITE(1,*)'SIGCROSS(1,1:NSYN),SIGCROSS(2,1:NSYN),SIGCROSS(3,1:NSYN)'
      WRITE(1,*)SIGCROSS(1,1:NSYN),SIGCROSS(2,1:NSYN),SIGCROSS(3,1:NSYN)

      WRITE(1,*)'SIGCLEAVE(1:NSYN)'
      WRITE(1,*)SIGCLEAVE(1:NSYN)

!      WRITE(1,*)'EDLOCKIN(1:NSYN)'
!      WRITE(1,*)EDLOCKIN(1:NSYN)

      WRITE(1,*)'EDCAP(1:NSYN)'
      WRITE(1,*)EDCAP(1:NSYN)

      WRITE(1,*)'EDHOLD(1,1:NSYN),EDHOLD(2,1:NSYN),EDHOLD(3,1:NSYN)'
      WRITE(1,*)EDHOLD(1,1:NSYN),EDHOLD(2,1:NSYN),EDHOLD(3,1:NSYN)

      WRITE(1,*)'CRLKAGE(1,1:NSYN),CRLKAGE(2,1:NSYN)'
      WRITE(1,*)CRLKAGE(1,1:NSYN),CRLKAGE(2,1:NSYN)

      WRITE(1,*)


   WRITE(1,*)'TEMLINK'
   write(1,*)ntemlink
   do n=1,ntemlink
      write(1,*)temlink(1:4,n)
   end do

   write(1,*)oldglytip(1,1:nsyn),oldglytip(2,1:nsyn),oldglysec(1,1:nsyn),oldglysec(2,1:nsyn)

      WRITE(1,*)

   WRITE(1,*)'REDIStributed point'
   write(1,*)jsyndir,na_red   

   CLOSE(1)

   END SUBROUTINE

!==========================================

   SUBROUTINE FIXPEPDIR(NPG,PGID,PGLEN,PEPDIR)

   IMPLICIT NONE

   INTEGER,VALUE::NPG

   INTEGER N,J,NA,JCHECK,J1,NA1,DIRCHECK,JDEL,NCOUNT

   INTEGER,ALLOCATABLE,DIMENSION(:)::PGLEN,PEPDIR
   INTEGER,ALLOCATABLE,DIMENSION(:,:)::PGID

   NCOUNT=0

   DO N=1,NPG

      DO J=1,PGLEN(N)

         NA=PGID(J,N)

         IF(PEPDIR(NA)==0)THEN

            NCOUNT=NCOUNT+1

            JCHECK=0

            DO J1=1,PGLEN(N)

               NA1=PGID(J1,N)

               IF(PEPDIR(NA1)/=0)THEN
                  DIRCHECK=PEPDIR(NA1)
                  JCHECK=J1
                  EXIT
               END IF

            END DO

            IF(JCHECK==0)THEN
               PRINT*,'ERROR AT FIX PEPDIR AT',N,NA
               STOP
            END IF

            JDEL=ABS(J-JCHECK)

            IF(MOD(JDEL,2)==0)THEN
               PEPDIR(NA)=DIRCHECK
            ELSE
               PEPDIR(NA)=-DIRCHECK
            END IF

         END IF

      END DO

   END DO

   PRINT*,'NUMBER OF PEPDIR FIXED:',NCOUNT

   END SUBROUTINE


!==========================================

END MODULE
