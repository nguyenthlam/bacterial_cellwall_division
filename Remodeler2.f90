PROGRAM PG_SYNTHESIS

USE DECLARE
USE MODS_GROW
!use ieee_arithmetic
IMPLICIT NONE
integer omp_get_thread_num,omp_get_num_threads,tid,jsig


   call paras(nthreads,nsyn,n_edhold,ntrap,lenmax,newlengthave,oldtypmax,joldbond,k_g,l_g,ktheta,theta_0,k_p,l_p, &
              LSTART,LREP,LSWITCH,KSWITCH,pres,pi,delta,invdelta,beta,ksur,kgtase,ktpase,ltpase, &
              dreact,dreact2,kedase,ledase,l_edcap,kgttp,lgttp13,lgttp2,kgted,lgted,kside,lside,ktped,ltped, &
              kwall,klead,kpair,lpair,invmpbp,wthick,pthick,kaxis,sigma,pconstrict0,delx,invdelx,lcons, &
              lambda,prmbond,growth,p_acti,p_deact,pterm0,p_glyin,p_glyout,p_transfast,nprint, &
              p_transslow,p_pepout,p_edhold,pcleave_spon,pedcap_spon,pedcap_inact,plyt,pgbelowbondmax)



   call omp_set_num_threads(nthreads)

!=============================================================================
!  TIMESTART is the clock time used to calculate the running time.
!  NSTART is the index of starting stage.
!  JFILE tells the program to write DCD files starting with index JFILE+1
!  LASTTIME is the total running time up to the stage of NSTART

   CALL system_clock(count_rate=cr)
   rate = REAL(cr)
   call system_clock(timestart)

   CALL INIT_RANDOM_SEED()

   OPEN(10,FILE='restart.dat')
   read(10,*)
   READ(10,*)NSTART,JFILE,jstart(nstart+1),LASTTIME,ratio(nstart+1)
   CLOSE(10)



!=============================================================================
!  SETUP THE INITIAL SYSTEM:

!  NATOM is the number of "atoms" in the system. Each atom is a PG unit.
!  X,Y,Z are the coordinates of atoms
!  PEPDIR tells the direction of peptide of the PG unit, it can be + or -
!  DNOR and ATOR tell donor and acceptor statuses
!  NPG is the number of residues. Each residue is a glycan strand
!  PGID denotes index of atoms on residues
!  PGLEN tells how long a strand is
!  PGDIR tells direction of strand
!  PGTYP tells if a strand is at 2 caps, old or new.
!  NBONDPEP is the number of peptide bonds, or peptide cross-links
!  BONDPEP denotes two PG unit indices associated with the bond
!  BONTYP tells if the bond is mature, immature or deleted
!  NBONDGLY is the number of glycosydic bonds, each connecting 2 PG units on the same strand
!  BONDGLY denotes the 2 PG unit indices
!  NLOOP is the number of "loop". Each loop is a complete circling of bonds
!  LOOP denotes the indices of atoms on a loop
!  LOOPLEN tells how many atoms the loop is consisted of
!  LOOPTYP tells if a loop is associated with a synthetic complex, non-associated, or deleted
!  SYNTHESIS tell status of a synthetic complex as active or inactive
!  SYNDIR dictates the direction of synthesis, or orientation of transglycosylase (PBP1)
!  GTLOAD tells if PBP1 (GTASE) is loaded with a precursor
!  SYNPG denotes the index of strand the complex is synthesizing. If value is 0: no strand
!  SYNLOOP denotes the active loop
!  GTASE denotes transglycosylase
!  TPASE denotes transpeptidase
!  EDASE denotes endopeptidase
!  SYNRAD denotes the local radius, or the distance from the surface to the cell axis at PBP1 position
!  ATOMRAD denotes the local radius, or the distance from the surface to the cell axis at an atom position
!  GLYTIP is the tip of the strand
!  GLYSEC is the second PG unit from the tip.
!  TPPEP is the units held by TPASE
!  EDPEP is the units held by EDASE





!  to speed up, only evaluate enzymes every 10 steps.

   JSKIP=10

   P_ACTI=P_ACTI*JSKIP
   PTERM=PTERM*JSKIP
   P_TRANSFAST=P_TRANSFAST*JSKIP
   P_TRANSSLOW=P_TRANSSLOW*JSKIP
   P_GLYIN=P_GLYIN*JSKIP
   P_GLYOUT=P_GLYOUT*JSKIP
   P_PEPOUT=P_PEPOUT*JSKIP
   P_EDHOLD=P_EDHOLD*JSKIP
   PCLEAVE_SPON=PCLEAVE_SPON*JSKIP
   PEDCAP_SPON=PEDCAP_SPON*JSKIP
   PEDCAP_INACT=PEDCAP_INACT*JSKIP


!------------------------------------------------------

!  LIMIT PG LENGTH NOT TO BE MORE THAN:


   CALL GETINFO(NSTART,NATOMOLD,NPGOLD,PGLENMAX,LOOPLENMAX,caplen)


   NATOMNEW=NATOMOLD*GROWTH+100



   NPGNEW=NATOMNEW/newlengthave!NPGOLD*GROWTH*1.1d0+10


   NATOMMAX=NATOMOLD+NATOMNEW

   NPEPMAX=NATOMMAX

   noldmax=npepmax/10

   NGLYMAX=NATOMMAX

   NLOOPMAX=NATOMMAX/2

   LOOPLENMAX=MAX(LOOPLENMAX*2,60)

   NTOTAL=NATOMold+2*natomnew+2*NPEPMAX+6*nsyn+ 4*nsyn + 12*nsyn+  2*nsyn    +4*nsyn+2*noldmax



!  the number of intervals for constriction pressure

   npress=lcons/delx
!-------------------------------------------------------

!  allocation of memory

   ALLOCATE(X(NATOMMAX),Y(NATOMMAX),Z(NATOMMAX),DNOR(NATOMMAX),ATOR(NATOMMAX),PEPDIR(NATOMMAX))
   ALLOCATE(PGID(PGLENMAX,NPGOLD),PGLEN(NPGOLD),newpgid(lenmax,npgnew),newpglen(npgnew))
   ALLOCATE(BONDGLY(2,NGLYMAX),BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX),PARTNER(2,NATOMMAX))
   ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX),LOOPLEN(NLOOPMAX),LOOPTYP(NLOOPMAX))


   ALLOCATE(SYNDIR(nsyn),SYNTHESIS(2,nsyn), &
            GTLOAD(nsyn),SYNPG(2,nsyn),SYNLOOP(nsyn))

   ALLOCATE(XGTASE(2,nsyn),YGTASE(2,nsyn),ZGTASE(2,nsyn))
   ALLOCATE(XTPASE(3,nsyn),YTPASE(3,nsyn),ZTPASE(3,nsyn))
   ALLOCATE(XEDASE(nsyn),YEDASE(nsyn),ZEDASE(nsyn))
   ALLOCATE(XEDASEOLD(nsyn),YEDASEOLD(nsyn),ZEDASEOLD(nsyn))

   ALLOCATE(GLYTIP(2,nsyn),GLYSEC(2,nsyn),TPPEP(2,3,nsyn),EDPEP(2,nsyn))

   allocate(oldglytip(2,nsyn),oldglysec(2,nsyn))

   ntemlinkmax=400/newlengthave
   print*,'hoop size is assumed 400 beads, change this step if otherwise'

   allocate(temlink(4,ntemlinkmax))



   allocate(oldbond(2,noldmax),oldtyp(noldmax),pgbelowbond(noldmax))


!-----------------------------------------------

   ALLOCATE(GTATRANS(nsyn))

   ALLOCATE(JDEACT(nsyn))

   ALLOCATE(JTERM(2,nsyn))

   ALLOCATE(SIGCROSS(3,nsyn))

   ALLOCATE(SIGCLEAVE(nsyn))


   ALLOCATE(EDCAP(nsyn))

   ALLOCATE(SIGBOND(2,200))


!  THIS IS USED TO TRACK PG UNITS AND BONDS IN ASSOCIATION WITH ENDOPEPTIDASE:

   ALLOCATE(EDHOLD(3,nsyn))

   ALLOCATE(J_EDHOLD(nsyn))

!  TRACK CROSS-LINKAGE OF STRAND ASSOCIATED WITH THE COMPLEXES:

   ALLOCATE(CRLKAGE(2,nsyn))

!  TO TELL IF A COMPLEX IS TRAPPED:

   ALLOCATE(JTRAP(nsyn))


!---------------------------------
!  FORCES ON ATOMS AND PBPS:

   ALLOCATE(FX(NATOMMAX),FY(NATOMMAX),FZ(NATOMMAX),FAMAG(NATOMMAX))
   ALLOCATE(FXREP(NATOMMAX),FYREP(NATOMMAX),FZREP(NATOMMAX))
   ALLOCATE(FXGLY(NATOMMAX),FYGLY(NATOMMAX),FZGLY(NATOMMAX))
   ALLOCATE(FXPEP(NATOMMAX),FYPEP(NATOMMAX),FZPEP(NATOMMAX))
   ALLOCATE(FXTHETA(NATOMMAX),FYTHETA(NATOMMAX),FZTHETA(NATOMMAX))
   ALLOCATE(FXPRES(NATOMMAX),FYPRES(NATOMMAX),FZPRES(NATOMMAX))

!   allocate(xsyn(nsyn))

   allocate(fycon(natommax),fzcon(natommax))


   ALLOCATE(FXGTASE(2,nsyn),FYGTASE(2,nsyn),FZGTASE(2,nsyn))
   ALLOCATE(FXTPASE(3,nsyn),FYTPASE(3,nsyn),FZTPASE(3,nsyn))
   ALLOCATE(FXEDASE(nsyn),FYEDASE(nsyn),FZEDASE(nsyn),FMAGPBP(nsyn))


   NRAND=1000000

   JRESET=NRAND/5

!   NSYNREP=MIN(NSYN+2,nsyn)
   ALLOCATE(RXGTASE(2*NRAND),RYGTASE(2*NRAND),RZGTASE(2*NRAND))
   ALLOCATE(RXTPASE(3*NRAND),RYTPASE(3*NRAND),RZTPASE(3*NRAND))
   ALLOCATE(RXEDASE(NRAND),RYEDASE(NRAND),RZEDASE(NRAND))


!  SETS OF RANDOM NUMBERS FOR CONVENIENCE:

   ALLOCATE(RANDS(NRAND))
!   JRANDS=NRAND

!-------------------

   ALLOCATE(LGTASE(2,nsyn))!,XLEAD(nsyn),YLEAD(nsyn),ZLEAD(nsyn))

   ALLOCATE(GLYNUM(nsyn),NEIGLY(2,100,nsyn))


!  constriction pressure

   allocate(pconstrict(npress))

!------------------------------------------------

!  reset arrays

   newPGID=0 ! this one is reset only once

35 PARTNER=0
   SYNTHESIS=0; GTLOAD=0; SYNPG=0; SYNLOOP=0
   GLYTIP=0; GLYSEC=0; TPPEP=0; EDPEP=0
   noldbond=0
   pgbelowbond=0
   GTATRANS=0
   JDEACT=0
   JTERM=0
   SIGCROSS=0
   SIGCLEAVE=0
   EDCAP=0
   SIGBOND=0
   EDHOLD=0
   J_EDHOLD=0
   CRLKAGE=0
   JTRAP=0

   FXREP=0.0D0; FYREP=0.0D0; FZREP=0.0D0
   FXGLY=0.0D0; FYGLY=0.0D0; FZGLY=0.0D0
   FXTHETA=0.0D0; FYTHETA=0.0D0; FZTHETA=0.0D0
   FXPEP=0.0D0; FYPEP=0.0D0; FZPEP=0.0D0
   FXPRES=0.0D0; FYPRES=0.0D0; FZPRES=0.0D0

   LGTASE=0.5D0
   JRANDS=NRAND
   JRFORCE=NRAND

   pconstrict=0.0!(1:1000)=pconstrictold(1:1000)

   ntemlink=0

   temlink=0

   oldglytip=0

   oldglysec=0

   jsyndir=0

   na_red=0

   CALL PG_INPUT(NSTART,NATOMOLD,NATOM,NATOMDEL,OLDNATOMDEL,DNOR,ATOR,PEPDIR,X,Y,Z, &
         NPGOLD,PGID,PGLEN,NPG,NEWPGID,NEWPGLEN,NBONDGLY,NBONDGLYOLD,NBONDPEP,NBONDDEL,caplen, &
          noldbond,oldbond,oldtyp,PARTNER,BONDGLY,BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
           jsyndir,na_red,NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
            GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDCAP,EDHOLD,CRLKAGE, &
             XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD, &
              npress,pconstrict,delx,pconstrict0,sigma,invdelta,oldglytip,oldglysec,temlink,ntemlink)


   IF(NSTART==0)THEN

      CALL VISUALPSF(npgold,NATOMOLD,NATOMnew,NPEPMAX,nsyn,noldmax,pglen,pgid)

   END IF



!   XEDASEOLD(1:NSYN)=XEDASE(1:NSYN)
!   YEDASEOLD(1:NSYN)=YEDASE(1:NSYN)
!   ZEDASEOLD(1:NSYN)=ZEDASE(1:NSYN)

!  NATOMCAP is the number of atoms on two caps
!  NATOMDEL is the number of atoms removed from the sacculus by lytic transglycosylases
!  PGCAP1, PGCAP2 are indices of 2 residues (hoops) connecting 2 caps with the cylinder
!  NATOMSTART is the number of atoms at NSTART=0, this is used for visualization purpose
!  NPGSTART is the number of residues (strands) at NSTART=0
!  L_G is the relaxed distance between 2 adjacent PG units on the same strand
!  K_G is the stretching constant of strands
!  L_P is the contour length of a peptide cross-link
!  K_P is the force constant in worm-like chain model for cross-links
!  THETA_0 is the relaxed angle at each unit on strands
!  KTHETA is the bending constant of strands
!  PRES is the turgor pressure

!  JUNIT is used to write DCD files
!  NTOTAL is the max number of atoms in visualization including glycan atoms, peptide atoms,
!  PBPs atoms, cross-linking atoms, deleted bond atoms, etc, anything that helps visualization
!  NATOMTOTAL is the number of glycan atoms in visualization
!  NPGTOTAL is the number of residues in visualization
!  MAXLEN is the max length of new strands in visualization
!  NPEPMAX is number of peptide cross-links in visualization
!  JFILE is the index of DCD files
!  NFRAME is the index of a frame in a DCD file


!-- WRITE DCD FILE:

   CALL RANDOM_NUMBER(R)
   JUNIT=80*R+11
   JFILE=JFILE+1

   CALL DCDHEADER(JUNIT,JFILE,NTOTAL)

   NFRAME=0


!============================================================


!=========================================================================
!  MISCELLANEOUS:
!  LGTASE denotes the relaxed distance from the end of new strand to GTASE
!  XLEAD,YLEAD,ZLEAD denote the prefered orientation of GTASE



!  NEIGHBOR BONDS FOR TOPOLOGICAL CONSTRAINT ON GTASE:
!  For each complex, count number of neighbor bonds to GTASE as NEINUM. NEIBOND denote the 
!  atoms on those bonds



   CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)


!=============================================================================
!=============================================================================
!=============================================================================


!  DYNAMICS OF GROWTH:

   PRINT*,'-------------------------------------------------------------'
   PRINT*,'DYNAMICS OF GROWTH'


   if(nstart==0) CALL WRITEDCD(JUNIT,NTOTAL,NATOM,NATOMOLD,NATOMnew,NPG,DNOR,ATOR,PEPDIR,X,Y,Z,GLYTIP, &
        newPGID,newPGLEN,NBONDPEP,BONDPEP,BONTYP,NPEPMAX,NSYN,TPPEP,noldbond,oldbond,oldtyp, &
        EDHOLD,EDPEP,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,NFRAME)


   NFRAME=0

   jrestart=0

   pterm=pterm0*ratio(nstart+1)


   DO JSTEP=jstart(nstart+1),5000000000  ! DO INSERTION

!     BEGINNING EACH STEP, CHECK AND UPDATE CONDITIONS

      JFORCE1=0
      JFORCE2=0


      IF(MOD(JSTEP-1,20)==0)THEN
         JFORCE2=1
      END IF

      if(mod(jstep,nprint)==0)then
         JPRINT=1
      else
         jprint=0
      end if

      SKIP=1

      IF(MOD(JSTEP-1,JSKIP)==0)THEN
         SKIP=0
         JFORCE1=1
      END IF


!-------------------------------------------------------------





!     UPDATE RANDOM FORCES:

      IF(JRFORCE+NSYN>NRAND)THEN
         CALL SETRAND(RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRFORCE,NRAND)

      END IF

      IF(MOD(JSTEP-1,1000)==0)THEN

         IF(MOD(JSTEP-1,10000)==0)THEN


            CALL RMPEPBOND(NATOM,NPG,newPGID,newPGLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP, &
                        DNOR,ATOR,NSYN,SYNPG,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOP)

            call updateoldbond(noldbond,oldbond,oldtyp,oldtypmax,joldbond,pgbelowbond, &
                 pgbelowbondmax,nsyn,x,y,z,xedase,yedase,zedase,prmbond,ator)


            CALL LYTGTASE(NATOM,NPG,newPGID,newPGLEN,DNOR,ATOR,PLYT,NATOMDEL,NSYN,SYNPG,SYNLOOP,LOOPLEN,LOOP, &
              noldbond,oldbond,oldtyp)


            NDEL=NATOMDEL-OLDNATOMDEL

            IF(NDEL>0.OR.NBONDDEL>0)THEN

               CALL UPDATESYS(natomold,NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,newPGID,newPGLEN, &
                 NSYN,SYNPG,GLYTIP,GLYSEC,TPPEP,EDPEP,EDHOLD,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                 NBONDGLY,NBONDGLYOLD,NBONDCAP,NBONDPEP,NBONDDEL,BONDGLY,BONDPEP,BONTYP,SYNLOOP,PARTNER, &
                 noldbond,oldbond,jpgnew)

               OLDNATOMDEL=NATOMDEL

               CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            END IF


         END IF

!        UPDATE NEIGHBOR BONDS:


         DO NS=1,NSYN
            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XEDASE,YEDASE,ZEDASE, &
                      NEIGLY,GLYNUM)

         END DO

      END IF



!====================== 

      IF(SKIP==1)THEN
         GOTO 41
      END IF


!     NOW WE CHECK AND UPDATE CONDITIONS OF ENZYMES



      DO NS=1,NSYN

!        UPDATE RANDOM NUMBERS:

         IF(JRANDS>NRAND-100)THEN
            CALL RANDOM_NUMBER(RANDS)
            JRANDS=0
         END IF

!        check if the complex is deactivated:

27       IF(JDEACT(NS)==1)THEN




            call deactivate(ns,nsyn,jstep,nbondpep,nbonddel,jdeact,newpglen,looplen,bontyp,dnor,ator,gtload,looptyp, &
              edcap,sigcleave,jtrap,synloop,synpg,newpgid,loop,bondpep,synthesis,glytip,glysec,crlkage, &
              sigcross,edhold,edpep,tppep)


            CYCLE

         END IF


!-----------------------------------------------------
!        PBP1 TRANSGLYCOSYLATION RULES:


         IF(SYNTHESIS(1,NS)+SYNTHESIS(2,NS)==2)THEN


!           CHECK THE FIRST GTASE:

            IF(GTLOAD(NS)==1)THEN


!              INITIATE NEW STRANDS:

               IF(SYNPG(1,NS)==0.and.SYNPG(2,NS)==0)THEN

                  CALL FIRSTPG(NS,NPG,newPGID,newPGLEN,NATOM,DNOR,ATOR,PEPDIR,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                         SYNPG,GTLOAD,SYNDIR,GLYTIP,GTATRANS)

                  JTRAP(NS)=JSTEP

               ELSE

!                 TO PREVENT ELONGATING FOREVER

                  IPG1=SYNPG(1,NS)

                  IPG2=SYNPG(2,NS)

                  jprevent=0

                  if(ipg1>0)then

                     IF(newPGLEN(IPG1)==LENMAX)then

                        GTLOAD(NS)=0

                        JTERM(1,NS)=1

                        jprevent=1

                     end if

                  end if


                  if(ipg2>0)then

                     if(newPGLEN(IPG2)==LENMAX)THEN

                        GTLOAD(NS)=0

                        JTERM(2,NS)=1

                        jprevent=1

                     end if



                  end if


                  if(jprevent==0)then



!                    UNLOADING PRECURSORS:

                     JRANDS=JRANDS+1


                     IF(P_GLYOUT>RANDS(JRANDS))THEN

                        GTLOAD(NS)=0
                        LGTASE(1:2,NS)=0.5D0


                     ELSE!IF(MOD(JSTEP-1,100)==0)THEN

!                       ELONGATING STRANDS:

                        CALL ELONGATE(NS,NPG,SYNPG,newPGID,newPGLEN,NBONDGLY,BONDGLY,X,Y,Z,XGTASE,YGTASE,ZGTASE,L_G,LGTASE, &
                             GLYTIP,GLYSEC,GTLOAD,NATOM,DNOR,ATOR,PEPDIR,GTATRANS,syndir, &
                             ntemlink,ntemlinkmax,oldglysec,oldglytip,temlink)

                        IF(GTLOAD(NS)==0)THEN

                           CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XEDASE,YEDASE,ZEDASE, &
                                NEIGLY,GLYNUM)

                           JTRAP(NS)=JSTEP


                        END IF

                     END IF

                  END IF

               END IF


            ELSE

!              TRANSLOCATION OF GTASES ON STRANDS:

               IF(GTATRANS(NS)==0)THEN

                  IF(SYNPG(1,NS)==0.and.SYNPG(2,NS)==0)THEN

                     GTATRANS(NS)=1

                  ELSE

!                    TRANSPEPTIDATION-FACILITATED:


                     p_trans=p_transslow

                     jfast1=0

                     ntip1=glytip(1,ns)

                     if(ntip1>0)then

                        if(dnor(ntip1)==0) jfast1=1

                        if(crlkage(1,ns)==0.and.pepdir(ntip1)==syndir(ns)) jfast1=1

                     else

                        jfast1=1

                     end if



                     jfast2=0

                     ntip2=glytip(2,ns)

                     if(ntip2>0)then

                        if(dnor(ntip2)==0.or.ator(ntip2)==0) jfast2=1

                        if(crlkage(2,ns)==0.and.pepdir(ntip2)==-syndir(ns)) jfast2=1

                     else

                        jfast2=1

                     end if


                     if(jfast1==1.and.jfast2==1) p_trans=p_transfast


!                    to prevent long, crosslinked double tails
                     if(crlkage(1,ns)==1.and.crlkage(2,ns)==1)then
                        if(dnor(ntip1)==1.or.(dnor(ntip2)==1.and.ator(ntip2)==1)) p_trans=0.0d0
                     end if

                     JRANDS=JRANDS+1

                     IF(P_TRANS>RANDS(JRANDS))THEN
                        GTATRANS(NS)=1
                     END IF



                  END IF

!              LOADING PRECURSORS:

               ELSE

                  ILOOP=SYNLOOP(NS)

                  AREA=0.01D0*LOOPLEN(ILOOP)**2

                  JRANDS=JRANDS+1

                  IF(P_GLYIN*AREA>RANDS(JRANDS))THEN
                     GTLOAD(NS)=1

!                    RELAXED DISTANCE FROM GTASE TO TIP:

                     IF(SYNPG(1,NS)>0)THEN
                        LGTASE(1,NS)=0.5D0+L_G


                     ELSE
                        LGTASE(1,NS)=0.5D0
                     END IF



                     IF(SYNPG(2,NS)>0)THEN
                        LGTASE(2,NS)=0.5D0+L_G
                     ELSE
                        LGTASE(2,NS)=0.5D0
                     END IF


                  END IF


               END IF


            END IF


!------------------------


!----------------------------------------------------------
!        check if the complex is activated:


         ELSE

            CALL ACTIVATE(NS,SYNTHESIS,SYNLOOP,XGTASE,YGTASE,ZGTASE,X,Y,Z,NEIGLY,GLYNUM, &
                          NLOOP,LOOP,LOOPLEN,LOOPTYP,P_ACTI,jredist)


            if(jredist==1) call REDISTRIBUTE(jsyndir,na_red,NATOM,X,Y,Z,NS,XGTASE,YGTASE,ZGTASE, &
                                XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,SYNDIR)



         END IF




!==================  TRANSPEPTIDATION RULES:


         IF(SYNTHESIS(1,NS)==0.OR.SYNTHESIS(2,NS)==0)THEN
            GOTO 30
         END IF





         IF(TPPEP(2,1,NS)==0)THEN



            IF(TPPEP(1,1,NS)==0)THEN

!              THE FIRST TPASE LOADING DONOR SITE:

               IF(SYNPG(1,NS)/=0.AND.CRLKAGE(1,NS)+CRLKAGE(2,NS)>0)THEN

                  NTIP=glytip(1,ns)!PGID(PGLEN(IPG),IPG)

                  IF(DNOR(NTIP)==1.AND.PEPDIR(NTIP)==SYNDIR(NS))THEN


                     DX=XTPASE(1,NS)-X(NTIP)
                     DY=YTPASE(1,NS)-Y(NTIP)
                     DZ=ZTPASE(1,NS)-Z(NTIP)

                     DIST=SQRT(DX**2+DY**2+DZ**2)

                     IF(DIST<DREACT)THEN

                        P_PEPIN=(1.0-DIST/DREACT)**2*JSKIP

                        JRANDS=JRANDS+1


                        IF(P_PEPIN>RANDS(JRANDS))THEN
                           TPPEP(1,1,NS)=NTIP
                        END IF

                     END IF


                  END IF

               END IF

!           THE FIRST TPASE CROSSLINKING:

            ELSE


               CALL TP1CRLK(NS,npgold,npg,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,newPGID,newPGLEN,BONDPEP,BONTYP, &
                    partner,edpep,PEPDIR,syndir,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,JRANDS,RANDS,DREACT,JSKIP)

!              THE FIRST TPASE WOULD RELEASE PEPTIDE DONOR:

               IF(TPPEP(2,1,NS)==0)THEN

                  JRANDS=JRANDS+1


                  IF(P_PEPOUT>RANDS(JRANDS))THEN
                     TPPEP(1,1,NS)=0

                  ELSE

                     NA=TPPEP(1,1,NS)

                     DX=XTPASE(1,NS)-X(NA)
                     DY=YTPASE(1,NS)-Y(NA)
                     DZ=ZTPASE(1,NS)-Z(NA)

                     IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
                        TPPEP(1,1,NS)=0
                     END IF

                  END IF

               END IF



            END IF


!        SIGNAL FOR THE FIRST TPASE TO CROSSLINK:

         ELSE

            JS=1

            CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,xtpase,ytpase,ztpase,SIGCROSS,DELTA)


!           THE FIRST TPASE WOULD RELEASE PEPTIDE ACCEPTOR:

            IF(SIGCROSS(1,NS)==0)THEN
               JRANDS=JRANDS+1


               IF(P_PEPOUT>RANDS(JRANDS))THEN
                  TPPEP(2,1,NS)=0
               END IF
            END IF



         END IF

!-------------



         IF(TPPEP(2,2,NS)==0)THEN



            IF(TPPEP(1,2,NS)==0)THEN


!              THE SECOND TPASE LOADING DONOR PEPTIDE:

               IF(SYNPG(1,NS)/=0)THEN


                  NTIP=glytip(1,ns)!PGID(PGLEN(IPG),IPG)

                  IF(DNOR(NTIP)==1.AND.PEPDIR(NTIP)==-syndir(ns))THEN


                     DX=XTPASE(2,NS)-X(NTIP)
                     DY=YTPASE(2,NS)-Y(NTIP)
                     DZ=ZTPASE(2,NS)-Z(NTIP)

                     DIST=SQRT(DX**2+DY**2+DZ**2)

                     IF(DIST<DREACT)THEN

                        P_PEPIN=(1.0-DIST/DREACT)**2*JSKIP

                        JRANDS=JRANDS+1


                        IF(P_PEPIN>RANDS(JRANDS))THEN
                           TPPEP(1,2,NS)=NTIP
                        END IF

                     END IF


                  END IF

               END IF


!           THE SECOND TRANSPEPTIDASE CROSSLINKING:

            ELSE


               CALL TP2CRLK(NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,glytip,syndir,PEPDIR, &
                    X,Y,Z,ATOR,JRANDS,RANDS,DREACT,JSKIP)



!              THE SECOND TPASE WOULD RELEASE PEPTIDE DONOR:


               IF(TPPEP(2,2,NS)==0)THEN
                  JRANDS=JRANDS+1


                  IF(P_PEPOUT>RANDS(JRANDS))THEN

                     TPPEP(1,2,NS)=0

                  ELSE

                     NA=TPPEP(1,2,NS)

                     DX=XTPASE(2,NS)-X(NA)
                     DY=YTPASE(2,NS)-Y(NA)
                     DZ=ZTPASE(2,NS)-Z(NA)

                     IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
                        TPPEP(1,2,NS)=0
                     END IF

                  END IF
               END IF



            END IF


!        SIGNAL FOR THE SECOND TPASE TO CROSSLINK:

         ELSE

            JS=2

            CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,xtpase,ytpase,ztpase,SIGCROSS,DELTA)

!           THE SECOND TPASE WOULD RELEASE PEPTIDE ACCEPTOR:

            IF(SIGCROSS(2,NS)==0)THEN
               JRANDS=JRANDS+1


               IF(P_PEPOUT>RANDS(JRANDS))THEN
                  TPPEP(2,2,NS)=0
               END IF
            END IF


         END IF



!-------------



         IF(TPPEP(2,3,NS)==0)THEN


            IF(TPPEP(1,3,NS)==0)THEN


!              THE THIRD TPASE LOADING DONOR SITE:

               IF(SYNPG(2,NS)/=0.AND.CRLKAGE(1,NS)+CRLKAGE(2,NS)>0)THEN

                  NTIP=glytip(2,ns)!PGID(PGLEN(IPG),IPG)

                  IF(DNOR(NTIP)==1.AND.PEPDIR(NTIP)==-SYNDIR(NS))THEN


                     DX=XTPASE(3,NS)-X(NTIP)
                     DY=YTPASE(3,NS)-Y(NTIP)
                     DZ=ZTPASE(3,NS)-Z(NTIP)

                     DIST=SQRT(DX**2+DY**2+DZ**2)

                     IF(DIST<DREACT)THEN

                        P_PEPIN=(1.0-DIST/DREACT)**2*JSKIP

                        JRANDS=JRANDS+1


                        IF(P_PEPIN>RANDS(JRANDS))THEN
                           TPPEP(1,3,NS)=NTIP
                        END IF

                     END IF


                  END IF


               END IF


!           THE THIRD TPASE CROSSLINKING:

            ELSE


               CALL TP3CRLK(NS,npgold,npg,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,newPGID,newPGLEN,BONDPEP,BONTYP, &
                    partner,edpep,PEPDIR,syndir,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,JRANDS,RANDS,DREACT,JSKIP)

!              THE THIRD TPASE WOULD RELEASE PEPTIDE DONOR:

               IF(TPPEP(2,3,NS)==0)THEN

                  JRANDS=JRANDS+1


                  IF(P_PEPOUT>RANDS(JRANDS))THEN
                     TPPEP(1,3,NS)=0

                  ELSE

                     NA=TPPEP(1,3,NS)

                     DX=XTPASE(3,NS)-X(NA)
                     DY=YTPASE(3,NS)-Y(NA)
                     DZ=ZTPASE(3,NS)-Z(NA)

                     IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
                        TPPEP(1,3,NS)=0
                     END IF

                  END IF

               END IF

            END IF


!        SIGNAL FOR THE THIRD TPASE TO CROSSLINK:

         ELSE

            JS=3

            CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,xtpase,ytpase,ztpase,SIGCROSS,DELTA)

!           THE THIRD TPASE WOULD RELEASE PEPTIDE ACCEPTOR:

            IF(SIGCROSS(3,NS)==0)THEN
               JRANDS=JRANDS+1


               IF(P_PEPOUT>RANDS(JRANDS))THEN
                  TPPEP(2,3,NS)=0
               END IF
            END IF



         END IF




!-------------------------------------


!----------------------------

!     POST CROSSLINKING:


         JCHANGE=0


         IF(SIGCROSS(1,NS)==1)THEN

            JS=1

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                       NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,EDPEP,JDEACT,PARTNER, &
                       X,Y,Z,XGTASE,YGTASE,ZGTASE,jrestart)

            if(jrestart==1) exit

            JCHANGE=1

         END IF

         IF(SIGCROSS(2,NS)==1)THEN

            JS=2

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                       NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,EDPEP,JDEACT,PARTNER, &
                       X,Y,Z,XGTASE,YGTASE,ZGTASE,jrestart)


            if(jrestart==1) exit

            JCHANGE=1

         END IF

         IF(SIGCROSS(3,NS)==1)THEN

            JS=3

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                       NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,EDPEP,JDEACT,PARTNER, &
                       X,Y,Z,XGTASE,YGTASE,ZGTASE,jrestart)

            if(jrestart==1) exit

            JCHANGE=1

         END IF

         IF(JCHANGE==1)THEN

            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XEDASE,YEDASE,ZEDASE, &
                      NEIGLY,GLYNUM)

            CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            JTRAP(NS)=JSTEP

         END IF


!============================  PBP4 ENDOPEPTIDATION RULES:

!      CHECK IF EDASE RELEASE PEPTIDES:

30       IF(EDPEP(1,NS)+EDPEP(2,NS)>0)THEN

            IF(J_EDHOLD(NS)+N_EDHOLD<JSTEP)THEN

               IF(EDPEP(1,NS)>0)THEN

                  JRANDS=JRANDS+1


                  IF(P_EDHOLD>RANDS(JRANDS))THEN
                     EDPEP(1,NS)=0
                  END IF

               END IF


               IF(EDPEP(2,NS)>0)THEN

                  JRANDS=JRANDS+1

                  IF(P_EDHOLD>RANDS(JRANDS))THEN
                     EDPEP(2,NS)=0
                  END IF

               END IF

            END IF


         ELSE


            IF(EDHOLD(1,NS)>0)THEN

               NB=EDHOLD(3,NS)

!              trimeric crosslink-triggered mechanism:

               IF(BONTYP(NB)==2)THEN

                     SIGCLEAVE(NS)=1

!              spontaneous mechanism:

               ELSE

                  JRANDS=JRANDS+1


                  IF(PCLEAVE_SPON>RANDS(JRANDS))THEN
                     SIGCLEAVE(NS)=1

!                 Edase might also release captured peptides:

                  ELSE

                     JRANDS=JRANDS+1


                     IF(P_EDHOLD>RANDS(JRANDS))THEN
                        EDHOLD(1:3,NS)=0
                     END IF


                  END IF


               END IF


!            CHECK IF EDASE CAPTURES PEPTIDES:

            ELSE


               IF(EDCAP(NS)>0)THEN

                  JRANDS=JRANDS+1

                  IF(PEDCAP_SPON>RANDS(JRANDS))THEN
                     EDCAP(NS)=0
                  END IF

               END IF


!               allowing cleavage when inactive:

               IF(SYNTHESIS(1,NS)+SYNTHESIS(2,NS)==0)THEN

                  JRANDS=JRANDS+1

                  IF(PEDCAP_INACT>RANDS(JRANDS))THEN
                     EDCAP(NS)=2
                  END IF

               ELSE

!                 trimeric crosslink-triggered mechanism:

                  IF(TPPEP(1,1,NS)>0.AND.TPPEP(2,1,NS)==0.AND. &
                     TPPEP(1,3,NS)>0.AND.TPPEP(2,3,NS)==0)THEN

                     EDCAP(NS)=1

                  ELSE

!                 spontaneous mechanism:

                     JRANDS=JRANDS+1

                     IF(PEDCAP_SPON>RANDS(JRANDS))THEN
                        EDCAP(NS)=3
                     END IF

                  END IF

               END IF

               IF(EDCAP(NS)>0)THEN

                  CALL ENDOPEP(NS,NSYN,SYNLOOP,LOOP,LOOPLEN,NBONDPEP,BONDPEP,BONTYP,SYNPG,newPGID,newPGLEN,PARTNER,EDHOLD, &
                  pepdir,syndir,EDCAP,X,Y,Z,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD,JRANDS,RANDS,JSKIP,glytip,l_edcap)

               END IF

            END IF

         END IF



!------------------------------------------------------------------------



!----------------------------------

!      POST CLEAVING:



         IF(SIGCLEAVE(NS)==1)THEN

            CALL POSTCLEAVE(NS,NSYN,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOP,NLOOPDEL,EDHOLD,SIGCLEAVE, &
                 JDEACT,X,Y,Z,NBONDDEL,NBONDPEP,BONDPEP,BONTYP,DNOR,ATOR,EDPEP,EDCAP,JSTEP,J_EDHOLD,PARTNER, &
                 noldbond,oldbond,oldtyp)


            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XEDASE,YEDASE,ZEDASE, &
                      NEIGLY,GLYNUM)

            CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            JTRAP(NS)=JSTEP

         END IF



!======================= TERMINATION


!        FIRST STRAND TERMINATION

         IF(JTERM(1,NS)==1)THEN

               SYNPG(1,NS)=0


               GTLOAD(NS)=0

               GLYTIP(1,NS)=0

               GLYSEC(1,NS)=0

               TPPEP(1:2,1,NS)=0

               TPPEP(1:2,2,NS)=0

               SIGCROSS(1:2,NS)=0

               CRLKAGE(1,NS)=0

               JTERM(1,NS)=0


         ELSEIF(CRLKAGE(1,NS)>2.AND.CRLKAGE(2,NS)>2.AND.TPPEP(1,1,NS)==0.AND.TPPEP(1,2,NS)==0)THEN

            IPG=SYNPG(1,NS)

            length=newpglen(ipg)

            ILOOP=SYNLOOP(NS)

            AREA=0.01D0*LOOPLEN(ILOOP)**2

            JRANDS=JRANDS+1

            IF(PTERM*CRLKAGE(1,NS)*lenmax/(lenmax-length)>RANDS(JRANDS)*AREA)THEN

               SYNPG(1,NS)=0


               GTLOAD(NS)=0

               GLYTIP(1,NS)=0

               GLYSEC(1,NS)=0

               SIGCROSS(1:2,NS)=0

               CRLKAGE(1,NS)=0



            END IF

         END IF

!        SECOND STRAND TERMINATION

         IF(JTERM(2,NS)==1)THEN

               SYNPG(2,NS)=0


               GTLOAD(NS)=0

               GLYTIP(2,NS)=0

               GLYSEC(2,NS)=0

               TPPEP(1:2,3,NS)=0

               TPPEP(1:2,2,NS)=0

               SIGCROSS(2:3,NS)=0

               CRLKAGE(2,NS)=0

               JTERM(2,NS)=0


         ELSEIF(CRLKAGE(2,NS)>2.AND.CRLKAGE(1,NS)>2.AND.TPPEP(1,3,NS)==0.AND.TPPEP(1,2,NS)==0)THEN

            IPG=SYNPG(2,NS)

            length=newpglen(ipg)

            ILOOP=SYNLOOP(NS)

            AREA=0.01D0*LOOPLEN(ILOOP)**2

            JRANDS=JRANDS+1


            IF(PTERM*CRLKAGE(2,NS)*lenmax/(lenmax-length)>RANDS(JRANDS)*AREA)THEN

               SYNPG(2,NS)=0


               GTLOAD(NS)=0

               GLYTIP(2,NS)=0

               GLYSEC(2,NS)=0

               SIGCROSS(2:3,NS)=0

               CRLKAGE(2,NS)=0



            END IF

         END IF

!======================= DEACTIVATION SIGNAL


!        signal for deactivation due to trapping:

         IF(JSTEP-JTRAP(NS)>NTRAP.AND.SYNPG(1,NS)>0.AND.SYNPG(2,NS)>0)THEN
            JDEACT(NS)=1
         END IF

!        signal for spontaneous deactivation:

         IF(JDEACT(NS)==0.AND.MOD(JSTEP-1,10000)==0.AND.SYNTHESIS(1,NS)==1.AND.SYNTHESIS(2,NS)==1)THEN

            JRANDS=JRANDS+1

            IF(P_DEACT>RANDS(JRANDS).or.abs(xedase(ns))>lambda)THEN
               JDEACT(NS)=1

               if(abs(xedase(ns))>lambda) then

                  call REDISTRIBUTE(jsyndir,na_red,NATOM,X,Y,Z,NS,XGTASE,YGTASE,ZGTASE, &
                       XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,SYNDIR)

                  call deactivate(ns,nsyn,jstep,nbondpep,nbonddel,jdeact,newpglen,looplen,bontyp,dnor,ator,gtload,looptyp, &
                       edcap,sigcleave,jtrap,synloop,synpg,newpgid,loop,bondpep,synthesis,glytip,glysec,crlkage, &
                       sigcross,edhold,edpep,tppep)

                  cycle


               end if

            END IF
         END IF

!        deactivation if two SYNCOMPs meet:

         if(ns<nsyn)then

            ns_ter=0

            do ns1=ns+1,nsyn

               if(synloop(ns)==synloop(ns1).and.synloop(ns)>0)then

                  jrands=jrands+1

                  if(RANDS(JRANDS)>0.5)then

                     ns_ter=ns

                  else

                     ns_ter=ns1

                  end if

                  exit

               end if

            end do

            if(ns_ter>0)then

               jdeact(ns_ter)=1

               call REDISTRIBUTE(jsyndir,na_red,NATOM,X,Y,Z,NS_ter,XGTASE,YGTASE,ZGTASE, &
                       XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,SYNDIR)

               call deactivate(ns_ter,nsyn,jstep,nbondpep,nbonddel,jdeact,newpglen,looplen,bontyp,dnor,ator,gtload,looptyp, &
                    edcap,sigcleave,jtrap,synloop,synpg,newpgid,loop,bondpep,synthesis,glytip,glysec,crlkage, &
                    sigcross,edhold,edpep,tppep)


            end if

         end if

      END DO





!================================================

!     this is to cope with large loop, and too many short strands:

      IF(MOD(JSTEP-1,10000)==0)THEN

         IF(MAXVAL(LOOPLEN(1:NLOOP))*2>LOOPLENMAX)THEN


            LOOPLENMAX=MAXVAL(LOOPLEN(1:NLOOP))*2

            ALLOCATE(NEWLOOP(MAXVAL(LOOPLEN(1:NLOOP)),NLOOP))

            DO NL=1,NLOOP
               NEWLOOP(1:LOOPLEN(NL),NL)=LOOP(1:LOOPLEN(NL),NL)
            END DO

            DEALLOCATE(LOOP)

            ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX))

            DO NL=1,NLOOP
               LOOP(1:LOOPLEN(NL),NL)=NEWLOOP(1:LOOPLEN(NL),NL)
            END DO

            DEALLOCATE(NEWLOOP)

            print*,'loop size expanded',LOOPLENMAX

         END IF


         IF(NPG>NPGnew-10)THEN

            print*,'reallocation of newPG',npg,npgnew

            ALLOCATE(NEWPGIDtemp(LENMAX,NPG),NEWPGLENtemp(NPG))
            NEWPGLENtemp(1:NPG)=newPGLEN(1:NPG)
            NEWPGIDtemp(1:LENMAX,1:NPG)=newPGID(1:LENMAX,1:NPG)

            DEALLOCATE(newPGID,newPGLEN)

            NPGnew=NPGnew+100

            ALLOCATE(newPGID(LENMAX,NPGnew),newPGLEN(NPGnew))

            newPGLEN(1:NPG)=NEWPGLENtemp(1:NPG)
            newPGID=0
            newPGID(1:LENMAX,1:NPG)=NEWPGIDtemp(1:LENMAX,1:NPG)

            DEALLOCATE(NEWPGIDtemp,NEWPGLENtemp)

         END IF

      END IF

!======================================================================



!     FORCE CALCULATION:

41    IF(JFORCE2==1)THEN

         CALL FPRES(NATOM,FXPRES,FYPRES,FZPRES,X,Y,Z,NLOOP,LOOP,LOOPLEN,LOOPTYP,PRES,DELTA,INVDELTA,caplen, &
          npress,pconstrict,invdelx,fycon,fzcon,nsyn)

      END IF



      IF(JFORCE1==1)THEN

         FXREP(1:NATOM)=FXPRES(1:NATOM)
         FYREP(1:NATOM)=FYPRES(1:NATOM)
         FZREP(1:NATOM)=FZPRES(1:NATOM)


         CALL FPEPBONDS(FXREP,FYREP,FZREP,X,Y,Z,NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND, &
                        L_P,K_P,LSTART,LREP,LSWITCH,KSWITCH,noldbond,oldbond,oldtyp)



         CALL FANGLES(NPGold,PGID,PGLEN,KTHETA,THETA_0,BETA,DELTA,INVDELTA,X,Y,Z,FXREP,FYREP,FZREP)


         if(npg>0) call FANGLES(NPG,newPGID,newPGLEN,KTHETA,THETA_0,BETA,DELTA,INVDELTA,X,Y,Z,FXREP,FYREP,FZREP)


         CALL SURFCONSTR(natom,NPG,newPGLEN,DNOR,ATOR,newPGID,KSUR,X,Y,Z,FYREP,FZREP,pthick)


         if(noldbond<128)then
            call oldbondpush(noldbond,npg,oldbond,oldtyp,pgbelowbond,newpgid,newpglen,wthick,kwall,x,y,z,fyrep,fzrep)

         else

            call oldbondpush1(nthreads,natomold,natom,noldbond,npg,oldbond,oldtyp, &
                 pgbelowbond,newpgid,newpglen,wthick,kwall,x,y,z,fyrep,fzrep)

         end if

      END IF



      FX(1:NATOM)=FXREP(1:NATOM)
      FY(1:NATOM)=FYREP(1:NATOM)
      FZ(1:NATOM)=FZREP(1:NATOM)


      CALL FGLYBONDS(NPGold,PGID,PGLEN,L_G,K_G,FX,FY,FZ,X,Y,Z)

      if(npg>0) call FGLYBONDS(NPG,newPGID,newPGLEN,L_G,K_G,FX,FY,FZ,X,Y,Z)



!     -----------------------------

      CALL RANDFORCES(NSYN,FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE, &
                      RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRFORCE, &
                      YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE)

!     -----------------------------

      CALL SYNCOMP(NSYN,SYNDIR,edpep,tppep,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                      KPAIR,LPAIR,KGTTP,LGTTP13,LGTTP2,KGTED,LGTED,KSIDE,LSIDE,DELTA,BETA,ktped,ltped, &
                      FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE)


!     -----------------------------


      CALL GTAHOLD(NSYN,GLYTIP,GLYSEC,KGTASE,LGTASE,KLEAD,KTHETA,THETA_0,DELTA,INVDELTA,BETA,X,Y,Z, &
                      FX,FY,FZ,XGTASE,YGTASE,ZGTASE,FXGTASE,FYGTASE,FZGTASE,syndir)




!     -----------------------------

      CALL TPAHOLD(NSYN,X,Y,Z,FX,FY,FZ,TPPEP,XTPASE,YTPASE,ZTPASE,FXTPASE,FYTPASE,FZTPASE,KTPASE,LTPASE)


!     -----------------------------

      CALL EDAHOLD(NSYN,EDHOLD,EDPEP,KEDASE,LEDASE,X,Y,Z,FX,FY,FZ, &
                   XEDASE,YEDASE,ZEDASE,FXEDASE,FYEDASE,FZEDASE)


!     -----------------------------------------


!     -----------------------------------------
!     CONSTRAIN PBPS AND UNLINKED PG TO THE SURFACE:


      CALL SURFLINK(NSYN,loop,looplen,synloop,y,z,YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE, &
              FYGTASE,FZGTASE,FYTPASE,FZTPASE,FYEDASE,FZEDASE,KSUR,pthick)


!     constrain the cell's axis to x-axis

      call axis(caplen,kaxis,y,z,fy,fz)



      call temlinkforce(ntemlink,temlink,KTHETA,THETA_0,BETA,INVDELTA,DELTA,l_g,k_g,X,Y,Z,FX,FY,FZ)





!======================================================================

!      IF(JFORCE1==1)THEN

         FAMAG(1:NATOM)=ABS(FX(1:NATOM))+ABS(FY(1:NATOM))+ABS(FZ(1:NATOM))

         FAMAX=MAXVAL(FAMAG(1:NATOM))

!      END IF

         FMAGPBP(1:NSYN)=ABS(FXGTASE(1,1:NSYN))+ABS(FYGTASE(1,1:NSYN))+ABS(FZGTASE(1,1:NSYN))

         FPMAX1=MAXVAL(FMAGPBP(1:NSYN))


         FMAGPBP(1:NSYN)=ABS(FXGTASE(2,1:NSYN))+ABS(FYGTASE(2,1:NSYN))+ABS(FZGTASE(2,1:NSYN))

         FPMAX2=MAXVAL(FMAGPBP(1:NSYN))



         FPBPMAX=MAX(FPMAX1,FPMAX2)

         GAM=MIN(0.01D0/FAMAX,0.1D0/FPBPMAX)

!      END IF

!---------------------------------
!-------------------------------------------------------

      IF(SKIP==0)THEN
         XEDASEOLD(1:NSYN)=XEDASE(1:NSYN)
         YEDASEOLD(1:NSYN)=YEDASE(1:NSYN)
         ZEDASEOLD(1:NSYN)=ZEDASE(1:NSYN)
      END IF

      CALL NEWCOOR(NATOM,X,Y,Z,FX,FY,FZ,NSYN,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                   FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE, &
                   FXEDASE,FYEDASE,FZEDASE,INVMPBP,GAM)

!======================================================================

!      IF(MOD(JSTEP,100000)==0)THEN


      if(jprint==1)then

         fconstrict=sum(sqrt(fycon(1:natom)**2+fzcon(1:natom)**2))

         WRITE(*,2)'STEP',JSTEP,'NATOM',NATOM,'MAX (nN)',0.01*max(FAMAX,FPBPMAX),'constrict (nN)',0.01*fconstrict

         CALL WRITEDCD(JUNIT,NTOTAL,NATOM,NATOMOLD,NATOMnew,NPG,DNOR,ATOR,PEPDIR,X,Y,Z,GLYTIP, &
              newPGID,newPGLEN,NBONDPEP,BONDPEP,BONTYP,NPEPMAX,NSYN,TPPEP,noldbond,oldbond,oldtyp, &
              EDHOLD,EDPEP,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,NFRAME)




         IF(NFRAME>=1000)THEN


            CALL PG_OUT(NSTART,NATOMOLD,NATOM,NATOMDEL,OLDNATOMDEL,DNOR,ATOR,PEPDIR,X,Y,Z, &
              NPG,NPGOLD,PGID,PGLEN,PGLENMAX,lenmax,newPGid,newpglen,NBONDPEP,NBONDDEL,noldbond,oldbond,oldtyp, &
               BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP,caplen, &
                jsyndir,na_red,NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDCAP,EDHOLD,CRLKAGE, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD, &
                  oldglytip,oldglysec,temlink,ntemlink)

            call system_clock(timerun)

            jstart(nstart+1)=jstep+1


!           update average length of new PG


               length=sum(newpglen(1:npg))

               newratio=1.0*length/npg/newlengthave

               newratio=newratio**4!newratio


               ratio(nstart+1)=ratio(nstart)*newratio

               pterm=pterm0*ratio(nstart+1)

               print*,'new PTERM, ratio',pterm,ratio(nstart+1)


            OPEN(10,FILE='restart.dat')
            write(10,*)'      nstart       jfile  jstart(nstart+1)   lasttime         ratio(nstart+1)'
            WRITE(10,*)NSTART,JFILE,jstart(nstart+1),(TIMERUN-TIMESTART)/rate+LASTTIME,ratio(nstart+1)
            CLOSE(10)


            CLOSE(JUNIT)
            JFILE=JFILE+1
            CALL DCDHEADER(JUNIT,JFILE,NTOTAL)
            NFRAME=0


         END IF


      END IF


!---------------------------------


!---------------------------------

      IF(NATOM>=NATOMMAX-100)THEN



            call system_clock(timerun)



            EXIT


      END IF

   END DO  ! END DO INSERTION


   if(jrestart==1)then

      print*,'restart job from configuration #',nstart

      CLOSE(JUNIT)


      OPEN(10,FILE='restart.dat')
      read(10,*)
      READ(10,*)NSTART,JFILE,jstart(nstart+1),LASTTIME,ratio(nstart+1)
      CLOSE(10)


      goto 35

   end if


!================================================

   call system_clock(timerun)

   TIMERUN=(TIMERUN-TIMESTART)/rate+LASTTIME

   days=timerun/86400
   timerun=timerun-86400*days

   hours=timerun/3600
   timerun=timerun-3600*hours

   mins=timerun/60
   secs=timerun-60*mins

   write(*,1)'RUNNING TIME =',DAYS,'days : ',HOURS,'hours : ',MINS,'mins : ',secs,'secs'

1  FORMAT(A14,2X,I2,1X,A7,I2,1X,A8,1X,I2,1X,A7,1x,I2,1x,A4)



2  FORMAT(A4,2X,I9,2X,A5,2X,I6,2X,A9,2X,F6.3,2x,a14,2x,F8.1)



END
