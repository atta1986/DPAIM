   program D20
   implicit real*8 (a-h,o-z)
!   implicit none
   Integer :: nsa,nsb,typ,nsmax,vers,maxat,max_elem,maxc,nn
   Real*8 :: aswitch,bswitch,conv,a0,zero,autoang,autokcal,c6conv,autoev
   parameter (nsmax=1000,aswitch=1.0,bswitch=8.8,conv=0.529177209,vers=3)

   parameter (maxat=20000)
   parameter (max_elem=94)
! maximum coordination number references per element
   parameter (maxc=5)
! coversion factors
   parameter (autoang=0.52917726d0)
   parameter (autokcal=627.509541d0)
   parameter (autoev= 27.21138505)

   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr,disper,rad2d,pi,AA0
!   Integer :: nsa,nsb,typ
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   Real*8,dimension(:,:) :: sita(10,nsmax),sitb(10,nsmax),sitaa(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax),tran(3,3)
   Real, dimension(:) :: aol(3), bol(3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(5,nsmax) :: latom(5,nsmax),atom(5,nsmax)
   Integer,dimension(:):: chara(nsmax),charb(nsmax)
   dimension:: bweight(nsmax,nsmax),sit_bnd(3)
   data a0 / 0.529177249d0/, zero /0.d0/
   DATA doit/'DOIT'/,done/'DONE'/,doitr6/'DOR6'/
   data printed /.false./
   Logical :: switch_on,onesystem
   Real*8, Dimension(98,98) :: AA,BB
!   pi = dacos(-1.d0)
!   rad2d = 180.d0/pi
!      implicit none             
!   integer :: maxat,max_elem,maxc                      
! conversion factors
!   real*8 :: autoang,autokcal,c6conv,autoev
!   parameter (maxat=20000)
!   parameter (max_elem=94)
! maximum coordination number references per element
!   parameter (maxc=5)
! coversion factors
!   parameter (autoang=0.52917726d0)
!   parameter (autokcal=627.509541d0)
!   parameter (autoev= 27.21138505)

! DFT-D version
   integer version
! number of atoms
!   integer n
! coordinates in au
   real*8,dimension(:,:), allocatable :: xyz,abc
! fixed atoms in geometry opt
   logical fix(maxat)
! lattice in au
   real*8 lat(3,3)
! gradient
   real*8,dimension(:,:), allocatable :: g      
   real*8 g_lat(3,3)
! cardinal numbers of elements
   integer,dimension(:), allocatable :: iz  
! cut-off radii for all element pairs
   real*8 r0ab(max_elem,max_elem)
! C6 for all element pairs 
   real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
! how many different C6 for one element
   integer mxc(max_elem)
! C6810 
   real*8 c6,c8,c10
! coordination numbers of the atoms
   real*8,dimension(:), allocatable :: cn  
! covalent radii
   real*8 rcov(max_elem)
! atomic <r^2>/<r^4> values
   real*8 r2r4(max_elem)
! energies
   real*8 e6, e8, e10, e12, disp, e6abc        
! THE PARAMETERS OF THE METHOD (not all a "free")
   real*8 rs6, rs8, rs10, s6, s18, alp6, alp8, alp10, s42, rs18, alp
! printout option
   logical echo
! grad ?
   logical grad
! analyse results ?
   logical anal
! third-order term?
   logical noabc
! gradient calctype
   logical numgrad
! special parameters
   logical tz
! periodic boundary conditions
   logical pbc
! repetitions of the unitcell to match the rthr and c_thr
   integer rep_vdw(3),rep_cn(3)
! R^2 distance neglect threshold (important for speed in case of large systems)
   real*8 rthr,rthr2
! R^2 distance to cutoff for CN_calculation
   real*8 cn_thr
! Integer for assigning only max/min cn C6 (0=normal, 1=min, 2=max)
! local and dummy variables
   character*80 atmp,btmp,ctmp,dtmp,etmp,ftmp,func
   character*2  esym 
   integer i,j,k,z,iat,jat,i1,i2,cnn
   integer ida(max_elem),ipot
   real*8  x,y,dispr,displ,gdsp,dum,xx(10),dum6(86)
   real*8  dum1,dum2,dum3(3),Veff3_a(nsmax),Veff3_b(nsmax),Veff5_a(nsmax),Veff5_b(nsmax)
   logical ex,pot,fdum
   logical minc6list(max_elem),maxc6list(max_elem),minc6,maxc6
! these new data are scaled with k2=4./3.  and converted a_0 via
! autoang=0.52917726d0
! C 1.88972601
   data rcov/&
   & 0.84628308, 1.15903197, 3.02356173, 2.36845659, 2.00011865,&
   & 1.95972601, 1.78894056, 1.58736983, 1.71256616, 1.68815527,&
   & 3.52748848, 3.14954334, 2.94718717, 2.72041997, 2.77159820,&
   & 2.57002732, 2.55443835, 2.41884923, 4.43455700, 3.88023730,&
   & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,&
   & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,&
   & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,&
   & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,&
   & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,&
   & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,&
   & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,&
   & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,&
   & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,&
   & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,&
   & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,&
   & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,&
   & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,&
   & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,&
   & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

! k1-k3
!   include 'param'
   echo=.true. 
   grad=.false.
   pot =.false.
   anal=.false.
   noabc=.true. 
   numgrad=.false.
   tz=.false.
   func=' none (read from parameter file)'
   version=0
   pbc=.false.
   fix=.false.
   minc6=.false.
   maxc6=.false.
   minc6list=.false.
   maxc6list=.false.
   fdum=.false.
! Cutoff r^2 thresholds for the gradient in bohr^2.
! rthr influences N^2 part of the gradient.
! rthr2 influences the N^3 part of the gradient. When using
! dftd3 in combination with semi-empirical methods or FFs, and large
! (>1000 atoms) systems, rthr2 is crucial for speed:
! Recommended values are 20^2 to 25^2 bohr.

   rthr=9000.0d0   ! UR, SE
   rthr2=1600.0d0
   cn_thr=1600.0d0

! J/mol nm^6 - > au
   c6conv=1.d-3/2625.4999d0/((0.052917726d0)**6)
   pi = dacos(-1.d0)
   rad2d = 180.d0/pi

! get coord filename
   call getarg(1,etmp)
   inquire(file=etmp,exist=ex)
!   if(.not.ex) call printoptions       
   ex=.false.
   ipot=0


   switch_on=.true.
!   onesystem=.true.
   onesystem=.false.
!   Read in the sites of monomer A
!   WRITE(fileA,'(a)') "fileA.dat"
   !open(unit=5,file='DPAIM_Dispersion_enegy.dat',form='formatted',status='old')

   open(unit=5,file='DPAIM_Dispersion_enegy.dat')
   if(onesystem) then
     open(unit=1,file='file.dat')
     read(1,*)nsa
     do i=1,nsa
         read(1,*)(sita(j,i),j=1,4)
         sita(1,i)=sita(1,i)/a0
         sita(2,i)=sita(2,i)/a0
         sita(3,i)=sita(3,i)/a0
     end do
     nsb=nsa
     sitaa=sita
     sitbb=sita
   else
!     print*,'else=',nsa
     open(unit=1,file='file.dat')
     read(1,*)nsa
!   read(1,*)
!   open(99,file='file.dat')
!   write(99,*)nsa

     do i=1,nsa
         read(1,*)(sita(j,i),j=1,6)
         sita(1,i)=sita(1,i) !/a0
         sita(2,i)=sita(2,i) !/a0
         sita(3,i)=sita(3,i) !/a0
         Veff3_a(i)=sita(5,i)
         Veff5_a(i)=sita(6,i)
!         print*,(sita(j,i),j=1,4)
!         print*,'****************'
!         print*,sita(:,1)
     end do

     read(1,*)nsb
   
     do i=1,nsb
         read(1,*)(sitb(j,i),j=1,6)
         sitb(1,i)=sitb(1,i) !/a0
         sitb(2,i)=sitb(2,i) !/a0
         sitb(3,i)=sitb(3,i) !/a0
         Veff3_b(i)=sitb(5,i)
         Veff5_b(i)=sitb(6,i)
!         print*,(sitb(j,i),j=1,4)
!         print*,ltempb(i)
     end do  
     close(1)
     sitaa=sita
     sitbb=sitb
   end if
   call analyze_hydrogens(sita,nsa,sitaa)
   sita=sitaa
!   do i=1,nsa
!      print*,'sita=',sita(:,i)
!   end do
!   print*,'connections in Monomer_A'
   call analyze_Boron(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_carbons2(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Nitrogen(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Oxygen(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Flourine(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Aluminium(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Silicon(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Phosphorus(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Sulfur(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Chlorine(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Bromine(sita,nsa,rcov,sitaa)
   sita=sitaa
   call analyze_Iodine(sita,nsa,rcov,sitaa)
   sita=sitaa
!   do i=1,nsa
!      print*,'sita=',sita(:,i)
!   end do
!   print*,'*****************'
!   do i=1,nsa      
!      print*,sita(5,i)
!   end do
   call analyze_hydrogens(sitb,nsb,sitbb)
   sitb=sitbb
!   do i=1,nsb
!      print*,sitb(:,i)
!   end do
!   print*,'connections in Monomer_B'
   call analyze_Boron(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_carbons2(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Nitrogen(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Oxygen(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Flourine(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Aluminium(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Silicon(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Phosphorus(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Sulfur(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Chlorine(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Bromine(sitb,nsb,rcov,sitbb)
   sitb=sitbb
   call analyze_Iodine(sitb,nsb,rcov,sitbb)
   sitb=sitbb
!   do i=1,nsb
!      print*,sitb(:,i)
!   end do
!   print*,sita(:,1)
!   disper=-10.5

!   call setr0ab(max_elem,autoang,r0ab)

   Call fitted_dispersion(sita,sitb,r0ab,nsa,nsb,Veff3_a,Veff3_b,Veff5_a,Veff5_b,disper)
   
!   BB=0.0
!   print*,'disper',disper
   write(5,'(A)')'**********'
   write(5,'(A)')' '
   write(5,'(A,F24.14)') 'DPAIM dispersion energy (kcal/mol)= ', 627.509474*disper
   write(5,'(A)')' '
   write(5,'(A)')'**********'
   close(5)
   open(27,file='aaa')
   write(27,'(F24.14)') 627.509474*disper
   close(27)
!
!   BB=AA
!   Do i=1,28
!     Do j=i+1,28
!       BB(i,j)=AA(i,j)+AA(j,i)
!     end do
!   end do

!   open(28,file='contribution.dat')
!   do i=1,28
!      do j=1,28
!         if (i.le.j) write(28,'(f18.12,1x)',advance='no') BB(i,j)
!      end do
!   end do
!   write(28,*) ''
!   close(28)

!   open(28,file='contribution.dat')
!   Do i=1,28
!     Do j=1,28
!       if (i.eq.j.or.i.lt.j) then
!          write(28,*) BB(i,j)
!       else
!         continue
!       end if
!     end do
!   end do
!   close(28)   


!   print*,'D20 Dispersion: ',disper
!
! For Periodic systems (from Grimme D3 dispersion function)
!
   if (pbc) then
     call pbcrdatomnumber(etmp,nn)
!   else
!     call rdatomnumber(etmp,n)
   end if
   allocate(xyz(3,nn))
   allocate(g(3,nn))
   allocate(iz(nn))
   allocate(cn(nn))
! reading coordinates and cell in VASP.5.2-format
! determing repetitions of unitcell
   if (pbc) then
           call pbcrdcoord(etmp,lat,nn,xyz,iz,autoang,sita)
           call set_criteria(rthr,lat,dum3)
           rep_vdw=int(dum3)+1
           call set_criteria(cn_thr,lat,dum3)
           rep_cn=int(dum3)+1
!   end if
     if(nn.lt.1)     call stoprun( 'no atoms' )
     if(nn.gt.maxat) call stoprun( 'too many atoms' )

!  Call PBC subroutine:
     call  pbcedisp(max_elem,maxc,nn,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
   &     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
   &     e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)

     e6   = e6   *s6

   ! e6abc= e6abc*s6 ! old and wrong !

     e8   = e8   *s18

     disp =-e6-e8-e6abc

   end if
   end program D20

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C compute energy
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
   subroutine pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
  &           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
  &           e6,e8,e10,e12,e63,lat,rthr,rep_vdw,cn_thr,rep_cn)
   implicit none  
   integer max_elem,maxc
   real*8 r2r4(max_elem),rcov(max_elem)
   real*8 rs6,rs8,rs10,alp6,alp8,alp10
   real*8 rthr,cn_thr,crit_cn
   integer rep_vdw(3),rep_cn(3)
   integer n,iz(*),version,mxc(max_elem)
!   integer rep_v(3)=rep_vdw!,rep_cn(3)
   real*8 xyz(3,*),r0ab(max_elem,max_elem),lat(3,3)!,r2r4(*)
!   real*8 rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
   real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
   real*8 e6, e8, e10, e12, e63!,crit_vdw,crit_cn
   logical noabc
 
   integer iat,jat,kat
   real*8 r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,c10,ang,rav,R0
   real*8 damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14
   real*8 cn(n),rxyz(3),dxyz(3)
   real*8 r2ab(n*n),cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,tau(3)
   integer ij,ik,jk !lin
   integer taux,tauy,tauz,counter
   real*8 a1,a2  !BJ-parameter
   real*8 bj_dmp6,bj_dmp8
   real*8 tmp1,tmp2


   e6 =0
   e8 =0
   e10=0
   e12=0
   e63=0
   tau=(/0.0,0.0,0.0/)
   counter=0
   crit_cn=cn_thr
!
   Do iat=1,n
      do jat=iat+1,n
! get C6
!         call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
!        &                              cn(iat),cn(jat),c6)

         rxyz=xyz(:,iat)-xyz(:,jat)
!         r42=r2r4(iz(iat))*r2r4(iz(jat))
!         bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
!         bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

!         if(.not.noabc)then
!           ij=lin(jat,iat)
! store C6 for C9, calc as sqrt
!           cc6ab(ij)=sqrt(c6)
!         endif
         do taux=-rep_vdw(1),rep_vdw(1)
         do tauy=-rep_vdw(2),rep_vdw(2)
         do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            
            dxyz=rxyz+tau

            r2=sum(dxyz*dxyz)
! cutoff
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r


            r6=r2**3      

            e6 =e6+c6/(r6+bj_dmp6)

! stored in main as sqrt
            c8 =3.0d0*c6*r42
            r8 =r6*r2

            e8 =e8+c8/(r8+bj_dmp8)

            counter=counter+1

         enddo !tauz
         enddo !tauy
         enddo !taux
      enddo !jat

! Now the self interaction
      jat=iat
! get C6
!      call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
!     &                                  cn(iat),cn(jat),c6)
!      r42=r2r4(iz(iat))*r2r4(iz(iat))
!      bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
!      bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8
           
!      if(.not.noabc)then
!         ij=lin(jat,iat)
! store C6 for C9, calc as sqrt
!         cc6ab(ij)=dsqrt(c6)
!      endif

      do taux=-rep_vdw(1),rep_vdw(1)
      do tauy=-rep_vdw(2),rep_vdw(2)
      do tauz=-rep_vdw(3),rep_vdw(3)
         if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            r2=sum(tau*tau)
! cutoff
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r


            r6=r2**3      

            e6 =e6+c6/(r6+bj_dmp6)*0.50d0

! stored in main as sqrt
            c8 =3.0d0*c6*r42
            r8 =r6*r2

            e8 =e8+c8/(r8+bj_dmp8)*0.50d0
            counter=counter+1

      enddo !tauz
      enddo !tauy
      enddo !taux
   enddo !iat
   end subroutine pbcedisp
!
   subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

      character(len=*),intent(inout) :: str,before
      character(len=*),intent(in) :: delims
      character,optional :: sep
      logical :: pres
      character :: ch,cha

      pres=present(sep)
      str=adjustl(str)
      call compact(str)
      lenstr=len_trim(str)
      if(lenstr == 0) return        ! string str is empty
      k=0
      ibsl=0                        ! backslash initially inactive
      before=' '
      do i=1,lenstr
         ch=str(i:i)
         if(ibsl == 1) then          ! backslash active
            k=k+1
            before(k:k)=ch
            ibsl=0
            cycle
         end if
         if(ch == '\') then          ! backslash with backslash inactive
            k=k+1
            before(k:k)=ch
            ibsl=1
            cycle
         end if
         ipos=index(delims,ch)         
         if(ipos == 0) then          ! character is not a delimiter
            k=k+1
            before(k:k)=ch
            cycle
         end if
         if(ch /= ' ') then          ! character is a delimiter that is not aspace
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
         cha=str(i+1:i+1)            ! character is a space delimiter
         iposa=index(delims,cha)
         if(iposa > 0) then          ! next character is a delimiter
            str=str(i+2:)
            if(pres) sep=cha
            exit
         else
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
      end do
      if(i >= lenstr) str=''
      str=adjustl(str)              ! remove initial spaces
      return

   end subroutine split
!
   subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control
! characters;
! removes initial spaces.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str)):: outstr
      
      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      isp=0
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
  
        select case(ich)
  
          case(9,32)     ! space or tab character
            if(isp==0) then
              k=k+1
              outstr(k:k)=' '
            end if
            isp=1
            
          case(33:)      ! not a space, quote, or control character
            k=k+1
            outstr(k:k)=ch
            isp=0
      
        end select
        
      end do

      str=adjustl(outstr)

   end subroutine compact

!**********************************************************************

   subroutine removesp(str)

      ! Removes spaces, tabs, and control characters in string str

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
        select case(ich)    
          case(0:32)  ! space, tab, or control character
               cycle       
          case(33:)  
            k=k+1
            outstr(k:k)=ch
        end select
      end do
      
      str=adjustl(outstr)
      
   end subroutine removesp
!
   subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0
      ibsl=0                        ! backslash initially inactive
      
      do i=1,lenstr
        ch=str(i:i)
        if(ibsl == 1) then          ! backslash active
         k=k+1
         outstr(k:k)=ch
         ibsl=0
         cycle
        end if
        if(ch == '\') then          ! backslash with backslash inactive
         ibsl=1
         cycle
        end if
        k=k+1
        outstr(k:k)=ch              ! non-backslash with backslash inactive
      end do
      
      str=adjustl(outstr)
      
   end subroutine removebksl



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       
!c            string pars procedures
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
   interface 
       subroutine split(str,delims,before,sep)
         character(len=*),intent(inout) :: str,before
         character(len=*),intent(in) :: delims
         character,optional,intent(inout) :: sep
       end subroutine split
   end interface

      character(len=*),intent(inout) :: str
      character(len=*),intent(in) :: delims
      character(len=len_trim(str)) :: strsav
      character(len=*),dimension(:),intent(inout) :: args
      integer, intent(out) :: nargs
      
      strsav=str
      call compact(str)
      na=size(args)
      do i=1,na
        args(i)=' '
      end do  
      nargs=0
      lenstr=len_trim(str)
      if(lenstr==0) return
      k=0

      do
         if(len_trim(str) == 0) exit
         nargs=nargs+1
         call split(str,delims,args(nargs))
         call removebksl(args(nargs))
      end do   
      str=strsav

   end subroutine parse




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c read coordinates in Angst and converts them to au 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   subroutine pbcrdcoord(fname,lattice,n,xyz,iat,autoang,sita)
   implicit none             
   interface
     subroutine parse(str,delims,args,nargs)
     character(len=*),intent(inout) :: str
     character(len=*),intent(in)  :: delims
     character(len=*),dimension(:),intent(inout) :: args
     integer, intent(out) :: nargs
     end subroutine parse
   end interface
   
   Real*8, Dimension(:,:):: sita(10,1000)
   real*8                :: xyz(3,*),atmnum
   real*8, INTENT(OUT)   ::lattice(3,3)
   integer, INTENT(out)               :: iat(*) 
   integer, INTENT(in)               :: n 
   character*(*), INTENT(IN)          :: fname
   logical              :: selective=.FALSE. ! Selective dynamics
   logical              :: cartesian=.TRUE.  ! Cartesian or direct
   real*8, INTENT(IN)   ::autoang

   real*8 xx(10),scalar
   character*200 line
   character*80 args(90),args2(90)
   
   integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck


   lattice=0
   
   ich=142
   open(unit=ich,file=fname)
   rewind(ich)
   ncheck=0
   ntype=0
   read(ich,'(a)',end=200)line !first line must contain Element Info
   call parse(line,' ',args,ntype)
   read(ich,'(a)',end=200)line !second line contains global scaling factor
   call readl(line,xx,nn)
   scalar=xx(1)/autoang        !the Ang->au conversion is included in the scaling factor
!  write(*,'(F8.6)')scalar
   DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call stoprun( 'Error reading unit cell vectors' )
     lattice(1,i)=xx(1)*scalar
     lattice(2,i)=xx(2)*scalar
     lattice(3,i)=xx(3)*scalar
   !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
   ENDDO
   read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
   line=adjustl(line)
   call readl(line,xx,nn)
   IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
   ENDIF
!       call elem(args(1),i_dummy2)
!       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!          args=args2
!       ENDIF
   IF (nn.NE.ntype ) THEN
     call stoprun( 'Error reading number of atomtypes')
   ENDIF
   ncheck=0
   DO i=1,nn
     i_dummy1=INT(xx(i))
     call elem(args(i),i_dummy2, atmnum)
     IF (i_dummy2<1 .OR. i_dummy2>94) then 
        call stoprun( 'Error: unknown element.')
     end IF
     DO j=1,i_dummy1
       ncheck=ncheck+1
       sita(4,ncheck)=atmnum
       iat(ncheck)=i_dummy2
     ENDDO
   ENDDO
   if (n.ne.ncheck) call stoprun('Error reading Number of Atoms')

   read(ich,'(a)',end=200)line
   line=adjustl(line)
   IF (line(:1).EQ.'s' .OR. line(:1).EQ.'S') THEN
     selective=.TRUE.
     read(ich,'(a)',end=200)line
     line=adjustl(line)
   ENDIF

!   write(*,*)line(:1)
   cartesian=(line(:1).EQ.'c' .OR. line(:1).EQ.'C' .OR. line(:1).EQ.'k' .OR. line(:1).EQ.'K')
   DO i=1,n
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn.NE.3) call stoprun( 'Error reading coordinates.')

     IF (cartesian) THEN
       xyz(1,i)=xx(1)*scalar
       xyz(2,i)=xx(2)*scalar
       xyz(3,i)=xx(3)*scalar
     ELSE
       xyz(1,i)=lattice(1,1)*xx(1)+lattice(1,2)*xx(2)+lattice(1,3)*xx(3)
       xyz(2,i)=lattice(2,1)*xx(1)+lattice(2,2)*xx(2)+lattice(2,3)*xx(3)
       xyz(3,i)=lattice(3,1)*xx(1)+lattice(3,2)*xx(2)+lattice(3,3)*xx(3)
     ENDIF
     sita(1,i)=xyz(1,i)
     sita(2,i)=xyz(2,i)
     sita(3,i)=xyz(3,i)
        
!      write(*,'(3F20.10,1X,I3)')xyz(:,i),iat(i)   !debug printout
      
   ENDDO
      
      
 200  continue

   close(ich)
   end subroutine pbcrdcoord


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Returns the number of a given element string (h-pu, 1-94)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   SUBROUTINE ELEM(KEY1, NATi, atmnum)
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   CHARACTER*(*) KEY1
   CHARACTER*2 ELEMNT(94),E
   Real*8 :: atmnum

   DATA ELEMNT/'h ','he',&
  & 'li','be','b ','c ','n ','o ','f ','ne',&
  & 'na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr',&
  & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
  & 'cd','in','sn','sb','te','i ','xe',&
  & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
  & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
  & 'au','hg','tl','pb','bi','po','at','rn',&
  & 'fr','ra','ac','th','pa','u ','np','pu'/

   nat=0
   e='  '
   k=1
   DO J=1,len(key1)
         if (k.gt.2)exit       
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
            k=k+1
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
   enddo

   DO I=1,94
      if(e.eq.elemnt(i))then
         NAT=I
         RETURN
      ENDIF
   ENDDO

   if (key1.eq.'H' .or. key1.eq.'h') then
     atmnum=1.00
   elseif (key1.eq.'He' .or. key1.eq.'HE' .or. key1.eq.'he') then
     atmnum=2.00
   elseif (key1.eq.'Li' .or. key1.eq.'LI' .or. key1.eq.'li') then
     atmnum=3.00
   elseif (key1.eq.'B' .or. key1.eq.'b') then
     atmnum=5.00
   elseif (key1.eq.'C' .or. key1.eq.'c') then
     atmnum=6.00
   elseif (key1.eq.'N' .or. key1.eq.'n') then
     atmnum=7.00
   elseif (key1.eq.'O' .or. key1.eq.'o') then
     atmnum=8.00
   elseif (key1.eq.'F' .or. key1.eq.'f') then
     atmnum=9.00
   elseif (key1.eq.'Ne' .or. key1.eq.'NE' .or. key1.eq.'ne') then
     atmnum=10.00
   elseif (key1.eq.'Al' .or. key1.eq.'AL' .or. key1.eq.'al') then
     atmnum=13.00
   elseif (key1.eq.'Si' .or. key1.eq.'SI' .or. key1.eq.'si') then
     atmnum=14.00
   elseif (key1.eq.'P' .or. key1.eq.'p') then
     atmnum=15.00
   elseif (key1.eq.'S' .or. key1.eq.'s') then
     atmnum=16.00
   elseif (key1.eq.'Cl' .or. key1.eq.'CL' .or. key1.eq.'cl') then
     atmnum=17.00
   elseif (key1.eq.'Ar' .or. key1.eq.'AR' .or. key1.eq.'ar') then
     atmnum=18.00
   elseif (key1.eq.'Br' .or. key1.eq.'BR' .or. key1.eq.'br') then
     atmnum=35.00
   elseif (key1.eq.'I' .or. key1.eq.'i') then
     atmnum=53.00
   end if
    
   end SUBROUTINE ELEM

!C     *****************************************************************         

   FUNCTION ESYM(I)
   CHARACTER*2 ESYM
   CHARACTER*2 ELEMNT(94)
   DATA ELEMNT/'h ','he',&
  & 'li','be','b ','c ','n ','o ','f ','ne',&
  & 'na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr',&
  & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
  & 'cd','in','sn','sb','te','i ','xe',&
  & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
  & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
  & 'au','hg','tl','pb','bi','po','at','rn',&
  & 'fr','ra','ac','th','pa','u ','np','pu'/
   ESYM=ELEMNT(I)
   RETURN
   END FUNCTION ESYM



   SUBROUTINE SET_CRITERIA(rthr,lat,tau_max)

     REAL*8 :: r_cutoff,rthr
     REAL*8 :: lat(3,3)
     REAL*8 :: tau_max(3)
     REAL*8 :: norm1(3),norm2(3),norm3(3)
     REAL*8 :: cos10,cos21,cos32
     real*8,external :: vectorsize

     r_cutoff=sqrt(rthr)
!       write(*,*) 'lat',lat
       !c find normal to the plane...
       ! kreuzprodukt means cross product
     call kreuzprodukt(lat(:,2),lat(:,3),norm1)
     call kreuzprodukt(lat(:,3),lat(:,1),norm2)
     call kreuzprodukt(lat(:,1),lat(:,2),norm3)
!     write(*,*) 'norm2',norm2
     !c ...normalize it...
     norm1=norm1/VECTORSIZE(norm1)
     norm2=norm2/VECTORSIZE(norm2)
     norm3=norm3/VECTORSIZE(norm3)
!     write(*,*) 'norm2_',norm2
       !c cos angles between normals and lattice vectors
     cos10=SUM(norm1*lat(:,1))
     cos21=SUM(norm2*lat(:,2))
     cos32=SUM(norm3*lat(:,3))
       !write(*,*) 'cos32',cos32
       !tau_max(1)=abs(2*r_cutoff/cos10)
       !tau_max(2)=abs(2*r_cutoff/cos21)
       !tau_max(3)=abs(2*r_cutoff/cos32)
       !write(*,*) 'r_cutoff',r_cutoff
     tau_max(1)=abs(r_cutoff/cos10)
     tau_max(2)=abs(r_cutoff/cos21)
     tau_max(3)=abs(r_cutoff/cos32)
!     write(*,'(3f8.4)')tau_max(1),tau_max(2),tau_max(3)
   END SUBROUTINE SET_CRITERIA


   SUBROUTINE kreuzprodukt(A,B,C)
     IMPLICIT NONE
  
     REAL*8 :: A(3),B(3)
     REAL*8 :: X,Y,Z
     REAL*8 :: C(3)
     
     X=A(2)*B(3)-B(2)*A(3)
     Y=A(3)*B(1)-B(3)*A(1)
     Z=A(1)*B(2)-B(1)*A(2)
     C=(/X,Y,Z/)
   END SUBROUTINE kreuzprodukt

   FUNCTION VECTORSIZE(VECT)

      REAL*8 :: VECT(3)
      REAL*8 :: SVECT(3)
      REAL*8 :: VECTORSIZE

      SVECT=VECT*VECT
      VECTORSIZE=SUM(SVECT)
      VECTORSIZE=VECTORSIZE**(0.5)
   END FUNCTION VECTORSIZE



  
!  These subroutines are copied and modified from Grimme D3 dispersion function
!  This subroutine reads POSCAR (of VASP package) file:
   subroutine pbcrdatomnumber(fname,n)
   implicit none
   interface
     subroutine parse(str,delims,args,nargs)
     character(len=*),intent(inout) :: str
     character(len=*),intent(in)  :: delims
     character(len=*),dimension(:),intent(inout) :: args
     integer, intent(out) :: nargs
     end subroutine parse
   end interface

   integer, INTENT(out)               :: n
   character*(*), INTENT(IN)          :: fname
   logical              :: selective=.FALSE. ! Selective dynamics
   logical              :: cartesian=.TRUE.  ! Cartesian or direct

   real*8 xx(10),scalar,fdum
   character*80 line,args(90),args2(90)

   integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2

   ich=142
   open(unit=ich,file=fname)
   n=0
   ntype=0
   read(ich,'(a)',end=200)line !first line must contain Element Info
   call parse(line,' ',args,ntype)
   read(ich,'(a)',end=200)line !second line contains global scaling factor
   call readl(line,xx,nn)
!   write(*,'(F8.6)')scalar
   DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call stoprun( 'Error reading unit cell vectors' )
   !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
   ENDDO
   read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
   line=adjustl(line)
   call readl(line,xx,nn)
   IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
   ENDIF
!    call elem(args(1),i_dummy2)
!    IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!       args=args2
!    ENDIF
   IF (nn.NE.ntype ) THEN
!      IF(nn.NE.ntype2) THEN
     call stoprun( 'Error reading number of atomtypes')
!      ELSE
!        ntype=ntype2
!      ENDIF
   ENDIF
   n=0
   DO i=1,nn
     i_dummy1=INT(xx(i))
       n=n+i_dummy1
   ENDDO

 200  continue

   close(ich)
   end subroutine pbcrdatomnumber
!
!     *****************************************************************
!     Reads a given line
!     *****************************************************************

   subroutine readl(A1,X,N)
   IMPLICIT REAL*8 (A-H,O-Z)
   CHARACTER*(*) A1
   DIMENSION X(*)
   I=0
   IS=1
10 I=I+1
   X(I)=READAA(A1,IS,IB,IE)
   IF(IB.GT.0 .AND. IE.GT.0) THEN
           IS=IE
           GOTO 10
   ENDIF
   N=I-1
   RETURN
   end subroutine readl

   FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
   IMPLICIT REAL*8 (A-H,O-Z)                                                 
   REAL*8 READAA                                                             
   CHARACTER*(*) A                                                      
   NINE=ICHAR('9')                                                           
   IZERO=ICHAR('0')                                                          
   MINUS=ICHAR('-')                                                          
   IDOT=ICHAR('.')                                                           
   ND=ICHAR('D')                                                             
   NE=ICHAR('E')                                                             
   IBL=ICHAR(' ')                                                            
   IEND=0                                                                    
   IEND2=0                                                                   
   IDIG=0                                                                    
   C1=0                                                                      
   C2=0                                                                      
   ONE=1.D0                                                                  
   X = 1.D0                                                                  
   NL=LEN(A) 
   DO 10 J=ISTART,NL-1                                                       
      N=ICHAR(A(J:J))                                                          
      M=ICHAR(A(J+1:J+1)) 
      IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
      IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20                                                 
10 CONTINUE                                                                  
   READAA=0.D0                                                               
   RETURN                                                                    
20 CONTINUE                                                                  
   IEND=J                                                                    
   DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
30 CONTINUE                                                                  
40 CONTINUE                                                                  
   IDIG=0                                                                    
   DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
50 CONTINUE                                                                  
!                                                                               
! PUT THE PIECES TOGETHER                                                       
!                                                                               
60 CONTINUE                                                                  
   READAA= ONE * ( C1 + C2 * X)                                              
   DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
   RETURN                                                                    
                                                                                
57 C1=0.0D0                                                                  
   ONE=1.0D0                                                                 
   DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
31 CONTINUE                                                                  
61 CONTINUE                                                                  
70 READAA=READAA*10**(ONE*C1)                                                
   RETURN                                                                    
   END FUNCTION READAA




   subroutine stoprun(s)
   character*(*) s
   write(*,*)'program stopped due to: ',s
   call system('touch dscf_problem')
   stop 'must stop!'
   end subroutine stoprun
!
!

!
!

   Subroutine AIM_polarizabilities(sita, mm, Veff3, Veff5, c11, c12, c22, aiw, aiwq)
   use dipole_polarizabilities_mod
   use quad_polarizabilities_mod
   integer, parameter :: vers = 3, nsmax = 1000
   Real*8 :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, &
   d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, &
   d29, d30, d31, dmin, tol, s1, s2, a, b, c, dd, e, f, c11, c12, c22, ratio_v3, ratio_v5
   !Real*8,intent(in) :: Veff33(nsmax), Veff55(nsmax)!, sita(10,nsmax)!, ratio_v3, ratio_v5
   Real*8,intent(in) :: Veff3, Veff5, sita(10,nsmax)
   Integer :: pick, i, j, mm
   Real :: aiw(1,23), aiwq(1,23)
   Real*8, dimension(20) :: AlF3v3, AlH3v3, Arv3, AlCl3v3, BCl3v3, BF3v3, &
   BH3v3, C2H2v3, C2H6v3, C6H6v3, CH3OHv3, CH4v3, CO2v3, H2v3, H2Ov3, H2Sv3, &
   HClv3, NH3v3, PCl3v3, PF3v3, PH3v3, SiCl4v3, SiF4v3, SiH4v3, C2H4v3, &
   HCONH2v3, HCOOHv3, Hev3, HFv3, N2v3, Nev3, C4H4N2v3, CH2Ov3, HBrv3, HIv3, &
   ClFv3, F2v3, Br2v3, CF3Brv3, CF3Iv3, CH3Brv3, CH3Iv3, mOHv3, I2v3, O2v3, &
   C5H12v3, HCNv3, SF2v3, SO2v3, SCl2v3, Si2H4v3, Si2H6v3, Cl2v3, CH3NH2v3, uracilv3
   Real*8, dimension(20) :: AlF3v5, AlH3v5, Arv5, AlCl3v5, BCl3v5, BF3v5, &
   BH3v5, C2H2v5, C2H6v5, C6H6v5, CH3OHv5, CH4v5, CO2v5, H2v5, H2Ov5, H2Sv5, &
   HClv5, NH3v5, PCl3v5, PF3v5, PH3v5, SiCl4v5, SiF4v5, SiH4v5, C2H4v5, &
   HCONH2v5, HCOOHv5, Hev5, HFv5, N2v5, Nev5, C4H4N2v5, CH2Ov5, HBrv5, HIv5, &
   ClFv5, F2v5, Br2v5, CF3Brv5, CF3Iv5, CH3Brv5, CH3Iv5, mOHv5, I2v5, O2v5, &
   C5H12v5, HCNv5, SF2v5, SO2v5, SCl2v5, Si2H4v5, Si2H6v5, Cl2v5, CH3NH2v5, uracilv5
   include 'volume_V3'
   include 'volume_V5'

   a=sita(4,mm)
   b=sita(5,mm)
   c=sita(6,mm)
   dd=sita(7,mm)
   e=sita(8,mm)
   f=sita(9,mm)
   i=1
   tol = 0.2
   !Veff3=Veff33(mm)
   !Veff5=Veff55(mm)
   !print*,'m,Veff3,Veff5',mm,Veff3,Veff5
   if (a.eq.1.0) then
     if (b.ne.0.0) then
                ratio_v3 = Veff3/H2v3(3)
                s1       = Veff3/H2v3(1)
                
                ratio_v5 = Veff5/H2v5(3)
                s2       = Veff5/H2v5(1)
                
                do j=1,23
                   aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * H2(j)
                   aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * H2q(j)
                end do
        ! -------------------- H atom polarizability (v3/v5) --------------------
       ! tol = 0.2
     elseif(.false.) then
        ! match against v3 "te" values (since aiw uses v3 scaling)
        print*,'Veff3,Veff5',Veff3,Veff5
        d1  = abs(Veff3 - H2v3(1))
        d2  = abs(Veff3 - BH3v3(2))
        d3  = abs(Veff3 - HCNv3(1))
        d4  = abs(Veff3 - C4H4N2v3(7))
        d5  = abs(Veff3 - NH3v3(2))
        d6  = abs(Veff3 - HCONH2v3(4))
        d7  = abs(Veff3 - HCOOHv3(2))
        d8  = abs(Veff3 - CH2Ov3(3))
        d9  = abs(Veff3 - C2H4v3(3))
        d10 = abs(Veff3 - C2H2v3(3))
        d11 = abs(Veff3 - C6H6v3(7))
        d12 = abs(Veff3 - CH4v3(2))
        d13 = abs(Veff3 - CH3OHv3(3))
        d14 = abs(Veff3 - C2H6v3(3))
        d15 = abs(Veff3 - C5H12v3(6))
        d16 = abs(Veff3 - H2Ov3(2))
        d17 = abs(Veff3 - HFv3(2))
        d18 = abs(Veff3 - AlH3v3(2))
        d19 = abs(Veff3 - SiH4v3(2))
        d20 = abs(Veff3 - Si2H4v3(3))
        d21 = abs(Veff3 - Si2H6v3(3))
        d22 = abs(Veff3 - PH3v3(2))
        d23 = abs(Veff3 - H2Sv3(2))
        d24 = abs(Veff3 - HBrv3(2))
        d25 = abs(Veff3 - CH3Brv3(3))
        d26 = abs(Veff3 - HIv3(2))
        d27 = abs(Veff3 - CH3Iv3(3))
        d28 = abs(Veff3 - CH3NH2v3(2))
        d29 = abs(Veff3 - CH3NH2v3(5))
        d30 = abs(Veff3 - uracilv3(2))
        d31 = abs(Veff3 - uracilv3(8))
 
        pick = 0
        dmin = huge(1.0)
        
        if (d1  <= tol .and. d1  < dmin) then; dmin = d1 ; pick =  1; end if
        if (d2  <= tol .and. d2  < dmin) then; dmin = d2 ; pick =  2; end if
        if (d3  <= tol .and. d3  < dmin) then; dmin = d3 ; pick =  3; end if
        if (d4  <= tol .and. d4  < dmin) then; dmin = d4 ; pick =  4; end if
        if (d5  <= tol .and. d5  < dmin) then; dmin = d5 ; pick =  5; end if
        if (d6  <= tol .and. d6  < dmin) then; dmin = d6 ; pick =  6; end if
        if (d7  <= tol .and. d7  < dmin) then; dmin = d7 ; pick =  7; end if
        if (d8  <= tol .and. d8  < dmin) then; dmin = d8 ; pick =  8; end if
        if (d9  <= tol .and. d9  < dmin) then; dmin = d9 ; pick =  9; end if
        if (d10 <= tol .and. d10 < dmin) then; dmin = d10; pick = 10; end if
        if (d11 <= tol .and. d11 < dmin) then; dmin = d11; pick = 11; end if
        if (d12 <= tol .and. d12 < dmin) then; dmin = d12; pick = 12; end if
        if (d13 <= tol .and. d13 < dmin) then; dmin = d13; pick = 13; end if
        if (d14 <= tol .and. d14 < dmin) then; dmin = d14; pick = 14; end if
        if (d15 <= tol .and. d15 < dmin) then; dmin = d15; pick = 15; end if
        if (d16 <= tol .and. d16 < dmin) then; dmin = d16; pick = 16; end if
        if (d17 <= tol .and. d17 < dmin) then; dmin = d17; pick = 17; end if
        if (d18 <= tol .and. d18 < dmin) then; dmin = d18; pick = 18; end if
        if (d19 <= tol .and. d19 < dmin) then; dmin = d19; pick = 19; end if
        if (d20 <= tol .and. d20 < dmin) then; dmin = d20; pick = 20; end if
        if (d21 <= tol .and. d21 < dmin) then; dmin = d21; pick = 21; end if
        if (d22 <= tol .and. d22 < dmin) then; dmin = d22; pick = 22; end if
        if (d23 <= tol .and. d23 < dmin) then; dmin = d23; pick = 23; end if
        if (d24 <= tol .and. d24 < dmin) then; dmin = d24; pick = 24; end if
        if (d25 <= tol .and. d25 < dmin) then; dmin = d25; pick = 25; end if
        if (d26 <= tol .and. d26 < dmin) then; dmin = d26; pick = 26; end if
        if (d27 <= tol .and. d27 < dmin) then; dmin = d27; pick = 27; end if
        if (d28 <= tol .and. d28 < dmin) then; dmin = d28; pick = 28; end if
        if (d29 <= tol .and. d29 < dmin) then; dmin = d29; pick = 29; end if
        if (d30 <= tol .and. d30 < dmin) then; dmin = d30; pick = 30; end if
        if (d31 <= tol .and. d31 < dmin) then; dmin = d31; pick = 31; end if
        select case (pick)
        
        case (1)   ! H2      vol(3) te(1)
           ratio_v3=Veff3/H2v3(3) ; s1 = Veff3/H2v3(1)
           ratio_v5=Veff5/H2v5(3) ; s2 = Veff5/H2v5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * H2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * H2q(j)
           end do
        
        case (2)   ! BH3     vol(5) te(1)
           ratio_v3=Veff3/BH3v3(5) ; s1 = Veff3/BH3v3(2)
           ratio_v5=Veff5/BH3v5(5) ; s2 = Veff5/BH3v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * BH3(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * BH3q(j)
           end do
        
        case (3)   ! HCN     vol(4) te(1)
           ratio_v3=Veff3/HCNv3(4) ; s1 = Veff3/HCNv3(1)
           ratio_v5=Veff5/HCNv5(4) ; s2 = Veff5/HCNv5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HCN(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HCNq(j)
           end do
        
        case (4)   ! C4H4N2  vol(11) te(1)
           ratio_v3=Veff3/C4H4N2v3(11) ; s1 = Veff3/C4H4N2v3(7)
           ratio_v5=Veff5/C4H4N2v5(11) ; s2 = Veff5/C4H4N2v5(7)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C4H4N2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C4H4N2q(j)
           end do
        
        case (5)   ! NH3     vol(5) te(1)
           ratio_v3=Veff3/NH3v3(5) ; s1 = Veff3/NH3v3(2)
           ratio_v5=Veff5/NH3v5(5) ; s2 = Veff5/NH3v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * NH3(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * NH3q(j)
           end do
        
        case (6)   ! HCONH2  vol(7) te(1)
           ratio_v3=Veff3/HCONH2v3(7) ; s1 = Veff3/HCONH2v3(4)
           ratio_v5=Veff5/HCONH2v5(7) ; s2 = Veff5/HCONH2v5(4)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HCONH2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HCONH2q(j)
           end do
        
        case (7)   ! HCOOH   vol(6) te(1)
           ratio_v3=Veff3/HCOOHv3(6) ; s1 = Veff3/HCOOHv3(2)
           ratio_v5=Veff5/HCOOHv5(6) ; s2 = Veff5/HCOOHv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HCOOH(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HCOOHq(j)
           end do
        
        case (8)   ! CH2O    vol(5) te(1)
           ratio_v3=Veff3/CH2Ov3(5) ; s1 = Veff3/CH2Ov3(3)
           ratio_v5=Veff5/CH2Ov5(5) ; s2 = Veff5/CH2Ov5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH2O(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH2Oq(j)
           end do
        
        case (9)   ! C2H4    vol(7) te(1)
           ratio_v3=Veff3/C2H4v3(7) ; s1 = Veff3/C2H4v3(3)
           ratio_v5=Veff5/C2H4v5(7) ; s2 = Veff5/C2H4v5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C2H4(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C2H4q(j)
           end do
        
        case (10)  ! C2H2    vol(5) te(1)
           ratio_v3=Veff3/C2H2v3(5) ; s1 = Veff3/C2H2v3(3)
           ratio_v5=Veff5/C2H2v5(5) ; s2 = Veff5/C2H2v5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C2H2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C2H2q(j)
           end do
        
        case (11)  ! C6H6    vol(13) te(1)
           ratio_v3=Veff3/C6H6v3(13) ; s1 = Veff3/C6H6v3(7)
           ratio_v5=Veff5/C6H6v5(13) ; s2 = Veff5/C6H6v5(7)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C6H6(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C6H6q(j)
           end do
        
        case (12)  ! CH4     vol(6) te(1)
           ratio_v3=Veff3/CH4v3(6) ; s1 = Veff3/CH4v3(2)
           ratio_v5=Veff5/CH4v5(6) ; s2 = Veff5/CH4v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH4(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH4q(j)
           end do
        
        case (13)  ! CH3OH   vol(7) te(1)
           ratio_v3=Veff3/CH3OHv3(7) ; s1 = Veff3/CH3OHv3(3)
           ratio_v5=Veff5/CH3OHv5(7) ; s2 = Veff5/CH3OHv5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH3OH(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH3OHq(j)
           end do
        
        case (14)  ! C2H6    vol(9) te(1)
           ratio_v3=Veff3/C2H6v3(9) ; s1 = Veff3/C2H6v3(3)
           ratio_v5=Veff5/C2H6v5(9) ; s2 = Veff5/C2H6v5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C2H6(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C2H6q(j)
           end do
        
        case (15)  ! C5H12   vol(15) te(1)
           ratio_v3=Veff3/C5H12v3(18) ; s1 = Veff3/C5H12v3(6)
           ratio_v5=Veff5/C5H12v5(18) ; s2 = Veff5/C5H12v5(6)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * C5H12(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * C5H12q(j)
           end do
        
        case (16)  ! H2O     vol(4) te(1)
           ratio_v3=Veff3/H2Ov3(4) ; s1 = Veff3/H2Ov3(2)
           ratio_v5=Veff5/H2Ov5(4) ; s2 = Veff5/H2Ov5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * H2O(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * H2Oq(j)
           end do
        
        case (17)  ! HF      vol(3) te(1)
           ratio_v3=Veff3/HFv3(3) ; s1 = Veff3/HFv3(2)
           ratio_v5=Veff5/HFv5(3) ; s2 = Veff5/HFv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HF(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HFq(j)
           end do
        
        case (18)  ! AlH3    vol(5) te(1)
           ratio_v3=Veff3/AlH3v3(5) ; s1 = Veff3/AlH3v3(2)
           ratio_v5=Veff5/AlH3v5(5) ; s2 = Veff5/AlH3v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * AlH3(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * AlH3q(j)
           end do
        
        case (19)  ! SiH4    vol(6) te(1)
           ratio_v3=Veff3/SiH4v3(6) ; s1 = Veff3/SiH4v3(2)
           ratio_v5=Veff5/SiH4v5(6) ; s2 = Veff5/SiH4v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * SiH4(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * SiH4q(j)
           end do

        case (20) ! Si2H4 vol(7)
           ratio_v3=Veff3/Si2H4v3(7); s1 = Veff3/Si2H4v3(3)
           ratio_v5=Veff5/Si2H4v5(7); s2 = Veff5/Si2H4v5(3)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*Si2H4(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*Si2H4q(j)
           end do

        case (21) ! Si2H6 vol(9)
           ratio_v3=Veff3/Si2H6v3(9); s1 = Veff3/Si2H6v3(3)
           ratio_v5=Veff5/Si2H6v5(9); s2 = Veff5/Si2H6v5(3)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*Si2H6(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*Si2H6q(j)
           end do

        case (22) ! PH3  vol(5)
           ratio_v3=Veff3/PH3v3(5); s1 = Veff3/PH3v3(2)
           ratio_v5=Veff5/PH3v5(5); s2 = Veff5/PH3v5(2)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*PH3(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*PH3q(j)
           end do
        
        case (23)  ! H2S     vol(4) te(1)   
           ratio_v3=Veff3/H2Sv3(4) ; s1 = Veff3/H2Sv3(2)
           ratio_v5=Veff5/H2Sv5(4) ; s2 = Veff5/H2Sv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * H2S(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * H2Sq(j)
           end do
        
        case (24)  ! HBr     vol(3) te(1)
           ratio_v3=Veff3/HBrv3(3) ; s1 = Veff3/HBrv3(2)
           ratio_v5=Veff5/HBrv5(3) ; s2 = Veff5/HBrv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HBr(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HBrq(j)
           end do
        
        case (25)  ! CH3Br   vol(6) te(1)
           ratio_v3=Veff3/CH3Brv3(6) ; s1 = Veff3/CH3Brv3(3)
           ratio_v5=Veff5/CH3Brv5(6) ; s2 = Veff5/CH3Brv5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH3Br(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH3Brq(j)
           end do
        
        case (26)  ! HI      vol(3) te(1)
           ratio_v3=Veff3/HIv3(3) ; s1 = Veff3/HIv3(2)
           ratio_v5=Veff5/HIv5(3) ; s2 = Veff5/HIv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * HI(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * HIq(j)
           end do
        
        case (27)  ! CH3I    vol(6) te(1)
           ratio_v3=Veff3/CH3Iv3(6) ; s1 = Veff3/CH3Iv3(3)
           ratio_v5=Veff5/CH3Iv5(6) ; s2 = Veff5/CH3Iv5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH3I(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH3Iq(j)
           end do

        case (28)  ! CH3NH2    vol(6) te(1)
           !print*,'pick CH3Iv3(2) hydrogen'
           ratio_v3=Veff3/CH3NH2v3(8) ; s1 = Veff3/CH3NH2v3(2)
           ratio_v5=Veff5/CH3NH2v5(8) ; s2 = Veff5/CH3NH2v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH3NH2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH3NH2q(j)
           end do

        case (29)  ! CH3NH2    vol(6) te(1)
           !print*,'pick CH3Iv3(5) hydrogen'
           ratio_v3=Veff3/CH3NH2v3(8) ; s1 = Veff3/CH3NH2v3(5)
           ratio_v5=Veff5/CH3NH2v5(8) ; s2 = Veff5/CH3NH2v5(5)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CH3NH2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CH3NH2q(j)
           end do

        case (30)  ! uracil    vol(2) 
           !print*,'pick CH3Iv3(5) hydrogen'
           ratio_v3=Veff3/uracilv3(13) ; s1 = Veff3/uracilv3(2)
           ratio_v5=Veff5/uracilv5(13) ; s2 = Veff5/uracilv5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * uracil(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * uracilq(j)
           end do

        case (31)  ! uracil    vol(8)
           !print*,'pick CH3Iv3(5) hydrogen'
           ratio_v3=Veff3/uracilv3(13) ; s1 = Veff3/uracilv3(8)
           ratio_v5=Veff5/uracilv5(13) ; s2 = Veff5/uracilv5(8)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * uracil(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * uracilq(j)
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(24), vvol3(24), vte3(24), vvol5(24), vte5(24)
           real :: pol(24,23), polq(24,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
           !print*,'picking avg. hydrogen pol.' 
           d = (/ d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24 /)
        
           vvol3 = (/ H2v3(3), BH3v3(5), HCNv3(4), C4H4N2v3(11), NH3v3(5), HCONH2v3(7), HCOOHv3(6), CH2Ov3(5), &
                      C2H4v3(7), C2H2v3(5), C6H6v3(13), CH4v3(6), CH3OHv3(7), C2H6v3(9), C5H12v3(18), H2Ov3(4), &
                      HFv3(3), AlH3v3(5), SiH4v3(6), H2Sv3(4), HBrv3(3), CH3Brv3(6), HIv3(3), CH3Iv3(6) /)
       
           vte3 = (/ H2v3(1), BH3v3(2), HCNv3(1), C4H4N2v3(7), NH3v3(2), HCONH2v3(4), HCOOHv3(2), CH2Ov3(3), &
                  C2H4v3(3), C2H2v3(3),  C6H6v3(7), CH4v3(1), CH3OHv3(3), C2H6v3(3), C5H12v3(6),  H2Ov3(2),  &
                  HFv3(2), AlH3v3(2), SiH4v3(2), H2Sv3(2), HBrv3(2), CH3Brv3(3), HIv3(2), CH3Iv3(3) /)
 
        
           vvol5 = (/ H2v5(3), BH3v5(5), HCNv5(4), C4H4N2v5(11), NH3v5(5), HCONH2v5(7), HCOOHv5(6), CH2Ov5(5), &
                      C2H4v5(7), C2H2v5(5), C6H6v5(13), CH4v5(6), CH3OHv5(7), C2H6v5(9), C5H12v5(18), H2Ov5(4), &
                      HFv5(3), AlH3v5(5), SiH4v5(6), H2Sv5(4), HBrv5(3), CH3Brv5(6), HIv5(3), CH3Iv5(6) /)

           vte5 = (/ H2v5(1), BH3v5(2), HCNv5(1), C4H4N2v5(7), NH3v5(2), HCONH2v5(4), HCOOHv5(2), CH2Ov5(3), &
                  C2H4v5(3), C2H2v5(3),  C6H6v5(7), CH4v5(1), CH3OHv5(3), C2H6v5(3), C5H12v5(6),  H2Ov5(2),  &
                  HFv5(2), AlH3v5(2), SiH4v5(2), H2Sv5(2), HBrv5(2), CH3Brv5(3), HIv5(2), CH3Iv5(3) /)        

           !vte5  = (/ H2v5(1), BH3v5(1), HCNv5(1), C4H4N2v5(1),  NH3v5(1), HCONH2v5(1), HCOOHv5(1), CH2Ov5(1), &
           !           C2H4v5(1), C2H2v5(1), C6H6v5(1), CH4v5(1), CH3OHv5(1), C2H6v5(1), C5H12v5(1), H2Ov5(1), &
           !           HFv5(1), AlH3v5(1), SiH4v5(1), H2Sv5(1), HBrv5(1), CH3Brv5(1), HIv5(1), CH3Iv5(1) /)
        
           do j=1,23
              pol(1,j)=H2(j);       polq(1,j)=H2q(j)
              pol(2,j)=BH3(j);      polq(2,j)=BH3q(j)
              pol(3,j)=HCN(j);      polq(3,j)=HCNq(j)
              pol(4,j)=C4H4N2(j);   polq(4,j)=C4H4N2q(j)
              pol(5,j)=NH3(j);      polq(5,j)=NH3q(j)
              pol(6,j)=HCONH2(j);   polq(6,j)=HCONH2q(j)
              pol(7,j)=HCOOH(j);    polq(7,j)=HCOOHq(j)
              pol(8,j)=CH2O(j);     polq(8,j)=CH2Oq(j)
              pol(9,j)=C2H4(j);     polq(9,j)=C2H4q(j)
              pol(10,j)=C2H2(j);    polq(10,j)=C2H2q(j)
              pol(11,j)=C6H6(j);    polq(11,j)=C6H6q(j)
              pol(12,j)=CH4(j);     polq(12,j)=CH4q(j)
              pol(13,j)=CH3OH(j);   polq(13,j)=CH3OHq(j)
              pol(14,j)=C2H6(j);    polq(14,j)=C2H6q(j)
              pol(15,j)=C5H12(j);   polq(15,j)=C5H12q(j)
              pol(16,j)=H2O(j);     polq(16,j)=H2Oq(j)
              pol(17,j)=HF(j);      polq(17,j)=HFq(j)
              pol(18,j)=AlH3(j);    polq(18,j)=AlH3q(j)
              pol(19,j)=SiH4(j);    polq(19,j)=SiH4q(j)
              pol(20,j)=H2S(j);     polq(20,j)=H2Sq(j)
              pol(21,j)=HBr(j);     polq(21,j)=HBrq(j)
              pol(22,j)=CH3Br(j);   polq(22,j)=CH3Brq(j)
              pol(23,j)=HI(j);      polq(23,j)=HIq(j)
              pol(24,j)=CH3I(j);    polq(24,j)=CH3Iq(j)
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)  = aiw_tmp(1:23)
           aiwq(i,1:23) = aiwq_tmp(1:23)
           end block
        
        end select
! ----------------------------------------------------------------------
     end if

   elseif (a.eq.2.0) then ! He atom
             Do j=1,23
               aiw(i,j)=He(j)
               aiwq(i,j)=Heq(j)
             end Do
   elseif (a.eq.5.0) then ! B atom

       ! tol = 0.2

        ! NOTE: use v3 "te" values for matching (since aiw uses v3)
        d1 = abs(Veff3 - BH3v3(1))
        d2 = abs(Veff3 - BF3v3(1))
        d3 = abs(Veff3 - BCl3v3(1))

        pick = 0
        dmin = huge(1.0)

        if (d1 <= tol .and. d1 < dmin) then
           dmin = d1
           pick = 1
        end if
        if (d2 <= tol .and. d2 < dmin) then
           dmin = d2
           pick = 2
        end if
        if (d3 <= tol .and. d3 < dmin) then
           dmin = d3
           pick = 3
        end if

        select case (pick)

        case (1)   ! BH3
           ratio_v3 = Veff3/(BH3v3(5)-3*BH3v3(2))
           s1       = Veff3/BH3v3(1)

           ratio_v5 = Veff5/(BH3v5(5)-3*BH3v5(2))
           s2       = Veff5/BH3v5(1)

           do j = 1, 23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3*(BH3(j)-1.5*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5*(BH3q(j)-1.5*H2q(j))
           end do

        case (2)   ! BF3
           ratio_v3 = Veff3/BF3v3(5)
           s1       = Veff3/BF3v3(1)

           ratio_v5 = Veff5/BF3v5(5)
           s2       = Veff5/BF3v5(1)

           do j = 1, 23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * BF3(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * BF3q(j)
           end do

        case (3)   ! BCl3
           ratio_v3 = Veff3/BCl3v3(5)
           s1       = Veff3/BCl3v3(1)

           ratio_v5 = Veff5/BCl3v5(5)
           s2       = Veff5/BCl3v5(1)

           do j = 1, 23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * BCl3(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * BCl3q(j)
           end do

        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j

           d     = (/ d1, d2, d3 /)

           vvol3 = (/ (BH3v3(5)-3*BH3v3(2)),  BF3v3(5),  BCl3v3(5) /)
           vte3  = (/ BH3v3(1),  BF3v3(1),  BCl3v3(1) /)

           vvol5 = (/ (BH3v5(5)-3*BH3v5(2)),  BF3v5(5),  BCl3v5(5) /)
           vte5  = (/ BH3v5(1),  BF3v5(1),  BCl3v5(1) /)

           do j = 1, 23
              pol(1,j)  = (BH3(j)-1.5*H2(j))   ; polq(1,j) = (BH3q(j)-1.5*H2q(j))
              pol(2,j)  = BF3(j)   ; polq(2,j) = BF3q(j)
              pol(3,j)  = BCl3(j)  ; polq(3,j) = BCl3q(j)
           end do

           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)

           aiw(i,1:23)  = aiw_tmp(1:23)
           aiwq(i,1:23) = aiwq_tmp(1:23)
           end block
        end select

   elseif (a.eq.6.0) then ! C2 atom
     !if (sitaa(5,i).eq.2.0) then
     if (b.ne.0.0.and.c.ne.0.0.and.dd.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then

        d1 = abs(Veff3 - C2H2v3(2))
        d2 = abs(Veff3 - CO2v3(1))
        d3 = abs(Veff3 - HCNv3(2))
        
        pick = 0
        dmin = huge(1.0)
        
        if (d1 <= tol .and. d1 < dmin) then; dmin = d1; pick = 1; end if
        if (d2 <= tol .and. d2 < dmin) then; dmin = d2; pick = 2; end if
        if (d3 <= tol .and. d3 < dmin) then; dmin = d3; pick = 3; end if
        
        select case (pick)
        
        case (1)   ! C2H2
           ratio_v3 = Veff3/(C2H2v3(5)-2*C2H2v3(3))
           s1       = Veff3/C2H2v3(2)
        
           ratio_v5 = Veff5/(C2H2v5(5)-2*C2H2v5(3))
           s2       = Veff5/C2H2v5(2)
        
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3*(C2H2(j)-H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5*(C2H2q(j)-H2q(j))
           end do
        
        case (2)   ! CO2
           ratio_v3 = Veff3/CO2v3(4)
           s1       = Veff3/CO2v3(1)
        
           ratio_v5 = Veff5/CO2v5(4)
           s2       = Veff5/CO2v5(1)
        
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CO2(j)
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CO2q(j)
           end do
        
        case (3)   ! HCN
           ratio_v3 = Veff3/(HCNv3(4)-HCNv3(1))
           s1       = Veff3/HCNv3(2)
        
           ratio_v5 = Veff5/(HCNv5(4)-HCNv5(1))
           s2       = Veff5/HCNv5(2)
        
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)            *ratio_v3*(HCN(j)-0.5*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCNq(j)-0.5*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1, d2, d3 /)
        
           vvol3 = (/ (C2H2v3(5)-2*C2H2v3(3)), CO2v3(4), (HCNv3(4)-HCNv3(1)) /)
           vte3  = (/ C2H2v3(2), CO2v3(1), HCNv3(2) /)
        
           vvol5 = (/ (C2H2v5(5)-2*C2H2v5(3)), CO2v5(4), (HCNv5(4)-HCNv5(1)) /)
           vte5  = (/ C2H2v5(2), CO2v5(1), HCNv5(2) /)
        
           do j=1,23
              pol(1,j)  = (C2H2(j)-H2(j)) ; polq(1,j) = (C2H2q(j)-H2q(j))
              pol(2,j)  = CO2(j)  ; polq(2,j) = CO2q(j)
              pol(3,j)  = (HCN(j)-0.5*H2(j))  ; polq(3,j) = (HCNq(j)-0.5*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)  = aiw_tmp(1:23)
           aiwq(i,1:23) = aiwq_tmp(1:23)
           end block
        
        end select
     !elseif (sitaa(5,i).eq.3.0) then C3
     elseif (b.ne.0.0.and.c.ne.0.0.and.dd.ne.0.0.and.e.eq.0.0.and.f.eq.0.0) then
        !tol = 0.2
        
        d1 = abs(Veff3 - C2H4v3(2))
        d2 = abs(Veff3 - C6H6v3(1))
        d3 = abs(Veff3 - HCONH2v3(2))
        d4 = abs(Veff3 - C4H4N2v3(1))
        d5 = abs(Veff3 - HCOOHv3(1))
        d6 = abs(Veff3 - CH2Ov3(2))
        d7 = 2.0 !abs(Veff3 - uracilv3(3))
        d8 = abs(Veff3 - uracilv3(7))
        d9 = abs(Veff3 - uracilv3(9))
        d10 = abs(Veff3 - uracilv3(11))
        
        pick = 0
        dmin = huge(1.0)
        
        if (d1 <= tol .and. d1 < dmin) then; dmin = d1; pick = 1; end if
        if (d2 <= tol .and. d2 < dmin) then; dmin = d2; pick = 2; end if
        if (d3 <= tol .and. d3 < dmin) then; dmin = d3; pick = 3; end if
        if (d4 <= tol .and. d4 < dmin) then; dmin = d4; pick = 4; end if
        if (d5 <= tol .and. d5 < dmin) then; dmin = d5; pick = 5; end if
        if (d6 <= tol .and. d6 < dmin) then; dmin = d6; pick = 6; end if
        if (d7 <= tol .and. d7 < dmin) then; dmin = d7; pick = 7; end if
        if (d8 <= tol .and. d8 < dmin) then; dmin = d8; pick = 8; end if
        if (d9 <= tol .and. d9 < dmin) then; dmin = d9; pick = 9; end if
        if (d10 <= tol .and. d10 < dmin) then; dmin = d10; pick = 10; end if
        !write(5,*)'pick',pick
        select case (pick)
        
        case (1)   ! C2H4  vol(7) te(2)
           ratio_v3 = Veff3/(C2H4v3(7)-4*C2H4v3(3));  s1 = Veff3/C2H4v3(2)
           ratio_v5 = Veff5/(C2H4v5(7)-4*C2H4v5(3));  s2 = Veff5/C2H4v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(C2H4(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(C2H4q(j)-2*H2q(j))
           end do
        
        case (2)   ! C6H6  vol(13) te(1)
           ratio_v3 = Veff3/(C6H6v3(13)-6*C6H6v3(7)); s1 = Veff3/C6H6v3(1)
           ratio_v5 = Veff5/(C6H6v5(13)-6*C6H6v5(7)); s2 = Veff5/C6H6v5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(C6H6(j)-3*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(C6H6q(j)-3*H2q(j))
           end do
        
        case (3)   ! HCONH2 vol(7) te(2)
           ratio_v3 = Veff3/(HCONH2v3(7)-3*HCONH2v3(4)); s1 = Veff3/HCONH2v3(2)
           ratio_v5 = Veff5/(HCONH2v5(7)-3*HCONH2v5(4)); s2 = Veff5/HCONH2v5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(HCONH2(j)-1.5*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCONH2q(j)-1.5*H2q(j))
           end do
        
        case (4)   ! C4H4N2 vol(11) te(1)
           ratio_v3 = Veff3/(C4H4N2v3(11)-4*C4H4N2v3(7)); s1 = Veff3/C4H4N2v3(1)
           ratio_v5 = Veff5/(C4H4N2v5(11)-4*C4H4N2v5(7)); s2 = Veff5/C4H4N2v5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(C4H4N2(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(C4H4N2q(j)-2*H2q(j))
           end do
        
        case (5)   ! HCOOH vol(6) te(1)
           ratio_v3 = Veff3/(HCOOHv3(6)-2*HCOOHv3(2)); s1 = Veff3/HCOOHv3(1)
           ratio_v5 = Veff5/(HCOOHv5(6)-2*HCOOHv5(2)); s2 = Veff5/HCOOHv5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (HCOOH(j)-H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (HCOOHq(j)-H2q(j))
           end do
        
        case (6)   ! CH2O vol(5) te(2)
           ratio_v3 = Veff3/(CH2Ov3(5)-2*CH2Ov3(3)); s1 = Veff3/CH2Ov3(2)
           ratio_v5 = Veff5/(CH2Ov5(5)-2*CH2Ov5(3)); s2 = Veff5/CH2Ov5(2)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (CH2O(j)-H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (CH2Oq(j)-H2q(j))
           end do

        case (7)   ! uracil vol(3) 
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(3)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(3)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (uracil(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (uracilq(j)-2*H2q(j))
           end do

        case (8)   ! uracil vol(7)
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(7)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(7)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (uracil(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (uracilq(j)-2*H2q(j))
           end do

        case (9)   ! uracil vol(9)
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(9)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(9)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (uracil(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (uracilq(j)-2*H2q(j))
           end do

        case (10)   ! uracil vol(11)
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(11)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(11)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * (uracil(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * (uracilq(j)-2*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(6), vvol3(6), vte3(6), vvol5(6), vte5(6)
           real :: pol(6,23), polq(6,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
       
           !write(5,*) 'picking avg. for C11' 
           d     = (/ d1, d2, d3, d4, d5, d6 /)
        
           vvol3 = (/ (C2H4v3(7)-4*C2H4v3(3)), (C6H6v3(13)-6*C6H6v3(7)), (HCONH2v3(7)-3*HCONH2v3(4)), &
                   (C4H4N2v3(11)-4*C4H4N2v3(7)), (HCOOHv3(6)-2*HCOOHv3(2)), (CH2Ov3(5)-2*CH2Ov3(3)) /)
           vte3  = (/ C2H4v3(2), C6H6v3(1),  HCONH2v3(2), C4H4N2v3(1),  HCOOHv3(1), CH2Ov3(2) /)
        
           vvol5 = (/ (C2H4v5(7)-4*C2H4v5(3)), (C6H6v5(13)-6*C6H6v5(7)), (HCONH2v5(7)-3*HCONH2v5(4)), &
                   (C4H4N2v5(11)-4*C4H4N2v5(7)), (HCOOHv5(6)-2*HCOOHv5(2)), (CH2Ov5(5)-2*CH2Ov5(3)) /)
           vte5  = (/ C2H4v5(2), C6H6v5(1),  HCONH2v5(2), C4H4N2v5(1),  HCOOHv5(1), CH2Ov5(2) /)
        
           do j=1,23
              pol(1,j)=(C2H4(j)-2*H2(j));   polq(1,j)=(C2H4q(j)-2*H2q(j))
              pol(2,j)=(C6H6(j)-3*H2(j));   polq(2,j)=(C6H6q(j)-3*H2q(j))
              pol(3,j)=(HCONH2(j)-1.5*H2(j)); polq(3,j)=(HCONH2q(j)-1.5*H2q(j))
              pol(4,j)=(C4H4N2(j)-2*H2(j)); polq(4,j)=(C4H4N2q(j)-2*H2q(j))
              pol(5,j)=(HCOOH(j)-H2(j));  polq(5,j)=(HCOOHq(j)-H2q(j))
              pol(6,j)=(CH2O(j)-H2(j));   polq(6,j)=(CH2Oq(j)-H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

     !elseif (sitaa(5,i).eq.4.0) then C4
     elseif (b.ne.0.0.and.c.ne.0.0.and.dd.ne.0.0.and.e.ne.0.0.and.f.eq.0.0) then

        !tol = 0.2
        
        d1 = abs(Veff3 - C2H6v3(2))
        d2 = abs(Veff3 - C5H12v3(1))
        d3 = abs(Veff3 - C5H12v3(2))
        d4 = abs(Veff3 - CH4v3(1))
        d5 = abs(Veff3 - CH3OHv3(1))
        d6 = abs(Veff3 - CF3Brv3(2))
        d7 = abs(Veff3 - CF3Iv3(2))
        d8 = abs(Veff3 - CH3Brv3(2))
        d9 = abs(Veff3 - CH3Iv3(2))
        d10 =abs(Veff3 - CH3NH2v3(4))
        
        pick = 0
        dmin = huge(1.0)
        
        if (d1 <= tol .and. d1 < dmin) then; dmin=d1; pick=1; end if
        if (d2 <= tol .and. d2 < dmin) then; dmin=d2; pick=2; end if
        if (d3 <= tol .and. d3 < dmin) then; dmin=d3; pick=3; end if
        if (d4 <= tol .and. d4 < dmin) then; dmin=d4; pick=4; end if
        if (d5 <= tol .and. d5 < dmin) then; dmin=d5; pick=5; end if
        if (d6 <= tol .and. d6 < dmin) then; dmin=d6; pick=6; end if
        if (d7 <= tol .and. d7 < dmin) then; dmin=d7; pick=7; end if
        if (d8 <= tol .and. d8 < dmin) then; dmin=d8; pick=8; end if
        if (d9 <= tol .and. d9 < dmin) then; dmin=d9; pick=9; end if
        if (d10 <= tol .and. d10 < dmin) then; dmin=d10; pick=10; end if

        select case (pick)
        
        case (1)   ! C2H6 vol(9) te(2)
           ratio_v3=Veff3/(C2H6v3(9)-6*C2H6v3(3)); s1 = Veff3/C2H6v3(2)
           ratio_v5=Veff5/(C2H6v5(9)-6*C2H6v5(3)); s2 = Veff5/C2H6v5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(C2H6(j)-3*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(C2H6q(j)-3*H2q(j))
           end do
        
        case (2)   ! C5H12 C1 vol(15) te(1)
           ratio_v3=Veff3/(C5H12v3(18)-12*C5H12v3(6)); s1 = Veff3/C5H12v3(1)
           ratio_v5=Veff5/(C5H12v5(18)-12*C5H12v5(6)); s2 = Veff5/C5H12v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(C5H12(j)-6*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(C5H12q(j)-6*H2q(j))
           end do
        
        case (3)   ! C5H12 C2 vol(15) te(2)
           ratio_v3=Veff3/(C5H12v3(18)-12*C5H12v3(6)); s1 = Veff3/C5H12v3(2)
           ratio_v5=Veff5/(C5H12v5(18)-12*C5H12v5(6)); s2 = Veff5/C5H12v5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(C5H12(j)-6*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(C5H12q(j)-6*H2q(j))
           end do
        
        case (4)   ! CH4 vol(6) te(1)
           ratio_v3=Veff3/(CH4v3(6)-4*CH4v3(2)); s1 = Veff3/CH4v3(1)
           ratio_v5=Veff5/(CH4v5(6)-4*CH4v5(2)); s2 = Veff5/CH4v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(CH4(j)-2*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(CH4q(j)-2*H2q(j))
           end do
        
        case (5)   ! CH3OH vol(7) te(1)
           ratio_v3=Veff3/(CH3OHv3(7)-4*CH3OHv3(3)); s1 = Veff3/CH3OHv3(1)
           ratio_v5=Veff5/(CH3OHv5(7)-4*CH3OHv5(3)); s2 = Veff5/CH3OHv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(CH3OH(j)-2*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(CH3OHq(j)-2*H2q(j))
           end do
        
        case (6)   ! CF3Br vol(6) te(2)
           ratio_v3=Veff3/CF3Brv3(6); s1 = Veff3/CF3Brv3(2)
           ratio_v5=Veff5/CF3Brv5(6); s2 = Veff5/CF3Brv5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*CF3Br(j)
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*CF3Brq(j)
           end do
        
        case (7)   ! CF3I vol(6) te(2)
           ratio_v3=Veff3/CF3Iv3(6); s1 = Veff3/CF3Iv3(2)
           ratio_v5=Veff5/CF3Iv5(6); s2 = Veff5/CF3Iv5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*CF3I(j)
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*CF3Iq(j)
           end do
        
        case (8)   ! CH3Br vol(6) te(2)
           ratio_v3=Veff3/(CH3Brv3(6)-3*CH3Brv3(3)); s1 = Veff3/CH3Brv3(2)
           ratio_v5=Veff5/(CH3Brv5(6)-3*CH3Brv5(3)); s2 = Veff5/CH3Brv5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(CH3Br(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(CH3Brq(j)-1.5*H2q(j))
           end do
        
        case (9)   ! CH3I vol(6) te(2)
           ratio_v3=Veff3/(CH3Iv3(6)-3*CH3Iv3(3)); s1 = Veff3/CH3Iv3(2)
           ratio_v5=Veff5/(CH3Iv5(6)-3*CH3Iv5(3)); s2 = Veff5/CH3Iv5(2)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(CH3I(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(CH3Iq(j)-1.5*H2q(j))
           end do

        case (10)   ! CH3NH2 vol(4) te(2)
           !print*,'pick CH3NH2v3(4) carbon'
           ratio_v3=Veff3/(CH3NH2v3(8)-5*CH3NH2v3(2)); s1 = Veff3/CH3NH2v3(4)
           ratio_v5=Veff5/(CH3NH2v5(8)-5*CH3NH2v5(2)); s2 = Veff5/CH3NH2v5(4)
           do j=1,23
              aiw(i,j)=abs(1.0-c11*s1)*ratio_v3*(CH3NH2(j)-2.5*H2(j))
              aiwq(i,j)=abs(1.0-c12*s2-c22*s2**2)*ratio_v5*(CH3NH2q(j)-2.5*H2q(j))
           end do

        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(9), vvol3(9), vte3(9), vvol5(9), vte5(9)
           real :: pol(9,23), polq(9,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
           !write(5,*),'pick avg. carbon pol.' 
           d     = (/ d1,d2,d3,d4,d5,d6,d7,d8,d9 /)
        
           vvol3 = (/ (C2H6v3(9)-6*C2H6v3(3)), (C5H12v3(18)-12*C5H12v3(6)), (C5H12v3(18)-12*C5H12v3(6)), (CH4v3(6)-4*CH4v3(2)),&
                   (CH3OHv3(7)-4*CH3OHv3(3)), CF3Brv3(6), CF3Iv3(6), (CH3Brv3(6)-3*CH3Brv3(3)), (CH3Iv3(6)-3*CH3Iv3(3)) /)
           vte3  = (/ C2H6v3(2), C5H12v3(1),  C5H12v3(2),  CH4v3(1), CH3OHv3(1), &
                      CF3Brv3(2), CF3Iv3(2),       CH3Brv3(2),       CH3Iv3(2) /)
        
           vvol5 = (/ (C2H6v5(9)-6*C2H6v5(3)), (C5H12v5(18)-12*C5H12v5(6)), (C5H12v5(18)-12*C5H12v5(6)), (CH4v5(6)-4*CH4v5(2)),&
                    (CH3OHv5(7)-4*CH3OHv5(3)), CF3Brv5(6), CF3Iv5(6), (CH3Brv5(6)-3*CH3Brv5(3)), (CH3Iv5(6)-3*CH3Iv5(3)) /)
           vte5  = (/ C2H6v5(2), C5H12v5(1),  C5H12v5(2),  CH4v5(1), CH3OHv5(1), &
                      CF3Brv5(2), CF3Iv5(2),       CH3Brv5(2),       CH3Iv5(2) /)
        
           do j=1,23
              pol(1,j)=(C2H6(j)-3*H2(j));       polq(1,j)=(C2H6q(j)-3*H2q(j))
              pol(2,j)=(C5H12(j)-6*H2(j));      polq(2,j)=(C5H12q(j)-6*H2q(j))
              pol(3,j)=(C5H12(j)-6*H2(j));      polq(3,j)=(C5H12q(j)-6*H2q(j))
              pol(4,j)=(CH4(j)-2*H2(j));        polq(4,j)=(CH4q(j)-2*H2q(j))
              pol(5,j)=(CH3OH(j)-2*H2(j));      polq(5,j)=(CH3OHq(j)-2*H2q(j))
              pol(6,j)=CF3Br(j);      polq(6,j)=CF3Brq(j)
              pol(7,j)=CF3I(j);       polq(7,j)=CF3Iq(j)
              pol(8,j)=(CH3Br(j)-1.5*H2(j)); polq(8,j)=(CH3Brq(j)-1.5*H2q(j))
              pol(9,j)=(CH3I(j)-1.5*H2(j));  polq(9,j)=(CH3Iq(j)-1.5*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select
     end if

   elseif (a.eq.7.0) then ! N1 and N2 atom
     !if (sitaa(5,i).eq.1.0) then
     if (b.ne.0.0.and.c.eq.0.0.and.dd.eq.0.0.and.e.eq.0.0.and.f.eq.0.0.or.b.ne.0.0.and.c.ne.0.0.and.dd.eq.0.0 &
        .and.e.eq.0.0.and.f.eq.0.0) then
       ! tol = 0.2
        
        d1 = abs(Veff3 - N2v3(1))
        d2 = abs(Veff3 - HCNv3(3))
        d3 = abs(Veff3 - C4H4N2v3(5))
        
        pick = 0
        dmin = huge(1.0)
        
        if (d1 <= tol .and. d1 < dmin) then; dmin=d1; pick=1; end if
        if (d2 <= tol .and. d2 < dmin) then; dmin=d2; pick=2; end if
        if (d3 <= tol .and. d3 < dmin) then; dmin=d3; pick=3; end if
        
        select case (pick)
        
        case (1)   ! N2 vol(3) te(1)
           ratio_v3=Veff3/N2v3(3); s1 = Veff3/N2v3(1)
           ratio_v5=Veff5/N2v5(3); s2 = Veff5/N2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*N2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*N2q(j)
           end do
        
        case (2)   ! HCN vol(4) te(3)
           ratio_v3=Veff3/(HCNv3(4)-HCNv3(1)); s1 = Veff3/HCNv3(3)
           ratio_v5=Veff5/(HCNv5(4)-HCNv5(1)); s2 = Veff5/HCNv5(3)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HCN(j)-0.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCNq(j)-0.5*H2q(j))
           end do

        case (3)   ! C4H4N2 vol(11) te(1)
           ratio_v3 = Veff3/(C4H4N2v3(11)-4*C4H4N2v3(7)); s1 = Veff3/C4H4N2v3(5)
           ratio_v5 = Veff5/(C4H4N2v5(11)-4*C4H4N2v5(7)); s2 = Veff5/C4H4N2v5(5)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(C4H4N2(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(C4H4N2q(j)-2*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
           !write(5,*) 'Picking Nitrogen avg. polarizability' 
           d     = (/ d1, d2, d3/)
        
           vvol3 = (/ N2v3(3), (HCNv3(4)-HCNv3(1)), (C4H4N2v3(11)-4*C4H4N2v3(7))/)
           vte3  = (/ N2v3(1), HCNv3(3), C4H4N2v3(5) /)
        
           vvol5 = (/ N2v5(3), (HCNv5(4)-HCNv5(1)), (C4H4N2v5(11)-4*C4H4N2v5(7))/)
           vte5  = (/ N2v5(1), HCNv5(3), C4H4N2v5(5) /)
        
           do j=1,23
              pol(1,j)=N2(j);  polq(1,j)=N2q(j)
              pol(2,j)=(HCN(j)-0.5*H2(j)); polq(2,j)=(HCNq(j)-0.5*H2q(j))
              pol(3,j)=(C4H4N2(j)-2*H2(j)); polq(3,j)=(C4H4N2q(j)-2*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           !write(5,*)'aiw(i,1),aiwq(i,1)',aiw(i,1),aiwq(i,1)
           end block
        
        end select

     !if (sitaa(5,i).eq.3.0) then N3
     elseif (b.ne.0.0.and.c.ne.0.0.and.dd.ne.0.0.and.e.eq.0.0.and.f.eq.0.0) then

        !tol = 0.2
        !print*,'picking Nitrogent polariz.'
        d1 = abs(Veff3 - NH3v3(1))
        d2 = abs(Veff3 - HCONH2v3(3))
        d3 = abs(Veff3 - CH3NH2v3(1))
        d4 = abs(Veff3 - uracilv3(1))
        d4 = abs(Veff3 - uracilv3(5))

        pick = 0
        dmin = huge(1.0)
        
        if (d1 <= tol .and. d1 < dmin) then; dmin=d1; pick=1; end if
        if (d2 <= tol .and. d2 < dmin) then; dmin=d2; pick=2; end if
        if (d3 <= tol .and. d3 < dmin) then; dmin=d3; pick=3; end if
        if (d4 <= tol .and. d4 < dmin) then; dmin=d4; pick=4; end if
        if (d5 <= tol .and. d5 < dmin) then; dmin=d5; pick=5; end if

        select case (pick)
        
        case (1)   ! NH3 vol(5) te(1)
           ratio_v3=Veff3/(NH3v3(5)-3*NH3v3(2)); s1 = Veff3/NH3v3(1)
           ratio_v5=Veff5/(NH3v5(5)-3*NH3v5(2)); s2 = Veff5/NH3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(NH3(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(NH3q(j)-1.5*H2q(j))
           end do
        
        case (2)   ! HCONH2 vol(7) te(3)
           ratio_v3=Veff3/(HCONH2v3(7)-3*HCONH2v3(4)); s1 = Veff3/HCONH2v3(3)
           ratio_v5=Veff5/(HCONH2v5(7)-3*HCONH2v5(4)); s2 = Veff5/HCONH2v5(3)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HCONH2(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCONH2q(j)-1.5*H2q(j))
           end do

        case (3)   ! CH3NH2 vol(7) te(3)
           !print*,'pick CH3NH2v3(1) Nitrogen'
           ratio_v3=Veff3/(CH3NH2v3(8)-5*CH3NH2v3(2)); s1 = Veff3/CH3NH2v3(1)
           ratio_v5=Veff5/(CH3NH2v5(8)-5*CH3NH2v5(2)); s2 = Veff5/CH3NH2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(CH3NH2(j)-2.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(CH3NH2q(j)-2.5*H2q(j))
           end do

        case (4)   ! uracil vol(1)
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(1)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(1)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(uracil(j)-2.0*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(uracilq(j)-2.0*H2q(j))
           end do

        case (5)   ! uracil vol(5)
           ratio_v3 = Veff3/(uracilv3(13)-4.0*uracilv3(2)); s1 = Veff3/uracilv3(5)
           ratio_v5 = Veff5/(uracilv5(13)-4.0*uracilv5(2)); s2 = Veff5/uracilv5(5)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(uracil(j)-2.0*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(uracilq(j)-2.0*H2q(j))
           end do

        case default

           block
           use polar_interp, only: polar_weighted_avg
           real :: d(2), vvol3(2), vte3(2), vvol5(2), vte5(2)
           real :: pol(2,23), polq(2,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
           !write(5,*)'picking Nitrogent average polariz.'
           d     = (/ d1, d2 /)
        
           vvol3 = (/ (NH3v3(5)-3*NH3v3(2)), (HCONH2v3(7)-3*HCONH2v3(4)) /)
           vte3  = (/ NH3v3(1), HCONH2v3(3) /)
        
           vvol5 = (/ (NH3v5(5)-3*NH3v5(2)), (HCONH2v5(7)-3*HCONH2v5(4)) /)
           vte5  = (/ NH3v5(1), HCONH2v5(3) /)
        
           do j=1,23
              pol(1,j)=(NH3(j)-1.5*H2(j));     polq(1,j)=(NH3q(j)-1.5*H2q(j))
              pol(2,j)=(HCONH2(j)-1.5*H2(j));  polq(2,j)=(HCONH2q(j)-1.5*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           !write(5,*)'aiw(i,1),aiwq(i,1)',aiw(i,1),aiwq(i,1)
           end block
        
        end select
     end if

   elseif (a.eq.8.0) then ! O1 atom
     !if (sitaa(5,i).eq.1.0) then
     if (b.ne.0.0.and.c.eq.0.0.and.dd.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then
        !tol = 0.2
        
        d1 = abs(Veff3 - CO2v3(2))
        d2 = abs(Veff3 - HCONH2v3(1))
        d3 = abs(Veff3 - HCOOHv3(3))
        d4 = abs(Veff3 - CH2Ov3(1))
        d5 = abs(Veff3 - uracilv3(4))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        if (d4<=tol .and. d4<dmin) then; dmin=d4; pick=4; end if
        if (d5<=tol .and. d5<dmin) then; dmin=d5; pick=5; end if
        
        select case (pick)
        
        case (1) ! CO2 vol(4) te(2)
           ratio_v3=Veff3/CO2v3(4); s1 = Veff3/CO2v3(2)
           ratio_v5=Veff5/CO2v5(4); s2 = Veff5/CO2v5(2)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*CO2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*CO2q(j)
           end do
        
        case (2) ! HCONH2 vol(7) te(1)
           ratio_v3=Veff3/(HCONH2v3(7)-5*CH3NH2v3(2)); s1 = Veff3/HCONH2v3(1)
           ratio_v5=Veff5/(HCONH2v5(7)-5*CH3NH2v5(2)); s2 = Veff5/HCONH2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HCONH2(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCONH2q(j)-1.5*H2q(j))
           end do
        
        case (3) ! HCOOH vol(6) te(3)
           ratio_v3=Veff3/(HCOOHv3(6)-2*HCOOHv3(2)); s1 = Veff3/HCOOHv3(3)
           ratio_v5=Veff5/(HCOOHv5(6)-2*HCOOHv5(2)); s2 = Veff5/HCOOHv5(3)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HCOOH(j)-H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCOOHq(j)-H2q(j))
           end do
        
        case (4) ! CH2O vol(5) te(1)
           ratio_v3=Veff3/(CH2Ov3(5)-2*CH2Ov3(3)); s1 = Veff3/CH2Ov3(1)
           ratio_v5=Veff5/(CH2Ov5(5)-2*CH2Ov5(3)); s2 = Veff5/CH2Ov5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(CH2O(j)-H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(CH2Oq(j)-H2q(j))
           end do

        case (5)   ! uracil vol(4)
           ratio_v3 = Veff3/(uracilv3(13)-4*uracilv3(2)); s1 = Veff3/uracilv3(4)
           ratio_v5 = Veff5/(uracilv5(13)-4*uracilv5(2)); s2 = Veff5/uracilv5(4)
           do j=1,23
              aiw(i,j)  = abs(1.0 - c11*s1)*ratio_v3*(uracil(j)-2*H2(j))
              aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(uracilq(j)-2*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(4), vvol3(4), vte3(4), vvol5(4), vte5(4)
           real :: pol(4,23), polq(4,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3,d4 /)
        
           vvol3 = (/ CO2v3(4), (HCONH2v3(7)-5*CH3NH2v3(2)), (HCOOHv3(6)-2*HCOOHv3(2)), (CH2Ov3(5)-2*CH2Ov3(3)) /)
           vte3  = (/ CO2v3(2), HCONH2v3(1), HCOOHv3(3), CH2Ov3(1) /)
        
           vvol5 = (/ CO2v5(4), (HCONH2v5(7)-5*CH3NH2v5(2)), (HCOOHv5(6)-2*HCOOHv5(2)), (CH2Ov5(5)-2*CH2Ov5(3)) /)
           vte5  = (/ CO2v5(2), HCONH2v5(1), HCOOHv5(3), CH2Ov5(1) /)
        
           do j=1,23
              pol(1,j)=CO2(j);     polq(1,j)=CO2q(j)
              pol(2,j)=(HCONH2(j)-1.5*H2(j));  polq(2,j)=(HCONH2q(j)-1.5*H2q(j))
              pol(3,j)=(HCOOH(j)-H2(j));   polq(3,j)=(HCOOHq(j)-H2q(j))
              pol(4,j)=(CH2O(j)-H2(j));    polq(4,j)=(CH2Oq(j)-H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select
        
     !elseif (sitaa(5,i).eq.2.0) then O2
     elseif (b.ne.0.0.and.c.ne.0.0.and.dd.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then
        
       ! tol = 0.2
        
        d1 = abs(Veff3 - H2Ov3(1))
        d2 = abs(Veff3 - CH3OHv3(2))
        d3 = abs(Veff3 - HCOOHv3(4))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        
        select case (pick)
        
        case (1) ! H2O vol(4) te(1)
           ratio_v3=Veff3/(H2Ov3(4)-2*H2Ov3(2)); s1 = Veff3/H2Ov3(1)
           ratio_v5=Veff5/(H2Ov5(4)-2*H2Ov5(2)); s2 = Veff5/H2Ov5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(H2O(j)-H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(H2Oq(j)-H2q(j))
           end do
        
        case (2) ! CH3OH vol(7) te(2)
           ratio_v3=Veff3/(CH3OHv3(7)-4*CH3OHv3(3)); s1 = Veff3/CH3OHv3(2)
           ratio_v5=Veff5/(CH3OHv5(7)-4*CH3OHv5(3)); s2 = Veff5/CH3OHv5(2)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(CH3OH(j)-2*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(CH3OHq(j)-2*H2q(j))
           end do
        
        case (3) ! HCOOH vol(6) te(4)
           ratio_v3=Veff3/(HCOOHv3(6)-2*HCOOHv3(2)); s1 = Veff3/HCOOHv3(4)
           ratio_v5=Veff5/(HCOOHv5(6)-2*HCOOHv5(2)); s2 = Veff5/HCOOHv5(4)
           !write(5,*)'ratio_v3,ratio_v5,s1,s2',ratio_v3,ratio_v5,s1,s2
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HCOOH(j)-H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HCOOHq(j)-H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
           !print*,'avg. value for O2'
           d     = (/ d1,d2,d3 /)
        
           vvol3 = (/ (H2Ov3(4)-2*H2Ov3(2)), (CH3OHv3(7)-4*CH3OHv3(3)), (HCOOHv3(6)-2*HCOOHv3(2)) /)
           vte3  = (/ H2Ov3(1), CH3OHv3(2), HCOOHv3(4) /)
        
           vvol5 = (/ (H2Ov5(4)-2*H2Ov5(2)), (CH3OHv5(7)-4*CH3OHv5(3)), (HCOOHv5(6)-2*HCOOHv5(2)) /)
           vte5  = (/ H2Ov5(1), CH3OHv5(2), HCOOHv5(4) /)
        
           do j=1,23
              pol(1,j)=(H2O(j)-H2(j));    polq(1,j)=(H2Oq(j)-H2q(j))
              pol(2,j)=(CH3OH(j)-2*H2(j));  polq(2,j)=(CH3OHq(j)-2*H2q(j))
              pol(3,j)=(HCOOH(j)-H2(j));  polq(3,j)=(HCOOHq(j)-H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select
     end if

   elseif (a.eq.9.0) then ! F atom
     if (b.eq.1.0) then

        ratio_v3 = Veff3/(HFv3(3)-HFv3(2))
        s1       = Veff3/HFv3(1)
        
        ratio_v5 = Veff5/(HFv5(3)-HFv5(2))
        s2       = Veff5/HFv5(1)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3*(HF(j)-0.5*H2(j))
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5*(HFq(j)-0.5*H2q(j))
        end do
        
     elseif (b.eq.5.0) then        
        
        ratio_v3 = Veff3/BF3v3(5)
        s1       = Veff3/BF3v3(2)
        
        ratio_v5 = Veff5/BF3v5(5)
        s2       = Veff5/BF3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * BF3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * BF3q(j)
        end do
        
     elseif (b.eq.6.0) then
        
        ratio_v3 = Veff3/CF3Brv3(6)
        s1       = Veff3/CF3Brv3(3)
        
        ratio_v5 = Veff5/CF3Brv5(6)
        s2       = Veff5/CF3Brv5(3)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * CF3Br(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * CF3Brq(j)
        end do
        
     elseif (b.eq.9.0) then
        
        ratio_v3 = Veff3/F2v3(3)
        s1       = Veff3/F2v3(1)
        
        ratio_v5 = Veff5/F2v5(3)
        s2       = Veff5/F2v5(1)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * F2(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * F2q(j)
        end do
        
     elseif (b.eq.13.0) then
        
        ratio_v3 = Veff3/AlF3v3(5)
        s1       = Veff3/AlF3v3(2)
        
        ratio_v5 = Veff5/AlF3v5(5)
        s2       = Veff5/AlF3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * AlF3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * AlF3q(j)
        end do
        
     elseif (b.eq.14.0) then
        
        ratio_v3 = Veff3/SiF4v3(6)
        s1       = Veff3/SiF4v3(2)
        
        ratio_v5 = Veff5/SiF4v5(6)
        s2       = Veff5/SiF4v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * SiF4(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * SiF4q(j)
        end do
        
     elseif (b.eq.15.0) then
        
        ratio_v3 = Veff3/PF3v3(5)
        s1       = Veff3/PF3v3(2)
        
        ratio_v5 = Veff5/PF3v5(5)
        s2       = Veff5/PF3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * PF3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * PF3q(j)
        end do
        
     elseif (b.eq.17.0) then
        
        ratio_v3 = Veff3/ClFv3(3)
        s1       = Veff3/ClFv3(2)
        
        ratio_v5 = Veff5/ClFv5(3)
        s2       = Veff5/ClFv5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * ClF(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * ClFq(j)
        end do
     end if


   elseif (a.eq.10.0) then ! Ne atom
             Do j=1,23
               aiw(i,j)=Ne(j)
               aiwq(i,j)=Neq(j)
             end Do

   elseif (a.eq.13.0) then ! Al atom

       ! tol = 0.2
        
        d1 = abs(Veff3 - AlH3v3(1))
        d2 = abs(Veff3 - AlF3v3(1))
        d3 = abs(Veff3 - AlCl3v3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        
        select case (pick)
        
        case (1) ! AlH3 vol(5) te(1)
           ratio_v3=Veff3/(AlH3v3(5)-3*AlH3v3(2)); s1 = Veff3/AlH3v3(1)
           ratio_v5=Veff5/(AlH3v5(5)-3*AlH3v5(2)); s2 = Veff5/AlH3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(AlH3(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(AlH3q(j)-1.5*H2q(j))
           end do
        
        case (2) ! AlF3
           ratio_v3=Veff3/AlF3v3(5); s1 = Veff3/AlF3v3(1)
           ratio_v5=Veff5/AlF3v5(5); s2 = Veff5/AlF3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*AlF3(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*AlF3q(j)
           end do
        
        case (3) ! AlCl3
           ratio_v3=Veff3/AlCl3v3(5); s1 = Veff3/AlCl3v3(1)
           ratio_v5=Veff5/AlCl3v5(5); s2 = Veff5/AlCl3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*AlCl3(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*AlCl3q(j)
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3 /)
        
           vvol3 = (/ (AlH3v3(5)-3*AlH3v3(2)), AlF3v3(5), AlCl3v3(5) /)
           vte3  = (/ AlH3v3(1), AlF3v3(1), AlCl3v3(1) /)
        
           vvol5 = (/ (AlH3v5(5)-3*AlH3v5(2)), AlF3v5(5), AlCl3v5(5) /)
           vte5  = (/ AlH3v5(1), AlF3v5(1), AlCl3v5(1) /)
        
           do j=1,23
              pol(1,j)=(AlH3(j)-1.5*H2(j));  polq(1,j)=(AlH3q(j)-1.5*H2q(j))
              pol(2,j)=AlF3(j);  polq(2,j)=AlF3q(j)
              pol(3,j)=AlCl3(j); polq(3,j)=AlCl3q(j)
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

   elseif (a.eq.14.0) then ! Si atom

       ! tol = 0.2
        
        d1 = abs(Veff3 - SiH4v3(1))
        d2 = abs(Veff3 - SiF4v3(1))
        d3 = abs(Veff3 - SiCl4v3(1))
        d4 = abs(Veff3 - Si2H4v3(1))
        d5 = abs(Veff3 - Si2H6v3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        if (d4<=tol .and. d4<dmin) then; dmin=d4; pick=4; end if
        if (d5<=tol .and. d5<dmin) then; dmin=d5; pick=5; end if
        
        select case (pick)
        
        case (1) ! SiH4 vol(6) te(1)
           ratio_v3=Veff3/(SiH4v3(6)-4*SiH4v3(2)); s1 = Veff3/SiH4v3(1)
           ratio_v5=Veff5/(SiH4v5(6)-4*SiH4v5(2)); s2 = Veff5/SiH4v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(SiH4(j)-2*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(SiH4q(j)-2*H2q(j))
           end do
        
        case (2) ! SiF4
           ratio_v3=Veff3/SiF4v3(6); s1 = Veff3/SiF4v3(1)
           ratio_v5=Veff5/SiF4v5(6); s2 = Veff5/SiF4v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*SiF4(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*SiF4q(j)
           end do
        
        case (3) ! SiCl4
           ratio_v3=Veff3/SiCl4v3(6); s1 = Veff3/SiCl4v3(1)
           ratio_v5=Veff5/SiCl4v5(6); s2 = Veff5/SiCl4v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*SiCl4(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*SiCl4q(j)
           end do

        case (4) ! Si2H4
           ratio_v3=Veff3/(Si2H4v3(7)-4*Si2H4v3(3)); s1 = Veff3/Si2H4v3(1)
           ratio_v5=Veff5/(Si2H4v5(7)-4*Si2H4v5(3)); s2 = Veff5/Si2H4v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(Si2H4(j)-2*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(Si2H4q(j)-2*H2q(j))
           end do

        case (5) ! Si2H6
           ratio_v3=Veff3/(Si2H6v3(9)-6*Si2H6v3(3)); s1 = Veff3/Si2H6v3(1)
           ratio_v5=Veff5/(Si2H6v5(9)-6*Si2H6v5(3)); s2 = Veff5/Si2H6v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(Si2H6(j)-3*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(Si2H6q(j)-3*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(5), vvol3(5), vte3(5), vvol5(5), vte5(5)
           real :: pol(5,23), polq(5,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3,d4,d5 /)
        
           vvol3 = (/ (SiH4v3(6)-4*SiH4v3(2)), SiF4v3(6), SiCl4v3(6), (Si2H4v3(7)-4*Si2H4v3(3)), (Si2H6v3(9)-6*Si2H6v3(3))/)
           vte3  = (/ SiH4v3(1), SiF4v3(1), SiCl4v3(1), Si2H4v3(1), Si2H6v3(1)/)
        
           vvol5 = (/ (SiH4v5(6)-4*SiH4v5(2)), SiF4v5(6), SiCl4v5(6), (Si2H4v5(7)-4*Si2H4v5(3)), (Si2H6v5(9)-6*Si2H6v5(3))/)
           vte5  = (/ SiH4v5(1), SiF4v5(1), SiCl4v5(1), Si2H4v5(1), Si2H6v5(1)/)
        
           do j=1,23
              pol(1,j)=(SiH4(j)-2*H2(j));  polq(1,j)=(SiH4q(j)-2*H2q(j))
              pol(2,j)=SiF4(j);  polq(2,j)=SiF4q(j)
              pol(3,j)=SiCl4(j); polq(3,j)=SiCl4q(j)
              pol(4,j)=(Si2H4(j)-2*H2(j)); polq(4,j)=(Si2H4q(j)-2*H2q(j))
              pol(5,j)=(Si2H6(j)-3*H2(j)); polq(5,j)=(Si2H6q(j)-3*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

   elseif (a.eq.15.0) then ! P atom

       ! tol = 0.2
        
        d1 = abs(Veff3 - PH3v3(1))
        d2 = abs(Veff3 - PCl3v3(1))
        d3 = abs(Veff3 - PF3v3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        
        select case (pick)
        
        case (1) ! PH3 vol(5) te(1)
           ratio_v3=Veff3/(PH3v3(5)-3*PH3v3(2)); s1 = Veff3/PH3v3(1)
           ratio_v5=Veff5/(PH3v5(5)-3*PH3v5(2)); s2 = Veff5/PH3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(PH3(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(PH3q(j)-1.5*H2q(j))
           end do
        
        case (2) ! PCl3
           ratio_v3=Veff3/PCl3v3(5); s1 = Veff3/PCl3v3(1)
           ratio_v5=Veff5/PCl3v5(5); s2 = Veff5/PCl3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*PCl3(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*PCl3q(j)
           end do
        
        case (3) ! PF3
           ratio_v3=Veff3/PF3v3(5); s1 = Veff3/PF3v3(1)
           ratio_v5=Veff5/PF3v5(5); s2 = Veff5/PF3v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*PF3(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*PF3q(j)
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(3), vvol3(3), vte3(3), vvol5(3), vte5(3)
           real :: pol(3,23), polq(3,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3 /)
        
           vvol3 = (/ (PH3v3(5)-3*PH3v3(2)), PCl3v3(5), PF3v3(5) /)
           vte3  = (/ PH3v3(1), PCl3v3(1), PF3v3(1) /)
        
           vvol5 = (/ (PH3v5(5)-3*PH3v5(2)), PCl3v5(5), PF3v5(5) /)
           vte5  = (/ PH3v5(1), PCl3v5(1), PF3v5(1) /)
        
           do j=1,23
              pol(1,j)=(PH3(j)-1.5*H2(j));   polq(1,j)=(PH3q(j)-1.5*H2q(j))
              pol(2,j)=PCl3(j);  polq(2,j)=PCl3q(j)
              pol(3,j)=PF3(j);   polq(3,j)=PF3q(j)
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

   elseif (a.eq.16.0) then ! S atom
  
       ! tol = 0.2
        
        d1 = abs(Veff3 - H2Sv3(1))
        d2 = abs(Veff3 - SO2v3(1))
        d3 = abs(Veff3 - SF2v3(1))
        d4 = abs(Veff3 - SCl2v3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        if (d4<=tol .and. d4<dmin) then; dmin=d4; pick=4; end if
        
        select case (pick)
        
        case (1) ! H2S vol(4) te(2)
           ratio_v3=Veff3/(H2Sv3(4)-2*H2Sv3(2)); s1 = Veff3/H2Sv3(1)
           ratio_v5=Veff5/(H2Sv5(4)-2*H2Sv5(2)); s2 = Veff5/H2Sv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(H2S(j)-H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(H2Sq(j)-H2q(j))
           end do
        
        case (2) ! SO2
           ratio_v3=Veff3/SO2v3(4); s1 = Veff3/SO2v3(1)
           ratio_v5=Veff5/SO2v5(4); s2 = Veff5/SO2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*SO2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*SO2q(j)
           end do
        
        case (3) ! SF2
           ratio_v3=Veff3/SF2v3(4); s1 = Veff3/SF2v3(1)
           ratio_v5=Veff5/SF2v5(4); s2 = Veff5/SF2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*SF2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*SF2q(j)
           end do
        
        case (4) ! SCl2
           ratio_v3=Veff3/SCl2v3(4); s1 = Veff3/SCl2v3(1)
           ratio_v5=Veff5/SCl2v5(4); s2 = Veff5/SCl2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*SCl2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*SCl2q(j)
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(4), vvol3(4), vte3(4), vvol5(4), vte5(4)
           real :: pol(4,23), polq(4,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3,d4 /)
        
           vvol3 = (/ (H2Sv3(4)-2*H2Sv3(2)), SO2v3(4), SF2v3(4), SCl2v3(4) /)
           vte3  = (/ H2Sv3(1), SO2v3(1), SF2v3(1), SCl2v3(1) /)
        
           vvol5 = (/ (H2Sv5(4)-2*H2Sv5(2)), SO2v5(4), SF2v5(4), SCl2v5(4) /)
           vte5  = (/ H2Sv5(1), SO2v5(1), SF2v5(1), SCl2v5(1) /)
        
           do j=1,23
              pol(1,j)=(H2S(j)-H2(j));  polq(1,j)=(H2Sq(j)-H2q(j))
              pol(2,j)=SO2(j);  polq(2,j)=SO2q(j)
              pol(3,j)=SF2(j);  polq(3,j)=SF2q(j)
              pol(4,j)=SCl2(j); polq(4,j)=SCl2q(j)
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

  
   elseif (a.eq.17.0) then ! Cl atom
     if (b.eq.1.0) then  
        ratio_v3 = Veff3/(HClv3(3)-HClv3(2))
        s1       = Veff3/HClv3(1)
        
        ratio_v5 = Veff5/(HClv5(3)-HClv5(2))
        s2       = Veff5/HClv5(1)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             *ratio_v3*(HCl(j)-0.5*H2(j))
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) *ratio_v5*(HClq(j)-0.5*H2q(j))
        end do
        
     elseif (b.eq.5.0) then
 
        ratio_v3 = Veff3/BCl3v3(5)
        s1       = Veff3/BCl3v3(2)
        
        ratio_v5 = Veff5/BCl3v5(5)
        s2       = Veff5/BCl3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * BCl3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * BCl3q(j)
        end do
        
     elseif (b.eq.9.0) then
        
        ratio_v3 = Veff3/ClFv3(3)
        s1       = Veff3/ClFv3(1)
        
        ratio_v5 = Veff5/ClFv5(3)
        s2       = Veff5/ClFv5(1)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * ClF(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * ClFq(j)
        end do
        
     elseif (b.eq.13.0) then
 
        ratio_v3 = Veff3/AlCl3v3(5)
        s1       = Veff3/AlCl3v3(2)
        
        ratio_v5 = Veff5/AlCl3v5(5)
        s2       = Veff5/AlCl3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * AlCl3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * AlCl3q(j)
        end do
        
     elseif (b.eq.14.0) then
 
        ratio_v3 = Veff3/SiCl4v3(6)
        s1       = Veff3/SiCl4v3(2)
        
        ratio_v5 = Veff5/SiCl4v5(6)
        s2       = Veff5/SiCl4v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * SiCl4(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * SiCl4q(j)
        end do
        
     elseif (b.eq.15.0) then
 
        ratio_v3 = Veff3/PCl3v3(5)
        s1       = Veff3/PCl3v3(2)
        
        ratio_v5 = Veff5/PCl3v5(5)
        s2       = Veff5/PCl3v5(2)
        
        do j=1,23
           aiw(i,j)  = abs(1.0 - c11*s1)             * ratio_v3 * PCl3(j)
           aiwq(i,j) = abs(1.0 - c12*s2 - c22*s2**2) * ratio_v5 * PCl3q(j)
        end do
        
     elseif (b.eq.17.0) then
 
        ratio_v3 = Veff3/Cl2v3(3)
        s1       = Veff3/Cl2v3(1)
        
        ratio_v5 = Veff5/Cl2v5(3)
        s2       = Veff5/Cl2v5(1)
        
        do j=1,23
           aiw(i,j)  = (1.0 - c11*s1)             * ratio_v3 * Cl2(j)
           aiwq(i,j) = (1.0 - c12*s2 - c22*s2**2) * ratio_v5 * Cl2q(j)
        end do
     end if


   elseif (a.eq.18.0) then ! Ar atom
  
             Do j=1,23
               aiw(i,j)=Ar(j)
               aiwq(i,j)=Arq(j)
             end Do

   elseif (a.eq.35.0) then ! Br atom
       ! tol = 0.2
        
        d1 = abs(Veff3 - HBrv3(1))
        d2 = abs(Veff3 - Br2v3(1))
        d3 = abs(Veff3 - CF3Brv3(1))
        d4 = abs(Veff3 - CH3Brv3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        if (d4<=tol .and. d4<dmin) then; dmin=d4; pick=4; end if
        
        select case (pick)
        
        case (1) ! HBr vol(3) te(1)
           ratio_v3=Veff3/(HBrv3(3)-HBrv3(2)); s1 = Veff3/HBrv3(1)
           ratio_v5=Veff5/(HBrv5(3)-HBrv5(2)); s2 = Veff5/HBrv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HBr(j)-0.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HBrq(j)-0.5*H2q(j))
           end do
        
        case (2) ! Br2 vol(3) te(1)
           ratio_v3=Veff3/Br2v3(3); s1 = Veff3/Br2v3(1)
           ratio_v5=Veff5/Br2v5(3); s2 = Veff5/Br2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*Br2(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*Br2q(j)
           end do
        
        case (3) ! CF3Br vol(6) te(1)
           ratio_v3=Veff3/CF3Brv3(6); s1 = Veff3/CF3Brv3(1)
           ratio_v5=Veff5/CF3Brv5(6); s2 = Veff5/CF3Brv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*CF3Br(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*CF3Brq(j)
           end do
        
        case (4) ! CH3Br vol(6) te(1)
           ratio_v3=Veff3/(CH3Brv3(6)-3*CH3Brv3(3)); s1 = Veff3/CH3Brv3(1)
           ratio_v5=Veff5/(CH3Brv5(6)-3*CH3Brv5(3)); s2 = Veff5/CH3Brv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(CH3Br(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(CH3Brq(j)-1.5*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(4), vvol3(4), vte3(4), vvol5(4), vte5(4)
           real :: pol(4,23), polq(4,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3,d4 /)
        
           vvol3 = (/ (HBrv3(3)-HBrv3(2)), Br2v3(3), CF3Brv3(6), (CH3Brv3(6)-3*CH3Brv3(3)) /)
           vte3  = (/ HBrv3(1), Br2v3(1), CF3Brv3(1), CH3Brv3(1) /)
        
           vvol5 = (/ (HBrv5(3)-HBrv5(2)), Br2v5(3), CF3Brv5(6), (CH3Brv5(6)-3*CH3Brv5(3)) /)
           vte5  = (/ HBrv5(1), Br2v5(1), CF3Brv5(1), CH3Brv5(1) /)
        
           do j=1,23
              pol(1,j)=(HBr(j)-0.5*H2(j));    polq(1,j)=(HBrq(j)-0.5*H2q(j))
              pol(2,j)=Br2(j);    polq(2,j)=Br2q(j)
              pol(3,j)=CF3Br(j);  polq(3,j)=CF3Brq(j)
              pol(4,j)=(CH3Br(j)-1.5*H2(j));  polq(4,j)=(CH3Brq(j)-1.5*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select
        
   elseif (a.eq.53.0) then ! I atom

      !  tol = 0.2
        
        d1 = abs(Veff3 - HIv3(1))
        d2 = abs(Veff3 - I2v3(1))
        d3 = abs(Veff3 - CF3Iv3(1))
        d4 = abs(Veff3 - CH3Iv3(1))
        
        pick=0; dmin=huge(1.0)
        if (d1<=tol .and. d1<dmin) then; dmin=d1; pick=1; end if
        if (d2<=tol .and. d2<dmin) then; dmin=d2; pick=2; end if
        if (d3<=tol .and. d3<dmin) then; dmin=d3; pick=3; end if
        if (d4<=tol .and. d4<dmin) then; dmin=d4; pick=4; end if
        
        select case (pick)
        
        case (1) ! HI vol(3) te(1)
           ratio_v3=Veff3/(HIv3(3)-HIv3(2)); s1 = Veff3/HIv3(1)
           ratio_v5=Veff5/(HIv5(3)-HIv5(2)); s2 = Veff5/HIv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(HI(j)-0.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(HIq(j)-0.5*H2q(j))
           end do
        
        case (2) ! I2 vol(3) te(1)
           ratio_v3=Veff3/I2v3(3); s1 = Veff3/I2v3(1)
           ratio_v5=Veff5/I2v5(3); s2 = Veff5/I2v5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*I22(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*I22q(j)
           end do
        
        case (3) ! CF3I vol(6) te(1)
           ratio_v3=Veff3/CF3Iv3(6); s1 = Veff3/CF3Iv3(1)
           ratio_v5=Veff5/CF3Iv5(6); s2 = Veff5/CF3Iv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*CF3I(j)
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*CF3Iq(j)
           end do
        
        case (4) ! CH3I vol(6) te(1)
           ratio_v3=Veff3/(CH3Iv3(6)-3*CH3Iv3(3)); s1 = Veff3/CH3Iv3(1)
           ratio_v5=Veff5/(CH3Iv5(6)-3*CH3Iv5(3)); s2 = Veff5/CH3Iv5(1)
           do j=1,23
              aiw(i,j)=abs(1.0 - c11*s1)*ratio_v3*(CH3I(j)-1.5*H2(j))
              aiwq(i,j)=abs(1.0 - c12*s2 - c22*s2**2)*ratio_v5*(CH3Iq(j)-1.5*H2q(j))
           end do
        
        case default
           block
           use polar_interp, only: polar_weighted_avg
           real :: d(4), vvol3(4), vte3(4), vvol5(4), vte5(4)
           real :: pol(4,23), polq(4,23)
           real :: aiw_tmp(23), aiwq_tmp(23)
           integer :: j
        
           d     = (/ d1,d2,d3,d4 /)
        
           vvol3 = (/ (HIv3(3)-HIv3(2)), I2v3(3), CF3Iv3(6), (CH3Iv3(6)-3*CH3Iv3(3)) /)
           vte3  = (/ HIv3(1), I2v3(1), CF3Iv3(1), CH3Iv3(1) /)
        
           vvol5 = (/ (HIv5(3)-HIv5(2)), I2v5(3), CF3Iv5(6), (CH3Iv5(6)-3*CH3Iv5(3)) /)
           vte5  = (/ HIv5(1), I2v5(1), CF3Iv5(1), CH3Iv5(1) /)
        
           do j=1,23
              pol(1,j)=(HI(j)-0.5*H2(j));    polq(1,j)=(HIq(j)-0.5*H2q(j))
              pol(2,j)=I22(j);   polq(2,j)=I22q(j)
              pol(3,j)=CF3I(j);  polq(3,j)=CF3Iq(j)
              pol(4,j)=(CH3I(j)-1.5*H2(j));  polq(4,j)=(CH3Iq(j)-1.5*H2q(j))
           end do
        
           call polar_weighted_avg(Veff3, Veff5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_tmp, aiwq_tmp)
        
           aiw(i,1:23)=aiw_tmp(1:23)
           aiwq(i,1:23)=aiwq_tmp(1:23)
           end block
        
        end select

   end if  
   end Subroutine AIM_polarizabilities



!
!
   Subroutine Coefficients(sita,coef6,coef8,damp,sph,dict,dict8,beta,alph,i,elem,vers)

   parameter (nsmax=100,conv=0.529177209)
   Real*8 :: a,b,c,d,e,f,sita
   Real, Intent(out) :: dict,dict8,beta,alph
   Integer :: i,elem
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax),dicta(nsmax),dicta8(nsmax),betaa(nsmax)
   Real*8, Dimension(99) :: coef6,coef8,damp,sph
   Data dicta/ 0.199730d0,0.625058d0,0.216319d0,1.5290d0,0.117956d0,0.155250d0,0.028773d0,&
         &  0.984668d0,0.957855d0,0.338996d0,0.405991d0,0.09588d0,0.185489d0,&
         &  0.185411d0,0.071774d0,0.003104d0,1.297139d0,0.914964d0,0.461576d0,&
         &  1.427069d0,0.72359d0,0.900101d0,0.261125d0,0.415583d0,0.009843d0,&
         &  0.890484d0,4.793738d0,6.254113d0,3.33151d0,9.995193d0,21.167203d0,69*0/
   Data dicta8/ 0.010225d0,0.004286d0,0.002841d0,0.1649d0,0.012842d0,0.003195d0,0.005962d0,&
           0.060477d0,0.011257d0,0.29648d0,0.011234d0,0.008433d0,0.020926d0,&
           0.045112d0,0.004783d0,0.264031d0,0.171713d0,0.268786d0,0.235496d0,&
           0.070811d0,0.052314d0,0.015629d0,0.028508d0,0.900131d0,0.894086d0,&
           2.415552d0,0.825826d0,0.691238d0,0.509268d0,1.169732d0,2.745814d0,69*0/
   Data betaa/ 1.898977d0,2.022739d0,1.737714d0,1.9509d0,1.456906d0,1.623847d0,1.618463d0,&
             1.828202d0,1.62796d0,0.589819d0,1.416878d0,1.66046d0,1.386093d0,&
             1.310652d0,2.225149d0,1.695239d0,1.951721d0,1.65896d0,1.780635d0,&
             2.396253d0,2.368908d0,2.458106d0,2.13224d0,1.133178d0,2.726511d0,&
             1.634346d0,2.279597d0,1.643112d0,1.843843d0,1.572187d0,1.463387d0,69*0/
   a=sita(4,i)
   b=sita(5,i)
   c=sita(6,i)
   d=sita(7,i)
   e=sita(8,i)
   f=sita(9,i)
!   print*,'vers',vers
!  version 3 corresponds to Das_2020
   if (vers==3.0) then
      if (a.eq.1.0.and.b.eq.1.0) then !opt 80
         dict=coef6(80)
         dict8=coef8(80)
         beta=0.32*damp(80)
         alph=sph(80)
         elem=1
      elseif (a.eq.1.0.and.b.eq.5.0) then !opt 81
         dict=coef6(81)
         dict8=coef8(81)
         beta=damp(81)
         alph=sph(81)
         elem=2
      elseif (a.eq.1.0.and.b.eq.6.0) then !opt 82
         dict=coef6(82) !0.2163190 !coef6(82)
         dict8=coef8(82) !0.002841 !coef8(82)
         beta=damp(82) !1.737714 !damp(82)
         alph=sph(82)
         elem=3
!      elseif (a.eq.6.0.and.b.eq.1.0) then
!         dict=dicta(4)
!         dict8=dicta8(4)
!         beta=betaa(4)
      elseif (a.eq.1.0.and.b.eq.7.0) then !opt 83
         dict=coef6(83)
         dict8=coef8(83)
         beta=damp(83)
         alph=sph(83)
         elem=4
      elseif (a.eq.1.0.and.b.eq.8.0) then !opt 84
         dict=coef6(84)
         dict8=coef8(84)
         beta=damp(84)
         alph=sph(84)
         elem=5
      elseif (a.eq.1.0.and.b.eq.9.0) then !opt 85
         dict=coef6(85)
         dict8=coef8(85)
         beta=damp(85)
         alph=sph(85)
         elem=6
        ! print*,'H-F',beta
      elseif (a.eq.1.0.and.b.eq.13.0) then !opt 86
         dict=coef6(86)
         dict8=coef8(86)
         beta=damp(86)
         alph=sph(86)
         elem=7
      elseif (a.eq.1.0.and.b.eq.14.0) then !opt 87
         dict=coef6(87)
         dict8=coef8(87)
         beta=0.32*damp(87)
         alph=sph(87)
         elem=8
      elseif (a.eq.1.0.and.b.eq.15.0) then !opt 88
         dict=coef6(88)
         dict8=coef8(88)
         beta=0.32*damp(88)
         alph=sph(88)
         elem=9
      elseif (a.eq.1.0.and.b.eq.16.0) then !opt 89
         dict=coef6(89)
         dict8=coef8(89)
         beta=damp(89)
         alph=sph(89)
         elem=10
      elseif (a.eq.1.0.and.b.eq.17.0) then !opt 90
         dict=coef6(90)
         dict8=coef8(90)
         beta=damp(90)
         alph=sph(90)
         elem=11
      elseif (a.eq.1.0.and.b.eq.35.0) then !opt 91
         dict=coef6(91)
         dict8=coef8(91)
         beta=damp(91)
         alph=sph(91)
         elem=12
      elseif (a.eq.1.0.and.b.eq.53.0) then
         dict=dicta(14)
         dict8=dicta8(14)
         beta=betaa(14)
      elseif (a.eq.2.0.and.b.eq.0.0) then !opt 96
         dict=coef6(96)
         dict8=coef8(96)
         beta=damp(96)
         alph=sph(96)
         elem=13
!      elseif (a.eq.5.0.and.b.eq.0.0) then ! Boron atom
!         dict=dicta(16)
!         dict8=dicta8(16)
!         beta=betaa(16)
      elseif (a.eq.5.0) then

        if (a.eq.5.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.0.0.and.f.eq.0.0) then ! HHHB type !opt 1
         dict=coef6(1)
         dict8=coef8(1)
         beta=damp(1)
         alph=sph(1)
         elem=14
        elseif (a.eq.5.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.0.0.and.f.eq.0.0) then ! FFFB type !opt 2
         dict=coef6(2)
         dict8=coef8(2)
         beta=0.32*damp(2)
         alph=sph(2)
         elem=14
        elseif (a.eq.5.0.and.b.eq.17.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then ! BClClCl !opt 3
         dict=coef6(3)
         dict8=coef8(3)
         beta=0.32*damp(3)
         alph=sph(3)
         elem=14
        else
         dict=(coef6(1)+coef6(2)+coef6(3))/3.0
         dict8=(coef8(1)+coef8(2)+coef8(3))/3.0
         beta=(damp(1)+0.32*damp(2)+0.32*damp(3))/3.0
         elem=14
        end if
!      elseif (a.eq.6.0.and.b.eq.2.0) then  ! Carbon atom
!         dict=dicta(17)
!         dict8=dicta8(17)
!         beta=betaa(17)
!      elseif (a.eq.6.0.and.b.eq.3.0) then
!         dict=dicta(18)
!         dict8=dicta8(18)
!         beta=betaa(18)
!      elseif (a.eq.6.0.and.b.eq.4.0) then
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
      elseif (a.eq.6.0.and.b.eq.7.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! N
         dict=dicta(4)
         dict8=dicta8(4)
         beta=betaa(4)
      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! H
         dict=dicta(4)
         dict8=dicta8(4)
         beta=betaa(4)
      elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! C
         dict=dicta(4)
         dict8=dicta8(4)
         beta=betaa(4)
      elseif (a.eq.6.0.and.b.eq.8.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! O
         dict=dicta(4)
         dict8=dicta8(4)
         beta=betaa(4)
      elseif (a.eq.6.0.and.b.ne.0.0.and.c.ne.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then
        if (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! C2H2 type !!!!opt 78
         dict=coef6(78)
         dict8=coef8(78)
         beta=damp(78)
         alph=sph(78)
         elem=15
!!
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.7.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HN !opt 4
         dict=coef6(4)     !1.297139 !coef6(4)
         dict8=coef8(4)    !0.1717130 !coef8(4)
         beta=damp(4)      !1.951721 !damp(4)
         alph=sph(4)
         elem=15
!         dict=(coef6(4)+coef6(5)+coef6(6)+coef6(7)+coef6(78))/5.0 !coef6(4)!1.297139 !coef6(4)
!         dict8=(coef8(4)+coef8(5)+coef8(6)+coef8(7)+coef8(78))/5.0 !coef8(4)!0.1717130 !coef8(4)
!         beta=(damp(4)+damp(5)+damp(6)+damp(7)+damp(78))/15.0 !damp(4)!1.951721 !damp(4)
!         alph=sph(4)
!         elem=15
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.9.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CF !opt 5
         dict=coef6(5)
         dict8=coef8(5)
         beta=damp(5)
         alph=sph(5)
         elem=15
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.7.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CN 
         dict=coef6(5)
         dict8=coef8(5)
         beta=damp(5)
         alph=sph(5)
         elem=15
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.17.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CCl !opt 6
         dict=coef6(6)
         dict8=coef8(6)
         beta=damp(6)
         alph=sph(6)
         elem=15
        elseif (a.eq.6.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !OO !opt 7
         dict=coef6(7)
         dict8=coef8(7)
         beta=0.35*damp(7)
         alph=sph(7)
         elem=15
!        elseif (a.eq.6.0.and.b.eq.8.0.and.c.eq.16.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !OS
!         dict=coef6(6) !dicta(17)
!         dict8=coef8(6) !dicta8(17)
!         beta=damp(6) !betaa(17)
!         alph=sph(6)
!         elem=15
!        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then! CC
!         dict=dicta(17)
!         dict8=dicta8(17)
!         beta=betaa(17)
!        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.15.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then! CP
!         dict=dicta(17)
!         dict8=dicta8(17)
!         beta=betaa(17)
        else !(a.eq.6.0.and.b.eq.1.0.and.c.eq.7.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HN !opt 4
         !dict=coef6(4)     !1.297139 !coef6(4)
         !dict8=coef8(4)    !0.1717130 !coef8(4)
         !beta=damp(4)      !1.951721 !damp(4)
         !alph=sph(4)
         !elem=15
         dict=(coef6(4)+coef6(5)+coef6(6)+coef6(7)+coef6(78))/5.0 !coef6(4)!1.297139 !coef6(4)
         dict8=(coef8(4)+coef8(5)+coef8(6)+coef8(7)+coef8(78))/5.0 !coef8(4)!0.1717130 !coef8(4)
         beta=(damp(4)+damp(5)+damp(6)+0.35*damp(7)+damp(78))/5.0 !damp(4)!1.951721 !damp(4)
         alph=sph(4)
         elem=15
        end if
!!
      elseif (a.eq.6.0.and.b.ne.0.0.and.c.ne.0.0.and.d.ne.0.0.and.e.eq.0.0.and.f.eq.0.0) then
        if (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !C2H4 type !!!!opt 79
         dict=coef6(79)
         dict8=coef8(79)
         beta=damp(79)
         alph=sph(79)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !HCC !opt 8
         dict=0.914964  !coef6(8)
         dict8=0.268786 !coef8(8)
         beta=damp(8) !1.65896  !damp(8)
         elem=16
        ! print*,'HCC',0.914964*10884.3258684610*(1.0/627.509474),0.268786*3886863.35457252*(1.0/627.509474)
!        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.16.0.and.e.eq.0.0.and.f.eq.0.0) then!CCS
!         dict=dicta(18)
!         dict8=dicta8(18)
!         beta=betaa(18)
!        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then!HCO
!         dict=dicta(18)
!         dict8=dicta8(18)
!         beta=betaa(18)
!        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.16.0.and.e.eq.0.0.and.f.eq.0.0) then!HCS
!         dict=dicta(17)
!         dict8=dicta8(17)
!         beta=betaa(17)
!        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.16.0.and.e.eq.0.0.and.f.eq.0.0) then!HHS
!         dict=dicta(18)
!         dict8=dicta8(18)
!         beta=betaa(18)
!!!
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.7.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !CNO opt 9
         dict=coef6(9)
         dict8=coef8(9)
         beta=damp(9)
         alph=sph(9)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then !CCCl opt 10
         dict=coef6(10)
         dict8=coef8(10)
         beta=0.5*damp(10)
         alph=sph(10)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then !CClCl opt 11
         dict=coef6(11)
         dict8=coef8(11)
         beta=damp(11)
         alph=sph(11)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.9.0.and.e.eq.0.0.and.f.eq.0.0) then !CCF opt 12
         dict=coef6(12)
         dict8=coef8(12)
         beta=damp(12)
         alph=sph(12)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.0.0.and.f.eq.0.0) then !CFF opt 13
         dict=coef6(13)
         dict8=coef8(13)
         beta=0.35*damp(13)
         alph=sph(13)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.35.0.and.e.eq.0.0.and.f.eq.0.0) then !CCBr opt 14
         dict=coef6(14)
         dict8=coef8(14)
         beta=0.35*damp(14)
         alph=sph(14)
         elem=16
        elseif (a.eq.6.0.and.b.eq.7.0.and.c.eq.7.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then !NNN opt 15
         dict=coef6(15)
         dict8=coef8(15)
         beta=0.35*damp(15)
         alph=sph(15)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then !HCN opt 16
         dict=coef6(16)
         dict8=coef8(16)
         beta=damp(16)
         alph=sph(16)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.7.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then !HNN opt 17
         dict=coef6(17)
         dict8=coef8(17)
         beta=damp(17)
         alph=sph(17)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !CCC opt 18
         dict=coef6(18) !0.914964  ! coef6(18)
         dict8=coef8(18) !0.268786 !coef8(18)
         beta=0.35*damp(18) !1.65896  !damp(18)
         alph=sph(18)
         elem=16
       !  print*,'CCC',dict*10884.3258684610*(1.0/627.509474),dict8*3886863.35457252*(1.0/627.509474)
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.7.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then !CNN opt 19
         dict=coef6(19)
         dict8=coef8(19)
         beta=damp(19)
         alph=sph(19)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then !CCN opt 20
         dict=coef6(20)
         dict8=coef8(20)
         beta=damp(20)
         alph=sph(20)
         elem=26
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !COO opt 21
         dict=coef6(21)
         dict8=coef8(21)
         beta=damp(21)
         alph=sph(21)
         elem=16
        elseif (a.eq.6.0.and.b.eq.7.0.and.c.eq.7.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !NNO opt 22
         dict=coef6(22)
         dict8=coef8(22)
         beta=damp(22)
         alph=sph(22)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !HHO opt 23
         dict=coef6(23)
         dict8=coef8(23)
         beta=0.5*damp(23)
         alph=sph(23)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !HOO opt 24
         dict=coef6(24)
         dict8=coef8(24)
         beta=damp(24)
         alph=sph(24)
         elem=16
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.7.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !HNO opt 25
         dict=coef6(25)
         dict8=coef8(25)
         beta=damp(25)
         alph=sph(25)
         elem=16
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !CCO opt 26
         dict=coef6(26)
         dict8=coef8(26)
         beta=damp(26)
         alph=sph(26)
         elem=16
        else
         dict=(coef6(9)+coef6(10)+coef6(11)+coef6(12)+coef6(13)+coef6(14)+coef6(15)+coef6(16)+&
     & coef6(17)+coef6(18)+coef6(19)+coef6(20)+coef6(21)+coef6(22)+coef6(23)+coef6(24)+coef6(25)+coef6(26))/18.0!coef6(31)
         dict8=(coef8(9)+coef8(10)+coef8(11)+coef8(12)+coef8(13)+coef8(14)+coef8(15)+coef8(16)+&
     & coef8(17)+coef8(18)+coef8(19)+coef8(20)+coef8(21)+coef8(22)+coef8(23)+coef8(24)+coef8(25)+coef8(26))/18.0!coef8(31)
         beta=(damp(9)+0.5*damp(10)+damp(11)+damp(12)+0.35*damp(13)+0.35*damp(14)+ &
     & 0.35*damp(15)+damp(16)+damp(17)+0.35*damp(18)+damp(19)+damp(20)+damp(21)+damp(22)+0.5*damp(23)+damp(24)+damp(25)+damp(26))/18.0!damp(31)
         elem=16
        end if
!!!
      elseif (a.eq.6.0.and.b.ne.0.0.and.c.ne.0.0.and.d.ne.0.0.and.e.ne.0.0.and.f.eq.0.0) then
        if (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.6.0.and.f.eq.0.0) then !HHCC opt 27
         dict=coef6(27)
         dict8=coef8(27)
         beta=damp(27)
         alph=sph(27)
         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.7.0.and.f.eq.0.0) then !HHCN
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.8.0.and.f.eq.0.0) then !HHCO
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.9.0.and.f.eq.0.0) then !CFFF
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.7.0.and.f.eq.0.0) then !CFFN
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.8.0.and.f.eq.0.0) then !HCCO
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.6.0.and.f.eq.0.0) then !HCCC
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
!      elseif (a.eq.6.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.9.0.and.f.eq.0.0) then !FFFF
!         dict=0.08*coef6(29) !dicta(19)
!         dict8=0.08*coef8(29) !coef8(29) !dicta8(19)
!         beta=0.9*damp(29) !damp(29) !betaa(19)
!         alph=sph(29)
!         elem=29
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.7.0.and.f.eq.0.0) then!HCCN
!         dict=0.01*coef6(29) !dicta(19)
!         dict8=0.01*coef8(29) !dicta8(19)
!         beta=damp(29) !betaa(19)
!         alph=sph(29)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.16.0.and.f.eq.0.0) then !HHCS
!         dict=dicta(19)
!         dict8=dicta8(19)
!         beta=betaa(19)
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.8.0.and.f.eq.0.0) then !HHHO opt 28
         dict=coef6(28)
         dict8=coef8(28)
         beta=0.5*damp(28)
         alph=sph(28)
         elem=17
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.9.0.and.f.eq.0.0) then !HHHF opt 29
         dict=coef6(29)
         dict8=coef8(29)
         beta=0.35*damp(29)
         alph=sph(29)
         elem=17
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.17.0.and.f.eq.0.0) then !HHHCl opt 30
         dict=coef6(30)
         dict8=coef8(30)
         beta=0.5*damp(30)
         alph=sph(30)
         elem=17
!      elseif (a.eq.6.0.and.b.eq.8.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.9.0.and.f.eq.0.0) then !OFFF opt 31
         !dict=coef6(31)
         !dict8=coef8(31)
         !beta=damp(31)
         !alph=sph(31)
         !elem=17
!         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
!     & coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(93)+coef6(94))/15.0!coef6(31)
!         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
!     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+ coef8(93)+coef8(94))/15.0!coef8(31)
!         beta=(damp(27)+damp(28) + damp(29)+damp(30) + damp(31)+damp(32) + &
!     & damp(33)+damp(34)+damp(35)+damp(36)+damp(37)+damp(38)+damp(39)+damp(93)+damp(94))/15.0!damp(31)
!         alph=sph(31)
!         elem=17
!       elseif (a.eq.6.0.and.b.eq.8.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.17.0.and.f.eq.0.0) then !OClClCl opt 32
         !dict=coef6(32)
         !dict8=coef8(32)
         !beta=damp(32)
         !alph=sph(32)
         !elem=17
!         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
!     &coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(93)+coef6(94))/15.0!coef6(31)
!         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
!     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+ coef8(93)+coef8(94))/15.0!coef8(31)
!         beta=(damp(27)+damp(28) + damp(29)+damp(30) + damp(31)+damp(32) + &
!     &damp(33)+damp(34)+damp(35)+damp(36)+damp(37)+damp(38)+damp(39)+damp(93)+damp(94))/15.0!damp(31)
!         alph=sph(31)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.17.0.and.f.eq.0.0) then !FFFCl opt 33
         !dict=coef6(33)
         !dict8=coef8(33)
         !beta=damp(33)
         !alph=sph(33)
         !elem=17
!         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
!     &coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(93)+coef6(94))/15.0!coef6(31)
!         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
!     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+ coef8(93)+coef8(94))/15.0!coef8(31)
!         beta=(damp(27)+damp(28) + damp(29)+damp(30) + damp(31)+damp(32) + &
!     &damp(33)+damp(34)+damp(35)+damp(36)+damp(37)+damp(38)+damp(39)+damp(93)+damp(94))/15.0!damp(31)
!         alph=sph(31)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.9.0.and.f.eq.0.0) then !HFFF opt 34
         !dict=coef6(34)
         !dict8=coef8(34)
         !beta=damp(34)
         !alph=sph(34)
         !elem=17
!         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
!     & coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(93)+coef6(94))/15.0!coef6(31)
!         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
!     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+ coef8(93)+coef8(94))/15.0!coef8(31)
!         beta=(damp(27)+damp(28) + damp(29)+damp(30) + damp(31)+damp(32) + &
!     & damp(33)+damp(34)+damp(35)+damp(36)+damp(37)+damp(38)+damp(39)+damp(93)+damp(94))/15.0!damp(31)
!         alph=sph(31)
!         elem=17
        elseif (a.eq.6.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.35.0.and.f.eq.0.0) then !FFFBr opt 35
         dict=coef6(35)
         dict8=coef8(35)
         beta=damp(35)
         alph=sph(35)
         elem=17
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.35.0.and.f.eq.0.0) then !HHHBr opt 36
         dict=coef6(36)
         dict8=coef8(36)
         beta=0.4*damp(36)
         alph=sph(36)
         elem=17
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.6.0.and.f.eq.0.0) then ! HHHC type !opt 37
         dict=coef6(37)
         dict8=coef8(37)
         beta=damp(37)
         alph=sph(37)
         elem=17
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.7.0.and.f.eq.0.0) then ! HHHN type !opt 38
         dict=coef6(38)
         dict8=coef8(38)
         beta=damp(38)
         alph=sph(38)
         elem=17
         !print*,'C-HHHN',beta
        elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.1.0.and.f.eq.0.0) then ! CH4 type !opt 39
         dict=coef6(39)
         dict8=coef8(39)
         beta=0.4*damp(39)
         alph=sph(39)
         elem=17
        elseif (a.eq.6.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.6.0.and.f.eq.0.0) then ! CCCC type !opt 92
         dict=coef6(92)
         dict8=coef8(92)
         beta=0.35*damp(92)
         elem=17
         alph=sph(92)
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.17.0.and.f.eq.0.0) then ! HClClCl type !opt 93
         !dict=coef6(93)
         !dict8=coef8(93)
         !beta=damp(93)
         !alph=sph(93)
         !elem=17
!         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
!     & coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(93)+coef6(94))/15.0!coef6(31)
!         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
!     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+ coef8(93)+coef8(94))/15.0!coef8(31)
!         beta=(damp(27)+damp(28) + damp(29)+damp(30) + damp(31)+damp(32) + &
!     & damp(33)+damp(34)+damp(35)+damp(36)+damp(37)+damp(38)+damp(39)+damp(93)+damp(94))/15.0!damp(31)
!         alph=sph(31)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.16.0.and.f.eq.0.0) then ! HHHS type !opt 94
         !dict=coef6(94)
         !dict8=coef8(94)
         !beta=damp(94)
         !alph=sph(94)
         !elem=17
        else
         dict=(coef6(27)+coef6(28)+coef6(29)+coef6(30)+coef6(31)+coef6(32)+coef6(33)+coef6(34)+&
     & coef6(35)+coef6(36)+coef6(37)+coef6(38)+coef6(39)+coef6(92)+coef6(93)+coef6(94))/16.0!coef6(31)
         dict8=(coef8(27)+coef8(28) + coef8(29)+coef8(30) + coef8(31)+coef8(32)&
     & + coef8(33)+coef8(34) + coef8(35)+coef8(36) + coef8(37)+coef8(38) +coef8(39)+coef8(92)+coef8(93)+coef8(94))/16.0!coef8(31)
         beta=(damp(27)+0.5*damp(28) + 0.35*damp(29)+0.5*damp(30) + damp(31)+damp(32) + &
     & damp(33)+damp(34)+damp(35)+0.4*damp(36)+damp(37)+damp(38)+0.4*damp(39)+0.4*damp(92)+0.35*damp(93)+0.35*damp(94))/16.0!damp(31)
         alph=sph(31)
         elem=17
        end if
!      elseif (a.eq.6.0.and.b.eq.17.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.17.0.and.f.eq.0.0) then ! ClClClCl
!         dict=coef6(93)
!         dict8=coef8(93)
!         beta=damp(93)
!         alph=sph(93)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.15.0.and.f.eq.0.0) then ! HHHP
!         dict=coef6(94)
!         dict8=coef8(94)
!         beta=damp(94)
!         alph=sph(94)
!         elem=17
!      elseif (a.eq.6.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.17.0.and.e.eq.17.0.and.f.eq.0.0) then ! HHClCl
!         dict=coef6(93)
!         dict8=coef8(93)
!         beta=damp(93)
!         alph=sph(93)
!      elseif (a.eq.7.0.and.b.eq.0.0) then ! N atom
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
      elseif (a.eq.7.0) then
        if (a.eq.7.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.1.0.and.f.eq.0.0) then ! NH4
         dict=coef6(99) !0.08*dicta(20)
         dict8=coef8(99) !0.08*dicta8(20)
         beta=damp(99) !betaa(20)
         elem=18
        elseif (a.eq.7.0.and.b.eq.7.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! N2
         dict=dicta(20)
         dict8=dicta8(20)
         beta=betaa(20)
         elem=18
!      elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.8.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! HO
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.7.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! CN
!         dict=1.1*dicta(20)
!         dict8=1.1*dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.7.0.and.c.eq.7.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then! NN
!         dict=0.9*dicta(20)
!         dict8=0.9*dicta8(20)
!         beta=betaa(20)
        elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.0.0.and.f.eq.0.0) then !HHH !opt 40
         dict=coef6(40)
         dict8=coef8(40)
         beta=1.2*damp(40)
         alph=sph(40)
         elem=18
        elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C opt 41
         dict=coef6(41) !0.65*1.427069 !0.65*1.427069 !coef6(41) ! 0.65 fact. added for TCCP_28 etc
         dict8=coef8(41) !0.65*0.070811003 !0.65*0.070811003 !coef8(41)
         beta=damp(41) !0.55*damp(41) !2.396253 !2.396253 !damp(41)
         alph=sph(41)
         elem=18
        elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !HHC opt 42
         dict=coef6(42)
         dict8=coef8(42)
         beta=damp(42)
         alph=sph(42)
         elem=18
        elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HC
         dict=coef6(42)
         dict8=coef8(42)
         beta=damp(42)
         alph=sph(42)
        elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !HCC opt 43
         dict=coef6(43)
         dict8=coef8(43)
         beta=1.3*damp(43) !1.5*damp(43)
         alph=sph(43)
         elem=18
         !print*,'N-HCC',beta
        elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !CCC !opt 44
         dict=coef6(44)
         dict8=coef8(44)
         beta=damp(44)
         alph=sph(44)
         elem=18
         !print*,'N-CCC',beta
        elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CC !opt 45
         dict=coef6(45)
         dict8=coef8(45)
         beta=damp(45)
         alph=sph(45)
         elem=18
        elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !COO
         dict=coef6(44)
         dict8=coef8(44)
         beta=damp(44)
         alph=sph(44)
         elem=18
        elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.6.0.and.f.eq.0.0) then ! HHHC
         dict=coef6(42) !dicta(20)
         dict8=coef8(42) !dicta8(20)
         beta=damp(42) !betaa(20)
         alph=sph(42)
         elem=18
        else
         dict=(coef6(40)+coef6(41)+coef6(42)+coef6(43)+coef6(44)+coef6(45))/6.0
         dict8=(coef8(40)+coef8(41)+coef8(42)+coef8(43)+coef8(44)+coef8(45))/6.0
         beta=(1.2*damp(40)+damp(41)+damp(42)+1.3*damp(43)+damp(44)+damp(45))/6.0
         !print*,beta
         elem=18
        end if
!      elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then ! CCCl
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then ! CCO
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.6.0.and.c.eq.16.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then ! CSCl
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.16.0.and.e.eq.0.0.and.f.eq.0.0) then ! HCS
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then ! HCN
!         dict=1.5*dicta(20)
!         dict8=1.5*dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.7.0.and.d.eq.7.0.and.e.eq.0.0.and.f.eq.0.0) then ! HNN
!         dict=dicta(20)
!         dict8=dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then ! OOO
!         dict=2.25*dicta(20)
!         dict8=2.25*dicta8(20)
!         beta=betaa(20)
!      elseif (a.eq.7.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !opt 19
!         dict=coef6(16)
!         dict8=coef8(16)
!         beta=damp(16)
!      elseif (a.eq.8.0.and.b.eq.0.0) then ! Oxygen atom
!         dict=dicta(21)
!         dict8=dicta8(21)
!         beta=betaa(21)
      elseif (a.eq.8.0) then
        if (a.eq.8.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HH !opt 46
         dict=coef6(46)
         dict8=coef8(46)
         beta=damp(46)
         alph=sph(46)
         elem=19
        elseif (a.eq.8.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C !opt 47
         dict=coef6(47)
         dict8=coef8(47)
         beta=damp(47)
         alph=sph(47)
         elem=19
        elseif (a.eq.8.0.and.b.eq.7.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !N
         dict=coef6(47)
         dict8=coef8(47)
         beta=damp(47)
         alph=sph(47)
         elem=19
!        elseif (a.eq.8.0.and.b.eq.16.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !S
!         dict=dicta(21) !coef6(47) ! dicta(21)
!         dict8=dicta8(21) !coef8(47) ! dicta8(21)
!         beta=betaa(21) !damp(47) ! betaa(21)
!      elseif (a.eq.8.0.and.b.eq.15.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !P
!         dict=dicta(21) !coef6(47) ! dicta(21)
!         dict8=dicta8(21) !coef8(47) ! dicta8(21)
!         beta=betaa(21) !damp(47) ! betaa(21)
        elseif (a.eq.8.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HC !opt 48
         dict=coef6(48)
         dict8=coef8(48)
         beta=damp(48)
         alph=sph(48)
         elem=19
        else
         dict=(coef6(46)+coef6(47)+coef6(48))/3.0
         dict8=(coef8(46)+coef8(47)+coef8(48))/3.0
         beta=(damp(46)+damp(47)+damp(48))/3.0
         elem=19
        end if
!      elseif (a.eq.8.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CC
!         dict=dicta(21)
!         dict8=dicta8(21)
!         beta=betaa(21)
!      elseif (a.eq.8.0.and.b.eq.1.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !H
!         dict=dicta(21)
!         dict8=dicta8(21)
!         beta=betaa(21)
!      elseif (a.eq.8.0.and.b.eq.17.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! Cl
!         dict=dicta(21)
!         dict8=dicta8(21)
!         beta=betaa(21)
!      elseif (a.eq.8.0.and.b.eq.1.0.and.c.eq.15.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! HP
!         dict=0.5*dicta(21)
!         dict8=0.5*dicta8(21)
!         beta=betaa(21)
!      elseif (a.eq.9.0.and.b.eq.0.0) then ! Florine
!         dict=dicta(22)
!         dict8=dicta8(22)
!         beta=betaa(22)
      elseif (a.eq.9.0) then
        if (a.eq.9.0.and.b.eq.1.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !H !opt 49
         dict=coef6(49)
         dict8=coef8(49)
         beta=damp(49)
         alph=sph(49)
         elem=20
         !print*,'F-H',beta
        elseif (a.eq.9.0.and.b.eq.5.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !B !opt 50
         dict=coef6(50)
         dict8=coef8(50)
         beta=0.5*damp(50)
         alph=sph(50)
         elem=20
        elseif (a.eq.9.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C !opt 51
         dict=coef6(51)   !1.2*dicta(22) !coef6(51)
         dict8=coef8(51)  !1.2*dicta8(22) !coef8(51)
         beta=2.0*damp(51)   !betaa(22) !damp(51)
         alph=sph(51)
         elem=51
        elseif (a.eq.9.0.and.b.eq.9.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !F !opt 52
         dict=coef6(52)
         dict8=coef8(52)
         beta=damp(52)
         alph=sph(52)
         elem=20
        elseif (a.eq.9.0.and.b.eq.13.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Al !opt 53
         dict=coef6(53)
         dict8=coef8(53)
         beta=damp(53)
         alph=sph(53)
         elem=20
        elseif (a.eq.9.0.and.b.eq.14.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Si !opt 54
         dict=coef6(54)
         dict8=coef8(54)
         beta=0.5*damp(54)
         alph=sph(54)
         elem=20
        elseif (a.eq.9.0.and.b.eq.15.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !P !opt 55
         dict=coef6(55)
         dict8=coef8(55)
         beta=0.5*damp(55)
         alph=sph(55)
         elem=20
        elseif (a.eq.9.0.and.b.eq.17.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Cl !opt 56
         dict=coef6(56)
         dict8=coef8(56)
         beta=damp(56)
         alph=sph(56)
         elem=20
        else
         dict=(coef6(49)+coef6(50)+coef6(51)+coef6(52)+coef6(53)+coef6(54)+coef6(55)+coef6(56))/8.0
         dict8=(coef8(49)+coef8(50)+coef8(51)+coef8(52)+coef8(53)+coef8(54)+coef8(55)+coef8(56))/8.0
         beta=(damp(49)+0.5*damp(50)+damp(51)+damp(52)+damp(53)+0.5*damp(54)+0.5*damp(55)+damp(56))/8.0
         elem=20
        end if
      elseif (a.eq.10.0.and.b.eq.0.0) then !opt 97
         dict=coef6(97)
         dict8=coef8(97)
         beta=damp(97)
         alph=sph(97)
         elem=21
!      elseif (a.eq.13.0.and.b.eq.0.0) then ! Aluminium atom
!         dict=dicta(24)
!         dict8=dicta8(24)
!         beta=betaa(24)
      elseif (a.eq.13.0) then
        if (a.eq.13.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.0.0.and.f.eq.0.0) then !HHH !opt 57
         dict=coef6(57)
         dict8=coef8(57)
         beta=damp(57)
         alph=sph(57)
         elem=22
        elseif (a.eq.13.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.0.0.and.f.eq.0.0) then !FFF !opt 58
         dict=coef6(58)
         dict8=coef8(58)
         beta=damp(58)
         alph=sph(58)
         elem=22
        elseif (a.eq.13.0.and.b.eq.17.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then !ClClCl !opt 59
         dict=coef6(59)
         dict8=coef8(59)
         beta=damp(59)
         alph=sph(59)
         elem=22
        else
         dict=(coef6(57)+coef6(58)+coef6(59))/3.0
         dict8=(coef8(57)+coef8(58)+coef8(59))/3.0
         beta=(damp(57)+damp(58)+damp(59))/3.0
         elem=22
        end if
!      elseif (a.eq.14.0.and.b.eq.0.0) then ! Silicon atom
!         dict=dicta(25)
!         dict8=dicta8(25)
!         beta=betaa(25)
      elseif (a.eq.14.0) then
        if (a.eq.14.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.1.0.and.f.eq.0.0) then !HHHH !opt 60
         dict=coef6(60)
         dict8=coef8(60)
         beta=0.35*damp(60)
         alph=sph(60)
         elem=23
        elseif (a.eq.14.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.9.0.and.f.eq.0.0) then !FFFF !opt 61
         dict=coef6(61)
         dict8=coef8(61)
         beta=0.35*damp(61)
         alph=sph(61)
         elem=23
        elseif (a.eq.14.0.and.b.eq.17.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.17.0.and.f.eq.0.0) then !ClClClCl !opt 62
         dict=coef6(62)
         dict8=coef8(62)
         beta=0.35*damp(62)
         alph=sph(62)
         elem=23
        else
         dict=(coef6(60)+coef6(61)+coef6(62))/3.0
         dict8=(coef8(60)+coef8(61)+coef8(62))/3.0
         beta=0.35*(damp(60)+damp(61)+damp(62))/3.0
         elem=23
        end if
!      elseif (a.eq.15.0.and.b.eq.0.0) then ! Phosphorus atom
!         dict=dicta(26)
!         dict8=dicta8(26)
!         beta=betaa(26)
      elseif (a.eq.15.0) then
        if (a.eq.15.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.0.0.and.f.eq.0.0) then !HHH !opt 63
         dict=coef6(63)
         dict8=coef8(63)
         beta=0.34*damp(63)
         alph=sph(63)
         elem=24
        elseif (a.eq.15.0.and.b.eq.9.0.and.c.eq.9.0.and.d.eq.9.0.and.e.eq.0.0.and.f.eq.0.0) then !FFF !opt 64
         dict=coef6(64)
         dict8=coef8(64)
         beta=0.34*damp(64)
         alph=sph(64)
         elem=24
        elseif (a.eq.15.0.and.b.eq.17.0.and.c.eq.17.0.and.d.eq.17.0.and.e.eq.0.0.and.f.eq.0.0) then !ClClCl !opt 65
         dict=coef6(65)
         dict8=coef8(65)
         beta=0.34*damp(65)
         alph=sph(65)
         elem=24
        else
         dict=(coef6(63)+coef6(64)+coef6(65))/3.0
         dict8=(coef8(63)+coef8(64)+coef8(65))/3.0
         beta=0.34*(damp(63)+damp(64)+damp(65))/3.0
         !print*,dict,dict8,beta
         elem=24
        end if
!      elseif (a.eq.15.0.and.b.eq.15.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! P
!         dict=1.4*dicta(26) !coef6(65)
!         dict8=1.4*dicta8(26) !coef8(65)
!         beta=betaa(26) !damp(65)
!      elseif (a.eq.15.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then ! C
!         dict=1.8*dicta(26) !coef6(65)
!         dict8=1.8*dicta8(26) !coef8(65)
!         beta=betaa(26) !damp(65)
!      elseif (a.eq.15.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !OOO
!         dict=coef6(65)
!         dict8=coef8(65)
!         beta=damp(65)
!         alph=sph(65)
!      elseif (a.eq.15.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.8.0.and.f.eq.0.0) then !OOOO
!         dict=0.5*coef6(65)
!         dict8=0.5*coef8(65)
!         beta=damp(65)
!         alph=sph(65)
!      elseif (a.eq.15.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.1.0.and.e.eq.8.0.and.f.eq.0.0) then !HHHO
!         dict=1.3*coef6(65)
!         dict8=1.3*coef8(65)
!         beta=damp(65)
!         alph=sph(65)
!      elseif (a.eq.15.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.8.0.and.f.eq.0.0) then !HHHO
!         dict= 1.3*coef6(65)
!         dict8=1.3*coef8(65)
!         beta=damp(65)
!         alph=sph(65)
!      elseif (a.eq.15.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.6.0.and.e.eq.0.0.and.f.eq.0.0) then !CCC
!         dict=0.5*coef6(65)
!         dict8=0.5*coef8(65)
!         beta=damp(65)
!         alph=sph(65)

!      elseif (a.eq.16.0.and.b.eq.0.0) then ! Sulfur atom
!         dict=dicta(27)
!         dict8=dicta8(27)
!         beta=betaa(27)
      elseif (a.eq.16.0) then
        if (a.eq.16.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HH !opt 66
         dict=coef6(66)
         dict8=coef8(66)
         beta=damp(66)
         alph=sph(66)
         elem=25
        elseif (a.eq.16.0.and.b.eq.1.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !HC !opt 95
         !dict=coef6(95) !0.45*coef6(95)
         !dict8=coef8(95)!0.45*coef8(95)
         !beta=damp(95)
         !alph=sph(95)
         !elem=25
         dict=coef6(66)
         dict8=coef8(66)
         beta=damp(66)
         alph=sph(66)
         elem=25
!        elseif (a.eq.16.0.and.b.eq.6.0.and.c.eq.7.0.and.d.eq.8.0.and.e.eq.8.0.and.f.eq.0.0) then !CNOO
!         dict=dicta(27) !coef6(95) !dicta(27)
!         dict8=dicta8(27) !coef8(95) !dicta8(27)
!         beta=betaa(27) !damp(95) !betaa(27)
!        elseif (a.eq.16.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C
!         dict=coef6(66) !dicta(27)
!         dict8=coef8(66) !dicta8(27)
!         beta=damp(66) !betaa(27)
!         alph=sph(66)
        elseif (a.eq.16.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !CC
         dict=coef6(95)
         dict8=coef8(95)
         beta=damp(66) !damp(95)
         alph=sph(95)
         elem=25
        elseif (a.eq.16.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !OO
         dict=coef6(95)
         dict8=coef8(95)
         beta=damp(66)!damp(95)
         alph=sph(95)
         elem=25
!      elseif (a.eq.16.0.and.b.eq.6.0.and.c.eq.6.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !CCO
!         dict=1.5*coef6(95) ! 1.5 fact. was added for TCCP_8
!         dict8=1.5*coef8(95)
!         beta=damp(95)
!         alph=sph(95)
        elseif (a.eq.16.0.and.b.eq.1.0.and.c.eq.1.0.and.d.eq.8.0.and.e.eq.0.0.and.f.eq.0.0) then !HHO
         dict=coef6(95)
         dict8=coef8(95)
         beta=damp(66) !damp(95)
         alph=sph(95)
        else
         dict=(coef6(66)+coef6(95))/2.0
         dict8=(coef8(66)+coef8(95))/2.0
         beta=(damp(66)+damp(66))/2.0
         elem=25
        end if
!      elseif (a.eq.17.0.and.b.eq.0.0) then ! Chlorine atom
!         dict=dicta(28)
!         dict8=dicta8(28)
!         beta=betaa(28)
     elseif (a.eq.17.0) then 
        if (a.eq.17.0.and.b.eq.1.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !H !opt 67
         dict=coef6(67)
         dict8=coef8(67)
         beta=damp(67)
         alph=sph(67)
         elem=26
        elseif (a.eq.17.0.and.b.eq.5.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !B !opt 68
         dict=coef6(68)
         dict8=coef8(68)
         beta=0.35*damp(68)
         alph=sph(68)
         elem=26
        elseif (a.eq.17.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C !opt 69
         dict=coef6(69) !dicta(28) !coef6(69)
         dict8=coef8(69) !dicta8(28) !coef8(69)
         beta=0.4*damp(69) !betaa(28) !damp(69)
         alph=sph(69)
         elem=26
        elseif (a.eq.17.0.and.b.eq.7.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !N
         dict=coef6(69) !dicta(28)
         dict8=coef8(69) !dicta8(28)
         beta=damp(69) !betaa(28)
         alph=sph(69)
         elem=26
        elseif (a.eq.17.0.and.b.eq.9.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !F !opt 70
         dict=coef6(70)
         dict8=coef8(70)
         beta=damp(70)
         alph=sph(70)
         elem=26
        elseif (a.eq.17.0.and.b.eq.13.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Al !opt 71
         dict=coef6(71)
         dict8=coef8(71)
         beta=damp(71)
         alph=sph(71)
         elem=26
        elseif (a.eq.17.0.and.b.eq.14.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Si !opt 72
         dict=coef6(72)
         dict8=coef8(72)
         beta=0.35*damp(72)
         alph=sph(72)
         elem=26
        elseif (a.eq.17.0.and.b.eq.15.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !P !opt 73
         dict=coef6(73)
         dict8=coef8(73)
         beta=0.33*damp(73)
         alph=sph(73)
         elem=26
        elseif (a.eq.17.0.and.b.eq.17.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Cl !opt 74
         dict=coef6(74)
         dict8=coef8(74)
         beta=1.2*damp(74)
         alph=sph(74)
         elem=26
         !print*,beta,a,b
        else
         dict=(coef6(67)+coef6(68)+coef6(69)+coef6(70)+coef6(71)+coef6(72)+coef6(73)+coef6(74))/8.0
         dict8=(coef8(67)+coef8(68)+coef8(69)+coef8(70)+coef8(71)+coef8(72)+coef8(73)+coef8(74))/8.0
         beta=(damp(67)+0.5*damp(68)+0.4*damp(69)+damp(70)+damp(71)+0.5*damp(72)+0.33*damp(73)+1.2*damp(74))/8.0
         elem=26
         !print*,'beta',beta
        end if
!      elseif (a.eq.17.0.and.b.eq.8.0.and.c.eq.8.0.and.d.eq.8.0.and.e.eq.8.0.and.f.eq.0.0) then !OOOO
!         dict=5.0*dicta(28) !coef6(69)
!         dict8=5.0*dicta8(28) !coef8(69)
!         beta=betaa(28) !damp(69)
!
      elseif (a.eq.18.0.and.b.eq.0.0) then !opt 98
         dict=coef6(98)
         dict8=coef8(98)
         beta=damp(98)
         alph=sph(98)
         elem=27
!      elseif (a.eq.35.0.and.b.eq.0.0) then ! Bromine atom
!         dict=dicta(30)
!         dict8=dicta8(30)
!         beta=betaa(30)
      elseif (a.eq.35.0) then
        if (a.eq.35.0.and.b.eq.1.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !H !opt 75
         dict=coef6(75)
         dict8=coef8(75)
         beta=0.8*damp(75)
         alph=sph(75)
         elem=28
        elseif (a.eq.35.0.and.b.eq.6.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !C !opt 76
         dict=coef6(76)
         dict8=coef8(76)
         beta=damp(76)
         alph=sph(76)
         elem=28
        elseif (a.eq.35.0.and.b.eq.35.0.and.c.eq.0.0.and.d.eq.0.0.and.e.eq.0.0.and.f.eq.0.0) then !Br !opt 77
         dict=coef6(77)
         dict8=coef8(77)
         beta=1.25*damp(77)
         alph=sph(77)
         elem=28
         !print*,beta
        else
         dict=(coef6(75)+coef6(76)+coef6(77))/3.0
         dict8=(coef8(75)+coef8(76)+coef8(77))/3.0
         beta=(damp(75)+0.8*damp(76)+1.25*damp(77))/3.0
         elem=28
        end if
      elseif (a.eq.53.0) then
         dict=dicta(31)
         dict8=dicta8(31)
         beta=1.0*betaa(31)
      else
        print*,'Coefficients not found:',a,b,c,d,e,f
        stop 'error'
      end if      
   end if
!   print*,beta
   end Subroutine

   Function catom(dict)
!  """Returns C6 coefficient for atom at, converted to kcal*bohr^6/mol."""
   Real :: dict,catom

   catom=dict*10884.3258684610

   End Function

   Function catomatom(dicta,dictb)

   Real*8 :: c6a,c6b,catomatom
   Real :: dicta,dictb

   c6a = catom(dicta)
   c6b = catom(dictb)
!   print*,c6a/627.509474,c6b/627.509474
   catomatom = sqrt(c6a*c6b)

   End Function catomatom 


   Function catom8(dict8)
!  """Returns C8 coefficient for atom at, converted to kcal*bohr^8/mol."""
   Real*8 :: catom8
   Real :: dict8

   catom8=dict8*3886863.35457252

   End Function

   Function catomatom8(dict8a,dict8b)
!  """Returns C8 coefficient between atoms with charges at1 and at2."""
   Real*8 :: c8a,c8b,catomatom8,catom8
   Real :: dict8a,dict8b
   c8a=catom8(dict8a)
   c8b=catom8(dict8b)
!   print*,c8a/627.509474,c8b/627.509474

   catomatom8=sqrt(c8a*c8b)

   End Function

   Function damp_TT(r,beta1,beta2)
!   """Returns the value of Tang-Toennies damping function at separation r."""
!  it is assumed below that the damping factor is the geometric mean of atomic
!  damping factors:
   Real*8 :: alpha,r,br,damp_TT,sm,term
   Real :: beta1,beta2

   alpha=sqrt(abs(beta1*beta2))
   br=alpha*r
!   print*,alpha
   sm=1.0
   term=1.0
   Do i=1,6
     term=term*br/i
     sm=sm+term
   end do
   damp_TT=1.0 - Exp(-br)*sm
   
   End Function


   Function damp_TT_8(r,beta1,beta2)
!  """Returns the value of Tang-Toennies damping function f_8 at separation r."""
! it is assumed below that the damping factor is the geometric mean of atomic
! damping factors:
   Real*8 :: alpha,sm,br,damp_TT_8,r,term
   Real :: beta1,beta2
!   Real :: beta1,beta2
!   print*,beta1,beta2
   alpha=sqrt(abs(beta1*beta2))
   br=alpha*r
!   print*,beta1,beta2
!   print*,alpha
   sm=1.0
   term=1.0
   Do i=1,8
     term=term*br/i
     sm=sm+term
   end do
!   print*,sm
   damp_TT8=1.0 - Exp(-br)*sm

   End Function
 

   Function covalent_radius(att)
!  """Returns the covalent (Bragg) radius of the atom at, converted to bohr.
!  The radii are taken from the Cambridge Structural Database."""
!  print('at[0]=',at)
   Real*8 :: r,covalent_radius
   Integer :: att,at
    at=att
!    print*,'covalent_radius=',at
!    at=att
    if (at==2) then
      r=1.50d0
    elseif (at==3) then
      r=1.28d0
    elseif (at==4) then
      r=0.96d0
    elseif (at==5) then
      r=0.83d0
    elseif (at==6) then
      r=0.68d0
    elseif (at==7) then
      r=0.68d0
    elseif (at==8) then
      r=0.68d0
    elseif (at==9) then
      r=0.64d0
    elseif (at==10) then
      r=1.50d0
    elseif (at==11) then
      r=1.66d0
    elseif (at==12) then
      r=1.41d0
    elseif (at==13) then
      r=1.21d0
    elseif (at==14) then
      r=1.20d0
    elseif (at==15) then
      r=1.05d0
    elseif (at==16) then
      r=1.02d0
    elseif (at==17) then
      r=0.99d0
    elseif (at==18) then
      r=1.51d0
    elseif (at==35) then
      r=1.21d0
    elseif (at==53) then
      r=1.40d0
    elseif (at==1) then
      r=0.23d0
!    else
!      print*,'Wrong atom type in covalent_radius'
    end if
    covalent_radius=r/0.529177209d0
   End Function

   Function switching_function(r,at1,at2,aswitch,bswitch)
!  """Returns the value of the dispersion switching (Fermi) function for
!  atoms at1,at2 at separation r."""
   Real*8 :: r,aswitch,bswitch,r0,switching_function,covalent_radius,esph
   Integer :: at1,at2
!   dimension(5,nsmax) :: at1,at2
!   print*,covalent_radius(at1),covalent_radius(at2)
   r0=aswitch*(covalent_radius(at1)+covalent_radius(at2))
   switching_function=1.0/(1.0 + Exp(-bswitch*(r/r0-1.0)))
   End Function

   Function spher_function(r,beta1,beta2,alph1,alph2)
!
   Real*8 :: r,spher_function
   Real :: beta1,beta2,alph1,alph2,alpha,aab
!
   alpha=sqrt(beta1*beta2)
   aab=40*sqrt(alph1*alph2)
   spher_function=aab*Exp(-alpha*r)
   End Function

   Function disp_total(r,at1,at2,dict1,dict81,beta1,alph1,dict2,dict82,beta2,alph2,AAcoef,switch_on,aswitch,bswitch)

!  """Returns the overall contribution to dispersion
!  energy (in kcal/mol) for atoms at1,at2 at separation r."""
!  if not(at1 in dict.keys()):
!    print ("Dispersion coefficients not defined for atom",at1)
!    raise "Error!"
!  if not(at2 in dict.keys()):
!    print ("Dispersion coefficients not defined for atom",at2)
!    raise "Error!"

   Real*8 :: disp_total,edisp,r,aswitch,bswitch, damp_TT_8, catomatom, catomatom8, &
    damp_TT,switching_function,spher_function,edisp68,AAcoef 
   Real :: dict1,dict81,beta1,dict2,dict82,beta2,alph1,alph2,alpha
   Logical :: switch_on
   Integer :: at1,at2
!   dimension(5,nsmax) :: at1,at2
   !print*,dict1,dict2,beta1,dict81,dict82,beta2
   alpha=sqrt(beta1*beta2)
   edisp=-(catomatom(dict1,dict2)/(r*r*r*r*r*r))*damp_TT(r,beta1,beta2)  &
       -  (catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r))*damp_TT_8(r,beta1,beta2) !&
!   print*,catomatom(dict1,dict2)/627.509474,catomatom8(dict81,dict82)/627.509474
!   & + (AAcoef*Exp(-alpha*r))
!   print*,AAcoef*Exp(-alpha*r)
!   edisp68=-catomatom(dict1,dict2)/(r*r*r*r*r*r)*damp_TT(r,beta1,beta2)  &
!   &    -  catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r)*damp_TT8(r,beta1,beta2)
!   exponn=Exp(-alpha*r)
!   print*,'c6=',-catomatom(dict1,dict2)/(r*r*r*r*r*r)*damp_TT(r,beta1,beta2)
!   print*,'c8=',-catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r)*damp_TT8(r,beta1,beta2)
   !print*,'edisp=',edisp
!   print*,damp_TT(r,beta1,beta2)
!   if (switch_on) then

   disp_total=edisp*switching_function(r,at1,at2,aswitch,bswitch)
!   print*,'switching_function ',switching_function(r,at1,at2,aswitch,bswitch)
!   end if

   End Function

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! numerical Casimir--Polder integration
      Function trapzd(pol)
!       USE ISO_FORTRAN_ENV  , only: wp => REAL64
       real, intent(in) :: pol(23)
       real*8 :: trapzd

       real*8, parameter :: freq(23) = [&
       0.000001, 0.050000, 0.100000,&
       0.200000, 0.300000, 0.400000,&
       0.500000, 0.600000, 0.700000,&
       0.800000, 0.900000, 1.000000,&
       1.200000, 1.400000, 1.600000,&
       1.800000, 2.000000, 2.500000,&
       3.000000, 4.000000, 5.000000,&
       7.500000, 10.00000]
      real*8, parameter :: weights(23) = 0.50d0 * [&
       ( freq (2) - freq (1) ),&
       ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),&
       ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),&
       ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),&
       ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),&
       ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),&
       ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),&
       ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),&
       ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),&
       ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),&
       ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),&
       ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),&
       ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),&
       ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),&
       ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),&
       ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),&
       ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),&
       ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),&
       ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),&
       ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),&
       ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),&
       ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),&
       ( freq(23) - freq(22) ) ]

      trapzd = sum(pol*weights)
      end function trapzd
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Function damp_TT6(r,bb)
!   """Returns the value of Tang-Toennies damping function at separation
!   r."""
!  it is assumed below that the damping factor is the geometric mean of
!  atomic damping factors:
      Real*8 :: alpha,r,br,damp_TT6,sm,term
      Real*8 :: beta1,beta2,cp1,bb
      Integer :: i

      !alpha=sqrt(abs(beta1*beta2))
      !br=cp1*alpha*r
      br=bb*r
!     print*,alpha
      sm=1.0
      term=1.0
      Do i=1,6
        term=term*br/i
        sm=sm+term
      end do
      damp_TT6=1.0 - Exp(-br)*sm

      End Function damp_TT6

      Function damp_TT8(r,bb)
!  """Returns the value of Tang-Toennies damping function f_8 at
!  separation r."""
! it is assumed below that the damping factor is the geometric mean of
! atomic
! damping factors:
      Real*8 :: alpha,sm,br,damp_TT8,r,term
      Real*8 :: beta1,beta2,cp2,bb
      Integer :: i
!     Real :: beta1,beta2
!     print*,beta1,beta2
      !alpha=sqrt(abs(beta1*beta2))
      !br=cp2*alpha*r
      br=r*bb
!     print*,beta1,beta2
!     print*,alpha
      sm=1.0
      term=1.0
      Do i=1,8
        term=term*br/i
        sm=sm+term
      end do
!
      damp_TT8=1.0 - Exp(-br)*sm

      End Function damp_TT8

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   Subroutine fitted_dispersion(laa,lab,r0ab,nsa,nsb,Veff3_a,Veff3_b,Veff5_a,Veff5_b,disper)
!  """Calculates the fitted dispersion energy, either within a molecule
!     or between two molecules."""
!  print('laa-fitted_dispersion',laa)
!  print('lab-fitted_dispersion',lab)
   integer, parameter :: vers = 3, nsmax = 1000
!   Parameter (vers=3,nsmax=100)
   Real*8 :: dispfit,x,y,disper,laa,lab,r,disp_total,AA0,expon,r2,r6,r8,&
   & edisp68,AAcoef,c6,c8,c6aa,c6bb,c8aa,c8bb,tmp,e6,e8,cp1,cp2,cp3,cp4, c11, c12, c22
   Real :: alpha_d(23),alpha_q(23)
   Real :: alpha_daa(23),alpha_qaa(23),alpha_dbb(23),alpha_qbb(23)
   
   Real :: aiw_a(1,23),aiw_b(1,23),aiwq_a(1,23),aiwq_b(1,23)
   Integer :: at1,at11,at2,at22,nsa,nsb,i,j,elem1,elem2
   Dimension :: x(10,nsmax),y(10,nsmax),laa(10,nsmax),lab(10,nsmax)
   Real*8, Dimension(99) :: coef6,coef8,damp,sph
   Real*8 :: AA(98,98),Veff3_a(nsmax),Veff3_b(nsmax),Veff5_a(nsmax),Veff5_b(nsmax)
   Real*8 :: Veff3_aa,Veff3_bb,Veff5_aa,Veff5_bb,damp_TT6,damp_TT8,r0ab(94,94),bb6,bb8
   Real*8, parameter :: pi =4.0d0*atan(1.0d0)
   Real*8, parameter :: thopi = 3.0d0/pi, fiftotwopi=15.0d0/(2.d0*pi)
   real(8), external :: trapzd
   logical :: onesystem,switch_on
!   onesystem=.true.
   onesystem=.false.

   open(81,file='param_damp.dat')
   read(81,*) cp1,cp2,cp3,cp4
   close(81)
   !write(*,*)'---Damping parameters---'
   !write(*,*)cp1,cp2,cp3,cp4
   !write(*,*)'------------------------'
   open(82,file='param_ckl.dat')
   read(82,*) c11,c12,c22
   close(82)

   x=laa
   y=lab
   e6 = 0.0
   e8 = 0.0
   !dispfit=0.0
!   c11=0.0
!   c12=0.0
!   c22=0.0

!   if(onesystem) then
!      Do i = 1,nsa
!         at1=x(4,i)
!         Call Coefficients(x,coef6,coef8,damp,sph,dict1,dict81,beta1,alph1,i,elem1,vers)
!         Do j = 1,nsa
!            at2=x(4,j)
!            Call Coefficients(x,coef6,coef8,damp,sph,dict2,dict82,beta2,alph2,j,elem2,vers)
!            r=sqrt((x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 + (x(3,i)-x(3,j))**2)
!            if(r.gt.0.0) then
!               dispfit=dispfit + disp_total(r,at1,at2,dict1,dict81,beta1,alph1,dict2,dict82,beta2,alph2,AAcoef,switch_on,aswitch,bswitch)
!
!            end if
!         end Do
!      end Do
!
!   else
   write(5,'(2A8,6A18)') "Atom 1","Atom 2","AIM_d_pol^{a}","AIM_q_pol^{a}","AIM_d_pol^{b}",&
   "AIM_q_pol^{b}","C6^{ab}","C8^{ab}"
   Do i = 1,nsa
      at1=x(4,i)
      Veff3_aa=Veff3_a(i)
      Veff5_aa=Veff5_a(i)
!      Call Coefficients(x,coef6,coef8,damp,sph,dict1,dict81,beta1,alph1,i,elem1,vers)
      !print*,'Veff3_a(i), Veff5_a(i)',Veff3_a(i), Veff5_a(i)
      !Call AIM_polarizabilities(x, i, Veff3_a(i), Veff5_a(i), c11, c12, c22, aiw_a, aiwq_a)
      Call AIM_polarizabilities(x, i, Veff3_aa, Veff5_aa, c11, c12, c22, aiw_a, aiwq_a)
      c6aa=0
      c8aa=0
      !alpha_daa = 0.0
      !alpha_qaa = 0.0
      !alpha_daa = aiw_a(1,:)*aiw_a(1,:)
      !alpha_qaa = 2*aiw_a(1,:)*aiwq_a(1,:)
      !c6aa = thopi*trapzd(alpha_daa)
      !c8aa = fiftotwopi*trapzd(alpha_qaa)
      Do j = 1,nsb
         at2=y(4,j)
         Veff3_bb=Veff3_b(j)
         Veff5_bb=Veff5_b(j)
!         Call Coefficients(y,coef6,coef8,damp,sph,dict2,dict82,beta2,alph2,j,elem2,vers)
         !Call AIM_polarizabilities(y, j, Veff3_b(j), Veff5_b(j), c11, c12, c22, aiw_b,aiwq_b)
         Call AIM_polarizabilities(y, j, Veff3_bb, Veff5_bb, c11, c12, c22, aiw_b,aiwq_b)
         r=sqrt((x(1,i)-y(1,j))**2 + (x(2,i)-y(2,j))**2 + (x(3,i)-y(3,j))**2)
         r2=r*r
         r6=r2**3
         r8=r6*r2
         !c6bb=0
         !c8bb=0
         c6=0.0
         c8=0.0
         alpha_d=0.0
         alpha_q=0.0
         !alpha_dbb=0.0
         !alpha_qbb=0.0
         alpha_d=aiw_a(1,:)*aiw_b(1,:)
         alpha_q=aiw_a(1,:)*aiwq_b(1,:) + aiwq_a(1,:)*aiw_b(1,:)

         !alpha_dbb=aiw_b(1,:)*aiw_b(1,:)
         !alpha_qbb=2*aiw_b(1,:)*aiwq_b(1,:)
         !c6bb = thopi*trapzd(alpha_dbb)
         !c8bb = fiftotwopi*trapzd(alpha_qbb)
         !print*,'i,j,Veff3_a(i),Veff3_b(j),Veff5_a(i),Veff5_b(j)',i,j,Veff3_a(i),Veff3_b(j),Veff5_a(i),Veff5_b(j) 
         !print*,at1,at2,alpha_d(1),alpha_q(1)
         !Do i=1,23
         !   alpha_q(i)=aiw(iat,i)*aiwq(jat,i)+aiwq(iat,i)*aiw(jat,i)
         !   print*,iat,aiw(iat,i),aiwq(iat,i)
         !end do
 
         c6 = thopi*trapzd(alpha_d)
         c8 = fiftotwopi*trapzd(alpha_q)
         
         write(5,'(2I8,6F18.3)') at1,at2,aiw_a(1,1),aiwq_a(1,1), aiw_b(1,1),aiwq_b(1,1),c6,c8
         
         
         tmp = sqrt(c8/c6)
         
         !bb6=(cp1*r0ab(at1,at2))+cp2
         !bb8=(cp3*r0ab(at1,at2))+cp4
         
         e6 = e6 +  c6/(r6+(cp1*tmp + cp2)**6) !((c6*damp_TT6(r,bb6))/r6) !c6/(r6+(cp1*tmp + cp2)**6)
         e8 = e8 +  c8/(r8+(cp3*tmp + cp4)**8) !((c8*damp_TT8(r,bb8))/r8) !c8/(r8+(cp3*tmp + cp4)**8)
 

         !dispfit=dispfit + disp_total(r,at1,at2,dict1,dict81,beta1,alph1,dict2,dict82,beta2,alph2,AAcoef,switch_on,aswitch,bswitch)
      
      end Do

   end Do
!   end if
   
   disper=-e6-e8 
   
   End subroutine 



   
   Subroutine analyze_carbons(latom,nat,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)


   Do i=1,nat
      atom=latom
      if (atom(4,i).ne.6.0) then
         cycle
      end if
      nneighbors=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.3.0) then
            nneighbors=nneighbors+1
         end if
      end Do
      sitout(5,i)=nneighbors
      !print*,'nneighbors=',nneighbors
   end Do

   end subroutine

   Subroutine analyze_Boron(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov

   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.5.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_carbons2(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Carbon neighbours'
   Do i=1,nat
      connected=0
      atom=latom
      if (atom(4,i).ne.6.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         print*,4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      !sitout(5,i)=nneighbors
      !print*,'nneighbors=',nneighbors
      ! Loading the atoms connected with carbon atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      !print*,'Before sorting'
      !write(*,*)(nonzeros(k),k=1,nonzero_count)
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'after sorting'
      !print*,size(connected)
      !write(*,'(6I)')(connected(k),k=1,size(connected))
      !write(*,*)(nonzeros(k),k=1,nonzero_count)
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do

   end subroutine


! Analyze Nitrogen i.e finding first neighbours of Nitrogen atom

   Subroutine analyze_Nitrogen(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Nitrogen neighbours'
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.7.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

! Analyze Oxygen i.e finding first neighbours of Oxygen atom

   Subroutine analyze_Oxygen(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Oxygen neighbours'
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.8.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

! Analyze Flourine i.e finding first neighbours of Flourine atom

   Subroutine analyze_Flourine(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Florine neighbours'
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.9.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_Aluminium(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Nitrogen neighbours'
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.13.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_Silicon(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.14.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_Phosphorus(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.15.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_Sulfur(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.16.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

! Analyze Chlorine i.e finding first neighbours of Chlorine atom

   Subroutine analyze_Chlorine(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   !print*,'Chlorine neighbours'
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.17.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         !print*,'rcov',rcov(int(atom(4,i))),rcov(int(latom(4,j))),4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            !print*,'rcov',4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine

   Subroutine analyze_Bromine(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.35.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine


   Subroutine analyze_Iodine(latom,nat,rcov,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ,cnt
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)
   Integer, dimension(6) :: connected
   integer, allocatable :: array(:), nonzeros(:)
   integer :: i, nonzero_count, total_size
   Real*8,Dimension(94) :: rcov
   Do i=1,nat
      connected=0.0
      atom=latom
      if (atom(4,i).ne.53.0) then
         cycle
      end if
      nneighbors=0
      cnt=0
      Do j=1,nat
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
         if (r.eq.0.0) then
            cycle
         elseif (r.lt.4.0d0/5.0d0*(rcov(int(atom(4,i)))+rcov(int(latom(4,j))))) then
            cnt=cnt+1
            nneighbors=nneighbors+1
            connected(cnt)=latom(4,j)
         end if
      end Do
      ! Loading the atoms connected with Nitrogen atom
      nonzero_count = count(connected /= 0)
      total_size = size(connected)
      allocate(nonzeros(nonzero_count))
      nonzeros = pack(connected, connected /= 0)  ! Extract nonzero values
      call sort_array(nonzeros,nonzero_count)
      ! Reconstruct the array: sorted nonzeros followed by zeros
      connected = 0.0  ! Initialize the array with zeros
      connected(1:nonzero_count) = nonzeros
      !print*,'Nitrogen neighbours'
      !write(*,'(6F6.2)')(connected(k),k=1,size(connected))
      Do k =1, 6
         sitout(k+4,i)=connected(k)
      end do
      deallocate(nonzeros)
   end Do
   end subroutine


   subroutine sort_array(arr,nn)
   integer , intent(inout) :: arr
   integer :: i, j, temp,nn
   dimension :: arr(nn)
   !print*,'subroutine sort_array',arr(1),arr(2),arr(3)
   do i = 1, size(arr) - 1
        do j = i + 1, size(arr)
            if (arr(i) > arr(j)) then
                 temp = arr(i)
                 arr(i) = arr(j)
                 arr(j) = temp
            end if
        end do
   end do
   end subroutine sort_array


   subroutine analyze_hydrogens(latom,nat,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=100,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr
   Integer :: nat,typ
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   dimension :: sita(10,nsmax),sitb(10,nsmax),sitout(10,nsmax),sitbb(10,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(10,nsmax) :: latom(10,nsmax),atom(10,nsmax)

!   print*,'atta',latom(:,1)
!   sita=latom
!   sitb=latom
   do i=1,nat
      atom=latom
      if (atom(4,i).ne.1.0) then
!         print*,'atom(1,i)=',atom(1,i)
         cycle
      end if
!      else
      closestr=100.0
!      print*,atom(:,i)
      Do j=1,nat
!         print*,'latom=',latom(1,j)
         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
         & (atom(3,i)-latom(3,j))**2)
!         print*,atom(1,i)
!         print*,'r=',r
!         print*,'closestr=',closestr
         if (r.eq.0.0) then
!            print*,atom(1,i)
            cycle
!         print*,'closestr=',closestr
!         print*,'r=',r
         elseif (r.lt.closestr) then
!            print*,'closestr=',closestr
!            print*,'r=',r
            closestr=r
            typ=latom(4,j)
!            print*,'type***************=',typ
         end if
!         print*,'closestr=',closestr
!         print*,'r=',r
      end do
      if (typ.eq.1) then
          !print*,"H connected to H"
          sitout(5,i)=typ
!          print*,'atom(5,i)*********************',sitaa(5,i)
      elseif (typ.eq.3) then
!     print "H connected to Li"
          sitout(5,i)=typ
      elseif (typ.eq.4) then
!     print "H connected to Be"
          sitout(5,i)=typ
      elseif (typ.eq.5) then
!     print "H connected to B"
          sitout(5,i)=typ
      elseif (typ.eq.6) then
!     print "H connected to C"
          sitout(5,i)=typ
!          sitaa(5,i)=int(sitaa(5,i))
!          print*,'a****',a
!          print*,'sitaa(5,i)**********************',sitaa(5,i)
      elseif (typ.eq.7) then
!     print "H connected to N"
          sitout(5,i)=typ
      elseif (typ.eq.8) then
          !print*,"H connected to O"
          sitout(5,i)=typ
!          print*,'sitaa(5,i)**********************',sitaa(5,i)
      elseif (typ.eq.9) then
!     print "H connected to F"
          sitout(5,i)=typ
      elseif (typ.eq.11) then
!     print "H connected to Na"
          sitout(5,i)=typ
      elseif (typ.eq.12) then
!     print "H connected to Mg"
          sitout(5,i)=typ
      elseif (typ.eq.13) then
!     print "H connected to Al"
          sitout(5,i)=typ
      elseif (typ.eq.14) then
!     print "H connected to Si"
          sitout(5,i)=typ
      elseif (typ.eq.15) then
!     print "H connected to P"
          sitout(5,i)=typ
      elseif (typ.eq.16) then
!     print "H connected to S"
          sitout(5,i)=typ
      elseif (typ.eq.17) then
!     print "H connected to Cl"
          sitout(5,i)=typ
      elseif (typ.eq.35) then
!     print "H connected to Br"
          sitout(5,i)=typ
      elseif (typ.eq.53) then
!     print "H connected to I"
          sitout(5,i)=typ
!      elseif (typ.eq.) 1,3,4,5,6,7,8,9,11,12,13,14,15,16,17 then
!          print*,'Wrong neighbour! ',typ
      end if
   end do
   end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C set cut-off radii
!C in parts due to INTEL compiler bug
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine setr0ab(max_elem,autoang,r)
      implicit none
      integer max_elem,i,j,k
      real*8 r(max_elem,max_elem),autoang
      real*8 r0ab(4465)
      r0ab(   1:  70)=(/ &
         2.1823,  1.8547,  1.7347,  2.9086,  2.5732,  3.4956,  2.3550 &
      ,  2.5095,  2.9802,  3.0982,  2.5141,  2.3917,  2.9977,  2.9484 &
      ,  3.2160,  2.4492,  2.2527,  3.1933,  3.0214,  2.9531,  2.9103 &
      ,  2.3667,  2.1328,  2.8784,  2.7660,  2.7776,  2.7063,  2.6225 &
      ,  2.1768,  2.0625,  2.6395,  2.6648,  2.6482,  2.5697,  2.4846 &
      ,  2.4817,  2.0646,  1.9891,  2.5086,  2.6908,  2.6233,  2.4770 &
      ,  2.3885,  2.3511,  2.2996,  1.9892,  1.9251,  2.4190,  2.5473 &
      ,  2.4994,  2.4091,  2.3176,  2.2571,  2.1946,  2.1374,  2.9898 &
      ,  2.6397,  3.6031,  3.1219,  3.7620,  3.2485,  2.9357,  2.7093 &
      ,  2.5781,  2.4839,  3.7082,  2.5129,  2.7321,  3.1052,  3.2962 &
      /)
      r0ab(  71: 140)=(/ &
         3.1331,  3.2000,  2.9586,  3.0822,  2.8582,  2.7120,  3.2570 &
      ,  3.4839,  2.8766,  2.7427,  3.2776,  3.2363,  3.5929,  3.2826 &
      ,  3.0911,  2.9369,  2.9030,  2.7789,  3.3921,  3.3970,  4.0106 &
      ,  2.8884,  2.6605,  3.7513,  3.1613,  3.3605,  3.3325,  3.0991 &
      ,  2.9297,  2.8674,  2.7571,  3.8129,  3.3266,  3.7105,  3.7917 &
      ,  2.8304,  2.5538,  3.3932,  3.1193,  3.1866,  3.1245,  3.0465 &
      ,  2.8727,  2.7664,  2.6926,  3.4608,  3.2984,  3.5142,  3.5418 &
      ,  3.5017,  2.6190,  2.4797,  3.1331,  3.0540,  3.0651,  2.9879 &
      ,  2.9054,  2.8805,  2.7330,  2.6331,  3.2096,  3.5668,  3.3684 &
      ,  3.3686,  3.3180,  3.3107,  2.4757,  2.4019,  2.9789,  3.1468 &
      /)
      r0ab( 141: 210)=(/ &
         2.9768,  2.8848,  2.7952,  2.7457,  2.6881,  2.5728,  3.0574 &
      ,  3.3264,  3.3562,  3.2529,  3.1916,  3.1523,  3.1046,  2.3725 &
      ,  2.3289,  2.8760,  2.9804,  2.9093,  2.8040,  2.7071,  2.6386 &
      ,  2.5720,  2.5139,  2.9517,  3.1606,  3.2085,  3.1692,  3.0982 &
      ,  3.0352,  2.9730,  2.9148,  3.2147,  2.8315,  3.8724,  3.4621 &
      ,  3.8823,  3.3760,  3.0746,  2.8817,  2.7552,  2.6605,  3.9740 &
      ,  3.6192,  3.6569,  3.9586,  3.6188,  3.3917,  3.2479,  3.1434 &
      ,  4.2411,  2.7597,  3.0588,  3.3474,  3.6214,  3.4353,  3.4729 &
      ,  3.2487,  3.3200,  3.0914,  2.9403,  3.4972,  3.7993,  3.6773 &
      ,  3.8678,  3.5808,  3.8243,  3.5826,  3.4156,  3.8765,  4.1035 &
      /)
      r0ab( 211: 280)=(/ &
         2.7361,  2.9765,  3.2475,  3.5004,  3.4185,  3.4378,  3.2084 &
      ,  3.2787,  3.0604,  2.9187,  3.4037,  3.6759,  3.6586,  3.8327 &
      ,  3.5372,  3.7665,  3.5310,  3.3700,  3.7788,  3.9804,  3.8903 &
      ,  2.6832,  2.9060,  3.2613,  3.4359,  3.3538,  3.3860,  3.1550 &
      ,  3.2300,  3.0133,  2.8736,  3.4024,  3.6142,  3.5979,  3.5295 &
      ,  3.4834,  3.7140,  3.4782,  3.3170,  3.7434,  3.9623,  3.8181 &
      ,  3.7642,  2.6379,  2.8494,  3.1840,  3.4225,  3.2771,  3.3401 &
      ,  3.1072,  3.1885,  2.9714,  2.8319,  3.3315,  3.5979,  3.5256 &
      ,  3.4980,  3.4376,  3.6714,  3.4346,  3.2723,  3.6859,  3.8985 &
      ,  3.7918,  3.7372,  3.7211,  2.9230,  2.6223,  3.4161,  2.8999 &
      /)
      r0ab( 281: 350)=(/ &
         3.0557,  3.3308,  3.0555,  2.8508,  2.7385,  2.6640,  3.5263 &
      ,  3.0277,  3.2990,  3.7721,  3.5017,  3.2751,  3.1368,  3.0435 &
      ,  3.7873,  3.2858,  3.2140,  3.1727,  3.2178,  3.4414,  2.5490 &
      ,  2.7623,  3.0991,  3.3252,  3.1836,  3.2428,  3.0259,  3.1225 &
      ,  2.9032,  2.7621,  3.2490,  3.5110,  3.4429,  3.3845,  3.3574 &
      ,  3.6045,  3.3658,  3.2013,  3.6110,  3.8241,  3.7090,  3.6496 &
      ,  3.6333,  3.0896,  3.5462,  2.4926,  2.7136,  3.0693,  3.2699 &
      ,  3.1272,  3.1893,  2.9658,  3.0972,  2.8778,  2.7358,  3.2206 &
      ,  3.4566,  3.3896,  3.3257,  3.2946,  3.5693,  3.3312,  3.1670 &
      ,  3.5805,  3.7711,  3.6536,  3.5927,  3.5775,  3.0411,  3.4885 &
      /)
      r0ab( 351: 420)=(/ &
         3.4421,  2.4667,  2.6709,  3.0575,  3.2357,  3.0908,  3.1537 &
      ,  2.9235,  3.0669,  2.8476,  2.7054,  3.2064,  3.4519,  3.3593 &
      ,  3.2921,  3.2577,  3.2161,  3.2982,  3.1339,  3.5606,  3.7582 &
      ,  3.6432,  3.5833,  3.5691,  3.0161,  3.4812,  3.4339,  3.4327 &
      ,  2.4515,  2.6338,  3.0511,  3.2229,  3.0630,  3.1265,  2.8909 &
      ,  3.0253,  2.8184,  2.6764,  3.1968,  3.4114,  3.3492,  3.2691 &
      ,  3.2320,  3.1786,  3.2680,  3.1036,  3.5453,  3.7259,  3.6090 &
      ,  3.5473,  3.5327,  3.0018,  3.4413,  3.3907,  3.3593,  3.3462 &
      ,  2.4413,  2.6006,  3.0540,  3.1987,  3.0490,  3.1058,  2.8643 &
      ,  2.9948,  2.7908,  2.6491,  3.1950,  3.3922,  3.3316,  3.2585 &
      /)
      r0ab( 421: 490)=(/ &
         3.2136,  3.1516,  3.2364,  3.0752,  3.5368,  3.7117,  3.5941 &
      ,  3.5313,  3.5164,  2.9962,  3.4225,  3.3699,  3.3370,  3.3234 &
      ,  3.3008,  2.4318,  2.5729,  3.0416,  3.1639,  3.0196,  3.0843 &
      ,  2.8413,  2.7436,  2.7608,  2.6271,  3.1811,  3.3591,  3.3045 &
      ,  3.2349,  3.1942,  3.1291,  3.2111,  3.0534,  3.5189,  3.6809 &
      ,  3.5635,  3.5001,  3.4854,  2.9857,  3.3897,  3.3363,  3.3027 &
      ,  3.2890,  3.2655,  3.2309,  2.8502,  2.6934,  3.2467,  3.1921 &
      ,  3.5663,  3.2541,  3.0571,  2.9048,  2.8657,  2.7438,  3.3547 &
      ,  3.3510,  3.9837,  3.6871,  3.4862,  3.3389,  3.2413,  3.1708 &
      ,  3.6096,  3.6280,  3.6860,  3.5568,  3.4836,  3.2868,  3.3994 &
      /)
      r0ab( 491: 560)=(/ &
         3.3476,  3.3170,  3.2950,  3.2874,  3.2606,  3.9579,  2.9226 &
      ,  2.6838,  3.7867,  3.1732,  3.3872,  3.3643,  3.1267,  2.9541 &
      ,  2.8505,  2.7781,  3.8475,  3.3336,  3.7359,  3.8266,  3.5733 &
      ,  3.3959,  3.2775,  3.1915,  3.9878,  3.8816,  3.5810,  3.5364 &
      ,  3.5060,  3.8097,  3.3925,  3.3348,  3.3019,  3.2796,  3.2662 &
      ,  3.2464,  3.7136,  3.8619,  2.9140,  2.6271,  3.4771,  3.1774 &
      ,  3.2560,  3.1970,  3.1207,  2.9406,  2.8322,  2.7571,  3.5455 &
      ,  3.3514,  3.5837,  3.6177,  3.5816,  3.3902,  3.2604,  3.1652 &
      ,  3.7037,  3.6283,  3.5858,  3.5330,  3.4884,  3.5789,  3.4094 &
      ,  3.3473,  3.3118,  3.2876,  3.2707,  3.2521,  3.5570,  3.6496 &
      /)
      r0ab( 561: 630)=(/ &
         3.6625,  2.7300,  2.5870,  3.2471,  3.1487,  3.1667,  3.0914 &
      ,  3.0107,  2.9812,  2.8300,  2.7284,  3.3259,  3.3182,  3.4707 &
      ,  3.4748,  3.4279,  3.4182,  3.2547,  3.1353,  3.5116,  3.9432 &
      ,  3.8828,  3.8303,  3.7880,  3.3760,  3.7218,  3.3408,  3.3059 &
      ,  3.2698,  3.2446,  3.2229,  3.4422,  3.5023,  3.5009,  3.5268 &
      ,  2.6026,  2.5355,  3.1129,  3.2863,  3.1029,  3.0108,  2.9227 &
      ,  2.8694,  2.8109,  2.6929,  3.1958,  3.4670,  3.4018,  3.3805 &
      ,  3.3218,  3.2815,  3.2346,  3.0994,  3.3937,  3.7266,  3.6697 &
      ,  3.6164,  3.5730,  3.2522,  3.5051,  3.4686,  3.4355,  3.4084 &
      ,  3.3748,  3.3496,  3.3692,  3.4052,  3.3910,  3.3849,  3.3662 &
      /)
      r0ab( 631: 700)=(/ &
         2.5087,  2.4814,  3.0239,  3.1312,  3.0535,  2.9457,  2.8496 &
      ,  2.7780,  2.7828,  2.6532,  3.1063,  3.3143,  3.3549,  3.3120 &
      ,  3.2421,  3.1787,  3.1176,  3.0613,  3.3082,  3.5755,  3.5222 &
      ,  3.4678,  3.4231,  3.1684,  3.3528,  3.3162,  3.2827,  3.2527 &
      ,  3.2308,  3.2029,  3.3173,  3.3343,  3.3092,  3.2795,  3.2452 &
      ,  3.2096,  3.2893,  2.8991,  4.0388,  3.6100,  3.9388,  3.4475 &
      ,  3.1590,  2.9812,  2.8586,  2.7683,  4.1428,  3.7911,  3.8225 &
      ,  4.0372,  3.7059,  3.4935,  3.3529,  3.2492,  4.4352,  4.0826 &
      ,  3.9733,  3.9254,  3.8646,  3.9315,  3.7837,  3.7465,  3.7211 &
      ,  3.7012,  3.6893,  3.6676,  3.7736,  4.0660,  3.7926,  3.6158 &
      /)
      r0ab( 701: 770)=(/ &
         3.5017,  3.4166,  4.6176,  2.8786,  3.1658,  3.5823,  3.7689 &
      ,  3.5762,  3.5789,  3.3552,  3.4004,  3.1722,  3.0212,  3.7241 &
      ,  3.9604,  3.8500,  3.9844,  3.7035,  3.9161,  3.6751,  3.5075 &
      ,  4.1151,  4.2877,  4.1579,  4.1247,  4.0617,  3.4874,  3.9848 &
      ,  3.9280,  3.9079,  3.8751,  3.8604,  3.8277,  3.8002,  3.9981 &
      ,  3.7544,  4.0371,  3.8225,  3.6718,  4.3092,  4.4764,  2.8997 &
      ,  3.0953,  3.4524,  3.6107,  3.6062,  3.5783,  3.3463,  3.3855 &
      ,  3.1746,  3.0381,  3.6019,  3.7938,  3.8697,  3.9781,  3.6877 &
      ,  3.8736,  3.6451,  3.4890,  3.9858,  4.1179,  4.0430,  3.9563 &
      ,  3.9182,  3.4002,  3.8310,  3.7716,  3.7543,  3.7203,  3.7053 &
      /)
      r0ab( 771: 840)=(/ &
         3.6742,  3.8318,  3.7631,  3.7392,  3.9892,  3.7832,  3.6406 &
      ,  4.1701,  4.3016,  4.2196,  2.8535,  3.0167,  3.3978,  3.5363 &
      ,  3.5393,  3.5301,  3.2960,  3.3352,  3.1287,  2.9967,  3.6659 &
      ,  3.7239,  3.8070,  3.7165,  3.6368,  3.8162,  3.5885,  3.4336 &
      ,  3.9829,  4.0529,  3.9584,  3.9025,  3.8607,  3.3673,  3.7658 &
      ,  3.7035,  3.6866,  3.6504,  3.6339,  3.6024,  3.7708,  3.7283 &
      ,  3.6896,  3.9315,  3.7250,  3.5819,  4.1457,  4.2280,  4.1130 &
      ,  4.0597,  3.0905,  2.7998,  3.6448,  3.0739,  3.2996,  3.5262 &
      ,  3.2559,  3.0518,  2.9394,  2.8658,  3.7514,  3.2295,  3.5643 &
      ,  3.7808,  3.6931,  3.4723,  3.3357,  3.2429,  4.0280,  3.5589 &
      /)
      r0ab( 841: 910)=(/ &
         3.4636,  3.4994,  3.4309,  3.6177,  3.2946,  3.2376,  3.2050 &
      ,  3.1847,  3.1715,  3.1599,  3.5555,  3.8111,  3.7693,  3.5718 &
      ,  3.4498,  3.3662,  4.1608,  3.7417,  3.6536,  3.6154,  3.8596 &
      ,  3.0301,  2.7312,  3.5821,  3.0473,  3.2137,  3.4679,  3.1975 &
      ,  2.9969,  2.8847,  2.8110,  3.6931,  3.2076,  3.4943,  3.5956 &
      ,  3.6379,  3.4190,  3.2808,  3.1860,  3.9850,  3.5105,  3.4330 &
      ,  3.3797,  3.4155,  3.6033,  3.2737,  3.2145,  3.1807,  3.1596 &
      ,  3.1461,  3.1337,  3.4812,  3.6251,  3.7152,  3.5201,  3.3966 &
      ,  3.3107,  4.1128,  3.6899,  3.6082,  3.5604,  3.7834,  3.7543 &
      ,  2.9189,  2.6777,  3.4925,  2.9648,  3.1216,  3.2940,  3.0975 &
      /)
      r0ab( 911: 980)=(/ &
         2.9757,  2.8493,  2.7638,  3.6085,  3.1214,  3.4006,  3.4793 &
      ,  3.5147,  3.3806,  3.2356,  3.1335,  3.9144,  3.4183,  3.3369 &
      ,  3.2803,  3.2679,  3.4871,  3.1714,  3.1521,  3.1101,  3.0843 &
      ,  3.0670,  3.0539,  3.3890,  3.5086,  3.5895,  3.4783,  3.3484 &
      ,  3.2559,  4.0422,  3.5967,  3.5113,  3.4576,  3.6594,  3.6313 &
      ,  3.5690,  2.8578,  2.6334,  3.4673,  2.9245,  3.0732,  3.2435 &
      ,  3.0338,  2.9462,  2.8143,  2.7240,  3.5832,  3.0789,  3.3617 &
      ,  3.4246,  3.4505,  3.3443,  3.1964,  3.0913,  3.8921,  3.3713 &
      ,  3.2873,  3.2281,  3.2165,  3.4386,  3.1164,  3.1220,  3.0761 &
      ,  3.0480,  3.0295,  3.0155,  3.3495,  3.4543,  3.5260,  3.4413 &
      /)
      r0ab( 981:1050)=(/ &
         3.3085,  3.2134,  4.0170,  3.5464,  3.4587,  3.4006,  3.6027 &
      ,  3.5730,  3.4945,  3.4623,  2.8240,  2.5960,  3.4635,  2.9032 &
      ,  3.0431,  3.2115,  2.9892,  2.9148,  2.7801,  2.6873,  3.5776 &
      ,  3.0568,  3.3433,  3.3949,  3.4132,  3.3116,  3.1616,  3.0548 &
      ,  3.8859,  3.3719,  3.2917,  3.2345,  3.2274,  3.4171,  3.1293 &
      ,  3.0567,  3.0565,  3.0274,  3.0087,  2.9939,  3.3293,  3.4249 &
      ,  3.4902,  3.4091,  3.2744,  3.1776,  4.0078,  3.5374,  3.4537 &
      ,  3.3956,  3.5747,  3.5430,  3.4522,  3.4160,  3.3975,  2.8004 &
      ,  2.5621,  3.4617,  2.9154,  3.0203,  3.1875,  2.9548,  2.8038 &
      ,  2.7472,  2.6530,  3.5736,  3.0584,  3.3304,  3.3748,  3.3871 &
      /)
      r0ab(1051:1120)=(/ &
         3.2028,  3.1296,  3.0214,  3.8796,  3.3337,  3.2492,  3.1883 &
      ,  3.1802,  3.4050,  3.0756,  3.0478,  3.0322,  3.0323,  3.0163 &
      ,  3.0019,  3.3145,  3.4050,  3.4656,  3.3021,  3.2433,  3.1453 &
      ,  3.9991,  3.5017,  3.4141,  3.3520,  3.5583,  3.5251,  3.4243 &
      ,  3.3851,  3.3662,  3.3525,  2.7846,  2.5324,  3.4652,  2.8759 &
      ,  3.0051,  3.1692,  2.9273,  2.7615,  2.7164,  2.6212,  3.5744 &
      ,  3.0275,  3.3249,  3.3627,  3.3686,  3.1669,  3.0584,  2.9915 &
      ,  3.8773,  3.3099,  3.2231,  3.1600,  3.1520,  3.4023,  3.0426 &
      ,  3.0099,  2.9920,  2.9809,  2.9800,  2.9646,  3.3068,  3.3930 &
      ,  3.4486,  3.2682,  3.1729,  3.1168,  3.9952,  3.4796,  3.3901 &
      /)
      r0ab(1121:1190)=(/ &
         3.3255,  3.5530,  3.5183,  3.4097,  3.3683,  3.3492,  3.3360 &
      ,  3.3308,  2.5424,  2.6601,  3.2555,  3.2807,  3.1384,  3.1737 &
      ,  2.9397,  2.8429,  2.8492,  2.7225,  3.3875,  3.4910,  3.4520 &
      ,  3.3608,  3.3036,  3.2345,  3.2999,  3.1487,  3.7409,  3.8392 &
      ,  3.7148,  3.6439,  3.6182,  3.1753,  3.5210,  3.4639,  3.4265 &
      ,  3.4075,  3.3828,  3.3474,  3.4071,  3.3754,  3.3646,  3.3308 &
      ,  3.4393,  3.2993,  3.8768,  3.9891,  3.8310,  3.7483,  3.3417 &
      ,  3.3019,  3.2250,  3.1832,  3.1578,  3.1564,  3.1224,  3.4620 &
      ,  2.9743,  2.8058,  3.4830,  3.3474,  3.6863,  3.3617,  3.1608 &
      ,  3.0069,  2.9640,  2.8427,  3.5885,  3.5219,  4.1314,  3.8120 &
      /)
      r0ab(1191:1260)=(/ &
         3.6015,  3.4502,  3.3498,  3.2777,  3.8635,  3.8232,  3.8486 &
      ,  3.7215,  3.6487,  3.4724,  3.5627,  3.5087,  3.4757,  3.4517 &
      ,  3.4423,  3.4139,  4.1028,  3.8388,  3.6745,  3.5562,  3.4806 &
      ,  3.4272,  4.0182,  3.9991,  4.0007,  3.9282,  3.7238,  3.6498 &
      ,  3.5605,  3.5211,  3.5009,  3.4859,  3.4785,  3.5621,  4.2623 &
      ,  3.0775,  2.8275,  4.0181,  3.3385,  3.5379,  3.5036,  3.2589 &
      ,  3.0804,  3.0094,  2.9003,  4.0869,  3.5088,  3.9105,  3.9833 &
      ,  3.7176,  3.5323,  3.4102,  3.3227,  4.2702,  4.0888,  3.7560 &
      ,  3.7687,  3.6681,  3.6405,  3.5569,  3.4990,  3.4659,  3.4433 &
      ,  3.4330,  3.4092,  3.8867,  4.0190,  3.7961,  3.6412,  3.5405 &
      /)
      r0ab(1261:1330)=(/ &
         3.4681,  4.3538,  4.2136,  3.9381,  3.8912,  3.9681,  3.7909 &
      ,  3.6774,  3.6262,  3.5999,  3.5823,  3.5727,  3.5419,  4.0245 &
      ,  4.1874,  3.0893,  2.7917,  3.7262,  3.3518,  3.4241,  3.5433 &
      ,  3.2773,  3.0890,  2.9775,  2.9010,  3.8048,  3.5362,  3.7746 &
      ,  3.7911,  3.7511,  3.5495,  3.4149,  3.3177,  4.0129,  3.8370 &
      ,  3.7739,  3.7125,  3.7152,  3.7701,  3.5813,  3.5187,  3.4835 &
      ,  3.4595,  3.4439,  3.4242,  3.7476,  3.8239,  3.8346,  3.6627 &
      ,  3.5479,  3.4639,  4.1026,  3.9733,  3.9292,  3.8667,  3.9513 &
      ,  3.8959,  3.7698,  3.7089,  3.6765,  3.6548,  3.6409,  3.5398 &
      ,  3.8759,  3.9804,  4.0150,  2.9091,  2.7638,  3.5066,  3.3377 &
      /)
      r0ab(1331:1400)=(/ &
         3.3481,  3.2633,  3.1810,  3.1428,  2.9872,  2.8837,  3.5929 &
      ,  3.5183,  3.6729,  3.6596,  3.6082,  3.5927,  3.4224,  3.2997 &
      ,  3.8190,  4.1865,  4.1114,  4.0540,  3.6325,  3.5697,  3.5561 &
      ,  3.5259,  3.4901,  3.4552,  3.4315,  3.4091,  3.6438,  3.6879 &
      ,  3.6832,  3.7043,  3.5557,  3.4466,  3.9203,  4.2919,  4.2196 &
      ,  4.1542,  3.7573,  3.7039,  3.6546,  3.6151,  3.5293,  3.4849 &
      ,  3.4552,  3.5192,  3.7673,  3.8359,  3.8525,  3.8901,  2.7806 &
      ,  2.7209,  3.3812,  3.4958,  3.2913,  3.1888,  3.0990,  3.0394 &
      ,  2.9789,  2.8582,  3.4716,  3.6883,  3.6105,  3.5704,  3.5059 &
      ,  3.4619,  3.4138,  3.2742,  3.7080,  3.9773,  3.9010,  3.8409 &
      /)
      r0ab(1401:1470)=(/ &
         3.7944,  3.4465,  3.7235,  3.6808,  3.6453,  3.6168,  3.5844 &
      ,  3.5576,  3.5772,  3.5959,  3.5768,  3.5678,  3.5486,  3.4228 &
      ,  3.8107,  4.0866,  4.0169,  3.9476,  3.6358,  3.5800,  3.5260 &
      ,  3.4838,  3.4501,  3.4204,  3.3553,  3.6487,  3.6973,  3.7398 &
      ,  3.7405,  3.7459,  3.7380,  2.6848,  2.6740,  3.2925,  3.3386 &
      ,  3.2473,  3.1284,  3.0301,  2.9531,  2.9602,  2.8272,  3.3830 &
      ,  3.5358,  3.5672,  3.5049,  3.4284,  3.3621,  3.3001,  3.2451 &
      ,  3.6209,  3.8299,  3.7543,  3.6920,  3.6436,  3.3598,  3.5701 &
      ,  3.5266,  3.4904,  3.4590,  3.4364,  3.4077,  3.5287,  3.5280 &
      ,  3.4969,  3.4650,  3.4304,  3.3963,  3.7229,  3.9402,  3.8753 &
      /)
      r0ab(1471:1540)=(/ &
         3.8035,  3.5499,  3.4913,  3.4319,  3.3873,  3.3520,  3.3209 &
      ,  3.2948,  3.5052,  3.6465,  3.6696,  3.6577,  3.6388,  3.6142 &
      ,  3.5889,  3.3968,  3.0122,  4.2241,  3.7887,  4.0049,  3.5384 &
      ,  3.2698,  3.1083,  2.9917,  2.9057,  4.3340,  3.9900,  4.6588 &
      ,  4.1278,  3.8125,  3.6189,  3.4851,  3.3859,  4.6531,  4.3134 &
      ,  4.2258,  4.1309,  4.0692,  4.0944,  3.9850,  3.9416,  3.9112 &
      ,  3.8873,  3.8736,  3.8473,  4.6027,  4.1538,  3.8994,  3.7419 &
      ,  3.6356,  3.5548,  4.8353,  4.5413,  4.3891,  4.3416,  4.3243 &
      ,  4.2753,  4.2053,  4.1790,  4.1685,  4.1585,  4.1536,  4.0579 &
      ,  4.1980,  4.4564,  4.2192,  4.0528,  3.9489,  3.8642,  5.0567 &
      /)
      r0ab(1541:1610)=(/ &
         3.0630,  3.3271,  4.0432,  4.0046,  4.1555,  3.7426,  3.5130 &
      ,  3.5174,  3.2884,  3.1378,  4.1894,  4.2321,  4.1725,  4.1833 &
      ,  3.8929,  4.0544,  3.8118,  3.6414,  4.6373,  4.6268,  4.4750 &
      ,  4.4134,  4.3458,  3.8582,  4.2583,  4.1898,  4.1562,  4.1191 &
      ,  4.1069,  4.0639,  4.1257,  4.1974,  3.9532,  4.1794,  3.9660 &
      ,  3.8130,  4.8160,  4.8272,  4.6294,  4.5840,  4.0770,  4.0088 &
      ,  3.9103,  3.8536,  3.8324,  3.7995,  3.7826,  4.2294,  4.3380 &
      ,  4.4352,  4.1933,  4.4580,  4.2554,  4.1072,  5.0454,  5.1814 &
      ,  3.0632,  3.2662,  3.6432,  3.8088,  3.7910,  3.7381,  3.5093 &
      ,  3.5155,  3.3047,  3.1681,  3.7871,  3.9924,  4.0637,  4.1382 &
      /)
      r0ab(1611:1680)=(/ &
         3.8591,  4.0164,  3.7878,  3.6316,  4.1741,  4.3166,  4.2395 &
      ,  4.1831,  4.1107,  3.5857,  4.0270,  3.9676,  3.9463,  3.9150 &
      ,  3.9021,  3.8708,  4.0240,  4.1551,  3.9108,  4.1337,  3.9289 &
      ,  3.7873,  4.3666,  4.5080,  4.4232,  4.3155,  3.8461,  3.8007 &
      ,  3.6991,  3.6447,  3.6308,  3.5959,  3.5749,  4.0359,  4.3124 &
      ,  4.3539,  4.1122,  4.3772,  4.1785,  4.0386,  4.7004,  4.8604 &
      ,  4.6261,  2.9455,  3.2470,  3.6108,  3.8522,  3.6625,  3.6598 &
      ,  3.4411,  3.4660,  3.2415,  3.0944,  3.7514,  4.0397,  3.9231 &
      ,  4.0561,  3.7860,  3.9845,  3.7454,  3.5802,  4.1366,  4.3581 &
      ,  4.2351,  4.2011,  4.1402,  3.5381,  4.0653,  4.0093,  3.9883 &
      /)
      r0ab(1681:1750)=(/ &
         3.9570,  3.9429,  3.9112,  3.8728,  4.0682,  3.8351,  4.1054 &
      ,  3.8928,  3.7445,  4.3415,  4.5497,  4.3833,  4.3122,  3.8051 &
      ,  3.7583,  3.6622,  3.6108,  3.5971,  3.5628,  3.5408,  4.0780 &
      ,  4.0727,  4.2836,  4.0553,  4.3647,  4.1622,  4.0178,  4.5802 &
      ,  4.9125,  4.5861,  4.6201,  2.9244,  3.2241,  3.5848,  3.8293 &
      ,  3.6395,  3.6400,  3.4204,  3.4499,  3.2253,  3.0779,  3.7257 &
      ,  4.0170,  3.9003,  4.0372,  3.7653,  3.9672,  3.7283,  3.5630 &
      ,  4.1092,  4.3347,  4.2117,  4.1793,  4.1179,  3.5139,  4.0426 &
      ,  3.9867,  3.9661,  3.9345,  3.9200,  3.8883,  3.8498,  4.0496 &
      ,  3.8145,  4.0881,  3.8756,  3.7271,  4.3128,  4.5242,  4.3578 &
      /)
      r0ab(1751:1820)=(/ &
         4.2870,  3.7796,  3.7318,  3.6364,  3.5854,  3.5726,  3.5378 &
      ,  3.5155,  4.0527,  4.0478,  4.2630,  4.0322,  4.3449,  4.1421 &
      ,  3.9975,  4.5499,  4.8825,  4.5601,  4.5950,  4.5702,  2.9046 &
      ,  3.2044,  3.5621,  3.8078,  3.6185,  3.6220,  3.4019,  3.4359 &
      ,  3.2110,  3.0635,  3.7037,  3.9958,  3.8792,  4.0194,  3.7460 &
      ,  3.9517,  3.7128,  3.5474,  4.0872,  4.3138,  4.1906,  4.1593 &
      ,  4.0973,  3.4919,  4.0216,  3.9657,  3.9454,  3.9134,  3.8986 &
      ,  3.8669,  3.8289,  4.0323,  3.7954,  4.0725,  3.8598,  3.7113 &
      ,  4.2896,  4.5021,  4.3325,  4.2645,  3.7571,  3.7083,  3.6136 &
      ,  3.5628,  3.5507,  3.5155,  3.4929,  4.0297,  4.0234,  4.2442 &
      /)
      r0ab(1821:1890)=(/ &
         4.0112,  4.3274,  4.1240,  3.9793,  4.5257,  4.8568,  4.5353 &
      ,  4.5733,  4.5485,  4.5271,  2.8878,  3.1890,  3.5412,  3.7908 &
      ,  3.5974,  3.6078,  3.3871,  3.4243,  3.1992,  3.0513,  3.6831 &
      ,  3.9784,  3.8579,  4.0049,  3.7304,  3.9392,  3.7002,  3.5347 &
      ,  4.0657,  4.2955,  4.1705,  4.1424,  4.0800,  3.4717,  4.0043 &
      ,  3.9485,  3.9286,  3.8965,  3.8815,  3.8500,  3.8073,  4.0180 &
      ,  3.7796,  4.0598,  3.8470,  3.6983,  4.2678,  4.4830,  4.3132 &
      ,  4.2444,  3.7370,  3.6876,  3.5935,  3.5428,  3.5314,  3.4958 &
      ,  3.4730,  4.0117,  4.0043,  4.2287,  3.9939,  4.3134,  4.1096 &
      ,  3.9646,  4.5032,  4.8356,  4.5156,  4.5544,  4.5297,  4.5083 &
      /)
      r0ab(1891:1960)=(/ &
         4.4896,  2.8709,  3.1737,  3.5199,  3.7734,  3.5802,  3.5934 &
      ,  3.3724,  3.4128,  3.1877,  3.0396,  3.6624,  3.9608,  3.8397 &
      ,  3.9893,  3.7145,  3.9266,  3.6877,  3.5222,  4.0448,  4.2771 &
      ,  4.1523,  4.1247,  4.0626,  3.4530,  3.9866,  3.9310,  3.9115 &
      ,  3.8792,  3.8641,  3.8326,  3.7892,  4.0025,  3.7636,  4.0471 &
      ,  3.8343,  3.6854,  4.2464,  4.4635,  4.2939,  4.2252,  3.7169 &
      ,  3.6675,  3.5739,  3.5235,  3.5126,  3.4768,  3.4537,  3.9932 &
      ,  3.9854,  4.2123,  3.9765,  4.2992,  4.0951,  3.9500,  4.4811 &
      ,  4.8135,  4.4959,  4.5351,  4.5105,  4.4891,  4.4705,  4.4515 &
      ,  2.8568,  3.1608,  3.5050,  3.7598,  3.5665,  3.5803,  3.3601 &
      /)
      r0ab(1961:2030)=(/ &
         3.4031,  3.1779,  3.0296,  3.6479,  3.9471,  3.8262,  3.9773 &
      ,  3.7015,  3.9162,  3.6771,  3.5115,  4.0306,  4.2634,  4.1385 &
      ,  4.1116,  4.0489,  3.4366,  3.9732,  3.9176,  3.8983,  3.8659 &
      ,  3.8507,  3.8191,  3.7757,  3.9907,  3.7506,  4.0365,  3.8235 &
      ,  3.6745,  4.2314,  4.4490,  4.2792,  4.2105,  3.7003,  3.6510 &
      ,  3.5578,  3.5075,  3.4971,  3.4609,  3.4377,  3.9788,  3.9712 &
      ,  4.1997,  3.9624,  4.2877,  4.0831,  3.9378,  4.4655,  4.7974 &
      ,  4.4813,  4.5209,  4.4964,  4.4750,  4.4565,  4.4375,  4.4234 &
      ,  2.6798,  3.0151,  3.2586,  3.5292,  3.5391,  3.4902,  3.2887 &
      ,  3.3322,  3.1228,  2.9888,  3.4012,  3.7145,  3.7830,  3.6665 &
      /)
      r0ab(2031:2100)=(/ &
         3.5898,  3.8077,  3.5810,  3.4265,  3.7726,  4.0307,  3.9763 &
      ,  3.8890,  3.8489,  3.2706,  3.7595,  3.6984,  3.6772,  3.6428 &
      ,  3.6243,  3.5951,  3.7497,  3.6775,  3.6364,  3.9203,  3.7157 &
      ,  3.5746,  3.9494,  4.2076,  4.1563,  4.0508,  3.5329,  3.4780 &
      ,  3.3731,  3.3126,  3.2846,  3.2426,  3.2135,  3.7491,  3.9006 &
      ,  3.8332,  3.8029,  4.1436,  3.9407,  3.7998,  4.1663,  4.5309 &
      ,  4.3481,  4.2911,  4.2671,  4.2415,  4.2230,  4.2047,  4.1908 &
      ,  4.1243,  2.5189,  2.9703,  3.3063,  3.6235,  3.4517,  3.3989 &
      ,  3.2107,  3.2434,  3.0094,  2.8580,  3.4253,  3.8157,  3.7258 &
      ,  3.6132,  3.5297,  3.7566,  3.5095,  3.3368,  3.7890,  4.1298 &
      /)
      r0ab(2101:2170)=(/ &
         4.0190,  3.9573,  3.9237,  3.2677,  3.8480,  3.8157,  3.7656 &
      ,  3.7317,  3.7126,  3.6814,  3.6793,  3.6218,  3.5788,  3.8763 &
      ,  3.6572,  3.5022,  3.9737,  4.3255,  4.1828,  4.1158,  3.5078 &
      ,  3.4595,  3.3600,  3.3088,  3.2575,  3.2164,  3.1856,  3.8522 &
      ,  3.8665,  3.8075,  3.7772,  4.1391,  3.9296,  3.7772,  4.2134 &
      ,  4.7308,  4.3787,  4.3894,  4.3649,  4.3441,  4.3257,  4.3073 &
      ,  4.2941,  4.1252,  4.2427,  3.0481,  2.9584,  3.6919,  3.5990 &
      ,  3.8881,  3.4209,  3.1606,  3.1938,  2.9975,  2.8646,  3.8138 &
      ,  3.7935,  3.7081,  3.9155,  3.5910,  3.4808,  3.4886,  3.3397 &
      ,  4.1336,  4.1122,  3.9888,  3.9543,  3.8917,  3.5894,  3.8131 &
      /)
      r0ab(2171:2240)=(/ &
         3.7635,  3.7419,  3.7071,  3.6880,  3.6574,  3.6546,  3.9375 &
      ,  3.6579,  3.5870,  3.6361,  3.5039,  4.3149,  4.2978,  4.1321 &
      ,  4.1298,  3.8164,  3.7680,  3.7154,  3.6858,  3.6709,  3.6666 &
      ,  3.6517,  3.8174,  3.8608,  4.1805,  3.9102,  3.8394,  3.8968 &
      ,  3.7673,  4.5274,  4.6682,  4.3344,  4.3639,  4.3384,  4.3162 &
      ,  4.2972,  4.2779,  4.2636,  4.0253,  4.1168,  4.1541,  2.8136 &
      ,  3.0951,  3.4635,  3.6875,  3.4987,  3.5183,  3.2937,  3.3580 &
      ,  3.1325,  2.9832,  3.6078,  3.8757,  3.7616,  3.9222,  3.6370 &
      ,  3.8647,  3.6256,  3.4595,  3.9874,  4.1938,  4.0679,  4.0430 &
      ,  3.9781,  3.3886,  3.9008,  3.8463,  3.8288,  3.7950,  3.7790 &
      /)
      r0ab(2241:2310)=(/ &
         3.7472,  3.7117,  3.9371,  3.6873,  3.9846,  3.7709,  3.6210 &
      ,  4.1812,  4.3750,  4.2044,  4.1340,  3.6459,  3.5929,  3.5036 &
      ,  3.4577,  3.4528,  3.4146,  3.3904,  3.9014,  3.9031,  4.1443 &
      ,  3.8961,  4.2295,  4.0227,  3.8763,  4.4086,  4.7097,  4.4064 &
      ,  4.4488,  4.4243,  4.4029,  4.3842,  4.3655,  4.3514,  4.1162 &
      ,  4.2205,  4.1953,  4.2794,  2.8032,  3.0805,  3.4519,  3.6700 &
      ,  3.4827,  3.5050,  3.2799,  3.3482,  3.1233,  2.9747,  3.5971 &
      ,  3.8586,  3.7461,  3.9100,  3.6228,  3.8535,  3.6147,  3.4490 &
      ,  3.9764,  4.1773,  4.0511,  4.0270,  3.9614,  3.3754,  3.8836 &
      ,  3.8291,  3.8121,  3.7780,  3.7619,  3.7300,  3.6965,  3.9253 &
      /)
      r0ab(2311:2380)=(/ &
         3.6734,  3.9733,  3.7597,  3.6099,  4.1683,  4.3572,  4.1862 &
      ,  4.1153,  3.6312,  3.5772,  3.4881,  3.4429,  3.4395,  3.4009 &
      ,  3.3766,  3.8827,  3.8868,  4.1316,  3.8807,  4.2164,  4.0092 &
      ,  3.8627,  4.3936,  4.6871,  4.3882,  4.4316,  4.4073,  4.3858 &
      ,  4.3672,  4.3485,  4.3344,  4.0984,  4.2036,  4.1791,  4.2622 &
      ,  4.2450,  2.7967,  3.0689,  3.4445,  3.6581,  3.4717,  3.4951 &
      ,  3.2694,  3.3397,  3.1147,  2.9661,  3.5898,  3.8468,  3.7358 &
      ,  3.9014,  3.6129,  3.8443,  3.6054,  3.4396,  3.9683,  4.1656 &
      ,  4.0394,  4.0158,  3.9498,  3.3677,  3.8718,  3.8164,  3.8005 &
      ,  3.7662,  3.7500,  3.7181,  3.6863,  3.9170,  3.6637,  3.9641 &
      /)
      r0ab(2381:2450)=(/ &
         3.7503,  3.6004,  4.1590,  4.3448,  4.1739,  4.1029,  3.6224 &
      ,  3.5677,  3.4785,  3.4314,  3.4313,  3.3923,  3.3680,  3.8698 &
      ,  3.8758,  4.1229,  3.8704,  4.2063,  3.9987,  3.8519,  4.3832 &
      ,  4.6728,  4.3759,  4.4195,  4.3952,  4.3737,  4.3551,  4.3364 &
      ,  4.3223,  4.0861,  4.1911,  4.1676,  4.2501,  4.2329,  4.2208 &
      ,  2.7897,  3.0636,  3.4344,  3.6480,  3.4626,  3.4892,  3.2626 &
      ,  3.3344,  3.1088,  2.9597,  3.5804,  3.8359,  3.7251,  3.8940 &
      ,  3.6047,  3.8375,  3.5990,  3.4329,  3.9597,  4.1542,  4.0278 &
      ,  4.0048,  3.9390,  3.3571,  3.8608,  3.8056,  3.7899,  3.7560 &
      ,  3.7400,  3.7081,  3.6758,  3.9095,  3.6552,  3.9572,  3.7436 &
      /)
      r0ab(2451:2520)=(/ &
         3.5933,  4.1508,  4.3337,  4.1624,  4.0916,  3.6126,  3.5582 &
      ,  3.4684,  3.4212,  3.4207,  3.3829,  3.3586,  3.8604,  3.8658 &
      ,  4.1156,  3.8620,  4.1994,  3.9917,  3.8446,  4.3750,  4.6617 &
      ,  4.3644,  4.4083,  4.3840,  4.3625,  4.3439,  4.3253,  4.3112 &
      ,  4.0745,  4.1807,  4.1578,  4.2390,  4.2218,  4.2097,  4.1986 &
      ,  2.8395,  3.0081,  3.3171,  3.4878,  3.5360,  3.5145,  3.2809 &
      ,  3.3307,  3.1260,  2.9940,  3.4741,  3.6675,  3.7832,  3.6787 &
      ,  3.6156,  3.8041,  3.5813,  3.4301,  3.8480,  3.9849,  3.9314 &
      ,  3.8405,  3.8029,  3.2962,  3.7104,  3.6515,  3.6378,  3.6020 &
      ,  3.5849,  3.5550,  3.7494,  3.6893,  3.6666,  3.9170,  3.7150 &
      /)
      r0ab(2521:2590)=(/ &
         3.5760,  4.0268,  4.1596,  4.1107,  3.9995,  3.5574,  3.5103 &
      ,  3.4163,  3.3655,  3.3677,  3.3243,  3.2975,  3.7071,  3.9047 &
      ,  3.8514,  3.8422,  3.8022,  3.9323,  3.7932,  4.2343,  4.4583 &
      ,  4.3115,  4.2457,  4.2213,  4.1945,  4.1756,  4.1569,  4.1424 &
      ,  4.0620,  4.0494,  3.9953,  4.0694,  4.0516,  4.0396,  4.0280 &
      ,  4.0130,  2.9007,  2.9674,  3.8174,  3.5856,  3.6486,  3.5339 &
      ,  3.2832,  3.3154,  3.1144,  2.9866,  3.9618,  3.8430,  3.9980 &
      ,  3.8134,  3.6652,  3.7985,  3.5756,  3.4207,  4.4061,  4.2817 &
      ,  4.1477,  4.0616,  3.9979,  3.6492,  3.8833,  3.8027,  3.7660 &
      ,  3.7183,  3.6954,  3.6525,  3.9669,  3.8371,  3.7325,  3.9160 &
      /)
      r0ab(2591:2660)=(/ &
         3.7156,  3.5714,  4.6036,  4.4620,  4.3092,  4.2122,  3.8478 &
      ,  3.7572,  3.6597,  3.5969,  3.5575,  3.5386,  3.5153,  3.7818 &
      ,  4.1335,  4.0153,  3.9177,  3.8603,  3.9365,  3.7906,  4.7936 &
      ,  4.7410,  4.5461,  4.5662,  4.5340,  4.5059,  4.4832,  4.4604 &
      ,  4.4429,  4.2346,  4.4204,  4.3119,  4.3450,  4.3193,  4.3035 &
      ,  4.2933,  4.1582,  4.2450,  2.8559,  2.9050,  3.8325,  3.5442 &
      ,  3.5077,  3.4905,  3.2396,  3.2720,  3.0726,  2.9467,  3.9644 &
      ,  3.8050,  3.8981,  3.7762,  3.6216,  3.7531,  3.5297,  3.3742 &
      ,  4.3814,  4.2818,  4.1026,  4.0294,  3.9640,  3.6208,  3.8464 &
      ,  3.7648,  3.7281,  3.6790,  3.6542,  3.6117,  3.8650,  3.8010 &
      /)
      r0ab(2661:2730)=(/ &
         3.6894,  3.8713,  3.6699,  3.5244,  4.5151,  4.4517,  4.2538 &
      ,  4.1483,  3.8641,  3.7244,  3.6243,  3.5589,  3.5172,  3.4973 &
      ,  3.4715,  3.7340,  4.0316,  3.9958,  3.8687,  3.8115,  3.8862 &
      ,  3.7379,  4.7091,  4.7156,  4.5199,  4.5542,  4.5230,  4.4959 &
      ,  4.4750,  4.4529,  4.4361,  4.1774,  4.3774,  4.2963,  4.3406 &
      ,  4.3159,  4.3006,  4.2910,  4.1008,  4.1568,  4.0980,  2.8110 &
      ,  2.8520,  3.7480,  3.5105,  3.4346,  3.3461,  3.1971,  3.2326 &
      ,  3.0329,  2.9070,  3.8823,  3.7928,  3.8264,  3.7006,  3.5797 &
      ,  3.7141,  3.4894,  3.3326,  4.3048,  4.2217,  4.0786,  3.9900 &
      ,  3.9357,  3.6331,  3.8333,  3.7317,  3.6957,  3.6460,  3.6197 &
      /)
      r0ab(2731:2800)=(/ &
         3.5779,  3.7909,  3.7257,  3.6476,  3.5729,  3.6304,  3.4834 &
      ,  4.4368,  4.3921,  4.2207,  4.1133,  3.8067,  3.7421,  3.6140 &
      ,  3.5491,  3.5077,  3.4887,  3.4623,  3.6956,  3.9568,  3.8976 &
      ,  3.8240,  3.7684,  3.8451,  3.6949,  4.6318,  4.6559,  4.4533 &
      ,  4.4956,  4.4641,  4.4366,  4.4155,  4.3936,  4.3764,  4.1302 &
      ,  4.3398,  4.2283,  4.2796,  4.2547,  4.2391,  4.2296,  4.0699 &
      ,  4.1083,  4.0319,  3.9855,  2.7676,  2.8078,  3.6725,  3.4804 &
      ,  3.3775,  3.2411,  3.1581,  3.1983,  2.9973,  2.8705,  3.8070 &
      ,  3.7392,  3.7668,  3.6263,  3.5402,  3.6807,  3.4545,  3.2962 &
      ,  4.2283,  4.1698,  4.0240,  3.9341,  3.8711,  3.5489,  3.7798 &
      /)
      r0ab(2801:2870)=(/ &
         3.7000,  3.6654,  3.6154,  3.5882,  3.5472,  3.7289,  3.6510 &
      ,  3.6078,  3.5355,  3.5963,  3.4480,  4.3587,  4.3390,  4.1635 &
      ,  4.0536,  3.7193,  3.6529,  3.5512,  3.4837,  3.4400,  3.4191 &
      ,  3.3891,  3.6622,  3.8934,  3.8235,  3.7823,  3.7292,  3.8106 &
      ,  3.6589,  4.5535,  4.6013,  4.3961,  4.4423,  4.4109,  4.3835 &
      ,  4.3625,  4.3407,  4.3237,  4.0863,  4.2835,  4.1675,  4.2272 &
      ,  4.2025,  4.1869,  4.1774,  4.0126,  4.0460,  3.9815,  3.9340 &
      ,  3.8955,  2.6912,  2.7604,  3.6037,  3.4194,  3.3094,  3.1710 &
      ,  3.0862,  3.1789,  2.9738,  2.8427,  3.7378,  3.6742,  3.6928 &
      ,  3.5512,  3.4614,  3.4087,  3.4201,  3.2607,  4.1527,  4.0977 &
      /)
      r0ab(2871:2940)=(/ &
         3.9523,  3.8628,  3.8002,  3.4759,  3.7102,  3.6466,  3.6106 &
      ,  3.5580,  3.5282,  3.4878,  3.6547,  3.5763,  3.5289,  3.5086 &
      ,  3.5593,  3.4099,  4.2788,  4.2624,  4.0873,  3.9770,  3.6407 &
      ,  3.5743,  3.5178,  3.4753,  3.3931,  3.3694,  3.3339,  3.6002 &
      ,  3.8164,  3.7478,  3.7028,  3.6952,  3.7669,  3.6137,  4.4698 &
      ,  4.5488,  4.3168,  4.3646,  4.3338,  4.3067,  4.2860,  4.2645 &
      ,  4.2478,  4.0067,  4.2349,  4.0958,  4.1543,  4.1302,  4.1141 &
      ,  4.1048,  3.9410,  3.9595,  3.8941,  3.8465,  3.8089,  3.7490 &
      ,  2.7895,  2.5849,  3.6484,  3.0162,  3.1267,  3.2125,  3.0043 &
      ,  2.9572,  2.8197,  2.7261,  3.7701,  3.2446,  3.5239,  3.4696 &
      /)
      r0ab(2941:3010)=(/ &
         3.4261,  3.3508,  3.1968,  3.0848,  4.1496,  3.6598,  3.5111 &
      ,  3.4199,  3.3809,  3.5382,  3.2572,  3.2100,  3.1917,  3.1519 &
      ,  3.1198,  3.1005,  3.5071,  3.5086,  3.5073,  3.4509,  3.3120 &
      ,  3.2082,  4.2611,  3.8117,  3.6988,  3.5646,  3.6925,  3.6295 &
      ,  3.5383,  3.4910,  3.4625,  3.4233,  3.4007,  3.2329,  3.6723 &
      ,  3.6845,  3.6876,  3.6197,  3.4799,  3.3737,  4.4341,  4.0525 &
      ,  3.9011,  3.8945,  3.8635,  3.8368,  3.8153,  3.7936,  3.7758 &
      ,  3.4944,  3.4873,  3.9040,  3.7110,  3.6922,  3.6799,  3.6724 &
      ,  3.5622,  3.6081,  3.5426,  3.4922,  3.4498,  3.3984,  3.4456 &
      ,  2.7522,  2.5524,  3.5742,  2.9508,  3.0751,  3.0158,  2.9644 &
      /)
      r0ab(3011:3080)=(/ &
         2.8338,  2.7891,  2.6933,  3.6926,  3.1814,  3.4528,  3.4186 &
      ,  3.3836,  3.2213,  3.1626,  3.0507,  4.0548,  3.5312,  3.4244 &
      ,  3.3409,  3.2810,  3.4782,  3.1905,  3.1494,  3.1221,  3.1128 &
      ,  3.0853,  3.0384,  3.4366,  3.4562,  3.4638,  3.3211,  3.2762 &
      ,  3.1730,  4.1632,  3.6825,  3.5822,  3.4870,  3.6325,  3.5740 &
      ,  3.4733,  3.4247,  3.3969,  3.3764,  3.3525,  3.1984,  3.5989 &
      ,  3.6299,  3.6433,  3.4937,  3.4417,  3.3365,  4.3304,  3.9242 &
      ,  3.7793,  3.7623,  3.7327,  3.7071,  3.6860,  3.6650,  3.6476 &
      ,  3.3849,  3.3534,  3.8216,  3.5870,  3.5695,  3.5584,  3.5508 &
      ,  3.4856,  3.5523,  3.4934,  3.4464,  3.4055,  3.3551,  3.3888 &
      /)
      r0ab(3081:3150)=(/ &
         3.3525,  2.7202,  2.5183,  3.4947,  2.8731,  3.0198,  3.1457 &
      ,  2.9276,  2.7826,  2.7574,  2.6606,  3.6090,  3.0581,  3.3747 &
      ,  3.3677,  3.3450,  3.1651,  3.1259,  3.0147,  3.9498,  3.3857 &
      ,  3.2917,  3.2154,  3.1604,  3.4174,  3.0735,  3.0342,  3.0096 &
      ,  3.0136,  2.9855,  2.9680,  3.3604,  3.4037,  3.4243,  3.2633 &
      ,  3.1810,  3.1351,  4.0557,  3.5368,  3.4526,  3.3699,  3.5707 &
      ,  3.5184,  3.4085,  3.3595,  3.3333,  3.3143,  3.3041,  3.1094 &
      ,  3.5193,  3.5745,  3.6025,  3.4338,  3.3448,  3.2952,  4.2158 &
      ,  3.7802,  3.6431,  3.6129,  3.5853,  3.5610,  3.5406,  3.5204 &
      ,  3.5036,  3.2679,  3.2162,  3.7068,  3.4483,  3.4323,  3.4221 &
      /)
      r0ab(3151:3220)=(/ &
         3.4138,  3.3652,  3.4576,  3.4053,  3.3618,  3.3224,  3.2711 &
      ,  3.3326,  3.2950,  3.2564,  2.5315,  2.6104,  3.2734,  3.2299 &
      ,  3.1090,  2.9942,  2.9159,  2.8324,  2.8350,  2.7216,  3.3994 &
      ,  3.4475,  3.4354,  3.3438,  3.2807,  3.2169,  3.2677,  3.1296 &
      ,  3.7493,  3.8075,  3.6846,  3.6104,  3.5577,  3.2052,  3.4803 &
      ,  3.4236,  3.3845,  3.3640,  3.3365,  3.3010,  3.3938,  3.3624 &
      ,  3.3440,  3.3132,  3.4035,  3.2754,  3.8701,  3.9523,  3.8018 &
      ,  3.7149,  3.3673,  3.3199,  3.2483,  3.2069,  3.1793,  3.1558 &
      ,  3.1395,  3.4097,  3.5410,  3.5228,  3.5116,  3.4921,  3.4781 &
      ,  3.4690,  4.0420,  4.1759,  4.0078,  4.0450,  4.0189,  3.9952 &
      /)
      r0ab(3221:3290)=(/ &
         3.9770,  3.9583,  3.9434,  3.7217,  3.8228,  3.7826,  3.8640 &
      ,  3.8446,  3.8314,  3.8225,  3.6817,  3.7068,  3.6555,  3.6159 &
      ,  3.5831,  3.5257,  3.2133,  3.1689,  3.1196,  3.3599,  2.9852 &
      ,  2.7881,  3.5284,  3.3493,  3.6958,  3.3642,  3.1568,  3.0055 &
      ,  2.9558,  2.8393,  3.6287,  3.5283,  4.1511,  3.8259,  3.6066 &
      ,  3.4527,  3.3480,  3.2713,  3.9037,  3.8361,  3.8579,  3.7311 &
      ,  3.6575,  3.5176,  3.5693,  3.5157,  3.4814,  3.4559,  3.4445 &
      ,  3.4160,  4.1231,  3.8543,  3.6816,  3.5602,  3.4798,  3.4208 &
      ,  4.0542,  4.0139,  4.0165,  3.9412,  3.7698,  3.6915,  3.6043 &
      ,  3.5639,  3.5416,  3.5247,  3.5153,  3.5654,  4.2862,  4.0437 &
      /)
      r0ab(3291:3360)=(/ &
         3.8871,  3.7741,  3.6985,  3.6413,  4.2345,  4.3663,  4.3257 &
      ,  4.0869,  4.0612,  4.0364,  4.0170,  3.9978,  3.9834,  3.9137 &
      ,  3.8825,  3.8758,  3.9143,  3.8976,  3.8864,  3.8768,  3.9190 &
      ,  4.1613,  4.0566,  3.9784,  3.9116,  3.8326,  3.7122,  3.6378 &
      ,  3.5576,  3.5457,  4.3127,  3.1160,  2.8482,  4.0739,  3.3599 &
      ,  3.5698,  3.5366,  3.2854,  3.1039,  2.9953,  2.9192,  4.1432 &
      ,  3.5320,  3.9478,  4.0231,  3.7509,  3.5604,  3.4340,  3.3426 &
      ,  4.3328,  3.8288,  3.7822,  3.7909,  3.6907,  3.6864,  3.5793 &
      ,  3.5221,  3.4883,  3.4649,  3.4514,  3.4301,  3.9256,  4.0596 &
      ,  3.8307,  3.6702,  3.5651,  3.4884,  4.4182,  4.2516,  3.9687 &
      /)
      r0ab(3361:3430)=(/ &
         3.9186,  3.9485,  3.8370,  3.7255,  3.6744,  3.6476,  3.6295 &
      ,  3.6193,  3.5659,  4.0663,  4.2309,  4.0183,  3.8680,  3.7672 &
      ,  3.6923,  4.5240,  4.4834,  4.1570,  4.3204,  4.2993,  4.2804 &
      ,  4.2647,  4.2481,  4.2354,  3.8626,  3.8448,  4.2267,  4.1799 &
      ,  4.1670,  3.8738,  3.8643,  3.8796,  4.0575,  4.0354,  3.9365 &
      ,  3.8611,  3.7847,  3.7388,  3.6826,  3.6251,  3.5492,  4.0889 &
      ,  4.2764,  3.1416,  2.8325,  3.7735,  3.3787,  3.4632,  3.5923 &
      ,  3.3214,  3.1285,  3.0147,  2.9366,  3.8527,  3.5602,  3.8131 &
      ,  3.8349,  3.7995,  3.5919,  3.4539,  3.3540,  4.0654,  3.8603 &
      ,  3.7972,  3.7358,  3.7392,  3.8157,  3.6055,  3.5438,  3.5089 &
      /)
      r0ab(3431:3500)=(/ &
         3.4853,  3.4698,  3.4508,  3.7882,  3.8682,  3.8837,  3.7055 &
      ,  3.5870,  3.5000,  4.1573,  4.0005,  3.9568,  3.8936,  3.9990 &
      ,  3.9433,  3.8172,  3.7566,  3.7246,  3.7033,  3.6900,  3.5697 &
      ,  3.9183,  4.0262,  4.0659,  3.8969,  3.7809,  3.6949,  4.2765 &
      ,  4.2312,  4.1401,  4.0815,  4.0580,  4.0369,  4.0194,  4.0017 &
      ,  3.9874,  3.8312,  3.8120,  3.9454,  3.9210,  3.9055,  3.8951 &
      ,  3.8866,  3.8689,  3.9603,  3.9109,  3.9122,  3.8233,  3.7438 &
      ,  3.7436,  3.6981,  3.6555,  3.5452,  3.9327,  4.0658,  4.1175 &
      ,  2.9664,  2.8209,  3.5547,  3.3796,  3.3985,  3.3164,  3.2364 &
      ,  3.1956,  3.0370,  2.9313,  3.6425,  3.5565,  3.7209,  3.7108 &
      /)
      r0ab(3501:3570)=(/ &
         3.6639,  3.6484,  3.4745,  3.3492,  3.8755,  4.2457,  3.7758 &
      ,  3.7161,  3.6693,  3.6155,  3.5941,  3.5643,  3.5292,  3.4950 &
      ,  3.4720,  3.4503,  3.6936,  3.7392,  3.7388,  3.7602,  3.6078 &
      ,  3.4960,  3.9800,  4.3518,  4.2802,  3.8580,  3.8056,  3.7527 &
      ,  3.7019,  3.6615,  3.5768,  3.5330,  3.5038,  3.5639,  3.8192 &
      ,  3.8883,  3.9092,  3.9478,  3.7995,  3.6896,  4.1165,  4.5232 &
      ,  4.4357,  4.4226,  4.4031,  4.3860,  4.3721,  4.3580,  4.3466 &
      ,  4.2036,  4.2037,  3.8867,  4.2895,  4.2766,  4.2662,  4.2598 &
      ,  3.8408,  3.9169,  3.8681,  3.8250,  3.7855,  3.7501,  3.6753 &
      ,  3.5499,  3.4872,  3.5401,  3.8288,  3.9217,  3.9538,  4.0054 &
      /)
      r0ab(3571:3640)=(/ &
         2.8388,  2.7890,  3.4329,  3.5593,  3.3488,  3.2486,  3.1615 &
      ,  3.1000,  3.0394,  2.9165,  3.5267,  3.7479,  3.6650,  3.6263 &
      ,  3.5658,  3.5224,  3.4762,  3.3342,  3.7738,  4.0333,  3.9568 &
      ,  3.8975,  3.8521,  3.4929,  3.7830,  3.7409,  3.7062,  3.6786 &
      ,  3.6471,  3.6208,  3.6337,  3.6519,  3.6363,  3.6278,  3.6110 &
      ,  3.4825,  3.8795,  4.1448,  4.0736,  4.0045,  3.6843,  3.6291 &
      ,  3.5741,  3.5312,  3.4974,  3.4472,  3.4034,  3.7131,  3.7557 &
      ,  3.7966,  3.8005,  3.8068,  3.8015,  3.6747,  4.0222,  4.3207 &
      ,  4.2347,  4.2191,  4.1990,  4.1811,  4.1666,  4.1521,  4.1401 &
      ,  3.9970,  3.9943,  3.9592,  4.0800,  4.0664,  4.0559,  4.0488 &
      /)
      r0ab(3641:3710)=(/ &
         3.9882,  4.0035,  3.9539,  3.9138,  3.8798,  3.8355,  3.5359 &
      ,  3.4954,  3.3962,  3.5339,  3.7595,  3.8250,  3.8408,  3.8600 &
      ,  3.8644,  2.7412,  2.7489,  3.3374,  3.3950,  3.3076,  3.1910 &
      ,  3.0961,  3.0175,  3.0280,  2.8929,  3.4328,  3.5883,  3.6227 &
      ,  3.5616,  3.4894,  3.4241,  3.3641,  3.3120,  3.6815,  3.8789 &
      ,  3.8031,  3.7413,  3.6939,  3.4010,  3.6225,  3.5797,  3.5443 &
      ,  3.5139,  3.4923,  3.4642,  3.5860,  3.5849,  3.5570,  3.5257 &
      ,  3.4936,  3.4628,  3.7874,  3.9916,  3.9249,  3.8530,  3.5932 &
      ,  3.5355,  3.4757,  3.4306,  3.3953,  3.3646,  3.3390,  3.5637 &
      ,  3.7053,  3.7266,  3.7177,  3.6996,  3.6775,  3.6558,  3.9331 &
      /)
      r0ab(3711:3780)=(/ &
         4.1655,  4.0879,  4.0681,  4.0479,  4.0299,  4.0152,  4.0006 &
      ,  3.9883,  3.8500,  3.8359,  3.8249,  3.9269,  3.9133,  3.9025 &
      ,  3.8948,  3.8422,  3.8509,  3.7990,  3.7570,  3.7219,  3.6762 &
      ,  3.4260,  3.3866,  3.3425,  3.5294,  3.7022,  3.7497,  3.7542 &
      ,  3.7494,  3.7370,  3.7216,  3.4155,  3.0522,  4.2541,  3.8218 &
      ,  4.0438,  3.5875,  3.3286,  3.1682,  3.0566,  2.9746,  4.3627 &
      ,  4.0249,  4.6947,  4.1718,  3.8639,  3.6735,  3.5435,  3.4479 &
      ,  4.6806,  4.3485,  4.2668,  4.1690,  4.1061,  4.1245,  4.0206 &
      ,  3.9765,  3.9458,  3.9217,  3.9075,  3.8813,  3.9947,  4.1989 &
      ,  3.9507,  3.7960,  3.6925,  3.6150,  4.8535,  4.5642,  4.4134 &
      /)
      r0ab(3781:3850)=(/ &
         4.3688,  4.3396,  4.2879,  4.2166,  4.1888,  4.1768,  4.1660 &
      ,  4.1608,  4.0745,  4.2289,  4.4863,  4.2513,  4.0897,  3.9876 &
      ,  3.9061,  5.0690,  5.0446,  4.6186,  4.6078,  4.5780,  4.5538 &
      ,  4.5319,  4.5101,  4.4945,  4.1912,  4.2315,  4.5534,  4.4373 &
      ,  4.4224,  4.4120,  4.4040,  4.2634,  4.7770,  4.6890,  4.6107 &
      ,  4.5331,  4.4496,  4.4082,  4.3095,  4.2023,  4.0501,  4.2595 &
      ,  4.5497,  4.3056,  4.1506,  4.0574,  3.9725,  5.0796,  3.0548 &
      ,  3.3206,  3.8132,  3.9720,  3.7675,  3.7351,  3.5167,  3.5274 &
      ,  3.3085,  3.1653,  3.9500,  4.1730,  4.0613,  4.1493,  3.8823 &
      ,  4.0537,  3.8200,  3.6582,  4.3422,  4.5111,  4.3795,  4.3362 &
      /)
      r0ab(3851:3920)=(/ &
         4.2751,  3.7103,  4.1973,  4.1385,  4.1129,  4.0800,  4.0647 &
      ,  4.0308,  4.0096,  4.1619,  3.9360,  4.1766,  3.9705,  3.8262 &
      ,  4.5348,  4.7025,  4.5268,  4.5076,  3.9562,  3.9065,  3.8119 &
      ,  3.7605,  3.7447,  3.7119,  3.6916,  4.1950,  4.2110,  4.3843 &
      ,  4.1631,  4.4427,  4.2463,  4.1054,  4.7693,  5.0649,  4.7365 &
      ,  4.7761,  4.7498,  4.7272,  4.7076,  4.6877,  4.6730,  4.4274 &
      ,  4.5473,  4.5169,  4.5975,  4.5793,  4.5667,  4.5559,  4.3804 &
      ,  4.6920,  4.6731,  4.6142,  4.5600,  4.4801,  4.0149,  3.8856 &
      ,  3.7407,  4.1545,  4.2253,  4.4229,  4.1923,  4.5022,  4.3059 &
      ,  4.1591,  4.7883,  4.9294,  3.3850,  3.4208,  3.7004,  3.8800 &
      /)
      r0ab(3921:3990)=(/ &
         3.9886,  3.9040,  3.6719,  3.6547,  3.4625,  3.3370,  3.8394 &
      ,  4.0335,  4.2373,  4.3023,  4.0306,  4.1408,  3.9297,  3.7857 &
      ,  4.1907,  4.3230,  4.2664,  4.2173,  4.1482,  3.6823,  4.0711 &
      ,  4.0180,  4.0017,  3.9747,  3.9634,  3.9383,  4.1993,  4.3205 &
      ,  4.0821,  4.2547,  4.0659,  3.9359,  4.3952,  4.5176,  4.3888 &
      ,  4.3607,  3.9583,  3.9280,  3.8390,  3.7971,  3.7955,  3.7674 &
      ,  3.7521,  4.1062,  4.3633,  4.2991,  4.2767,  4.4857,  4.3039 &
      ,  4.1762,  4.6197,  4.8654,  4.6633,  4.5878,  4.5640,  4.5422 &
      ,  4.5231,  4.5042,  4.4901,  4.3282,  4.3978,  4.3483,  4.4202 &
      ,  4.4039,  4.3926,  4.3807,  4.2649,  4.6135,  4.5605,  4.5232 &
      /)
      r0ab(3991:4060)=(/ &
         4.4676,  4.3948,  4.0989,  3.9864,  3.8596,  4.0942,  4.2720 &
      ,  4.3270,  4.3022,  4.5410,  4.3576,  4.2235,  4.6545,  4.7447 &
      ,  4.7043,  3.0942,  3.2075,  3.5152,  3.6659,  3.8289,  3.7459 &
      ,  3.5156,  3.5197,  3.3290,  3.2069,  3.6702,  3.8448,  4.0340 &
      ,  3.9509,  3.8585,  3.9894,  3.7787,  3.6365,  4.1425,  4.1618 &
      ,  4.0940,  4.0466,  3.9941,  3.5426,  3.8952,  3.8327,  3.8126 &
      ,  3.7796,  3.7635,  3.7356,  4.0047,  3.9655,  3.9116,  4.1010 &
      ,  3.9102,  3.7800,  4.2964,  4.3330,  4.2622,  4.2254,  3.8195 &
      ,  3.7560,  3.6513,  3.5941,  3.5810,  3.5420,  3.5178,  3.8861 &
      ,  4.1459,  4.1147,  4.0772,  4.3120,  4.1207,  3.9900,  4.4733 &
      /)
      r0ab(4061:4130)=(/ &
         4.6157,  4.4580,  4.4194,  4.3954,  4.3739,  4.3531,  4.3343 &
      ,  4.3196,  4.2140,  4.2339,  4.1738,  4.2458,  4.2278,  4.2158 &
      ,  4.2039,  4.1658,  4.3595,  4.2857,  4.2444,  4.1855,  4.1122 &
      ,  3.7839,  3.6879,  3.5816,  3.8633,  4.1585,  4.1402,  4.1036 &
      ,  4.3694,  4.1735,  4.0368,  4.5095,  4.5538,  4.5240,  4.4252 &
      ,  3.0187,  3.1918,  3.5127,  3.6875,  3.7404,  3.6943,  3.4702 &
      ,  3.4888,  3.2914,  3.1643,  3.6669,  3.8724,  3.9940,  4.0816 &
      ,  3.8054,  3.9661,  3.7492,  3.6024,  4.0428,  4.1951,  4.1466 &
      ,  4.0515,  4.0075,  3.5020,  3.9158,  3.8546,  3.8342,  3.8008 &
      ,  3.7845,  3.7549,  3.9602,  3.8872,  3.8564,  4.0793,  3.8835 &
      /)
      r0ab(4131:4200)=(/ &
         3.7495,  4.2213,  4.3704,  4.3300,  4.2121,  3.7643,  3.7130 &
      ,  3.6144,  3.5599,  3.5474,  3.5093,  3.4853,  3.9075,  4.1115 &
      ,  4.0473,  4.0318,  4.2999,  4.1050,  3.9710,  4.4320,  4.6706 &
      ,  4.5273,  4.4581,  4.4332,  4.4064,  4.3873,  4.3684,  4.3537 &
      ,  4.2728,  4.2549,  4.2032,  4.2794,  4.2613,  4.2491,  4.2375 &
      ,  4.2322,  4.3665,  4.3061,  4.2714,  4.2155,  4.1416,  3.7660 &
      ,  3.6628,  3.5476,  3.8790,  4.1233,  4.0738,  4.0575,  4.3575 &
      ,  4.1586,  4.0183,  4.4593,  4.5927,  4.4865,  4.3813,  4.4594 &
      ,  2.9875,  3.1674,  3.4971,  3.6715,  3.7114,  3.6692,  3.4446 &
      ,  3.4676,  3.2685,  3.1405,  3.6546,  3.8579,  3.9637,  4.0581 &
      /)
      r0ab(4201:4270)=(/ &
         3.7796,  3.9463,  3.7275,  3.5792,  4.0295,  4.1824,  4.1247 &
      ,  4.0357,  3.9926,  3.4827,  3.9007,  3.8392,  3.8191,  3.7851 &
      ,  3.7687,  3.7387,  3.9290,  3.8606,  3.8306,  4.0601,  3.8625 &
      ,  3.7269,  4.2062,  4.3566,  4.3022,  4.1929,  3.7401,  3.6888 &
      ,  3.5900,  3.5350,  3.5226,  3.4838,  3.4594,  3.8888,  4.0813 &
      ,  4.0209,  4.0059,  4.2810,  4.0843,  3.9486,  4.4162,  4.6542 &
      ,  4.5005,  4.4444,  4.4196,  4.3933,  4.3741,  4.3552,  4.3406 &
      ,  4.2484,  4.2413,  4.1907,  4.2656,  4.2474,  4.2352,  4.2236 &
      ,  4.2068,  4.3410,  4.2817,  4.2479,  4.1921,  4.1182,  3.7346 &
      ,  3.6314,  3.5168,  3.8582,  4.0927,  4.0469,  4.0313,  4.3391 &
      /)
      r0ab(4271:4340)=(/ &
         4.1381,  3.9962,  4.4429,  4.5787,  4.4731,  4.3588,  4.4270 &
      ,  4.3957,  2.9659,  3.1442,  3.4795,  3.6503,  3.6814,  3.6476 &
      ,  3.4222,  3.4491,  3.2494,  3.1209,  3.6324,  3.8375,  3.9397 &
      ,  3.8311,  3.7581,  3.9274,  3.7085,  3.5598,  4.0080,  4.1641 &
      ,  4.1057,  4.0158,  3.9726,  3.4667,  3.8802,  3.8188,  3.7989 &
      ,  3.7644,  3.7474,  3.7173,  3.9049,  3.8424,  3.8095,  4.0412 &
      ,  3.8436,  3.7077,  4.1837,  4.3366,  4.2816,  4.1686,  3.7293 &
      ,  3.6709,  3.5700,  3.5153,  3.5039,  3.4684,  3.4437,  3.8663 &
      ,  4.0575,  4.0020,  3.9842,  4.2612,  4.0643,  3.9285,  4.3928 &
      ,  4.6308,  4.4799,  4.4244,  4.3996,  4.3737,  4.3547,  4.3358 &
      /)
      r0ab(4341:4410)=(/ &
         4.3212,  4.2275,  4.2216,  4.1676,  4.2465,  4.2283,  4.2161 &
      ,  4.2045,  4.1841,  4.3135,  4.2562,  4.2226,  4.1667,  4.0932 &
      ,  3.7134,  3.6109,  3.4962,  3.8352,  4.0688,  4.0281,  4.0099 &
      ,  4.3199,  4.1188,  3.9768,  4.4192,  4.5577,  4.4516,  4.3365 &
      ,  4.4058,  4.3745,  4.3539,  2.8763,  3.1294,  3.5598,  3.7465 &
      ,  3.5659,  3.5816,  3.3599,  3.4024,  3.1877,  3.0484,  3.7009 &
      ,  3.9451,  3.8465,  3.9873,  3.7079,  3.9083,  3.6756,  3.5150 &
      ,  4.0829,  4.2780,  4.1511,  4.1260,  4.0571,  3.4865,  3.9744 &
      ,  3.9150,  3.8930,  3.8578,  3.8402,  3.8073,  3.7977,  4.0036 &
      ,  3.7604,  4.0288,  3.8210,  3.6757,  4.2646,  4.4558,  4.2862 &
      /)
      r0ab(4411:4465)=(/ &
         4.2122,  3.7088,  3.6729,  3.5800,  3.5276,  3.5165,  3.4783 &
      ,  3.4539,  3.9553,  3.9818,  4.2040,  3.9604,  4.2718,  4.0689 &
      ,  3.9253,  4.4869,  4.7792,  4.4918,  4.5342,  4.5090,  4.4868 &
      ,  4.4680,  4.4486,  4.4341,  4.2023,  4.3122,  4.2710,  4.3587 &
      ,  4.3407,  4.3281,  4.3174,  4.1499,  4.3940,  4.3895,  4.3260 &
      ,  4.2725,  4.1961,  3.7361,  3.6193,  3.4916,  3.9115,  3.9914 &
      ,  3.9809,  3.9866,  4.3329,  4.1276,  3.9782,  4.5097,  4.6769 &
      ,  4.5158,  4.3291,  4.3609,  4.3462,  4.3265,  4.4341 &
      /)

      k=0
      do i=1,max_elem
         do j=1,i
            k=k+1
            r(i,j)=r0ab(k)/autoang
            r(j,i)=r0ab(k)/autoang
         enddo
      enddo

      end subroutine setr0ab
