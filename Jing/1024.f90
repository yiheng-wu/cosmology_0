program jing2gadget
implicit none

include './inc/inc6519.f90'

real(4),allocatable::x(:,:)
real(4),allocatable::v(:,:)
integer(4),allocatable::id(:)

! Jing's variables
integer(4)::np,ips
real(4)::ztp,omegat,lambdat,rLbox,xscale,vscale
real(4)::scalep,scalepf,vfact2
real(4)::a,H

! conversion variables
real(4),parameter::rhocrit=2.77e11 ! (Msun/h) / (Mpc/h)^3
real(4)::massp

! multifile variables
integer(4)::npfile
integer(4)::list(2,nfile)
integer(4)::nplist(nfile)

! Gadget's variables
integer(4)::npart(6)
real(8)::massarr(6)
real(8)::Ga
real(8)::redshift
integer(4)::flag_sfr,flag_feedback
integer(4)::nall(6)
integer(4)::flag_cooling
integer(4)::numfiles
real(8)::boxsize
real(8)::omegam,omegal,hubblepara
integer(4)::flag_age,flag_metals
integer(4)::nallhw(6)
integer(4)::flag_entr_ics
real(4)::unused(15)

! miscs
character(128)::str
character(4)::step_string ! Jing's step
character(3)::gstep ! Gadget's step
character(3)::fstr ! Gadget's file id
integer(4)::step
character(512)::filename
integer(4)::i,j,k

! debug
logical(4)::debug=.false.

call getarg(1,str)
read(str,*) step
write(step_string,'(i4.4)') step
write(*,*) "Jing's step: ", step_string

call getarg(2,str)
read(str,*) step
write(gstep,'(i3.3)') step
write(*,*) "Gadget's step: ", gstep


filename=trim(indir)//'pos'//simucode//'.'//step_string
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',form='unformatted',convert='big_endian')
read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
close(31)
write(*,*) '========== JING HEAD =========='
write(*,*) 'particle number:',np
write(*,*) 'ips:',ips
write(*,*) 'ztp:',ztp
write(*,*) 'omegat:',omegat
write(*,*) 'lambdat:',lambdat
write(*,*) 'boxsize:',rLbox
write(*,*) 'xscale:',xscale
write(*,*) 'vscale:',vscale
write(*,*) 'a:',(ips*dp+1)/(nstepf*dp+1)
write(*,*) 'z:',(nstepf*dp+1)/(ips*dp+1)-1
write(*,*) '========== JING HEAD =========='

allocate(x(3,np))
allocate(v(3,np))
allocate(id(np))

filename=trim(indir)//'pos'//simucode//'.'//step_string
open(31,file=filename,status='old',form='unformatted',convert='big_endian')
read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
read(31) x(1,:),x(2,:),x(3,:)
close(31)
if (debug) then
write(*,*) 'min,max:',minval(x(1,:)),maxval(x(1,:))
write(*,*) 'min,max:',minval(x(2,:)),maxval(x(2,:))
write(*,*) 'min,max:',minval(x(3,:)),maxval(x(3,:))
endif

filename=trim(indir)//'vel'//simucode//'.'//step_string
open(31,file=filename,status='old',form='unformatted',convert='big_endian')
read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
read(31) v(1,:),v(2,:),v(3,:)
close(31)
if (debug) then
write(*,*) 'min,max:',minval(v(1,:)),maxval(v(1,:))
write(*,*) 'min,max:',minval(v(2,:)),maxval(v(2,:))
write(*,*) 'min,max:',minval(v(3,:)),maxval(v(3,:))
endif

filename=trim(indir)//'id'//simucode//'.'//step_string
open(31,file=filename,status='old',form='unformatted',convert='big_endian')
read(31) id(1:np)
close(31)
if (debug) write(*,*) 'min,max:',minval(id(:)),maxval(id(:))

! convert position to gadget
x=x*rLbox*1000.  ! kpc/h

! convert velocity to gadget
scalepf=nstepf*dp+1.
scalep=ips*dp+1.
a=scalep/scalepf
H=100.*sqrt(omega0/a**3.+1-omega0)
vfact2=alpha*scalep
!!v=vfact2*v*rLbox*H*a !! convert Jing's velocity to physical
!!v=vfact2*v*rLbox*H   !! convert Jing's velocity to comoving
v=vfact2*v*rLbox*H*sqrt(a)  !! convert Jing's velocity to Gadget unit

! convert id to gadget
id=id-1

! convert mass
massp=rhocrit*rLbox**3*omega0/np
write(*,*) 'particle mass:',massp,'Msun/h'
massp=massp/1.e10 

! multifile setting
npfile=np/nfile
do i=1,nfile-1
  list(1,i)=(i-1)*npfile+1
  list(2,i)=list(1,i)+npfile-1
enddo
list(1,nfile)=(nfile-1)*npfile+1
list(2,nfile)=np
do i=1,nfile
  nplist(i)=list(2,i)-list(1,i)+1
enddo
write(*,*) 'multifile information'
do i=1,nfile
  write(*,*) nplist(i),list(1,i),list(2,i)
enddo

! define Gadget head
massarr(1:6)=[0.,massp,0.,0.,0.,0.]
Ga=(ips*dp+1)/(nstepf*dp+1)
redshift=(nstepf*dp+1)/(ips*dp+1)-1
flag_sfr=0
flag_feedback=0
nall(1:6)=[0,np,0,0,0,0]
flag_cooling=0
numfiles=nfile
boxsize=rLbox*1000.
omegam=omega0
omegal=1-omega0
hubblepara=Jinghubble
flag_age=0
flag_metals=0
nallhw(:)=0
flag_entr_ics=0
unused(:)=0

write(*,*) '========== GADGET HEAD =========='
write(*,*) 'scale factor=',a
write(*,*) 'redshift=',redshift
write(*,*) 'flag_sfr=',flag_sfr
write(*,*) 'flag_feedback=',flag_feedback
write(*,*) 'nall=',nall
write(*,*) 'flag_cooling=',flag_cooling
write(*,*) 'numfiles=',numfiles
write(*,*) 'boxsize=',boxsize
write(*,*) 'omega0',omega0
write(*,*) 'omegal=',omegal
write(*,*) 'hubblepara=',hubblepara
write(*,*) 'flag_age=',flag_age
write(*,*) 'flag_metals=',flag_metals
write(*,*) 'nallhw=',nallhw
write(*,*) 'flag_entr_ics=',flag_entr_ics
write(*,*) '========== GADGET HEAD =========='

! write Gadget

do i=1,numfiles

  npart(1:6)=[0,nplist(i),0,0,0,0]

  write(fstr,'(i3)') i-1
  fstr=adjustl(fstr)
  filename=trim(outdir)//'snapshot_'//gstep//'.'//trim(fstr)
  write(*,*) 'writing: ',trim(filename)
  open(61,file=filename,form='unformatted',status='replace')
  write(61) npart, massarr, Ga, redshift, flag_sfr, flag_feedback, nall, flag_cooling, &
           numfiles, boxsize, omegam, omegal, hubblepara, flag_age, flag_metals, &
           nallhw, flag_entr_ics, unused
  write(61) x(1:3,list(1,i):list(2,i))
  write(61) v(1:3,list(1,i):list(2,i))
  write(61) id(list(1,i):list(2,i))
  close(61)

enddo

deallocate(x,v,id)

endprogram jing2gadget
