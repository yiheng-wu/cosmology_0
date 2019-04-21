program p2d
implicit none

character(4)::simucode='6620'
character(4)::step_string='2448'
character(100)::indir='/data/s1/simu/Jing6620/'

real(4),allocatable::x(:,:)
real(4),allocatable::v(:,:)

integer(8)::np,ips
real(4)::ztp,omegat,lambdat,rLbox,xscale,vscale
real(4)::alpha=1.,dp=0.0288,nstepf=5000.,scalep,scalepf,vfact2
real(4)::a,H,omega0=0.268
integer(4),parameter::node=4

logical(4)::rsd=.false.

integer(4)::i,j,k
integer(8)::pid
real(4)::t1,t2

call memo(1)
call readposition
if (rsd) then 
  call readvelocity
  call rsdshift
endif
call memo(0)

contains

subroutine memo(command)
implicit none
integer(4)::command
character(100)::fname
if (command.eq.1) then
  fname=trim(indir)//'pos'//simucode//'.'//step_string//'.'//'01'
  write(*,*) 'reading:',trim(fname) 
  open(31,file=fname,status='old',form='unformatted')
  read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
  close(31)
  write(*,*) 'particle number:',np
  write(*,*) 'ips:',ips
  write(*,*) 'ztp:',ztp
  write(*,*) 'a:',1./(1.+ztp)
  write(*,*) 'omegat:',omegat
  write(*,*) 'lambdat:',lambdat
  write(*,*) 'boxsize:',rLbox
  write(*,*) 'xscale:',xscale
  write(*,*) 'vscale:',vscale
  write(*,*) 'memo for x:',(3*float(np))*(4)/(1024.**3),'G'
  allocate(x(3,np))
  if (rsd) then
    write(*,*) 'memo for v:',(3*float(np))*(4)/(1024.**3),'G'
    allocate(v(3,np))
  endif
elseif (command.eq.0) then
  deallocate(x)
  if (rsd) deallocate(v)
  write(*,*) 'WORK DONE!'
  write(*,*) ' '
endif
endsubroutine memo

subroutine readposition
implicit none
integer(4)::fid,uid
character(2)::f_string
character(100)::fname
real(4)::tt1,tt2
real(8)::li
call cpu_time(t1)
do fid=1,node
  write(f_string,'(i2.2)') fid
  fname=trim(indir)//'pos'//simucode//'.'//step_string//'.'//f_string
  write(*,*) 'reading:',trim(fname)
  uid=30+fid
  open(uid,file=fname,status='old',form='unformatted')
  if (fid.eq.1) read(uid) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
  write(*,*) 'reading particle from',(fid-1)*np/node+1,'to',fid*np/node
  call cpu_time(tt1)
  read(uid) x(1:3,(fid-1)*np/node+1:fid*np/node)
  call cpu_time(tt2)
  write(*,*) 'time consumed for reading command:',tt2-tt1,'seconds'
  close(uid)
enddo
call cpu_time(t2)
write(*,*) 'time consumed for reading position:',t2-t1,'seconds'
endsubroutine readposition

subroutine readvelocity
implicit none
integer(4)::fid,uid
character(2)::f_string
character(100)::fname
call cpu_time(t1)
do fid=1,node
  write(f_string,'(i2.2)') fid
  fname=trim(indir)//'vel'//simucode//'.'//step_string//'.'//f_string
  write(*,*) 'reading:',trim(fname)
  uid=30+fid
  open(uid,file=fname,status='old',form='unformatted')
  if (fid.eq.1) read(uid) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
  write(*,*) 'reading particle from',(fid-1)*np/node+1,'to',fid*np/node
  read(uid) v(1:3,(fid-1)*np/node+1:fid*np/node)
  close(uid)
enddo
call cpu_time(t2)
write(*,*) 'time consumed for reading velocity:',t2-t1,'seconds'
endsubroutine readvelocity

subroutine rsdshift
implicit none
call cpu_time(t1)
!H=100.*sqrt(omega0*(1.+ztp)**3+1-omega0)
write(*,*) 'processing RSD...'
write(*,*) 'redshift:',ztp
write(*,*) 'scale factor:',1./(1.+ztp)
!write(*,*) 'Hubble parameter:',H
scalep=ips*dp+1.
write(*,*) 'scalep:',scalep
vfact2=alpha*scalep
write(*,*) 'vfact2:',vfact2
x(3,:)=x(3,:)+v(3,:)*vfact2
where(x(3,:).gt.1.) x(3,:)=x(3,:)-1.
where(x(3,:).le.0.) x(3,:)=x(3,:)+1.
call cpu_time(t2)
write(*,*) 'time consumed for rsdshift:',t2-t1,'seconds'
endsubroutine rsdshift

endprogram p2d
