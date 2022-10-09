! Mean square displacement calculation for interfacial water molecules
! Water molecules near 6 angstorm from protein backbone atom are considerd only
! Written by Abhik Ghosh Moulick, Department of Physics of Complex Systems
! SNBNCBS
	
integer,parameter:: atomno = 369   ! change,check last atomno of complex from pdb file
integer,parameter:: solno =  13118 ! change, check lst solno from pdb file
integer,parameter:: maxframe =  1005 ! change, maximum frame to consider
real :: xc(maxframe,1:atomno,3),xw(maxframe,atomno+1:solno,3) 
real :: num(atomno+1:solno), distmatrix(maxframe,atomno+1:solno) 
integer :: ichain,atnum,gap,p
character :: line*80
character :: atnam*4, junk*1
integer :: resnum,j,i,steps,nbin,resinum,resn
real :: distance(0:100),counter,dist1,dist2
real :: dist, dist_x, dist_y, dist_z, dist_gap
real :: dist_a, dist_b, dist_c, Pt
integer :: delay, origin
data delt/0.1/

open(10, file='Input.pdb') ! Change file name
open(20, file='msd-water-around-protein.dat')! Change file name
	
! Read pdb file
ichain=1
do while (ichain.ge.1)	
read (10,'(a)') line
if (line(1:4).eq.'ATOM') then
read (line,1) resnum,atnam,junk,resn,x,y,z
!write (22,1) resnum,atnam,junk,resn,x,y,z
1 format (6x,i5,2x,a4,3x,a2,1x,i4,3x,3f8.3)
if(atnam.ne.'OW') then
	xc(ichain,resnum,1) = x
	xc(ichain,resnum,2) = y
	xc(ichain,resnum,3) = z
else if(atnam.eq.'OW') then
	xw(ichain,resnum,1) = x
	xw(ichain,resnum,2) = y
	xw(ichain,resnum,3) = z
end if 
end if
				
if (line(1:4).eq.'TER ') THEN
ichain = ichain + 1
maxchain = ichain
!write(*,*) ichain
if (maxchain .eq. 1002) exit
end if
end do
 close(10)

! MSD calculation      	
! Initialize array to control interfacial water number
! Any water molecule which is at 6 angstorm distance from backbone atom
! of any residue should not consider for further calculations

steps = maxchain-1
origin = 950  ! Change
delay = 50    ! Change
do ichain = 1,steps
	do j = atomno+1, solno
	Pt = 0.0
		do i = 1, atomno
			dist_x = xc(ichain,i,1)-xw(ichain,j,1)
      			dist_y = xc(ichain,i,2)-xw(ichain,j,2)
      			dist_z = xc(ichain,i,3)-xw(ichain,j,3)
      			dist = sqrt(dist_x**2+dist_y**2+dist_z**2)
      			if (dist .le. 6.0) then
      				Pt = Pt + 1.0
      			end if
      		end do
      	if (Pt .gt. 0.0) then
		distmatrix(ichain,j) = 1.0
	else 
		distmatrix(ichain,j) = 0.0
	end if
	end do
end do
do gap = 1,delay
distance(gap) = 0.0
	do ichain = 1, origin 
		counter = 0.0
		dist2 = 0.0
		do j = atomno+1, solno
			dist1 = 0.0
			do p = 1,gap
				dist1 = dist1 + distmatrix(ichain+p,j)
			end do
			if (dist1 .eq. gap) then
				counter = counter + 1.0
				dist_a = xw(ichain,j,1)-xw(ichain+gap,j,1)
      				dist_b = xw(ichain,j,2)-xw(ichain+gap,j,2)
      				dist_c = xw(ichain,j,3)-xw(ichain+gap,j,3)
      				dist_gap = (dist_a**2)+(dist_b**2)+(dist_c**2)
      				dist2 = dist2 + dist_gap
      			end if
      		end do
	if (dist2 .eq. 0.0 .or. counter .eq. 0.0) then
		distance(gap) = distance(gap)
	else
		distance(gap) = distance(gap)+(dist2/counter)
	end if
	end do
	distance(gap) = distance(gap)/(float(origin))
	
write(20,35)float(gap),distance(gap)
35 	format (f10.4,'	',f16.5)
end do
end
