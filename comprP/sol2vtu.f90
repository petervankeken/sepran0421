! sol2vtu
! output solution to vtu file
! PvK October 30 2017
subroutine sol2vtu(ichoice,isol,basename,idegfd,namedof)
use sepmodulecomio
use sepmodulecpack
use sepmodulemain
use sepmodulekmesh
use sepmodulemeshinf
use sepmodulemesh
use sepmodulekprob
use sepmoduleprobinf
use coeff
use control
implicit none
interface 
  subroutine sol2vtu01(ichoice,npoint,npelm,nelem,nphys,ndim,idegfd,usol,nusol,coor,kmeshc,kprobp_copy,filename,  &
        & namedof,num_dofs_is_constant)
    integer, intent(in) :: ichoice,npoint,npelm,nelem,nphys,ndim,idegfd,kmeshc(:),nusol
    integer, intent(inout) :: kprobp_copy(:,:)
    real(kind=8), intent(in) :: usol(:),coor(:,:)
    character(len=*), intent(in) :: filename,namedof(:)
    logical, intent(in) :: num_dofs_is_constant
  end subroutine sol2vtu01
end interface
integer,intent(in) :: ichoice ! 1=single scalar; 2=(u,v,p) and output (u,v) as vector, p as scalar
integer,intent(in) :: isol  ! standard sepran array
character(len=*),intent(in) :: basename,namedof(:)
character(len=120) :: filename
integer :: idegfd  ! 0=output all dof's, otherwise output only idegfdth dof
integer :: idof,inode_here,ip
integer :: i,j,kproba_min,kproba_max,allocate_status,kproba_here
character(len=4) :: ext
logical :: lu_used,num_dofs_is_constant
integer,dimension(:,:),allocatable :: kprobp_copy

call sepactsolbf1(isol)

if (print_node) write(irefwr,*) 'output isol to file with basename ',basename(1:len_trim(basename))

num_dofs_is_constant=.true.
if (ichoice==2 .and. npelm==6) then
   ! Stokes solution on P2 element 
   if (itype_stokes==903) then
     num_dofs_is_constant=.false.
     if (kprobloc(14)==0) then
        if (print_node) write(irefwr,*) 'PERROR(sol2vtu): expecting P2 903 but kprob(14)==0: ',kprobloc(14)
        call instop
     endif
     if (indprpi==0) then
        if (print_node) write(irefwr,*) 'PERROR(sol2vtu): expecting P2 903 but indprpi==0: ',indprpi
        call instop
     endif
   endif
endif
   
if (idegfd>nphys) then
   ! one shouldn't ask to output a dof beyond nphys
   if (print_node) write(irefwr,*) 'PERROR(sol2vtu): idegfd > nphys'
   if (print_node) write(irefwr,*) 'idegfd = ',idegfd,' nphys = ',nphys
   call instop
endif

if (.not.num_dofs_is_constant) then
   ! We want to interpolate pressure in the points on the sides of the elements
   ! and locally modify kprobp
   allocate(kprobp_copy(npoint,nphys),stat=allocate_status)
   if (allocate_status /= 0) then
      if (print_node) write(irefwr,*) 'PERROR(isol2vtu): problem allocating kprobp_copy'
      call instop
   endif
   kprobp_copy=kprobpi
endif

!write(irefwr,*) 'kmeshc: '
!do i=1,nelem
!   ip=(i-1)*npelm
!   write(irefwr,'(''element '',i5,'' has nodals: '',6i5)') i,(kmeshc(ip+j),j=1,npelm)
!enddo

!if (indprpi>0) then
!  write(irefwr,*) 'in first triangle: '
!  do i=1,npelm
!     ip=kmeshc(i)
!     write(irefwr,'(''nodal point '',i5,'' has kprobp: '',3i5)') ip,kprobpi(ip,1),kprobpi(ip,2),kprobpi(ip,3)
!  enddo
!endif
     

if (parallel) then
   ! Sepran parallelism with domain decomposition. Report info and errors for all nodes

   write(irefwr,*) 'output isol to file with basename ',basename(1:len_trim(basename))
   write(irefwr,*) 'in parallel for node : ',inodenr

   ! check on limits of fixed strings
   if (iacnodes>999) then
      write(irefwr,*) 'PERROR(sol2vtu): number of processes >= 1000: ',iacnodes
      write(irefwr,*) 'change the number of digits appended to basename'
      call instop
   endif
   if (len_trim(basename) + 8 > 120) then
      write(irefwr,*) 'PERROR(sol2vtu) :: length of basename is too much'
      write(irefwr,*) 'length should be less than 116 but is: ',len_trim(basename)
   endif

   ! set file name to basename + inodenr + '.vtu'
   write(ext,'(i4.4)') inodenr
   filename=basename(1:len_trim(basename)) // ext // '.vtu'

   ! write out vtu file for this node
   call sol2vtu01(ichoice,npoint,npelm,nelem,nphysi,ndim,idegfd,ks(isol)%sol,nusoli, & 
        & coor,kmeshc,kprobp_copy,filename,namedof,num_dofs_is_constant)

   ! inodenr goes from 1 to iacnodes
   if (inodenr == 1) then
      ! produce pvtu file 
      filename=basename(1:len_trim(basename)) // '.pvtu'
      inquire(unit=11,opened=lu_used)
      if (lu_used) then
         write(irefwr,*) 'PERROR(isol2vtu) : unit 11 is already used by another file'
         call instop 
      endif
      open(11,file=filename) 
      write(11,'(''<?xml version="1.0"?>'')') 
      write(11,'(''<VTKFile type="PUnstructuredGrid">'')') 
      write(11,'(''  <PUnstructuredGrid GhostLevel="0">'')')
      write(11,'(''    <PPointData>'')') 
      if (idegfd==0) then
         ! output all dofs
         if (ichoice==1) then
           do idof=1,nphys
              write(11,*) '      <PDataArray type="Float64" Name="',namedof(idof)(1:len_trim(namedof(idof))), &
                 & '" format="ascii"/>'
           enddo
         else if (ichoice==2) then
           ! special case for (u,v),p
           write(11,'(''        <PDataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii"/>'')')
           write(11,'(''        <PDataArray type="Float64" Name="Pressure" NumberOfComponents="1" format="ascii"/>'')')
         else
           write(irefwr,*) 'PERROR(sol2vtu01): ichoice is invalid: ',ichoice
           call instop
         endif
      else
         idof=idegfd
         write(11,*) '      <PDataArray type="Float64" Name="',namedof(1)(1:len_trim(namedof(1))), &
               & '" format="ascii"/>'
      endif
      write(11,'(''    </PPointData>'')') 
      write(11,'(''    <PPoints>'')') 
      write(11,'(''      <PDataArray type="Float32" NumberOfComponents="3" format="ascii"/>'')')
      write(11,'(''    </PPoints>'')') 
      do inode_here=1,iacnodes
         write(ext,'(i4.4)') inode_here
         filename=basename(1:len_trim(basename)) // ext // '.vtu'
         write(11,*) '    <Piece Source="',filename(1:len_trim(filename)),'"/>'
      enddo
      write(11,'(''  </PUnstructuredGrid>'')')
      write(11,'(''</VTKFile>'')')
      close(11)
   endif

else

   ! write(irefwr,*) 'output isol to file with basename ',basename(1:len_trim(basename))
   ! write(irefwr,*) 'in serial'
   if (len_trim(basename) + 4 > 120) then
      write(irefwr,*) 'PERROR(sol2vtu) :: length of basename is too much'
      write(irefwr,*) 'length should be less than 116 but is: ',len_trim(basename)
      call instop
   endif
   filename=basename(1:len_trim(basename)) // '.vtu'
   call sol2vtu01(ichoice,npoint,npelm,nelem,nphysi,ndim,idegfd,ks(isol)%sol,nusoli, & 
        & coor,kmeshc,kprobp_copy,filename,namedof,num_dofs_is_constant)
endif

if (allocated(kprobp_copy)) deallocate(kprobp_copy)
   
end subroutine sol2vtu

! Actual writing of the .vtu file 
! Write out coordinates coor(ndim,1:npoint)
! Write out nodal point numbers of the elements kmeshc(npelm,nelem)
! Write out offsets (npelm,2*npelm,etc.)
! Write out element type (VTK_TETRA=10; VTK_TRIANGE=5, VTK_QUADRATIC_TRIANGLE=22)
! Write out values of dofs per node usol(nphysi,npoint)
!     This latter part needs to be adjusted for variables dofs in the elements (quadratic TH element)
! PvK October 30 2017
! PvK October 31 2017: split out output of single scalar field (ichoice=1) and
! vector (u,v) + scalar (p) field
! PvK November 20 2017 Expansion to P2P1 element with interpolated pressure
subroutine sol2vtu01(ichoice,npoint,npelm,nelem,nphys,ndim,idegfd,usol,nusol,coor,kmeshc,kprobp,filename, & 
      & namedof,num_dofs_is_constant)
use sepmodulecomio
implicit none
integer, intent(in) :: ichoice,npoint,npelm,nelem,nphys,ndim,idegfd,kmeshc(:),nusol
integer, intent(inout) :: kprobp(:,:)
real(kind=8), intent(in) :: usol(:),coor(:,:)
character(len=*), intent(in) :: filename,namedof(:)
logical, intent(in) :: num_dofs_is_constant
integer :: i,j,ip,element_type,idof,ipp,ivelx,ively,ipres,ip1,ip2,ipp1,ipp2
character(len=120) :: namedofhere
real(kind=8), allocatable :: edge_node_pressure(:)
integer :: ip_node_pressure,allocate_status
real(kind=8) :: pressure_vertex1,pressure_vertex2

! Select VTK element type based on number of points in the element
element_type=0
if (npelm==3 .and. ndim==2) then
   ! linear triangle
   element_type=5
else if (npelm==4 .and. ndim==3) then
   ! linear tetrahedron
   element_type=10
else if (npelm==6 .and. ndim==2) then 
   ! quadratic triangle
   element_type=22
else 
   write(irefwr,*) 'PERROR(sol2vtu01): npelm / ndim combination invalid'
   write(irefwr,*) 'npelm = ',npelm, 'ndim = ',ndim
   call instop
endif

if (ichoice==2 .and. idegfd /= 0) then
   write(irefwr,*) 'PERROR(sol2vtu01): inconsistent input: ichoice=2 but idegfd=',idegfd 
   call instop
endif
!write(irefwr,*) 'num_dofs_is_constant: ',num_dofs_is_constant
if (.not.num_dofs_is_constant) then
   ! set up an array to store the pressure in the edges of the elements
   if (npoint*nphys-nusol <= 0) then
      write(irefwr,*) 'PERROR(sol2vtu01): in setting up kprobp_copy: '
      write(irefwr,*) 'PERROR(sol2vtu01): expecting npoint*nphys > nusol but: '
      write(irefwr,*) 'PERROR(sol2vtu01): npoint*nphys = ',npoint*nphys
      write(irefwr,*) 'PERROR(sol2vtu01): nusol        = ',nusol
      call instop
   endif
   !write(irefwr,*) 'set up edge_node_pressure array of length: ',npoint*nphys-nusol
   allocate(edge_node_pressure(npoint*nphys-nusol))
   ! kprobp is laid out as referring to u1,u2,u3,...,v1,v2,v3...,p1,p2,p3...
   ! where kprobp(:,3) is zero when a node is on the edge. We will compute the 
   ! pressure in those nodes by averaging the values in the surrounding vertices
   ! and store sequentially in edge_node_pressure; we'll store the negative of
   ! the index to edge_node_pressure  in the corresponding entries of kprobp(:,3).
   edge_node_pressure=1e7_8
   ! Loop over elements (redundant but algorithmically easy)
   ip=0
   ip_node_pressure=0
   do i=1,nelem
      !write(irefwr,'(''element : '',i5,'' has nodal point numbers: '',6i5)') i,(kmeshc((i-1)*npelm+j),j=1,npelm)
      !write(irefwr,'(''          '',5x,'' with pressure stored in  '',6i5)') (kprobp(kmeshc((i-1)*npelm+j),3),j=1,npelm)
      do j=2,npelm,2 ! loop over edge nodes
         ip = kmeshc((i-1)*npelm+j) ! global nodal point number of local node j in this element
!        write(irefwr,*) 'kprobp in point ',ip,kprobp(ip,1),kprobp(ip,2),kprobp(ip,3)
         ! sepran nodal point order in the element is vertex1, edge1, vertex2,
         ! edge2, vertex3, edge3
         if (kprobp(ip,3) == 0) then
            ! edge point has not yet been filled with pressure
            ip_node_pressure=ip_node_pressure+1
            if (ip_node_pressure > npoint*nphys - nusol) then
               write(irefwr,*) 'PERROR(sol2vtu01) ip_node_pressure is too large : ',ip_node_pressure
               write(irefwr,*) 'but should be <= ',npoint*nphys-nusol
               call instop
            endif
            ip1 = kmeshc((i-1)*npelm+mod(j-1,npelm))  ! nodal point number of vertex before edge 
            ip2 = kmeshc((i-1)*npelm+mod(j+1,npelm))  ! nodal point number of vertex after edge 
            !write(irefwr,'('' pressure solution in vertices is in usol at : '',2i5)') kprobp(ip1,3),kprobp(ip2,3)
            if (kprobp(ip1,3) == 0) then
               write(irefwr,*) 'PERROR: vertex pressure is not filled'
               write(irefwr,'(''ielem, j, ip, ipp1 = '',4i5)') i,j,ip,ipp1
               write(irefwr,'(''kprobp(ipp1) =       '',i5)') kprobp(ip1,3)
               call instop
            endif
            if (kprobp(ip2,3) == 0) then
               write(irefwr,*) 'PERROR: vertex pressure is not filled'
               write(irefwr,'(''ielem, j, ip, ipp2 = '',4i5)') i,j,ip,ipp2
               write(irefwr,'(''kprobp(ipp2) =       '',i5)') kprobp(ip2,3)
               !call instop
            endif
               
            pressure_vertex1=usol(kprobp(ip1,3))
            pressure_vertex2=usol(kprobp(ip2,3))
            edge_node_pressure(ip_node_pressure)=0.5*(pressure_vertex1+pressure_vertex2)
            ! write(irefwr,*) 'P: ',pressure_vertex1,pressure_vertex2,edge_node_pressure(ip_node_pressure)
            kprobp(ip,3) = - ip_node_pressure ! store negative index into kprobp
            !write(irefwr,'('' in elem '',i7,'' j= '',i7,'' ip/ipp/ip_node_pressure = '',9i7)') & 
            !   & i,j,ip,ipp,ip_node_pressure,nelem,npoint*nphys-nusol,nusol
         endif
      enddo
      ! write(irefwr,*) 
   enddo
   if (ip_node_pressure /= npoint*nphys - nusol) then
      write(irefwr,*) 'PERROR(sol2vtu01) after interpolation of pressure to edge nodes'
      write(irefwr,*) 'the total number of filled edges should be ',npoint*nphys-nusol
      write(irefwr,*) 'but is                                     ',ip_node_pressure
      call instop
   endif

!  do i=1,npoint
!     ivelx = kprobp(i,1)
!     ively = kprobp(i,2)
!     ipres = kprobp(i,3)
!     if (ipres<0) then
!        write(irefwr,'(i5,3i7,3f12.3)') i,ivelx,ively,ipres,usol(ivelx),usol(ively),edge_node_pressure(-ipres)
!     else
!        write(irefwr,'(i5,3i7,3f12.3)') i,ivelx,ively,ipres,usol(ivelx),usol(ively),usol(ipres)
!     endif
!  enddo
!  do i=1,npoint
!     ipres = kprobp(i,3)
!     if (ipres<0) then
!        write(irefwr,*) edge_node_pressure(-ipres)
!     else
!        write(irefwr,*) usol(ipres)
!     endif
!  enddo
endif


! XML header for vtu file
open(11,file=filename)
write(11,'(''<?xml version="1.0"?>'')')  ! note that using free format (11,*) adds a space before "<" which causes paraview to croak
write(11,'(''<VTKFile type="UnstructuredGrid" version="0.1">'')')
write(11,*) '    <UnstructuredGrid>'
write(11,'(''      <Piece NumberOfPoints="'',i10,''" NumberOfCells="'',i10,''">'')') npoint,nelem
! Write out coordinates coor(ndim,1:npoint)
write(11,*) '        <Points>'
write(11,*) '          <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
ip=0
if (ndim == 2) then
   do i=1,npoint
      write(11,'(3e15.7)') coor(1,i),coor(2,i),0e0
   enddo
else 
   do i=1,npoint
      write(11,'(3e15.7)') (coor(j,i),j=1,3)
   enddo
endif
write(11,*) '          </DataArray>'
write(11,*) '        </Points>'


! Now write out nodal point numbers of the elements kmeshc(npelm,nelem)
write(11,*) '        <Cells>'
write(11,*) '          <DataArray type="UInt32" Name="connectivity" format="ascii">' 
ip=0
if (npoint > int(1e7)) then
   write(irefwr,*) 'PERROR(sol2vtu01): npoint is too large for formatting in connectivity: ',npoint
   write(irefwr,*) 'should be less than 1e7'
   call instop
endif

if (element_type==22) then

  ! quadratic triangle with 6 nodal points.
  do i=1,nelem 
     ! sepran nodal numbers are: vertex 1, center 2, vertex 3, center 4, vertex 5, center 6 (counter clock wise)
     ! VTK numbers are: vertex 0, vertex 1, vertex 2, center 3, center 4, center 5 (counter clock wise)
     write(11,'(3i8,$)') (kmeshc(ip+j)-1,j=1,5,2)  ! subtract 1 to conform to C-style indexing
     write(11,'(3i8)') (kmeshc(ip+j)-1,j=2,6,2)  ! subtract 1 to conform to C-style indexing
     ip=ip+npelm
  enddo

else if (element_type==5) then

  ! linear triangle
  do i=1,nelem
     write(11,'(4i8:)') (kmeshc(ip+j)-1,j=1,npelm)  ! subtract 1 to conform to C-style indexing
     ip=ip+npelm
  enddo

else 
   
   write(irefwr,*) 'PERROR(sol2vtu01): not yet suited for element_type = ',element_type
   call instop

endif

write(11,*) '          </DataArray>'
! Write out offsets (npelm,2*npelm,etc.)
write(11,*) '          <DataArray type="UInt32" Name="offsets" format="ascii">'
write(11,'(100i10:)') (npelm*i,i=1,nelem) 
write(11,*) '          </DataArray>'
! Write out element type (VTK_TETRA=10; VTK_TRIANGE=5, VTK_QUADRATIC_TRIANGLE=22)
write(11,*) '          <DataArray type="UInt8" Name="types" format="ascii">'
write(11,'(100i3)') (element_type,i=1,nelem)
write(11,*) '          </DataArray>'
write(11,*) '        </Cells>'

! Write out actual solution; idegfd controls which dofs are output
write(11,*) '        <PointData>'

if (idegfd==0) then
   ! output all dofs
   if (ichoice==1) then

     ! output all dofs starting with 1 as scalar fields
     do idof=1,nphys
        namedofhere=namedof(idof)
        write(11,*) '        <DataArray type="Float64" Name="',namedofhere(1:len_trim(namedofhere)), &
                              & '" format="ascii">' 
        if (num_dofs_is_constant) then
           ! number of dofs per point in standard element is constant
           write(11,'(20e15.7,:)') (usol(idof+(i-1)*nphys),i=1,npoint)
        else 
           ! P2P1 values are in usol; interpolate pressure in edge nodes is in
           ! edge_node_pressure
           do i=1,npoint
              ip=kprobp(i,idof) 
              if (ip>0) then
                 write(11,'(e17.7,$)') usol(ip)
              else
                 write(11,'(e17.7,$)') edge_node_pressure(-ip)
              endif
           enddo
           write(11,*)
        endif
        write(11,*) '        </DataArray>'
     enddo

   else if (ichoice==2) then

     ! output dofs (1,2) as vector and dofs (3) as scalar
     write(11,'(''        <DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'')')
     if (num_dofs_is_constant) then
       ! number of dofs per point in standard element is constant
       do i=1,npoint
          write(11,'(3e15.7)') (usol(idof+(i-1)*nphys),idof=1,2),0.0_8
       enddo
     else 
       ! P2P1 values are in usol; interpolate pressure in edge nodes is in
       ! edge_node_pressure
       do i=1,npoint
          ip1=kprobp(i,1)
          ip2=kprobp(i,2)
          !write(irefwr,'(''about to write vector velocity for point '',i5,'' with solution in usol(ip1/ip2) where ip1/2='',2i5)') i,ip1,ip2
          if (ip1>0.and.ip2>0) then ! well, ip's shouldn't be <= 0....
             write(11,'(3e15.7,$)') usol(ip1),usol(ip2),0.0_8   ! (u,v,0)
          endif
       enddo
       write(11,*)
     endif
     write(11,*) '        </DataArray>'
     idof=3
     write(11,'(''        <DataArray type="Float64" Name="Pressure" NumberOfComponents="1" format="ascii">'')')
     if (num_dofs_is_constant) then
        write(11,'(20e15.7)') (usol(idof+(i-1)*nphys),i=1,npoint)
     else 
        ! P2P1 values are in usol; interpolate pressure in edge nodes is in
        ! edge_node_pressure
        do i=1,npoint
           ip=kprobp(i,idof)
           if (ip>0) then
              write(11,'(e15.7,$)') usol(ip)
              ! write(irefwr,*) 'pressure in point: ',i,ip,usol(ip)
           else
              write(11,'(e15.7,$)') edge_node_pressure(-ip)
              ! write(irefwr,*) 'pressure in point: ',i,ip,edge_node_pressure(-ip)
           endif
        enddo
        write(11,*)
     endif
     write(11,*) '        </DataArray>'

   else

     write(irefwr,*) 'PERROR(sol2vtu01): invalid choice for ichoice: ',ichoice
     call instop 

   endif

else 

   ! output only idegfd'th dof
   idof=idegfd
   namedofhere=namedof(1)
   write(11,*) '        <DataArray type="Float64" Name="',namedofhere(1:len_trim(namedofhere)), &
                            & '" format="ascii">' 
   if (num_dofs_is_constant) then
      write(11,'(20e15.7,:)') (usol(idof+(i-1)*nphys),i=1,npoint)
   else
      ! P2P1 values are in usol; interpolate pressure in edge nodes is in
      ! edge_node_pressure
      do i=1,npoint
         ip=kprobp(i,idof) 
         if (ip>0) then
            write(11,'(e15.7,$)') usol(ip)
         else
            write(11,'(e15.7,$)') edge_node_pressure(-ip)
         endif
      enddo
      write(11,*)
   endif
      
   write(11,*) '        </DataArray>'

endif
write(11,*) '        </PointData>'
write(11,'(''      </Piece>'')')
write(11,*) '    </UnstructuredGrid>'
write(11,'(''</VTKFile>'')')
close(11)

if (allocated(edge_node_pressure)) deallocate(edge_node_pressure)

end subroutine sol2vtu01
