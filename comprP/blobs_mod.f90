module blobs
implicit none
integer :: iblob,nblob
integer,parameter :: NBLOB_MAX=4,LU_BLOBS(1:NBLOB_MAX)=(/ (201+iblob,iblob=1,NBLOB_MAX) /)
real(kind=8),dimension(NBLOB_MAX) :: x_blob,y_blob,r_blob,delta_rhop
integer :: icenterblob(NBLOB_MAX)
real(kind=8) :: eta_blob
logical ::  CH94_blob
end module blobs
        
