!Contains info pertaining file I/O
module filedata
    implicit none

integer,parameter::BASISFILE=28,&!File containing information parameters required to run program BASISFILE = basis.in
    POTFILE=29,&!File containing potential information CONFILE = nadvibs.in
    OUTFILE=30,&!Majority of output information gets printed to OUTFILE = output.dat
    RESTARTFILE=32,&!File required to initiate a restart of lanczos algorithm
    ROOTINFO=33,&!Created only if idroots > 0
    ARCHIVE=34,&!DRA handle for archive file
    QRESTART=35,&!DRA handle for restart lanczos vectors
    SOINFO=36,&!File containing spinorbit parameters
    EIGVECS=37!File containing eigenvectors of H
character*100::outdir!Director to output text files to

contains
subroutine ShowTime()!Show date hour minute second
    integer,dimension(8)::time
    call date_and_time(values=time)
    write(unit=OUTFILE,fmt='(I4,A1,I2,A1,1x,I2,2x,I2,A1,I2,A1,I2)')time(1),'-',time(2),'-',time(3),'-',time(5),':',time(6),':',time(7)
end subroutine ShowTime

end module filedata
