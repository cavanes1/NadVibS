module progdata
    implicit none

integer,parameter::zero=0,one=1,two=2,three=3
real*8,parameter::zerodp=0d0,ELECMASSU=1822.888484929d0,&!Mass of a proton in atomic units
    MWAU2WAVE=5140.49227189d0,&!Conversion factor for mass weighted a.u. to cm-1, specifically, -> 219474.63*sqrt(1/ELECMASSU)
    AU2WAVE=219474.6313705d0,&!Conversion factor for a.u. to cm-1
    AU2EV=27.211384523d0!Conversion factor for a.u. to eV
logical::neworigin,&!Shift reference geometry away from initial state
    orthog,&!TRUE if re-orthogonalization is requested
    orthogexact,&!TRUE if exact dot products are to be used to monitor lanczos vector orthogonality
    firststep,&!Global orthogonalization, when performed, must be done on subsequent steps
    restartrun,&!TRUE if current run is restarted from previous results
    savevecs!TRUE if lanczos vectors are to be saved
integer::dimen,&!The product of the number of basis functions for each mode
    ordr,&!Order of the expansion
    nstblks,&!The number of unique state blocks
    nproc,&!Number of processors employed
    myid,&!Unique name of the current process
    nstates,&!Number of electronic states
    nmodes,&!Number of vibrational modes
    nmax,&!Maximum number of basis function in any given mode
    natoms,&!Number of atoms in molecule
    niter,&!Number of lanczos iterations to be performed
    iiter,&!Initial iteration, = 1 if new job, = iteration to restart from otherwise
    totstart,totend,iterstart,iterend,&!Timing variable
    maxstor,&!The maximum number of lanczos vectors that can be stored on disk
    northog,&!The number of orthogonalization batches to perform
    chkorthog,&!If using recurrence, compute exact dot products every chkorthog iterations
    idroots,&!How many of the lead contributors to a vibronic level to identify
    soroots,&!How many roots to compute spin orbit parameters for
    nseg,&!Number of segements to divide lanczos vector into
    nirreps,&!Number of irreps in point group symmetry
    q1,q2,q3,&!Global array handles
    ndiag,&!Number of diagonal terms
    nmatel!Length of hamiltonian mat. el. array
integer,allocatable,dimension(:)::nfunc,&!Array holding the number of harmonic oscillator functions for each node
    noterms,&!Number of unique coupling terms per order
    shft,&!Offset array for indexing basis sets
    loindex,&!Index to start orthogonalizing when doing partial orthogonalization
    roindex,&!Index to stop orthogonalizing lanczos vector when doing P.O.
    istate,&!Initial state of the reference (i.e. anion ground state)
    nalltrms,&!Total number of terms in a processor batch
    nztrms!Number of non-zero unique terms per order (unique means this term consists of all different dimensions)
integer,allocatable,dimension(:,:)::nzindx,&!Index array of non-zero potential terms
    nzblks!Block locations of non-zero potential terms
integer::ndblks!Number of blocks that have diagonal contributions, excluding (i,i), which do by definition
integer,allocatable,dimension(:)::dblks!Identity of those blocks
integer,dimension(:,:,:),allocatable::basind,bstart,rcshft!hterms: data for doing matrix vector product
integer,dimension(:,:),allocatable::nelems!Number of non-zero elements in a particular section of Hamiltonian
integer,dimension(:),allocatable::nsteps!Number of different matrix element types for a potential term
integer,dimension(:,:),allocatable::basdif,stype,&!Delta values between harmonic oscillator matrix elements
    bmax,&!b value at which counter resets. Does that help?  
    iordr!Array of indicies for location of h.o. matrix elements in homatel
integer,dimension(:),allocatable::dordr,&!Order of diagonal term i
    dmap,&!Location of term information for diagonal term i
    npirrep!Number of coordinates per irrep
integer,dimension(:,:),allocatable::vbounds,&!Batch sizes
    cpindx!Coupling block indices
integer,dimension(2)::lo,hi!Indices of the Lanczos vector a particular process is responsible for
real*8,dimension(:),allocatable::concoef!dcoef: constant terms for each state block
real*8,dimension(:,:),allocatable::nzcoef!Array of non-zero potential coefficients
real*8::bjiconv,&!The degree of convergence before eigenvalues are printed -> 1e-(bjiconv)
    epsilon,&!Machine precision of the architecture
    eta,&!Epsilon^(3/4)
    maxdisk,&!Maximum amount of disk available for use (in MB)
    ztoler!Any coefficient less than 1e-(ztoler) is set to 0
real*8,dimension(:),allocatable::homatel,&!Harmonic oscillator matrix elements
    alpha,beta,&!Hold the alpha and beta constants used to construct tridiagonal matrix
    aomega,&!Value of omega for harmonic oscillator functions
    bomega!Alpha coefficients for determining overlap with reference state
real*8,dimension(:,:),allocatable::omega!Monitors lanczos vector orthogonality, used only if re-orthog is requested
real*8,dimension(:),allocatable::statew,&!Initial weights of the vibronic states
    cocoef,&!Constant coefficients
    dvec,&!d vector for determining overlap with reference state
    tmat!T matrix for determining overlap with reference state
real*8,dimension(:,:),allocatable::scr1,scr2!Scratch buffers

end module progdata