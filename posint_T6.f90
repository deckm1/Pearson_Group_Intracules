!---------------------------------------------------------------------------------------------
!	Calculate a Position Intracule (using data from a .wfn file)  
!
!		Written by Dalton Mackenzie (2013) 
!		with much thanks to Adam, Brendan, Zosia and Dr. Pearson
!
!		Much of the wavefunction reading code comes from Zosia, a previous member of the group.
!
!      This is the version of the Code with resubmit/block capability.
!      Tolerance has been included, and the threshold value is 1.0E-6.
!       
!       Calling method:
!       ./outfile wfnFile.wfn intFile.int COEFS.coef Molecular_Orbital_# chlineStart chlineEnd
!
!---------------------------------------------------------------------------------------------
MODULE SHELL_PAIR_DETERMINATION
	integer sSetCount, pSetCount, dSetCount
	integer shellPairListSize, shellPairCallListSize
	integer curShellPairIndex, curShellPairCallIndex
	integer shellPairSSoffset, shellPairPSoffset, shellPairPPoffset
	integer shellPairDSoffset, shellPairDPoffset, shellPairDDoffset
	integer SSshellPairCount, PSshellPairCount, PPshellPairCount
	integer DSshellPairCount, DPshellPairCount, DDshellPairCount
	integer size_ssss, size_psss, size_ppss, size_psps, size_ppps
	integer size_pppp, size_dsss, size_dpss, size_dsps, size_ddss
	integer size_dsds, size_dpps, size_dspp, size_ddps, size_dpds
	integer size_ddds, size_dppp, size_ddpp, size_dpdp, size_dddp
	integer size_dddd
	
!	Arrays for building our shell pair list, and consequently call lists
	integer, dimension(:), allocatable :: sOptions
	integer, dimension(:), allocatable :: pOptions
	integer, dimension(:), allocatable :: dOptions
	
!	Arrays for Shell Pair lists and to call Josh's positionRR.f90
	integer*8, dimension(:), allocatable :: shellPairList
	integer*8, dimension(:), allocatable :: shellPairCallList	

END MODULE

MODULE HF2PDM
!Here we will store the 2d array used for the 2 particle density matrix.
	double precision, dimension(:,:), allocatable :: pMatrixFromCoefs						!2d array used for 2 particle density matrix
END MODULE

	program position_intracule

	USE SHELL_PAIR_DETERMINATION
	USE HF2PDM
!	USE DIMENSIONS

	implicit none
	integer i, j,ia,ib,ic,id, a, b, c, d, FileStat, MO, lineStart, lineEnd, curLine			
	character g*8, gaussian*8, input*50, output*50, localized*50, orbital*50				! These characters are for error messages, and input/output files.
	character chlineStart*50, chlineEnd*50													! lineStart and lineEnd as read from command line	
	integer NMO, oNPRIM, NNUC																! Number of molecular orbitals, primitives, and nuclei in the system
  integer breakFromCallLoopCount
	double precision, dimension(:,:), allocatable :: COOR									! Coordinate array for each nucleus of size(NNUC x 3) (for x,y,z)
  integer, dimension(:), allocatable :: indexedZeroes               ! This array will have length oNprim and start with a 1 at every index.
                                                                    ! As zeroes are determined, a zero will replace every index belonging
                                                                    ! to a set where every coef falls below TOLERANCE.

	integer, dimension(:), allocatable :: oCENT												! Old Centres (including 0 COEFS)
	integer, dimension(:), allocatable :: oTYPE												! Old Types (including 0 COEFS)
	double precision, dimension(:), allocatable :: oEXPO									! Old Exponents (including 0 COEFS)
	double precision, dimension(:,:), allocatable :: oCOEF									! Old COEFS (including 0s)

	integer, dimension(:), allocatable :: classSizeList										! array to track class sizes (for more accurate addition) 
	integer classSizeCounter																! a counter for the array above

	double precision cor_a(3), cor_b(3), cor_c(3), cor_d(3)									! Gives the coordinates of a,b,c,d
	integer getAngMom, getGaussSum1to, getPrefactorCode, ang_a, ang_b, ang_c, ang_d			! Gives the angular momentum of a,b,c,d
  integer curAngMom, bfsetIrrelevant, curBF, otherCurBF, otherCurAngMom, lcode               !Tracks the angmom of a bf in question
	double precision PHFtotal, PHFab, PHFcd, PHFad, PHFbc, HFT, P, curClassSum				! 2 particle density matrix elements
	double precision expo_a, expo_b, expo_c, expo_d											! Exponents of centre's a,b,c,d for the current call to positionRR.f90
	integer psInCall, dsInCall, posint_size, posint_index									! Prefactor values
	integer aLoopIndex, bLoopIndex, cLoopIndex, dLoopIndex									! Prefactor values
	integer tempa, tempb, tempc, tempd														! Prefactor values
	integer curPrefactor, verbose 
  integer*8 bfsequencesAccountedFor,bfsequencesExpected 														! Prefactor values
	double precision, dimension(1296) :: Pos_int 											! Position integral array returned from positionRR.f90
																							! it has size 1296 because dddd would need 6^4 = 1296 indices
!-----------------------------------------------------------------------
!	MEGA ARRAY OFFSETS
!-----------------------------------------------------------------------

	integer sCount,pCount,dCount															! Count the number of each type of basis function s,p,d
	double precision, dimension(:), allocatable :: WFN_MegaArray							! WFN_MegaArray - a megaarray for the overall wavefunction
	double precision, dimension(:), allocatable :: sArray									! sArray - an array for all of the s type atomic orbitals
	double precision, dimension(:), allocatable :: pArray									! pArray - an array for all of the p type atomic orbitals
	double precision, dimension(:), allocatable :: dArray									! dArray - an array for all of the d type atomic orbitals
	integer megaArraySize, curInt
	integer typeOffset, exponentOffset, centreOffset, coefOffset, angMomOffset				! Note that only the coefs for that one orbital will be used.
	integer sTypeOffset, sExponentOffset, sCentreOffset, sCoefOffset, sindexOffset			! Note that only the coefs for this particular MO will be used.
	integer pTypeOffset, pExponentOffset, pCentreOffset, pCoefOffset, pindexOffset			! Note that only the coefs for this particular MO will be used.
	integer dTypeOffset, dExponentOffset, dCentreOffset, dCoefOffset, dindexOffset			! Note that only the coefs for this particular MO will be used.
	integer NPRIM, curNewRef, sIndex, pIndex, dIndex,zeroCount							! Values to track the removal of zeroes, and splitting of mega array.

!-----------------------------------------------------------------------
! ADDITIONAL PARAMETERS
!-----------------------------------------------------------------------

	double precision pi, XU, U, MAXLINES, R, TOLERANCE										! XU, MAXLINES, R scale the output plot (U vs P(u))
	integer*8 MAXNPRIM
	parameter (pi = 3.14159265359d0)
	parameter (MAXLINES = 250.0d0)															! Parameters for Mura-Knowles quadrature
	parameter (R = 3.00d0)																	! Can set up such that IMAX and R are entered by User.
	parameter (TOLERANCE = 1.0D-6) 															! SET COEF TOLERANCE TO 1 x 10 ^ -6
	parameter (MAXNPRIM = 10000)															! MAXNPRIM is the maximum number of primitives, and it MUST be a power of 10.
																							! it is used to separate basis function numbers in the shellPairCallList.
!-----------------------------------------------------------------------
! Getting the arguments from the command line.
!-----------------------------------------------------------------------

	call GETARG(1, input)																	! First argument from the command line is the input .wfn file
	call GETARG(2, output)																	! Second argument is the output .dat file
	call GETARG(3, localized)																! Third argument is the localized coefficient file (USE 'NA' FOR NON LOCALIZED)
	call GETARG(4, orbital)																	! Fourth argument is the orbital being calculated (0 is total intracule)
	call GETARG(5, chlineStart)																! Fifth argument is the line to begin calculating the intracule
	call GETARG(6, chlineEnd)																! Sixth argument is the line to end calculating the intracule
	read (orbital,*) MO																		! converts the string "orbital" into the integer, "MO"
	read (chlineStart,*) lineStart															! converts the string "chlineStart" into the integer "lineStart"
	read (chlineEnd,*) lineEnd																! converts the string "chlineEnd" into the integer "lineEnd"

!-----------------------------------------------------------------------
! Attempting to open the .wfn file for reading
!-----------------------------------------------------------------------

	open(unit = 7,file = input,IOSTAT=FileStat,status='old')								! unit 7 is the .wfn file
	  if (FileStat > 0) then																
	    print *, 'Error opening input file:'												
	    print *, ' -> Please ensure the .wfn file is correct, ' 
	    print *, '    and the filename does not contain illegal characters.'
	    stop
	  endif

	open(unit = 8,file = output)															! unit 8 is the output file that u and P(u) are written to

	open(unit = 10,file = localized,IOSTAT=FileStat,status='old')							! unit 10 is the localized .coef file
	  if (FileStat > 0) then
	    print *, 'Error opening localized orbital coefficients file:'
	    print *, ' -> Please ensure the file is correct, ' 
	    print *, '    and the filename does not contain illegal characters.'
	    stop
	  endif

!----------------------------------------------------------------------
!	Read .wfn file
!----------------------------------------------------------------------

	read (unit = 7, FMT = 100) G,NMO,oNPRIM,NNUC											! reads NMO, NPRIM, NNUC from input
	gaussian = 'GAUSSIAN'
	if (g.ne.gaussian) then																	!  Possible Error: if 'GAUSSIAN' is not present on the second line of the wfn file
	  print *, 'ERROR: Invalid input file type(no GAUSSIAN statement)'
	  print *, '	    -> Ensure the input is a .wfn file.'
	  stop
	endif

	ALLOCATE(COOR(1:nnuc,1:3))
	j = 3
	do i = 1, nnuc
	read (7,"(25x,f11.8,1x,f11.8,1x,f11.8)") coor(i,1:j)									! reads the coordinates of each nucleus from .wfn
	enddo						

	ALLOCATE(oCENT(1:onprim))
	read (7, 101) (ocent(i),i=1,onprim)														! reads the centre assignments

	ALLOCATE(oTYPE(1:onprim))
	read (7, 101) (otype(i),i=1,onprim)														! reads the type assignments

	ALLOCATE(oEXPO(1:onprim))
	read (7, 102) (oexpo(i),i=1,onprim)														! reads the exponents

!	Input formats:
100 format (/a8,11x,I4,15x,I5,16x,I4) 														! 100 is for G,NMO, oNPRIM, NNUC
																							! Read 8 characters ('GAUSSIAN'), skip 11 characters
																							! read 4 digits, skip 15, read 5, skip 16, read 4
101	format (20x,20i3)																		! 101 is for centre and type assignments
																							! skip 20 characters [20x], then read 20 3 character long integers
102	format (11x,5e14.8)																		! 102 is the specific format used to read the exponents
																							! skip 11 spaces and then read in 5 doubles of the form e14.8
103	format (5(e15.9))																		! 103 is the specific format used to read in coefficients
																							! Read in 5 doubles of the format e15.9
!----------------------------------------------------------------------
!	DECREASE ARRAY SIZES AND NPRIM VALUES.
!----------------------------------------------------------------------
!	Build a coefficient array of size nmo by number of primitives
	ALLOCATE(oCOEF(1:nmo,1:onprim))

!  Read localized coefficients, if necessary:
	if(localized/='NA')then
	do i = 1, NMO
	    read (10, 103) (ocoef(i,j),j=1,onprim)
	enddo
	endif

!  Normalize localized coefficients, if necessary:
	if(localized/='NA')then
		do i = 1, onprim
			if(otype(i)==1) then
	   			do j = 1, nmo
					ocoef(j,i) = ocoef(j,i)*(((2.0d0/pi)**(3.0d0/4.0d0))*((oexpo(i))**(3.0d0/4.0d0)))
	    		enddo
	  		else if((otype(i)==2).or.(otype(i)==3).or.(otype(i)==4))then
	    		do j = 1, nmo
					ocoef(j,i) = ocoef(j,i)*2*(((2.0d0/pi)**(3.0d0/4.0d0))*((oexpo(i))**(5.0d0/4.0d0)))
	    		enddo
	  		else if((otype(i)==5).or.(otype(i)==6).or.(otype(i)==7))then
	    		do j = 1, nmo
					ocoef(j,i) = ocoef(j,i)*4*((((2.0d0/pi)**(3.0d0/4.0d0))*((oexpo(i))**(7.0d0/4.0d0)))/(dsqrt(3.0d0)))		
	    		enddo
	  		else if((otype(i)==8).or.(otype(i)==9).or.(otype(i)==10))then
	    		do j = 1, nmo
					ocoef(j,i) = ocoef(j,i)*4*(((2.0d0/pi)**(3.0d0/4.0d0))*((oexpo(i))**(7.0d0/4.0d0)))
	    		enddo
	  		else
	    		print *, 'ERROR: Are you trying to do f-orbitals? ... You crazy fool!'
	    		stop
	  		endif
		enddo
	endif

!-------------------------------------------------------------------------------
!	This is the section where tolerance has been implemented.
!
!
!	It never hurts to be ambitious. Unless you're McBeth
!
!	*Note, for tolerance to work here, we need to work in sets of orbitals.
!	An S orbital can be removed if its coefficient is below the given
!	tolerance value, but a p-orbital set can only be removed if all 3
!	that is, x,y and z p-orbitals have coefficients below the tolerance
!	value. Similarly with a d-orbital set, the xx,yy,zz,xy,xz and yz orbitals
!	must all have coefficient values that are effectively zero.
!
!-------------------------------------------------------------------------------
!
! Algorithm for using tolerance with basis function sets.
! The print lines that have been commented out can be used to determine which basis functions
! have been cut by tolerance.

  ALLOCATE(indexedZeroes(1:onprim))

  do i = 1,onprim
    indexedZeroes(i) = 1
  enddo
  if((MO.ge.1).and.(MO.le.NMO)) then  !MO is a valid MO and tolerance can be used.
    zeroCount = 0
    i = 1
    do while(i.le.onprim)
      curAngMom = getAngMom(oType(i))
      if(curAngMom.eq.0) then
        !S Case
        if(abs(ocoef(MO,i)).le.TOLERANCE) then
!             print *, "There was an s basis function to be cut at bf: ", i, "with coef ", ocoef(MO,i)
          indexedZeroes(i) = 0
          zeroCount = zeroCount + 1
        endif
        i = i + 1    
      elseif(curAngMom.eq.1) then
        !P Case
          bfsetIrrelevant = 1 !Assume the set of p orbitals is useless.
          j = 0
          do while((j.lt.3).and.(bfsetIrrelevant.eq.1))
            curBF = i+j
            if(abs(ocoef(MO,curBF)).gt.TOLERANCE) then
              bfsetIrrelevant = 0
            endif
            j = j+1
          enddo
          if(bfsetIrrelevant.eq.1) then
!             print *, "There were p basis functions to be cut when at bf: ", i
            indexedZeroes(i) = 0
            indexedZeroes(i+1) = 0
            indexedZeroes(i+2) = 0
            zeroCount = zeroCount + 3
          endif
        i = i + 3
      elseif(curAngMom.eq.2) then
        !D Case
          bfsetIrrelevant = 1 !Assume the set of d orbitals is useless.
          j = 0
          do while((j.lt.6).and.(bfsetIrrelevant.eq.1))
            curBF = i+j
            if(abs(ocoef(MO,curBF)).gt.TOLERANCE) then
              bfsetIrrelevant = 0
            endif
            j = j+1
          enddo
          if(bfsetIrrelevant.eq.1) then
 !            print *, "There were d basis functions to be cut when at bf: ", i
            indexedZeroes(i) = 0
            indexedZeroes(i+1) = 0
            indexedZeroes(i+2) = 0
            indexedZeroes(i+3) = 0
            indexedZeroes(i+4) = 0
            indexedZeroes(i+5) = 0
            zeroCount = zeroCount + 6
          endif
        i = i + 6
      endif  !end ang mom if-else-if    
  
    enddo ! enddo while(i.le.onprim)

  endif !end is MO between 1 and NMO if statement

NPRIM = onprim - zeroCount

!    print *,""
!    print *,"Old NPRIM = ",onprim
!    print *,"zeroCount = ",zeroCount
!    print *,"New NPRIM = ",nprim
!    print *,""

!``````````````````````````````````````````````````````````````````````````````````````````````````````
!	ALLOCATE THE MEGA ARRAY NOW
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!	The array will have 3*NNUC indices for x,y,z coordinates, plus 5*NPRIM indices for the centre, type, exponent, coef and ang mom of each basis function
	megaArraySize = 3*NNUC + 5*NPRIM !+ 1
	ALLOCATE(WFN_MEGAARRAY(1:megaArraySize))

! Offsets to reference back to the data later
	centreOffset = 3*NNUC 																	! The centres are located immediately after the coordinates (3*NNUC)
	typeOffset = 3*NNUC + NPRIM 															! The types are located NPRIM spaces after the centres begin
	exponentOffset = 3*NNUC + 2*NPRIM 														! The exponents are located NPRIM spaces after the types begin
	coefOffset = 3*NNUC + 3*NPRIM 															! The coefs are located NPRIM spaces after the exponents begin
	angMomOffset = 3*NNUC + 4*NPRIM															! The ang moms are located NPRIM spaces after the coefs for that MO begin.

	scount = 0																				!A count for the basis functions of angular momentum 0
	pcount = 0																				!A count for the basis functions of angular momentum 1
	dcount = 0																				!A count for the basis functions of angular momentum 2

!``````````````````````````````````````````````````````````````````````````````````````````````````````
!	Let's begin by copying over the coordinates to the mega array
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

	do j = 1, NNUC		!For every nucleus:
		wfn_megaarray(j*3 - 2) = coor(j,1)													! Add the x coordinate for the given nucleus to the mega array
    
		wfn_megaarray(j*3 - 1) = coor(j,2)													! Add the y coordinate for the given nucleus to the mega array

		wfn_megaarray(j*3) = coor(j,3)														! Add the z coordinate for the given nucleus to the mega array

	enddo
	
!``````````````````````````````````````````````````````````````````````````````````````````````````````
!	Now add the other assignments to the megaarray
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!	*Please note, if you want to save yourself from some reading, the next little section is only
!	to help people understand why curNewRef is needed in the loop to cut values which are smaller
!	than the tolerance value.

!	Begin with a curNewRef value of 1, this will be our tracker for our current position in the new
!	array. We will need it if some of the values are not to be copied over from the old values to
!	the new tolerance clipped values.

!	A simple example of this would be if our coefficient array (originally) was:
!	{ 0.5 , 0.6 , 0 , 1.2, 0 }
!	We would start with a counter, say j, and scan from 1 to 5 through the array.
!	Our new output array starts like this: {}
!	with a curNewRef value of 1.
!	We look at 0.5, and see that it has an important value, we add it to the new array and increment j, and now we have:
!	{ 0.5 } with a curNewRef of 2 and a j value of 2
!	We look at the next value oldCoefs[j] and see that it's a 0.6. Now we add that to the new array, increment j, and now we have:
!	{ 0.5 , 0.6 } with a curNewRef of 3 and a j value of 3
!	Next we look at oldCoefs[j] and see that it's a zero. Since we don't want zeroes, we DO NOT add this to the new array, we increment j,
!	but curNewRef remains at 3.
!	Looking at oldCoefs[4] now, we see 1.2, and add it to the new array. We increment j and curNewRef, and we have:
!	{ 0.5 , 0.6 , 1.2 } , a curNewRef of 4 and a j value of 5.
!	Finally, we will end up skipping the last zero, finishing with a curNewRef of 4, and a j value of 6
!	(both of these counters point to position SIZE + 1, as they are to reference the NEXT element in each of their respective arrays.

! For each basis function, add the centre, type, exponent, coefficient and angular momentum to the mega array


	curNewRef = 1
	do j = 1, oNPRIM																		!Scan through all of the old basis functions
		if (indexedZeroes(j).eq.1) then	!THIS LINE USED TO BE: .gt.TOLERANCE, but for now we want ALL of the basis functions
			wfn_megaarray(centreOffset + curNewRef) = oCENT(j)
			wfn_megaarray(typeOffset + curNewRef) = oTYPE(j)
			wfn_megaarray(exponentOffset + curNewRef) = oEXPO(j)
			wfn_megaarray(coefOffset + curNewRef) = oCOEF(MO,j)
			wfn_megaarray(angMomOffset + curNewRef) = getAngMom(oType(j))	
		curNewRef = curNewRef + 1														!Increment as long as we have actually added an entry to the new array.
		endif		!End the if checking whether or not the coefficient has passed the tolerance test. (No tolerance in initial version of this code)
	enddo

!``````````````````````````````````````````````````````````````````````````````````````````````````````
! So now that the angular momentums are in the megaarray, we should scan through and count them.
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

	do i = 1, NPRIM																			! For all of the basis functions
		curInt = nint(wfn_megaarray(angMomOffset + i))		! Turn ang mom from double to integer		
		
    select case (curInt)
		case (0)																			! If ang mom (L) = 0, then s orbital
			scount = scount + 1
		case (1)																			! If ang mom (L) = 1, then p orbital
			pcount = pcount + 1
		case (2)																			! If ang mom (L) = 2, then d orbital
			dcount = dcount + 1
		end select	
	enddo

!``````````````````````````````````````````````````````````````````````````````````````````
! 	We will also build 3 subarrays for s p and d angular momentums.
!	Remember, these subarrays need to remember their previous index as well.
!	But it's kind of silly to keep track of their angular momentums anymore,
!	as they are inherent to the array we choose to store them in. So for each
!	one, we will no longer store the angular momentum, but now we will store
!	the old index of the basis function within the larger wavefunction mega array.
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

	if(scount.gt.0) then
		ALLOCATE(sarray(1:(3*NNUC + 5*scount))) !Only allocate the array if there is information for it.
		do j = 1, NNUC
!	Put the coordinates in the s array
			sarray(j*3 - 2) = coor(j,1)
			sarray(j*3 - 1) = coor(j,2)
			sarray(j*3) = coor(j,3)
		enddo	
	endif

	if(pcount.gt.0) then
		ALLOCATE(parray(1:(3*NNUC + 5*pcount))) !Only allocate the array if there is information for it.
		do j = 1, NNUC
!	Put the coordinates in the p array
			parray(j*3 - 2) = coor(j,1)
			parray(j*3 - 1) = coor(j,2)
			parray(j*3) = coor(j,3)
		enddo
	endif

	if(dcount.gt.0) then
		ALLOCATE(darray(1:(3*NNUC + 5*dcount))) !Only allocate the array if there is information for it.
		do j = 1, NNUC
!	Put the coordinates in the d array
			darray(j*3 - 2) = coor(j,1)
			darray(j*3 - 1) = coor(j,2)
			darray(j*3) = coor(j,3)
		enddo
	endif

	!Now let's fill the arrays with the appropriate information from the mega array.	
	!This uses the same idea with curNewRef underneath the tolerance section from earlier.	
	sindex = 1
	pindex = 1
	dindex = 1

    if(scount.gt.0) then																	! sArray offsets
        sCentreOffset = 3*NNUC
        sTypeOffset = sCentreOffset + scount
        sExponentOffset = sTypeOffset + scount
        sCoefOffset = sExponentOffset + scount
        sIndexOffset = sCoefOffset + scount
    endif

    if(pcount.gt.0) then																	! pArray offsets
       pCentreOffset = 3*NNUC
       pTypeOffset = pCentreOffset + pcount
       pExponentOffset = pTypeOffset + pcount
       pCoefOffset = pExponentOffset + pcount
       pIndexOffset = pCoefOffset + pcount
    endif

    if(dcount.gt.0) then																	! dArray offsets
        dCentreoffset = 3*NNUC
        dTypeOffset = dCentreOffset + dcount
        dExponentOffset = dTypeOffset + dcount
        dCoefOffset = dExponentOffset + dcount
        dIndexOffset = dCoefOffset + dcount
   endif

	do i = 1, NPRIM
		curInt = nint(wfn_megaarray(angMomOffset + i))
		select case (curInt)
		case (0) 																			! Add to s array
			sarray(scentreoffset + sindex) = wfn_megaarray(centreOffset + i)
			sarray(stypeOffset + sindex) = wfn_megaarray(typeOffset + i)
			sarray(sexponentOffset + sindex) = wfn_megaarray(exponentOffset + i)
			sarray(scoefOffset + sindex) = wfn_megaarray(coefOffset + i)
			sarray(sindex) = i       										! this means basis function i is an s orbital
			sindex = sindex + 1

		case (1)																			! Add to p array
			parray(pcentreoffset + pindex) = wfn_megaarray(centreOffset + i)
			parray(ptypeOffset + pindex) = wfn_megaarray(typeOffset + i)
			parray(pexponentOffset + pindex) = wfn_megaarray(exponentOffset + i)
			parray(pcoefOffset + pindex) = wfn_megaarray(coefOffset + i)
			parray(pindex) = i       										! this means basis function i is a p orbital
			pindex = pindex + 1

		case (2)																			! Add to d array
			darray(dcentreoffset + dindex) = wfn_megaarray(centreOffset + i)
			darray(dtypeOffset + dindex) = wfn_megaarray(typeOffset + i)
			darray(dexponentOffset + dindex) = wfn_megaarray(exponentOffset + i)
			darray(dcoefOffset + dindex) = wfn_megaarray(coefOffset + i)
			darray(dindex) = i       										! this means basis function i is a d orbital
			dindex = dindex + 1
		end select	
	enddo

	if(MOD(pCount , 3).ne.0) then
		print *, 'The wfn file is flawed, there is a set of p orbitals with less than 3 p bfs.'
		stop
	endif

	if(MOD(dCount, 6).ne.0) then
		print *, 'The wfn file is flawed, there is a set of d orbitals with less than 6 d bfs.'
		stop
	endif

	sSetCount = scount
	pSetCount = pcount / 3
	dSetCount = dcount / 6	

!----------------------------------------------------------------------
!	WE WILL NOW COMPUTE THE 2D ARRAY USED FOR THE HF2PDM	
!----------------------------------------------------------------------
!First we shall allocate the array.
	if(MO.ne.0)then
		ALLOCATE(pMatrixFromCoefs(1:nprim,1:nprim))

!Now we shall use a loop to fill it
		do i = 1,nprim
			do j = 1, nprim
				pMatrixFromCoefs(i,j) = (wfn_megaarray(coefOffset + i)) * (wfn_megaarray(coefOffset + j)) * 2
			enddo
		enddo
	endif


!----------------------------------------------------------------------
!	MAKE THE SHELL PAIR ARRAY FOR THE SYSTEM
!----------------------------------------------------------------------
!	To do this we will first need to determine the size of this
!	array. Then we can set up a few offsets to work with it.
!	From there we can work through 6 double loops to fill it.

!	Start out by getting the sizes of the different sections of the shell pair array,
!	and building the offsets.

	SSshellPairCount = getGaussSum1to(sSetCount)
	PSshellPairCount = pSetCount * sSetCount
	PPshellPairCount = getGaussSum1to(pSetCount)
	DSshellPairCount = dSetCount * sSetCount
	DPshellPairCount = dSetCount * pSetCount
	DDshellPairCount = getGaussSum1to(dSetCount)

	shellPairSSoffset = 0
	shellPairPSoffset = shellPairSSoffset + SSshellPairCount
	shellPairPPoffset = shellPairPSoffset + PSshellPairCount
	shellPairDSoffset = shellPairPPoffset + PPshellPairCount
	shellPairDPoffset = shellPairDSoffset + DSshellPaircount
	shellPairDDoffset = shellPairDPoffset + DPshellPaircount

	!Allocate shell pair list
	shellPairListSize = SSshellPairCount + PSshellPairCount + PPshellPairCount + DSshellPairCount + DPshellPairCount + DDshellPairCount
	ALLOCATE(shellPairList(1:shellPairListSize))
	curShellPairIndex = 1

!````````````````````````````````````````````````````````````````````````````````````````
!	NOW LETS FILL THE SHELLPAIR ARRAY, USING 6 LOOPS (OPTIONAL BASED ON IF)
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
! [ss]
	if(SSshellPairCount.gt.0)then
		do i = 1, sSetCount
			do j = 1, i
				
				shellPairList(curShellPairIndex) = MAXNPRIM * sarray(i) + sarray(j)

				curShellPairIndex = curShellPairIndex + 1
			enddo			
		enddo
	endif

! [ps]
	if(PSshellPairCount.gt.0)then
		do i = 1, pSetCount
			do j = 1, sSetCount
				
				shellPairList(curShellPairIndex) = MAXNPRIM * parray((3*i) - 2) + sarray(j)
	
			curShellPairIndex = curShellPairIndex + 1
			enddo
		enddo
	endif

! [pp]
	if(PPshellPairCount.gt.0)then

		do i = 1, pSetCount
			do j = 1, i
				
				shellPairList(curShellPairIndex) = MAXNPRIM * parray((3*i) - 2) + parray((3*j) - 2)

				curShellPairIndex = curShellPairIndex + 1
			enddo			
		enddo
	endif

! [ds]
	if(DSshellPairCount.gt.0)then
		do i = 1, dSetCount
			do j = 1, sSetCount
				
				shellPairList(curShellPairIndex) = MAXNPRIM * darray((6*i) - 5) + sarray(j)

				curShellPairIndex = curShellPairIndex + 1
			enddo			
		enddo
	endif

! [dp]
	if(DPshellPairCount.gt.0)then
		do i = 1, dSetCount
			do j = 1, pSetCount

				shellPairList(curShellPairIndex) = MAXNPRIM * darray((6*i) - 5) + parray((3*j) - 2)

				curShellPairIndex = curShellPairIndex + 1
			enddo			
		enddo
	endif

! [dd]
	if(DDshellPairCount.gt.0)then
		do i = 1, dSetCount
			do j = 1, i

				shellPairList(curShellPairIndex) = MAXNPRIM * darray((6*i) - 5) + darray((6*j) - 5)

				curShellPairIndex = curShellPairIndex + 1
			enddo			
		enddo
	endif

! We need to build the proper lists of shell pairs to work with.
! Since [ssss] is symmetrical between the two atoms, we include the ab <= cd rule for [abcd] shell pair generation.
! So we need to build a list of shell pairs using units from 1 to scount (i = i <= scount) 
!	The shell pair list generated here should be : {11,12,13,22,23,33}, which in turn should be used to generate:
!	{1111,1112,1113,1122,1123,1133,1212,1213,1222,1223,1233,1313,1322,1323,1333,2222,2223,2233,2323,2333,3333}
!	Two rules are applied here, a <= b and c <= d, and as mentioned previously, because the atoms are symmetrical,
!	We add the rule whereby ab (represented as a*10 + b) <= cd (represented as c*10 + d) so ab must be <= cd to avoid doubles.
!	We will now generate these lists, so that we have the proper indices to use when calling Josh's code.
!	We only have to work with S basis functions here, so first let's figure out the shell
!	It makes sense to build the total shellPairArray here and now.
!	This array will have size: 
!	gaussSum1to(sSetCount) + gaussSum1to(pSetCount) + gaussSum1to(dSetCount) 
!	+ (sSetCount * pSetCount) + (sSetCount * dSetCount) + (pSetCount * dSetCount)
!	We need to perform the following work:
!	ss, ps, pp, ds, dp, dd
!	So to build the shell call list, we will add class groups one at a time.

!----------------------------------------------------------------------
!	Calculate Position Intracule (output P(u) and u)
!----------------------------------------------------------------------
!	First build the shellPairCallList, and then work with it at the end.
!	We must allocate shellPairCallList but we need the size of shellPairCallList first:
  size_ssss = getGaussSum1to(ssshellpaircount)											  !size_ssss
  size_psss = psshellpaircount * ssshellpaircount											!size_psss
  size_ppss = ppshellpaircount * ssshellpaircount											!size_ppss
  size_psps = getGaussSum1to(psshellpaircount)											  !size_psps
  size_ppps = ppshellpaircount * psshellpaircount											!size_ppps
  size_pppp = getGaussSum1to(ppshellpaircount)											  !size_pppp
  size_dsss = dsshellpaircount * ssshellpaircount											!size_dsss
  size_dpss = dpshellpaircount * ssshellpaircount											!size_dpss
  size_dsps = dsshellpaircount * psshellpaircount											!size_dsps
  size_ddss = ddshellpaircount * ssshellpaircount											!size_ddss
  size_dsds = getGaussSum1to(dsshellpaircount)											  !size_dsds
  size_dpps = dpshellpaircount * psshellpaircount											!size_dpps
  size_dspp = dsshellpaircount * ppshellpaircount											!size_dspp
  size_ddps = ddshellpaircount * psshellpaircount											!size_ddps
  size_dpds = dpshellpaircount * dsshellpaircount											!size_dpds
  size_ddds = ddshellpaircount * dsshellpaircount											!size_ddds
  size_dppp = dpshellpaircount * ppshellpaircount											!size_dppp
  size_ddpp = ddshellpaircount * ppshellpaircount											!size_ddpp
  size_dpdp = getGaussSum1to(dpshellpaircount)											  !size_dpdp
  size_dddp = ddshellpaircount * dpshellpaircount											!size_dddp
  size_dddd = getGaussSum1to(ddshellpaircount)											  !size_dddd

	shellPairCallListSize = size_ssss + size_psss + size_ppss + size_psps + size_ppps + size_pppp
	shellPairCallListSize = shellPairCallListSize + size_dsss + size_dpss + size_dsps + size_ddss
	shellPairCallListSize = shellPairCallListSize + size_dsds + size_dpps + size_dspp + size_ddps
	shellPairCallListSize = shellPairCallListSize + size_dpds + size_ddds + size_dppp + size_ddpp
	shellPairCallListSize = shellPairCallListSize + size_dpdp + size_dddp + size_dddd

	allocate(shellPairCallList(1:shellPairCallListSize))
		
	allocate(classSizeList(1:21))															! An array to add a full class to P(u) at a time to improve accuracy! 
	
	classSizeList(1) = size_ssss
	classSizeList(2) = size_psss + classSizeList(1)
	classSizeList(3) = size_ppss + classSizeList(2)
	classSizeList(4) = size_psps + classSizeList(3)
	classSizeList(5) = size_ppps + classSizeList(4)
	classSizeList(6) = size_pppp + classSizeList(5)
	classSizeList(7) = size_dsss + classSizeList(6)
	classSizeList(8) = size_dpss + classSizeList(7)
	classSizeList(9) = size_dsps + classSizeList(8)
	classSizeList(10) = size_ddss + classSizeList(9)
	classSizeList(11) = size_dsds + classSizeList(10)
	classSizeList(12) = size_dpps + classSizeList(11)
	classSizeList(13) = size_dspp + classSizeList(12)
	classSizeList(14) = size_ddps + classSizeList(13)
	classSizeList(15) = size_dpds + classSizeList(14)
	classSizeList(16) = size_ddds + classSizeList(15)
	classSizeList(17) = size_dppp + classSizeList(16)
	classSizeList(18) = size_ddpp + classSizeList(17)
	classSizeList(19) = size_dpdp + classSizeList(18)
	classSizeList(20) = size_dddp + classSizeList(19)
	classSizeList(21) = size_dddd + classSizeList(20)

	curShellPairCallIndex = 1
![ss ss]	
	do i = 1, SSshellPairCount
		do j = 1, i
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(i) + shellPairList(j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo

![ps ss]
	do i = 1, PSshellPairCount
		do j = 1, SSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairPSoffset+i) + shellPairList(shellPairSSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![pp ss]
	do i = 1, PPshellPairCount
		do j = 1, SSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairPPoffset+i) + shellPairList(shellPairSSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![ps ps]
	do i = 1, PSshellPairCount
		do j = 1, i	
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairPSoffset+i) + shellPairList(shellPairPSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![pp ps]
	do i = 1, PPshellPairCount

		do j = 1, PSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairPPoffset+i) + shellPairList(shellPairPSoffset+j)

			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![pp pp]
	do i = 1, PPshellPairCount
		do j = 1, i
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairPPoffset+i) + shellPairList(shellPairPPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo

![ds ss]
	do i = 1, DSshellPairCount
		do j = 1, SSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDSoffset+i) + shellPairList(shellPairSSoffset+j)			
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dp ss]
	do i = 1, DPshellPairCount
		do j = 1, SSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDPoffset+i) + shellPairList(shellPairSSoffset+j)	
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![ds ps]
	do i = 1, DSshellPairCount
		do j = 1, PSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDSoffset+i) + shellPairList(shellPairPSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd ss]
	do i = 1, DDshellPairCount
		do j = 1, SSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairSSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![ds ds]
	do i = 1, DSshellPairCount
		do j = 1, i
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDSoffset+i) + shellPairList(shellPairDSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dp ps]
	do i = 1, DPshellPairCount
		do j = 1, PSshellPairCount			
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDPoffset+i) + shellPairList(shellPairPSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![ds pp]
	do i = 1, DSshellPairCount
		do j = 1, PPshellPairCount			
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDSoffset+i) + shellPairList(shellPairPPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd ps]
	do i = 1, DDshellPairCount
		do j = 1, PSshellPairCount			
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairPSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dp ds]
	do i = 1, DPshellPairCount
		do j = 1, DSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDPoffset+i) + shellPairList(shellPairDSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd ds]
	do i = 1, DDshellPairCount
		do j = 1, DSshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairDSoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dp pp]
	do i = 1, DPshellPairCount
		do j = 1, PPshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDPoffset+i) + shellPairList(shellPairPPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd pp]
	do i = 1, DDshellPairCount
		do j = 1, PPshellPairCount
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairPPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dp dp]
	do i = 1, DPshellPairCount
		do j = 1, i
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDPoffset+i) + shellPairList(shellPairDPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd dp]
	do i = 1, DDshellPairCount
		do j = 1, DPshellPairCount	
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairDPoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo
![dd dd]
	do i = 1, DDshellPairCount
		do j = 1, i
			shellPairCallList(curShellPairCallIndex) = MAXNPRIM ** 2 * shellPairList(shellPairDDoffset+i) + shellPairList(shellPairDDoffset+j)
			curShellPairCallIndex = curShellPairCallIndex + 1
		enddo
	enddo

! Start using the shellPairCallList to call positionRR.f90
	do curLine = lineStart, lineEnd
	  classSizeCounter = 1
	  P = 0.0D0
	  curClassSum = 0.0D0
	  xu = curLine/(MAXLINES + 1)
	  u = -R * (log(1 - (xu**3)))

    bfsequencesExpected = nprim ** 4
    bfsequencesAccountedFor = 0

	  do i = 1, shellPairCallListSize		

! First get a,b,c,d the basis function numbers
	    a = shellPairCallList(i) / MAXNPRIM ** 3
	    b = MOD((shellPairCallList(i) / MAXNPRIM ** 2), MAXNPRIM)
	    c = MOD((shellPairCallList(i) / MAXNPRIM), MAXNPRIM)
	    d = MOD((shellPairCallList(i)), MAXNPRIM)

!       print *,a,b,c,d

! Now get everything else.
! To figure out the number of values to pull from Josh's array, we will need to 
! know the number of pSets and dSets within the call. To determine that,
! reference the angmom values of these a,b,c and d basis functions within
! the megaarray. From there, you should determine a curPSetCount and a curDSetCount.
! The number of values to pull from Josh's array = 3 ** pSetCount * 6 ** dSetCount ( * 1 ** sSetCount = 1)
! Now we should fill the values needed to call Josh's code
		
!*****************************************************************************************************************
	    cor_a(1) = wfn_megaarray((wfn_megaarray(centreOffset + a) - 1)*3 + 1)
      cor_a(2) = wfn_megaarray((wfn_megaarray(centreOffset + a) - 1)*3 + 2)
      cor_a(3) = wfn_megaarray((wfn_megaarray(centreOffset + a) - 1)*3 + 3)

	    cor_b(1) = wfn_megaarray((wfn_megaarray(centreOffset + b) - 1)*3 + 1)
      cor_b(2) = wfn_megaarray((wfn_megaarray(centreOffset + b) - 1)*3 + 2)
      cor_b(3) = wfn_megaarray((wfn_megaarray(centreOffset + b) - 1)*3 + 3)

	    cor_c(1) = wfn_megaarray((wfn_megaarray(centreOffset + c) - 1)*3 + 1)
      cor_c(2) = wfn_megaarray((wfn_megaarray(centreOffset + c) - 1)*3 + 2)
      cor_c(3) = wfn_megaarray((wfn_megaarray(centreOffset + c) - 1)*3 + 3)

	    cor_d(1) = wfn_megaarray((wfn_megaarray(centreOffset + d) - 1)*3 + 1)
      cor_d(2) = wfn_megaarray((wfn_megaarray(centreOffset + d) - 1)*3 + 2)
      cor_d(3) = wfn_megaarray((wfn_megaarray(centreOffset + d) - 1)*3 + 3)

	    ang_a = wfn_megaarray(angMomOffset + a)
	    ang_b = wfn_megaarray(angMomOffset + b)
	    ang_c = wfn_megaarray(angMomOffset + c)
	    ang_d = wfn_megaarray(angMomOffset + d)

!*****************************************************************************************************************

	    !With these ang values, we can determine the number of ps and ds in this current call.
	    !With that knowledge, we can determine the size of the array.

	    psInCall = 0
	    dsInCall = 0

	    posint_size = 0

	    aLoopIndex = 0
	    bLoopIndex = 0
	    cLoopIndex = 0
	    dLoopIndex = 0

	    call determinePsAndDsInCall(ang_a,ang_b,ang_c,ang_d,aLoopIndex,bLoopIndex,cLoopIndex,dLoopIndex,psInCall,dsInCall,posint_size)

	    !So now that we have the 4 loop indices, we are about to call Josh's code.

	    expo_a = wfn_megaarray(exponentOffset + a)
	    expo_b = wfn_megaarray(exponentOffset + b)
	    expo_c = wfn_megaarray(exponentOffset + c)
	    expo_d = wfn_megaarray(exponentOffset + d)

	    !We have gotten the exponent values, and we are now calling Josh's code.
		
    call PositionIntegral(cor_a,cor_b,cor_c,cor_d,expo_a,expo_b,expo_c,expo_d,ang_a,ang_b,ang_c,ang_d,u,verbose,Pos_int)

	    posint_index = 1

		!We should now do a loop through the 4 loop indices to account for all of the posint
		!values returned by Josh's posint array.

		curPrefactor = getPrefactorCode(a,b,c,d,MAXNPRIM)	!We have 5 possible values here, 1 through 5. 

		!They mean the following:

		!1: group 1
		!2: group 1
		!3: group 1

		!No flip required, multiply by 2**(getPrefactorCode-1) and proceed.

		!4: group 2
		!5: group 2

		!A flip is required. For 4, multiply by 2, flip and run again, also multiplying by 2
		!For 5, multiply by 4, flip and run again, also multiplying by 4.

		do ia = 0,(aLoopIndex - 1)
			do ib = 0,(bLoopIndex - 1)
				do ic = 0,(cLoopIndex - 1)
					do id = 0,(dLoopIndex - 1)

						if(curPrefactor.ge.1.and.curPrefactor.le.3) then

							!Within here, we will iterate through the posint values in Josh's array,
							!all the while determining the new coefficient values to multiply by,
							!and multiplying by the prefactor for the given sequence call to 
							!Josh's code, contracting to add the integrals to our overall value.

							tempa = a + ia	
							tempb = b + ib
							tempc = c + ic
							tempd = d + id

							if (MO == 0) then				! Calculate ALL molecular orbitals
								PHFab = PHFtotal(tempa,tempb,nmo,nprim,ocoef)	
		  						PHFcd = PHFtotal(tempc,tempd,nmo,nprim,ocoef)
		  						PHFad = PHFtotal(tempa,tempd,nmo,nprim,ocoef)
		  						PHFbc = PHFtotal(tempb,tempc,nmo,nprim,ocoef)
		  						HFT = (((2*PHFab*PHFcd)-(PHFad*PHFbc))/4.0d0)

							else if ((MO > 0).and.(MO <= NMO)) then		! Calculate ONE molecular orbital
								PHFab = pMatrixFromCoefs(tempa,tempb)
			  					PHFcd = pMatrixFromCoefs(tempc,tempd)
		  						PHFad = pMatrixFromCoefs(tempa,tempd)
		  						PHFbc = pMatrixFromCoefs(tempb,tempc)
		  						HFT = ((0.5d0*PHFab*PHFcd)-(0.250d0*PHFad*PHFbc))
							else
								print*, 'ERROR: You have entered an invalid value for molecular orbital.'
								stop
							endif

							curClassSum = curClassSum + HFT * pos_int(posint_index) * (2**(curPrefactor-1))	

              bfsequencesAccountedFor = bfsequencesAccountedFor + (2**(curPrefactor-1))

						elseif(curPrefactor.eq.4.or.curPrefactor.eq.5) then

							!Within here, we will iterate through the posint values in Josh's array,
							!all the while determining the new coefficient values to multiply by,
							!and multiplying by the prefactor for the given sequence call to 
							!Josh's code, contracting to add the integrals to our overall value.

							tempa = a + ia	
							tempb = b + ib
							tempc = c + ic
							tempd = d + id

							if (MO == 0) then				! Calculate ALL molecular orbitals
								PHFab = PHFtotal(tempa,tempb,nmo,nprim,ocoef)	
		  						PHFcd = PHFtotal(tempc,tempd,nmo,nprim,ocoef)
		  						PHFad = PHFtotal(tempa,tempd,nmo,nprim,ocoef)
		  						PHFbc = PHFtotal(tempb,tempc,nmo,nprim,ocoef)		  					
		  						HFT = (((2*PHFab*PHFcd)-(PHFad*PHFbc))/4.0d0)
							else if ((MO > 0).and.(MO <= NMO)) then		! Calculate ONE molecular orbital
								PHFab = pMatrixFromCoefs(tempa,tempb)
			  					PHFcd = pMatrixFromCoefs(tempc,tempd)
		  						PHFad = pMatrixFromCoefs(tempa,tempd)
		  						PHFbc = pMatrixFromCoefs(tempb,tempc)
		  						HFT = ((0.5d0*PHFab*PHFcd)-(0.250d0*PHFad*PHFbc))
							else
								print*, 'ERROR: You have entered an invalid value for molecular orbital.'
								stop
							endif

							curClassSum = curClassSum + HFT * pos_int(posint_index) * (2**(curPrefactor-3))							
              bfsequencesAccountedFor = bfsequencesAccountedFor + (2**(curPrefactor-3))
!--------------------------------------------------------------
!	THE FOLLOWING IS NOT A MISTAKE!
!	WE NEED TO FLIP ONE OF THE SETS, AND GET THE NEW HFT VALUE,
! SO WE'LL FLIP A AND B.
!--------------------------------------------------------------						
							tempa = b + ib
							tempb = a + ia
							tempc = c + ic
							tempd = d + id

							if (MO == 0) then				! Calculate ALL molecular orbitals
								PHFab = PHFtotal(tempa,tempb,nmo,nprim,ocoef)	
		  						PHFcd = PHFtotal(tempc,tempd,nmo,nprim,ocoef)
		  						PHFad = PHFtotal(tempa,tempd,nmo,nprim,ocoef)
		  						PHFbc = PHFtotal(tempb,tempc,nmo,nprim,ocoef)			
		  						HFT = (((2*PHFab*PHFcd)-(PHFad*PHFbc))/4.0d0)
							else if ((MO > 0).and.(MO <= NMO)) then		! Calculate ONE molecular orbital
								PHFab = pMatrixFromCoefs(tempa,tempb)
			  					PHFcd = pMatrixFromCoefs(tempc,tempd)
		  						PHFad = pMatrixFromCoefs(tempa,tempd)
		  						PHFbc = pMatrixFromCoefs(tempb,tempc)
		  						HFT = ((0.5d0*PHFab*PHFcd)-(0.250d0*PHFad*PHFbc))
							else
								print*, 'ERROR: You have entered an invalid value for molecular orbital.'
								stop
							endif
							curClassSum = curClassSum + HFT * pos_int(posint_index) * (2**(curPrefactor-3))
              bfsequencesAccountedFor = bfsequencesAccountedFor + (2**(curPrefactor-3))			
						endif	!curPrefactor loop
						posint_index = posint_index + 1
					enddo ! end id loop
				enddo ! end ic loop
			enddo ! end ib loop
		enddo ! end ia loop
!If we have finished all of the sequence calls for this class, then we can work on the next class

			if(i == classSizeList(classSizeCounter))then			
			P = P + curClassSum
			curClassSum = 0.0D0
			classSizeCounter = classSizeCounter + 1
			endif
		enddo !end do i = 1 -> shellPairCallListSize

		write (8,110) u,P
	enddo	!end curLine loop
	 110 format	(2(e15.8,8x))	! To write u and P(u) to the output .dat file
  end		!end program

!---------------------------------------------------------------------------
!	FUNCTIONS
!---------------------------------------------------------------------------
!```````````````````````````````````````````````````````````````````````````
!	getGaussSum1to(n) gives the sum from 1 to n going by increments of 1
!	using the gaussian formula. 
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	integer function getGaussSum1to(n)
	implicit none
	integer n, output
	output = 0
	output = (n * (n + 1)) / 2
	getGaussSum1to = output
	return  
	end
!```````````````````````````````````````````````````````````````````````````
!	getAngMom(x) takes a type assignment, and converts it to an angular momentum
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	integer function getAngMom(x)
        implicit none
        integer x
        if (x == 1) then 
        	getAngMom = 0 
        else if ((x == 2).OR.(x == 3).OR.(x == 4)) then
        	getAngMom = 1 
        else if (x >= 5) then 
        	getAngMom = 2 
        endif   
        return  
	end
!````````````````````````````````````````````````````````````````````````````
!	determinePrefactor(x) takes a 4 digit input (representing an integral,
!	and determines the symmetry prefactor for it.
!	x is of the form [abcd] where abcd = 1000a + 100b + 10c + d
!	and a,b,c,d represent specific basis functions.
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	integer function getPrefactorCode(valA,valB,valC,valD,MAXNPRIM)
        implicit none
        integer valA, valB, valC, valD
        integer degreeOfSymmetry1, degreeOfSymmetry2, degreeOfSymmetry3
        integer posintPrefactor, oneAtomConsistent
		integer*8 MAXNPRIM
		degreeOfSymmetry1 = 0  !atm 1
		degreeOfSymmetry2 = 0  !atm 2
		degreeOfSymmetry3 = 0  !both atms
!	is a == b?
        if (vala == valb) then 
          degreeOfSymmetry1 = 1
        else
        	degreeOfSymmetry1 = 2
        endif
!	is c == d         
         if( valc == vald ) then
         	degreeOfSymmetry2 = 1
         else
         	degreeOfSymmetry2 = 2
         endif
!	is ab == cd
	if((vala*MAXNPRIM+valb) == (valc*MAXNPRIM+vald))then
		degreeOfSymmetry3 = 1
	else
		degreeOfSymmetry3 = 2
	endif
	
	posintPrefactor = degreeOfSymmetry1 * degreeOfSymmetry2 * degreeOfSymmetry3

	if(degreeOfSymmetry1 == 1 .or. degreeOfSymmetry2 == 1) then
		oneAtomConsistent = 1
	else
		oneAtomConsistent = 0
	endif
!If oneAtomConsistent = 1, then the situation is in group 1, and falls under one of the following:
!aaaa
!bbaa
!baaa
!bcaa	
	if(oneAtomConsistent == 1)then
		if(posintPrefactor == 1) then
			getPrefactorCode = 1
		elseif(posintPrefactor == 2) then
			getPrefactorCode = 2
		elseif(posintPrefactor == 4) then
			getPrefactorCode = 3
		endif	!End posintPrefactor if statement
	endif	!End if oneAtomConsistent == 1
!Otherwise, the situation is in group 2, and falls under one of the following:
!baba
!baca
!abcd
	if(oneAtomConsistent == 0) then
		if(posintPrefactor == 4) then
			getPrefactorCode = 4
		elseif(posintPrefactor == 8) then
			getPrefactorCode = 5
		endif !End posintPrefactor if statement
	endif !End if oneAtomConsistent == 0
	return
	end
!``````````````````````````````````````````````````````````````````````````````````````````````
!	PHFTotal is a function that determines the proper scaling coefficient to 
!	multiply each integral by for a total molecular intracule.
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	double precision function PHFtotal(x,y,nmo,nprim,coef)
	implicit none
	integer x,y,i,nmo,nprim
	double precision coef(nmo,nprim)
	PHFtotal = 0
	do i = 1, nmo
	  PHFtotal = PHFtotal + (coef(i,x)*coef(i,y))
	enddo
	PHFtotal = PHFtotal*2
	return
	end
!``````````````````````````````````````````````````````````````````````````````````````````````
	subroutine determinePsAndDsInCall(ang_a,ang_b,ang_c,ang_d,aLI,bLI,cLI,dLI,psInCall,dsInCall,posint_size)
!``````````````````````````````````````````````````````````````````````````````````````````````
! LI here stands for Loop Index
! aLoopIndex represents the number of dimensions to consider for the ath basis function
! The same is true for bLI,cLI, and dLI for bth,cth, and dth basis function respectively
	implicit none	
	integer ang_a,ang_b,ang_c,ang_d
	integer psInCall, dsInCall
	integer aLI,bLI,cLI,dLI
	integer posint_size
! a Loop Index
	select case (ang_a)
	case(0)
		aLI = 1	
	case(1)
		aLI = 3
		psInCall = psInCall + 1
	case(2)
		aLI = 6
		dsInCall = dsInCall + 1
	end select
! b Loop Index
	select case (ang_b)
	case(0)
		bLI = 1	
	case(1)
		bLI = 3
		psInCall = psInCall + 1
	case(2)
		bLI = 6
		dsInCall = dsInCall + 1
	end select
! c Loop Index
	select case (ang_c)
	case(0)
		cLI = 1	
	case(1)
		cLI = 3
		psInCall = psInCall + 1
	case(2)
		cLI = 6
		dsInCall = dsInCall + 1
	end select
! d Loop Index
	select case (ang_d)
	case(0)
		dLI = 1	
	case(1)
		dLI = 3
		psInCall = psInCall + 1
	case(2)
		dLI = 6
		dsInCall = dsInCall + 1
	end select	
	posint_size = 3 ** psInCall * 6 ** dsInCall
	end
