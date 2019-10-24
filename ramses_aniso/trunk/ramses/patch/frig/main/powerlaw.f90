subroutine sample_powerlaw(x, a, b, alpha, n)
    ! Sample from a power-law between a and b, with an index of alpha (for the PDF)
    use amr_commons
    use pm_commons
    use random


    implicit none


    real(8), dimension(1:n), intent(out):: x
    real(8), intent(in):: a, b, alpha
    integer, intent(in):: n

    integer ,dimension(1:ncpu,1:IRandNumSize)::allseed

    real(8):: u, p, q
    integer:: i

    p = alpha + 1.0_8
    q = 1.0_8 / p


    ! If necessary, initialize random number generator
    if(localseed(1)==-1)then
      call rans(ncpu,iseed,allseed)
      localseed=allseed(myid,1:IRandNumSize)
    end if


    do i = 1, n
        call Ranf( localseed, u )

        write(*,*) 'random number generated ', u
        
!        call random_number(u)
        ! u follows an uniform law between 0 and 1
        ! Scale it to b^p..a^p
        u = b**p + (a**p - b**p) * u

        ! Calculate x(i)
        x(i) = u**q
    end do
end subroutine sample_powerlaw

