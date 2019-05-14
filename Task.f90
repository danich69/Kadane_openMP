module Task
    use omp_lib
    implicit none
    contains
 
    subroutine GetMaxCoordinates(A, x1, y1, x2, y2)                     ! A - input matrix; x1, y1 - lower coordinates of submatrix; x2, y2 - upper coordinates of submatrix

        real(8), intent(in), dimension(:,:) :: A
        integer(4), intent(out) :: x1, y1, x2, y2
        
        real(8), dimension(size(A, dim = 2)) :: maximum_S                 ! array of maximum sums started at i-th row
        real(8), dimension(size(A, dim = 1)):: p                         ! array of column inner sums
        real(8) :: maximum, CurrentSum, t2, t1, T
        
        integer(4), dimension(2) :: coords
        integer(4), dimension(size(A, dim = 2)) :: maximum_L, maximum_R, maximum_B ! array of left, right, bottom coordinates of submatrixes with maximum sums started at i-th row
        integer(4) :: left, right                                       ! coordinates of maximum subarray of p 
        integer(4) :: i, j, m, n, k        

        if(maxval(A) < 0) then
            coords = maxloc(A)
            x2 = coords(1)
            x1 = coords(1)
            y2 = coords(2)
            y1 = coords(2)           
            return
        endif

        m = size(A, dim = 1)
        n = size(A, dim = 2)
        
        call omp_set_num_threads(1)
        
        !$omp parallel shared(a, m, n, maximum_S, maximum_L, maximum_R, maximum_B) private(i, j, left, right, CurrentSum, p)

        !$omp do schedule(dynamic)
        do i = 1, n                                                     ! main body of algorythm
            p = 0
            do j = i, n
                p = p + A(:, j)
                
                call Kande(p, left, right, CurrentSum)
                
                if (CurrentSum  >=  maximum_S(i) .or. i == j) then
                    maximum_S(i) = CurrentSum
                    maximum_L(i) = left
                    maximum_R(i) = right
                    maximum_B(i) = j
                endif
            
            enddo
        enddo
        !$omp end do

        !$omp end parallel
        
        y1 = maxloc(maximum_S, dim = 1)
        x1 = maximum_L(y1)
        y2 = maximum_B(y1)
        x2 = maximum_R(y1)

    end subroutine

    subroutine Kande(a, x1, x2, summary)                                ! 1D Kadane algorythm
    
        real(8), intent(in), dimension(:) :: a
        real(8), intent(out) :: summary
        real(8) :: Max_End, possible_1, possible_2
        
        integer(4) :: i, leftIndex, n
        integer(4), intent(out) :: x1, x2    

        n = size(a)

        summary = a(1)
        x1 = 1
        x2 = 1
        Max_End = a(1)
        leftIndex = 1
        
        do i = 2, n
            possible_1 = a(i)
            possible_2 = Max_End + a(i)

            if (possible_1 > possible_2) then
                Max_End = possible_1
                leftIndex = i
            else
                    Max_End = possible_2
            endif

            if (Max_End >= summary) then
                summary = Max_End
                x1 = leftIndex
                x2 = i
            endif
        enddo
        
    end subroutine

end module
