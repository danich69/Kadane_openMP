module omp_kadane_second_smirnov
    implicit none

    contains

        subroutine GetMaxCoordinates(A,x1,y1,x2,y2, thrdnmb)
            use omp_lib

            integer(4),intent(out) :: x1, y1, x2, y2
            integer(4) :: minor_matrix_size, maijor_matrix_size, i, j, k
            integer(4) :: left_edge, right_edge, top_edge, bottom_edge
            integer(4) :: buffer, minimal_posible, thread_num, thrdnmb
            real(8) :: current_sum, resulting_sum                               !maximum_sum

            real(8),dimension(:,:), intent(in) :: A
            real(8),allocatable :: current_column(:),B(:,:)
            real(8),allocatable :: result_sums(:)
            integer, allocatable :: coords(:,:)

            logical :: transp

            maijor_matrix_size = size( A, dim=1 )
            minor_matrix_size = size( A, dim=2 )
            transp = .FALSE.

            if ( maijor_matrix_size < minor_matrix_size ) then

                transp = .TRUE.
                allocate ( B ( size( A, dim = 2 ), size( A, dim = 1 ) ) )
                B = transpose( A )
                maijor_matrix_size = size( B, dim=1 )
                minor_matrix_size = size( B, dim=2 )

            else

                allocate ( B ( maijor_matrix_size, minor_matrix_size ) )
                B = A

            endif

            allocate( current_column ( maijor_matrix_size ) )

            !maximum_sum = B(1,1)
            x1 = 1
            y1 = 1
            x2 = 1
            y2 = 1

            call omp_set_num_threads(thrdnmb)

            !$omp parallel private ( current_column, right_edge, left_edge, top_edge, bottom_edge, resulting_sum, minimal_posible, current_sum ) default ( shared )

            !$omp single

                thread_num = omp_get_num_threads() - 1

                allocate( result_sums (0 : thread_num ) )
                allocate( coords ( 0 : thread_num, 4 ) )

                coords = 1
                result_sums = B(1,1)

            !$omp end single

            !$omp do schedule (dynamic)

            do left_edge = 1, minor_matrix_size

                current_column = B ( :, left_edge )

                do right_edge = left_edge, minor_matrix_size

                    if (right_edge > left_edge) then

                        current_column = current_column + B ( :, right_edge )

                    endif

                    current_sum = current_column( 1 )
                    top_edge = 1
                    bottom_edge = 1
                    resulting_sum = 0
                    minimal_posible = 0;

                    do i = 1, size ( current_column )

                        resulting_sum = resulting_sum+current_column( i )

                        if ( resulting_sum > current_sum ) then

                            current_sum = resulting_sum
                            top_edge = minimal_posible+1
                            bottom_edge = i

                        endif

                        if ( resulting_sum < 0 ) then

                            resulting_sum = 0
                            minimal_posible = i

                        endif

                    enddo

                    k = omp_get_thread_num()

                    if ( current_sum > result_sums( k ) ) then

                        result_sums( k ) = current_sum
                        coords( k, 1 ) = top_edge
                        coords( k, 2 ) = bottom_edge
                        coords( k, 3 ) = left_edge
                        coords( k, 4 ) = right_edge

                    end if


                enddo

            enddo

            !$omp end do

            !$omp end parallel

            k = maxloc ( result_sums(:), dim = 1 )

            x1 = coords( k, 1 )
            x2 = coords( k, 2 )
            y1 = coords( k, 3 )
            y2 = coords( k, 4 )

            deallocate( result_sums, coords )
            deallocate( current_column )

            if (transp) then

                buffer = x1
                x1 = y1
                y1 = buffer

                buffer = y2
                y2 = x2
                x2 = buffer

            endif

        end subroutine GetMaxCoordinates
end module
