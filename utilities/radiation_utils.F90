module radiation_utils
    implicit none

    private :: transpose_real, transpose_block_real
    integer, parameter :: block_size = 64

    contains

    pure subroutine reorder_dims_1_3(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output
      integer :: n(1:3), i

      n = shape(input)
      do i = 1, n(2)
        call transpose_real(input(:,i,:), output(:,i,:))
      end do
      !output = reshape(input, shape(output), order=(/3, 2, 1/))
    end subroutine reorder_dims_1_3

    pure subroutine reorder_dims_1_2(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output
      integer :: n(1:3), i

      n = shape(input)
      do i = 1, n(3)
        call transpose_real(input(:,:,i), output(:,:,i))
      end do
      !output = reshape(input, shape(output), order=(/2, 1, 3/))
    end subroutine reorder_dims_1_2

    !This subroutine performs cache-friendly transposition by processing the matrix in "small" blocks of 64x64
    pure subroutine transpose_real(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:), intent(in) :: input
      real(jprb), dimension(:,:), intent(out) :: output
      real(jprb), dimension(block_size,block_size) :: from, to
      integer :: ri ,ci, rbi, cbi, n_r_blocks, n_c_blocks, n_rows, n_cols, i, rem_rows, rem_cols

      n_rows = size(input, 2)
      n_cols = size(input, 1)
      n_r_blocks = n_rows / block_size
      rem_rows = mod(n_rows, block_size)
      n_c_blocks = n_cols / block_size
      rem_cols = mod(n_cols, block_size)
      ri = 1
      do rbi = 1, n_r_blocks
        ci = 1
        do cbi = 1, n_c_blocks
          from = input(ci:ci+block_size-1, ri:ri+block_size-1)
          to = transpose(from)
          output(ri:ri+block_size-1, ci:ci+block_size-1) = to
          ci = ci + block_size
        end do
        if (rem_cols > 0) then
            from(1:rem_cols, :) = input(ci:n_cols, ri:ri+block_size-1)
            to(:, 1:rem_cols) = transpose(from(1:rem_cols, :))
            output(ri:ri+block_size-1, ci:n_cols) = to(:, 1:rem_cols)
        end if
        ri = ri + block_size
      end do
      if (rem_rows > 0) then
        ci = 1
        do cbi = 1, n_c_blocks
          from(:, 1:rem_rows) = input(ci:ci+block_size-1, ri:n_rows)
          to(1:rem_rows, :) = transpose(from(:, 1:rem_rows))
          output(ri:n_rows, ci:ci+block_size-1) = to(1:rem_rows, :)
          ci = ci + block_size
        end do
        if (rem_cols > 0) then
            from(1:rem_cols, 1:rem_rows) = input(ci:n_cols, ri:n_rows)
            to(1:rem_rows, 1:rem_cols) = transpose(from(1:rem_cols, :))
            output(ri:n_rows, ci:n_cols) = to(1:rem_rows, 1:rem_cols)
        end if
      end if

!      do ri = 1, n_rows
!        do ci = 1, n_cols
!          output(ri, ci) = input(ci, ri)
!        end do
!      end do
      !output = transpose(input)
    end subroutine transpose_real

end module radiation_utils
