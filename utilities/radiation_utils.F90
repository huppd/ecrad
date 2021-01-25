module radiation_utils
    contains

    pure subroutine reorder_dims_1_3(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output

      output = reshape(input, shape(output), order=(/3, 2, 1/))
    end subroutine reorder_dims_1_3

    pure subroutine reorder_dims_1_2(input, output)
      use parkind1, only : jprb
      real(jprb), dimension(:,:,:), intent(in) :: input
      real(jprb), dimension(:,:,:), intent(out) :: output

      output = reshape(input, shape(output), order=(/2, 1, 3/))
    end subroutine reorder_dims_1_2

end module radiation_utils
