! Created on 02.02.21.
! Author: Mikhail Zhigun

module m_memory_pool
   use, intrinsic :: iso_c_binding

   implicit none

    private

    type, public  :: stack_memory_pool
        integer(kind=c_int8_t), dimension(:), pointer, private :: mem => null()
        integer(kind=c_intptr_t), private        :: start_ptr = 0
        integer(kind=c_intptr_t), private        :: end_ptr = 0
    contains
        procedure :: create => stack_memory_pool_create
        !procedure :: destroy => stack_memory_pool_destroy
        procedure :: allocate => stack_memory_pool_allocate
        procedure :: deallocate => stack_memory_pool_deallocate
        final :: stack_memory_pool_destroy
    end type stack_memory_pool

    interface sallocate_array
        module procedure sallocate_int_1d
    end interface sallocate_array

    interface sdeallocate_array
        module procedure sdeallocate_int_1d
    end interface sdeallocate_array

    public :: sallocate_array, sdeallocate_array
  contains
    subroutine stack_memory_pool_create(this, byte_size)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_int32_t), intent(in) :: byte_size

        if (this%start_ptr == 0) then
            allocate(this%mem(byte_size))
            this%start_ptr =  transfer(c_loc(this%mem(1)), this%start_ptr)
            this%end_ptr = this%start_ptr
        else
            stop "stack_memory_pool already initialized"
        end if
    end subroutine stack_memory_pool_create

    subroutine stack_memory_pool_destroy(this)
        type(stack_memory_pool) :: this

        if (this%start_ptr .ne. 0) then
            deallocate(this%mem)
            this%start_ptr = 0
            this%end_ptr = 0
        end if
    end subroutine stack_memory_pool_destroy

    function stack_memory_pool_allocate(this, byte_size) result(allocated_mem_ptr)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_int32_t), intent(in) :: byte_size
        integer(kind=c_intptr_t) :: allocated_mem_ptr
        integer(kind=c_int32_t) :: used_mem_size
        integer(kind=c_int32_t) :: free_mem_size
        integer(kind=c_int32_t) :: total_mem_size

        if (this%start_ptr .ne. 0) then
            total_mem_size = size(this%mem, 1, kind=c_int32_t)
            used_mem_size = this%end_ptr - this%start_ptr
            free_mem_size = total_mem_size - used_mem_size
            if (free_mem_size >=  byte_size) then
                allocated_mem_ptr = this%end_ptr
                this%end_ptr = this%end_ptr + byte_size
            else
                stop "not enough memory"
            end if
        else
            stop "stack_memory_pool not initialized"
        end if
    end function stack_memory_pool_allocate

    subroutine stack_memory_pool_deallocate(this, allocated_mem_ptr)
        class(stack_memory_pool), intent(inout) :: this
        integer(kind=c_intptr_t), intent(in) :: allocated_mem_ptr

        if (this%start_ptr .ne. 0) then
            if (allocated_mem_ptr >=this%start_ptr) then
                if (allocated_mem_ptr < this%end_ptr) then
                    this%end_ptr = allocated_mem_ptr
                else
                    !Unfortunately due to separated declaration and initialization
                    !stack deallocations can happen out of order in Fortran
                    !stop "address is within deallocated block"
                end if
            else
                stop "allocated_mem_ptr not in the memory pool"
            end if
        else
            stop "stack_memory_pool not initialized"
        end if
    end subroutine stack_memory_pool_deallocate

    subroutine sallocate_int_1d(stack_mem_pool, arr_ptr, size)
        class(stack_memory_pool), intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), dimension(:), pointer :: arr_ptr
        integer(kind=c_int32_t), intent(in) :: size
        type(c_ptr) :: ptr
        ptr = transfer(stack_mem_pool%allocate(size * c_int32_t), ptr)
        call c_f_pointer(ptr, arr_ptr, shape=[size])
    end subroutine sallocate_int_1d

    subroutine sdeallocate_int_1d(stack_mem_pool, arr_ptr)
        class(stack_memory_pool), intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), dimension(:), pointer, intent(in) :: arr_ptr
        integer(kind=c_intptr_t) :: ptr
        if (size(arr_ptr, 1) > 0) then
            ptr = transfer(c_loc(arr_ptr(1)), ptr)
            call stack_mem_pool%deallocate(ptr)
        end if
    end subroutine sdeallocate_int_1d

end module m_memory_pool

module m_array_int_1d
   use, intrinsic :: iso_c_binding
   use m_memory_pool

   implicit none

   private

    type, public :: array_int_1d
         integer(kind=c_int32_t), dimension(:), pointer, private :: arr_ptr => null()
         type(stack_memory_pool), pointer, private :: stack_mem_pool => null()
    contains
        procedure :: create => array_int_1d_create
        final :: array_int_1d_destroy
    end type array_int_1d

   contains

    function array_int_1d_create(this, stack_mem_pool, size) result(arr_ptr)
        class(array_int_1d), intent(inout) :: this
        class(stack_memory_pool), target, intent(inout) :: stack_mem_pool
        integer(kind=c_int32_t), intent(in) :: size
        integer(kind=c_int32_t), dimension(:), pointer :: arr_ptr

        if (.not. associated(this%arr_ptr)) then
            call sallocate_array(stack_mem_pool, this%arr_ptr, size)
            arr_ptr => this%arr_ptr
            this%stack_mem_pool => stack_mem_pool
        else
            stop "Array already allocated"
        end if
    end function array_int_1d_create

    subroutine array_int_1d_destroy(this)
        type(array_int_1d) :: this

        if (associated(this%arr_ptr)) then
            call sdeallocate_array(this%stack_mem_pool, this%arr_ptr)
        else
            stop "Array not allocated"
        end if
    end subroutine array_int_1d_destroy

end module m_array_int_1d

module m_mempool_tests
use, intrinsic :: iso_c_binding

implicit none

integer(kind=c_int32_t), parameter :: KB = 1024, ARR_SIZE = 10, N_REPEATS = 100

contains

subroutine test_auto_arrays()
    ! If compiler decides to allocate on the heap, there will be big performance penalty
    integer(kind=c_int32_t), dimension(ARR_SIZE):: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_auto_arrays start"

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c
    print *, "test_auto_arrays end"
end subroutine test_auto_arrays

subroutine test_stack_mempool_manual_dealloc(mem_pool)
    use m_memory_pool
    use m_array_int_1d

    type(stack_memory_pool), intent(inout) :: mem_pool
    integer(kind=c_int32_t), dimension(:), pointer :: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_stack_mempool_manual_dealloc start"

    call sallocate_array(mem_pool, a, ARR_SIZE)
    call sallocate_array(mem_pool, b, ARR_SIZE)
    call sallocate_array(mem_pool, c, ARR_SIZE)

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c
    ! It is dangerous to leave deallocation up to the programmer,
    ! as it could be forgotten
    call sdeallocate_array(mem_pool, c)
    call sdeallocate_array(mem_pool, b)
    call sdeallocate_array(mem_pool, a)

    print *, "test_stack_mempool_manual_dealloc end"
end subroutine test_stack_mempool_manual_dealloc

subroutine test_stack_mempool_auto_dealloc(mem_pool)
    use m_memory_pool
    use m_array_int_1d
    type(stack_memory_pool), intent(inout) :: mem_pool
    type(array_int_1d) :: a_, b_, c_
    integer(kind=c_int32_t), dimension(:), pointer :: a, b, c
    integer(kind=c_int32_t) :: i, j

    print *, "test_stack_mempool_auto_dealloc start"

    a => a_%create(mem_pool, ARR_SIZE)
    b => b_%create(mem_pool, ARR_SIZE)
    c => c_%create(mem_pool, ARR_SIZE)

    do i = 1, ARR_SIZE
        a(i) = i
        b(i) = ARR_SIZE - i
    end do

    c = a + b

    print *, a
    print *, b
    print *, c

    print *, "test_stack_mempool_auto_dealloc end"
end subroutine test_stack_mempool_auto_dealloc

end module m_mempool_tests


program mempool_test

use m_mempool_tests
use m_memory_pool

implicit none

type(stack_memory_pool) :: mem_pool
integer(kind=c_int32_t) :: i, j

!------------------------
do j = 1, N_REPEATS
    call test_auto_arrays()
end do
!------------------------
call mem_pool%create(1 * KB)
!------------------------
do j = 1, N_REPEATS
    call test_stack_mempool_manual_dealloc(mem_pool)
end do
!------------------------
do j = 1, N_REPEATS
    call test_stack_mempool_auto_dealloc(mem_pool)
end do

end program mempool_test