! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! DR_HOOK is a profiling and debugging system for the IFS, and should
! be called at the beginning and end of each subroutine. This is a
! dummy implementation for offline packages.

module omptimer

  use omp_lib

  public

  ! allocate two large arrays
  DOUBLE PRECISION, allocatable :: total_time(:,:)
  DOUBLE PRECISION, allocatable :: tstart(:,:)
  integer, allocatable :: ncalls(:)
  CHARACTER(len=80), allocatable :: names(:)
  integer, parameter :: hash_size = 1000
  integer :: omp_num_threads

contains

  SUBROUTINE char_to_hash(c, a)
    CHARACTER(len=*), INTENT(in) :: c
    INTEGER, INTENT(out) :: a
    INTEGER :: i

    INTEGER :: p = 31
    INTEGER :: m = hash_size - 3
    INTEGER :: p_pow

    a = 0
    p_pow = 1
    DO i=1,len_trim(c)
        a = MOD(a + (ichar(c(i:i)) + 1) * p_pow, m)
        p_pow = MOD(p_pow * p, m)
    END DO

  END SUBROUTINE

  subroutine omptimer_mark(proc_name, iswitch, proc_key)

    use parkind1, only : jprb

    character(len=*), intent(in)    :: proc_name
    integer,          intent(in)    :: iswitch
    real(jprb),       intent(inout) :: proc_key
    integer :: idx, thread_id
    double precision :: tstop 

   thread_id = omp_get_thread_num()

   if (iswitch == 0) then
      call char_to_hash(proc_name, idx)
      proc_key = real(idx)
      names(idx) = proc_name
      ncalls(idx) = ncalls(idx) + 1
      tstart(idx,thread_id) = omp_get_wtime()

    else if (iswitch == 1) then
      idx = INT(proc_key)
      tstop = omp_get_wtime()

      total_time(idx,thread_id) = total_time(idx,thread_id) + (tstop - tstart(idx, thread_id))
    endif

  end subroutine omptimer_mark

  subroutine omptimer_initialize()

    CHARACTER(len=255) :: ntstr
    CALL get_environment_variable("OMP_NUM_THREADS", ntstr)
    READ(ntstr, '(I2)') omp_num_threads
!    omp_num_threads = omp_get_num_threads()
 
    allocate(total_time(hash_size, 0:omp_num_threads-1))
    allocate(tstart(hash_size, 0:omp_num_threads-1))
    allocate(ncalls(hash_size))
    allocate(names(hash_size))

    total_time = 0
    ncalls = 0

  end subroutine omptimer_initialize

  subroutine omptimer_finalize()

    integer :: idx
    double precision :: ttime

    open(1, file="omptiming.txt", action="write")

    do idx = 1,1000
      ttime = 0
      if(total_time(idx,0) > 0.0) then
        ttime = sum(total_time(idx,:)) / omp_num_threads
        rms = sqrt( sum( (total_time(idx,:) - ttime)**2) / omp_num_threads)
        write(1, '(A80,E10.3,E10.3,I10)'), names(idx), ttime, rms, ncalls(idx)
      end if
    end do

    close(1)

  end subroutine omptimer_finalize

end module omptimer
