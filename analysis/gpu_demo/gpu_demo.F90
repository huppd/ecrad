! Created by  on 28.01.21.
! Author: Mikhail Zhigun

program gpu_demo
#ifdef _OPENACC
  use openacc
  use cudafor
implicit none
  call run()
#else
  print *, "Error: compiled without OpenACC"
  stop
#endif
contains
  subroutine run()
    integer, parameter :: DEV_NUM = 0 ! Always pick up first found device
    integer :: num_dev, istat
    character*100 :: dev_name
    integer(8) :: gpu_mem_size, gpu_mem_size_mb, num_mp, max_threads_per_mp, max_threads_per_block, warp_size, &
                  & max_threads_per_mp
    type (cudaDeviceProp) :: prop

    num_dev = acc_get_num_devices(ACC_DEVICE_NVIDIA)
    if (num_dev .lt. 1) then
      print *, 'Error: There are no devices available on this host'
      stop
    endif
    write(*,"(' Number of Nvidia GPU devices: ',i0)") num_dev
    call acc_get_property_string(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_NAME, dev_name)
    write(*,"(' Selecting device:             ',i0,' ""',a,'""')") DEV_NUM, dev_name
    call acc_set_device_num(DEV_NUM, acc_device_nvidia)
    gpu_mem_size = acc_get_property(DEV_NUM, ACC_DEVICE_NVIDIA, ACC_PROPERTY_MEMORY)
    gpu_mem_size_mb = gpu_mem_size / (1024*1024)
    write(*,"(' Device memory size:           ',i0,'mb')") gpu_mem_size_mb
    istat = cudaGetDeviceProperties(prop, DEV_NUM)
    num_mp = prop%multiProcessorCount
    max_threads_per_block = prop%maxThreadsPerBlock
    warp_size = prop%warpSize
    max_threads_per_mp = prop%maxThreadsPerMultiProcessor
    write(*,"(' Number of multiprocessors:    ',i0)") num_mp
    write(*,"(' Max threads per MP:           ',i0)") max_threads_per_mp
    write(*,"(' Max threads per block:        ',i0)") max_threads_per_block
    write(*,"(' Warp (OpenACC vector) size:   ',i0)") warp_size
  end subroutine run

end program gpu_demo