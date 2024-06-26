! Generated automatically.  DO NOT EDIT!

! integers

! unsigned
  integer(C_INT), parameter :: PFFT_TRANSPOSED_NONE = 0
  integer(C_INT), parameter :: PFFT_SHIFTED_NONE = 0
  integer(C_INT), parameter :: PFFT_MEASURE = 0
  integer(C_INT), parameter :: PFFT_NO_TUNE = 0
  integer(C_INT), parameter :: FPFFT_DEFAULT_BLOCKS = -1
  integer(C_INT), parameter :: FPFFT_NO_GCELLS = -1
  integer(C_INT), parameter :: PFFT_INT = 1
  integer(C_INT), parameter :: PFFT_PTRDIFF_T = 2
  integer(C_INT), parameter :: PFFT_FLOAT = 3
  integer(C_INT), parameter :: PFFT_DOUBLE = 4
  integer(C_INT), parameter :: PFFT_LDOUBLE = 5
  integer(C_INT), parameter :: PFFT_UNSIGNED = 6
  integer(C_INT), parameter :: PFFT_SWITCH = 7
  integer(C_INT), parameter :: PFFT_GC_TRANSPOSED_NONE = 0

! shifted unsigned
  integer(C_INT), parameter :: PFFT_TRANSPOSED_IN = 1
  integer(C_INT), parameter :: PFFT_TRANSPOSED_OUT = 2
  integer(C_INT), parameter :: PFFT_SHIFTED_IN = 4
  integer(C_INT), parameter :: PFFT_SHIFTED_OUT = 8
  integer(C_INT), parameter :: PFFT_ESTIMATE = 16
  integer(C_INT), parameter :: PFFT_PATIENT = 32
  integer(C_INT), parameter :: PFFT_EXHAUSTIVE = 64
  integer(C_INT), parameter :: PFFT_TUNE = 128
  integer(C_INT), parameter :: PFFT_PRESERVE_INPUT = 256
  integer(C_INT), parameter :: PFFT_DESTROY_INPUT = 512
  integer(C_INT), parameter :: PFFT_BUFFERED_INPLACE = 1024
  integer(C_INT), parameter :: PFFT_PADDED_R2C = 2048
  integer(C_INT), parameter :: PFFT_GC_TRANSPOSED = 1
  integer(C_INT), parameter :: PFFT_GC_SENDRECV = 2
  integer(C_INT), parameter :: PFFT_GC_RMA = 4
  integer(C_INT), parameter :: PFFT_GC_R2C = 8
  integer(C_INT), parameter :: PFFT_GC_PADDED = 16

! redirections
  integer(C_INT), parameter :: PFFT_PADDED_C2R = PFFT_PADDED_R2C
  integer(C_INT), parameter :: PFFT_GC_C2R = PFFT_GC_R2C
  integer(C_INT), parameter :: PFFT_GC_PADDED_R2C = PFFT_GC_R2C  &
      + PFFT_GC_PADDED
  integer(C_INT), parameter :: PFFT_GC_PADDED_C2R = PFFT_GC_C2R  &
      + PFFT_GC_PADDED
  integer(C_INT), parameter :: PFFT_R2HC = FFTW_R2HC
  integer(C_INT), parameter :: PFFT_HC2R = FFTW_HC2R
  integer(C_INT), parameter :: PFFT_DHT = FFTW_DHT
  integer(C_INT), parameter :: PFFT_REDFT00 = FFTW_REDFT00
  integer(C_INT), parameter :: PFFT_REDFT01 = FFTW_REDFT01
  integer(C_INT), parameter :: PFFT_REDFT10 = FFTW_REDFT10
  integer(C_INT), parameter :: PFFT_REDFT11 = FFTW_REDFT11
  integer(C_INT), parameter :: PFFT_RODFT00 = FFTW_RODFT00
  integer(C_INT), parameter :: PFFT_RODFT01 = FFTW_RODFT01
  integer(C_INT), parameter :: PFFT_RODFT10 = FFTW_RODFT10
  integer(C_INT), parameter :: PFFT_RODFT11 = FFTW_RODFT11
  integer(C_INT), parameter :: PFFT_FORWARD = FFTW_FORWARD
  integer(C_INT), parameter :: PFFT_BACKWARD = FFTW_BACKWARD
  integer(C_INT), parameter :: PFFT_DEFAULT_BLOCK = FFTW_MPI_DEFAULT_BLOCK


  interface
    subroutine pfft_init() bind(C, name='pfft_init')
      import
    end subroutine pfft_init
    
    subroutine pfft_cleanup() bind(C, name='pfft_cleanup')
      import
    end subroutine pfft_cleanup
    
    type(C_PTR) function pfft_malloc(Nos) bind(C, name='pfft_malloc')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfft_malloc
    
    type(C_PTR) function pfft_alloc_real(Nos) bind(C, name='pfft_alloc_real')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfft_alloc_real
    
    type(C_PTR) function pfft_alloc_complex(Nos) bind(C, name='pfft_alloc_complex')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfft_alloc_complex
    
    subroutine pfft_free(p) bind(C, name='pfft_free')
      import
      type(C_PTR), value :: p
    end subroutine pfft_free
    
    subroutine pfft_execute(plan) bind(C, name='pfft_execute')
      import
      type(C_PTR), value :: plan
    end subroutine pfft_execute
    
    subroutine pfft_execute_dft(plan,in,out) bind(C, name='pfft_execute_dft')
      import
      type(C_PTR), value :: plan
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine pfft_execute_dft
    
    subroutine pfft_execute_dft_r2c(plan,in,out) bind(C, name='pfft_execute_dft_r2c')
      import
      type(C_PTR), value :: plan
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine pfft_execute_dft_r2c
    
    subroutine pfft_execute_dft_c2r(plan,in,out) bind(C, name='pfft_execute_dft_c2r')
      import
      type(C_PTR), value :: plan
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine pfft_execute_dft_c2r
    
    subroutine pfft_execute_r2r(plan,in,out) bind(C, name='pfft_execute_r2r')
      import
      type(C_PTR), value :: plan
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine pfft_execute_r2r
    
    subroutine pfft_destroy_plan(plan) bind(C, name='pfft_destroy_plan')
      import
      type(C_PTR), value :: plan
    end subroutine pfft_destroy_plan
    
    subroutine pfft_init_input_complex_3d(Nos,local_n,local_n_start,data) bind(C, name='pfft_init_input_complex_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_complex_3d
    
    subroutine pfft_init_input_complex(rnk_n,Nos,local_n,local_start,data) bind(C, name='pfft_init_input_complex')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_complex
    
    subroutine pfft_init_input_complex_hermitian_3d(Nos,local_n,local_n_start,data) &
               bind(C, name='pfft_init_input_complex_hermitian_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_complex_hermitian_3d
    
    subroutine pfft_init_input_complex_hermitian(rnk_n,Nos,local_n,local_start,data) &
               bind(C, name='pfft_init_input_complex_hermitian')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_complex_hermitian
    
    subroutine pfft_init_input_real_3d(Nos,local_n,local_n_start,data) bind(C, name='pfft_init_input_real_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      real(C_DOUBLE), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_real_3d
    
    subroutine pfft_init_input_real(rnk_n,Nos,local_n,local_start,data) bind(C, name='pfft_init_input_real')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      real(C_DOUBLE), dimension(*), intent(out) :: data
    end subroutine pfft_init_input_real
    
    real(C_DOUBLE) function pfft_check_output_complex_3d(Nos,local_n,local_n_start,data,comm) &
                            bind(C, name='pfft_check_output_complex_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_complex_3d
    
    real(C_DOUBLE) function pfft_check_output_complex(rnk_n,Nos,local_n,local_start,data,comm) &
                            bind(C, name='pfft_check_output_complex_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_complex
    
    real(C_DOUBLE) function pfft_check_output_complex_hermitian_3d(Nos,local_n,local_n_start,data,comm) &
                            bind(C, name='pfft_check_output_complex_hermitian_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_complex_hermitian_3d
    
    real(C_DOUBLE) function pfft_check_output_complex_hermitian(rnk_n,Nos,local_n,local_start,data,comm) &
                            bind(C, name='pfft_check_output_complex_hermitian_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_complex_hermitian
    
    real(C_DOUBLE) function pfft_check_output_real_3d(Nos,local_n,local_n_start,data,comm) &
                            bind(C, name='pfft_check_output_real_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      real(C_DOUBLE), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_real_3d
    
    real(C_DOUBLE) function pfft_check_output_real(rnk_n,Nos,local_n,local_start,data,comm) &
                            bind(C, name='pfft_check_output_real_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      real(C_DOUBLE), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfft_check_output_real
    
    integer(C_INTPTR_T) function pfft_local_size_dft_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft_3d
    
    integer(C_INTPTR_T) function pfft_local_size_dft_r2c_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                            local_o_start) bind(C, name='pfft_local_size_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft_r2c_3d
    
    integer(C_INTPTR_T) function pfft_local_size_dft_c2r_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                            local_o_start) bind(C, name='pfft_local_size_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft_c2r_3d
    
    integer(C_INTPTR_T) function pfft_local_size_r2r_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_r2r_3d
    
    integer(C_INTPTR_T) function pfft_local_size_dft(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft
    
    integer(C_INTPTR_T) function pfft_local_size_dft_r2c(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                         local_o_start) bind(C, name='pfft_local_size_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft_r2c
    
    integer(C_INTPTR_T) function pfft_local_size_dft_c2r(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                         local_o_start) bind(C, name='pfft_local_size_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_dft_c2r
    
    integer(C_INTPTR_T) function pfft_local_size_r2r(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_r2r
    
    integer(C_INTPTR_T) function pfft_local_size_many_dft(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                          local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_many_dft
    
    integer(C_INTPTR_T) function pfft_local_size_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                              local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_many_dft_r2c
    
    integer(C_INTPTR_T) function pfft_local_size_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                              local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_many_dft_c2r
    
    integer(C_INTPTR_T) function pfft_local_size_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                          local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfft_local_size_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfft_local_size_many_r2r
    
    subroutine pfft_local_block_dft_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft_3d
    
    subroutine pfft_local_block_dft_r2c_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft_r2c_3d
    
    subroutine pfft_local_block_dft_c2r_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft_c2r_3d
    
    subroutine pfft_local_block_r2r_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_r2r_3d
    
    subroutine pfft_local_block_dft(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft
    
    subroutine pfft_local_block_dft_r2c(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft_r2c
    
    subroutine pfft_local_block_dft_c2r(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_dft_c2r
    
    subroutine pfft_local_block_r2r(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfft_local_block_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_r2r
    
    subroutine pfft_local_block_many_dft(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                         local_o_start) bind(C, name='pfft_local_block_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_many_dft
    
    subroutine pfft_local_block_many_dft_r2c(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                             local_o_start) bind(C, name='pfft_local_block_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_many_dft_r2c
    
    subroutine pfft_local_block_many_dft_c2r(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                             local_o_start) bind(C, name='pfft_local_block_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_many_dft_c2r
    
    subroutine pfft_local_block_many_r2r(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                         local_o_start) bind(C, name='pfft_local_block_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfft_local_block_many_r2r
    
    type(C_PTR) function pfft_plan_dft_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft_3d
    
    type(C_PTR) function pfft_plan_dft_r2c_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft_r2c_3d
    
    type(C_PTR) function pfft_plan_dft_c2r_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft_c2r_3d
    
    type(C_PTR) function pfft_plan_r2r_3d(Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfft_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_r2r_3d
    
    type(C_PTR) function pfft_plan_dft(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft
    
    type(C_PTR) function pfft_plan_dft_r2c(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft_r2c
    
    type(C_PTR) function pfft_plan_dft_c2r(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_dft_c2r
    
    type(C_PTR) function pfft_plan_r2r(rnk_n,Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfft_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_r2r
    
    type(C_PTR) function pfft_plan_many_dft(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfft_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft
    
    type(C_PTR) function pfft_plan_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfft_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft_r2c
    
    type(C_PTR) function pfft_plan_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfft_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft_c2r
    
    type(C_PTR) function pfft_plan_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,kinds,pfft_flags) &
                         bind(C, name='pfft_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_r2r
    
    type(C_PTR) function pfft_plan_many_dft_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                    pfft_flags) bind(C, name='pfft_plan_many_dft_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft_skipped
    
    type(C_PTR) function pfft_plan_many_dft_r2c_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                        pfft_flags) bind(C, name='pfft_plan_many_dft_r2c_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft_r2c_skipped
    
    type(C_PTR) function pfft_plan_many_dft_c2r_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                        pfft_flags) bind(C, name='pfft_plan_many_dft_c2r_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_dft_c2r_skipped
    
    type(C_PTR) function pfft_plan_many_r2r_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,kinds, &
                                                    pfft_flags) bind(C, name='pfft_plan_many_r2r_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_many_r2r_skipped
    
    integer(C_INTPTR_T) function pfft_prod_INT(d,vec) bind(C, name='pfft_prod_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec
    end function pfft_prod_INT
    
    integer(C_INTPTR_T) function pfft_sum_INT(d,vec) bind(C, name='pfft_sum_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec
    end function pfft_sum_INT
    
    integer(C_INT) function pfft_equal_INT(d,vec1,vec2) bind(C, name='pfft_equal_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
    end function pfft_equal_INT
    
    subroutine pfft_vcopy_INT(d,vec1,vec2) bind(C, name='pfft_vcopy_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(out) :: vec2
    end subroutine pfft_vcopy_INT
    
    subroutine pfft_vadd_INT(d,vec1,vec2,sum) bind(C, name='pfft_vadd_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
      integer(C_INTPTR_T), dimension(*), intent(out) :: sum
    end subroutine pfft_vadd_INT
    
    subroutine pfft_vsub_INT(d,vec1,vec2,sum) bind(C, name='pfft_vsub_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
      integer(C_INTPTR_T), dimension(*), intent(out) :: sum
    end subroutine pfft_vsub_INT
    
    subroutine pfft_apr_complex_3d(data,local_n,local_start,name,comm) bind(C, name='pfft_apr_complex_3d_f03')
      import
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_apr_complex_3d
    
    subroutine pfft_apr_complex_permuted_3d(data,local_n,local_start,perm0,perm1,perm2,name,comm) &
               bind(C, name='pfft_apr_complex_permuted_3d_f03')
      import
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INT), value :: perm0
      integer(C_INT), value :: perm1
      integer(C_INT), value :: perm2
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_apr_complex_permuted_3d
    
    subroutine pfft_apr_real_3d(data,local_n,local_start,name,comm) bind(C, name='pfft_apr_real_3d_f03')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_apr_real_3d
    
    subroutine pfft_apr_real_permuted_3d(data,local_n,local_start,perm0,perm1,perm2,name,comm) &
               bind(C, name='pfft_apr_real_permuted_3d_f03')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INT), value :: perm0
      integer(C_INT), value :: perm1
      integer(C_INT), value :: perm2
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_apr_real_permuted_3d
    
    subroutine pfft_reset_timer(ths) bind(C, name='pfft_reset_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_reset_timer
    
    type(C_PTR) function pfft_get_timer(ths) bind(C, name='pfft_get_timer')
      import
      type(C_PTR), value :: ths
    end function pfft_get_timer
    
    subroutine pfft_print_average_timer(ths,comm) bind(C, name='pfft_print_average_timer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfft_print_average_timer
    
    subroutine pfft_print_average_timer_adv(ths,comm) bind(C, name='pfft_print_average_timer_adv_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfft_print_average_timer_adv
    
    subroutine pfft_write_average_timer(ths,name,comm) bind(C, name='pfft_write_average_timer_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_write_average_timer
    
    subroutine pfft_write_average_timer_adv(ths,name,comm) bind(C, name='pfft_write_average_timer_adv_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_write_average_timer_adv
    
    type(C_PTR) function pfft_copy_timer(orig) bind(C, name='pfft_copy_timer')
      import
      type(C_PTR), value :: orig
    end function pfft_copy_timer
    
    subroutine pfft_average_timer(ths) bind(C, name='pfft_average_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_average_timer
    
    type(C_PTR) function pfft_add_timers(sum1,sum2) bind(C, name='pfft_add_timers')
      import
      type(C_PTR), value :: sum1
      type(C_PTR), value :: sum2
    end function pfft_add_timers
    
    type(C_PTR) function pfft_reduce_max_timer(ths,comm) bind(C, name='pfft_reduce_max_timer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end function pfft_reduce_max_timer
    
    type(C_PTR) function pfft_convert_timer2vec(ths) bind(C, name='pfft_convert_timer2vec')
      import
      type(C_PTR), value :: ths
    end function pfft_convert_timer2vec
    
    type(C_PTR) function pfft_convert_vec2timer(times) bind(C, name='pfft_convert_vec2timer')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: times
    end function pfft_convert_vec2timer
    
    subroutine pfft_destroy_timer(ths) bind(C, name='pfft_destroy_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_destroy_timer
    
    integer(C_INT) function pfft_create_procmesh(rnk_n,comm,np,comm_cart) bind(C, name='pfft_create_procmesh_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INT32_T), value :: comm
      integer(C_INT), dimension(*), intent(in) :: np
      integer(C_INT32_T), intent(out) :: comm_cart
    end function pfft_create_procmesh
    
    integer(C_INT) function pfft_create_procmesh_1d(comm,np0,comm_cart_1d) bind(C, name='pfft_create_procmesh_1d_f03')
      import
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: np0
      integer(C_INT32_T), intent(out) :: comm_cart_1d
    end function pfft_create_procmesh_1d
    
    integer(C_INT) function pfft_create_procmesh_2d(comm,np0,np1,comm_cart_2d) bind(C, name='pfft_create_procmesh_2d_f03')
      import
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: np0
      integer(C_INT), value :: np1
      integer(C_INT32_T), intent(out) :: comm_cart_2d
    end function pfft_create_procmesh_2d
    
    integer(C_INTPTR_T) function pfft_local_size_gc_3d(local_n,local_start,gc_below,gc_above,local_ngc,local_gc_start) &
                                 bind(C, name='pfft_local_size_gc_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfft_local_size_gc_3d
    
    integer(C_INTPTR_T) function pfft_local_size_gc(rnk_n,local_n,local_start,gc_below,gc_above,local_ngc,local_gc_start) &
                                 bind(C, name='pfft_local_size_gc')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfft_local_size_gc
    
    integer(C_INTPTR_T) function pfft_local_size_many_gc(rnk_n,local_n,local_start,howmany,gc_below,gc_above,local_ngc, &
                                                         local_gc_start) bind(C, name='pfft_local_size_many_gc')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfft_local_size_many_gc
    
    type(C_PTR) function pfft_plan_rgc_3d(Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfft_plan_rgc_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_DOUBLE), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_rgc_3d
    
    type(C_PTR) function pfft_plan_cgc_3d(Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfft_plan_cgc_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_cgc_3d
    
    type(C_PTR) function pfft_plan_rgc(rnk_n,Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfft_plan_rgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_DOUBLE), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_rgc
    
    type(C_PTR) function pfft_plan_cgc(rnk_n,Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfft_plan_cgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_cgc
    
    type(C_PTR) function pfft_plan_many_rgc(rnk_n,Nos,howmany,block,gc_below,gc_above,data,comm_cart,gc_flags) &
                         bind(C, name='pfft_plan_many_rgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: block
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_DOUBLE), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_many_rgc
    
    type(C_PTR) function pfft_plan_many_cgc(rnk_n,Nos,howmany,block,gc_below,gc_above,data,comm_cart,gc_flags) &
                         bind(C, name='pfft_plan_many_cgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: block
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfft_plan_many_cgc
    
    subroutine pfft_exchange(ths) bind(C, name='pfft_exchange')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_exchange
    
    subroutine pfft_reduce(ths) bind(C, name='pfft_reduce')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_reduce
    
    subroutine pfft_destroy_gcplan(ths) bind(C, name='pfft_destroy_gcplan')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_destroy_gcplan
    
    subroutine pfft_reset_gctimers(ths) bind(C, name='pfft_reset_gctimers')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_reset_gctimers
    
    type(C_PTR) function pfft_get_gctimer_exg(ths) bind(C, name='pfft_get_gctimer_exg')
      import
      type(C_PTR), value :: ths
    end function pfft_get_gctimer_exg
    
    type(C_PTR) function pfft_get_gctimer_red(ths) bind(C, name='pfft_get_gctimer_red')
      import
      type(C_PTR), value :: ths
    end function pfft_get_gctimer_red
    
    subroutine pfft_print_average_gctimer(ths,comm) bind(C, name='pfft_print_average_gctimer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfft_print_average_gctimer
    
    subroutine pfft_print_average_gctimer_adv(ths,comm) bind(C, name='pfft_print_average_gctimer_adv_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfft_print_average_gctimer_adv
    
    subroutine pfft_write_average_gctimer(ths,name,comm) bind(C, name='pfft_write_average_gctimer_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_write_average_gctimer
    
    subroutine pfft_write_average_gctimer_adv(ths,name,comm) bind(C, name='pfft_write_average_gctimer_adv_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfft_write_average_gctimer_adv
    
    type(C_PTR) function pfft_copy_gctimer(orig) bind(C, name='pfft_copy_gctimer')
      import
      type(C_PTR), value :: orig
    end function pfft_copy_gctimer
    
    subroutine pfft_average_gctimer(ths) bind(C, name='pfft_average_gctimer')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_average_gctimer
    
    type(C_PTR) function pfft_add_gctimers(sum1,sum2) bind(C, name='pfft_add_gctimers')
      import
      type(C_PTR), value :: sum1
      type(C_PTR), value :: sum2
    end function pfft_add_gctimers
    
    type(C_PTR) function pfft_reduce_max_gctimer(ths,comm) bind(C, name='pfft_reduce_max_gctimer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end function pfft_reduce_max_gctimer
    
    subroutine pfft_convert_gctimer2vec(ths,times) bind(C, name='pfft_convert_gctimer2vec')
      import
      type(C_PTR), value :: ths
      real(C_DOUBLE), dimension(*), intent(out) :: times
    end subroutine pfft_convert_gctimer2vec
    
    type(C_PTR) function pfft_convert_vec2gctimer(times) bind(C, name='pfft_convert_vec2gctimer')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: times
    end function pfft_convert_vec2gctimer
    
    subroutine pfft_destroy_gctimer(ths) bind(C, name='pfft_destroy_gctimer')
      import
      type(C_PTR), value :: ths
    end subroutine pfft_destroy_gctimer
    
  end interface


  interface
    subroutine pfftf_init() bind(C, name='pfftf_init')
      import
    end subroutine pfftf_init
    
    subroutine pfftf_cleanup() bind(C, name='pfftf_cleanup')
      import
    end subroutine pfftf_cleanup
    
    type(C_PTR) function pfftf_malloc(Nos) bind(C, name='pfftf_malloc')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfftf_malloc
    
    type(C_PTR) function pfftf_alloc_real(Nos) bind(C, name='pfftf_alloc_real')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfftf_alloc_real
    
    type(C_PTR) function pfftf_alloc_complex(Nos) bind(C, name='pfftf_alloc_complex')
      import
      integer(C_SIZE_T), value :: Nos
    end function pfftf_alloc_complex
    
    subroutine pfftf_free(p) bind(C, name='pfftf_free')
      import
      type(C_PTR), value :: p
    end subroutine pfftf_free
    
    subroutine pfftf_execute(plan) bind(C, name='pfftf_execute')
      import
      type(C_PTR), value :: plan
    end subroutine pfftf_execute
    
    subroutine pfftf_execute_dft(plan,in,out) bind(C, name='pfftf_execute_dft')
      import
      type(C_PTR), value :: plan
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine pfftf_execute_dft
    
    subroutine pfftf_execute_dft_r2c(plan,in,out) bind(C, name='pfftf_execute_dft_r2c')
      import
      type(C_PTR), value :: plan
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine pfftf_execute_dft_r2c
    
    subroutine pfftf_execute_dft_c2r(plan,in,out) bind(C, name='pfftf_execute_dft_c2r')
      import
      type(C_PTR), value :: plan
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine pfftf_execute_dft_c2r
    
    subroutine pfftf_execute_r2r(plan,in,out) bind(C, name='pfftf_execute_r2r')
      import
      type(C_PTR), value :: plan
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine pfftf_execute_r2r
    
    subroutine pfftf_destroy_plan(plan) bind(C, name='pfftf_destroy_plan')
      import
      type(C_PTR), value :: plan
    end subroutine pfftf_destroy_plan
    
    subroutine pfftf_init_input_complex_3d(Nos,local_n,local_n_start,data) bind(C, name='pfftf_init_input_complex_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_complex_3d
    
    subroutine pfftf_init_input_complex(rnk_n,Nos,local_n,local_start,data) bind(C, name='pfftf_init_input_complex')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_complex
    
    subroutine pfftf_init_input_complex_hermitian_3d(Nos,local_n,local_n_start,data) &
               bind(C, name='pfftf_init_input_complex_hermitian_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_complex_hermitian_3d
    
    subroutine pfftf_init_input_complex_hermitian(rnk_n,Nos,local_n,local_start,data) &
               bind(C, name='pfftf_init_input_complex_hermitian')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_complex_hermitian
    
    subroutine pfftf_init_input_real_3d(Nos,local_n,local_n_start,data) bind(C, name='pfftf_init_input_real_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      real(C_FLOAT), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_real_3d
    
    subroutine pfftf_init_input_real(rnk_n,Nos,local_n,local_start,data) bind(C, name='pfftf_init_input_real')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      real(C_FLOAT), dimension(*), intent(out) :: data
    end subroutine pfftf_init_input_real
    
    real(C_FLOAT) function pfftf_check_output_complex_3d(Nos,local_n,local_n_start,data,comm) &
                           bind(C, name='pfftf_check_output_complex_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_complex_3d
    
    real(C_FLOAT) function pfftf_check_output_complex(rnk_n,Nos,local_n,local_start,data,comm) &
                           bind(C, name='pfftf_check_output_complex_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_complex
    
    real(C_FLOAT) function pfftf_check_output_complex_hermitian_3d(Nos,local_n,local_n_start,data,comm) &
                           bind(C, name='pfftf_check_output_complex_hermitian_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_complex_hermitian_3d
    
    real(C_FLOAT) function pfftf_check_output_complex_hermitian(rnk_n,Nos,local_n,local_start,data,comm) &
                           bind(C, name='pfftf_check_output_complex_hermitian_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_complex_hermitian
    
    real(C_FLOAT) function pfftf_check_output_real_3d(Nos,local_n,local_n_start,data,comm) &
                           bind(C, name='pfftf_check_output_real_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n_start
      real(C_FLOAT), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_real_3d
    
    real(C_FLOAT) function pfftf_check_output_real(rnk_n,Nos,local_n,local_start,data,comm) &
                           bind(C, name='pfftf_check_output_real_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      real(C_FLOAT), dimension(*), intent(in) :: data
      integer(C_INT32_T), value :: comm
    end function pfftf_check_output_real
    
    integer(C_INTPTR_T) function pfftf_local_size_dft_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft_3d
    
    integer(C_INTPTR_T) function pfftf_local_size_dft_r2c_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                             local_o_start) bind(C, name='pfftf_local_size_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft_r2c_3d
    
    integer(C_INTPTR_T) function pfftf_local_size_dft_c2r_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                             local_o_start) bind(C, name='pfftf_local_size_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft_c2r_3d
    
    integer(C_INTPTR_T) function pfftf_local_size_r2r_3d(Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_r2r_3d
    
    integer(C_INTPTR_T) function pfftf_local_size_dft(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                      local_o_start) bind(C, name='pfftf_local_size_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft
    
    integer(C_INTPTR_T) function pfftf_local_size_dft_r2c(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                          local_o_start) bind(C, name='pfftf_local_size_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft_r2c
    
    integer(C_INTPTR_T) function pfftf_local_size_dft_c2r(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                          local_o_start) bind(C, name='pfftf_local_size_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_dft_c2r
    
    integer(C_INTPTR_T) function pfftf_local_size_r2r(rnk_n,Nos,comm_cart,pfft_flags,local_ni,local_i_start,local_no, &
                                                      local_o_start) bind(C, name='pfftf_local_size_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_r2r
    
    integer(C_INTPTR_T) function pfftf_local_size_many_dft(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                           local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_many_dft
    
    integer(C_INTPTR_T) function pfftf_local_size_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags, &
                                                               local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_many_dft_r2c
    
    integer(C_INTPTR_T) function pfftf_local_size_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags, &
                                                               local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_many_dft_c2r
    
    integer(C_INTPTR_T) function pfftf_local_size_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,comm_cart,pfft_flags,local_ni, &
                                                           local_i_start,local_no,local_o_start) &
                                 bind(C, name='pfftf_local_size_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end function pfftf_local_size_many_r2r
    
    subroutine pfftf_local_block_dft_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft_3d
    
    subroutine pfftf_local_block_dft_r2c_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft_r2c_3d
    
    subroutine pfftf_local_block_dft_c2r_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft_c2r_3d
    
    subroutine pfftf_local_block_r2r_3d(Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_r2r_3d
    
    subroutine pfftf_local_block_dft(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft
    
    subroutine pfftf_local_block_dft_r2c(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft_r2c
    
    subroutine pfftf_local_block_dft_c2r(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_dft_c2r
    
    subroutine pfftf_local_block_r2r(rnk_n,Nos,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no,local_o_start) &
               bind(C, name='pfftf_local_block_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_r2r
    
    subroutine pfftf_local_block_many_dft(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                          local_o_start) bind(C, name='pfftf_local_block_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_many_dft
    
    subroutine pfftf_local_block_many_dft_r2c(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                              local_o_start) bind(C, name='pfftf_local_block_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_many_dft_r2c
    
    subroutine pfftf_local_block_many_dft_c2r(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                              local_o_start) bind(C, name='pfftf_local_block_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_many_dft_c2r
    
    subroutine pfftf_local_block_many_r2r(rnk_n,ni,no,iblock,oblock,comm_cart,pid,pfft_flags,local_ni,local_i_start,local_no, &
                                          local_o_start) bind(C, name='pfftf_local_block_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: pid
      integer(C_INT), value :: pfft_flags
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ni
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_i_start
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_no
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_o_start
    end subroutine pfftf_local_block_many_r2r
    
    type(C_PTR) function pfftf_plan_dft_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft_3d
    
    type(C_PTR) function pfftf_plan_dft_r2c_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft_r2c_3d
    
    type(C_PTR) function pfftf_plan_dft_c2r_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft_c2r_3d
    
    type(C_PTR) function pfftf_plan_r2r_3d(Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfftf_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_r2r_3d
    
    type(C_PTR) function pfftf_plan_dft(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft
    
    type(C_PTR) function pfftf_plan_dft_r2c(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft_r2c
    
    type(C_PTR) function pfftf_plan_dft_c2r(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_dft_c2r
    
    type(C_PTR) function pfftf_plan_r2r(rnk_n,Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfftf_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_r2r
    
    type(C_PTR) function pfftf_plan_many_dft(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfftf_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft
    
    type(C_PTR) function pfftf_plan_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfftf_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft_r2c
    
    type(C_PTR) function pfftf_plan_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
                         bind(C, name='pfftf_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft_c2r
    
    type(C_PTR) function pfftf_plan_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,kinds,pfft_flags) &
                         bind(C, name='pfftf_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_r2r
    
    type(C_PTR) function pfftf_plan_many_dft_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                     pfft_flags) bind(C, name='pfftf_plan_many_dft_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft_skipped
    
    type(C_PTR) function pfftf_plan_many_dft_r2c_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                         pfft_flags) bind(C, name='pfftf_plan_many_dft_r2c_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft_r2c_skipped
    
    type(C_PTR) function pfftf_plan_many_dft_c2r_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
                                                         pfft_flags) bind(C, name='pfftf_plan_many_dft_c2r_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: sign
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_dft_c2r_skipped
    
    type(C_PTR) function pfftf_plan_many_r2r_skipped(rnk_n,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,kinds, &
                                                     pfft_flags) bind(C, name='pfftf_plan_many_r2r_skipped_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: ni
      integer(C_INTPTR_T), dimension(*), intent(in) :: no
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
      integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
      integer(C_INT), dimension(*), intent(in) :: skip_trafos
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_many_r2r_skipped
    
    integer(C_INTPTR_T) function pfftf_prod_INT(d,vec) bind(C, name='pfftf_prod_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec
    end function pfftf_prod_INT
    
    integer(C_INTPTR_T) function pfftf_sum_INT(d,vec) bind(C, name='pfftf_sum_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec
    end function pfftf_sum_INT
    
    integer(C_INT) function pfftf_equal_INT(d,vec1,vec2) bind(C, name='pfftf_equal_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
    end function pfftf_equal_INT
    
    subroutine pfftf_vcopy_INT(d,vec1,vec2) bind(C, name='pfftf_vcopy_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(out) :: vec2
    end subroutine pfftf_vcopy_INT
    
    subroutine pfftf_vadd_INT(d,vec1,vec2,sum) bind(C, name='pfftf_vadd_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
      integer(C_INTPTR_T), dimension(*), intent(out) :: sum
    end subroutine pfftf_vadd_INT
    
    subroutine pfftf_vsub_INT(d,vec1,vec2,sum) bind(C, name='pfftf_vsub_INT')
      import
      integer(C_INT), value :: d
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec1
      integer(C_INTPTR_T), dimension(*), intent(in) :: vec2
      integer(C_INTPTR_T), dimension(*), intent(out) :: sum
    end subroutine pfftf_vsub_INT
    
    subroutine pfftf_apr_complex_3d(data,local_n,local_start,name,comm) bind(C, name='pfftf_apr_complex_3d_f03')
      import
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_apr_complex_3d
    
    subroutine pfftf_apr_complex_permuted_3d(data,local_n,local_start,perm0,perm1,perm2,name,comm) &
               bind(C, name='pfftf_apr_complex_permuted_3d_f03')
      import
      complex(C_FLOAT_COMPLEX), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INT), value :: perm0
      integer(C_INT), value :: perm1
      integer(C_INT), value :: perm2
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_apr_complex_permuted_3d
    
    subroutine pfftf_apr_real_3d(data,local_n,local_start,name,comm) bind(C, name='pfftf_apr_real_3d_f03')
      import
      real(C_FLOAT), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_apr_real_3d
    
    subroutine pfftf_apr_real_permuted_3d(data,local_n,local_start,perm0,perm1,perm2,name,comm) &
               bind(C, name='pfftf_apr_real_permuted_3d_f03')
      import
      real(C_FLOAT), dimension(*), intent(in) :: data
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INT), value :: perm0
      integer(C_INT), value :: perm1
      integer(C_INT), value :: perm2
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_apr_real_permuted_3d
    
    subroutine pfftf_reset_timer(ths) bind(C, name='pfftf_reset_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_reset_timer
    
    type(C_PTR) function pfftf_get_timer(ths) bind(C, name='pfftf_get_timer')
      import
      type(C_PTR), value :: ths
    end function pfftf_get_timer
    
    subroutine pfftf_print_average_timer(ths,comm) bind(C, name='pfftf_print_average_timer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_print_average_timer
    
    subroutine pfftf_print_average_timer_adv(ths,comm) bind(C, name='pfftf_print_average_timer_adv_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_print_average_timer_adv
    
    subroutine pfftf_write_average_timer(ths,name,comm) bind(C, name='pfftf_write_average_timer_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_write_average_timer
    
    subroutine pfftf_write_average_timer_adv(ths,name,comm) bind(C, name='pfftf_write_average_timer_adv_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_write_average_timer_adv
    
    type(C_PTR) function pfftf_copy_timer(orig) bind(C, name='pfftf_copy_timer')
      import
      type(C_PTR), value :: orig
    end function pfftf_copy_timer
    
    subroutine pfftf_average_timer(ths) bind(C, name='pfftf_average_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_average_timer
    
    type(C_PTR) function pfftf_add_timers(sum1,sum2) bind(C, name='pfftf_add_timers')
      import
      type(C_PTR), value :: sum1
      type(C_PTR), value :: sum2
    end function pfftf_add_timers
    
    type(C_PTR) function pfftf_reduce_max_timer(ths,comm) bind(C, name='pfftf_reduce_max_timer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end function pfftf_reduce_max_timer
    
    type(C_PTR) function pfftf_convert_timer2vec(ths) bind(C, name='pfftf_convert_timer2vec')
      import
      type(C_PTR), value :: ths
    end function pfftf_convert_timer2vec
    
    type(C_PTR) function pfftf_convert_vec2timer(times) bind(C, name='pfftf_convert_vec2timer')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: times
    end function pfftf_convert_vec2timer
    
    subroutine pfftf_destroy_timer(ths) bind(C, name='pfftf_destroy_timer')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_destroy_timer
    
    integer(C_INT) function pfftf_create_procmesh(rnk_n,comm,np,comm_cart) bind(C, name='pfftf_create_procmesh_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INT32_T), value :: comm
      integer(C_INT), dimension(*), intent(in) :: np
      integer(C_INT32_T), intent(out) :: comm_cart
    end function pfftf_create_procmesh
    
    integer(C_INT) function pfftf_create_procmesh_1d(comm,np0,comm_cart_1d) bind(C, name='pfftf_create_procmesh_1d_f03')
      import
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: np0
      integer(C_INT32_T), intent(out) :: comm_cart_1d
    end function pfftf_create_procmesh_1d
    
    integer(C_INT) function pfftf_create_procmesh_2d(comm,np0,np1,comm_cart_2d) bind(C, name='pfftf_create_procmesh_2d_f03')
      import
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: np0
      integer(C_INT), value :: np1
      integer(C_INT32_T), intent(out) :: comm_cart_2d
    end function pfftf_create_procmesh_2d
    
    integer(C_INTPTR_T) function pfftf_local_size_gc_3d(local_n,local_start,gc_below,gc_above,local_ngc,local_gc_start) &
                                 bind(C, name='pfftf_local_size_gc_3d')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfftf_local_size_gc_3d
    
    integer(C_INTPTR_T) function pfftf_local_size_gc(rnk_n,local_n,local_start,gc_below,gc_above,local_ngc,local_gc_start) &
                                 bind(C, name='pfftf_local_size_gc')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfftf_local_size_gc
    
    integer(C_INTPTR_T) function pfftf_local_size_many_gc(rnk_n,local_n,local_start,howmany,gc_below,gc_above,local_ngc, &
                                                          local_gc_start) bind(C, name='pfftf_local_size_many_gc')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: local_start
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_ngc
      integer(C_INTPTR_T), dimension(*), intent(out) :: local_gc_start
    end function pfftf_local_size_many_gc
    
    type(C_PTR) function pfftf_plan_rgc_3d(Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfftf_plan_rgc_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_FLOAT), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_rgc_3d
    
    type(C_PTR) function pfftf_plan_cgc_3d(Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfftf_plan_cgc_3d_f03')
      import
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_cgc_3d
    
    type(C_PTR) function pfftf_plan_rgc(rnk_n,Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfftf_plan_rgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_FLOAT), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_rgc
    
    type(C_PTR) function pfftf_plan_cgc(rnk_n,Nos,gc_below,gc_above,data,comm_cart,gc_flags) bind(C, name='pfftf_plan_cgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_cgc
    
    type(C_PTR) function pfftf_plan_many_rgc(rnk_n,Nos,howmany,block,gc_below,gc_above,data,comm_cart,gc_flags) &
                         bind(C, name='pfftf_plan_many_rgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: block
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      real(C_FLOAT), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_many_rgc
    
    type(C_PTR) function pfftf_plan_many_cgc(rnk_n,Nos,howmany,block,gc_below,gc_above,data,comm_cart,gc_flags) &
                         bind(C, name='pfftf_plan_many_cgc_f03')
      import
      integer(C_INT), value :: rnk_n
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), dimension(*), intent(in) :: block
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_below
      integer(C_INTPTR_T), dimension(*), intent(in) :: gc_above
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: data
      integer(C_INT32_T), value :: comm_cart
      integer(C_INT), value :: gc_flags
    end function pfftf_plan_many_cgc
    
    subroutine pfftf_exchange(ths) bind(C, name='pfftf_exchange')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_exchange
    
    subroutine pfftf_reduce(ths) bind(C, name='pfftf_reduce')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_reduce
    
    subroutine pfftf_destroy_gcplan(ths) bind(C, name='pfftf_destroy_gcplan')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_destroy_gcplan
    
    subroutine pfftf_reset_gctimers(ths) bind(C, name='pfftf_reset_gctimers')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_reset_gctimers
    
    type(C_PTR) function pfftf_get_gctimer_exg(ths) bind(C, name='pfftf_get_gctimer_exg')
      import
      type(C_PTR), value :: ths
    end function pfftf_get_gctimer_exg
    
    type(C_PTR) function pfftf_get_gctimer_red(ths) bind(C, name='pfftf_get_gctimer_red')
      import
      type(C_PTR), value :: ths
    end function pfftf_get_gctimer_red
    
    subroutine pfftf_print_average_gctimer(ths,comm) bind(C, name='pfftf_print_average_gctimer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_print_average_gctimer
    
    subroutine pfftf_print_average_gctimer_adv(ths,comm) bind(C, name='pfftf_print_average_gctimer_adv_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_print_average_gctimer_adv
    
    subroutine pfftf_write_average_gctimer(ths,name,comm) bind(C, name='pfftf_write_average_gctimer_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_write_average_gctimer
    
    subroutine pfftf_write_average_gctimer_adv(ths,name,comm) bind(C, name='pfftf_write_average_gctimer_adv_f03')
      import
      type(C_PTR), value :: ths
      character(C_CHAR), dimension(*), intent(in) :: name
      integer(C_INT32_T), value :: comm
    end subroutine pfftf_write_average_gctimer_adv
    
    type(C_PTR) function pfftf_copy_gctimer(orig) bind(C, name='pfftf_copy_gctimer')
      import
      type(C_PTR), value :: orig
    end function pfftf_copy_gctimer
    
    subroutine pfftf_average_gctimer(ths) bind(C, name='pfftf_average_gctimer')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_average_gctimer
    
    type(C_PTR) function pfftf_add_gctimers(sum1,sum2) bind(C, name='pfftf_add_gctimers')
      import
      type(C_PTR), value :: sum1
      type(C_PTR), value :: sum2
    end function pfftf_add_gctimers
    
    type(C_PTR) function pfftf_reduce_max_gctimer(ths,comm) bind(C, name='pfftf_reduce_max_gctimer_f03')
      import
      type(C_PTR), value :: ths
      integer(C_INT32_T), value :: comm
    end function pfftf_reduce_max_gctimer
    
    subroutine pfftf_convert_gctimer2vec(ths,times) bind(C, name='pfftf_convert_gctimer2vec')
      import
      type(C_PTR), value :: ths
      real(C_DOUBLE), dimension(*), intent(out) :: times
    end subroutine pfftf_convert_gctimer2vec
    
    type(C_PTR) function pfftf_convert_vec2gctimer(times) bind(C, name='pfftf_convert_vec2gctimer')
      import
      real(C_DOUBLE), dimension(*), intent(in) :: times
    end function pfftf_convert_vec2gctimer
    
    subroutine pfftf_destroy_gctimer(ths) bind(C, name='pfftf_destroy_gctimer')
      import
      type(C_PTR), value :: ths
    end subroutine pfftf_destroy_gctimer
    
  end interface
