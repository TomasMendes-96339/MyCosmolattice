      integer, parameter :: PFFT_REDFT00 = FFTW_REDFT00
      integer, parameter :: PFFT_REDFT01 = FFTW_REDFT01
      integer, parameter :: PFFT_REDFT10 = FFTW_REDFT10
      integer, parameter :: PFFT_REDFT11 = FFTW_REDFT11
      integer, parameter :: PFFT_RODFT00 = FFTW_RODFT00
      integer, parameter :: PFFT_RODFT01 = FFTW_RODFT01
      integer, parameter :: PFFT_RODFT10 = FFTW_RODFT10
      integer, parameter :: PFFT_RODFT11 = FFTW_RODFT11

      integer, parameter :: PFFT_TRANSPOSED_NONE   = 0
      integer, parameter :: PFFT_TRANSPOSED_IN     = 1
      integer, parameter :: PFFT_TRANSPOSED_OUT    = 2
      integer, parameter :: PFFT_SHIFTED_IN        = 4
      integer, parameter :: PFFT_SHIFTED_OUT       = 8
      integer, parameter :: PFFT_MEASURE           = 0
      integer, parameter :: PFFT_ESTIMATE          = 16
      integer, parameter :: PFFT_PATIENT           = 32
      integer, parameter :: PFFT_EXHAUSTIVE        = 64
      integer, parameter :: PFFT_NO_TUNE           = 0
      integer, parameter :: PFFT_TUNE              = 128
      integer, parameter :: PFFT_PRESERVE_INPUT    = 256
      integer, parameter :: PFFT_DESTROY_INPUT     = 512
      integer, parameter :: PFFT_BUFFERED_INPLACE  = 1024
      
      integer, parameter :: PFFT_GC_NONTRANSPOSED  = 0
      integer, parameter :: PFFT_GC_TRANSPOSED     = 1
      integer, parameter :: PFFT_GC_SENDRECV       = 2
      integer, parameter :: PFFT_GC_RMA            = 4
      
      integer, parameter :: PFFT_INT               = 1
      integer, parameter :: PFFT_PTRDIFF_T         = 2
      integer, parameter :: PFFT_FLOAT             = 3
      integer, parameter :: PFFT_DOUBLE            = 4
      integer, parameter :: PFFT_UNSIGNED          = 5

      integer, parameter :: PFFT_FORWARD           = -1
      integer, parameter :: PFFT_BACKWARD          = +1

      integer, parameter ::ptrdiff_t_kind = 8

      integer(ptrdiff_t_kind), parameter :: PFFT_DEFAULT_BLOCK = 0
      integer(ptrdiff_t_kind), parameter :: PFFT_DEFAULT_BLOCKS = -1
      integer(ptrdiff_t_kind), parameter :: PFFT_NO_GCELLS = -1







