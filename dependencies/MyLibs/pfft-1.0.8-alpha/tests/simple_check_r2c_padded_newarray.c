#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err;
  double *planned_in, *executed_in;
  pfft_complex *planned_out, *executed_out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE| PFFT_PADDED_R2C,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory for planning */
  planned_in  = pfft_alloc_real (2 * alloc_local);
  planned_out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft_r2c_3d(
      n, planned_in, planned_out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_PADDED_R2C);
  
  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft_c2r_3d(
      n, planned_out, planned_in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_PADDED_C2R);

  /* Free planning arrays since we use other arrays for execution */
  pfft_free(planned_in); pfft_free(planned_out);

  /* Allocate memory for execution */
  executed_in  = pfft_alloc_real(2 * alloc_local);
  executed_out = pfft_alloc_complex(alloc_local);

  /* Initialize input with random numbers */
  pfft_init_input_real(3, n, local_ni, local_i_start,
      executed_in);
  
  /* execute parallel forward FFT */
  pfft_execute_dft_r2c(plan_forw, executed_in, executed_out);

  /* clear the old input */
  pfft_clear_input_real(3, n, local_ni, local_i_start,
      executed_in);

  /* execute parallel backward FFT */
  pfft_execute_dft_c2r(plan_back, executed_out, executed_in);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    executed_in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = pfft_check_output_real(3, n, local_ni, local_i_start, executed_in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(executed_in); pfft_free(executed_out);
  MPI_Finalize();
  return 0;
}
