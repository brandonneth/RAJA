//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-20, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "RAJA/RAJA.hpp"
#include "camp/resource.hpp"


/*
 * RAJA Teams Example: Upper Triangular Pattern + Shared Memory
 *
 * Teams introduces hierarchal parallism through the concept of
 * teams and threads.  Computation is executed in a pre-defined grid
 * composed of threads and grouped into teams. The teams model enables
 * developers to express parallism through loops over teams, and sub loops
 * over threads. Team loops are executed in parallel and
 * threads within a team should be treated as sub-parallel regions.
 *
 * Team shared memory is allocated between team and thread loops.
 * Memory allocated within thread loops are thread private.
 * The example below demonstrates composing an upper triangular
 * loop pattern, and shared memory.
 *
 */

/*
 * Define potential launch policies
 * Currenly we only support two:
 * Host and Device.
 */
using launch_policy = RAJA::expt::LaunchPolicy<
#if defined(RAJA_ENABLE_OPENMP)
    RAJA::expt::omp_launch_t
#else
    RAJA::expt::seq_launch_t
#endif
#if defined(RAJA_ENABLE_CUDA)
    ,
    RAJA::expt::cuda_launch_t<false>
#endif
#if defined(RAJA_ENABLE_HIP)
    ,
    RAJA::expt::hip_launch_t<false>
#endif
    >;

/*
 * Define team policies.
 * Up to 3 dimension are supported: x,y,z
 */
using teams_x = RAJA::expt::LoopPolicy<RAJA::loop_exec
#if defined(RAJA_ENABLE_CUDA)
                                       ,
                                       RAJA::cuda_block_x_direct
#endif
#if defined(RAJA_ENABLE_HIP)
                                       ,
                                       RAJA::hip_block_x_direct
#endif
                                       >;
/*
 * Define thread policies.
 * Up to 3 dimension are supported: x,y,z
 */
using threads_x = RAJA::expt::LoopPolicy<RAJA::loop_exec
#if defined(RAJA_ENABLE_CUDA)
                                         ,
                                         RAJA::cuda_thread_x_loop
#endif
#if defined(RAJA_ENABLE_HIP)
                                         ,
                                         RAJA::hip_thread_x_loop
#endif
                                         >;


int main(int RAJA_UNUSED_ARG(argc), char **RAJA_UNUSED_ARG(argv[]))
{

  // Resource object for host
  camp::resources::Host host_res;

  // Resource objects for CUDA or HIP
#if defined(RAJA_ENABLE_CUDA)
  camp::resources::Cuda device_res;
#endif

#if defined(RAJA_ENABLE_HIP)
  camp::resources::Hip device_res;
#endif

  std::cout << "\n Running RAJA-Teams examples...\n";

  for (int exec_place = 0; exec_place < (int)RAJA::expt::NUM_PLACES;
       ++exec_place) {
    RAJA::expt::ExecPlace select_cpu_or_gpu = (RAJA::expt::ExecPlace)exec_place;

    // auto select_cpu_or_gpu = RAJA::HOST;
    // auto select_cpu_or_gpu = RAJA::DEVICE;

    int N_tri = 5;
    int *Ddat;
    if (select_cpu_or_gpu == RAJA::expt::HOST)
      Ddat = host_res.allocate<int>(N_tri * N_tri);
    if (select_cpu_or_gpu == RAJA::expt::DEVICE)
      Ddat = device_res.allocate<int>(N_tri * N_tri);

    /*
     * launch just starts a "kernel" it's doesn't provide any looping.
     *
     * The first argument determines which policy should be executed,
     *
     * The second argument is the number of teams+threads needed for each of the
     * policies.
     *
     * Third argument is the lambda for the policy.
     *
     *
     * The lambda takes a "resource" object, which has the teams+threads
     */


    std::cout << "\n Running Upper triangular pattern example...\n";


    RAJA::View<int, RAJA::Layout<2>> D(Ddat, N_tri, N_tri);

    RAJA::expt::launch<launch_policy>(
        select_cpu_or_gpu,
        RAJA::expt::Resources(RAJA::expt::Teams(N_tri),
                              RAJA::expt::Threads(N_tri)),
        [=] RAJA_HOST_DEVICE(RAJA::expt::LaunchContext ctx) {
          RAJA::expt::loop<teams_x>(ctx, RAJA::RangeSegment(0, N_tri), [&](int r) {

                // Array shared within threads of the same team
                TEAM_SHARED int s_A[1];

                RAJA::expt::loop<threads_x>(ctx, RAJA::RangeSegment(r, N_tri), [&](int c) {
                    if (c == r) s_A[0] = r;
                    D(r, c) = r * N_tri + c;
                });  // loop j

                RAJA::expt::loop<threads_x>(ctx, RAJA::RangeSegment(r, N_tri), [&](int c) {

                    printf("r=%d, c=%d : D=%d : s_A = %d \n", r, c, D(r, c), s_A[0]);

                    });  // loop c
              });        // loop r
        });              // outer lambda


    if (select_cpu_or_gpu == RAJA::expt::HOST) {
      host_res.deallocate(Ddat);
    }
    if (select_cpu_or_gpu == RAJA::expt::DEVICE) {
      device_res.deallocate(Ddat);
    }

  }  // Execution places loop


}  // Main
