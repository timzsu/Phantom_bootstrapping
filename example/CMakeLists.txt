# test examples
add_executable(example_context
  1_bfv.cu
  2_bgv.cu
  3_ckks.cu
  4_kernel_fusing.cu
  5_bootstrapping.cu
  6_nn.cu
  example.cu
)

target_link_libraries(example_context PRIVATE Phantom PhantomBootstrapping)
target_include_directories(example_context PUBLIC ${CMAKE_SOURCE_DIR}/src)