set(SRCFILES
  src/CacheManager.cpp
  )
set(HEADERS
  include/CacheManager.h
  include/CacheParameters.h
  include/CacheWeights.h
  include/WeightNormalization.h
  include/WeightMonotonicSpline.h
  include/WeightUniformSpline.h
  include/WeightGeneralSpline.h
  include/WeightBase.h
  include/CacheIndexedSums.h
  )

if(CMAKE_CUDA_COMPILER)
  cmessage(STATUS "CUDA utilities being compiled")
  set(SRCFILES ${SRCFILES} src/WeightNormalization.cu)
  set(SRCFILES ${SRCFILES} src/WeightMonotonicSpline.cu)
  set(SRCFILES ${SRCFILES} src/WeightUniformSpline.cu)
  set(SRCFILES ${SRCFILES} src/WeightGeneralSpline.cu)
  set(SRCFILES ${SRCFILES} src/CacheParameters.cu)
  set(SRCFILES ${SRCFILES} src/CacheWeights.cu)
  set(SRCFILES ${SRCFILES} src/CacheIndexedSums.cu)
else(CMAKE_CUDA_COMPILER)
  cmessage(STATUS "CUDA utilities are not being compiled")
  set(SRCFILES ${SRCFILES} src/WeightNormalization.cpp)
  set(SRCFILES ${SRCFILES} src/WeightMonotonicSpline.cpp)
  set(SRCFILES ${SRCFILES} src/WeightUniformSpline.cpp)
  set(SRCFILES ${SRCFILES} src/WeightGeneralSpline.cpp)
  set(SRCFILES ${SRCFILES} src/CacheParameters.cpp)
  set(SRCFILES ${SRCFILES} src/CacheWeights.cpp)
  set(SRCFILES ${SRCFILES} src/CacheIndexedSums.cpp)
endif(CMAKE_CUDA_COMPILER)

if( USE_STATIC_LINKS )
  add_library(GundamCache STATIC ${SRCFILES})
else()
  add_library(GundamCache SHARED ${SRCFILES})
endif()

# Make sure the current directories are available for the later
# compilation.
target_include_directories( GundamCache PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries( GundamCache PUBLIC GundamFitParameters GundamFitSamples ${ROOT_LIBRARIES} )

#  set_target_properties(GundamCache PROPERTIES VERSION "${GUNDAM_VERSION_STRING}")

install(TARGETS GundamCache DESTINATION lib)
