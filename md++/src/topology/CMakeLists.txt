set(TOP_CC
    topology.cc
    topology_global.cc
    exclusions.cc
)

add_library(grotopology OBJECT ${TOP_CC})

target_include_directories(grotopology PUBLIC
    ..
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
)
target_link_libraries(grotopology PRIVATE ${EXTERNAL_LIBRARIES})