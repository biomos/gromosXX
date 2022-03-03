set(MATH_CC
    math.cc
    volume.cc
    random.cc
    transformation.cc
    boundary_checks.cc
)

add_library(gromath OBJECT ${MATH_CC})

target_include_directories(gromath PUBLIC
    ..
    ${PROJECT_BINARY_DIR}
)