set(CONFIG_CC
    configuration_global.cc
    configuration.cc
    energy.cc
    average.cc
    state_properties.cc
    mesh.cc
    kspace.cc
    influence_function.cc
)

add_library(groconfiguration OBJECT ${CONFIG_CC})

target_include_directories(groconfiguration PUBLIC
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ..
)