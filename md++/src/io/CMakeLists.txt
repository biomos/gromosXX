set(IO_CC
    io.cc
    argument.cc
    blockinput.cc
    instream.cc
    message.cc
    print_block.cc
    read_input.cc
    read_special.cc
    gzstream.cc
    parameter/in_parameter.cc
    parameter/check_parameter.cc 
    topology/in_topology.cc
    topology/in_perturbation.cc
    configuration/inframe.cc
    configuration/in_configuration.cc
    configuration/out_configuration.cc
    topology/in_posres.cc
    topology/in_distanceres.cc
    topology/in_angrest.cc
    topology/in_dihrest.cc
    topology/in_jvalue.cc
    topology/in_friction.cc
    topology/in_xray.cc
    topology/in_gamd.cc
    topology/in_order.cc
    topology/in_leus.cc
    topology/in_bsleus.cc
    topology/in_qmmm.cc
    topology/in_symrest.cc
    topology/in_rdc.cc
    topology/in_tf_rdc.cc
    topology/in_zaxisoribias.cc
)

add_library(groio OBJECT ${IO_CC})

target_include_directories(groio PUBLIC
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    .
    ..
    configuration
    parameter
    topology
)
