set(REPEX_CC
    replicaExchange.cpp
    replica_graph_control.cpp
    replica/_replica_Interface.cpp
    replica/replica.cc
    replica/replica_MPI_master.cpp
    replica/replica_MPI_slave.cpp
    replica_exchangers/replica_exchange_base_interface.cc
    replica_exchangers/replica_exchange_master_interface.cc
    replica_exchangers/replica_exchange_slave_interface.cc
    replica_exchangers/2D_T_lambda_REPEX/re_base_2d_l_T_HREMD.cc
    replica_exchangers/2D_T_lambda_REPEX/re_master_2d_l_T_HREMD.cc
    replica_exchangers/2D_T_lambda_REPEX/re_slave_2d_l_T_HREMD.cc
    replica_exchangers/1D_S_RE_EDS/re_base_eds.cc
    replica_exchangers/1D_S_RE_EDS/re_master_eds.cc
    replica_exchangers/1D_S_RE_EDS/re_slave_eds.cc
    replica_exchangers/2D_S_Eoff_RE_EDS/re_base_2d_s_eoff_eds.cc
    replica_exchangers/2D_S_Eoff_RE_EDS/re_master_2d_s_eoff_eds.cc
    replica_exchangers/2D_S_Eoff_RE_EDS/re_slave_2d_s_eoff_eds.cc
)

add_library(grorepex OBJECT ${REPEX_CC})

target_include_directories(grorepex PUBLIC
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    .
    ..
)
