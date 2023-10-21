AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/replicaExchange

xxincludedir = $(includedir)/md++/replicaExchange

xxinclude_HEADERS = repex_mpi.h\
    replica_graph_control.h\
    replica_data.h\
	replica/_replica_Interface.h\
	replica/replica.h\
	replica/replica_MPI_master.h\
	replica/replica_MPI_slave.h\
	replica_exchangers/replica_exchange_base_interface.h\
	replica_exchangers/replica_exchange_master_interface.h\
	replica_exchangers/replica_exchange_slave_interface.h\
	replica_exchangers/2D_T_lambda_REPEX/re_base_2d_l_T_HREMD.h\
	replica_exchangers/2D_T_lambda_REPEX/re_master_2d_l_T_HREMD.h\
	replica_exchangers/2D_T_lambda_REPEX/re_slave_2d_l_T_HREMD.h\
	replica_exchangers/1D_S_RE_EDS/re_base_eds.h\
	replica_exchangers/1D_S_RE_EDS/re_master_eds.h\
	replica_exchangers/1D_S_RE_EDS/re_slave_eds.h\
	replica_exchangers/2D_S_Eoff_RE_EDS/re_base_2d_s_eoff_eds.h\
	replica_exchangers/2D_S_Eoff_RE_EDS/re_master_2d_s_eoff_eds.h\
	replica_exchangers/2D_S_Eoff_RE_EDS/re_slave_2d_s_eoff_eds.h

noinst_LTLIBRARIES = libgrorepex.la

libgrorepex_la_SOURCES = replicaExchange.cpp\
    replica_graph_control.cpp\
    replica/_replica_Interface.cpp\
    replica/replica.cc\
    replica/replica_MPI_master.cpp\
    replica/replica_MPI_slave.cpp\
    replica_exchangers/replica_exchange_base_interface.cc\
    replica_exchangers/replica_exchange_master_interface.cc\
    replica_exchangers/replica_exchange_slave_interface.cc\
    replica_exchangers/2D_T_lambda_REPEX/re_base_2d_l_T_HREMD.cc\
    replica_exchangers/2D_T_lambda_REPEX/re_master_2d_l_T_HREMD.cc\
    replica_exchangers/2D_T_lambda_REPEX/re_slave_2d_l_T_HREMD.cc\
    replica_exchangers/1D_S_RE_EDS/re_base_eds.cc\
    replica_exchangers/1D_S_RE_EDS/re_master_eds.cc\
    replica_exchangers/1D_S_RE_EDS/re_slave_eds.cc\
    replica_exchangers/2D_S_Eoff_RE_EDS/re_base_2d_s_eoff_eds.cc\
    replica_exchangers/2D_S_Eoff_RE_EDS/re_master_2d_s_eoff_eds.cc\
    replica_exchangers/2D_S_Eoff_RE_EDS/re_slave_2d_s_eoff_eds.cc


libgrorepex_la_LIBADD =

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
