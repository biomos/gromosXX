#!/troll2/markus/programs/bin/python

import os

_name = "default"
_nr = 1

program = "/troll2/markus/programs/gromosXX/LINUX_GCC_RELEASE/program/md"
simdir = os.getcwd()
workdir = "/scrloc/${USER}"
python = "python"
pert = 0

def set_name(name, nr=1):
    global _name, _nr, _script_name
    _name = name
    _nr = nr
    return None

def write(rep):
    global _name, _nr, program, simdir, workdir

    repdir = simdir + "/" + str(rep) + "/"
    full_name = _name + "_" + str(_nr)
    old_name = _name + "_" + str(_nr-1)
    script_name = repdir + "/j" + full_name + ".sh"

    f = open(script_name, "w")
    f.write("#!/bin/sh\n")
    f.write("# job script for " + _name + " (" + str(_nr) + ") \n\n")
    
    f.write("PROG=" + program + "\n\n")
    
    f.write("WORKDIR=" + workdir + "\n")
    f.write("SIMDIR=" + repdir + "\n")
    f.write("mkdir -p ${WORKDIR}\n")
    f.write("cd ${WORKDIR}\n\n")

    f.write("TOPO=${SIMDIR}/../" + _name + ".topo\n")
    f.write("INPUT=${SIMDIR}/" + full_name + ".in\n")
    if crd == "":
        f.write("CONF=${SIMDIR}/" + old_name + ".fin\n\n")
    else:
        f.write("CONF="+crd+"\n\n")

    f.write("OUTPUT=" + full_name + ".out\n")
    f.write("FIN=" + full_name + ".fin\n")
    f.write("TRJ=" + full_name + ".trj\n")
    f.write("TRV=" + full_name + ".trv\n")
    f.write("TRE=" + full_name + ".tre\n")
    f.write("BAE=" + full_name + ".bae\n")
    if pert != 0:
        f.write("TRG=" + full_name + ".trg\n")
        f.write("BAG=" + full_name + ".bag\n")

    f.write("REP=" + full_name + ".rep\n")

    f.write("\n\n")

    f.write("export LD_LIBRARY_PATH=/troll2/markus/programs/lib\n")

    f.write("\n\n")

    f.write("${PROG} \\\n")
    f.write("\t@topo\t${TOPO} \\\n")
    f.write("\t@conf\t${CONF} \\\n")
    f.write("\t@input\t${INPUT} \\\n")
    f.write("\t@fin\t${FIN} \\\n")
    f.write("\t@trj\t${TRJ} \\\n")
    f.write("\t@trv\t${TRV} \\\n")
    f.write("\t@tre\t${TRE} \\\n")
    f.write("\t@bae\t${BAE} \\\n")
    if pert != 0:
        f.write("\t@trg\t${TRG} \\\n")
        f.write("\t@bag\t${BAG} \\\n")
    f.write("\t@rep\t${REP} \\\n")
    f.write("\t>\t${OUTPUT}\n")

    f.write("\n\n")
    f.write("uname -a >> ${OUTPUT}\n\n")

    f.write("gzip ${TRJ}\n")
    f.write("gzip ${TRV}\n")
    f.write("gzip ${TRE}\n")
    f.write("gzip ${BAE}\n")
    if pert != 0:
        f.write("gzip ${TRG}\n")
        f.write("gzip ${BAG}\n")

    f.write("\n\n")

    f.write("cp ${OUTPUT}\t\t${SIMDIR}\n")
    f.write("cp ${FIN}\t\t${SIMDIR}\n")
    f.write("cp ${TRJ}.gz\t\t${SIMDIR}\n")
    f.write("cp ${TRV}.gz\t\t${SIMDIR}\n")
    f.write("cp ${TRE}.gz\t\t${SIMDIR}\n")
    f.write("cp ${BAE}.gz\t\t${SIMDIR}\n")
    if pert != 0:
        f.write("cp ${TRG}.gz\t\t${SIMDIR}\n")
        f.write("cp ${BAG}.gz\t\t${SIMDIR}\n")

    f.write("cp ${REP}\t\t${SIMDIR}\n")
    
    f.write("\n\n")

    f.write("cd ${SIMDIR}\n")
    f.write("rm -r ${WORKDIR}\n")

    f.write("\n\n")
    
    f.close()

    os.chmod(script_name, 0755)
    
    return None


def write_qmaster(rep1, rep2, T):
    global _name, _nr, program, simdir, workdir

    repdir1 = simdir + "/" + str(T[rep1])
    repdir2 = simdir + "/" + str(T[rep2])

    full_name = _name + "_" + str(_nr)
    script_name = simdir + "/mj" + full_name + "_" + str(T[rep1]) + "_" + str(T[rep2]) +".sh"

    job1 = repdir1 + "/j" + full_name + ".sh"
    job2 = repdir2 + "/j" + full_name + ".sh"
    
    f = open(script_name, "w")
    f.write("#!/bin/sh\n")
    f.write("# master job script for " + _name + " (" + str(_nr) + ") \n\n")

    f.write("# first replica\n")
    f.write(job1 + " &\n\n")

    f.write("# second replica\n")
    f.write(job2 + " &\n\n")

    f.write("# wait for jobs to finish\n")
    f.write("wait\n\n")

    f.write("# check switch\n")
    f.write(python + " repex.py " + str(rep1) + " " + str(rep1 + 1) + " " + str(_nr + 1) + " >> repex.out\n")
    f.write(python + " repex.py " + str(rep2 + 1) + " " + str(rep2 + 2) + " " + str(_nr + 1) + " >> repex.out\n")    

    f.write("\n\n")

    os.chmod(script_name, 0755)
    return None


def write_smaster(rep1, rep2):
    global _name, _nr, program, simdir, workdir

    repdir1 = simdir + "/rep" + str(rep1)
    repdir2 = simdir + "/rep" + str(rep2)

    full_name = _name + "_" + str(_nr)
    script_name = simdir + "/mj" + full_name + "_" + str(rep1) + "_" + str(rep2) +".sh"

    job1 = "j" + full_name + ".sh"
    job2 = "j" + full_name + ".sh"
    
    f = open(script_name, "w")
    f.write("#!/bin/sh\n")
    f.write("# master job script for " + _name + " (" + str(_nr) + ") \n\n")
    f.write("echo " + _name + " run " + str(_nr) + " rep1 = " + str(rep1) + " rep2 = " + str(rep2) + "\n\n")

    f.write("# first replica\n")
    f.write("cd " + repdir1 + "\n")
    f.write(job1 + "\n\n")

    f.write("# second replica\n")
    f.write("cd " + repdir2 + "\n")    
    f.write(job2 + "\n\n")

    f.write("cd " + simdir + "\n\n")

    f.write("# check switch\n")
    f.write(python + " repex.py " + str(rep1 - 1) + " " + str(rep1) + " " + str(_nr + 1) + "\n")
    f.write(python + " repex.py " + str(rep2) + " " + str(rep2 + 1) + " " + str(_nr + 1) + "\n")    

    f.write("\n\n")

    os.chmod(script_name, 0755)
    return None



