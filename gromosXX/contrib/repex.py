#!/troll2/markus/programs/bin/python

################################################################################
##
## REPLICA EXCHANGE SIMULATIONS
##
################################################################################

import sys
import locale
import random
import shutil
import os
import script
import math

################################################################################
# run parameter
#
# number of MD simulations followed by a MC temperature exchange trial move
max_run = 1000
# system name
name = "my_system"
# path (simulation directory)
path = os.getcwd()
# Temperatures to run replicas at
T = [300, 305, 310, 315, 321, 326, 332, 337, 343, 349, 354, 360]
#
# md binary
program = "/troll2/markus/programs/gromosXX/LINUX_GCC_RELEASE/bin/md"
# temporary (run) directory (temperature will be added)
scr = "/scrloc/${USER}_"
# submit command (script name will be added)
submit = "psub -s penguin 2 "
# single run (don't use the queue, run only on 1 processor)
single = False
#
################################################################################

################################################################################
# constants
# Boltzmann constant
kB = 0.00831441
#
################################################################################

################################################################################
# private
_conf = ["", ""]
################################################################################

################################################################################
# get_name
# produce a file name from temperature index, run number and postfix
################################################################################
def get_name(t, run, post):
    global path, name, T
    return path + "/" + str(T[t]) + "/" + name + "_" + str(run) + post
#
# end get_name
################################################################################

################################################################################
# write_ts
# write out time series of exchanges, exchange probablities
################################################################################
def write_ts(ID1, T01, ID2, T02):
    # write time series of temperatures per replica ...
    nr1 = "replica_" + str(ID1) + ".dat"
    nr2 = "replica_" + str(ID2) + ".dat"

    f1 = open(nr1, "a")
    f2 = open(nr2, "a")

    f1.write("\t" + str(run) + "\t\t" + str(T01) + "\n")
    f2.write("\t" + str(run) + "\t\t" + str(T02) + "\n")

    f1.close()
    f2.close()

    # ... and replicas per temperature
    nt1 = "temp_" + str(T01) + ".dat"
    nt2 = "temp_" + str(T02) + ".dat"

    f1 = open(nt1, "a")
    f2 = open(nt2, "a")

    f1.write("\t" + str(run) + "\t\t" + str(ID1) + "\n")
    f2.write("\t" + str(run) + "\t\t" + str(ID2) + "\n")

    f1.close()
    f2.close()
#
# end write_ts
################################################################################

################################################################################
#
# Replica Exchange
# calculate exchange probability, exchange starting coordinates if necessary
# submit next jobs
#
################################################################################
def repex(r1, r2, run):
    global max_run, path, name, program, kB, _conf

    b = 1
    if r1 == 0:
        r1 = len(T)
        b = 0
    if r2 > len(T):
        r2 = 1
        b = 0

    # 0 base indices (but command line arguments taken from 1..N)
    r1 -= 1
    r2 -= 1

    if run > max_run + 1:
        # this is the end...
        print "maximum number of runs reached for ", T[r1], " and ", T[r2], "...\n\n"
        sys.exit()
    
    name1 = get_name(r1, run-1, ".rep")
    try:
        f1 = open(name1, "r")
    except:
        print "\twait: ", T[r1], " - ", T[r2], " run ", run, " -- ", T[r1], " not finished"
        return False

    name2 = get_name(r2, run-1, ".rep")
    try:        
        f2 = open(name2, "r")
    except:
        print "\twait: ", T[r1], " - ", T[r2], " run ", run, " -- ", T[r2], " not finished"
        return False

    # this is ugly: parse the files...
    # ID Temperature(ref) Energy Temperature(current)
    data = 0
    for line in f1:
        if line == "END\n":
            data = 0
        if data == 1:
            spec = line
        
        if line == "REPEXDATA\n":
            data = 1

    l = spec.split()
    ID1 = locale.atoi(l[0])
    T01 = locale.atof(l[1])
    E1 = locale.atof(l[2])
    T1 = locale.atof(l[3])
    s1 = 0

    f1.close()

    data = 0
    for line in f2:
        if line == "END\n":
            data = 0
        if data == 1:
            spec = line
        
        if line == "REPEXDATA\n":
            data = 1

    l = spec.split()
    ID2 = locale.atoi(l[0])
    T02 = locale.atof(l[1])
    E2 = locale.atof(l[2])
    T2 = locale.atof(l[3])
    s2 = 0
    
    f2.close()
    
    # calculate exchange probability
    p = math.exp(-(1.0 / (kB * T1) - 1.0 / (kB * T2)) * (E2 - E1))
    
    # draw a random number
    r = random.random()

    print "preparing run ", run
    print "\tReplica ", ID1, " : ", T[r1], " (", T1,"K): ", E1, "kJ/mol"
    print "\tReplica ", ID2, " : ", T[r2], " (", T2,"K): ", E2, "kJ/mol"
    print "\t=> P(switch) : ", p, "\t\trandom : ", r

    # set the startup files correctly
    # these will be used when writing the script
    _conf[0] = get_name(r1, run - 1, ".fin")
    _conf[1] = get_name(r2, run - 1, ".fin")
    
    # don't switch the 1 - max_run pair
    if p > r and b == 1:
        print "\t!!! switch " + str(T[r1]) + " - " + str(T[r2]) + " for run " + str(run)
        dummy = _conf[0]
        _conf[0] = _conf[1]
        _conf[1] = dummy
        dummy = ID1
        ID1 = ID2
        ID2 = dummy
        s1 = 1
        s2 = 1
    else:
        print "\t!!! keep " + str(T[r1]) + " - " + str(T[r2]) + " for run " + str(run)

    # prepare the input files
    path1 = path + "/" + str(T[r1]) + "/"
    path2 = path + "/" + str(T[r2]) + "/"

    t1 = path1 + name + "_" + str(run) + ".in"
    t2 = path2 + name + "_" + str(run) + ".in"

    # check if file already exist, then somebody else was faster!!!
    if os.path.exists(t1):
        print "--- synchronisation problem: ", r1+1, " - ", r2+1, " run ", run
        sys.exit()

    if os.path.exists(t2):
        print "--- synchronisation problem: ", r1+1, " - ", r2+1, " run ", run
        sys.exit()

    shutil.copy(path1 + name + ".in", t1)
    shutil.copy(path2 + name + ".in", t2)

    # append the REPLICA blocks with reference temperatures
    f1 = open(t1, "a")
    f2 = open(t2, "a")

    f1.write("REPLICA\n")
    f1.write("#\tID\tT0\tscale\n")
    f1.write("\t" + str(ID1) + "\t" + str(T01) + "\t" + str(s1)+ "\n")
    f1.write("END\n")
    
    f2.write("REPLICA\n")
    f2.write("#\tID\tT0\tscale\n")
    f2.write("\t" + str(ID2) + "\t" + str(T02) + "\t" + str(s2)+ "\n")
    f2.write("END\n")

    f1.close()
    f2.close()

    # write out the time - series
    write_ts(ID1, T01, ID2, T02, E1, E2, p, s1)
    
    return True
#
# end repex
################################################################################


################################################################################
#
# run over the queue
#
# uses _conf[0] and _conf[1], set by repex!
# writes the scripts and the master script
################################################################################
#
def qrun(r1, r2, run, start=1):
    global path, name, T, _conf, scr, single

    if r1 == 0:
        r1 = len(T)
    if r2 > len(T):
        r2 = 1

    r1 -= 1
    r2 -= 1

    script.workdir = scr + name + "_" + str(T[r1])
    script.crd = _conf[0]
    
    script.write(T[r1])

    script.workdir = scr + name + "_" + str(T[r2])
    script.crd = _conf[1]
    
    script.write(T[r2])

    script.workdir = scr + name
    script.crd = ""

    if not single:
        script.write_qmaster(r1, r2, T)

        if start==1:
            script_name = "mj" + name + "_" + str(run) + "_" + str(T[r1]) + "_" + str(T[r2]) + ".sh"
            e = os.system(submit + script_name)
            if e != 0:
                print "submitting master job(" + str(T[r1]) + ", " + str(T[r2]) + ") failed!"
            else:
                print "job ", T[r1], " - ", T[r2], " run ", run, " submitted\n"
        else:
            print "queue master script written\n"

    return None
#
# end qrun
################################################################################


################################################################################
#
# single run
#
# keeps a stack of jobs it still has to run
################################################################################
#
def srun(run):
    global path, T, program, name, _conf, scr
    class item:
        def __init__(self, T, run):
            self.T = T
            self.run = run
    stack = []

    for t in T:
        stack.append(item(t, run))

    while len(stack) > 0:
        i = stack.pop()
        print "starting ", i.T, " run ", i.run

        scriptname = path + "/" + str(i.T) + "/j" + name + "_" + str(i.run) + ".sh"
        if not os.path.exists(scriptname):
            print "script not found!"
            sys.exit()

        e = os.system(scriptname)
        if e != 0:
            print "running script failed!"
            sys.exit()

        # try exchange
        if (i.run % 2) + ((T.index(i.T)+1) % 2) == 1:
            r2 = T.index(i.T) + 2
        else:
            r2 = T.index(i.T) + 1

        r1 = r2 - 1

        if repex(r1, r2, i.run + 1):
            # continue simulations
            if r1 < 1:
                r1 = len(T)
            if r2 > len(T):
                r2 = 1

            r1 -= 1
            r2 -= 1

            script.set_name(name, i.run + 1)

            script.workdir = scr + name + "_" + str(T[r1])
            script.crd = _conf[0]
    
            script.write(T[r1])

            script.workdir = scr + name + "_" + str(T[r2])
            script.crd = _conf[1]

            script.write(T[r2])

            script.workdir = scr + name
            script.crd = ""

            # add the jobs to the stack...
            stack.append(item(T[r1], i.run + 1))
            stack.append(item(T[r2], i.run + 1))    

    print "single run done\n\n"

    return None
            

################################################################################
#
# setup the system
# if run == 1, it creates directories, sym-links and master scripts
# otherwise, it just creates the input files
#
################################################################################
def setup(run):
    global path, name, program, T

    simdir = path + "/"
    i = 0
    for t in T:
        i += 1
        print "\npreparing replica ", i, " at T ", t, " (run ", run, ")\n"
    
        # prepare the input files
        path1 = path + "/" + str(t) + "/"
        if run == 1:
            os.mkdir(path1)
            shutil.copy(simdir + name + ".in", path1 + name + ".in")
            shutil.copy(simdir + name + ".crd", path1 + name + "_" + str(run-1) + ".fin")
        
        shutil.copy(path1 + name + ".in", path1 + name + "_" + str(run) + ".in")

        f1 = open(path1 + name + "_" + str(run) + ".in", "a")

        f1.write("REPLICA\n")
        f1.write("#\tID\tT0\tscale\n")
        f1.write("\t" + str(i) + "\t" + str(t) + "\t" + str(1)+ "\n")
        f1.write("END\n")

        f1.close()

        # write time series of temperatures per replica ...
        nr1 = "replica_" + str(i) + ".dat"

        f1 = open(nr1, "a")
        f1.write("\t" + str(run) + "\t\t" + str(t) + "\n")
        f1.close()

        # ... and replicas per temperature
        nt1 = "temp_" + str(t) + ".dat"
        f1 = open(nt1, "a")
        f1.write("\t" + str(run) + "\t\t" + str(i) + "\n")
        f1.close()

    # write scripts and master scripts
    if run == 1:
        for i in range(1, len(T) + 1, 2):
            qrun(i, i+1, run, 0)
        
    return None

################################################################################
#
# usage
#
################################################################################
def print_usage():
    print "usage:\n"
    print "\trepex rep1 rep2 run"
    print "\trepex setup run"
    print "\trepex srun rep1 rep2 run"
    print "\trepex qrun rep1 rep2 run\n"
    print "0.\tprepare replica directories (repN) and put"
    print "\tcoordinates (name_0.fin) and"
    print "\tinput files (name.in) there"
    print "1.\trun setup on all replicas"
    print "2.\trun qrun (or srun) on each pair of simulations"
    print "3.\tstart the first n/2 master scripts manually or"
    print "\twrite a script that does it for you"
    print "4.\thope...\n\n"
    return None


################################################################################
#
# main
#
################################################################################
if len(sys.argv) < 3:
    print_usage()
    sys.exit()
    
if sys.argv[1] == "setup":
    run = locale.atoi(sys.argv[2])
    
    script.set_name(name, run)
    script.program = program
    script.simdir = path

    setup(run)
    sys.exit()
        
elif sys.argv[1] == "qrun":
    if len(sys.argv) < 5:
        print_usage()
        sys.exit()
        
    r1 = locale.atoi(sys.argv[2])
    r2 = locale.atoi(sys.argv[3])
    run = locale.atoi(sys.argv[4])
    # only write but don't run scripts
    script.set_name(name, run)
    script.program = program
    script.simdir = path
	
    qrun(r1, r2, run, 0)
    sys.exit()    

elif sys.argv[1] == "srun":
    run = locale.atoi(sys.argv[2])
    
    script.set_name(name, run)
    script.program = program
    script.simdir = path

    srun(run)
    sys.exit()
    
else:
    r1 = locale.atoi(sys.argv[1])
    r2 = locale.atoi(sys.argv[2])
    run = locale.atoi(sys.argv[3])

    script.set_name(name, run)
    script.program = program
    script.simdir = path
    if repex(r1, r2, run):
        qrun(r1, r2, run)
    else:
        sys.exit()
