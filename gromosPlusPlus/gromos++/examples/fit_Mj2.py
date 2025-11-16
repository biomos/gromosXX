"""
This file is part of GROMOS.

Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
See <https://www.gromos.net> for details.

GROMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

"""
author dp

fit_Mj2.py script fits the Mj2.out data from epsilon program as described in [Schroeder & Steinhauser, J. Chem. Phys. 132, 244109 (2010); doi. 10.1063/1.3432620]
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

def plot_data(data2pl,mul=False,show=True):
	if mul:
		for i in range(len(data2pl)-1):
			plt.plot(data2pl[0],data2pl[i+1])
	else:
		plt.plot(data2pl[0],data2pl[1])
	if show:
		plt.show()
		plt.close()

def f_n(t,*p):
	data_new=p[0]*t
	for j in xrange(len(p)/2):
		data_new+=p[2*j+1]*(np.exp(-t/p[2*j+2])-1)
	return data_new

def f_l(t,k,a):
	data_new=k*t
	data_new+=a
	return data_new

def fit_Mj2(data,n=3):
	if n%2!=1:
		raise Exception("wrong number of parameters given")
	if 50.0/data[0][1] not in data[1]:
		params=[data[1][-1]/data[0][-1],1,10]
	else:
		params=[data[1][-1]/data[0][-1],500,0.1]
	params.extend([1]*(n-3))
	return curve_fit(f_n,data[0],data[1],params,maxfev = 100000000)

def read_data(f,p=0.625):
	data=[]
	for i in f:
		if not i.startswith("#"):
			data.append([])
			temp=i.split()
			for j in temp:
				data[-1].append(float(j))
	N_fr=int(len(data)*p)+1
	return np.array(zip(*data[:N_fr]))

if __name__=="__main__":
	par = argparse.ArgumentParser()
	par.add_argument("-d",dest="data_file", help="data files (usually Mj2.out)", type=file,required=True)
	par.add_argument("-e",dest="e_data_file", help="epsilon data files", type=file)
	par.add_argument("-n",dest="n", help="number of parameters",type=int, choices=range(3, 12, 2),default=3)
	par.add_argument("-f",dest="f", help="percent of data",type=float,default=0.625)
	par.add_argument("-p",dest="p", help="enable ploting",action='store_true',default=False)
	par.add_argument("-s",dest="s", help="save the plot",action='store_true',default=False)
	args = par.parse_args()
	if args.f<0 or args.p>1:
		print "-f must be [0,1]"
		par.print_help()
	data=read_data(args.data_file,args.f)
	fit=fit_Mj2(data,args.n)
	print "k\t",fit[0][0]
	A_sum=0
	for j in xrange(len(fit[0])/2):
		print "A"+str(j+1)+"\t",fit[0][2*j+1]
		A_sum-=fit[0][2*j+1]/2.0
		print "tau"+str(j+1)+"\t",fit[0][2*j+2]

	if args.e_data_file:
		l=args.e_data_file.readlines()
		temp=l[-4].split()
		eps_Md2=temp[-1]
		temp=l[-1].split()
		fac1=float(temp[1])#2epsilon_rf
		Md2=float(temp[2])
		fac2=float(temp[4])
		M_tot=Md2+A_sum
		M_tot-=float(temp[3])
		a=fac1*M_tot+fac2
		b=fac2-M_tot
		print "\n\n<Md2> = ",Md2
		print "<Md>2 = ",temp[3]
		print "im[Mjj] = ",A_sum,"\n"
		print "e_stat_Md2 = ",eps_Md2
		print "e_stat = ",a/b
		

	if args.p:
		plot_data([data[0],data[1],f_n(data[0],*fit[0])],True,False)
		plt.xlabel("time ps")
		plt.ylabel(r"Mj$^2$ nm$^2$e$^2$")
		if args.s:
			plt.savefig("Mj2.png")
		else:
			plt.show()
		plt.close()

