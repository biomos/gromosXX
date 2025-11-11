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

fit_JJ.py script fits the JJ.out data from epsilon program as described in [Schroeder & Steinhauser, J. Chem. Phys. 132, 244109 (2010); doi. 10.1063/1.3432620]
"""

import numpy as np
import math
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
	if len(p)%2!=0:
		raise Exception("wrong number of parameters given")
	data_new=0
	for j in xrange(len(p)/4):
		data_new+=p[4*j]*np.cos(p[4*j+1]*t+p[4*j+2])*np.exp(-t/p[4*j+3])
	if len(p)%4!=0:
		data_new+=p[-2]*np.exp(-t/p[-1])
	return data_new

def fit_Mjj(data,n=10):
	if n%2!=0:
		raise Exception("wrong number of parameters given - has to be an even number")
	params=[1]*n
	return curve_fit(f_n,data[0],data[1],params,maxfev = 100000000)

def im_Mjj(fit_params,pr=False):
	if len(fit_params)%4!=0:
		params=fit_params[:-1]
		params=np.append(params,[0]*2)
		params=np.append(params,fit_params[-1])
	else:	
		params=list(fit_params)
	im_Mjj=0
	for i in xrange(len(params)/4):
		tau2=params[4*i+3]**2
		omega2=params[4*i+1]**2
		to=params[4*i+3]*params[4*i+1]
		to2=tau2*omega2
		fac1=params[4*i]*tau2/((to2+1)**2)
		fac2=math.cos(params[4*i+2])*(to2-1)+2*to*math.sin(params[4*i+2])
		im_Mjj+=fac1*fac2
	return im_Mjj

def read_data(f,p=1):
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
	par.add_argument("-d",dest="data_file", help="data files (usually JJ.out)", type=file,required=True)
	par.add_argument("-e",dest="e_data_file", help="epsilon data files", type=file)
	par.add_argument("-n",dest="n", help="number of parameters",type=int, choices=range(4, 19,2),default=10)
	par.add_argument("-f",dest="f", help="percent of data",type=float,default=1)
	par.add_argument("-p",dest="p", help="enable ploting",action='store_true',default=False)
	par.add_argument("-s",dest="s", help="save the plot",action='store_true',default=False)
	args = par.parse_args()
	if args.f<0 or args.p>1:
		print "-f must be [0,1]"
		par.print_help()
	data=read_data(args.data_file,args.f)
	fit=fit_Mjj(data,args.n)

	A_sum=0
	for j in xrange(len(fit[0])/4):
		print "A"+str(j+1)+"\t",fit[0][4*j]
		A_sum-=fit[0][4*j]
		print "omega"+str(j+1)+"\t",fit[0][4*j+1]
		print "delta"+str(j+1)+"\t",fit[0][4*j+2]
		print "tau"+str(j+1)+"\t",fit[0][4*j+3]
	if len(fit[0])%4!=0:
		print "A"+str(j+2)+"\t",fit[0][-2]
		print "tau"+str(j+2)+"\t",fit[0][-1]
	
	iMjj=im_Mjj(fit[0])

	if args.e_data_file:
		l=args.e_data_file.readlines()
		temp=l[-4].split()
		eps_Md2=temp[-1]
		temp=l[-1].split()
		fac1=float(temp[1])#2epsilon_rf
		Md2=float(temp[2])
		fac2=float(temp[4])
		M_tot=Md2+iMjj
		M_tot-=float(temp[3])
		a=fac1*M_tot+fac2
		b=fac2-M_tot
		print "\n\nMd2 = ",Md2
		print "<Md>2 = ",temp[3]
		print "im[Mjj] = ",iMjj,"\n"
		print "\n"

		a=fac1*M_tot+fac2
		b=fac2-M_tot
		print "e_stat_Md2 = ",eps_Md2
		print "e_stat = ",a/b
		
		
	if args.p:
		plot_data([data[0],data[1],f_n(data[0],*fit[0])],True,False)
		plt.xlabel("time ps")
		plt.ylabel(r"Mjj e$^2$nm$^2$ps$^{-2}$")
		if args.s:
			plt.savefig("Mjj.png")
		else:
			plt.show()
		plt.close()

