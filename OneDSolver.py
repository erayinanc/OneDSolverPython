#!/usr/bin/python
# 1d scalar transport
# by Eray Inanc (eray.inanc@uni-due.de)
# transport two predefined profiles (sin and top-hat)
# three numerical schemes for convective fluxes are considered, easily extendable
# a profile for velocity can also be given instead of a fixed value 
# default settings: CFL=0.25, D=1e-5, delta=0.2e-3, u=1
# you need numpy library for computation, matplotlib library for post-proc 

# libraries
import numpy as np

# options
scheme  = 3 		# 1:cds / 2:uds / 3:TVD-blend
profile = 2 		# 1:sinus profile / 2:top-hat profile

# computation parameters
CFL     = 0.25       	# timestep condition
tmax    = 0.01       	# computation time [s]
rho 	= 1.0 		# density [kg/m3] 
D       = 1.0e-5 	# kin. diffusivity [m2/2] 
nx 	= 0.01	  	# domain lenght [m] 
delta 	= 0.2e-3    	# cell size [m] 
u 	= 1.0         	# Flow velocity (can be a profile) 
fac 	= 0.5         	# Factor between CSD and UDS for scheme 3 

# set boundaries
def initialise(phiU,psi):
	global exact,dt,dom
	# initialise domain
	dom = np.linspace(0,int(nx/delta),len(psi))
	# initialise scalar field
	if profile == 1:
		psi[:] = 0.5 * u * np.sin(np.linspace(0,2*np.pi,len(psi))) + 0.5 * u
	elif profile == 2:
		psi[midX-int(nI/10):midX+int(nI/10)] = u
	# initialise velocity field
	phiU[:] = u 
	# initialise time step width 
	dt = CFL*delta/np.amax(u)
	# exact solution for comparison
	exact = psi[:]
	return phiU,psi

# generate aliases and helper parameters
def genAlias():
	global nI,midX,Ifi,Ifip,Ifim,Ila,Ilap,Ilam
	# helper parameters 
	nI = len(np.arange(int(nx/delta)))
	midX = int(nI/2)
	Ifi  = 1
	Ifip = Ifi + 1
	Ifim = Ifi - 1
	Ila  = nI - 1
	Ilam = Ila - 1
	Ilap = Ila + 1
	return nI,midX,Ifi,Ifip,Ifim,Ila,Ilap,Ilam

# velocity at faces 
def center2surface(phiU):
	fluxU = np.zeros((int(nx/delta)))
	fluxU[Ifim:Ila] = 0.5 * (phiU[Ifi:Ilap] + phiU[Ifim:Ila])
	return fluxU

# convective fluxes
def calcFluxCon(fluxU,psi):
	fluxConX = np.zeros((int(nx/delta)))
	# cds
	if scheme == 1:
		fluxConX[Ifim:Ila] = 0.5 * fluxU[Ifim:Ila] * (psi[Ifim:Ila] + psi[Ifi:Ilap])
	# uds
	elif scheme == 2:
		fluxConX[Ifim:Ila] = fluxU[Ifim:Ila] * psi[Ifim:Ila] 
	# tvd-blend
	elif scheme == 3:
		weight = fac * ( u * dt / delta + 1 )
		fluxConX[Ifim:Ila] = fluxU[Ifim:Ila] * ( \
				( weight      ) * 0.5 * (psi[Ifim:Ila] + psi[Ifi:Ilap]) + \
				( 1.0 - weight) * psi[Ifim:Ila] )
	return fluxConX

# diffusive fluxes
def calcFluxDif(psi):
	fluxDifX = np.zeros((int(nx/delta)))
	fluxDifX[Ifim:Ila] = (D/delta) * (psi[Ifi:Ilap]-psi[Ifim:Ila])
	return fluxDifX

# plot windows
def outputIni(psi):
	global plt,fig,line
	# plot library
	import matplotlib.pyplot as plt
	# plot properties
	plt.close('all')
	plt.ion()
	fig = plt.figure() 
	plt.plot(dom, exact, 'k-', label='exact')
	line, = plt.plot(dom, psi, 'r-', label='computed')
	plt.legend(loc=1)
	plt.xlabel(r'$[mm]$')
	plt.ylabel(r'$\psi$')
	plt.ylim((-0.1,1.1))
	return plt,fig,line

# plot final
def outputFinal(psi):
	plt.ioff()
	plt.plot(dom, exact, 'k-')
	plt.plot(dom, psi, 'r-')
	plt.show()

# main program
def main():
	global dt
	# allocation
	phiU 	= np.zeros((int(nx/delta)))
	psi 	= np.zeros((int(nx/delta)))
	psiP 	= np.zeros((int(nx/delta)))
	exact 	= np.zeros((int(nx/delta)))
	# initialise
	genAlias() 			# helper parameters
	phiU,psi = initialise(phiU,psi) # flow field
	outputIni(psi) 			# output window
	# time integration
	t=0
	for n in range(int(tmax/dt)):
		print 'n: ', n
		print 'time [ms]: ', t*1000
		print '-------------------'
		# velocity fluxes
		fluxU = center2surface(phiU)
		# convective fluxes
		fluxConX = calcFluxCon(fluxU,psi) 
		# diffusive fluxes
		fluxDifX = calcFluxDif(psi) 
	 	# apply fluxes 
		psiP[Ifi:Ila] = psi[Ifi:Ila] + rho*dt/delta * ( \
				- (fluxConX[Ifi:Ila] - fluxConX[Ifim:Ilam]) + \
				  (fluxDifX[Ifi:Ila] - fluxDifX[Ifim:Ilam]) )
	     	# apply periodic boundaries, careful with indices for numpy!
	    	psiP[Ifim] = psiP[-2]
	     	psiP[-1] = psiP[Ifi]
	     	# update t and phi
	     	t = t + dt
		dt = CFL*delta/np.amax(u)
		psi = np.copy(psiP)
	     	# post-proc
		line.set_ydata(psiP)
		fig.canvas.draw()
	# final plot
	outputFinal(psi)
	print 'finished!'
''' call the main function only within the main Thread'''
if __name__ == "__main__":
	import sys
	main()


   	
