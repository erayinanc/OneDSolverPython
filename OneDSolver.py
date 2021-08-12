# 1d scalar transport in serial
# v.210812a
# contact eray.inanc@uni-due.de
# transport two predefined profiles (sin and top-hat)
# three numerical schemes for convective fluxes
# top-hat velocity profile 
# default settings: CFL=0.25, D=1e-5 m2/s, delta=0.2e-3 m, u=1 m/s

# libraries
import numpy as np, time

# options
scheme  = 3    # 1:cds / 2:uds / 3:TVD-blend
profile = 2    # 1:sinus profile / 2:top-hat profile
postp = False  # post-processing?

# computation parameters
CFL   = 0.25   # timestep condition
tmax  = 0.0001 # computation time [s]
rho   = 1.0    # density [kg/m3] 
D     = 1.0e-5 # kin. diffusivity [m2/2] 
nx    = 0.01   # domain lenght [m] 
delta = 0.2e-3 # cell size [m] 
u     = 1.0    # Flow velocity (can be a profile) 
fac   = 0.5    # Factor between CSD and UDS for scheme 3 

# set boundaries
def initialise():
    # num of elements in subdomain, assign globally
    global nI,sT
    nI = int(nx/delta)

    # starting time
    sT = time.time()

    # allocation
    phiU  = np.zeros(nI)
    psi   = np.zeros(nI)

    # initialise domain
    dom = np.linspace(0,nx,nI)

    # middle of subdomain
    midX = int(nI/2)

    # initialise scalar field
    if profile == 1:
        psi[:] = 0.5 * u * np.sin(np.linspace(0,2*np.pi,nI)) + 0.5 * u
    elif profile == 2:
        psi[midX-int(nI/10):midX+int(nI/10)] = u

    # initialise velocity field
    phiU[:] = u

    # exact solution for comparison
    exact = np.copy(psi)

    # output window
    if postp:
        outputIni(psi,dom,exact)

    return phiU,psi,exact,dom

# velocity at faces 
def center2surface(phiU):
    flux = np.zeros(nI)
    flux[0:-1] = 0.5 * (phiU[1::] + phiU[0:-1])

    return flux

# convective fluxes
def calcFluxC(fluxU,psi,dt):
    flux = np.zeros(nI)

    # cds
    if scheme == 1:
        flux[0:-1] = 0.5 * fluxU[0:-1] * (psi[0:-1] + psi[1::])

    # uds
    elif scheme == 2:
        flux[0:-1] = fluxU[0:-1] * psi[0:-1]

    # blend
    elif scheme == 3:
        weight = fac * ( u * dt / delta + 1 )
        flux[0:-1] = fluxU[0:-1] * ( \
                ( weight      ) * 0.5 * (psi[0:-1] + psi[1::]) + \
                ( 1.0 - weight) * psi[0:-1] )

    return flux

# diffusive fluxes
def calcFluxD(psi):
    flux = np.zeros(nI)
    flux[0:-1] = (D/delta) * (psi[1::]-psi[0:-1])

    return flux

# plot windows
def outputIni(psi,dom,exact):
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

# plot final
def outputFinal(psi,exact,dom):
    plt.ioff()
    plt.plot(dom, exact, 'k-')
    plt.plot(dom, psi, 'r-')
    plt.show()

# main program
def main():
    # initialise
    phiU,psi,exact,dom = initialise()

    # time integration
    t=0
    n=0
    print(f'------------------------')
    print(f'starting')
    print(f'------------------------\n')
    while t < tmax:
        print(f'------------------------')
        print(f'n: {n} \t '\
              f'time [ms]: {t*1e3}')
        print(f'------------------------')

        # time step width
        dt = CFL*delta/np.amax(u)

        # velocity fluxes
        flux = center2surface(phiU)

        # convective fluxes
        fluxC = calcFluxC(flux,psi,dt)

        # diffusive fluxes
        fluxD = calcFluxD(psi)

        # apply fluxes 
        psi[1:-1] += rho*dt/delta * ( \
                - (fluxC[1:-1] - fluxC[0:-2]) + \
                  (fluxD[1:-1] - fluxD[0:-2]) )

        # apply periodic boundaries, careful with indices for numpy!
        psi[0] = psi[-2]
        psi[-1] = psi[1]

        # update t, n
        t+=dt
        n+=1

        # # post-proc
        if postp:
            line.set_ydata(psi)
            fig.canvas.draw()

    # final plot
    if postp:
        outputFinal(psi,exact,dom)
    print(f'\nfinished in {time.time()-sT:2.5f} seconds!')

if __name__ == "__main__":
    main()

#eof
