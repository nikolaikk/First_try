import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

def solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no):
    
    x = np.arange(xo,xend,x_mesh)
    z = np.arange(zo,zend,z_mesh)
    Nz = len(z)
    Nx = len(x)
    
    epsilon = 8.854187817e-12
    c = 299792458

    ko = 2*np.pi*no/Lambda
    n = 1
    k = 2*np.pi*n*np.ones(len(x))
    f = np.zeros([Nx,Nz])
    f[:,0] = np.exp(-(x/wo)**2)
    f[0,0] = 0
    f[-1,0] = 0
    
    ## Absorption
    
    region = 1/2
    k_grad = ((min(x)+max(x)*region-x[0:int(np.round(len(x)*region/2))])*region/2)**2
    k_abs = np.zeros(len(x))
    k_abs[0:len(k_grad)]=k_grad
    k_abs[-len(k_grad)-1:-1] = np.flip(k_grad,0)
    
    n = n + 0.0013j*k_abs
#    k = 2*np.pi/Lambda*n;

    ## Definition of L matrix
    
    L = np.zeros([Nx,Nx])
    for i in range(Nx):
        L[i,i]=-2+(ko**2-k[i]**2)*x_mesh**2    # in free space^2 - in material^2
    
    for i in range(1,Nx):
        L[i,i-1]=1;
        L[i-1,i]=1;
        
    L = (1/x_mesh**2)*L
    
    ## Psi evolution
    
    E = np.eye(Nx)
    const_to_simplify = (np.linalg.pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)))
    #const_to_simplify(find(const_to_simplify<1e-14))=0
    
    for i in range(Nz-1):
        f[:,i+1]= np.dot(const_to_simplify,f[:,i])
    
        f[0,i+1] = 0
        f[-1,i+1] = 0
    
        
    periodic = np.array([np.exp(-1j*ko*z),]*Nx)
    f=np.multiply(f,periodic);
    
    I = epsilon*no*c*0.5*np.abs(f)**2
    
    return f,I,x,z

def analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no):
    x = np.arange(xo,xend,x_mesh)
    z = np.arange(zo,zend,z_mesh)
    Nz = len(z)
    Nx = len(x)
    E = 1
    k = 2*np.pi*no/Lambda    
    zR = np.pi*wo**2/Lambda
    Ez = np.zeros([Nx,Nz])
    epsilon = 8.854187817e-12
    c = 299792458
    
    for i in range(Nz):
        
        w = wo*np.sqrt(1+(z[i]/zR)**2)
        r = z[i]+(zR**2/z[i])
        gouy = np.arctan(z[i]/zR)

        # Amplitude Normalized
        Ez[:,i] = E * (wo/w)**0.5*np.exp(-(x**2)/w**2)*np.exp(1j*(k*z[i]+k*(x**2)/(2*r)-gouy))   
        
        I = epsilon*no*c*0.5*np.abs(Ez)**2
    
    return Ez,I
    
# =============================================================================
# =============================================================================

## Initialization
    
t = time.time()

zo = 0                    # micrometer
zend = 600                # micrometer
z_mesh = 1                # micrometer
ratio = 1
xo = -400                 # micrometer
xend = -xo                # micrometer
x_mesh = z_mesh/ratio     # micrometer
Lambda = 10               # micrometer
wo = 25                   # micrometer
no = 1

## Start Calculation

part_of_program = 0;
part_of_program = 1;

## Convergence

if part_of_program == 0:
    step = 0.1
    from_ = 4
    to_ = 0.5
#    range_z = np.arange(from_,to_,-step)
#    range_z = np.array([2,1*np.sqrt(2),1,np.sqrt(2)*0.5,0.5,0.25*np.sqrt(2),0.25])
    range_z = np.array([1,0.75,0.5,0.375,0.25])
    
    err_l2 = np.zeros(len(range_z))
    err_max = err_l2
    err_point = err_l2
    Number_Nz = err_l2
    Number_Nz_Nx = err_l2
    ii = 0;
    
    
    Ez,I,x,z = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
    Ez_analytic, I_analytic = analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
    
    
else:
    Ez,I,x,z = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
    Ez_analytic, I_analytic = analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
    
    #       PLOT 3D
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure()
#    ax = fig.add_subplot(111,projection='3d')
    Z, X = np.meshgrid(z, x)

    plt.subplot(211)
    plt.contour(Z, X, np.real(Ez_analytic))
    plt.subplot(212)
    plt.contour(Z, X, np.real(Ez))

    
    plt.show()
    

#    for i_z_mesh in range_z:
#        x_mesh = i_z_mesh
#        
#        Ez,I,x,z = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
#        
#        Ez_analytic, I_analytic = analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
#        
#        Number_Nz_Nx[ii] = len(z)*len(x)
#        
#        Ez_diff = np.abs(np.abs(np.real(Ez))-np.abs(np.real(Ez_analytic)))
#        err_l2[ii] = 1/Number_Nz_Nx[ii] * np.sqrt(np.sum(Ez_diff**2))
#        err_max[ii] = np.max(Ez_diff)
#        ii = ii+1

        
        
        
##        print(i_z_mesh)
#    
#    
#    
#    
#    
#    zo = 0  # meter
#    zend = 0.3  # meter
#    Nz = 200
#    delta_z = (zend-zo)/Nz
#    #z=zo:delta_z:zend
#    z = np.linspace(zo,zend,Nz)
#    
#    Lambda = 500e-9
#    no=1
#    n=1
#    k=2*np.pi*n/Lambda
#    ko=2*np.pi*no/Lambda
#    
#    xo=-0.001
#    xend=-xo
#    delta_x=(xend-xo)/(Nz-1)
#    #x=np.arange(xo,xend,delta_x)
#    x = np.linspace(xo,xend,Nz)
#    wo=0.0001
#    f=np.zeros([Nz,Nz])
#    f[:,0] = np.exp(-(x/wo)**2)
#
# ## Define L matrix
#
#L = np.zeros([Nz,Nz]);
#for i in range(0,Nz):
#     L[i,i]=-2+(ko**2-k**2)*delta_x**2 # in free space^2 - in material^2
#
#for i in range (1,Nz):
#     L[i,i-1]=1
#     L[i-1,i]=1
#
#L = (1/delta_x**2)*L
#
#
#E = np.eye(Nz)
#
## #   Psi evolution
##import scipy.io
##import sys
##mat = scipy.io.loadmat('const_to_simplify.mat')
##const_to_simplify = mat['const_to_simplify']
#
#const_to_simplify = np.dot((E+(delta_z*L/(4*1j*ko))), np.linalg.inv (E-(delta_z*L/(4*1j*ko))))
#const_to_simplify[np.where(const_to_simplify<1e-15)] = 0
#
#for i in range(0, Nz-1):
##for i in range(0, 2):
#    f[:,i+1] = np.float128(np.dot(const_to_simplify,  f[:,i]))
#
#elapsed = time.time()-t
#print('Elapsed time : ',elapsed)
#
##f = scipy.io.loadmat('f')['f']
#I = np.multiply(np.conjugate(f),f)
#
#
###============================================================================
##       PLOT 3D
###============================================================================
#import matplotlib.gridspec as gridspec
#from mpl_toolkits.mplot3d import Axes3D
#
#gs = gridspec.GridSpec(1,5)
#fig1 = plt.subplot(gs[0,0:2])
#plt.imshow(np.real(I))
#plt.xlabel('z (m)')
#plt.ylabel('x (mm)')
#fig1.set_xticklabels(   fig1.get_xticks()*delta_z   )
#fig1.set_yticklabels(   fig1.get_yticks()*delta_x*1000   )
#
#fig2 = plt.subplot(gs[0,2:6],projection='3d')
#u = np.linspace(0,  Nz, Nz)
#x = np.outer(np.ones(np.size(u)), u)
#y = x.T
#z = np.multiply(np.conjugate(f),f)
#fig2.plot_surface(x, y, z, cmap = 'RdBu_r')
#
#plt.show()