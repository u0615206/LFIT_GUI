import numpy as np
import scipy.signal as signal
from scipy.stats import norm
import scipy.sparse.linalg as splinalg
from scipy.sparse import csr_matrix
from scipy.sparse import diags
import scipy.special as sf
import cosmolopy.distance as cd
from scipy.spatial import Delaunay
import scipy.linalg as linalg

def xy_transform(x, y, x_cen, y_cen, phi):
	xnew=(x-x_cen)*np.cos(np.pi*phi/180.0)+(y-y_cen)*np.sin(np.pi*phi/180.0)
	ynew=-(x-x_cen)*np.sin(np.pi*phi/180.0)+(y-y_cen)*np.cos(np.pi*phi/180.0)
	return (xnew, ynew)
	
def gauss_2d(x, y, spar):
	(xnew, ynew)=xy_transform(x, y, spar[1], spar[2], spar[4])
	r=np.sqrt(spar[5]*xnew**2.+ynew**2./spar[5])
	return spar[0]*np.exp(-0.5*(r/spar[3])**spar[6])

#def sersic_2d(x, y, spar):
#    (xnew, ynew)=xy_transform(x, y, spar[1], spar[2], spar[4])
#    r=np.sqrt(spar[5]*xnew**2.+ynew**2./spar[5])
#    return spar[0]*np.exp(-0.5*(r/spar[3])**spar[6])

def sersic_2d(x, y, spar):
# note the position angle in the sersic_2d is off by 90 degrees 
# compared to sersic_phot
	(xnew, ynew)=xy_transform(x, y, spar[1], spar[2], spar[4])
	n=spar[6]
	if n >= 0.36: # from Ciotti & Bertin 1999, truncated to n^-3
		k=2.0*n-1./3+4./(405.*n)+46./(25515.*n**2.)+131./(1148175.*n**3.)
	else: # from MacArthur et al. 2003
		k=0.01945-0.8902*n+10.95*n**2.-19.67*n**3.+13.43*n**4.
	r=np.sqrt(spar[5]*xnew**2.+ynew**2./spar[5])
	return spar[0]*np.exp(-k*(r/spar[3])**(1./n))

def sersic_phot(x, y, par):
	(xnew, ynew)=xy_transform(x, y, par[1], par[2], par[4])
	n=par[6]
	if n >= 0.36: # from Ciotti & Bertin 1999, truncated to n^-3
		k=2.0*n-1./3+4./(405.*n)+46./(25515.*n**2.)+131./(1148175.*n**3.)
	else: # from MacArthur et al. 2003
		k=0.01945-0.8902*n+10.95*n**2.-19.67*n**3.+13.43*n**4.
	r=np.sqrt(xnew**2./par[5]+par[5]*ynew**2.)
	return par[7]+par[0]*np.exp(-k*(r/par[3])**(1./n))

def csersic_phot(x, y, par):
	n=par[6]
	(xnew, ynew)=xy_transform(x, y, par[1], par[2], par[4])
	if n >= 0.36: # from Ciotti & Bertin 1999, truncated to n^-3
		k=2.0*n-1./3+4./(405.*n)+46./(25515.*n**2.)+131./(1148175.*n**3.)
	else: # from MacArthur et al. 2003
		k=0.01945-0.8902*n+10.95*n**2.-19.67*n**3.+13.43*n**4.
	r=np.sqrt(xnew**2./par[5]+par[5]*ynew**2.)
	r_c=par[7]
	r_e=par[3]
	return par[10]+par[0]*(1.+(r_c/r)**par[8])**(par[9]/par[8])*np.exp(-k*((r**par[8]+r_c**par[8])/r_e**par[8])**(1./(par[6]*par[8])))

def hernquist_phot(x, y, par):
	(xnew, ynew)=xy_transform(x, y, par[1], par[2], par[4])
	r=np.sqrt(xnew**2./par[5]+par[5]*ynew**2.)
	r_s=par[3]
	return par[6]+0.5*par[0]*r_s/np.pi/r/(r+r_s)**3.
		
def def_pm(x, y, lpar):
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], 0.0)
	r_sie=np.sqrt(xnew**2.+ynew**2.)
	alpha_x=xnew/(r_sie**2.+(r_sie==0))
	alpha_y=ynew/(r_sie**2.+(r_sie==0))
	return (lpar[0]**2.*alpha_x, lpar[0]**2.*alpha_y)

def def_sie(x, y, lpar):
# Calculating the deflection angle of an SIE mass profile following Kormann 1993
# The convergence has the form of kappa(x, y)=0.5*sqrt(q)*b_sie/sqrt(x^2+q^2*y^2)
# In this form, b_sie is the Einstein radius in the intermediate-axis convention
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	if lpar[3] > 180.0:
		lpar[3]=lpar[3]-180.0
	elif lpar[3] < 0.0:
		lpar[3]=lpar[3]+180.0   
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], lpar[3])
	r_sie=np.sqrt(xnew**2.+ynew**2.)
	qfact=np.sqrt((1.0/lpar[4]-lpar[4]))
	eps=10.**(-8)
	if np.abs(qfact) <= eps:
		alpha_x=xnew/(r_sie+(r_sie==0))
		alpha_y=ynew/(r_sie+(r_sie==0))
	else:
		#delta=np.sqrt(xnew**2.+(lpar[4]*ynew)**2.)/r_sie
		alpha_x=np.arcsinh(np.sqrt(1.0/lpar[4]**2.0-1.0)*xnew/(r_sie+(r_sie==0)))/qfact
		alpha_y=np.arcsin(np.sqrt(1.0-lpar[4]**2.0)*ynew/(r_sie+(r_sie==0)))/qfact   
	(alpha_x_new, alpha_y_new)=xy_transform(alpha_x, alpha_y, 0.0, 0.0, -lpar[3])
	return (lpar[0]*alpha_x_new, lpar[0]*alpha_y_new)

#def def_sple_old(x, y, lpar):
## Calculating the deflection angle of an SPLE mass profile 
## following Tessore & Metcalf 2015 (arXiv:1507.01819)
## The convergence has the form of kappa(x, y)=0.5*(2-t)*(b/sqrt(q^2*x^2+y^2))^t
## In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention
#	if lpar[4] > 1.0:
#		lpar[4]=1.0/lpar[4]
#		lpar[3]=lpar[3]+90.0    
#	phi_sple=lpar[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
#	if phi_sple > 180.0:
#		phi_sple=phi_sple-180.0
#	elif phi_sple < 0.0:
#		phi_sple=phi_sple+180.0
#	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], phi_sple)
#
#	if np.size(x)==1:
#		n_x=1
#		n_y=1
#		q=lpar[4]
#		t=lpar[5]
#		f=(1.0-q)/(1.0+q)
#		b=lpar[0]*np.sqrt(q)
#		
#		if t==0:
#			print 'WARNING: the profile corresponds to a constant mass sheet'
#			print 'WARNING: the deflection angle is not usable'
#			alpha_sple_x=2.0*q/(1.0+q)*xnew
#			alpha_sple_y=2.0/(1.0+q)*ynew
#		elif t==2:
#			r=np.sqrt(xnew**2.+ynew**2.)
#			alpha_sple_x=lpar[0]**2.*xnew/(r**2.+(r==0.))
#			alpha_sple_y=lpar[0]**2.*ynew/(r**2.+(r==0.))
#		else:
#			alpha_sple_x=0.0*x
#			alpha_sple_y=0.0*y		
#			phi=np.arctan2(ynew, q*xnew)
#			R=np.sqrt(q**2.*xnew**2.+ynew**2.)
#			z=complex(np.cos(phi), np.sin(phi))
#			tmp=2.0*b/(1.0+q)*(b/R)**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
#			alpha_sple_x=tmp.real
#			alpha_sple_y=tmp.imag
#	else:
#		n_x=x.shape[0]
#		n_y=x.shape[1]
#
#		q=lpar[4]
#		t=lpar[5]
#		f=(1.0-q)/(1.0+q)
#		b=lpar[0]*np.sqrt(q)
#		
#		if t==0:
#			print 'WARNING: the profile corresponds to a constant mass sheet'
#			print 'WARNING: the deflection angle is not usable'
#			alpha_sple_x=2.0*q/(1.0+q)*xnew
#			alpha_sple_y=2.0/(1.0+q)*ynew
#		elif t==2:
#			r=np.sqrt(xnew**2.+ynew**2.)
#			alpha_sple_x=lpar[0]**2.*xnew/(r**2.+(r==0.))
#			alpha_sple_y=lpar[0]**2.*ynew/(r**2.+(r==0.))
#		else:
#			alpha_sple_x=0.0*x
#			alpha_sple_y=0.0*y		
#			phi=np.arctan2(ynew, q*xnew)
#			R=np.sqrt(q**2.*xnew**2.+ynew**2.)
#			for i in np.arange(n_x):
#				for j in np.arange(n_y):
#					z=complex(np.cos(phi[i, j]), np.sin(phi[i, j]))
#					tmp=2.0*b/(1.0+q)*(b/R[i, j])**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
#					alpha_sple_x[i, j]=tmp.real
#					alpha_sple_y[i, j]=tmp.imag
#	
#	(alpha_sple_xnew, alpha_sple_ynew)=xy_transform(alpha_sple_x, alpha_sple_y, 0.0, 0.0, -phi_sple)
#	return (alpha_sple_xnew, alpha_sple_ynew)

def def_sple(x, y, lpar):
# Calculating the deflection angle of an SPLE mass profile 
# following Tessore & Metcalf 2015 (arXiv:1507.01819)
# The convergence has the form of kappa(x, y)=0.5*(2-t)*(b/sqrt(q^2*x^2+y^2))^t
# In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	phi_sple=lpar[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
	if phi_sple > 180.0:
		phi_sple=phi_sple-180.0
	elif phi_sple < 0.0:
		phi_sple=phi_sple+180.0
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], phi_sple)

	q=lpar[4]
	t=lpar[5]
	f=(1.0-q)/(1.0+q)
	b=lpar[0]*np.sqrt(q)
	
	if t==0:
		print 'WARNING: the profile corresponds to a constant mass sheet'
		print 'WARNING: the deflection angle is not usable'
		alpha_sple_x=2.0*q/(1.0+q)*xnew
		alpha_sple_y=2.0/(1.0+q)*ynew
	elif t==2:
		r=np.sqrt(xnew**2.+ynew**2.)
		alpha_sple_x=lpar[0]**2.*xnew/(r**2.+(r==0.))
		alpha_sple_y=lpar[0]**2.*ynew/(r**2.+(r==0.))
	else:
		phi=np.arctan2(ynew, q*xnew)
		R=np.sqrt(q**2.*xnew**2.+ynew**2.)
		z=np.cos(phi)*(1+0j)+np.sin(phi)*(0+1j)
		tmp=2.0*b/(1.0+q)*(b/R)**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
		alpha_sple_x=tmp.real
		alpha_sple_y=tmp.imag
		
	(alpha_sple_xnew, alpha_sple_ynew)=xy_transform(alpha_sple_x, alpha_sple_y, 0.0, 0.0, -phi_sple)
	return (alpha_sple_xnew, alpha_sple_ynew)
	
def def_sple_ppn(x, y, lpar):
# Calculating the deflection angle of an SPLE mass profile 
# following Tessore & Metcalf 2015 (arXiv:1507.01819)
# The convergence has the form of kappa(x, y)=0.5*(2-t)*(1+gamma_ppn)/2*(b/sqrt(q^2*x^2+y^2))^t
# In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	phi_sple=lpar[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
	if phi_sple > 180.0:
		phi_sple=phi_sple-180.0
	elif phi_sple < 0.0:
		phi_sple=phi_sple+180.0
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], phi_sple)

	if np.size(x)==1:
		n_x=1
		n_y=1
		q=lpar[4]
		t=lpar[5]
		gamma_ppn=lpar[6]
		f=(1.0-q)/(1.0+q)
		b=lpar[0]*np.sqrt(q)
		
		if t==0:
			print 'WARNING: the profile corresponds to a constant mass sheet'
			print 'WARNING: the deflection angle is not usable'
			alpha_sple_x=2.0*q/(1.0+q)*xnew
			alpha_sple_y=2.0/(1.0+q)*ynew
		elif t==2:
			r=np.sqrt(xnew**2.+ynew**2.)
			alpha_sple_x=0.5*(1+gamma_ppn)*lpar[0]**2.*xnew/(r**2.+(r==0.))
			alpha_sple_y=0.5*(1+gamma_ppn)*lpar[0]**2.*ynew/(r**2.+(r==0.))
		else:
			alpha_sple_x=0.0*x
			alpha_sple_y=0.0*y		
			phi=np.arctan2(ynew, q*xnew)
			R=np.sqrt(q**2.*xnew**2.+ynew**2.)
			z=complex(np.cos(phi), np.sin(phi))
			tmp=2.0*b/(1.0+q)*(b/R)**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
			alpha_sple_x=0.5*(1+gamma_ppn)*tmp.real
			alpha_sple_y=0.5*(1+gamma_ppn)*tmp.imag
	else:
		n_x=x.shape[0]
		n_y=x.shape[1]

		q=lpar[4]
		t=lpar[5]
		f=(1.0-q)/(1.0+q)
		gamma_ppn=lpar[6]
		b=lpar[0]*np.sqrt(q)
		
		if t==0:
			print 'WARNING: the profile corresponds to a constant mass sheet'
			print 'WARNING: the deflection angle is not usable'
			alpha_sple_x=2.0*q/(1.0+q)*xnew
			alpha_sple_y=2.0/(1.0+q)*ynew
		elif t==2:
			r=np.sqrt(xnew**2.+ynew**2.)
			alpha_sple_x=0.5*(1+gamma_ppn)*lpar[0]**2.*xnew/(r**2.+(r==0.))
			alpha_sple_y=0.5*(1+gamma_ppn)*lpar[0]**2.*ynew/(r**2.+(r==0.))
		else:
			alpha_sple_x=0.0*x
			alpha_sple_y=0.0*y		
			phi=np.arctan2(ynew, q*xnew)
			R=np.sqrt(q**2.*xnew**2.+ynew**2.)
			for i in np.arange(n_x):
				for j in np.arange(n_y):
					z=complex(np.cos(phi[i, j]), np.sin(phi[i, j]))
					tmp=2.0*b/(1.0+q)*(b/R[i, j])**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
					alpha_sple_x[i, j]=0.5*(1+gamma_ppn)*tmp.real
					alpha_sple_y[i, j]=0.5*(1+gamma_ppn)*tmp.imag
	
	(alpha_sple_xnew, alpha_sple_ynew)=xy_transform(alpha_sple_x, alpha_sple_y, 0.0, 0.0, -phi_sple)
	return (alpha_sple_xnew, alpha_sple_ynew)

def def_softie(x, y, lpar):
# Calculating the deflection angle of an softened isothermal ellipsoid mass profile 
# following Keeton & Kochanek 1998 (ApJ 495, 157)
# The convergence in Keeton's paper has the form of kappa(x, y)=0.5*b/sqrt(q^2*s^2+q^2*x^2+y^2)
# However, because of the extra 90 degree in the coordinate transformation, 
# the convergence associated with this deflection angle calculation is 
# kappa(x, y)=0.5*b/sqrt(q^2*s^2+x^2+q^2*y^2)

# q*s is the core radius
# In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention,  
# in the limit of s->0
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	if lpar[3] > 180.0:
		lpar[3]=lpar[3]-180.0
	elif lpar[3] < 0.0:
		lpar[3]=lpar[3]+180.0
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], lpar[3]+90.) # the extra 90 is due to the axis-flip between Keeton & Kochanek 1998 and Kormann 1993

	q=lpar[4]
	b=lpar[0]*np.sqrt(q)
	s=lpar[5]
	eps=10.**(-8)
	psi=np.sqrt(q**2*(s**2+xnew**2)+ynew**2)
	qfact=np.sqrt(1.0-q**2.)
	if  (qfact <= eps):
		alpha_x=b*xnew/(psi+s+(psi+s == 0))
		alpha_y=b*ynew/(psi+s+(psi+s == 0))
	else:
		alpha_x=b/qfact*np.arctan(qfact*xnew/(psi+s+(psi+s == 0)))
		alpha_y=b/qfact*np.arctanh(qfact*ynew/(psi+q**2.*s+(psi+q**2.*s == 0)))

	(alpha_x_new, alpha_y_new)=xy_transform(alpha_x, alpha_y, 0.0, 0.0, -lpar[3]-90.)
	return (alpha_x_new, alpha_y_new)

def mathcal_F(x):
# calculte the special function F(x) defined in Bartelmann 1996, AA, 313, 697
# and Golse & Kneib 2002, AA, 390, 821, especially in the limit of x==1, 
# take the correct 1/3 from Golse & Kneib 2002. 
	if np.size(x)==1:
		if x==1.:
			return 1./3.
		elif x > 1.:
			return (1.0-2.*np.arctan(np.sqrt((x-1.)/(1.+x)))/np.sqrt(x**2.-1.))/(x**2.-1.)
		elif x < 1:
			return (1.0-2.*np.arctanh(np.sqrt((1.-x)/(1.+x)))/np.sqrt(1-x**2.))/(x**2.-1.)
	elif np.size(x) > 1:
		f=np.zeros_like(x)
		wh_eq1=np.where(x==1.)
		wh_gt1=np.where(x>1.)
		wh_lt1=np.where(x<1.)
		f[wh_eq1]=1./3.
		f[wh_gt1]=(1.0-2.*np.arctan(np.sqrt((x[wh_gt1]-1.)/(1.+x[wh_gt1])))/np.sqrt(x[wh_gt1]**2.-1.))/(x[wh_gt1]**2.-1.)
		f[wh_lt1]=(1.0-2.*np.arctanh(np.sqrt((1.-x[wh_lt1])/(1.+x[wh_lt1])))/np.sqrt(1-x[wh_lt1]**2.))/(x[wh_lt1]**2.-1.)
		return f

def g_snfw(x):
	if np.size(x)==1:
		if x==1.:
			return np.log(0.5*x)+1.
		elif x > 1.:
			return np.log(0.5*x)+2.*np.arctan(np.sqrt((x-1.)/(1.+x)))/np.sqrt(x**2.-1.)
		elif x < 1:
			return np.log(0.5*x)+2.*np.arctanh(np.sqrt((1.-x)/(1.+x)))/np.sqrt(1-x**2.)
	elif np.size(x) > 1:
		g_return=np.zeros_like(x)
		wh_eq1=np.where(x==1.)
		wh_gt1=np.where(x>1.)
		wh_lt1=np.where(x<1.)
		g_return[wh_eq1]=np.log(0.5*x[wh_eq1])+1.
		g_return[wh_gt1]=np.log(0.5*x[wh_gt1])+2.*np.arctan(np.sqrt((x[wh_gt1]-1.)/(1.+x[wh_gt1])))/np.sqrt(x[wh_gt1]**2.-1.)
		g_return[wh_lt1]=np.log(0.5*x[wh_lt1])+2.*np.arctanh(np.sqrt((1.-x[wh_lt1])/(1.+x[wh_lt1])))/np.sqrt(1-x[wh_lt1]**2.)
		return g_return

def def_snfw(x, y, lpar):
# Calculating the deflection angle of a SPHERICAL NFW mass profile 
# following Bartelmann 1996, AA, 313, 697 and Golse & Kneib 2002, AA, 390, 821
# The NFW profile has the form of rho(r)=lpar[0]*sigma_crit/(r*(1+r/rs)^2.)
	rs=lpar[5]
	if lpar[3] > 180.0:
		lpar[3]=lpar[3]-180.0
	elif lpar[3] < 0.0:
		lpar[3]=lpar[3]+180.0
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], lpar[3])
	phicoord=np.arctan2(ynew, xnew)
	r_new=np.sqrt(xnew**2.+ynew**2.)
	reduced_r=r_new/rs
	alpha_x=4.0*lpar[0]*r_new*g_snfw(reduced_r)/reduced_r**2.*np.cos(phicoord)
	alpha_y=4.0*lpar[0]*r_new*g_snfw(reduced_r)/reduced_r**2.*np.sin(phicoord)
	return (alpha_x, alpha_y)

def def_shear(x, y, shear_par):
# The potential of the external shear has the form of 
# psi=-0.5*shear*r^2*cos2(phi-phi_shear)
# Note that the shear position angle is off by 90 deg 
# from others'
    shear=shear_par[0]
    phig=np.deg2rad(shear_par[1])
    phicoord=np.arctan2(y, x)
    rcoord=np.sqrt(x**2.+y**2.)
    # dpsi/dr
    dpsi_1=-rcoord*shear*np.cos(2.*(phicoord-phig))
    # 1/r * dpis/dphi
    dpsi_2=rcoord*shear*np.sin(2.*(phicoord-phig))

    alpha_x=dpsi_1*np.cos(phicoord)-dpsi_2*np.sin(phicoord)
    alpha_y=dpsi_1*np.sin(phicoord)+dpsi_2*np.cos(phicoord)
    return (alpha_x, alpha_y)

def def_pix(x, y, kappa):
   """
   Function to take a convergence map and to generate
   a lensing potential gradient function from it, with
   one free scaling parameter.

   Return arrays have deflections in units of pixels.
   """
   # Create centered kernels with which to convolve:
   nx, ny = kappa.shape
   dpix_x=x[0, 1]-x[0, 0]
   dpix_y=y[1, 0]-y[0, 0]
   ximage=np.outer(np.ones(2*ny+1), np.arange(2*nx+1)-nx)*dpix_x
   yimage=np.outer(np.arange(2*ny+1)-ny, np.ones(2*nx+1))*dpix_y  
   xkernel = (1./np.pi) * ximage / ((ximage**2 + yimage**2) + (ximage == 0))
   ykernel = (1./np.pi) * yimage / ((ximage**2 + yimage**2) + (yimage == 0))
   # make a padded version of kappa:
   kappa_pad = np.zeros((3*nx, 3*ny))
   kappa_pad[nx:2*nx,ny:2*ny] = kappa
   # Do convolutions for deflections:
   xdefl = signal.fftconvolve(kappa_pad, xkernel, mode='same')*(dpix_x*dpix_y)
   ydefl = signal.fftconvolve(kappa_pad, ykernel, mode='same')*(dpix_x*dpix_y)
   # For debugging, the divergence of the deflection field:
   #div_defl = (xdefl[ny:2*ny,nx+1:2*nx+1] - xdefl[ny:2*ny,nx-1:2*nx-1] +
   #            ydefl[ny+1:2*ny+1,nx:2*nx] - ydefl[ny-1:2*ny-1,nx:2*nx]) / 2.
   alpha_x=xdefl[nx:2*nx,ny:2*ny]
   alpha_y=ydefl[nx:2*nx,ny:2*ny]
   return (alpha_x, alpha_y)
   # Need to make local copies, since they can be changed at the caller level:
   #xbase_loc = copy.deepcopy(x)
   #ybase_loc = copy.deepcopy(y)
   #def ltm_grad(x, y, par):
   #    xg_out = par[0] * bilinear(x, y, xdefl, xbase=xbase_loc, ybase=ybase_loc)
   #    yg_out = par[0] * bilinear(x, y, ydefl, xbase=xbase_loc, ybase=ybase_loc)
   #    return (xg_out, yg_out)
   #return ltm_grad

def def_total(x, y, lpar, lens_model, n_lens=1):
# p is a list of dictionaries containing the lensing parameters
	p=make_dict_lens(lpar, lens_model, n_lens)
	alpha_x=0.0*x
	alpha_y=0.0*y
	for i_lens in np.arange(n_lens):
		lens_model_i=lens_model[i_lens]
		lpar_i_dict=p[i_lens]
		if lens_model_i == 'sie':
			lpar_i=[lpar_i_dict['b_SIE'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen'], lpar_i_dict['lens_pa'], lpar_i_dict['lens_q'], lpar_i_dict['lens_gamma']]
			(alpha_x_i, alpha_y_i)=def_sie(x, y, lpar_i)
		elif lens_model_i == 'sple':
			lpar_i=[lpar_i_dict['b_SIE'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen'], lpar_i_dict['lens_pa'], lpar_i_dict['lens_q'], lpar_i_dict['lens_gamma']]
			(alpha_x_i, alpha_y_i)=def_sple(x, y, lpar_i)
		elif lens_model_i == 'softie':
			lpar_i=[lpar_i_dict['b_SIE'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen'], lpar_i_dict['lens_pa'], lpar_i_dict['lens_q'], lpar_i_dict['lens_rc']]		
			(alpha_x_i, alpha_y_i)=def_softie(x, y, lpar_i)
		elif lens_model_i == 'pm':
			lpar_i=[lpar_i_dict['theta_e'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen']]
			(alpha_x_i, alpha_y_i)=def_pm(x, y, lpar_i)
		elif lens_model_i == 'snfw':
			lpar_i=[lpar_i_dict['kappa_s'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen'], lpar_i_dict['lens_pa'], lpar_i_dict['lens_q'], lpar_i_dict['lens_rs']]
			(alpha_x_i, alpha_y_i)=def_snfw(x, y, lpar_i)
		elif lens_model_i == 'ext. shear':
			lpar_i=[lpar_i_dict['shear_amp'], lpar_i_dict['shear_pa']]
			(alpha_x_i, alpha_y_i)=def_shear(x, y, lpar_i)
		elif lens_model_i == 'sple_ppn':
			lpar_i=[lpar_i_dict['b_SIE'], lpar_i_dict['lens_xcen'], lpar_i_dict['lens_ycen'], lpar_i_dict['lens_pa'], lpar_i_dict['lens_q'], lpar_i_dict['lens_gamma'], lpar_i_dict['lens_gamma_ppn']]
			(alpha_x_i, alpha_y_i)=def_sple_ppn(x, y, lpar_i)			
		alpha_x=alpha_x+alpha_x_i
		alpha_y=alpha_y+alpha_y_i
	return (alpha_x, alpha_y)

def make_dict_lens(lpar, lens_model, n_lens):
	n_lens_par_dict={'sie': 6, 'sple': 6, 'softie':6, 'pm': 3, 'snfw':6, 'ext. shear': 2, 'sple_ppn':7}
	p=[]
	current_index=0
	for i_lens in np.arange(n_lens):
		lens_model_i=lens_model[i_lens]
		n_par_i=n_lens_par_dict[lens_model_i]
		lpar_i=lpar[current_index:current_index+n_par_i]
		if lens_model_i == 'sie' or lens_model_i == 'sple':
			p_i=dict([('b_SIE', lpar_i[0]), ('lens_xcen', lpar_i[1]), ('lens_ycen', lpar_i[2]), ('lens_pa', lpar_i[3]), ('lens_q', lpar_i[4]), ('lens_gamma', lpar_i[5])])
		elif lens_model_i == 'softie':
			p_i=dict([('b_SIE', lpar_i[0]), ('lens_xcen', lpar_i[1]), ('lens_ycen', lpar_i[2]), ('lens_pa', lpar_i[3]), ('lens_q', lpar_i[4]), ('lens_rc', lpar_i[5])])
		elif lens_model_i == 'pm':
			p_i=dict([('theta_e', lpar_i[0]), ('lens_xcen', lpar_i[1]), ('lens_ycen', lpar_i[2])])
		elif lens_model_i == 'snfw':
			p_i=dict([('kappa_s', lpar_i[0]), ('lens_xcen', lpar_i[1]), ('lens_ycen', lpar_i[2]), ('lens_pa', lpar_i[3]), ('lens_q', lpar_i[4]), ('lens_rs', lpar_i[5])])
		elif lens_model_i == 'ext. shear':
			p_i=dict([('shear_amp', lpar_i[0]), ('shear_pa', lpar_i[1])])
		elif lens_model_i == 'sple_ppn':
			p_i=dict([('b_SIE', lpar_i[0]), ('lens_xcen', lpar_i[1]), ('lens_ycen', lpar_i[2]), ('lens_pa', lpar_i[3]), ('lens_q', lpar_i[4]), ('lens_gamma', lpar_i[5]), ('lens_gamma_ppn', lpar_i[6])])
		p.append(p_i)		
		current_index=current_index+n_par_i
	return p

def make_dict_source(spar, source_model, n_source):
	n_source_par_dict={'sersic': 7, 'gaussian': 7}
	p=[]
	current_index=0
	for i_source in np.arange(n_source):
		source_model_i=source_model[i_source]
		n_par_i=n_source_par_dict[source_model_i]
		spar_i=spar[current_index:current_index+n_par_i]
		if source_model_i == 'sersic' or source_model_i == 'gaussian':
			p_i=dict([('source_amp', spar_i[0]), ('source_xcen', spar_i[1]), ('source_ycen', spar_i[2]), ('source_sigma', spar_i[3]), ('source_pa', spar_i[4]), ('source_q', spar_i[5]), ('source_n', spar_i[6])])
		p.append(p_i)		
		current_index=current_index+n_par_i
	return p

def model(x, y, lpar, spar, lens_model, source_model, PSF=np.array([0]), n_source=1, n_lens=1):
# Again, p_source is a list of dictionaries containing the source parameters
	p_source=make_dict_source(spar, source_model, n_source)

	I_fit=0.0*x
	(alpha_x, alpha_y)=def_total(x, y, lpar, lens_model, n_lens=n_lens)

	for i_s in np.arange(n_source):
		if source_model[i_s] == 'gaussian':
			spar_i_dict=p_source[i_s]
			spar_i=[spar_i_dict['source_amp'], spar_i_dict['source_xcen'], spar_i_dict['source_ycen'], spar_i_dict['source_sigma'], spar_i_dict['source_pa'], spar_i_dict['source_q'], spar_i_dict['source_n']]
			I_fit_tmp=gauss_2d(x-alpha_x, y-alpha_y, spar_i)
		elif source_model[i_s] == 'sersic':
			spar_i_dict=p_source[i_s]
			spar_i=[spar_i_dict['source_amp'], spar_i_dict['source_xcen'], spar_i_dict['source_ycen'], spar_i_dict['source_sigma'], spar_i_dict['source_pa'], spar_i_dict['source_q'], spar_i_dict['source_n']]		
			I_fit_tmp=sersic_2d(x-alpha_x, y-alpha_y, spar_i)   
		I_fit=I_fit+I_fit_tmp    
	if (PSF.any()):
		I_model=signal.fftconvolve(I_fit, PSF, mode='same')
	else:
		I_model=I_fit
	return I_model

def phot_model(x, y, ppar, photmodel, PSF=np.array([0]), n_phot=1):
	n_phot_par_dict={'sersic': 8, 'csersic': 11, 'hernquist': 7, 'bspline': 7}
	I_phot=0.0*x
	current_phot_index=0
	for i_phot in np.arange(n_phot):
		ppar_i=ppar[current_phot_index:current_phot_index+n_phot_par_dict[photmodel[i_phot]]]
		if photmodel[i_phot] == 'sersic':
			I_phot+=sersic_phot(x, y, ppar_i)
		elif photmodel[i_phot] == 'csersic':
			I_phot+=csersic_phot(x, y, ppar_i)	
		elif photmodel[i_phot] == 'hernquist':
			I_phot+=hernquist_phot(x, y, ppar_i)				
		current_phot_index+=n_phot_par_dict[photmodel[i_phot]]
	I_phot=signal.fftconvolve(I_phot, PSF, mode='same')
	return I_phot

def source_model(x, y, spar, sourcemodel, n_source=1):
	n_source_par_dict={'sersic': 7, 'gaussian': 7}
	I_source=0.0*x
	current_source_index=0
	for i_source in np.arange(n_source):
		spar_i=spar[current_source_index:current_source_index+n_source_par_dict[sourcemodel[i_source]]]
		if sourcemodel[i_source] == 'sersic':
			I_source+=sersic_2d(x, y, spar_i)
		elif sourcemodel[i_source] == 'gaussian':
			I_source+=gaussian_2d(x, y, spar_i)
		current_source_index+=n_source_par_dict[sourcemodel[i_source]]
	return I_source

def n_max(arr, n):
	indices_mask=arr.ravel().argsort()[-n:]
	indices=np.unravel_index(indices_mask, arr.shape)
	return indices

def select_used_pixels(data, mask, n, scheme=0):
	if scheme == 0: # select the m boundary pixel and n-m brightest pixels
		mask_2d=np.zeros_like(data)
		mask_2d[mask]=1.
		used_pixels_mask=np.zeros(mask[0].shape[0], dtype=int)
		# identify all the pixels in the boundary for inversion
		bnd_pixels_x=np.array([], dtype=int)
		bnd_pixels_y=np.array([], dtype=int)
		for i in range(mask[0].shape[0]):
			try:
				tmp=mask_2d[mask[0][i]-1, mask[1][i]]*mask_2d[mask[0][i], mask[1][i]+1]*mask_2d[mask[0][i]+1, mask[1][i]]*mask_2d[mask[0][i], mask[1][i]-1]
			except:
				tmp=0
			if tmp == 0:
				used_pixels_mask[i]=1
				bnd_pixels_x=np.append(bnd_pixels_x, mask[0][i])
				bnd_pixels_y=np.append(bnd_pixels_y, mask[1][i])
		
		# check to see if the code above works
		#mask_2d[bnd_pixels_x, bnd_pixels_y]=2.0
		#plt.imshow(mask_2d, **myargs)
		#plt.show()
		
		# select the other brightest n-bnd_pixels_x.shape[0] pixels used for inversion	
		junk=data.copy()
		junk[bnd_pixels_x, bnd_pixels_y]=0.0
		bright_pixels=n_max(junk*mask_2d, n-bnd_pixels_x.shape[0])
		used_pixels_mask[(junk*mask_2d)[mask].ravel().argsort()[-(n-bnd_pixels_x.shape[0]):]]=1
	
		# check to see if the code above works
		#mask_2d[bright_pixels]=3.0
		#plt.imshow(mask_2d, **myargs)
		#plt.show()
		
		used_pixels=(np.append(bnd_pixels_x, bright_pixels[0]), np.append(bnd_pixels_y, bright_pixels[1]))
		# The ultimate check to see if all the pixels needed are selected
		#mask_2d[used_pixels]=2.0
		#plt.imshow(mask_2d, **myargs)
		#plt.show()
		#mask_2d=np.zeros_like(data)
		#for i in range(used_pixels_mask.shape[0]):
		#	if used_pixels_mask[i]==1:
		#		mask_2d[lens_mask_1d[0][i], lens_mask_1d[1][i]]=1
		#
		#plt.imshow(mask_2d, **myargs)
		#plt.show()	
	elif scheme == 1: # select every other pixel
		used_pixels_mask=np.zeros(mask[0].shape[0], dtype=int)
		used_pixels_x=np.array([], dtype=int)	
		used_pixels_y=np.array([], dtype=int)	
		for i in range(mask[0].shape[0]):
			if (i % 2 == 0):
				used_pixels_x=np.append(used_pixels_x, mask[0][i])
				used_pixels_y=np.append(used_pixels_y, mask[1][i])
				used_pixels_mask[i]=1
		used_pixels=(used_pixels_x, used_pixels_y)
	return (used_pixels, used_pixels_mask)

def construct_source_grid(theta_x, theta_y, alpha_x, alpha_y, used_pixels_mask):
	A_m=np.nonzero(used_pixels_mask)[0].shape[0]
	used_pixels_indices=np.nonzero(used_pixels_mask)[0]
	points=np.zeros((A_m, 2))
	points[:, 0]=theta_x[used_pixels_indices]-alpha_x[used_pixels_indices]
	points[:, 1]=theta_y[used_pixels_indices]-alpha_y[used_pixels_indices]
#	for i in range(A_m):
#		points[i, 0]=theta_x[used_pixels_indices[i]]-alpha_x[used_pixels_indices[i]]
#		points[i, 1]=theta_y[used_pixels_indices[i]]-alpha_y[used_pixels_indices[i]]
	
	tri=Delaunay(points)
	return tri

def lens_mapping(theta_x, theta_y, alpha_x, alpha_y, lens_mask, used_pixels, used_pixels_mask, PSF_csr=np.array([0])):
	A_n=lens_mask[0].shape[0]
	A_m=used_pixels[0].size
	used_pixels_indices=np.nonzero(used_pixels_mask)[0]
	not_used_pixels_indices=np.where(used_pixels_mask == 0)[0]
	points=np.zeros((A_m, 2))
	points[:, 0]=theta_x[used_pixels_indices]-alpha_x[used_pixels_indices]
	points[:, 1]=theta_y[used_pixels_indices]-alpha_y[used_pixels_indices]
	
	tri=Delaunay(points)
	
	# plot the source pixels
	#fig, ax=plt.subplots()
	#ax.triplot(points[:,0], points[:,1], tri.simplices.copy(), 'k-')
	#ax.plot(points[:,0], points[:,1], 'k.')
	#ax.imshow(I_source, vmax=0.3, extent=ext, **myargs)
	#ax.set_xlim(source_xbase.min(), source_xbase.max(), emit=False)
	#ax.set_ylim(source_ybase.min(), source_ybase.max(), emit=False)
	#plt.show()
	
	A_mapping=np.zeros((A_n, A_m))
	for i_column in range(A_m):
		A_mapping[used_pixels_indices[i_column], i_column]=1.0
	
	for i_row in not_used_pixels_indices:
		point_p=[theta_x[i_row]-alpha_x[i_row], theta_y[i_row]-alpha_y[i_row]]
		column_a, column_b, column_c=tri.simplices[Delaunay.find_simplex(tri, point_p)]
		pa=(((points[column_a]-point_p)**2.).sum())
		pb=(((points[column_b]-point_p)**2.).sum())
		pc=(((points[column_c]-point_p)**2.).sum())
		weight_a=pb*pc/(pb*pc+pa*pc+pa*pb)
		weight_b=pa*pc/(pb*pc+pa*pc+pa*pb)
		weight_c=pa*pb/(pb*pc+pa*pc+pa*pb)
		A_mapping[i_row, column_a]=weight_a
		A_mapping[i_row, column_b]=weight_b
		A_mapping[i_row, column_c]=weight_c
	
	if (PSF_csr.ndim == 2):
		A_mapping=(PSF_csr.dot(csr_matrix(A_mapping))).toarray()

	return (A_mapping, tri)

def generate_psf_csr(lens_mask, nx, ny, PSF):
	A_n=lens_mask[0].shape[0]
	all_pixels_x=lens_mask[0]
	all_pixels_y=lens_mask[1]
	PSF_hw=(PSF.shape[0]-1)/2
	PSF_operator=np.zeros((A_n, A_n))
	for i_column in range(A_n):
		indx_xcen=all_pixels_x[i_column]
		indx_ycen=all_pixels_y[i_column]
		junk=np.zeros((nx, ny))
		junk[indx_xcen, indx_ycen]=1.0
		junk_conv=signal.fftconvolve(junk, PSF, mode='same')
		PSF_operator[:, i_column]=junk_conv[lens_mask]
	wh=np.where(PSF_operator < 10.**(-8)*PSF_operator.max())
	PSF_operator[wh]=0.0
	PSF_csr=csr_matrix(PSF_operator)
	return PSF_csr

def pix_source(data, invvar, fmask, used_pixels, used_pixels_mask, x, y, lpar, lens_model, n_lens=1, PSF_csr=np.array([0]), lam=0.0, reg_scheme=0, return_cov=0):
#	import time
#	t0=time.time()
	lens_mask=np.where(fmask == 1)
	theta_x=x[lens_mask].flatten()
	theta_y=y[lens_mask].flatten()
	(alpha_x, alpha_y)=def_total(x, y, lpar, lens_model, n_lens=n_lens)
	alpha_x_1d=alpha_x[lens_mask].flatten()
	alpha_y_1d=alpha_y[lens_mask].flatten()
	(A_mapping, tri)=lens_mapping(theta_x, theta_y, alpha_x_1d, alpha_y_1d, lens_mask, used_pixels, used_pixels_mask, PSF_csr=PSF_csr)
#	print 'Total time elapsed at 1 is', time.time()-t0, 'seconds.'
	data_1d=data[lens_mask].flatten()
	weight_1d=invvar[lens_mask].flatten()
	data_1d_norm=data_1d*np.sqrt(weight_1d)
	A_mapping_norm=A_mapping*np.sqrt(weight_1d[:, None])
#   use low-level routines scipy.linalg.blas for multiplications 
#	for better performance. Need to convert array to Fortran order 
#	use F.flags to check if it's C order or Fortran order
	F=linalg.blas.dgemm(alpha=1.0, a=A_mapping_norm.T, b=A_mapping_norm.T, trans_b=True)
	D=linalg.blas.dgemm(alpha=1.0, a=A_mapping_norm.T, b=data_1d_norm)
#	F=np.dot(A_mapping_norm.transpose(), A_mapping_norm)
#	D=np.dot(A_mapping_norm.transpose(), data_1d_norm)
	if lam !=0.0:
		if reg_scheme==0: # zeroth-order regularization
			H=np.identity(F.shape[0])
		elif reg_scheme==1: # gradient form
			H=np.identity(F.shape[0])
			(indices, indptr)=tri.vertex_neighbor_vertices
			for i_column in range(H.shape[0]):
				H[i_column, i_column]=(indices[i_column+1]-indices[i_column])
				H[indptr[indices[i_column]:indices[i_column+1]], i_column]=-1.0
		#Vec_S=splinalg.bicgstab(F+lam*H, D)
		Vec_S=linalg.solve(F+lam*H, D) # a faster solver
	else:
		#Vec_S=splinalg.bicgstab(F, D)
		Vec_S=linalg.solve(F, D, sym_pos=True)
	#Vec_S=Vec_S[0] # specifically needed if using bicgstab
	#fit_1d=np.dot(A_mapping, Vec_S)
	fit_1d=linalg.blas.dgemm(alpha=1.0, a=A_mapping.T, b=Vec_S, trans_a=True)
	fit=np.zeros_like(data)
	fit[lens_mask]=fit_1d.reshape(-1)

	if return_cov:
		if lam == 0:
			Cov_matrix=linalg.inv(F)
		else:
			R=linalg.inv(F+lam*H)
			RH=linalg.blas.dgemm(alpha=1.0, a=R, b=H.T, trans_b=True)
			Cov_matrix=R-lam*linalg.blas.dgemm(alpha=1.0, a=R, b=RH.T)
			#RH=np.dot(R, H)
			#Cov_matrix=R-lam*np.dot(R, RH.transpose())			
	else:
		Cov_matrix=None
	output=dict([('fit', fit), ('solution', (Vec_S).reshape(-1)), ('Covariance', Cov_matrix)])
#	print 'Total time elapsed at 3 is', time.time()-t0, 'seconds.'
	return output

def psi_sie(x, y, lpar):
# Analytical effective lensing potential for an SIE lens
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	if lpar[3] > 180.0:
		lpar[3]=lpar[3]-180.0
	elif lpar[3] < 0.0:
		lpar[3]=lpar[3]+180.0   
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], lpar[3])
	r_sie=np.sqrt(xnew**2.+ynew**2.)
	qfact=np.sqrt((1.0/lpar[4]-lpar[4]))
	eps=10.**(-8)
	cos_phi=xnew/(r_sie+(r_sie == 0))
	sin_phi=ynew/(r_sie+(r_sie == 0))
	if np.abs(qfact) <= eps:
		print 'q=1'
		psi=lpar[0]*r_sie
	else:
		psi=lpar[0]*r_sie*(sin_phi*np.arcsin(np.sqrt(1.0-lpar[4]**2.0)*sin_phi)+cos_phi*np.arcsinh(np.sqrt(1./lpar[4]**2.-1.)*cos_phi))/qfact
	return psi

def psi_sple(x, y, lpar):
# Calculating the lensing potential of an SPLE mass profile 
# following Tessore & Metcalf 2015 (arXiv:1507.01819)
# The convergence has the form of kappa(x, y)=0.5*(2-t)*(b/sqrt(q^2*x^2+y^2))^t
# In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention
	t=lpar[5]
	if t==0:
		print 'WARNING: the profile corresponds to a constant mass sheet'
		print 'WARNING: the lensing potential is not usable'
		return np.zeros_like(x)
	elif t==2: # point mass
		return psi_pm(x, y, lpar) # Caution: might not be correct
	else:
		if lpar[4] > 1.0:
			lpar[4]=1.0/lpar[4]
			lpar[3]=lpar[3]+90.0    
		phi_sple=lpar[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
		if phi_sple > 180.0:
			phi_sple=phi_sple-180.0
		elif phi_sple < 0.0:
			phi_sple=phi_sple+180.0
		(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], phi_sple)
		(alpha_sple_x, alpha_sple_y)=def_sple(x, y, lpar)
		(alpha_sple_xnew, alpha_sple_ynew)=xy_transform(alpha_sple_x, alpha_sple_y, 0.0, 0.0, phi_sple) # transforming back to the rotated lens plane for scalars
		return (xnew*alpha_sple_xnew+ynew*alpha_sple_ynew)/(2.0-t)

def psi_pm(x, y, lpar):
	psi=0.5*lpar[0]**2.*np.log((x-lpar[1])**2.+(y-lpar[2])**2.)
	return psi

def psi_shear(x, y, shear_par):
	r=np.sqrt(x**2.+y**2.)
	cos_phi=x/(r+(r == 0))
	sin_phi=y/(r+(r == 0))
	return 0.5*shear_par[0]*r**2.*((2.0*cos_phi**2.-1.)*np.cos(np.deg2rad(2.0*shear_par[1]))+2.0*sin_phi*cos_phi*np.sin(np.deg2rad(2.0*shear_par[1])))
	
def det_Jacobian_matrix(x, y, lpar, lens_model, n_lens=1):
	n_lens_par_dict={'sie': 6, 'sple': 6, 'pm': 6, 'softie': 6, 'snfw': 6, 'ext. shear': 2}
	psi_11=np.zeros_like(x)
	psi_12=np.zeros_like(x)
	psi_21=np.zeros_like(x)
	psi_22=np.zeros_like(x)
	ind=0
	for i_lens in range(n_lens):
		lpar_i=lpar[ind:ind+n_lens_par_dict[lens_model[i_lens]]]
		ind+=n_lens_par_dict[lens_model[i_lens]]
		if lens_model[i_lens] == 'sie':
			psi_11_tmp, psi_12_tmp, psi_21_tmp, psi_22_tmp=grad_psi_sie(x, y, lpar_i)
			psi_11+=psi_11_tmp
			psi_12+=psi_12_tmp
			psi_21+=psi_21_tmp
			psi_22+=psi_22_tmp
		elif lens_model[i_lens] == 'pm':
			psi_11_tmp, psi_12_tmp, psi_21_tmp, psi_22_tmp=grad_psi_pm(x, y, lpar_i)
			psi_11+=psi_11_tmp
			psi_12+=psi_12_tmp
			psi_21+=psi_21_tmp
			psi_22+=psi_22_tmp
		elif lens_model[i_lens] == 'ext. shear':
			psi_11_tmp, psi_12_tmp, psi_21_tmp, psi_22_tmp=grad_psi_shear(x, y, lpar_i)
			psi_11+=psi_11_tmp
			psi_12+=psi_12_tmp
			psi_21+=psi_21_tmp
			psi_22+=psi_22_tmp
	return (1.0-psi_11)*(1.0-psi_22)-psi_12*psi_21

def grad_psi_sie(x, y, lpar):
	if lpar[4] > 1.0:
		lpar[4]=1.0/lpar[4]
		lpar[3]=lpar[3]+90.0    
	if lpar[3] > 180.0:
		lpar[3]=lpar[3]-180.0
	elif lpar[3] < 0.0:
		lpar[3]=lpar[3]+180.0   
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], lpar[3])
	r_sie=np.sqrt(xnew**2.+ynew**2.)
	q=lpar[4]
	cos_phi=xnew/(r_sie+(r_sie == 0))
	sin_phi=ynew/(r_sie+(r_sie == 0))
	kappa_sie=0.5*lpar[0]*np.sqrt(q/(xnew**2.+q**2.*ynew**2.))
	psi_11=2.0*kappa_sie*sin_phi**2.
	psi_12=-2.0*kappa_sie*sin_phi*cos_phi
	psi_21=psi_12
	psi_22=2.0*kappa_sie*cos_phi**2.
	psi_11_new=np.cos(lpar[3])**2.*psi_11-2.0*np.sin(lpar[3])*np.cos(lpar[3])*psi_12+np.sin(lpar[3])**2.*psi_22
	psi_12_new=np.sin(lpar[3])*np.cos(lpar[3])*psi_11+np.cos(2.0*lpar[3])*psi_12-np.sin(lpar[3])*np.cos(lpar[3])*psi_22
	psi_21_new=psi_12_new
	psi_22_new=np.sin(lpar[3])**2.*psi_11+np.sin(2.0*lpar[3])*psi_12+np.cos(lpar[3])**2.*psi_22
	return (psi_11_new, psi_12_new, psi_21_new, psi_22_new)

def grad_psi_sple(x, y, lpar):
# Calculating the gradients of the lensing potential of an SPLE mass profile 
# following Tessore & Metcalf 2015 (arXiv:1507.01819)
# The convergence has the form of kappa(x, y)=0.5*(2-t)*(b/sqrt(q^2*x^2+y^2))^t
# In this form, b/sqrt(q) is the Einstein radius in the intermediate-axis convention
	t=lpar[5]
	if t==0:
		print 'WARNING: the profile corresponds to a constant mass sheet'
		print 'WARNING: the lensing potential is not usable'
		return np.zeros_like(x)
	elif t==2:
		return grad_psi_pm(x, y, lpar)
	else:
		if lpar[4] > 1.0:
			lpar[4]=1.0/lpar[4]
			lpar[3]=lpar[3]+90.0
		phi_sple=lpar[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
		if phi_sple > 180.0:
			phi_sple=phi_sple-180.0
		elif phi_sple < 0.0:
			phi_sple=phi_sple+180.0
		q=lpar[4]
		t=lpar[5]
		f=(1.0-q)/(1.0+q)
		b=lpar[0]*np.sqrt(q)
		(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], phi_sple) 		
		phi=np.arctan2(ynew, q*xnew)
		R=np.sqrt(q**2.*xnew**2.+ynew**2.)
		r=np.sqrt(xnew**2.+ynew**2.)
		theta=np.arctan2(ynew, xnew)
		
		gamma_sple_x=np.zeros_like(x)
		gamma_sple_y=np.zeros_like(y)		
		kappa_sple=0.5*(2-t)*(b/R)**t
		n_x=x.shape[0]
		n_y=x.shape[1]
		for i in np.arange(n_x):
			for j in np.arange(n_y):
				z=complex(np.cos(phi[i, j]), np.sin(phi[i, j]))
				z_polar=complex(np.cos(theta[i, j]), np.sin(theta[i, j]))
				alpha_tmp=2.0*b/(1.0+q)*(b/R[i, j])**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
				gamma_tmp=-z_polar**2.*kappa_sple[i, j]+(1-t)*z_polar*alpha_tmp/r[i, j]
				gamma_sple_x[i, j]=gamma_tmp.real
				gamma_sple_y[i, j]=gamma_tmp.imag

		psi_11=kappa_sple+gamma_sple_x
		psi_22=kappa_sple-gamma_sple_x
		psi_12=gamma_sple_y
		psi_21=psi_12
		psi_11_new=np.cos(phi_sple)**2.*psi_11-2.0*np.sin(phi_sple)*np.cos(phi_sple)*psi_12+np.sin(phi_sple)**2.*psi_22
		psi_12_new=np.sin(phi_sple)*np.cos(phi_sple)*psi_11+np.cos(2.0*phi_sple)*psi_12-np.sin(phi_sple)*np.cos(phi_sple)*psi_22
		psi_21_new=psi_12_new
		psi_22_new=np.sin(phi_sple)**2.*psi_11+np.sin(2.0*phi_sple)*psi_12+np.cos(phi_sple)**2.*psi_22
		return (psi_11_new, psi_12_new, psi_21_new, psi_22_new)

def grad_psi_pm(x, y, lpar):
	(xnew, ynew)=xy_transform(x, y, lpar[1], lpar[2], 0.0)
	r_sie_2=xnew**2.+ynew**2.
	psi_11=lpar[0]**2.*(ynew**2.-xnew**2.)/(r_sie_2**2.+(r_sie_2 == 0))
	psi_12=-lpar[0]**2.*2.*xnew*ynew/(r_sie_2**2.+(r_sie_2 == 0))
	psi_21=psi_12
	psi_22=lpar[0]**2.*(xnew**2.-ynew**2.)/(r_sie_2**2.+(r_sie_2 == 0))
	return (psi_11, psi_12, psi_21, psi_22)

def grad_psi_shear(x, y, shear_par):
	shear=shear_par[0]
	phi_shear=np.deg2rad(shear_par[1])
	psi_11=-shear*np.cos(2.0*phi_shear)
	psi_12=-shear*np.sin(2.0*phi_shear)
	psi_21=psi_12
	psi_22=shear*np.cos(2.0*phi_shear)
	return (psi_11, psi_12, psi_21, psi_22)
	
def tdelay(x, y, beta_x, beta_y, lpar, n_lens, lens_model, z_l, z_s):
	n_lens_par_dict={'sie': 6, 'sple': 6, 'pm': 6, 'softie': 6, 'snfw': 6, 'ext. shear': 2}
	psi=np.zeros_like(x)
	ind=0
	for i_lens in np.arange(n_lens):
		lens_model_i=lens_model[i_lens]
		lpar_i=lpar[ind:ind+n_lens_par_dict[lens_model[i_lens]]]
		ind+=n_lens_par_dict[lens_model[i_lens]]
		if lens_model_i == 'sie':
			psi+=psi_sie(x, y, lpar_i)
		elif lens_model_i == 'pm':
			psi+=psi_pm(x, y, lpar_i)
		elif lens_model_i == 'ext. shear':
			psi+=psi_shear(x, y, lpar_i)			
	cosmo={'omega_M_0':0.274, 'omega_lambda_0':0.726, 'h':0.7}
	cosmo=cd.set_omega_k_0(cosmo)
	d_l=cd.comoving_distance(z_l, **cosmo) # in Mpc
	d_s=cd.comoving_distance(z_s, **cosmo) # in Mpc
#	does NOT need a (1+z_l) because d_l is the comoving distance
	time_delay=(365.*10.**6*3.26*np.pi**2./(180.*3600.)**2.)*(d_s*d_l)/(d_s-d_l)*(0.5*(x-beta_x)**2.+0.5*(y-beta_y)**2.-psi)
	return time_delay # in days

def inverse_magnification_map(x, y, lpar, lens_model, n_lens=1):
	psi_11=0.0
	psi_12=0.0
	psi_22=0.0
	for i_lens in range(n_lens):
		lpar_i=lpar[i_lens*6:i_lens*6+6]
		lens_model_i=lens_model[i_lens]
		if lens_model_i == 'sie':
			(xnew, ynew)=xy_transform(x, y, lpar_i[1], lpar_i[2], lpar_i[3])
			r=np.sqrt(xnew**2.+ynew**2.)
			sin_phi=ynew/(r+(r==0))
			cos_phi=xnew/(r+(r==0))
			kappa_i=0.5*lpar_i[0]/np.sqrt(xnew**2./lpar_i[4]+lpar_i[4]*ynew**2.)
			psi_11+=2.0*kappa_i*sin_phi**2.
			psi_12+=-2.0*kappa_i*sin_phi*cos_phi
			psi_22+=2.0*kappa_i*cos_phi**2.
		elif lens_model_i == 'sple':
			t=lpar_i[5]
			if t==0:
				print 'WARNING: the profile corresponds to a constant mass sheet'
				print 'WARNING: the lensing potential is not usable'
				return np.zeros_like(x)
			elif t==2:
				(xnew, ynew)=xy_transform(x, y, lpar_i[1], lpar_i[2], lpar_i[3])
				r=np.sqrt(xnew**2.+ynew**2.)
				psi_11+=lpar_i[0]**2.*(ynew**2.-xnew**2.)/(r**4.+(r==0))
				psi_12+=-2.0*lpar_i[0]**2.*xnew*ynew/(r**4.+(r==0))
				psi_22+=lpar_i[0]**2.*(xnew**2.-ynew**2.)/(r**4.+(r==0))
			else:
				phi_sple=lpar_i[3]+90.0# the extra 90 is due to the axis-flip between Tessore & Metcalf 2015 and Kormann 1993
				if phi_sple > 180.0:
					phi_sple=phi_sple-180.0
				elif phi_sple < 0.0:
					phi_sple=phi_sple+180.0
				q=lpar_i[4]
				t=lpar_i[5]
				f=(1.0-q)/(1.0+q)
				b=lpar_i[0]*np.sqrt(q)
				(xnew, ynew)=xy_transform(x, y, lpar_i[1], lpar_i[2], phi_sple)
				phi=np.arctan2(ynew, q*xnew)
				R=np.sqrt(q**2.*xnew**2.+ynew**2.)
				r=np.sqrt(xnew**2.+ynew**2.)
				theta=np.arctan2(ynew, xnew)
				
				gamma_sple_x=np.zeros_like(x)
				gamma_sple_y=np.zeros_like(y)		
				kappa_sple=0.5*(2-t)*(b/R)**t
				n_x=x.shape[0]
				n_y=x.shape[1]
				for i in np.arange(n_x):
					for j in np.arange(n_y):
						z=complex(np.cos(phi[i, j]), np.sin(phi[i, j]))
						z_polar=complex(np.cos(theta[i, j]), np.sin(theta[i, j]))
						alpha_tmp=2.0*b/(1.0+q)*(b/R[i, j])**(t-1.)*z*sf.hyp2f1(1., 0.5*t, 2.-0.5*t, -f*z**2.)
						gamma_tmp=-z_polar**2.*kappa_sple[i, j]+(1-t)*z_polar*alpha_tmp/r[i, j]
						gamma_sple_x[i, j]=gamma_tmp.real
						gamma_sple_y[i, j]=gamma_tmp.imag

				psi_11=kappa_sple+gamma_sple_x
				psi_22=kappa_sple-gamma_sple_x
				psi_12=gamma_sple_y
		elif lens_model_i == 'pm':
			(xnew, ynew)=xy_transform(x, y, lpar_i[1], lpar_i[2], 0.0)
			r=np.sqrt(xnew**2.+ynew**2.)
			psi_11+=lpar_i[0]**2.*(ynew**2.-xnew**2.)/(r**4.+(r==0))
			psi_12+=-2.0*lpar_i[0]**2.*xnew*ynew/(r**4.+(r==0))
			psi_22+=lpar_i[0]**2.*(xnew**2.-ynew**2.)/(r**4.+(r==0))
		elif lens_model_i == 'ext. shear':
			psi_11+=-lpar_i[0]*np.cos(2.0*lpar_i[1])
			psi_12+=-lpar_i[0]*np.sin(2.0*lpar_i[1])
			psi_22+=lpar_i[0]*np.cos(2.0*lpar_i[1])
			
	A11=1.0-psi_11
	A12=-psi_12
	A21=-psi_12
	A22=1.0-psi_22
	return A11*A22-A12*A21

def gauss_prior(x, par1, par2):
	return par1+norm.ppf(x)*par2

def uniform_prior(x, par1, par2):
	return par1+(par2-par1)*x

def prior_A(cube, ndim, nparams):
	for i_source in np.arange(n_source):
		cube[0+i_source*7]=uniform_prior(cube[0+i_source*7], priorpars[0+i_source*7, 0], priorpars[0+i_source*7, 1])
		cube[1+i_source*7]=uniform_prior(cube[1+i_source*7], priorpars[1+i_source*7, 0], priorpars[1+i_source*7, 1])
		cube[2+i_source*7]=uniform_prior(cube[2+i_source*7], priorpars[2+i_source*7, 0], priorpars[2+i_source*7, 1])
		cube[3+i_source*7]=uniform_prior(cube[3+i_source*7], priorpars[3+i_source*7, 0], priorpars[3+i_source*7, 1])
		cube[4+i_source*7]=uniform_prior(cube[4+i_source*7], priorpars[4+i_source*7, 0], priorpars[4+i_source*7, 1])
		cube[5+i_source*7]=uniform_prior(cube[5+i_source*7], priorpars[5+i_source*7, 0], priorpars[5+i_source*7, 1])
		cube[6+i_source*7]=uniform_prior(cube[6+i_source*7], priorpars[6+i_source*7, 0], priorpars[6+i_source*7, 1])
	#endfor
	for i_lens in np.arange(n_lens):
		cube[0+n_source*7+i_lens*6]=gauss_prior(cube[0+n_source*7+i_lens*6], priorpars[0+n_source*7+i_lens*6, 0], priorpars[0+n_source*7+i_lens*6, 1])
		cube[1+n_source*7+i_lens*6]=uniform_prior(cube[1+n_source*7+i_lens*6], priorpars[1+n_source*7+i_lens*6, 0], priorpars[1+n_source*7+i_lens*6, 1])
		cube[2+n_source*7+i_lens*6]=uniform_prior(cube[2+n_source*7+i_lens*6], priorpars[2+n_source*7+i_lens*6, 0], priorpars[2+n_source*7+i_lens*6, 1])
		cube[3+n_source*7+i_lens*6]=uniform_prior(cube[3+n_source*7+i_lens*6], priorpars[3+n_source*7+i_lens*6, 0], priorpars[3+n_source*7+i_lens*6, 1])
		cube[4+n_source*7+i_lens*6]=uniform_prior(cube[4+n_source*7+i_lens*6], priorpars[4+n_source*7+i_lens*6, 0], priorpars[4+n_source*7+i_lens*6, 1])
		cube[5+n_source*7+i_lens*6]=uniform_prior(cube[5+n_source*7+i_lens*6], priorpars[5+n_source*7+i_lens*6, 0], priorpars[5+n_source*7+i_lens*6, 1])
		if cube[4+n_source*7+i_lens*6] > 1.0:
			cube[4+n_source*7+i_lens*6]=1.0/cube[4+n_source*7+i_lens*6]
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+90.0

		if cube[3+n_source*7+i_lens*6] < 0.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+180.0

		if cube[3+n_source*7+i_lens*6] > 180.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]-180.0		

	if has_shear:
		cube[n_source*7+n_lens*6]=gauss_prior(cube[n_source*7+n_lens*6], priorpars[n_source*7+n_lens*6, 0], priorpars[n_source*7+n_lens*6, 1])
		cube[n_source*7+n_lens*6+1]=uniform_prior(cube[n_source*7+n_lens*6+1], priorpars[n_source*7+n_lens*6+1, 0], priorpars[n_source*7+n_lens*6+1, 1])
		
def prior_C(cube, ndim, nparams):
	for i_source in np.arange(n_source):
		cube[0+i_source*7]=uniform_prior(cube[0+i_source*7], priorpars[0+i_source*7, 0], priorpars[0+i_source*7, 1])
		cube[1+i_source*7]=uniform_prior(cube[1+i_source*7], priorpars[1+i_source*7, 0], priorpars[1+i_source*7, 1])
		cube[2+i_source*7]=uniform_prior(cube[2+i_source*7], priorpars[2+i_source*7, 0], priorpars[2+i_source*7, 1])
		cube[3+i_source*7]=uniform_prior(cube[3+i_source*7], priorpars[3+i_source*7, 0], priorpars[3+i_source*7, 1])
		cube[4+i_source*7]=uniform_prior(cube[4+i_source*7], priorpars[4+i_source*7, 0], priorpars[4+i_source*7, 1])
		cube[5+i_source*7]=uniform_prior(cube[5+i_source*7], priorpars[5+i_source*7, 0], priorpars[5+i_source*7, 1])
		cube[6+i_source*7]=uniform_prior(cube[6+i_source*7], priorpars[6+i_source*7, 0], priorpars[6+i_source*7, 1])			
	#endfor
	for i_lens in np.arange(n_lens):
		cube[0+n_source*7+i_lens*6]=uniform_prior(cube[0+n_source*7+i_lens*6], priorpars[0+n_source*7+i_lens*6, 0], priorpars[0+n_source*7+i_lens*6, 1])
		cube[1+n_source*7+i_lens*6]=uniform_prior(cube[1+n_source*7+i_lens*6], priorpars[1+n_source*7+i_lens*6, 0], priorpars[1+n_source*7+i_lens*6, 1])
		cube[2+n_source*7+i_lens*6]=uniform_prior(cube[2+n_source*7+i_lens*6], priorpars[2+n_source*7+i_lens*6, 0], priorpars[2+n_source*7+i_lens*6, 1])
		cube[3+n_source*7+i_lens*6]=uniform_prior(cube[3+n_source*7+i_lens*6], priorpars[3+n_source*7+i_lens*6, 0], priorpars[3+n_source*7+i_lens*6, 1])
		cube[4+n_source*7+i_lens*6]=uniform_prior(cube[4+n_source*7+i_lens*6], priorpars[4+n_source*7+i_lens*6, 0], priorpars[4+n_source*7+i_lens*6, 1])
		cube[5+n_source*7+i_lens*6]=uniform_prior(cube[5+n_source*7+i_lens*6], priorpars[5+n_source*7+i_lens*6, 0], priorpars[5+n_source*7+i_lens*6, 1])
		if cube[4+n_source*7+i_lens*6] > 1.0:
			cube[4+n_source*7+i_lens*6]=1.0/cube[4+n_source*7+i_lens*6]
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+90.0

		if cube[3+n_source*7+i_lens*6] < 0.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+180.0

		if cube[3+n_source*7+i_lens*6] > 180.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]-180.0		
		if has_shear:
			cube[n_source*7+n_lens*6]=uniform_prior(cube[n_source*7+n_lens*6], priorpars[n_source*7+n_lens*6, 0], priorpars[n_source*7+n_lens*6, 1])
			cube[n_source*7+n_lens*6+1]=uniform_prior(cube[n_source*7+n_lens*6+1], priorpars[n_source*7+n_lens*6+1, 0], priorpars[n_source*7+n_lens*6+1, 1])

def prior_A_pixelized(cube, ndim, nparams):
	for i_lens in np.arange(n_lens):
		cube[0+i_lens*6]=gauss_prior(cube[0+i_lens*6], priorpars[0+i_lens*6, 0], priorpars[0+i_lens*6, 1])
		cube[1+i_lens*6]=uniform_prior(cube[1+i_lens*6], priorpars[1+i_lens*6, 0], priorpars[1+i_lens*6, 1])
		cube[2+i_lens*6]=uniform_prior(cube[2+i_lens*6], priorpars[2+i_lens*6, 0], priorpars[2+i_lens*6, 1])
		cube[3+i_lens*6]=uniform_prior(cube[3+i_lens*6], priorpars[3+i_lens*6, 0], priorpars[3+i_lens*6, 1])
		cube[4+i_lens*6]=uniform_prior(cube[4+i_lens*6], priorpars[4+i_lens*6, 0], priorpars[4+i_lens*6, 1])
		cube[5+i_lens*6]=uniform_prior(cube[5+i_lens*6], priorpars[5+i_lens*6, 0], priorpars[5+i_lens*6, 1])
		if cube[4+i_lens*6] > 1.0:
			cube[4+i_lens*6]=1.0/cube[4+i_lens*6]
			cube[3+i_lens*6]=cube[3+i_lens*6]+90.0

		if cube[3+i_lens*6] < 0.0:
			cube[3+i_lens*6]=cube[3+i_lens*6]+180.0

		if cube[3+i_lens*6] > 180.0:
			cube[3+i_lens*6]=cube[3+i_lens*6]-180.0		

	if has_shear:
		cube[n_lens*6]=gauss_prior(cube[n_lens*6], priorpars[n_lens*6, 0], priorpars[n_lens*6, 1])
		cube[n_lens*6+1]=uniform_prior(cube[n_lens*6+1], priorpars[n_lens*6+1, 0], priorpars[n_lens*6+1, 1])
		
	cube[-1]=uniform_prior(cube[-1], priorpars[-1, 0], priorpars[-1, 1])

def loglike(cube, ndim, nparams):
	p=np.zeros(ndim)
	for i in range(ndim):
		p[i]=cube[i]
	for i_lens in range(n_lens):
		if (p[0+n_source*7+i_lens*6] < 0 or p[4+n_source*7+i_lens*6] < 0):
			return np.log(0.0)
	I_model=model(x1, y1, p, lens_model, source_model, PSF=PSF, n_source=n_source, n_lens=n_lens, has_shear=has_shear)
	chi2=np.sum((I_image1-I_model)**2.*I_invvar1)
	return -0.5*chi2

def loglike_1011(cube, ndim, nparams):
	p=np.zeros(ndim)
	for i in range(ndim):
		p[i]=cube[i]
	for i_source in range(n_source):
		p[3+i_source*7]=10.**(p[3+i_source*7])
		p[5+i_source*7]=10.**(p[5+i_source*7])
		p[6+i_source*7]=10.**(p[6+i_source*7])
	for i_lens in range(n_lens):
		if (p[0+n_source*7+i_lens*6] < 0 or p[4+n_source*7+i_lens*6] < 0):
			return np.log(0.0)
	I_model=model(x1, y1, p[n_source*7:], p[0:n_source*7], lens_model, source_model, PSF=PSF, n_source=n_source, n_lens=n_lens)
	chi2=np.sum((I_image1-I_model)**2.*I_invvar1)
	return -0.5*chi2

def prior_pixelized(cube, ndim, nparams):
	for i in range(ndim):
		cube[i]=uniform_prior(cube[i], priorpars[i, 0], priorpars[i, 1])

def loglike_pixelized(cube, ndim, nparams):
	p=np.zeros(ndim)
	for i in range(ndim):
		p[i]=cube[i]
	lpar=p[0:-1]
	lam=p[-1]
	output=pix_source(I_image1, I_invvar1, lens_mask_2d, used_pixels, used_pixels_mask, x1, y1, lpar, lens_model, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=0)
	I_fit_pix=output['fit']
	Vec_S=output['solution']	
	chi2=np.sum((I_image1-I_fit_pix)**2.*I_invvar1)
	if reg_scheme == 0:
		chi2+=lam*(Vec_S**2.).sum()
	return -0.5*chi2

def prior_phot(cube, ndim, nparams):
	for i_source in np.arange(n_source):
		cube[0+i_source*7]=uniform_prior(cube[0+i_source*7], priorpars[0+i_source*7, 0], priorpars[0+i_source*7, 1])
		cube[1+i_source*7]=uniform_prior(cube[1+i_source*7], priorpars[1+i_source*7, 0], priorpars[1+i_source*7, 1])
		cube[2+i_source*7]=uniform_prior(cube[2+i_source*7], priorpars[2+i_source*7, 0], priorpars[2+i_source*7, 1])
		cube[3+i_source*7]=uniform_prior(cube[3+i_source*7], priorpars[3+i_source*7, 0], priorpars[3+i_source*7, 1])
		cube[4+i_source*7]=uniform_prior(cube[4+i_source*7], priorpars[4+i_source*7, 0], priorpars[4+i_source*7, 1])
		cube[5+i_source*7]=uniform_prior(cube[5+i_source*7], priorpars[5+i_source*7, 0], priorpars[5+i_source*7, 1])
		cube[6+i_source*7]=uniform_prior(cube[6+i_source*7], priorpars[6+i_source*7, 0], priorpars[6+i_source*7, 1])
	#endfor
	for i_lens in np.arange(n_lens):
		cube[0+n_source*7+i_lens*6]=uniform_prior(cube[0+n_source*7+i_lens*6], priorpars[0+n_source*7+i_lens*6, 0], priorpars[0+n_source*7+i_lens*6, 1])
		cube[1+n_source*7+i_lens*6]=uniform_prior(cube[1+n_source*7+i_lens*6], priorpars[1+n_source*7+i_lens*6, 0], priorpars[1+n_source*7+i_lens*6, 1])
		cube[2+n_source*7+i_lens*6]=uniform_prior(cube[2+n_source*7+i_lens*6], priorpars[2+n_source*7+i_lens*6, 0], priorpars[2+n_source*7+i_lens*6, 1])
		cube[3+n_source*7+i_lens*6]=uniform_prior(cube[3+n_source*7+i_lens*6], priorpars[3+n_source*7+i_lens*6, 0], priorpars[3+n_source*7+i_lens*6, 1])
		cube[4+n_source*7+i_lens*6]=uniform_prior(cube[4+n_source*7+i_lens*6], priorpars[4+n_source*7+i_lens*6, 0], priorpars[4+n_source*7+i_lens*6, 1])
		cube[5+n_source*7+i_lens*6]=uniform_prior(cube[5+n_source*7+i_lens*6], priorpars[5+n_source*7+i_lens*6, 0], priorpars[5+n_source*7+i_lens*6, 1])
		if cube[4+n_source*7+i_lens*6] > 1.0:
			cube[4+n_source*7+i_lens*6]=1.0/cube[4+n_source*7+i_lens*6]
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+90.0

		if cube[3+n_source*7+i_lens*6] < 0.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]+180.0

		if cube[3+n_source*7+i_lens*6] > 180.0:
			cube[3+n_source*7+i_lens*6]=cube[3+n_source*7+i_lens*6]-180.0		

	if has_shear:
		cube[n_source*7+n_lens*6]=uniform_prior(cube[n_source*7+n_lens*6], priorpars[n_source*7+n_lens*6, 0], priorpars[n_source*7+n_lens*6, 1])
		cube[n_source*7+n_lens*6+1]=uniform_prior(cube[n_source*7+n_lens*6+1], priorpars[n_source*7+n_lens*6+1, 0], priorpars[n_source*7+n_lens*6+1, 1])

	for i_lens in np.arange(n_lens):
		cube[0+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[0+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[0+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[0+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[1+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[1+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[1+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[1+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[2+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[2+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[2+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[2+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[3+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[3+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[3+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[3+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[4+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[4+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[4+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[4+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[5+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[5+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[5+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[5+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])
		cube[6+n_source*7+n_lens*6+has_shear*2+i_lens*7]=uniform_prior(cube[6+n_source*7+n_lens*6+has_shear*2+i_lens*7], priorpars[6+n_source*7+n_lens*6+has_shear*2+i_lens*7, 0], priorpars[6+n_source*7+n_lens*6+has_shear*2+i_lens*7, 1])

def loglike_phot(cube, ndim, nparams):
	p=np.zeros(ndim)
	for i in range(ndim):
		p[i]=cube[i]
	for i_lens in range(n_lens):
		if (p[0+n_source*7+i_lens*6] < 0 or p[4+n_source*7+i_lens*6] < 0):
			return np.log(0.0)
	p_lens=p[0:-7*n_lens]
	p_phot=p[-7*n_lens:]
	I_model=model(x1, y1, p_lens, lens_model, source_model, PSF=PSF, n_source=n_source, n_lens=n_lens, has_shear=has_shear)
	I_lens=0.0*x1
	for i_lens in np.arange(n_lens_phot):
		p_phot_i=p_phot[i_lens*7:i_lens*7+7]
		I_lens+=sersic_phot(x1, y1, p_phot_i)
	I_lens=signal.fftconvolve(I_lens, PSF, mode='same')
	chi2=np.sum((I_data1-I_model-I_lens)**2.*I_invvar1)
	return -0.5*chi2

def ln_gaussian_prior(x, par1, par2):
	return -0.5*((x-par1)/par2)**2.-np.log(np.sqrt(2.*np.pi)*par2)

def ln_uniform_prior(x, par1, par2):
	if par1 == par2:
		return 0.0
	if x < par1 or x > par2:
#		print 'parameter x out of bound'
#		print x, par1, par2
		return -np.inf
	else:
		return -1.0*np.log(par2-par1)

def lnprior(theta, priorpars):
	lnprior=0.
	for i in np.arange(len(theta)):
		if priorpars[i, 0] == 0: # uniform prior
			lnprior+=ln_uniform_prior(theta[i], priorpars[i, 1], priorpars[i, 2])
		elif priorpars[i, 0] == 1: # gaussian prior
			lnprior+=ln_gaussian_prior(theta[i], priorpars[i, 1], priorpars[i, 2])
	return lnprior

def lnlike_phot(theta, x, y, PSF, n_phot, photmodel, I_data, I_invvar):
	I_phot=phot_model(x, y, theta, photmodel, PSF=PSF, n_phot=n_phot)
	chi2=np.sum((I_data-I_phot)**2.*I_invvar)
	return -0.5*chi2

def lnprob_phot(theta, x, y, PSF, n_phot, photmodel, I_data, I_invvar, priorpars):
	lp=lnprior(theta, priorpars)
	if not np.isfinite(lp):
		return -np.inf
	return lp+lnlike_phot(theta, x, y, PSF, n_phot, photmodel, I_data, I_invvar)

def lnlike(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF):
	n_lens_par_dict={'sie':6, 'sple':6, 'pm':3, 'snfw':6, 'ext. shear':2, 'sple_ppn':7}
	n_source_par_dict={'sersic': 7, 'gaussian': 7, 'pix':2}
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]
	spar=theta[0:spar_index]
	lpar=theta[spar_index:spar_index+lpar_index]
	ppar=theta[spar_index+lpar_index:]

	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	I_model=model(x, y, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)
	chi2=np.sum((I_data-I_phot-I_model)**2.0*I_invvar)
	return -0.5*chi2
		
def lnprob(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF, priorpars):
	lp=lnprior(theta, priorpars)
	if not np.isfinite(lp):
		return -np.inf
	return lp+lnlike(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF) 

def f(alpha, delta, beta):
	psi=alpha+delta-2.0
	return 2.0*np.sqrt(np.pi)*(delta-3.0)/(psi-3.0)/(psi-2.0*beta)*(sf.gamma(0.5*(psi-1))/sf.gamma(0.5*psi)-beta*sf.gamma(0.5*(psi+1.0))/sf.gamma(0.5*(psi+2.0)))*(sf.gamma(0.5*delta)*sf.gamma(0.5*alpha)/sf.gamma(0.5*(delta-1.))/sf.gamma(0.5*(alpha-1.)))

def g(alpha, delta, beta):
	psi=alpha+delta-2.0
	return (sf.gamma(0.5*(psi-1))/sf.gamma(0.5*psi)-beta*sf.gamma(0.5*(psi+1.0))/sf.gamma(0.5*(psi+2.0)))*(sf.gamma(0.5*delta)*sf.gamma(0.5*alpha)/sf.gamma(0.5*(delta-1.))/sf.gamma(0.5*(alpha-1.)))*(sf.gamma(0.5*(5.0-delta-alpha))/sf.gamma(0.5*(3.-delta)))/(psi-2.0*beta)

#def lnlike_kinematics(theta_kinematics, sigma_obs, sigma_err, theta_e, D_s, D_ds):
#	alpha, delta, beta, gamma = theta_kinematics
#	alpha+=1.0
#	c=2.997924*10.0**5.
#	theta_a=2.8/2.355
#	sigma_predicted=c/(360.0*np.sqrt(20.0))*2.0**(0.5-0.25*alpha)*np.sqrt((5.0-delta-alpha)/(3.0-delta)*(theta_e*D_s/D_ds*(2./(1.+gamma))*f(alpha, delta, beta)*(theta_a/theta_e)**(2.0-alpha))*(sf.gamma(0.5*(5.0-delta-alpha))/sf.gamma(0.5*(3.-delta))))
#	if not np.isfinite(sigma_predicted):
#		return -np.inf
#	return -0.5*((sigma_obs-sigma_predicted)/sigma_err)**2.0

def lnlike_kinematics_ppn(theta_kinematics, sigma_obs, sigma_err, D_s, D_ds):
	alpha, delta, beta, theta_e = theta_kinematics
	alpha+=1.0
	c=2.997924*10.0**5.
	theta_a=2.8/2.355
	sigma_predicted=(c/(360.0*np.sqrt(10.0)))*(np.pi)**0.25*2.0**(0.5-0.25*alpha)*np.sqrt(D_s/D_ds)*np.sqrt(theta_e)*(theta_a/theta_e)**(1-0.5*alpha)*np.sqrt(g(alpha, delta, beta))
	if not np.isfinite(sigma_predicted):
		return -np.inf
	return -0.5*((sigma_obs-sigma_predicted)/sigma_err)**2.0

def lnprob_kinematics(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF, sigma_obs, sigma_err, theta_e, D_s, D_ds, priorpars):
	lp=lnprior(theta, priorpars)
	if not np.isfinite(lp):
		return -np.inf
	theta_lens=theta[0:-3]
	theta_kinematics=theta[-4:]
	lnlike_kine=lnlike_kinematics(theta_kinematics, sigma_obs, sigma_err, theta_e, D_s, D_ds)
	if not np.isfinite(lnlike_kine):
		return -np.inf
	return lp+lnlike(theta_lens, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF)+lnlike_kine

def lnprob_kinematics_ppn(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF, sigma_obs, sigma_err, D_s, D_ds, priorpars):
	lp=lnprior(theta, priorpars)
	if not np.isfinite(lp):
		return -np.inf
	theta_lens=theta[0:-2]
	theta_kinematics=np.zeros(4)
	theta_kinematics[0]=theta[-4]
	theta_kinematics[1]=theta[-2]
	theta_kinematics[2]=theta[-1]
	theta_kinematics[3]=theta[-9]
	lnlike_kine=lnlike_kinematics_ppn(theta_kinematics, sigma_obs, sigma_err, D_s, D_ds)
	if not np.isfinite(lnlike_kine):
		return -np.inf
	return lp+lnlike(theta_lens, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, PSF)+lnlike_kine

def lnlike_pix(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, PSF, PSF_csr):
	n_lens_par_dict={'sie':6, 'sple':6, 'pm':3, 'snfw':6, 'ext. shear':2}
	n_source_par_dict={'sersic': 7, 'gaussian': 7, 'pix':2}
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]
	spar=theta[0:spar_index]
	lpar=theta[spar_index:spar_index+lpar_index]
	ppar=theta[spar_index+lpar_index:]
	lam=spar[0]
	reg_scheme=spar[1]

	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	output=pix_source(I_data-I_phot, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, x, y, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=0)
	I_model=output['fit']
	Vec_S=output['solution']

	if reg_scheme == 0:
		G_L=np.sum(Vec_S**2.0)
	elif reg_scheme == 1:
		theta_x=x[lens_mask_2d].flatten()
		theta_y=y[lens_mask_2d].flatten()
		(alpha_x, alpha_y)=def_total(x, y, lpar, lensmodel, n_lens=n_lens)
		alpha_x_1d=alpha_x[lens_mask_2d].flatten()
		alpha_y_1d=alpha_y[lens_mask_2d].flatten()
		tri=construct_source_grid(theta_x, theta_y, alpha_x_1d, alpha_y_1d, used_pixels_mask)
		(indices, indptr)=tri.vertex_neighbor_vertices
		G_L=np.array([0])
		for i in range(Vec_S.shape[0]):
			G_L[0]+=((Vec_S[indptr[indices[i]:indices[i+1]]]-Vec_S[i])**2.0).sum()

	chi2=np.sum((I_data-I_phot-I_model)**2.0*I_invvar)+lam*G_L
	return -0.5*chi2
		
def lnprob_pix(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, PSF, PSF_csr, priorpars):
	lp=lnprior(theta, priorpars)
	if not np.isfinite(lp):
		return -np.inf
	return lp+lnlike_pix(theta, x, y, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, PSF, PSF_csr)


