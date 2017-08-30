#!/usr/bin/env python
import numpy as np
import os, sys, commands, pyfits
import lens_funcs_gui as lf
from lmfit import Parameters, minimize, fit_report, conf_interval, printfuncs
import string
from scipy.interpolate import griddata
import emcee
from emcee.utils import MPIPool

class lfit:
	def __init__(self):
		self.data_file=''
		self.data_ext=0
		self.invvar_file=''
		self.invvar_ext=0
		self.psf_file=''
		self.psf_ext=0
		self.jmask_file=''
		self.jmask_ext=0

		self.hw=0
		self.dpix=0.0
		self.n_lens=0
		self.lensmodel=[]
		self.n_source=0
		self.sourcemodel=[]
		self.n_phot=0
		self.photmodel=[]
		self.fit_mode=0

		self.PSF_csr=None

		self.n_lens_par_dict={'sie':6, 'sple':6, 'pm':3, 'snfw':6, 'ext. shear':2}
		self.n_source_par_dict={'sersic': 7, 'gaussian': 7, 'pix':2}
		self.n_phot_par_dict={'sersic': 8, 'csersic': 11, 'hernquist': 7, 'bspline': 0}

	def params_setting(self, var_string):
		setting={'value': 0.0, 'vary': 1, 'min': None, 'max': None, 'expr': None}
		try:
			float(var_string)
		except:
			pass
		else:
			setting['value']=float(var_string)
			return setting
		try:
			low=string.index(var_string, 'value=')
		except: 
			pass
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				setting['value']=float(var_string[low+6:high])
			else:
				setting['value']=float(var_string[low+6:])
		try:
			low=string.index(var_string, 'vary=')
		except: 
			pass
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				setting['vary']=int(var_string[low+5:high])
			else:
				setting['vary']=int(var_string[low+5:])
	
		try:
			low=string.index(var_string, 'min=')
		except: 
			pass
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				setting['min']=float(var_string[low+4:high])
			else:
				setting['min']=float(var_string[low+4:])
	
		try:
			low=string.index(var_string, 'max=')
		except: 
			pass
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				setting['max']=float(var_string[low+4:high])
			else:
				setting['max']=float(var_string[low+4:])
	
		try:
			low=string.index(var_string, 'expr=')
		except: 
			pass
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				setting['expr']=var_string[low+5:high]
			else:
				setting['expr']=var_string[low+5:]
		return setting

	def residual(self, pars, x, y, PSF, n_lens, n_source, lensmodel, sourcemodel, data=None, I_invvar=None):
		vals=pars.valuesdict()
		p=vals.values()
		spar_index=0
		for i_source in np.arange(n_source):
			spar_index+=self.n_source_par_dict[sourcemodel[i_source]]
		lpar_index=0
		for i_lens in np.arange(n_lens):
			lpar_index+=self.n_lens_par_dict[lensmodel[i_lens]]
		
		spar=p[0:spar_index]
		lpar=p[spar_index:spar_index+lpar_index]
		ppar=p[spar_index+lpar_index:]
		I_model=lf.model(x, y, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)
		if (self.n_phot==0 or (self.n_phot==1 and self.photmodel[0] == 'bspline')):
			I_phot=self.I_bspline1.copy()
		else:
			I_phot=lf.phot_model(x, y, ppar, self.photmodel, PSF=PSF, n_phot=self.n_phot)
	
		if data is None:
	 		return (I_model+I_phot).reshape(-1)
		if I_invvar is None:
			return (data-I_model-I_phot).reshape(-1)
		else:
			return ((data-I_model-I_phot)*np.sqrt(I_invvar)).reshape(-1)
		
	# Displaying the residual after every iteration
	# Works ok, but significantly slows down the program
	#	ax2=fig2.add_subplot(111)
	#	im2=ax2.imshow(I_image1-I_model, vmin=0.0, vmax=np.amax((I_image1-I_model)*(1-gmask1)), **myargs)
	#	canvas2.draw()
	
	def residual_pix(self, pars, x, y, n_lens, lensmodel, lens_mask_2d, used_pixels, used_pixels_mask, data=None, I_invvar=None, PSF=np.array([0]), PSF_csr=np.array([0])):
		vals=pars.valuesdict()
		p=vals.values()
		spar_index=0
		for i_source in np.arange(self.n_source):
			spar_index+=self.n_source_par_dict[self.sourcemodel[i_source]]
		lpar_index=0
		for i_lens in np.arange(self.n_lens):
			lpar_index+=self.n_lens_par_dict[lensmodel[i_lens]]
		spar=p[0:spar_index]
		lpar=p[spar_index:spar_index+lpar_index]
		ppar=p[spar_index+lpar_index:]
		lam=spar[0]
		reg_scheme=spar[1]
		if (self.n_phot==0 or (self.n_phot==1 and self.photmodel[0] == 'bspline')):
			I_phot=self.I_bspline1.copy()
		else:
			I_phot=lf.phot_model(x, y, ppar, self.photmodel, PSF=PSF, n_phot=self.n_phot)
	
		output=lf.pix_source(data-I_phot, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, x, y, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=0)
		I_model=output['fit']
		Vec_S=output['solution']
	
		if reg_scheme == 0:
			G_L=Vec_S
		elif reg_scheme == 1:
			theta_x=x[lens_mask_2d].flatten()
			theta_y=y[lens_mask_2d].flatten()
			(alpha_x, alpha_y)=lf.def_total(x, y, lpar, lensmodel, n_lens=n_lens)
			alpha_x_1d=alpha_x[lens_mask_2d].flatten()
			alpha_y_1d=alpha_y[lens_mask_2d].flatten()
			tri=lf.construct_source_grid(theta_x, theta_y, alpha_x_1d, alpha_y_1d, used_pixels_mask)
			(indices, indptr)=tri.vertex_neighbor_vertices
			G_L=np.array([0])
			for i in range(Vec_S.shape[0]):
				G_L[0]+=((Vec_S[indptr[indices[i]:indices[i+1]]]-Vec_S[i])**2.0).sum()
			G_L=np.sqrt(G_L)
		if data is None:
	 		return np.concatenate(((I_model+I_phot).reshape(-1), np.sqrt(lam)*G_L))
		if I_invvar is None:
			return np.concatenate(((data-I_model-I_phot).reshape(-1), np.sqrt(lam)*G_L))
		else:
			return np.concatenate((((data-I_model-I_phot)*np.sqrt(I_invvar)).reshape(-1), np.sqrt(lam)*G_L))
	
	def messg(self, params, iter, resid, *args, **kws):
		print '# of Iteration: %i' %iter
		print 'Chi2= %10.4f' %np.sum(resid**2.0)
		print fit_report(params, show_correl=False)

	def generate_models(self, n_lens, lpar, lensmodel, n_source, spar, sourcemodel, n_phot, ppar, photmodel, PSF, PSF_csr=np.array([0])):
		if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
			I_phot=self.I_bspline1
		else:
			I_phot=lf.phot_model(self.x1, self.y1, ppar, photmodel, PSF=PSF, n_phot=n_phot)

		if (len(sourcemodel) == 1 and sourcemodel[0] == 'pix'):
			lam=spar[0]
			reg_scheme=spar[1]
			output=lf.pix_source(self.I_data1-I_phot, self.I_invvar1, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask, self.x1, self.y1, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=1)
			I_fit=output['fit']
			Vec_S=output['solution']
			Cov_matrix=output['Covariance']
			Vec_S_err=np.sqrt(np.diagonal(Cov_matrix))
		else:
			I_fit=lf.model(self.x1, self.y1, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)

		n_x_source=8*100+1
		n_y_source=8*100+1
		dpix_source=0.01
		self.dpix_source=dpix_source
		source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2)*dpix_source
		source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2, np.ones(n_x_source))*dpix_source
		if len(sourcemodel) ==1 and sourcemodel[0] == 'pix':
			theta_x=self.x1[self.lens_mask_1d].flatten()
			theta_y=self.y1[self.lens_mask_1d].flatten()
			(alpha_x, alpha_y)=lf.def_total(self.x1, self.y1, lpar, lensmodel, n_lens=n_lens)
			alpha_x_1d=alpha_x[self.lens_mask_1d].flatten()
			alpha_y_1d=alpha_y[self.lens_mask_1d].flatten()
			
			tri=lf.construct_source_grid(theta_x, theta_y, alpha_x_1d, alpha_y_1d, self.used_pixels_mask)
			I_source=griddata(tri.points, Vec_S, (source_xbase, source_ybase), method='linear', fill_value=0.0)
			I_source_err=griddata(tri.points, Vec_S_err, (source_xbase, source_ybase), method='linear', fill_value=0.0)
		else:
			I_source=0.0*source_xbase
			current_source_index=0
			for i_source in np.arange(n_source):
				spar_i=spar[current_source_index:current_source_index+self.n_source_par_dict[self.sourcemodel[i_source]]]
				I_source+=lf.sersic_2d(source_xbase, source_ybase, spar_i)
				current_source_index+=self.n_source_par_dict[sourcemodel[i_source]]

		self.I_fit=I_fit
		self.I_phot=I_phot
		self.I_source=I_source
		if len(self.sourcemodel) ==1 and self.sourcemodel[0] == 'pix':
			self.I_source_err=I_source_err
	
	def save_fit(self):
		fname=(self.output_fname).strip()
		col1=pyfits.Column(name='spar', format='D', array=self.spar)
		col1_err=pyfits.Column(name='spar_err', format='D', array=self.spar_err)
		col2=pyfits.Column(name='lpar', format='D', array=self.lpar)
		col2_err=pyfits.Column(name='lpar_err', format='D', array=self.lpar_err)
		col3=pyfits.Column(name='p_phot', format='D', array=self.ppar)
		col3_err=pyfits.Column(name='p_phot_err', format='D', array=self.ppar_err)
		tbhdu1=pyfits.BinTableHDU.from_columns([col1, col1_err])
		tbhdu2=pyfits.BinTableHDU.from_columns([col2, col2_err])
		tbhdu3=pyfits.BinTableHDU.from_columns([col3, col3_err])
		
		if (self.fit_mode == 1):
			hdu=pyfits.PrimaryHDU()
		elif self.fit_mode == 2:
			hdu=pyfits.PrimaryHDU(self.mcmc_out)
		hdu.header['target']=self.obj_name
		hdu.header['N_Lens']=self.n_lens
		for i_lens in np.arange(self.n_lens):
			hdu.header['lmodel_'+str(i_lens+1)]=self.lensmodel[i_lens]
		hdu.header['N_Source']=self.n_source
		for i_source in np.arange(self.n_source):
			hdu.header['smodel_'+str(i_source+1)]=self.sourcemodel[i_source]
		hdu.header['N_Phot']=self.n_phot
		for i_phot in np.arange(self.n_phot):
			hdu.header['pmodel_'+str(i_phot+1)]=self.photmodel[i_phot]
		try:
			hdu.header['DOF']=self.dof
		except:
			hdu.header['DOF']=self.out.nfree
		try:
			hdu.header['Chi2']=self.chi2
		except:
			hdu.header['Chi2']=self.out.chisqr
		hdu.header['dpix_l']=self.dpix
		try:
			hdu.header['dpix_s']=self.dpix_source
		except:
			hdu.header['dpix_s']=0.01
		hdu.header.add_comment('1st extension is source parameters')
		hdu.header.add_comment('2nd extension is lens model parameters')
		hdu.header.add_comment('3rd extension is photometric parameters')
		hdu.header.add_comment('4th extension is the HST image')
		hdu.header.add_comment('5th extension is the inverse variance')
		hdu.header.add_comment('6th extension is the photometric fit')
		hdu.header.add_comment('7th extension is the lens model')
		hdu.header.add_comment('8th extension is the source model')
		hdulist=pyfits.HDUList([hdu, tbhdu1, tbhdu2, tbhdu3])
		hdulist.writeto(fname, clobber=True)
		pyfits.append(fname, self.I_data1)
		pyfits.append(fname, self.I_invvar1)
		pyfits.append(fname, self.I_phot)
		pyfits.append(fname, self.I_fit)
		pyfits.append(fname, self.I_source)
		if (len(self.sourcemodel) == 1 and self.sourcemodel[0] == 'pix'):
			hdu.header.add_comment('9th extension is the source error')
		 	pyfits.append(fname, self.I_source_err)
		print 'Saved to the file '+fname

	def generate_bizfile(self, parfile):
		# read in the par file
		setup_file=open(parfile, 'r')
		data_file, data_ext=(setup_file.readline()).split()
		data_ext=int(data_ext)
		invvar_file, invvar_ext=(setup_file.readline()).split()
		invvar_ext=int(invvar_ext)
		psf_file, psf_ext=(setup_file.readline()).split()
		psf_ext=int(psf_ext)
		jmask_file, jmask_ext=(setup_file.readline()).split()
		jmask_ext=int(jmask_ext)
		gmask_file, gmask_ext=(setup_file.readline()).split()
		gmask_ext=int(gmask_ext)
		phot_file, phot_ext=(setup_file.readline()).split()
		phot_ext=int(phot_ext)
		xcen, ycen, q, pa=(setup_file.readline()).split()
		xcen=float(xcen)
		ycen=float(ycen)
		q=float(q)
		pa=float(pa)
		cycle, program, visit=(setup_file.readline()).split()
		cycle=cycle.strip()
		program=program.strip()
		visit=visit.strip()
		output_fname=setup_file.readline()
	
		setup_file.close()

		if not os.path.exists(os.getenv('HST_DATAROOT')+'/'+cycle):
			os.makedirs(os.getenv('HST_DATAROOT')+'/'+cycle)
		if not os.path.exists(os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program):
			os.makedirs(os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program)
		if not os.path.exists(os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program+'/'+visit):
			os.makedirs(os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program+'/'+visit)		
		hdulist=pyfits.open(data_file)
		I_data=hdulist[data_ext].data
		hdulist.close()
		hdulist=pyfits.open(invvar_file)
		I_invvar=hdulist[invvar_ext].data
		hdulist.close()
		hdulist=pyfits.open(psf_file)
		PSF=hdulist[psf_ext].data
		PSF=PSF.astype(float) # the original PSF is float32, while the default is float64 in python
		PSF=PSF/PSF.sum()
		hdulist.close()
		hdulist=pyfits.open(jmask_file)
		jmask=hdulist[jmask_ext].data
		hdulist.close()
		hdulist=pyfits.open(gmask_file)
		gmask=hdulist[gmask_ext].data
		hdulist.close()
		hdulist=pyfits.open(phot_file)
		I_bspline=hdulist[phot_ext].data
		hdulist.close()
			
		phw=18
		nx_psf=PSF.shape[0]
		if phw > nx_psf/2: 
			phw=nx_psf/2
		wmax=np.argmax(PSF)
		pxc=np.mod(wmax, nx_psf)
		pyc=wmax/nx_psf
		tpsf=PSF[pxc-phw:pxc+phw+1,pyc-phw:pyc+phw+1]
		tpsf=tpsf/tpsf.sum()

		fname=(output_fname).strip()
		hdu=pyfits.PrimaryHDU(I_data)
		hdulist=pyfits.HDUList([hdu])
		hdulist.append(pyfits.ImageHDU(I_invvar))
		hdulist.append(pyfits.ImageHDU(0*I_data)) # masks for negative pixels and cosmic rays in the true biz file
		hdulist.append(pyfits.ImageHDU(PSF))
		hdulist.append(pyfits.ImageHDU(gmask))
		hdulist.append(pyfits.ImageHDU(jmask))
		hdulist.append(pyfits.ImageHDU(gmask))
		header=pyfits.Header()
		header['xcenter']=xcen
		header['ycenter']=ycen
		header['axisrat']=q
		header['axisangl']=pa
		hdulist.append(pyfits.ImageHDU(data=I_bspline, header=header))
		hdulist.append(pyfits.ImageHDU(I_bspline))
		hdulist.append(pyfits.ImageHDU(tpsf))
		hdulist.writeto(os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program+'/'+visit+'/'+fname, clobber=True)
		print 'Saved to '+os.getenv('HST_DATAROOT')+'/'+cycle+'/'+program+'/'+visit+'/'+fname

	def setup_lfit(self, parfile):
		# read in the par file
		setup_file=open(parfile, 'r')
		self.setup_file=setup_file
		data_file, data_ext=(setup_file.readline()).split()
		data_ext=int(data_ext)
		self.obj_name=data_file[data_file.find('SLACSJ'):data_file.find('SLACSJ')+14+1]
		invvar_file, invvar_ext=(setup_file.readline()).split()
		invvar_ext=int(invvar_ext)
		psf_file, psf_ext=(setup_file.readline()).split()
		psf_ext=int(psf_ext)
		jmask_file, jmask_ext=(setup_file.readline()).split()
		jmask_ext=int(jmask_ext)
		gmask_file, gmask_ext=(setup_file.readline()).split()
		gmask_ext=int(gmask_ext)
		phot_file, phot_ext=(setup_file.readline()).split()
		phot_ext=int(phot_ext)

		self.output_fname=setup_file.readline()

		hw, dpix=(setup_file.readline()).split()
		hw=int(hw)
		dpix=float(dpix)
		self.dpix=dpix

		n_source, sourcemodel_temp=(setup_file.readline()).split()
		n_source=int(n_source)
		self.n_source=n_source
		sourcemodel=[model for model in sourcemodel_temp.split(',')]
		self.sourcemodel=sourcemodel
		n_lens, lensmodel_temp=(setup_file.readline()).split()
		n_lens=int(n_lens)
		self.n_lens=n_lens
		lensmodel=[model for model in lensmodel_temp.split(',')]
		self.lensmodel=lensmodel
		n_phot, photmodel_temp=(setup_file.readline()).split()
		n_phot=int(n_phot)
		self.n_phot=n_phot
		photmodel=[model for model in photmodel_temp.split(',')]
		self.photmodel=photmodel
		
		self.fit_mode=int(setup_file.readline())
		
		hdulist=pyfits.open(data_file)
		I_data=hdulist[data_ext].data
		hdulist.close()
		hdulist=pyfits.open(invvar_file)
		I_invvar=hdulist[invvar_ext].data
		hdulist.close()
		hdulist=pyfits.open(psf_file)
		PSF=hdulist[psf_ext].data
		PSF=PSF.astype(float) # the original PSF is float32, while the default is float64 in python
		self.PSF=PSF
		hdulist.close()
		hdulist=pyfits.open(jmask_file)
		jmask=hdulist[jmask_ext].data
		self.jmask=jmask
		hdulist.close()
		hdulist=pyfits.open(gmask_file)
		gmask=hdulist[gmask_ext].data
		self.gmask=gmask
		hdulist.close()
		hdulist=pyfits.open(phot_file)
		I_bspline=hdulist[phot_ext].data
		self.I_bspline=I_bspline
		hdulist.close()
		
		I_invvar=I_invvar*jmask
		I_image=I_data-I_bspline
		n_x, n_y=I_data.shape
		x=np.outer(np.ones(n_y), np.arange(n_x)-n_x/2)*dpix
		y=np.outer(np.arange(n_y)-n_y/2, np.ones(n_x))*dpix
		
		scale_flag=1
		if scale_flag:
			width_shell=15
			wh=np.where(((x <= width_shell*dpix-(n_x/2)*dpix) +(x >= (n_x/2)*dpix-width_shell*dpix)+(y <= width_shell*dpix-(n_y/2)*dpix)+(y >= (n_y/2)*dpix-width_shell*dpix))*(jmask > 0))
			scale=1.0/(np.mean(np.sqrt(I_invvar[wh]))*np.std(I_data[wh]*jmask[wh]))
			I_invvar=I_invvar*scale**2.

		phw=18
		nx_psf=PSF.shape[0]
		if phw > nx_psf/2: 
			phw=nx_psf/2
		wmax=np.argmax(PSF)
		pxc=np.mod(wmax, nx_psf)
		pyc=wmax/nx_psf
		tpsf=PSF[pxc-phw:pxc+phw+1,pyc-phw:pyc+phw+1]
		tpsf=tpsf/tpsf.sum()
		self.tpsf=tpsf

		x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		I_image1=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		I_bspline1=I_bspline[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		jmask1=jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
		gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]

		self.x1=x1
		self.y1=y1
		self.I_data1=I_data1
		self.I_image1=I_image1
		self.I_bspline1=I_bspline1
		self.I_invvar1=I_invvar1
		self.jmask1=jmask1
		self.gmask1=gmask1

		if len(sourcemodel) ==1 and sourcemodel[0] == 'pix':
			lens_mask_2d=1-gmask1
			lens_mask_1d=np.where(lens_mask_2d == 1)
			(used_pixels, used_pixels_mask)=lf.select_used_pixels(I_image1, lens_mask_1d, lens_mask_1d[0].shape[0]/2, scheme=0)
			used_pixels_2d=np.zeros_like(x1, dtype=int)
			used_pixels_2d[used_pixels]=1
			PSF_csr=lf.generate_psf_csr(lens_mask_1d, x1.shape[0], x1.shape[1], tpsf)
			self.lens_mask_1d=lens_mask_1d
			self.lens_mask_2d=lens_mask_2d
			self.used_pixels=used_pixels
			self.used_pixels_mask=used_pixels_mask
			self.used_pixels_2d=used_pixels_2d
			self.PSF_csr=PSF_csr

	def run_lmfit(self):

		leastsq_kws={'maxfev':8000}

		fit_params=Parameters()
		for i_source in np.arange(self.n_source):
			if self.sourcemodel[i_source] == 'gaussian' or self.sourcemodel[i_source] == 'sersic':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('source'+str(i_source)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif self.sourcemodel[i_source] == 'pix':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lam', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('reg_scheme', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		for i_lens in np.arange(self.n_lens):
			if self.lensmodel[i_lens]=='ext. shear':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('shear_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('shear_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif (self.lensmodel[i_lens]=='sie' or self.lensmodel[i_lens]=='sple' or self.lensmodel[i_lens]=='softie'):
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_bsie', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_gamma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif self.lensmodel[i_lens]=='pm':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_theta_e', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif (self.lensmodel[i_lens]=='snfw'):
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_kappa_s', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('lens'+str(i_lens)+'_r_s', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		for i_phot in np.arange(self.n_phot):
			if self.photmodel[i_phot] == 'sersic':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])	
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif self.photmodel[i_phot] == 'csersic':
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_rc', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_alpha', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_gamma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])	
			elif self.photmodel[i_phot] == 'hernquist':		
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_rs', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				setting=self.params_setting(self.setup_file.readline())
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif self.photmodel[i_phot] == 'bspline':
				pass

		self.setup_file.close()

		if (len(self.sourcemodel)==1 and self.sourcemodel[0] == 'pix'):
			out=minimize(self.residual_pix, fit_params, args=(self.x1, self.y1, self.n_lens, self.lensmodel, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask), kws={'data':self.I_data1, 'I_invvar':self.I_invvar1, 'PSF':self.PSF, 'PSF_csr':self.PSF_csr}, iter_cb=self.messg, **leastsq_kws)
		else:	
			out=minimize(self.residual, fit_params, args=(self.x1, self.y1, self.PSF, self.n_lens, self.n_source, self.lensmodel, self.sourcemodel), kws={'data':self.I_data1, 'I_invvar':self.I_invvar1}, iter_cb=self.messg, **leastsq_kws)
		dof=out.nfree-(np.where(self.I_invvar1==0.0))[0].size
		self.out=out
		print "==============================================================="
		print "==============================================================="
		print "Fitting finished"
		print (fit_report(fit_params))
#		ci, trace=conf_interval(out, p_names=['lens0_bsie', 'lens0_xcen', 'lens0_ycen', 'lens1_bsie', 'lens1_xcen', 'lens1_ycen'], sigmas=[0.68,0.95], trace=True, verbose=False)
#		printfuncs.report_ci(ci)
		print "==============================================================="

		fit_par=fit_params.valuesdict()
		pars=fit_par.values()
		spar_index=0
		for i_source in np.arange(self.n_source):
			spar_index+=self.n_source_par_dict[self.sourcemodel[i_source]]
		lpar_index=0
		for i_lens in np.arange(self.n_lens):
			lpar_index+=self.n_lens_par_dict[self.lensmodel[i_lens]]

		spar=pars[0:spar_index]
		lpar=pars[spar_index:spar_index+lpar_index]
		ppar=pars[spar_index+lpar_index:]

		spar_err=np.zeros_like(spar)
		current_source_index=0
		for i_source in np.arange(self.n_source):
			if self.sourcemodel[i_source] == 'sersic' or self.sourcemodel[i_source] == 'gaussian':
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_amp'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_xcen'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_ycen'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_sigma'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_pa'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_q'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['source'+str(i_source)+'_n'].stderr
				current_source_index+=1
			elif self.sourcemodel[i_source] == 'pix':
				spar_err[current_source_index]=fit_params['lam'].stderr
				current_source_index+=1
				spar_err[current_source_index]=fit_params['reg_scheme'].stderr
				current_source_index+=1

		lpar_err=np.zeros_like(lpar)
		current_lens_index=0
		for i_lens in np.arange(self.n_lens):
			if self.lensmodel[i_lens] == 'sie' or self.lensmodel[i_lens] == 'sple' or self.lensmodel[i_lens] == 'softie':
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_bsie'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_xcen'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_ycen'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_pa'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_q'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_gamma'].stderr
				current_lens_index+=1
			elif self.lensmodel[i_lens] == 'pm':
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_theta_e'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_xcen'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_ycen'].stderr
				current_lens_index+=1
			elif self.lensmodel[i_lens] == 'ext. shear':
				lpar_err[current_lens_index]=fit_params['shear_amp'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['shear_pa'].stderr
				current_lens_index+=1	
			elif self.lensmodel[i_lens] == 'snfw':
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_kappa_s'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_xcen'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_ycen'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_pa'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_q'].stderr
				current_lens_index+=1
				lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_r_s'].stderr
				current_lens_index+=1

		ppar_err=np.zeros_like(ppar)
		current_phot_index=0
		for i_phot in np.arange(self.n_phot):
			if self.photmodel[i_phot] == 'sersic':
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_amp'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_xcen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_ycen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_sigma'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_pa'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_q'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_n'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_I0'].stderr
				current_phot_index+=1
			elif self.photmodel[i_phot] == 'csersic':
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_amp'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_xcen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_ycen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_sigma'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_pa'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_q'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_n'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_rc'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_alpha'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_gamma'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_I0'].stderr
				current_phot_index+=1
			elif self.photmodel[i_phot] == 'hernquist':
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_amp'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_xcen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_ycen'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_rs'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_pa'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_q'].stderr
				current_phot_index+=1
				ppar_err[current_phot_index]=fit_params['phot'+str(i_phot)+'_I0'].stderr
				current_phot_index+=1
			elif self.photmodel[i_phot] == 'bspline':
				ppar_err=[]

		self.lpar=lpar
		self.lpar_err=lpar_err
		self.spar=spar
		self.spar_err=spar_err
		self.ppar=ppar
		self.ppar_err=ppar_err

		self.generate_models(self.n_lens, lpar, self.lensmodel, self.n_source, spar, self.sourcemodel, self.n_phot, ppar, self.photmodel, self.PSF, self.PSF_csr)

	def get_source_priors(self, n_source, sourcemodel):
		n_par_source=0
		for i_source in range(self.n_source):
			n_par_source+=self.n_source_par_dict[self.sourcemodel[i_source]]
		scr_priorpars=np.zeros((n_par_source, 3))
		for i_par in range(n_par_source):
			ptype, par1, par2=(self.setup_file.readline()).split()
			scr_priorpars[i_par, 0]=int(ptype)
			scr_priorpars[i_par, 1]=float(par1)
			scr_priorpars[i_par, 2]=float(par2)
		return scr_priorpars

	def get_lens_priors(self, n_lens, lensmodel):
		n_par_lens=0
		for i_lens in range(self.n_lens):
			n_par_lens+=self.n_lens_par_dict[self.lensmodel[i_lens]]
		lens_priorpars=np.zeros((n_par_lens, 3))
		for i_par in range(n_par_lens):
			ptype, par1, par2=(self.setup_file.readline()).split()
			lens_priorpars[i_par, 0]=int(ptype)
			lens_priorpars[i_par, 1]=float(par1)
			lens_priorpars[i_par, 2]=float(par2)
		return lens_priorpars

	def get_phot_priors(self, n_phot, photmodel):
		n_par_phot=0
		for i_phot in range(self.n_phot):
			n_par_phot+=self.n_phot_par_dict[self.photmodel[i_phot]]
		phot_priorpars=np.zeros((n_par_phot, 3))
		for i_par in range(n_par_phot):
			ptype, par1, par2=(self.setup_file.readline()).split()
			phot_priorpars[i_par, 0]=int(ptype)
			phot_priorpars[i_par, 1]=float(par1)
			phot_priorpars[i_par, 2]=float(par2)
		return phot_priorpars

	def init_pos(self, nwalkers, priorpars):
		ndim=(priorpars.shape)[0]
		pos=[]
		for i in range(nwalkers):
			p_i=[]
			for j in range(ndim):
				if priorpars[j, 0] == 0:
					tmp=(priorpars[j, 2]-priorpars[j, 1])*np.random.random()+priorpars[j, 1]
				elif priorpars[j, 0] == 1:
					tmp=priorpars[j, 1]+priorpars[j, 2]*np.random.randn()
				p_i.append(tmp)
			pos.append(p_i)
		return pos

	def run_mcmc(self):
		src_priorpars=self.get_source_priors(self.n_source, self.sourcemodel)
		lens_priorpars=self.get_lens_priors(self.n_lens, self.lensmodel)
		phot_priorpars=self.get_phot_priors(self.n_phot, self.photmodel)
		priorpars=np.vstack((src_priorpars, lens_priorpars, phot_priorpars))
		ndim=(priorpars.shape)[0]
		nthreads, nwalkers, burnin, post_burnin=(self.setup_file.readline()).split()
		nthreads=int(nthreads)
		nwalkers=int(nwalkers)
		burnin=int(burnin)
		post_burnin=int(post_burnin)
		print 'Setting up the sampler...............'
		lf.I_bspline=self.I_bspline1
		if self.n_source == 1 and self.sourcemodel[0] == 'pix':
			sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob_pix, args=(self.x1, self.y1, self.n_lens, self.lensmodel, self.n_source, self.sourcemodel, self.n_phot, self.photmodel, self.I_data1, self.I_invvar1, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask, self.PSF, self.PSF_csr, priorpars), threads=nthreads) 
		else:
			sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob, args=(self.x1, self.y1, self.n_lens, self.lensmodel, self.n_source, self.sourcemodel, self.n_phot, self.photmodel, self.I_data1, self.I_invvar1, self.PSF, priorpars), threads=nthreads) 
		print 'Setting up the sampler...Completed...'
		
		pos=self.init_pos(nwalkers, priorpars)
#		pos=[np.concatenate((spar, lpar, ppar))+0.01*np.random.randn(ndim) for i in range(nwalkers)]
#		for i in range(nwalkers):
#			pos[i][26]=1.0
#			pos[i][32]=1.0
#			pos[i][40]=0.0
#			pos[i][48]=0.0
#			pos[i][56]=0.0
		print 'Running MCMC.........................'
		# run xxx steps as a burn-in
		pos, prob, state=sampler.run_mcmc(pos, burnin)
		# reset the chain to remove the burn-in samples
		sampler.reset()
		print 'Running MCMC....burn-in completed....'
		# starting from the final position in the burn-in chain, sample for xxx
		sampler.run_mcmc(pos, post_burnin, rstate0=state)
		print 'Running MCMC........Completed........'
		print("Mean acceptance fraction: {0:.6f}".format(np.mean(sampler.acceptance_fraction)))
		print 'Autocorrelation time is ', sampler.acor
	
		samples=sampler.chain
		lnprobability=sampler.lnprobability
		lnprobability=lnprobability.reshape(((lnprobability.shape)[0], (lnprobability.shape)[1], 1))
		mcmc_out=np.concatenate((samples, lnprobability), axis=2)
		self.mcmc_out=mcmc_out

		samples=samples[:, :, :].reshape((-1, ndim))
		res=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))	
		p_mcmc=[]
		p_err_mcmc=[]
		for i in range(ndim):
			p_mcmc.append(res[i][0])
			p_err_mcmc.append(0.5*(res[i][1]+res[i][2]))
		
		spar_index=0
		for i_source in np.arange(self.n_source):
			spar_index+=self.n_source_par_dict[self.sourcemodel[i_source]]
		lpar_index=0
		for i_lens in np.arange(self.n_lens):
			lpar_index+=self.n_lens_par_dict[self.lensmodel[i_lens]]

		spar=p_mcmc[0:spar_index]
		lpar=p_mcmc[spar_index:spar_index+lpar_index]
		ppar=p_mcmc[spar_index+lpar_index:]

		spar_err=p_err_mcmc[0:spar_index]
		lpar_err=p_err_mcmc[spar_index:spar_index+lpar_index]
		ppar_err=p_err_mcmc[spar_index+lpar_index:]

		print 'lpar ', lpar
		print 'lpar_err ', lpar_err
		print 'spar ', spar
		print 'spar_err ', spar_err
		print 'ppar ', ppar
		print 'ppar_err ', ppar_err

		self.lpar=lpar
		self.lpar_err=lpar_err
		self.spar=spar
		self.spar_err=spar_err
		self.ppar=ppar
		self.ppar_err=ppar_err

		self.generate_models(self.n_lens, lpar, self.lensmodel, self.n_source, spar, self.sourcemodel, self.n_phot, ppar, self.photmodel, self.PSF, self.PSF_csr)

		self.samples=samples
		self.dof=self.I_data1.size-(np.where(self.I_invvar1==0.0))[0].size
		if (self.n_source == 1 and self.sourcemodel[0] == 'pix'):
			pars=np.concatenate((spar, lpar, ppar))
			self.chi2=self.residual_pix(pars, self.x1, self.y1, self.n_lens, self.lensmodel, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask, data=self.I_data1, I_invvar=self.I_invvar1, PSF=self.PSF, PSF_csr=self.PSF_csr)
		else:
			self.chi2=((self.I_data1-self.I_fit-self.I_phot)**2.0*self.I_invvar1).sum()

	def run_mcmc_mpi(self):
		pool=MPIPool(loadbalance=True)
		if not pool.is_master():
			pool.wait()
			sys.exit(0)

		src_priorpars=self.get_source_priors(self.n_source, self.sourcemodel)
		lens_priorpars=self.get_lens_priors(self.n_lens, self.lensmodel)
		phot_priorpars=self.get_phot_priors(self.n_phot, self.photmodel)
		priorpars=np.vstack((src_priorpars, lens_priorpars, phot_priorpars))
		ndim=(priorpars.shape)[0]
		nwalkers, burnin, post_burnin=(self.setup_file.readline()).split()
		nwalkers=int(nwalkers)
		burnin=int(burnin)
		post_burnin=int(post_burnin)
		print 'Setting up the sampler...............'
		lf.I_bspline=self.I_bspline1
		if self.n_source == 1 and self.sourcemodel[0] == 'pix':
			sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob_pix, args=(self.x1, self.y1, self.n_lens, self.lensmodel, self.n_source, self.sourcemodel, self.n_phot, self.photmodel, self.I_data1, self.I_invvar1, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask, self.PSF, self.PSF_csr, priorpars), pool=pool) 
		else:
			sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob, args=(self.x1, self.y1, self.n_lens, self.lensmodel, self.n_source, self.sourcemodel, self.n_phot, self.photmodel, self.I_data1, self.I_invvar1, self.PSF, priorpars), pool=pool) 
		print 'Setting up the sampler...Completed...'
		
		pos=self.init_pos(nwalkers, priorpars)
#		pos=[np.concatenate((spar, lpar, ppar))+0.01*np.random.randn(ndim) for i in range(nwalkers)]
#		for i in range(nwalkers):
#			pos[i][26]=1.0
#			pos[i][32]=1.0
#			pos[i][40]=0.0
#			pos[i][48]=0.0
#			pos[i][56]=0.0
		print 'Running MCMC.........................'
		# run xxx steps as a burn-in
		pos, prob, state=sampler.run_mcmc(pos, burnin)
		# reset the chain to remove the burn-in samples
		sampler.reset()
		print 'Running MCMC....burn-in completed....'
		# starting from the final position in the burn-in chain, sample for xxx
		sampler.run_mcmc(pos, post_burnin, rstate0=state)
		print 'Running MCMC........Completed........'
		print("Mean acceptance fraction: {0:.6f}".format(np.mean(sampler.acceptance_fraction)))
		print 'Autocorrelation time is ', sampler.acor
		pool.close()
		samples=sampler.chain[:, :, :].reshape((-1, ndim))
		
		res=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))	
		p_mcmc=[]
		p_err_mcmc=[]
		for i in range(ndim):
			p_mcmc.append(res[i][0])
			p_err_mcmc.append(0.5*(res[i][1]+res[i][2]))
		
		spar_index=0
		for i_source in np.arange(self.n_source):
			spar_index+=self.n_source_par_dict[self.sourcemodel[i_source]]
		lpar_index=0
		for i_lens in np.arange(self.n_lens):
			lpar_index+=self.n_lens_par_dict[self.lensmodel[i_lens]]

		spar=p_mcmc[0:spar_index]
		lpar=p_mcmc[spar_index:spar_index+lpar_index]
		ppar=p_mcmc[spar_index+lpar_index:]

		spar_err=p_err_mcmc[0:spar_index]
		lpar_err=p_err_mcmc[spar_index:spar_index+lpar_index]
		ppar_err=p_err_mcmc[spar_index+lpar_index:]

		print 'lpar ', lpar
		print 'lpar_err ', lpar_err
		print 'spar ', spar
		print 'spar_err ', spar_err
		print 'ppar ', ppar
		print 'ppar_err ', ppar_err

		self.lpar=lpar
		self.lpar_err=lpar_err
		self.spar=spar
		self.spar_err=spar_err
		self.ppar=ppar
		self.ppar_err=ppar_err

		self.generate_models(self.n_lens, lpar, self.lensmodel, self.n_source, spar, self.sourcemodel, self.n_phot, ppar, self.photmodel, self.PSF, self.PSF_csr)

		self.samples=samples
		self.dof=self.I_data1.size-(np.where(self.I_invvar1==0.0))[0].size
		if (self.n_source == 1 and self.sourcemodel[0] == 'pix'):
			pars=np.concatenate((spar, lpar, ppar))
			self.chi2=self.residual_pix(pars, self.x1, self.y1, self.n_lens, self.lensmodel, self.lens_mask_2d, self.used_pixels, self.used_pixels_mask, data=self.I_data1, I_invvar=self.I_invvar1, PSF=self.PSF, PSF_csr=self.PSF_csr)
		else:
			self.chi2=((self.I_data1-self.I_fit-self.I_phot)**2.0*self.I_invvar1).sum()

	def run_lfit(self, parfile):
		self.setup_lfit(parfile)
		if (self.fit_mode == 1):
			self.run_lmfit()
		elif (self.fit_mode == 2):
			self.run_mcmc()
		elif (self.fit_mode == 3):
			self.run_mcmc_mpi()

		self.save_fit()

"""
Usage

junk=lfit()

if len(sys.argv) !=2:
   print 'Error!'
   print 'Usage: python lfit_script.py <parfile>'
   sys.exit(1)

junk.run_lfit(sys.argv[1])
"""
