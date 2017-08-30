import os, commands, pyfits
import numpy as np
from lmfit import Parameters, minimize, fit_report, conf_interval, printfuncs
import lens_funcs_gui as lf
import tkMessageBox
from Tkinter import *
from tkFileDialog import *
import ttk
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import emcee
import string
from scipy.interpolate import griddata
import scipy.special as sf

myargs = {'origin': 'lower'}

def get_cycle_list():
# search directories under HST_DIR to get the cycle list
	command='pwd'
	current_path=commands.getoutput(command)
	flag=1
	fpath=os.getenv('HST_DIR')
	if fpath is None: 
		print "Please set up $HST_DIR."
		exit()

	command="cd "+fpath+"/ && ls -d */ | sed 's/\/$//'"
	output=commands.getoutput(command)
	dir_list=output.split( )
	command='cd '+current_path
	cycle_list=[]
	for i in np.arange(len(dir_list)):
		try:
			cycle_list.append(int(dir_list[i]))
		except Exception:
			continue
	return cycle_list
	
def get_program_list(cycle):
# dig into a particular cycle to get the program list
	command='pwd'
	current_path=commands.getoutput(command)
	fpath=os.getenv('HST_DIR')
	command="cd "+fpath+"/"+str(cycle)+"/ && ls -d */ | sed 's/\/$//'"
	output=commands.getoutput(command)
	dir_list=output.split()
	command='cd '+current_path
	program_list=[]
	for i in np.arange(len(dir_list)):
		try:
			program_list.append(int(dir_list[i]))
		except Exception:
			continue
	return program_list

def update_program_list(*args):
# update the program list once a new cycle is selected
	clear_parameter_panel()
	program_list=get_program_list(var_cycle.get())
	var_program.set(program_list[0])
	menu=program_menu['menu']
	menu.delete(0, END)
	for program in program_list:
		menu.add_command(label=program, command=lambda program1=program: var_program.set(program1))
#	getvalues()

def get_target_list(cycle, program):
# get the target list in a particular cycle, program combination
	command='pwd'
	current_path=commands.getoutput(command)
	fpath=os.getenv('HST_DIR')
	command="cd "+fpath+"/"+str(cycle)+"/"+str(program)+"/ && ls */*biz.fits"
	visit_list=[]
	target_list=[]
	status, output=commands.getstatusoutput(command)
	if status==0:
		file_list=output.split( )
		for i in np.arange(len(file_list)):
			tmp=file_list[i].split('/')
			visit_list.append(tmp[0])
			if program == 13000:
				target_list.append(tmp[1][0:14])
			elif program == 14189:
				obj_name='SDSSJ'+tmp[1][tmp[1].find('SLACSJ')+6:tmp[1].find('SLACSJ')+24]
				target_list.append(obj_name)
			else:
				target_list.append(tmp[1][0:15])
				
	command='cd '+current_path
	return (visit_list, target_list)

def update_target_list(*args):
# update the target list when a new cycle or program is selected
	global visit
	clear_parameter_panel()
	visit_list, target_list=get_target_list(var_cycle.get(), var_program.get())
	target_tree.delete(*target_tree.get_children())
	for i in np.arange(len(visit_list)):
		target_tree.insert("", "end", text=str(i+1), values=(visit_list[i], target_list[i]))
	try:
		visit=visit_list[0]
	except Exception:
		print 'No targets in this program'
#	getvalues()

def select_visit(event):
# select a particular visit from the menu
	global visit, n_lens, lensmodel, current_row_lens, current_column_lens, frame_lens, n_source, sourcemodel, current_row_source, current_column_source, frame_source
	clear_parameter_panel()
	var_show_residual_bspline.set(0)
	var_show_residual_sersic.set(0)
	select=target_tree.selection()
	select_info=target_tree.set(select)
	visit=select_info['visitID']
	getvalues()

def update_comments():
	comments=(str(comment_text.get(1.0, END))).replace("\n", "") # convert unicode to string and chop off the new line "\n"
	hdulist=pyfits.open(fname, mode='update')
	primary_header=hdulist[0].header
	primary_header['sys_cmts']=comments
	hdulist.flush()
	hdulist.close()
	print 'comments updated for '+fname

def master_on_closing():
# closing the main window dialog
	if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
		master.quit()

def on_closing(popup):
	popup.destroy()
	plt.close('all')

def clear_parameter_panel():
# clear the lens and source parameter panels
	global n_lens, lensmodel, current_row_lens, current_column_lens, lens_frame_list, drop_lens_button_list, lpar_entry_list, lpar_fix_list, n_phot, photmodel, current_row_phot, current_column_phot, phot_frame_list, drop_phot_button_list, ppar_entry_list, ppar_fix_list, n_source, sourcemodel, current_row_source, current_column_source, source_frame_list, drop_source_button_list, spar_entry_list, spar_fix_list
	try:
		n_lens=0
		lensmodel=[]
		current_row_lens=0
		current_column_lens=0
		lens_frame_list=[]
		drop_lens_button_list=[]
		lpar_entry_list=[]
		lpar_fix_list=[]

		n_phot=0
		photmodel=[]
		current_row_phot=0
		current_column_phot=0
		phot_frame_list=[]
		drop_phot_button_list=[]
		ppar_entry_list=[]
		ppar_fix_list=[]

		n_source=0
		sourcemodel=[]
		current_row_source=0
		current_column_source=0
		source_frame_list=[]
		drop_source_button_list=[]
		spar_entry_list=[]
		spar_fix_list=[]

		for widget in lens_panel.winfo_children():
			widget.destroy()
		for widget in phot_panel.winfo_children():
			widget.destroy()
		for widget in source_panel.winfo_children():
			widget.destroy()
	except Exception:
		pass

def getvalues():
# when a visit is selected, get the imaging data from the biz file
	global cycle, program, visit, fname, I_data, I_lens, I_image, I_invvar, jmask, gmask, n_x, n_y, x, y, lens_xcen, lens_ycen, lens_q0, lens_pa0, PSF, tpsf, obj_name, dpix, fpath
	cycle=var_cycle.get()
	program=var_program.get()
	fpath=os.getenv('HST_DATAROOT')+'/'+str(cycle)+'/'+str(program)+'/'+visit
	command='ls '+fpath+'/*biz.fits'
	fname=commands.getstatusoutput(command)[1]
	if program<13000:
		obj_name=fname[fname.find('SLACSJ'):fname.find('SLACSJ')+14+1]
	elif program==13000:
		obj_name=fname[fname.find('SHARP'):fname.find('SHARP')+14]
	elif program==14189:
		obj_name='SDSSJ'+fname[fname.find('SLACSJ')+6:fname.find('SLACSJ')+24]
	hdulist=pyfits.open(fname)
	I_data=hdulist[0].data
	primary_header=hdulist[0].header
	I_invvar=hdulist[1].data
	PSF=hdulist[3].data
	PSF=PSF.astype(float) # the original PSF is float32, while the default is float64 in python
	mask=hdulist[4].data
	jmask=hdulist[5].data
	gmask=hdulist[6].data
	I_lens=hdulist[7].data
	tpsf=hdulist[9].data
	I_image=I_data-I_lens
#	I_invvar=I_invvar*jmask

	hdr=hdulist[7].header
	n_x=hdr['NAXIS1']
	n_y=hdr['NAXIS2']
	if program==13000:
		dpix=0.01
	elif program==13100:
		dpix=0.396
	elif program==14189 or program == 11602:
		dpix=0.04
	elif program==10000:
		dpix=0.185	
	elif program==11000:
		dpix=1
	else:
		dpix=0.05
		
	x=np.outer(np.ones(n_y), np.arange(n_x)-n_x/2)*dpix
	y=np.outer(np.arange(n_y)-n_y/2, np.ones(n_x))*dpix
	
	lens_xcen=hdr['xcenter']*dpix-(n_x/2)*dpix
	lens_ycen=hdr['ycenter']*dpix-(n_y/2)*dpix
	lens_q0=hdr['axisrat']
	lens_pa0=hdr['axisangl']
	if lens_q0 > 1.0:
		lens_q0=1.0/lens_q0
		lens_pa0=lens_pa0+90.0
	
	if lens_pa0 < 0.0:
		lens_pa0=lens_pa0+180.0
	
	if lens_pa0 > 180.0:
		lens_pa0=lens_pa0-180.0

	scale_flag=1
	if program == 10000:
		scale_flag=0
	if scale_flag:
		width_shell=15
		wh=np.where(((x <= width_shell*dpix-(n_x/2)*dpix) +(x >= (n_x/2)*dpix-width_shell*dpix)+(y <= width_shell*dpix-(n_y/2)*dpix)+(y >= (n_y/2)*dpix-width_shell*dpix))*(jmask > 0))
		scale=1.0/(np.mean(np.sqrt(I_invvar[wh]))*np.std(I_data[wh]*jmask[wh]))
		I_invvar=I_invvar*scale**2.

	if 'sys_cmts' in primary_header:
		comment_text.delete(1.0, END)
		comment_text.insert(END, primary_header['sys_cmts'])
	else:
		comment_text.delete(1.0, END)
		comment_text.insert(END, 'Insert comments here. ')
	hdulist.close()

def show_image():
# show the image of a selected system
	global popup, frame_popup, canvas, I_data1, I_image1, gmask1, hw_entry, fig, ax, current_cmap, hw, x_entry, y_entry, val_entry, remove_foreground_status, mouse_mode_text, mask_mode_text, brush_text, mouse_mode_button, mask_mode_button, brush_button, vmin, vmax, draw_mask_start

	popup=Toplevel()
	popup.wm_title("Preview")
	popup.protocol("WM_DELETE_WINDOW", lambda: popup.destroy())
	frame_popup=Frame(popup)
	for i in range(30):
		frame_popup.grid_rowconfigure(i, weight=1)
		frame_popup.grid_columnconfigure(i, weight=1)
	default_cmap=plt.get_cmap('jet')
	current_cmap=default_cmap
	frame_popup.bind("<B1-Motion>", lambda event: adjust_image(event))
	frame_popup.pack()

	mouse_mode_label=Label(frame_popup, text="Mouse", width=12)
	mouse_mode_label.grid(row=9, column=0)
	mouse_mode_list=['color', 'mask']
	mouse_mode_text=StringVar()
	mouse_mode_text.set('color')
	mouse_mode_button=ttk.Combobox(frame_popup, textvariable=mouse_mode_text, state='readonly')
	mouse_mode_button.grid(row=9, column=1, sticky=W)
	mouse_mode_button['values']=tuple(mouse_mode_list)
	mouse_mode_button.config(width=10)

	mask_mode_label=Label(frame_popup, text="Mask", width=12)
	mask_mode_label.grid(row=9, column=2)
	mask_mode_list=['none', 'feature', 'junk', 'all']
	mask_mode_text=StringVar()
	mask_mode_text.set('none')
	mask_mode_button=ttk.Combobox(frame_popup, textvariable=mask_mode_text, state='readonly')
	mask_mode_button.grid(row=9, column=3, sticky=W)
	mask_mode_button['values']=tuple(mask_mode_list)
	mask_mode_button.config(width=10)
	mask_mode_button.bind("<<ComboboxSelected>>", lambda event: show_mask(event))

	brush_label=Label(frame_popup, text="Brush", width=12)
	brush_label.grid(row=9, column=4)
	brush_list=['1', '2', '4', '5', '10', '20']
	brush_text=StringVar()
	brush_text.set('1')
	brush_button=ttk.Combobox(frame_popup, textvariable=brush_text, state='readonly')
	brush_button.grid(row=9, column=5, sticky=W)
	brush_button['values']=tuple(brush_list)
	brush_button.config(width=10)

	save_button=Button(frame_popup, text='Save', width=10, command=save_masks)
	save_button.grid(row=9, column=7, sticky=W)

	Label(frame_popup, text="HW", width=12).grid(row=10, column=0)
	hw_entry=Entry(frame_popup, width=12)
	hw_entry.grid(row=10, column=1, sticky=W)
	if 'hw' in globals():
		hw_entry.insert(0, hw)
	else:
		hw_entry.insert(0, 80)
	hw_entry.bind('<Return>', resize_image)
#	blank_label=Label(frame_popup, width=12, text=' ').grid(row=0, column=10)

	x_label=Label(frame_popup, text="x", width=12)
	x_label.grid(row=11, column=0)
	x_entry=Entry(frame_popup, width=12)
	x_entry.grid(row=11, column=1, sticky=W)
	x_entry.insert(0, 0.0)
	Label(frame_popup, text="y", width=12).grid(row=11, column=2)
	y_entry=Entry(frame_popup, width=12)
	y_entry.grid(row=11, column=3, sticky=W)
	y_entry.insert(0, 0.0)
	Label(frame_popup, text="value", width=12).grid(row=11, column=4)
	val_entry=Entry(frame_popup, width=12)
	val_entry.grid(row=11, column=5, sticky=W)

	show_residual_bspline_button=Checkbutton(frame_popup, variable=var_show_residual_bspline, justify=LEFT, text="Remove foreground (bspline)", command=show_residual)
	show_residual_bspline_button.grid(row=12, column=0, columnspan=2, sticky=W)
	show_residual_bspline_button.deselect()

	show_residual_sersic_button=Checkbutton(frame_popup, variable=var_show_residual_sersic, justify=LEFT, text="Remove foreground (sersic)", command=show_residual)
	show_residual_sersic_button.grid(row=13, column=0, columnspan=2, sticky=W)
	show_residual_sersic_button.deselect()

	fig=plt.figure(figsize=(8, 8))
	canvas=FigureCanvasTkAgg(fig, master=frame_popup)
	NavigationToolbar2TkAgg(canvas, popup)
	canvas.get_tk_widget().grid(row=0, column=0, rowspan=8, columnspan=10, sticky='NW', padx=10, pady=10)	
	ax=fig.add_subplot(111)

	hw=int(hw_entry.get())
	if var_show_residual_bspline.get() == 1:
		remove_foreground_status=1
	elif var_show_residual_sersic.get() == 1:
		remove_foreground_status=2
	else:
		remove_foreground_status=0	

	vmin=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1].min()
	vmax=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1].max()
	update_image()
	canvas.get_tk_widget().focus_set()

	draw_mask_start=0
	button1_press=canvas.mpl_connect('button_press_event', button_press_event)
	button1_release=canvas.mpl_connect('button_release_event', button_release_event)
	mouse_position=canvas.mpl_connect('motion_notify_event', mouse_position_event)
	
def show_mask(event):
	update_image()
	canvas.get_tk_widget().focus_set()

def button_press_event(event):
# get the corresponding source position of an arbitrary position in the image plane
# that is control+click'ed based on the current lens model
	global draw_mask_start
	canvas.get_tk_widget().focus_set()
	if ((event.key == 'control' or event.key == 'ctrl+') and mouse_mode_button.get() == 'color'):
		theta_x=event.xdata
		theta_y=event.ydata
		if n_lens == 0:
			print "Null lens model"
			pass
		else:
			lpar=[]
			for i_lens in np.arange(n_lens):
				if lensmodel[i_lens] =='sie' or lensmodel[i_lens] == 'sple' or lensmodel[i_lens] == 'softie':
					setting=params_setting((lpar_entry_list[i_lens]['b_SIE_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_pa_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_q_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['gamma_entry']).get())	
					lpar.append(setting['value'])		
				elif lensmodel[i_lens] == 'pm':
					setting=params_setting((lpar_entry_list[i_lens]['theta_e_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
					lpar.append(setting['value'])	
				elif lensmodel[i_lens] == 'ext. shear':
					setting=params_setting((lpar_entry_list[i_lens]['shear_amp_entry']).get())										
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['shear_pa_entry']).get())										
					lpar.append(setting['value'])
				elif lensmodel[i_lens] =='snfw':

					setting=params_setting((lpar_entry_list[i_lens]['kappa_s_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_pa_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['lens_q_entry']).get())	
					lpar.append(setting['value'])
					setting=params_setting((lpar_entry_list[i_lens]['r_s_entry']).get())	
					lpar.append(setting['value'])		
		try:	
			(alpha_x, alpha_y)=lf.def_total(theta_x, theta_y, lpar, lensmodel, n_lens=n_lens)
			update_parameter_value(spar_entry_list[-1], 'source_xcen_entry', theta_x-alpha_x)
			update_parameter_value(spar_entry_list[-1], 'source_ycen_entry', theta_y-alpha_y)
		except Exception:
			pass
	elif event.key == 'shift':
		draw_mask_start=2
	else:
		draw_mask_start=1

def button_release_event(event):
	global draw_mask_start
	draw_mask_start=0

def mouse_position_event(event):
	global gmask1, jmask1
	canvas.get_tk_widget().focus_set()
	if draw_mask_start == 0:
		if event.xdata == None or event.ydata == None:
			pass
		else:
			ind=np.unravel_index(np.argmin((x1-event.xdata)**2.+(y1-event.ydata)**2.), x1.shape)
			x_entry.delete(0, END)
			x_entry.insert(0, x1[ind])
			y_entry.delete(0, END)
			y_entry.insert(0, y1[ind])
			val_entry.delete(0, END)
			if remove_foreground_status == 0:
				val_entry.insert(0, I_data1[ind])
			elif remove_foreground_status == 1:
				val_entry.insert(0, I_image1[ind])	
			elif remove_foreground_status == 2:
				val_entry.insert(0, (I_data1-I_phot)[ind])
	elif draw_mask_start > 0:
		if mouse_mode_button.get() == 'color':
			pass
		elif mouse_mode_button.get() == 'mask':
			if mask_mode_button.get() == 'all' or mask_mode_button.get() == 'none':
				pass
			elif mask_mode_button.get() == 'feature':
				ind=np.unravel_index(np.argmin((x1-event.xdata)**2.+(y1-event.ydata)**2.), x1.shape)
				if draw_mask_start == 1:
					gmask1[ind[0]-int(brush_button.get()):ind[0]+int(brush_button.get())+1, ind[1]-int(brush_button.get()):ind[1]+int(brush_button.get())+1]=0
				elif draw_mask_start == 2:
					gmask1[ind[0]-int(brush_button.get()):ind[0]+int(brush_button.get())+1, ind[1]-int(brush_button.get()):ind[1]+int(brush_button.get())+1]=1
				update_image()
			elif mask_mode_button.get() == 'junk':
				ind=np.unravel_index(np.argmin((x1-event.xdata)**2.+(y1-event.ydata)**2.), x1.shape)
				if draw_mask_start == 1:
					jmask1[ind[0]-int(brush_button.get()):ind[0]+int(brush_button.get())+1, ind[1]-int(brush_button.get()):ind[1]+int(brush_button.get())+1]=0
				elif draw_mask_start == 2:
					jmask1[ind[0]-int(brush_button.get()):ind[0]+int(brush_button.get())+1, ind[1]-int(brush_button.get()):ind[1]+int(brush_button.get())+1]=1
				update_image()

def save_masks():
	jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]=jmask1
	gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]=gmask1
	temp=pyfits.open(fname, mode='update')
	temp[5].data=jmask
	temp[6].data=gmask
	temp.flush()
	temp.close()
	print "Saved masks to "+fname

def adjust_image(event):
# adjust the color scheme of the image by moving the mouse with left button being held down
	global current_cmap, hw_entry, hw, vmin, vmax
	if var_show_residual_bspline.get() == 1:
		remove_foreground_status=1
	elif var_show_residual_sersic.get() == 1:
		remove_foreground_status=2
	else:
		remove_foreground_status=0
	if remove_foreground_status == 0:
		I_display=I_data1
	elif remove_foreground_status == 1:
		I_display=I_image1
	elif remove_foreground_status == 2:
		I_display=I_phot
	med=np.median(I_display)
	std=np.std(I_display)
	imax=min(med+10.*std, I_display.max())
	imin=max(med-2.*std, I_display.min())
	x=min(event.x/1000., 2.)
	y=min(event.y/1000., 2.)
	vmax=(imax-imin)*x+imin
	vmin=(imax-imin)*y+imin
	if vmin > vmax:
		vmin=vmax

	update_image()
	canvas.get_tk_widget().focus_set()

def resize_image(event):
# resize the image based on the input hw
	global hw, remove_foreground_status
	hw=int(hw_entry.get())
	update_image()
	canvas.get_tk_widget().focus_set()

def show_residual():
# remove the foreground light
	global hw, remove_foreground_status
	hw=int(hw_entry.get())
	if var_show_residual_bspline.get() == 1:
		remove_foreground_status=1
	elif var_show_residual_sersic.get() == 1:
		remove_foreground_status=2
	else:
		remove_foreground_status=0
	update_image()
	canvas.get_tk_widget().focus_set()
	
def update_image():
	global I_data1, I_image1, x1, y1, jmask1, gmask1
	I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_image1=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	jmask1=jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	ext_lens=[x1.min(), x1.max(), y1.min(), y1.max()]
	
	ax.clear()
	if remove_foreground_status == 0:
		im=ax.imshow(I_data1, extent=ext_lens, vmin=vmin, vmax=vmax, cmap=current_cmap, **myargs)
	elif remove_foreground_status == 1:
		im=ax.imshow(I_image1, extent=ext_lens, vmin=vmin, vmax=vmax, cmap=current_cmap, **myargs)	
	elif remove_foreground_status == 2:
		im=ax.imshow(I_data1-I_phot, extent=ext_lens, vmin=vmin, vmax=vmax, cmap=current_cmap, **myargs)
	ax.autoscale(False)
	if (mask_mode_button.get() == 'junk'):
		wh_junk=np.where(jmask1 == 0.)
		ax.plot(x1[wh_junk], y1[wh_junk], 'sg', markerfacecolor='none', markeredgecolor='g')
	elif (mask_mode_button.get() == 'feature'):
		wh_feature=np.where(gmask1 == 0.)
		ax.plot(x1[wh_feature], y1[wh_feature], 'sr', markerfacecolor='none', markeredgecolor='r')
	elif (mask_mode_button.get() == 'all'):
		wh_junk=np.where(jmask1 == 0.)
		ax.plot(x1[wh_junk], y1[wh_junk], 'sg', markerfacecolor='none', markeredgecolor='g')
		wh_feature=np.where(gmask1 == 0.)
		ax.plot(x1[wh_feature], y1[wh_feature], 'sr', markerfacecolor='none', markeredgecolor='r')		
	ax.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	plt.subplots_adjust(left=0.05, right=0.98, bottom=0.15, wspace=0.0)	
	canvas.draw()
	canvas.get_tk_widget().focus_set()
	
def add_a_lens(event):
	global n_lens, lensmodel, current_row_lens, lens_model_label
	n_lens+=1
	lensmodel.append(lensmodel_list[add_lens_model_button.current()])
	add_lens_model_text.set('Add a lens')
	if n_lens==1:
		lens_model_label=Label(lens_panel, text='Lens Model')
		lens_model_label.config(relief=RAISED, bd=2)
		lens_model_label.grid(row=current_row_lens, column=0, sticky='w')
		current_row_lens=current_row_lens+1
	add_a_lens_entry(current_row_lens, 0)
	current_row_lens=current_row_lens+1

def add_a_lens_entry(row, column):
# add a lens entry
# the parameter entries are saved as a dictionary for each lens

	current_lens_frame=Frame(lens_panel)
	for i in range(30):
		current_lens_frame.grid_rowconfigure(i, weight=1)
		current_lens_frame.grid_columnconfigure(i, weight=1)
	current_lens_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	lens_frame_list.append(current_lens_frame)
	current_drop_lens_button=Button(current_lens_frame, text="Drop", width=12, command=lambda: drop_a_lens(row))
	current_drop_lens_button.grid(row=0, column=column+12)
	drop_lens_button_list.append(current_drop_lens_button)
	if lensmodel[-1] == 'sie' or lensmodel[-1] == 'sple':
		Label(current_lens_frame, text="b_SIE", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		Label(current_lens_frame, text="lens_pa", width=12).grid(row=0, column=column+6)
		Label(current_lens_frame, text="lens_q", width=12).grid(row=0, column=column+8)
		Label(current_lens_frame, text="gamma", width=12).grid(row=0, column=column+10)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar3_entry=Entry(current_lens_frame, width=12)
		current_lpar4_entry=Entry(current_lens_frame, width=12)
		current_lpar5_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)
		current_lpar3_entry.grid(row=0, column=column+7)
		current_lpar4_entry.grid(row=0, column=column+9)
		current_lpar5_entry.grid(row=0, column=column+11)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)
		current_lpar3_entry.delete(0, END)
		current_lpar4_entry.delete(0, END)
		current_lpar5_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=1.0, min=0.0")
		current_lpar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_lpar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		current_lpar3_entry.insert(0, "value=%10.6f, min=0.01, max=180.0" %lens_pa0)
		current_lpar4_entry.insert(0, "value=%10.6f, min=0.01, max=1.0" %lens_q0)
		current_lpar5_entry.insert(0, "value=1.0")
		if lensmodel[-1] =='sie':
			current_lpar5_entry.insert(END, ", vary=0")
		lpar_entry_i={'b_SIE_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry, 'lens_pa_entry': current_lpar3_entry, 'lens_q_entry': current_lpar4_entry, 'gamma_entry': current_lpar5_entry}
		lpar_entry_list.append(lpar_entry_i)
	elif lensmodel[-1] == 'pm':
		Label(current_lens_frame, text="theta_e", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=1.0, min=0.0")
		current_lpar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_lpar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		lpar_entry_i={'theta_e_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry}
		lpar_entry_list.append(lpar_entry_i)
	elif lensmodel[-1] == 'ext. shear':
		Label(current_lens_frame, text="shear_amp", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="shear_pa", width=12).grid(row=0, column=column+2)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=0.05")
		current_lpar1_entry.insert(0, "value=30.0, min=0.01, max=180.0")
		
		lpar_entry_i={'shear_amp_entry':current_lpar0_entry, 'shear_pa_entry': current_lpar1_entry}
		lpar_entry_list.append(lpar_entry_i)
	elif lensmodel[-1] == 'snfw':
		Label(current_lens_frame, text="kappa_s", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		Label(current_lens_frame, text="lens_pa", width=12).grid(row=0, column=column+6)
		Label(current_lens_frame, text="lens_q", width=12).grid(row=0, column=column+8)
		Label(current_lens_frame, text="r_s", width=12).grid(row=0, column=column+10)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar3_entry=Entry(current_lens_frame, width=12)
		current_lpar4_entry=Entry(current_lens_frame, width=12)
		current_lpar5_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)
		current_lpar3_entry.grid(row=0, column=column+7)
		current_lpar4_entry.grid(row=0, column=column+9)
		current_lpar5_entry.grid(row=0, column=column+11)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)
		current_lpar3_entry.delete(0, END)
		current_lpar4_entry.delete(0, END)
		current_lpar5_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=1.0, min=0.0")
		current_lpar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_lpar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		current_lpar3_entry.insert(0, "value=0.0, vary=0")
		current_lpar4_entry.insert(0, "value=1.0, vary=0")
		current_lpar5_entry.insert(0, "value=5.0, min=0.0")
		lpar_entry_i={'kappa_s_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry, 'lens_pa_entry': current_lpar3_entry, 'lens_q_entry': current_lpar4_entry, 'r_s_entry': current_lpar5_entry}
		lpar_entry_list.append(lpar_entry_i)

def drop_a_lens(row):
	global n_lens, current_row_lens
	print 'Dropping lens %i' %(row)
	for widget in lens_frame_list[row-1].winfo_children():
		widget.destroy()
	n_lens-=1
	lensmodel.pop(row-1)
	lens_frame_list.pop(row-1)
	drop_lens_button_list.pop(row-1)
	lpar_entry_list.pop(row-1)
	drop_lens_number=row
	for i_lens in range(drop_lens_number-1, n_lens):
		lens_frame_list[i_lens].grid_forget()
		this_row=(i_lens+1)
		lens_frame_list[i_lens].grid(row=this_row, column=0, rowspan=1, sticky='ew')
		drop_lens_button_list[i_lens].config(text="Drop", command=lambda this_row=this_row: drop_a_lens(this_row))
	if n_lens==0:
		lens_model_label.destroy()
		current_row_lens=current_row_lens-1
	current_row_lens=current_row_lens-1

def add_a_phot(event):
	global n_phot, photmodel, current_row_phot, phot_model_label
	n_phot+=1
	photmodel.append(photmodel_list[add_phot_model_button.current()])
	add_phot_model_text.set('Add a phot')
	if n_phot==1:
		phot_model_label=Label(phot_panel, text='Phot. Model')
		phot_model_label.config(relief=RAISED, bd=2)
		phot_model_label.grid(row=current_row_phot, column=0, sticky='w')
		current_row_phot=current_row_phot+1
	add_a_phot_entry(current_row_phot, 0)
	current_row_phot=current_row_phot+1

def add_a_phot_entry(row, column):
# add a phot entry
# the parameter entries are saved as a dictionary for each phot

	current_phot_frame=Frame(phot_panel)
	for i in range(30):
		current_phot_frame.grid_rowconfigure(i, weight=1)
		current_phot_frame.grid_columnconfigure(i, weight=1)
	current_phot_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	phot_frame_list.append(current_phot_frame)
	if photmodel[-1] == 'sersic':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="sigma", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="n", width=12).grid(row=0, column=column+12)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+14)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar7_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)
		current_ppar7_entry.grid(row=0, column=column+15)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)
		current_ppar7_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=10.0")
		current_ppar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_ppar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		current_ppar3_entry.insert(0, "value=1.0")
		current_ppar4_entry.insert(0, "value=%10.6f, min=1, max=180.0" %lens_pa0)
		current_ppar5_entry.insert(0, "value=%10.6f, min=0.1, max=1.0" %lens_q0)
		current_ppar6_entry.insert(0, "value=4.0, min=0.01")
		current_ppar7_entry.insert(0, "value=0.01, min=-0.02, max=0.02")
		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_sigma_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_n_entry': current_ppar6_entry, 'phot_I0_entry': current_ppar7_entry}
		ppar_entry_list.append(ppar_entry_i)

		current_drop_phot_button=Button(current_phot_frame, text="Drop", width=12, command=lambda: drop_a_phot(row))
		current_drop_phot_button.grid(row=0, column=column+22)
		drop_phot_button_list.append(current_drop_phot_button)
	elif photmodel[-1] == 'csersic':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="sigma", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="n", width=12).grid(row=0, column=column+12)
		Label(current_phot_frame, text="rc", width=12).grid(row=0, column=column+14)
		Label(current_phot_frame, text="alpha", width=12).grid(row=0, column=column+16)
		Label(current_phot_frame, text="gamma", width=12).grid(row=0, column=column+18)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+20)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar7_entry=Entry(current_phot_frame, width=12)
		current_ppar8_entry=Entry(current_phot_frame, width=12)
		current_ppar9_entry=Entry(current_phot_frame, width=12)
		current_ppar10_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)
		current_ppar7_entry.grid(row=0, column=column+15)
		current_ppar8_entry.grid(row=0, column=column+17)
		current_ppar9_entry.grid(row=0, column=column+19)
		current_ppar10_entry.grid(row=0, column=column+21)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)
		current_ppar7_entry.delete(0, END)
		current_ppar8_entry.delete(0, END)
		current_ppar9_entry.delete(0, END)
		current_ppar10_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=10.0")
		current_ppar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_ppar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		current_ppar3_entry.insert(0, "value=1.0")
		current_ppar4_entry.insert(0, "value=%10.6f, min=1, max=180.0" %lens_pa0)
		current_ppar5_entry.insert(0, "value=%10.6f, min=0.1, max=1.0" %lens_q0)
		current_ppar6_entry.insert(0, "value=1.0, min=0.0")
		current_ppar7_entry.insert(0, "value=0.1, min=0.0")
		current_ppar8_entry.insert(0, "value=10.0, min=0.0")
		current_ppar9_entry.insert(0, "value=1.0, min=0.0")
		current_ppar10_entry.insert(0, "value=0.01, min=-0.02, max=0.02")

		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_sigma_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_n_entry': current_ppar6_entry, 'phot_rc_entry': current_ppar7_entry, 'phot_alpha_entry': current_ppar8_entry, 'phot_gamma_entry': current_ppar9_entry, 'phot_I0_entry': current_ppar10_entry}
		ppar_entry_list.append(ppar_entry_i)

		current_drop_phot_button=Button(current_phot_frame, text="Drop", width=12, command=lambda: drop_a_phot(row))
		current_drop_phot_button.grid(row=0, column=column+22)
		drop_phot_button_list.append(current_drop_phot_button)
	elif photmodel[-1] == 'hernquist':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="rs", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+12)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=10.0")
		current_ppar1_entry.insert(0, "value=%10.6f" %lens_xcen)
		current_ppar2_entry.insert(0, "value=%10.6f" %lens_ycen)
		current_ppar3_entry.insert(0, "value=0.5, min=0.0")
		current_ppar4_entry.insert(0, "value=%10.6f, min=1, max=180.0" %lens_pa0)
		current_ppar5_entry.insert(0, "value=%10.6f, min=0.1, max=1.0" %lens_q0)
		current_ppar6_entry.insert(0, "value=0.01, min=-0.02, max=0.02")
		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_rs_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_I0_entry': current_ppar6_entry}
		ppar_entry_list.append(ppar_entry_i)

		current_drop_phot_button=Button(current_phot_frame, text="Drop", width=12, command=lambda: drop_a_phot(row))
		current_drop_phot_button.grid(row=0, column=column+22)
		drop_phot_button_list.append(current_drop_phot_button)

def drop_a_phot(row):
	global n_phot, current_row_phot
	print 'Dropping phot %i' %(row)
	for widget in phot_frame_list[row-1].winfo_children():
		widget.destroy()
	n_phot-=1
	photmodel.pop(row-1)
	phot_frame_list.pop(row-1)
	drop_phot_button_list.pop(row-1)
	ppar_entry_list.pop(row-1)
	drop_phot_number=row
	for i_phot in range(drop_phot_number-1, n_phot):
		phot_frame_list[i_phot].grid_forget()
		this_row=(i_phot+1)
		phot_frame_list[i_phot].grid(row=this_row, column=0, rowspan=1, sticky='ew')
		drop_phot_button_list[i_phot].config(text="Drop", command=lambda this_row=this_row: drop_a_phot(this_row))
	if n_phot==0:
		phot_model_label.destroy()
		current_row_phot=current_row_phot-1
	current_row_phot=current_row_phot-1

def add_a_source(event):
	global n_source, sourcemodel, current_row_source, source_model_label
	n_source+=1
	sourcemodel.append(sourcemodel_list[add_source_model_button.current()])
	add_source_model_text.set('Add a source')
	if n_source==1:
		source_model_label=Label(source_panel, text='Source Model')
		source_model_label.config(relief=RAISED, bd=2)
		source_model_label.grid(row=current_row_source, column=0, sticky='w')
		current_row_source=current_row_source+1
	add_a_source_entry(current_row_source, 0)
	current_row_source=current_row_source+1

def add_a_source_entry(row, column):
# add a source entry
# the parameter entries are saved as a dictionary for each source
	current_source_frame=Frame(source_panel)
	current_source_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	for i in range(30):
		current_source_frame.grid_rowconfigure(i, weight=1)
		current_source_frame.grid_columnconfigure(i, weight=1)
	source_frame_list.append(current_source_frame)
	current_drop_source_button=Button(current_source_frame, text="Drop", width=12, command=lambda: drop_a_source(row))
	current_drop_source_button.grid(row=0, column=column+14)
	drop_source_button_list.append(current_drop_source_button)
	if sourcemodel[-1] == 'sersic' or sourcemodel[-1] == 'gaussian':
		Label(current_source_frame, text="src_amp", width=12).grid(row=0, column=column)
		Label(current_source_frame, text="src_xcen", width=12).grid(row=0, column=column+2)
		Label(current_source_frame, text="src_ycen", width=12).grid(row=0, column=column+4)
		Label(current_source_frame, text="src_sigma", width=12).grid(row=0, column=column+6)
		Label(current_source_frame, text="src_pa", width=12).grid(row=0, column=column+8)
		Label(current_source_frame, text="src_q", width=12).grid(row=0, column=column+10)
		Label(current_source_frame, text="src_n", width=12).grid(row=0, column=column+12)

		current_source_amp_entry=Entry(current_source_frame, width=12)
		current_source_xcen_entry=Entry(current_source_frame, width=12)
		current_source_ycen_entry=Entry(current_source_frame, width=12)
		current_source_sigma_entry=Entry(current_source_frame, width=12)
		current_source_pa_entry=Entry(current_source_frame, width=12)
		current_source_q_entry=Entry(current_source_frame, width=12)
		current_source_n_entry=Entry(current_source_frame, width=12)

		current_source_amp_entry.grid(row=0, column=column+1)
		current_source_xcen_entry.grid(row=0, column=column+3)
		current_source_ycen_entry.grid(row=0, column=column+5)
		current_source_sigma_entry.grid(row=0, column=column+7)
		current_source_pa_entry.grid(row=0, column=column+9)
		current_source_q_entry.grid(row=0, column=column+11)
		current_source_n_entry.grid(row=0, column=column+13)

		current_source_amp_entry.delete(0, END)
		current_source_xcen_entry.delete(0, END)
		current_source_ycen_entry.delete(0, END)
		current_source_sigma_entry.delete(0, END)
		current_source_pa_entry.delete(0, END)
		current_source_q_entry.delete(0, END)
		current_source_n_entry.delete(0, END)

		current_source_amp_entry.insert(0, "value=1., min=0.0")
		current_source_xcen_entry.insert(0, "value=0.")
		current_source_ycen_entry.insert(0, "value=0.")
		current_source_sigma_entry.insert(0, "value=0.05, min=0.000001")
		current_source_pa_entry.insert(0, "value=20.0, min=1, max=180.0")
		current_source_q_entry.insert(0, "value=0.6, min=0.1, max=1.0")
		current_source_n_entry.insert(0, "value=2, min=0.01")

		spar_entry_i={'source_amp_entry': current_source_amp_entry, 'source_xcen_entry': current_source_xcen_entry, 'source_ycen_entry': current_source_ycen_entry, 'source_sigma_entry': current_source_sigma_entry, 'source_pa_entry': current_source_pa_entry, 'source_q_entry': current_source_q_entry, 'source_n_entry': current_source_n_entry}
		spar_entry_list.append(spar_entry_i)
	elif sourcemodel[-1] == 'pix':
		Label(current_source_frame, text="lambda", width=12).grid(row=0, column=column)
		Label(current_source_frame, text="reg. scheme", width=12).grid(row=0, column=column+2)
		current_source_lambda_entry=Entry(current_source_frame, width=12)
		current_source_reg_scheme_entry=Entry(current_source_frame, width=12)
		current_source_lambda_entry.grid(row=0, column=column+1)
		current_source_reg_scheme_entry.grid(row=0, column=column+3)
		current_source_lambda_entry.delete(0, END)
		current_source_reg_scheme_entry.delete(0, END)
		current_source_lambda_entry.insert(0, "value=10.0, min=0., vary=0")
		current_source_reg_scheme_entry.insert(0, "value=0, vary=0")
		spar_entry_i={'lambda_entry': current_source_lambda_entry, 'reg_scheme_entry': current_source_reg_scheme_entry}
		spar_entry_list.append(spar_entry_i)
		
def drop_a_source(row):
	global n_source, current_row_source
	print 'Dropping source %i' %(row)
	for widget in source_frame_list[row-1].winfo_children():
		widget.destroy()
	n_source-=1
	source_frame_list.pop(row-1)
	drop_source_button_list.pop(row-1)
	sourcemodel.pop(row-1)
	spar_entry_list.pop(row-1)
	drop_source_number=row
	for i_source in range(drop_source_number-1, n_source):
		source_frame_list[i_source].grid_forget()
		this_row=(i_source+1)
		source_frame_list[i_source].grid(row=this_row, column=0, rowspan=1, sticky='ew')
		drop_source_button_list[i_source].config(text="Drop", command=lambda this_row=this_row: drop_a_source(this_row))
	if n_source==0:
		source_model_label.destroy()
		current_row_source=current_row_source-1
	current_row_source=current_row_source-1

def residual(pars, x, y, PSF, n_lens, n_source, lensmodel, sourcemodel, data=None, I_invvar=None):
	global I_model
	vals=pars.valuesdict()
	p=vals.values()
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]

	spar=p[0:spar_index]
	lpar=p[spar_index:spar_index+lpar_index]
	ppar=p[spar_index+lpar_index:]
	I_model=lf.model(x, y, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)
	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=lf.phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

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

def residual_pix(pars, x, y, n_lens, lensmodel, lens_mask_2d, used_pixels, used_pixels_mask, data=None, I_invvar=None, PSF=np.array([0]), PSF_csr=np.array([0])):
	global I_model
	vals=pars.valuesdict()
	p=vals.values()
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]
	spar=p[0:spar_index]
	lpar=p[spar_index:spar_index+lpar_index]
	ppar=p[spar_index+lpar_index:]
	lam=spar[0]
	reg_scheme=spar[1]
	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=lf.phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

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

def residual_phot(pars, x, y, PSF, n_lens, n_phot, n_source, lensmodel, photmodel, sourcemodel, data=None, I_invvar=None):
	global I_model, I_phot
	vals=pars.valuesdict()
	p=vals.values()
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	ppar_index=0
	for i_phot in np.arange(n_phot):
		ppar_index+=n_phot_par_dict[photmodel[i_phot]]

	spar=p[0:spar_index]
	lpar=p[spar_index:-ppar_index]
	ppar=p[-ppar_index:]

	I_model=lf.model(x, y, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)
	
	I_phot=lf.phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	if data is None:
 		return (I_model+I_phot).reshape(-1)
	if I_invvar is None:
		return (data-I_model-I_phot).reshape(-1)
	else:
		return ((data-I_model-I_phot)*np.sqrt(I_invvar)).reshape(-1)

def residual_phot_pix(pars, x, y, n_lens, n_phot, lensmodel, photmodel, lens_mask_2d, used_pixels, used_pixels_mask, data=None, I_invvar=None, PSF_csr=np.array([0]), PSF=np.array([0])):
	global I_model, I_phot
	vals=pars.valuesdict()
	p=vals.values()
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	ppar_index=0
	for i_phot in np.arange(n_phot):
		ppar_index+=n_phot_par_dict[photmodel[i_phot]]

	spar=p[0:spar_index]
	lpar=p[spar_index:-ppar_index]
	ppar=p[-ppar_index:]
	lam=spar[0]
	reg_scheme=spar[1]

	I_phot=lf.phot_model(x, y, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	output=lf.pix_source(data-I_phot, I_invvar, lens_mask_2d, used_pixels, used_pixels_mask, x, y, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=0)
	I_model=output['fit']

	if data is None:
 		return (I_model+I_phot).reshape(-1)
	if I_invvar is None:
		return (data-I_model-I_phot).reshape(-1)
	else:
		return ((data-I_model-I_phot)*np.sqrt(I_invvar)).reshape(-1)

def initial_phot_fit(pars, x, y, PSF, n_phot, photmodel, data=None, I_invvar=None):
	global I_phot
	vals=pars.valuesdict()
	p=vals.values()
	I_phot=lf.phot_model(x, y, p, photmodel, PSF=PSF, n_phot=n_phot)
	if data is None:
 		return np.sqrt(I_phot).reshape(-1)
	if I_invvar is None:
		return (data-I_phot).reshape(-1)
	else:
		return ((data-I_phot)*np.sqrt(I_invvar)).reshape(-1)
	
def messg(params, iter, resid, *args, **kws):
	print '# of Iteration: %i' %iter
	print 'Chi2= %10.4f' %np.sum(resid**2.0)
	print fit_report(params, show_correl=False)

def params_setting(var_string):
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

def update_parameter_value(dictionary, key, value):
	var_string=dictionary[key].get()
	try:
		float(var_string)
	except: # the entry is a string
		low=string.find(var_string, 'value=')
		if low < 0:
			pass # "value" is missing, no need to update its value
		else:
			high=string.find(var_string, ',', low)
			if high > 0:
				dictionary[key].delete(low+6, high)
				dictionary[key].insert(low+6, value)
			else:
				dictionary[key].delete(low+6, END)
				dictionary[key].insert(low+6, value)
	else: # the entry is a pure number, which should be the value of this parameter
		dictionary[key].delete(0, END)
		dictionary[key].insert(0, value)

def phot_fit():
	global I_data1, x1, y1, I_invvar1, I_phot, ppar, ppar_err
	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	leastsq_kws={'maxfev':8000}
	fit_params=Parameters()
	if (n_phot > 0 and not ('bspline' in photmodel)):
		for i_phot in np.arange(n_phot):
			if photmodel[i_phot] == 'sersic':
				setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])

				setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())					
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
					
				setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_sigma_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
					
				setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
					
				setting=params_setting((ppar_entry_list[i_phot]['phot_n_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
					
				setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			elif photmodel[i_phot] == 'csersic':
				setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_sigma_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
					
				setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_n_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_rc_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_rc', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_alpha_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_alpha', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_gamma_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_gamma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])	
			elif photmodel[i_phot] == 'hernquist':		
				setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_rs_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_rs', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
				setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
				fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])

		out=minimize(initial_phot_fit, fit_params, args=(x1, y1, PSF, n_phot, photmodel), kws={'data':I_data1, 'I_invvar':I_invvar1*gmask1}, iter_cb=messg, **leastsq_kws)	
	else:
		pass
	print "==============================================================="
	print "==============================================================="
	print "Fitting finished"
	print (fit_report(fit_params))
	print "==============================================================="
	fit_par=fit_params.valuesdict()
	ppar=fit_par.values()

	ppar_err=np.zeros_like(ppar)
	current_phot_index=0
	for i_phot in np.arange(n_phot):
		if photmodel[i_phot] == 'sersic':
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
		elif photmodel[i_phot] == 'csersic':
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
		elif photmodel[i_phot] == 'hernquist':
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
		elif photmodel[i_phot] == 'bspline':
			ppar_err=[]

	current_phot_index=0
	for i_phot in range(n_phot):
		ppar_entry_i=ppar_entry_list[i_phot]
		if photmodel[i_phot] == 'sersic':
			update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_sigma_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_n_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
			current_phot_index+=1
		elif photmodel[i_phot] == 'csersic':
			update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_sigma_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_n_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_rc_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_alpha_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_gamma_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
			current_phot_index+=1
		elif photmodel[i_phot] == 'hernquist':
			update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_rs_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
			current_phot_index+=1
			update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
			current_phot_index+=1

	I_phot=lf.phot_model(x1, y1, ppar, photmodel, PSF=PSF, n_phot=n_phot)

def phot_fit_mcmc():
	global frame_mcmc, phot_priorpars_entry_list, phot_priorpars
	popup_mcmc=Toplevel()
	popup_mcmc.wm_title("MCMC Fitting")
	popup_mcmc.protocol("WM_DELETE_WINDOW", lambda: popup_mcmc.destroy())
	frame_mcmc=Frame(popup_mcmc)
	for i in range(30):
		frame_mcmc.grid_rowconfigure(i, weight=1)
		frame_mcmc.grid_columnconfigure(i, weight=1)
	frame_mcmc.pack(fill=BOTH, expand=YES)
	
	n_par_phot=0
	phot_priorpars_entry_list=[]
	current_row=1
	current_column=0
	for i_phot in range(n_phot):
		if i_phot == 0:
			Label(frame_mcmc, text="phot model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_phot+=n_phot_par_dict[photmodel[i_phot]]
		add_a_phot_entry_mcmc(photmodel[i_phot], current_row, current_column)
		current_row+=n_phot_par_dict[photmodel[i_phot]]
	current_column+=4
	phot_priorpars=np.zeros((n_par_phot, 3))

	priorpars=phot_priorpars
	Button(frame_mcmc, text='Run', width=12, command=run_phot_fit_mcmc).grid(row=30, column=current_column+1, pady=4)
	Label(frame_mcmc, text="* 0-uniform prior; 1-gaussian prior").grid(row=current_row, column=current_column+3)
	
def init_pos(nwalkers, priorpars):
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

def run_phot_fit_mcmc():
	get_phot_priors()

	global I_data1, x1, y1, I_invvar1, I_phot
	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]

	ndim=(phot_priorpars.shape)[0]
	nwalkers=100
	print 'Setting up the sampler...............'
	sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob_phot, args=(x1, y1, PSF, n_phot, photmodel, I_data1, I_invvar1, phot_priorpars)) 
	print 'Setting up the sampler...Completed...'
	
	pos=init_pos(nwalkers, phot_priorpars)
	#pos=[p_phot0+np.random.randn(ndim) for i in range(nwalkers)]
	print 'Running MCMC.........................'
	# run 1000 steps as a burn-in
	pos, prob, state=sampler.run_mcmc(pos, 100)
	# reset the chain to remove the burn-in samples
	sampler.reset()
	print 'Running MCMC....burn-in completed....'
	# starting from the final position in the burn-in chain, sample for 5000
	sampler.run_mcmc(pos, 500, rstate0=state)
	print 'Running MCMC........Completed........'
	
	samples=sampler.chain[:, :, :].reshape((-1, ndim))
	
	res=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))	
	p_mcmc=[]
	for i in range(ndim):
		p_mcmc.append(res[i][0])
	
	print p_mcmc
	I_phot=lf.phot_model(x1, y1, p_mcmc, photmodel, PSF=PSF, n_phot=n_phot)
	update_phot_parameters(p_mcmc, n_phot, photmodel)
#	import triangle
#	fig=triangle.corner(samples[0:14000, :], quantiles=[0.16, 0.5, 0.84])
#	plt.show()

def lfit():
	global I_image1, I_fit, I_bspline, I_phot, I_source, I_source_err, I_invvar1, gmask1, jmask1, x1, y1, source_xbase, source_ybase, n_lens, n_phot, n_source, lensmodel, photmodel, sourcemodel, spar, lpar, ppar, spar_err, lpar_err, ppar_err, out, dpix_source, dof, chi2

	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_image1=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_bspline=I_lens[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	jmask1=jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]

	leastsq_kws={'maxfev':8000}

	if len(sourcemodel) ==1 and sourcemodel[0] == 'pix':
		lens_mask_2d=(1-gmask1.astype(int))
		lens_mask_1d=np.where(lens_mask_2d == 1)
		(used_pixels, used_pixels_mask)=lf.select_used_pixels(I_image1, lens_mask_1d, lens_mask_1d[0].shape[0]/2, scheme=0)
		used_pixels_2d=np.zeros_like(x1, dtype=int)
		used_pixels_2d[used_pixels]=1
		PSF_csr=lf.generate_psf_csr(lens_mask_1d, x1.shape[0], x1.shape[1], tpsf)

	fit_params=Parameters()
	for i_source in np.arange(n_source):
		if sourcemodel[i_source] == 'gaussian' or sourcemodel[i_source] == 'sersic':
			setting=params_setting((spar_entry_list[i_source]['source_amp_entry']).get())										
			fit_params.add('source'+str(i_source)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_xcen_entry']).get())										
			fit_params.add('source'+str(i_source)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_ycen_entry']).get())										
			fit_params.add('source'+str(i_source)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_sigma_entry']).get())										
			fit_params.add('source'+str(i_source)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_pa_entry']).get())										
			fit_params.add('source'+str(i_source)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_q_entry']).get())										
			fit_params.add('source'+str(i_source)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['source_n_entry']).get())										
			fit_params.add('source'+str(i_source)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif sourcemodel[i_source] == 'pix':
			setting=params_setting((spar_entry_list[i_source]['lambda_entry']).get())
			fit_params.add('lam', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((spar_entry_list[i_source]['reg_scheme_entry']).get())
			fit_params.add('reg_scheme', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
	for i_lens in np.arange(n_lens):
		if lensmodel[i_lens]=='ext. shear':
			setting=params_setting((lpar_entry_list[i_lens]['shear_amp_entry']).get())										
			fit_params.add('shear_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['shear_pa_entry']).get())										
			fit_params.add('shear_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif (lensmodel[i_lens]=='sie' or lensmodel[i_lens]=='sple' or lensmodel[i_lens]=='softie'):
			setting=params_setting((lpar_entry_list[i_lens]['b_SIE_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_bsie', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_pa_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_q_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['gamma_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_gamma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif lensmodel[i_lens]=='pm':
			setting=params_setting((lpar_entry_list[i_lens]['theta_e_entry']).get())
			fit_params.add('lens'+str(i_lens)+'_theta_e', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif (lensmodel[i_lens]=='snfw'):
			setting=params_setting((lpar_entry_list[i_lens]['kappa_s_entry']).get())
			fit_params.add('lens'+str(i_lens)+'_kappa_s', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_xcen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_ycen_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_pa_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['lens_q_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			setting=params_setting((lpar_entry_list[i_lens]['r_s_entry']).get())	
			fit_params.add('lens'+str(i_lens)+'_r_s', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])

	for i_phot in np.arange(n_phot):
		if photmodel[i_phot] == 'sersic':
			setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())
			fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])

			setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())					
			fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
			setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_sigma_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
			setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
			setting=params_setting((ppar_entry_list[i_phot]['phot_n_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
			setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif photmodel[i_phot] == 'csersic':
			setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_sigma_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_sigma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
				
			setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_n_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_n', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_rc_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_rc', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_alpha_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_alpha', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_gamma_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_gamma', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])	
		elif photmodel[i_phot] == 'hernquist':		
			setting=params_setting((ppar_entry_list[i_phot]['phot_amp_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_amp', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_xcen_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_xcen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_ycen_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_ycen', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_rs_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_rs', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_pa_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_pa', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_q_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_q', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
			
			setting=params_setting((ppar_entry_list[i_phot]['phot_I0_entry']).get())										
			fit_params.add('phot'+str(i_phot)+'_I0', value=setting['value'], vary=setting['vary'], min=setting['min'], max=setting['max'], expr=setting['expr'])
		elif photmodel[i_phot] == 'bspline':
			pass

	if (len(sourcemodel)==1 and sourcemodel[0] == 'pix'):
		out=minimize(residual_pix, fit_params, args=(x1, y1, n_lens, lensmodel, lens_mask_2d, used_pixels, used_pixels_mask), kws={'data':I_data1, 'I_invvar':I_invvar1*jmask1, 'PSF':PSF, 'PSF_csr':PSF_csr}, iter_cb=messg, **leastsq_kws)
	else:	
#		mini=lmfit.Minimizer(residual, fit_params, fcn_args=(x1, y1, PSF, n_lens, n_source, lensmodel, sourcemodel), fcn_kws={'data':I_data1, 'I_invvar':I_invvar1}, iter_cb=messg, **leastsq_kws)
#		out=mini.minimize(method='leastsq')
#		print 'yes'
		out=minimize(residual, fit_params, args=(x1, y1, PSF, n_lens, n_source, lensmodel, sourcemodel), kws={'data':I_data1, 'I_invvar':I_invvar1*jmask1}, iter_cb=messg, **leastsq_kws)
	dof=out.nfree-(np.where(I_invvar1*jmask1==0.0))[0].size
	chi2=out.chisqr
	
	print "==============================================================="
	print "==============================================================="
	print "Fitting finished"
	print (fit_report(fit_params))
#	ci, trace=conf_interval(out, p_names=['lens0_bsie', 'lens0_xcen', 'lens0_ycen', 'lens1_bsie', 'lens1_xcen', 'lens1_ycen'], sigmas=[0.68,0.95], trace=True, verbose=False)
#	printfuncs.report_ci(ci)
	print "==============================================================="

	fit_par=fit_params.valuesdict()
	pars=fit_par.values()
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]

	spar=pars[0:spar_index]
	lpar=pars[spar_index:spar_index+lpar_index]
	ppar=pars[spar_index+lpar_index:]

	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=lf.phot_model(x1, y1, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	if (len(sourcemodel) == 1 and sourcemodel[0] == 'pix'):
		lam=spar[0]
		reg_scheme=spar[1]
		output=lf.pix_source(I_data1-I_phot, I_invvar1*jmask1, lens_mask_2d, used_pixels, used_pixels_mask, x1, y1, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=1)
		I_fit=output['fit']
		Vec_S=output['solution']
		Cov_matrix=output['Covariance']
		Vec_S_err=np.sqrt(np.diagonal(Cov_matrix))
	else:
		I_fit=lf.model(x1, y1, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)

	spar_err=np.zeros_like(spar)
	current_source_index=0
	for i_source in np.arange(n_source):
		if sourcemodel[i_source] == 'sersic' or sourcemodel[i_source] == 'gaussian':
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
		elif sourcemodel[i_source] == 'pix':
			spar_err[current_source_index]=fit_params['lam'].stderr
			current_source_index+=1
			spar_err[current_source_index]=fit_params['reg_scheme'].stderr
			current_source_index+=1

	lpar_err=np.zeros_like(lpar)
	current_lens_index=0
	for i_lens in np.arange(n_lens):
		if lensmodel[i_lens] == 'sie' or lensmodel[i_lens] == 'sple' or lensmodel[i_lens] == 'softie':
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
		elif lensmodel[i_lens] == 'pm':
			lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_theta_e'].stderr
			current_lens_index+=1
			lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_xcen'].stderr
			current_lens_index+=1
			lpar_err[current_lens_index]=fit_params['lens'+str(i_lens)+'_ycen'].stderr
			current_lens_index+=1
		elif lensmodel[i_lens] == 'ext. shear':
			lpar_err[current_lens_index]=fit_params['shear_amp'].stderr
			current_lens_index+=1
			lpar_err[current_lens_index]=fit_params['shear_pa'].stderr
			current_lens_index+=1	
		elif lensmodel[i_lens] == 'snfw':
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
	for i_phot in np.arange(n_phot):
		if photmodel[i_phot] == 'sersic':
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
		elif photmodel[i_phot] == 'csersic':
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
		elif photmodel[i_phot] == 'hernquist':
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
		elif photmodel[i_phot] == 'bspline':
			ppar_err=[]

	n_x_source=10*100+1
	n_y_source=10*100+1
	dpix_source=0.1*dpix
	source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2)*dpix_source
	source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2, np.ones(n_x_source))*dpix_source
	if len(sourcemodel) ==1 and sourcemodel[0] == 'pix':
		theta_x=x1[lens_mask_1d].flatten()
		theta_y=y1[lens_mask_1d].flatten()
		(alpha_x, alpha_y)=lf.def_total(x1, y1, lpar, lensmodel, n_lens=n_lens)
		alpha_x_1d=alpha_x[lens_mask_1d].flatten()
		alpha_y_1d=alpha_y[lens_mask_1d].flatten()
		
		tri=lf.construct_source_grid(theta_x, theta_y, alpha_x_1d, alpha_y_1d, used_pixels_mask)
		I_source=griddata(tri.points, Vec_S, (source_xbase, source_ybase), method='linear', fill_value=0.0)
		I_source_err=griddata(tri.points, Vec_S_err, (source_xbase, source_ybase), method='linear', fill_value=0.0)
	else:
		I_source=0.0*source_xbase
		current_source_index=0
		for i_source in np.arange(n_source):
			spar_i=spar[current_source_index:current_source_index+n_source_par_dict[sourcemodel[i_source]]]
			I_source+=lf.sersic_2d(source_xbase, source_ybase, spar_i)
			current_source_index+=n_source_par_dict[sourcemodel[i_source]]

#	update the parameter entries
	update_phot_parameters(ppar, n_phot, photmodel)
	update_lens_parameters(lpar, n_lens, lensmodel)
	update_source_parameters(spar, n_source, sourcemodel)

def update_lens_parameters(lpar, n_lens, lensmodel):
	current_lens_index=0
	for i_lens in range(n_lens):
		lpar_entry_i=lpar_entry_list[i_lens]
		if lensmodel[i_lens] == 'sie' or lensmodel[i_lens] == 'sple' or lensmodel[i_lens] == 'softie':
			update_parameter_value(lpar_entry_i, 'b_SIE_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_xcen_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_ycen_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_pa_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_q_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'gamma_entry', lpar[current_lens_index])
			current_lens_index+=1
		elif lensmodel[i_lens] == 'pm':
			update_parameter_value(lpar_entry_i, 'theta_e_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_xcen_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_ycen_entry', lpar[current_lens_index])
			current_lens_index+=1
		elif lensmodel[i_lens] == 'ext. shear':
			update_parameter_value(lpar_entry_i, 'shear_amp_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'shear_pa_entry', lpar[current_lens_index])
			current_lens_index+=1
		elif lensmodel[i_lens] == 'snfw':
			update_parameter_value(lpar_entry_i, 'kappa_s_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_xcen_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_ycen_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_pa_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'lens_q_entry', lpar[current_lens_index])
			current_lens_index+=1
			update_parameter_value(lpar_entry_i, 'r_s_entry', lpar[current_lens_index])
			current_lens_index+=1

def update_source_parameters(spar, n_source, sourcemodel):
	current_source_index=0
	for i_source in range(n_source):
		spar_entry_i=spar_entry_list[i_source]
		if sourcemodel[i_source] == 'sersic' or sourcemodel[i_source] == 'gaussian':
			update_parameter_value(spar_entry_i, 'source_amp_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_xcen_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_ycen_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_sigma_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_pa_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_q_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'source_n_entry', spar[current_source_index])
			current_source_index+=1
		elif sourcemodel[i_source] == 'pix':
			update_parameter_value(spar_entry_i, 'lambda_entry', spar[current_source_index])
			current_source_index+=1
			update_parameter_value(spar_entry_i, 'reg_scheme_entry', spar[current_source_index])
			current_source_index+=1

def update_phot_parameters(ppar, n_phot, photmodel):
	current_phot_index=0
	for i_phot in range(n_phot):
		if ((n_phot == 0) or (n_phot == 1 and photmodel[0] == 'bspline')): 
			continue
		else:
			ppar_entry_i=ppar_entry_list[i_phot]
			if photmodel[i_phot] == 'sersic':
				update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_sigma_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_n_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
				current_phot_index+=1
			elif photmodel[i_phot] == 'csersic':
				update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_sigma_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_n_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_rc_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_alpha_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_gamma_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
				current_phot_index+=1
			elif photmodel[i_phot] == 'hernquist':
				update_parameter_value(ppar_entry_i, 'phot_amp_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_xcen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_ycen_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_rs_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_pa_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_q_entry', ppar[current_phot_index])
				current_phot_index+=1
				update_parameter_value(ppar_entry_i, 'phot_I0_entry', ppar[current_phot_index])
				current_phot_index+=1

def add_a_phot_entry_mcmc(model, row, column):
	global phot_priorpars_entry_list
	current_row=row
	current_column=column
	if model == 'sersic':
		Label(frame_mcmc, text="amp", width=12).grid(row=current_row, column=current_column)
		current_phot_amp_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_amp_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_amp_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_amp_ptype_entry.delete(0, END)
		current_phot_amp_ptype_entry.insert(0, 0)	
		current_phot_amp_p1_entry.delete(0, END)
		current_phot_amp_p1_entry.insert(0, 1)	
		current_phot_amp_p2_entry.delete(0, END)
		current_phot_amp_p2_entry.insert(0, 500)	
		current_row+=1
		Label(frame_mcmc, text="xcen", width=12).grid(row=current_row, column=current_column)
		current_phot_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_xcen_ptype_entry.delete(0, END)
		current_phot_xcen_ptype_entry.insert(0, 0)	
		current_phot_xcen_p1_entry.delete(0, END)
		current_phot_xcen_p1_entry.insert(0, -0.1)	
		current_phot_xcen_p2_entry.delete(0, END)
		current_phot_xcen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="ycen", width=12).grid(row=current_row, column=current_column)
		current_phot_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_ycen_ptype_entry.delete(0, END)
		current_phot_ycen_ptype_entry.insert(0, 0)	
		current_phot_ycen_p1_entry.delete(0, END)
		current_phot_ycen_p1_entry.insert(0, -0.1)	
		current_phot_ycen_p2_entry.delete(0, END)
		current_phot_ycen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="sigma", width=12).grid(row=current_row, column=current_column)
		current_phot_sigma_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_sigma_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_sigma_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_sigma_ptype_entry.delete(0, END)
		current_phot_sigma_ptype_entry.insert(0, 0)	
		current_phot_sigma_p1_entry.delete(0, END)
		current_phot_sigma_p1_entry.insert(0, 0.1)	
		current_phot_sigma_p2_entry.delete(0, END)
		current_phot_sigma_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="pa", width=12).grid(row=current_row, column=current_column)
		current_phot_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_pa_ptype_entry.delete(0, END)
		current_phot_pa_ptype_entry.insert(0, 0)	
		current_phot_pa_p1_entry.delete(0, END)
		current_phot_pa_p1_entry.insert(0, 1.0)	
		current_phot_pa_p2_entry.delete(0, END)
		current_phot_pa_p2_entry.insert(0, 180.0)	
		current_row+=1
		Label(frame_mcmc, text="q", width=12).grid(row=current_row, column=current_column)
		current_phot_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_q_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_q_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_q_ptype_entry.delete(0, END)
		current_phot_q_ptype_entry.insert(0, 0)	
		current_phot_q_p1_entry.delete(0, END)
		current_phot_q_p1_entry.insert(0, 0.1)	
		current_phot_q_p2_entry.delete(0, END)
		current_phot_q_p2_entry.insert(0, 1.0)	
		current_row+=1
		Label(frame_mcmc, text="n", width=12).grid(row=current_row, column=current_column)
		current_phot_n_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_n_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_n_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_n_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_n_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_n_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_n_ptype_entry.delete(0, END)
		current_phot_n_ptype_entry.insert(0, 0)	
		current_phot_n_p1_entry.delete(0, END)
		current_phot_n_p1_entry.insert(0, 0.1)	
		current_phot_n_p2_entry.delete(0, END)
		current_phot_n_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="bg", width=12).grid(row=current_row, column=current_column)
		current_phot_bg_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_bg_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_bg_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_bg_ptype_entry.delete(0, END)
		current_phot_bg_ptype_entry.insert(0, 0)	
		current_phot_bg_p1_entry.delete(0, END)
		current_phot_bg_p1_entry.insert(0, -0.02)	
		current_phot_bg_p2_entry.delete(0, END)
		current_phot_bg_p2_entry.insert(0, 0.02)	
		current_row+=1
		phot_priorpars_entry_i={
			'phot_amp_ptype': current_phot_amp_ptype_entry, 
			'phot_amp_p1': current_phot_amp_p1_entry, 
			'phot_amp_p2': current_phot_amp_p2_entry, 
			'phot_xcen_ptype': current_phot_xcen_ptype_entry, 
			'phot_xcen_p1': current_phot_xcen_p1_entry, 
			'phot_xcen_p2': current_phot_xcen_p2_entry, 
			'phot_ycen_ptype': current_phot_ycen_ptype_entry, 
			'phot_ycen_p1': current_phot_ycen_p1_entry, 
			'phot_ycen_p2': current_phot_ycen_p2_entry, 
			'phot_sigma_ptype': current_phot_sigma_ptype_entry, 
			'phot_sigma_p1': current_phot_sigma_p1_entry, 
			'phot_sigma_p2': current_phot_sigma_p2_entry, 
			'phot_pa_ptype': current_phot_pa_ptype_entry, 
			'phot_pa_p1': current_phot_pa_p1_entry, 
			'phot_pa_p2': current_phot_pa_p2_entry, 
			'phot_q_ptype': current_phot_q_ptype_entry, 
			'phot_q_p1': current_phot_q_p1_entry, 
			'phot_q_p2': current_phot_q_p2_entry, 
			'phot_n_ptype': current_phot_n_ptype_entry, 
			'phot_n_p1': current_phot_n_p1_entry, 
			'phot_n_p2': current_phot_n_p2_entry, 
			'phot_bg_ptype': current_phot_bg_ptype_entry, 
			'phot_bg_p1': current_phot_bg_p1_entry, 
			'phot_bg_p2': current_phot_bg_p2_entry}
	elif model == 'csersic':
		Label(frame_mcmc, text="amp", width=12).grid(row=current_row, column=current_column)
		current_phot_amp_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_amp_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_amp_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_amp_ptype_entry.delete(0, END)
		current_phot_amp_ptype_entry.insert(0, 0)	
		current_phot_amp_p1_entry.delete(0, END)
		current_phot_amp_p1_entry.insert(0, 0)	
		current_phot_amp_p2_entry.delete(0, END)
		current_phot_amp_p2_entry.insert(0, 500)	
		current_row+=1
		Label(frame_mcmc, text="xcen", width=12).grid(row=current_row, column=current_column)
		current_phot_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_xcen_ptype_entry.delete(0, END)
		current_phot_xcen_ptype_entry.insert(0, 0)	
		current_phot_xcen_p1_entry.delete(0, END)
		current_phot_xcen_p1_entry.insert(0, -0.1)	
		current_phot_xcen_p2_entry.delete(0, END)
		current_phot_xcen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="ycen", width=12).grid(row=current_row, column=current_column)
		current_phot_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_ycen_ptype_entry.delete(0, END)
		current_phot_ycen_ptype_entry.insert(0, 0)	
		current_phot_ycen_p1_entry.delete(0, END)
		current_phot_ycen_p1_entry.insert(0, -0.1)	
		current_phot_ycen_p2_entry.delete(0, END)
		current_phot_ycen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="sigma", width=12).grid(row=current_row, column=current_column)
		current_phot_sigma_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_sigma_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_sigma_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_sigma_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_sigma_ptype_entry.delete(0, END)
		current_phot_sigma_ptype_entry.insert(0, 0)	
		current_phot_sigma_p1_entry.delete(0, END)
		current_phot_sigma_p1_entry.insert(0, 0.0)	
		current_phot_sigma_p2_entry.delete(0, END)
		current_phot_sigma_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="pa", width=12).grid(row=current_row, column=current_column)
		current_phot_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_pa_ptype_entry.delete(0, END)
		current_phot_pa_ptype_entry.insert(0, 0)	
		current_phot_pa_p1_entry.delete(0, END)
		current_phot_pa_p1_entry.insert(0, 0.0)	
		current_phot_pa_p2_entry.delete(0, END)
		current_phot_pa_p2_entry.insert(0, 180.0)	
		current_row+=1
		Label(frame_mcmc, text="q", width=12).grid(row=current_row, column=current_column)
		current_phot_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_q_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_q_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_q_ptype_entry.delete(0, END)
		current_phot_q_ptype_entry.insert(0, 0)	
		current_phot_q_p1_entry.delete(0, END)
		current_phot_q_p1_entry.insert(0, 0.0)	
		current_phot_q_p2_entry.delete(0, END)
		current_phot_q_p2_entry.insert(0, 1.0)	
		current_row+=1
		Label(frame_mcmc, text="n", width=12).grid(row=current_row, column=current_column)
		current_phot_n_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_n_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_n_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_n_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_n_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_n_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_n_ptype_entry.delete(0, END)
		current_phot_n_ptype_entry.insert(0, 0)	
		current_phot_n_p1_entry.delete(0, END)
		current_phot_n_p1_entry.insert(0, 0.0)	
		current_phot_n_p2_entry.delete(0, END)
		current_phot_n_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="rc", width=12).grid(row=current_row, column=current_column)
		current_phot_rc_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_rc_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_rc_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_rc_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_rc_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_rc_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_rc_ptype_entry.delete(0, END)
		current_phot_rc_ptype_entry.insert(0, 0)	
		current_phot_rc_p1_entry.delete(0, END)
		current_phot_rc_p1_entry.insert(0, 0.0)	
		current_phot_rc_p2_entry.delete(0, END)
		current_phot_rc_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="alpha", width=12).grid(row=current_row, column=current_column)
		current_phot_alpha_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_alpha_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_alpha_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_alpha_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_alpha_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_alpha_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_alpha_ptype_entry.delete(0, END)
		current_phot_alpha_ptype_entry.insert(0, 0)	
		current_phot_alpha_p1_entry.delete(0, END)
		current_phot_alpha_p1_entry.insert(0, 0.0)	
		current_phot_alpha_p2_entry.delete(0, END)
		current_phot_alpha_p2_entry.insert(0, 200.0)	
		current_row+=1
		Label(frame_mcmc, text="gamma", width=12).grid(row=current_row, column=current_column)
		current_phot_gamma_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_gamma_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_gamma_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_gamma_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_gamma_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_gamma_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_gamma_ptype_entry.delete(0, END)
		current_phot_gamma_ptype_entry.insert(0, 0)	
		current_phot_gamma_p1_entry.delete(0, END)
		current_phot_gamma_p1_entry.insert(0, 0.0)	
		current_phot_gamma_p2_entry.delete(0, END)
		current_phot_gamma_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="bg", width=12).grid(row=current_row, column=current_column)
		current_phot_bg_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_bg_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_bg_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_bg_ptype_entry.delete(0, END)
		current_phot_bg_ptype_entry.insert(0, 0)	
		current_phot_bg_p1_entry.delete(0, END)
		current_phot_bg_p1_entry.insert(0, -0.02)	
		current_phot_bg_p2_entry.delete(0, END)
		current_phot_bg_p2_entry.insert(0, 0.02)	
		current_row+=1
		phot_priorpars_entry_i={
			'phot_amp_ptype': current_phot_amp_ptype_entry, 
			'phot_amp_p1': current_phot_amp_p1_entry, 
			'phot_amp_p2': current_phot_amp_p2_entry, 
			'phot_xcen_ptype': current_phot_xcen_ptype_entry, 
			'phot_xcen_p1': current_phot_xcen_p1_entry, 
			'phot_xcen_p2': current_phot_xcen_p2_entry, 
			'phot_ycen_ptype': current_phot_ycen_ptype_entry, 
			'phot_ycen_p1': current_phot_ycen_p1_entry, 
			'phot_ycen_p2': current_phot_ycen_p2_entry, 
			'phot_sigma_ptype': current_phot_sigma_ptype_entry, 
			'phot_sigma_p1': current_phot_sigma_p1_entry, 
			'phot_sigma_p2': current_phot_sigma_p2_entry, 
			'phot_pa_ptype': current_phot_pa_ptype_entry, 
			'phot_pa_p1': current_phot_pa_p1_entry, 
			'phot_pa_p2': current_phot_pa_p2_entry, 
			'phot_q_ptype': current_phot_q_ptype_entry, 
			'phot_q_p1': current_phot_q_p1_entry, 
			'phot_q_p2': current_phot_q_p2_entry, 
			'phot_n_ptype': current_phot_n_ptype_entry, 
			'phot_n_p1': current_phot_n_p1_entry, 
			'phot_n_p2': current_phot_n_p2_entry, 
			'phot_rc_ptype': current_phot_rc_ptype_entry, 
			'phot_rc_p1': current_phot_rc_p1_entry, 
			'phot_rc_p2': current_phot_rc_p2_entry, 
			'phot_alpha_ptype': current_phot_alpha_ptype_entry, 
			'phot_alpha_p1': current_phot_alpha_p1_entry, 
			'phot_alpha_p2': current_phot_alpha_p2_entry, 
			'phot_gamma_ptype': current_phot_gamma_ptype_entry, 
			'phot_gamma_p1': current_phot_gamma_p1_entry, 
			'phot_gamma_p2': current_phot_gamma_p2_entry, 
			'phot_bg_ptype': current_phot_bg_ptype_entry, 
			'phot_bg_p1': current_phot_bg_p1_entry, 
			'phot_bg_p2': current_phot_bg_p2_entry}
	elif model == 'hernquist':
		Label(frame_mcmc, text="amp", width=12).grid(row=current_row, column=current_column)
		current_phot_amp_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_amp_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_amp_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_amp_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_amp_ptype_entry.delete(0, END)
		current_phot_amp_ptype_entry.insert(0, 0)	
		current_phot_amp_p1_entry.delete(0, END)
		current_phot_amp_p1_entry.insert(0, 0)	
		current_phot_amp_p2_entry.delete(0, END)
		current_phot_amp_p2_entry.insert(0, 500)	
		current_row+=1
		Label(frame_mcmc, text="xcen", width=12).grid(row=current_row, column=current_column)
		current_phot_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_xcen_ptype_entry.delete(0, END)
		current_phot_xcen_ptype_entry.insert(0, 0)	
		current_phot_xcen_p1_entry.delete(0, END)
		current_phot_xcen_p1_entry.insert(0, -0.1)	
		current_phot_xcen_p2_entry.delete(0, END)
		current_phot_xcen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="ycen", width=12).grid(row=current_row, column=current_column)
		current_phot_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_ycen_ptype_entry.delete(0, END)
		current_phot_ycen_ptype_entry.insert(0, 0)	
		current_phot_ycen_p1_entry.delete(0, END)
		current_phot_ycen_p1_entry.insert(0, -0.1)	
		current_phot_ycen_p2_entry.delete(0, END)
		current_phot_ycen_p2_entry.insert(0, 0.1)	
		current_row+=1
		Label(frame_mcmc, text="rs", width=12).grid(row=current_row, column=current_column)
		current_phot_rs_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_rs_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_rs_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_rs_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_rs_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_rs_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_rs_ptype_entry.delete(0, END)
		current_phot_rs_ptype_entry.insert(0, 0)	
		current_phot_rs_p1_entry.delete(0, END)
		current_phot_rs_p1_entry.insert(0, 0.0)	
		current_phot_rs_p2_entry.delete(0, END)
		current_phot_rs_p2_entry.insert(0, 10.0)	
		current_row+=1
		Label(frame_mcmc, text="pa", width=12).grid(row=current_row, column=current_column)
		current_phot_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_pa_ptype_entry.delete(0, END)
		current_phot_pa_ptype_entry.insert(0, 0)	
		current_phot_pa_p1_entry.delete(0, END)
		current_phot_pa_p1_entry.insert(0, 0.0)	
		current_phot_pa_p2_entry.delete(0, END)
		current_phot_pa_p2_entry.insert(0, 180.0)	
		current_row+=1
		Label(frame_mcmc, text="q", width=12).grid(row=current_row, column=current_column)
		current_phot_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_q_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_q_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_q_ptype_entry.delete(0, END)
		current_phot_q_ptype_entry.insert(0, 0)	
		current_phot_q_p1_entry.delete(0, END)
		current_phot_q_p1_entry.insert(0, 0.0)	
		current_phot_q_p2_entry.delete(0, END)
		current_phot_q_p2_entry.insert(0, 1.0)	
		current_row+=1
		Label(frame_mcmc, text="bg", width=12).grid(row=current_row, column=current_column)
		current_phot_bg_ptype_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_ptype_entry.grid(row=current_row, column=current_column+1)
		current_phot_bg_p1_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p1_entry.grid(row=current_row, column=current_column+2)
		current_phot_bg_p2_entry=Entry(frame_mcmc, width=12)
		current_phot_bg_p2_entry.grid(row=current_row, column=current_column+3)
		current_phot_bg_ptype_entry.delete(0, END)
		current_phot_bg_ptype_entry.insert(0, 0)	
		current_phot_bg_p1_entry.delete(0, END)
		current_phot_bg_p1_entry.insert(0, -0.02)	
		current_phot_bg_p2_entry.delete(0, END)
		current_phot_bg_p2_entry.insert(0, 0.02)	
		current_row+=1
		phot_priorpars_entry_i={
			'phot_amp_ptype': current_phot_amp_ptype_entry, 
			'phot_amp_p1': current_phot_amp_p1_entry, 
			'phot_amp_p2': current_phot_amp_p2_entry, 
			'phot_xcen_ptype': current_phot_xcen_ptype_entry, 
			'phot_xcen_p1': current_phot_xcen_p1_entry, 
			'phot_xcen_p2': current_phot_xcen_p2_entry, 
			'phot_ycen_ptype': current_phot_ycen_ptype_entry, 
			'phot_ycen_p1': current_phot_ycen_p1_entry, 
			'phot_ycen_p2': current_phot_ycen_p2_entry, 
			'phot_rs_ptype': current_phot_rs_ptype_entry, 
			'phot_rs_p1': current_phot_rs_p1_entry, 
			'phot_rs_p2': current_phot_rs_p2_entry, 
			'phot_pa_ptype': current_phot_pa_ptype_entry, 
			'phot_pa_p1': current_phot_pa_p1_entry, 
			'phot_pa_p2': current_phot_pa_p2_entry, 
			'phot_q_ptype': current_phot_q_ptype_entry, 
			'phot_q_p1': current_phot_q_p1_entry, 
			'phot_q_p2': current_phot_q_p2_entry, 
			'phot_bg_ptype': current_phot_bg_ptype_entry, 
			'phot_bg_p1': current_phot_bg_p1_entry, 
			'phot_bg_p2': current_phot_bg_p2_entry}
	phot_priorpars_entry_list.append(phot_priorpars_entry_i)

def add_a_lens_entry_mcmc(model, row, column):
	global lens_priorpars_entry_list
	current_row=row
	current_column=column
	if model == 'sie' or model == 'sple':
		Label(frame_mcmc, text="b_SIE", width=12).grid(row=current_row, column=current_column)
		current_lens_b_SIE_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_b_SIE_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_b_SIE_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_b_SIE_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_b_SIE_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_b_SIE_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_b_SIE_ptype_entry.delete(0, END)
		current_lens_b_SIE_ptype_entry.insert(0, 0)	
		current_lens_b_SIE_p1_entry.delete(0, END)
		current_lens_b_SIE_p1_entry.insert(0, 0.5)	
		current_lens_b_SIE_p2_entry.delete(0, END)
		current_lens_b_SIE_p2_entry.insert(0, 2.0)	

		current_row+=1
		Label(frame_mcmc, text="lens_xcen", width=12).grid(row=current_row, column=current_column)
		current_lens_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_xcen_ptype_entry.delete(0, END)
		current_lens_xcen_ptype_entry.insert(0, 0)	
		current_lens_xcen_p1_entry.delete(0, END)
		current_lens_xcen_p1_entry.insert(0, -0.2)	
		current_lens_xcen_p2_entry.delete(0, END)
		current_lens_xcen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="lens_ycen", width=12).grid(row=current_row, column=current_column)
		current_lens_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_ycen_ptype_entry.delete(0, END)
		current_lens_ycen_ptype_entry.insert(0, 0)	
		current_lens_ycen_p1_entry.delete(0, END)
		current_lens_ycen_p1_entry.insert(0, -0.2)	
		current_lens_ycen_p2_entry.delete(0, END)
		current_lens_ycen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="lens_pa", width=12).grid(row=current_row, column=current_column)
		current_lens_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_pa_ptype_entry.delete(0, END)
		current_lens_pa_ptype_entry.insert(0, 0)	
		current_lens_pa_p1_entry.delete(0, END)
		current_lens_pa_p1_entry.insert(0, 1)	
		current_lens_pa_p2_entry.delete(0, END)
		current_lens_pa_p2_entry.insert(0, 180)	

		current_row+=1
		Label(frame_mcmc, text="lens_q", width=12).grid(row=current_row, column=current_column)
		current_lens_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_q_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_q_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_q_ptype_entry.delete(0, END)
		current_lens_q_ptype_entry.insert(0, 0)	
		current_lens_q_p1_entry.delete(0, END)
		current_lens_q_p1_entry.insert(0, 0.1)	
		current_lens_q_p2_entry.delete(0, END)
		current_lens_q_p2_entry.insert(0, 1.0)	

		current_row+=1
		Label(frame_mcmc, text="gamma", width=12).grid(row=current_row, column=current_column)
		current_lens_gamma_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_gamma_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_gamma_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_gamma_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_gamma_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_gamma_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_gamma_ptype_entry.delete(0, END)
		current_lens_gamma_ptype_entry.insert(0, 0)	
		current_lens_gamma_p1_entry.delete(0, END)
		current_lens_gamma_p1_entry.insert(0, 1.0)	
		current_lens_gamma_p2_entry.delete(0, END)
		current_lens_gamma_p2_entry.insert(0, 1.0)	
		
		lens_priorpars_entry_i={
			'lens_b_SIE_ptype': current_lens_b_SIE_ptype_entry, 
			'lens_b_SIE_p1': current_lens_b_SIE_p1_entry, 
			'lens_b_SIE_p2': current_lens_b_SIE_p2_entry, 
			'lens_xcen_ptype': current_lens_xcen_ptype_entry, 
			'lens_xcen_p1': current_lens_xcen_p1_entry, 
			'lens_xcen_p2': current_lens_xcen_p2_entry, 
			'lens_ycen_ptype': current_lens_ycen_ptype_entry, 
			'lens_ycen_p1': current_lens_ycen_p1_entry, 
			'lens_ycen_p2': current_lens_ycen_p2_entry, 
			'lens_pa_ptype': current_lens_pa_ptype_entry, 
			'lens_pa_p1': current_lens_pa_p1_entry, 
			'lens_pa_p2': current_lens_pa_p2_entry, 
			'lens_q_ptype': current_lens_q_ptype_entry, 
			'lens_q_p1': current_lens_q_p1_entry, 
			'lens_q_p2': current_lens_q_p2_entry, 
			'lens_gamma_ptype': current_lens_gamma_ptype_entry, 
			'lens_gamma_p1': current_lens_gamma_p1_entry, 
			'lens_gamma_p2': current_lens_gamma_p2_entry}
	elif model == 'pm':
		Label(frame_mcmc, text="theta_e", width=12).grid(row=current_row, column=current_column)
		current_lens_theta_e_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_theta_e_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_theta_e_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_theta_e_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_theta_e_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_theta_e_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_theta_e_ptype_entry.delete(0, END)
		current_lens_theta_e_ptype_entry.insert(0, 0)	
		current_lens_theta_e_p1_entry.delete(0, END)
		current_lens_theta_e_p1_entry.insert(0, 0.5)	
		current_lens_theta_e_p2_entry.delete(0, END)
		current_lens_theta_e_p2_entry.insert(0, 2.0)	

		current_row+=1
		Label(frame_mcmc, text="lens_xcen", width=12).grid(row=current_row, column=current_column)
		current_lens_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_xcen_ptype_entry.delete(0, END)
		current_lens_xcen_ptype_entry.insert(0, 0)	
		current_lens_xcen_p1_entry.delete(0, END)
		current_lens_xcen_p1_entry.insert(0, -0.2)	
		current_lens_xcen_p2_entry.delete(0, END)
		current_lens_xcen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="lens_ycen", width=12).grid(row=current_row, column=current_column)
		current_lens_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_ycen_ptype_entry.delete(0, END)
		current_lens_ycen_ptype_entry.insert(0, 0)	
		current_lens_ycen_p1_entry.delete(0, END)
		current_lens_ycen_p1_entry.insert(0, -0.2)	
		current_lens_ycen_p2_entry.delete(0, END)
		current_lens_ycen_p2_entry.insert(0, 0.2)	

		lens_priorpars_entry_i={
			'lens_theta_e_ptype': current_lens_theta_e_ptype_entry, 
			'lens_theta_e_p1': current_lens_theta_e_p1_entry, 
			'lens_theta_e_p2': current_lens_theta_e_p2_entry, 
			'lens_xcen_ptype': current_lens_xcen_ptype_entry, 
			'lens_xcen_p1': current_lens_xcen_p1_entry, 
			'lens_xcen_p2': current_lens_xcen_p2_entry, 
			'lens_ycen_ptype': current_lens_ycen_ptype_entry, 
			'lens_ycen_p1': current_lens_ycen_p1_entry, 
			'lens_ycen_p2': current_lens_ycen_p2_entry}
	elif model == 'snfw':
		Label(frame_mcmc, text="kappa_s", width=12).grid(row=current_row, column=current_column)
		current_lens_kappa_s_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_kappa_s_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_kappa_s_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_kappa_s_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_kappa_s_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_kappa_s_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_kappa_s_ptype_entry.delete(0, END)
		current_lens_kappa_s_ptype_entry.insert(0, 0)	
		current_lens_kappa_s_p1_entry.delete(0, END)
		current_lens_kappa_s_p1_entry.insert(0, 0.5)	
		current_lens_kappa_s_p2_entry.delete(0, END)
		current_lens_kappa_s_p2_entry.insert(0, 2.0)	

		current_row+=1
		Label(frame_mcmc, text="lens_xcen", width=12).grid(row=current_row, column=current_column)
		current_lens_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_xcen_ptype_entry.delete(0, END)
		current_lens_xcen_ptype_entry.insert(0, 0)	
		current_lens_xcen_p1_entry.delete(0, END)
		current_lens_xcen_p1_entry.insert(0, -0.2)	
		current_lens_xcen_p2_entry.delete(0, END)
		current_lens_xcen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="lens_ycen", width=12).grid(row=current_row, column=current_column)
		current_lens_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_ycen_ptype_entry.delete(0, END)
		current_lens_ycen_ptype_entry.insert(0, 0)	
		current_lens_ycen_p1_entry.delete(0, END)
		current_lens_ycen_p1_entry.insert(0, -0.2)	
		current_lens_ycen_p2_entry.delete(0, END)
		current_lens_ycen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="lens_pa", width=12).grid(row=current_row, column=current_column)
		current_lens_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_pa_ptype_entry.delete(0, END)
		current_lens_pa_ptype_entry.insert(0, 0)	
		current_lens_pa_p1_entry.delete(0, END)
		current_lens_pa_p1_entry.insert(0, 0.0)	
		current_lens_pa_p2_entry.delete(0, END)
		current_lens_pa_p2_entry.insert(0, 0.0)	

		current_row+=1
		Label(frame_mcmc, text="lens_q", width=12).grid(row=current_row, column=current_column)
		current_lens_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_q_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_q_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_q_ptype_entry.delete(0, END)
		current_lens_q_ptype_entry.insert(0, 0)	
		current_lens_q_p1_entry.delete(0, END)
		current_lens_q_p1_entry.insert(0, 1.0)	
		current_lens_q_p2_entry.delete(0, END)
		current_lens_q_p2_entry.insert(0, 1.0)	

		current_row+=1
		Label(frame_mcmc, text="r_s", width=12).grid(row=current_row, column=current_column)
		current_lens_r_s_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_r_s_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_r_s_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_r_s_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_r_s_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_r_s_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_r_s_ptype_entry.delete(0, END)
		current_lens_r_s_ptype_entry.insert(0, 0)	
		current_lens_r_s_p1_entry.delete(0, END)
		current_lens_r_s_p1_entry.insert(0, 1.0)	
		current_lens_r_s_p2_entry.delete(0, END)
		current_lens_r_s_p2_entry.insert(0, 50.0)	
		
		lens_priorpars_entry_i={
			'lens_kappa_s_ptype': current_lens_kappa_s_ptype_entry, 
			'lens_kappa_s_p1': current_lens_kappa_s_p1_entry, 
			'lens_kappa_s_p2': current_lens_kappa_s_p2_entry, 
			'lens_xcen_ptype': current_lens_xcen_ptype_entry, 
			'lens_xcen_p1': current_lens_xcen_p1_entry, 
			'lens_xcen_p2': current_lens_xcen_p2_entry, 
			'lens_ycen_ptype': current_lens_ycen_ptype_entry, 
			'lens_ycen_p1': current_lens_ycen_p1_entry, 
			'lens_ycen_p2': current_lens_ycen_p2_entry, 
			'lens_pa_ptype': current_lens_pa_ptype_entry, 
			'lens_pa_p1': current_lens_pa_p1_entry, 
			'lens_pa_p2': current_lens_pa_p2_entry, 
			'lens_q_ptype': current_lens_q_ptype_entry, 
			'lens_q_p1': current_lens_q_p1_entry, 
			'lens_q_p2': current_lens_q_p2_entry, 
			'lens_r_s_ptype': current_lens_r_s_ptype_entry, 
			'lens_r_s_p1': current_lens_r_s_p1_entry, 
			'lens_r_s_p2': current_lens_r_s_p2_entry}
	elif model == 'ext. shear':
		Label(frame_mcmc, text="shear_amp", width=12).grid(row=current_row, column=current_column)
		current_lens_shear_amp_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_amp_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_shear_amp_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_amp_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_shear_amp_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_amp_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_shear_amp_ptype_entry.delete(0, END)
		current_lens_shear_amp_ptype_entry.insert(0, 0)	
		current_lens_shear_amp_p1_entry.delete(0, END)
		current_lens_shear_amp_p1_entry.insert(0, 0.0)	
		current_lens_shear_amp_p2_entry.delete(0, END)
		current_lens_shear_amp_p2_entry.insert(0, 0.20)	

		current_row+=1
		Label(frame_mcmc, text="shear_pa", width=12).grid(row=current_row, column=current_column)
		current_lens_shear_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_lens_shear_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_lens_shear_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_lens_shear_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_lens_shear_pa_ptype_entry.delete(0, END)
		current_lens_shear_pa_ptype_entry.insert(0, 0)	
		current_lens_shear_pa_p1_entry.delete(0, END)
		current_lens_shear_pa_p1_entry.insert(0, 0.1)	
		current_lens_shear_pa_p2_entry.delete(0, END)
		current_lens_shear_pa_p2_entry.insert(0, 180.0)	

		lens_priorpars_entry_i={
			'lens_shear_amp_ptype': current_lens_shear_amp_ptype_entry, 
			'lens_shear_amp_p1': current_lens_shear_amp_p1_entry, 
			'lens_shear_amp_p2': current_lens_shear_amp_p2_entry, 
			'lens_shear_pa_ptype': current_lens_shear_pa_ptype_entry, 
			'lens_shear_pa_p1': current_lens_shear_pa_p1_entry, 
			'lens_shear_pa_p2': current_lens_shear_pa_p2_entry}

	lens_priorpars_entry_list.append(lens_priorpars_entry_i)

def add_a_source_entry_mcmc(model, row, column):
	global source_priorpars_entry_list
	current_row=row
	current_column=column
	if model == 'sersic' or model == 'gaussian':
		Label(frame_mcmc, text="src_amp", width=12).grid(row=current_row, column=current_column)
		current_src_amp_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_amp_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_amp_p1_entry=Entry(frame_mcmc, width=12)
		current_src_amp_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_amp_p2_entry=Entry(frame_mcmc, width=12)
		current_src_amp_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_amp_ptype_entry.delete(0, END)
		current_src_amp_ptype_entry.insert(0, 0)	
		current_src_amp_p1_entry.delete(0, END)
		current_src_amp_p1_entry.insert(0, 0.01)	
		current_src_amp_p2_entry.delete(0, END)
		current_src_amp_p2_entry.insert(0, 2.0)	

		current_row+=1
		Label(frame_mcmc, text="src_xcen", width=12).grid(row=current_row, column=current_column)
		current_src_xcen_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_xcen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_xcen_p1_entry=Entry(frame_mcmc, width=12)
		current_src_xcen_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_xcen_p2_entry=Entry(frame_mcmc, width=12)
		current_src_xcen_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_xcen_ptype_entry.delete(0, END)
		current_src_xcen_ptype_entry.insert(0, 0)	
		current_src_xcen_p1_entry.delete(0, END)
		current_src_xcen_p1_entry.insert(0, -0.2)	
		current_src_xcen_p2_entry.delete(0, END)
		current_src_xcen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="src_ycen", width=12).grid(row=current_row, column=current_column)
		current_src_ycen_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_ycen_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_ycen_p1_entry=Entry(frame_mcmc, width=12)
		current_src_ycen_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_ycen_p2_entry=Entry(frame_mcmc, width=12)
		current_src_ycen_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_ycen_ptype_entry.delete(0, END)
		current_src_ycen_ptype_entry.insert(0, 0)	
		current_src_ycen_p1_entry.delete(0, END)
		current_src_ycen_p1_entry.insert(0, -0.2)	
		current_src_ycen_p2_entry.delete(0, END)
		current_src_ycen_p2_entry.insert(0, 0.2)	

		current_row+=1
		Label(frame_mcmc, text="src_sigma", width=12).grid(row=current_row, column=current_column)
		current_src_sigma_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_sigma_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_sigma_p1_entry=Entry(frame_mcmc, width=12)
		current_src_sigma_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_sigma_p2_entry=Entry(frame_mcmc, width=12)
		current_src_sigma_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_sigma_ptype_entry.delete(0, END)
		current_src_sigma_ptype_entry.insert(0, 0)	
		current_src_sigma_p1_entry.delete(0, END)
		current_src_sigma_p1_entry.insert(0, 0.0001)	
		current_src_sigma_p2_entry.delete(0, END)
		current_src_sigma_p2_entry.insert(0, 3.0)	

		current_row+=1
		Label(frame_mcmc, text="src_pa", width=12).grid(row=current_row, column=current_column)
		current_src_pa_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_pa_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_pa_p1_entry=Entry(frame_mcmc, width=12)
		current_src_pa_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_pa_p2_entry=Entry(frame_mcmc, width=12)
		current_src_pa_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_pa_ptype_entry.delete(0, END)
		current_src_pa_ptype_entry.insert(0, 0)	
		current_src_pa_p1_entry.delete(0, END)
		current_src_pa_p1_entry.insert(0, 1)	
		current_src_pa_p2_entry.delete(0, END)
		current_src_pa_p2_entry.insert(0, 180)	

		current_row+=1
		Label(frame_mcmc, text="src_q", width=12).grid(row=current_row, column=current_column)
		current_src_q_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_q_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_q_p1_entry=Entry(frame_mcmc, width=12)
		current_src_q_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_q_p2_entry=Entry(frame_mcmc, width=12)
		current_src_q_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_q_ptype_entry.delete(0, END)
		current_src_q_ptype_entry.insert(0, 0)	
		current_src_q_p1_entry.delete(0, END)
		current_src_q_p1_entry.insert(0, 0.1)	
		current_src_q_p2_entry.delete(0, END)
		current_src_q_p2_entry.insert(0, 1.0)	

		current_row+=1
		Label(frame_mcmc, text="src_n", width=12).grid(row=current_row, column=current_column)
		current_src_n_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_n_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_n_p1_entry=Entry(frame_mcmc, width=12)
		current_src_n_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_n_p2_entry=Entry(frame_mcmc, width=12)
		current_src_n_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_n_ptype_entry.delete(0, END)
		current_src_n_ptype_entry.insert(0, 0)	
		current_src_n_p1_entry.delete(0, END)
		current_src_n_p1_entry.insert(0, 0.1)	
		current_src_n_p2_entry.delete(0, END)
		current_src_n_p2_entry.insert(0, 10.0)	
		
		src_priorpars_entry_i={
			'src_amp_ptype': current_src_amp_ptype_entry, 
			'src_amp_p1': current_src_amp_p1_entry, 
			'src_amp_p2': current_src_amp_p2_entry, 
			'src_xcen_ptype': current_src_xcen_ptype_entry, 
			'src_xcen_p1': current_src_xcen_p1_entry, 
			'src_xcen_p2': current_src_xcen_p2_entry, 
			'src_ycen_ptype': current_src_ycen_ptype_entry, 
			'src_ycen_p1': current_src_ycen_p1_entry, 
			'src_ycen_p2': current_src_ycen_p2_entry, 
			'src_sigma_ptype': current_src_sigma_ptype_entry, 
			'src_sigma_p1': current_src_sigma_p1_entry, 
			'src_sigma_p2': current_src_sigma_p2_entry, 
			'src_pa_ptype': current_src_pa_ptype_entry, 
			'src_pa_p1': current_src_pa_p1_entry, 
			'src_pa_p2': current_src_pa_p2_entry, 
			'src_q_ptype': current_src_q_ptype_entry, 
			'src_q_p1': current_src_q_p1_entry, 
			'src_q_p2': current_src_q_p2_entry, 
			'src_n_ptype': current_src_n_ptype_entry, 
			'src_n_p1': current_src_n_p1_entry, 
			'src_n_p2': current_src_n_p2_entry}
	elif model == 'pix':
		Label(frame_mcmc, text="lambda", width=12).grid(row=current_row, column=current_column)
		current_src_lambda_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_lambda_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_lambda_p1_entry=Entry(frame_mcmc, width=12)
		current_src_lambda_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_lambda_p2_entry=Entry(frame_mcmc, width=12)
		current_src_lambda_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_lambda_ptype_entry.delete(0, END)
		current_src_lambda_ptype_entry.insert(0, 0)	
		current_src_lambda_p1_entry.delete(0, END)
		current_src_lambda_p1_entry.insert(0, 0.)	
		current_src_lambda_p2_entry.delete(0, END)
		current_src_lambda_p2_entry.insert(0, 2000.0)	

		current_row+=1
		Label(frame_mcmc, text="reg. scheme", width=12).grid(row=current_row, column=current_column)
		current_src_scheme_ptype_entry=Entry(frame_mcmc, width=12)
		current_src_scheme_ptype_entry.grid(row=current_row, column=current_column+1)
		current_src_scheme_p1_entry=Entry(frame_mcmc, width=12)
		current_src_scheme_p1_entry.grid(row=current_row, column=current_column+2)
		current_src_scheme_p2_entry=Entry(frame_mcmc, width=12)
		current_src_scheme_p2_entry.grid(row=current_row, column=current_column+3)
		current_src_scheme_ptype_entry.delete(0, END)
		current_src_scheme_ptype_entry.insert(0, 0)	
		current_src_scheme_p1_entry.delete(0, END)
		current_src_scheme_p1_entry.insert(0, 1)	
		current_src_scheme_p2_entry.delete(0, END)
		current_src_scheme_p2_entry.insert(0, 1)	
		
		src_priorpars_entry_i={
			'src_lambda_ptype': current_src_lambda_ptype_entry, 
			'src_lambda_p1': current_src_lambda_p1_entry, 
			'src_lambda_p2': current_src_lambda_p2_entry, 
			'src_scheme_ptype': current_src_scheme_ptype_entry, 
			'src_scheme_p1': current_src_scheme_p1_entry, 
			'src_scheme_p2': current_src_scheme_p2_entry}

	source_priorpars_entry_list.append(src_priorpars_entry_i)

def lfit_mcmc():
	global frame_mcmc, lens_priorpars_entry_list, source_priorpars_entry_list, phot_priorpars_entry_list, lens_priorpars, src_priorpars, phot_priorpars, nthreads_entry, nwalker_entry, burnin_entry, post_burnin_entry
	popup_mcmc=Toplevel()
	popup_mcmc.wm_title("MCMC Fitting")
	popup_mcmc.protocol("WM_DELETE_WINDOW", lambda: popup_mcmc.destroy())
	frame_mcmc=Frame(popup_mcmc)
	for i in range(40):
		frame_mcmc.grid_rowconfigure(i, weight=1)
		frame_mcmc.grid_columnconfigure(i, weight=1)
	frame_mcmc.pack(fill=BOTH, expand=YES)

	n_par_lens=0
	lens_priorpars_entry_list=[]
	current_row=1
	current_column=0
	for i_lens in range(n_lens):
		if i_lens == 0:
			Label(frame_mcmc, text="lens model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_lens+=n_lens_par_dict[lensmodel[i_lens]]
		add_a_lens_entry_mcmc(lensmodel[i_lens], current_row, current_column)
		current_row+=n_lens_par_dict[lensmodel[i_lens]]
	lens_priorpars=np.zeros((n_par_lens, 3))
	current_column+=5

	n_par_src=0
	source_priorpars_entry_list=[]
	current_row=1
	for i_src in range(n_source):
		if i_src == 0:
			Label(frame_mcmc, text="source model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_src+=n_source_par_dict[sourcemodel[i_src]]
		add_a_source_entry_mcmc(sourcemodel[i_src], current_row, current_column)
		current_row+=n_source_par_dict[sourcemodel[i_src]]
	src_priorpars=np.zeros((n_par_src, 3))
	current_column+=5

	n_par_phot=0
	phot_priorpars_entry_list=[]
	current_row=1
	for i_phot in range(n_phot):
		if i_phot == 0:
			Label(frame_mcmc, text="phot model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_phot+=n_phot_par_dict[photmodel[i_phot]]
		add_a_phot_entry_mcmc(photmodel[i_phot], current_row, current_column)
		current_row+=n_phot_par_dict[photmodel[i_phot]]
	phot_priorpars=np.zeros((n_par_phot, 3))
	current_column+=5

	Label(frame_mcmc, text="nthreads").grid(row=36, column=current_column+1)
	nthreads_entry=Entry(frame_mcmc, width=12)
	nthreads_entry.grid(row=36, column=current_column+2, sticky=W)
	nthreads_entry.insert(0, 8)

	Label(frame_mcmc, text="nwalker").grid(row=37, column=current_column+1)
	nwalker_entry=Entry(frame_mcmc, width=12)
	nwalker_entry.grid(row=37, column=current_column+2, sticky=W)
	nwalker_entry.insert(0, 100)

	Label(frame_mcmc, text="burn-in").grid(row=38, column=current_column+1)
	burnin_entry=Entry(frame_mcmc, width=12)
	burnin_entry.grid(row=38, column=current_column+2, sticky=W)
	burnin_entry.insert(0, 1000)

	Label(frame_mcmc, text="post burn-in").grid(row=39, column=current_column+1)
	post_burnin_entry=Entry(frame_mcmc, width=12)
	post_burnin_entry.grid(row=39, column=current_column+2, sticky=W)
	post_burnin_entry.insert(0, 5000)

	Button(frame_mcmc, text='Run', width=12, command=run_lfit_mcmc).grid(row=40, column=current_column+1, pady=4)
	Label(frame_mcmc, text="* 0-uniform prior; 1-gaussian prior").grid(row=0, column=current_column+2)

def lfit_multinest():
	global frame_mcmc, lens_priorpars_entry_list, source_priorpars_entry_list, phot_priorpars_entry_list, lens_priorpars, src_priorpars, phot_priorpars, nlive_entry, ins_entry, resume_entry
	popup_mcmc=Toplevel()
	popup_mcmc.wm_title("MultiNest Fitting")
	popup_mcmc.protocol("WM_DELETE_WINDOW", lambda: popup_mcmc.destroy())
	frame_mcmc=Frame(popup_mcmc)
	for i in range(40):
		frame_mcmc.grid_rowconfigure(i, weight=1)
		frame_mcmc.grid_columnconfigure(i, weight=1)
	frame_mcmc.pack(fill=BOTH, expand=YES)

	n_par_lens=0
	lens_priorpars_entry_list=[]
	current_row=1
	current_column=0
	for i_lens in range(n_lens):
		if i_lens == 0:
			Label(frame_mcmc, text="lens model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_lens+=n_lens_par_dict[lensmodel[i_lens]]
		add_a_lens_entry_mcmc(lensmodel[i_lens], current_row, current_column)
		current_row+=n_lens_par_dict[lensmodel[i_lens]]
	lens_priorpars=np.zeros((n_par_lens, 3))
	current_column+=5

	n_par_src=0
	source_priorpars_entry_list=[]
	current_row=1
	for i_src in range(n_source):
		if i_src == 0:
			Label(frame_mcmc, text="source model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_src+=n_source_par_dict[sourcemodel[i_src]]
		add_a_source_entry_mcmc(sourcemodel[i_src], current_row, current_column)
		current_row+=n_source_par_dict[sourcemodel[i_src]]
	src_priorpars=np.zeros((n_par_src, 3))
	current_column+=5

	n_par_phot=0
	phot_priorpars_entry_list=[]
	current_row=1
	for i_phot in range(n_phot):
		if i_phot == 0:
			Label(frame_mcmc, text="phot model", width=12).grid(row=0, column=current_column)
			Label(frame_mcmc, text="prior type*", width=12).grid(row=0, column=current_column+1)
			Label(frame_mcmc, text="prior par1", width=12).grid(row=0, column=current_column+2)
			Label(frame_mcmc, text="prior par2", width=12).grid(row=0, column=current_column+3)
		n_par_phot+=n_phot_par_dict[photmodel[i_phot]]
		add_a_phot_entry_mcmc(photmodel[i_phot], current_row, current_column)
		current_row+=n_phot_par_dict[photmodel[i_phot]]
	phot_priorpars=np.zeros((n_par_phot, 3))
	current_column+=5

	Label(frame_mcmc, text="nlive").grid(row=36, column=current_column+1)
	nlive_entry=Entry(frame_mcmc, width=12)
	nlive_entry.grid(row=36, column=current_column+2, sticky=W)
	nlive_entry.insert(0, 300)

	Label(frame_mcmc, text="ins").grid(row=37, column=current_column+1)
	ins_entry=Entry(frame_mcmc, width=12)
	ins_entry.grid(row=37, column=current_column+2, sticky=W)
	ins_entry.insert(0, 'true')

	Label(frame_mcmc, text="resume").grid(row=38, column=current_column+1)
	resume_entry=Entry(frame_mcmc, width=12)
	resume_entry.grid(row=38, column=current_column+2, sticky=W)
	resume_entry.insert(0, 'false')

	Button(frame_mcmc, text='Run', width=12, command=run_lfit_multinest).grid(row=40, column=current_column+1, pady=4)
	Label(frame_mcmc, text="* 0-uniform prior; 1-gaussian prior").grid(row=0, column=current_column+2)

def get_source_priors():
	global scr_priorpars
	current_index=0
	for i_src in range(n_source):
		if sourcemodel[i_src] == 'sersic' or sourcemodel[i_src] == 'gaussian':
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_amp_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_amp_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_amp_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_xcen_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_xcen_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_xcen_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_ycen_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_ycen_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_ycen_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_sigma_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_sigma_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_sigma_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_pa_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_pa_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_pa_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_q_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_q_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_q_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_n_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_n_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_n_p2'].get())
			current_index+=1
		elif sourcemodel[i_src] == 'pix':
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_lambda_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_lambda_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_lambda_p2'].get())
			current_index+=1
			src_priorpars[current_index, 0]=float(source_priorpars_entry_list[i_src]['src_scheme_ptype'].get())
			src_priorpars[current_index, 1]=float(source_priorpars_entry_list[i_src]['src_scheme_p1'].get())
			src_priorpars[current_index, 2]=float(source_priorpars_entry_list[i_src]['src_scheme_p2'].get())
			current_index+=1

def get_lens_priors():
	global lens_priorpars
	current_index=0
	for i_lens in range(n_lens):
		if lensmodel[i_lens] == 'sie' or lensmodel[i_lens] == 'sple':
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_pa_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_pa_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_pa_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_q_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_q_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_q_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_gamma_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_gamma_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_gamma_p2'].get())
			current_index+=1
		elif lensmodel[i_lens] == 'pm':
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_theta_e_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_theta_e_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_theta_e_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())
			current_index+=1
		elif lensmodel[i_lens] == 'snfw':
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_kappa_s_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_kappa_s_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_kappa_s_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_pa_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_pa_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_pa_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_q_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_q_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_q_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_r_s_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_r_s_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_r_s_p2'].get())
			current_index+=1
		elif lensmodel[i_lens] == 'ext. shear':
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_shear_amp_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_shear_amp_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_shear_amp_p2'].get())
			current_index+=1
			lens_priorpars[current_index, 0]=float(lens_priorpars_entry_list[i_lens]['lens_shear_pa_ptype'].get())
			lens_priorpars[current_index, 1]=float(lens_priorpars_entry_list[i_lens]['lens_shear_pa_p1'].get())
			lens_priorpars[current_index, 2]=float(lens_priorpars_entry_list[i_lens]['lens_shear_pa_p2'].get())
			current_index+=1

def get_phot_priors():
	global phot_priorpars
	current_index=0
	for i_phot in range(n_phot):
		if photmodel[i_phot] == 'sersic':
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_amp_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_pa_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_q_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_q_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_q_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_n_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_n_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_n_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_bg_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p2'].get())
			current_index+=1
		elif photmodel[i_phot] == 'csersic':
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_amp_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_pa_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_q_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_q_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_q_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_n_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_n_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_n_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_rc_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_rc_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_rc_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_alpha_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_alpha_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_alpha_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_gamma_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_gamma_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_gamma_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_bg_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p2'].get())
			current_index+=1
		elif photmodel[i_phot] == 'hernquist':
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_amp_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_amp_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_xcen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_ycen_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_rs_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_rs_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_rs_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_pa_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_pa_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_q_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_q_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_q_p2'].get())
			current_index+=1
			phot_priorpars[current_index, 0]=float(phot_priorpars_entry_list[i_phot]['phot_bg_ptype'].get())
			phot_priorpars[current_index, 1]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p1'].get())
			phot_priorpars[current_index, 2]=float(phot_priorpars_entry_list[i_phot]['phot_bg_p2'].get())
			current_index+=1

def config_sky_priors():
	global npar_sky
	prior_type_dict={'0': 'unif', '1': 'norm'}
	prior_type=prior_type_dict[phot_priorpars_entry_list[0]['phot_bg_ptype'].get()]
	arg1=float(phot_priorpars_entry_list[0]['phot_bg_p1'].get())
	arg2=float(phot_priorpars_entry_list[0]['phot_bg_p2'].get())
	if (prior_type == 0 and arg1 == arg2):
		lensed_ini.write('foregr.bg = '+str(arg1)+'\n')
	else:
		lensed_ini.write('foregr.bg = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
	lensed_ini.write('foregr.dx = 0\n')
	lensed_ini.write('foregr.dy = 0\n')
	npar_sky=3

def config_phot_priors():
	global npar_phot
	prior_type_dict={'0': 'unif', '1': 'norm'}
	for i_phot in range(n_phot):
		if photmodel[i_phot] == 'sersic':
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_amp_ptype'].get()]
			arg1=-2.5*np.log10(float(phot_priorpars_entry_list[i_phot]['phot_amp_p2'].get()))
			arg2=-2.5*np.log10(float(phot_priorpars_entry_list[i_phot]['phot_amp_p1'].get()))
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.mag = '+str(arg1)+'\n')
			else:
				lensed_ini.write('host'+str(i_phot+1)+'.mag = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_xcen_ptype'].get()]
			arg1=hw+float(phot_priorpars_entry_list[i_phot]['phot_xcen_p1'].get())/dpix+1
			arg2=hw+float(phot_priorpars_entry_list[i_phot]['phot_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.x = '+str(arg1)+'\n')
			else: 
				lensed_ini.write('host'+str(i_phot+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_ycen_ptype'].get()]
			arg1=hw+float(phot_priorpars_entry_list[i_phot]['phot_ycen_p1'].get())/dpix+1
			arg2=hw+float(phot_priorpars_entry_list[i_phot]['phot_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('host'+str(i_phot+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_sigma_ptype'].get()]
			arg1=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p1'].get())/dpix
			arg2=float(phot_priorpars_entry_list[i_phot]['phot_sigma_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.r = '+str(arg1)+'\n')
			else:					
				lensed_ini.write('host'+str(i_phot+1)+'.r = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_pa_ptype'].get()]
			arg1=float(phot_priorpars_entry_list[i_phot]['phot_pa_p1'].get())
			arg2=float(phot_priorpars_entry_list[i_phot]['phot_pa_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.pa = '+str(arg1)+'\n')
			else:
				lensed_ini.write('host'+str(i_phot+1)+'.pa = wrap '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_q_ptype'].get()]
			arg1=float(phot_priorpars_entry_list[i_phot]['phot_q_p1'].get())
			arg2=float(phot_priorpars_entry_list[i_phot]['phot_q_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.q = '+str(arg1)+'\n')
			else:
				lensed_ini.write('host'+str(i_phot+1)+'.q = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1
			prior_type=prior_type_dict[phot_priorpars_entry_list[i_phot]['phot_n_ptype'].get()]
			arg1=float(phot_priorpars_entry_list[i_phot]['phot_n_p1'].get())
			arg2=float(phot_priorpars_entry_list[i_phot]['phot_n_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('host'+str(i_phot+1)+'.n = '+str(arg1)+'\n')
			else:
				lensed_ini.write('host'+str(i_phot+1)+'.n = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_phot+=1

def config_lens_priors():
	global npar_lens
	prior_type_dict={'0': 'unif', '1': 'norm'}
	for i_lens in range(n_lens):
		if lensmodel[i_lens] == 'sie':
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_b_SIE_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p1'].get())/dpix
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_pa_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_pa_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_pa_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.pa = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.pa = wrap '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_q_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_q_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_q_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.q = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.q = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
		elif lensmodel[i_lens] == 'sple':
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_b_SIE_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p1'].get())/dpix
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_b_SIE_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_pa_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_pa_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_pa_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.pa = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.pa = wrap '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_q_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_q_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_q_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.q = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.q = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_gamma_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_gamma_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_gamma_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.t = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.t = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
		elif lensmodel[i_lens] == 'pm':
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_theta_e_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_theta_e_p1'].get())/dpix
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_theta_e_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.r = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
		elif lensmodel[i_lens] == 'ext. shear':
			lensed_ini.write('lens'+str(i_lens+1)+'.x = '+str(hw+1.0)+'\n')
			lensed_ini.write('lens'+str(i_lens+1)+'.y = '+str(hw+1.0)+'\n')
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_shear_amp_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_shear_amp_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_shear_amp_p2'].get())
			lensed_ini.write('lens'+str(i_lens+1)+'.g1 = '+prior_type+' '+str(-arg2)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_shear_pa_ptype'].get()]
			lensed_ini.write('lens'+str(i_lens+1)+'.g2 = '+prior_type+' '+str(-arg2)+' '+str(arg2)+'\n')
			npar_lens+=1
		elif lensmodel[i_lens] == 'snfw':
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_kappa_s_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_kappa_s_p1'].get())
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_kappa_s_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.kappa = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.kappa = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_xcen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_ycen_ptype'].get()]
			arg1=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p1'].get())/dpix+1
			arg2=hw+float(lens_priorpars_entry_list[i_lens]['lens_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1
			prior_type=prior_type_dict[lens_priorpars_entry_list[i_lens]['lens_r_s_ptype'].get()]
			arg1=float(lens_priorpars_entry_list[i_lens]['lens_r_s_p1'].get())/dpix
			arg2=float(lens_priorpars_entry_list[i_lens]['lens_r_s_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('lens'+str(i_lens+1)+'.rs = '+str(arg1)+'\n')
			else:
				lensed_ini.write('lens'+str(i_lens+1)+'.rs = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_lens+=1

def config_source_priors():
	global npar_source
	prior_type_dict={'0': 'unif', '1': 'norm'}
	for i_source in range(n_source):
		if sourcemodel[i_source] == 'sersic':
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_amp_ptype'].get()]
			arg1=-2.5*np.log10(float(source_priorpars_entry_list[i_source]['src_amp_p2'].get()))
			arg2=-2.5*np.log10(float(source_priorpars_entry_list[i_source]['src_amp_p1'].get()))
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.mag = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.mag = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_xcen_ptype'].get()]
			arg1=hw+float(source_priorpars_entry_list[i_source]['src_xcen_p1'].get())/dpix+1
			arg2=hw+float(source_priorpars_entry_list[i_source]['src_xcen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.x = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.x = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_ycen_ptype'].get()]
			arg1=hw+float(source_priorpars_entry_list[i_source]['src_ycen_p1'].get())/dpix+1
			arg2=hw+float(source_priorpars_entry_list[i_source]['src_ycen_p2'].get())/dpix+1
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.y = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.y = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_sigma_ptype'].get()]
			arg1=float(source_priorpars_entry_list[i_source]['src_sigma_p1'].get())/dpix
			arg2=float(source_priorpars_entry_list[i_source]['src_sigma_p2'].get())/dpix
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.r = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.r = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_pa_ptype'].get()]
			arg1=float(source_priorpars_entry_list[i_source]['src_pa_p1'].get())
			arg2=float(source_priorpars_entry_list[i_source]['src_pa_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.pa = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.pa = wrap '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_q_ptype'].get()]
			arg1=float(source_priorpars_entry_list[i_source]['src_q_p1'].get())
			arg2=float(source_priorpars_entry_list[i_source]['src_q_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.q = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.q = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1
			prior_type=prior_type_dict[source_priorpars_entry_list[i_source]['src_n_ptype'].get()]
			arg1=float(source_priorpars_entry_list[i_source]['src_n_p1'].get())
			arg2=float(source_priorpars_entry_list[i_source]['src_n_p2'].get())
			if (prior_type == 'unif' and arg1 == arg2):
				lensed_ini.write('source'+str(i_source+1)+'.n = '+str(arg1)+'\n')
			else:
				lensed_ini.write('source'+str(i_source+1)+'.n = '+prior_type+' '+str(arg1)+' '+str(arg2)+'\n')
			npar_source+=1

def run_lfit_multinest():
	global I_data1, x1, y1, I_invvar1, jmask1, lensed_ini, I_phot, I_fit, I_source, I_source_err, dpix_source, source_xbase, source_ybase, lpar, spar, ppar, lpar_err, spar_err, ppar_err, npar_sky, npar_phot, npar_lens, npar_source

	To_lensed_dict={'sie': 'sie', 'sple': 'epl', 'ext. shear': 'shear', 'pm': 'point_mass', 'snfw': 'snfw', 'sersic': 'sersic_lfit'}

	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_image1=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_bspline=I_lens[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	jmask1=jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	lf.I_bspline=I_bspline

	os.system("mkdir multinest_output_dir")
	lensed_ini=open(obj_name+'_lensed.ini', 'w')
	lensed_ini.write('; Configuration file for '+obj_name+'\n')
	lensed_ini.write('gain = 1.0\n')
	lensed_ini.write('offset = 0.0\n')
	lensed_ini.write('rule = point\n')
	lensed_ini.write('output = true\n')
	lensed_ini.write('root = multinest_output_dir/'+obj_name+'\n')

	fname='multinest_output_dir/'+obj_name+'_image_temp.fits'
	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		hdu=pyfits.PrimaryHDU(I_data1-I_bspline)
	else:
		hdu=pyfits.PrimaryHDU(I_data1)
	hdu.header['target']=obj_name
	hdu.header['dpix_l']=dpix
	hdulist=pyfits.HDUList([hdu])
	hdulist.writeto(fname, clobber=True)
	lensed_ini.write('image = '+fname+'\n')

	fname='multinest_output_dir/'+obj_name+'_weight_temp.fits'	
	hdu=pyfits.PrimaryHDU(I_invvar1)
	hdu.header['target']=obj_name
	hdu.header['dpix_l']=dpix
	hdulist=pyfits.HDUList([hdu])
	hdulist.writeto(fname, clobber=True)
	lensed_ini.write('weight = '+fname+'\n')

	fname='multinest_output_dir/'+obj_name+'_psf_temp.fits'	
	hdu=pyfits.PrimaryHDU(tpsf)
	hdulist=pyfits.HDUList([hdu])
	hdulist.writeto(fname, clobber=True)
	lensed_ini.write('psf = '+fname+'\n')

#	The mask seems not to be working correctly. Need an idea.
	fname='multinest_output_dir/'+obj_name+'_mask_temp.fits'	
	hdu=pyfits.PrimaryHDU(1-jmask1)
	hdulist=pyfits.HDUList([hdu])
	hdulist.writeto(fname, clobber=True)
	lensed_ini.write('mask = '+fname+'\n')

	nlive=int(nlive_entry.get())
	lensed_ini.write('nlive = '+str(nlive)+'\n')
	ins=ins_entry.get()
	lensed_ini.write('ins = '+ins+'\n')
	resume=resume_entry.get()
	lensed_ini.write('resume = '+resume+'\n')

	lensed_ini.write('\n')
	lensed_ini.write('[objects]\n')
	
	for i_phot in np.arange(n_phot):
		if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
			continue
		else:
			if i_phot == 0:
				lensed_ini.write('foregr = sky\n')			
			lensed_ini.write('host'+str(i_phot+1)+' = '+To_lensed_dict[photmodel[i_phot]]+'\n')
	for i_lens in np.arange(n_lens):
		lensed_ini.write('lens'+str(i_lens+1)+' = '+To_lensed_dict[lensmodel[i_lens]]+'\n')
	for i_source in np.arange(n_source):
		lensed_ini.write('source'+str(i_source+1)+' = '+To_lensed_dict[sourcemodel[i_source]]+'\n')

	lensed_ini.write('\n')
	lensed_ini.write('[priors]\n')
	npar_sky=0
	npar_phot=0
	npar_lens=0
	npar_source=0
	if (n_phot > 0) and photmodel[0] != 'bspline': 
		config_sky_priors()
		config_phot_priors()
	config_lens_priors()
	config_source_priors()
	lensed_ini.close()

	os.system("lensed "+obj_name+'_lensed.ini')

	npar=npar_sky+npar_phot+npar_lens+npar_source
	mapping_dict={}
	paraname=open('./multinest_output_dir/'+obj_name+'.paramnames', 'r')
	res=open('./multinest_output_dir/'+obj_name+'stats.dat', 'r')
# 	skip the first two lines in the stats.dat file
	for i in np.arange(4): 
		junk=res.readline()
	for i in np.arange(npar):
		key=(paraname.readline()).rstrip()
		junk, value, err=(res.readline()).split()
		value=float(value)
		err=float(err)
		mapping_dict[key]=value
		mapping_dict[key+'_err']=err
	paraname.close()
	res.close()

	ppar=[]
	ppar_err=[]
	for i_phot in np.arange(n_phot):
		if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
			continue
		else: 
			if (photmodel[i_phot] == 'sersic'):
				# conversion from L_tot defined in Lensed to my I_0 
				L_tot=10.0**(-0.4*mapping_dict['host'+str(i_phot+1)+'.mag'])
				r=mapping_dict['host'+str(i_phot+1)+'.r']
				q=mapping_dict['host'+str(i_phot+1)+'.q']
				n=mapping_dict['host'+str(i_phot+1)+'.n']
				ppar.append(L_tot)
				ppar.append((mapping_dict['host'+str(i_phot+1)+'.x']-hw-1)*dpix)
				ppar.append((mapping_dict['host'+str(i_phot+1)+'.y']-hw-1)*dpix)
				ppar.append(r*dpix)
				ppar.append((mapping_dict['host'+str(i_phot+1)+'.pa']+90.0) % 180.0)
				ppar.append(q)
				ppar.append(n)
				if i_phot==0:
					ppar.append(mapping_dict['foregr.bg'])
				else: 
					ppar.append(0.0)
				ppar_err.append(0.4*np.log(10.0)*mapping_dict['host'+str(i_phot+1)+'.mag_err']*L_tot)
				ppar_err.append((mapping_dict['host'+str(i_phot+1)+'.x_err'])*dpix)
				ppar_err.append((mapping_dict['host'+str(i_phot+1)+'.y_err'])*dpix)
				ppar_err.append(mapping_dict['host'+str(i_phot+1)+'.r_err']*dpix)
				ppar_err.append(mapping_dict['host'+str(i_phot+1)+'.pa_err'])
				ppar_err.append(mapping_dict['host'+str(i_phot+1)+'.q_err'])
				ppar_err.append(mapping_dict['host'+str(i_phot+1)+'.n_err'])
				if i_phot==0:
					ppar_err.append(mapping_dict['foregr.bg_err'])
				else: 
					ppar_err.append(0.0)

	lpar=[]
	lpar_err=[]
	for i_lens in np.arange(n_lens):
		if (lensmodel[i_lens] == 'sie'):
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.r']*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.x']-hw-1)*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.y']-hw-1)*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.pa']+90.0) % 180.0)
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.q'])
			lpar.append(1.0)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.r_err']*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.x_err'])*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.y_err'])*dpix)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.pa_err'])
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.q_err'])
			lpar_err.append(0.0)
		elif (lensmodel[i_lens] == 'sple'):
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.r']*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.x']-hw-1)*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.y']-hw-1)*dpix)
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.pa'])
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.q'])
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.t'])
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.r_err']*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.x_err'])*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.y_err'])*dpix)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.pa_err'])
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.q_err'])
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.t_err'])
		elif (lensmodel[i_lens] == 'pm'):
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.r']*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.x']-hw-1)*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.y']-hw-1)*dpix)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.r_err']*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.x_err'])*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.y_err'])*dpix)
		elif (lensmodel[i_lens] == 'ext. shear'):
			g1=mapping_dict['lens'+str(i_lens+1)+'.g1']
			g2=mapping_dict['lens'+str(i_lens+1)+'.g2']
			g1_err=mapping_dict['lens'+str(i_lens+1)+'.g1_err']
			g2_err=mapping_dict['lens'+str(i_lens+1)+'.g2_err']
			amp=np.sqrt(g1**2.0+g2**2.0)
			amp_err=np.sqrt(g1**2.0*g1_err**2.0+g2**2.0*g2_err**2.0)/amp
			pa=(np.arctan2(g2, g1)+np.pi)/2.0/np.pi*180.0 # the extra pi is used to shift the range from [-pi, pi] to [0, 2pi]
			pa_err=0.5*np.sqrt(g2**2.0*g1_err**2.0+g1**2.0*g2_err**2.0)/amp**2.0/np.pi*180.0
			lpar.append(amp)
			lpar.append(pa)
			lpar_err.append(amp_err)
			lpar_err.append(pa_err)
		elif (lensmodel[i_lens] == 'snfw'):
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.kappa'])
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.x']-hw-1)*dpix)
			lpar.append((mapping_dict['lens'+str(i_lens+1)+'.y']-hw-1)*dpix)
			lpar.append(0.0)
			lpar.append(1.0)
			lpar.append(mapping_dict['lens'+str(i_lens+1)+'.rs']*dpix)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.kappa_err'])
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.x_err'])*dpix)
			lpar_err.append((mapping_dict['lens'+str(i_lens+1)+'.y_err'])*dpix)
			lpar_err.append(0.0)
			lpar_err.append(0.0)
			lpar_err.append(mapping_dict['lens'+str(i_lens+1)+'.rs_err']*dpix)	

	spar=[]
	spar_err=[]
	for i_source in np.arange(n_source):
		if (sourcemodel[i_source] == 'sersic'):
			L_tot=10.0**(-0.4*mapping_dict['source'+str(i_source+1)+'.mag'])
			r=mapping_dict['source'+str(i_source+1)+'.r']
			q=mapping_dict['source'+str(i_source+1)+'.q']
			n=mapping_dict['source'+str(i_source+1)+'.n']
			spar.append(L_tot)
			spar.append((mapping_dict['source'+str(i_source+1)+'.x']-hw-1)*dpix)
			spar.append((mapping_dict['source'+str(i_source+1)+'.y']-hw-1)*dpix)
			spar.append(r*dpix)
			spar.append((mapping_dict['source'+str(i_source+1)+'.pa']) % 180.0)
			spar.append(q)
			spar.append(n)
			spar_err.append(0.4*np.log(10.0)*mapping_dict['source'+str(i_source+1)+'.mag_err']*L_tot)
			spar_err.append((mapping_dict['source'+str(i_source+1)+'.x_err'])*dpix)
			spar_err.append((mapping_dict['source'+str(i_source+1)+'.y_err'])*dpix)
			spar_err.append(mapping_dict['source'+str(i_source+1)+'.r_err']*dpix)
			spar_err.append(mapping_dict['source'+str(i_source+1)+'.pa_err'])
			spar_err.append(mapping_dict['source'+str(i_source+1)+'.q_err'])
			spar_err.append(mapping_dict['source'+str(i_source+1)+'.n_err'])

	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=lf.phot_model(x1, y1, ppar, photmodel, PSF=tpsf, n_phot=n_phot)

	I_fit=lf.model(x1, y1, lpar, spar, lensmodel, sourcemodel, PSF=tpsf, n_source=n_source, n_lens=n_lens)

	n_x_source=10*100+1
	n_y_source=10*100+1
	dpix_source=0.1*dpix
	source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2)*dpix_source
	source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2, np.ones(n_x_source))*dpix_source
	I_source=0.0*source_xbase
	current_source_index=0
	for i_source in np.arange(n_source):
		spar_i=spar[current_source_index:current_source_index+n_source_par_dict[sourcemodel[i_source]]]
		I_source+=lf.sersic_2d(source_xbase, source_ybase, spar_i)
		current_source_index+=n_source_par_dict[sourcemodel[i_source]]

	update_phot_parameters(ppar, n_phot, photmodel)
	update_lens_parameters(lpar, n_lens, lensmodel)
	update_source_parameters(spar, n_source, sourcemodel)

def run_lfit_mcmc():
	global I_data1, x1, y1, I_invvar1, I_phot, I_fit, I_source, I_source_err, dpix_source, source_xbase, source_ybase, lpar, spar, ppar, lpar_err, spar_err, ppar_err

	get_lens_priors()
	get_source_priors()
	get_phot_priors()
	priorpars=np.vstack((src_priorpars, lens_priorpars, phot_priorpars))

	x1=x[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	y1=y[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_data1=I_data[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_image1=I_image[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_bspline=I_lens[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	I_invvar1=I_invvar[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	jmask1=jmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	gmask1=gmask[(n_x-1)/2-hw:(n_x-1)/2+hw+1, (n_y-1)/2-hw:(n_y-1)/2+hw+1]
	lf.I_bspline=I_bspline
	
	ndim=(priorpars.shape)[0]
	nwalkers=int(nwalker_entry.get())
	print 'Setting up the sampler...............'
	if n_source == 1 and sourcemodel[0] == 'pix':
		lens_mask_2d=1-gmask1
		lens_mask_1d=np.where(lens_mask_2d == 1)
		(used_pixels, used_pixels_mask)=lf.select_used_pixels(I_image1, lens_mask_1d, lens_mask_1d[0].shape[0]/2)
		used_pixels_2d=np.zeros_like(x1, dtype=int)
		used_pixels_2d[used_pixels]=1
		PSF_csr=lf.generate_psf_csr(lens_mask_1d, x1.shape[0], x1.shape[1], tpsf)
		sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob_pix, args=(x1, y1, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data1, I_invvar1*jmask1, lens_mask_2d, used_pixels, used_pixels_mask, PSF, PSF_csr, priorpars), threads=int(nthreads_entry.get())) 
	else:
		print int(nthreads_entry.get())
		sampler=emcee.EnsembleSampler(nwalkers, ndim, lf.lnprob, args=(x1, y1, n_lens, lensmodel, n_source, sourcemodel, n_phot, photmodel, I_data1, I_invvar1*jmask1, PSF, priorpars), threads=int(nthreads_entry.get())) 
	print 'Setting up the sampler...Completed...'
	
	pos=init_pos(nwalkers, priorpars)
#	pos=[np.concatenate((spar, lpar, ppar))+0.01*np.concatenate((spar, lpar, ppar))*np.random.randn(ndim) for i in range(nwalkers)]
#	for i in range(nwalkers):
#		pos[i][12]=1.0
##		pos[i][0]=1000.0
##		pos[i][6]=4.0
	print 'Running MCMC.........................'
	# run xxx steps as a burn-in
	pos, prob, state=sampler.run_mcmc(pos, int(burnin_entry.get()))
	# reset the chain to remove the burn-in samples
	sampler.reset()
	print 'Running MCMC....burn-in completed....'
	# starting from the final position in the burn-in chain, sample for xxx
	sampler.run_mcmc(pos, int(post_burnin_entry.get()), rstate0=state)
	print 'Running MCMC........Completed........'
	print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
	print 'Autocorrelation time is ', sampler.acor

#	samples=sampler.chain[:, :, :].reshape((-1, ndim))
#	lnprobability=sampler.lnprobability.flatten()
#	lnprobability=lnprobability.reshape((lnprobability.size, 1))
#	mcmc_out=np.hstack((samples, lnprobability))

	samples=sampler.chain
	lnprobability=sampler.lnprobability
	lnprobability=lnprobability.reshape(((lnprobability.shape)[0], (lnprobability.shape)[1], 1))
	mcmc_out=np.concatenate((samples, lnprobability), axis=2)

	samples=samples[:, :, :].reshape((-1, ndim))
	res=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))	
	p_mcmc=[]
	p_err_mcmc=[]
	for i in range(ndim):
		p_mcmc.append(res[i][0])
		p_err_mcmc.append(0.5*(res[i][1]+res[i][2]))
	
	spar_index=0
	for i_source in np.arange(n_source):
		spar_index+=n_source_par_dict[sourcemodel[i_source]]
	lpar_index=0
	for i_lens in np.arange(n_lens):
		lpar_index+=n_lens_par_dict[lensmodel[i_lens]]

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

	if (n_phot==0 or (n_phot==1 and photmodel[0] == 'bspline')):
		I_phot=I_bspline
	else:
		I_phot=lf.phot_model(x1, y1, ppar, photmodel, PSF=PSF, n_phot=n_phot)

	if (len(sourcemodel) == 1 and sourcemodel[0] == 'pix'):
		lam=spar[0]
		reg_scheme=spar[1]
		output=lf.pix_source(I_data1-I_phot, I_invvar1*jmask1, lens_mask_2d, used_pixels, used_pixels_mask, x1, y1, lpar, lensmodel, n_lens=n_lens, PSF_csr=PSF_csr, lam=lam, reg_scheme=reg_scheme, return_cov=1)
		I_fit=output['fit']
		Vec_S=output['solution']
		Cov_matrix=output['Covariance']
		Vec_S_err=np.sqrt(np.diagonal(Cov_matrix))
	else:
		I_fit=lf.model(x1, y1, lpar, spar, lensmodel, sourcemodel, PSF=PSF, n_source=n_source, n_lens=n_lens)

	n_x_source=10*100+1
	n_y_source=10*100+1
	dpix_source=0.1*dpix
	source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2)*dpix_source
	source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2, np.ones(n_x_source))*dpix_source
	if len(sourcemodel) ==1 and sourcemodel[0] == 'pix':
		theta_x=x1[lens_mask_1d].flatten()
		theta_y=y1[lens_mask_1d].flatten()
		(alpha_x, alpha_y)=lf.def_total(x1, y1, lpar, lensmodel, n_lens=n_lens)
		alpha_x_1d=alpha_x[lens_mask_1d].flatten()
		alpha_y_1d=alpha_y[lens_mask_1d].flatten()
		
		tri=lf.construct_source_grid(theta_x, theta_y, alpha_x_1d, alpha_y_1d, used_pixels_mask)
		I_source=griddata(tri.points, Vec_S, (source_xbase, source_ybase), method='linear', fill_value=0.0)
		I_source_err=griddata(tri.points, Vec_S_err, (source_xbase, source_ybase), method='linear', fill_value=0.0)
	else:
		I_source=0.0*source_xbase
		current_source_index=0
		for i_source in np.arange(n_source):
			spar_i=spar[current_source_index:current_source_index+n_source_par_dict[sourcemodel[i_source]]]
			I_source+=lf.sersic_2d(source_xbase, source_ybase, spar_i)
			current_source_index+=n_source_par_dict[sourcemodel[i_source]]

	update_phot_parameters(ppar, n_phot, photmodel)
	update_lens_parameters(lpar, n_lens, lensmodel)
	update_source_parameters(spar, n_source, sourcemodel)

	fname=obj_name+'_SIE_mcmc_F814W.fits'
	col1=pyfits.Column(name='spar', format='D', array=spar)
	col1_err=pyfits.Column(name='spar_err', format='D', array=spar_err)
	col2=pyfits.Column(name='lpar', format='D', array=lpar)
	col2_err=pyfits.Column(name='lpar_err', format='D', array=lpar_err)
	col3=pyfits.Column(name='p_phot', format='D', array=ppar)
	col3_err=pyfits.Column(name='p_phot_err', format='D', array=ppar_err)
	tbhdu1=pyfits.BinTableHDU.from_columns([col1, col1_err])
	tbhdu2=pyfits.BinTableHDU.from_columns([col2, col2_err])
	tbhdu3=pyfits.BinTableHDU.from_columns([col3, col3_err])
	
	hdu=pyfits.PrimaryHDU(mcmc_out)
	hdu.header['target']=obj_name
	hdu.header['N_Lens']=n_lens
	for i_lens in np.arange(n_lens):
		hdu.header['lmodel_'+str(i_lens+1)]=lensmodel[i_lens]
	hdu.header['N_Source']=n_source
	for i_source in np.arange(n_source):
		hdu.header['smodel_'+str(i_source+1)]=sourcemodel[i_source]
	hdu.header['N_Phot']=n_phot
	for i_phot in np.arange(n_phot):
		hdu.header['pmodel_'+str(i_phot+1)]=photmodel[i_phot]
	try:
		hdu.header['DOF']=dof
	except:
		hdu.header['DOF']=I_data1.size
	try:
		hdu.header['Chi2']=chi2
	except:
		hdu.header['Chi2']=((I_data1-I_phot-I_fit)**2.0*I_invvar1*jmask1).sum()
	hdu.header['dpix_l']=dpix
	try:
		hdu.header['dpix_s']=dpix_source
	except:
		hdu.header['dpix_s']=0.1*dpix
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
	pyfits.append(fname, I_data1)
	pyfits.append(fname, I_invvar1*jmask1)
	pyfits.append(fname, I_phot)
	pyfits.append(fname, I_fit)
	pyfits.append(fname, I_source)
	if (len(sourcemodel) == 1 and sourcemodel[0] == 'pix'):
		hdu.header.add_comment('9th extension is the source error')
	 	pyfits.append(fname, I_source_err)
	print 'Saved the samples to '+fname
   
#	import corner
#	fig=corner.corner(samples[:, 2:7], quantiles=[0.16, 0.5, 0.84])
#	fig.savefig(obj_name+"_mcmc.png")
#	print "Saved the corner plot to "+obj_name+"_mcmc.png"

def show_phot_fit():
	popup=Toplevel()
	popup.wm_title("Foreground light fitting results")
	popup.protocol("WM_DELETE_WINDOW", lambda: on_closing(popup))
	frame=Frame(popup)
	for i in range(30):
		frame.grid_rowconfigure(i, weight=1)
		frame.grid_columnconfigure(i, weight=1)

	frame.pack(fill=BOTH, expand=YES)

	fig_phot, (ax1_phot, ax2_phot, ax3_phot)=plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
	canvas_phot=FigureCanvasTkAgg(fig_phot, master=frame)
	NavigationToolbar2TkAgg(canvas_phot, popup)	
	canvas_phot.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=15, sticky='NW', padx=4, pady=4)
	vmax=(I_data1-I_phot).max()
	vmin=(I_data1-I_phot).min()
	ext_lens=[x1.min(), x1.max(), y1.min(), y1.max()]
	ax1_phot.imshow(I_data1, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax1_phot.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax1_phot.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	ax2_phot.imshow(I_phot, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax2_phot.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax2_phot.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)
	ax3_phot.imshow(I_data1-I_phot, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax3_phot.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax3_phot.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	plt.subplots_adjust(left=0.05, right=0.98, bottom=0.15, wspace=0.16)
	canvas_phot.draw()
	canvas_phot.get_tk_widget().focus_set()
	
def show_fit():
	global save_fit_button, save_fit_button_text, save_fit_options, fig1, ax1, ax2, ax3, ax4, ax5, ax6, canvas_fit, zf_entry, vmin, vmax
	popup=Toplevel()
	popup.wm_title("Fitting results")
#	popup.protocol("WM_DELETE_WINDOW", lambda: popup.destroy())
	popup.protocol("WM_DELETE_WINDOW", lambda: on_closing(popup))
	frame=Frame(popup)
	for i in range(30):
		frame.grid_rowconfigure(i, weight=1)
		frame.grid_columnconfigure(i, weight=1)
	frame.bind("<B1-Motion>", lambda event: adjust_result(event))
	frame.pack(fill=BOTH, expand=YES)

	Label(frame, text="zoom factor").grid(row=0, column=21)
	zf_entry=Entry(frame, width=12)
	zf_entry.grid(row=0, column=22, sticky=W)
	zf_entry.insert(0, 1.0)
	zf_entry.bind('<Return>', adjust_source)

	fig1, (ax1, ax2, ax3, ax4, ax5, ax6)=plt.subplots(nrows=1, ncols=6, figsize=(24, 4))
	canvas_fit=FigureCanvasTkAgg(fig1, master=frame)
	NavigationToolbar2TkAgg(canvas_fit, popup)	
	canvas_fit.get_tk_widget().grid(row=2, column=0, rowspan=4, columnspan=24, sticky='NW', padx=4, pady=4)
	vmax=(I_data1-I_phot-I_fit).max()
	vmin=(I_data1-I_phot-I_fit).min()
	show_result(vmin=vmin, vmax=vmax, source_zoom=float(zf_entry.get()))

	save_fit_button_text=StringVar()
	save_fit_options=['fits', 'eps']
	save_fit_button_text.set('Export as')
	save_fit_button=ttk.Combobox(frame, textvariable=save_fit_button_text)
	save_fit_button.grid(row=0, column=0, padx=1, pady=1)
	save_fit_button['values']=tuple(save_fit_options)
	save_fit_button.config(width=15)
	save_fit_button.bind("<<ComboboxSelected>>", save_fit)

def adjust_result(event):
# adjust the color scheme of the fit by moving the mouse with left button being held down	
	global vmin, vmax
	imax=I_image1.max()
	imin=I_image1.min()
	x=min(event.x/800., 2.)
	y=min(event.y/800., 2.)
	vmax=(imax-imin)*x+imin
	vmin=(imax-imin)*y+imin
	if vmin > vmax:
		vmin=vmax
	show_result(vmin=vmin, vmax=vmax, source_zoom=float(zf_entry.get()))

def adjust_source(event):
# adjust the source size	
	show_result(vmin=vmin, vmax=vmax, source_zoom=float(zf_entry.get()))

def show_result(vmin=0., vmax=1., source_zoom=1.0):
	xsize, ysize=I_source.shape
	hw_x=(xsize-1)/2
	hw_y=(ysize-1)/2
	I_source_sub=I_source[hw_x-hw_x/source_zoom:hw_x+hw_x/source_zoom+1, hw_y-hw_y/source_zoom:hw_y+hw_y/source_zoom+1]

	ax1.clear()
	ax2.clear()
	ax3.clear()
	ax4.clear()
	ax5.clear()
	ax6.clear()
	
	ext_lens=[x1.min(), x1.max(), y1.min(), y1.max()]
	ext_source=[source_xbase[hw_x-hw_x/source_zoom:hw_x+hw_x/source_zoom+1, hw_y-hw_y/source_zoom:hw_y+hw_y/source_zoom+1].min(), source_xbase[hw_x-hw_x/source_zoom:hw_x+hw_x/source_zoom+1, hw_y-hw_y/source_zoom:hw_y+hw_y/source_zoom+1].max(), source_ybase[hw_x-hw_x/source_zoom:hw_x+hw_x/source_zoom+1, hw_y-hw_y/source_zoom:hw_y+hw_y/source_zoom+1].min(), source_ybase[hw_x-hw_x/source_zoom:hw_x+hw_x/source_zoom+1, hw_y-hw_y/source_zoom:hw_y+hw_y/source_zoom+1].max()]
	ax1.imshow(I_data1, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax1.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax1.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	ax2.imshow(I_phot, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax2.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax2.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)
	ax3.imshow(I_data1-I_phot, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax3.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax3.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	ax4.imshow(I_fit, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax4.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax4.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	ax5.imshow(I_data1-I_phot-I_fit, vmin=vmin, vmax=vmax, extent=ext_lens, **myargs)
	ax5.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax5.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)	
	ax6.imshow(I_source_sub, vmin=0.0, vmax=5.0*I_source_sub.std(), extent=ext_source, **myargs)
	ax6.set_xlabel(r'$\theta_x$ [arcsec]', fontsize=12)
	ax6.set_ylabel(r'$\theta_y$ [arcsec]', fontsize=12)
	plt.subplots_adjust(left=0.05, right=0.95, bottom=0.15, wspace=0.16)
	canvas_fit.draw()
	canvas_fit.get_tk_widget().focus_set()
	
def save_fit(event):
	save_fit_option=save_fit_options[save_fit_button.current()]
	if save_fit_option=='fits':
		fname=asksaveasfilename(defaultextension='.fits', initialdir=fpath, initialfile=obj_name+'_SIE.fits')
		if not fname:
			return
		col1=pyfits.Column(name='spar', format='D', array=spar)
		col1_err=pyfits.Column(name='spar_err', format='D', array=spar_err)
		col2=pyfits.Column(name='lpar', format='D', array=lpar)
		col2_err=pyfits.Column(name='lpar_err', format='D', array=lpar_err)
		col3=pyfits.Column(name='p_phot', format='D', array=ppar)
		col3_err=pyfits.Column(name='p_phot_err', format='D', array=ppar_err)
		tbhdu1=pyfits.BinTableHDU.from_columns([col1, col1_err])
		tbhdu2=pyfits.BinTableHDU.from_columns([col2, col2_err])
		tbhdu3=pyfits.BinTableHDU.from_columns([col3, col3_err])
		
		hdu=pyfits.PrimaryHDU()
		hdu.header['target']=obj_name
		hdu.header['N_Lens']=n_lens
		for i_lens in np.arange(n_lens):
			hdu.header['lmodel_'+str(i_lens+1)]=lensmodel[i_lens]
		hdu.header['N_Source']=n_source
		for i_source in np.arange(n_source):
			hdu.header['smodel_'+str(i_source+1)]=sourcemodel[i_source]
		hdu.header['N_Phot']=n_phot
		for i_phot in np.arange(n_phot):
			hdu.header['pmodel_'+str(i_phot+1)]=photmodel[i_phot]
		try:
			hdu.header['DOF']=dof
		except:
			hdu.header['DOF']=out.nfree
		try:
			hdu.header['Chi2']=chi2
		except:
			hdu.header['Chi2']=out.chisqr
		hdu.header['dpix_l']=dpix
		try:
			hdu.header['dpix_s']=dpix_source
		except:
			hdu.header['dpix_s']=0.1*dpix
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
		pyfits.append(fname, I_data1)
		pyfits.append(fname, I_invvar1*jmask1)
		pyfits.append(fname, I_phot)
		pyfits.append(fname, I_fit)
		pyfits.append(fname, I_source)
		if (len(sourcemodel) == 1 and sourcemodel[0] == 'pix'):
			hdu.header.add_comment('9th extension is the source error')
		 	pyfits.append(fname, I_source_err)
		print 'Saved to the file '+fname+'...'
	else:
		fname=asksaveasfilename(defaultextension='.eps', initialfile=obj_name+'_SIE.eps')
		if not fname:
			return
		fig1.savefig(fname)
		print 'Saved to the file '+fname+'...'

	save_fit_button_text.set('Export as')

def recover_a_lens_entry(row, column, ind):
# recover a lens entry from file
# the parameter entries are saved as a dictionary for each lens
	global current_lens_index
	current_lens_frame=Frame(lens_panel)
	current_lens_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	for i in range(30):
		current_lens_frame.grid_rowconfigure(i, weight=1)
		current_lens_frame.grid_columnconfigure(i, weight=1)
	lens_frame_list.append(current_lens_frame)
	current_drop_lens_button=Button(current_lens_frame, text="Drop", width=12, command=lambda: drop_a_lens(row))
	current_drop_lens_button.grid(row=0, column=column+12)
	drop_lens_button_list.append(current_drop_lens_button)
	if lensmodel[ind] == 'sie' or lensmodel[ind] == 'sple':
		Label(current_lens_frame, text="b_SIE", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		Label(current_lens_frame, text="lens_pa", width=12).grid(row=0, column=column+6)
		Label(current_lens_frame, text="lens_q", width=12).grid(row=0, column=column+8)
		Label(current_lens_frame, text="gamma", width=12).grid(row=0, column=column+10)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar3_entry=Entry(current_lens_frame, width=12)
		current_lpar4_entry=Entry(current_lens_frame, width=12)
		current_lpar5_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)
		current_lpar3_entry.grid(row=0, column=column+7)
		current_lpar4_entry.grid(row=0, column=column+9)
		current_lpar5_entry.grid(row=0, column=column+11)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)
		current_lpar3_entry.delete(0, END)
		current_lpar4_entry.delete(0, END)
		current_lpar5_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar1_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar2_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar3_entry.insert(0, "value=%16.10f, min=0.0, max=180." %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar4_entry.insert(0, "value=%16.10f, min=0.01, max=1.0" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar5_entry.insert(0, "value=%16.10f, vary=0" %lpar[current_lens_index])
		current_lens_index+=1
		lpar_entry_i={'b_SIE_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry, 'lens_pa_entry': current_lpar3_entry, 'lens_q_entry': current_lpar4_entry, 'gamma_entry': current_lpar5_entry}
		lpar_entry_list.append(lpar_entry_i)	
	elif lensmodel[ind] == 'pm':
		Label(current_lens_frame, text="theta_e", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=%16.10f, min=0.0" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar1_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar2_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		lpar_entry_i={'theta_e_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry}
		lpar_entry_list.append(lpar_entry_i)
	elif lensmodel[ind] == 'ext. shear':
		Label(current_lens_frame, text="shear_amp", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="shear_pa", width=12).grid(row=0, column=column+2)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=%16.10f, min=0.0" %lpar[current_lens_index])
		current_lens_index+=1		
		current_lpar1_entry.insert(0, "value=%16.10f, min=0.0, max=180." %lpar[current_lens_index])
		current_lens_index+=1
		
		lpar_entry_i={'shear_amp_entry':current_lpar0_entry, 'shear_pa_entry': current_lpar1_entry}
		lpar_entry_list.append(lpar_entry_i)
	elif lensmodel[ind] == 'snfw':
		Label(current_lens_frame, text="kappa_s", width=12).grid(row=0, column=column)
		Label(current_lens_frame, text="lens_xcen", width=12).grid(row=0, column=column+2)
		Label(current_lens_frame, text="lens_ycen", width=12).grid(row=0, column=column+4)
		Label(current_lens_frame, text="lens_pa", width=12).grid(row=0, column=column+6)
		Label(current_lens_frame, text="lens_q", width=12).grid(row=0, column=column+8)
		Label(current_lens_frame, text="r_s", width=12).grid(row=0, column=column+10)
		current_lpar0_entry=Entry(current_lens_frame, width=12)
		current_lpar1_entry=Entry(current_lens_frame, width=12)
		current_lpar2_entry=Entry(current_lens_frame, width=12)
		current_lpar3_entry=Entry(current_lens_frame, width=12)
		current_lpar4_entry=Entry(current_lens_frame, width=12)
		current_lpar5_entry=Entry(current_lens_frame, width=12)
		current_lpar0_entry.grid(row=0, column=column+1)
		current_lpar1_entry.grid(row=0, column=column+3)
		current_lpar2_entry.grid(row=0, column=column+5)
		current_lpar3_entry.grid(row=0, column=column+7)
		current_lpar4_entry.grid(row=0, column=column+9)
		current_lpar5_entry.grid(row=0, column=column+11)

		current_lpar0_entry.delete(0, END)
		current_lpar1_entry.delete(0, END)
		current_lpar2_entry.delete(0, END)
		current_lpar3_entry.delete(0, END)
		current_lpar4_entry.delete(0, END)
		current_lpar5_entry.delete(0, END)

		current_lpar0_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar1_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar2_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar3_entry.insert(0, "value=%16.10f, vary=0" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar4_entry.insert(0, "value=%16.10f, vary=0" %lpar[current_lens_index])
		current_lens_index+=1
		current_lpar5_entry.insert(0, "value=%16.10f" %lpar[current_lens_index])
		current_lens_index+=1
		lpar_entry_i={'kappa_s_entry':current_lpar0_entry, 'lens_xcen_entry': current_lpar1_entry, 'lens_ycen_entry': current_lpar2_entry, 'lens_pa_entry': current_lpar3_entry, 'lens_q_entry': current_lpar4_entry, 'r_s_entry': current_lpar5_entry}
		lpar_entry_list.append(lpar_entry_i)

def recover_a_source_entry(row, column, ind):
# recover a source entry from file
# the parameter entries are saved as a dictionary for each source
	global current_source_index
	current_source_frame=Frame(source_panel)
	current_source_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	for i in range(30):
		current_source_frame.grid_rowconfigure(i, weight=1)
		current_source_frame.grid_columnconfigure(i, weight=1)
	source_frame_list.append(current_source_frame)
	current_drop_source_button=Button(current_source_frame, text="Drop", width=12, command=lambda: drop_a_source(row))
	current_drop_source_button.grid(row=0, column=column+14)
	drop_source_button_list.append(current_drop_source_button)
	if sourcemodel[ind] == 'sersic' or sourcemodel[ind] == 'gaussian':
		Label(current_source_frame, text="src_amp", width=12).grid(row=0, column=column)
		Label(current_source_frame, text="src_xcen", width=12).grid(row=0, column=column+2)
		Label(current_source_frame, text="src_ycen", width=12).grid(row=0, column=column+4)
		Label(current_source_frame, text="src_sigma", width=12).grid(row=0, column=column+6)
		Label(current_source_frame, text="src_pa", width=12).grid(row=0, column=column+8)
		Label(current_source_frame, text="src_q", width=12).grid(row=0, column=column+10)
		Label(current_source_frame, text="src_n", width=12).grid(row=0, column=column+12)

		current_source_amp_entry=Entry(current_source_frame, width=12)
		current_source_xcen_entry=Entry(current_source_frame, width=12)
		current_source_ycen_entry=Entry(current_source_frame, width=12)
		current_source_sigma_entry=Entry(current_source_frame, width=12)
		current_source_pa_entry=Entry(current_source_frame, width=12)
		current_source_q_entry=Entry(current_source_frame, width=12)
		current_source_n_entry=Entry(current_source_frame, width=12)

		current_source_amp_entry.grid(row=0, column=column+1)
		current_source_xcen_entry.grid(row=0, column=column+3)
		current_source_ycen_entry.grid(row=0, column=column+5)
		current_source_sigma_entry.grid(row=0, column=column+7)
		current_source_pa_entry.grid(row=0, column=column+9)
		current_source_q_entry.grid(row=0, column=column+11)
		current_source_n_entry.grid(row=0, column=column+13)

		current_source_amp_entry.delete(0, END)
		current_source_xcen_entry.delete(0, END)
		current_source_ycen_entry.delete(0, END)
		current_source_sigma_entry.delete(0, END)
		current_source_pa_entry.delete(0, END)
		current_source_q_entry.delete(0, END)
		current_source_n_entry.delete(0, END)

		current_source_amp_entry.insert(0, "value=%16.10f, min=0.0" %spar[current_source_index])
		current_source_index+=1
		current_source_xcen_entry.insert(0, "value=%16.10f" %spar[current_source_index])
		current_source_index+=1
		current_source_ycen_entry.insert(0, "value=%16.10f" %spar[current_source_index])
		current_source_index+=1
		current_source_sigma_entry.insert(0, "value=%16.10f" %spar[current_source_index])
		current_source_index+=1
		current_source_pa_entry.insert(0, "value=%16.10f, min=0.0, max=180.0" %spar[current_source_index])
		current_source_index+=1
		current_source_q_entry.insert(0, "value=%16.10f, min=0.0, max=1.0" %spar[current_source_index])
		current_source_index+=1
		current_source_n_entry.insert(0, "value=%16.10f, min=0.01" %spar[current_source_index])
		current_source_index+=1

		spar_entry_i={'source_amp_entry': current_source_amp_entry, 'source_xcen_entry': current_source_xcen_entry, 'source_ycen_entry': current_source_ycen_entry, 'source_sigma_entry': current_source_sigma_entry, 'source_pa_entry': current_source_pa_entry, 'source_q_entry': current_source_q_entry, 'source_n_entry': current_source_n_entry}
		spar_entry_list.append(spar_entry_i)
	elif sourcemodel[ind] == 'pix':
		Label(current_source_frame, text="lambda", width=12).grid(row=0, column=column)
		Label(current_source_frame, text="reg. scheme", width=12).grid(row=0, column=column+2)
		current_source_lambda_entry=Entry(current_source_frame, width=12)
		current_source_reg_scheme_entry=Entry(current_source_frame, width=12)
		current_source_lambda_entry.grid(row=0, column=column+1)
		current_source_reg_scheme_entry.grid(row=0, column=column+3)
		current_source_lambda_entry.delete(0, END)
		current_source_reg_scheme_entry.delete(0, END)
		current_source_lambda_entry.insert(0, "value=%16.10f" %spar[current_source_index])
		current_source_index+=1		
		current_source_reg_scheme_entry.insert(0, "value=%16.10f, vary=0" %spar[current_source_index])
		current_source_index+=1		
		spar_entry_i={'lambda_entry': current_source_lambda_entry, 'reg_scheme_entry': current_source_reg_scheme_entry}
		spar_entry_list.append(spar_entry_i)
		
def recover_a_phot_entry(row, column, ind):
# recover a phot entry from file
# the parameter entries are saved as a dictionary for each phot
	global current_phot_index
	current_phot_frame=Frame(phot_panel)
	current_phot_frame.grid(row=row, column=column, rowspan=1, sticky='ew')
	for i in range(30):
		current_phot_frame.grid_rowconfigure(i, weight=1)
		current_phot_frame.grid_columnconfigure(i, weight=1)
	phot_frame_list.append(current_phot_frame)
	current_drop_phot_button=Button(current_phot_frame, text="Drop", width=12, command=lambda: drop_a_phot(row))
	current_drop_phot_button.grid(row=0, column=column+22)
	drop_phot_button_list.append(current_drop_phot_button)
	if photmodel[ind] == 'sersic':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="sigma", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="n", width=12).grid(row=0, column=column+12)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+14)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar7_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)
		current_ppar7_entry.grid(row=0, column=column+15)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)
		current_ppar7_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar1_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar2_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar3_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar4_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar5_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar6_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar7_entry.insert(0, "value=%16.10f, min=-0.02, max=0.02" %ppar[current_phot_index])
		current_phot_index+=1
		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_sigma_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_n_entry': current_ppar6_entry, 'phot_I0_entry': current_ppar7_entry}
		ppar_entry_list.append(ppar_entry_i)
	elif photmodel[ind] == 'csersic':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="sigma", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="n", width=12).grid(row=0, column=column+12)
		Label(current_phot_frame, text="rc", width=12).grid(row=0, column=column+14)
		Label(current_phot_frame, text="alpha", width=12).grid(row=0, column=column+16)
		Label(current_phot_frame, text="gamma", width=12).grid(row=0, column=column+18)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+20)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar7_entry=Entry(current_phot_frame, width=12)
		current_ppar8_entry=Entry(current_phot_frame, width=12)
		current_ppar9_entry=Entry(current_phot_frame, width=12)
		current_ppar10_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)
		current_ppar7_entry.grid(row=0, column=column+15)
		current_ppar8_entry.grid(row=0, column=column+17)
		current_ppar9_entry.grid(row=0, column=column+19)
		current_ppar10_entry.grid(row=0, column=column+21)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)
		current_ppar7_entry.delete(0, END)
		current_ppar8_entry.delete(0, END)
		current_ppar9_entry.delete(0, END)
		current_ppar10_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar1_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar2_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar3_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar4_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar5_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar6_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar7_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar8_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar9_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar10_entry.insert(0, "value=%16.10f, min=-0.02, max=0.02" %ppar[current_phot_index])
		current_phot_index+=1

		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_sigma_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_n_entry': current_ppar6_entry, 'phot_rc_entry': current_ppar7_entry, 'phot_alpha_entry': current_ppar8_entry, 'phot_gamma_entry': current_ppar9_entry, 'phot_I0_entry': current_ppar10_entry}
		ppar_entry_list.append(ppar_entry_i)
	elif photmodel[ind] == 'hernquist':
		Label(current_phot_frame, text="amp", width=12).grid(row=0, column=column)
		Label(current_phot_frame, text="xcen", width=12).grid(row=0, column=column+2)
		Label(current_phot_frame, text="ycen", width=12).grid(row=0, column=column+4)
		Label(current_phot_frame, text="rs", width=12).grid(row=0, column=column+6)
		Label(current_phot_frame, text="pa", width=12).grid(row=0, column=column+8)
		Label(current_phot_frame, text="q", width=12).grid(row=0, column=column+10)
		Label(current_phot_frame, text="bg", width=12).grid(row=0, column=column+12)
		current_ppar0_entry=Entry(current_phot_frame, width=12)
		current_ppar1_entry=Entry(current_phot_frame, width=12)
		current_ppar2_entry=Entry(current_phot_frame, width=12)
		current_ppar3_entry=Entry(current_phot_frame, width=12)
		current_ppar4_entry=Entry(current_phot_frame, width=12)
		current_ppar5_entry=Entry(current_phot_frame, width=12)
		current_ppar6_entry=Entry(current_phot_frame, width=12)
		current_ppar0_entry.grid(row=0, column=column+1)
		current_ppar1_entry.grid(row=0, column=column+3)
		current_ppar2_entry.grid(row=0, column=column+5)
		current_ppar3_entry.grid(row=0, column=column+7)
		current_ppar4_entry.grid(row=0, column=column+9)
		current_ppar5_entry.grid(row=0, column=column+11)
		current_ppar6_entry.grid(row=0, column=column+13)

		current_ppar0_entry.delete(0, END)
		current_ppar1_entry.delete(0, END)
		current_ppar2_entry.delete(0, END)
		current_ppar3_entry.delete(0, END)
		current_ppar4_entry.delete(0, END)
		current_ppar5_entry.delete(0, END)
		current_ppar6_entry.delete(0, END)

		current_ppar0_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar1_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar2_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar3_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar4_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar5_entry.insert(0, "value=%16.10f" %ppar[current_phot_index])
		current_phot_index+=1
		current_ppar6_entry.insert(0, "value=%16.10f, min=-0.02, max=0.02" %ppar[current_phot_index])
		current_phot_index+=1
		ppar_entry_i={'phot_amp_entry':current_ppar0_entry, 'phot_xcen_entry': current_ppar1_entry, 'phot_ycen_entry': current_ppar2_entry, 'phot_rs_entry': current_ppar3_entry, 'phot_pa_entry': current_ppar4_entry, 'phot_q_entry': current_ppar5_entry, 'phot_I0_entry': current_ppar6_entry}
		ppar_entry_list.append(ppar_entry_i)

def import_lens_model():
	global n_lens, n_source, n_phot, lensmodel, sourcemodel, photmodel, lpar, spar, ppar, lpar_err, spar_err, ppar_err, dof, chi2, current_row_lens, current_row_source, current_row_phot, lens_model_label, source_model_label, phot_model_label, current_lens_index, current_source_index, current_phot_index, I_data1, I_image1, I_invvar1, I_phot, I_fit, I_source, hw, x1, y1, source_xbase, source_ybase, dpix, dpix_source

	popup=Toplevel()
	popup.wm_title("Import Lens Model")
	popup.protocol("WM_DELETE_WINDOW", lambda: popup.destroy())
	frame_import=Frame(popup)
	frame_import.pack()
	fname=askopenfilename(defaultextension='.fits', initialdir=fpath, initialfile=obj_name+'_SIE.fits')
	if fname:
		print 'Importing lens model from %s' %fname
	
		clear_parameter_panel()
		current_lens_index=0
		current_source_index=0
		current_phot_index=0
		
		hdulist=pyfits.open(fname)
		n_lens=hdulist[0].header['N_Lens']
		n_source=hdulist[0].header['N_Source']
		n_phot=hdulist[0].header['N_Phot']
		lensmodel=[]
		for i_lens in np.arange(n_lens):
			lensmodel.append(hdulist[0].header['lmodel_'+str(i_lens+1)])
		sourcemodel=[]
		for i_source in np.arange(n_source):
			sourcemodel.append(hdulist[0].header['smodel_'+str(i_source+1)])
		photmodel=[]
		for i_phot in np.arange(n_phot):
			photmodel.append(hdulist[0].header['pmodel_'+str(i_phot+1)])
		#dof=hdulist[0].header['DOF']
		chi2=hdulist[0].header['Chi2']
		dpix=hdulist[0].header['dpix_l']
		dpix_source=hdulist[0].header['dpix_s']
		spar=((hdulist[1].data)['spar']).copy()
		lpar=((hdulist[2].data)['lpar']).copy()
		ppar=((hdulist[3].data)['p_phot']).copy()
		spar_err=((hdulist[1].data)['spar_err']).copy()
		lpar_err=((hdulist[2].data)['lpar_err']).copy()
		ppar_err=((hdulist[3].data)['p_phot_err']).copy()
		I_data1=(hdulist[4].data).copy()
		I_invvar1=(hdulist[5].data).copy()
		I_phot=(hdulist[6].data).copy()
		I_fit=(hdulist[7].data).copy()
		I_source=(hdulist[8].data).copy()
		I_image1=I_data1-I_phot
		hdulist.close()
	
		n_x=I_data1.shape[0]
		n_y=I_data1.shape[1]
		hw=(n_x-1)/2
		x1=np.outer(np.ones(n_y), np.arange(n_x)-n_x/2)*dpix
		y1=np.outer(np.arange(n_y)-n_y/2, np.ones(n_x))*dpix
		n_x_source=I_source.shape[0]
		n_y_source=I_source.shape[1]
		source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2)*dpix_source
		source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2, np.ones(n_x_source))*dpix_source
		
		lens_model_label=Label(lens_panel, text='Lens Model')
		lens_model_label.config(relief=RAISED, bd=2)
		lens_model_label.grid(row=current_row_lens, column=0, sticky='w')
		current_row_lens=current_row_lens+1
		for i_lens in np.arange(n_lens):
			recover_a_lens_entry(current_row_lens, 0, i_lens)
			current_row_lens=current_row_lens+1
	
		source_model_label=Label(source_panel, text='Source Model')
		source_model_label.config(relief=RAISED, bd=2)
		source_model_label.grid(row=current_row_source, column=0, sticky='w')
		current_row_source=current_row_source+1
		for i_source in np.arange(n_source):
			recover_a_source_entry(current_row_source, 0, i_source)
			current_row_source=current_row_source+1
	
		if n_phot > 0:
			phot_model_label=Label(phot_panel, text='Phot. Model')
			phot_model_label.config(relief=RAISED, bd=2)
			phot_model_label.grid(row=current_row_phot, column=0, sticky='w')
			current_row_phot=current_row_phot+1
			for i_phot in np.arange(n_phot):
				recover_a_phot_entry(current_row_phot, 0, i_phot)
				current_row_phot=current_row_phot+1
	else:
		pass

