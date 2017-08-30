#!/usr/bin/env python
import os, commands, pyfits
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from Tkinter import *
import ttk
import lfit_gui_funcs as lfit
from numpy import arange

class Application:
	master=Tk()
	master.wm_title("Lens fitting tool")
	master.geometry("1600x800+10+10")
	for i in range(50):
	#	master.rowconfigure(i, weight=1)
	#	master.columnconfigure(i, weight=1)
		master.grid_rowconfigure(i, weight=1)
		master.grid_columnconfigure(i, weight=1)
	master.resizable(width=True, height=True)
	lfit.master=master
	
	main_panel=PanedWindow(master)
	main_panel.pack(fill=BOTH, expand=1)
	for i in range(50):
		main_panel.grid_rowconfigure(i, weight=1)
		main_panel.grid_columnconfigure(i, weight=1)
	
	left_panel=PanedWindow(main_panel, orient=VERTICAL)
	main_panel.add(left_panel)
	left_panel.config(width=350)
	for i in range(50):
		left_panel.grid_rowconfigure(i, weight=1)
		left_panel.grid_columnconfigure(i, weight=1)

	target_panel=Frame(left_panel, relief=SUNKEN, borderwidth=1)
	left_panel.add(target_panel)
	target_panel.config(height=600)
	for i in range(20):
		target_panel.grid_rowconfigure(i, weight=1)
		target_panel.grid_columnconfigure(i, weight=1)

	comment_panel=Frame(left_panel, relief=SUNKEN, borderwidth=1)
	left_panel.add(comment_panel)
	comment_panel.config(height=200)
	for i in range(20):
		comment_panel.grid_rowconfigure(i, weight=1)
		comment_panel.grid_columnconfigure(i, weight=1)

	comment_text=Text(comment_panel, borderwidth=1, relief=SUNKEN)
	comment_text.grid(row=1, column=0, columnspan=100)

	Button(comment_panel, text='Save comments', width=15, command=lfit.update_comments).grid(row=100, column=99)

	right_panel=PanedWindow(main_panel, orient=VERTICAL)
	main_panel.add(right_panel)
	for i in range(20):
		right_panel.grid_rowconfigure(i, weight=1)
		right_panel.grid_columnconfigure(i, weight=1)
	
	control_panel=Frame(right_panel, relief=SUNKEN, borderwidth=1)
	right_panel.add(control_panel)
	control_panel.config(height=200)
	for i in range(20):
		control_panel.grid_rowconfigure(i, weight=1)
		control_panel.grid_columnconfigure(i, weight=1)
	
	phot_panel=Frame(right_panel, relief=SUNKEN, borderwidth=1)
	right_panel.add(phot_panel)
	phot_panel.config(height=250)
	for i in range(20):
		phot_panel.grid_rowconfigure(i, weight=1)
		phot_panel.grid_columnconfigure(i, weight=1)
	
	lens_panel=Frame(right_panel, relief=SUNKEN, borderwidth=1)
	right_panel.add(lens_panel)
	lens_panel.config(height=250)
	for i in range(20):
		lens_panel.grid_rowconfigure(i, weight=1)
		lens_panel.grid_columnconfigure(i, weight=1)
	
	source_panel=Frame(right_panel, relief=SUNKEN, borderwidth=1)
	right_panel.add(source_panel)
	source_panel.config(height=300)
	for i in range(20):
		source_panel.grid_rowconfigure(i, weight=1)
		source_panel.grid_columnconfigure(i, weight=1)
	
	lfit.left_panel=left_panel
	lfit.target_panel=target_panel
	lfit.comment_panel=comment_panel
	lfit.comment_text=comment_text
	lfit.contro_panel=control_panel
	lfit.lens_panel=lens_panel
	lfit.phot_panel=phot_panel
	lfit.source_panel=source_panel
	
	# Find all avaiable cycles under $HST_DIR
	# Only needs to be done once when initiating the GUI
	cycle_list=lfit.get_cycle_list()
	
	cycle_label=Label(target_panel, text="Cycle", width=8)
	cycle_label.grid(row=0, column=0, padx=2, pady=2, sticky='nw')
	var_cycle=IntVar()
	cycle_menu=apply(OptionMenu, (target_panel, var_cycle)+tuple(cycle_list))
	cycle_menu.grid(row=0, column=1, padx=2, pady=2, sticky='nw')
	cycle_menu.config(width=15)
	var_cycle.set(23)
	lfit.var_cycle=var_cycle
	
	program_list=lfit.get_program_list(var_cycle.get())
	program_label=Label(target_panel, text="Program", width=8)
	program_label.grid(row=0, column=2, padx=2, pady=2, sticky='nw')
	var_program=IntVar()
	program_menu=apply(OptionMenu, (target_panel, var_program)+tuple(program_list))
	program_menu.grid(row=0, column=3, padx=2, pady=2, sticky='nw')
	program_menu.config(width=15)
	var_program.set(program_list[0])
	lfit.var_program=var_program
	lfit.program_menu=program_menu
	
	targetColumns=("visitID", "Target")
	target_tree=ttk.Treeview(target_panel, columns=targetColumns)
	for i in range(2):
		target_tree.column(i, anchor='center')
	targetColumn="#0"
	target_tree.heading(targetColumn, text="Number")
	target_tree.column(targetColumn, width=40)
	target_tree.column(targetColumns[0], width=40)
	target_tree.column(targetColumn[1], width=40)
	target_tree.heading(targetColumns[0], text="Visit #")
	target_tree.heading(targetColumns[1], text="Target")
	target_scrollbar=ttk.Scrollbar(target_panel, orient="vertical", style='blue.Vertical.TScrollbar', command=target_tree.yview)
	target_tree.configure(yscrollcommand=target_scrollbar.set)
	for i in range(2):
		target_tree.grid_rowconfigure(i, weight=1)
		target_tree.grid_columnconfigure(i, weight=1)
	target_tree.grid(column=0, row=1, columnspan=4, rowspan=25, padx=2, pady=2, sticky="news")
	target_tree.configure(height=30)
	target_scrollbar.grid(column=4, row=1, rowspan=25, sticky='ns')
	lfit.target_tree=target_tree
	
	visit_list, target_list=lfit.get_target_list(var_cycle.get(), var_program.get())
	for i in arange(len(visit_list)):
		target_tree.insert("", "end", text=str(i+1), values=(visit_list[i], target_list[i]))
	
	visit=visit_list[0]
	lfit.visit=visit
	
	target_tree.bind("<<TreeviewSelect>>", lfit.select_visit)
	
	n_lens=0
	lensmodel=[]
	current_row_lens=0
	current_column_lens=0
	lfit.n_lens=n_lens
	lfit.lensmodel=lensmodel
	lens_frame_list=[]
	drop_lens_button_list=[]
	lpar_entry_list=[]
	lpar_fix_list=[]
	n_lens_par_dict={'sie':6, 'sple':6, 'pm':3, 'snfw':6, 'ext. shear':2}
	
	lfit.current_row_lens=current_row_lens
	lfit.current_column_lens=current_column_lens
	lfit.lens_frame_list=lens_frame_list
	lfit.drop_lens_button_list=drop_lens_button_list
	lfit.lpar_entry_list=lpar_entry_list
	lfit.lpar_fix_list=lpar_fix_list
	lfit.n_lens_par_dict=n_lens_par_dict
	
	lensmodel_list=['sie', 'sple', 'pm', 'snfw', 'ext. shear']
	add_lens_model_text=StringVar()
	add_lens_model_text.set('Add a lens')
	add_lens_model_button=ttk.Combobox(control_panel, textvariable=add_lens_model_text, state='readonly')
	add_lens_model_button.grid(row=1, column=3, pady=4)
	add_lens_model_button['values']=tuple(lensmodel_list)
	add_lens_model_button.config(width=15)
	add_lens_model_button.bind("<<ComboboxSelected>>", lfit.add_a_lens)
	lfit.lensmodel_list=lensmodel_list
	lfit.add_lens_model_button=add_lens_model_button
	lfit.add_lens_model_text=add_lens_model_text
	
	n_phot=0
	photmodel=[]
	current_row_phot=0
	current_column_phot=0
	lfit.n_phot=n_phot
	lfit.photmodel=photmodel
	phot_frame_list=[]
	drop_phot_button_list=[]
	ppar_entry_list=[]
	ppar_fix_list=[]
	n_phot_par_dict={'sersic': 8, 'csersic': 11, 'hernquist': 7, 'bspline': 0}
	
	lfit.current_row_phot=current_row_phot
	lfit.current_column_phot=current_column_phot
	lfit.phot_frame_list=phot_frame_list
	lfit.drop_phot_button_list=drop_phot_button_list
	lfit.ppar_entry_list=ppar_entry_list
	lfit.ppar_fix_list=ppar_fix_list
	lfit.n_phot_par_dict=n_phot_par_dict
	
	photmodel_list=['sersic', 'csersic', 'hernquist', 'bspline']
	add_phot_model_text=StringVar()
	add_phot_model_text.set('Add a phot')
	add_phot_model_button=ttk.Combobox(control_panel, textvariable=add_phot_model_text, state='readonly')
	add_phot_model_button.grid(row=1, column=2, pady=4)
	add_phot_model_button['values']=tuple(photmodel_list)
	add_phot_model_button.config(width=15)
	add_phot_model_button.bind("<<ComboboxSelected>>", lfit.add_a_phot)
	lfit.photmodel_list=photmodel_list
	lfit.add_phot_model_button=add_phot_model_button
	lfit.add_phot_model_text=add_phot_model_text
	
	n_source=0
	sourcemodel=[]
	current_row_source=0
	current_column_source=0
	lfit.n_source=n_source
	lfit.sourcemodel=sourcemodel
	
	source_frame_list=[]
	drop_source_button_list=[]
	source_amp_entry_list=[]
	source_xcen_entry_list=[]
	source_ycen_entry_list=[]
	source_sigma_entry_list=[]
	source_pa_entry_list=[]
	source_q_entry_list=[]
	source_n_entry_list=[]
	source_amp_fix_list=[]
	source_xcen_fix_list=[]
	source_ycen_fix_list=[]
	source_sigma_fix_list=[]
	source_pa_fix_list=[]
	source_q_fix_list=[]
	source_n_fix_list=[]
	spar_entry_list=[]
	spar_fix_list=[]
	n_source_par_dict={'sersic': 7, 'gaussian': 7, 'pix':2}
	
	lfit.current_row_source=current_row_source
	lfit.current_column_source=current_column_source
	lfit.source_frame_list=source_frame_list
	lfit.drop_source_button_list=drop_source_button_list
	lfit.spar_entry_list=spar_entry_list
	lfit.spar_fix_list=spar_fix_list
	lfit.n_source_par_dict=n_source_par_dict
	
	sourcemodel_list=['sersic', 'pix', 'gaussian']
	add_source_model_text=StringVar()
	add_source_model_text.set('Add a source')
	add_source_model_button=ttk.Combobox(control_panel, textvariable=add_source_model_text, state='readonly')
	add_source_model_button.grid(row=1, column=4, pady=4)
	add_source_model_button['values']=tuple(sourcemodel_list)
	add_source_model_button.config(width=15)
	add_source_model_button.bind("<<ComboboxSelected>>", lfit.add_a_source)
	lfit.sourcemodel_list=sourcemodel_list
	lfit.add_source_model_button=add_source_model_button
	lfit.add_source_model_text=add_source_model_text
	
	# Need to initiate show_residual on the top level 
	# in order for that to work
	var_show_residual_bspline=IntVar()
	var_show_residual_sersic=IntVar()
	lfit.var_show_residual_bspline=var_show_residual_bspline
	lfit.var_show_residual_sersic=var_show_residual_sersic
	
	# Passing the functions to lfit
	
	lfit.getvalues()
	var_cycle.trace('w', lfit.update_program_list)
	var_program.trace('w', lfit.update_target_list)
	
	Button(control_panel, text='Quit', width=15, command=master.quit).grid(row=1, column=9, pady=4)
	Button(control_panel, text='Preview', width=15, command=lfit.show_image).grid(row=1, column=0, pady=4, sticky='nw')
	Button(control_panel, text='Import lens model from file', width=30, command=lfit.import_lens_model).grid(row=2, column=0, columnspan=3, pady=4, sticky='nw')
	Button(control_panel, text='L-M Phot. fit', width=15, command=lfit.phot_fit).grid(row=1, column=6, pady=4)
	Button(control_panel, text='MCMC Phot. fit', width=15, command=lfit.phot_fit_mcmc).grid(row=1, column=7, pady=4)
	Button(control_panel, text='Show Phot. fit', width=15, command=lfit.show_phot_fit).grid(row=1, column=8, pady=4)
	Button(control_panel, text='L-M Optimize', width=15, command=lfit.lfit).grid(row=2, column=6, pady=4)
	Button(control_panel, text='MCMC lfit', width=15, command=lfit.lfit_mcmc).grid(row=2, column=7, pady=4)
	Button(control_panel, text='MultiNest', width=15, command=lfit.lfit_multinest).grid(row=2, column=8, pady=4)
	Button(control_panel, text='Show fit', width=15, command=lfit.show_fit).grid(row=2, column=9, pady=4)
	
	master.protocol("WM_DELETE_WINDOW", lfit.master_on_closing)
	master.mainloop()	
