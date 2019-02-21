'''
GUI for using the command line mapping code I ahve written so far

InitMapper is opened first as a Tk GUI interface. This interface asks you to choose what code to run.
Based on what you select, certain fields are added to the next window that pops up.

The second window, the input frame, provides a run, show figure and save file button as defaults
'''


import sys, re, os

if sys.version_info[0] < 3:
    vers=2
    print("Using python version 2.x")
else:
    vers=3
    print("Using python version 3.x")

if vers == 2:
	from Tkinter import *
	import tkFileDialog as filedialog
elif vers == 3:
	from tkinter import *
	import tkinter.filedialog as filedialog

class Chooser:
	'''
	To add new code to the GUI:
	1) Add button for the mapping/plotting code in __init__
	2) Map the button click to a new method (keep the same naming convention as other buttons)
	3) Write the new method at the bottom of the other methods
		3a) Add widgets in the method based on what the code needs. If there is no widget for what 
		the code needs yet, write the new widget in the Input class (see notes for that class for 
		instructions on how to do this)
	'''
	def __init__(self, parent):
		self.begin = parent
		self.begin.title("Choose what to plot...")
		# Main frame
		self.frame = Frame(self.begin,bd=2,relief=GROOVE)
		self.frame.pack(pady=2,fill=X)
		# Control button and text
		self.DiffMap_Button = Button(self.frame, text="DiffMap", command = self.runDiffMap, height = 1)		
		self.DiffMap_Button.pack(anchor = "center")
		self.DiffZonalMean_Button = Button(self.frame, text="DiffZonalMean", command = self.runDiffZonalMean, height = 1)		
		self.DiffZonalMean_Button.pack(anchor = "center")
		self.DiffMeridionalMean_Button = Button(self.frame, text="DiffMeridionalMean", command = self.runDiffMeridionalMean, height = 1)		
		self.DiffMeridionalMean_Button.pack(anchor = "center")
		self.DiffPressureOverBox_Button = Button(self.frame, text="DiffPressureOverBox", command = self.runDiffPressureOverBox, height = 1)		
		self.DiffPressureOverBox_Button.pack(anchor = "center")
		self.DiffPressureVsLong_Button = Button(self.frame, text="DiffPressureVsLong", command = self.runDiffPressureVsLong, height = 1)		
		self.DiffPressureVsLong_Button.pack(anchor = "center")
		self.DiffPressureVsLat_Button = Button(self.frame, text="DiffPressureVsLat", command = self.runDiffPressureVsLat, height = 1)		
		self.DiffPressureVsLat_Button.pack(anchor = "center")
		self.DiffCorrelationMap_Button = Button(self.frame, text="DiffCorrelationMap", command = self.runDiffCorrelationMap, height = 1)		
		self.DiffCorrelationMap_Button.pack(anchor = "center")
		self.DiffEOFmap_Button = Button(self.frame, text="DiffEOFmap", command = self.runDiffEOFmap, height = 1)		
		self.DiffEOFmap_Button.pack(anchor = "center")
		self.CorrelationMap_Button = Button(self.frame, text="CorrelationMap", command = self.runCorrelationMap, height = 1)		
		self.CorrelationMap_Button.pack(anchor = "center")
		# Add a quit button to shut it all down
		self.bottomframe = Frame(self.begin,bd = 2, relief = GROOVE)
		self.bottomframe.pack(pady=2,fill=X)
		self.Quit_Button = Button(self.bottomframe, text = "QUIT", command = self.close_all_windows, height = 1)
		self.Quit_Button.pack(anchor = "center")
	#
	def close_all_windows(self):
		try:
			self.root2.destroy()
		except (AttributeError, TclError):
			pass
		root.destroy()
	#
	def runDiffMap(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffMap")
		window.AddChooseFile("c", "Select control file")
		window.AddChooseFile("t", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-40 40 335 333')
		window.AddText('v', 'Enter variable name', default = "PRECT")
		window.AddText('clev', 'Enter clev (min max nlev)', default = '')
		window.AddText('diffclev', 'Enter diffclev (min max nlev)', default = '')
		self.root2.mainloop()
	#
	def runDiffZonalMean(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffZonalMean")
		window.AddChooseFile("c", "Select control file")
		window.AddChooseFile("t", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "PRECT")
		self.root2.mainloop()
	#
	def runDiffMeridionalMean(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffMeridionalMean")
		window.AddChooseFile("control", "Select control file")
		window.AddChooseFile("test", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "PRECT")
		self.root2.mainloop()
	#
	def runDiffPressureOverBox(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffPressureOverBox")
		window.AddChooseFile("control", "Select control file")
		window.AddChooseFile("test", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "T")
		self.root2.mainloop()
	#
	def runDiffPressureVsLong(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffPressureVsLong")
		window.AddChooseFile("control", "Select control file")
		window.AddChooseFile("test", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "T")
		window.AddText('clev', 'Enter clev (min max nlev)', default = '')
		window.AddText('diffclev', 'Enter diffclev (min max nlev)', default = '')
		self.root2.mainloop()
	#
	def runDiffPressureVsLat(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffPressureVsLat")
		window.AddChooseFile("control", "Select control file")
		window.AddChooseFile("test", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "T")
		window.AddText('clev', 'Enter clev (min max nlev)', default = '')
		window.AddText('diffclev', 'Enter diffclev (min max nlev)', default = '')
		self.root2.mainloop()
	#
	def runDiffPressureVsLat(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffCorrelationMap")
		window.AddChooseFile("control", "Select control file")
		window.AddChooseFile("test", "Select test file")
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-90 90 0 360')
		window.AddText('v', 'Enter variable name', default = "T")
		window.AddText('clev', 'Enter clev (min max nlev)', default = '')
		window.AddText('diffclev', 'Enter diffclev (min max nlev)', default = '')
		self.root2.mainloop()
	
	def runDiffCorrelationMap(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffCorrelationMap")
		window.AddChooseDir("cdir", "Set control directory")
		window.AddChooseDir("tdir", "Set test directory")
		window.AddText('field_var', 'Enter global variable name', default = 'PRECT')
		window.AddText('loc_var', 'Enter local (point) variable name', default = "PRECT_d18O")
		window.AddText('latlon', 'Enter lat/lon coordinate of the point (2 numbers)', default = '0 0')
		window.AddText('del', 'Enter delta - the number of degrees to make a box around the lat/lon coordinate', default = '2')
		window.AddText('grep_pre', 'Enter expression for grep before date string', default = '*cam.h0')
		window.AddText('grep_post', 'Enter expression for grep after date string', default = '.nc')
		window.AddDate(years = True, months = True, days = False)
		self.root2.mainloop()

	def runDiffEOFmap(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "DiffEOFmap")
		window.AddChooseDir("cdir", "Set control directory")
		window.AddChooseDir("tdir", "Set test directory")
		window.AddText('v', 'Enter variable name', default = 'PRECT')
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-50 50 0 360')
		window.AddText('n', 'Enter number of EOFs to plot', default = '1')
		window.AddText('grep_pre', 'Enter expression for grep before date string', default = 'kn')
		window.AddText('grep_post', 'Enter expression for grep after date string', default = 'climo*nc')
		window.AddDate(years = True, months = True, days = False)
		self.root2.mainloop()
	def runCorrelationMap(self):
		self.root2 = Toplevel()
		window = Input(self.root2, "CorrelationMap")
		window.AddChooseDir("dir", "Set directory")
		window.AddText('field_var', 'Enter global variable name', default = 'PRECT')
		window.AddText('loc_var', 'Enter local (point) variable name', default = "PRECT_d18O")
		window.AddText('latlon', 'Enter lat/lon coordinate of the point (2 numbers)', default = '0 0')
		window.AddText('del', 'Enter delta - the number of degrees to make a box around the lat/lon coordinate', default = '2')
		window.AddText('box', 'Enter box bounds (bottom, top, left, right)', default = '-50 50 0 360')
		window.AddText('grep_pre', 'Enter expression for grep before date string', default = ' ')
		window.AddText('grep_post', 'Enter expression for grep after date string', default = 'climo*nc')
		window.AddDate(years = True, months = True, days = False)
		self.root2.mainloop()

class Input:
	'''
	How to add new widgets:
	1) Do not change __init__, it should be good to go
	2) Add the method for the new widget using the existing naming conventions
		2a) First line of the method should be widgets.append("widgetName")
	3) Add the appropriate behavior in the OnRunClick method for the new widget by using 
	the existing if statements as examples (e.g. "if "widgetName" in widgets:")
	4) Append the correct string to the execute string so that the code runs correctly
		4a) This is based on the argparse section of the new code being added.
	'''
	def __init__(self, parent, code):
		self.code = code+".py"
		self.begin = parent
		self.begin.title("Choose input data for " + code)
		# Add stuff to click on and type in
		self.topframe = Frame(self.begin,bd=2,relief=GROOVE)
		self.topframe.pack(pady=2,fill=X)
		# Frame for save and show figure checkbuttons
		self.saveshowframe = Frame(self.begin, bd = 2, relief = GROOVE)
		self.saveshowframe.pack(pady=2,fill=Y)
		self.savefigbool = IntVar() 
		self.savefigbool.set(False) # set to false by default
		self.Save_Checkbutton = Checkbutton(self.saveshowframe, text = "Save file", variable = self.savefigbool)
		self.Save_Checkbutton.pack(anchor = 'center')
		self.showfigbool = IntVar()
		self.showfigbool.set(True) # set to true by default
		self.Show_Checkbutton = Checkbutton(self.saveshowframe, text = "Show image", variable = self.showfigbool)
		self.Show_Checkbutton.pack(anchor = 'center')
		self.Show_Checkbutton.select() # Actually not necessary, but this just puts a check on the Show image box for reassurance
		# Frame for run button
		self.runframe = Frame(self.begin,bd = 2, relief = GROOVE)
		self.runframe.pack(pady=2,fill=Y)
		self.Run_Button = Button(self.runframe, text = "RUN", command = self.OnRunClick, height = 1)
		self.Run_Button.pack(anchor = "center",side = LEFT)
		self.Quit_Button = Button(self.runframe, text = "QUIT", command = self.close_all_windows, height = 1)
		self.Quit_Button.pack(anchor = "center",side = LEFT)
		# Set a bunch of empty dicts to keep track of what's been added (is there a better way to do this? Probably...)
		self.doDates = False
		self.widgets = []
		self.Buttons = {}
		self.Labels = {}
		self.Texts = {}

	def AddChooseFile(self, name, label):
		self.widgets.append(name)
		self.Buttons[name] = Button(self.topframe, text=label, command = lambda: self.OnFileChooseClick(name, self.Texts[name]), height = 1)		
		self.Buttons[name].pack(anchor = "center")
		self.Texts[name] = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.Texts[name].pack(anchor = "center")

	def AddChooseDir(self, name, label):
		self.widgets.append(name)
		self.Buttons[name] = Button(self.topframe, text=label, command = lambda: self.OnDirChooseClick(name, self.Texts[name]), height = 1)		
		self.Buttons[name].pack(anchor = "center")
		self.Texts[name] = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.Texts[name].pack(anchor = "center")

	def AddText(self, name, label, default=''):
		self.widgets.append(name)
		self.Labels[name] = Label(self.topframe, text = label,height=1)
		self.Labels[name].pack(pady=4, anchor = "center")
		self.Texts[name] = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.Texts[name].insert(END, default)
		self.Texts[name].pack(anchor = "center")

	def AddDate(self, years = True, months = True, days = True):
		self.doDates = True
		self.datesToAdd = []
		if years:
			self.datesToAdd.append('years')
			self.Labels['years'] = Label(self.topframe, text = 'Set the years (2 numbers: start_year end_year) (enter -1 -1 for no years)', height=1)
			self.Labels['years'].pack(anchor = "center")
			self.Texts['years'] = Text(self.topframe,height=1,bd=2,relief=RIDGE)
			self.Texts['years'].pack(anchor = "center")
		if months:
			self.datesToAdd.append('months')
			self.MonthCheckButtons = {}
			self.monthframe = Frame(self.topframe, bd = 2, relief = GROOVE)
			self.monthframe.pack(pady=2,fill=Y)
			months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
			self.monthbools = []
			for i, month in enumerate(months):
				self.monthbools.append(IntVar())
				self.MonthCheckButtons[i] = Checkbutton(self.monthframe, text = month, variable = self.monthbools[i])
				self.MonthCheckButtons[i].pack(anchor = 'center', side = LEFT)
				self.MonthCheckButtons[i].select()
			self.AllCheckbutton = Checkbutton(self.topframe, text = "ANN", command = lambda: checkMonths(range(12)))
			self.AllCheckbutton.pack(side = LEFT)
			self.JJACheckbutton = Checkbutton(self.topframe, text = "JJA", command = lambda: checkMonths([5,6,7]))
			self.JJACheckbutton.pack(side = LEFT)
			self.DJFCheckbutton = Checkbutton(self.topframe, text = "DJF", command = lambda: checkMonths([0,1,11]))
			self.DJFCheckbutton.pack(side = LEFT)
			self.JJASCheckbutton = Checkbutton(self.topframe, text = "JJAS", command = lambda: checkMonths([5,6,7,8]))
			self.JJASCheckbutton.pack(side = LEFT)
			def checkMonths(months = range(12)): # defaults all months
				for i in range(12):
					self.MonthCheckButtons[i].deselect()
				for i in months:
					self.MonthCheckButtons[i].select()
		if days:
			self.widgets.append('days')


	def close_all_windows(self):
		try:
			self.begin.destroy()
		except (AttributeError, TclError):
			None
		root.destroy()
	#
	def OnFileChooseClick(self, name, textbox):
		title = "Select file for -" + name
		dlg=filedialog.Open(title = title)
		filepath = dlg.show()
		textbox.delete(1.0, END)
		textbox.insert(END, filepath)
	
	def OnDirChooseClick(self, name, textbox):
		title = "Select directory for -" + name
		dlg=filedialog.askdirectory(title = title)
		textbox.delete(1.0, END)
		textbox.insert(END, dlg)

	def OnRunClick(self):
		# Build the execute string line by line based on what info is provided
		# Not the prettiest way to do things
		# stackoverflow would really make fun of me but it works I guess
		pythonloc = "python" # Forces user to set shell env variable python to the correct interpreter
		codeloc = os.path.expanduser("~/python/bin/pygoda/") + self.code # Works for cheyenne and climate
		execute = pythonloc + " " + codeloc
		for widget in self.widgets:
			execute += ' -' + widget + " " + self.Texts[widget].get("1.0",'end-1c')
		if not self.savefigbool.get():
			execute += " -nosave "
		if self.showfigbool.get():
			execute += " -show "
		if self.doDates:
			if 'years' in self.datesToAdd:
				execute += " -years " + str(self.Texts['years'].get("1.0",'end-1c'))
			if 'months' in self.datesToAdd:
				execute += " -months "
				months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
				for i, month in enumerate(months):
					if self.monthbools[i].get():
						execute += " " + month
		print(execute)
		os.system(execute)


root = Tk()
first = Chooser(root)
root.mainloop()

