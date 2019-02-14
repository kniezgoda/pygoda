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
		# Add a quit button to shut it all down
		self.bottomframe = Frame(self.begin,bd = 2, relief = GROOVE)
		self.bottomframe.pack(pady=2,fill=X)
		self.Quit_Button = Button(self.bottomframe, text = "QUIT", command = self.close_all_windows, height = 1)
		self.Quit_Button.pack(anchor = "center")
	#
	def close_all_windows(self):
		try:
			self.root2.destroy()
		except NameError:
			None
		root.destroy()
	#
	def runDiffMap(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffMap")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-40 40 335 333')
		window.AddVariable("PRECT")
		window.AddClev()
		window.AddDiffClev()
		self.root2.mainloop()
	#
	def runDiffZonalMean(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffZonalMean")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-90 90 0 360')
		window.AddVariable("PRECT")
		self.root2.mainloop()
	#
	def runDiffMeridionalMean(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffMeridionalMean")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-90 90 0 360')
		window.AddVariable("PRECT")
		self.root2.mainloop()
	#
	def runDiffPressureOverBox(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffPressureOverBox")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-90 90 0 360')
		window.AddVariable("T")
		self.root2.mainloop()
	#
	def runDiffPressureVsLong(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffPressureVsLong")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-90 90 0 360')
		window.AddVariable("T")
		window.AddClev()
		window.AddDiffClev()
		self.root2.mainloop()
	#
	def runDiffPressureVsLat(self):
		self.root2 = Tk()
		window = Input(self.root2, "DiffPressureVsLat")
		window.AddControlFile()
		window.AddTestFile()
		window.AddBox('-90 90 0 360')
		window.AddVariable("T")
		window.AddClev()
		window.AddDiffClev()
		self.root2.mainloop()
	#

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
		self.savefigbool = BooleanVar() 
		self.savefigbool.set(False) # set to false by default
		self.Save_Checkbutton = Checkbutton(self.saveshowframe, text = "Save file", variable = self.savefigbool)
		self.Save_Checkbutton.pack(anchor = 'center')
		self.showfigbool = BooleanVar()
		self.showfigbool.set(True) # set to true by default
		self.Show_Checkbutton = Checkbutton(self.saveshowframe, text = "Show file", variable = self.showfigbool)
		self.Show_Checkbutton.pack(anchor = 'center')
		self.Show_Checkbutton.select()
		# Frame for run button
		self.runframe = Frame(self.begin,bd = 2, relief = GROOVE)
		self.runframe.pack(pady=2,fill=Y)
		self.Run_Button = Button(self.runframe, text = "RUN", command = self.OnRunClick, height = 1)
		self.Run_Button.pack(anchor = "center",side = LEFT)
		self.Quit_Button = Button(self.runframe, text = "QUIT", command = self.close_all_windows, height = 1)
		self.Quit_Button.pack(anchor = "center",side = LEFT)
		self.widgets = []
	#
	def AddControlFile(self):
		self.widgets.append("control")
		self.ControlFileChoose_Button = Button(self.topframe, text="Choose control file", command = lambda: self.OnFileChooseClick(0), height = 1)		
		self.ControlFileChoose_Button.pack(anchor = "center")
		self.ControlFileChoose_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.ControlFileChoose_Text.pack(anchor = "center")
	#
	def AddTestFile(self):
		self.widgets.append("test")
		self.TestFileChoose_Button = Button(self.topframe, text="Choose test file", command = lambda: self.OnFileChooseClick(1), height = 1)		
		self.TestFileChoose_Button.pack(anchor = "center")
		self.TestFileChoose_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.TestFileChoose_Text.pack(anchor = "center")
	#
	def AddBox(self, box):
		self.widgets.append("box")
		self.Box_Label = Label(self.topframe, text = "Box (bottom, top, left, right):",height=1)
		self.Box_Label.pack(pady=4, anchor = "center")
		self.Box_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.Box_Text.insert(END, box)
		self.Box_Text.pack(anchor = "center")
	#
	def AddVariable(self, variable):
		self.widgets.append("variable")
		self.Variable_Label = Label(self.topframe, text = "Variable:",height=1)
		self.Variable_Label.pack(pady=4, anchor = "center")
		self.Variable_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.Variable_Text.insert(END, variable)
		self.Variable_Text.pack(anchor = "center")
	#
	def AddClev(self):
		self.widgets.append("clev")
		self.ColorLevel_Label = Label(self.topframe, text = "Color level for control and test plots (min max nlev)\nLeave blank for default:",height=2)
		self.ColorLevel_Label.pack(pady=4, anchor = "center")
		self.ColorLevel_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.ColorLevel_Text.pack(anchor = "center")
	#
	def AddDiffClev(self):
		self.widgets.append("diffclev")
		self.DiffColorLevel_Label = Label(self.topframe, text = "Color level for difference plot (min max nlev)\nLeave blank for default:",height=2)
		self.DiffColorLevel_Label.pack(pady=4, anchor = "center")
		self.DiffColorLevel_Text = Text(self.topframe,height=1,bd=2,relief=RIDGE)
		self.DiffColorLevel_Text.pack(anchor = "center")
	#
	def close_all_windows(self):
		try:
			self.begin.destroy()
		except NameError:
			None
		root.destroy()
	#
	def OnFileChooseClick(self, which):
		title = "Select"
		if which == 0: 
			title += " control file"
		if which == 1: 
			title += " test file"
		dlg=filedialog.Open(title = title)
		#assign the path name of the selected file to temp
		temp = dlg.show()
		if temp is '':
			return None
		#assign the file name only to hold
		hold=re.split('/',temp)
		hold=hold[len(hold)-1]
		#insert text to the next avaiable input text box 
		if which == 0:
			self.control_filepath = temp
			self.control_filename = hold
			self.ControlFileChoose_Text.delete(1.0, END)
			self.ControlFileChoose_Text.insert(END, self.control_filename)
		if which == 1:
			self.test_filepath = temp
			self.test_filename = hold
			self.TestFileChoose_Text.delete(1.0, END)
			self.TestFileChoose_Text.insert(END, self.test_filename)
	#
	def OnRunClick(self):
		# Build the execute string line by line based on what info is provided
		# Not the prettiest way to do things
		# stackoverflow would really make fun of me but it works I guess
		pythonloc = "python" # Forces user to set shell env variable python to the correct interpreter
		codeloc = os.path.expanduser("~/python/bin/pygoda/") + self.code # Works for cheyenne and climate
		execute = pythonloc + " " + codeloc
		if "control" in self.widgets:
			execute += " -c " + self.control_filepath
		if "test" in self.widgets:
			execute += " -t " + self.test_filepath
		if "box" in self.widgets:
			getbox = self.Box_Text.get("1.0",'end-1c') # end-1c deletes the newline char at the end of text boxes
			self.box = getbox
			execute += " -box " + self.box
		if "variable" in self.widgets:
			self.variable = self.Variable_Text.get("1.0", "end-1c")
			execute += " -v " + self.variable
		if "clev" in self.widgets:
			self.clev = self.ColorLevel_Text.get("1.0", "end-1c")
			if self.clev is not '':
				execute += " -clev " + self.clev
		if "diffclev" in self.widgets:
			self.diffclev = self.DiffColorLevel_Text.get("1.0", "end-1c")
			if self.diffclev is not '':
				execute += " -diffclev " + self.diffclev
		self.savefig = self.savefigbool.get()
		self.showfig = self.showfigbool.get()
		if not self.savefig:
			execute += " -nosave "
		if self.showfig:
			execute += " -show "
		print(execute)
		os.system(execute)


root = Tk()
first = Chooser(root)
root.mainloop()


