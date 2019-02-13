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
	
	
class InitMapper:
	def __init__(self, parent):
		self.begin = parent
		self.begin.title("Choose what to plot...")
		# Main frame
		self.frame = Frame(self.begin,bd=2,relief=GROOVE)
		self.frame.pack(pady=2,fill=X)
		# Code buttons
		# DiffMap
		self.DiffMap_Button = Button(self.frame, text="DiffMap", command = self.runDiffMap, height = 1)		
		self.DiffMap_Button.pack(anchor = "center")
		# DiffZonalMean
		self.DiffZonalMean_Button = Button(self.frame, text="DiffZonalMean", command = self.runDiffZonalMean, height = 1)		
		self.DiffZonalMean_Button.pack(anchor = "center")
		
	def runDiffMap(self):
		root2 = Tk()
		diffmap = DiffMap(root2)
		root2.mainloop()
		
	def runDiffZonalMean(self):
		root2 = Tk()
		diffzonalmean = DiffZonalMean(root2)
		root2.mainloop()
		
class DiffZonalMean:
	def __init__(self, parent, control_file = '', test_file = '', box = [-40,40,335,333], var = "PRECT", clev = '', diffclev = ''):
		self.begin = parent
		self.begin.title("Choose input data...")
		# Things that need defaults
		self.control_filepath = control_file
		self.control_filename = ''
		if self.control_filepath is not '':
			self.control_filename = self.control_filepath.split("/")[-1]
		self.test_filepath = test_file
		self.test_filename = ''
		if self.test_filepath is not '':
			self.test_filename = self.test_filepath.split("/")[-1]
		self.box = str(box[0]) + " " + str(box[1]) + " " + str(box[2]) + " " + str(box[3])
		self.variable = var
		self.clev = clev
		self.diffclev = diffclev
		# Add stuff to click on and type in
		self.frame = Frame(self.begin,bd=2,relief=GROOVE)
		self.frame.pack(pady=2,fill=X)
		# Control button and text
		self.ControlFileChoose_Button = Button(self.frame, text="Choose control file", command = lambda: self.OnFileChooseClick(0), height = 1)		
		self.ControlFileChoose_Button.pack(anchor = "center")
		self.ControlFileChoose_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.ControlFileChoose_Text.insert(END, self.control_filename)
		self.ControlFileChoose_Text.pack(anchor = "center")
		# Test button and text
		self.TestFileChoose_Button = Button(self.frame, text="Choose test file", command = lambda: self.OnFileChooseClick(1), height = 1)		
		self.TestFileChoose_Button.pack(anchor = "center")
		self.TestFileChoose_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.TestFileChoose_Text.insert(END, self.test_filename)
		self.TestFileChoose_Text.pack(anchor = "center")
		# Box text
		self.Box_Label = Label(self.frame, text = "Box (bottom, top, left, right):",height=1)
		self.Box_Label.pack(pady=4, anchor = "center")
		self.Box_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.Box_Text.insert(END, self.box)
		self.Box_Text.pack(anchor = "center")
		# Variable text
		self.Variable_Label = Label(self.frame, text = "Variable:",height=1)
		self.Variable_Label.pack(pady=4, anchor = "center")
		self.Variable_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.Variable_Text.insert(END, self.variable)
		self.Variable_Text.pack(anchor = "center")
		# clev text
		self.ColorLevel_Label = Label(self.frame, text = "Color level for control and test plots (min max nlev)\nLeave blank for default:",height=2)
		self.ColorLevel_Label.pack(pady=4, anchor = "center")
		self.ColorLevel_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.ColorLevel_Text.insert(END, self.clev)
		self.ColorLevel_Text.pack(anchor = "center")
		# diffclev text
		self.DiffColorLevel_Label = Label(self.frame, text = "Color level for difference plot (min max nlev)\nLeave blank for default:",height=2)
		self.DiffColorLevel_Label.pack(pady=4, anchor = "center")
		self.DiffColorLevel_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.DiffColorLevel_Text.insert(END, self.diffclev)
		self.DiffColorLevel_Text.pack(anchor = "center")
		# Save and show check boxes
		self.savefigbool = BooleanVar() 
		self.savefigbool.set(False) # set to false by default
		self.showfigbool = BooleanVar()
		self.showfigbool.set(True) # set to true by default
		self.Save_Checkbutton = Checkbutton(self.frame, text = "Save file", variable = self.savefigbool)
		self.Save_Checkbutton.pack(anchor = 'center')
		self.Show_Checkbutton = Checkbutton(self.frame, text = "Show file", variable = self.showfigbool)
		self.Show_Checkbutton.pack(anchor = 'center')
		# Run button
		self.Run_Button = Button(self.frame, text = "RUN", command = self.runDiffMap, height = 1)
		self.Run_Button.pack(anchor = "center")
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
	def runDiffMap(self):
		self.box = self.Box_Text.get("1.0",'end-1c') # end-1c deletes the newline char at the end of text boxes
		self.variable = self.Variable_Text.get("1.0", "end-1c")
		self.clev = self.ColorLevel_Text.get("1.0", "end-1c")
		self.diffclev = self.DiffColorLevel_Text.get("1.0", "end-1c")
		self.savefig = self.savefigbool.get()
		self.showfig = self.showfigbool.get()
		# Build the string line by line to pass into os.system()
		pythonloc = '/home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python'
		codeloc = '/home/server/student/homes/kniezgod/python/bin/pygoda/DiffMap.py'
		execute = pythonloc + " " + codeloc
		execute += " -c " + self.control_filepath
		execute += " -t " + self.test_filepath
		execute += " -box " + self.box
		execute += " -v " + self.variable
		if self.clev is not '':
			execute += " -clev " + self.clev
		if self.diffclev is not '':
			execute += " -diffclev " + self.diffclev
		if not self.savefig:
			execute += " -nosave "
		if self.showfig:
			execute += " -show "
		print(execute)
		os.system(execute)

class DiffMap:
	def __init__(self, parent, control_file = '', test_file = '', box = [-40,40,335,333], var = "PRECT", clev = '', diffclev = ''):
		self.begin = parent
		self.begin.title("Choose input data...")
		# Things that need defaults
		self.control_filepath = control_file
		self.control_filename = ''
		if self.control_filepath is not '':
			self.control_filename = self.control_filepath.split("/")[-1]
		self.test_filepath = test_file
		self.test_filename = ''
		if self.test_filepath is not '':
			self.test_filename = self.test_filepath.split("/")[-1]
		self.box = str(box[0]) + " " + str(box[1]) + " " + str(box[2]) + " " + str(box[3])
		self.variable = var
		self.clev = clev
		self.diffclev = diffclev
		# Add stuff to click on and type in
		self.frame = Frame(self.begin,bd=2,relief=GROOVE)
		self.frame.pack(pady=2,fill=X)
		# Control button and text
		self.ControlFileChoose_Button = Button(self.frame, text="Choose control file", command = lambda: self.OnFileChooseClick(0), height = 1)		
		self.ControlFileChoose_Button.pack(anchor = "center")
		self.ControlFileChoose_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.ControlFileChoose_Text.insert(END, self.control_filename)
		self.ControlFileChoose_Text.pack(anchor = "center")
		# Test button and text
		self.TestFileChoose_Button = Button(self.frame, text="Choose test file", command = lambda: self.OnFileChooseClick(1), height = 1)		
		self.TestFileChoose_Button.pack(anchor = "center")
		self.TestFileChoose_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.TestFileChoose_Text.insert(END, self.test_filename)
		self.TestFileChoose_Text.pack(anchor = "center")
		# Box text
		self.Box_Label = Label(self.frame, text = "Box (bottom, top, left, right):",height=1)
		self.Box_Label.pack(pady=4, anchor = "center")
		self.Box_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.Box_Text.insert(END, self.box)
		self.Box_Text.pack(anchor = "center")
		# Variable text
		self.Variable_Label = Label(self.frame, text = "Variable:",height=1)
		self.Variable_Label.pack(pady=4, anchor = "center")
		self.Variable_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.Variable_Text.insert(END, self.variable)
		self.Variable_Text.pack(anchor = "center")
		# clev text
		self.ColorLevel_Label = Label(self.frame, text = "Color level for control and test plots (min max nlev)\nLeave blank for default:",height=2)
		self.ColorLevel_Label.pack(pady=4, anchor = "center")
		self.ColorLevel_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.ColorLevel_Text.insert(END, self.clev)
		self.ColorLevel_Text.pack(anchor = "center")
		# diffclev text
		self.DiffColorLevel_Label = Label(self.frame, text = "Color level for difference plot (min max nlev)\nLeave blank for default:",height=2)
		self.DiffColorLevel_Label.pack(pady=4, anchor = "center")
		self.DiffColorLevel_Text = Text(self.frame,height=1,bd=2,relief=RIDGE)
		self.DiffColorLevel_Text.insert(END, self.diffclev)
		self.DiffColorLevel_Text.pack(anchor = "center")
		# Save and show check boxes
		self.savefigbool = BooleanVar() 
		self.savefigbool.set(False) # set to false by default
		self.showfigbool = BooleanVar()
		self.showfigbool.set(True) # set to true by default
		self.Save_Checkbutton = Checkbutton(self.frame, text = "Save file", variable = self.savefigbool)
		self.Save_Checkbutton.pack(anchor = 'center')
		self.Show_Checkbutton = Checkbutton(self.frame, text = "Show file", variable = self.showfigbool)
		self.Show_Checkbutton.pack(anchor = 'center')
		# Run button
		self.Run_Button = Button(self.frame, text = "RUN", command = self.runDiffMap, height = 1)
		self.Run_Button.pack(anchor = "center")
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
	def runDiffMap(self):
		self.box = self.Box_Text.get("1.0",'end-1c') # end-1c deletes the newline char at the end of text boxes
		self.variable = self.Variable_Text.get("1.0", "end-1c")
		self.clev = self.ColorLevel_Text.get("1.0", "end-1c")
		self.diffclev = self.DiffColorLevel_Text.get("1.0", "end-1c")
		self.savefig = self.savefigbool.get()
		self.showfig = self.showfigbool.get()
		# Build the string line by line to pass into os.system()
		pythonloc = '/home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python'
		codeloc = '/home/server/student/homes/kniezgod/python/bin/pygoda/DiffMap.py'
		execute = pythonloc + " " + codeloc
		execute += " -c " + self.control_filepath
		execute += " -t " + self.test_filepath
		execute += " -box " + self.box
		execute += " -v " + self.variable
		if self.clev is not '':
			execute += " -clev " + self.clev
		if self.diffclev is not '':
			execute += " -diffclev " + self.diffclev
		if not self.savefig:
			execute += " -nosave "
		if self.showfig:
			execute += " -show "
		print(execute)
		os.system(execute)



root = Tk()
first = InitMapper(root)
root.mainloop()


