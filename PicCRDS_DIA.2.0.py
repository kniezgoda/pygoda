#!/usr/bin/env python -W ignore::DeprecationWarning

'''
PicCRDS_DIA: Picarro Cavity Ring-down Spectrometer Data Identifier Algorithm
EOS Version 2 (Feb 1, 2016):
    -Changed function of program:
        -now, the FirstWidget instance, named 'first' will stay open even after the LabelWidget and FigureWidget are closed
        -this allows for quick changing of files and parameters
    -In terms of code, nothing is done in the global body of the code
        -everything is now in 'run' method of FirstWidget or in the EnterStds method of LabelWidget
        -functions are defined below the classes, and the only global body code is the code that initiates the FirstWidget instance
    -Allows user to choose what is plotter (d18O or dD)
    -Drift correction techinque now uses inner third method


Author: Kyle Niezgoda
kniezgoda@ntu.edu.sg (work)
kniezgo@gmail.com (personal)


The following non-standard python packages are required to run this program:

matplotlib, pandas, numpy, pylab, tkinter/Tkinter, scipy

All of the above packages can be downloaded by installing Anaconda for Python.


This app is optimized for use on Mac operating systems, but will work on Windows and Linux as well.
Some slight changes in GUI aesthetics will occur, but the data calibration and identification methods do not change.

App is currently not available for python version 3.x on Windows using the Anaconda installer, as the way Anaconda 
accesses Windows Tcl/Tk software is distrupting some of the printing procedures in the code. 
Windows users are encouraged to install Python version 2.x in order to use this software. 


Data files from the CRDS may be in .dat or .csv format.
'''


#define python version
import sys, re, os
if sys.version_info[0] < 3:
    vers=2
else:
    vers=3

import matplotlib
matplotlib.use('TkAgg')
import pandas as pd
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
if vers == 2:
    from Tkinter import *
    import tkFileDialog as filedialog
elif vers == 3:
    from tkinter import *
    import tkinter.filedialog as filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
from os.path import expanduser
from math import *
from scipy import stats
import time
import datetime


if vers==2:
    print("Using python version 2.x")
if vers==3:
    print("Using python version 3.x")

home=os.path.expanduser("~")

#turn off pandas warnings
pd.options.mode.chained_assignment = None

'''
FirstWidget is the first tkinter parent
It contains multiple frames:
Frame 1: file selection, working directory selection, file out name
Frame 2: calibration standard data selection
Frame 3: climate data selection (only for EOS verison)
Frame 4: variable initiation (memory effect value and average by value)
'''
class FirstWidget:
    def __init__(self, parent):
        #initiate parent
        self.FirstTk = parent
        self.FirstTk.title("Choose files and set variables")

        #initiate frames
        #frames are packed from top to bottom (default packing order)

        #file selection
        self.Frame1 = Frame(self.FirstTk,bd=2,relief=GROOVE)
        self.Frame1.pack(pady=2,fill=X)

        #calibration data selection
        self.Frame2 = Frame(self.FirstTk,bd=2,relief=GROOVE)
        self.Frame2.pack(pady=2,fill=X)

        # #climate data selection
        # self.Frame3 = Frame(self.FirstTk,bd=2,relief=GROOVE)
        # self.Frame3.pack(pady=2,fill=X)

        #define parameters
        self.Frame4 = Frame(self.FirstTk,bd=2,relief=GROOVE)
        self.Frame4.pack(pady=2,fill=X)

        #plot what?
        self.Frame5 = Frame(self.FirstTk,bd=2,relief=GROOVE)
        self.Frame5.pack(pady=2,fill=X)

        #run/quit
        self.Frame6 = Frame(self.FirstTk)
        self.Frame6.pack(pady=2,fill=X)

        #####FRAME INFO#####

        ###Frame 1: set WD, input files, and output file names
        #left side:
        self.Frame1_left = Frame(self.Frame1)
        self.Frame1_left.pack(side='left', fill=Y)
        self.cd_button = Button(self.Frame1_left,text='Set Working Directory', command=self.OnCdButtonClick,height=1)
        self.cd_button.pack(anchor='w',pady=5)
        self.filechoose_frame1_left = Frame(self.Frame1_left)
        self.filechoose_frame1_left.pack(pady=20,fill=X)
        self.filechoose_button = Button(self.filechoose_frame1_left,text="Choose Input file", command=self.OnFilechooseClick,height=1)
        self.filechoose_button.pack(anchor='center')
        self.filechoose_delete_button = Button(self.filechoose_frame1_left, text='Remove entries', command=self.DeleteFiles)
        self.filechoose_delete_button.pack(anchor='center')
        self.fileout_label = Label(self.Frame1_left,text='File output name:',height=1)
        self.fileout_label.pack(side='bottom',pady=4,anchor='center')
        #right side
        self.Frame1_right = Frame(self.Frame1)
        self.Frame1_right.pack(side='left', fill=X)
        self.cd_label = Text(self.Frame1_right,height=1,bd=2,relief=RIDGE)
        self.cd_label.insert(END,'Please set the working directory')
        self.cd_label.pack(anchor='center')
        self.filechoose_frame1_right = Frame(self.Frame1_right)
        self.filechoose_frame1_right.pack(pady=5)
        self.filechoose_text1 = Text(self.filechoose_frame1_right,height=1,bd=2,relief=RIDGE)
        self.filechoose_text1.pack(anchor='center')
        self.filechoose_text1.insert(END, 'Input File 1')
        self.filechoose_text2 = Text(self.filechoose_frame1_right,height=1,bd=2,relief=RIDGE)
        self.filechoose_text2.pack(anchor='center')
        self.filechoose_text2.insert(END, 'Input File 2 (optional)')
        self.filechoose_text3 = Text(self.filechoose_frame1_right,height=1,bd=2,relief=RIDGE)
        self.filechoose_text3.pack(anchor='center')
        self.filechoose_text3.insert(END, 'Input File 3 (optional)')
        self.fileout_text = Text(self.Frame1_right,height=1,bd=2,relief=SUNKEN)
        self.fileout_text.pack(side='bottom',anchor='center')

        ###Frame 2: Standard data selection
        #left side
        self.Frame2_left = Frame(self.Frame2)
        self.Frame2_left.pack(side=LEFT,fill=Y)
        self.std_button = Button(self.Frame2_left,text='Choose Calibration Standards', command=self.OnStdButtonClick)
        self.std_button.pack(pady=4)
        self.standard_method_frame = Frame(self.Frame2_left)
        self.standard_method_frame.pack(pady=10)
        self.standard_method = StringVar()
        self.standard_method.trace('w', self.standard_method_varchange)
        self.std_useDefault_checkbutton = Radiobutton(self.standard_method_frame,text='Set/use default standards', variable=self.standard_method, value='def')
        self.std_useDefault_checkbutton.pack(anchor='w')
        self.or_label = Label(self.standard_method_frame,text='---OR---')
        self.or_label.pack(anchor='center')
        self.std_maninput_checkbutton = Radiobutton(self.standard_method_frame,text='Manual Standard Entry', variable=self.standard_method, value='man')
        self.std_maninput_checkbutton.pack(anchor='w')
        #right side
        self.Frame2_right = Frame(self.Frame2)
        self.Frame2_right.pack(side=LEFT,fill=BOTH)
        self.stddata_text = Text(self.Frame2_right,height=1,bd=2,relief=RIDGE)
        self.stddata_text.pack(anchor='n')
        self.useDefaultStd_text_Frame = Frame(self.Frame2_right)
        self.useDefaultStd_text_Frame.pack(pady = 25)
        self.scrollbar = Scrollbar(self.useDefaultStd_text_Frame, orient=HORIZONTAL)
        self.useDefaultStd_text = Text(self.useDefaultStd_text_Frame,height=1, wrap=NONE, xscrollcommand=self.scrollbar.set)
        self.useDefaultStd_text.pack(anchor='n',expand=YES)
        self.useDefaultStd_text.insert(END, 'Choose this box to set your default standards file to the file choosen above.')
        if os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
            f = open(home+'/.default_standard_data_PICARROALG.txt', 'r')
            path = f.read()
            self.useDefaultStd_text.delete(1.0, END)
            self.useDefaultStd_text.insert(END, 'Default file is set to '+path)
        self.scrollbar.pack(fill=X,anchor='n')
        self.scrollbar.config(command=self.useDefaultStd_text.xview)
        #cancel button
        self.Frame2_cancel = Frame(self.Frame2)
        self.Frame2_cancel.pack(side=RIGHT,fill=Y)
        self.cancel_button = Button(self.Frame2_cancel, text='Remove', command=self.RemoveStdData,height=1)
        self.cancel_button.pack(pady=4,padx=2)

        # ###Frame3: Climate data selection
        # #left side
        # self.climatedata = ''
        # self.Frame3_left = Frame(self.Frame3)
        # self.Frame3_left.pack(side=LEFT,fill=Y)
        # self.climatedata_button = Button(self.Frame3_left,text='Choose Climate Data', command=self.ChooseClimateData, height=1)
        # self.climatedata_button.pack(pady=4,fill=X)
        # self.useDefaultClimateData = IntVar()
        # self.useDefaultClimateData.trace('w', self.useDefaultClimateData_varchange)
        # self.useDefaultClimate_button = Checkbutton(self.Frame3_left, text='Use Default Climate Data', variable=self.useDefaultClimateData)
        # self.useDefaultClimate_button.pack(pady=4)
        # #right side
        # self.Frame3_right = Frame(self.Frame3)
        # self.Frame3_right.pack(side=LEFT,fill=BOTH)
        # self.climatedata_text = Text(self.Frame3_right,height=1,bd=2,relief=RIDGE)
        # self.climatedata_text.pack(anchor='n')

        ###Frame4: variable assignment (avgby and memeffval)
        #avg by
        self.Frame4_avgby = Frame(self.Frame4)
        self.Frame4_avgby.pack(anchor='w')
        self.avgby_label = Label(self.Frame4_avgby,text="Average By Value:")
        self.avgby_label.pack(side=LEFT,anchor='s')
        self.avgby_scale = Scale(self.Frame4_avgby,from_=2,to=60,orient=HORIZONTAL)
        self.avgby_scale.set(30)
        self.avgby_scale.pack(side=LEFT)
        self.use_avgby = IntVar()
        self.avgby_button = Checkbutton(self.Frame4_avgby, text='Disable average-by', variable=self.use_avgby)
        self.avgby_button.pack(side=LEFT, anchor='s')
        #memory effect 
        self.Frame4_memeff = Frame(self.Frame4)
        self.Frame4_memeff.pack(anchor='w') #packs this frame below the avg by frame
        self.memeffval_label = Label(self.Frame4_memeff,text="Memory Effect Value (permil dD):")
        self.memeffval_label.pack(side=LEFT,anchor='s')
        self.memeffval_scale = Scale(self.Frame4_memeff,from_=0,to=2,orient=HORIZONTAL,resolution=.05)
        self.memeffval_scale.set(.75)
        self.memeffval_scale.pack(side=LEFT)


        #plot what? (d18O or dD)
        self.Frame5_plotwhat = Frame(self.Frame5)
        self.Frame5_plotwhat.pack(anchor='w')
        self.plotwhat_label = Label(self.Frame5_plotwhat,text='Choose what to plot (choose one):')
        self.plotwhat_label.pack(anchor='center')
        self.plotwhat = StringVar()
        self.d18O_radiobutton = Radiobutton(self.Frame5_plotwhat, text='d18O', variable=self.plotwhat, value='raw_d18O')
        self.d18O_radiobutton.pack(anchor='w')
        self.dD_radiobutton = Radiobutton(self.Frame5_plotwhat, text='dD', variable=self.plotwhat, value='raw_dD')
        self.dD_radiobutton.pack(anchor='w')

        ###Frame5: run and quit buttons
        self.run_button = Button(self.Frame6,text="Run",command=self.run, width=5)
        self.run_button.pack(side=LEFT,padx=3)
        self.cancel_button = Button(self.Frame6,text="Quit",command=sys.exit, width=5)
        self.cancel_button.pack(side=LEFT,padx=3)

        #class variables:
        self.std_data_path = ''

    ###instructions for clicking the CD button
    def OnCdButtonClick(self):
        self.cd_label.grid_forget()
        ###user enters new directory
        self.dir_=filedialog.askdirectory()+'/'
        if self.dir_ != '':
            os.chdir(self.dir_)
            self.cd_label.delete(1.0, END)
            self.cd_label.insert(END, self.dir_)
        else:
            self.cd_label.delete(1.0, END)
            self.cd_label.insert(END, os.getcwd())

    ###instructions for clicking the filechoose button
    def OnFilechooseClick(self):
        dlg=filedialog.Open()
        #assign the path name of the selected file to temp
        temp = dlg.show()
        if temp is '':
            return None
        #assign the file name only to hold
        hold=re.split('/',temp)
        hold=hold[len(hold)-1]
        #insert text to the next avaiable input text box 
        if self.filechoose_text1.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text1.delete(1.0, END)
            self.filechoose_text1.insert(END, hold)
            self.fi1=temp
        elif self.filechoose_text2.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text2.delete(1.0, END)
            self.filechoose_text2.insert(END, hold)
            self.fi2=temp
        elif self.filechoose_text3.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text3.delete(1.0, END)
            self.filechoose_text3.insert(END, hold)
            self.fi3=temp
        else:
            print("All available file slots are taken. Remove input files using the 'Remove entries' button before selecting new files.")

    def DeleteFiles(self):
        if not self.filechoose_text3.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text3.delete(1.0, END)
            self.filechoose_text3.insert(END, 'Input file 3 (optional)')
        elif not self.filechoose_text2.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text2.delete(1.0, END)
            self.filechoose_text2.insert(END, 'Input file 2 (optional)')
        elif not self.filechoose_text1.get("1.0",'end-1c').split(" ")[0] == "Input":
            self.filechoose_text1.delete(1.0, END)
            self.filechoose_text1.insert(END, 'Input file 1')
        else:
            print('All files have been removed.')

    def OnStdButtonClick(self):
        dlg=filedialog.Open()
        self.std_data_path=dlg.show()
        ###get the name of the file only
        hold=re.split('/',self.std_data_path)
        hold=hold[len(hold)-1]
        ###create a label showing the new file name
        if not hold == '':
            self.stddata_text.delete(1.0, END)
            self.stddata_text.insert(END, hold)
            if self.standard_method.get() == 'def':
                self.useDefaultStd_text.delete(1.0, END)
                self.useDefaultStd_text.insert(END, 'Default file will be set to '+hold)

    def standard_method_varchange(self, *args):
        if self.standard_method.get() == 'man':
            self.std_data_path = ''
            self.useDefaultStd_text.delete(1.0, END)
            self.stddata_text.delete(1.0, END)
        if self.standard_method.get() == 'def':
            if self.std_data_path == '':
                if os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                    f = open(home+'/.default_standard_data_PICARROALG.txt', 'r')
                    path = f.read()
                    self.useDefaultStd_text.delete(1.0, END)
                    self.useDefaultStd_text.insert(END, 'Default file is set to '+path)
                if not os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                    self.useDefaultStd_text.delete(1.0, END)
                    self.useDefaultStd_text.insert(END, 'Default file will be set to the file choosen above')
            if not self.std_data_path == '':
                hold=re.split('/',self.std_data_path)
                hold=hold[len(hold)-1]
                self.useDefaultStd_text.delete(1.0, END)
                self.useDefaultStd_text.insert(END, 'Default file will be set to '+hold)
                    
    def RemoveStdData(self):
        self.std_data_path = ''
        self.stddata_text.delete(1.0, END)
        if self.standard_method.get() == 'def':
            self.useDefaultStd_text.delete(1.0, END)
            if os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                f = open(home+'/.default_standard_data_PICARROALG.txt', 'r')
                path = f.read()
                self.useDefaultStd_text.insert(END, 'Default file is set to '+path)
            if not os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                    self.useDefaultStd_text.insert(END, 'Default file will be set to the file choosen above')

    def ChooseClimateData(self):
        dlg=filedialog.Open()
        #assign the path name of the selected file to temp
        temp = dlg.show()
        if temp is '':
            return None
        #assign the file name only to hold
        hold=re.split('/',temp)
        hold=hold[len(hold)-1]
        self.climatedata_text.delete(1.0, END)
        self.climatedata_text.insert(END, hold)
        self.climatedata=temp
        f = open(home+'/.default_climate_data_PICARROALG.txt', 'w')
        f.write(temp)


    def useDefaultClimateData_varchange(self, *args):
        if self.useDefaultClimateData.get():
            if os.path.isfile(home+'/.default_climate_data_PICARROALG.txt'):
                f = open(home+'/.default_climate_data_PICARROALG.txt', 'r')
                self.climatedata = f.read()
                self.climatedata_text.delete(1.0, END)
                self.climatedata_text.insert(END, 'Default file is set to '+self.climatedata)
            else:
                self.climatedata_text.delete(1.0, END)
                self.climatedata_text.insert(END, 'No default climate data is available. Choose data manually first.')
                self.climatedata = ''
        else: 
            self.climatedata_text.delete(1.0, END)
            self.climatedata_text.insert(END, '')
            self.climatedata = ''

    ##instructions on clicking the run button:
    def run(self): 
        '''
        variables without self. infront of them are local to the run method
        this is ok as long as said variables don't need to be referenced elsewhere 
        or don't need to be pulled from another method inside the FirstWidget class
        '''

        #set the name for the output file
        self.fo=self.fileout_text.get(1.0,END)[:-1]

        #set standard methods and check for errors using self std_data_path and the std_method chosen
        #check over three conditions (standard_method.get())
        #1)
        if self.standard_method.get() == '':
            if not self.std_data_path == '':
                self.calibrate_bool = True
            if self.std_data_path == '':
                self.calibrate_bool = False
        #2)
        if self.standard_method.get() == 'man':
            if self.std_data_path == '':
                self.calibrate_bool = True
            if not self.std_data_path == '':
                print('Cannot have manual standard entry and standard data chosen at the same time! \n \
                Please chose either manual standard entry or select a standard data file.')
                #ends the method, returning user to FirstWidget
                return
        #3)
        if self.standard_method.get() == 'def':
            self.calibrate_bool = True
            if not self.std_data_path == '':
                #write the path to the standard data as the first line in the text file
                f = open(home+'/.default_standard_data_PICARROALG.txt', 'w')
                f.write(self.std_data_path)
                f.close()
            if self.std_data_path == '':
                #if the default text file does not exist yet print an error and exit the program
                if not os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                    print('ERROR: Default file ~/.default_standard_data_PICARROALG.txt does not exist \n \
                    Please choose a standard data file first, then use the "Set default standards" option.')
                    #ends the method, sending the user back to the FirstWidget
                    return
                #if the text files exists, read the text from the file as the path to the standard data
                if os.path.isfile(home+'/.default_standard_data_PICARROALG.txt'):
                    f = open(home+'/.default_standard_data_PICARROALG.txt', 'r')
                    #first chance for the program to set the std_data variable here
                    self.std_data_path = f.read()

        #set how many files will be used
        if not self.filechoose_text3.get("1.0",'end-1c').split(" ")[0] == "Input":
            use=2
        elif not self.filechoose_text2.get("1.0",'end-1c').split(" ")[0] == "Input":
            use=1
        else:
            use=0

        #determine whether or not data will be averaged
        avgby_bool = self.use_avgby.get()
        if avgby_bool == 1:
            avgby_bool = False
        else:
            avgby_bool = True 

        #set file input name and cancel run if there is no input file selected
        filename=self.filechoose_text1.get("1.0",'end-1c')
        if filename.split(" ")[0] == 'Input':
            print('Please select at least one file')
            return

        if self.plotwhat.get() == '':
            print('Please choose what to plot!')
            return

        #set fi to the input files that will be read
        if use==0:
            fi=[self.fi1]
        if use==1:
            fi=[self.fi1, self.fi2]
        if use==2:
            fi=[self.fi1, self.fi2, self.fi3]

        #set avgby
        avgby=self.avgby_scale.get()
        #set memeffval
        memeffval=self.memeffval_scale.get()
            
        #Run the manual standard entry widget if need be
        if self.standard_method.get() == 'man':
            root2 = Tk()
            ManStdEntry = ManualStdEntryWidget(root2)
            root2.mainloop()
        elif (self.standard_method.get() == '') and (self.std_data_path == ''):
            print('\nWARNING: No standard data selected or no standard entry method choosen! The time-series \
        graph will be displayed but data synthesis is not possible with calibration standards. \n')
        else:
            #standard date is read here for most cases
            self.std_data = pd.read_csv(self.std_data_path)

        #set data_inc_avg to the return value of the MainAlgorithm
        global data_inc_avg
        data_inc_avg = MainAlgorithm(avgby = avgby, memeffval = memeffval, fi = fi, fo = self.fo, dir_ = self.dir_, avgby_bool = avgby_bool)
        
        #create the initial figure with no highlights or thickens
        #this sets the initial value for f, col, num, t, and numticks
        fig = CreateFigure(TITLE=self.fo, COLOR='', NUMBER=0, NUMTICKS=16)

        # #creates the label widget
        #these two widgets are canceled by buttons inside the widgets
        self.root3 = Tk()
        self.label = LabelWidget(self.root3)
        #creates the figure widget inside the mainloop of the label widget
        self.root4 = Tk()

        #initializes a FigureWidget instance with the global variable f, created by CreateFigure() method above
        self.figure = FigureWidget(self.root4, fig)

        self.root4.mainloop() #inner loop = figure
        self.root3.mainloop() #outter loop = labels

#manual standard entry class
class ManualStdEntryWidget:
    def __init__(self, parent):
        self.ManInput = parent
        self.ManInput.title("Set standard data manually")

        self.frame1 = Frame(self.ManInput,bd=2,relief=GROOVE)
        self.frame1.grid(row=0,column=0,padx=3)
        self.frame2 = Frame(self.ManInput,bd=2,relief=GROOVE)
        self.frame2.grid(row=0,column=1,padx=3)
        self.frame3 = Frame(self.ManInput,bd=2,relief=GROOVE)
        self.frame3.grid(row=0,column=2,padx=3)
        self.button_frame = Frame(self.ManInput)
        self.button_frame.grid(row=1,column=1,)

        self.std_label = Label(self.frame1, text='Enter standard names below separated by commas')
        self.std_label.grid(row=0,column=0)
        self.v1=StringVar()
        self.std_entry = Entry(self.frame1,borderwidth=2,relief=SUNKEN,textvariable=self.v1)
        self.std_entry.grid(row=1,column=0)

        self.d18O_label = Label(self.frame2, text='Enter true delta-18O values below separated by commas')
        self.d18O_label.grid(row=0,column=0)
        self.v2=StringVar()
        self.d18O_entry = Entry(self.frame2,borderwidth=2,relief=SUNKEN,textvariable=self.v2)
        self.d18O_entry.grid(row=1,column=0)

        self.dD_label = Label(self.frame3, text='Enter true delta-deuterium values below separated by commas')
        self.dD_label.grid(row=0,column=0)
        self.v3=StringVar()
        self.dD_entry = Entry(self.frame3,borderwidth=2,relief=SUNKEN,textvariable=self.v3)
        self.dD_entry.grid(row=1,column=0)

        self.std_enter_button = Button(self.button_frame,text='Enter',command=lambda: self.EnterStds(Standards=self.v1,true_d18O=self.v2,true_dD=self.v3))
        self.std_enter_button.grid(row=0,column=0,columnspan=1,padx=3)
        self.cancel_button = Button(self.button_frame,text='Cancel',command=self.ManInput.destroy)
        self.cancel_button.grid(row=0,column=1,columnspan=1,padx=3)

    def EnterStds(self, Standards, true_d18O, true_dD):
        global std
        self.standards=pd.Series(re.split(",",Standards.get()), dtype='object')
        self.standards=list(map(str.strip, self.standards))
        self.d18O=pd.Series(re.split(",",true_d18O.get()), dtype='float64')
        self.dD=pd.Series(re.split(",",true_dD.get()), dtype='float64')
        std = pd.DataFrame({'standard' : self.standards, 'real.d18O' : self.d18O, 'real.D_H' : self.dD})
        self.ManInput.destroy()

class FigureWidget:
    def __init__(self, parent, fig):
        self.FigureRoot = parent
        self.FigureRoot.title("Time Series Plot")
        self.canvas = FigureCanvasTkAgg(fig, master=self.FigureRoot)
        self.canvas = self.canvas.get_tk_widget()
        self.canvas.pack(fill=BOTH,expand=YES)

    def UpdatePlot(self, fig):
        self.canvas.destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.FigureRoot)
        self.canvas = self.canvas.get_tk_widget()
        self.canvas.pack(fill=BOTH,expand=YES)

    def destroy(self):
        self.FigureRoot.destroy()
    
class LabelWidget:
    global num, col, numticks, t
    def __init__(self, parent):
        self.Label_Root = parent
        self.Label_Root.title("Enter event descriptors")
        self.scrollbar = Scrollbar(parent)
        self.scrollbar.grid(row=0, column=0)
        self.textbox_list = [] #holds identifiers for individual text boxes
        self.answer_list = pd.DataFrame({'color_group' : pd.Series(color_list), 'event' : pd.Series(repeat('event', len(color_list)))})

        self.Top_Label_Root = Frame(self.Label_Root,borderwidth=2,relief=GROOVE)
        self.Top_Label_Root.grid(row=0,column=0,padx=3,pady=3)
        self.Bottom_Label_Root = Frame(self.Label_Root,borderwidth=2,relief=GROOVE)
        self.Bottom_Label_Root.grid(row=1,column=0,padx=3,pady=3)

        #labels for the colors going across the top row
        self.redbutton = Button(self.Top_Label_Root,text='red',command=lambda: self.TriggerPlot(COL='red',NUM=num,T=first.fo,NT=numticks),width=4)
        self.redbutton.grid(row=0, column=1, sticky=S, pady=2)
        self.bluebutton = Button(self.Top_Label_Root,text='blue',command=lambda: self.TriggerPlot(COL='blue',NUM=num,T=first.fo,NT=numticks),width=5)
        self.bluebutton.grid(row=0, column=2, sticky=S, pady=2)
        self.greenbutton = Button(self.Top_Label_Root,text='green',command=lambda: self.TriggerPlot(COL='green',NUM=num,T=first.fo,NT=numticks),width=6)
        self.greenbutton.grid(row=0, column=3, sticky=S, pady=2)
        self.brownbutton = Button(self.Top_Label_Root,text='brown',command=lambda: self.TriggerPlot(COL='brown',NUM=num,T=first.fo,NT=numticks),width=6)
        self.brownbutton.grid(row=0, column=4, sticky=S, pady=2)

        #print the numbers vertically in column 0
        num_colorgroups = int(ceil(len(color_list)/4.0))+1
        button_list = []
        for I in range(1,num_colorgroups):
            Button(self.Top_Label_Root,text=I,command=lambda I=I: self.TriggerPlot(COL=col,NUM=I,T=first.fo,NT=numticks),width=2).grid(row=I, column=0, sticky=E, padx=2)

        #print text boxes for entry in the matrix
        for i in range(0,len(color_list)):
            self.textbox_list.append(Text(self.Top_Label_Root,height=1,width=25,borderwidth=2,relief=SUNKEN))
            self.textbox_list[i].bind("<Tab>", self.focus_next_window)
            self.textbox_list[i].grid(row=int(floor((i/4)+1)),column=(int(i%4)+1))

        #create enter button for the Label_Root
        self.highliteall_button = Button(self.Bottom_Label_Root, text='Highlight all', command=lambda: self.TriggerPlot(HL='all',T=first.fo,NT=numticks), width=18)
        self.highliteall_button.grid(row=0, column=0,sticky=E,padx=2)
        self.removehighlites_button = Button(self.Bottom_Label_Root, text='Remove all highlights', command=lambda: self.TriggerPlot(T=first.fo,NT=numticks), width=18)
        self.removehighlites_button.grid(row=0,column=1,sticky=W,padx=2)

        #create increase and decrease grid lines buttons
        self.increasegridlines_button = Button(self.Bottom_Label_Root, text='Increase gridlines', command=lambda: self.TriggerPlot(HL=hl,NUM=num,COL=col,T=first.fo,NT=numticks+5), width=18)
        self.increasegridlines_button.grid(row=1, column=0,sticky=E,padx=2)
        self.decreasegridlines_button = Button(self.Bottom_Label_Root, text='Decrease gridlines', command=lambda: self.TriggerPlot(HL=hl,NUM=num,COL=col,T=first.fo,NT=numticks-5), width=18)
        self.decreasegridlines_button.grid(row=1, column=1,sticky=W,padx=2)

        #create reference water label and entry widget for drift correction
        self.ref_label = Label(self.Bottom_Label_Root,text='Reference standard for drift correction:')
        self.ref_label.grid(row=2,column=0,sticky='NESW')
        self.ref_entry = Text(self.Bottom_Label_Root,height=1,width=25,borderwidth=2,relief=SUNKEN)
        self.ref_entry.grid(row=2,column=1,sticky='NESW')
        
        self.enter_button = Button(self.Bottom_Label_Root, text='OK', command=self.Frame2_LabelEntry, width=5)
        self.enter_button.grid(row=3, column=0,sticky=E,padx=2,pady=10)  
        self.quit_button = Button(self.Bottom_Label_Root, text='Quit', command=self.quit, width=5)
        self.quit_button.grid(row=3,column=1,sticky=W,padx=2,pady=10)
    
    def TriggerPlot(self,HL=None,NUM=0,COL='',NT=16,T=''):
        #creates a new f
        fig = CreateFigure(HIGHLIGHT=HL,NUMBER=NUM,COLOR=COL,NUMTICKS=NT,TITLE=T)
        #updates the figure widget with the new f
        first.figure.UpdatePlot(fig)
    
    def focus_next_window(self, event):
        event.widget.tk_focusNext().focus()
        return('break')

    def quit(self):
        #clear the working directory
        first.cd_label.delete(1.0, END)
        first.cd_label.insert(END, 'Set Working Directory')
        #clear the input files
        if not first.filechoose_text3.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text3.delete(1.0, END)
            first.filechoose_text3.insert(END, 'Input file 3 (optional)')
        if not first.filechoose_text2.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text2.delete(1.0, END)
            first.filechoose_text2.insert(END, 'Input file 2 (optional)')
        if not first.filechoose_text1.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text1.delete(1.0, END)
            first.filechoose_text1.insert(END, 'Input file 1')
        #clear the file output name
        first.fileout_text.delete(1.0, END)
        first.fo = ''
        first.root3.destroy()
        first.root4.destroy()
    
    def Frame2_LabelEntry(self):
        global data_inc_avg
        #foor loop through the cells to find the event values
        for i in range(len(self.answer_list['event'])):
            textbox_list_entry = self.textbox_list[i].get("1.0",'end-1c')
            #if cells values are empty, use the pervious cell's value
            if textbox_list_entry == '':
                #if the first cell is empty, make the first cell a blank string
                if i == 0:
                    self.answer_list['event'][i] = ''
                else:
                    textbox_list_entry = self.answer_list['event'][i-1]
            self.answer_list['event'][i] = textbox_list_entry
        #if all the cells are empty, tell the program to just save the plot figure and end the program
        self.PlotFigOnly_bool = False
        if not any(self.answer_list['event'] != ''):
            self.PlotFigOnly_bool = True
        self.drift_entry = self.ref_entry.get("1.0", 'end-1c')
        self.drift_bool = True
        if self.drift_entry == '':
            self.drift_bool = False

        #prints out the time series plot
        plt.clf()
        f = CreateFigure(HIGHLIGHT='all', TITLE = first.fo)
        self.c=FigureCanvas(f)
        self.c.print_figure(first.dir_+first.fo+'.pdf')
        print('Time series plot printed to working directory.')

        #enter in the event descriptors to the data_inc_avg dataframe
        EventDescriptor()

        #drift correct if instructed to
        #need to change this to the inner thirds method
        if self.drift_bool:
            #initiate the boolean as false to start
            self.UseBasicDriftCorrection = False
            #try to do IT drift correction
            IT_DriftCorrection()
            #if it doesn't work, the boolean variable initiated above will trigger basic drift correction
            if self.UseBasicDriftCorrection:
                DriftCorrection()
        else:
            print('Data will not be drift corrected!')
            #these columns will be removed before the data is printed out
            #they need to be in for the time being so that the data works with the following sections
            data_inc_avg['drift_corrected_d18O'] = data_inc_avg['d18O']
            data_inc_avg['drift_corrected_dD'] = data_inc_avg['dD']


        #do calibrations if instructed to, otherwise, just leave data as is
        if first.calibrate_bool:
            Calibration()
        else:
            print('Date will not be calibrated!')
            #these will be deleted later, but need to be included for now
            data_inc_avg['calculated_d18O'] = data_inc_avg['d18O']
            data_inc_avg['calculated_dD'] = data_inc_avg['dD']

        # data_inc_avg.to_csv(home+"/Desktop/data_inc_avg.csv")
        # sys.exit()
        #print out the lmwl
        LMWL()

        #correct the data based on drift_bool and calibrate_bool
        if not first.label.drift_bool:
            data_inc_avg = data_inc_avg.drop('drift_corrected_d18O', 1)
            data_inc_avg = data_inc_avg.drop('drift_corrected_dD', 1)
        if not first.calibrate_bool:
            data_inc_avg = data_inc_avg.drop('calculated_d18O', 1)
            data_inc_avg = data_inc_avg.drop('calculated_dD', 1)

        #add in climate data if need be
        if first.climatedata is not '':
            ClimateData()

        #write out the data
        data_inc_avg.to_csv(first.dir_+first.fo+'_analysis.csv')
        print('Analysis data written to working directory.')

        ###compute averages and standard deviations for all non-mem events
        grouped = data_inc_avg.groupby('event', as_index = False)
        iso_aggs = grouped['calculated_d18O', 'calculated_dD'].agg([np.mean, np.std]) 
        iso_aggs.to_csv(first.dir_+first.fo+'_details.csv')
        print('Details data written to working directory.')

        #clear the working directory
        first.cd_label.delete(1.0, END)
        first.cd_label.insert(END, 'Set Working Directory')
        #clear the input files
        if not first.filechoose_text3.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text3.delete(1.0, END)
            first.filechoose_text3.insert(END, 'Input file 3 (optional)')
        if not first.filechoose_text2.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text2.delete(1.0, END)
            first.filechoose_text2.insert(END, 'Input file 2 (optional)')
        if not first.filechoose_text1.get("1.0",'end-1c').split(" ")[0] == "Input":
            first.filechoose_text1.delete(1.0, END)
            first.filechoose_text1.insert(END, 'Input file 1')
        #clear the file output name
        first.fileout_text.delete(1.0, END)
        first.fo = ''


        first.figure.destroy()
        self.Label_Root.destroy()


###GLOBAL METHODS
#following method modifies global variable f to represent the figure you intend to show
#use global variable f when creating an instance of FigureWidget to tell the instance what image to display
#also can use the class method FigureWidget.UpdatePlot(f) to change the figure on the widget
def CreateFigure(HIGHLIGHT=None, COLOR='', NUMBER=0, TITLE='', NUMTICKS=16):
    #set global variables for continuity between graphs
    global a,f,num,col,color_list,numticks,t,hl
    num=NUMBER
    col=COLOR
    numticks=NUMTICKS
    hl=HIGHLIGHT
    t=TITLE
    f = plt.Figure()
    a = f.add_subplot(111)
    color_list = []

    #print out details about which color/number combo or highlight description is active
    if (COLOR=='') & (NUMBER==0):
        print('highlight: ' + str(HIGHLIGHT))
    else:
        if COLOR=='':
            c='none selected'
        else:
            c=COLOR
        if NUMBER==0:
            n='none selected'
        else:
            n=NUMBER
        print('color: ' + str(c) + ', number: ' + str(n))

    max_group = max(data_inc_avg['group'])
    #print each color section 1 by 1
    for i in range(1,(max_group+2)):
        if i == max_group + 1:
            hold = data_inc_avg[data_inc_avg['group'] == 0]
            a.plot(hold['index'], hold[first.plotwhat.get()], color = 'black', marker = '.', linestyle = 'None', alpha = .1) #this line plots memory effects
            continue
        if i != max_group + 1:
            hold = data_inc_avg[data_inc_avg['group'] == i]

            #plot all highlighted and thickened
            if HIGHLIGHT is 'all':
                if i % 4 == 1:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'red')
                    s='red '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 2:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'blue')
                    s='blue '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 3:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'green')
                    s='green '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 0:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'saddlebrown')
                    s='brown '+str(int(round(i/4)))
                    color_list.append(s)
                    continue

            #plot nothing thickened or highlighted
            if (COLOR is '') and (NUMBER is 0):
                if i % 4 == 1:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'red')
                    s='red '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 2:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'blue')
                    s='blue '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 3:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'green')
                    s='green '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 0:
                    a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'saddlebrown')
                    s='brown '+str(int(round(i/4)))
                    color_list.append(s)
                    continue

            #plot a certain color highlighted and thickened
            if (col is not '') and (num is 0):
                if i % 4 == 1:
                    if col is 'red':
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'red')
                    else:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'red')
                    s='red '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 2:
                    if col is 'blue':
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'blue')
                    else:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'blue')
                    s='blue '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 3:
                    if col is 'green':
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'green')
                    else:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'green')
                    s='green '+str(int((round(i/4))+1))
                    color_list.append(s)
                    continue
                if i % 4 == 0:
                    if col is 'brown':
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'saddlebrown')
                    else:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'saddlebrown')
                    s='brown '+str(int(round(i/4)))
                    color_list.append(s)
                    continue

            #plot a certain number highlighted and thickened
            if (col is '') and (num is not 0):
                if i in [x+((num-1)*4) for x in range(1,5)]:
                    if i % 4 == 1:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'red')
                        s='red '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 2:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'blue')
                        s='blue '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 3:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'green')
                        s='green '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 0:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'saddlebrown')
                        s='brown '+str(int(round(i/4)))
                        color_list.append(s)
                        continue
                else:
                    if i % 4 == 1:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'red')
                        s='red '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 2:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'blue')
                        s='blue '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 3:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'green')
                        s='green '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 0:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'saddlebrown')
                        s='brown '+str(int(round(i/4)))
                        color_list.append(s)
                        continue
            #plot a certain color/number combination as highlighted and thickened
            if (col is not '') and (num is not 0):
                if i in [x+((num-1)*4) for x in range(1,5)]:
                    if i % 4 == 1:
                        if col is 'red':
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'red')
                        else:
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'red')
                        s='red '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 2:
                        if col is 'blue':
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'blue')
                        else:
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'blue')
                        s='blue '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 3:
                        if col is 'green':
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'green')
                        else:
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'green')
                        s='green '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 0:
                        if col is 'brown':
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 4.5, color = 'saddlebrown')
                        else:
                            a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'saddlebrown')
                        s='brown '+str(int(round(i/4)))
                        color_list.append(s)
                        continue
                else:
                    if i % 4 == 1:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'red')
                        s='red '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 2:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'blue')
                        s='blue '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 3:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'green')
                        s='green '+str(int((round(i/4))+1))
                        color_list.append(s)
                        continue
                    if i % 4 == 0:
                        a.plot(hold['index'], hold[first.plotwhat.get()], linewidth = 3.0, alpha = .25, color = 'saddlebrown')
                        s='brown '+str(int(round(i/4)))
                        color_list.append(s)
                        continue

    by = len(data_inc_avg.index)/numticks 
    at = [x*by for x in range(1, (numticks+1))]
    l = []
    for i in at:
        l.append(data_inc_avg['time'][i])

    a.set_xticks(at, minor=False)
    a.set_xticklabels(l, fontdict=None, minor=False, rotation=70)
    a.set_xlabel('time')
    a.set_ylabel(first.plotwhat.get())
    a.set_title(t)
    a.grid(which='major')
    return f


'''
All classes and global methods are now defined.

The following code runs the actual program 
by calling the widgets in order and assigning 
values to globa variables that need to be used 
by future widgets and code.
'''



def MainAlgorithm(avgby, memeffval, fi, fo, dir_, avgby_bool):
    print('Averaging by: %d' % avgby)
    print('Memory Effect Value: %f' % float(memeffval))
    print('Running Program on: %s' % fi)
    print('Files named '+fo+' will be written out to '+dir_)
    os.chdir(home)

    print("Preparing data...")

    ###read in the data
    frames = []
    for f in fi:
        filetype= f.split(".")
        if filetype[len(filetype)-1] == 'dat':
            temp = pd.read_csv(f, delim_whitespace=True)
        elif filetype[len(filetype)-1] == 'csv':
            temp = pd.read_csv(f) 
        else:
            print("ERROR: input file must either be a .csv or a .dat file")
            sys.exit()
        frames.append(temp)
    master_data = pd.concat(frames, ignore_index=True)
    #removes empty rows at the end of the file
    master_data = master_data[np.isfinite(master_data['EPOCH_TIME'])]

    #add time to the master data using the epoch time column
    master_data['seconds'] = map(lambda e: int(time.strftime('%S', time.localtime(e))), master_data['EPOCH_TIME'])
    master_data['hourminutes'] = map(lambda e: time.strftime('%H:%M', time.localtime(e)), master_data['EPOCH_TIME'])
    master_data['date'] = map(lambda e: time.strftime('%d/%m/%Y', time.localtime(e)), master_data['EPOCH_TIME'])

    nrows = len(master_data.axes[0])
    ncols = len(master_data.axes[1])
    delta_18O = master_data['Delta_18_16']

    data_inc_avg = pd.DataFrame({'date' : pd.Series([NaN]), 'time' : pd.Series([NaN]), 'raw_d18O' : pd.Series([NaN]), 'raw_dD' : pd.Series([NaN]), 'H2O' : pd.Series([NaN])})
    #if avgby_bool is true, avg by avgby
    if avgby_bool:
        for i in range(0, nrows-1):
            if i == 0:
                start = i
                continue
            #find what second we are at in the current iteration of the for loop    
            s = master_data['seconds'][i]
            #if we are at a seconds value equal to 0 mod avgby, average the preceding rows
            if s % avgby != 0:
                continue
            elif s % avgby == 0:
                #case for if we are at the last row
                if i == nrows - 1:
                    temp = pd.DataFrame({'date' : pd.Series([master_data['date'][i]]),
                        'time' : pd.Series([master_data['hourminutes'][i]]),
                        'raw_d18O' : pd.Series([master_data['Delta_18_16'][start:i].mean()]), 
                        'raw_dD' : pd.Series([master_data['Delta_D_H'][start:i].mean()]),
                        'H2O' : pd.Series([master_data['H2O'][start:i].mean()])})
                    data_inc_avg = data_inc_avg.append(temp, ignore_index = True)
                    break

                s_next = master_data['seconds'][i+1]
                #case for if the next line is also 0 mod avgby
                if s_next % avgby == 0:
                    continue
                #case if the next line is not 0 mod avgby
                if s_next % avgby != 0:
                    temp = pd.DataFrame({'date' : pd.Series([master_data['date'][start]]),
                        'time' : pd.Series([master_data['hourminutes'][start]]),
                        'raw_d18O' : pd.Series([master_data['Delta_18_16'][start:i].mean()]), 
                        'raw_dD' : pd.Series([master_data['Delta_D_H'][start:i].mean()]),
                        'H2O' : pd.Series([master_data['H2O'][start:i].mean()])})
                    data_inc_avg = data_inc_avg.append(temp, ignore_index = True)
                    start = i + 1
                    continue
        #delete the first row of NaN's
        data_inc_avg = data_inc_avg.ix[1:]
        data_inc_avg['index'] = list(range(1, len(data_inc_avg.axes[0])+1))
    #if label.useavgby is false, then just use the master data as data_inc_avg
    else:
        data_inc_avg = data_inc_avg.append(pd.DataFrame({'date' : pd.Series([master_data['date']]), \
            'time' : pd.Series(master_data['hourminutes']), \
            'raw_d18O' : pd.Series(master_data['Delta_18_16']), \
            'raw_dD' : pd.Series(master_data['Delta_D_H']), \
            'H2O' : pd.Series(master_data['H2O'])}), ignore_index = True)
        data_inc_avg = data_inc_avg.ix[1:]
        data_inc_avg['index'] = list(range(1, len(data_inc_avg.axes[0])+1))

    #create event column and enter in 'mem' for wherever the change is > memeffval
    mem_series = ['']
    for i in data_inc_avg['index']:
        if i == 1:
            continue
        if abs(data_inc_avg['raw_dD'][i] - data_inc_avg['raw_dD'][i-1]) > memeffval:
            mem_series.append('mem')
        else:
            mem_series.append('')
    data_inc_avg['event'] = mem_series

    #fill in mem gaps 
    for i in data_inc_avg['index']:
        if (i == 1) or (i == 2):
            continue
        if (data_inc_avg['event'][i] == 'mem') and \
        (data_inc_avg['event'][i-2] == 'mem') and \
        (abs(data_inc_avg['raw_dD'][i] - data_inc_avg['raw_dD'][i-1]) > (memeffval - .15)):
            data_inc_avg['event'][i-1] = 'mem'
            continue

    #give group numbers based on separation by mems
    group = []
    for i in data_inc_avg['index']:
        if i == 1:
            tracker = 1
        if i == len(data_inc_avg['event']):
            if data_inc_avg['event'][i] == 'mem':
                group.append(0)
                break
            if data_inc_avg['event'][i] != 'mem':
                group.append(tracker)
                break
        if data_inc_avg['event'][i] != 'mem':
            group.append(tracker)
            continue
        if data_inc_avg['event'][i] == 'mem':
            group.append(0)
            if data_inc_avg['event'][i+1] == 'mem':
                continue
            if data_inc_avg['event'][i+1] != 'mem':
                tracker = tracker + 1
                continue
    data_inc_avg['group'] = group

    #add col_group for later matching
    col_group = []
    for i in data_inc_avg['group']:
        if i!=0:
            if i%4==1:
                col_group.append('red '+str(int((round(i/4))+1)))
            if i%4==2:
                col_group.append('blue '+str(int((round(i/4))+1)))
            if i%4==3:
                col_group.append('green '+str(int((round(i/4))+1)))
            if i%4==0:
                col_group.append('brown '+str(int(round(i/4))))
        if i==0:
            col_group.append('mem')
    data_inc_avg['col_group'] = col_group

    return data_inc_avg


def EventDescriptor():
    print('Filling things in...')
    event = []
    for i in data_inc_avg['col_group']:
        if i == 'mem':
            event.append('mem')
            continue
        if i != 'mem':
            if vers == 2:
                event.append(str([first.label.answer_list['event'][j] for j in range(len(first.label.answer_list['color_group'])) if first.label.answer_list['color_group'][j] == i])[3:-2])
            elif vers == 3:
                event.append(str([first.label.answer_list['event'][j] for j in range(len(first.label.answer_list['color_group'])) if first.label.answer_list['color_group'][j] == i])[2:-2])
    data_inc_avg['event'] = event

def DriftCorrection():
    print("Use basic drift correction...")
    ###Drift Correction
    ###keep only reference water data
    drift_data = data_inc_avg[data_inc_avg['event'] == first.label.drift_entry]
    ###get rid of the first and last 10 minutes
    drift_data = drift_data[21:-20]
    ###keep the remaining first and last 10% of the data
    keep = int(max(drift_data['index'])*.1)
    keep_d18O = [np.mean(drift_data['raw_d18O'][0:keep]), np.mean(drift_data['raw_d18O'][-keep:max(drift_data['index'])])]
    keep_dD = [np.mean(drift_data['raw_dD'][0:keep]), np.mean(drift_data['raw_dD'][-keep:max(drift_data['index'])])]
    drift_d18O = (keep_d18O[1]-keep_d18O[0])/max(drift_data['index'])
    drift_dD = (keep_dD[1]-keep_dD[0])/max(drift_data['index'])
    print('Drift correction d18O slopes by segment: \n'+ str(drift_d18O))
    DCd18O = []
    DCdD = []
    for n, d in enumerate(data_inc_avg['raw_d18O']):
        DCd18O.append(d - n*drift_d18O)
    for n, d in enumerate(data_inc_avg['raw_dD']):
        DCdD.append(d - n*drift_dD)
    data_inc_avg['drift_corrected_d18O'] = DCd18O
    data_inc_avg['drift_corrected_dD'] = DCdD

def IT_DriftCorrection():
    data_inc_avg['drift_corrected_d18O'] = data_inc_avg['raw_d18O']
    data_inc_avg['drift_corrected_dD'] = data_inc_avg['raw_dD']
    opt1 = None
    drift_entry = first.label.drift_entry

    #find TW group numbers
    TWgroups = [data_inc_avg['group'].iloc[n] for n, i in enumerate(data_inc_avg['event']) if i == drift_entry]
    #keep TW group numbers that are TW segments longer than 21 minutes
    useable_twg = []
    for twg in pd.unique(TWgroups):
        dt = [datetime.datetime.strptime(data_inc_avg['date'].iloc[n] + data_inc_avg['time'].iloc[n], '%d/%m/%Y%H:%M') for n, i in enumerate(data_inc_avg['group']) if i == twg]
        dt_diff = dt[len(dt)-1]-dt[0]
        if (dt_diff.seconds / 60) >= 21:
            useable_twg.append(twg)

    #find non-TW and non-mem group numbers
    CALgroups = [data_inc_avg['group'].iloc[n] for n, i in enumerate(data_inc_avg['event']) if (i != 'mem') & (i != drift_entry) & (i != opt1)]

    #find which TW segments will be needed
    need_twg = []
    can_continue = True
    for cg in pd.unique(CALgroups):
        #if the one of the Calibration groups is group number 1, there is no TW data before it, so IT method cannot be used
        if cg == 1:
            first.label.UseBasicDriftCorrection = True
            print("IT drift correction failed")
            return
        #if cg != 1, find the above and below TW groups
        if cg > 1:
            twg_below = [t for t in useable_twg if t < cg]
            twg_above = [t for t in useable_twg if t > cg]
        #if the twg_below or _above are empty, we have to use basic drift correction
        if (len(twg_below) == 0) | (len(twg_above) == 0):
            first.label.UseBasicDriftCorrection = True
            print("IT drift correction failed")
            return
        #otherwise, append twg_below and _above to need_twg
        need_twg.append(twg_below[-1])
        need_twg.append(twg_above[0])

    need_twg = pd.unique(need_twg)
    print("Using inner-thirds drift correction method...")

    #find the color groups that were used - useful for identification of TW segments used in drift correction for future reference
    need_twcolgroups = [data_inc_avg[data_inc_avg.group == x]['col_group'].iloc[0] for x in need_twg]

    #average the inner thirds of the needed TW groups
    TWsegment_center = []
    IT_d18O = []
    IT_d2H = []
    for ntwg in need_twg:
        #extract rows 
        temp = data_inc_avg[data_inc_avg.group == ntwg]

        #calculate inner third row sizes
        nrows = temp.shape[0]
        inner_third_rows = range(nrows/3, 2*nrows/3)

        #find the center row
        TWsegment_center.append((nrows/2) + temp.index.values[0])

        #slice out the inner third rows and average
        ITtemp = temp.iloc[inner_third_rows]
        IT_d18O.append(np.mean(ITtemp['raw_d18O']))
        IT_d2H.append(np.mean(ITtemp['raw_dD']))

    #compute the slopes
    slope_d18O = []
    slope_d2H = []
    for i in range(len(IT_d18O)):
        if i == (len(IT_d18O)-1):
            continue
        else:
            run = TWsegment_center[i+1] - TWsegment_center[i]
            rise_d18O = IT_d18O[i+1] - IT_d18O[i]
            rise_d2H = IT_d2H[i+1] - IT_d2H[i]
            slope_d18O.append(rise_d18O / run)
            slope_d2H.append(rise_d2H / run)

    #apply the slopes to the rows that are contained in the drift correction regions
    for n, i in enumerate(TWsegment_center):
        if n == (len(TWsegment_center)-1):
            continue
        else:
            temp_d18O = data_inc_avg.iloc[TWsegment_center[n]:TWsegment_center[n+1]]['raw_d18O']
            temp_dD = data_inc_avg.iloc[TWsegment_center[n]:TWsegment_center[n+1]]['raw_dD']
            hold_d18O = []
            hold_dD = []
            for n2, t in enumerate(temp_d18O):
                hold_d18O.append(t - (n2*slope_d18O[n]))
            for n2, t in enumerate(temp_dD):
                hold_dD.append(t - (n2*slope_d2H[n]))
            data_inc_avg['drift_corrected_d18O'][TWsegment_center[n]:TWsegment_center[n+1]] = hold_d18O
            data_inc_avg['drift_corrected_dD'][TWsegment_center[n]:TWsegment_center[n+1]] = hold_dD
    print('Used TW color groups for drift correction: \n'+ str(need_twcolgroups))
    print('Drift correction d18O slopes by segment: \n'+ str(slope_d18O))


def Calibration():
    print("Calibrating...")
    ###Calibration
    # data_inc_avg = pd.read_csv(home+'/Desktop/test.csv')
    #std_data = pd.read_csv(home+'/Documents/EOS/data/NTU.reference.standards.csv')
    grouped = data_inc_avg.groupby('event', as_index = False)
    # event_means = data_inc_avg.groupby('event', as_index = False).mean()
    event_means = grouped.agg(np.mean)
    # event_stddev = data_inc_avg.groupby('event', as_index = False).std()
    measured_stds = event_means[event_means['event'].isin(first.std_data['standard'])]
    real_d18O = []
    real_dD = []
    for i in measured_stds['event']:
        real_d18O.append(float(first.std_data['real.d18O'][first.std_data['standard']==i]))
        real_dD.append(float(first.std_data['real.D_H'][first.std_data['standard']==i]))
    ###compute linear regression stats
    d18O_calibration_plot = plt.plot(measured_stds['drift_corrected_d18O'], real_d18O)
    dD_calibration_plot = plt.plot(measured_stds['drift_corrected_dD'], real_dD)
    slope_d18O, intercept_d18O, r_value_d18O, p_value_d18O, std_err_d18O = stats.linregress(measured_stds['drift_corrected_d18O'], real_d18O)
    slope_dD, intercept_dD, r_value_dD, p_value_dD, std_err_dD = stats.linregress(measured_stds['drift_corrected_dD'], real_dD)
    ###predict calculated values using the linear regression info
    data_inc_avg['calculated_d18O'] = [slope_d18O*x + intercept_d18O for x in data_inc_avg['drift_corrected_d18O']]
    data_inc_avg['calculated_dD'] = [slope_dD*x + intercept_dD for x in data_inc_avg['drift_corrected_dD']]

#function that checks a vector for 'rain' entries
def checkForRain(inputString):
    return bool(re.search('rain', inputString))

def LMWL():
    plt.clf()
    d = data_inc_avg    
    #check for rain
    keep = list(map(checkForRain, d['event']))
    keep = [i for i in range(len(keep)) if keep[i]]
    #if there is rain, print the LMWL
    if len(keep) > 0:
        keep = d.iloc[keep,:]
        x = keep['calculated_d18O']
        y = keep['calculated_dD']
        m, b = np.polyfit(x, y, 1)
        plt.plot(x, y, '.')
        y_hat = [m*X + b for X in x]
        plt.plot(x, y_hat, '-')
        plt.ylabel('dD')
        plt.xlabel('d18O')
        plt.title('LWML')
        plt.savefig(first.dir_+first.fo+'_LMWL.pdf')
        print('')
        print('LWML printed')
        print('')



#add in the climate data:
#requires the compiled climate data in csv format. 
#compiled data is created using new.climate.data.aggregation.R in ~/R/bin
#read in climate data

# def ClimateData():
#     print("Adding Climate Data")
#     global data_inc_avg
#     climate_data = pd.read_csv(first.climatedata)
#     #add in datetime info
#     climate_data_datetime = map(lambda y,mo,d,h,mi: datetime.datetime(year=y,month=mo,day=d,hour=h,minute=mi), climate_data['year'],climate_data['month'],climate_data['day'],climate_data['hour'],climate_data['minute'])
#     climate_data['datetime'] = climate_data_datetime
#     #add datetime info to data_inc_Avg
#     data_inc_avg['datetime'] = map(lambda d, t: datetime.datetime.strptime(d + t, '%d/%m/%Y%H:%M'), data_inc_avg['date'],data_inc_avg['time'])

#     #remove entries from climate_data we don't need
#     first_dt = data_inc_avg['datetime'].irow(0)
#     last_dt = data_inc_avg['datetime'].irow(-1)
#     hold = [n for n, x in enumerate(climate_data['datetime']) if first_dt <= x <= last_dt]
#     hold = climate_data.iloc[hold]

#     #match the indices of the two dataframes
#     hold = hold.set_index('datetime')
#     hold.index.name = None
#     data_inc_avg = data_inc_avg.set_index('datetime')
#     data_inc_avg.index.name = None

#     #transfer columns to data_inc_avg
#     data_inc_avg['mm.rain'] = hold['mm.rain']
#     data_inc_avg['temp'] = hold['temp']
#     data_inc_avg['RH'] = hold['RH']

#     #remove the cluttered index names from data_inc_avg
#     data_inc_avg = data_inc_avg.reset_index()



#Run the first widget
root = Tk()
first = FirstWidget(root)
root.mainloop()

sys.exit()

