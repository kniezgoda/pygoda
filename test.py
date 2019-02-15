
from Tkinter import *

class First:
	def __init__(self):
		self.root = Tk()
		self.frame = Frame(self.root)
		self.frame.pack()
		self.Run_Button = Button(self.frame, text = "RUN", command = self.onRunClicked)
		self.Run_Button.pack()
		self.fun = dict()
	
	def onRunClicked(self):
		self.window = Second()
		self.window.root.mainloop()

class Second:
	def __init__(self):
		self.root = Toplevel()

		self.frame = Frame(self.root,bd=2,relief=GROOVE)
		self.frame.pack()

		self.savefigbool = IntVar() 

		self.checkbutton = Checkbutton(self.frame, variable=self.savefigbool)
		self.checkbutton.pack()
		
		self.Run_Button = Button(self.root, text = "RUN", command = self.onRunClicked, height = 1)
		self.Run_Button.pack()

	def onRunClicked(self):
		print(self.savefigbool.get())


def call_second_window():
    second_window = Second()


app = First()
app.fun = call_second_window
app.root.mainloop()