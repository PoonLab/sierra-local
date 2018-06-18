import tkinter as tk
from tkinter import filedialog
import os
import subprocess

class App():
   def __init__(self, parent):
      self.parent = parent

      self.parent.winfo_toplevel().title("sierralocal")

      # INPUT FRAME

      self.inputframe = tk.Frame(self.parent)
      self.inputframe.pack()

      #select file frame

      self.selectframe = tk.Frame(self.inputframe)
      self.selectframe.pack(side = tk.TOP)

      self.selectInputfileButton = tk.Button(self.selectframe, text="Select...", fg="black", command=self.getInputFile)
      self.selectInputfileButton.pack( side = tk.RIGHT)

      self.inputFileLabelVar = tk.StringVar()
      self.inputFileLabel = tk.Label(self.selectframe, textvariable=self.inputFileLabelVar, font = "Verdana 10 bold")
      self.inputFileLabelVar.set("Input file")
      self.inputFileLabel.pack(side=tk.LEFT)

      self.inputFileNameVar = tk.StringVar()
      self.inputFileName = tk.Label(self.selectframe, textvariable=self.inputFileNameVar)
      self.inputFileName.pack(side=tk.LEFT)

      # paste text frame

      self.pasteframe = tk.Frame(self.inputframe)
      self.pasteframe.pack(side = tk.RIGHT)

      self.inputFileLabelVar = tk.StringVar()
      self.inputFileLabel = tk.Label(self.pasteframe, textvariable=self.inputFileLabelVar, font = "Verdana 10 bold")
      self.inputFileLabelVar.set(" or Paste Sequences")
      self.inputFileLabel.pack()

      self.pastetext = tk.Text(self.pasteframe)
      self.pastetext.pack()

      # OUTPUT FRAME

      self.outputframe = tk.Frame(self.parent)
      self.outputframe.pack()

      self.selectOutputfileButton = tk.Button(self.outputframe, text="Select...", fg="black", command=self.setOutputFile)
      self.selectOutputfileButton.pack( side = tk.RIGHT)

      self.outputFileLabelVar = tk.StringVar()
      self.outputFileLabel = tk.Label(self.outputframe, textvariable=self.outputFileLabelVar, font = "Verdana 10 bold")
      self.outputFileLabelVar.set("Output file")
      self.outputFileLabel.pack(side=tk.LEFT)

      self.outputfileNameVar = tk.StringVar()
      self.outputFileName = tk.Label(self.outputframe, textvariable=self.outputfileNameVar)
      self.outputFileName.pack(side=tk.LEFT)

      # OPTIONS FRAME

      self.optionframe = tk.Frame(self.parent)
      self.optionframe.pack()

      self.skipalign = tk.IntVar()
      self.skipalignbutton = tk.Checkbutton(self.optionframe, text="Skip NucAmino alignment", variable=self.skipalign)
      self.skipalignbutton.pack(side=tk.LEFT)

      self.cleanup = tk.IntVar()
      self.cleanupbutton = tk.Checkbutton(self.optionframe, text="Preserve NucAmino file.", variable=self.cleanup)
      self.cleanupbutton.select()
      self.cleanupbutton.pack(side=tk.RIGHT)

      # CONSOLE FRAME

      consoleFrame = tk.Frame(self.parent)
      consoleFrame.pack()

      self.text = tk.Text(consoleFrame)
      self.text.pack()

      # RUN FRAME

      self.runframe = tk.Frame(self.parent)
      self.runframe.pack(side=tk.BOTTOM)

      self.runbutton = tk.Button(self.runframe, text="RUN", fg="black", command=lambda: self.runSierra(self.inputFileNameVar.get(), self.outputfileNameVar.get(), self.skipalign.get(), self.cleanup.get()))
      self.runbutton.pack( side = tk.LEFT)

      self.exitbutton = tk.Button(self.runframe, text="EXIT", fg="black", command=self.parent.destroy)
      self.exitbutton.pack(side=tk.RIGHT)

   def getInputFile(self):
      inputfile = filedialog.askopenfilename(initialdir = "~/git/sierra-local",title = "Select file",filetypes = (("FASTA files","*.fa"),("all files","*.*")))
      self.inputFileNameVar.set(inputfile)

   def setOutputFile(self):
      outputfile = filedialog.asksaveasfilename(initialdir = "~/git/sierra-local",title = "Select file",filetypes = (("JSON file","*.json"),("all files","*.*")))
      self.outputfileNameVar.set(outputfile)

   def runSierra(self, inputfile, outputfile, skipalign, cleanup):
      pastedtext = self.pastetext.get("1.0",'end-1c').strip()
      if len(pastedtext) < 10: #no pasted sequence
         self.text.delete(1.0,tk.END)
         skip = lambda x:' -skipalign' if x==1 else ''
         clean = lambda x:' -cleanup' if x==0 else ''
         p = subprocess.Popen('python3 sierralocal.py {} -o {}{}{}'.format(inputfile, outputfile, skip(skipalign), clean(cleanup)),stdout=subprocess.PIPE, shell=True)
         output, err = p.communicate()
         self.text.insert(tk.END, output)
      else: #pasted sequence
         with open("temp.fa", "w") as temp:
            temp.write(pastedtext)
         self.text.delete(1.0,tk.END)
         skip = lambda x:' -skipalign' if x==1 else ''
         clean = lambda x:' -cleanup' if x==0 else ''
         p = subprocess.Popen('python3 sierralocal.py {} -o {}{}{}'.format("temp.fa", outputfile, skip(skipalign), clean(cleanup)),stdout=subprocess.PIPE, shell=True)
         output, err = p.communicate()
         self.text.insert(tk.END, output)
         os.remove("temp.fa")



def main():
   root = tk.Tk()
   app = App(root)
   root.mainloop()

if __name__ == "__main__":
    main()
    