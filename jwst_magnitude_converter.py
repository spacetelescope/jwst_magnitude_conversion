#! /usr/bin/env python
#
"""
This python program is for reading in a set of observed magnitude/colour 
values and converting them to expected JWST magnitudes.  The conversion 
is based on simulated stellar colours from magnitudes.py.  One can also 
transform from one JWST instrument magnitude to another.

To use the code non-interactively one has to specify parameters in a configobj
file.  See the example file hst_to_nis.cfg for an example.  The magnitude 
names and other parameters are specified in parameter=value pairs in the file 
whereupon the code can produce an output file without bringing up the interface. 
Note that if the non-interactive calculation fails there is no indication of 
why it fails.

"""

from __future__ import print_function
from __future__ import division
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
import numpy 
import numpy.polynomial.legendre as legendre
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import astropy.io.fits as fits
from configobj import ConfigObj
import sys, os

if sys.version_info[0] == 2:
    import Tkinter as Tk
    import ttk 
    from ScrolledText import ScrolledText
    import tkFileDialog
    import tkMessageBox
    import tkColorChooser
elif sys.version_info[0] == 3:
    import tkinter as Tk
    import tkinter.ttk as ttk
    from tkinter.scrolledtext import ScrolledText
    import tkinter.filedialog as tkFileDialog
    from tkinter.messagebox import Message as tkMessageBox
    from tkinter.colorchooser import Chooser as tkColorChooser

path = os.path.dirname(__file__)
kuruczfilternames=["Sloan u  ","Sloan g  ","Sloan r  ","Sloan i  ",
"Sloan z  ","Bessel U  ","Bessel B  ","Bessel V  ",
"Bessel/Cousins R  ","Bessel/Cousins I  ","Johnson U  ",
"Johnson B  ","Johnson V  ","Johnson R  ",
"Johnson I  ","Johnson J  ","Johnson H  ",
"Johnson K  ","Johnson L  ","Johnson M  ",
"Johnson N  ","IRAS [12] ","IRAS [25] ","IRAS [60] ",
"IRAS [100] ","2MASS J ","2MASS H ","2MASS Ks ","IRAC [3.6] ","IRAC [4.4] ",
"IRAC [5.7] ","IRAC [8.0] ","MIPS [24] ","MIPS [70] ","MIPS [160] ",
"DENIS i ","DENIS J ","DENIS Ks ","WISE W1 ","WISE W2 ","WISE W3 ","WISE W4 ",
"UKIDSS Z ","UKIDSS Y ","UKIDSS J ","UKIDSS H ",
"UKIDSS K ","Pan-STARRS1 g ","Pan-STARRS1 r ",
"Pan-STARRS1 i ","Pan-STARRS1 z ","Pan-STARRS1 y ",
"Pan-STARRS1 w ","HST ACS F220W ","HST ACS F250W ","HST ACS F330W ",
"HST ACS F435W ","HST ACS F475W ","HST ACS F555W ","HST ACS F606W ",
"HST ACS F625W ","HST ACS F775W ","HST ACS F814W ","HST WFC3 F218W ",
"HST WFC3 F225W ","HST WFC3 F275W ","HST WFC3 F336W ","HST WFC3 F390W ",
"HST WFC3 F438W ","HST WFC3 F475W ","HST WFC3 F555W ","HST WFC3 F606W ",
"HST WFC3 F625W ","HST WFC3 F775W ","HST WFC3 F814W ","HST WFC3 F105W ",
"HST WFC3 F110W ","HST WFC3 F125W ","HST WFC3 F140W ","HST WFC3 F160W ",
"GAIA g","GAIA gbp","GAIA grp",
"NIRCam F070W ","NIRCam F090W ","NIRCam F115W ","NIRCam F140M ",
"NIRCam F150W ","NIRCam F150W2 ","NIRCam F162M ","NIRCam F164N ",
"NIRCam F182M ","NIRCam F187M ","NIRCam F200W ","NIRCam F210M ",
"NIRCam F212N ","NIRCam F250M ","NIRCam F277W ","NIRCam F300M ",
"NIRCam F322W2 ","NIRCam F323N ","NIRCam F335M ","NIRCam F356W ",
"NIRCam F360M ","NIRCam F405N ","NIRCam F410M ","NIRCam F430M ",
"NIRCam F444W ","NIRCam F460M ","NIRCam F466N ","NIRCam F470N ",
"NIRCam F480M ","MIRI F560W ","MIRI F770W ","MIRI F1000W ",
"MIRI F1280W ","MIRI F1130W ","MIRI F1500W ","MIRI F1800W ",
"MIRI F2100W ","MIRI F2550W ",'NIRSpec F110W','NIRSpec F140X',
'NIRSpec Clear','NIRSpec F070LP','NIRSec F100LP','NIRSpec F170LP',
'NIRSPec F290LP','Guider 1','Guider 2']

phoenixfilternames=['Sloan u ','Sloan g ','Sloan r ','Sloan i ','Sloan z ',
'Bessel U ','Bessel B ','Bessel V ','Bessel/Cousins R ','Bessel/Cousins I ',
'Johnson U ','Johnson B ','Johnson V ','Johnson R ','Johnson I ','Johnson J ',
'Johnson H ','Johnson K ','Johnson L ','2MASS J ','2MASS H ','2MASS Ks ',
'IRAC [3.6] ','IRAC [4.4] ','DENIS i ','DENIS J ','DENIS Ks ','WISE W1 ',
'WISE W2 ','UKIDSS Z ','UKIDSS Y ','UKIDSS J ','UKIDSS H ','UKIDSS K ',
'Pan-STARRS1 g ','Pan-STARRS1 r ','Pan-STARRS1 i ','Pan-STARRS1 z ',
'Pan-STARRS1 y ','Pan-STARRS1 w ','HST ACS F220W ','HST ACS F250W ',
"GAIA g","GAIA gbp","GAIA grp",
'HST ACS F330W ','HST ACS F435W ','HST ACS F475W ','HST ACS F555W ',
'HST ACS F606W ','HST ACS F625W ','HST ACS F775W ','HST ACS F814W ',
'HST WFC3 F218W ','HST WFC3 F225W ','HST WFC3 F275W ','HST WFC3 F336W ',
'HST WFC3 F390W ','HST WFC3 F438W ','HST WFC3 F475W ','HST WFC3 F555W ',
'HST WFC3 F606W ','HST WFC3 F625W ','HST WFC3 F775W ','HST WFC3 F814W ',
'HST WFC3 F105W ','HST WFC3 F110W ','HST WFC3 F125W ','HST WFC3 F140W ',
'HST WFC3 F160W '
"NIRCam F070W ","NIRCam F090W ","NIRCam F115W ","NIRCam F140M ",
"NIRCam F150W ","NIRCam F150W2 ","NIRCam F162M ","NIRCam F164N ",
"NIRCam F182M ","NIRCam F187M ","NIRCam F200W ","NIRCam F210M ",
"NIRCam F212N ","NIRCam F250M ","NIRCam F277W ","NIRCam F300M ",
"NIRCam F322W2 ","NIRCam F323N ","NIRCam F335M ","NIRCam F356W ",
"NIRCam F360M ","NIRCam F405N ","NIRCam F410M ","NIRCam F430M ",
"NIRCam F444W ","NIRCam F460M ","NIRCam F466N ","NIRCam F470N ",
"NIRCam F480M ","MIRI F560W ","MIRI F770W ","MIRI F1000W ",
"MIRI F1280W ","MIRI F1130W ","MIRI F1500W ","MIRI F1800W ",
"MIRI F2100W ","MIRI F2550W ",'NIRSpec F110W','NIRSpec F140X',
'NIRSpec Clear','NIRSpec F070LP','NIRSec F100LP','NIRSpec F170LP',
'NIRSPec F290LP','Guider 1','Guider 2']

jwstnames=["NIRISS F090W ","NIRISS F115W ","NIRISS F140M ","NIRISS F150W ",
"NIRISS F158M ","NIRISS F200W ","NIRISS F277W ","NIRISS F356W ",
"NIRISS F380M ","NIRISS F430M ","NIRISS F444W ","NIRISS F480M ",
'Guider 1','Guider 2',
"NIRCam F070W ","NIRCam F090W ","NIRCam F115W ","NIRCam F140M ",
"NIRCam F150W ","NIRCam F150W2 ","NIRCam F162M ","NIRCam F164N ",
"NIRCam F182M ","NIRCam F187M ","NIRCam F200W ","NIRCam F210M ",
"NIRCam F212N ","NIRCam F250M ","NIRCam F277W ","NIRCam F300M ",
"NIRCam F322W2 ","NIRCam F323N ","NIRCam F335M ","NIRCam F356W ",
"NIRCam F360M ","NIRCam F405N ","NIRCam F410M ","NIRCam F430M ",
"NIRCam F444W ","NIRCam F460M ","NIRCam F466N ","NIRCam F470N ",
"NIRCam F480M ","MIRI F560W ","MIRI F770W ","MIRI F1000W ",
"MIRI F1280W ","MIRI F1130W ","MIRI F1500W ","MIRI F1800W ",
"MIRI F2100W ","MIRI F2550W ",
'NIRSpec F110W', 'NIRSpec F140X', 'NIRSpec Clear',
'NIRSpec F070LP','NIRSec F100LP','NIRSpec F170LP', 'NIRSPec F290LP']

class magConGUI():
  def __init__(self,parent=None,**args):
    self.haveData=False
    self.plotOption=0
    self.yinvert=False
    self.haveTransformation=False
    self.plotFrames=[None,None,None]
    self.canvasList=[None,None,None]
    self.subplotList=[None,None,None]
    self.mplfigList=[None,None,None]
    self.posLabelText=[None,None,None]
    self.posLabel=[None,None,None]
    self.xdata=[None,None,None]
    self.ydata=[None,None,None]
    self.xlabel=[None,None,None]
    self.ylabel=[None,None,None]
    self.rangeEntry=[None,None,None,None]
# for each plot, xmin, xmax, ymin, ymax values saved here....
    self.plotLimits=numpy.zeros((4,3),dtype=numpy.float32)
# for each plot, xmin, xmax, ymin, ymax data values saved here...
    self.dataLimits=numpy.zeros((4,3),dtype=numpy.float32)
# self.fitResults will hold the fit parameters for from 1 to 59 JWST filters for the
# different instruments (12 NIRISS, 2 Guider, 29 NIRCam, 9 MIRI, 7 NIRSpec)
    self.fitResults=[None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
# self.fitRange holds the colour range of validity of the fitting
    self.fitRange=numpy.zeros((2),dtype=numpy.float32)

  def runGUI(self,root):
    if root is None:
      pass
    else:
#      root.__init__(self,parent,args)
      self.makeWidgets(root)
#      self.pack(expand=Tk.YES,fill=Tk.BOTH)
    
  def roundFloat(self,value,flag):
    if value == 0.:
      return value
    else:
      if value < 0.:
        sign=-1.
        value= -value
      else:
        sign=1.
      power=math.log10(value)
      if power < 0.:
        exp=int(power-1.)
      else:
        exp=0.0
      shift=10.**exp
      x=value/shift
      if flag == 0 and sign > 0.:
        x=math.floor(x)
      if flag == 0 and sign < 0.:
        x=math.ceil(x)
      if flag == 1 and sign > 0.:
        x=math.ceil(x)
      if flag == 1 and sign < 0.:
        x=math.floor(x)
    return x*shift*sign

  def sepLine(self,parent,w1,h1,pad):
    lincanvas=Tk.Canvas(parent,height=h1,width=w1)
    try:
      lincanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=Tk.YES)
    except:
      pass
    lincanvas.create_line(pad,h1/2,w1-pad,h1/2)
    return lincanvas

  def makeWidgets(self,root):
    self.root=root
    self.areaLabel=Tk.Label(root,text="Message Area")
    self.areaLabel.pack(side=Tk.TOP)
    self.messageText=ScrolledText(root,height=6,width=50,bd=1,relief=Tk.RIDGE,wrap=Tk.NONE)
    self.messageText.config(font=('courier',12,'bold'))
    self.messageText.pack()
    sep1=self.sepLine(root,400,16,10)
    selectFrame=Tk.Frame(root)
    selectFrame.pack(side=Tk.TOP)
    label1=Tk.Label(selectFrame,text='Input Magnitude 1')
    label1.grid(row=0,column=0)
    label2=Tk.Label(selectFrame,text='Input Magnitude 2')
    label2.grid(row=0,column=1)
    self.mag1box=ttk.Combobox(selectFrame,width=20)
    self.mag1box.grid(row=1,column=0)
    self.mag1box['values']=kuruczfilternames
    self.mag1box.current(11)
    self.mag2box=ttk.Combobox(selectFrame,width=20)
    self.mag2box.grid(row=1,column=1)
    self.mag2box['values']=kuruczfilternames
    self.mag2box.current(12)
    label1=Tk.Label(selectFrame,text='Column in File')
    label1.grid(row=2,column=0)
    label2=Tk.Label(selectFrame,text='Column in File')
    label2.grid(row=2,column=1)
    self.mag1column=Tk.Entry(selectFrame,width=5)
    self.mag1column.grid(row=3,column=0)
    self.mag1column.insert(0,'1')
    self.mag2column=Tk.Entry(selectFrame,width=5) 
    self.mag2column.grid(row=3,column=1)
    self.mag2column.insert(0,'2')
    sel1=Tk.Frame(selectFrame)
    sel1.grid(row=4,column=0)
    self.mag1option=Tk.IntVar()
    opt1=Tk.Radiobutton(sel1,text='Magnitude',variable=self.mag1option,value=0)
    opt1.pack(side=Tk.LEFT)
    opt2=Tk.Radiobutton(sel1,text='Colour',variable=self.mag1option,value=1)
    opt2.pack(side=Tk.LEFT)
    sel2=Tk.Frame(selectFrame)
    sel2.grid(row=4,column=1)
    self.mag2option=Tk.IntVar()
    opt1=Tk.Radiobutton(sel2,text='Magnitude',variable=self.mag2option,value=0)
    opt1.pack(side=Tk.LEFT)
    opt2=Tk.Radiobutton(sel2,text='Colour',variable=self.mag2option,value=1)
    opt2.pack(side=Tk.LEFT)
    self.yAxisValue=Tk.IntVar()
    opt1=Tk.Radiobutton(selectFrame,text='Y axis?',variable=self.yAxisValue,value=0)
    opt1.grid(row=5,column=0)
    opt2=Tk.Radiobutton(selectFrame,text='Y axis?',variable=self.yAxisValue,value=1)
    opt2.grid(row=5,column=1)
    label3=Tk.Label(selectFrame,text='Optional RA Column:')
    label3.grid(row=6,column=0)
    label4=Tk.Label(selectFrame,text='Optional Dec Column:')
    label4.grid(row=6,column=1)
    self.racolumn=Tk.Entry(selectFrame,width=5)
    self.racolumn.grid(row=7,column=0)
    self.racolumn.insert(0,'0')
    self.deccolumn=Tk.Entry(selectFrame,width=5) 
    self.deccolumn.grid(row=7,column=1)
    self.deccolumn.insert(0,'0')
    Tk.Button(root,text="Read Data",command=self.readData).pack(side=Tk.TOP)
    sep1=self.sepLine(root,400,16,10)
    outFrame=Tk.Frame(root)
    outFrame.pack(side=Tk.TOP)
    label1=Tk.Label(outFrame,text="Output JWST Filter(s):")
    label1.pack(side=Tk.TOP)
    self.outOption=Tk.IntVar()
    sel3=Tk.Frame(outFrame)
    sel3.pack(side=Tk.TOP)
    opt1=Tk.Radiobutton(sel3,text='All Filters       ',variable=self.outOption,value=0)
    opt1.grid(row=0,column=0,sticky=Tk.W)
    opt2=Tk.Radiobutton(sel3,text='One Filter        ',variable=self.outOption,value=1)
    opt2.grid(row=1,column=0,sticky=Tk.W)
    opt3=Tk.Radiobutton(sel3,text='Two Filters       ',variable=self.outOption,value=2)
    opt3.grid(row=2,column=0,sticky=Tk.W)
    setFrame=Tk.Frame(root)
    setFrame.pack(side=Tk.TOP)
    self.setOption=Tk.IntVar()
    label1=Tk.Label(setFrame,text="Use magnitude set:")
    label1.pack(side=Tk.LEFT)
    mopt1=Tk.Radiobutton(setFrame,text='Kurucz',variable=self.setOption,value=0,command=self.setLabels)
    mopt1.pack(side=Tk.LEFT)
    mopt2=Tk.Radiobutton(setFrame,text='Phoenix',variable=self.setOption,value=1,command=self.setLabels)
    mopt2.pack(side=Tk.LEFT)
    mopt3=Tk.Radiobutton(setFrame,text='BOSZ',variable=self.setOption,value=3,command=self.setLabels)
    mopt3.pack(side=Tk.LEFT)
    mopt4=Tk.Radiobutton(setFrame,text='Blackbodies',variable=self.setOption,value=2,command=self.setLabels)
    mopt4.pack(side=Tk.LEFT)
    l1=Tk.Frame(root)
    l1.pack()
    label1=Tk.Label(l1,text="Order of fit: ")
    label1.pack(side=Tk.LEFT)
    self.fitOrder=Tk.Entry(l1,width=3)
    self.fitOrder.pack()
    self.fitOrder.insert(0,'4')
    Tk.Button(root,text="Transform",command=self.doTransformation).pack()
    Tk.Button(root,text="Write Out Transformation",command=self.writeTransformation).pack()
    self.outOption.set(2)
    self.mag3box=ttk.Combobox(sel3,width=15)
    self.mag3box.grid(row=1,column=1)
    self.mag3box['values']=jwstnames
    self.mag3box.current(0)
    self.mag4box=ttk.Combobox(sel3,width=15)
    self.mag4box.grid(row=2,column=1)
    self.mag4box['values']=jwstnames
    self.mag4box.current(1)
    sep1=self.sepLine(root,400,16,10)
    buttonFrame=Tk.Frame(root)
    buttonFrame.pack(side=Tk.TOP)
    b1=Tk.Button(buttonFrame,text="Plot Input Colour-Magnitude Diagram",command=self.plotInputData)
    b1.pack(side=Tk.TOP)
    b2=Tk.Button(buttonFrame,text="Plot Colour Transformation",command=self.plotTransform)
    b2.pack(side=Tk.TOP)
    b3=Tk.Button(buttonFrame,text="Plot Output Colour-Magnitude Diagram",command=self.plotOutputData)
    b3.pack(side=Tk.TOP)
    sep1=self.sepLine(root,400,16,10)
    Tk.Button(root,text="Write Magnitude Values",command=self.writeValues).pack(side=Tk.TOP)
    Tk.Button(root,text="Close Widget",command=self.onExit).pack(side=Tk.TOP)
# To use the full blackbody range change the comment on the first "status=" line
# over to the second line.
#
#    status=self.readModelValues('magslist_old_kurucz.new','magslist_phoenix_grid.new','magslist_blackbody_fullrange.new')
    status=self.readModelValues('magslist_old_kurucz.new','magslist_phoenix_grid.new','magslist_blackbody.new','magslist_bosz_normal.new')
    if not status:
      self.putMessage('Error reading in the model magnitude values.  No transformations can be done.\n',self.messageText)
      self.haveModelMags=False
    else:
      self.putMessage('Have read in the Kurucz, BOSZ, Phoenix, and blackbody \nmodel magnitude values.\n',self.messageText)
      self.haveModelMags=True

  def setLabels(self):
    setoption=self.setOption.get()
    if setoption == 1:
      self.mag1box['values']=phoenixfilternames
      self.mag2box['values']=phoenixfilternames
    else:
      self.mag1box['values']=kuruczfilternames
      self.mag2box['values']=kuruczfilternames
    if self.haveModelMags and self.haveData:
        self.doTransformation()
      
  def writeValues(self):
    if not self.haveData:
      return
    outfilename=tkFileDialog.asksaveasfilename(filetypes=[('data files','*.data'),('data files','*.dat'),('data files','*.txt'),('data files','*.tbl'),('All files','*')])
    if isinstance(outfilename,basestring):
      self.writeMags(outfilename)

  def writeMags(self,outfilename):
    try:
      value=self.jwstMags[0,0]
    except:
      self.putMessage('Error: no transformed values to write out',self.messageText)
      return
    outfile=open(outfilename,'w')
    s1='# %s | %s ' % (self.xlabel[0],self.ylabel[0])
    nout=[]
    for n in range(59):
      values=self.jwstMags[:,n]
      amin=numpy.min(values)
      amax=numpy.max(values)
      dela=amax-amin
      if dela > 0.:
        s1=s1+' | '+jwstnames[n]
        nout.append(n)
    s1=s1.rstrip(' ,')
    if not self.ravalues is None:
      s1=s1+' | RA | Dec '
    print(s1,file=outfile)
    npoints=len(self.xdata[0])
    for n in range(npoints):
      s1='%10.5f %10.5f ' % (self.xdata[0][n],self.ydata[0][n])
      for m in range(len(nout)):
        s1=s1+'%10.5f ' % (self.jwstMags[n,nout[m]])
      if not self.ravalues is None:
        s1=s1+'%s %s' % (self.ravalues[n],self.decvalues[n])
      print(s1,file=outfile)
    outfile.close()
      
  def readModelValues(self,filename1,filename2,filename3,filename4):
    try:
      self.kuruczMagValues,self.kuruczModelMagLabels,self.kuruczFilterPars=self.readMagslist(os.path.join(path,filename1),142)
      self.phoenixMagValues,self.phoenixModelMagLabels,self.phoenixFilterPars=self.readMagslist(os.path.join(path,filename2),121)
      self.blackbodyMagValues,self.blackbodyModelMagLabels,self.blackbodyFilterPars=self.readMagslist(os.path.join(path,filename3),142)
      self.boszMagValues,self.boszModelMagLabels,self.boszFilterPars=self.readMagslist(os.path.join(path,filename4),142)
      if self.boszMagValues is None or self.kuruczMagValues is None or self.phoenixMagValues is None or self.blackbodyMagValues is None:
        return False
      return True
    except:
      return False

  def gettext(self,filename,n):
    if n < 0:
      return None
    try:
      infile=open(filename,'r')
      lines=infile.readlines()
      infile.close()
      outvalues=[]
      for line in lines:
        if '#' in line[0:1] or '\\' in line[0:1] or '|' in line[0:1]:
          pass
        else:
          values=line.split()
          outvalues.append(values[n])
      return outvalues
    except:
      return None
  
  def readMagslist(self,filename,nfilters):
    try:
      infile=open(filename,'r')
      modelMagLabels=[]
      filterPars=numpy.zeros((nfilters,3),dtype=numpy.float32)
      line=infile.readline()
      line=infile.readline()
      line=infile.readline()
      line=infile.readline()
      line=infile.readline()
      for i in range(nfilters):
        line=infile.readline()
        values=line.split()
        m=len(values)
        str=''
        for n in range(2,m-3):
          str=str+values[n]+' '
        str=str.replace('filter','')
        str=str.replace('Filter','')
        modelMagLabels.append(str)
        filterPars[i,0]=float(values[-3])
        filterPars[i,1]=float(values[-2])
        filterPars[i,2]=float(values[-1])
      infile.close()
      modelMagValues=numpy.loadtxt(filename)
      return modelMagValues,modelMagLabels,filterPars
    except:
      return None,None,None
      
  def set1(self,event):
    self.currentPlot=1
    
  def set2(self,event):
    self.currentPlot=2
    
  def set3(self,event):
    self.currentPlot=3
    
  def plotInputData(self):
    self.plotOption=1
    self.makePlot(self.plotOption)

  def plotTransform(self):
    self.plotOption=2
    self.makePlot(self.plotOption)

  def doTransformation(self):
    if not self.haveModelMags:
      self.putMessage('Error: no transformation can be calculated.\n',self.messageText)
      return
    if not self.haveData:
      self.putMessage('Error: read in data first.\n',self.messageText)
      return
    self.fitResults=[None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
    filter1=self.mag1box.get()
    filter2=self.mag2box.get()
    filter3=self.mag3box.get()
    filter4=self.mag4box.get()
    magopt=self.outOption.get()
    setopt=self.setOption.get()
    mopt1,mopt2,mopt3,mopt4=self.matchFilter(setopt,filter1,filter2,filter3,filter4)
    if mopt1 < 0 or mopt2 < 0 or mopt3 < 0 or mopt4 < 0:
      self.putMessage('Error: could not match the specfied filters, so no transformation can be calculated.\n',self.messageText)
      return
    self.xlabel[1]=filter1+' - '+filter2
    yopt=self.yAxisValue.get()
    if yopt == 0:
      self.ylabel[1]=filter1+' - '+filter3
    else:
      self.ylabel[1]=filter2+' - '+filter3
    ndatapoints=len(self.xdata[0])
    self.jwstMags=numpy.zeros((ndatapoints,59),dtype=numpy.float32)
    self.yopt=yopt
    norder=int(self.fitOrder.get())
    self.fit1(magopt,setopt,yopt,norder,mopt1,mopt2,mopt3,mopt4,True)
    xmin,xmax,ymin,ymax=self.getRange(self.xdata[1],self.ydata[1])
    self.dataLimits[0,1]=xmin
    self.dataLimits[1,1]=xmax
    self.dataLimits[2,1]=ymin
    self.dataLimits[3,1]=ymax
    self.haveTransformation=True
    self.transformedJWSTMags=numpy.zeros((ndatapoints,50),dtype=numpy.float32)
    self.jwstInds=[mopt1,mopt2,mopt3,mopt4]
    if yopt == 1:
      self.xdata[2]=numpy.copy(self.xdata[0])
      self.ydata[2]=numpy.copy(self.jwstMags[:,mopt3])
    else:
      self.xdata[2]=self.jwstMags[:,mopt3]-self.jwstMags[:,mopt4]
      self.ydata[2]=self.jwstMags[:,mopt3]
    xmin,xmax,ymin,ymax=self.getRange(self.xdata[2],self.ydata[2])
    self.dataLimits[0,2]=xmin
    self.dataLimits[1,2]=xmax
    self.dataLimits[2,2]=ymin
    self.dataLimits[3,2]=ymax
    self.xlabel[2]=filter3+' - '+filter4
    self.ylabel[2]=filter3
    self.putMessage('The transformation has been calculated.\n',self.messageText)

  def matchFilter(self,setopt,filter1,filter2,filter3,filter4):
    mopt1=-10
    mopt2=-10
    mopt3=-10
    mopt4=-10
    if setopt == 0:
      for i in range(len(self.kuruczModelMagLabels)):
        if filter1 in self.kuruczModelMagLabels[i]:
          mopt1=i
        if filter2 in self.kuruczModelMagLabels[i]:
          mopt2=i
        if filter3 in self.kuruczModelMagLabels[i]:
          mopt3=i
        if filter4 in self.kuruczModelMagLabels[i]:
          mopt4=i
    if setopt == 1:
      for i in range(len(self.phoenixModelMagLabels)):
        if filter1 in self.phoenixModelMagLabels[i]:
          mopt1=i
        if filter2 in self.phoenixModelMagLabels[i]:
          mopt2=i
        if filter3 in self.phoenixModelMagLabels[i]:
          mopt3=i
        if filter4 in self.phoenixModelMagLabels[i]:
          mopt4=i
    if setopt == 2:
      for i in range(len(self.blackbodyModelMagLabels)):
        if filter1 in self.blackbodyModelMagLabels[i]:
          mopt1=i
        if filter2 in self.blackbodyModelMagLabels[i]:
          mopt2=i
        if filter3 in self.blackbodyModelMagLabels[i]:
          mopt3=i
        if filter4 in self.blackbodyModelMagLabels[i]:
          mopt4=i
    if setopt == 3:
      for i in range(len(self.boszModelMagLabels)):
        if filter1 in self.boszModelMagLabels[i]:
          mopt1=i
        if filter2 in self.boszModelMagLabels[i]:
          mopt2=i
        if filter3 in self.boszModelMagLabels[i]:
          mopt3=i
        if filter4 in self.boszModelMagLabels[i]:
          mopt4=i
    return mopt1,mopt2,mopt3,mopt4

  def doFit(self,x1,y1,y2,n,yopt,norder,interactive):
# case 1:  something like J - K vs J, do J -K to J - JWST filter
      if yopt == 0:
        self.modelCol1=x1-y1
        self.modelCol2=x1-y2
# case2: something like B - V vs V, do B - V to V - JWST filter
      else:
        self.modelCol1=x1-y1
        self.modelCol2=y1-y2
      try:
        if norder < 2:
          norder=4
          if interactive:
            self.putMessage('Error: bad fit order.  Putting the value to 4.',self.messageText)
            self.fitOrder.delete(0,Tk.END)
            self.fitOrder.insert(0,str(norder))
      except:
        norder=4
        if interactive:
          self.putMessage('Error: bad fit order.  Putting the value to 4.',self.messageText)
          self.fitOrder.delete(0,Tk.END)
          self.fitOrder.insert(0,str(norder))
      self.fitResults[n]=legendre.legfit(self.modelCol1,self.modelCol2,norder)
      fitvalues=legendre.legval(self.modelCol1,self.fitResults[n])
      rms=numpy.sqrt(numpy.sum((fitvalues-self.modelCol2)*(fitvalues-self.modelCol2))/len(fitvalues))
      mindev=numpy.min(numpy.abs(fitvalues-self.modelCol2))
      maxdev=numpy.max(numpy.abs(fitvalues-self.modelCol2))
      if interactive:
        self.putMessage('Transformation RMS value for filter\n %s: %.4f\nRange: %.4f to %.4f\n\n' % (jwstnames[n],rms,mindev,maxdev),self.messageText)
      else:
        print('Transformation RMS value for filter\n %s: %.4f\nRange: %.4f to %.4f\n\n' % (jwstnames[n],rms,mindev,maxdev))
      self.fitRange[0]=numpy.min(self.modelCol1)    
      self.fitRange[1]=numpy.max(self.modelCol1)
      xmin=numpy.min(self.modelCol1)-0.5
      xmax=numpy.max(self.modelCol1)+0.5
      newcol=legendre.legval(self.xdata[0],self.fitResults[n])
      ymin=legendre.legval(xmin,self.fitResults[n])
      ymax=legendre.legval(xmax,self.fitResults[n])
      newcol[self.xdata[0] < xmin]=ymin
      newcol[self.xdata[0] > xmax]=ymax
      if yopt == 0:
        self.jwstMags[:,n]=self.ydata[0]-newcol
      else:
        self.jwstMags[:,n]=self.ydata[0]+self.xdata[0]-newcol
    
  def plotOutputData(self):
    self.plotOption=3
    self.makePlot(self.plotOption)

  def putMessage(self,s1,messagearea):
    messagearea.insert(Tk.END,s1)
    messagearea.see(Tk.END)
    
  def makePlot(self,plotOption):
    if not self.haveData:
      self.putMessage('Error: no data values have been read in yet.\n',self.messageText)
      return
    if plotOption == 3 and not self.haveTransformation:
      self.putMessage('Error: no transformation has been produced in yet.\n',self.messageText)
      return
    ind=plotOption-1
    if self.plotFrames[ind] is None:
      self.plotFrames[ind]=Tk.Toplevel()
      if ind == 0:
        self.plotFrames[ind].bind('<Enter>',self.set1)
      if ind == 1:
        self.plotFrames[ind].bind('<Enter>',self.set2)
      if ind == 2:
        self.plotFrames[ind].bind('<Enter>',self.set3)
      self.posLabelText[ind]=Tk.StringVar()
      self.posLabel[ind]=Tk.Label(self.plotFrames[ind],textvariable=self.posLabelText[ind])
      self.posLabel[ind].pack()
      self.posLabelText[ind].set("Position: ")
      self.mplfigList[ind]=Figure(figsize=(5,5),dpi=100)
      self.subplotList[ind]=self.mplfigList[ind].add_subplot(1,1,1)
      self.canvasList[ind]=FigureCanvasTkAgg(self.mplfigList[ind],master=self.plotFrames[ind])
      self.canvasList[ind].show()
      self.canvasList[ind].get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=Tk.YES)
      self.canvasList[ind].mpl_connect("motion_notify_event",self.setPlotPosition)
      rangeFrame=Tk.Frame(self.plotFrames[ind])
      rangeFrame.pack(side=Tk.TOP)
      lab1=Tk.Label(rangeFrame,text='x min')
      lab1.grid(row=0,column=0)
      self.rangeEntry[0]=Tk.Entry(rangeFrame,width=10)
      self.rangeEntry[0].grid(row=0,column=1)
      lab1=Tk.Label(rangeFrame,text='y min')
      lab1.grid(row=1,column=0)
      self.rangeEntry[2]=Tk.Entry(rangeFrame,width=10)
      self.rangeEntry[2].grid(row=1,column=1)
      lab1=Tk.Label(rangeFrame,text='x max')
      lab1.grid(row=0,column=2)
      self.rangeEntry[1]=Tk.Entry(rangeFrame,width=10)
      self.rangeEntry[1].grid(row=0,column=3)
      lab1=Tk.Label(rangeFrame,text='y max')
      lab1.grid(row=1,column=2)
      self.rangeEntry[3]=Tk.Entry(rangeFrame,width=10)
      self.rangeEntry[3].grid(row=1,column=3)
      Tk.Button(rangeFrame,text="Apply Range",command=self.applyRange).grid(row=0,column=4)
      Tk.Button(rangeFrame,text="Default Range",command=self.defaultRange).grid(row=1,column=4)
      Tk.Button(rangeFrame,text="Invert y",command=self.inverty).grid(row=0,column=5)
      buttonFrame=Tk.Frame(self.plotFrames[ind])
      buttonFrame.pack(side=Tk.TOP)
      Tk.Button(buttonFrame,text="Save as PS",command=self.makePS).pack(side=Tk.LEFT)
      Tk.Button(buttonFrame,text="Save as PNG",command=self.makePNG).pack(side=Tk.LEFT)
      Tk.Button(buttonFrame,text="Write Plot Values",command=self.onPrint).pack(side=Tk.LEFT)
      if ind == 1:
        Tk.Button(buttonFrame,text="Plot Other Transformation",command=self.altPlot).pack(side=Tk.LEFT)
      Tk.Button(self.plotFrames[ind],text="Close Window",command=self.plotFrames[ind].withdraw).pack()
    else:
      self.plotFrames[ind].deiconify()
      self.subplotList[ind].clear()
    if self.xdata[ind] is None:
      self.putMessage('Error: no data values to plot.\n',self.messageText)
      return
    elif self.ydata[ind] is None:
      self.putMessage('Error: no data values to plot.\n',self.messageText)
      return
    try:
      self.rangeEntry[0].delete(0,Tk.END)
      self.rangeEntry[1].delete(0,Tk.END)
      self.rangeEntry[2].delete(0,Tk.END)
      self.rangeEntry[3].delete(0,Tk.END)
      self.rangeEntry[0].insert(0,str(self.dataLimits[0,ind]))
      self.rangeEntry[1].insert(0,str(self.dataLimits[1,ind]))
      self.rangeEntry[2].insert(0,str(self.dataLimits[2,ind]))
      self.rangeEntry[3].insert(0,str(self.dataLimits[3,ind]))
    except:
      pass
    if ind == 0:
      self.subplotList[ind].plot(self.xdata[ind],self.ydata[ind],'o',markersize=1.,color='blue')
      self.subplotList[ind].invert_yaxis()
    if ind == 1:
      self.subplotList[ind].plot(self.xdata[ind],self.ydata[ind],'o',markersize=5.,color='red')
      xmin=numpy.min(self.xdata[ind])-0.5
      xmax=numpy.max(self.xdata[ind])+0.5
      delx=(xmax-xmin)/1000.
      xmodel=numpy.arange(xmin,xmax,delx)
      ymodel=legendre.legval(xmodel,self.fitResults[self.jwstInds[2]])
      self.subplotList[ind].plot(xmodel,ymodel,color='cyan')
    if ind == 2:
      self.subplotList[ind].plot(self.xdata[ind],self.ydata[ind],'o',markersize=1.,color='blue')
      self.subplotList[ind].invert_yaxis()
    self.plotLimits[0,ind],self.plotLimits[1,ind]=self.subplotList[ind].get_xlim()
    self.plotLimits[2,ind],self.plotLimits[3,ind]=self.subplotList[ind].get_ylim()
    self.subplotList[ind].set_xlabel(self.xlabel[ind])
    self.subplotList[ind].set_ylabel(self.ylabel[ind])
    self.canvasList[ind].show()
    return

  def altPlot(self):
    self.subplotList[1].clear()
    ind=self.setOption.get()
    if ind == 0:
      xd1=self.kuruczMagValues[:,self.jwstInds[0]]-self.kuruczMagValues[:,self.jwstInds[1]]
      if self.yopt == 0:
        yd1=self.kuruczMagValues[:,self.jwstInds[0]]-self.kuruczMagValues[:,self.jwstInds[3]]
      else:
        yd1=self.kuruczMagValues[:,self.jwstInds[1]]-self.kuruczMagValues[:,self.jwstInds[3]]
    if ind == 1:
      xd1=self.phoenixMagValues[:,self.jwstInds[0]]-self.phoenixMagValues[:,self.jwstInds[1]]
      if self.yopt == 0:
        yd1=self.phoenixMagValues[:,self.jwstInds[0]]-self.phoenixMagValues[:,self.jwstInds[3]]
      else:
        yd1=self.phoenixMagValues[:,self.jwstInds[1]]-self.phoenixMagValues[:,self.jwstInds[3]]
    if ind == 2:
      xd1=self.blackbodyMagValues[:,self.jwstInds[0]]-self.blackbodyMagValues[:,self.jwstInds[1]]
      if self.yopt == 0:
        yd1=self.blackbodyMagValues[:,self.jwstInds[0]]-self.blackbodyMagValues[:,self.jwstInds[3]]
      else:
        yd1=self.blackbodyMagValues[:,self.jwstInds[1]]-self.blackbodyMagValues[:,self.jwstInds[3]]
    if ind == 3:
      xd1=self.boszMagValues[:,self.nirissInds[0]]-self.boszMagValues[:,self.nirissInds[1]]
      if self.yopt == 0:
        yd1=self.boszMagValues[:,self.nirissInds[0]]-self.boszMagValues[:,self.nirissInds[3]]
      else:
        yd1=self.boszMagValues[:,self.nirissInds[1]]-self.boszMagValues[:,self.nirissInds[3]]
    self.subplotList[1].plot(xd1,yd1,'o',markersize=5.,color='red')
    xmin=numpy.min(xd1)-0.5
    xmax=numpy.max(xd1)+0.5
    delx=(xmax-xmin)/1000.
    xmodel=numpy.arange(xmin,xmax,delx)
    ymodel=legendre.legval(xmodel,self.fitResults[self.jwstInds[3]])
    self.subplotList[1].plot(xmodel,ymodel,color='cyan')
    self.subplotList[1].set_xlabel(self.xlabel[1])
    filter1=self.mag1box.get()
    filter2=self.mag2box.get()
    filter3=self.mag3box.get()
    filter4=self.mag4box.get()
    yl1=filter1+' - '+filter4
    self.subplotList[1].set_ylabel(yl1)
    self.canvasList[1].show()
  
  def inverty(self):
    self.subplotList[self.plotOption-1].invert_yaxis()
    self.canvasList[self.plotOption-1].show()
  
  def applyRange(self):
    try:
      xmin=float(self.rangeEntry[0].get())
      xmax=float(self.rangeEntry[1].get())
      ymin=float(self.rangeEntry[2].get())
      ymax=float(self.rangeEntry[3].get())
      self.subplotList[self.plotOption-1].set_xlim(xmin,xmax)
      self.subplotList[self.plotOption-1].set_ylim(ymin,ymax)
      if self.yinvert:
        self.subplotList[self.plotOption-1].invert_yaxis()
      self.canvasList[self.plotOption-1].show()
    except:
      pass

  def defaultRange(self):
    try:
      xmin=self.plotLimits[0,self.plotOption-1]
      xmax=self.plotLimits[1,self.plotOption-1]
      ymin=self.plotLimits[2,self.plotOption-1]
      ymax=self.plotLimits[3,self.plotOption-1]
      self.subplotList[self.plotOption-1].set_xlim(xmin,xmax)
      self.subplotList[self.plotOption-1].set_ylim(ymin,ymax)
      self.canvasList[self.plotOption-1].show()
      self.rangeEntry[0].delete(0,Tk.END)
      self.rangeEntry[1].delete(0,Tk.END)
      self.rangeEntry[2].delete(0,Tk.END)
      self.rangeEntry[3].delete(0,Tk.END)
      self.rangeEntry[0].insert(0,str(self.plotLimits[0,self.plotOption-1]))
      self.rangeEntry[1].insert(0,str(self.plotLimits[1,self.plotOption-1]))
      if self.plotLimits[2,self.plotOption-1] < self.plotLimits[3,self.plotOption-1]:
        self.rangeEntry[2].insert(0,str(self.plotLimits[2,self.plotOption-1]))
        self.rangeEntry[3].insert(0,str(self.plotLimits[3,self.plotOption-1]))
      else:
        self.rangeEntry[2].insert(0,str(self.plotLimits[3,self.plotOption-1]))
        self.rangeEntry[3].insert(0,str(self.plotLimits[2,self.plotOption-1]))
    except:
      pass
  
  def onPrint(self):
    if not self.haveData:
      return
    outfilename=tkFileDialog.asksaveasfilename(filetypes=[('data files','*.data'),('data files','*.dat'),('data files','*.txt'),('All files','*')])
    if isinstance(outfilename,basestring):
      try:
        outfile=open(outfilename,"w")
        if outfile != None and outfile != '':
          print("# x: %s" % self.xlabel[self.plotOption-1], file=outfile)
          print("# y: %s" % self.ylabel[self.plotOption-1], file=outfile)
          for i in range(len(self.xdata[self.plotOption-1])):
            print('%12.5f %12.5f' % (self.xdata[self.plotOption-1][i],self.ydata[self.plotOption-1][i]),file=outfile)
        outfile.close()
      except:
        self.messageText.insert(Tk.END,"Error writing file: "+outfilename+"\n")
        self.messageText.see(Tk.END)

  def makePS(self):
    if not self.haveData:
      return
    outfile=tkFileDialog.asksaveasfilename(filetypes=[('postscript','*.ps')])
    if isinstance(outfile,basestring) and outfile != '':
      self.mplfigList[self.plotOption-1].savefig(outfile,format="PS")

  def makePNG(self):
    if not self.haveData:
      return
    outfile=tkFileDialog.asksaveasfilename(filetypes=[('PNG','*.png')])
    if isinstance(outfile,basestring) and outfile != '':
      self.mplfigList[self.plotOption-1].savefig(outfile,format="PNG")

#  def clearPlot(self):
#    self.subplotList[self.plotOption-1].clear()
#    self.canvasList[self.plotOption-1].show()
  
  def setPlotPosition(self,event):
    try:
      s1="Position: %.6g %.6g" % (event.xdata,event.ydata)
      self.posLabelText[self.plotOption-1].set(s1)
    except:
      pass
    
  def readData(self):
    try:
      xindex=int(self.mag1column.get())-1
      yindex=int(self.mag2column.get())-1
      raindex=int(self.racolumn.get())-1
      decindex=int(self.deccolumn.get())-1
      if xindex == yindex or xindex < 0 or yindex < 0:
        self.putMessage('Error: bad column numbers defined.  Check the inputs.\n',self.messageText)
        return
    except:
      self.putMessage('Error: the magnitude columns need to be defined properly.  Check the inputs.\n',self.messageText)
      return
    opt1=self.mag1option.get()
    opt2=self.mag2option.get()
    if opt1 == 1 and opt2 == 1:
      self.putMessage('Error: the magnitude columns need to be defined properly.  Two colours will not work.  Check the inputs.\n',self.messageText)
      return
    opt3=self.yAxisValue.get()
    self.haveTransformation=False
    ravalues=None
    decvalues=None
    filename=tkFileDialog.askopenfilename(filetypes=[('data files','*.data'),('data files','*.dat'),('data files','*.txt'),('data files','*.out'),('data files','*.tbl'),('fits files','*.fits'),('All files','*')])
    if isinstance(filename,type('string')):
      if '.fits' in filename:
        hdulist=fits.open(filename)
        data1=hdulist[1].data
        try:
          indata1=data1.field(xindex)
          indata2=data1.field(yindex)
          if raindex >= 0 and decindex >= 0:
            try:
              ravalues=data1.field(raindex)
              decvalues=data1.field(decindex)
              print('Have read in the RA/Dec values.')
            except:
              ravalues=None
              decvalues=None
              s1='Error trying to read the RA/Dec values from columns %d and %d of file %s.  These will not be used.' % (raindex,decindex,filename)
              self.putMessage(s1,self.messageText)
          hdulist.close()
        except:
          hdulist.close()
          s1='Error trying to read columns %i and %i in file %s.\nPlease check your inputs.\n' % (xindex,yindex,filename)
          self.putMessage(s1,self.messageText)
          return
      else:
        try:
          indata1=numpy.loadtxt(filename,usecols=(xindex,))
          indata2=numpy.loadtxt(filename,usecols=(yindex,))
        except:
          try:
            indata1=numpy.loadtxt(filename,usecols=(xindex,),comments=['#','\\','|'])
            indata2=numpy.loadtxt(filename,usecols=(yindex,),comments=['#','\\','|'])
          except:
            s1='Error trying to read columns %i and %i in file %s.\nPlease check your inputs.\n' % (xindex,yindex,filename)
            self.putMessage(s1,self.messageText)
            return
        if raindex >= 0 and decindex >= 0:
          try:
            ravalues=self.gettext(filename,raindex)
            decvalues=self.gettext(filename,decindex)
          except:
            ravalues=None
            decvalues=None
            s1='Error trying to read the RA/Dec values from columns %d and %d of file %s.  These will not be used.' % (raindex,decindex,filename)
            self.putMessage(s1,self.messageText)
      if not ravalues is None:
        n1=len(ravalues)
        n2=len(decvalues)
        n3=len(indata1)
        if n1 != n2 or n1 != n3:
          ravalues=None
          decvalues=None
          s1='Error: mis-match of RA/Dec values with photometry.  RA and Dec will not be used.'
          self.putMessage(s1,self.messageText)
      label1=self.mag1box.get()
      label2=self.mag2box.get()
      xdata,ydata,xlabel,ylabel=self.makexy(opt1,opt2,opt3,label1,label2,indata1,indata2)
      self.haveData=True
      xmin,xmax,ymin,ymax=self.getRange(xdata,ydata)
      self.dataLimits[0,0]=xmin
      self.dataLimits[1,0]=xmax
      self.dataLimits[2,0]=ymin
      self.dataLimits[3,0]=ymax
      npoints=len(xdata)
      s1='Have read in %d points from file %s\n' % (npoints,filename)
      self.putMessage(s1,self.messageText)
      self.xdata[0]=numpy.copy(xdata)
      self.ydata[0]=numpy.copy(ydata)
      self.xlabel[0]=xlabel
      self.ylabel[0]=ylabel
      self.ravalues=ravalues
      self.decvalues=decvalues
      for i in range(3):
        if self.subplotList[i] is None:
          pass
        else:
          self.subplotList[i].clear()

  def makexy(self,opt1,opt2,opt3,label1,label2,indata1,indata2):
    if opt1 == 0 and opt2 == 0:
      xdata=indata1-indata2
      xlabel=label1+' - '+label2
      if opt3 == 0:
        ydata=numpy.copy(indata1)
        ylabel=label1
      else:
        ydata=numpy.copy(indata2)
        ylabel=label2
    if opt1 == 0 and opt2 == 1:
      if opt3 == 0:
        ydata=numpy.copy(indata1)
        xdata=numpy.copy(indata2)
        xlabel=label1+' - '+label2
        ylabel=label1
      else:
        ydata=numpy.copy(indata2+indata1)
        xlabel=label1+' - '+label2
        ylabel=label2
        xdata=numpy.copy(indata2)
    if opt1 == 1 and opt2 == 0:
      if opt3 == 1:
        ydata=numpy.copy(indata2)
        xdata=numpy.copy(indata1)
        xlabel=label2+' - '+label1
        ylabel=label2
      else:
        ydata=numpy.copy(indata1+indata2)
        xdata=numpy.copy(indata1)
        xlabel=label2+' - '+label1
        ylabel=label1
    return xdata,ydata,xlabel,ylabel
          
  def writeTransformation(self):
    if not self.haveTransformation:
      return
    outfilename=tkFileDialog.asksaveasfilename(filetypes=[('data files','*.data'),('data files','*.dat'),('data files','*.txt'),('All files','*')])
    if isinstance(outfilename,basestring):
      try:
        outfile=open(outfilename,'w')
        if self.yopt == 0:
          for n in range(59):
            self.writeParams(self.fitResults[n],jwstnames[n],n,outfile)
        if self.yopt == 1:
          self.writeParams(self.fitResults[self.jwstInds[2]],jwstnames[self.jwstInds[2]],self.jwstInds[n],outfile)
        if self.yopt == 2:
          self.writeParams(self.fitResults[self.jwstInds[2]],jwstnames[self.jwstInds[2]],self.jwstInds[2],outfile)
          self.writeParams(self.fitResults[self.jwstInds[3]],jwstnames[self.jwstInds[3]],self.jwstInds[3],outfile)
        outfile.close()
      except:
        self.messageText.insert(Tk.END,"Error writing file: "+outfilename+"\n")
        self.messageText.see(Tk.END)

  def writeParams(self,fitvalues,filtername,n,outfile):
    print("# Legendre polynomial fit for filter %s" % (filtername),file=outfile)
    if self.setOption.get() == 0 | self.setOption.get() == 2:
      label=self.kuruczModelMagLabels[self.jwstInds[0]]+' - '+filtername
    else:
      label=self.phoenixModelMagLabels[self.jwstInds[0]]+' - '+filtername
    print("# colour 1 (x): %s " % (self.xlabel[1]), file=outfile)
    print("# colour 2 (y): %s " % (label), file=outfile)
    print("# ",file=outfile)
    for n in range(len(fitvalues)):
      print("# order %d coefficient: %f" % (n,fitvalues[n]),file=outfile)
    print("# ",file=outfile)
    print("# range of values in fit (x): %f to %f " % (numpy.min(self.xdata[1]),numpy.max(self.xdata[1])),file=outfile)
    print('# \n# Fitted values:',file=outfile)
    xmin=numpy.min(self.xdata[1])-0.5
    xmax=numpy.max(self.xdata[1])+0.5
    delx=(xmax-xmin)/1000.
    xmodel=numpy.arange(xmin,xmax,delx)
    ymodel=legendre.legval(xmodel,fitvalues)
    for i in range(len(xmodel)):
      print('%10.5f %10.5f' % (xmodel[i],ymodel[i]),file=outfile)
        
  def getRange(self,xdata,ydata):
    xmin=numpy.min(xdata)
    xmax=numpy.max(xdata)
    ymin=numpy.min(ydata)
    ymax=numpy.max(ydata)
    xmin=self.roundFloat(xmin,0)
    xmax=self.roundFloat(xmax,1)
    ymin=self.roundFloat(ymin,0)
    ymax=self.roundFloat(ymax,1)
    return xmin,xmax,ymin,ymax
        
  def onExit(self):
    self.root.quit()

  def fit1(self,magopt,setopt,yopt,norder,mopt1,mopt2,mopt3,mopt4,interactive):
    ndatapoints=len(self.xdata[0])
    self.jwstMags=numpy.zeros((ndatapoints,59),dtype=numpy.float32)
    self.yopt=yopt
    if setopt == 0:
      if magopt == 0:
        for n in range(59):
          self.doFit(self.kuruczMagValues[:,mopt1],self.kuruczMagValues[:,mopt2],self.kuruczMagValues[:,n],n,yopt,norder,interactive)
      if magopt == 1:
        self.doFit(self.kuruczMagValues[:,mopt1],self.kuruczMagValues[:,mopt2],self.kuruczMagValues[:,mopt3],mopt3,yopt,norder,interactive)
      if magopt == 2:
        self.doFit(self.kuruczMagValues[:,mopt1],self.kuruczMagValues[:,mopt2],self.kuruczMagValues[:,mopt3],mopt3,yopt,norder,interactive)
        self.doFit(self.kuruczMagValues[:,mopt1],self.kuruczMagValues[:,mopt2],self.kuruczMagValues[:,mopt4],mopt4,yopt,norder,interactive)
      self.xdata[1]=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt2]
      if yopt == 0:
        self.ydata[1]=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt3]
      else:
        self.ydata[1]=self.kuruczMagValues[:,mopt2]-self.kuruczMagValues[:,mopt3]
    if setopt == 1:
      if magopt == 0:
        for n in range(59):
          self.doFit(self.phoenixMagValues[:,mopt1],self.phoenixMagValues[:,mopt2],self.phoenixMagValues[:,n],n,yopt,norder,interactive)
      if magopt == 1:
        self.doFit(self.phoenixMagValues[:,mopt1],self.phoenixMagValues[:,mopt2],self.phoenixMagValues[:,mopt3],mopt3,yopt,norder,interactive)
      if magopt == 2:
        self.doFit(self.phoenixMagValues[:,mopt1],self.phoenixMagValues[:,mopt2],self.phoenixMagValues[:,mopt3],mopt3,yopt,norder,interactive)
        self.doFit(self.phoenixMagValues[:,mopt1],self.phoenixMagValues[:,mopt2],self.phoenixMagValues[:,mopt4],mopt4,yopt,norder,interactive)
      self.xdata[1]=self.phoenixMagValues[:,mopt1]-self.phoenixMagValues[:,mopt2]
      self.ydata[1]=self.phoenixMagValues[:,mopt1]-self.phoenixMagValues[:,mopt3]
    if setopt == 2:
      if magopt == 0:
        for n in range(59):
          self.doFit(self.blackbodyMagValues[:,mopt1],self.blackbodyMagValues[:,mopt2],self.blackbodyMagValues[:,n],n,yopt,norder,interactive)
      if magopt == 1:
        self.doFit(self.blackbodyMagValues[:,mopt1],self.blackbodyMagValues[:,mopt2],self.blackbodyMagValues[:,mopt3],mopt3,yopt,norder,interactive)
      if magopt == 2:
        self.doFit(self.blackbodyMagValues[:,mopt1],self.blackbodyMagValues[:,mopt2],self.blackbodyMagValues[:,mopt3],mopt3,yopt,norder,interactive)
        self.doFit(self.blackbodyMagValues[:,mopt1],self.blackbodyMagValues[:,mopt2],self.blackbodyMagValues[:,mopt4],mopt4,yopt,norder,interactive)
      self.xdata[1]=self.blackbodyMagValues[:,mopt1]-self.blackbodyMagValues[:,mopt2]
      self.ydata[1]=self.blackbodyMagValues[:,mopt1]-self.blackbodyMagValues[:,mopt3]
    if setopt == 3:
      if magopt == 0:
        for n in range(12):
          self.doFit(self.boszMagValues[:,mopt1],self.boszMagValues[:,mopt2],self.boszMagValues[:,n],n,yopt,norder,interactive)
      if magopt == 1:
        self.doFit(self.boszMagValues[:,mopt1],self.boszMagValues[:,mopt2],self.boszMagValues[:,mopt3],mopt3,yopt,norder,interactive)
      if magopt == 2:
        self.doFit(self.boszMagValues[:,mopt1],self.boszMagValues[:,mopt2],self.boszMagValues[:,mopt3],mopt3,yopt,norder,interactive)
        self.doFit(self.boszMagValues[:,mopt1],self.boszMagValues[:,mopt2],self.boszMagValues[:,mopt4],mopt4,yopt,norder,interactive)
      self.xdata[1]=self.boszMagValues[:,mopt1]-self.boszMagValues[:,mopt2]
      self.ydata[1]=self.boszMagValues[:,mopt1]-self.boszMagValues[:,mopt3]
    
  def autoTransform(self,parameters):
    try:
      if parameters['Output_Filter_Values']['modelset'] == 'Kurucz':
        setopt=0
        self.kuruczMagValues,self.kuruczModelMagLabels,self.kuruczFilterPars=self.readMagslist(os.path.join(path,'magslist_old_kurucz.new'),142)
        if self.kuruczMagValues is None:
          self.exit()
      if parameters['Output_Filter_Values']['modelset'] == 'Phoenix':
        setopt=1
        self.phoenixMagValues,self.phoenixModelMagLabels,self.phoenixFilterPars=self.readMagslist(os.path.join(path,'magslist_phoenix_grid.new'),121)
        if self.phoenixMagValues is None:
          self.exit()
      if parameters['Output_Filter_Values']['modelset'] == 'blackbody':
        setopt=2
        self.blackbodyMagValues,self.blackbodyMagLabels,self.blackbodyFilterPars=self.readMagslist(os.path.join(path,'magslist_backbody.new'),142)
        if self.blackbodyMagValues is None:
          self.exit()
      if parameters['Output_Filter_Values']['modelset'] == 'BOSZ':
        setopt=3
        self.boszMagValues,self.boszModelMagLabels,self.boszFilterPars=self.readMagslist(os.path.join(path,'magslist_bosz_normal.new'),142)
        if self.boszMagValues is None:
          self.exit()
      if parameters['Input_Magnitude_Parameters']['column1type'] == 'magnitude':
        opt1=0
      else:
        opt1=1
      if parameters['Input_Magnitude_Parameters']['column2type'] == 'magnitude':
        opt2=0
      else:
        opt2=1
      filter1=parameters['Input_Magnitude_Parameters']['filter1'].replace('_',' ')
      filter2=parameters['Input_Magnitude_Parameters']['filter2'].replace('_',' ')
      filter3=parameters['Output_Filter_Values']['jwst1'].replace('_',' ')
      filter4=parameters['Output_Filter_Values']['jwst2'].replace('_',' ')
      norder=int(parameters['Output_Filter_Values']['fitorder'])
      if norder < 2:
        sys.exit()
      mopt1,mopt2,mopt3,mopt4=self.matchFilter(setopt,filter1,filter2,filter3,filter4)
      if mopt1 < 0 or mopt2 < 0 or mopt3 < 0 or mopt4 < 0:
        sys.exit()
      self.jwstInds=[mopt1,mopt2,mopt3,mopt4]
      xindex=int(parameters['Input_Magnitude_Parameters']['column1'])-1
      yindex=int(parameters['Input_Magnitude_Parameters']['column2'])-1
      raindex=int(parameters['Input_Magnitude_Parameters']['racolumn'])-1
      decindex=int(parameters['Input_Magnitude_Parameters']['deccolumn'])-1
      if xindex == yindex or xindex < 0 or yindex < 0:
        sys.exit()
      if opt1 == 1 & opt2 == 1:
        sys.exit()
      yaxis=int(parameters['Input_Magnitude_Parameters']['yvalue'])-1
      filename=parameters['Input_Magnitude_Parameters']['datafile']
      if '.fits' in filename[-5:]:
        hdulist=fits.open(filename)
        data1=hdulist[1].data
        indata1=data1.field(xindex)
        indata2=data1.field(yindex)
        if raindex >= 0 and decindex >= 0:
          ravalues=data1.field(raindex)
          decvalues=data1.field(decindex)
        else:
          ravalues=None
          decvalues=None
      else:
        try:
          indata1=numpy.loadtxt(filename,usecols=(xindex,))
          indata2=numpy.loadtxt(filename,usecols=(yindex,))
        except:
          indata1=numpy.loadtxt(filename,usecols=(xindex,),comments=['#','\\','|'])
          indata2=numpy.loadtxt(filename,usecols=(yindex,),comments=['#','\\','|'])
        if raindex >= 0 and decindex >= 0:
          ravalues=self.gettext(filename,raindex)
          decvalues=self.gettext(filename,decindex)
        else:
          ravalues=None
          decvalues=None
      opt3=0
      xdata,ydata,xlabel,ylabel=self.makexy(opt1,opt2,opt3,filter1,filter2,indata1,indata2)
      self.xdata[0]=numpy.copy(xdata)
      self.ydata[0]=numpy.copy(ydata)
      self.xlabel[0]=xlabel
      self.ylabel[0]=ylabel
      if filter1 == filter2:
        magopt=1
      else:
        magopt=2
      ndatapoints=len(self.xdata[0])
      self.jwstMags=numpy.zeros((ndatapoints,59),dtype=numpy.float32)
      self.fit1(magopt,setopt,yaxis,norder,mopt1,mopt2,mopt3,mopt4,False)
      outfilename=parameters['Output_Filter_Values']['outfilename']
      self.ravalues=ravalues
      self.decvalues=decvalues
      self.writeMags(outfilename)
      print('Output file %s has been written.' % (outfilename))
      xmin=numpy.min(self.xdata[1])-0.5
      xmax=numpy.max(self.xdata[1])+0.5
      delx=(xmax-xmin)/1000.
      xmodel=numpy.arange(xmin,xmax,delx)
      setnames=['kurucz','phoenix','blackbody','bosz']
      for n in range(2):
        if n == yaxis:
          ymodel=legendre.legval(xmodel,self.fitResults[self.jwstInds[2]])
          matplotlib.pyplot.plot(self.xdata[1],self.ydata[1],'o',markersize=1.,color='red')
        else:
          if setopt == 0:
            xd1=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt2]
            if yaxis == 0:
              yd1=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt4]
            else:
              yd1=self.kuruczMagValues[:,mopt2]-self.kuruczMagValues[:,mopt4]
            matplotlib.pyplot.plot(xd1,yd1,'o',markersize=1.,color='red')
          if setopt == 1:
            xd1=self.phoenixMagValues[:,mopt1]-self.phoenixMagValues[:,mopt2]
            if yaxis == 0:
              yd1=self.phoenixMagValues[:,mopt1]-self.phoenixMagValues[:,mopt4]
            else:
              yd1=self.phoenixMagValues[:,mopt2]-self.phoenixMagValues[:,mopt4]
            matplotlib.pyplot.plot(xd1,yd1,'o',markersize=1.,color='red')
          if setopt == 2:
            xd1=self.blackbodyMagValues[:,mopt1]-self.blackbodyMagValues[:,mopt2]
            if yaxis == 0:
              yd1=self.blackbodyMagValues[:,mopt1]-self.blackbodyMagValues[:,mopt4]
            else:
              yd1=self.blackbodyMagValues[:,mopt2]-self.blackbodyMagValues[:,mopt4]
            matplotlib.pyplot.plot(xd1,yd1,'o',markersize=1.,color='red')
          if setopt == 4:
            xd1=self.boszMagValues[:,mopt1]-self.boszMagValues[:,mopt2]
            if yaxis == 0:
              yd1=self.boszMagValues[:,mopt1]-self.boszMagValues[:,mopt4]
            else:
              yd1=self.boszMagValues[:,mopt2]-self.boszMagValues[:,mopt4]
            matplotlib.pyplot.plot(xd1,yd1,'o',markersize=1.,color='red')
          ymodel=legendre.legval(xmodel,self.fitResults[self.jwstInds[3]])
          xd1=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt2]
          if yaxis == 0:
            yd1=self.kuruczMagValues[:,mopt1]-self.kuruczMagValues[:,mopt4]
          else:
            yd1=self.kuruczMagValues[:,mopt2]-self.kuruczMagValues[:,mopt4]
        matplotlib.pyplot.plot(xmodel,ymodel,color='cyan')
        if n == 0:
          matplotlib.pyplot.ylabel(ylabel+' - JWST '+filter3)
          figname='transformation '+xlabel+' to '+ylabel+' - '+filter3+' '+setnames[setopt]+'.png'
        else:
          matplotlib.pyplot.ylabel(ylabel+' - JWST '+filter4)
          figname='transformation '+xlabel+' to '+ylabel+' - '+filter4+' '+setnames[setopt]+'.png'
        matplotlib.pyplot.xlabel(xlabel)
        figname=figname.replace(' ','_')
        figname=figname.replace('-','minus')
        matplotlib.pyplot.savefig(figname)
        matplotlib.pyplot.close()
    except:
      print('Error trying to do the transformation.')
      sys.exit()
    
def parseArguments(argv):
  sectionlist=['Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Input_Magnitude_Parameters','Output_Filter_Values','Output_Filter_Values','Output_Filter_Values','Output_Filter_Values','Output_Filter_Values']
  keylist=['filter1','filter2','column1','column2','column1type','column2type','yvalue','datafile','racolumn','deccolumn','jwst1','jwst2','modelset','fitorder','outfilename']
  try:
    filename=argv[-1]
    if '.cfg' in filename: 
      print('Reading cfg file %s' % filename)   
      config=ConfigObj(filename)
      for i in range(len(sectionlist)):
        value=config[sectionlist[i]][keylist[i]]
      return config
    else:
      return None
  except KeyError:
    raise RuntimeError('Failed to parse {}'.format(filename))

def main(argv):
  # print('Calling main function')
  parameters=parseArguments(argv)
  print('Parameters', parameters)
  if parameters is None:
    root=Tk.Tk()
    root.title("JWST Magnitude Simulation Tool")
    x=magConGUI()
    x.runGUI(root)
    root.mainloop()
  else:
    root=None
    x=magConGUI()
    x.autoTransform(parameters)
      
if __name__ == "__main__":
  main(sys.argv)
