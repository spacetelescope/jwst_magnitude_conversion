#! /usr/bin/env python
#
"""
This program is for looking at colour-colour relatiions based on the simulated
magnitude files from magnitudes.py, which are input to jwst_magnitude_converter.py.
The code here allows approximate extinction effects to be included, just based 
on a mean extinction relative to V in each filter.

Some of the code is taken from jwst_magnitude_converter.py.  This code only 
works with the simulated magnitudes and not with actual data.
"""
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
import numpy 
import numpy.polynomial.legendre as legendre
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import astropy.io.fits as fits
from configobj import ConfigObj
import sys
import os
import tkinter as Tk
import tkinter.ttk as ttk
from tkinter.scrolledtext import ScrolledText
import tkinter.filedialog as tkFileDialog
from tkinter.messagebox import Message as tkMessageBox
from tkinter.colorchooser import Chooser as tkColorChooser
kurucz_filter_names = [
"NIRISS F090W ","NIRISS F115W ","NIRISS F140M ","NIRISS F150W ",
"NIRISS F158M ","NIRISS F200W ","NIRISS F277W ","NIRISS F356W ",
"NIRISS F380M ","NIRISS F430M ","NIRISS F444W ","NIRISS F480M ",
"Guider 1","Guider 2",
"NIRCam F070W ","NIRCam F090W ","NIRCam F115W ","NIRCam F140M ",
"NIRCam F150W ","NIRCam F150W2 ","NIRCam F162M ","NIRCam F164N ",
"NIRCam F182M ","NIRCam F187M ","NIRCam F200W ","NIRCam F210M ",
"NIRCam F212N ","NIRCam F250M ","NIRCam F277W ","NIRCam F300M ",
"NIRCam F322W2 ","NIRCam F323N ","NIRCam F335M ","NIRCam F356W ",
"NIRCam F360M ","NIRCam F405N ","NIRCam F410M ","NIRCam F430M ",
"NIRCam F444W ","NIRCam F460M ","NIRCam F466N ","NIRCam F470N ",
"NIRCam F480M ","MIRI F560W ","MIRI F770W ","MIRI F1000W ",
"MIRI F1280W ","MIRI F1130W ","MIRI F1500W ","MIRI F1800W ",
"MIRI F2100W ","MIRI F2550W ","NIRSpec F110W","NIRSpec F140X",
"NIRSpec Clear","NIRSpec F070LP",
"NIRSec F100LP","NIRSpec F170LP","NIRSPec F290LP",
"Sloan u  ","Sloan g  ","Sloan r  ","Sloan i  ",
"Sloan z  ","Bessell U  ","Bessell B  ","Bessell V  ",
"Bessell/Cousins R  ","Bessell/Cousins I  ","Johnson U  ",
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
"GAIA g","GAIA gbp","GAIA grp"]
phoenix_filter_names = [
'NIRISS F090W','NIRISS F115W','NIRISS F140M','NIRISS F150W','NIRISS F158M',
'NIRISS F200W','NIRISS F277W','NIRISS F356W','NIRISS F380M','NIRISS F430M',
'NIRISS F444W','NIRISS F480M','Guider 1','Guider 2','NIRCam F070W',
'NIRCam F090W','NIRCam F115W','NIRCam F140M','NIRCam F150W','NIRCam F150W2',
'NIRCam F162M','NIRCam F164N','NIRCam F182M','NIRCam F187M','NIRCam F200W',
'NIRCam F210M','NIRCam F212N','NIRCam F250M','NIRCam F277W','NIRCam F300M',
'NIRCam F322W2','NIRCam F323N','NIRCam F335M','NIRCam F356W','NIRCam F360M',
'NIRCam F405N','NIRCam F410M','NIRCam F430M','NIRCam F444W','NIRCam F460M',
'NIRCam F466N','NIRCam F470N','NIRCam F480M','NIRSpec F110W','NIRSpec F140X',
'NIRSpec Clear','NIRSpec F070LP','NIRSpec F100LP','NIRSpec F170LP',
'NIRSpec F290LP','Sloan u','Sloan g','Sloan r',
'Sloan i','Sloan z','Bessell U','Bessell B',
'Bessell V','Bessell/Cousins R','Bessell/Cousins I',
'Johnson U','Johnson B','Johnson V','Johnson R',
'Johnson I','Johnson J','Johnson H','Johnson K',
'Johnson L','Johnson M','2MASS J','2MASS H','2MASS Ks',
'IRAC [3.6]','IRAC [4.4]','DENIS i','DENIS J','DENIS Ks','WISE W1','WISE W2',
'UKIDSS Z','UKIDSS Y','UKIDSS J','UKIDSS H',
'UKIDSS K','Pan-STARRS1 g','Pan-STARRS1 r',
'Pan-STARRS1 i','Pan-STARRS1 z','Pan-STARRS1 y',
'Pan-STARRS1 w','HST ACS F220W','HST ACS F250W','HST ACS F330W',
'HST ACS F435W','HST ACS F475W','HST ACS F555W','HST ACS F606W','HST ACS F625W',
'HST ACS F775W','HST ACS F814W','HST WFC3 F218W','HST WFC3 F225W',
'HST WFC3 F275W','HST WFC3 F336W','HST WFC3 F390W','HST WFC3 F438W',
'HST WFC3 F475W','HST WFC3 F555W','HST WFC3 F606W','HST WFC3 F625W',
'HST WFC3 F775W','HST WFC3 F814W','HST WFC3 F105W','HST WFC3 F110W',
'HST WFC3 F125W','HST WFC3 F140W','HST WFC3 F160W','GAIA g',
'GAIA gbp','GAIA grp']

class magnitudeTransform():
    """
This the driver class that brings up the user interface for the plotting 
and fitting of the simulated colour-colour relations.  It has no parameters 
and does not return anything.
    """
    def __init__(self,parent=None,**args):
        """
This routine initiallizes the object with some values.  One needs to use the 
run_gui routine to start the actual user interface.  The optional parameters 
parent and args are not actually used.
        """
        self.range_entry = [None,None,None,None]
# for the plot, xmin, xmax, ymin, ymax values saved here....
        self.plot_limits = numpy.zeros((4),dtype=numpy.float32)
# for the plot, xmin, xmax, ymin, ymax data values saved here...
        self.data_limits = numpy.zeros((4),dtype=numpy.float32)
# self.fitRange holds the x colour range of validity of the fitting for a given plot.
        self.fit_range = numpy.zeros((2),dtype=numpy.float32)
        self.plot_frame = None

    def run_gui(self,root):
        """
This routine starts the user interface for the program.  It requires as a 
parameter the parent "root" Tk window.  No values are returned.  If root is 
None then the code does not start the interface.  This provides the possibility
of using the code non-interactively.  For this particular program that option 
is not currently used.

Parameter

root  :  a Tk top-level window object, or None

Returns

No value is returned by the routine.

        """
        if root is None:
            pass
        else:
            self.make_widgets(root)
    
    def round_float(self,value,minimum_flag):
        """
This routine performs rounding of floating point values to convenient values 
for such things as calculating plot ranges.  

It takes the following parameters:

value  :            A floating point number to be rounded off

minimum_flag  :     A boolean value value to determine whether the value should
                    be rounded up or down.  If True values are rounded down as 
                    one wishes for the minimum value in a plot, and if False 
                    values are rounded up as one wishes for the maximum value 
                    in a plot.

There is a single return value:

rounded_value  :    The rounded off floating point number calculated from 
                    `value'.

The code differs from use of the floor and ceil functions in that it tries to 
round to the nearest significant figure for a given value.  So passing in 
a value of 1.2e-08, for exmaple, with minimum_flag = False returns 2.0e-08, 
and if the flag is True it returns 1.0e-08.

The code uses the math log10, floor, and ceil functions.

        """
        if value == 0.:
            return value
        else:
            if value < 0.:
                sign = -1.
                value = -value
            else:
                sign = 1.
            power = math.log10(value)
            if power < 0.:
                exp = int(power-1.)
            else:
                exp = 0.0
            shift = 10.**exp
            x = value/shift
            if x < 1.7:
                x=x*10.
                shift=shift/10.
            elif x < 2.5:
                x=x*5.
                shift=shift/5.
            if (minimum_flag)  and sign > 0.:
                x = math.floor(x)
            elif (minimum_flag)  and sign < 0.:
                x = math.ceil(x)
            elif (not minimum_flag) and sign > 0.:
                x = math.ceil(x)
            elif (not minimum_flag) and sign < 0.:
                x = math.floor(x)
            rounded_value = x*shift*sign
            # If the rounded value is very close to the input value, offset
            # by 5%...not normally part of the routine, but needed for
            # matplotlib plots.
            ratio = abs(value/rounded_value)
            if ratio > 0.97:
                if ((minimum_flag == 0) and rounded_value < 0.) or ((minimum_flag == 1) and rounded_value > 0.):
                    rounded_value = rounded_value*1.05
                else:
                    rounded_value = rounded_value/1.05
            return rounded_value

    def sep_line(self,parent,width_value,height_value,pad_value):
        """
This routine makes a separator line in Tkinter, and returns the Canvas object.
It is a convience routine if one is making lines as dividers in an interface, 
and is using "pack" to position items from top to bottom within a frame.  It 
is NOT a general routine for making lines in Tkinter.

Usage:
    
    line=self.sep_line(parent,width,height,pad)

Parameters:

    parent  :  Tk.Frame or Tk.Toplevel object

       The Tk object that will contain the line.

    width_value  :  int

       The desired width of the line in pixels.

    height_value  :  int 

       The desired height of the line in pixels.

    pad_value  :  int

       The desired padding of the line in pixels.

Returns:

    lincanvas  :  Tk.Canvas object

       The object created by the call to Tk.Canvas is returned.

`width_value', `height_value', and `pad_value' must be positive non-zero 
integers.

Note that the code does not check these values, it just passes them to the 
Tk.Canvas widget call to create the line object.

Also note that the code calls the "pack" function with side=Tk.TOP to make 
a separator line.  This is assuming that items are being positioned using 
pack not grid when the Tk window is made.  If one is using grid in the 
parent then this will cause the code to hang. 

        """
        line_canvas = Tk.Canvas(parent,height=height_value,width=width_value)
        try:
            line_canvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=Tk.YES)
        except:
            pass
        line_canvas.create_line(pad_value,height_value/2,width_value-pad_value,height_value/2)
        return line_canvas

    def make_widgets(self,root):
        """
This routine reads in the model magnitude values that are needed for the 
code to operate.  Assuming that this works, the code then makes the main 
Tkinter window used by the code for the interface.   Otherwise the code 
simply exits as nothing can be done without the model magnitude values.

Parameters:

    root  :  Tkinter Toplevel or main window object

Returns:

    No values are returned by the routine.



        """
        # Read in the simulated magnitude values needed for the calculations.
        # If this fails, the program cannot do anything so it exits.
        status = self.read_model_values('magslist_bosz_normal.new','magslist_old_kurucz.new','magslist_phoenix_grid.new','magslist_blackbody.new')
        if not status:
            print('Error reading in the model magnitude values.  Exiting.\n')
            sys.exit(1)
        self.root = root
        self.area_label = Tk.Label(root,text="Message Area")
        self.area_label.pack(side=Tk.TOP)
        self.message_text = ScrolledText(root,height=6,width=60,bd=1,relief=Tk.RIDGE,wrap=Tk.NONE)
        self.message_text.config(font=('courier',12,'bold'))
        self.message_text.pack()
        sep1=self.sep_line(root,450,16,10)
        label1=Tk.Label(root,text='Input filter names for the colour-colour plot:')
        label1.pack(side=Tk.TOP)
        select_frame = Tk.Frame(root)
        select_frame.pack(side=Tk.TOP)
        label1 = Tk.Label(select_frame,text='X Filter Magnitude 1')
        label1.grid(row=0,column=0)
        self.mag_boxes = []
        mag1box = ttk.Combobox(select_frame,width=20)
        mag1box.grid(row=1,column=0)
        mag1box['values'] = kurucz_filter_names
        mag1box.current(70)
        self.mag_boxes.append(mag1box)
        label2 = Tk.Label(select_frame,text='X Filter Magnitude 2')
        label2.grid(row=0,column=1)
        mag1box = ttk.Combobox(select_frame,width=20)
        mag1box.grid(row=1,column=1)
        mag1box['values'] = kurucz_filter_names
        mag1box.current(71)
        self.mag_boxes.append(mag1box)
        label3 = Tk.Label(select_frame,text = 'Y Filter Magnitude 1')
        label3.grid(row=2,column=0)
        mag1box = ttk.Combobox(select_frame,width=20)
        mag1box.grid(row=3,column=0)
        mag1box['values'] = kurucz_filter_names
        mag1box.current(0)
        self.mag_boxes.append(mag1box)
        label4 = Tk.Label(select_frame,text='Y Filter Magnitude 2')
        label4.grid(row=2,column=1)
        mag1box = ttk.Combobox(select_frame,width=20)
        mag1box.grid(row=3,column=1)
        mag1box['values'] = kurucz_filter_names
        mag1box.current(1)
        self.mag_boxes.append(mag1box)
        sep1 = self.sep_line(root,450,16,10)
        l1=Tk.Label(root,text='Values to use in plotting:')
        l1.pack()
        b1 = Tk.Frame(root)
        b1.pack()
        self.model_flags = []
        text = ['BOSZ','Kurucz','Phoenix','blackbody']
        for loop in range(4):
            ivar = Tk.IntVar()
            self.model_flags.append(ivar)
            cb = Tk.Checkbutton(b1,text=text[loop],var=ivar)
            cb.pack(side=Tk.LEFT)
            if loop < 3:
                ivar.set(1)
        self.plot_symbol_size = Tk.DoubleVar()
        self.symbol_size = Tk.Scale(root,orient=Tk.HORIZONTAL,to=6.0,from_=0.0,tickinterval=1.0,resolution=0.01,variable=self.plot_symbol_size,command=self.replot,length=350,label="Plot Symbol Size")
        self.symbol_size.pack()
        self.plot_symbol_size.set(3.0)
        self.fit_flag = Tk.IntVar()
        cb = Tk.Checkbutton(root,text='Plot Fit Values: ',var=self.fit_flag)
        cb.pack()
        l1 = Tk.Frame(root)
        l1.pack()
        label1 = Tk.Label(l1,text="Order of fit: ")
        label1.pack(side=Tk.LEFT)
        self.fit_order = Tk.Entry(l1,width=3)
        self.fit_order.pack()
        self.fit_order.insert(0,'4')
        button_frame = Tk.Frame(root)
        button_frame.pack(side=Tk.TOP)
        b1 = Tk.Button(button_frame,text="Plot Colour-Colour Diagram",command=self.plot_colour_colour)
        b1.pack(side=Tk.TOP)
        sep1 = self.sep_line(root,450,16,10)
        Tk.Button(root,text="Close Widget",command=self.root.quit).pack(side=Tk.TOP)
        self.put_message('Have read in the BOSZ Kurucz, Phoenix, and blackbody \nmodel magnitude values.\n',self.message_text)

    def replot(self,event):
        """
This routine is called when the marker size slider is changed.  It checks to see 
if the plot window is active, and if so it replots with the limits currently 
in the entry fields.  One has to do this indirectly because one can change the 
point size on the main window before the plot is made.

Parameters:
    
    event : TK event variable

        This is the event from the slider.  It is not currently used.

Returns:

    No values are returned from this routine.

        """
        if self.plot_frame is None:
            return
        else:
            self.make_plot(True)
        
    def put_message(self,string1,message_area):
        """
This routine posts a message string to a scrolled text message area.

Parameters:

    string1  :  string variable

        This is the string that is appended to the message area

   message_area  :   Tk.ScrolledText object

        This is the ScrolledText area that receives the text.

Returns:

   No values are returned by this routine.


        """
        message_area.insert(Tk.END,string1)
        message_area.see(Tk.END)
    
    def read_model_values(self,filename1,filename2,filename3,filename4):
        """
This routine reads in four sets of magnitude values simulated for a set of 
spectral shapes.  The first three are the BOSZ, Kurucz, and Phoenix stellar 
atmosphere model sets.  The last one is a blackbody model set.  The code 
places these magnitude values into variables within the object, so no values
are returned.

Parameters:

    filename1  :  string 

        The BOSZ magnitude list file name (without the path...)

    filename2  :  string 

        The Kurucz magnitude list file name (without the path...)

    filename3  :  string 

        The Phoenix grid magnitude list file name (without the path...)

    filename4  :  string 

        The blackbody magnitude list file name (without the path...)

All these files are produced by the magnitudes.py code.  Contact Kevin Volk 
if details about this are needed.

        """
        # The path variable is used in reading in the simulated magnitude values.
        # If is either the path to the python code or the value of $SIMULATED_MAGNITUDES_PATH
        try:
            path = os.environ['SIMULATED_MAGNITUDES_PATH']
        except:
            path  =  os.path.dirname(__file__)
        if path[-1:] != '/':
            path = path+'/'
        try:
            bosz_magnitude_values,bosz_magnitude_labels,bosz_filter_parameters = self.read_magnitude_list(os.path.join(path,filename1),len(kurucz_filter_names))
            kurucz_magnitude_values,kurucz_magnitude_labels,kurucz_filter_parameters = self.read_magnitude_list(os.path.join(path,filename2),len(kurucz_filter_names))
            phoenix_magnitude_values,phoenix_magnitude_labels,phoenix_filter_parameters = self.read_magnitude_list(os.path.join(path,filename3),len(phoenix_filter_names))
            blackbody_magnitude_values,blackbody_magnitude_labels,blackbody_filter_parameters = self.read_magnitude_list(os.path.join(path,filename4),len(kurucz_filter_names))
            if kurucz_magnitude_values is None or phoenix_magnitude_values is None or blackbody_magnitude_values is None:
                self.magnitudes_values = [None,None,None,None]
                self.magnitudes_labels = [None,None,None,None]
                self.filter_parameters = [None,None,None,None]
                return False
            else:
                self.magnitude_values = [bosz_magnitude_values,kurucz_magnitude_values,phoenix_magnitude_values,blackbody_magnitude_values]
                self.magnitude_labels = [bosz_magnitude_labels,kurucz_magnitude_labels,phoenix_magnitude_labels,blackbody_magnitude_labels]
                self.filter_parameters = [bosz_filter_parameters,kurucz_filter_parameters,phoenix_filter_parameters,blackbody_filter_parameters]
                return True
        except:
            return False

    def get_ascii_text_column(self,input_file_name,column_number):
        """
This routine reads the numerical values from a specific column from an input 
ascii file and returns these as a numpy arrry.  If there is an issue reading 
the file, the it returns None.  The code uses numpy.loadtxt to read the input 
file specified; characters #, |, and \ are assumed to be comment characters.
Hence the code can read IPAC table files, for example.

Parameters:

    input_file_name  :  string 

        The name of the ascii file name to read from.

    column_number  :  int

        The column number from which to read values.  It must be a postive 
        integer value.  Number is 0-based as normal in python/numpy.

Returns:

    out_values  : numpy float array  or  None

        The array of values read from the file is returned as a numpy array if
        reading proceeds successfully.  Any error produces a refurn value of 
        None. 

        """
        if column_number < 0:
            return None
        try:
            infile = open(input_file_name,'r')
            lines = infile.readlines()
            infile.close()
            out_values = []
            for line in lines:
                if '#' in line[0:1] or '\\' in line[0:1] or '|' in line[0:1]:
                    pass
                else:
                    line = line.strip('\n')
                    values = line.split()
                    out_values.append(values[column_number])
            out_values = numpy.asarray(out_values)
            return out_values
        except:
            return None
  
    def read_magnitude_list(self,input_file_name,number_of_filters):
        """
The code here reads in a "standard" magslist.out file from magnitudes.py.  The 
file is assumed to have a header listing the run parameters on lines 1 to 5, 
then a list of the filter properties line by line, followed by the grid of 
output magnitudes.  The code first reads in the filter names and the filter 
properties (filter wavelength, zero magnitude flux density in W/m^2/micron, 
zero magnitude flux density in Jansky).  These are assumed to be given at the 
top of the file in comment lines starting with #.  The number_of_filters value 
tells the routine how many filters are listed, as this varies between the 
models with full wavelength coverage (BOSZ, Kurucz, and blackbody models) and 
the phoenix grid models that extend only to 5 microns and hence which are not 
suitable for use with MIRI filters and some other mid-infrared filters.  For 
that case, the number of filters is reduced to remove the long wavelength 
filters from the list.

Parameters:

   input_file_name  :  string

        The input ascii file name.  This is either a full path name or the 
        name of a file in the current directory.

   number_of_filters  :  integer

        The number of filter magnitudes listed in the file.  One expects
        that the value matches either the length of the kurucz_filter_names 
        list or the length of the phoenix_filter_names list.

Returns:

   model_magnitude_values  :  numpy floating point array, or None

      A two-dimensional array of model magnitude values with nominal size 
      `number_of_filters' times the number of lines of magnitudes in the
      file.  The latter value is not known in advance.  Note that 
      there is no check of the size of this array.

   model_magnitude_labels  :  list of strings, or None

       The list of names of the filters.

   filter_parameters  :  numpy floating point array, or None

       A two-demensional array of size `number_of_filters' x 3 with the filter
       effective wavelength and zero magnitude flux density values as read from 
       the magnitude list file.

If there is any error in trying to read the file the three values are None.

        """
        if (number_of_filters != len(phoenix_filter_names)) & (number_of_filters != len(kurucz_filter_names)):
            return None,None,None
        try:
            infile = open(input_file_name,'r')
            model_magnitude_labels = []
            filter_parameters = numpy.zeros((number_of_filters,3),dtype=numpy.float32)
            line = infile.readline()
            line = infile.readline()
            line = infile.readline()
            line = infile.readline()
            line = infile.readline()
            for i in range(number_of_filters):
                line = infile.readline()
                values = line.split()
                m = len(values)
                str = ''
                for n in range(2,m-3):
                    str = str+values[n]+' '
                    str = str.replace('filter','')
                    str = str.replace('Filter','')
                model_magnitude_labels.append(str)
                filter_parameters[i,0] = float(values[-3])
                filter_parameters[i,1] = float(values[-2])
                filter_parameters[i,2] = float(values[-1])
            infile.close()
            model_magnitude_values = numpy.loadtxt(input_file_name)
            return model_magnitude_values,model_magnitude_labels,filter_parameters
        except:
            return None,None,None

    def plot_colour_colour(self):
        """
This subroutine produces a stand-alone plot window and then calls the 
colour-colour plot routine to make the plot.  If the window exists, it just 
beings it back to visibility and renews the plot.

The plotting is a separate routine because one may want to change the range 
and re-plot, so in that event the code calls the plot subroutine again.

Parameters:

    There are no parameters.

Returns:

    No values are returned.

        """
        if self.plot_frame is None:
            self.plot_frame = Tk.Toplevel()
            self.plot_position_label_text = Tk.StringVar()
            self.plot_position_label = Tk.Label(self.plot_frame,textvariable=self.plot_position_label_text)
            self.plot_position_label.pack(side=Tk.TOP)
            self.plot_position_label_text.set('Position: ')
            self.plot_figure = Figure(figsize=(5,5),dpi=100)
            self.subplot = self.plot_figure.add_subplot(1,1,1)
            self.plot_canvas = FigureCanvasTkAgg(self.plot_figure,master=self.plot_frame)
            self.plot_canvas.show()
            self.plot_canvas.get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=Tk.YES)
            self.plot_canvas.mpl_connect("motion_notify_event",self.set_plot_position)
            range_frame = Tk.Frame(self.plot_frame)
            range_frame.pack(side=Tk.TOP)
            lab1 = Tk.Label(range_frame,text='x min')
            lab1.grid(row=0,column=0)
            self.range_entry[0] = Tk.Entry(range_frame,width=10)
            self.range_entry[0].grid(row=0,column=1)
            lab1 = Tk.Label(range_frame,text='y min')
            lab1.grid(row=1,column=0)
            self.range_entry[2]=Tk.Entry(range_frame,width=10)
            self.range_entry[2].grid(row=1,column=1)
            lab1 = Tk.Label(range_frame,text='x max')
            lab1.grid(row=0,column=2)
            self.range_entry[1] = Tk.Entry(range_frame,width=10)
            self.range_entry[1].grid(row=0,column=3)
            lab1=Tk.Label(range_frame,text='y max')
            lab1.grid(row=1,column=2)
            self.range_entry[3] = Tk.Entry(range_frame,width=10)
            self.range_entry[3].grid(row=1,column=3)
            Tk.Button(range_frame,text="Apply Range",command=lambda: self.make_plot(True)).grid(row=0,column=4)
            Tk.Button(range_frame,text="Default Range",command=lambda: self.make_plot(False)).grid(row=1,column=4)
            button_frame=Tk.Frame(self.plot_frame)
            button_frame.pack(side=Tk.TOP)
            Tk.Button(button_frame,text="Save as PS",command=self.makePS).pack(side=Tk.LEFT)
            Tk.Button(button_frame,text="Save as PNG",command=self.makePNG).pack(side=Tk.LEFT)
            Tk.Button(button_frame,text="Write Plot Values",command=self.on_print).pack(side=Tk.LEFT)
            Tk.Button(self.plot_frame,text="Close Window",command=self.plot_frame.withdraw).pack()
        else:
            self.plot_frame.deiconify()
            self.subplot.clear()
        try:
            self.range_entry[0].delete(0,Tk.END)
            self.range_entry[1].delete(0,Tk.END)
            self.range_entry[2].delete(0,Tk.END)
            self.range_entry[3].delete(0,Tk.END)
            self.range_entry[0].insert(0,str(self.data_limits[0,ind]))
            self.range_entry[1].insert(0,str(self.data_limits[1,ind]))
            self.range_entry[2].insert(0,str(self.data_limits[2,ind]))
            self.range_entry[3].insert(0,str(self.data_limits[3,ind]))
        except:
            pass
        self.make_plot(False)

    def make_plot(self,rangeflag):
        """
This is the main plot routine for the colour-colour window.  It will plot the 
simulated magnitudes and the fits.

Parameters:

    rangeflag  :  Boolean

        If equal to True, read the range values from the entry fields and 
        apply them.  If equal to False, get the values from the plot and 
        write them to the entry fields for the range.

Returns:

    There are no return values from this routine.

        """
        self.subplot.clear()
        pointsize = self.symbol_size.get()
        fit_flag = self.fit_flag.get()
        phoenix_inds = []
        kurucz_inds = []
        names = []
        flag = True
        # The following loop gets the selected magnitude index values and sees 
        # if they are in range.  One needs to check for indexes that produce
        # zero colours as well.  If flag is set to False, the routine will
        # print an error message and return.
        for loop in range(4):
            kurucz_inds.append(self.mag_boxes[loop].current())
            input_name = self.mag_boxes[loop].get()
            input_name = input_name.rstrip(' ')
            names.append(input_name)
            if kurucz_inds[-1] == -1:
                flag = False
            else:
                try:
                    indvalue = phoenix_filter_names.index(names[-1])
                except:
                    indvalue = -1
                phoenix_inds.append(indvalue)
        if (kurucz_inds[0] == kurucz_inds[1]) | (kurucz_inds[2] == kurucz_inds[3]):
            flag = False
        if not flag:
            self.put_message('Error: the requested filters are not set properly.]nPleae check your inputs.')
            return
        xlabel = kurucz_filter_names[kurucz_inds[0]]+' - '+kurucz_filter_names[kurucz_inds[1]]
        ylabel = kurucz_filter_names[kurucz_inds[2]]+' - '+kurucz_filter_names[kurucz_inds[3]]
        xlabel = xlabel.replace('filter','')
        xlabel = xlabel.replace('Filter','')
        ylabel = ylabel.replace('filter','')
        ylabel = ylabel.replace('Filter','')
        colours = ['black','blue','green','red']
        fit_colours = ['#708090','cyan','#00FF00','orange']
        self.data_limits = numpy.zeros((4),dtype=numpy.float32)
        self.plot_limits = numpy.zeros((4),dtype=numpy.float32)
        fit_order = 0
        modelfits = []
        if fit_flag == 1:
            try:
                fit_order = int(self.fit_order.get())
                if (fit_order < 1) | (fit_order > 20):
                    self.put_message('Error in the fit order (range 1 to 20).\nNo fits will be plotted.')
            except:
                fit_order = 0
                self.put_message('Error in the fit order (range 1 to 20).\nNo fits will be plotted.')
        for loop in range(4):
            xdata = None
            ydata = None
            value = self.model_flags[loop].get()
            if loop != 2:
                inds = kurucz_inds
            else:
                inds = phoenix_inds
            if value == 1:
                if min(inds) < 0:
                    modelfits.append(None)
                else:
                    xdata = self.magnitude_values[loop][:,inds[0]]-self.magnitude_values[loop][:,inds[1]]
                    ydata = self.magnitude_values[loop][:,inds[2]]-self.magnitude_values[loop][:,inds[3]]
                    plotrange = self.get_range(xdata,ydata)
                    if loop == 0:
                        self.data_limits = numpy.copy(plotrange[:,0])
                        self.plot_limits = numpy.copy(plotrange[:,1])
                    else:
                        self.data_limits[0] = min(self.data_limits[0],plotrange[0,0])
                        self.data_limits[1] = max(self.data_limits[1],plotrange[1,0])
                        self.data_limits[2] = min(self.data_limits[2],plotrange[2,0])
                        self.data_limits[3] = max(self.data_limits[3],plotrange[3,0])
                        self.plot_limits[0] = min(self.plot_limits[0],plotrange[0,1])
                        self.plot_limits[1] = max(self.plot_limits[1],plotrange[1,1])
                        self.plot_limits[2] = min(self.plot_limits[2],plotrange[2,1])
                        self.plot_limits[3] = max(self.plot_limits[3],plotrange[3,1])
                    self.subplot.plot(xdata,ydata,'o',color=colours[loop],ms=pointsize)
                    if fit_order > 0:
                        fit_results = legendre.legfit(xdata,ydata,fit_order)
                        modelfits.append(fit_results)
                    else:
                        modelfits.append(None)
            else:
                modelfits.append(None)
        delx = (self.plot_limits[1]-self.plot_limits[0])/10.
        if delx > 0.5:
            delx = 0.5
        model_xmin = self.plot_limits[0]-delx
        model_xmax = self.plot_limits[1]+delx
        delx = (model_xmax-model_xmin)/1000
        xmodel = numpy.arange(model_xmin,model_xmax,delx)
        for loop in range(4):
            if fit_order > 0:
                if modelfits[loop] is None:
                    pass
                else:
                    ymodel = legendre.legval(xmodel,modelfits[loop])
                    self.subplot.plot(xmodel,ymodel,color=fit_colours[loop])
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        if rangeflag:
            try:
                xmin = float(self.range_entry[0].get())
                xmax = float(self.range_entry[1].get())
                ymin = float(self.range_entry[2].get())
                ymax = float(self.range_entry[3].get())
                self.subplot.set_xlim(xmin,xmax)
                self.subplot.set_ylim(ymin,ymax)
            except:
                pass
        else:
            self.subplot.set_xlim(self.plot_limits[0],self.plot_limits[1])
            self.subplot.set_ylim(self.plot_limits[2],self.plot_limits[3])
            self.range_entry[0].delete(0,Tk.END)
            self.range_entry[1].delete(0,Tk.END)
            self.range_entry[2].delete(0,Tk.END)
            self.range_entry[3].delete(0,Tk.END)
            self.range_entry[0].insert(0,str(self.plot_limits[0]))
            self.range_entry[1].insert(0,str(self.plot_limits[1]))
            self.range_entry[2].insert(0,str(self.plot_limits[2]))
            self.range_entry[3].insert(0,str(self.plot_limits[3]))
        self.plot_canvas.show()

    def get_range(self,xdata,ydata):
        """
The routine takes as input (x,y) values assumed to be for a plot and returns 
the range values as a numpy array of length 4x2.  These are the raw values 
and then rounded values.

Parameters:
    
    xdata  :  numpy array

        The set of x values for a plot.

    ydata  :  numpy array

        The set of y values for a plot.

Returns:

    plotrange  :  numpy array, length 4,2

        plotrange [:,0] -> The values raw values xmin, xmax, ymin, ymax
        plotrange [:,1] -> The values rounded values xmin, xmax, ymin, ymax

        """
        plotrange = numpy.zeros((4,2),dtype=numpy.float32)
        plotrange[0,0] = numpy.min(xdata)
        plotrange[1,0] = numpy.max(xdata)
        plotrange[2,0] = numpy.min(ydata)
        plotrange[3,0] = numpy.max(ydata)
        plotrange[0,1] = self.round_float(plotrange[0,0],True)
        plotrange[1,1] = self.round_float(plotrange[1,0],False)
        plotrange[2,1] = self.round_float(plotrange[2,0],True)
        plotrange[3,1] = self.round_float(plotrange[3,0],False)
        return plotrange
        
    def on_print(self):
        pass
            
    def makePS(self):
        """
This saves the current plot as a positscript output file.  There are no 
parameters and no return values.
        """
        outfilename = tkFileDialog.asksaveasfilename(filetypes=[('postscript','*.ps')])
        if isinstance(outfilename,basestring) and outfilename != '':
            self.plot_figure.savefig(outfile,format="PS")

    def makePNG(self):
        """
This saves the current plot as a PNG output file.  There are not parameters and
no return values.
        """
        outfilename = tkFileDialog.asksaveasfilename(filetypes=[('PNG','*.png')])
        if isinstance(outfilename,basestring) and outfilename != '':
            self.plot_figure.savefig(outfile,format="PNG")

    def set_plot_position(self,event):
        """
This is call-back routine for motion in the colour-colour plot window.  When
an event comes in it writes the (x,y) position to the plot label area.

Parameters:

    event  :  Tk event variable

Returns:

    No values are returned from the routine.

        """
        try:
            string1="Position: %.6g %.6g" % (event.xdata,event.ydata)
            self.plot_position_label_text.set(string1)
        except:
            pass
    
            
def main(argv):
    """
This is the main program for the code when called from the command line.  It 
currently assumes that the code is to be run interactively. 

Parameters:

    argv  :  list of strings (sys.argv by assumption)

        This is not used currently.  It is present in case the code needs to 
        be able to be used non-interactively in the future.

Returns:

    No values are returned by the routine.

    """
    root = Tk.Tk()
    root.title("Colour-Colour Transformation Tool")
    x = magnitudeTransform()
    x.run_gui(root)
    root.mainloop()

if __name__ == "__main__":
    main(sys.argv)

