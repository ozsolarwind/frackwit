# frackwit

FrackWIT
========

(C) Copyright 2022 Angus King

FrackWIT is a program to generate a theoretical grid of fracking wellheads. It's aim is to
help visualise the impact of such a grid on the landscape. This is a very early prototype.
The program uses and creates a KML (Google Earth) file.

The program accepts a number of inputs, set out below and generates a KML file containing
a grid of wellheads (pads) and connecting tracks. The output file is stored in the same
folder as the input file and the filename is based upon the input filename with either a
prefix and/or suffix added to the filename. The inputs are:

o  Area of Interest: A KML (Google Earth) file with a LineString defining the area to be
   examined. This is the input file and is required. To change the file click the 'Change'
   button
   
o  Named Area: Choose a particular named area (Placemark) from a file with many areas 
   (sets of coordinates). If left empty the first set of coordinates is chosen
   
o  Initial roads/track: A KML file with a number of LineStrings defining existing
   tracks/roads that the fracking grid can connect to. This field is optional but helps to
   provide some variability to the network of roads/tracks produced. To change the file
   click the 'Change' button.
   
o  Track width: The width of the tracks in metres.

o  Track Colour: The colour the generated tracks are to appear in. To change the colour click
   the colour field.
   
o  Output file prefix: A prefix for the output filename

o  Output file suffix. A suffix for the output filename

o  Pad area: The size of the wellpad in thousands of square metres. The resulting pad is
   square in shape
   
o  Pad spacing: The east-west and north-south spacing between pads in metres

o  Pad Colour: The colour the generated pads are to appear in. To change the colour click
   the colour field
   
o  No of wells: The number of wells to be generated in hundreds

o  Wells per Pad: How many wells are accomodated for each pad.

Action buttons are:
o  Quit. Exit the program

o  Create grid. Create the output grid as a KML file. The program will start in the centre
   of the chosen area and radiate outwards placing pads where they fall within the area. A
   track will be generated to connect from the pad to the nearest (existing) track.
   
o  Help. Show this help
