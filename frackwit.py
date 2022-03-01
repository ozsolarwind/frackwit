#!/usr/bin/python3
#
#  Copyright (C) 2022 Angus King
#
#  frackwit.py - This file is (part of) frackwit.
#
#  frackwit is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  frackwit is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General
#  Public License along with catalogue.  If not, see
#  <http://www.gnu.org/licenses/>.
#

import configparser  # decode .ini file
import io
from math import sin, cos, radians, asin, acos, atan2, sqrt, degrees, pi
import os
import sys
import tkinter as tk
from tkinter import colorchooser
from tkinter import filedialog
from tkinter import ttk
from xml.etree.ElementTree import ElementTree, fromstring
import zipfile

RADIUS = 6367.   # km is the radius of the Earth

class Line:
    def __init__(self, name, coordinates):
        self.name = name
        self.coordinates = []
        for i in range(len(coordinates)):
            self.coordinates.append([])
            for j in range(len(coordinates[i])):
                self.coordinates[-1].append(round(coordinates[i][j], 6))

def dust(pyd, pxd, y1d, x1d, y2d, x2d):   # debug
    if y1d == y2d and x1d == x2d:
        return [-1]
    px = radians(pxd)
    py = radians(pyd)
    x1 = radians(x1d)
    y1 = radians(y1d)
    x2 = radians(x2d)
    y2 = radians(y2d)
    p_x = x2 - x1
    p_y = y2 - y1
    something = p_x * p_x + p_y * p_y
    u = ((px - x1) * p_x + (py - y1) * p_y) / float(something)
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    x = x1 + u * p_x
    y = y1 + u * p_y
    dx = x - px
    dy = y - py
    dist = sqrt(dx * dx + dy * dy)
    return [round(abs(dist) * RADIUS, 3), round(degrees(y), 6), round(degrees(x), 6)]

def gridConnect(lat, lon, lines, ignore=[], dummy_fix=False):
    shortest = [99999, -1., -1., -1]
    for l in range(len(lines)):
        if l in ignore:
            continue
        for i in range(len(lines[l].coordinates) - 1):
            if dummy_fix:
                dist = dust(lat, lon, lines[l].coordinates[i][0], lines[l].coordinates[i][1],
                       lines[l].coordinates[i + 1][0], lines[l].coordinates[i + 1][1])
            else:
                dist = DistancePointLine(lat, lon, lines[l].coordinates[i][0], lines[l].coordinates[i][1],
                       lines[l].coordinates[i + 1][0], lines[l].coordinates[i + 1][1])
            if dist[0] >= 0 and dist[0] < shortest[0]:
                shortest = dist[:]
                shortest.append(l)
    if shortest[0] == 99999:
         shortest[0] = -1
    return shortest   # length, lat, lon, line#

def DistancePointLine(pyd, pxd, y1d, x1d, y2d, x2d):
# px,py is the point to test.
# x1,y1,x2,y2 is the line to check distance.
# Returns distance from the line, or if the intersecting point on the line nearest
# the point tested is outside the endpoints of the line, the distance to the
# nearest endpoint.
# Returns -1 on 0 denominator conditions.
    px = radians(pxd)
    py = radians(pyd)
    x1 = radians(x1d)
    y1 = radians(y1d)
    x2 = radians(x2d)
    y2 = radians(y2d)
    b13 = Bearing(y1, x1, py, px)
    b12 = Bearing(y1, x1, y2, x2)
    d13 = Distance(y1, x1, py, px)
    d23 = Distance(y2, x2, py, px)
    dxt = asin(sin(d13) * sin(b13 - b12))
    dat = acos(cos(d13) / cos(dxt))
    iy = asin(sin(y1) * cos(dat) + cos(y1) * sin(dat) * cos(b12))
    ix = x1 + atan2(sin(b12) * sin(dat) * cos(y1), cos(dat) - sin(y1) * sin(iy))
    if abs(ix - x1) > abs(x1 - x2) or abs(iy - y1) > abs(y1 - y2):
        dst = d13
        ix = x1
        iy = y1
    else:
        dst = Distance(iy, ix, py, px)
        if d13 < dst:   # must be another way but this'll do for now
            dst = d13
            ix = x1
            iy = y1
    if d23 < dst:   # must be another way but this'll do for now
        dst = d23
        ix = x2
        iy = y2
    return [round(abs(dst) * RADIUS, 3), round(degrees(iy), 6), round(degrees(ix), 6)]

def Bearing(y1, x1, y2, x2):
# find the bearing between the coordinates
    return atan2(sin(x2 - x1) * cos(y2), cos(y1) * sin(y2) - sin(y1) * cos(y2) * cos(x2 - x1))

def destinationxy(lat1, lon1, bearing, distance):
    """
    Given a start point, initial bearing, and distance, calculate
    the destination point and final bearing travelling along a
    (shortest distance) great circle arc
    """
# convert decimal degrees to radians
    ln1, lt1, baring = list(map(radians, [lon1, lat1, bearing]))
# "reverse" haversine formula
    lat2 = asin(sin(lt1) * cos(distance / RADIUS) + cos(lt1) * sin(distance / RADIUS) * cos(baring))
    lon2 = ln1 + atan2(sin(baring) * sin(distance / RADIUS) * cos(lt1), cos(distance / RADIUS) - sin(lt1) * sin(lat2))
    return degrees(lat2), degrees(lon2)

def Distance(y1, x1, y2, x2):
# find the differences between the coordinates
    dy = y2 - y1
    dx = x2 - x1
    ra13 = pow(sin(dy / 2.), 2) + cos(y1) * cos(y2) * pow(sin(dx / 2.), 2)
    return 2 * asin(min(1, sqrt(ra13)))

def actualDistance(coord1, coord2, full=False):
    y1d = coord1[0]
    x1d = coord1[1]
    y2d = coord2[0]
    x2d = coord2[1]
    x1 = radians(x1d)
    y1 = radians(y1d)
    x2 = radians(x2d)
    y2 = radians(y2d)
    dst = Distance(y1, x1, y2, x2)
    if full:
        return abs(dst) * RADIUS
    else:
        return round(abs(dst) * RADIUS, 3)

def area_of_polygon(coords):
    """Calculates the area of an arbitrary polygon given its verticies"""
    def reproject(latitude, longitude):
        """Returns the x & y coordinates in meters using a sinusoidal projection"""
        earth_radius = RADIUS # in meters
        lat_dist = pi * earth_radius / 180.0
        y = [lat * lat_dist for lat in latitude]
        x = [int * lat_dist * cos(radians(lat)) for lat, int in zip(latitude, longitude)]
        return x, y

    latitude = []
    longitude = []
    for coord in coords:
        latitude.append(coord[0])
        longitude.append(coord[1])
    x, y = reproject(latitude, longitude)
    area = 0.0
    for i in range(-1, len(x)-1):
        area += x[i] * (y[i+1] - y[i-1])
    return abs(area) / 2.0


def within_map(x, y, poly):
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

def get_order(cells):
    midl = cells[0] // 2 * cells[1] + cells[1] // 2
    cr = 0
    cc = 0
    arry = []
    ctr = 0
    for y in range(cells[0]):
        arry.append([])
        for x in range(cells[1]):
            arry[-1].append(ctr)
            if ctr == midl:
                cr = y
                cc = x
            ctr += 1
    done = False
    lp = 1
    ctr = 1
    order = []
    order.append(arry[cr][cc])
    for c in range(cells[0] * cells[1] - 1):
        for col in range(lp):
            cc += 1 # east
            if cc >= cells[1] or cr < 0 or cc < 0:
                continue
            try:
                ctr += 1
                if ctr > (cells[0] * cells[1]):
                    done = True
                    break
                order.append(arry[cr][cc])
            except:
                ctr -= 1
        if done:
            break
        for col in range(lp):
            cr += 1 # north
            if cr >= cells[0] or cr < 0 or cc < 0:
                continue
            try:
                ctr += 1
                if ctr > (cells[0] * cells[1]):
                    done = True
                    break
                order.append(arry[cr][cc])
            except:
                ctr -= 1
        if done:
            break
        for col in range(lp + 1):
            cc -= 1 # west
            if cc < 0 or cr < 0:
                continue
            try:
                ctr += 1
                if ctr > (cells[0] * cells[1]):
                    done = True
                    break
                order.append(arry[cr][cc])
            except:
                ctr -= 1
        if done:
            break
        for col in range(lp + 1):
            cr -= 1 # south
            if cr < 0 or cc < 0:
                continue
            try:
                ctr += 1
                if ctr > (cells[0] * cells[1]):
                    done = True
                    break
                order.append(arry[cr][cc])
            except:
                ctr -= 1
        if done:
            break
        lp += 2
    return order

def getRoads(kml_file):
    zipped = False
    lines = []
    if kml_file[-4:] == '.kmz': # zipped file?
        zipped = True
        zf = zipfile.ZipFile(kml_file, 'r')
        inner_file = ''
        for name in zf.namelist():
            if name[-4:] == '.kml':
                inner_file = name
                break
        if inner_file == '':
            return
        memory_file = io.BytesIO()
        memory_file.write(zf.open(inner_file).read())
        root = ElementTree(fromstring(memory_file.getvalue()))
    else:
        kml_data = open(kml_file, 'rb')
        root = ElementTree(fromstring(kml_data.read()))
     # Create an iterator
    if sys.version_info[1] < 9: # before python 3.9
        iterat = root.getiterator()
    else:
        iterat = root.iter()
    placemark_id = ''
    line_names = []
    for element in iterat:
        elem = element.tag[element.tag.find('}') + 1:]
        if elem == 'name':
            line_name = element.text
            if placemark_id != '':
                line_name += placemark_id
                placemark_id = ''
        elif elem == 'Placemark':
            for key, value in list(element.items()):
                if key == 'id':
                    if value[:4] == 'kml_':
                        placemark_id = value[3:]
                    else:
                        placemark_id = value
        elif elem == 'coordinates':
            coords = []
            coordinates = ' '.join(element.text.split()).split()
            for i in range(len(coordinates)):
                coords.append([float(coordinates[i].split(',')[1]), float(coordinates[i].split(',')[0])])
            if line_name in line_names:
                i = 2
                while line_name + '#' + str(i) in line_names:
                    i += 1
                line_name += '#' + str(i)
            line_names.append(line_name)
            lines.append(Line(line_name, coords))
    if zipped:
        memory_file.close()
        zf.close()
    else:
        kml_data.close()
    return lines

def quit():
    sys.exit()

def get_kml_file():
    global kml_var, kml_widget
    file = kml_var.get()
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        olddir = file[:file.rfind('/')]
    else:
        olddir = file[:file.rfind('/')]
    # Open and return file path
    newfile = filedialog.askopenfilename(title='Select Area File', filetypes=(('KML', '*.kml'), ('KMZ', '*.kmz')),
                                         initialdir=olddir)
    if newfile != '':
        kml_var.set(newfile)
        if len(newfile) < 50:
            kml_widget.config(anchor=tk.W)
        else:
            kml_widget.config(anchor=tk.E)

def get_trk_file():
    global trk_var
    file = trk_var.get()
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        olddir = file[:file.rfind('/')]
    else:
        olddir = file[:file.rfind('/')]
    newfile = filedialog.askopenfilename(title='Select Tracks File', filetypes=(('KML', '*.kml'), ('KMZ', '*.kmz')),
                                         initialdir=olddir)
    trk_var.set(newfile)

def pick_pad_colour():
    global b4, pad_colour
    color = colorchooser.askcolor(pad_colr.get(), title ='Choose Well Pads Colour')
    b4.config(bg=color[1])
    pad_colr.set(color[1].upper())

def pick_trk_colour():
    global b3, trk_colour
    color = colorchooser.askcolor(trk_colr.get(), title ='Choose Tracks Colour')
    b3.config(bg=color[1])
    trk_colr.set(color[1].upper())

def show_help():
    root = tk.Toplevel(window)
    root.title('FrackWIT - Help')
    container = ttk.Frame(root)
    canvas = tk.Canvas(container, width=600)
    scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
    scrollable_frame = ttk.Frame(canvas)

    scrollable_frame.bind(
       "<Configure>",
       lambda e: canvas.configure(
           scrollregion=canvas.bbox("all")
       )
    )
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)
    help_lines = ['FrackWIT',
              "FrackWIT is a program to generate a theoretical grid of fracking wellheads. It's aim is to",
              'help visualise the impact of such a grid on the landscape.',
              'The program accepts a number of inputs, set out below and generates a KML file containing',
              'a grid of wellheads (pads) and connecting tracks. The output file is stored in the same',
              '   folder as the input file and the filename is based upon the input filename with either',
              '   a prefix and/or suffix added to the filename. The inputs are:',
              "o Area of Interest: A KML (Google Earth) file with a LineString defining the area to be",
              "   examined. This is the input file and is required. To change the file click the 'Change'",
              '   button',
              "o Initial roads/track: A KML file with a number of LineStrings defining existing",
              "   tracks/roads that the fracking grid can connect to. This field is optional but helps to",
              "   provide some variability to the network of roads/tracks produced. To change the file",
              "   click the 'Change' button.",
              'o Track width: The width of the tracks in metres.',
              'o Track Colour: The colour the generated tracks are to appear in. To change the colour click',
              '    the colour field.',
              'o Output file prefix: A prefix for the output filename',
              'o Output file suffix. A suffix for the output filename',
              'o Pad area: The size of the wellpad in thousands of square metres. The resulting pad is',
              '   square in shape',
              'o Pad spacing: The east-west and north-south spacing between pads in metres',
              'o Pad Colour: The colour the generated pads are to appear in. To change the colour click',
              '   the colour field',
              'o No of wells: The number of wells to be generated in hundreds',
              'o Wells per Pad: How many wells are accomodated for each pad.',
              'Action buttons are:',
              'o Quit. Exit the program',
              'o Create grid. Create the output grid as a KML file. The program will start in the centre',
              '   of the chosen area and radiate outwards placing pads where they fall within the area. A',
              '   track will be generated to connect from the pad to the nearest (existing) track.',
              'o Help. Show this help']
    tk.ttk.Label(scrollable_frame, text=help_lines[0], font=("Helvetica", 14, "bold"), wraplength=100, justify='left').pack(anchor='w')
    for line in help_lines[1:]:
        ttk.Label(scrollable_frame, text=line).pack(anchor='w')
    container.pack()
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

def create_grid():
    def add_kml_line(ctr, cell, grid_point, point):
        lline.append('\t<Placemark>')
        lline.append('\t\t<name>' + str(cell) + '</name><styleUrl>#line</styleUrl>')
        lline.append('\t\t<LineString>\n\t\t\t<coordinates>\n\t\t\t\t' + str(round(grid_point[2], 6)) + ',' + str(round(grid_point[1], 6)) + \
                     ',0 ' + str(round(point[1], 6)) + ',' + str(round(point[0], 6)) + ',0\n\t\t\t</coordinates>\n\t\t</LineString>')
        lline.append('\t\t<Polygon>\n\t\t\t<outerBoundaryIs><LinearRing>\n\t\t\t<coordinates>')
        dirn = degrees(Bearing(grid_point[1], grid_point[2], point[0], point[1]))
        road_len = actualDistance(grid_point[1 : 3], point, full=True)
        e = destinationxy(grid_point[1], grid_point[2], dirn + 90, minor_road_width_in_m / 2000.)
        n = destinationxy(grid_point[1], grid_point[2], dirn - 90, minor_road_width_in_m / 2000.)
        s = destinationxy(point[0], point[1], dirn + 90, minor_road_width_in_m / 2000.)
        w = destinationxy(point[0], point[1], dirn - 90, minor_road_width_in_m / 2000.)
        coords = '\t\t\t\t%0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 ' % (n[1], n[0], e[1], e[0],
                 s[1], s[0], w[1], w[0], n[1], n[0])
        lline.append(coords)
        lline.append('\t\t\t</coordinates>\n\t\t\t</LinearRing></outerBoundaryIs>\n\t\t</Polygon>')
        lline.append('\t</Placemark>')
        lines.append(Line(str(cell),[[point[0], point[1]], [grid_point[1], grid_point[2]]]))
        return road_len

    kml_file = kml_var.get()
    if not os.path.exists(kml_file):
        msg_widget.config(fg='red')
        msg_var.set(kml_file + ' not found')
        return
    tgt_name = tgt_var.get()
    trk_file = trk_var.get()
    if trk_file != '':
        if not os.path.exists(trk_file):
            msg_widget.config(fg='red')
            msg_var.set(trk_file + ' not found')
            return
    out_prefix = prefix.get()
    if len(out_prefix) > 0:
        if out_prefix[-1] != '_':
            out_prefix += '_'
    out_suffix = suffix.get()
    if len(out_suffix) > 0:
        if out_suffix[0] != '_':
            out_suffix = '_' + out_suffix
    if out_prefix == '' and out_suffix == '':
        out_prefix = 'fracking_grid_'
    pad_area_in_sqm = tk_area.get() * 1000
    pad_spacing_in_m = [tk_space1.get(), tk_space2.get()]
    no_of_wells = tk_wells.get() * 100
    wells_per_pad = tk_perpad.get()
    clr = pad_colr.get()
    colour2 = 'ff' + clr[5:] + clr[3:5] + clr[1:3]
    clr = trk_colr.get()
    colour = 'ff' + clr[5:] + clr[3:5] + clr[1:3]
    minor_road_width_in_m = tk_wdth.get()
    size = sqrt(pad_area_in_sqm) / 1000
    area_top_left = [None, None]
    area_bot_right = [None, None]
    cells = [0, 0]
    area_coords = None
    zipped = False
    if kml_file[-4:] == '.kmz': # zipped file?
        zipped = True
        zf = zipfile.ZipFile(kml_file, 'r')
        inner_file = ''
        for name in zf.namelist():
            if name[-4:] == '.kml':
                inner_file = name
                break
        if inner_file == '':
            return
        memory_file = io.BytesIO()
        memory_file.write(zf.open(inner_file).read())
        root = ElementTree(fromstring(memory_file.getvalue()))
    else:
        kml_data = open(kml_file, 'rb')
        root = ElementTree(fromstring(kml_data.read()))
     # Create an iterator
    if sys.version_info[1] < 9: # before python 3.9
        iterat = root.getiterator()
    else:
        iterat = root.iter()
    do_coords = False
    for element in iterat:
        elem = element.tag[element.tag.find('}') + 1:]
        if tgt_name != '':
            if elem == 'name':
                name = element.text.strip()
            elif elem == 'Placemark':
                try:
                    name = element[id]
                except:
                    pass
        if elem == 'Point':
            do_coords = False
        elif elem == 'LineString' or elem == 'LinearRing':
            do_coords = True
        elif elem == 'coordinates':
            if not do_coords:
                continue
            if tgt_name != '' and name != tgt_name:
                continue
            area_coords = []
            coordinates = ' '.join(element.text.split()).split()
            for i in range(len(coordinates)):
                area_coords.append([round(float(coordinates[i].split(',')[1]), 6),
                  round(float(coordinates[i].split(',')[0]), 6)])
            break
    if zipped:
        memory_file.close()
        zf.close()
    else:
        kml_data.close()
    if area_coords is None:
        msg_widget.config(fg='red')
        msg = tgt_name + ' coordinates not found in ' + kml_file
        msg = msg.strip()
        msg = msg[:60] + '\n' + msg[60:120] + '\n' + msg[120:]
        msg_var.set(msg)
        return
    for coord in area_coords:
        for i in range(2):
            if area_top_left[i] is None:
                area_top_left[i] = coord[i]
            elif coord[i] > area_top_left[i]:
                area_top_left[i] = coord[i]
            if area_bot_right[i] is None:
                area_bot_right[i] = coord[i]
            elif coord[i] < area_bot_right[i]:
                area_bot_right[i] = coord[i]
    ew_distance = actualDistance(area_top_left, [area_top_left[0], area_bot_right[1]])
    ns_distance = actualDistance(area_top_left, [area_bot_right[0], area_top_left[1]])
    cells[1] = int(ew_distance / (pad_spacing_in_m[1] / 1000.))
    cells[0] = int(ns_distance / (pad_spacing_in_m[0] / 1000.))
    max_no_of_pads = cells[0] * cells[1]
    if max_no_of_pads == 0:
        msg_widget.config(fg='red')
        msg_var.set('Chosen area too small')
        return
    if no_of_wells == 0:
        no_of_wells = max_no_of_pads * wells_per_pad
    else:
        no_of_wells = min(max_no_of_pads * wells_per_pad, no_of_wells)
  #  msg = 'potential pads = {:,d} (need {:,d} for {:,d} wells)'.format(max_no_of_pads, int(no_of_wells / wells_per_pad), no_of_wells)
   # msg_var.set(msg)
    centre = [area_top_left[0] + (area_bot_right[0] - area_top_left[0]) / 2.,
              area_top_left[1] + (area_bot_right[1] - area_top_left[1]) / 2]
    tlla = destinationxy(centre[0], centre[1], 0., (pad_spacing_in_m[0] * cells[0] / 2) / 1000.)
    tllo = destinationxy(centre[0], centre[1], 270., (pad_spacing_in_m[1] * cells[1] / 2) / 1000.)
    top_left = [tlla[0], tllo[1]]
    trlo = destinationxy(centre[0], centre[1], 90., (pad_spacing_in_m[1] * cells[1] / 2) / 1000.)
    top_right = [tlla[0], trlo[1]]
    blla = destinationxy(centre[0], centre[1], 180., (pad_spacing_in_m[0] * cells[0] / 2) / 1000.)
    bot_left = [blla[0], tllo[1]]
    bot_right = [blla[0], trlo[1]]
    ew_step = destinationxy(centre[0], centre[1], 90., pad_spacing_in_m[1] / 1000.)
    ns_step = destinationxy(centre[0], centre[1], 0., pad_spacing_in_m[0] / 1000.)
    if trk_file == '':
        lines = []
        point1 = None
    else:
        lines = getRoads(trk_file)
    left_delta = []
    for i in range(2):
        left_delta.append((top_left[i] - bot_left[i]) / (cells[0] - 1))
    top_delta = []
    for i in range(2):
        top_delta.append((top_right[i] - top_left[i]) / (cells[1] - 1))
    points = []
    cells_2 = [int(cells[0] / 2), int(cells[1] / 2)]
    r1 = cells_2[0]
    c1 = cells_2[1]
    for row in range(cells[0]):
        lat1 = bot_left[0] + left_delta[0] * row
        lon1 = bot_left[1] + left_delta[1] * row
        for col in range(cells[1]):
            points.append([lat1 + top_delta[0] * col, lon1 + top_delta[1] * col])
    ctr = 0
    road_len = 0
    roads = 0
    order = get_order(cells)
    pline = ['<?xml version="1.0" encoding="UTF-8"?>',
                     '<kml xmlns="http://www.opengis.net/kml/2.2">',
                     '<Document>']
    lline = ['<Folder>\n<name>Roads</name>']
    pline.append('<name>Fracking Grid</name>')
    cmt_line = len(pline)
    pline.append('<!-- Inputs:\nArea of interest: {}\nPad area: {:,d}\nWells per Pad: {:d}\nOutput file prefix: {}\nPad spacing: {:d} x {:d}\nOutput file suffix: {}\nInitial roads/tracks: {}\nNo of wells: {:d}\nMinor road width: {:d}\n     Outputs:\n'.format(kml_file, pad_area_in_sqm, wells_per_pad, out_prefix.strip('_'), pad_spacing_in_m[0], pad_spacing_in_m[1], out_suffix.strip('_'), trk_file,
       no_of_wells, minor_road_width_in_m))
    pline.append('<Style id="well">\n\t<LineStyle>\n\t\t<color>' + colour + '</color>\n\t\t<width>1</width>\n\t</LineStyle>' \
                 + '\n\t<PolyStyle>\n\t\t<color>' + colour + '</color>\n\t</PolyStyle>\n</Style>')
    pline.append('<Style id="line">\n\t<LineStyle>\n\t\t<color>' + colour2 + '</color>\n\t\t<width>1</width>\n\t</LineStyle>' \
                 + '\n\t<PolyStyle>\n\t\t<color>' + colour2 + '</color>\n\t</PolyStyle>\n</Style>')
    descr_line = len(pline)
    pline.append('')
    pline.append('<Folder>')
    pline.append('<name>Well pads</name>')
    last_point = []
    tilt = -45.
    check_ctr = int(no_of_wells / wells_per_pad)
    for cell in order:
        point = points[cell]
        if not within_map(point[0], point[1], area_coords):
            continue
        ctr += 1
        if ctr > check_ctr:
            break
        pline.append('\t<Placemark>')
        pline.append('\t\t<name>%s</name>\n\t\t<styleUrl>#well</styleUrl>\n\t\t<Point>\n\t\t\t<coordinates>\n\t\t\t\t%0.6f,%0.6f,0\n\t\t\t</coordinates>\n\t\t</Point>' % (cell, point[1], point[0]))
        pline.append('\t\t<Polygon>\n\t\t\t<outerBoundaryIs><LinearRing>\n\t\t\t<coordinates>')
        e = destinationxy(point[0], point[1], 90. + tilt, size)
        w = destinationxy(point[0], point[1], 270. + tilt, size)
        n = destinationxy(point[0], point[1], 0. + tilt, size)
        s = destinationxy(point[0], point[1], 180. + tilt, size)
        coords = '\t\t\t\t%0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 %0.6f,%0.6f,0 ' % (n[1], n[0], e[1], e[0],
                 s[1], s[0], w[1], w[0], n[1], n[0])
        pline.append(coords)
        pline.append('\t\t\t</coordinates>\n\t\t\t</LinearRing></outerBoundaryIs>\n\t\t</Polygon>')
        pline.append('\t</Placemark>')
        if len(lines) == 0:
            if point1 is None:
                point1 = [point[:]]
                continue
            l1 = [pad_spacing_in_m[0] / 1000.] + point1[0]
            road_len += add_kml_line(ctr, 0, l1, point)
            roads += 1
            point1.append(point)
            lines.append(Line('1', point1))
        grid_point = gridConnect(point[0], point[1], lines)
        road_len += add_kml_line(ctr, cell, grid_point, point)
        roads += 1
        lines.append(Line(str(cell),[[point[0], point[1]], [grid_point[1], grid_point[2]]]))
        if len(last_point) < 0:
            if last_point[1] == point[1] or last_point[0] == point[0]:
                road_len += add_kml_line(ctr, cell, grid_point, point)
                roads += 1
        last_point[:] = point
        last_cell = cell
    road_len -= roads * size / 2. # reduce by part in the pad
    pline.append('</Folder>')
    pline[cmt_line] += 'Area of interest is: {:0,.2f} sq km\n'.format(area_of_polygon(area_coords))
    pline[cmt_line] += 'Total road length is: {:0,.1f} km. Road area is: {:,} sqm ({:0.2f} sq km)\n'.format(road_len,
                       int(round(road_len * minor_road_width_in_m / 1000., 3) * 1000000), road_len * minor_road_width_in_m / 1000.)
    pline[cmt_line] += 'Total wellpad area is: {:,} sqm ({:0.2f} sq km); {:,} wells\n -->'.format(ctr * pad_area_in_sqm, ctr * pad_area_in_sqm / 1000000., ctr * wells_per_pad)
    pline[descr_line] = '<description><![CDATA[This KML file shows a fracking grid for ' + kml_file[:-4].title() + \
                        '. It has '+  str(ctr) + ' wellpads.]]></description>'
    lline.append('</Folder>\n</Document>\n</kml>')
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        i = kml_file.rfind('/')
    else:
        i = kml_file.rfind('/')
    kml_out = kml_file[:i + 1] + out_prefix + kml_file.lower()[i + 1:-4]
    if tgt_name != '':
        if kml_out.find(tgt_name.lower()) < 0 and out_suffix.lower() != tgt_name.lower():
            kml_out += '_' + tgt_name.lower()
    kml_out += out_suffix + '.kml'
    k_file = open(kml_out, 'w')
    for line in pline:
        k_file.write(line + '\n')
    for line in lline:
        k_file.write(line + '\n')
    k_file.close()
#   save config
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        try:
            ini_file = '~\\frackwit.ini'.replace('~', os.environ['HOME'])
        except:
            ini_file = 'frackwit.ini'
    else:
        ini_file = '~/.user/frackwit.ini'.replace('~', os.path.expanduser('~'))
    if os.path.exists(ini_file):
        if os.path.exists(ini_file + '~'):
            os.remove(ini_file + '~')
        if os.path.exists(ini_file):
            os.rename(ini_file, ini_file + '~')
    sou = open(ini_file, 'w')
    lines = '[Variables]\nkml_file={}\ntgt_name={}\npad_area_in_sqm={:d}\nwells_per_pad={:d}\nout_prefix={}\npad_spacing_in_m={:d},{:d}\nout_suffix={}\ntrk_file={}\nno_of_wells={:d}\nminor_road_width_in_m={:d}\ntrk_colour={}\npad_colour={}'.format(kml_file, tgt_name, pad_area_in_sqm, wells_per_pad, out_prefix.strip('_'), pad_spacing_in_m[0], pad_spacing_in_m[1], out_suffix.strip('_'), trk_file,
       no_of_wells, minor_road_width_in_m, trk_colr.get(), pad_colr.get())
    sou.write(lines)
    sou.close()
    msg = str(ctr) + ' Wellpads saved to ' + kml_out
    msg = msg[:60] + '\n' + msg[60:120] + '\n' + msg[120:]
    msg_widget.config(fg='green')
    msg_var.set(msg)

file_dir = '' #'/home/' + os.getlogin() + '/Documents/Gus/Duplicated/SEN/Tech Stuff/KCER 2021/fracking/'
pad_spacing_in_m = [800, 800]
pad_area_in_sqm = 15000
wells_per_pad = 1
no_of_wells = 0
kml_file = 'area.kml'
tgt_name = ''
trk_file = 'tracks.kml'
out_prefix = 'fracking_grid'
out_suffix = ''
trk_colour = '#FFFF00'
pad_colour = '#FFFF00'
minor_road_width_in_m = 12
config = configparser.RawConfigParser()
if sys.platform == 'win32' or sys.platform == 'cygwin':
    try:
        config_file = '~\\frackwit.ini'.replace('~', os.environ['HOME'])
    except:
        config_file = 'frackwit.ini'
else:
    config_file = '~/.user/frackwit.ini'.replace('~', os.path.expanduser('~'))
if os.path.exists(config_file):
    config.read(config_file)
    items = config.items('Variables')
    for key, value in items:
        if key == 'kml_file':
             kml_file = value
        elif key == 'tgt_name':
             tgt_name = value
        elif key == 'trk_file':
             trk_file = value
        elif key == 'pad_spacing_in_m':
             try:
                 spc = value.split(',')
                 pad_spacing_in_m = [int(spc[0]), int(spc[1])]
             except:
                 pass
        elif key == 'pad_area_in_sqm':
             try:
                 pad_area_in_sqm = int(value) / 1000
             except:
                 pass
        elif key == 'wells_per_pad':
             try:
                 wells_per_pad = int(value)
             except:
                 pass
        elif key == 'out_prefix':
             out_prefix = value
        elif key == 'out_suffix':
             out_suffix = value
        elif key == 'trk_colour':
             trk_colour = value
        elif key == 'pad_colour':
             pad_colour = value
        elif key == 'minor_road_width_in_m':
             try:
                 minor_road_width_in_m = int(value)
             except:
                 pass
        elif key == 'no_of_wells':
             try:
                 no_of_wells = int(value) / 100
             except:
                 pass
window = tk.Tk()
window.geometry('640x520')
window.title('frackwit - Create a theoretical grid of Fracking wells')
row = 0
tk.Label(window, text='Area of interest:').grid(row=row, sticky=tk.W, padx=5, pady=5)
kml_var = tk.StringVar()
kml_var.set(file_dir + kml_file)
kml_widget = tk.Label(window, textvariable=kml_var, width=50, anchor=tk.E, bg='white', relief=tk.SUNKEN)
kml_widget.grid(row=row, column=1, columnspan=3, pady=5, padx=5, sticky=tk.W)
if len(kml_file) < 50:
    kml_widget.config(anchor=tk.W)
else:
    kml_widget.config(anchor=tk.E)
b1 = tk.Button(window, text='Change', command=get_kml_file)
b1.grid(row=row, column=3)
row += 1
tk.Label(window, text='(KML file with single LineString)').grid(row=row, column=1, sticky=tk.W, padx=5)
row += 1
tk.Label(window, text='Named Area:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tgt_var = tk.StringVar()
tgt_var.set(tgt_name)
tgt_widget = tk.Entry(window, textvariable=tgt_var)
tgt_widget.grid(row=row, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
tk.Label(window, text='(Case sensitive)').grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Initial roads/tracks:').grid(row=row, sticky=tk.W, padx=5, pady=5)
trk_var = tk.StringVar()
trk_var.set(trk_file)
trk_widget = tk.Label(window, textvariable=trk_var, width=50, anchor=tk.E, bg='white', relief=tk.SUNKEN)
trk_widget.grid(row=row, column=1, columnspan=3, pady=5, padx=5, sticky=tk.W)
if len(trk_file) < 50:
    trk_widget.config(anchor=tk.W)
else:
    trk_widget.config(anchor=tk.E)
b2 = tk.Button(window, text='Change', command=get_trk_file)
b2.grid(row=row, column=3)
row += 1
tk.Label(window, text='(KML file with multiple LineStrings)').grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Track width:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tk_wdth = tk.IntVar()
tk_wdth.set(minor_road_width_in_m)
trk_wdth_widget = tk.Spinbox(window, from_=4, to=50, width=10, textvariable=tk_wdth)
trk_wdth_widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
tk.Label(window, text='(m)').grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Track Colour:').grid(row=row, sticky=tk.W, padx=5, pady=5)
trk_colr = tk.StringVar()
trk_colr.set(trk_colour)
b3 = tk.Button(window, textvariable=trk_colr, command=pick_trk_colour, padx=5, pady=5, bg=trk_colour)
b3.grid(row=row, column=1, sticky=tk.W)
tk.Label(window, text='(Click colour to change)').grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Output file prefix:').grid(row=row, sticky=tk.W, padx=5, pady=5)
prefix = tk.StringVar()
prefix.set(out_prefix)
out_prefix_widget = tk.Entry(window, textvariable=prefix)
out_prefix_widget.grid(row=row, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Output file suffix:').grid(row=row, sticky=tk.W, padx=5, pady=5)
suffix = tk.StringVar()
suffix.set(out_suffix)
out_suffix_widget = tk.Entry(window, textvariable=suffix)
out_suffix_widget.grid(row=row, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Pad area:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tk_area = tk.IntVar()
tk_area.set(pad_area_in_sqm)
pad_area_widget = tk.Spinbox(window, from_=1, to=50, width=10, textvariable=tk_area)
pad_area_widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
tk.Label(window, text="('000 sqm)").grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Pad spacing:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tk_space1 = tk.IntVar()
tk_space1.set(pad_spacing_in_m[0])
pad_space1_widget = tk.Spinbox(window, from_=500, to=5000, width=10, textvariable=tk_space1)
pad_space1_widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
tk_space2 = tk.IntVar()
tk_space2.set(pad_spacing_in_m[1])
pad_space2_widget = tk.Spinbox(window, from_=500, to=5000, width=10, textvariable=tk_space2)
pad_space2_widget.grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
tk.Label(window, text="(m x m)").grid(row=row, column=3, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Pad Colour:').grid(row=row, sticky=tk.W, padx=5, pady=5)
pad_colr = tk.StringVar()
pad_colr.set(pad_colour)
b4 = tk.Button(window, textvariable=pad_colr, command=pick_pad_colour, padx=5, pady=5, bg=pad_colour)
b4.grid(row=row, column=1, sticky=tk.W)
tk.Label(window, text='(Click colour to change)').grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='No of wells:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tk_wells = tk.IntVar()
tk_wells.set(no_of_wells)
wells_widget = tk.Spinbox(window, from_=1, to=1000, width=10, textvariable=tk_wells)
wells_widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
tk.Label(window, text="(100's)").grid(row=row, column=2, sticky=tk.W, padx=5, pady=5)
row += 1
tk.Label(window, text='Wells per Pad:').grid(row=row, sticky=tk.W, padx=5, pady=5)
tk_perpad = tk.IntVar()
tk_perpad.set(wells_per_pad)
wells_per_widget = tk.Spinbox(window, from_=1, to=50, width=10, textvariable=tk_perpad)
wells_per_widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
row += 1
q = tk.Button(window, text='Quit', command=quit)
q.grid(row=row, column=0, padx=5, pady=5)
b = tk.Button(window, text='Create grid', command=create_grid)
b.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
h = tk.Button(window, text='Help', command=show_help)
h.grid(row=row, column=3, sticky=tk.W, padx=5, pady=5)
row += 1
msg_var = tk.StringVar()
msg_var.set('')
msg_widget = tk.Label(window, textvariable=msg_var, width=60, anchor='w', fg='green',
                      justify=tk.LEFT)
msg_widget.grid(row=row, rowspan=3, column=1, columnspan=3, pady=5, padx=5, sticky=tk.NW)
tk.mainloop()
sys.exit()
