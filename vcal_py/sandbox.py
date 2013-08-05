import site
site.addsitedir('')
import VC
import Utils
import Converter
import os
import cPickle
import quakelib
import h5py

#import cProfile
import math

import numpy as np
#import Chaos
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import matplotlib.lines as mlines

#import eqsim
import subprocess
#import glob
#from PIL import Image
#import shutil
import profile


print Utils.magnitude(1.17, 7.05e8)

''' Print the total slip of each event 
data_file = '/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10.h5'
sys_name = os.path.basename(data_file).split('.')[0]
vc_sys = VC.VCSys(sys_name, data_file)

start_event = vc_sys.eventForYear(10000.0)
end_event = vc_sys.eventForYear(40000.0)

f = h5py.File(vc_sys.data_file, 'r')
events = f['event_table']
sweeps = f['event_sweep_table']

for event in events[start_event:end_event]:
    slip = 0.0
    ruptured_elements = {}
    for sweep in sweeps[event[8]:event[9]]:
        slip += sweep[3]
        try:
            tmp = ruptured_elements[sweep[2]]
        except KeyError:
            ruptured_elements[sweep[2]] = True

    if len(ruptured_elements) > 10:
        print event[0], slip, slip/float(len(ruptured_elements))

f.close()
'''

''' Print the surface areas of all the sections 
data_file = '/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10.h5'
sys_name = os.path.basename(data_file).split('.')[0]
vc_sys = VC.VCSys(sys_name, data_file)

for sid, section in vc_sys.sections.iteritems():
    print sid, section.sname, section.surface_area
'''

''' Given an event find the x year interval around it 
data_file = '/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10.h5'
sys_name = os.path.basename(data_file).split('.')[0]
vc_sys = VC.VCSys(sys_name, data_file)

target_evid = 487847
interval = 500.0

event_info = vc_sys.webgui_getEventDetail(target_evid)

start_year = event_info['year'] - interval/2.0
end_year = event_info['year'] + interval/2.0

print start_year, end_year
'''



''' create the webgui caches 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_webgui/media/vc_data/ALLCAL2_no-creep_dt-08_st-10.h5'
sys_name = os.path.basename(data_file).split('.')[0]
vc_sys = VC.VCSys(sys_name, data_file)

#print vc_sys.eventForYear(10000)

vc_sys.webgui_checkCaches()
'''

''' webgui profiling 

def doThing():
    data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_webgui/media/vc_data/ALLCAL2_no-creep_dt-08_st-10.h5'
    sys_name = os.path.basename(data_file).split('.')[0]
    vc_sys = VC.VCSys(sys_name, data_file)
    print vc_sys.webgui_getEventList(0,1000)

profile.run('doThing()')
'''

''' examine properties of outputted EQSim files
# parse the eqsim file
f = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_110912-122869_Events_slip-map-5.5.dat','r')
eqsim_events = {}
eqsim_slip_maps = {}
for line in f:
    if line.startswith('200'):
        linedat = line.split(' ')
        
        event_id        = int(linedat[1])
        magnitude       = float(linedat[2])
        time            = float(linedat[3])
        duration        = float(linedat[4])
        sid             = int(linedat[5])
        depth_lo        = float(linedat[6])
        depth_hi        = float(linedat[7])
        das_lo          = float(linedat[8])
        das_hi          = float(linedat[9])
        hypo_depth      = float(linedat[10])
        hypo_das        = float(linedat[11])
        area            = float(linedat[12])
        mean_slip       = float(linedat[13])
        moment          = float(linedat[14])
        shear_before    = float(linedat[15])
        shear_after     = float(linedat[16])
        normal_before   = float(linedat[17])
        normal_after    = float(linedat[18])
        
        event = [event_id, magnitude, time, duration, sid, depth_lo, depth_hi, das_lo, das_hi, hypo_depth, hypo_das, area, mean_slip, moment, shear_before, shear_after, normal_before, normal_after]
        
        try:
            eqsim_events[event_id].append(event)
        except KeyError:
            eqsim_events[event_id] = [event]
    
    elif line.startswith('201'):
        linedat = line.split(' ')
        
        depth_lo        = float(linedat[1])
        depth_hi        = float(linedat[2])
        das_lo          = float(linedat[3])
        das_hi          = float(linedat[4])
        area            = float(linedat[5])
        mean_slip       = float(linedat[6])
        moment          = float(linedat[7])
        shear_before    = float(linedat[8])
        shear_after     = float(linedat[9])
        normal_before   = float(linedat[10])
        normal_after    = float(linedat[11])
        element_id      = int(linedat[12])
        
        ele_slip = [depth_lo, depth_hi, das_lo, das_hi, area, mean_slip, moment, shear_before, shear_after, normal_before, normal_after, element_id]
        
        try:
            eqsim_slip_maps[event_id][sid].append(ele_slip)
        except KeyError:
            try:
                eqsim_slip_maps[event_id][sid] = [ele_slip]
            except KeyError:
                eqsim_slip_maps[event_id] = {sid: [ele_slip]}

# get the element id mapping file
eid_remap_file = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_element-remap.dat','r')
eid_remap = {}
for ln, line in enumerate(eid_remap_file):
    if ln != 0:
        orig_id, vc_id = map(int, line.split())
        eid_remap[vc_id] = orig_id

# parse the hdf5 file
fh5 = h5py.File('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10.h5', 'r')
h5_events = fh5['event_table']
h5_sweeps = fh5['event_sweep_table']

#event_h5_mags = []
#calc_h5_mags = []

#print 'h5_eid', 'eqsim_eid', 'calc_h5_mag', 'calc_eqsim_mag', 'event_h5_mag', 'event_eqsim_mag', 'num_of_sections', 'ave_slips_per_block'
for eqsim_eid in range(1,len(eqsim_events)):
    # do analysis
    #eqsim_eid = 10
    #print eqsim_eid
    # find the h5 slips and areas for the given event
    the_eid = 110912 - 1 + eqsim_eid
    h5_event = h5_events[the_eid]
    h5_event_sweeps = h5_sweeps[h5_event[8]:h5_event[9]]
    h5_areas = {}
    h5_slips = {}
    h5_num_of_slips = {}
    for sweep in h5_event_sweeps:
        #print sweep[2], sweep[4]
        h5_areas[sweep[2]] = float(sweep[4])
        try:
            h5_num_of_slips[sweep[2]] += 1
            h5_slips[sweep[2]] += float(sweep[3])
        except KeyError:
            h5_num_of_slips[sweep[2]] = 1
            h5_slips[sweep[2]] = float(sweep[3])

    try:
        # find the eqsim slips and areas for the given event
        eqsim_areas = {}
        eqsim_slips = {}
        eqsim_event = eqsim_events[eqsim_eid]
        eqsim_event_slipmap = eqsim_slip_maps[eqsim_eid]
        for sid in eqsim_event_slipmap.keys():
            for ele in eqsim_event_slipmap[sid]:
                eqsim_areas[ele[11]] = ele[4]
                eqsim_slips[ele[11]] = ele[5]
    except KeyError:
        eqsim_areas = None
        eqsim_slips = None


    lame_mu = 3.000000e+010
    #print the_eid, eqsim_eid
    h5_moment = 0.0
    eqsim_moment = 0.0
    for eid in h5_areas.keys():
        h5_moment += h5_areas[eid] * h5_slips[eid] * lame_mu
        if eqsim_areas is None and eqsim_slips is None:
            eqsim_moment = None
        else:
            eqsim_moment += eqsim_areas[eid_remap[eid]] * eqsim_slips[eid_remap[eid]] * lame_mu
        #print eid, eid_remap[eid], h5_areas[eid], eqsim_areas[eid_remap[eid]], h5_slips[eid], eqsim_slips[eid_remap[eid]]
    #event_h5_mags.append(h5_event[3])
    #calc_h5_mags.append((2.0/3.0)*math.log(1.0e7 * h5_moment,10) - 10.7)


    print the_eid, eqsim_eid, (2.0/3.0)*math.log(1.0e7 * h5_moment,10) - 10.7, (2.0/3.0)*math.log(1.0e7 * eqsim_moment,10) - 10.7 if eqsim_moment is not None else eqsim_event[0][1], h5_event[3], eqsim_event[0][1], len(eqsim_events[eqsim_eid]), np.mean(h5_num_of_slips.values())
'''

''' output scaling relations with mags from the sweeps instead of the reported event mags 
# parse the hdf5 file
fh5 = h5py.File('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/ALLCAL2_no-creep_dt-08_st-10.h5', 'r')
h5_events = fh5['event_table']
h5_sweeps = fh5['event_sweep_table']
h5_blocks = fh5['block_info_table']

start_event = 110912
end_event = 471207

event_h5_mags = []
calc_h5_mags = []
h5_rupture_areas = []
h5_average_slips = []
#h5_rupture_length = []

lame_mu = 3.000000e+010

number_of_years = h5_events[end_event][1] - h5_events[start_event][1]

print number_of_years

# get values from data
for h5_event in h5_events[start_event:end_event]:
    h5_event_sweeps = h5_sweeps[h5_event[8]:h5_event[9]]
    h5_areas = {}
    h5_slips = {}
    h5_surface_lengths = {}
    for sweep in h5_event_sweeps:
        #h5_surface_lengths[sweep[2]] = math.sqrt(float(sweep[4])) if int(h5_blocks[sweep[2]][7]) != 0 else 0.0
        h5_areas[sweep[2]] = float(sweep[4])
        try:
            h5_slips[sweep[2]] += float(sweep[3])
        except KeyError:
            h5_slips[sweep[2]] = float(sweep[3])

    h5_moment = 0.0
    h5_rupture_area = sum(h5_areas.values())
    h5_average_slip = np.mean(h5_slips.values())
    #h5_rupture_length = sum(h5_surface_lengths.values())

    for eid in h5_areas.keys():
        h5_moment += h5_areas[eid] * h5_slips[eid] * lame_mu

    event_h5_mags.append(h5_event[3])
    calc_h5_mags.append((2.0/3.0)*math.log(1.0e7 * h5_moment,10) - 10.7)

    h5_rupture_areas.append(h5_rupture_area)
    h5_average_slips.append(h5_average_slip)
    #h5_rupture_length.append(h5_rupture_length)

total_events = len(event_h5_mags)
print total_events

# frequency - magnitude
f1 = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/fm_event.dat','w')
f2 = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/fm_calc.dat','w')

event_h5_cum_freq = {}
calc_h5_cum_freq = {}

for num, magnitude in enumerate(sorted(event_h5_mags)):
    event_h5_cum_freq['%.10f'%magnitude] = total_events - (num + 1)

for num, magnitude in enumerate(sorted(calc_h5_mags)):
    calc_h5_cum_freq['%.10f'%magnitude] = total_events - (num + 1)

f1.write('event_mag event_freq\n')

for magnitude in sorted(event_h5_cum_freq.iterkeys()):
    f1.write('%f %f\n'%(float(magnitude), float(event_h5_cum_freq[magnitude])/number_of_years))

f2.write('calc_mag calc_freq\n')

for magnitude in sorted(calc_h5_cum_freq.iterkeys()):
    f2.write('%f %f\n'%(float(magnitude), float(calc_h5_cum_freq[magnitude])/number_of_years))

f1.close()
f2.close()

# magnitude - rupture area
f1 = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/mra.dat','w')

f1.write('rupture_area event_mag calc_mag\n')

for e, rupture_area in enumerate(h5_rupture_areas):
    f1.write('%f %f %f\n'%(rupture_area, event_h5_mags[e], calc_h5_mags[e]))

f1.close()

# magnitude - average slip
f1 = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/mas.dat','w')

f1.write('average_slip event_mag calc_mag\n')

for e, average_slip in enumerate(h5_average_slips):
    f1.write('%f %f %f\n'%(average_slip, event_h5_mags[e], calc_h5_mags[e]))

f1.close()

'''

''' examine properties of outputted EQSim files 

f = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/ALLCAL2_no-creep_dt-08_st-10_110912-112005_Events_slip-map-5.5.dat','r')

lame_mu = 3.000000e+010
events = []
slip_maps = {}

for line in f:
    if line.startswith('200'):
        linedat = line.split(' ')
        
        event_id        = int(linedat[1])
        magnitude       = float(linedat[2])
        time            = float(linedat[3])
        duration        = float(linedat[4])
        sid             = int(linedat[5])
        depth_lo        = float(linedat[6])
        depth_hi        = float(linedat[7])
        das_lo          = float(linedat[8])
        das_hi          = float(linedat[9])
        hypo_depth      = float(linedat[10])
        hypo_das        = float(linedat[11])
        area            = float(linedat[12])
        mean_slip       = float(linedat[13])
        moment          = float(linedat[14])
        shear_before    = float(linedat[15])
        shear_after     = float(linedat[16])
        normal_before   = float(linedat[17])
        normal_after    = float(linedat[18])
    
        events.append([event_id,magnitude,time,duration,sid,depth_lo,depth_hi,das_lo,das_hi,hypo_depth,hypo_das,area,mean_slip,moment,shear_before,shear_after,normal_before,normal_after])
    
    elif line.startswith('201'):
        linedat = line.split(' ')
        
        depth_lo        = float(linedat[1])
        depth_hi        = float(linedat[2])
        das_lo          = float(linedat[3])
        das_hi          = float(linedat[4])
        area            = float(linedat[5])
        mean_slip       = float(linedat[6])
        moment          = float(linedat[7])
        shear_before    = float(linedat[8])
        shear_after     = float(linedat[9])
        normal_before   = float(linedat[10])
        normal_after    = float(linedat[11])
        element_id      = int(linedat[12])
        
        slip_map = [depth_lo,depth_hi,das_lo,das_hi,area,mean_slip,moment,shear_before,shear_after,normal_before,normal_after,element_id]
        
        try:
            slip_maps['%i %i'%(event_id, sid)].append(slip_map)
        except KeyError:
            slip_maps['%i %i'%(event_id, sid)] = [slip_map]

previous_event_id = None
for event in events:
    try:
        if previous_event_id == event[0]:
            event_slip_map += slip_maps['%i %i'%(event[0], event[4])]
        else:
            event_slip_map = slip_maps['%i %i'%(event[0], event[4])]
            event_slips = []
            event_areas = []
    except KeyError:
        event_slip_map = None
    
    if event_slip_map is not None:
        for ele_slip_map in event_slip_map:
            event_slips.append(ele_slip_map[5])
            event_areas.append(ele_slip_map[4])
        moment = 0.0
        for i, slip in enumerate(event_slips):
            moment += lame_mu * slip * event_areas[i]
            
        print event[0], event[1], (2.0/3.0)*math.log(1.0e7 * moment,10) - 10.7, (2.0/3.0)*math.log(1.0e7 * np.mean(event_slips) * sum(event_areas) * lame_mu,10) - 10.7, event[12], np.mean(event_slips), event[11], sum(event_areas)
    previous_event_id == event[0]
'''

''' look at properties of the parsed vcsys 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2/ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

for sid in vc_sys.geometry.layered_sections.iterkeys():
    print sid, vc_sys.geometry.layered_sections[sid]
'''

''' create new models with a specified number of elements from the ALLCAL2 model 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2/ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

max_elements = 10000
num_elements = 0
sections_to_delete = []

for sid, section in vc_sys.geometry.sections.iteritems():
    num_elements += len(section.selement_ids)

    if num_elements > max_elements:
        sections_to_delete.append(sid)
        

for sid in sections_to_delete:
    vc_sys.deleteSection(sid)

vc_sys.name = vc_sys.name + '_Trunc-%i'%max_elements

vc_sys.plotFaultMap()
vc_sys.exportKMLGeometry(vc_sys.name+'.kml')
vc_sys.exportEqsimGeometry(vc_sys.name+'_Geometry.dat', vc_sys.name+'_Friction.dat')
'''

''' look at slips based on element area and stress drop 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run4/ALLCAL2_no-creep_dt-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

for eid, ele in vc_sys.geometry.elements.iteritems():
    print eid, (ele.area()**0.5 * ele.static_strength) * (1.0/ele.lame_mu), ele.static_strength, ele.area()**0.5, ele.lame_mu, vc_sys.geometry.sections[ele.sid].sname
'''

''' find out which section an element is in 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run3/ALLCAL2_no-creep_dt-05_st-00_long.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

print vc_sys.geometry.elements[4870].sid
'''

''' plot greens functions

# plot parameters
imw = 1024.0 # the full image width
lm = 40.0
rm = 50.0
tm = 40.0
bm = 40.0
res = 72.0
cbh = 20.0
cbs = 40.0
vss = 50.0 # vertical section spacing

arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
        
# data files
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/new_vcal/build/Debug/B-Lagoon-et-al.h5'
greens_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/new_vcal/build/Debug/B-Lagoon-et-al_Greens.h5'


os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

f = h5py.File(greens_file, 'r')

greens_normal = f['greens_normal']
greens_shear = f['greens_shear']

greens_normal_np = np.array(greens_normal)
greens_shear_np = np.array(greens_shear)

# create colormaps so that zero values map to white
shear_white_pt = -greens_shear_np.min()/(greens_shear_np.max() - greens_shear_np.min())
shear_color_dict = {'red':      (   (0.0, 0.0, 1.0),
                                    (shear_white_pt, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
    
                    'green':    (   (0.0, 0.0, 0.0),
                                    (shear_white_pt, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
    
                    'blue':     (   (0.0, 0.0, 0.0),
                                    (shear_white_pt, 1.0, 1.0),
                                    (1.0, 1.0, 0.0))
                    }

normal_white_pt = -greens_normal_np.min()/(greens_normal_np.max() - greens_normal_np.min())
normal_color_dict = {'red':      (   (0.0, 0.0, 1.0),
                                    (normal_white_pt, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
    
                    'green':    (   (0.0, 0.0, 0.0),
                                    (normal_white_pt, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
    
                    'blue':     (   (0.0, 0.0, 0.0),
                                    (normal_white_pt, 1.0, 1.0),
                                    (1.0, 1.0, 0.0))
                    }
    

shear_cmap = LinearSegmentedColormap('shear_cmap', shear_color_dict, N=256, gamma=1.0)
normal_cmap = LinearSegmentedColormap('normal_cmap', normal_color_dict, N=256, gamma=1.0)

imh = imw/2.0 + cbh + cbs
imwi = imw/res
imhi = imh/res
fig = plt.figure(figsize=(imwi, imhi), dpi=res)
ph = imh - tm - bm - cbh - cbs # the height for both matricies
pw = ph
shear_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
normal_ax = fig.add_axes(((imw - pw - rm)/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))

#patch_size =  pw/float(greens_shear_np.shape[0])

#maxWeight = 2**np.ceil(np.log(np.abs(greens_shear_np).max())/np.log(2))


#print greens_shear_np + abs(greens_shear_np.min())

#, norm=LogNorm(vmin=greens_shear_np.min(), vmax=greens_shear_np.max())
shear_ax.pcolor(greens_shear_np, cmap=shear_cmap )
normal_ax.pcolor(greens_normal_np, cmap=normal_cmap  )

shear_ax.invert_yaxis()
normal_ax.invert_yaxis()

shear_ax.axis('tight')
normal_ax.axis('tight')

for label in shear_ax.xaxis.get_ticklabels() + shear_ax.yaxis.get_ticklabels() + normal_ax.xaxis.get_ticklabels() + normal_ax.yaxis.get_ticklabels():
    label.set_fontproperties(arial12)

for tick in shear_ax.xaxis.get_major_ticks() + normal_ax.xaxis.get_major_ticks():
    tick.label1On = False
    tick.label2On = True

shear_cb_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, cbh/imh))
normal_cb_ax = fig.add_axes(((imw - pw - rm)/imw, bm/imh, pw/imw, cbh/imh))

norm = mpl.colors.Normalize(vmin=greens_shear_np.min(), vmax=greens_shear_np.max())
cb = mpl.colorbar.ColorbarBase(shear_cb_ax, cmap=shear_cmap, norm=norm, orientation='horizontal')

norm = mpl.colors.Normalize(vmin=greens_normal_np.min(), vmax=greens_normal_np.max())
cb = mpl.colorbar.ColorbarBase(normal_cb_ax, cmap=normal_cmap, norm=norm, orientation='horizontal')

for label in shear_cb_ax.xaxis.get_ticklabels() + normal_cb_ax.xaxis.get_ticklabels():
    label.set_fontproperties(arial10)
for line in shear_cb_ax.xaxis.get_ticklines() + normal_cb_ax.xaxis.get_ticklines():
    line.set_alpha(0)

#patches = []
#for (x,y),w in np.ndenumerate(greens_shear_np):
    #print x,y,w
    #if w > 0: color = 'white'
    #else:     color = 'black'
    #size = np.sqrt(np.abs(w))
    #ssmr_patches.append(mpatches.Rectangle(((hi*ppe)/pw, vsl - (vi*ppe)/ph), ppe/pw, ppe/ph, fc=element_color))
    #fig_ax.add_patch(rect)
#fig_ax.autoscale_view()

# Reverse the yaxis limits
#fig_ax.set_ylim(*fig_ax.get_ylim()[::-1])

# mark what rows correspond to which sections
for sid, section in vc_sys.geometry.sections.iteritems():
    min_id = min(section.selement_ids)
    max_id = max(section.selement_ids)
    
    shear_ax.add_line(mlines.Line2D((0,greens_shear_np.shape[0] + 1), (max_id + 1,max_id + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
    shear_ax.add_line(mlines.Line2D((max_id + 1,max_id + 1), (0,greens_shear_np.shape[0] + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
    
    shear_ax.text((max_id + min_id)/2, max_id, '%i %s'%(sid, section.sname), ha='center', va='bottom', fontproperties=arial7_light)

    normal_ax.add_line(mlines.Line2D((0,greens_shear_np.shape[0] + 1), (max_id + 1,max_id + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
    normal_ax.add_line(mlines.Line2D((max_id + 1,max_id + 1), (0,greens_shear_np.shape[0] + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))

    normal_ax.text((max_id + min_id)/2, max_id, '%i %s'%(sid, section.sname), ha='center', va='bottom', fontproperties=arial7_light)

print np.unravel_index(np.argmax(greens_shear_np), greens_shear_np.shape), np.unravel_index(np.argmin(greens_shear_np), greens_shear_np.shape)

for (x,y),w in np.ndenumerate(greens_shear_np):
    print x,y,w

fig.savefig('%s_Greens.png'%(sys_name), format='png')
'''

''' plot a matrix of the sections and the FM spectrum for each section 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/ALLCAL2_no-creep_dt-08_st-10.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

start_event = vc_sys.eventForYear(float(10000))
end_event = vc_sys.eventForYear(float(10010))

print start_event, end_event
f = h5py.File(data_file, 'r')
        
events = f['event_table']
sweeps = f['event_sweep_table']

mags = []

mags_in_section = {}

for event in events[start_event:end_event]:
    mag = event[3]
    mags.append(mag)
    involved_sections = {}
    event_sweeps = sweeps[event[8]:event[9]]

    for sweep in event_sweeps:
        eid = sweep[2]
        sid = vc_sys.geometry.elements[eid].sid
    
        involved_sections[sid] = sid

    for sid in involved_sections.keys():
        try:
            mags_in_section[sid].append(mag)
        except KeyError:
            mags_in_section[sid] = [mag]

print mags_in_section


    #event_sections.append({'id':event[0], 'year':event[1], 'mag':event[3], 'sections':sorted(involved_sections.keys())})


hist, mag_bins = np.histogram(mags, bins = 100)

section_vs_mag = np.zeros( (mag_bins.size, len(vc_sys.geometry.sections)) )

section_ids = np.array(sorted(vc_sys.geometry.sections.keys()))

section_id_index_mapping = {}
section_names = []

for i, sec_id in enumerate(section_ids):
    section_id_index_mapping[sec_id] = i
    section_names.append(vc_sys.geometry.sections[sec_id].sname)

print mag_bins
print mag_bins.size



for event in events[start_event:end_event]:
    involved_sections = {}
    event_sweeps = sweeps[event[8]:event[9]]
    
    for sweep in event_sweeps:
        eid = sweep[2]
        sid = vc_sys.geometry.elements[eid].sid
        involved_sections[sid] = sid

    event_mag = event[3]
    involved_section_ids = involved_sections.keys()

    mag_bin_num = np.digitize(np.array([event_mag]), mag_bins) - 1

    #print event[0], event_mag, mag_bin_num,
    for sec_id in involved_section_ids:
        #print vc_sys.geometry.sections[sec_id].sname,
        section_vs_mag[mag_bin_num,section_id_index_mapping[sec_id]] += 1
    #print


#print section_vs_mag

 # plot parameters
imw = 1280.0 # the full image width
lm = 60.0
rm = 20.0
tm = 20.0
bm = 30.0
res = 72.0
cbh = 10.0
cbs = 80.0
vss = 50.0 # vertical section spacing

arial7 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7)

cmap = cm.Purples

imh = imw * (float(len(mag_bins))/float(len(section_ids))) + cbh + cbs
imwi = imw/res
imhi = imh/res
fig = plt.figure(figsize=(imwi, imhi), dpi=res)
ph = imh - tm - bm - cbh - cbs # the height for the matrix
pw = imw - lm - rm
fig_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))

print section_vs_mag.min(), section_vs_mag.max()
norm = mpl.colors.LogNorm(vmin=1, vmax=section_vs_mag.max())

fig_ax.pcolor(section_vs_mag, cmap=cmap, norm=norm )

fig_ax.axis('tight')
fig_ax.set_xticks(range(0, len(section_names)))
for i in range(0, len(section_names)):
    fig_ax.axvline(x=float(i), lw=0.2, c='0.2')
fig_ax.set_xticklabels(section_names)

fig_ax.set_yticks(range(0, len(mag_bins)))
mag_bin_labels = []
for i in range(0, len(mag_bins)):
    fig_ax.axhline(y=float(i)-0.1, lw=0.2, c='0.2')
    mag_bin_labels.append('%f'%(mag_bins[i]))
fig_ax.set_yticklabels(mag_bin_labels)

for line in fig_ax.xaxis.get_ticklines() + fig_ax.yaxis.get_ticklines():
    line.set_alpha(0)

for label in fig_ax.xaxis.get_ticklabels():
    label.set_fontproperties(arial7)
    label.set_rotation(90)
    label.set_ha('left')
for label in fig_ax.yaxis.get_ticklabels():
    label.set_fontproperties(arial7)

cb_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, cbh/imh))
cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='horizontal')

for label in cb_ax.xaxis.get_ticklabels():
    label.set_fontproperties(arial7)
for line in cb_ax.xaxis.get_ticklines():
    line.set_alpha(0)

fig.savefig('%s_Mag_By_Section.png'%(vc_sys.name), format='png')

f.close()

'''

''' plot a matrix of the sections and the binned number of events in each section 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run4/ALLCAL2_no-creep_dt-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

start_event = vc_sys.eventForYear(float(10000))
end_event = vc_sys.eventForYear(float(40000))

print start_event, end_event
f = h5py.File(data_file, 'r')
        
events = f['event_table']
sweeps = f['event_sweep_table']

mags = []

for event in events[start_event:end_event]:
    mags.append(event[3])
    #involved_sections = {}
    #event_sweeps = sweeps[event[8]:event[9]]

    #for sweep in event_sweeps:
    #    eid = sweep[2]
    #    sid = vc_sys.geometry.elements[eid].sid
    
    #    involved_sections[sid] = sid
    
    #event_sections.append({'id':event[0], 'year':event[1], 'mag':event[3], 'sections':sorted(involved_sections.keys())})


hist, mag_bins = np.histogram(mags, bins = 100)

section_vs_mag = np.zeros( (mag_bins.size, len(vc_sys.geometry.sections)) )

section_ids = np.array(sorted(vc_sys.geometry.sections.keys()))

section_id_index_mapping = {}
section_names = []

for i, sec_id in enumerate(section_ids):
    section_id_index_mapping[sec_id] = i
    section_names.append(vc_sys.geometry.sections[sec_id].sname)

print mag_bins
print mag_bins.size



for event in events[start_event:end_event]:
    involved_sections = {}
    event_sweeps = sweeps[event[8]:event[9]]
    
    for sweep in event_sweeps:
        eid = sweep[2]
        sid = vc_sys.geometry.elements[eid].sid
        involved_sections[sid] = sid

    event_mag = event[3]
    involved_section_ids = involved_sections.keys()

    mag_bin_num = np.digitize(np.array([event_mag]), mag_bins) - 1

    #print event[0], event_mag, mag_bin_num,
    for sec_id in involved_section_ids:
        #print vc_sys.geometry.sections[sec_id].sname,
        section_vs_mag[mag_bin_num,section_id_index_mapping[sec_id]] += 1
    #print


#print section_vs_mag

 # plot parameters
imw = 1280.0 # the full image width
lm = 60.0
rm = 20.0
tm = 20.0
bm = 30.0
res = 72.0
cbh = 10.0
cbs = 80.0
vss = 50.0 # vertical section spacing

arial7 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7)

cmap = cm.Purples

imh = imw * (float(len(mag_bins))/float(len(section_ids))) + cbh + cbs
imwi = imw/res
imhi = imh/res
fig = plt.figure(figsize=(imwi, imhi), dpi=res)
ph = imh - tm - bm - cbh - cbs # the height for the matrix
pw = imw - lm - rm
fig_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))

print section_vs_mag.min(), section_vs_mag.max()
norm = mpl.colors.LogNorm(vmin=1, vmax=section_vs_mag.max())

fig_ax.pcolor(section_vs_mag, cmap=cmap, norm=norm )

fig_ax.axis('tight')
fig_ax.set_xticks(range(0, len(section_names)))
for i in range(0, len(section_names)):
    fig_ax.axvline(x=float(i), lw=0.2, c='0.2')
fig_ax.set_xticklabels(section_names)

fig_ax.set_yticks(range(0, len(mag_bins)))
mag_bin_labels = []
for i in range(0, len(mag_bins)):
    fig_ax.axhline(y=float(i)-0.1, lw=0.2, c='0.2')
    mag_bin_labels.append('%f'%(mag_bins[i]))
fig_ax.set_yticklabels(mag_bin_labels)

for line in fig_ax.xaxis.get_ticklines() + fig_ax.yaxis.get_ticklines():
    line.set_alpha(0)

for label in fig_ax.xaxis.get_ticklabels():
    label.set_fontproperties(arial7)
    label.set_rotation(90)
    label.set_ha('left')
for label in fig_ax.yaxis.get_ticklabels():
    label.set_fontproperties(arial7)

cb_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, cbh/imh))
cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='horizontal')

for label in cb_ax.xaxis.get_ticklabels():
    label.set_fontproperties(arial7)
for line in cb_ax.xaxis.get_ticklines():
    line.set_alpha(0)

fig.savefig('%s_Mag_By_Section.png'%(vc_sys.name), format='png')

f.close()
'''

''' custom colormap tests 
cdict1 = {  'red':  (   (0.0, 0.0, 1.0),
                        (0.3, 1.0, 1.0),
                        (1.0, 0.0, 0.0)),
    
            'green': (  (0.0, 0.0, 0.0),
                        (0.3, 1.0, 1.0),
                        (1.0, 0.0, 0.0)),
    
            'blue': (   (0.0, 0.0, 0.0),
                        (0.3, 1.0, 1.0),
                        (1.0, 1.0, 0.0))
        }
    


#blue_red1 = LinearSegmentedColormap.from_list('cmap', ['#86b8f5','#ffffff'], N=256, gamma=1.0)
blue_red1 = LinearSegmentedColormap('cmap', cdict1, N=256, gamma=1.0)

a = np.linspace(0, 1, 256).reshape(1,-1)
a = np.vstack((a,a))

fig = plt.figure(figsize=(5,10))

plt.axis("off")
plt.imshow(a, aspect='auto', cmap=blue_red1, origin='lower')

fig.savefig("disp-map-test/cmap.png", format='png')
'''

''' plot a series of event rupture maps 

os.chdir('/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_py/')
for eid in range(234137, 234214):
    proc_args = 'python vcal_py.py -d /Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_no-creep_mod-stress-drop/run1/ALLCAL2_no-creep_mod-stress-drop_dt-08_st-10.h5 --erm %i'%eid
    proc = subprocess.Popen(proc_args, shell=True)
    proc.wait()
'''

''' changing model defined stress drops 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep/run5/ALLCAL2_no-creep_dt-08_st-10.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

for eid, element in vc_sys.geometry.elements.iteritems():
    element.static_strength *= 0.9

new_name = 'ALLCAL2_no-creep_mod-stress-drop'
vc_sys.exportEqsimGeometry('%s_Geometry.dat'%new_name, '%s_Friction.dat'%new_name)
'''

''' changing model defined slip rates 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/10x4/ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))

for eid, element in vc_sys.geometry.elements.iteritems():
    element.slip *= 2.5

new_name = 'ALLCAL2_no-creep_mod-slip'
vc_sys.exportEqsimGeometry('%s_Geometry.dat'%new_name, '%s_Friction.dat'%new_name)
'''

''' create webgui test files for different list sizes and sorts 
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_webgui/media/vc_data/ALLCAL2_no-creep_dt-08_st-10.h5'

os.chdir(os.path.dirname(data_file))
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))


sizes = [500]
#sort_fields = { 'evid':             0,
#                'evyear':           1,
#                'evmag':            2,
#                'evtriggersec':     3,
#                'evtriggerele':     4,
#                'evinvolvedsec':    5,
#                'evinvolvedele':    6,
#                'evaveslip':        7
#                }

sort_fields = { 'evid':             0}
sort_directions = ['ascending']

for size in sizes:
    for fname, fid in sort_fields.iteritems():
        for direction in sort_directions:
            print '%i %s %s'%(size, fname, direction)
            event_list = vc_sys.webgui_getEventList(data_file, 0, size, fid, direction)
            #filename = '/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_webgui/design/templates/event_list_test_%i_dat_%s_%s.json'%(size, fname, direction)
            #f = open(filename,'w')
            #f.write(event_list)
            #f.close()
'''


''' find the lat lon of two points a certian distance apart 
converter = Converter.Converter()

converter.setLatLon0( 38.536693, -121.751596)

print converter.lat0, converter.lon0
print converter.xy_latlon(3000, 0)
'''

'''
plot okada stresses for a single element
c           = 0.0
dip         = (math.pi/180) * 90.0
L           = 1000.0
W           = 1000.0
US          = -10
UD          = 0.0
UT          = 0.0
lame_lambda = 3.2e10
lame_mu     = 3.0e10

plot_width = 1.993
plot_height = 1.993
resolution = 300.0
h_pix = int(plot_width * resolution)
v_pix = int(plot_height * resolution)

cmap = plt.get_cmap('rainbow')

action = 'top'
do_plot = True

if action is 'top':
    x_l = -5000
    x_r = 6000
    y_b = -5500
    y_t = 5500

    Xs = np.linspace(x_l,x_r,h_pix)
    Ys = np.linspace(y_b,y_t,v_pix)

    Xs,Ys = np.meshgrid(Xs,Ys)

    Ds = np.empty(Xs.shape)
    Ds_log = np.empty(Ds.shape)
    D_colors = np.empty((Ds.shape[0],Ds.shape[1],4))

    okada = eqsim.Okada()
    okada.reset()
    okada.precalc(dip, lame_lambda, lame_mu)

    it = np.nditer(Xs, flags=['multi_index'])
    while not it.finished:
        Ds[it.multi_index] = okada.get_exy(Xs[it.multi_index], Ys[it.multi_index], 0, c, L, W, US, UD, UT)
        it.iternext()
        
    it = np.nditer(Ds, flags=['multi_index'])
    while not it.finished:
        val = Ds[it.multi_index]
        if val <= 0:
            try:
                Ds_log[it.multi_index] = -math.log(-val,10)
            except:
                print val
        else:
            Ds_log[it.multi_index] =  math.log(val,10)
        it.iternext()
    
    if do_plot:
        it = np.nditer(Ds_log, flags=['multi_index'])
        while not it.finished:
            val = (Ds_log[it.multi_index] - Ds_log.min())/(Ds_log.max() - Ds_log.min())
            r,g,b,a = cmap(val)
            D_colors[it.multi_index[0], it.multi_index[1], 0] = r
            D_colors[it.multi_index[0], it.multi_index[1], 1] = g
            D_colors[it.multi_index[0], it.multi_index[1], 2] = b
            D_colors[it.multi_index[0], it.multi_index[1], 3] = a
            it.iternext()
        
        fig = plt.figure(figsize=[plot_width, plot_height],dpi=resolution)
        the_ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])

        the_ax.plot((0,L), (0,0), color='k', linewidth=1, solid_capstyle='round', solid_joinstyle='round')
        the_ax.arrow( 0, 0, L, 0, color='k', linewidth=1, shape='right', head_width=300, head_length=300, length_includes_head=True)
        the_ax.arrow( L, -0, -L, 0, color='k', linewidth=1, shape='right', head_width=300, head_length=300, length_includes_head=True)
        the_ax.imshow(D_colors, interpolation='none',extent=(x_l, x_r, y_b, y_t))

        for label in the_ax.xaxis.get_ticklabels() + the_ax.yaxis.get_ticklabels():
            label.set_alpha(0)
                
        for line in the_ax.xaxis.get_ticklines() + the_ax.yaxis.get_ticklines():
            line.set_alpha(0)

        plt.savefig('/Users/sachs/Dropbox/ORALSMK2/prospectus/graphics/okada_stress_top.pdf', format='pdf')
    
    else:
        print Ds.min() * 2.0 * lame_mu, Ds.max() * 2.0 *lame_mu
elif action is 'side':
    x_l = -5000
    x_r = 6000
    z_b = -10000
    z_t = 1000

    Xs = np.linspace(x_l,x_r,h_pix)
    Zs = np.linspace(z_t,z_b,v_pix)

    Xs,Zs = np.meshgrid(Xs,Zs)
    
    Ds = np.empty(Xs.shape)
    Ds_log = np.empty(Ds.shape)
    D_colors = np.empty((Ds.shape[0],Ds.shape[1],4))

    okada = eqsim.Okada()
    okada.reset()
    okada.precalc(dip, lame_lambda, lame_mu)

    it = np.nditer(Xs, flags=['multi_index'])
    while not it.finished:
        Ds[it.multi_index] = okada.get_exy(Xs[it.multi_index], 0, Zs[it.multi_index], c, L, W, US, UD, UT)
        it.iternext()

    
    it = np.nditer(Ds, flags=['multi_index'])
    while not it.finished:
        val = Ds[it.multi_index]
        if val <= 0:
            try:
                Ds_log[it.multi_index] = -math.log(-val,10)
            except:
                pass
        else:
            Ds_log[it.multi_index] =  math.log(val,10)
        it.iternext()
        
    if do_plot:
        it = np.nditer(Ds_log, flags=['multi_index'])
        while not it.finished:
            val = (Ds_log[it.multi_index] - Ds_log.min())/(Ds_log.max() - Ds_log.min())
            r,g,b,a = cmap(val)
            D_colors[it.multi_index[0], it.multi_index[1], 0] = r
            D_colors[it.multi_index[0], it.multi_index[1], 1] = g
            D_colors[it.multi_index[0], it.multi_index[1], 2] = b
            D_colors[it.multi_index[0], it.multi_index[1], 3] = a
            it.iternext()


        fig = plt.figure(figsize=[plot_width, plot_height],dpi=resolution)
        the_ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        
        f_ele = mpatches.Rectangle((0, -W), L, W, ec='none', fc='k')
        f_sky = mpatches.Rectangle((x_l, 0), abs(x_r) + abs(x_l), z_t, ec='k', fc='none', hatch='////')
        the_ax.add_patch(f_ele)
        the_ax.add_patch(f_sky)
        the_ax.plot((x_l,x_r), (0,0), color='k', linewidth=0.5, solid_capstyle='round', solid_joinstyle='round')
        the_ax.arrow( L, -W*0.5, L*0.5, 0, color='k', linewidth=1, shape='right', head_width=300, head_length=300, length_includes_head=True)
        the_ax.arrow( 0, -W*0.5, -L*0.5, 0, color='k', linewidth=1, shape='right', head_width=300, head_length=300, length_includes_head=True)
        the_ax.plot((L*0.5,0), (-W*0.5,-W*0.5), color='white', linewidth=1)
        the_ax.imshow(D_colors, interpolation='none',extent=(x_l, x_r, z_b, z_t))

        for label in the_ax.xaxis.get_ticklabels() + the_ax.yaxis.get_ticklabels():
            label.set_alpha(0)
                
        for line in the_ax.xaxis.get_ticklines() + the_ax.yaxis.get_ticklines():
            line.set_alpha(0)

        plt.savefig('/Users/sachs/Dropbox/ORALSMK2/prospectus/graphics/okada_stress_side.pdf', format='pdf')
        
    else:
        print Ds.min() * 2.0 * lame_mu, Ds.max() * 2.0 *lame_mu
elif action is 'cbar':
    a = np.linspace(0, 1, 256).reshape(1,-1)
    a = np.vstack((a,a))
    
    fig = plt.figure(figsize=[plot_width, .5],dpi=resolution)
    the_ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    plt.axis("off")
    the_ax.imshow(a, aspect='auto', cmap=cmap, origin='lower')
    
    plt.savefig('/Users/sachs/Dropbox/ORALSMK2/prospectus/graphics/okada_stress_cbar.pdf', format='pdf')
'''

'''
dir = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep_param-sweep/vc_05_20/'
sys_name = '/ALLCAL2_1-7-11_no-creep_dyn-05_st-20'

NE_pt = MyVector.LatLonPoint( 36.764708,-120.447200)
SW_pt = MyVector.LatLonPoint( 36.356748,-120.888881)

vc_sys = cPickle.load(open('%s%s_cache/%s.pkl'%(dir,sys_name,sys_name), 'rb'))

#vc_sys = VC.VCSys()
#vc_sys.loadHdf5Geometry(data_file)
print vc_sys.sectionsInArea(NE_pt, SW_pt)
'''

'''
dir = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep_param-sweep/vc_05_20/anim_test/'

for imagename in glob.glob('%s/*.png'%dir):
    #print os.path.dirname(imagename)
    im = Image.open(imagename)
    imr = im.resize((614,874))
    imr.save('%s/r_%s'%(os.path.dirname(imagename),os.path.basename(imagename)), "PNG")


print 'Creating movie'
proc_args = 'ffmpeg -y -r 15 -b:v 64k -i %sr_%s.png anim.mp4'%(dir, '%d')
proc = subprocess.Popen(proc_args, shell=True)
proc.wait()
'''

'''
if not os.path.exists(dir_new):
    os.makedirs(dir_new)

for imagename in glob.glob('%s/*.pdf'%dir_old):
    index = int(os.path.basename(imagename).split('.')[0])
    shutil.copyfile(imagename, '%s%i.pdf'%(dir_new,index+1))
    #im = Image.open(imagename)
    #im.save('%s%i.png'%(dir_new,index+1), "PNG")
'''

'''
#108722
total = len(range(108722,109895))
for num, evid in enumerate(range(108722,109895)):
#    print evid
    print '##################### %i of %i #####################'%(num+1,total)
    proc_args = 'python vcal_py.py -d /Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep_param-sweep/vc_05_20/ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5 --edm %i'%(evid)
    proc = subprocess.Popen(proc_args, shell=True)
    proc.wait()


'''

'''
This code creates a movie of simple deformation

okada = eqsim.Okada()

for n,d in enumerate(np.arange(-5,5.1,0.1)):

    X,Y = np.meshgrid(np.arange(-20,21,1),np.arange(-20,21,1) )

    U = np.empty(X.shape)
    V = np.empty(Y.shape)

    it = np.nditer(X, flags=['multi_index'])

    while not it.finished:
        
        loc = eqsim.Vec3(float(X[it.multi_index]),float(Y[it.multi_index]),0)
        disp = okada.calc_displacement_vector(loc, 1, math.pi/2, 1, 1, d,0,0,1,1)
        
        U[it.multi_index] = disp[0]
        V[it.multi_index] = disp[1]
        
        it.iternext()

    M = np.zeros(U.shape, dtype='bool')
    M[2*U.shape[0]/5:3*U.shape[0]/5,2*U.shape[1]/5:3*U.shape[1]/5] = True

    U = np.ma.masked_array(U, mask=M)
    V = np.ma.masked_array(V, mask=M)

    fig = plt.figure(figsize=(6,6))

    ax = fig.add_axes((.1,.1,.8,.8),projection='rectilinear')

    ax.quiver(X,Y,U,V,angles='xy', scale_units='xy', scale = 0.01)

    fig.savefig("disp_test/disp_%i.png"%n, format='png')

print "Creating movie"
proc_args = "ffmpeg -y -r 30 -b 64k -i disp_test/disp_%s.png disp.mp4"%("%d")
proc = subprocess.Popen(proc_args, shell=True)
proc.wait()
'''

'''
print X
print Y

loc = eqsim.Vec3(10,10,0)
print okada.calc_displacement_vector(loc, 1000, 90, 1000, 1000, 1,0,0,3e10,3e10)
'''

'''

f = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep_param-sweep/vc_05_20/ALLCAL2_1-7-11_no-creep_dyn-05_st-20_108764-277803_Events_slip-map-5.5.dat','r')

slip_map = False
event_header = False
ele_shear_sum = 0
shear_sum = 0
ele_slip_sum = 0
slip_sum = 0
enum = None
for line in f:
    
    
    if line.startswith('201'):
        slip_map = True
    if line.startswith('200'):
        event_header = True
    
    if slip_map or event_header:
        dat = line.split()
    
    if event_header:
        enum = int(dat[1])
    
    #print enum
    
    if enum == 21:
        if event_header:
            slip_sum += float(dat[13])
            shear_sum +=  float(dat[16]) - float(dat[15])
        if slip_map:
            ele_slip_sum += float(dat[6])
            ele_shear_sum += float(dat[9]) - float(dat[8])

    slip_map = False
    event_header = False

print shear_sum, ele_shear_sum
print slip_sum, ele_slip_sum
        
'''

'''
#rv_X = []
#rv_Y = []
print 'Random Variable'
for L in range(2, 10):
    #rv_X.append(L)
    #rv_Y.append(Chaos.permutationEntropy(rv,L))
    print L, Chaos.permutationEntropy(rv,L)
print

#cf_X = []
#cf_Y = []
print 'Coin Flip'
for L in range(2, 10):
    #cf_X.append(L)
    #cf_Y.append(Chaos.permutationEntropy(cf,L))
    print L, Chaos.permutationEntropy(cf,L)
'''

'''

data = np.zeros(100000)
for x in range(100000):
    if x%2 == 0:
        data[x] = -500000.0
    else:
        data[x] = -2.42144e-08

print data

cf_X = []
cf_Y = []
for L in range(2, 10):
    cf_X.append(L)
    cf_Y.append(Chaos.permutationEntropy(data,L))
    
fig = plt.figure()

ax = fig.add_axes((.1,.1,.8,.8),projection='rectilinear')

ax.plot(cf_X, cf_Y)

ax.set_xlim((0,10))
ax.set_ylim((0,10))

ax.set_ylabel('H(L)')
ax.set_xlabel('L')

ax.legend(('random walk','random var', 'coin flip'), loc=4)

fig.savefig("test.png", format='png')
'''

''' average velocities

data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_no-creep_param-sweep/vc_05_20/ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'

os.chdir(os.path.dirname(data_file))
sys_name = os.path.basename(data_file).split('.')[0]
vc_sys = cPickle.load(open('%s.pkl'%sys_name, 'rb'))

slip_velocities = []

for element in vc_sys.geometry.elements.itervalues():
    print element.slip
    slip_velocities.append( element.slip )
print
print sum(slip_velocities)/float(len(slip_velocities))
'''

''' year filter?
f = open('/Users/sachs/Documents/VirtualCalifornia/trunk/models/ALLCAL2_1-7-11_SAF-no-creep/ALLCAL2_1-7-11_SAF-no-creep_SAF-Parkfield_stresses_dyn-0-9_st-20.dat','r')


for line_num, line in enumerate(f):
    data = line.split()
    if line_num == 0:
        labels = data
        labels_filtered = [labels[1],labels[2]]
        label_indicies_filtered = [1,2]
        
        for label_num, label in enumerate(labels):
            if label.endswith('_cff'):
                labels_filtered.append(labels[label_num])
                label_indicies_filtered.append(label_num)
        
        last_year = None
        
        print ','.join(labels_filtered)
    else:
        current_year = data[1]
        if current_year != last_year:
            data_filtered = []
            for i, d in enumerate(data):
                if i in label_indicies_filtered:
                    data_filtered.append(d)
            print ','.join(data_filtered)
        
        last_year = current_year
'''