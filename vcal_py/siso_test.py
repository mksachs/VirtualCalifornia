#!/usr/bin/env python
import site

site.addsitedir('')
import Converter

#events_f = open('/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_eqsim/110912_471207_events_slip-map-5.5.dat','r')
events_f = open('/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_eqsim/110912_471207_events_slip-map-5.5_eid-fix.dat','r')
model_f = open('/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/ALLCAL2_1-7-11_no-creep_Geometry.dat','r')
element_remap_f = open('/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_element-remap.dat','r')

vc_to_eqsim_eid_remap = {}
for ln, line in enumerate(element_remap_f):
    if ln != 0:
        eqsim_eid, vc_eid = map(int, line.split())
        vc_to_eqsim_eid_remap[vc_eid] = eqsim_eid
        
''' test slip in slip out '''
si = {}

for line in model_f:
    if line.startswith('204'):
        dat = line.split()
        si[int(dat[1])] = float(dat[7])

event_times = []
so = {}
for line in events_f:
    if line.startswith('200'):
        event_times.append(float(line.split()[3]))
    elif line.startswith('201'):
        dat = line.split()
        #eid = vc_to_eqsim_eid_remap[int(dat[12])]
        eid = int(dat[12])
        slip = float(dat[6])
        try:
            so[eid] += slip
        except KeyError:
            so[eid] = slip

sorted_event_times = sorted(event_times)
total_time_sec = sorted_event_times[-1] - sorted_event_times[0]

print 'eid', 'so', 'si'

for eid in sorted(so.iterkeys()):
    print eid, so[eid]/total_time_sec, si[eid]
    #try:
    #    print eid, so[eid]/total_time_sec, si[eid]
    #except:
    #    print eid, 0.0, si[eid]


''' fix the element id issue 
events_fix_f = open('/Users/sachs/Documents/VirtualCaliforniaSVN/trunk/models/ALLCAL2_1-7-11_no-creep/run7/ALLCAL2_no-creep_dt-08_st-10_eqsim/110912_471207_events_slip-map-5.5_eid-fix.dat','w')

for line in events_f:
    if line.startswith('201'):
        dat = line.split()
        eid = int(dat[12])
        events_fix_f.write('%s %s %s %s %s %s %s %s %s %s %s %s %i\r'%(dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9], dat[10], dat[11], vc_to_eqsim_eid_remap[eid]))
    else:
        events_fix_f.write('%s'%line)

events_fix_f.close()
'''