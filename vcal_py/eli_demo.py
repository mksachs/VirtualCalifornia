import site
site.addsitedir('')

# the VC modules
import VC

# 'standard' modules
import os
import cPickle

# install via macports
import h5py

# should install by doing make on new_vcal
#import quakelib


''' testing webgui output '''
data_file = '/Users/sachs/Documents/VirtualCalifornia/trunk/vcal_webgui/media/vc_data/ALLCAL2_no-creep_dt-08_st-10.h5'

os.chdir(os.path.dirname(data_file))

# grab the system name from the file name
data_file_split = os.path.basename(data_file).split('.')
tmp_sys_name = data_file_split[0]
end_index = len(tmp_sys_name)
sys_name = tmp_sys_name[0:end_index]

vc_sys = VC.VCSys(sys_name, data_file)

#event_detail = vc_sys.webgui_getEventDetail(data_file, 156)

#print event_detail

event_list = vc_sys.webgui_getEventList(0, 1000)

print event_list

#traces = vc_sys.webgui_getFaultTraces()
#center = vc_sys.webgui_getSystemCenter()

#print traces
#print
#print center

