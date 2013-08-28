from django.db import models
import os
from django.conf import settings
import cPickle
import h5py

class VCModel(models.Model):
    name = models.CharField(max_length=200)
    date_added = models.DateTimeField('date added')
    hdf5_output_file = models.FileField(upload_to='vc_data/')
    events_per_page_options = models.CommaSeparatedIntegerField(max_length=100)
    default_events_per_page = 1000
    
    def system_name(self):
        return os.path.basename(self.hdf5_output_file.path).split('.')[0]
    
    def cache_path(self):
        return '%svc_data/%s_cache/'%(settings.MEDIA_ROOT, self.system_name())
    
    def event_list(self, vc_sys, page_number=0, number_per_page=500, sort_field=0, sort_direction='ascending'):
        event_list = vc_sys.webgui_getEventList(page_number, number_per_page, sort_field, sort_direction)
        return event_list
    
    def event_detail(self, vc_sys, evid):
        try:
            event_detail = cPickle.load(open('%swebgui/event-detail/event-detail_%s.pkl'%(self.cache_path(), evid),'rb'))
        except IOError:
            #vc_sys = cPickle.load(open('%s%s.pkl'%(self.cache_path(),self.system_name()), 'rb'))
    
            event_detail = vc_sys.webgui_getEventDetail(evid)
            
            if not os.path.exists('%swebgui/event-detail'%(self.cache_path())):
                os.makedirs('%swebgui/event-detail'%(self.cache_path()))
            
            cPickle.dump(event_detail, open('%swebgui/event-detail/event-detail_%s.pkl'%(self.cache_path(), evid), 'wb'))
        return event_detail
    
    def fault_map_data(self, vc_sys):
        try:
            fault_map_data = cPickle.load(open('%swebgui/fault-map-data.pkl'%(self.cache_path()),'rb'))
        except IOError:
            #vc_sys = cPickle.load(open('%s%s.pkl'%(self.cache_path(),self.system_name()), 'rb'))

            fault_map_data = vc_sys.webgui_getFaultTraces()
            
            if not os.path.exists('%swebgui'%(self.cache_path())):
                os.makedirs('%swebgui'%(self.cache_path()))
            
            cPickle.dump(fault_map_data, open('%swebgui/fault-map-data.pkl'%(self.cache_path()), 'wb'))

        return fault_map_data
    
    def system_info(self, vc_sys):
        try:
            system_info = cPickle.load(open('%swebgui/system-info.pkl'%(self.cache_path()),'rb'))
        except IOError:
            
            f = h5py.File(self.hdf5_output_file.path, 'r')
            
            events = f['event_table']
            
            system_info = {'sections': len(vc_sys.geometry.sections), 'elements': len(vc_sys.geometry.elements), 'events':len(events), 'years':events[-1][1], }
            
            if not os.path.exists('%swebgui'%(self.cache_path())):
                os.makedirs('%swebgui'%(self.cache_path()))
            
            cPickle.dump(system_info, open('%swebgui/system-info.pkl'%(self.cache_path()), 'wb'))
            
            f.close()
        
        return system_info
    
    def get_action(self, action_params, vc_sys):
        #vc_sys = cPickle.load(open('%s%s.pkl'%(self.cache_path(),self.system_name()), 'rb'))
        
        start_year = action_params['start_year']
        end_year = action_params['end_year']
        
        if action_params['server_data']['from_view'] == 'fault_map':
            selected_sections = action_params['server_data']['selected']
            
            section_filter = []
            
            for sid in sorted(map(int, selected_sections.keys())):
                if selected_sections['%i'%sid]:
                    section_filter.append(sid)
        
            vc_sys.section_filter = section_filter
        
            if start_year == 0:
                start_evid = 0
            else:
                start_evid  = vc_sys.eventForYear(float(start_year))
            if end_year is not None:
                end_evid  = vc_sys.eventForYear(float(end_year))
            else:
                end_evid = None
            
        action_plots = vc_sys.webgui_getActionPlots(start_evid=start_evid, end_evid=end_evid,
                                                        fm=action_params['fm'],
                                                        mra=action_params['mra'],
                                                        mas=action_params['mas'],
                                                        asrl=action_params['asrl'],
                                                        st=action_params['st']
                                                    )
        
        return action_plots

    def __unicode__(self):
        return self.name

