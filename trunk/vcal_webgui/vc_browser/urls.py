from django.conf.urls import patterns, include, url

urlpatterns = patterns('vc_browser.views',
    url(r'^$', 'fault_map', name='fault_map'),
    url(r'event_list.html$', 'event_list_init', name='event_list_init'),
    url(r'get_event_list.json$', 'event_list', name='event_list'),
    url(r'event_detail.html$', 'event_detail', name='event_detail'),
    
    url(r'action.html$', 'action', name='action')
)