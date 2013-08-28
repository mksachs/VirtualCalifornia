from vc_browser.models import VCModel
from django.shortcuts import render, render_to_response, get_object_or_404
from django.shortcuts import HttpResponseRedirect
from django.template import RequestContext
from django.conf import settings
from django.http import HttpResponse

import site
import sys
import os
import cPickle
import json

site.addsitedir('')

import VC

import cProfile
import pstats

# this handles ajax calls for additional event lists. used for paging, sorting etc. sends a json object to the ui for display.
def event_list(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    sort_fields = { 'evid':             0,
                    'evyear':           1,
                    'evmag':            2,
                    'evtriggersec':     3,
                    'evtriggerele':     4,
                    'evinvolvedsec':    5,
                    'evinvolvedele':    6,
                    'evaveslip':        7
                    }
    
    vc_sys = getVCSys(request)
    
    return HttpResponse(
        model_list[0].event_list(
            vc_sys,
            number_per_page=int(request.GET['items_to_display']),
            sort_field=sort_fields[request.GET['sort_column']],
            sort_direction=request.GET['sort_direction'],
            page_number=int(request.GET['page_number'])
            ),
        mimetype="application/json")

# this loads the event list template. used when the event list is first requested.
def event_list_init(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    options = []
    for option in map( int , model_list[0].events_per_page_options.split(',') ):
        options.append( option )
    
    vc_sys = getVCSys(request)

    return render_to_response('vc_browser/event_list.html',
        {   'system_info': model_list[0].system_info(vc_sys),
            'system_name':model_list[0].name,
            'event_list':model_list[0].event_list(vc_sys, number_per_page=model_list[0].default_events_per_page),
            'events_per_page_options': options,
            'default_events_per_page':model_list[0].default_events_per_page
        })


''' this is for profiling
def event_list_init(request):
    cProfile.runctx('doStuff(request)',globals(),{'request':request},'profile')
    stats = pstats.Stats('profile')
    # Clean up filenames for the report
    stats.strip_dirs()

    # Sort the statistics by the cumulative time spent in the function
    stats.sort_stats('cumulative')

    stats.print_stats()


def doStuff(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    options = []
    for option in map( int , model_list[0].events_per_page_options.split(',') ):
        options.append( option )
    
    vc_sys = getVCSys(request)

    return render_to_response('vc_browser/event_list.html',
        {   'system_info': model_list[0].system_info(vc_sys),
            'system_name':model_list[0].name,
            'event_list':model_list[0].event_list(vc_sys, number_per_page=model_list[0].default_events_per_page),
            'events_per_page_options': options,
            'default_events_per_page':model_list[0].default_events_per_page
        })
'''

# this handles the ajax call for the event detail. unlike the event list view, this returns html from the given template
# to be loaded into an element by the ui.
def event_detail(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    vc_sys = getVCSys(request)
    
    event_detail = model_list[0].event_detail(vc_sys, int(request.GET['event_num']))
    return render_to_response('vc_browser/event_detail.html',
        {   'event_detail': event_detail
        },
        context_instance=RequestContext(request))

def fault_map(request):
    # right now we only have one model
    
    model_list = VCModel.objects.all()
    
    #print model_list[0].fault_map_data()
    
    vc_sys = getVCSys(request)
    
    return render_to_response('vc_browser/fault_map.html',
        {   'system_info': model_list[0].system_info(vc_sys),
            'system_name':model_list[0].name,
            'fault_map_data': model_list[0].fault_map_data(vc_sys)
        })

def action(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    try:
        fm = True if request.GET['frequency_magnitude'] == 'true' else False
    except KeyError:
        fm = False

    try:
        mra = True if request.GET['magnitude_rupture_area'] == 'true' else False
    except KeyError:
        mra = False

    try:
        mas = True if request.GET['magnitude_average_slip'] == 'true' else False
    except KeyError:
        mas = False

    try:
        asrl = True if request.GET['average_slip_rupture_length'] == 'true' else False
    except KeyError:
        asrl = False

    try:
        st = True if request.GET['space_time'] == 'true' else False
    except KeyError:
        st = False
    
    try:
        start_year = float(request.GET['start_year'])
    except KeyError:
        start_year = 0

    try:
        end_year = float(request.GET['end_year'])
    except KeyError:
        end_year = None

    try:
        server_data = json.loads(request.GET['server_data'])
    except KeyError:
        server_data = None

    action_params = { 'fm':fm,
                        'mra':mra,
                        'mas':mas,
                        'asrl':asrl,
                        'st':st,
                        'start_year':start_year,
                        'end_year':end_year,
                        'server_data':server_data
                    }

    vc_sys = getVCSys(request)
        
    #action_detail = None
    #cProfile.runctx('action_detail = model_list[0].get_action(action_params,vc_sys)',globals(),{'action_params':action_params, 'action_detail':action_detail, 'model_list':model_list, 'vc_sys':vc_sys})

    action_detail = model_list[0].get_action(action_params, vc_sys)
    
    return render_to_response('vc_browser/action.html',
        {   'action_detail':action_detail,
            'system_info': model_list[0].system_info(vc_sys),
            'system_name':model_list[0].name,
            'from_view':server_data['from_view']
        }, context_instance=RequestContext(request))

def getVCSys(request):
    # right now we only have one model
    model_list = VCModel.objects.all()
    
    #del request.session['vc_sys']
    
    '''
    try:
        vc_sys = request.session['vc_sys']
    except KeyError:
        print '*** load vc_sys ***'
        request.session['vc_sys'] = VC.VCSys(model_list[0].system_name(), model_list[0].hdf5_output_file.path)
        #request.session['vc_sys'] = cPickle.load(open('%s%s.pkl'%(model_list[0].cache_path(),model_list[0].system_name()), 'rb'))
        vc_sys = request.session['vc_sys']
    '''
    vc_sys = VC.VCSys(model_list[0].system_name(), model_list[0].hdf5_output_file.path)
    return vc_sys
