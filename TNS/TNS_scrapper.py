#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:38:11 2017

Developed and tested in : 
- Python version 3.6.3
- Linux Ubuntu version 16.04 LTS (64-bit)

@author: Nikola Knezevic
"""


import os
import sys
import numpy as np
import requests
import json
from collections import OrderedDict
from time import sleep

from astropy.time import Time

sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from camera_pointings import cam_pointings
from catalogs.catalog import sector_times

#new header needed as of 2021 May
header = {'User-Agent':'tns_marker{"tns_id":54047,"type": "bot", "name":"tess1"}'}
#header = {'User-Agent':'tns_marker{"tns_id":870,"type": "user", "name":"mmfausnaugh"}'}

############################# PARAMETERS #############################
# API key for Bot                                                    #
with open('api_key.txt','r') as f:
    api_key = f.read().strip()
# list that represents json file for search obj                      #
search_obj=[("ra",""), ("dec",""), ("radius",""), ("units",""),      #
            ("objname",""), ("internal_name","")]                    #
# list that represents json file for get obj                         #
get_obj=[("objname",""), ("photometry","0"), ("spectra","1")]        #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
#old url before dec 6 2020
#url_tns_api="https://wis-tns.weizmann.ac.il/api/get"                 #
#url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"     #
url_tns_api="https://www.wis-tns.org/api/get"                 #
url_tns_sandbox_api="https://sandbox.wis-tns.org/api/get"     #
######################################################################

############################# DIRECTORIES ############################
# current working directory                                          #
cwd=os.getcwd()                                                      #
# directory for downloaded files                                     #
download_dir=os.path.join(cwd,'downloaded_files')                    #
######################################################################

########################## API FUNCTIONS #############################
# function for changing data to json format                          #
def format_to_json(source):                                          #
    # change data to json format and return                          #
    parsed=json.loads(source,object_pairs_hook=OrderedDict)          #
    result=json.dumps(parsed,indent=4)                               #
    return result                                                    #
#--------------------------------------------------------------------#
# function for search obj                                            #
def search(url,json_list):                                           #
  try:                                                               #
    # url for search obj                                             #
    search_url=url+'/search'                                         #
    # change json_list to json format                                #
    json_file=OrderedDict(json_list)                                 #
    # construct the list of (key,value) pairs                        #
    search_data={'api_key': api_key,                         #
                 'data':json.dumps(json_file)}              #
    # search obj using request module      #
    response=requests.post(search_url, 
                           headers=header,
                           data=search_data)            #
#    search_data=[('api_key',(None, api_key)),                        #
#                 ('data',(None,json.dumps(json_file)))]              #
#    # search obj using request module      #
#    response=requests.post(search_url, 
#                           headers=header,
#                           files=search_data)            #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for get obj                                               #
def get(url,json_list):                                              #
  try:                                                               #
    # url for get obj                                                #
    get_url=url+'/object'                                            #
    # change json_list to json format                                #
    json_file=OrderedDict(json_list)                                 #
    # construct the list of (key,value) pairs                        #

    get_data={'api_key':api_key,'data': json.dumps(json_file)}
    response=requests.post(get_url, 
                           headers=header,
                           data=get_data)                  #
#    get_data=[('api_key',(None, api_key)),                           #
#                 ('data',(None,json.dumps(json_file)))]              #
    # get obj using request module                                   #
#    response=requests.post(get_url, 
#                           headers=header,
#                           files=get_data)                  #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for downloading file                                      #
def get_file(url):                                                   #
  try:                                                               #
    # take filename                                                  #
    filename=os.path.basename(url)                                   #
    # downloading file using request module                          #
    response=requests.post(url, files=[('api_key',(None, api_key))],
                           header=headers, #
                           stream=True)                              #
    # saving file                                                    #
    path=os.path.join(download_dir,filename)                         #
    if response.status_code == 200:                                  #
        with open(path, 'wb') as f:                                  #
            for chunk in response:                                   #
                f.write(chunk)                                       #
        print ('File : '+filename+' is successfully downloaded.')    #
    else:                                                            #
        print ('File : '+filename+' was not downloaded.')            #
        print ('Please check what went wrong.')                      #
  except Exception as e:                                             #
    print ('Error message : \n'+str(e))                              #
######################################################################


			


#gotta do s2cam4, s15 cam4, s16 3&4, s18 2,3,4, s20+4, s22 3+4, s27 3+4, s40 3+4, ,s42 4, 46 2,3,4
#active_sectors = np.r_[49,50,51]

#starting in S52, modified to only pull transients from within the last 3 months of sector start
active_sectors = np.r_[57:60]
#active_sectors = [58]

#these are imported from catalog2tess_px/camera_pointings/cam_pointings.py
cams = [cam_pointings.cam1, 
        cam_pointings.cam2, 
        cam_pointings.cam3, 
        cam_pointings.cam4]

for s in active_sectors:
    for ii,cam in enumerate(cams):
        #pick out individual cameras
        sleep(3)
        print('searching sector {}, camera {}'.format(s, ii+1))

        catfile = 's{:02d}/sector{}_cam{}_transients.txt'.format(s, s, ii+1)
        if os.path.isfile(catfile):
            catname = np.genfromtxt(catfile,usecols=(1),dtype=str, comments='@')
        else:
            catname = np.array([])
            with open(catfile,'w') as fout:
                fout.write('#%-6s %-8s %-15s %-25s %-14s %-14s %-22s %-15s %-7s %-22s %-12s %-40s %-12s\n'%(
                    'prefix','name','group','internal_name',
                    'RA','DEC',
                    'Discovery Time','Type','Mag','Filter','Redshift',
                    'Host Galaxy','Host Redshift'))

        fout = open(catfile, 'a')

        search_obj=[("ra", "{}".format(cam[s-1][0])),
                    ("dec","{}".format(cam[s-1][1])),
                    ("radius","17"),
                    ("units","degrees"),
                    ("objname",""),
                    ("internal_name",""),
                    ("discoverydate","")]                    
        print('doing cone search')
        response=search(url_tns_api,search_obj)
        print('cone search done')
        
        if None not in response:
            # Here we just display the full json data as the response
            json_data=json.loads(response.text)
            objs = np.array([ elem['objname'] for elem in json_data['data']['reply'] ])
            times_sort = np.array([ elem['objid'] for elem in json_data['data']['reply'] ])
            idx = np.argsort(times_sort)
            objs = objs[idx]
            times_sort = times_sort[idx]
            #np.savetxt('tmp.txt',np.c_[objs, times_sort],fmt='%s')
            for obj in objs:
                print(obj)
                if '2022' not in obj and '2023' not in obj:
                    continue
                
                if any(np.in1d(catname, str(obj))):
                    continue
                if len(obj) >8:
                    continue
                get_obj = [("objname",obj)]
                sleep(3.0)
                print('doing obj query')
                response=get(url_tns_api, get_obj)
                print('done with obj query')
                
                json_data2 = json.loads(response.text)
                try:
                    json_data2 = json_data2['data']['reply']
                except:
                    print(json_data2)
                    raise
                try:
                    times = json_data2['discoverydate']
                    t = Time(times.replace(' ','T'),  scale='utc',format='isot').jd - 2457000.0
                    #print(times, t, sector_times['s{:d}'.format(s)][1] + 30.0)
                    if t > sector_times['s{:d}'.format(s)][1] + 30.0:
                        print('object {} is more than 30 days after end of sector'.format(obj) )
                        break
                    print(t, sector_times['s{:d}'.format(s)][0], t < sector_times['s{:d}'.format(s)][0] - 90.0)
                    if t < sector_times['s{:d}'.format(s)][0] - 90.0:
                        continue
                    mags  = float(json_data2['discoverymag'])
                except Exception as e:
                    #print('error on {}, skipping'.format(obj))
                    print(e)
                    #raise
                    continue

                if json_data2['discmagfilter']['name'] and json_data2['discmagfilter']['family'] is not None:
                    fband = str(json_data2['discmagfilter']['name']) + '_'+str(json_data2['discmagfilter']['family'])
                else:
                    fband = 'None'
                RAs   = json_data2['ra']
                DECs  = json_data2['dec']
                if json_data2['name_prefix'] is None:
                    prefix = 'AT'
                if 'SN' not in json_data2['name_prefix']:
                    prefix  = 'AT'
                else:
                    prefix = json_data2['name_prefix']

                if json_data2['hostname'] is not None:
                    print(json_data2['hostname'], json_data2['internal_names'])
                    if isinstance(json_data2['hostname'], str):
                        host  = json_data2['hostname'].replace(' ','')
                        if len(host) == 0:
                            host = 'None'
                    else:
                        host = 'None'
                else:
                    host = 'None'

                if json_data2['host_redshift'] is not None:                    
                    host_redshift = float(json_data2['host_redshift'])
                else:
                    host_redshift = -99

                if json_data2['redshift'] is not None:                    
                    redshift = float(json_data2['redshift'])
                else:
                    redshift = -99

                group = json_data2['reporting_group']['group_name']

                if json_data2['internal_names'] is None:
                    internal_name = 'None'
                elif json_data2['discoverer_internal_name'] is None:
                    internal_name = 'None'
                elif len(json_data2['discoverer_internal_name']) == 0:
                    internal_name = 'None'
#                elif json_data2['internal_name'] is None:
#                    internal_name = 'None'
#                elif len(json_data2['internal_name']) == 0:
#                    internal_name = 'None'
                else:
                    try:
                        internal_name = json_data2['internal_name']
                        internal_name = internal_name.replace(' ','-') 
                    except:
                        internal_name = json_data2['discoverer_internal_name']
                        internal_name = internal_name.replace(' ','-') 
                        

                if json_data2['object_type']['name'] is not None:
                    obj_type = json_data2['object_type']['name'].replace(' ','')
                else:
                    obj_type = 'Unclassified'
                    
                print('  ',internal_name,  mags, obj, times)

                fout.write('%6s %8s %15s %25s %14s %14s %22s %15s %7s %22s %12s %40s %12s\n'%(
                    prefix, obj, group, internal_name,
                    RAs, DECs,
                    times, obj_type, mags, fband, redshift,
                    host,host_redshift))
                fout.flush()
                
        else:
            print (response[1])

        fout.close()
