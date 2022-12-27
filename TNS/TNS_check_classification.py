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
import numpy as np
import requests
import json
from collections import OrderedDict
from astropy.time import Time
import datetime
import sys
from time import sleep

#new header needed as of 2021 May
header = {'User-Agent':'tns_marker{"tns_id":54047,"type": "bot", "name":"tess1"}'}
#header={'user-agent':'tns_marker{"tns_id":870,"type": "user", "name":"mmfausnaugh"}'}

############################# PARAMETERS #############################
# API key for Bot
with open('api_key.txt','r') as f:
    api_key = f.read().strip()#
# list that represents json file for search obj                      #
search_obj=[("ra",""), ("dec",""), ("radius",""), ("units",""),      #
            ("objname",""), ("internal_name","")]                    #
# list that represents json file for get obj                         #
get_obj=[("objname",""), ("photometry","0"), ("spectra","1")]        #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
#old, before Dec 6 2020
#url_tns_api="https://wis-tns.weizmann.ac.il/api/get"                 #
#url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"     #
#updated for their cloud migration
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
    search_data=[('api_key',(None, api_key)),                        #
                 ('data',(None,json.dumps(json_file)))]              #
    # search obj using request module                                #
    response=requests.post(search_url, 
                           headers=header,
                           files=search_data)            #
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
    get_data=[('api_key',(None, api_key)),                           #
                 ('data',(None,json.dumps(json_file)))]              #
    # get obj using request module                                   #
    response=requests.post(get_url, 
                           headers = header,
                           files=get_data)                  #
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
    response=requests.post(url, files=[('api_key',(None, api_key))], #
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



active_sectors = np.r_[58:60]
#active_sectors = [56]


if len(sys.argv) > 2:
    print('\n\ncan only specify one argument, the time_offset in days\n\n')
    raise RuntimeError                
elif len(sys.argv) == 2:
    time_offset = float(sys.argv[1])
else:
    time_offset = 30.0

print('\n\nchecking for changes {} days in the past...\n\n'.format(time_offset))


for s in active_sectors:
    #loop over cameras
    for ii in range(4):
        sleep(60)
        print('searching sector {}, camera {}'.format(s, ii+1))

        catfile = 's{:02d}/sector{}_cam{}_transients.txt'.format(s, s,ii+1)
        if os.path.isfile(catfile):
            catname = np.genfromtxt(catfile,usecols=(1),dtype=str)
            prefix  = np.genfromtxt(catfile,usecols=(0),dtype=str)
            classification = np.genfromtxt(catfile,usecols=(8),dtype=str)
            t1,t2 = np.genfromtxt(catfile,usecols=(6,7),unpack=1,dtype=str)
            t = np.array([ z[0]+'T'+z[1] for z in zip(t1,t2)])
            #            print(Time(t[0:20], format='isot',scale='utc'))

            try:
                t = Time(t, format='isot',scale='utc')
            except:
                for et in t:
                    try:
                        Time(et, format='isot',scale='utc')
                    except:
                        print(et)
                        raise
                        

            tjd = t.jd - 2457000.0

            moment = datetime.datetime.now().isoformat()
            tjd_moment = Time(moment, format='isot', scale='utc').jd - 2457000.0

            m_use = tjd > tjd_moment - time_offset

            for obj in catname[m_use]:
                print('object: {}'.format(obj),file=sys.stderr)
                get_obj = [("objname",obj)]
                sleep(3.0)
                response=get(url_tns_api, get_obj)
                try:
                    json_data2 = json.loads(response.text)
                    #print(json_data2)
                    json_data2 = json_data2['data']['reply']
                    m = np.in1d(catname, obj)

                    if json_data2['object_type']['name']  is not  None:
                        obj_type = json_data2['object_type']['name'].replace(' ','')
                    else:
                        obj_type = 'Unclassified'

                    if (prefix[m][0] != json_data2['name_prefix']) and (classification[m][0] != obj_type):
                        print('Update for {}:  {} {} --> {} {}'.format(obj,
                                                                       prefix[m][0],
                                                                       classification[m][0],
                                                                       json_data2['name_prefix'],
                                                                       obj_type))
                    

                except:
                    print('obj',obj)
                    print('response',response)
                    #continue
                    raise
