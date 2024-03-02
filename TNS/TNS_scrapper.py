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
url_tns_api_search="https://www.wis-tns.org/search?&"                 #
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
#active_sectors = np.r_[69:74]
active_sectors = [72]

#these are imported from catalog2tess_px/camera_pointings/cam_pointings.py
cams = [cam_pointings.cam1, 
        cam_pointings.cam2, 
        cam_pointings.cam3, 
        cam_pointings.cam4]

for s in active_sectors:
    for ii,cam in enumerate(cams):

        #params to access TNS api
        #Shri had these as 500 in his script,

        #need to reset inorder to loop through pages
        NMAX = "1000"
        PAGENO = 0
        tmpfile = "tmp_tns_file_{}_cam{}.tsv".format(s, ii + 1)
        savefile = "save_tns_file_{}_cam{}.tsv".format(s, ii + 1)

        #pick out individual cameras
        sleep(61.0)
        print('searching sector {}, camera {}'.format(s, ii+1))

        catfile = 's{:02d}/sector{}_cam{}_transients.txt'.format(s, s, ii+1)
        if os.path.isfile(catfile):
            catname = np.genfromtxt(catfile,usecols=(1),dtype=str, comments='@')
        else:
            catname = np.array([])
            with open(catfile,'w') as fout:
                fout.write('#%-6s %-8s %-15s %-25s %-14s %-14s %-22s %-15s %-7s %-22s %-12s %-40s %-12s %-25s\n'%(
                    'prefix','name','group','internal_name',
                    'RA','DEC',
                    'Discovery Time','Type','Mag','Filter','Redshift',
                    'Host Galaxy','Host Redshift','bibcode'))

        fout = open(catfile, 'a')

        #anything discovered 5 days prior to sector star for followup
        sector_time_start = Time( sector_times['s{:d}'.format(s)][0] + 2457000 - 5,
                                  format='jd').isot
        #anything discovered 30 days after TESS, for precovery
        sector_time_end   = Time( sector_times['s{:d}'.format(s)][1] + 2457000 + 30,
                                  format='jd').isot

        ra_search = cam[s-1][0]
        dec_search = cam[s-1][1]
        
        print(cam[s-1][0], cam[s-1][1], sector_time_start, sector_time_end.split('T')[0])


        #update on 2024-03-01; to make search
        #managable, I need to limit by the few months around the times of TESS
        #observations.  Shri Kulkarni's qTNSm.sh tool is the only thing
        #I have found that can do this.  Need to construct the URL by hand
        #and submit by curl

        start_time = sector_time_start.split('T')[0]
        end_time = sector_time_end.split('T')[0]


        #need to loop over PAGENO in query, until no records are returned
        n_rows = 1000
        while n_rows > 1:
            sleep(5.0)
            s2="date_start%5Bdate%5D=" + start_time + "&date_end%5Bdate%5D=" + end_time
            s3="&ra={}&decl={}&radius=20&coords_unit=deg&".format(ra_search, dec_search)
            s4="&num_page=" + NMAX+ "&page=" + str(PAGENO) + "&format=tsv"

            curl_string = "curl"+  " -s"+  " -o " +  tmpfile +\
                          " -H " + "'User-Agent: " + json.dumps(header["User-Agent"]) + "'"+\
                          " -d"+ " 'api_key=" + api_key + "' " +\
                            "'" + url_tns_api_search  +  s2 + s3 + s4 + "'"
                
            print(curl_string)
            os.system(curl_string)
            #will still have a header
            n_rows = 1
            with open(tmpfile,'r') as fin:
                for row in fin:
                    #check for html 
                    row.strip()
                    #weird white space remains, seems to be ^M character
                    #this is more reliable, I guess?
                    if len(row) > 10 and '<' not in row and \
                       '@' not in row and '{' not in row and \
                       'http' not in row and 'ID' not in row:
                            n_rows += 1
                            #split and write to fout, which is the catfile
                            data = row.split("\t")

                            id_num = data[0].replace('"','')
                            obj = data[1].replace('"','')
                            prefix, obj = obj.split()
                            #do not write if it is already in the catalog
                            if any(np.in1d(catname, str(obj))):
                                continue

                            times = data[20].replace('"','')
                            t = Time(times.replace(' ','T'),  scale='utc',format='isot').jd - 2457000.0
                            mags  = data[18].replace('"','')

                            fband = data[19].replace('"','')
                            if len(fband) == 0:
                                fband='None'

                            RAs   = data[2].replace('"','')
                            DECs  = data[3].replace('"','')

                            host = data[6].replace('"','')
                            if len(host) == 0:
                                host = 'None'
                            else:
                                host = host.replace(' ','')

                            host_redshift = data[7].replace('"','')
                            if len(host_redshift) == 0:
                                host_redshift= -99
                            else:
                                host_redshift = float(host_redshift)


                            redshift = data[5].replace('"','')
                            if len(redshift) == 0:
                                redshift = -99
                            else:
                                redshift = float(redshift)


                            group = data[8].replace('"','')
                            group = group.replace(', ',',')

                            internal_name = data[12].replace('"','')
                            if len(internal_name) == 0:
                                internal_name = 'None'
                            internal_name = internal_name.replace(' ','-')

                            obj_type = data[4].replace('"','')
                            if len(obj_type) == 0:
                                obj_type = 'Unclassified'
                            else:
                                obj_type = obj_type.replace(' ','')


                            bibcode = data[24].replace('"','')
                            
                            print('  ',internal_name,  group, mags, obj, times,bibcode)


                            fout.write('%6s %8s %35s %25s %14s %14s %22s %15s %7s %22s %12s %40s %12s %25s\n'%(
                                prefix, obj, group, internal_name,
                                RAs, DECs,
                                times, obj_type, mags, fband, redshift,
                                host,host_redshift,bibcode))


            PAGENO += 1
            print(PAGENO, n_rows)


        print('exited loop')
        #os.remove(tmpfile)
              
        fout.close()
