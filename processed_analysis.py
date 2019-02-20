#!/usr/bin/python


#processed_analysis.py compiles processed data from different experiments,
# groups the data by experiment type, and performs summary analyses

#Python 2.7

##USAGE 
#python processed_analysis.py processed_date_range.txt

##A text file (e.g. 'processed_date_range.txt' specifies the dates, or range of dates,
#over which to analyze data.
#A single date specified as YYYY.MM.DD
#A range of dates specified as YYYY.MM.DD-YYYY.MM.DD




import numpy as np
import os
import pylab
import pdb
import pickle
from py_utilities import circ
import scipy.stats.morestats as st
from py_utilities import tw_plot_library3 as plt
import datetime
from py_utilities import tw_filehandling as fh
from py_utilities import tw_calc_library as calc
import datetime
DOWNSAMPLE_VEC_STRENGTH=True
DOWNSAMPLE_SAMPLE_RATE=10.0
SAVE_VEC_STRENGTH=True
ADD_START_ANGLE=True
ADD_TIME_OF_DAY=True
HEAT_ANALYSIS=False
REMOVE_STRIPE_OUTLIERS=True
NUM_DATA_TYPES=5
base_path='/users/tim/'


input_data_path=base_path+'Dropbox/flystuff/uo_data/'
sum_data_path=base_path+'Dropbox/flystuff/'


EXP_KEYS=['es_uaskirj','ss96_uaskirj','rna_eyeless_flies','rna_control_flies','ss96toy_lc_flies','estoy_lc_flies','mcherryhs_flies','unclassified']
SUMSTATS_KEYS_TO_REMOVE=['lftwng','gain_coefficient','raw_time','mot_position','pol_sensor','time_in_min','mot_deg','rtwng','mot_rad','normhst','mn_dot_product_lst','heat_map','heat_map_vls','mot_vel', 'xrad','vechst','mot_position_in_rad']
SPLITSTATS_KEYS_TO_REMOVE=['realigned_norm_heat_map_vls','norm_heat_map_vls','r','mot_deg','theta','rtwng','mot_rad','vec_time_lst','len_vector_lst','mn_vector_lst','normhst','mn_dot_product_lst','heat_map','heat_map_vls','mot_vel','vechst','xrad']
ADD_ALL_SEGMENTS=True
RESTRICT_TO_SHARED_FIRST_AND_SECOND_SEGMENT=True

class Closed_Anal():

    #calls prep_data_fxn for each file
    def __init__ (self,input_args):
        self.exec_file=os.path.abspath(__file__)
        
        self.day_array,self.exp_type,self.anal_type=fh.read_days(input_args[1])

        if len(input_args)==2:
            self.save_dt_name=[]
        else:
            self.save_dt_name=input_args[2]


    def initialize_variables(self):
        tst=1

    def run_analysis(self):


        indata_str=['_pair_sumdata.pck']

        savedata_str=['vecstrength_unpaired_sumsundata.pck']


        for indvl,foo in enumerate(indata_str):

            self.cr_indata_str=indata_str[indvl]

            if not self.save_dt_name:
                cr_savedata_str=savedata_str[indvl]
            else:
                cr_savedata_str=self.save_dt_name
            fig=[]
            outdt={}
            self.sumdt={}
            self.adt={}
            self.initialize_sumdt_keys()
            self.flyct=0

            for crday in self.day_array:

                self.crday=pylab.num2date(crday)

                self.get_data_by_day()

            self.parse_summary_data()
            
            pckfname=sum_data_path+'summary_data/'+cr_savedata_str

\
            outdt['sumdt']=self.reduce_size(self.sumdt)
            outdt['adt']=self.adt
            

            fh.save_to_file(pckfname,outdt,save_executable=self.exec_file)
            

    def initialize_sumdt_keys(self):
        for crkey in EXP_KEYS:
            self.sumdt[crkey]=[]
            self.adt[crkey]={}
            self.adt[crkey]['vec_across_pair']=[]
            for stim_type in np.arange(NUM_DATA_TYPES):
                self.adt[crkey][stim_type]={}

                for dat_key in ['vec_strength','mnrad','flynum']:
                    self.adt[crkey][stim_type][dat_key]=[]
                for dat_key in ['paired_rad','paired_vec']:
                    self.adt[crkey][stim_type][dat_key]={}
                    for ind in [0,1,'flynum']:
                      self.adt[crkey][stim_type][dat_key][ind]=[]  

    def get_data_by_day(self):

        cr_data_path=input_data_path
        self.datapath=fh.make_data_path(cr_data_path,self.crday)
        dirlist=fh.get_directories(self.datapath)
        print 'getting data from %s'%self.datapath

        if dirlist:
            for self.crflyname in dirlist:

                self.get_data_by_fly()
    
    #for each fly append to self.sumdt

    def reduce_size(self,sumdt):

        outdt=sumdt

        flylist=[x for x in self.sumdt.keys() if type(x)  is int]
        for crflynum in flylist:
            crfly=self.sumdt[crflynum]
            explist=[x for x in crfly.keys() if type(x)  is int]
            for crexpnum in explist:
                try:
                    crexp=crfly[crexpnum]
                except:
                    print('error opening data')

                
                for key in outdt[crflynum][crexpnum]['sumstats'].keys():
                    for crkey_to_remove in SUMSTATS_KEYS_TO_REMOVE:    
                        outdt[crflynum][crexpnum]['sumstats'][key][crkey_to_remove]=[]
        return outdt

    def get_vec_dat(self,crdt):
        
        return vecdt

    def get_data_by_fly(self):
        #print
        EXPCT=0
        self.suffix=self.crflyname+'/cloop/'

        files = sorted(os.listdir(self.datapath+self.suffix))
        crdt={}

        pck_files = [f for f in files if f.find(self.cr_indata_str)>0]


        for self.crfile in pck_files:

            crdata=fh.open_pickle(self.datapath+self.suffix+self.crfile)
            
            if crdata:
                crdt[EXPCT]=crdata
                if ADD_START_ANGLE:
                    skip_flag=False

                    if 'exp' in crdata['sumstats'].keys():
                        crkey='exp'
                    elif 'ctrl' in crdata['sumstats'].keys():
                        crkey='ctrl'
                    else:
                        skip_flag=True
                    if not skip_flag:
                        cr_rad=crdata['sumstats'][crkey]['mot_rad']
                        mnvl=circ.circmean(cr_rad[0:2])
                        crdt[EXPCT]['sumstats'][crkey]['start_rad']=circ.circmean(mnvl)

                        if mnvl>100:
                            print('mnvl greater than 100')
                if ADD_TIME_OF_DAY:

                    if 'exp' in crdata['sumstats'].keys():
                        crkey='exp'
                    elif 'ctrl' in crdata['sumstats'].keys():
                        crkey='ctrl'
                    else:
                        skip_flag=True
                    if not skip_flag:
                        tst=1

                EXPCT=EXPCT+1

        if crdt:
            crdt=self.add_total_vec_strength_across_paired(crdt)
           
            self.sumdt[self.flyct]=crdt

            self.sumdt[self.flyct]['crday']=self.crday

            self.flyct=self.flyct+1


    def add_total_vec_strength_across_paired(self,crdt):

        deg_per_bin=10
        degbins = np.arange(0,370,deg_per_bin)

        min_time=4
        

        for crind in crdt.keys():
            try:
                if crdt[crind]['sumstats'][0]['time_in_min'][-1]-crdt[crind]['sumstats'][0]['time_in_min'][0]:
                    if crdt[crind]['sumstats'][1]['time_in_min'][-1]-crdt[crind]['sumstats'][1]['time_in_min'][0]:
                        COMBINED_DATA_FLAG=True
                    else:
                        COMBINED_DATA_FLAG=False
                else:
                    COMBINED_DATA_FLAG=False
            except:
                COMBINED_DATA_FLAG=False
            
            #check if there  is valid, combined_data
            if COMBINED_DATA_FLAG:
                combined_deg=np.concatenate((crdt[crind]['sumstats'][0]['mot_deg'],crdt[crind]['sumstats'][1]['mot_deg']))


                hstout_dt=calc.make_hist_calculations(combined_deg,degbins)
                
                crdt[crind]['sumstats'][0]['vec_across_pair']=hstout_dt['circvar_mod360']

        return crdt
    def heat_analysis(self):
        self.heatdt['average_heat_map']={}
        self.heatdt['split_average_heat_map']={}
        self.heatdt['gain_values']={}
        self.heatdt['num_flies']={}
        for vlind,self.crvl_type in enumerate(EXP_KEYS):
            CT=0
            self.initial_flag=1

            self.heatdt['average_heat_map'][self.crvl_type]={}
            self.heatdt['num_flies'][self.crvl_type]={}


            for crflyind in self.sumdt[self.crvl_type]:
                
                crheatmap=self.sumdt[crflyind][0]['sumstats']['stripe'][0]['mnrad']
                self.heatdt['average_heat_map'][self.crvl_type]=self.add_heat_map_value(self.heatdt['average_heat_map'][self.crvl_type],crheatmap,crflyind)
                self.heatdt['num_flies'][self.crvl_type]=CT
                CT=CT+1



    def add_heat_map_value(self,cr_ave_dt,crheatmap,crflyind):
        if np.isnan(np.sum(crheatmap['norm_heat_map_vls'])):
            print('no heat map value')

        if self.initial_flag:

             cr_ave_dt=crheatmap

             cr_ave_dt['num_vls']=0
             self.initial_flag=0
             cr_ave_dt['ratio']=[]
             cr_ave_dt['flyind']=[]
             cr_ave_dt['num_vls']=cr_ave_dt['num_vls']+1
        else:

            for crkey in ['realigned_norm_heat_map_vls','norm_heat_map_vls']:

                cr_ave_dt[crkey]=cr_ave_dt[crkey]+crheatmap[crkey]
            cr_ave_dt['num_vls']=cr_ave_dt['num_vls']+1
            if crkey is 'norm_heat_map_vls':
                cr_ratio=np.sum(crheatmap['norm_heat_map_vls'][0:37,:])/np.sum(crheatmap['norm_heat_map_vls'][37:72,:])
                cr_ave_dt['ratio'].append(cr_ratio)
                if self.crvl_type is 'all_linear_polarizer_flies':
                    if np.isinf(cr_ratio):
                        print('ratio error')
                cr_ave_dt['flyind'].append(crflyind)



        return cr_ave_dt



    def parse_summary_data(self):

        self.initialize_variables()
        self.flylist=[x for x in self.sumdt.keys() if type(x)  is int]
        
        for self.crflynum in self.flylist:

            self.gal_tnt_ind=[]
            self.rna_control_ind=[]
            self.tb_ind=[]

            self.crfly=self.sumdt[self.crflynum]


            self.assign_data_to_experimental_groups()

    def assign_data_to_experimental_groups(self):
        minimum_duration={}
        for ind in [0,1,2]:
            minimum_duration[ind]=4
        for ind in [3,4]:
            minimum_duration[ind]=0.9
        try:
            exp_type=self.crfly[0]['params']['fly_type']
        except:
            exp_type='unclassified'
  
        if 'ss96es_rnai' in exp_type:
            cr_exp_type='toyes_rnai_flies'
        

        
        elif 'es_uaskirj' in exp_type:
            cr_exp_type='es_uaskirj'
        elif 'ss96_uaskirj' in exp_type:
            cr_exp_type='ss96_uaskirj'
       

        elif 'estoy_lc' in exp_type:
            cr_exp_type='estoy_lc_flies'

        elif 'ss96toy_lc' in exp_type:
            cr_exp_type='ss96toy_lc_flies'

       

        elif 'eyrnai_con' in exp_type:
            cr_exp_type='ey_rnai_con_flies'
        elif 'mcherryhs' in exp_type:
            cr_exp_type='mcherryhs_flies'

        elif 'eyrnai_exp' in exp_type:
            cr_exp_type='eyrnai_exp_flies'
        else:
            cr_exp_type='unclassified'

        num_exps_added={}
        crheading={}
        crvec_strength={}
        data_added_flag=False
        for stim_type in np.arange(NUM_DATA_TYPES):
            num_exps_added[stim_type]=0
            crheading[stim_type]=[]
            crvec_strength[stim_type]=[]
        self.sumdt[cr_exp_type].append(self.crflynum)
        exps=[f for f in self.sumdt[self.crflynum].keys() if type(f) is int]

        for crexp in exps:
            for crtype in self.sumdt[self.crflynum][crexp]['sumstats'].keys():
                try:
                    crdt=self.sumdt[self.crflynum][crexp]['sumstats'][crtype]
                except:
                    tst=1

                try:
                    max_time=crdt['time_in_min'][-1]-crdt['time_in_min'][0]
                except:
                    max_time=0

                #'vec_across_pair' only saved at 0th index
                try:
                    if crtype==0:

                        crvec_across_pair=crdt['vec_across_pair']
                    
                except:
                    
                    crvec_across_pair=np.nan
                if max_time>minimum_duration[crtype]:
                #                         if num_exps_added[crtype]<2:

                    crvec=1-crdt['circvar_mod360']
                    crmn=crdt['mnrad_360']
                    if REMOVE_STRIPE_OUTLIERS:
                        if crvec>0.99:
                            crvec=np.nan
                            crmn=np.nan
                    crvec_strength[crtype].append(crvec)
                    crheading[crtype].append(crmn)

                    num_exps_added[crtype]=num_exps_added[crtype]+1
                    data_added_flag=True
                else:
                    tst=1

            if data_added_flag:

                self.fill_in_with_nans(crheading,crvec_strength)
                self.add_to_master_list(crheading,crvec_strength,cr_exp_type,crvec_across_pair)

            else:
                tst=1
    def fill_in_with_nans(self,crheading,crvec_strength):
        for crtype in np.arange(NUM_DATA_TYPES):
            heading_len=len(crheading[crtype])
            vec_len=len(crvec_strength[crtype])
            add_heading_num=1-heading_len
            add_vec_num=1-vec_len

            for crvl in np.arange(add_heading_num):
                crheading[crtype].append(np.nan)
            for crvl in np.arange(add_vec_num):

                crvec_strength[crtype].append(np.nan)

    def add_to_master_list(self,crheading,crvec_strength,cr_exp_type,crvec_across_pair):
        for crtype in np.arange(NUM_DATA_TYPES):
            
            crvec=crvec_strength[crtype]
            crhead=crheading[crtype]

            if len(crvec)>1:
                #this is a problem b/c there should only be one value
                
                not_nan_ind=np.where(~np.isnan(crvec))[0]
                if not_nan_ind:
                    crvec_strength[crtype]=[crvec[not_nan_ind[0]]]
            if len(crhead)>1:
                #this is a problem b/c there should only be one value
                
                not_nan_ind=np.where(~np.isnan(crhead))[0]
                if not_nan_ind:
                    crheading[crtype]=[crhead[not_nan_ind[0]]]
            
            if not RESTRICT_TO_SHARED_FIRST_AND_SECOND_SEGMENT:
                
                if crtype==1:
                    if np.isnan(crheading[1]):
                        PROCEED_FLAG=False
                    else:
                        PROCEED_FLAG=True
                else:
                    PROCEED_FLAG=True
                


                if PROCEED_FLAG:
                    
                    if not np.isnan(crheading[crtype]):
                        
                        self.adt[cr_exp_type][crtype]['mnrad'].append(crheading[crtype])


                        self.adt[cr_exp_type][crtype]['vec_strength'].append(crvec_strength[crtype])
                        self.adt[cr_exp_type][crtype]['flynum'].append(self.crflynum)
                        if crtype==1:
                            self.adt[cr_exp_type]['vec_across_pair'].append(crvec_across_pair)
                    #add both first and second value to paired data structure
                    #if crtype is 1 agnd proceed_flag is true
                    if crtype==1:
                        
                        for crind in [0,1]:
                            self.adt[cr_exp_type][crind]['paired_rad'][0].append(crheading[0])
                            self.adt[cr_exp_type][crind]['paired_vec'][0].append(crvec_strength[0])
                            self.adt[cr_exp_type][crind]['paired_rad']['flynum'].append(self.crflynum)
                            self.adt[cr_exp_type][crind]['paired_rad'][1].append(crheading[1])
                            self.adt[cr_exp_type][crind]['paired_vec'][1].append(crvec_strength[1])
            else:
                
                if np.isnan(crheading[1]):
                        tst=1
                        #DO NOT ADD value because missing paired value
                else:
                        self.adt[cr_exp_type][crtype]['mnrad'].append(crheading[crtype])


                        self.adt[cr_exp_type][crtype]['vec_strength'].append(crvec_strength[crtype])
                        self.adt[cr_exp_type][crtype]['flynum'].append(self.crflynum) 
              

if __name__ == '__main__':
    import sys
    print sys.argv[1]
    closed_anal=Closed_Anal(sys.argv)
    closed_anal.run_analysis()
