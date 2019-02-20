#write_elife_repo_csv.py used to create the individual csv files in the Dryad Repository
#for Sullivan, Warren, Doe 2019.

#Python 2.7

#DEPENDENCIES
#various libraries from https://github.com/timothylwarren/py_utilities

##USAGE 
#python write_elife_repo_csv.py

##FLY_TYPE_TO_SAVE variable sets the data type for csv file



import numpy as np
import os
import pylab
import pdb
import pickle
import datetime
import csv
from py_utilities import set_up_figure as setup
from py_utilities import tw_filehandling as fh
from py_utilities import tw_calc_library as calc



inpath='/users/tim/'
save_file_folder=inpath+'Dropbox/flystuff/elife_repository_data/'
dropbox_path=inpath
path=setup.set_up_path(inpath)
vl_types_to_plot=['mcherryhs_flies','ey_rnai_con_flies','eyrnai_exp_flies','estoy_lc_flies','ss96toy_lc_flies','ss96_uaskirj','es_uaskirj']
#inds_to_plot=['uninterrupted_paired_inds','interrupted_paired_inds','all_paired_inds']
COLVLS={}
COLVLS['uninterrupted']='k'
COLVLS['interrupted']='c'
#COLVLS=['k','c','r','g']
combine_colvls=['k','r', 'c','b','m','g','k','k','k']
PLOT_START_ANGLE=False
PLOT_VECTOR_STRENGTH=True
SHORT_PAIRED_DELAY=[14,35]
NUM_EXAMPLE_FLIES=5
OMIT_FIRST_HIST=True
AXISPAD=2
PLOT_POSTER=False
DO_EXTRA_ANALYSIS=False
PLOT_PAIRED_DIFF=True
ADD_EXP_FIT=True

FLY_TYPE_TO_SAVE=['es_uaskirj']
IND_TYPE_TO_SAVE=['es_uaskirj']
OFFSET_TO_SUBTRACT=[0]

class Plot_Sum():
    
    def make_plots(self):
        self.exec_file=os.path.abspath(__file__)
        self.python_directory=os.path.split(self.exec_file)[0]

        



        self.datapath=path['sum_data_path']+'summary_data/'
        
        self.pck_open_name='vecstrength_unpaired_sumsundata.pck'
        self.datapath=path['input_data_path']+'summary_data/'

        self.sumdt=fh.open_pickle(self.datapath+self.pck_open_name)['sumdt']

        for data_type_ind,crflytype in enumerate(FLY_TYPE_TO_SAVE):
            motor_offset=OFFSET_TO_SUBTRACT[data_type_ind]
            self.crindtype=IND_TYPE_TO_SAVE[data_type_ind]

            self.csv_ind=0
            flyinds=self.sumdt[crflytype]
            for self.crflynum in flyinds:
                
                self.get_fly_data(self.crflynum)
                if self.crdt:
                    for crind in np.arange(len(self.crdt)):
                        #0 is index for first sun flight
                        #1 is index for second sun flight
                        #2 is index for stripe flight
                        self.add_to_cr_csv(crind=crind)



    def get_fly_data(self,crflynum):
        try:
            fname_to_plot=self.sumdt[crflynum][0]['sumstats'][0]['fname']
        except:
            self.crdt={}
            return
        split_name=fname_to_plot.split('Dropbox')
        
        datpath=dropbox_path+'Dropbox'+split_name[1].strip('.txt')+'_rawdata.pck'

        self.crdt=fh.open_pickle(datpath)
        
        print 'pickle file opened'

    def add_to_cr_csv(self, **kwargs):
        
        outfile = save_file_folder+self.crindtype +'_'+str(self.csv_ind)+'_.csv'
        outsize = 0.6*(1024 * 1024 * 1024) # 600MB
        file_exists = os.path.isfile(outfile)

        make_new_file=False
        if file_exists:
            if os.path.getsize(outfile)>outsize:
                if not paired_flag:
                    make_new_file=True
                else:
                    if not ctrval:
                        make_new_file=True

        if make_new_file:
            self.csv_ind += 1
            outfile = save_file_folder+self.crindtype +'_'+str(self.csv_ind)+'_.csv'
                #rewrite file_exists to False so that new header is written
            file_exists=False
        else:

            outfile = save_file_folder+self.crindtype +'_'+str(self.csv_ind)+'_.csv'
        #make matric to add
        csv_header,cr_data_matrix,nanline=self.prepare_data_matrix(**kwargs)

        with open(outfile, 'a+') as csvfile:


            writer=csv.writer(csvfile)
            if not file_exists:

                writer.writerow(csv_header)

            np.savetxt(csvfile,cr_data_matrix,delimiter=',')
            #write a nan line after first flight bout of paired flight
            

    def prepare_data_matrix(self,**kwargs):
        crind=kwargs['crind']
        ln=len(self.crdt[crind]['raw_time'])
            #name_text_array=np.repeat(self.crdt['fname'].split('data/')[1].split('/cloop')[0],ln)
        name_num_array=np.repeat(self.crflynum,ln)

        tmp_mot=self.crdt[crind]['mot_rad']
        final_mot=calc.force_angle_to_range(tmp_mot,force_range='0_2pi')
        lenvls=len(tmp_mot)
        stim_type_array=np.full((lenvls),crind)
        data_to_save=(name_num_array,self.crdt[crind]['raw_time'],self.crdt[crind]['lftwng'],self.crdt[crind]['rtwng'],final_mot,stim_type_array)

        header=['fly_id','raw_time','lftwng','rtwng','stimulus_in_radians','stim_type']
        nanline=np.repeat(np.nan,len(header))
        
        return header,np.column_stack(data_to_save),nanline





if __name__ == '__main__':

    plot_anal=Plot_Sum()
    dat=plot_anal.make_plots()
