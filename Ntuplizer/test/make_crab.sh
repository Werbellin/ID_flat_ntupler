#!/bin/sh

### C. Ochando / 2013
#
#
# Usage : dbs.sh <list> 
#

# ========================================================================
# Configuration Options
# ========================================================================
list=$1

#prod="PAT_CMG_V5_15_0"
prod="phys14fix"

mkdir /data_CMS/cms/ochando/CrabJobs/${prod}

echo "---> Working with list: " $list

# process
#all_process=`echo $list | awk -F/ '{print $7}' | awk -F. '{print $1}' | awk -F_ '{print $2}'` #from folder DBS
all_process=`echo $list | awk -F/ '{print $2}' | awk -F. '{print $1}' | awk -F_ '{print $2}'`
echo "---> Working with process: " $all_process

# ======================
# Create multicrab.cfg
# ======================

multicrab_file=multicrab_${all_process}_${prod}.cfg

echo "---> Building " $multicrab_file

echo "[MULTICRAB]" >> ${multicrab_file}
echo "cfg=crab_LLR.cfg" >> ${multicrab_file}
echo " " >>  ${multicrab_file} 
echo "[COMMON]"  >> ${multicrab_file}
echo " " >>  ${multicrab_file} 


# ======================
#  Loop over datasets
# ======================
for datasetpath in `less $list`;do 

process=`echo $datasetpath | awk -F/ '{print $2}' | sed "s?_??g" | sed "s?-??g"`
PUprofile=`echo $datasetpath | awk -F/ '{print $3}' | awk -F- '{print $2}' | sed "s?_??g" | sed "s?-??g"`

datasample="DATA"
isMCflag=0
ispythia6flag=1

isDATA=`echo $datasetpath | grep Double`
if [ $isDATA = ""]
    then
    datasample="MC"
    isMCflag=1
fi

isPYTHIA6=`echo $datasetpath | grep pythia6`
if [ $isPYTHIA6 = ""]
    then
    ispythia6flag=0
fi

echo "[${process}${PUprofile}]" >>  ${multicrab_file} 
echo "CMSSW.datasetpath=${datasetpath}" >>  ${multicrab_file} 
echo "CMSSW.pycfg_params= isMC=${isMCflag} ispythia6=${ispythia6flag}" >>  ${multicrab_file} 
echo "USER.user_remote_dir=/eID/${datasample}/${prod}/${datasetpath}" >>  ${multicrab_file} 
echo "USER.ui_working_dir=/data_CMS/cms/ochando/CrabJobs/${prod}/${process}${PUprofile}" >> ${multicrab_file} 

echo " " >>  ${multicrab_file} 

done
