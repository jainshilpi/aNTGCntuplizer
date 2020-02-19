#!/usr/bin/env python

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys

config = config()


#**************************submit function***********************
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
def submit(config):
	try:
		crabCommand('submit', config = config)
	except HTTPException as hte:
		print "Failed submitting task: %s" % (hte.headers)
	except ClientException as cle:
		print "Failed submitting task: %s" % (cle)
#****************************************************************


workarea='/afs/hep.wisc.edu/home/shilpi/private/forMETStudy/CMSSW_9_4_13/src/ggAnalysis/ggNtuplizer/test/'

mainOutputDir = '/store/user/shilpi/aNTGC/ggNtuplizerSkim_AOD/BHStudy/'


config.General.requestName = 'BH2p5TeV'
config.General.transferLogs = True
config.General.workArea = '%s' % workarea


config.Site.storageSite = 'T2_US_Wisconsin'
#config.Site.whitelist = [#whitelist]
#config.Site.blacklist = [#blacklist]


config.JobType.psetName  = 'run_MC2017AOD_94X.py'
config.JobType.pluginName  = 'Analysis'


config.Data.inputDBS = 'phys03'
config.Data.inputDataset = "/BeamHalo_2017/shilpi-BeamHalo_AODSIM_v1_lessMemory-8e0901deea2974d15dce1d0266e4c5f9/USER"
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '%s' % mainOutputDir
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 100
config.JobType.allowUndistributedCMSSW = True
#config.Data.ignoreLocality = True
#config.Data.totalUnits = #totalUnits
#config.Data.lumiMask = '#lumiMaskFile'
submit(config)
