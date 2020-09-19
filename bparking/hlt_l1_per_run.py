#! /bin/env python

import os
import json
import collections
from argparse import ArgumentParser
import FWCore.ParameterSet.Config as cms

def syscall(cmd):
	print 'Executing: %s' % cmd
	retval = os.system(cmd)
	if retval != 0:
		raise RuntimeError('Command failed!')

def get_cmssw_version():
        cmssw_version = os.getenv("CMSSW_VERSION")
        if cmssw_version is not None:
                cmssw_version = cmssw_version.split('_')
                cmssw_version = cmssw_version[1:]
                cmssw_version = [int(ver) for ver in cmssw_version]

        return tuple(cmssw_version)

def validate_environment():
        host_name = os.getenv("HOSTNAME")
        if "lxplus" not in host_name:
                raise RuntimeError('Code relies on hltConfigFromDB, which needs to be ran at lxplus!')
        min_version = (10, 6, 2)
        for v1, v2 in zip(get_cmssw_version(), min_version):
                if v1 > v2:
                        return True
                if v1 < v2:
                        return False

parser = ArgumentParser()
parser.add_argument('runs')
parser.add_argument('datasets')
args = parser.parse_args()

if not validate_environment():
        raise RuntimeError('Code relies on hltConfigFromDB, which needs to be ran from CMSSW > 10_6_2!!')

runs = args.runs.split(',')
datasets = args.datasets.split(',')

hlt_l1 = {}
for run in runs:
        print run
        syscall('hltConfigFromDB --v2 --runNumber {} > hlttmp.py'.format(run))
        import hlttmp
        reload(hlttmp)
        from hlttmp import process

        for dataset in datasets:
                if hasattr(process.datasets, dataset):
                        paths = getattr(process.datasets, dataset).value()

                        for path in paths:
                                modules = getattr(process, path).moduleNames()
                                for module in modules:
                                        module_attr = getattr(process, module)
                                        if module_attr.type_() is 'HLTL1TSeed':
                                                seeds = module_attr.L1SeedsLogicalExpression.configValue()
                                                seeds = seeds.replace("'", '')
                                                seeds = seeds.split(' OR ')
                                                if not path in hlt_l1:
                                                        hlt_l1[path] = seeds
                                                else:
                                                        if set(seeds) != set(hlt_l1[path]):
                                                                print 'Stored seed list:'
                                                                print hlt_l1[path]
                                                                print 'New seed list:'
                                                                print seeds
                                                                raise RuntimeError('Path {} has conflicting L1 seed lists!'.format(path))
        syscall('rm hlttmp.py')

with open('hlt_l1.json', 'w') as out:
        json.dump(collections.OrderedDict(sorted(hlt_l1.items())), out, indent = 2)


# python hlt_l1_per_run.py $(cat gold_json_hlt_2018 | grep -v '#' | awk '{print $1}' | tr '\n' ',') ParkingBPH1,ParkingBPH2,ParkingBPH3,ParkingBPH4,ParkingBPH5,ParkingBPH6
