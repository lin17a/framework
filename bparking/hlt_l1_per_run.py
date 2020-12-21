#! /bin/env python
# python hlt_l1_per_run.py --runs run1,run2... --datasets dataset1,dataset2... <--keep/--ignore path1,path2... --output_name output>
# python hlt_l1_per_run.py --runs $(cat gold_json_hlt_2018 | grep -v '#' | awk '{print $1}' | tr '\n' ',') --datasets ParkingBPH1,ParkingBPH2,ParkingBPH3,ParkingBPH4,ParkingBPH5,ParkingBPH6

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
parser.add_argument('--runs', help = 'runs to consider, comma separated', required = True)
parser.add_argument('--datasets', help = 'datasets to consider, comma separated', required = True)
parser.add_argument('--output_name', help = '', default = 'hlt_l1.json', required = False)
parser.add_argument('--keep', help = 'keep the following HLT paths, comma separated. partial name keeps all matches', default = '', required = False)
parser.add_argument('--ignore', help = 'ignore the following HLT paths, comma separated. partial name ignores all matches', default = '', required = False)
args = parser.parse_args()

if not validate_environment():
        raise RuntimeError('Code relies on hltConfigFromDB, which needs to be ran from CMSSW > 10_6_2!!')

runs = args.runs.split(',')
datasets = args.datasets.split(',')
keep_only = args.keep.split(',')
ignore_only = args.ignore.split(',')

if (keep_only != [''] and ignore_only != ['']):
        raise RuntimeError("--keep and --ignore aren't meant to be used together!!")

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
                                skip_path = False
                                if (ignore_only != ['']):
                                        for ignore in ignore_only:
                                                skip_path = skip_path or (path.find(ignore) != -1)
                                        if (skip_path):
                                                continue

                                keep_path = False
                                if (keep_only != ['']):
                                        for keep in keep_only:
                                                keep_path = keep_path or (path.find(keep) != -1)
                                        if (not keep_path):
                                                continue

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
        syscall('rm hlttmp.pyc')

with open(args.output_name, 'w') as out:
        json.dump(collections.OrderedDict(sorted(hlt_l1.items())), out, indent = 2)

