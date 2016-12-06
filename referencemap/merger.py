#!/usr/bin/env python

##############################################################################
# Merges raw files to complete configuration ones (ex: all of IT73)
##############################################################################

import numpy as np
import healpy as hp
import glob
import os
import itertools
import argparse

import dataFunctions as df


def merger(params, nmaps=3, **opts):

    # Option for overwriting existing merged file
    # prefix = '/data/user/jbourbeau/anisotropy/maps/'
    outBase = '_'.join(params) + '.fits'
    outFile = opts['outdir'] + '/' + outBase
    if os.path.isfile(outFile):
        if not opts['overwrite']:
            return
        os.remove(outFile)

    # Collect all fits files with every parameter in params
    config = params[0]
    # * gets all possible parameter options
    masterList = glob.glob(opts['mapdir'] + '/' + config +
                           '/' + config + '_*_????-??-??.fits')
    masterList.sort()
    fileList = []
    for f in masterList:
        f_params = os.path.basename(f)[:-5].split('_')
        if f_params[:-1] == params:
            fileList.append(f)

    print('Reading files with {}'.format(params))
    print('{} files found...'.format(len(fileList)))
    if len(fileList) == 0:
        return

    # Merge files
    for i, file in enumerate(fileList):
        print('Working on {}'.format(file))
        temp = hp.read_map(file, range(nmaps), verbose=False)
        if i == 0:
            combined_map = [np.zeros(temp[j].shape) for j in range(nmaps)]
        for k in range(nmaps):
            # print temp[k]
            combined_map[k] += temp[k]

    # Write to file
    hp.write_map(outFile, combined_map)


def projectMerge(detector, nmaps=3):

    # Clean existing merged files
    prefix = '/data/user/jbourbeau/anisotropy/maps/merged'
    rmList = glob.glob('%s/%s_*.fits' % (prefix, detector))
    print rmList
    for f in rmList:
        print(f)
        os.remove(f)

    # Get fileList for given detector
    fileList = glob.glob(prefix + '/' + detector + '?*_*H_*.fits')
    fileList.sort()

    # Get list of parameters
    paramList = []
    for file in fileList:
        params = os.path.basename(file)[:-5].split('_')
        if params[1:] not in paramList:
            paramList += [params[1:]]       # ignore detector configuration

    # Merge maps with given parameters
    for params in paramList:

        tempList = []
        for file in fileList:
            f_params = os.path.basename(file)[:-5].split('_')
            if f_params[1:] == params:
                tempList.append(file)

        if len(tempList) == 1:
            continue

        print 'Working on', detector, params
        print len(tempList), 'files found...'

        for i, file in enumerate(tempList):
            temp = hp.read_map(file, range(nmaps), verbose=False)
            if i == 0:
                combined_map = [np.zeros(temp[j].shape) for j in range(nmaps)]
            for k in range(nmaps):
                combined_map[k] += temp[k]

        # Write to file (overwrite if necessary)
        outFile = '%s/%s_%s.fits' % (prefix, detector, '_'.join(params))
        # print combined_map.title
        hp.write_map(outFile, combined_map)


if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Creates merged map files')
    p.add_argument('-c', '--config', dest='config', nargs='*',
                   help='Detector configuration to be merged.')
    p.add_argument('--mapdir', dest='mapdir',
                   default='/data/user/jbourbeau/anisotropy/maps',
                   help='Directory containing maps to be merged.')
    p.add_argument('--outdir', dest='outdir',
                   default='/data/user/jbourbeau/anisotropy/maps/merged/',
                   help='Output directory where merged maps will be saved.')
    p.add_argument('--overwrite', dest='overwrite',
                   default=False, action='store_true',
                   help='Option to overwrite existing merged maps')
    args = p.parse_args()
    opts = vars(args).copy()

    if not args.outdir:
        opts['outdir'] = opts['mapdir']
    # Collect all raw map files
    #prefix = '/data/user/fmcnally/anisotropy/maps/raw/'
    #prefix = '/data/user/jbourbeau/anisotropy/maps/raw/'
    # prefix = '/data/user/jbourbeau/anisotropy/maps/'
    if not args.config:
        configs = df.getICconfigs()
    else:
        configs = args.config
    #configs += ['IT59','IT73','IT81','IT81-II','IT81-III','IT81-IV']
    masterList = []
    for config in configs:
        masterList += glob.glob(args.mapdir + '/' +
                                config + '/*_????-??-??.fits')
    masterList.sort()

    # Calculate comprehensive list of all unique parameters
    paramList = []
    # Removes '.fits' from the entries in masterList
    testList = [os.path.basename(f)[:-5] for f in masterList]
    for f in testList:
        # Need to be split for the next section of code--to get rid of
        # duplicate parameters
        f_params = f.split('_')
        f_params = f_params[:-1]        # exclude date
        paramList.append(f_params)
    paramList.sort()

    # Python black magic eliminates duplicates
    paramList = list(a for a, _ in itertools.groupby(paramList))

    # Run merger on each set of parameters
    for params in paramList:
        merger(params, **opts)

    # Run project merger
    detectorList = ['IT', 'IC']
    for detector in detectorList:
        projectMerge(detector)
