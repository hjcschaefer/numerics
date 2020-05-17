#!/usr/bin/env python

import os
import sys

DATADIR="/home/hschaefe/code/python/getnss/data"
RESDIR="analysis"
CLUSTER="/home/hschaefe/code/kmeans/mm/cluster2"

if __name__ == '__main__':

    funds = [f for f in os.listdir(DATADIR) if os.path.isfile(os.path.join(DATADIR, f))]
    for fund in funds:
        print(fund)
        try:
            os.makedirs(os.path.join(RESDIR, fund))
        except FileExistsError as e:
            pass

        os.system("cd %s; %s 4 %s 20"%(os.path.join(RESDIR, fund), CLUSTER, os.path.join(DATADIR, fund)))




