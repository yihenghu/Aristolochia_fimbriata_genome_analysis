#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    -------------------------------------------------
    File Name:       get_mostly_single_copy_orthogroup.py
    Description:     This script is for select mostly single copy orthogroup from OrthoMCL
    Author:          Yiheng Hu
    Date:            2020/3/10
    Email:           yihenghu@yeah.net
    -------------------------------------------------
    Change Activity: 2020/3/10
    -------------------------------------------------
"""
__author__ = 'Yiheng Hu'

import sys
import re
import collections

# read orthogroup of orthoMCL
def readGroups(groups):
    orthogroups = {}
    with open(groups, "r") as OGs:
        for og in OGs:
            if not og.startswith("\t"):
                orthogroups0 = re.split(r'\s*', re.sub(',', ' ', og.rstrip()))
                orthogroups[orthogroups0[0]] = orthogroups0[1:]
    return orthogroups  # {OG1:[geneA1,geneB1,geneC1],OG2:[geneA2,geneB2]}

# stat orthogroup
def statGroups(orthogroups):
    orthogroups_stat = {}
    for og_id in orthogroups.keys():
        species = re.findall(r"(\w+)\|", "\t".join(orthogroups[og_id]))
        species_copy = dict(collections.Counter(species))
        orthogroups_stat[og_id] = species_copy
    return orthogroups_stat

# multiple of list
def mult(list0):
    result = 1
    for i in list0:
        result *= i
    return result

# stat monocots, eudicots, magnoliids, basal angiosperm per OGs cotain number of genes
division = {
    "Smo": "Out",
    "Pab": "Out",
    "Gbi": "Out",
    "Atr": "Out",
    "Nad": "Out",
    "Spo": "Mon",
    "Peq": "Mon",
    "Mac": "Mon",
    "Sbi": "Mon",
    "Osa": "Mon",
    "Afi": "Mog",
    "Cmi": "Mog",
    "Lir": "Mog",
    "Pni": "Mog",
    "Pam": "Mog",
    "Aco": "Eud",
    "Pso": "Eud",
    "Nnu": "Eud",
    "Vvi": "Eud",
    "Sly": "Eud",
    "Ptr": "Eud",
    "Ath": "Eud"}


def main():
    orthogroups = readGroups(sys.argv[1])
    orthogroups_stat = statGroups(orthogroups)
    new_orthogroups_stat = {}
    for og,stats in orthogroups_stat.items():
        new_stats = {}
        for sp in stats:
            new_stats[(sp, division[sp])] = stats[sp]
        new_orthogroups_stat[og] = new_stats
    for og, stats in new_orthogroups_stat.items():
        out, mon, mog, eud = 0, 0, 0, 0
        out0, mon0, mog0, eud0 = 0, 0, 0, 0
        for sp_di in stats:
            if sp_di[1] == "Out":
                if stats[sp_di] == 1:
                    out += 1
                else:
                    out0 += 1
            elif sp_di[1] == "Mon":
                if stats[sp_di] == 1:
                    mon += 1
                else:
                    mon0 += 1
            elif sp_di[1] == "Mog":
                if sp_di[0] == "Pni":
                    if 0 < stats[sp_di] <= 2:
                        mog += 1
                    else:
                        mog0 += 1
                else:
                    if stats[sp_di] == 1:
                        mog += 1
                    else:
                        mog0 += 1
            elif sp_di[1] == "Eud":
                if sp_di[0] == "Pso":
                    if 0 < stats[sp_di] <= 2:
                        eud += 1
                    else:
                        eud0 += 1
                else:
                    if stats[sp_di] == 1:
                        eud += 1
                    else:
                        eud0 += 1
        if out0 == 0 and mon0 == 0 and mog0 == 0 and eud0 == 0:
            if out >= 3 and mon >= 3 and mog >= 3 and eud >= 4:
                print og

if __name__ == "__main__":
    main()

