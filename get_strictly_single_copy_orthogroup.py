#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    -------------------------------------------------
    File Name:       get_strictly_single_copy_orthogroup.py
    Description:     This script is for select strictly single copy orthogroup form OrthoMCL
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

# read orthogroup of orthoMCL or orthoFinder
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


def main():
    orthogroups = readGroups(sys.argv[1])
    orthogroups_stat = statGroups(orthogroups)

    # select "Pso==2，Pni==2", and other species are single copy
    for og, stats in orthogroups_stat.items():
        if len(stats) == 22:
            if stats["Pso"] <= 2 and stats["Pni"] <= 2:
                sp_num = {}
                for sp, number in stats.items():
                    if sp == "Pso" or sp == "Pni":
                        pass
                    else:
                        if number == 1:
                            sp_num[sp] = number
                if sum(sp_num.values()) == 20 and mult(sp_num.values())  == 1:
                    print og


if __name__ == "__main__":
    main()
