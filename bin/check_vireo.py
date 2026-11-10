#!/usr/bin/env python3
from statistics import median
import argparse
import sys

def check_donors(path, min_cells):
    counts = {}
    fragment={}
    with open(path) as f:
        next(f)
        for line in f:
            donor=line.split("\t")[1]
            if (donor != "doublet") and (donor != "unassigned"):
                counts[donor]=counts.get(donor,0)+1
                counts["all"]=counts.get("all",0)+1
        for key in counts.keys():
            if key != "all":
                if min_cells>0:
                    if counts[key]<100:
                        fragment[key]="fail"
                    else:
                        fragment[key]="pass"
                else:
                    if counts[key]<(0.1*counts["all"]):
                        fragment[key]="fail"
                    else:
                        fragment[key]="pass"
    return fragment

def check_pairs(path, min_sites):
    gt_close={}
    usable = {}
    with open(path) as f:
        next(f)
        for line in f:
            d1= line.split("\t")[0]
            d2= line.split("\t")[1]
            donors=(d1, d2)
            sites = int(line.split("\t")[2])
            if sites > min_sites:
                frac = int(line.split("\t")[3])/sites
                usable[donors]=frac
    if len (usable) >=3:
        for key in usable.keys():
            disc=usable[key]
            if disc<0.45:
                other_disc=[]
                for key2 in usable.keys():
                    if key2 != key:
                        other_disc.append(usable[key2])
                if disc<0.6*median(other_disc):
                    gt_close[key[0]]=gt_close.get(key[0],0)+1
                    gt_close[key[1]]=gt_close.get(key[1],0)+1
            else:
                gt_close[key[0]]=gt_close.get(key[0],0)+0
                gt_close[key[1]]=gt_close.get(key[1],0)+0
    return gt_close

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("donor_ids", help="Vireo donor_ids.tsv")
    p.add_argument("pairwise", help="pairwise_gt_diff.tsv")
    p.add_argument("--min_cells", type=int, default=0,
        help="absolute minimum number of cells per donor (default: %(default)s = off)")
    p.add_argument("--min_sites", type=int, default=200, help="min sites for a usable pair (default: %(default)s)")

    ##p.add_argument("-o", "--output", default="donor_qc.tsv", help="output TSV (default: %(default)s)")
    # p.add_argument("--cell_fraction", type=float, default=0.1,
    #                help="fail donors with cells < this fraction of the mean (default: %(default)s)")
    # p.add_argument("--min_pairs", type=int, default=3, help="min usable pairs to run the pair check (default: %(default)s)")
    # p.add_argument("--rel", type=float, default=0.60, help="relative cutoff vs median of other pairs (default: %(default)s)")
    # p.add_argument("--abs_cut", type=float, default=0.45, help="absolute fraction cutoff (default: %(default)s)")
    a = p.parse_args()

    gt_close=check_pairs(a.pairwise, a.min_sites)
    fragment=check_donors(a.donor_ids, a.min_cells)
    check=0
    text=""
    for key in fragment.keys():
        if gt_close.get(key,0)>0 or fragment[key]=="fail":
            text=text+key+" fail\n"
            check+=1
        else:
            text=text+key+" pass\n"
    if check > 0:
        with open("failed_donors_check", "w", newline="") as out:
            out.write(text)
    else:
        with open("passed_donors_check", "w", newline="") as out:
            out.write(text)

if __name__ == "__main__":
    main()