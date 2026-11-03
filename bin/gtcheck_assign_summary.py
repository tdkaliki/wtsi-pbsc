#!/usr/bin/env python3

## make donor assignment across panels

import sys
import csv
import re
import statistics
import pandas as pd

VERBOSE = True


nargs = len(sys.argv)
if nargs < 3:
    sys.exit("usage: {:s} <summary output> <assignment output panel 1> [<assignment output panel 2> ....]"
        .format(sys.argv[0]))

oufn = sys.argv[1] 
ZSCORE_THRESH = float(sys.argv[2]) # Default is 8 - cardinal has been run with this threshold
ZSCORE_DIST_THRESH = float(sys.argv[3]) # Default is 8 - cardinal has been run with this threshold
infns = sys.argv[4:]

sys.stdout.write("compare scores from {:d} files and write overall assignment to file {:s}\n".
    format(len(infns), oufn))


SMVAL = float(1e-9)
# ZSCORE_THRESH = float(3)
# ZSCORE_DIST_THRESH = float(3)

DONOR_CELLINE_MATCHSTR = re.compile(".*celline_.*")
FNAM_MATCHSTR = re.compile("^pool_(\S+)_panel_(\S+)_gtcheck_donor_assignments")

class AssignmentTables:
    def __init__(self):
        self.donors = []
        self.cell_lines = {}
        self.panels = {}
        self.cell_line_panel = []

    def parse_assignment(self, infn):
        tabd = {}
        ctr = 0
        fm = FNAM_MATCHSTR.match(infn)
        if not fm:
            sys.exit("ERROR: unexpected input file name '{:s}'.".format(infn))
        panel = fm.group(2)
        infh = open(infn, 'r')
        for row in csv.DictReader(infh):
            ctr += 1
            donor_id = row['donor_query']
            if donor_id not in self.donors:
                self.donors.append(donor_id)
                # donor_query,donor_gt,score0,score1,score_n,n,mean,sd,z0,z1
            tabd[row['donor_query']] = (row['donor_gt'], row['score0'], row['score1'], row['z0'], row['z1'], row['score_n'],row['n'],row['mean'],row['sd'])
        infh.close()
        self.panels[panel] = tabd
        return ctr

    def parse_assignment_output_files(self, infns = []):
        ctr = 0
        for fn in infns:
            ctr += 1
            if VERBOSE:
                sys.stderr.write("# input file {:d}: {:s}\n".format(ctr, fn))
            self.parse_assignment(fn)
        spike_in_panel = self.identify_cell_line_panel()
        if VERBOSE:
            try:
                sys.stderr.write("# cell line panel is '{:s}'\n".format(spike_in_panel))
            except:
                _ = 'spikeins_dont exist'
        return ctr

    def identify_cell_line_panel(self):
        fnam_cell_line = None
        for nam in self.panels:
            panel = self.panels[nam]
            is_cell_line = None
            is_error = False
            is_ambiguous= False
            z_dist = {}
            for donor in panel:
                dtup = panel[donor]
                assigned_donor = dtup[0]
                score = float(dtup[1])
                score1 = float(dtup[2])
                # Its not only celline panels bu in general if there is only a binary donor panel.
                z1 = float(dtup[3])
                z0 = float(dtup[4])
                z_dist[donor]={'z1':z1,'z0':z0}
            z_dist = pd.DataFrame(z_dist).T
            if (all(z_dist['z1'].isin([0, 1])) or all(z_dist['z0'].isin([0, 1]))) and all((z_dist['z1']-z_dist['z0']).isin([0, 1])):
                
                fnam_cell_line = nam
                self.cell_line_panel.append(fnam_cell_line)
                for donor in panel:
                    dtup = panel[donor]
                    assigned_donor = dtup[0]
                    cell_line = assigned_donor
                    score = float(dtup[1])
                    score1 = float(dtup[2])
                    
                    if fnam_cell_line not in self.cell_lines:
                        self.cell_lines[fnam_cell_line]={}
                    if cell_line in self.cell_lines[fnam_cell_line]:
                        self.cell_lines[fnam_cell_line][cell_line][donor] = (score, score1)
                    else:
                        self.cell_lines[fnam_cell_line][cell_line] = {donor: (score, score1)}
            else:
                _='not a celline'
                # cm = DONOR_CELLINE_MATCHSTR.match(assigned_donor)
                # if ((z1 ==0 or z0==0) and z1-z0==0 and not is_cell_line==None):
                #     continue
                # if (((z1 ==1 or z0==1) and z1-z0==1) or ((z1 ==0 or z0==0) and z1-z0==0)):
                        
                #     if is_cell_line is None:
                #         is_cell_line = True
                #         # if fnam_cell_line is not None and fnam_cell_line != nam:
                #         #     sys.exit("ERROR: there are multiple panels with cell line names in them, e.g. {:s} and {:s}."
                #         #         .format(fnam_cell_line, fnam))
                #         fnam_cell_line = nam
                #         self.cell_line_panel.append(fnam_cell_line)
                #     else:
                #         is_error = not is_cell_line
                #     cell_line = assigned_donor
                #     if fnam_cell_line not in self.cell_lines:
                #         self.cell_lines[fnam_cell_line]={}
                #     if cell_line in self.cell_lines[fnam_cell_line]:
                #         self.cell_lines[fnam_cell_line][cell_line][donor] = (score, score1)
                #     else:
                #         self.cell_lines[fnam_cell_line][cell_line] = {donor: (score, score1)}
                # elif is_cell_line is None:
                #     is_cell_line = False
                # else:
                #     is_error = is_cell_line
                # if is_error:
                    
                #     sys.exit(f"ERROR: panel {panel} has cell lines mixed with non-cell line labels.")
                    # in this particular situation i would like to select it as not being a celline.
        
        return fnam_cell_line

    def make_cell_line_assignment(self, donor,panel):
        # calculate mean of scores across donors with 1 donor removed
        assignment = None # dictionary of cell lines, values are true/false for a confiden assignment outcome
        n_confident = 0
        for cl in self.cell_lines[panel]:
            print(cl)
            scores = []
            donor_score = 0
            donor_score1 = 0
            scd = self.cell_lines[panel][cl]
            is_present = False
            for d in self.donors:
                if d not in scd:
                    continue
                sc, sc1 = scd[d]
                if d == donor:
                    donor_score = float(sc)
                    donor_score1 = float(sc1)
                    is_present = True
                else:
                    scores.append(float(sc))

            ns = len(scores)
            if ns > 0:
                m = statistics.mean(scores)
                if ns > 1:
                    sd = statistics.stdev(scores)
                    if (sd==0):
                        sd=0.0000001
                    z = (m - donor_score)/sd
                    is_confident_assignment = z > ZSCORE_THRESH
                    if donor == 'donor8':
                        print("donor8", donor_score, z, sd, m, is_confident_assignment)
                else:
                    is_confident_assignment = 2*donor_score < m
            else:
                is_confident_assignment = donor_score*3 < donor_score1

            if is_confident_assignment and donor_score > SMVAL:
                n_confident += 1
                assignment = cl
        # check that there is only one confident assingment
        if donor == 'donor8':
            print("cell line assignment donor8", donor_score, donor_score1,is_confident_assignment, n_confident)

        if n_confident != 1:
            assignment = None
            z=None
        return z,assignment

    def make_panel_assignment(self, panel, donor):
        donor_assigned, score0, score1, z0, z1 ,score_n, n, mean, sd= self.panels[panel][donor]
        is_confident_assignment = float(z0) > ZSCORE_THRESH and float(z0)-float(z1) > ZSCORE_DIST_THRESH
        # is_confident_assignment = True
        if is_confident_assignment:
            assignment = donor_assigned
        else:
            assignment = None
        return z0,z1,assignment

    def assign_donors(self, oufh):
        oufh.write("donor_query,donor_gt,panel\n")
        df ={}
        if (len(self.cell_line_panel) ==0):
            self.identify_cell_line_panel()
        for donor in self.donors:
            final_assignment = None
            panel_assignment = None
            assigned_donors=[]
            assigned_cellines=[]
            cell_line_assignment = None
            panel_ass_panel = None
            cell_line_panel = None
            final_panel = None
            n_assignments = 0
            for panel in self.panels:
                if panel in self.cell_line_panel:
                    z,cell_line_assignment = self.make_cell_line_assignment(donor,panel)
                    if cell_line_assignment is not None:
                        cell_line_panel = cell_line_assignment
                        print("donor:", donor,"panel:", panel, "cell_line_assignment", cell_line_assignment)
                        assigned_cellines.append({'panel':panel,'z':z,'assignment':cell_line_assignment})
                else:
                    z0,z1,assignment = self.make_panel_assignment(panel, donor)
                    print("donor:", donor, ",panel:", panel, "assignment:", assignment)
                    if assignment is not None:
                        panel_assignment = assignment
                        panel_ass_panel = panel
                        n_assignments += 1
                        assigned_donors.append({'panel':panel_ass_panel,'assignment':panel_assignment,'z0':z0,'z1':z1})
            if n_assignments == 1:
                final_assignment = panel_assignment
                final_panel = panel_ass_panel
            elif (n_assignments >1):
                # Here we pick the best of 2 matches
                assignments = pd.DataFrame(assigned_donors)
                assignments['z_diff']= assignments['z0'].astype(float)-assignments['z1'].astype(float)
                final_assignment =assignments[assignments['z_diff']==max(assignments['z_diff'])]['assignment'].values[0]
                final_panel=assignments[assignments['z_diff']==max(assignments['z_diff'])]['panel'].values[0]
            elif n_assignments == 0:
                DF = pd.DataFrame(assigned_cellines)
                if (len(DF))>0:
                    Best_Match = DF[DF['z'].astype(float)==max(DF['z'].astype(float))]
                    final_assignment = Best_Match['assignment'].values[0]
                    final_panel = Best_Match['panel'].values[0]
                else:
                    final_assignment=cell_line_assignment
                    final_panel=None
            if final_assignment is None:
                final_assignment = 'NONE'
                final_panel = 'NONE'
            if final_panel is None:
                final_panel = 'NONE'
            oustr = "{:s},{:s},{:s}".format(donor, final_assignment, final_panel)

            if(final_panel=='NONE'):
                print('panel is none')
                df[donor]={'donor_query':donor,'donor_gt':'NONE','score0':0,'score1':0,'score_n':0,'n':0,'mean':0,'sd':0,'z0':0,'z1':0,'final_panel':final_panel}
            else:
                # df[donor]={'final_assignment':final_assignment,'final_panel':final_panel}
                # 0 row['donor_gt'], 1 row['score0'], 2 row['score1'], 3 row['z0'], 4 row['z1'], 5 row['score_n'], 6 row['n'],7 row['mean'],8 row['sd']
                df[donor]={'donor_query':donor,'donor_gt':self.panels[final_panel][donor][0],
                'score0':self.panels[final_panel][donor][1],
                'score1':self.panels[final_panel][donor][2],
                'score_n':self.panels[final_panel][donor][5],
                'n':self.panels[final_panel][donor][6],
                'mean':self.panels[final_panel][donor][7],
                'sd':self.panels[final_panel][donor][8],
                'z0':self.panels[final_panel][donor][3],
                'z1':self.panels[final_panel][donor][4],
                'final_panel':final_panel}


            oufh.write(oustr + '\n')
            print(oustr)
        # build from records with explicit columns: pd.DataFrame(dict_of_dicts).T sorts columns alphabetically when there is a single donor
        columns = ['donor_query','donor_gt','score0','score1','score_n','n','mean','sd','z0','z1','final_panel']
        return pd.DataFrame(list(df.values()), columns=columns)

if __name__ == '__main__':

    asstab = AssignmentTables()
    fctr = asstab.parse_assignment_output_files(infns)
    if VERBOSE:
        sys.stderr.write("# {:d} assignment files read.\n".format(fctr))
    print(asstab.cell_lines)
    oufh = open(oufn, 'w')
    donor_panels_and_stats = asstab.assign_donors(oufh)
    donor_panels_and_stats.to_csv(f'stats_{oufn}',sep=',',index=False)
    oufh.close()

    sys.exit()
