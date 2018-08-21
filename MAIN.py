#!/usr/bin/python3

"""
Main secondary metabolism pipeline. Imports modules and processes data to write a dataframe which can later be processed with R.
"""

import os
import csv
import sys
import argparse
import pandas as pd
import shutil

from smModule.smServerSide import tmpSmBiTable, createBidirSmurf, mysqlSmChecker
from smModule.aspSMDl import dlSMdata
import smModule.bioSlim3 as bio # tupleToFasta
from smModule.processMibig3 import processMibig, dlSmurfProteins, writeMibigFormatted, processBlastResult, dlMibig

from smModule.misc import readConfig



parser=argparse.ArgumentParser(description='''
Script to execute the seconary metabolite analysis pipeline. In case you want to leave out some analysis you have to modify the script.\n

Example:

python3 MAIN.py -o nigri_orgs.txt -bibase biblast_table -biFinal smurf_bidir_hits_test -t species_tree.nwk -l run.log -od out_dir
''')
parser.add_argument("--orgs", "-o",
					dest="filename",
					required=True,
					help="Input file with jgi names (one per row) of organisms", metavar="FILE")
parser.add_argument("--treeFile", "-t",
					dest="tree",
					required=True,
					help="Tree file in newick format",
					metavar="FILE")
parser.add_argument("--biblastTable", "-bibase",
					dest="bibase",
					required=True,
					help="Specify the original biblast table to use as basis for the smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--clusterBiblast", "-biFinal",
					dest="biFinal",
					required=True,
					help="Specify name for smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--log", "-l",
					dest="logFile",
					required=True,
					help="Specify the name for log file", metavar="CHAR")
parser.add_argument("--outdir", "-od",
					dest="sn",
					required=True,
					help="Output directory, preferably the name of your set", metavar="CHAR")
parser.add_argument("--otherBlast", "-ob",
					dest="otherBlast",
					action="store_true",
					default=False,
					help="If problems are encountered with blast use this command for another query")
parser.add_argument("--fetch_raw_data", "-fr",
					dest="fetch_raw_data",
					action="store_true",
					default=False,
					help="This command uses non subsetted tables")
parser.add_argument("--cluster_and_visualize", "-cv",
					dest="cluster_and_visualize",
					action="store_true",
					default=False,
					help="If you only want to rerun clustering and visualization, choose this option")
args=parser.parse_args()




with open("configNew.txt") as c:
	config = readConfig(c.readlines())

print("Queries to format smurf databases")
if args.otherBlast:
	print("%s -in smurf.fasta -dbtype prot" % config['makeblastdbPath'])
else:
	print("%s -i smurf.fasta" % config['formatDbPath'])
	


########
# LOADING ORGS

filename = args.filename 

biblastBaseTable = args.bibase
smurfBidirHitsName = args.biFinal 
testLogName = args.logFile 
treeFile = args.tree 
setName = args.sn 

with open(filename, "r") as tmp:
	orgSet = [item.strip() for item in tmp.readlines()]

if setName not in os.listdir():
	os.mkdir(setName)

else:
	input("A folder for the specified set is already available. If you want to rerun the analysis press Enter, else ctrl+c/d\n")

shutil.copy2(treeFile, setName)

def download_and_processing():
	os.chdir(setName)
	# Checking data
	mysqlSmChecker(orgSet, testLogName)

	# Creating smurf bidir hits table for dataset.
	tmpSmBiTable(smtable = biblastBaseTable) # Creating a temporary table


	# THIS needs to be written to disk
	createBidirSmurf(smtable = smurfBidirHitsName) # Creating a sm cluster table
	# 

	# DOWNLOADING SM DATA
	input("Did you create a smurfbiblast table? If yes, press Enter, if not, ctrl+c/d")


	smIpGf = dlSMdata(orgSet)

	# FILE WRITING HERE
	with open("data.tsv", "w") as handle:
		tsv_out = csv.writer(handle, delimiter="\t")
		tsv_out.writerows(smIpGf)
	# FILE WRITING END

	##########
	# MIBIG PART

	input("Press Enter to continue with dereplication of mibig entries in your dataset")

	print("Creating directory for mibig")

	if "mibig" not in os.listdir():
		os.mkdir("mibig")

	os.chdir("mibig")
	dlMibig()

	# this also needs to be written to file

	sp = dlSmurfProteins(orgSet)
	bio.tupleToFasta(sp,"smurf.fasta")
	print("Wrote smurf proteins to disk")

	#  File writing end

	try:
		print("Formatting smurf database")
		if args.otherblast:
			os.system("%s -i smurf.fasta" % config['formatDbPath'])
		else:
			os.system("%s -in smurf.fasta -dbtype prot" % config['makeblastdbPath'])
	except Exception as e:
		print("Could not format smurf proteins to blast database")

	print(os.getcwd())
	print("processing mibig data")
	geneClusters = processMibig("mibigGbkFiles")
	print("Writing mibig data to disk")
	writeMibigFormatted(geneClusters)


	input("""Press Enter to perform blast: \n
	%s -query mibigDb.fasta -subject smurf.fasta -max_target_seqs 25 -evalue 1e-50 -num_threads 3 -outfmt '6 std qlen slen' -out mibigVsSmurf.txt\n WILL TAKE APPROX 30MIN TO 1H.""" % config['blastpPath'])
	# os.getcwd()
	print("performing blast")
	os.system("%s -query mibigDb.fasta -subject smurf.fasta -max_target_seqs 25 -evalue 1e-50 -num_threads 3 -outfmt '6 std qlen slen' -out mibigVsSmurf.txt" % config['blastpPath'])
	print("done")

	print("Processing mibig blast results")
	mibigBlast = processBlastResult(blastFile = "mibigVsSmurf.txt", translationFile = "translateBgc.txt")

	mibigBlast['org_id'] = mibigBlast['org_id'].astype('int64')

	mibigBlast['protein_id'] = mibigBlast['protein_id'].astype('int64')

	os.chdir("..")


	try:
	    print("Joining mibig annotation on sm data")
	    smData = pd.merge(smIpGf, mibigBlast[['org_id', 'protein_id', 'compound']], how = "left", left_on=['org_id','protein_id'], right_on = ['org_id','protein_id'])
	except Exception as e:
	    raise
	    print("Cannot merge mibig results and sm data")

	smData.to_csv("sm_data_"+setName+".tsv", sep = '\t', encoding = "UTF-8", index = False, na_rep='none')


	# Using r scipt to create cluster families


def clustering_and_output():
	print("Current working directory")
	print(os.getcwd())

	cluster_blast_all =bio.dbFetch("""SELECT *, ROUND(
                 ( COALESCE( ( pident_tailoring /(clust_size-q_max_bb)  )*0.35,0)+
                 COALESCE( (pident_bb/q_max_bb)*0.65, 0)
                 ),2) AS pident_score FROM (

                 SELECT bidir.q_org, bidir.q_clust_id, bidir.h_org, bidir.h_clust_id,
                 SUM(CASE WHEN sm_short != 'none' then pident else 0 end) AS pident_bb,
                 smurf.clust_size,
                 COALESCE(SUM(sm_short != 'none'),0) AS count_bb,
                 tqmax.q_max_bb,
                 SUM(CASE WHEN sm_short = 'none' then pident else 0 end) AS pident_tailoring,
                 COALESCE(SUM(sm_short = 'none'),0) AS count_tailoring

                 FROM (
                 SELECT bia.* FROM (
                 SELECT * FROM %s) bia
                 JOIN (
                 SELECT * FROM %s ) bib
                 ON bia.h_clust_id = bib.q_clust_id
                 WHERE bia.q_clust_id = bib.h_clust_id
                 GROUP BY bia.q_clust_id, bia.q_protein_id, bia.h_clust_id, bia.h_protein_id
                 ) AS bidir

                 JOIN smurf
                 ON bidir.q_org = smurf.org_id AND bidir.q_protein_id = smurf.sm_protein_id AND bidir.q_clust_id != bidir.h_clust_id

                 JOIN (SELECT CONCAT(org_id, '_' , clust_backbone,'_', clust_size) AS q_clust_id, SUM(sm_short != 'none') AS q_max_bb FROM smurf
                 GROUP BY q_clust_id) tqmax
                 ON bidir.q_clust_id = tqmax.q_clust_id

                 GROUP BY q_clust_id, h_clust_id ) ta;""" %  (smurfBidirHitsName, smurfBidirHitsName))
	with open("clusterBlastAll.csv", "w") as handle:
		cluster_writer = csv.writer(handle, delimiter = ",")
		cluster_writer.writerow(["q_org", "q_clust_id", "h_org", "h_clust_id", "pident_bb", "clust_size", "count_bb", "q_max_bb", "pident_tailoring", "count_tailoring", "pident_score"])
		cluster_writer.writerows(cluster_blast_all)
	smFile = "sm_data_"+setName+".tsv"


	os.chdir("..")

	cmdString = "Rscript clusterData.R %s %s" % (smFile, setName) # File was clustered before, so it will be renamed with a _c extension
	print(cmdString)
	print("Calculating cluster families")

	os.system(cmdString)


	print("Searching for unique secondary metabolic gene clusters in tree")
	cmdString = "Rscript uniquesAtNodes.R %s %s %s" % (smFile.replace(".tsv","")+"_c.tsv", setName, treeFile)

	os.system(cmdString)

	cmdString = "Rscript visualizeSm.R %s %s %s" % (smFile.replace(".tsv","")+"_c_n.tsv", setName, treeFile)
	print("Preparing visualization")
	os.system(cmdString)

if __name__ == '__main__' and not args.cluster_and_visualize:
	download_and_processing()
	clustering_and_output()

if __name__ == '__main__' and args.cluster_and_visualize:
	clustering_and_output()
