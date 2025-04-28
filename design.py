# !/usr/bin/env python
# coding: utf-8

import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import util as uti
#from sklearn.ensemble import RandomForestClassifier



import subprocess
import os
working_dir = './'
indivfilename="individual_list.txt"

#$ FAUST 
#foldx='/home/pkekeneshuskey/sources/foldxLinux64/foldx_20251231'
foldx='/Users/alecloftus/Desktop/parvalbumin_foldx-main/foldx_20251231'

def run_foldx(command, input_files=None, output_prefix="foldx_output", working_dir=None):
    """
    Runs FoldX command via subprocess.

    Args:
        command (str): The FoldX command string (e.g., "RepairPDB", "BuildModel", etc.).
        input_files (list, optional): List of input file paths. Defaults to None.
        output_prefix (str, optional): Prefix for output files. Defaults to "foldx_output".
        working_dir (str, optional): The directory to run FoldX in. Defaults to the current working directory.

    Returns:
        tuple: A tuple containing the standard output and standard error of FoldX.
               Returns (None, None) if an error occurs.
    """

    if working_dir is None:
        working_dir = os.getcwd()

    if input_files:
        input_str = " ".join(input_files)  # Concatenate input files into a space-separated string
        full_command = f"{foldx} --command={command} --pdb={input_str} --output-file={output_prefix}"
    else:
        full_command = f"{foldx} --command={command} --output-file={output_prefix}"

    try:
        print(full_command)
        process = subprocess.Popen(
            full_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=working_dir,
            text=True  # Important for getting text output, not bytes
        )
        stdout, stderr = process.communicate()
        return stdout, stderr

    except FileNotFoundError:
        print("FoldX executable not found. Make sure it's in your PATH.")
        return None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None

#action='getpdb'  # get pdb structure 
#action='repair'  # fix it via foldx 
#action='ptraj'    # find residues at interface 
#action='build'   # introduce mutations via foldx [see mutList below]
#action='process' # process scores from foldx:w
pdbCode ='2NLN' #'1RWY'
working_directory = "./" # or wherever your files are.

# selecting interface residues 
chain="A"
residueAB= (8,33) # interface 1
residueCD= (50,89) #modified interface 2
#residueCD= (60,89)# interface 2
distance = 4.0   

# use the blast notebook for this             
mutList = [
 'E81K',
 'E81R',
 'S78K',
 """
 'D45K',
 'D59E',
 'G60E',
 'G98D',
 'K69G',
 'L58I'
 """
]
# obtained from running code with action ptraj
pairsList = ["D20-E81","P21-E81","P21-S80","P21-T78","P21-R75","D22-S80","D22-E81","F24-E81","F24-S84","E25-E81","E25-R75","E25-F66","P26-L85","P26-F66","P26-E81","P26-S84"]

#1RWY pairslist
#pairsList = ["I11-F70","I15-L67","F18-R75","T19-R75","A20-R75","A21-R75","D22-S78","D22-E81","D22-K80","S23-E81","S23-R75","F24-R75","F24-E81","F24-L67","F24-L85","H26-L85","H26-A88","F29-F70","F29-L67","V33-F70"]
#1RWY_mutlist = [
# 'S23K,R75D',
#'S78K',
# 'S23K',
# 'F29W',
# 'V33F',
# 'I15T L67T'
#"""
#]
rotabase_file = "rotabase.txt" #Make sure rotabase.txt is in the working directory or provide full path.

fasta="""\
>NP_037127.1 oncomodulin [Rattus norvegicus]
MSITDILSAEDIAAALQECQDPDTFEPQKFFQTSGLSKMSASQVKDIFRFIDNDQSGYLDGDELKYFLQKFQSDARELTE
SETKSLMDAADNDGDGKIGADEFQEMVHS"""

#1RWY model
#>XP_006241991.1 parvalbumin alpha isoform X1 [Rattus norvegicus]
#MSMTDLLSAEDIKKAIGAFTAADSFDHKKFFQMVGLKKKSADDVKKVFHILDKDKSGFIEEDELGSILKG
#FSSDARDLSAKETKTLMAAGDKDGDGKIGVEEFSTLVAES"""

def doit(
  action=None,
  dryRun=False) :         
  """
  FoldX calls
  """
  print(f"Performing {action}")
  
  ## processes 
  if 'getpdb' in action: 
    # Get file 
    print(f"Looking for {pdbCode}") 
    cmd = f"wget https://files.rcsb.org/view/{pdbCode}.pdb"
    process = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=working_dir,
                text=True  # Important for getting text output, not bytes
            )
    stdout, stderr = process.communicate()
    print(stderr[:-1])

    if 'Bad Request' in stderr:
        raise RuntimeError()
    elif 'Saving to' in stderr:
        print(pdbCode + '.pdb has been saved.')
    else:
        raise RuntimeError('Error not recognized. Please review PDB code and try again.')
    #repeating this task creates another pdb file saved as '1RWY.pdb.1'
    #   and so forth. Consider rewriting to overwrite prev versions?
    
  if 'repair' in action: 
    # Fix file
    # - takes a few minutes ; add progress bar ?
    pdb_file = f"{pdbCode}.pdb"
    output_prefix = "repaired"
    print("Might have to run manually....") 
    stdout, stderr = run_foldx("RepairPDB", input_files=[pdb_file], output_prefix=output_prefix, working_dir=working_directory)
    print(stderr)
  
  #import util as ut
  #make sure to remove the top 3 lines of the {pdbCode}_Repair.pdb
  import pytraj as pt
  if 'ptraj' in action:
    pdb_file = f"{pdbCode}_Repair.pdb"
    print(f"Getting range ids") 
    aaRange1="8-33"
    aaRange2="50-89" #"60-89"
    print('checkpoint 1')
    frag1=uti.get_aa(pdb_file,aaRange=[8,33])
    print(frag1)
    frag2=uti.get_aa(pdb_file,aaRange=[50,89])
    print(f"Getting interface pairs {pdb_file}") 
    uti.get_interface_pairs(pdb_file,distance,
      aaRange1=aaRange1,
      aaRange2=aaRange2) 
  
  
  if 'build' in action:
    result = []
    pdb_file = f"{pdbCode}_Repair.pdb"
    print(pdb_file)
    print("Using first few from pairsList")
    proposedMuts = uti.get_muts_pairs(pairsList)
    print('Here are the proposed mutants')
    print(proposedMuts)
    proposedMuts=['Y57F,D59E,G60E','Y57F,D59E,G60E,K69G','Y57F,K69G']#['D45K','N52K','Q54K','S55D','Y57F','L58I','D59E','G60E','K69G','G98D']
    for line in proposedMuts:
      print("Prepping mutations")
      
      muts = line.strip().split(',')
      mut_converted = []
  
      for mut in muts:
          #print(f" {mut}")
          
          words = mut.strip().split()
          modified_words = [
              f"{word[0]}{chain}{word[1:]}" for word in words
          ]
          
          out_line = ",".join(modified_words)
          mut_converted.append(out_line)
  
      print(f" ==> {mut_converted}")
      
      result_line = ",".join(mut_converted) + ";\n"
      result.append(result_line)
  
    with open(indivfilename, 'w') as f:
      delimiter=""
      f.write(delimiter.join(map(str, result))) #map to string in case of numbers
    #print(f"List successfully written to {filename}")
    if dryRun:
      print("Exiting early") 
      exit()
  
    if os.path.exists(pdb_file) and os.path.exists(rotabase_file):
        stdout2,stderr2 = run_foldx("BuildModel", 
          input_files=[pdb_file, rotabase_file,f"--mutant-file={indivfilename}"], 
          output_prefix="build_model", working_dir=working_directory)
        print(stdout2)
  
  if 'process' in action or 'build' in action:
    # get last N lines of file 
    filename = f"Average_build_model_{pdbCode}_Repair.fxout"                
    with open(indivfilename, "r") as file:
      lines = file.readlines()  # Reads each line into a list
    muts = [line.replace(";\n","") for line in lines]

    n = len(muts)
    lines = uti.get_last_n_lines(filename, n)
    for i in range(n): 
      line = lines[i]
      z = line.split("\t") 
      print(f"{muts[i]} {z[0]} {z[2]}") 
    
  if 'binding' in action:
      pdb_file= f"{pdbCode}_Repair.pdb"
      
      if os.path.exists(pdb_file) and os.path.exists(rotabase_file):
          stdout3,stderr3= run_foldx("MetalBinding",input_files=[pdb_file,rotabase_file],
                                     output_prefix="predict_output",working_dir=working_directory)
          print(stdout3)
          print(stderr3)

#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#

#
# Message printed when program run without arguments 
#
def helpmsg():
    scriptName = sys.argv[0]
    msg = f"""
Purpose: 
 
Usage:
  {scriptName} -action <build/process/ptraj>
  action='getpdb'  # get pdb structure 
  action='repair'  # fix it via foldx 
  action='ptraj'    # find residues at interface 
  action='build'   # introduce mutations via foldx [see mutList below]
  action='process' # process scores from foldx:w
  
Notes:

"""
    return msg 
#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-action"):
      print(f"{arg}") 
      action=sys.argv[i+1] 
      print(f"{arg} {action}") 
      doit(action)
      quit()
  





  raise RuntimeError("Arguments not understood")




