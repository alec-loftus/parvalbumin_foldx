#!/usr/bin/env python
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
#$ FAUST 
#foldx='/home/pkekeneshuskey/sources/foldxLinux64/foldx_20251231'
foldx='/home/pkekeneshuskey/source/foldxLinux64/foldx_20251231'

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
action='ptraj'    # find residues at interface 
#action='build'   # introduce mutations via foldx
#action='process' # process scores from foldx:w
pdbCode = '1RWY'
working_directory = "./" # or wherever your files are.

# selecting interface residues 
chain="A"
residueAB= (8,33) # interface 1
residueCD= (60,89)# interface 2
distance = 4.0   

# use the blast notebook for this             
mutList = [
 'I15V',
 'I15L',
 'T19S',
 'T19K',
 'T19Q',
 'A20D',
 'A21T',
 'S23T',
 'H26Y',
 'V33S',
# NOT FOUND  'G67K',
 'S78T',
 'K80E',
 'K83E'
]
#'S78K',
# 'S23K',
# 'F29W',
# 'V33F',
# 'I15T L67T'
rotabase_file = "rotabase.txt" #Make sure rotabase.txt is in the working directory or provide full path.

fasta="""\
>XP_006241991.1 parvalbumin alpha isoform X1 [Rattus norvegicus]
MSMTDLLSAEDIKKAIGAFTAADSFDHKKFFQMVGLKKKSADDVKKVFHILDKDKSGFIEEDELGSILKG
FSSDARDLSAKETKTLMAAGDKDGDGKIGVEEFSTLVAES"""


## processes 
if 'getpdb' in action: 
  # Get file 
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

if 'repair' in action: 
  # Fix file
  # - takes a few minutes 
  pdb_file = f"{pdbCode}.pdb"
  output_prefix = "repaired"
  stdout, stderr = run_foldx("RepairPDB", input_files=[pdb_file], output_prefix=output_prefix, working_dir=working_directory)
  

import util as ut
import pytraj as pt
if 'ptraj' in action:
  pdb_file = f"{pdbCode}_Repair.pdb"
  uti.get_interface_pairs(pdb_file,distance) 


if 'build' in action:
  result = []
  pdb_file = f"{pdbCode}_Repair.pdb"
  filename="individual_list.txt"
  for line in mutList:
    print("Prepping mutations") 
    words = line.split()
    modified_words = []
    for word in words:
        modified_word = word[0] + chain + word[1:]
        modified_words.append(modified_word)
    outLine = ",".join(modified_words) 
    result.append(outLine+";\n") 

  with open(filename, 'w') as f:
    delimiter=""
    f.write(delimiter.join(map(str, result))) #map to string in case of numbers
  #print(f"List successfully written to {filename}")

  if os.path.exists(pdb_file) and os.path.exists(rotabase_file):
      stdout2,stderr2 = run_foldx("BuildModel", 
        input_files=[pdb_file, rotabase_file,f"--mutant-file={filename}"], 
        output_prefix="build_model", working_dir=working_directory)
      print(stdout2)

if 'process' in action:
  # get last N lines of file 
  filename = f"Average_build_model_{pdbCode}_Repair.fxout"                
  n = len(mutList) 
  lines = uti.get_last_n_lines(filename, n)
  for i in range(n): 
    line = lines[i]
    z = line.split("\t") 
    print(f"{mutList[i]} {z[0]} {z[2]}") 
  
