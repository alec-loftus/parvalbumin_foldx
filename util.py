import pytraj as pt
import re 


aa1= {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

import random 

amino_acid_map = {
    'A': ['G', 'V'], 'R': ['K', 'H'], 'N': ['D', 'Q'], 'D': ['E', 'N'],
    'C': ['S', 'T'], 'Q': ['N', 'E'], 'E': ['D', 'Q'], 'G': ['A', 'P'],
    'H': ['R', 'K'], 'I': ['L', 'V'], 'L': ['I', 'M'], 'K': ['R', 'H'],
    'M': ['L', 'F'], 'F': ['Y', 'W'], 'P': ['G', 'A'], 'S': ['T', 'C'],
    'T': ['S', 'V'], 'W': ['F', 'Y'], 'Y': ['F', 'W'], 'V': ['I', 'A']
}

# Create a hash (dictionary) mapping each amino acid to a random choice


def mutate(aaName):
  randomized_aa_hash = {aa: random.choice(choices) for aa, choices in amino_acid_map.items()}
    
  aaNew = randomized_aa_hash[aaName]
  return aaNew

  
def get_muts_pairs(pairsList):
  keep=5

  proposedMuts=[]

  for pair in random.sample(pairsList,keep):
    parts = pair.split('-')
    
    # Randomly select one of the parts
    chosen_part = random.choice(parts)
    
    # Extract the first letter
    letter = chosen_part[0]
    
    # Call the function with the letter
    mut = mutate(letter)
    newPair = f"{chosen_part}{mut}"
    print(f"{letter} --> {newPair}")  
    proposedMuts.append(newPair) 

  #proposedMuts = get_muts_pairs(pairsList)  
  return proposedMuts


def getaainfo(traj,resnum,oneLetter=False):
  atom_index= pt.select_atoms(traj.top, f":{resnum}@CA")
  resname = traj.top.atom(atom_index[-1]).resname 
  
  if oneLetter:
    resname = aa1[resname.upper()]
  return f"{resname}{resnum}"

import re 
def get_aa(pdb_file,
  aaRange=[8,33]
):
  traj = pt.load(pdb_file)      
  # Extract residue numbers and names
  daList=""
  reses=[res for res in traj.topology.residues]          
  for res in reses[aaRange[0]+0:aaRange[1]+0]:
  #for res in traj.topology.select("resid 8-10"):
    res_name = aa1[res.name.upper()]
    res_num = int(res.original_resid)               
    daList+=f"{res_name}{res_num}"
    
  print(daList)
  return daList
    

  # Print the mapping (Residue Number -> Amino Acid Name)

def get_interface_pairs(pdb_file,distance,
  aaRange1="8-33",
  aaRange2="60-89"
):
  print("WARNING: number is off!!!!!!!!! by 1 ") 
  print("Manually had to delete first few lines before ptraj would read") 
  print(f"HARDCODED AARANGE RIGHT NOW, {aaRange1}") 

  # get closest pairs 
  traj = pt.load(pdb_file)      
  # get calpha   
  #calpha1 = uti.get_calpha_coordinates(traj,residue_range=residueAB)
  #calpha2 = uti.get_calpha_coordinates(traj,residue_range=residueCD)


  # native contacts
  #data = pt.native_contacts(traj, 
  #   mask=[1,2], mask2=[20,30], ref=traj, distance=998.0)
  #command='nativecontacts :1-2 :20-30 distance 200.0 series savenative resseries present'
  print(f"Searching for all contacts within {distance}") 
  command=f"""\
nativecontacts :{aaRange1}&!@H= :{aaRange2}&!@H=  \
writecontacts native-contacts.dat \
resout resout.dat \
distance {distance} \
byresidue out all-residues.dat mindist maxdist \
map mapout gnu \
series seriesout native-contacts-series.dat 
"""
  #`print(command)
  contactData=pt.compute(command,traj)
  #print(contactData)

  # this will return a list of keys of the form
  # :33@O_:89@O
  pairs = dict()
  for s in contactData.keys():
    #print(contactData.keys())
    #print(s)
    match = re.findall(r":(\d+)@[A-Z]+", s)
    try: 
      key = f"{match[0]}-{match[1]}"
      pairs[key] = match # pair  # contactData[s]
    except:
      1 


  print("Printing interfacial pairs") 
  edges=""
  for s in pairs.keys():
    pair = []
    lint,rint = pairs[s]
    pair.append(getaainfo(traj,lint,oneLetter=True))          
    pair.append(getaainfo(traj,rint,oneLetter=True))          
    #print(pair) 
    edges+=f"{pair[0]}-{pair[1]},"
    #print(indices)
    
  print(edges)

def get_calpha_coordinates(traj, residue_range=(1, 10)):
    """
    Extracts coordinates of C-alpha atoms for a specified residue range from a trajectory.

    Args:
        traj (pytraj.Trajectory): The trajectory object.
        residue_range (tuple, optional): A tuple (start, end) defining the residue range.
                                         Defaults to (1, 10).

    Returns:
        numpy.ndarray: A NumPy array containing the coordinates of C-alpha atoms.
                       Shape: (n_frames, n_residues, 3).
                       Returns None if no C-alpha atoms are found.
    """

    start_res, end_res = residue_range
    # chain doesn't seem to work 
    #calpha_mask = f"::{chain}:{start_res}-{end_res}@CA"
    calpha_mask = f":{start_res}-{end_res}@CA"
    print(calpha_mask)

    try:
        #calpha_coords = pt.coords(traj, mask=calpha_mask)
        indices = pt.select_atoms(traj.top, calpha_mask)
        print(indices) 
        calpha_coords = traj.xyz[0:,indices,:] 
        #Reshape to (n_frames, n_residues, 3)
        n_residues = end_res - start_res + 1
        n_frames = traj.n_frames
        calpha_coords = calpha_coords.reshape((n_frames, n_residues, 3))
        return calpha_coords

    except ValueError: #No C-alpha atoms found.
        print(f"No C-alpha atoms found in residues {start_res}-{end_res}.")
        return None


def get_last_n_lines(filename, n):
  """
  Opens a file and returns a list containing the last n lines.

  Args:
    filename: The name of the file to open.
    n: The number of lines to read from the end.

  Returns:
    A list of strings, where each string is a line from the file.
    Returns an empty list if the file cannot be opened or if n is invalid.
  """
  if n <= 0:
    return []

  try:
    with open(filename, 'r') as file:
      lines = file.readlines()
      if len(lines) < n:
        return [line.rstrip('\n') for line in lines] #return all if file has less lines than n.
      else:
        return [line.rstrip('\n') for line in lines[-n:]] #slice and return last n.
  except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    return []
  except Exception as e:
    print(f"An error occurred: {e}")
    return []

