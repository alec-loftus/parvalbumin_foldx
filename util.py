import pytraj as pt
import re 
def get_interface_pairs(pdb_file,distance):
  print("WARNING: number is off!!!!!!!!! by 1 ") 
  print("MWARNUALLY DELECTED CHAINS B C ") 
  # get closest pairs 
  print(pdb_file) 
  traj = pt.load(pdb_file)      
  # get calpha   
  #calpha1 = uti.get_calpha_coordinates(traj,residue_range=residueAB)
  #calpha2 = uti.get_calpha_coordinates(traj,residue_range=residueCD)


  # native contacts
  #data = pt.native_contacts(traj, 
  #   mask=[1,2], mask2=[20,30], ref=traj, distance=998.0)
  #command='nativecontacts :1-2 :20-30 distance 200.0 series savenative resseries present'
  print(f"Searching for all contacts within {distance}") 
  print("HNARDCODED RTANGE RIGHT NOW") 
  command=f"""\
nativecontacts :8-33&!@H= :60-89&!@H=  \
writecontacts native-contacts.dat \
resout resout.dat \
distance {distance} \
byresidue out all-residues.dat mindist maxdist \
map mapout gnu \
series seriesout native-contacts-series.dat 
"""
  contactData=pt.compute(command,traj)
  #print(contactData)

  # this will return a list of keys of the form
  # :33@O_:89@O
  pairs = dict()
  for s in contactData.keys():
    #print(contactData.keys())
    print(s)
    match = re.findall(r":(\d+)@[A-Z]+", s)
    try: 
      key = f"{match[0]}-{match[1]}"
      pairs[key] = match # pair  # contactData[s]
    except:
      1 

  def getinfo(traj,resnum):
    atom_index= pt.select_atoms(traj.top, f":{resnum}@CA")
    resname = traj.top.atom(atom_index[-1]).resname 
    return f"{resname}{resnum}"

  print("Printing interfacial pairs") 
  for s in pairs.keys():
    pair = []
    lint,rint = pairs[s]
    pair.append(getinfo(traj,lint))          
    pair.append(getinfo(traj,rint))          
    print(pair) 
    #print(indices)
    
    
    #print( pt.get_residue_info(traj,lint) ) 
    #print(s, pairs[s])
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

