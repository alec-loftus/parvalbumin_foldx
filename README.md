# Setup

## Dependencies

Install the required software:

* **FoldX:** Follow the installation instructions provided on the [FoldX website](http://foldxsuite.crg.eu/).
* **PyTraj:**
    ```bash
    conda install -c ambermd pytraj
    ```
    *Note: PyTraj is often included with AmberTools.*
* **BioBlast:**
    ```bash
    conda install -c bioconda bioblast # Or follow specific installation guide if different
    ```

* need to make a symbolic link to the foldx rotatebase.txt file 

## Environment Configuration

1.  **Source Configuration File:**
    Before running, make sure to source the necessary environment variables or settings from your configuration file:
    ```bash
    source config.bash
    ```

2.  **Activate Conda Environment:**
    Activate the Conda environment where AmberTools (and potentially the other dependencies) are installed. This example uses `AmberTools23`:
    ```bash
    conda activate AmberTools23
    ```

## Running the Design

Execute the design script 
* To get your protein of interest *assumes 1RwY*
    ```bash
    # Example command structure (replace with your actual command)
    python design.py -action getpdb
    ```
* To repair it using foldx
    ```bash
    python design.py -action repair 
    ```
    Note, i manually prepped the structure by removing its other chains and heteroatoms 
* To identify the closest neighboring residues on the interface, use the `-ptraj` option when running the design script:
    ```bash
    python design.py -action ptraj
    ```


  action='getpdb'  # get pdb structure
  
  action='repair'  # fix it via foldx
  
  action='ptraj'    # find residues at interface
  
  action='build'   # introduce mutations via foldx [see mutList below]
  
  action='process' # process scores from foldx:w

