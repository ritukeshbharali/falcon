# Python script to create POD basis

import ast
import argparse
import os.path
import numpy as np
import re
import shutil
import vtk
import subprocess
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt

def file_to_wip(fname):

    # Define a VTK reader
    reader = vtk.vtkXMLUnstructuredGridReader()

    # Parse files
    print('   Reading', fname, '...')
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()

    # Get some data
    nNodes     = reader.GetNumberOfPoints()
    nElems     = reader.GetNumberOfCells()
    nCellData  = reader.GetNumberOfCellArrays()
    nPointData = reader.GetNumberOfPointArrays()

    # Get integration point weights (if present)
    dofIdx = -1
    for idx in range(reader.GetNumberOfCellArrays()):
        print(reader.GetCellArrayName(idx))
        if reader.GetCellArrayName(idx) == "wip":
            dofIdx = idx
            print(f'   - wip found! ...')
            wip = vtk_to_numpy(data.GetCellData().GetArray(dofIdx))
        else:
            print(f'   - {reader.GetCellArrayName(idx)} found! ...')

    return np.squeeze(wip)

def file_to_slip(fname):

    # Define a VTK reader
    reader = vtk.vtkXMLUnstructuredGridReader()

    # Parse files
    print('   Reading', fname, '...')
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()

    # Get some data
    nNodes     = reader.GetNumberOfPoints()
    nElems     = reader.GetNumberOfCells()
    nCellData  = reader.GetNumberOfCellArrays()
    nPointData = reader.GetNumberOfPointArrays()

    # Get integration point weights (if present)
    dofIdx = -1
    for idx in range(reader.GetNumberOfCellArrays()):
        if reader.GetCellArrayName(idx) == "slip":
            dofIdx = idx
            print(f'   - slip found! ...')
            slip = vtk_to_numpy(data.GetCellData().GetArray(dofIdx))
        else:
            RuntimeError(f'   - slip not found! ...')

        rows, cols = slip.shape
        slip = slip.reshape(rows*cols,1)

    return np.squeeze(slip)

class POD():
    """ Initialize the POD class """
    def __init__(self, props:dict):

        # Get some props
        if 'basis_dof' in props:
            self.basisDof   = props['basis_dof']
        else:
            RuntimeError('basis_dof must be provided!')

        if 'sens_dof' in props:
            self.sensDof   = props['sens_dof']
        else:
            RuntimeError('sens_dof must be provided!')

        if 'cutoff_tol' in props:
            self.cutOffTol   = props['cutoff_tol']
        else:
            print('cutoff_tol is not provided!')
            print('cutoff_tol is set to default: 1e-6!')
            self.cutOffTol   = 1.e-6

        if 'nsteps_online' in props:
            self.nStepsOnline   = props['nsteps_online']
        else:
            RuntimeError('nsteps_online is not provided!')

        if 'jem_export' in props:
            self.jem_export   = props['jem_export']
        else:
            RuntimeError('jem_env must be provided!')

        if 'jive_export' in props:
            self.jive_export  = props['jive_export']
        else:
            RuntimeError('jive_env must be provided!')

        if 'run_falcon' in props:
            self.run_falcon   = props['run_falcon']
        else:
            RuntimeError('run_falcon must be provided!')

        # Set some paths to required directories
        self.basePath  = os.getcwd()
        self.offPath   = os.path.join(self.basePath,'offline')
        self.onPath    = os.path.join(self.basePath,'online')
        self.basisPath = os.path.join(self.basePath,'basis')
        self.sensPath  = os.path.join(self.basePath,'sensitivity')
        self.unitPath = os.path.join(self.sensPath,'unitLoad')

        # Make the additional directories
        if not os.path.exists(self.onPath):
            os.makedirs(self.onPath, exist_ok=True)
            print(f'Directory created: {self.onPath}')
        else:
            print(f'Directory already exists: {self.onPath}')

        if not os.path.exists(self.basisPath):
            os.makedirs(self.basisPath, exist_ok=True)
            print(f'Directory created: {self.basisPath}')
        else:
            print(f'Directory already exists: {self.basisPath}')

        if not os.path.exists(self.sensPath):
            os.makedirs(self.sensPath, exist_ok=True)
            print(f'Directory created: {self.sensPath}')
        else:
            print(f'Directory already exists: {self.sensPath}')

        if not os.path.exists(self.unitPath):
            os.makedirs(self.unitPath, exist_ok=True)
            print(f'Directory created: {self.unitPath}')
        else:
            print(f'Directory already exists: {self.unitPath}')

        # Set some default parameters
        self.nSteps     = -1    # simulation steps
        self.nNodes     = -1    # number of nodes
        self.nElems     = -1    # number of elements
        self.nCellData  = -1    # number of cell data
        self.nPointData = -1    # number of point data
        self.wip        = None  # integration point weights
        self.snapMat    = None  # snapshot matrix
        self.C          = None  # correlation matrix
        self.eigvals    = None  # eigen values
        self.eigvecs    = None  # eigen vectors
        self.podBasis   = None  # POD basis matrix
        self.nBasis     = -1    # number of POD basis


    def run_offline(self):
        """ Run the offline simulation """

        # Prepare the environment
        os.system(self.jem_export)
        os.system(self.jive_export)

        # Change directory
        os.chdir(self.offPath)

        # Run simulation
        os.system(self.run_falcon)

        # Change directory
        os.chdir(self.basePath)


    def __create_snapshot_matrix(self):
        """ Creates a snapshot matrix from VTU files """

        # Change directory
        os.chdir(self.offPath)

        # Read the pvd file
        with open('problem_out.pvd','r') as file:
            lines = file.readlines()

        # Extract the number of files (should be contiguous)
        self.nSteps = sum(1 for line in lines if '<DataSet file=' in line)

        # File initial part
        fname0    = 'problem_out_p0_'

        # Define a VTK reader
        reader = vtk.vtkXMLUnstructuredGridReader()

        # Parse the first file
        fname = fname0+str(0)+'.vtu'
        print('   Reading', fname, '...')
        reader.SetFileName("problem_out_p0_0.vtu")
        reader.Update()
        data = reader.GetOutput()

        # Get some mesh data
        self.nNodes     = reader.GetNumberOfPoints()
        self.nElems     = reader.GetNumberOfCells()
        self.nCellData  = reader.GetNumberOfCellArrays()
        self.nPointData = reader.GetNumberOfPointArrays()

        # Extract the integration point weights (if present)
        dofIdx = -1
        for idx in range(reader.GetNumberOfCellArrays()):
            if reader.GetCellArrayName(idx) == "wip":
                dofIdx = idx
                print(f'   - wip found! ...')
                self.wip = vtk_to_numpy(data.GetCellData().GetArray(dofIdx))
            else:
                print(f'   - wip not found! Skipping ...')

        # Extract the basis dof as CellData
        dofIdx = -1
        for idx in range(reader.GetNumberOfCellArrays()):
            if reader.GetCellArrayName(idx) == self.basisDof:
                dofIdx = idx
                print(f'   - dof {self.basisDof} found! ...')
                basisField = np.matrix(vtk_to_numpy(data.GetCellData().GetArray(dofIdx)))
            
        if dofIdx == -1:
            RuntimeError(f'   - dof {self.basisDof} not found!')

        # Get the shape of the basis field (CellData)
        # Columns will give the number of components
        rows, cols = basisField.shape
        self.basisNComp = cols

        # Print some info to terminal
        print(f'   - Found {cols} components of {self.basisDof}! ...')
        print('   - Components will be stored in block format!')

        # Check for cell data consistency
        assert rows == self.nElems

        # Fill in the snapshot matrix
        self.snapMat = np.matrix(np.zeros((rows*cols, self.nSteps)))

        for i in range(cols):
            startIdx = i * rows
            endIdx   = (i+1) * rows
            self.snapMat[startIdx:endIdx, 0] = basisField[:,i]

        # Parse the remaining VTU files and fill the snapshot
        # matrix
        for i in range(1,self.nSteps):

            # Define a VTK reader
            reader = vtk.vtkXMLUnstructuredGridReader()

            # Get data from vtu file
            fname = fname0+str(i)+'.vtu'
            print('   Reading', fname, '...')
            reader.SetFileName(fname)
            reader.Update()
            data = reader.GetOutput()

            # Get the basis dof
            basisField = np.matrix(vtk_to_numpy(data.GetCellData().GetArray(dofIdx)))

            for j in range(cols):
                startIdx = j * rows
                endIdx   = (j+1) * rows
                self.snapMat[startIdx:endIdx, i] = basisField[:,j]
        
        # Change directory
        os.chdir(self.basePath)


    def __create_correlation_matrix(self):
        """ Creates a correction matrix """
        self.C = np.matmul(self.snapMat.T,self.snapMat)


    def __do_eigendecomposition(self):
        """ Perform eigen decomposition """

        eigvals, eigvecs = np.linalg.eig(self.C)

        # Store only real parts
        self.eigvals = np.real(eigvals)
        self.eigvecs = np.real(eigvecs)

        self.eigvals /= np.max(self.eigvals)

        plt.figure(figsize=(10, 6))
        plt.plot(self.eigvals, 'o-', markersize=5)
        plt.xlabel('Elements [-]')
        plt.ylabel('Eigen Values [-]')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig('pod_values.png')
        plt.show(block=False)
        plt.pause(2)
        plt.close()

    #-------------------------------------------------------------------------
    # Create POD basis (create snapshot matrix from offline simulation, then
    # compute the correlation matrix, followed by eigen decomposition, and 
    # finally selection of basis based on a cut off tolerance.)
    #-------------------------------------------------------------------------

    def create_pod_basis(self):

        # Create the snapshot matrix
        self.__create_snapshot_matrix()

        # Create the correlation matrix
        self.__create_correlation_matrix()

        # Perform the eigenvalue decomposition
        self.__do_eigendecomposition()

        # Get index of eigvals corresponding to cutoff tol
        self.cutOffIdx = np.where(self.eigvals \
                                   < self.cutOffTol)[0][0]

        # Compute the pod basis
        self.podBasis = np.matrix(np.matmul(self.snapMat,\
                                  self.eigvecs[:,0:self.cutOffIdx]))

        # Get the shape of POD basis
        rows, cols  = self.podBasis.shape
        self.nBasis = cols

        # Change directory
        os.chdir(self.offPath)

        # Read the iPoint-Element mapping file
        with open('problem.ipElemMap','r') as file:
            lines = file.read()

        # Convert lines into an ipMap dictionary
        ipMap = ast.literal_eval(lines)

        # Change to basis directory
        os.chdir(self.basisPath)

        # Store basis in Jive format constraints file for
        # sensitivity computations. These constraints will
        # be enforced on the dummy integration point nodes.

        # Loop over the basis, create a new file for each basis.
        for ibasis in range(cols):

            filename = 'pod_'+self.basisDof+'.cons'+str(ibasis)
            print(f'   Creating {filename}')

            with open(filename,'w') as fout:
                fout.write("<NodeConstraints>\n")

                # Loop over the components (if any!)
                for icomp in range(self.basisNComp):

                    # Extract the pod basis block for this component
                    startIdx = (icomp) * self.nElems
                    endIdx   = (icomp+1) * self.nElems
                    podComp  = self.podBasis[startIdx:endIdx,ibasis]

                    for idx in range(self.nElems):
                        fout.write(f"{self.basisDof}{str(icomp)}[{int(ipMap[idx+1])}]={podComp[idx,0]}; \n")
                fout.write("</NodeConstraints>\n")

        # Store the basis files as plain text for the online
        # simulation. Loop over the basis, create a new file
        # for each basis.
        for ibasis in range(cols):

            filename = 'pod_'+self.basisDof+'.basis'+str(ibasis)
            print(f'   Creating {filename}') 

            # Create a new basis file
            with open('pod_'+self.basisDof+'.basis'+str(ibasis),'w') as fout:

                # Loop over the elements
                for idx in range(self.nElems):
                    for icol in range(self.basisNComp):
                        fout.write(f"{self.podBasis[idx+icol*self.nElems,ibasis]} ")
                    fout.write(f"\n")

        # Change to base directory
        os.chdir(self.basePath)

    #-------------------------------------------------------------------------
    # Prepare sensitivity directories (files copied from offline and 
    # modified with regex)
    #-------------------------------------------------------------------------

    def prepare_sensitivity_directories(self):

        # Get files to be copied from offline studies
        proFile  = 'problem.pro'
        proPath  = os.path.join(self.offPath,proFile)

        dataFile = 'problem.data'
        dataPath = os.path.join(self.offPath,dataFile)

        meshFile = 'problem.mesh'
        meshPath = os.path.join(self.offPath,meshFile)

        # Copy the files
        print(f'Copying .pro,.data,.mesh files to {self.unitPath}')
        shutil.copy(proPath,  self.unitPath)
        shutil.copy(dataPath, self.unitPath)
        shutil.copy(meshPath, self.unitPath)

        # Change directory to unitLoad
        os.chdir(self.unitPath)

        # Modify the .pro file
        with open(proFile,'r') as file:
            text = file.read()

            # Set number of steps to 1
            text = re.sub(r'(runWhile\s*=\s*\"i\s*<\s*)\d+(\s*\";)', r'\g<1>1\2', text)

            # fgMode is switched off
            text = re.sub(r'(fgMode\s*=\s*)true(;)', r'\1false\2', text)

            # loadIncr is set to 1
            text = re.sub(r'^\s*loadIncr\s*=\s*\d+(\.\d+)?;', '         loadIncr = 1.0;', text, flags=re.MULTILINE)

            # Remove lodi module
            text = re.sub(r'\"lodi\",\s*', '', text)

            # Remove graph module
            text = re.sub(r'\"graph\",\s*', '', text)

            # FemView module is switched off
            text = re.sub(r'\"view\",\s*', '', text)
            text = re.sub(r'\"view2\",\s*', '', text)

            # Re-write the file
            with open(proFile, 'w') as file:
                file.write(text)

        # Extract the domains in the mesh (they should be 
        # named as Domain1Elems, Domain2Elems) for grains 1,2 ...
        cons = ''
        with open(proFile,'r') as file:
            text = file.read()

            # Bulk section patterns
            bulk_pattern = re.compile(r'bulk\d+\s*=\s*"GradientCrystalPlasticity"\s*{[^}]*?elements\s*=\s*"([^"]+)"', re.DOTALL)
            
            # Find all matches
            matches = bulk_pattern.findall(text)

            # Create a bulk elements dict
            bulk_elems = dict()

            # Find each match
            for i, match in enumerate(matches, start=1):
                domain_match = re.match(r'(Domain\d+)', match)
                if domain_match:
                    bulk_elems[f'bulk{i}'] = domain_match.group(1)

            # Create constraint string
            for bulk, elements in bulk_elems.items():
                for i in range(self.basisNComp):
                    cons += f'  {self.basisDof}{i}[{elements}IPNodes] = 0.0;\n'

        # Add Node constraints for the basis dofs
        with open(dataFile,'r') as file:
            lines = file.readlines()

        # Find insertion point
        for j, line in enumerate(lines):
            if '</NodeConstraints>' in line:
                insertion_point = j - 1
                lines.insert(insertion_point, cons)
                break

        with open(dataFile,'w') as file:
            file.writelines(lines)

        # Change to base directory
        os.chdir(self.basePath)

        # Loop over the basis
        for ibasis in range(self.nBasis):

            # Set up a new directory
            dirPath = os.path.join(self.sensPath,self.basisDof+str(ibasis))

            if not os.path.exists(dirPath):
                os.makedirs(dirPath, exist_ok=True)
                print(f'Directory created: {dirPath}')
            else:
                print(f'Directory already exists: {dirPath}')

            # Get the corresponding basis file and path
            basisFile      = f'pod_{self.basisDof}.cons{ibasis}'
            basisFilePath  = os.path.join(self.basisPath,basisFile)

            # Copy simulation input files
            shutil.copy(basisFilePath,  dirPath)
            shutil.copy(proPath,   dirPath)
            shutil.copy(dataPath,  dirPath)
            shutil.copy(meshPath,  dirPath)

            # Change directory
            os.chdir(dirPath)

            # Modify the problem.pro simulation file
            with open(proFile,'r') as file:
                text = file.read()

                # Set number of steps to 1
                text = re.sub(r'(runWhile\s*=\s*\"i\s*<\s*)\d+(\s*\";)', r'\g<1>1\2', text)

                # fgMode is switched off
                text = re.sub(r'(fgMode\s*=\s*)true(;)', r'\1false\2', text)

                # Neumann constraints are removed
                text = re.sub(r'"cons",?', '', text)
                text = re.sub(r'cons\s*=\s*\{[^}]*\{[^}]*\{[^}]*\}[^}]*\}[^}]*\};', '', text, flags=re.DOTALL)

                # Remove lodi module
                text = re.sub(r'\"lodi\",\s*', '', text)

                # Remove graph module
                text = re.sub(r'\"graph\",\s*', '', text)

                # FemView module is switched off
                text = re.sub(r'\"view\",\s*', '', text)
                text = re.sub(r'\"view2\",\s*', '', text)

                with open(proFile, 'w') as file:
                    file.write(text)

            # Insert pod cons file to *.data
            line_insert = '<Include source="pod_'+self.basisDof \
                           +'.cons'+str(ibasis)+'"/>\n'

            with open(dataFile,'r') as file:
                lines = file.readlines()

            for j, line in enumerate(lines):
                if '<Include source="problem.mesh"/>' in line:
                    insertion_point = j + 1
                    lines.insert(insertion_point, line_insert)
                    break

            with open(dataFile,'w') as file:
                file.writelines(lines)

        # Change to base directory
        os.chdir(self.basePath)

    #-------------------------------------------------------------------------
    # Run sensitivity simulations (w.r.t to each basis and unit load)
    #-------------------------------------------------------------------------

    def run_sensitivity_simulations(self):

        # Prepare the environment for Jive simulation
        os.system(self.jem_export)
        os.system(self.jive_export)

        # Change directory to unitLoad
        os.chdir(self.unitPath)

        # Run sensitivity simulation for the unit load case
        os.system(self.run_falcon)

        # Run sensitivity simulations for each basis
        for ibasis in range(self.nBasis):

            # Get the corresponding directory path
            dirPath = os.path.join(self.sensPath,self.basisDof+str(ibasis))

            # Change to basis directory
            os.chdir(dirPath)

            # Run the simulation
            os.system(self.run_falcon)
        
        # Change to base directory
        os.chdir(self.basePath)

    #-------------------------------------------------------------------------
    # Prepare online directories (files copied from offline and 
    # modified with regex)
    #-------------------------------------------------------------------------

    def prepare_online_directory(self):

        # Prepare the files to be copied from offline studies
        proFile  = 'problem.pro'
        proPath  = os.path.join(self.offPath,proFile)

        dataFile = 'problem.data'
        dataPath = os.path.join(self.offPath,dataFile)

        meshFile = 'problem.mesh'
        meshPath = os.path.join(self.offPath,meshFile)

        # Copy the simulation input files
        print(f'Copying .pro,.data,.mesh files to {self.onPath}')
        shutil.copy(proPath,  self.onPath)
        shutil.copy(dataPath, self.onPath)
        shutil.copy(meshPath, self.onPath)

        # Copy all the basis files
        for ibasis in range(self.nBasis):

            # Get the corresponding basis file path
            basisFile      = f'pod_{self.basisDof}.basis{ibasis}'
            basisFilePath  = os.path.join(self.basisPath,basisFile)

            # Copy
            shutil.copy(basisFilePath, self.onPath)

        # Get the solution fields sensitivity files
        for ibasis in range(self.nBasis):

            dirPath = os.path.join(self.sensPath,self.basisDof+str(ibasis))

            if not os.path.exists(dirPath):
                RuntimeError(f'Directory does not exist: {dirPath}')
            else:
                print(f'Directory exist: {dirPath}')

            # Change directory
            os.chdir(dirPath)

            # Define a VTK reader
            reader = vtk.vtkXMLUnstructuredGridReader()

            # Parse the vtu file
            fname = 'problem_out_p0_0.vtu'
            print('Reading', fname, 'in', dirPath)
            reader.SetFileName("problem_out_p0_0.vtu")
            reader.Update()
            data = reader.GetOutput()

            # Get index of sensitivity dof (this dof is required
            # to assemble the POD online system of equations)
            sdofIdx = -1
            for idx in range(reader.GetNumberOfCellArrays()):
                if reader.GetCellArrayName(idx) == self.sensDof:
                    sdofIdx = idx

            if sdofIdx == -1:
                RuntimeError(f'dof {self.sensDof} does not exist in vtu files')

            # Get sensitivity dof field
            sensField = np.matrix(vtk_to_numpy(data.GetCellData().GetArray(sdofIdx)))

            rows, cols = sensField.shape

            # Print some info to terminal
            print(f'   - Found {cols} components of {self.sensDof}! ...')

            # Print data to file
            basisFile     = f'pod_{self.sensDof}Hat.basis{ibasis}'
            basisFilePath = os.path.join(self.onPath,basisFile)

            # Create the file
            with open(basisFilePath,'w') as fout:
                # Loop over rows
                for irow in range(rows):
                    for jcol in range(cols):
                        fout.write(f"{sensField[irow,jcol]} ")
                    fout.write(f"\n")

        # Copy data from unit load sensitivity computations
        os.chdir(self.unitPath)

        # Define a VTK reader
        reader = vtk.vtkXMLUnstructuredGridReader()

        # Parse the vtu file
        fname = 'problem_out_p0_0.vtu'
        print('Reading', fname, 'in', dirPath)
        reader.SetFileName("problem_out_p0_0.vtu")
        reader.Update()
        data = reader.GetOutput()

        # Get index of sdof
        sdofIdx = -1
        for idx in range(reader.GetNumberOfCellArrays()):
            if reader.GetCellArrayName(idx) == self.sensDof:
                sdofIdx = idx

        if sdofIdx == -1:
            RuntimeError(f'dof {self.sensDof} does not exist in vtu files')

        # Get sensitivity dof field
        sensField = np.matrix(vtk_to_numpy(data.GetCellData().GetArray(sdofIdx)))

        rows, cols = sensField.shape

        # Print some info to terminal
        print(f'   - Found {cols} components of {self.sensDof}! ...')

        # Print data to file
        basisFile     = f'pod_{self.sensDof}Hat0.basis'
        basisFilePath = os.path.join(self.onPath,basisFile)

        # Create the file
        with open(basisFilePath,'w') as fout:
            # Loop over rows
            for irow in range(rows):
                for jcol in range(cols):
                    fout.write(f"{sensField[irow,jcol]} ")
                fout.write(f"\n")

        # Go to online directory
        os.chdir(self.onPath)

        # Modify the problem.pro simulation file
        with open(proFile,'r') as file:
            text = file.read()

        # Set number of steps to 1
        text = re.sub(r'(runWhile\s*=\s*\"i\s*<\s*)\d+(\s*\";)', rf'\g<1>{self.nStepOnline}\2', text)

        # Replace models with a single model
        text = re.sub(r'models\s*=\s*\[[^\]]+\];', 'models = ["bulk1"];', text, flags=re.MULTILINE)

        # Replace model type
        text = re.sub(r'bulk1\s*=\s*"GradientCrystalPlasticity"', 'bulk1 = "GradientCrystalPlasticityROM"', text)

        # Replace elements for bulk1 with full domain
        text = re.sub(r'elements\s*=\s*"Domain1Elems";', 'elements = "elems2DElems";', text)

        # Replace 'ipNodes  = "Domain1IPNodes";' with 'modeNodes = "modeNodes";'
        text = re.sub(r'ipNodes\s*=\s*"Domain1IPNodes";', 'modeNodes = "modeNodes";', text)

        # Remove the 'material' block
        material_pattern = re.compile(r'material\s*=\s*\{[^}]*\};', re.DOTALL)
        text = material_pattern.sub('', text)

        # Remove lines with 'rotation = some number;'
        text = re.sub(r'rotation\s*=\s*[0-9.e+-]+\s*;', '', text)

        # Remove the 'cons' block
        #cons_pattern = re.compile(r'cons\s*=\s*\{[^}]*\};', re.DOTALL)
        #text = cons_pattern.sub('', text)

        # Remove the 'lodi' block
        #lodi_pattern = re.compile(r'lodi\s*=\s*\{[^}]*\};', re.DOTALL)
        #text = lodi_pattern.sub('', text)

        # Get lines about the material properties
        lines_to_insert = set()
        matches = re.findall(r'\b(tstar|tauY|n)\b\s*=\s*([0-9.e+-]+);', text)
        for match in matches:
            lines_to_insert.add(f'      {match[0]} = {match[1]};')

        # Remove all slip systems
        slips_pattern = re.compile(r'slips\s*=\s*\[([^\]]+)\];')
        slips_match = slips_pattern.search(text)

        if slips_match:
            slips_content = slips_match.group(1)
            slip_names = re.findall(r'slip\d+', slips_content)

            # Remove the slip blocks
            for slip_name in slip_names:
                slip_pattern = re.compile(rf'{slip_name}\s*=\s*\{{[^}}]*\}};', re.DOTALL)
                text = slip_pattern.sub('', text)
            text = slips_pattern.sub('', text)

        # Find the line with dtime
        match_dtime = re.search(r'dtime\s*=\s*[^\n]+;', text)
        if match_dtime:
            dtime_line = match_dtime.group(0)
            index_dtime = match_dtime.end()+1

            # Insert stored lines after 'dtime  ='
            text = text[:index_dtime] + '\n'.join(lines_to_insert) + '\n\n' + text[index_dtime:]
        
        # Insert more info
        # Find lines with 'tstar = some number;'
        matches = re.finditer(r'tstar\s*=\s*([0-9.e+-]+);', text)

        # Prepare lines to insert after 'tstar = some number;'
        lines_to_insert = []
        for match in matches:
            tstar_value = match.group(1)
            lines_to_insert.append(f'      tstar = {tstar_value};\n'
                                   f'      nmodes = {self.nBasis};\n'
                                   f'      nslips = {self.basisNComp};\n'
                                   f'      slipHat0File = "pod_slipHat0.basis";\n'
                                   f'      slipHatFile = "pod_slipHat.basis";\n'
                                   f'      tauHatFile = "pod_tau.basis";\n'
                                   f'      basis = {list(range(self.nBasis))};\n')

        # Insert the constructed lines after each 'tstar = some number;'
        for insert_line in reversed(lines_to_insert):  # Insert in reverse to maintain line order
            text = re.sub(r'tstar\s*=\s*[0-9.e+-]+;', insert_line, text, count=1)

        # Setup the required modules
        text = re.sub(r'modules\s*=\s*\[[^\]]+\];', 'modules = ["solver","vtkout"];', text, flags=re.MULTILINE)

        # Setup the solver
        text = re.sub(r'^\s*parFactorize\s*=\s*true;', '        parFactorize = false;', text, flags=re.MULTILINE)
        text = re.sub(r'^\s*numThreads\s*=\s*4;', '        numThreads = 1;', text, flags=re.MULTILINE)


        # Setup the post-processing data
        text = re.sub(r'cellData\s*=\s*\[[^\]]+\];', 'cellData = ["slip","tau"];', text, flags=re.MULTILINE)

        # Re-write the pro file
        with open(proFile,'w') as file:
            file.write(text)

        # Modify problem.data
        with open(dataFile, 'r') as file:
            text = file.read()

            match = re.search(re.escape('<Include source="problem.mesh"/>'), text)

            if match:
                # Truncate the content at the position of the marker
                text = text[:match.end()]

            # Write back the truncated content
            with open(dataFile, 'w') as file:
                file.write(text)

        # Modify mesh file
        with open(meshFile, 'r') as file:
            text = file.read()

            pattern = r'<NodeGroup name=Domain1IPNodes>\s*{\s*([^}]*)\s*}\s*</NodeGroup>'
            match = re.search(pattern, text, re.DOTALL)

            if match:
                numbers_str = match.group(1).strip()
                numbers = re.findall(r'\d+', numbers_str)  # Find all numbers in the string
                numbers = list(map(int, numbers))  # Convert strings to integers

            # Generate output for <NodeGroup name=modeNodes>
            mode_nodes_content = f'<NodeGroup name=modeNodes>\n{{{" ".join(str(num) for num in numbers[0:self.nBasis])}}}\n</NodeGroup>\n'

        # Append or overwrite the content at the end of the file
        with open(meshFile, 'a') as file:
            file.write(mode_nodes_content)

        # Create directories based on number of modes selected

        # Prepare the files to be copied from offline studies
        proFile  = 'problem.pro'
        proPath  = os.path.join(self.onPath,proFile)

        dataFile = 'problem.data'
        dataPath = os.path.join(self.onPath,dataFile)

        meshFile = 'problem.mesh'
        meshPath = os.path.join(self.offPath,meshFile)

        for ibasis in range(self.nBasis):

            dirPath = os.path.join(self.onPath,'basis'+str(ibasis+1))

            if not os.path.exists(dirPath):
                os.makedirs(dirPath, exist_ok=True)
                print(f'Directory created: {dirPath}')
            else:
                print(f'Directory already exists: {dirPath}')

            # Copy the files
            print(f'Copying .pro,.data,.mesh files to {dirPath}')
            shutil.copy(proPath,  dirPath)
            shutil.copy(dataPath, dirPath)
            shutil.copy(meshPath, dirPath)

            # Copy unit load sensitivity files
            sens0File       = f'pod_{self.sensDof}Hat0.basis'
            sens0FilePath   = os.path.join(self.onPath,sens0File)
            shutil.copy(sens0FilePath, dirPath)

            for jbasis in range(ibasis+1):

                # Get the corresponding basis file path
                basisFile      = f'pod_{self.basisDof}.basis{jbasis}'
                basisFilePath  = os.path.join(self.onPath,basisFile)

                sensFile       = f'pod_{self.sensDof}Hat.basis{jbasis}'
                sensFilePath   = os.path.join(self.onPath,sensFile)

                # Copy
                shutil.copy(basisFilePath, dirPath)
                shutil.copy(sensFilePath,  dirPath)

            # Change to dirPath
            os.chdir(dirPath)

            # Modify the problem.pro simulation file
            with open(proFile,'r') as file:
                text = file.read()

            text = re.sub(r'nmodes\s*=\s*\d+\s*;', f'nmodes = {ibasis+1};', text)

            # Generate the new basis array based on the value of n
            newBasisArray   = list(range(ibasis+1))
            newBasisContent = f'basis = [{",".join(map(str, newBasisArray))}];'

            # Define the pattern to find the basis array
            pattern = r'basis\s*=\s*\[.*?\]\s*;'

            # Replace the existing basis array with the new one
            text = re.sub(pattern, newBasisContent, text, flags=re.DOTALL)

            # Re-write the pro file
            with open(proFile,'w') as file:
                file.write(text)

            # Modify mesh file
            with open(meshFile, 'r') as file:
                text = file.read()

            pattern = r'<NodeGroup name=Domain1IPNodes>\s*{\s*([^}]*)\s*}\s*</NodeGroup>'
            match = re.search(pattern, text, re.DOTALL)

            if match:
                numbers_str = match.group(1).strip()
                numbers = re.findall(r'\d+', numbers_str)  # Find all numbers in the string
                numbers = list(map(int, numbers))  # Convert strings to integers

            # Generate output for <NodeGroup name=modeNodes>
            mode_nodes_content = f'<NodeGroup name=modeNodes>\n{{{" ,".join(str(num) for num in numbers[0:ibasis+1])},}}\n</NodeGroup>\n'

            # Append or overwrite the content at the end of the file
            with open(meshFile, 'a') as file:
                file.write(mode_nodes_content)

    #--------------------------------------------------------------------------
    # Run online simulation (for number of basis chosen)
    #--------------------------------------------------------------------------

    def run_online_simulation(self):

        # Prepare the environment for Jive simulation
        os.system(self.jem_export)
        os.system(self.jive_export)

        # Loop over the directories with varying number of basis selected
        for ibasis in range(self.nBasis):

            # Get the directory path
            dirPath = os.path.join(self.onPath,'basis'+str(ibasis+1))

            # Change directory
            os.chdir(dirPath)

            # Run online simulation
            os.system(self.run_falcon)

        # Print completion
        print(f'\n\nAll online simulations complete!\n\n')

    #--------------------------------------------------------------------------
    # Post-process data (compute error w.r.t number of basis chosen)
    #--------------------------------------------------------------------------

    def post_process_data(self,finalStep=-1):

        # Override the number of steps (if required)
        if finalStep > 0:
            self.nStepsOnline = finalStep

        # First get the error from the simulation log files
        # Initialize some variables
        error       = []
        modes       = []

        # Define the regex pattern to match "error: some float number"
        pattern = r'error: ([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)'

        # Loop over the online sub-directories to extract the error
        for ibasis in range(self.nBasis):

            # Get the directory path
            dirPath = os.path.join(self.onPath,'basis'+str(ibasis+1))

            # Change directory
            print(f'Entering {dirPath}')
            os.chdir(dirPath)

            # Initialize an empty array to errors from all steps
            all_errors = []

            # Read the file line by line
            with open('problem.log', 'r') as file:
                for line in file:
                    # Use regex to find matches
                    match = re.search(pattern, line)
                    if match:
                        # Extract the float number from the match
                        float_number = float(match.group(1))
                        # Add the float number to the array
                        all_errors.append(float_number)

            np.savetxt(f'errorMeasureWith{str(ibasis+1)}Basis.txt',all_errors)
            modes.append(int(ibasis+1))
            error.append(all_errors[self.nStepsOnline-1])

        plt.figure(figsize=(8, 6))  # Specify figure size
        plt.plot(modes, error, marker='o', linestyle='-', color='black')
        plt.xlabel('# of modes [-]', fontsize=14)
        plt.ylabel('Error [-]', fontsize=14)
        plt.grid(True)
        #plt.legend()
        #plt.xlim(0, 10)  # Adjust limits based on your data
        #plt.ylim(0, 1e-1)  # Adjust limits based on your data
        plt.yscale("log")   
        plt.tick_params(axis='both', which='major', labelsize=12)
        os.chdir(self.basePath)
        plt.savefig('plot_errormeasure.png', dpi=400)  # Adjust dpi as needed
        plt.show(block=False)
        plt.pause(3)
        plt.close()

        # Second, compute the relative L2 error in slip computed from
        # online simulations
        # Initialize some variables
        offSlip = np.empty((self.nElems*self.basisNComp,self.nStepsOnline))
        errorL2 = []

        # Change directory
        os.chdir(self.offPath)

        # First get the weight of integration points
        fname = f'problem_out_p0_0.vtu'
        wip   = file_to_wip(fname)

        # Expand wip for each component
        wipC = np.empty((self.basisNComp * len(wip),), dtype=wip.dtype)
        for i in range(self.basisNComp):
            wipC[i::self.basisNComp] = wip

        # Loop over the vtu files 
        for i in range(self.nStepsOnline):

            # Get filename
            fname = f'problem_out_p0_{i}.vtu'

            # Get the slip
            offSlip[:,i] = file_to_slip(fname)

        # Loop over the basis directories
        for ibasis in range(self.nBasis):

            dirPath = os.path.join(self.onPath,f'basis{ibasis+1}')
            
            print(f'Entering {dirPath}')
            os.chdir(dirPath)
            
            # Initialize error for this basis
            berrN = 0.0
            berrD = 0.0

            # Loop over the vtu files
            for j in range(self.nStepsOnline):

                # Get filename
                fname = f'problem_out_p0_{j}.vtu'

                # Get the slip
                slip = file_to_slip(fname)

                for k in range(self.basisNComp * len(wip)):

                    # Compute error
                    berrN += ((slip[k] - offSlip[k,j])**2) \
                             * wipC[k]
                    berrD += ((offSlip[k,j])**2) \
                             * wipC[k]

            berr = berrN/berrD
            errorL2.append(berr)

        plt.figure(figsize=(8, 6))  # Specify figure size
        plt.plot(modes, errorL2, marker='o', linestyle='-', color='black')
        plt.xlabel('# of modes [-]', fontsize=14)
        plt.ylabel('L2 Error (Slip) [-]', fontsize=14)
        plt.grid(True)
        #plt.legend()
        #plt.xlim(0, 10)  # Adjust limits based on your data
        #plt.ylim(0, 1e-1)  # Adjust limits based on your data
        plt.yscale("log")   
        plt.tick_params(axis='both', which='major', labelsize=12)
        os.chdir(self.basePath)
        plt.savefig('plot_l2error.png', dpi=400)  # Adjust dpi as needed
        plt.show(block=False)
        plt.pause(3)
        plt.close()    

#==========================================================
# main
#==========================================================

def main():

    # Create props and globdat dictionaries
    props   = dict()
    globdat = dict()

    # ------- USER MODIFICATION BEGINS HERE! ----------------

    props['basis_dof']      = 'tau'   # dof used for pod basis
    props['sens_dof']       = 'slip'  # dof used for sensitivity
    props['cutoff_tol']     = 1e-3    # cut off tol to choose pod basis
    props['nsteps_online']  = 500     # number of steps in online simulation
    props['jem_export']  = "export JEMDIR='/home/ritukesh/codes/remote/jemjive/jem-3.0'"
    props['jive_export'] = "export JIVEDIR='/home/ritukesh/codes/remote/jemjive/jive-3.0'"
    props['run_falcon']  = 'bash -c "source /opt/intel/oneapi/setvars.sh && /home/ritukesh/codes/github/falcon/build/falcon-opt problem.pro"'

    # ------- USER MODIFICATION ENDS HERE! ------------------

    # Initialize an instance of POD class
    pod = POD(props)

    # Perform offline simulation
    print(f'Running the offline simulation ... \n\n')
    pod.run_offline()

    # Create the pod basis
    print(f'Creating POD basis ... \n\n')
    pod.create_pod_basis()

    # Prepare sensitivity directories
    print(f'Preparing directories to compute sensitivities ... \n\n')
    pod.prepare_sensitivity_directories()

    # Compute sensitivities
    print(f'Computing sensitivities ... \n\n')
    pod.run_sensitivity_simulations()

    # Prepare online directory
    print(f'Preparing online simulation directories ... \n\n')
    pod.prepare_online_directory()

    # Run online simulation
    pod.run_online_simulation()

    # Post-process data
    pod.post_process_data()


if __name__ == "__main__":
    main()