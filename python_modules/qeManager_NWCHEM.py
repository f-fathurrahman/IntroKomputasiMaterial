from ase import Atoms
import numpy as np

ANG2BOHR = 1.889725989

def read_nwchem_out_version(filename):
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        # Search for string that describes the version
        if "Northwest Computational Chemistry Package (NWChem)" in line:
            version = line.split()[-1]
            break
    f.close()
    return version

def read_nwchem_out_Natoms(filename):
    Natoms = 0
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        # Search for string that describes the version
        if "XYZ format geometry" in line:
            line = f.readline() # skip
            Natoms = int( f.readline() )
            break
    f.close()
    return Natoms


def print_traj_xyz_forces(traj, traj_forces, filename):
    f = open(filename, "w")
    Ntraj = len(traj)
    for i in range(Ntraj):
        atoms = traj[i]
        forces = traj_forces[i]
        Natoms = len(atoms)
        symbs = atoms.symbols
        r = atoms.get_positions()
        f.write("%d\n" % Natoms)
        f.write("\n")
        for ia in range(Natoms):
            f.write("%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n" % \
                (symbs[ia], r[ia,0], r[ia,1], r[ia,2], \
                    forces[ia][0], forces[ia][1], forces[ia][2]))
    f.close()


def read_nwchem_out_geoopt(filename):

    Natoms = read_nwchem_out_Natoms(filename)
    
    f = open(filename, "r")
    
    first_time = True

    traj = []
    traj_forces = []

    ifound_geom = 0
    ifound_grad = 0

    while True:
        
        line = f.readline()
        
        if not line:
            break

        if "Geometry \"geometry\" -> \"geometry\"" in line:
            
            ifound_geom = ifound_geom + 1
            
            #print("Found geometry")
            #print(line)

            # Skip lines
            for i in range(6):
                line = f.readline()

            symbs = []
            coords = []
            for ia in range(Natoms):
                line = f.readline()
                #print(line)
                ll = line.split()
                symbs.append(ll[1])
                # FIXME: Check this
                x = float( ll[3] )
                y = float( ll[4] )
                z = float( ll[5] )
                coords.append([x,y,z])
            traj.append( Atoms(symbs, coords) )

        if "DFT ENERGY GRADIENTS" in line:
            ifound_grad = ifound_grad + 1
            # Skip lines
            for i in range(3):
                line = f.readline()
            forces = []
            for ia in range(Natoms):
                line = f.readline()
                #print(line)
                ll = line.split()
                fx = float( ll[5] )/ANG2BOHR
                fy = float( ll[6] )/ANG2BOHR
                fz = float( ll[7] )/ANG2BOHR
                forces.append([fx,fy,fz])
            traj_forces.append(forces)

    # Append zeros forces for the last atoms
    forces = []
    for ia in range(Natoms):
        forces.append([0.0, 0.0, 0.0])
    traj_forces.append(forces)

    f.close()
    print("ifound_geom = ", ifound_geom)
    print("ifound_grad = ", ifound_grad)
    print_traj_xyz_forces(traj, traj_forces, "TEMP_GEOOPT.xyz")
