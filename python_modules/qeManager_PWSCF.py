#
# This file is originally part of qeManager by me (Fadjar Fathurrahman).
# Now the project seems to be abandoned in favor of other projects.
#
# This file is intended to help managing and automating PWSCF calculations.
# To avoid installing qeManager, I combined several files into one big file
# like this which is meant to be imported directly. 
#
# TODO: Documentation
#
# `PWSCFInput` is meant to be a low level PWSCF input manipulator. It consists
# of lower "subclasses" such as `ControlNameList`, etc. All units conform
# with the units used in PWSCF.


import sys
import os

import numpy as np
import matplotlib.pyplot as plt

from ase import Atom, Atoms
from ase.units import Bohr, Ry
from ase.dft.kpoints import *
import ase.io


class ControlNameList:
    """
    Values equal to `None` will not be written in the input file.
    """
    def __init__(self):

        self.calculation = "scf"
        self.nstep = 1000
        self.tprnfor = True
        self.outdir = "./tmp"
        self.prefix = "pwscf"
        self.wf_collect = True
        self.restart_mode = "from_scratch"
        self.pseudo_dir = "./"
        self.verbosity = "high"


        self.title = None
        self.iprint = None
        self.tstress = None
        self.dt = None
        self.wfcdir = None
        self.lkpoint_dir = None
        self.max_seconds = None
        self.etot_conv_thr = None
        self.forc_conv_thr = None
        self.disk_io = None
        self.tefield = None
        self.dipfield = None
        self.lelfield = None
        self.nberrycyc = None
        self.lorbm = None
        self.lberry = None
        self.gdir = None
        self.nppstr = None
        self.lfcpopt = None
        self.monopole = None


    def write(self,f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&CONTROL\n")
        f.write("  prefix = \"%s\"\n" % self.prefix)
        f.write("  outdir = \"%s\"\"\n" % self.outdir)
        f.write("  pseudo_dir = \"%s\"\n" % self.outdir)
        f.write("  tprnfor = %s\n" % self.tprnfor)
        f.write("  nstep = %d\n" % self.nstep)
        f.write("/\n\n")

    def write_all(self,f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&CONTROL\n")
        sdict = self.__dict__
        for k in sdict:
            if not( sdict[k] == None ):
                if type(sdict[k]) == str:
                    f.write("  %s = \"%s\"\n" % (k,sdict[k]))
                else:
                    f.write("  %s = %s\n" % (k,sdict[k]))
        f.write("/\n\n")






class SystemNameList:

    """
    A class to represent SYSTEM namelist in PWSCF input.
    Atomic structures will be specified using in CELL_PARAMETERS.
    """
    def __init__(self, atoms):

        self.ibrav = 0
        self.celldm = None
        self.A = None
        self.B = None
        self.C = None
        self.cosAB = None
        self.cosAC = None
        self.cosBC = None
        self.nat = len(atoms)

        self.ntyp = len( np.unique( atoms.get_atomic_numbers() ) )

        self.nbnd = None
        self.tot_charge = None
        self.tot_magnetization = None
        self.starting_magnetization = None
        self.ecutwfc = 30.0
        self.ecutrho = 4.*self.ecutwfc
        self.ecutfock = None
        self.nr1 = None
        self.nr2 = None
        self.nr3 = None
        self.nr1s = None
        self.nr2s = None
        self.nr3s = None
        self.nosym = None
        self.nosym_evc = None

        self.noinv = None
        self.no_t_rev = None
        self.force_symmorphic = None
        self.use_all_frac = None
        self.occupations = None
        self.one_atom_occupations = None
        self.starting_spin_angle = None
        self.degauss = None
        self.smearing = None
        self.nspin = None
        self.noncolin = None
        self.ecfixed = None
        self.qcutz = None
        self.q2sigma = None
        self.input_dft = None

        self.exx_fraction = None
        self.screening_parameter = None
        self.exxdiv_treatment = None
        self.x_gamma_extrapolation = None
        self.ecutvcut = None
        self.nqx1 = None
        self.nqx2 = None
        self.nqx3 = None
        self.lda_plus_u = None
        self.lda_plus_u_kind = None
        self.Hubbard_U = None
        self.Hubbard_J0 = None
        self.Hubbard_alpha = None
        self.Hubbard_beta = None
        #Hubbard_J(i,ityp)
        #starting_ns_eigenvalue(m,ispin,I)
        self.U_projection_type = None
        self.edir = None
        self.emaxpos = None
        self.eopreg = None
        self.eamp = None
        self.angle1 = None
        self.angle2 = None
        self.constrained_magnetization = None
        self.fixed_magnetization = None
        #lambda
        self.report = None
        self.lspinorb = None
        self.assume_isolated = None
        self.esm_bc = None
        self.esm_w = None
        self.esm_efield = None
        self.esm_nfit = None
        self.fcp_mu = None
        self.vdw_corr = None
        self.london = None
        self.london_s6 = None
        self.london_c6 = None
        self.london_rvdw = None
        self.london_rcut = None
        self.ts_vdw_econv_thr = None
        self.ts_vdw_isolated = None
        self.xdm = None
        self.xdm_a1 = None
        self.xdm_a2 = None
        self.space_group = None
        self.uniqueb = None
        self.origin_choice = None
        self.rhombohedral = None
        self.zmon = None
        self.realxz = None
        self.block = None
        self.block_1 = None
        self.block_2 = None
        self.block_height = None

    def set_ecutwfc(self,ecutwfc):
        self.ecutwfc = ecutwfc
        self.ecutrho = 4.0*ecutwfc

    def set_spinpolarized(self):
        self.nspin = 2
        self.starting_magnetization = []
        for i in range(self.ntyp):
            self.starting_magnetization.append(0.5)

    def set_smearing(self, smearing_type="mv", degauss=0.01):
        self.occupations = "smearing"
        self.smearing = smearing_type
        self.degauss = degauss

    def write(self, f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&SYSTEM\n")
        f.write("  ibrav = %d\n" % self.ibrav)
        f.write("  nat = %d\n" % self.nat)
        f.write("  ntyp = %d\n" % self.ntyp)
        f.write("  ecutwfc = %f\n" % self.ecutwfc)
        f.write("  ecutrho = %f\n" % self.ecutrho)
        f.write("/\n\n")

    def write_all(self,f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&SYSTEM\n")
        sdict = self.__dict__
        for k in sdict:
            if not( sdict[k] == None ):
                if type(sdict[k]) == str:
                    f.write("  %s = \"%s\"\n" % (k,sdict[k]))
                elif type(sdict[k]) == list:
                    Nlist = len(sdict[k])
                    for i in range(Nlist):
                        f.write("  %s(%d) = %f\n" % (k,i+1,sdict[k][i]))
                else:
                    f.write("  %s = %s\n" % (k,sdict[k]))
        f.write("/\n\n")




class ElectronsNameList:
    """
    A class to represent ELECTRONS namelist in PWSCF input.
    """

    def __init__(self):
        self.electron_maxstep = 150
        self.conv_thr = 1e-6
        self.mixing_mode = "plain"
        self.mixing_beta = 0.7
        self.mixing_ndim = 8
        self.diagonalization = "david"
        self.scf_must_converge = None
        self.adaptive_thr = None
        self.conv_thr_init = None
        self.conv_thr_multi = None
        self.mixing_fixed_ns = None
        self.ortho_para = None
        self.diago_thr_init = None
        self.diago_cg_maxiter = None
        self.diago_david_ndim = None
        self.diago_full_acc = None
        self.efield = None
        self.efield_cart = None
        self.efield_phase = None
        self.startingpot = None
        self.startingwfc = None
        self.tqr = None

    def write_all(self,f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&ELECTRONS\n")
        sdict = self.__dict__
        for k in sdict:
            if not( sdict[k] == None ):
                if type(sdict[k]) == str:
                    f.write("  %s = \"%s\"\n" % (k,sdict[k]))
                else:
                    f.write("  %s = %s\n" % (k,sdict[k]))
        f.write("/\n\n")



class IonsNameList:
    """
    A class to represent IONS namelist in PWSCF
    """

    def __init__(self):
         self.ion_dynamics = None
         self.ion_positions = None
         self.pot_extrapolation = None
         self.wfc_extrapolation = None
         self.remove_rigid_rot = None
         self.ion_temperature = None
         self.tempw = None
         self.tolp = None
         self.delta_t = None
         self.nraise = None
         self.refold_pos = None
         self.upscale = None
         self.bfgs_ndim = None
         self.trust_radius_max = None
         self.trust_radius_min = None
         self.trust_radius_ini = None
         self.w_1 = None
         self.w_2 = None

    def write_all(self,f=None):
        if f == None:
            f = sys.stdout
        #
        f.write("&IONS\n")
        sdict = self.__dict__
        for k in sdict:
            if not( sdict[k] == None ):
                if type(sdict[k]) == str:
                    f.write("  %s = \"%s\"\n" % (k,sdict[k]))
                else:
                    f.write("  %s = %s\n" % (k,sdict[k]))
        f.write("/\n\n")


class PWSCFBandstructure:

    def __init__(self,pwinput):
        self.pwinput = pwinput


    def write_scf(self):
        self.pwinput.filename = "PWINPUT_scf"
        self.pwinput.write()

    def write_bands(self):
        self.pwinput.set_calc_bands("fcc",Nkpts=100)
        self.pwinput.filename = "PWINPUT_bands"
        self.pwinput.write()

    def write_collect_bands(self):
        write_bands_inp()

    def run(self):
        pass

    def plot_bands(self, use_collect_bands_data=True,
                   shift_to_efermi=False, ):
        #
        Nbands, Nkpts = get_Nbands_Nkpts("LOG_bands")
        print("Nbands = ", Nbands, "Nkpts = ", Nkpts)

        if not use_collect_bands_data:
            xcoords = pwinput.bands_xcoords
            ebands = read_bandstructure("LOG_bands")
        else:
            xcoords, ebands = read_bands_gnu("bands.out.gnu", Nbands, Nkpts)

        Efermi = read_fermi_level("LOG_scf")

        if shift_to_efermi:
            ebands[:,:] = ebands[:,:] - Efermi

        #Emin = np.min(ebands)
        #Emax = np.max(ebands)
        Emin = 10
        Emax = 30

        print("Emin   = ", Emin)
        print("Emax   = ", Emax)
        print("Efermi = ", Efermi)

        plt.clf()
        for ib in range(Nbands):
            plt.plot( xcoords[ib,:], ebands[ib,:], marker="o" )

        plt.grid()
        plt.ylim(Emin,Emax)

        #specialk_xcoords = pwinput.specialk_xcoords
        specialk_xcoords = read_specialk_xcoords("LOG_collect_bands")

        kmin = np.min(specialk_xcoords)
        kmax = np.max(specialk_xcoords)

        plt.xlim(kmin,kmax)

        for p in specialk_xcoords:
          plt.plot([p, p], [Emin, Emax], "k-")

        plt.plot([kmin,kmax],[Efermi,Efermi], "k-", linewidth=2)

        plt.savefig("bands.png", dpi=300)


class PWSCFInput:

    """
    A simple Python class to generate PWSCF input
    """

    def __init__(self, atoms, pspFiles, filename=None, move_atoms=False,
                 gamma_only=False, kpt_automatic=False, Nk=[1,1,1], nkshift=[0,0,0]):

        self.atoms = atoms
        self.pspFiles = pspFiles
        if filename == None:
            self.filename = sys.stdout
        else:
            self.filename = filename
        self.move_atoms = move_atoms
        self.gamma_only = gamma_only
        self.kpt_automatic = kpt_automatic
        self.Nk = Nk
        self.nkshift = nkshift
        #
        self.CONTROL = ControlNameList()
        self.SYSTEM = SystemNameList(atoms)
        self.ELECTRONS = ElectronsNameList()
        self.IONS = IonsNameList()
        # for bandstructure
        self.bandstructure_mode = False
        self.kpts = None
        self.kweights = None
        self.special_kcoords = None


    def set_spinpolarized(self):
        self.SYSTEM.set_spinpolarized()

    def set_ecutwfc(self,ecutwfc):
        self.SYSTEM.set_ecutwfc(ecutwfc)

    def set_smearing(self,smearing_type="mv",degauss=0.01):
        self.SYSTEM.set_smearing(smearing_type=smearing_type,degauss=degauss)

    def set_calc_bands(self, lattice, Nkpts=60):
        self.bandstructure_mode = True
        self.kpt_automatic = False
        self.gamma_only = False
        self.CONTROL.calculation = "bands"
        self.CONTROL.verbosity = "high"
        #
        kpts, x, Xkpt = gen_kpath( self.atoms, lattice, Nkpts=Nkpts )
        Nkpts = len(x)
        print("Nkpts = ", Nkpts)
        for ik in range(Nkpts):
            sys.stdout.write("%.8f %.8f %.8f %.8f\n" % (kpts[ik,0],kpts[ik,1],kpts[ik,2],x[ik]))
        print(Xkpt)
        #
        self.bands_xcoords = x
        self.kpts = kpts
        self.kweights = np.ones(len(x))
        self.specialk_xcoords = Xkpt


    def write(self):
        #
        inpFile = open(self.filename,"w")
        #
        self.CONTROL.write_all(f=inpFile)
        self.SYSTEM.write_all(f=inpFile)
        self.ELECTRONS.write_all(f=inpFile)
        if self.move_atoms:
            self.IONS.write_all(f=inpFile)
        write_atomic_species( self.atoms, pspFiles=self.pspFiles, f=inpFile )
        write_atomic_positions( self.atoms, f=inpFile )
        #
        if self.bandstructure_mode:
            write_kpoints( f=inpFile, gamma_only=self.gamma_only, automatic=self.kpt_automatic,
                       Nk=self.Nk, nkshift=self.nkshift,
                       bandstructure=True,
                       kpts=self.kpts, weights=self.kweights )
        else:
            write_kpoints( f=inpFile, gamma_only=self.gamma_only, automatic=self.kpt_automatic,
                       Nk=self.Nk, nkshift=self.nkshift )
        #
        write_cell( self.atoms, f=inpFile )
        #
        inpFile.close()


class PWSCFOutput:

    def __init__(self, fname):
        self.filename = fname
        self.file = open(fname, "r")

    def parse(self):
        f = self.file
        line = f.readline()
        while line:
            #
            if "Program PWSCF" in line:
                self.version = line.split()[2]
                self.start_date = line.split()[8]
                # work-around for Python3
                self.start_time = "".join( line.split()[10:] ).replace(" ","")
                print("Version = %s" % self.version)
                print("Start date = %s" % self.start_date)
                print("Start time = %s" % self.start_time)
            #
            line = f.readline()

    def close(self):
        self.file.close()




def read_fermi_level(filename):
    f = open(filename,"r")
    while True:
        line = f.readline()
        if not line:
            break
        #
        if("the Fermi energy" in line):
            Efermi = float( line.split()[4] )
    return Efermi

def read_specialk_xcoords(filename):
    f = open(filename,"r")
    xcoords = []
    while True:
        line = f.readline()
        if not line:
            break
        #
        if("high-symmetry" in line):
            xcoords.append( float(line.split()[7]) )
    f.close()
    return np.array(xcoords)


def read_bands_gnu(filband, Nbands, Nkpts ):
    databands = np.loadtxt(filband)

    ebands = np.zeros( (Nbands, Nkpts) )
    kvec   = np.zeros( (Nbands, Nkpts) )

    for ib in range(Nbands):
        idx1 = (ib)*Nkpts
        idx2 = (ib+1)*Nkpts
        # For QE-6 no need to convert the band energy from Ry to eV
        ebands[ib,:] = databands[idx1:idx2,1]
        kvec[ib,:]   = databands[idx1:idx2,0]
    #
    return kvec, ebands


def get_Nbands_Nkpts(logfile):
    f = open(logfile,"r")
    #
    while True:
        line = f.readline()
        if not line:
            break
        # Read number of bands
        if("number of Kohn-Sham states" in line):
            Nbands = int( line.split()[4] )
        # read number of k-points
        if("number of k points" in line):
            Nkpts = int( line.split()[4] )
    #
    f.close()
    return Nbands, Nkpts



def read_bandstructure(logfile):
    #
    Nbands, Nkpts = get_Nbands_Nkpts(logfile)
    #
    ebands = np.zeros((Nbands,Nkpts))
    #
    ENE_PER_LINE = 8
    #
    f = open(logfile,"r")
    #
    while True:
        line = f.readline()
        if not line:
            break
        #
        Nline = Nbands/ENE_PER_LINE
        Nextra = Nbands%ENE_PER_LINE
        #
        if("End of band structure calculation" in line):
            #ikpt = 0
            for ikpt in range(Nkpts):
                f.readline()
                f.readline()
                f.readline()
                ibnd = 0
                for i in range(Nline):
                    line = f.readline()
                    #print(line,end="")
                    for i in range(8):
                        ebands[ibnd,ikpt] = float(line.split()[i])
                        ibnd = ibnd + 1
                #
                if Nextra > 0:
                    line = f.readline()
                    #print(line,end="")
                    for i in range(Nextra):
                        ebands[ibnd,ikpt] = float(line.split()[i])
                        ibnd  = ibnd + 1
                #
                #print(ebands[:,ikpt])
    #
    print("Nkpts = ", Nkpts, " Nbands = ", Nbands)
    f.close()
    return ebands


# If unit=1 no conversion is done
def read_pwscf_energy(logfile,last=True,unit=1):
    """
    Read PWSCF energy
    """
    f = open(logfile,"r")
    energy = []
    while True:
        line = f.readline()
        if not line:
            break
        if("!    total energy" in line):
            #sys.stdout.write(line)
            energy.append( float(line.split()[4]) )
    #
    f.close()
    #
    if last:
        return energy[-1]*unit
    else:
        return energy*unit


#
# FIXME: This is not working when ControlNamelist.verbose = "high"
#
def read_pwscf_force(logfile,last=True):
    print("Reading force for file: ", logfile)
    f = open(logfile,"r")
    all_force = []

    while True:
        line = f.readline()
        # break if EOF
        if not line:
            break
        # Read number of atoms
        if("number of atoms" in line):
            #print line
            Natoms = int( line.split()[4] )
        # Read force
        if("Forces acting" in line):
            line = f.readline()
            #
            force = np.zeros( (Natoms,3) )
            for ia in range(Natoms):
                line = f.readline()
                #print line
                force[ia,0] = float( line.split()[6] )
                force[ia,1] = float( line.split()[7] )
                force[ia,2] = float( line.split()[8] )
                #
                #print()
                #print(force*Ry/Bohr)
                all_force.append(force)
    #
    f.close()
    #
    if last:
        return all_force[-1]*(Ry/Bohr)
    else:
        return all_force*(Ry/Bohr)


def read_pwscf_total_force(logfile,last=True):

    f = open(logfile,"r")
    
    total_force = None

    while True:
        line = f.readline()
        # break if EOF
        if not line:
            break
        # Read total force
        if("Total force" in line):
            print("line = ", line)
            total_force = float( line.split()[3] )
    f.close()
    #
    return total_force


def write_bands_inp( filename="bands.inp", tmpdir="./tmp"):
    #
    f = open(filename, "w")
    f.write("&BANDS\n")
    f.write("  outdir = \"%s\"\n" % tmpdir)
    f.write("/\n")
    f.close()


def write_kpoints(f=None, kpts=None, weights=None,
                  Nk=[1,1,1], nkshift=[0,0,0],
                  gamma_only=False,
                  automatic=False,
                  bandstructure=False):
    if gamma_only:
        f.write("K_POINTS gamma\n")
        f.write("\n")
    elif automatic:
        f.write("K_POINTS automatic\n")
        f.write("%d %d %d %d %d %d\n\n" % (Nk[0], Nk[1], Nk[2],
                 nkshift[0], nkshift[1], nkshift[2]))
    elif bandstructure:
        f.write("K_POINTS crystal\n")
        Nkpts = kpts.shape[0]
        f.write("%d\n" % Nkpts)
        for ik in range(Nkpts):
            f.write("%.9f %.9f %.9f %.9f\n" % (kpts[ik,0], kpts[ik,1], kpts[ik,2], weights[ik]))
        f.write("\n")
    else:
        f.write("K_POINTS\n")
        Nkpts = kpts.shape[0]
        f.write("%d\n" % Nkpts)
        for ik in range(Nkpts):
            f.write("%.9f %.9f %.9f %.9f\n" % (kpts[ik,0], kpts[ik,1], kpts[ik,2], weights[ik]))
        f.write("\n")


def write_atomic_species(atoms, f=None, pspFiles=None, masses=None):
    if f==None:
        f = sys.stdout
    #
    unique_species = np.unique(atoms.get_chemical_symbols())

    f.write("ATOMIC_SPECIES\n")
    #
    if pspFiles is not None:
        isp = 0
        for s in unique_species:
            f.write("%5s %8.2f  %s\n" % (s, Atom(s).mass, pspFiles[isp]))
            isp = isp + 1
    else:
        raise RuntimeError("pspFiles needs to be specified")
    #
    f.write("\n")




def write_atomic_positions(atoms, f=None):
    if f == None:
        f = sys.stdout
    #
    f.write("ATOMIC_POSITIONS angstrom\n")
    for a in atoms:
        f.write("%5s %18.10f %18.10f %18.10f\n" %
                (a.symbol, a.position[0], a.position[1], a.position[2]))
    f.write("\n")

def write_cell(atoms,f=None):
    cell = atoms.get_cell()/Bohr # convert to Bohr
    if f == None:
        f = sys.stdout
    #
    f.write("CELL_PARAMETERS bohr\n")
    f.write("%18.10f %18.10f %18.10f\n" % (cell[0,0],cell[0,1],cell[0,2]) )
    f.write("%18.10f %18.10f %18.10f\n" % (cell[1,0],cell[1,1],cell[1,2]) )
    f.write("%18.10f %18.10f %18.10f\n" % (cell[2,0],cell[2,1],cell[2,2]) )
    f.write("\n")


def find_ntyp(atoms):
    """
    Return number of unique atom types (species)
    Not really needed, he he he :-)
    np.unique can be used instead
    """
    symbs = atoms.get_chemical_symbols()
    unique_symbols = []
    for s in symbs:
        if not (s in unique_symbols):
            unique_symbols.append(s)
    #
    return len(unique_symbols)


def gen_kpath( atoms, lattice, Nkpts=60 ):
    """
    Automatically generate k-points for a band structure calculation.
    """
    #
    points = get_special_points(atoms.cell, lattice)
    paths = parse_path_string(special_paths[lattice])
    #print(paths[0])
    kpts_spec = [points[k] for k in paths[0]]
    kpts, x, Xkpt = get_bandpath(kpts_spec,atoms.cell,Nkpts)
    #
    # TODO: also return string for special k-points" symbol
    # probably using variable `paths`.
    return kpts, x, Xkpt



class ConvergenceTest:
    """
    A class for doing 1D convergence test.
    """

    def __init__(self, pwinput, what=None, values=None, prefixInp='TEMP_PWINPUT_', prefixOut='LOG_'):
        """
        `what` can be one of ecutwfc, kpts
        `values`
        """
        #
        if what == None and values==None:
            raise RuntimeWarning('Test convergence is set to ecutwfc')
            self.what = 'ecutwfc'
            self.values = np.arange(20,80,10)
        #
        self.pwinput = pwinput
        self.what = what
        self.values = values
        self.Ndata = len(values)
        #
        self.energies = None
        #
        self.prefixInp = prefixInp  # XXX need to save this ?
        self.prefixOut = prefixOut
        self.inpFiles = []
        self.outFiles = []
        if self.what == 'kpoints':
            for v in values:
                strv = str(v).strip('[]').replace(',','').replace(' ','_')
                self.inpFiles.append(self.prefixInp + what + '_' + strv)
                self.outFiles.append(self.prefixOut + what + '_' + strv)
        else:
            for v in values:
                self.inpFiles.append(self.prefixInp + what + '_' + str(v))
                self.outFiles.append(self.prefixOut + what + '_' + str(v))
        #
        self.inputs_have_been_written = False


    def run(self):
        """
        one-time run
        """
        #
        if not self.inputs_have_been_written:
            self.write()
        #
        for i in range(self.Ndata):
            os.system('pw.x < ' + self.inpFiles[i] + ' > ' + self.outFiles[i])


    def write(self):
        """
        writes only the required input files
        """
        #
        if self.what == 'ecutwfc':
            for i in range(self.Ndata):
                self.pwinput.filename = self.inpFiles[i]
                self.pwinput.SYSTEM.set_ecutwfc(self.values[i])
                self.pwinput.write()
        #
        elif self.what == 'ecutrho':
            for i in range(self.Ndata):
                self.pwinput.filename = self.inpFiles[i]
                self.pwinput.SYSTEM.ecutrho = self.values[i]
                self.pwinput.write()
        elif self.what == 'kpoints':
            for i in range(self.Ndata):
                self.pwinput.filename = self.inpFiles[i]
                self.pwinput.Nk = self.values[i]
                self.pwinput.write()
        #
        else:
            raise RuntimeError('what = %s is not implemented yet' % (self.what))
        #
        self.inputs_have_been_written = True


    def read(self):
        """
        Read data
        """
        energies = []
        total_forces = []
        #forces = []
        for i in range(self.Ndata):
            energies.append( read_pwscf_energy(self.outFiles[i]) )
            #forces.append( read_pwscf_force(self.outFiles[i]) )
            total_forces.append( read_pwscf_total_force(self.outFiles[i]) )
        self.energies = np.array(energies)
        #self.forces = forces
        self.total_forces = total_forces
        return self.values, self.energies, self.total_forces
