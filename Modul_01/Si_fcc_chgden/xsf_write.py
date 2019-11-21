import numpy as np

def make_data_periodic(x):
    
    Nx, Ny, Nz = x.shape
    
    x_p = np.zeros((Nx+1,Ny+1,Nz+1))

    for i in range(Nx+1):
        if i == Nx:
            ii = 0
        else:
            ii = i
        for j in range(Ny+1):
            if j == Ny:
                jj = 0
            else:
                jj = j
            for k in range(Nz+1):
                if k == Nz:
                    kk = 0
                else:
                    kk = k
                #
                x_p[i,j,k] = x[ii,jj,kk]

    return x_p





def xsf_write(filename, atoms, data=None):
    
    f = open(filename, "w")

    f.write("CRYSTAL\n")
    f.write("PRIMVEC\n")

    # XXX Check this!
    cell = atoms.get_cell().array
    for i in range(3):
        f.write(" {:.14f}  {:.14f}  {:.14f}\n".format(cell[i,0], cell[i,1], cell[i,2]))

    f.write("PRIMCOORD\n")
    Natoms = len(atoms)
    atpos = atoms.get_positions()
    f.write(" {} 1\n".format(Natoms))
    for i in range(Natoms):
        symb = atoms.symbols[0]
        f.write(" {:s}  {:.14f}  {:.14f}  {:.14f}\n".format(symb,atpos[i,0],atpos[i,1],atpos[i,2]))

    if data is None:
        f.close()
        return

    # Assume the system is periodic and the data is only given at Fourier grid
    # We need to make it periodic first (generic grid)
    data_p = make_data_periodic(data)
    Nx, Ny, Nz = data_p.shape

    f.write("BEGIN_BLOCK_DATAGRID_3D\n")
    f.write("data\n")
    f.write("BEGIN_DATAGRID_3Dgrid#1\n")
    
    f.write("{} {} {}\n".format(Nx,Ny,Nz))
    
    f.write("0.0 0.0 0.0\n")  # center
    # write the cell again
    for i in range(3):
        f.write(" {:.14f}  {:.14f}  {:.14f}\n".format(cell[i,0], cell[i,1], cell[i,2]))

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                f.write(" {:18.10f}".format(data_p[i,j,k]))
            f.write("\n")


    f.write("END_DATAGRID_3D\n")
    f.write("END_BLOCK_DATAGRID_3D\n")

    f.close()

    return

