from typing import no_type_check_decorator
import numpy as np

def read_scat_spec(filename="scat_spec.dat", Nth_obs=None, Ne=None):
    """
    Reads the scat_spec.dat file and returns Ispec, Qspec, Uspec, and sorted versions.

    Parameters:
        filename (str): The file to read.
        Nth_obs (int): Number of theta observation points.
        Ne (int): Number of energy bins.

    Returns:
        tuple: (Ispec, Qspec, Uspec, Ispec_s, Qspec_s, Uspec_s)
    """
    if Nth_obs is None or Ne is None:
        raise ValueError("Nth_obs and Ne must be specified.")

    shape = (Nth_obs + 1, Ne + 1)
    shape_s = (6, Nth_obs + 1, Ne + 1)

    print(f"Reading {filename}...")
    print(f"Expecting Ispec, Qspec, Uspec with shape {shape}")
    print(f"Expecting Ispec_s, Qspec_s, Uspec_s with shape {shape_s}")

    with open(filename, "r") as f:
        lines = f.readlines()

    print(f"Total lines in file: {len(lines)}")

    # Read Ispec
    index = 0
    print(f"Reading Ispec from lines {index} to {index + shape[0]}")
    Ispec = np.loadtxt(lines[index:index + shape[0]], dtype=np.float64)
    index += shape[0]

    # Read Qspec
    print(f"Reading Qspec from lines {index} to {index + shape[0]}")
    Qspec = np.loadtxt(lines[index:index + shape[0]], dtype=np.float64)
    index += shape[0]

    # Read Uspec
    print(f"Reading Uspec from lines {index} to {index + shape[0]}")
    Uspec = np.loadtxt(lines[index:index + shape[0]], dtype=np.float64)
    index += shape[0]

    # Read Ispec_s
    print(f"Reading Ispec_s (sorted spectra) from lines {index} onwards")
    Ispec_s = np.array([np.loadtxt(lines[index + i * shape[0]: index + (i + 1) * shape[0]], dtype=np.float64)
                         for i in range(6)])
    index += shape_s[0] * shape[0]

    # Read Qspec_s
    print(f"Reading Qspec_s (sorted spectra) from lines {index} onwards")
    Qspec_s = np.array([np.loadtxt(lines[index + i * shape[0]: index + (i + 1) * shape[0]], dtype=np.float64)
                         for i in range(6)])
    index += shape_s[0] * shape[0]

    # Read Uspec_s
    print(f"Reading Uspec_s (sorted spectra) from lines {index} onwards")
    Uspec_s = np.array([np.loadtxt(lines[index + i * shape[0]: index + (i + 1) * shape[0]], dtype=np.float64)
                         for i in range(6)])

    print("Finished reading file!")
    print(f"Ispec shape: {Ispec.shape}, Qspec shape: {Qspec.shape}, Uspec shape: {Uspec.shape}")
    print(f"Ispec_s shape: {Ispec_s.shape}, Qspec_s shape: {Qspec_s.shape}, Uspec_s shape: {Uspec_s.shape}")

    return Ispec, Qspec, Uspec, Ispec_s, Qspec_s, Uspec_s

def read_scat_line(filename="scat_line.0000.dat", Nth_obs=None, Ne=None):
    """
    Reads the scat_line.0000.dat file and returns Lspec.

    Parameters:
        filename (str): The file to read.
        Nth_obs (int): Number of theta observation points.
        Ne (int): Number of energy bins.

    Returns:
        np.ndarray: Lspec array of shape (Nth_obs + 1, Ne + 1).
    """
    if Nth_obs is None or Ne is None:
        raise ValueError("Nth_obs and Ne must be specified.")

    shape = (Nth_obs + 1, Ne + 1)
    
    print(f"Reading {filename}...")
    print(f"Expecting Lspec with shape {shape}")

    with open(filename, "r") as f:
        lines = f.readlines()

    print(f"Total lines in file: {len(lines)}")

    # Read Lspec
    print(f"Reading Lspec from lines 0 to {shape[0]}")
    Lspec = np.loadtxt(lines[:shape[0]], dtype=np.float64)

    print("Finished reading file!")
    print(f"Lspec shape: {Lspec.shape}")

    return Lspec

def read_scat_disk(filename="scat_disk.0000.dat", Nr=None, Ne=None):
    """
    Reads the scat_disk.0000.dat file and returns Rspec.

    Parameters:
        filename (str): The file to read.
        Nr (int): Number of radial bins.
        Ne (int): Number of energy bins.

    Returns:
        np.ndarray: Rspec array of shape (Nr + 1, Ne + 1).
    """
    if Nr is None or Ne is None:
        raise ValueError("Nr and Ne must be specified.")

    shape = (Nr + 1, Ne + 1)
    
    print(f"Reading {filename}...")
    print(f"Expecting Rspec with shape {shape}")

    with open(filename, "r") as f:
        lines = f.readlines()

    print(f"Total lines in file: {len(lines)}")

    # Read Rspec
    print(f"Reading Rspec from first {shape[0]} lines")
    Rspec = np.loadtxt(lines[:shape[0]], dtype=np.float64)

    print("Finished reading file!")
    print(f"Rspec shape: {Rspec.shape}")

    return Rspec

def read_scat_specr(filename="scat_spcr.0000.dat", Nr=None, Nth_obs=None, Ne=None):
    """
    Reads the scat_spcr.0000.dat file and extracts the spectral components.

    Parameters:
        filename (str): The file to read.
        Nr (int): Number of radial bins.
        Nth_obs (int): Number of observer angles.
        Ne (int): Number of energy bins.

    Returns:
        tuple: (Ispecr, (rr, L_factr, G_factr), Qspecr, Uspecr, Rspecr, Cspecr)
    """
    if Nr is None or Nth_obs is None or Ne is None:
        raise ValueError("Nr, Nth_obs, and Ne must be specified.")

    print(f"Reading {filename}...")
    with open(filename, "r") as f:
        lines = f.readlines()

    line_index = 0

    # Read Ispecr (shape: Nr+1, Nth_obs+1, Ne+1)
    Ispecr = np.zeros((Nr+1, Nth_obs+1, Ne+1))
    print(f"Reading Ispecr with shape {(Nr+1, Nth_obs+1, Ne+1)}")
    
    for jr in range(Nr+1):
        for jth in range(Nth_obs+1):
            Ispecr[jr, jth] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1
        line_index += 1  # Skip blank line

    # Read rr, L_factr, G_factr (shape: Nr+1)
    print("Reading rr, L_factr, G_factr")
    rr = np.zeros(Nr+1)
    L_factr = np.zeros(Nr+1)
    G_factr = np.zeros(Nr+1)

    for jr in range(Nr+1):
        rr[jr], L_factr[jr], G_factr[jr] = map(float, lines[line_index].split())
        # print(rr[jr], L_factr[jr], G_factr[jr])
        line_index += 1
    line_index += 1  # Skip blank line

    # Read Qspecr (shape: Nr+1, Nth_obs+1, Ne+1)
    print(f"Reading Qspecr with shape {(Nr+1, Nth_obs+1, Ne+1)}")
    Qspecr = np.zeros((Nr+1, Nth_obs+1, Ne+1))

    for jr in range(Nr+1):
        for jth in range(Nth_obs+1):
            #print(line_index, lines[line_index])
            Qspecr[jr, jth] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1
        line_index += 1  # Skip blank line
    line_index += 1  # Skip blank line

    # Read Uspecr (shape: Nr+1, Nth_obs+1, Ne+1)
    print(f"Reading Uspecr with shape {(Nr+1, Nth_obs+1, Ne+1)}")
    Uspecr = np.zeros((Nr+1, Nth_obs+1, Ne+1))

    for jr in range(Nr+1):
        for jth in range(Nth_obs+1):
            #print(f"Line {line_index}: {lines[line_index].strip()}")
            Uspecr[jr, jth] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1
        line_index += 1  # Skip blank line
    line_index += 1  # Skip blank line
    # Read Rspecr (shape: Nr+1, Nth_obs+1, Ne+1)
    print(f"Reading Rspecr with shape {(Nr+1, Nth_obs+1, Ne+1)}")
    Rspecr = np.zeros((Nr+1, Nth_obs+1, Ne+1))

    for jr in range(Nr+1):
        for jth in range(Nth_obs+1):
            Rspecr[jr, jth] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1
        line_index += 1  # Skip blank line

    # Read Cspecr (shape: Nr+1, Ne+1)
    print(f"Reading Cspecr with shape {(Nr+1, Ne+1)}")
    Cspecr = np.zeros((Nr+1, Ne+1))

    for jr in range(Nr+1):
        Cspecr[jr] = np.array(list(map(float, lines[line_index].split())))
        line_index += 1

    print("Finished reading file!")
    return Ispecr, (rr, L_factr, G_factr), Qspecr, Uspecr, Rspecr, Cspecr

def read_scat_imag(filename, Nth_obs, Ni, Ne_obs):
    """
    Reads the `scat_imag.XXXX.dat` file and extracts `image` and `spcimage` arrays.

    Parameters:
        filename (str): Path to the data file.
        Nth_obs (int): Number of observation angles.
        Ni (int): Image grid size.
        Ne_obs (int): Number of observed energy bins.

    Returns:
        tuple: (image, spcimage)
    """
    print(f"Reading {filename}...")
    with open(filename, "r") as f:
        lines = f.readlines()
    
    line_index = 0
    total_lines = len(lines)
    print(f"Total lines in file: {total_lines}")
    
    # Read image array (Nth_obs+1, Ni+1, Ni+1)
    image = np.zeros((Nth_obs+1, Ni+1, Ni+1))
    print(f"Reading image array with shape {(Nth_obs+1, Ni+1, Ni+1)}...")
    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            image[it, iy] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1
    
    # Read spcimage array (Nth_obs+1, Ni+1, Ni+1, Ne_obs+1)
    spcimage = np.zeros((Nth_obs+1, Ni+1, Ni+1, Ne_obs+1))
    print(f"Reading spcimage array with shape {(Nth_obs+1, Ni+1, Ni+1, Ne_obs+1)}...")
    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            for ix in range(Ni+1):
                spcimage[it, ix, iy] = np.array(list(map(float, lines[line_index].split())))
                line_index += 1

    print("Finished reading data.")
    return image, spcimage

def read_scat_ipol(filename, Nth_obs, Ni, Ne_obs):
    """
    Reads the `scat_ipol.XXXX.dat` file and extracts `imagex`, `imagey`, `spcimagex`, and `spcimagey`.

    Parameters:
        filename (str): Path to the data file.
        Nth_obs (int): Number of observation angles.
        Ni (int): Image grid size.
        Ne_obs (int): Number of observed energy bins.

    Returns:
        tuple: (imagex, imagey, spcimagex, spcimagey)
    """
    print(f"Reading {filename}...")

    with open(filename, "r") as f:
        lines = f.readlines()

    line_index = 0
    total_lines = len(lines)
    print(f"Total lines in file: {total_lines}")

    # Read imagex array
    imagex = np.zeros((Nth_obs+1, Ni+1, Ni+1))
    print(f"Reading imagex with shape {(Nth_obs+1, Ni+1, Ni+1)}...")

    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            # print(f"Reading line {line_index} for imagex[it={it}, iy={iy}]")
            imagex[it, iy] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1

    # Read imagey array
    imagey = np.zeros((Nth_obs+1, Ni+1, Ni+1))
    print(f"Reading imagey with shape {(Nth_obs+1, Ni+1, Ni+1)}...")

    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            # print(f"Reading line {line_index} for imagey[it={it}, iy={iy}]")
            imagey[it, iy] = np.array(list(map(float, lines[line_index].split())))
            line_index += 1

    # Read spcimagex array
    spcimagex = np.zeros((Nth_obs+1, Ni+1, Ni+1, Ne_obs+1))
    print(f"Reading spcimagex with shape {(Nth_obs+1, Ni+1, Ni+1, Ne_obs+1)}...")

    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            for ix in range(Ni+1):
                # print(f"Reading line {line_index} for spcimagex[it={it}, ix={ix}, iy={iy}]")
                spcimagex[it, ix, iy] = np.array(list(map(float, lines[line_index].split())))
                line_index += 1

    # Read spcimagey array
    spcimagey = np.zeros((Nth_obs+1, Ni+1, Ni+1, Ne_obs+1))
    print(f"Reading spcimagey with shape {(Nth_obs+1, Ni+1, Ni+1, Ne_obs+1)}...")

    for it in range(Nth_obs+1):
        for iy in range(Ni+1):
            for ix in range(Ni+1):
                # print(f"Reading line {line_index} for spcimagey[it={it}, ix={ix}, iy={iy}]")
                spcimagey[it, ix, iy] = np.array(list(map(float, lines[line_index].split())))
                line_index += 1

    print("Finished reading data.")
    return imagex, imagey, spcimagex, spcimagey

def read_scat_cpow(filename, Nr, Nth, Nph):
    """
    Reads the `scat_cpow.XXXX.dat` file and extracts `corpow_ijk`.

    Parameters:
        filename (str): Path to the data file.
        Nr (int): Number of radial points.
        Nth (int): Number of theta (angular) points.
        Nph (int): Number of phi (azimuthal) points.

    Returns:
        corpow_ijk (numpy.ndarray): 3D array of correlation power values.
    """
    print(f"Reading {filename}...")

    with open(filename, "r") as f:
        lines = f.readlines()

    line_index = 0
    total_lines = len(lines)
    print(f"Total lines in file: {total_lines}")

    # Initialize the array
    corpow_ijk = np.zeros((Nr+1, Nth+1, Nph+1))
    print(f"Reading corpow_ijk with shape {(Nr+1, Nth+1, Nph+1)}...")

    for ir in range(Nr+1):
        for ith in range(Nth+1):
            # print(f"Reading line {line_index} for corpow_ijk[ir={ir}, ith={ith}]")
            values = list(map(float, lines[line_index].split()))
            corpow_ijk[ir, ith] = np.array(values)
            line_index += 1

    print("Finished reading data.")
    return corpow_ijk

def read_scat_ithp(filename, Nth_obs, Nph_obs, Ni):
    """
    Reads the `scat_ithp.XXXX.dat` file and extracts `phimage`, `phimagex`, `phimagey`.

    Parameters:
        filename (str): Path to the data file.
        Nth_obs (int): Number of theta observation points.
        Nph_obs (int): Number of phi observation points.
        Ni (int): Image resolution (assumed square).

    Returns:
        tuple: (phimage, phimagex, phimagey) as 4D numpy arrays.
    """
    print(f"Reading {filename}...")

    with open(filename, "r") as f:
        lines = f.readlines()

    line_index = 0
    total_lines = len(lines)
    print(f"Total lines in file: {total_lines}")

    # Initialize arrays
    phimage = np.zeros((Nth_obs+1, Nph_obs+1, Ni+1, Ni+1))
    phimagex = np.zeros((Nth_obs+1, Nph_obs+1, Ni+1, Ni+1))
    phimagey = np.zeros((Nth_obs+1, Nph_obs+1, Ni+1, Ni+1))
    
    print(f"Reading phimage, phimagex, phimagey with shape {(Nth_obs+1, Nph_obs+1, Ni+1, Ni+1)}...")

    def read_image(array_name, array):
        nonlocal line_index
        for it in range(Nth_obs+1):
            for ip in range(Nph_obs+1):
                for ix in range(Ni+1):
                    values = list(map(float, lines[line_index].split()))
                    array[it, ip, ix] = np.array(values)
                    line_index += 1

                # Extra newline after each ip loop
                if line_index < total_lines and lines[line_index].strip() == "":
                    line_index += 1

            # Extra newline after each it loop
            if line_index < total_lines and lines[line_index].strip() == "":
                line_index += 1

        print(f"Finished reading {array_name}.")

    read_image("phimage", phimage)
    read_image("phimagex", phimagex)
    read_image("phimagey", phimagey)

    print("Finished reading all data.")
    return phimage, phimagex, phimagey

def read_scat_spcp(filename, Nth_obs, Nph_obs, Ne):
    """
    Reads the `scat_spcp.XXXX.dat` file and extracts `Ispecp`.

    Parameters:
        filename (str): Path to the data file.
        Nth_obs (int): Number of theta observation points.
        Nph_obs (int): Number of phi observation points.
        Ne (int): Number of energy bins.

    Returns:
        np.ndarray: Ispecp array with shape (Nth_obs+1, Nph_obs+1, Ne+1).
    """
    print(f"Reading {filename}...")

    with open(filename, "r") as f:
        lines = f.readlines()

    line_index = 0
    total_lines = len(lines)
    print(f"Total lines in file: {total_lines}")

    # Initialize array
    Ispecp = np.zeros((Nth_obs+1, Nph_obs+1, Ne+1))

    print(f"Reading Ispecp with shape {(Nth_obs+1, Nph_obs+1, Ne+1)}...")

    for jth in range(Nth_obs+1):
        for jph in range(Nph_obs+1):
            values = list(map(float, lines[line_index].split()))
            Ispecp[jth, jph, :] = np.array(values)
            line_index += 1

        # Extra newline after each jph loop
        if line_index < total_lines and lines[line_index].strip() == "":
            line_index += 1

    print("Finished reading Ispecp.")
    return Ispecp




if __name__ == "__main__":
    Nth_obs = 40
    Ne = 100
    Nr = 95
    Ni = 80
    Ne_obs = 10
    Nth=63
    Nph=47
    Nth_obs=40
    Nph_obs=0
    print('====')
    Ispec, Qspec, Uspec, Ispec_s, Qspec_s, Uspec_s = read_scat_spec(filename='scat_spec.0624.dat', Nth_obs=40, Ne=Ne)
    print('====')
    Lspec = read_scat_line(filename='scat_line.0624.dat', Nth_obs=Nth_obs, Ne=Ne)
    print('====')
    Rspec = read_scat_disk(filename='scat_disk.0624.dat', Nr=Nr, Ne=Ne)
    print('====')
    Ispecr, (rr, L_factr, G_factr), Qspecr, Uspecr, Rspecr, Cspecr = read_scat_specr(filename="scat_spcr.0624.dat", Nr=Nr, Nth_obs=Nth_obs, Ne=Ne)
    print('====')
    image, spcimage = read_scat_imag(filename='scat_imag.0624.dat', Nth_obs=Nth_obs, Ni=Ni, Ne_obs=Ne_obs)
    print('====')
    read_scat_ipol(filename='scat_ipol.0624.dat', Nth_obs=Nth_obs, Ni=Ni, Ne_obs=Ne_obs)
    print('====')
    read_scat_cpow(filename='scat_cpow.0624.dat', Nr=Nr, Nth=Nth, Nph=Nph)
    print('====')
    phimage, phimagex, phimagey = read_scat_ithp(filename='scat_ithp.0624.dat', Nth_obs=Nth_obs, Nph_obs=Nph_obs, Ni=Ni)
    print('====')
    Ispecp = read_scat_spcp(filename='scat_spcp.0624.dat', Nth_obs=Nth_obs, Nph_obs=Nph_obs, Ne=Ne)
