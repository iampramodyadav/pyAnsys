import sys
import numpy as np
import pyvista as pv
import pandas as pd
from ansys.mapdl.core import launch_mapdl
mapdl = launch_mapdl(start_timeout=120,)
# print(mapdl)
PYTHONTOOLSDIR = r"C:\Users\pramod.kumar\OneDrive - SIEMENSGAMESA\Work\ContactDetect"
sys.path.append(PYTHONTOOLSDIR)
import ContactDetection as cd

def paraimport(paraname, parapath, cmname):
    """Import a parametric solid into a APDL and create component.

    Args:
        paraname (str): The name of the parametric solid.
        parapath (str): The path to the parametric solid file.
        cmname (str): The name of the component to create in the APDL

    Returns:
        None.

    """
    mapdl.allsel(labt='ALL', entity='ALL')
    vlist=mapdl.geometry.vnum
    vmin = 0
    if len(vlist)!=0:
        vmin=mapdl.geometry.vnum.max()

    mapdl.parain(name=paraname,
                 extension='x_t',
                 path=parapath,
                 entity='SOLIDS',
                 fmt='0')

    vmax=mapdl.geometry.vnum.max()
    mapdl.vsel(type_='S',vmin=vmin+1, vmax=vmax)
    mapdl.cm(cname=cmname, entity='VOLU')
    mapdl.cmsel(type_='S', name=cmname, entity='VOLU' )
    vlist=mapdl.geometry.vnum
    print(f"{paraname} imported as component {cmname} with vlum num {vlist}")


def area_cent(a_num):
    """This function takes a list of area numbers as input and returns the centroid of the areas.

    The centroid is calculated by first summing the x, y, and z coordinates of all the areas in the list.
    The average of these sums is then returned as the centroid.

    Args:
        a_num (list):list of area numbers

    Returns:
        List of the centroid coordinate x, y, and z 

    """
    for i in range(len(a_num)):
        if i==0:
            mapdl.asel('S','AREA', '', a_num[i])
        else:
            mapdl.asel('A','AREA', '', a_num[i])

    mapdl.asum()
    
    cent_x = mapdl.get(entity='AREA', item1='CENT', it1num='X')
    cent_y = mapdl.get(entity='AREA', item1='CENT', it1num='Y')
    cent_z = mapdl.get(entity='AREA', item1='CENT', it1num='Z')
    
    return np.array([cent_x,cent_y,cent_z])


def clocal_area(a_num):
    """
    This function is to find the centroid and rotation angle in ansys apdl using pyMAPDL library and python library.

    Parameters
    ----------
    a_num : int
        The area number.

    Returns
    -------
    centroid : np.array
        The centroid of the area.
    r_z : float
        The rotation angle about the z-axis.
    r_x : float
        The rotation angle about the x-axis.
    r_y : float
        The rotation angle about the y-axis.

    """
    centroid = np.array(area_cent([a_num]))
    ## selecting kp from area and geneartng kp if required
    # mapdl.csys(kcn=0)
    mapdl.asel("S", "AREA", "", a_num)
    mapdl.allsel(labt="BELOW ", entity="AREA")
    mapdl.allsel(labt="BELOW ", entity="LINE")
    llist = mapdl.geometry.lnum
    # print(llist)
    mapdl.allsel(labt="BELOW ", entity="KP")
    klist = mapdl.geometry.knum
    # print(klist)

    if len(klist) <= 4:
        for lnum in llist:
            mapdl.ldiv(lnum, ndiv=5)
        llist = mapdl.geometry.lnum
        klist = mapdl.geometry.knum

    # print(llist)
    # print(klist)

    ## key point coordinate in dictionary key as kp
    kp_coord = {}
    for kpoiint in klist:
        k_x = mapdl.queries.kx(kpoiint)
        k_y = mapdl.queries.ky(kpoiint)
        k_z = mapdl.queries.kz(kpoiint)
        k1 = np.array([k_x, k_y, k_z])
        kp_coord[kpoiint] = k1

    # k1 = kp_coord[klist[0]]
    # k2 = kp_coord[klist[1]]
    k1 = centroid
    k2 = kp_coord[klist[1]]
    

    for i in range(2, len(klist)):
        k3 = kp_coord[klist[i]]
        if cd.check_collinear(k1, k2, k3) > 1.0e-4:
            # print(klist[i])
            # print(k3)
            break

    x_dcm, y_dcm, z_dcm = cd.coordinate_dcm(k1, k2, k3)
    # R = np.column_stack((x_dcm, y_dcm, z_dcm))
    R = np.row_stack((x_dcm, y_dcm, z_dcm))
    r_z, r_x, r_y = cd.dcm2angleZXY(R)

    return centroid, r_z, r_x, r_y

def contact_between_cm(cm1, cm2, cs_num=1000):
    """
       - This function takes two volume components and returns a dictionary of area pairs.
       - The dictionary key is the area number of the first volume component, and the value is
       - the list of area numbers of the second volume component that are on the same plane as the area number of the first volume component.
        - pyansys and numpy libarary used

    Args:
        cm1 (str): Name of the first volume component.
        cm2 (str): Name of the second volume component.
        cs_num (int, optional): The number of the coordinate system. Defaults to 1000.

    Returns:
        dict: A dictionary of area pairs.
    """
    .

    mapdl.csys(kcn=0)
    mapdl.cmsel(type_="S", name=cm1, entity="VOLU")
    mapdl.allsel(labt="BELOW ", entity="VOLU")
    mapdl.allsel(labt="BELOW ", entity="AREA")
    alist = mapdl.geometry.anum
    cs_num = cs_num

    area_pair = {}
    for a_num in alist:
        mapdl.csys(kcn=0)
        # mapdl.wpcsys(kcn=0)
        centroid, r_z, r_x, r_y = clocal_area(a_num)
        # print(centroid, r_z, r_x, r_y)
        # print(a_num)

        cs_num = cs_num + 1
        cs_typ = 0
        cent_x = centroid[0]
        cent_y = centroid[1]
        cent_z = centroid[2]
        rxy = r_z
        ryz = r_x
        rzx = r_y
        mapdl.clocal(
            kcn=cs_num,
            kcs=cs_typ,
            xl=cent_x,
            yl=cent_y,
            zl=cent_z,
            thxy=rxy,
            thyz=ryz,
            thzx=rzx,
        )

        mapdl.cmsel(type_="S", name=cm2, entity="VOLU")
        mapdl.allsel(labt="BELOW ", entity="VOLU")
        mapdl.allsel(labt="BELOW ", entity="AREA")
        alist1 = mapdl.geometry.anum
        # print(alist1)
        mapdl.wpcsys(kcn=cs_num)
        mapdl.wpstyl(wpctyp=0)
        mapdl.csys(kcn=4)
        mapdl.asel("R", "LOC", "Z", -1.0e-3, 1.0e-3)
        # mapdl.allsel(labt="BELOW ", entity="AREA")
        alist_cm = mapdl.geometry.anum
        area_pair[a_num]=alist_cm
        # print(alist_cm)
        # print()
        # mapdl.csys(kcn=0)
    return area_pair


if __name__ == '__main__':
    paraName = 'BBSPACERS_D4118297-001_SG5XGR01_TEST'
    path = 'K:\\Users\\pramod.kumar\\Sandbox\\PyAnsys'
    cmname= "part2"
    paraimport(paraName, path, cmname)
    mapdl.vplot(vtk=True, show_lines=True, show_axes=True, smooth_shading=True)
