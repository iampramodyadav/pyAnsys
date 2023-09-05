import sys
import numpy as np
import pyvista as pv
import pandas as pd
from ansys.mapdl.core import launch_mapdl

mapdl = launch_mapdl(start_timeout=120,)
print(mapdl)

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
    """This function takes a list of area numbers as input and returns the centroid of the areas and area val.

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

    area_val = mapdl.get(entity='AREA', item1='AREA')
    
    return np.array([cent_x,cent_y,cent_z]), area_val


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
    centroid = np.array(area_cent([a_num])[0])
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
        mapdl.asel("R", "LOC", "Z", -1.0e-6, 1.0e-6)
        # mapdl.allsel(labt="BELOW ", entity="AREA")
        alist_cm = mapdl.geometry.anum
        area_pair[a_num]=alist_cm
        # print(alist_cm)
        # print()
        # mapdl.csys(kcn=0)
    return area_pair

def area21_kp_vector(anum1, anum2):
    """
    Calculates the vector from a keypoint in area 1 to the nearest keypoint in area 2.

    Args:
        anum1 (int): The number of area 1.
        anum2 (int): The number of area 2.

    Returns:
        list: A list of vectors from the nearest keypoint in area 2 to the keypoints in area 1.
        max length vector at first position
    """
    
    mapdl.csys(kcn=0)

    # ------- Area-1 & area-2 centroid -------
    centroid1 = area_cent([anum1])[0]
    centroid2 = area_cent([anum2])[0]
    
    # ------- Area-1 keypoints coordinate -------
    mapdl.asel("S", "AREA", "", anum1)
    mapdl.allsel(labt="BELOW ", entity="AREA")
    mapdl.allsel(labt="BELOW ", entity="LINE")
    llist1 = mapdl.geometry.lnum
    mapdl.allsel(labt="BELOW ", entity="KP")
    klist1 = mapdl.geometry.knum
    len(klist1)

    if len(klist1) <= 40:
        for lnum in llist1:
            mapdl.ldiv(lnum, ndiv=4)
        llist1 = mapdl.geometry.lnum
        klist1 = mapdl.geometry.knum

    kp_coord1 = {}
    for kpoiint in klist1:
        k_x = mapdl.queries.kx(kpoiint)
        k_y = mapdl.queries.ky(kpoiint)
        k_z = mapdl.queries.kz(kpoiint)
        k1 = np.array([k_x, k_y, k_z])
        kp_coord1[kpoiint] = k1

    # ------- Area-2 keypoints coordinate -------
    mapdl.asel("S", "AREA", "", anum2)
    mapdl.allsel(labt="BELOW ", entity="AREA")
    mapdl.allsel(labt="BELOW ", entity="LINE")
    llist2 = mapdl.geometry.lnum
    mapdl.allsel(labt="BELOW ", entity="KP")
    klist2 = mapdl.geometry.knum
    len(klist2)

    if len(klist2) <= 40:
        for lnum in llist2:
            mapdl.ldiv(lnum, ndiv=8)
        llist2 = mapdl.geometry.lnum
        klist2 = mapdl.geometry.knum

    kp_coord2 = {}
    for kpoiint in klist2:
        k_x = mapdl.queries.kx(kpoiint)
        k_y = mapdl.queries.ky(kpoiint)
        k_z = mapdl.queries.kz(kpoiint)
        k1 = np.array([k_x, k_y, k_z])
        kp_coord2[kpoiint] = k1

    # ------- Distance between centroid of area-1 and keypoints of area-2 -------
    kp_list2 = list(kp_coord2.keys())
    dist_list = [np.linalg.norm(centroid1 - kp_coord2[kpoint]) for kpoint in kp_list2]
    # Kp of area-2 nearest to area-1
    min_value = min(dist_list)
    # print(min_value)
    min_index = dist_list.index(min_value)
    # print(min_index)
    
    # ------- Coordinate of nearest Kp -------
    v0 = kp_coord2[kp_list2[min_index]]
    # print(v0)
    
    # ------- Create vector from nearest Kp of area-2 to all Kps of area-1 -------
    kplist1 = list(kp_coord1.keys())
    vector_21_list = [kp_coord1[kpn] - v0 for kpn in kplist1]
    # print(len(vector_21_list))

    # -------max length vector at first -------
    len_vector_21_list = [np.linalg.norm(vec) for vec in vector_21_list]
    max_value = max(len_vector_21_list)
    max_index = len_vector_21_list.index(max_value)
    # len_vector_21_list = len_vector_21_list[max_index: ] + len_vector_21_list[ :max_index]
    # # print(len_vector_21_list)
    # vector_21_list = vector_21_list[max_index: ] + vector_21_list[ :max_index]
    
    return vector_21_list, max_index

def check_overlap(anum1, anum2):
    overlap = 0
    eq_area = 0
    
    centroid1, area_val1 = area_cent([anum1])
    centroid2, area_val2 = area_cent([anum2])
    centroid1 = np.array(centroid1)
    centroid2 = np.array(centroid2)
    
    dist = np.linalg.norm(centroid1 - centroid2)
    
    if dist <= 1.0e-3:
        # print(f"area overlapping, ccentroid dist: {dist}")
        overlap = 1

    if np.isclose(area_val1, area_val2, rtol=1e-05, atol=1e-06, equal_nan=False):
        eq_area = 1

    return overlap, eq_area
    
if __name__ == '__main__':
    paraName = 'BBSPACERS_D4118297-001_SG5XGR01_TEST'
    path = 'K:\\Users\\pramod.kumar\\Sandbox\\PyAnsys'
    cmname= "part2"
    paraimport(paraName, path, cmname)
    mapdl.vplot(vtk=True, show_lines=True, show_axes=True, smooth_shading=True)
