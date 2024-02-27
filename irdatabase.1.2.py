#! python3

"""
Script to fill a csv file for the SBNAF IR database from all datafiles.

With this setup it only works if the id is the id number of the asteroid.

The script does not check if a parameter is valid or not, so in this state
it is vulnerable to user error.

All data and a header is written to the output file in csv format.
R. Szakats, 2018-2021, Konkoly Obs.
V1.3
"""
import sys
import os
import csv
from numpy import array
import numpy as np
import datetime
import time

# import subprocess
from BB_cc_Wright2010_call import ccorrect_wise
from BB_cc_AKARI_call import ccorrect_akari
from BB_cc_MSX_call import ccorrect_msx
from BB_cc_IRAS_call import ccorrect_simps
import wget
from decimal import Decimal

# import math
from astropy.table import Table, setdiff, Column
from progress import progress
import warnings
import urllib3

# import re
import getopt
import logging

warnings.simplefilter(action="ignore", category=FutureWarning)
# warnings.simplefilter(action='ignore', category=ssl-warnings)

urllib3.disable_warnings()


###############################################################################


def get_jplh(target, center, tlist, ismoon):
    # API documentation:
    # https://ssd-api.jpl.nasa.gov/doc/horizons.html
    # urlbase = str("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1")
    urlbase = str("https://ssd.jpl.nasa.gov/api/horizons.api?format=text")
    if os.path.exists("/tmp/tmp.csv"):
        os.remove("/tmp/tmp.csv")
    good = False
    while good is False:
        # url="https://ssd.jpl.nasa.gov/api/horizons.api?format=text"
        # url=${url}"&COMMAND='${target}%3B'&OBJ_DATA='YES'&MAKE_EPHEM='YES'"
        # url=${url}"&EPHEM_TYPE='OBSERVER'&CENTER='${center}'&TLIST='${TLIST}'"
        # url=${url}"&QUANTITIES='1'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'"
        if isasteroid:
            url = urlbase + str("&COMMAND='") + str(target) + str("%3B'")
            url = url + str("&CENTER='") + str(center) + str("'")
        if ismoon:
            url = urlbase + str("&COMMAND=") + str("'") + str(target) + \
                            str("'")
            url = url + str("&CENTER=coord@399&COORD_TYPE=GEODETIC")
            url += str("&SITE_COORD='157.500000%20 -11.950000 852.000000'")

        url = url + str("&OBJ_DATA='YES'&MAKE_EPHEM='YES'")
        url = url + str("&EPHEM_TYPE='OBSERVER'&TLIST='" + str(tlist) + "'")
        url = url + str("&QUANTITIES=''&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
        url = url + str("&CAL_FORMAT=BOTH")
        # url = url+str("&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'")
        # url = url+str("&TLIST='")+str(tlist)+str("'")
        # url = url+str("&QUANTITIES=''&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
        # url = url+str("&CAL_FORMAT=BOTH")
        # if ismoon:
        #     url = str(urlbase+str(&TABLE_TYPE=OBSERVER&QUANTITIES=''
        # &COMMAND=')+target+'&CSV_FORMAT=YES&CAL_FORMAT=BOTH&ANG_FORMAT=DEG
        # &APPARENT=AIRLESS&REF_SYSTEM=J2000&EXTRA_PREC=NO
        # &CENTER=coord@399&COORD_TYPE=GEODETIC&SITE_COORD=%27157.500000
        # %20-11.950000%20852.000000%27&TLIST=2449317.23542&SKIP_DAYLT=NO)
        # else:
        #     url = str(urlbase+"&COMMAND=")+str("'")+str(target)+str
        # (":'&CENTER='")+str(center)+str("'&MAKE_EPHEM='YES'&TABLE_TYPE=
        # 'OBSERVER'&TLIST='")+str(tlist)+str("'&QUANTITIES=''&CSV_FORMAT=
        # 'YES'&ANG_FORMAT='DEG'")
        # print(url)
        try:
            http = urllib3.PoolManager()
            r = http.request("GET", url)
            f = open("/tmp/tmp.csv", "wb")
            f.write(r.data)
            f.close()
        except Exception as ex:
            logging.exception(f"Caught an error at {time.strftime('%Y:%m:%d %H:%M:%S', time.localtime())}")
            time.sleep(20)
            good = False
        if os.path.exists("/tmp/tmp.csv"):
            statinfo = os.stat("/tmp/tmp.csv")
            if statinfo.st_size > 2000:
                good = True
            else:
                time.sleep(10)
                print(f"File size is too small! Trying again...")
        else:
            time.sleep(10)

    F = open("/tmp/tmp.csv", "r")
    for line in F:
        if "Date__" in line:
            header = line.rstrip()
            break
    F.close()

    i = 0
    F = open("/tmp/tmp.csv", "r")
    for line in F:
        while i == 1:
            data = line.rstrip()
            # print(line)
            i = i + 1
            break
        if "$$SOE" in line:
            i = i + 1
    F.close()

    header = header.replace(", , ,", ",empty,empty1,")
    header = header.replace("T-O-I/IB_Illu%", "T-O-I,IB_Illu%,")
    header = header[:-1]
    data = data[:-1]

    file = open("/tmp/tmp.dat", "w")
    file.write(header)
    file.write("\n")
    file.write(data)
    file.close()

    inp = Table.read("/tmp/tmp.dat", format="ascii.csv", header_start=0, data_start=1)
    date = inp["Date__(UT)__HR:MN:SC.fff"]
    # ra=inp['R.A._(ICRF/J2000.0)']
    # dec=inp['DEC_(ICRF/J2000.0)']
    ra = inp["R.A._(ICRF)"]
    dec = inp["DEC_(ICRF)"]
    RA_rate = inp["dRA*cosD"]
    DEC_rate = inp["d(DEC)/dt"]
    try:
        V = inp["APmag"]
    except Exception as e:
        V = inp["T-mag"]
    if V == 'n.a.':
        V = [-99.9]
    r = inp["r"]
    delta = inp["delta"]
    # lighttime=inp['1-way_LT']
    lighttime = inp["1-way_down_LT"]
    elong = inp["S-O-T"]
    elongFlag = inp["/r"]
    alpha = inp["S-T-O"]
    ObsEclLon = inp["ObsEcLon"]
    ObsEclLat = inp["ObsEcLat"]
    
    t = Table(
        [
            date,
            ra,
            dec,
            RA_rate,
            DEC_rate,
            V,
            r,
            delta,
            lighttime,
            elong,
            elongFlag,
            alpha,
            ObsEclLon,
            ObsEclLat,
        ],
        names=(
            "date",
            "ra",
            "dec",
            "RA_rate",
            "DEC_rate",
            "V",
            "r",
            "delta",
            "lighttime",
            "elong",
            "elongFlag",
            "alpha",
            "ObsEclLon",
            "ObsEclLat",
        ),
        meta={"name": "first table"},
    )
    return t


###############################################################################


def get_jplhelements(target, tlist, ismoon):
    # urlbase = str("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&")
    if os.path.exists("/tmp/tmp2.csv"):
        os.remove("/tmp/tmp2.csv")

    good = False
    while good is False:
        urlbase = str("https://ssd.jpl.nasa.gov/api/horizons.api?format=text")
        if ismoon:
            # url = str(urlbase+"&COMMAND='")+str(target)+str(":'&CENTER='500@3'&MAKE_EPHEM='YES'&TABLE_TYPE='elements'&TLIST='")+str(tlist)+str("'&QUANTITIES=''&OUT_UNITS='AU-D'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
            url = urlbase + str("&COMMAND='" + str(target) + "'")
            url += str("&CENTER='500@3'")
            url += str("&MAKE_EPHEM='YES'&TABLE_TYPE='elements'")
            url += str("&TLIST='" + str(tlist) + "'&QUANTITIES=''")
            url += str("&OUT_UNITS='AU-D'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
            # url = https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=%27301%27&CENTER=%27500@3%27&MAKE_EPHEM=%27YES%27&TABLE_TYPE=%27elements%27&TLIST=%272449317.23542%27&QUANTITIES=%27%27&OUT_UNITS=%27AU-D%27&CSV_FORMAT=%27YES%27&ANG_FORMAT=%27DEG%27
        else:
            # url = str(urlbase+"&COMMAND='")+str(target)+str(":'&CENTER='500@10'&MAKE_EPHEM='YES'&TABLE_TYPE='elements'&TLIST='")+str(tlist)+str("'&QUANTITIES=''&OUT_UNITS='AU-D'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
            url = urlbase + str("&COMMAND='" + str(target) + ":'")
            url += str("&CENTER='500@10'")
            url += str("&MAKE_EPHEM='YES'&TABLE_TYPE='elements'")
            url += str("&TLIST='" + str(tlist) + "'&QUANTITIES=''")
            url += str("&OUT_UNITS='AU-D'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'")
        # print(url)
        try:
            http = urllib3.PoolManager()
            r = http.request("GET", url)
            f = open("/tmp/tmp2.csv", "wb")
            f.write(r.data)
            f.close()
        except Exception as ex:
            logging.exception(f"Caught an error at {time.strftime('%Y:%m:%d %H:%M:%S', time.localtime())}")
            time.sleep(20)
            good = False
        if os.path.exists("/tmp/tmp2.csv"):
            statinfo = os.stat("/tmp/tmp2.csv")
            if statinfo.st_size > 2000:
                good = True
            else:
                time.sleep(10)
                print(f"File size is too small! Trying again...")
        else:
            time.sleep(10)

    F = open("/tmp/tmp2.csv", "r")
    for line in F:
        if "JDTDB" in line:
            header = line.rstrip()
            break
    F.close()

    i = 0
    F = open("/tmp/tmp2.csv", "r")
    for line in F:
        while i == 1:
            data = line.rstrip()
            i = i + 1
            break
        if "$$SOE" in line:
            i = i + 1
    F.close()
    file = open("/tmp/tmp2.dat", "w")
    file.write(header)
    file.write("\n")
    file.write(data)
    file.close()

    inp = Table.read("/tmp/tmp2.dat", format="ascii.csv")
    a = inp["A"]
    e = inp["EC"]
    incl = inp["IN"]
    Omega = inp["OM"]
    w = inp["W"]
    M = inp["MA"]
    t = Table(
        [a, e, incl, Omega, w, M],
        names=("a", "e", "incl", "Omega", "w", "M"),
        meta={"name": "first table"},
    )
    return t


###############################################################################


def getvecs1(target, tlist):
    if os.path.exists("/tmp/tmp0.csv"):
        os.remove("/tmp/tmp0.csv")
    good = False
    while good is False:
        # url=str("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&TABLE_TYPE='VECTORS'&CENTER='500@10'&COMMAND='")+str(target)+":'"+str("&TLIST='")+str(tlist)+str("'&CAL_FORMAT='JD'&VEC_TABLE='1'&OUT_UNITS='AU-D'&CSV_FORMAT='YES'")
        urlbase = str("https://ssd.jpl.nasa.gov/api/horizons.api?format=text")
        url = urlbase + str("&COMMAND='" + str(target) + ":'")
        url += str("&TLIST='" + str(tlist) + "'")
        url += str("&CAL_FORMAT='JD'&VEC_TABLE='1'&OUT_UNITS='AU-D'")
        url += str("&CSV_FORMAT='YES'&EPHEM_TYPE='VECTORS'")
        url += str("&CENTER='500@10'")
        try:
            http = urllib3.PoolManager()
            r = http.request("GET", url)
            f = open("/tmp/tmp0.csv", "wb")
            f.write(r.data)
            f.close()
        except Exception as ex:
            logging.exception(f"Caught an error at {time.strftime('%Y:%m:%d %H:%M:%S', time.localtime())}")
            time.sleep(20)
            good = False
        if os.path.exists("/tmp/tmp0.csv"):
            statinfo = os.stat("/tmp/tmp0.csv")
            if statinfo.st_size > 2000:
                good = True
            else:
                time.sleep(10)
                print(f"File size is too small! Trying again...")
        else:
            time.sleep(10)
    i = 0
    F = open("/tmp/tmp0.csv", "r")
    for line in F:
        while i == 1:
            data = line.rstrip().split()
            # print(line)
            i = i + 1
            break
        if "$$SOE" in line:
            i = i + 1
    F.close()
    return data[4].replace(",", ""), data[5].replace(",", ""), data[6].replace(",", "")


###############################################################################


def getvecs2(target, center, tlist):
    if os.path.exists("/tmp/tmp00.csv"):
        os.remove("/tmp/tmp00.csv")
    good = False
    while good is False:
        # url=str("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&TABLE_TYPE='VECTORS'&CENTER='")+str(center)+str("'&COMMAND='")+str(target)+":'"+str("&TLIST='")+str(tlist)+str("'&CAL_FORMAT='JD'&VEC_TABLE='1'&OUT_UNITS='AU-D'&CSV_FORMAT='YES'")
        # https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='33342:'&TLIST='2457338.48845'&CAL_FORMAT='JD'&VEC_TABLE='1'&OUT_UNITS='AU-D'&CSV_FORMAT='YES'&EPHEM_TYPE='VECTORS'
        urlbase = str("https://ssd.jpl.nasa.gov/api/horizons.api?format=text")
        url = urlbase + str("&COMMAND='" + str(target) + ":'")
        url += str("&TLIST='" + str(tlist) + "'")
        url += str("&CAL_FORMAT='JD'&VEC_TABLE='1'&OUT_UNITS='AU-D'")
        url += str("&CSV_FORMAT='YES'&EPHEM_TYPE='VECTORS'")
        url += str("&CENTER='" + str(center) + "'")
        try:
            http = urllib3.PoolManager()
            r = http.request("GET", url)
            f = open("/tmp/tmp00.csv", "wb")
            f.write(r.data)
            f.close()
        except Exception as ex:
            logging.exception(f"Caught an error at {time.strftime('%Y:%m:%d %H:%M:%S', time.localtime())}")
            time.sleep(20)
            good = False
        if os.path.exists("/tmp/tmp00.csv"):
            statinfo = os.stat("/tmp/tmp00.csv")
            if statinfo.st_size > 2000:
                good = True
            else:
                time.sleep(10)
                print(f"File size is too small! Trying again...")
        else:
            time.sleep(10)
    i = 0
    F = open("/tmp/tmp00.csv", "r")
    for line in F:
        while i == 1:
            data = line.rstrip().split()
            # print(line)
            i = i + 1
            break
        if "$$SOE" in line:
            i = i + 1
    F.close()
    return data[4].replace(",", ""), data[5].replace(",", ""), data[6].replace(",", "")


################################################################################
def cont_check(infile, outfile):
    tmpfile = "/tmp/input.tmp"
    if os.path.exists("/tmp/input.tmp"):
        os.remove("/tmp/input.tmp")

    num_lines_in = sum(1 for line in open(infile)) - 1

    if os.path.exists(outfile):
        num_lines_out = sum(1 for line in open(outfile)) - 2
    else:
        num_lines_out = 0
    diff = num_lines_in - num_lines_out - 1
    # print(num_lines_in,num_lines_out)

    if num_lines_in == num_lines_out:
        return tmpfile, diff, True
    else:
        file = open(infile, "r")
        content = file.readlines()
        file = open(tmpfile, "w")
        for l in content[num_lines_out:num_lines_in]:
            file.write(str(l))
        file.close()
        return tmpfile, diff, False


################################################################################


def getname(designations, id, did, dnames):
    tname = ""
    id = str(id)
    try:
        id = id.lstrip("0")
    except Exception as e:
        error = e
    try:
        id = int(id)
    except Exception as e:
        error = e
    # id(array(id))
    # print(type(id))
    for d in designations:
        try:
            # id=int(row[0])
            idx = did == id
            slen = len(str(dnames[idx]).split())
        except Exception as e:
            error = e
        if dnames[idx].size > 0:
            # print("Processing",str(dnames[idx]).replace('[','').replace(']','').replace("'",'').strip())
            if desig0s[idx]:
                tname = (
                    str(dnames[idx])
                    .replace("[", "")
                    .replace("]", "")
                    .strip("'")
                    .strip('"')
                    .strip()
                    .split()[1:]
                )
                altname = (
                    str(id)
                    + "#"
                    + str(desig0s[idx])
                    .replace("[", "")
                    .replace("]", "")
                    .strip("'")
                    .strip('"')
                    .strip()
                    + "#"
                    + str(naifids[idx]).replace("[", "").replace("]", "")
                )
            else:
                if desig1s[idx]:
                    tname = (
                        str(dnames[idx])
                        .replace("[", "")
                        .replace("]", "")
                        .strip("'")
                        .strip('"')
                        .strip()
                        .split()[1:]
                    )
                    altname = (
                        str(id)
                        + "#"
                        + str(desig1s[idx])
                        .replace("[", "")
                        .replace("]", "")
                        .strip("'")
                        .strip('"')
                        .strip()
                        + "#"
                        + str(naifids[idx]).replace("[", "").replace("]", "")
                    )
            break
        else:
            # id=str(row[0])
            idx = d == id
            if dnames[idx]:
                if desig1s[idx]:
                    tname = (
                        str(dnames[idx])
                        .replace("[", "")
                        .replace("]", "")
                        .strip("'")
                        .strip('"')
                        .strip()
                        .split()[1:]
                    )
                    altname = (
                        str(id)
                        + "#"
                        + str(desig1s[idx])
                        .replace("[", "")
                        .replace("]", "")
                        .strip("'")
                        .strip('"')
                        .strip()
                        + "#"
                        + str(naifids[idx]).replace("[", "").replace("]", "")
                    )
                else:
                    if desig2s[idx]:
                        tname = (
                            str(dnames[idx])
                            .replace("[", "")
                            .replace("]", "")
                            .strip("'")
                            .strip('"')
                            .strip()
                            .split()[1:]
                        )
                        altname = (
                            str(id)
                            + "#"
                            + str(desig2s[idx])
                            .replace("[", "")
                            .replace("]", "")
                            .strip("'")
                            .strip('"')
                            .strip()
                            + "#"
                            + str(naifids[idx]).replace("[", "").replace("]", "")
                        )
                # print("Processing",str(dnames[idx]).replace('[','').replace(']','').replace("'",'').strip())
                break
    namestring = []

    if tname:

        for i in range(0, len(tname)):
            if i == 0:
                namestring = str(tname[i])
            else:
                # print("tname[i]",tname[i])
                namestring = str(namestring) + " " + tname[i]
    else:
        namestring = (
            str(str(dnames[idx]).strip().split(" ")[1])
            .replace("[", "")
            .replace("]", "")
            .strip("'")
            .strip('"')
            .strip()
        )
        tname = namestring
        altname = str(id) + "#" + str(naifids[idx]).replace("[", "").replace("]", "")
    # print("id,idx",id,idx)
    # print(tname,id)
    naifid = str(naifids[idx]).replace("[", "").replace("]", "")
    return (
        str(tname).replace("[", "").replace("]", "").replace("'", "").replace(",", ""),
        altname,
        naifid,
        idx,
    )


################################################################################


def process(
    observatory_project,
    obscode,
    id,
    tname,
    jd,
    band,
    wvl,
    influx,
    einflux,
    quality_flags,
    comments_remarks,
    inp,
    i,
    documents_references,
    instrument_detector,
    idx,
    obsmode,
):
    influx = float(influx)
    einflux = float(einflux)
    # print(len(oinp['#naifid']))

    print("Horizons query")
    # print(id)
    dq = get_jplh(id, obscode, jd, ismoon)
    dq2 = get_jplhelements(id, jd, ismoon)
    # time.sleep(1)

    target_X = 0
    target_Y = 0
    target_Z = 0
    observer_X = 0
    observer_Y = 0
    observer_Z = 0
    jpl_obj_radius = 0
    jpl_obj_albedo = 0
    # print("H,G",dq['H'],dq['G'])
    print("XYZ vectors query")
    target_X, target_Y, target_Z = getvecs1(id, jd)
    observer_X, observer_Y, observer_Z = getvecs2(id, obscode, jd)

    x1 = float(target_X)
    y1 = float(target_Y)
    z1 = float(target_Z)
    x2 = float(observer_X)
    y2 = float(observer_Y)
    z2 = float(observer_Z)
    # Calculating the X-Y-Z coordinates of observer@sun
    x3 = "{:.15E}".format(Decimal(x1 - x2))
    y3 = "{:.15E}".format(Decimal(y1 - y2))
    z3 = "{:.15E}".format(Decimal(z1 - z2))

    # phase angle with sign
    r = np.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    D = np.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    cosa = (x1 * x2 + y1 * y2 + z1 * z2) / (r * D)
    zcom = x1 * y2 - y1 * x2
    if zcom > 0.0:
        # phase_angle=(atan2(sqrt(1-cosa*cosa),cosa)*180./3.1415926536)
        phase_angle_alpha = "%+5.2f" % (
            np.arctan2(np.sqrt(1 - cosa * cosa), cosa) * 180.0 / np.pi
        )
    if zcom < 0.0:
        phase_angle_alpha = "%+5.2f" % (
            -1.0 * np.arctan2(np.sqrt(1 - cosa * cosa), cosa) * 180.0 / np.pi
        )
    if zcom == 0.0:
        phase_angle_alpha = 0.0
    # print("zcom", zcom)
    # print("phase_angle_alpha", phase_angle_alpha)

    # Lighttime corrected epoch
    LTcorrected_epoch = (
        float(observation_mid_time) - (dq["lighttime"][0] * 60.0) / 3600.0 / 24.0
    )

    comet = []
    F = open("/tmp/tmp.csv", "r")
    for line in F:
        if "Comet" in line:
            comet = line.split()
            # print(line.split())
            break
    F.close()

    if comet:
        comments_remarks = (
            comments_remarks
            + " Object flagged as comet in JPL Horizons. M1 value was used instead of H for absolute magnitude. "
        )
        F = open("/tmp/tmp00.csv", "r")
        for line in F:
            if " M1=" in line:
                dq["M1"] = line.split()[1]
                # print(line.split()[1])
                break
        F.close()

        i = 0
        F = open("/tmp/tmp00.csv", "r")
        for line in F:
            while i == 1:
                jpl_obj_radius = line.split()[3]
                i = i + 1
                break
            if "RAD" in line:
                i = i + 1
                # break
        F.close()

        dq["H"] = dq["M1"]
        dq["G"] = ""

    else:
        F = open("/tmp/tmp.csv", "r")
        for line in F:
            if " H=" in line:
                dq["H"] = line.split()[1]
                break
        F.close()

        F = open("/tmp/tmp.csv", "r")
        for line in F:
            if " H=" in line:
                dq["G"] = line.split()[3]
                break
        F.close()

        F = open("/tmp/tmp2.csv", "r")
        for line in F:
            if "RAD" in line:
                jpl_obj_radius = line.split()[3]
                break
        F.close()
        F = open("/tmp/tmp2.csv", "r")
        for line in F:
            if "ALBEDO" in line:
                jpl_obj_albedo = line.split()[1]
                break
        F.close()

    altflag = 0
    if jpl_obj_albedo == "n.a.":
        # print("naifids[idx]",int(naifids[idx]))
        idx2 = naifids[idx] == anaifid
        if alb[idx2].size > 0:
            jpl_obj_albedo = (
                str(alb[idx2]).replace("[", "").replace("]", "").replace("'", "")
            )
            comments_remarks = (
                comments_remarks
                + " Albedo is taken from Ali-Lagoa & Delbo' (2017) and Ali-Lagoa et al. (2018. in press). "
            )
        elif id in albedos2["name"]:
            jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["name"] == id)][0]
            comments_remarks = (
                str(comments_remarks)
                + " Albedo is taken from "
                + str(albedos2["ref"][np.where(albedos2["name"] == id)][0])
            )
            print("Replacing albedo value!", jpl_obj_albedo)
            altflag = 1
        elif tname in albedos2["name"]:
            jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["name"] == tname)][0]
            comments_remarks = (
                str(comments_remarks)
                + " Albedo is taken from "
                + str(albedos2["ref"][np.where(albedos2["name"] == tname)][0])
            )
            print("Replacing albedo value!", jpl_obj_albedo)
            altflag = 1
        elif id in albedos2["number"]:
            jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["number"] == id)][0]
            comments_remarks = (
                str(comments_remarks)
                + " Albedo is taken from "
                + str(albedos2["ref"][np.where(albedos2["number"] == id)][0])
            )
            print("Replacing albedo value!", jpl_obj_albedo)
            altflag = 1
        else:

            jpl_obj_albedo = 0.10
            comments_remarks = (
                str(comments_remarks)
                + " An assumed geometric albedo of 0.1 was used to calculate the colour-correction factor. "
            )
            # print("jpl_obj_albedo:",jpl_obj_albedo)

    if float(jpl_obj_albedo) > 0.5:
        jpl_obj_albedo = 0.5
    if float(jpl_obj_albedo) < 0.03:
        jpl_obj_albedo = 0.03

    if naifid == "2162173":
        jpl_obj_radius = 0.45

    if jpl_obj_radius == "n.a.":
        # print("naifids[idx]",int(naifids[idx]))
        idx2 = naifids[idx] == anaifid
        if diam[idx2].size > 0:
            jpl_obj_radius = (
                float(
                    str(diam[idx2]).replace("[", "").replace("]", "").replace("'", "")
                )
                / 2.0
            )
            comments_remarks = (
                comments_remarks
                + " Radius is taken from Ali-Lagoa & Delbo' (2017) and Ali-Lagoa et al. (2018. in press). "
            )
        else:

            jpl_obj_radius = 0.0
            comments_remarks = comments_remarks + " No radius data for this object! "
    try:
        jplG = dq["G"][0]
    except Exception as e:
        jplG = ""

    if jplG == "n.a.":
        comments_remarks = (
            comments_remarks
            + " No G value for this object in JPL Horizons. An assumed slope parameter of 0.15 was used to calculate the colour-correction factor. "
        )
        jplG = 0.15

    if jplG == "":
        comments_remarks = (
            comments_remarks
            + " No G value for this object in JPL Horizons. An assumed slope parameter of 0.15 was used to calculate the colour-correction factor. "
        )
        jplG = 0.15

    try:
        if dq["H"][0] == "":
            comments_remarks = (
                comments_remarks + " No H value for this object in JPL Horizons. "
            )
    except Exception as e:
        comments_remarks = (
            comments_remarks + " No H value for this object in JPL Horizons. "
        )
        dq["H"] = ""

    try:
        V = dq["V"][0]
    except Exception as e:
        V = dq["Tmag"][0]

    # alt. albedo
    # print(id)
    if id in albedos2["name"] and altflag == 0:
        jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["name"] == id)][0]
        comments_remarks = (
            comments_remarks
            + " Albedo is taken from "
            + str(albedos2["ref"][np.where(albedos2["name"] == id)][0])
        )
        print("Replacing albedo value!", jpl_obj_albedo)
        altflag = 1
    if id in albedos2["number"] and altflag == 0:
        jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["number"] == id)][0]
        comments_remarks = (
            comments_remarks
            + " Albedo is taken from "
            + str(albedos2["ref"][np.where(albedos2["number"] == id)][0])
        )
        print("Replacing albedo value!", jpl_obj_albedo)
        altflag = 1
    if tname in albedos2["name"] and altflag == 0:
        jpl_obj_albedo = albedos2["albedo"][np.where(albedos2["name"] == tname)][0]
        comments_remarks = (
            str(comments_remarks)
            + " Albedo is taken from "
            + str(albedos2["ref"][np.where(albedos2["name"] == tname)][0])
        )
        print("Replacing albedo value!", jpl_obj_albedo)
        altflag = 1
    # print(comments_remarks)
    # calc. T_eff
    # T_eff = (392.0/np.sqrt(dq['r'][0]))*(1-float(jpl_obj_albedo)*(0.290 + 0.684*float(jplG)))**0.25
    print(dq["r"][0], jpl_obj_albedo, jplG)
    T_eff = (393.6 / np.sqrt(dq["r"][0])) * (
        1 - float(jpl_obj_albedo) * (0.290 + 0.684 * float(jplG))
    ) ** 0.25

    if observatory_project == "WISE":
        w1, w2, w3, w4 = ccorrect_wise(T_eff)
        # corrected flux and error
        if wvl == 3.35:
            abs_err_f = 0.05
            colour_correction_factor = float(("%5.3f" % float(w1)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            colour_corrected_flux_density = "%5.6f" % (float(influx) / float(w1))
            absolute_flux_error = "%5.6f" % (
                np.sqrt(
                    float(einflux) ** 2
                    + (abs_err_f * float(influx)) ** 2
                    + (cerr * float(influx)) ** 2
                )
                / float(w2)
            )
            # W2 data in general might contain reflected sunlight!
            # NotSaturated; QFr10; Warning! W1 data in general might contain reflected sunlight!
            comments_remarks = str(
                comments_remarks
                + ";"
                + warning
                + "; Warning! W1 data in general might contain reflected sunlight!"
            )

        if wvl == 4.60:
            abs_err_f = 0.05
            colour_correction_factor = float(("%5.3f" % float(w2)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            colour_corrected_flux_density = "%5.6f" % (float(influx) / float(w2))
            absolute_flux_error = "%5.6f" % (
                np.sqrt(
                    float(einflux) ** 2
                    + (abs_err_f * float(influx)) ** 2
                    + (cerr * float(influx)) ** 2
                )
                / float(w2)
            )
            comments_remarks = str(
                comments_remarks
                + ";"
                + warning
                + " ;Warning! W2 data in general might contain reflected sunlight!"
            )

        if wvl == 11.10:
            abs_err_f = 0.045
            colour_correction_factor = float(("%5.3f" % float(w3)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            colour_corrected_flux_density = "%5.6f" % (float(influx) / float(w3))
            absolute_flux_error = "%5.6f" % (
                np.sqrt(
                    float(einflux) ** 2
                    + (abs_err_f * float(influx)) ** 2
                    + (cerr * float(influx)) ** 2
                )
                / float(w3)
            )
            comments_remarks = str(comments_remarks + ";" + warning)

        elif wvl == 22.64:
            abs_err_f = 0.057
            colour_correction_factor = float(("%5.3f" % float(w4)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            colour_corrected_flux_density = "%5.6f" % (float(influx) / float(w4))
            absolute_flux_error = "%5.6f" % (
                np.sqrt(
                    float(einflux) ** 2
                    + (abs_err_f * float(influx)) ** 2
                    + (cerr * float(influx)) ** 2
                )
                / float(w4)
            )
            comments_remarks = str(comments_remarks + ";" + warning)

        documents_references = "WISEASD"
        obsmode = "survey"
        comments_remarks = str(comments_remarks)

    ##################################
    if observatory_project == "AKARI":
        # colour correction
        cc4, cc7, cc9, cc11, cc15, cc18, cc24 = ccorrect_akari(T_eff)
        # corrected flux and error
        abs_err_f = 0.035
        if wvl == 4:
            colour_correction_factor = float(("%5.3f" % float(cc9)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc4)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc4)
                )
            instrument_detector = "IRC-NIR"
            band = "N4"
            wvl = "4.0"

        if wvl == 7:
            colour_correction_factor = float(("%5.3f" % float(cc9)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03

            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc7)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc7)
                )
            instrument_detector = "IRC-MIR-S"
            band = "S7"
            wvl = "7.0"

        if wvl == 9:
            colour_correction_factor = float(("%5.3f" % float(cc9)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03

            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc9)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc9)
                )
            instrument_detector = "IRC-MIR-S"
            band = "S9W"
            wvl = "9.0"

        if wvl == 11:
            colour_correction_factor = float(("%5.3f" % float(cc9)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03

            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc11)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc11)
                )
            instrument_detector = "IRC-MIR-S"
            band = "S11"
            wvl = "11.0"

        elif wvl == 18:
            colour_correction_factor = float(("%5.3f" % float(cc18)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc18)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc18)
                )
            instrument_detector = "IRC-MIR-L"
            band = "L18W"
            wvl = "18.0"

        if wvl == 15:
            colour_correction_factor = float(("%5.3f" % float(cc15)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc15)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc15)
                )
            instrument_detector = "IRC-MIR-L"
            band = "L15"
            wvl = "15.0"

        elif wvl == 24:
            colour_correction_factor = float(("%5.3f" % float(cc24)))
            if colour_correction_factor >= 0.95 and colour_correction_factor <= 1.05:
                cerr = 0.01
            elif colour_correction_factor >= 0.9 and colour_correction_factor < 0.95:
                cerr = 0.02
            elif colour_correction_factor < 0.9 or colour_correction_factor > 1.05:
                cerr = 0.03
            if tname == "Itokawa":
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc24)
                )
            else:
                colour_corrected_flux_density = "%5.6f" % (float(influx) / float(cc9))
                absolute_flux_error = "%5.6f" % (
                    np.sqrt(
                        float(einflux) ** 2
                        + (abs_err_f * float(influx)) ** 2
                        + (cerr * float(influx)) ** 2
                    )
                    / float(cc24)
                )
            instrument_detector = "IRC-MIR-L"
            band = "L24"
            wvl = "24.0"

        if documents_references == "U11":
            documents_references = "Usui et al. 2011"

        if documents_references == "H13":
            documents_references = "Hasegawa et al. 2013"

        if documents_references == "M16":
            documents_references = "Mueller T. G. et al. 2017"

        if documents_references == "M14":
            documents_references = "Mueller T. G. et al. 2014"

        documents_references = str(documents_references + "#" + "AKARIAFC")

    ##################################
    if observatory_project == "HSO":
        abs_err = 0.05
        cc_err = 0.01
        rel_err = 0.01
        if od_num == "1441" and id == 1:
            rel_err = 0.0

        oids = [
            1342202078,
            1342221723,
            1342205032,
            1342202076,
            1342202077,
            1342202079,
            1342186132,
            1342186133,
            1342221724,
            1342221725,
            1342205033,
            1342203465,
            1342203466,
            1342232383,
            1342232384,
        ]

        if obsmode == "Scanmap":
            if band == "blue":
                noise = 0.004
                colour_correction_factor = 1.005
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.015
            if band == "green":
                noise = 0.005
                colour_correction_factor = 1.025
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.02
            if band == "red":
                noise = 0.013
                colour_correction_factor = 1.06
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.04

        if obsmode == "ChopNod":
            if band == "blue":
                noise = 0.021
                colour_correction_factor = 1.005
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.015
            if band == "green":
                noise = 0.020
                colour_correction_factor = 1.025
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.02
            if band == "red":
                noise = 0.052
                colour_correction_factor = 1.06
                einflux = np.sqrt((0.01 * influx) ** 2 + noise ** 2)
                if obs_id in oids:
                    rel_err = 0.04

        if od_num == "1441" and id == 1:
            einflux = 0.01 * influx

        rel_error = np.sqrt(noise ** 2 + (rel_err * influx) ** 2)
        flux_cc = influx / colour_correction_factor
        abs_error = np.sqrt(
            (rel_error / colour_correction_factor) ** 2
            + (abs_err * flux_cc) ** 2
            + (cc_err * flux_cc) ** 2
        )
        colour_corrected_flux_density = "%5.6f" % (flux_cc)
        absolute_flux_error = "%5.6f" % (abs_error)

        if documents_references == "1":
            documents_references = "M10"
        if documents_references == "2":
            documents_references = "L10"
        if documents_references == "3":
            documents_references = "LIM10"
        if documents_references == "4":
            documents_references = "SS12"
        if documents_references == "5":
            documents_references = "MM12"
        if documents_references == "6":
            documents_references = "V12"
        if documents_references == "7":
            documents_references = "P12"
        if documents_references == "8":
            documents_references = "F13"
        if documents_references == "9":
            documents_references = "L13"
        if documents_references == "10":
            documents_references = "V14"
        if documents_references == "11":
            documents_references = "D14"
        if documents_references == "12":
            documents_references = "SS17"
        if documents_references == "13":
            documents_references = "M19A"
        if documents_references == "14":
            documents_references = "V18"
        if documents_references == "15":
            documents_references = "FT20"
        if documents_references == "16":
            documents_references = "K13"
        if documents_references == "17":
            documents_references = "L16"

        if documents_references == "18":
            documents_references = "P15"
        if documents_references == "19":
            documents_references = "P16"
        if documents_references == "20":
            documents_references = "M17"
        if documents_references == "21":
            documents_references = "M12"
        if documents_references == "22":
            documents_references = "M13"
        if documents_references == "23":
            documents_references = "M14"
        if documents_references == "M19B":
            documents_references = "M19B"
        # if documents_references == '17':
        #    documents_references="L16"
        # if documents_references == '17':
        #    documents_references="L16"

        if "color-corrected" in comments_remarks:
            colour_corrected_flux_density = "%5.4f" % (influx)
            absolute_flux_error = "%5.4f" % (einflux)
            # influx=0.0
            # einflux=0.0

        # documents_references="M18"

    ##################################
    if observatory_project == "MSX":
        # colour correction
        ccA, ccC, ccD, ccE = ccorrect_msx(T_eff)

        # corrected flux and error
        if wvl == 8.28:
            colour_correction_factor = "%5.3f" % float(ccA)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(ccA))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.05 * float(influx)) ** 2) / float(ccA)
            )

        elif wvl == 12.13:
            colour_correction_factor = "%5.3f" % float(ccC)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(ccC))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.05 * float(influx)) ** 2) / float(ccC)
            )

        elif wvl == 14.65:
            colour_correction_factor = "%5.3f" % float(ccD)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(ccD))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.06 * float(influx)) ** 2) / float(ccD)
            )

        elif wvl == 21.34:
            colour_correction_factor = "%5.3f" % float(ccE)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(ccE))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.06 * float(influx)) ** 2) / float(ccE)
            )
        documents_references = "T02MSX"
        obsmode = "survey"
    ##################################

    if observatory_project == "IRAS":
        # colour correction
        cc12, cc25, cc60, cc100 = ccorrect_simps(T_eff)

        # corrected flux and error
        if wvl == 12.0:
            colour_correction_factor = "%5.3f" % float(cc12)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(cc12))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.1 * float(influx)) ** 2) / float(cc12)
            )

        elif wvl == 25.0:
            colour_correction_factor = "%5.3f" % float(cc25)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(cc25))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.1 * float(influx)) ** 2) / float(cc25)
            )

        elif wvl == 60.0:
            colour_correction_factor = "%5.3f" % float(cc60)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(cc60))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.15 * float(influx)) ** 2) / float(cc60)
            )

        elif wvl == 100.0:
            colour_correction_factor = "%5.3f" % float(cc100)
            colour_corrected_flux_density = "%5.3f" % (float(influx) / float(cc100))
            absolute_flux_error = "%5.3f" % (
                np.sqrt(float(einflux) ** 2 + (0.15 * float(influx)) ** 2)
                / float(cc100)
            )
        documents_references = "T02IRAS"
        obsmode = "survey"

    if observatory_project == "HSO":
        observation_IDs = obs_id
    else:
        observation_IDs = ""
        # str(id)+str(observatory_project)
    # granule_uid=str(observation_IDs)+"_proc"

    # print("wvl,colour_correction_factor,colour_corrected_flux_density:",wvl,colour_correction_factor,colour_corrected_flux_density)
    print("Writing results to file", ofile)
    # print("-------------------------------\n")

    f = open((ofile), "a")
    towrite = (
        str(naifid).replace("[", "").replace("]", "")
        + ","
        + str(tname)
        + ","
        + str(observatory_project)
        + ","
        + str(obscode)
        + ","
        + str(instrument_detector)
        + ","
        + str(obsmode)
        + ","
        + str(observation_IDs)
        + ","
        + str("%5.5f" % float(observation_start_time))
        + ","
        + str("%5.5f" % float(observation_mid_time))
        + ","
        + str("%5.5f" % float(observation_end_time))
        + ","
        + str(dq["date"][0]).replace("[", "").replace("]", "").replace("'", "")
        + ","
        + str(band)
        + ","
        + str(influx)
        + ","
        + str(einflux)
        + ","
        + str(quality_flags)
        + ","
        + str(dq2["a"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq2["e"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq2["incl"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq2["Omega"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq2["w"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq2["M"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["H"][0]).replace("[", "").replace("]", "")
        + ","
        + str(str(jplG).replace("[", "").replace("]", ""))
        + ","
        + str(str(jpl_obj_radius)).replace("[", "").replace("]", "")
        + ","
        + str(jpl_obj_albedo)
        + ","
        + str(dq["ra"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["dec"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["RA_rate"][0] / 3600.0).replace("[", "").replace("]", "")
        + ","
        + str(dq["DEC_rate"][0] / 3600.0).replace("[", "").replace("]", "")
        + ","
        + str(dq["V"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["r"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["delta"][0]).replace("[", "").replace("]", "")
        + ", "
        + str(dq["lighttime"][0] * 60.0).replace("[", "").replace("]", "")
        + ","
        + str(dq["elong"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["elongFlag"][0]).replace("[", "").replace("]", "")
        + ","
        + str(phase_angle_alpha)
        + ","
        + str(dq["ObsEclLon"][0]).replace("[", "").replace("]", "")
        + ","
        + str(dq["ObsEclLat"][0]).replace("[", "").replace("]", "")
        + ","
        + str(target_X)
        + ","
        + str(target_Y)
        + ","
        + str(target_Z)
        + ","
        + str(observer_X)
        + ","
        + str(observer_Y)
        + ","
        + str(observer_Z)
        + ","
        + str(x3)
        + ","
        + str(y3)
        + ","
        + str(z3)
        + ","
        + str(wvl)
        + ","
        + str(colour_correction_factor)
        + ","
        + str(colour_corrected_flux_density)
        + ","
        + str(absolute_flux_error)
        + ","
        + str(comments_remarks)
        + ", "
        + (str(LTcorrected_epoch).replace("[", "").replace("]", ""))[0:14]
        + ", "
        + str(documents_references)
        + ","
        + str(ofname)
        + ","
        + str(nowdate)
        + ","
        + str(altname)
    )
    f.write(towrite)
    f.write("\n")
    f.close
    # print(documents_references,ofname)


################################################################################
################################################################################


def getwisedata(inp, i):
    jd = inp["JD"][i]
    observation_start_time = str(jd).strip()
    observation_mid_time = str(jd).strip()
    observation_end_time = str(jd).strip()
    band = inp["Band"][i]
    wvl = inp["L(um)"][i]
    influx = "%5.6f" % (inp["F(Jy)"][i])
    einflux = "%5.6f" % (inp["eF(Jy)"][i])
    quality_flags = inp["AQ"][i]
    comments_remarks = str(inp["Sat_flag"][i])
    try:
        warning = str(inp["Warning"][i])
    except Exception as e:
        print(e)
        warning = str("")

    return (
        jd,
        observation_start_time,
        observation_mid_time,
        observation_end_time,
        band,
        wvl,
        influx,
        einflux,
        quality_flags,
        comments_remarks,
        warning,
    )


################################################################################
def getakaridata(inp, i):
    jd = inp["JDmid-time"][i]
    observation_start_time = str(jd).strip()
    observation_mid_time = str(jd).strip()
    observation_end_time = str(jd).strip()
    if tname == "Itokawa":
        influx = "%5.6f" % (inp["in-bandflux"][i])
        einflux = "%5.6f" % (inp["Ein-bandflux"][i])
    else:
        influx = "%5.6f" % (inp["in-bandflux"][i])
        einflux = "%5.6f" % (inp["Ein-bandflux"][i])
    wvl = inp["wavelength"][i]
    if wvl == 4:
        band = "N4"
    if wvl == 7:
        band = "S7"
    if wvl == 9:
        band = "S9W"
    if wvl == 11:
        band = "S11"
    if wvl == 18:
        band = "L18W"
    if wvl == 15:
        band = "L15"
    if wvl == 24:
        band = "L24"
    obsmode = inp["Observationmode"][i]
    documents_references = str(inp["Reference"][i]).replace(" ", "")
    return (
        jd,
        observation_start_time,
        observation_mid_time,
        observation_end_time,
        band,
        wvl,
        influx,
        einflux,
        obsmode,
        documents_references,
    )


################################################################################
def gethsodata(inp, i):
    obs_id = inp["obs_id"][i]
    od_num = inp["od_num"][i]
    instrument = inp["instrument"][i]
    obs_mode = inp["obs_mode"][i]
    observation_start_time = str(inp["jd_start"][i]).strip()
    observation_end_time = str(inp["jd_end"][i]).strip()
    observation_mid_time = "{0:.5f}".format(
        ((float(observation_start_time)) + float((observation_end_time))) / 2
    )
    band = inp["band"][i]
    if obs_mode == "PS_Photom":
        obs_mode = "ChopNod"
    jd = observation_mid_time
    if band == "blue":
        wvl = 70
    if band == "green":
        wvl = 100
    if band == "red":
        wvl = 160
    influx = "%5.4f" % (inp["flux_eef_corr"][i])
    einflux = "%5.4f" % (inp["unc_flux_eef_corr"][i])
    """
    com=str(inp['comment'][i])
    m=(re.search(r"\(([0-9])+\)",com))
    documents_references=(m.group(0).replace("(","")).replace(")","")
    mm=m.group(0)
    print("documents_references",documents_references)
    instrument_detector=instrument
    obsmode=obs_mode
    comments_remarks=com.replace(mm,"")
    """
    documents_references = (str(inp["ref"][i])).strip()
    obsmode = obs_mode
    instrument_detector = instrument
    comments_remarks = inp["comment"][i]
    return (
        jd,
        observation_start_time,
        observation_mid_time,
        observation_end_time,
        band,
        wvl,
        influx,
        einflux,
        obsmode,
        obs_id,
        od_num,
        instrument_detector,
        documents_references,
        comments_remarks,
    )


################################################################################
def getmsxdata(inp, i):
    jd = "{0:.5f}".format(float(inp["JD"][i]))
    observation_start_time = jd
    observation_end_time = jd
    observation_mid_time = jd
    band = inp["Band"][i]
    wvl = inp["Wvl"][i]
    influx = "%5.3f" % (inp["f(Jy)"][i])
    einflux = "%5.3f" % (inp["ef(Jy)"][i])
    return (
        jd,
        observation_start_time,
        observation_mid_time,
        observation_end_time,
        band,
        wvl,
        influx,
        einflux,
    )


################################################################################
def getirasdata(inp, i):
    jd = "{0:.5f}".format(float(inp["JD"][i]))
    observation_start_time = jd
    observation_end_time = jd
    observation_mid_time = jd
    band = inp["Band"][i]
    wvl = inp["Wvl"][i]
    influx = "%5.3f" % (inp["f(Jy)"][i])
    einflux = "%5.3f" % (inp["ef(Jy)"][i])
    return (
        jd,
        observation_start_time,
        observation_mid_time,
        observation_end_time,
        band,
        wvl,
        influx,
        einflux,
    )


################################################################################
fullCmdArguments = sys.argv
argumentList = fullCmdArguments[1:]
unixOptions = "h:wd:i:"
gnuOptions = ["help", "workdir=", "inst="]

try:
    arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

for currentArgument, currentValue in arguments:
    if currentArgument in ("-h", "--help"):
        print("Usage:")
        print(
            """python3ginop_main.py --help --workdir=wd --inst=inst
              """
        )
        sys.exit(2)
    elif currentArgument in ("-wd", "--workdir"):
        workdir = currentValue
    elif currentArgument in ("-i", "--inst"):
        inst = currentValue

# entering working directory
os.chdir(workdir)
scriptdir = "/data/munka/git/SBNAF-IRDB/"
namedatabase = "DASTCOM.IDX"  # datafile for names and designations
albedofile = "albedos.csv"  # datafile for missing albedos
prog = True  # show progressbar
ismoon = False
isasteroid = True

if os.path.exists(namedatabase) is False:
    print("Downloading DASTCOM.IDX.")
    url = "ftp://ssd.jpl.nasa.gov/pub/xfr/DASTCOM.IDX"
    wget.download(url)

print("\nReading DASTCOM.IDX")
with open(namedatabase) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=",")
    did = []
    dnames = []
    naifids = []
    desig0s = []
    desig1s = []
    desig2s = []
    eee = []
    for row in readCSV:
        dname = row[0].strip()
        id = int((dname.strip()).split(" ")[0])
        try:
            naifid = int(row[1])
        except Exception as e:
            error = e
            naifid = 0
        desig0 = row[2]
        desig1 = row[3]
        try:
            desig2 = row[4]
        except Exception as e:
            error = e
            desig2 = ""
        dnames.append(dname)
        naifids.append(naifid)
        desig0s.append(desig0)
        desig1s.append(desig1)
        desig2s.append(desig2)
        did.append(id)

did = array(did)
naifids = array(naifids)
desig0s = array(desig0s)
desig1s = array(desig1s)
desig2s = array(desig2s)
dnames = array(dnames)
designations = [desig0s, desig1s, desig2s]
################################################################################
print("Reading albedos.csv")
with open(albedofile) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=",")
    anaifid = []
    astnum = []
    H = []
    G = []
    diam = []
    alb = []
    for row in readCSV:
        nid = int(row[0])
        anum = row[1]
        aH = row[2]
        aG = row[3]
        adiam = row[4]
        aalb = row[5]
        anaifid.append(nid)
        astnum.append(anum)
        H.append(aH)
        G.append(aG)
        diam.append(adiam)
        alb.append(aalb)

anaifid = array(anaifid)
astnum = array(astnum)
H = array(H)
G = array(G)
diam = array(diam)
alb = array(alb)
################################################################################
print("Reading albedos2.csv")
albedos2 = Table.read("obj_pv_hv_r_v3.csv", format="ascii.csv", data_start=1)


################################################################################
observatory_projects = [inst]
print(observatory_projects)
# observatory_projects=['AKARI','HSO','MSX','IRAS','WISE']

for observatory_project in observatory_projects:
    print(observatory_project)
    if observatory_project == "AKARI":
        fname = "AKARI_allsky_wheader.csv"  # input filename
        # fname='akari_test.csv'
        ofname = fname
        observatory_project = "AKARI"
        obscode = "500@399"
        ofile = "AstIrDbTbl_akari.csv"
        # ofile='akaritest_out.csv'

    if observatory_project == "HSO":
        # fname='hso_pacs_mba_phot_all_sorted_wheader.csv'
        # fname='hso_pacs_tno_wheader.csv'
        fname = "hso_wheader.csv"
        ofname = fname
        observatory_project = "HSO"
        obscode = "500@-486"
        # ofile='AstIrDbTbl_hso_test.csv' # output filename
        ofile = "hso_out.csv"

    if observatory_project == "MSX":
        fname = "MIMPS_MSX_datafile2_forPostgerSQL_singleflux_nozeros_wheader.csv"  # input filename
        # fname='msxtest.csv'
        ofname = fname
        observatory_project = "MSX"
        obscode = "500@399"
        ofile = "AstIrDbTbl_mimps.csv"  # output filename
        # ofile='msxtest_out.csv'

    if observatory_project == "IRAS":
        fname = "SIMPS.FP208B_forPostgerSQL_singleflux_filtered_wheader.csv"  # input filename
        # fname='irastest.csv'
        ofname = fname
        observatory_project = "IRAS"
        obscode = "500@399"
        ofile = "AstIrDbTbl_simps.csv"  # output filename
        # ofile='irastest_out.csv'

    if observatory_project == "WISE":
        fname='wise_akari_wheader.csv' # input filename
        # fname = '03200_SBNAF.csv'
        # fname = "33342_MovObjSearch.csv"
        ofname = fname
        observatory_project = "WISE"
        obscode = "500@-163"
        ofile='AstIrDbTbl_wise_akari.csv' # output filename
        # ofile = "33342_MovObjSearch_out.csv"

    ################################################################################

    finname, diff, cond = cont_check(fname, ofile)

    if os.path.exists(ofile):
        exists = 1
    else:
        exists = 0

    if exists == 0:
        # writing csv header
        f = open((ofile), "w")
        towrite = "#naifid,targetname,observatory_project,observatory_code,instrument_detector,obsmode,observation_IDs,observation_start_time,observation_mid_time,observation_end_time,datetime,band_filter,calibrated_inband_flux_Jy,inband_flux_error_Jy,quality_flags,orbital_param_A,orbital_param_EC,orbital_param_IN,orbital_param_OM,orbital_param_W,orbital_param_MA,absolute_magnitude_H,slope_parameter_G,jpl_obj_radius,jpl_obj_albedo,Right_Ascension_RA,Declination_DEC,RA_rate,DEC_rate,apparent_magnitude_V,heliocentric_distance_r,obscentric_distance_delta,lighttime,solar_elongation_elong,before_after_opposition,phase_angle_alpha,ObsEclLon,ObsEclLat,target_X_sun,target_Y_sun,target_Z_sun,target_X_observer,target_Y_observer,target_Z_observer,observer_X_sun,observer_Y_sun,observer_Z_sun,reference_wavelengths_micron,colour_correction_factor,colour_corrected_flux_density,absolute_flux_error,comments_remarks,LTcorrected_epoch,documents_references,input_table_source,data_last_modification,alt_target_name"
        f.write(towrite)
        f.write("\n")
        towrite = "#LONG,STRING,STRING,STRING,STRING,STRING,STRING,DOUBLE,DOUBLE,DOUBLE,STRING,STRING,DOUBLE,DOUBLE,STRING,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,DOUBLE,DOUBLE,FLOAT,FLOAT,STRING,FLOAT,FLOAT,FLOAT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,FLOAT,FLOAT,DOUBLE,DOUBLE,STRING,DOUBLE,STRING,STRING,STRING,STRING"
        f.write(towrite)
        f.write("\n")
        f.close()

    if sum(1 for line in open(ofile)) < 2:
        # writing csv header
        f = open((ofile), "w")
        towrite = "#naifid,targetname,observatory_project,observatory_code,instrument_detector,obsmode,observation_IDs,observation_start_time,observation_mid_time,observation_end_time,datetime,band_filter,calibrated_inband_flux_Jy,inband_flux_error_Jy,quality_flags,orbital_param_A,orbital_param_EC,orbital_param_IN,orbital_param_OM,orbital_param_W,orbital_param_MA,absolute_magnitude_H,slope_parameter_G,jpl_obj_radius,jpl_obj_albedo,Right_Ascension_RA,Declination_DEC,RA_rate,DEC_rate,apparent_magnitude_V,heliocentric_distance_r,obscentric_distance_delta,lighttime,solar_elongation_elong,before_after_opposition,phase_angle_alpha,ObsEclLon,ObsEclLat,target_X_sun,target_Y_sun,target_Z_sun,target_X_observer,target_Y_observer,target_Z_observer,observer_X_sun,observer_Y_sun,observer_Z_sun,reference_wavelengths_micron,colour_correction_factor,colour_corrected_flux_density,absolute_flux_error,comments_remarks,LTcorrected_epoch,documents_references,input_table_source,data_last_modification,alt_target_name"
        f.write(towrite)
        f.write("\n")
        towrite = "#LONG,STRING,STRING,STRING,STRING,STRING,STRING,DOUBLE,DOUBLE,DOUBLE,STRING,STRING,DOUBLE,DOUBLE,STRING,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,DOUBLE,DOUBLE,FLOAT,FLOAT,STRING,FLOAT,FLOAT,FLOAT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,FLOAT,FLOAT,DOUBLE,DOUBLE,STRING,DOUBLE,STRING,STRING,STRING,STRING"
        f.write(towrite)
        f.write("\n")
        f.close()

    now = datetime.datetime.now()
    nowdate = now.strftime("%Y-%m-%d %H:%M:%S")
    print("Reading input file", fname)
    inp = Table.read(fname, format="ascii.csv", data_start=1)

    print("Reading output file", ofile)
    oinp = Table.read(ofile, format="ascii.csv", data_start=2)
    ################################################################################

    if len(oinp["#naifid"]) == 0:
        for i in range(0, len(inp["Desig"])):
            id = " "
            jd = " "
            instrument_detector = " "
            observation_IDs = " "
            quality_flags = " "
            target_X = " "
            target_Y = " "
            target_Z = " "
            observer_X = " "
            observer_Y = " "
            observer_Z = " "
            colour_correction_factor = " "
            colour_corrected_flux_density = " "
            absolute_flux_error = " "
            LTcorrected_epoch = " "
            documents_references = " "
            comments_remarks = " "
            band = " "
            altname = " "
            # granule_uid=' '
            obsmode = " "
            id = inp["Desig"][i]
            # print("id",id)
            tname, altname, naifid, idx = getname(designations, id, did, dnames)

            if observatory_project == "WISE":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    quality_flags,
                    comments_remarks,
                    warning,
                ) = getwisedata(inp, i)

            if observatory_project == "AKARI":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    documents_references,
                ) = getakaridata(inp, i)

            if observatory_project == "HSO":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    obs_id,
                    od_num,
                    instrument_detector,
                    documents_references,
                    comments_remarks,
                ) = gethsodata(inp, i)

            if observatory_project == "MSX":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getmsxdata(inp, i)

            if observatory_project == "IRAS":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getirasdata(inp, i)
            if prog is True:
                id = str(id)
                try:
                    id = id.lstrip("0")
                except Exception as e:
                    error = e
                stat = "Processing " + str(id) + " " + str(tname)
                progress(i + 1, len(inp["Desig"]), status=stat)
                print(" ")
            else:
                print("Processing", id, tname)
            process(
                observatory_project,
                obscode,
                id,
                tname,
                jd,
                band,
                wvl,
                influx,
                einflux,
                quality_flags,
                comments_remarks,
                inp,
                i,
                documents_references,
                instrument_detector,
                idx,
                obsmode,
            )

    ################################################################################################################################################################################################################################################

    else:
        # ts1 = time.time()
        print("Continuing.")
        toprocess = []
        toremove = []
        total = len(oinp["#naifid"]) * len(inp["Desig"])
        ii = 0
        jj = 0
        naifs = []
        omts = []
        omts = array(omts)
        bands = []
        for i in range(0, len(inp["Desig"])):
            ii = ii + 1
            id = inp["Desig"][i]
            tname, altname, naifid, idx = getname(designations, id, did, dnames)
            if observatory_project == "WISE":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    quality_flags,
                    comments_remarks,
                    warning,
                ) = getwisedata(inp, i)

            if observatory_project == "AKARI":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    documents_references,
                ) = getakaridata(inp, i)

            if observatory_project == "HSO":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    obs_id,
                    od_num,
                    instrument_detector,
                    documents_references,
                    comments_remarks,
                ) = gethsodata(inp, i)

            if observatory_project == "MSX":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getmsxdata(inp, i)

            if observatory_project == "IRAS":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getirasdata(inp, i)
            naifs.append(naifid)
            omts = np.append(omts, "%5.5f" % float(str(observation_mid_time).strip()))
            bands.append(band)
        omts.astype(float)
        inp["#naifids"] = naifs
        try:
            inp["#naifids"] = inp["#naifids"].astype(str)
        except Exception as e:
            error = e
        newcol = Column(name="observation_mid_time", data=omts, dtype=float)
        inp.add_column(newcol)
        # inp.info()
        inp["observation_mid_time"] = inp["observation_mid_time"].astype(float)
        inp["band_filter"] = bands
        oinp["#naifids"] = oinp["#naifid"].astype(str)
        sdiff = setdiff(
            inp, oinp, keys=["#naifids", "observation_mid_time", "band_filter"]
        )
        inp = sdiff

        if len(inp["Desig"]) == 0:
            print("Input file already processed!")

        for i in range(0, len(inp["Desig"])):
            id = " "
            jd = " "
            instrument_detector = " "
            observation_IDs = " "
            quality_flags = " "
            target_X = " "
            target_Y = " "
            target_Z = " "
            observer_X = " "
            observer_Y = " "
            observer_Z = " "
            colour_correction_factor = " "
            colour_corrected_flux_density = " "
            absolute_flux_error = " "
            LTcorrected_epoch = " "
            documents_references = " "
            comments_remarks = " "
            band = " "
            altname = " "
            # granule_uid=' '
            obsmode = " "

            id = inp["Desig"][i]
            # print(id)
            tname, altname, naifid, idx = getname(designations, id, did,
                                                  dnames)
            # print(tname, altname, naifid)
            if observatory_project == "WISE":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    quality_flags,
                    comments_remarks,
                    warning,
                ) = getwisedata(inp, i)

            if observatory_project == "AKARI":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    documents_references,
                ) = getakaridata(inp, i)

            if observatory_project == "HSO":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                    obsmode,
                    obs_id,
                    od_num,
                    instrument_detector,
                    documents_references,
                    comments_remarks,
                ) = gethsodata(inp, i)

            if observatory_project == "MSX":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getmsxdata(inp, i)

            if observatory_project == "IRAS":
                (
                    jd,
                    observation_start_time,
                    observation_mid_time,
                    observation_end_time,
                    band,
                    wvl,
                    influx,
                    einflux,
                ) = getirasdata(inp, i)

            if prog is True:
                id = str(id)
                try:
                    id = id.lstrip("0")
                except Exception as e:
                    error = e
                stat = "Processing " + str(id) + " " + str(tname)
                progress(len(oinp["#naifid"]) + i + 1, ii, status=stat)
                print(" ")
            else:
                print("Processing", id, tname)
            process(
                observatory_project,
                obscode,
                id,
                tname,
                jd,
                band,
                wvl,
                influx,
                einflux,
                quality_flags,
                comments_remarks,
                inp,
                i,
                documents_references,
                instrument_detector,
                idx,
                obsmode,
            )

print("All done!")
