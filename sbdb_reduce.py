#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 08:31:09 2021

@author: iggymatheson
"""
#%%
def sbdb_des_to_packed_mpc_des(des1): # bring in SBDB designations, get out MPC designations
    des1 = str(des1)
    if len(des1) == 5:
        des2 = des1 # leave it unchanged, example is 15760 for Albion
    if len(des1) == 6: # change the first two numerical digits to a letter
        first_two_digits = des1[0:2]
        last_four_digits = des1[2:6]
        if first_two_digits == '10':
            one_digit = 'A'
        elif first_two_digits == '11':
            one_digit = 'B'
        elif first_two_digits == '12':
            one_digit = 'C'
        elif first_two_digits == '13':
            one_digit = 'D'
        elif first_two_digits == '14':
            one_digit = 'E'
        elif first_two_digits == '15':
            one_digit = 'F'
        elif first_two_digits == '16':
            one_digit = 'G'
        elif first_two_digits == '17':
            one_digit = 'H'
        elif first_two_digits == '18':
            one_digit = 'I'
        elif first_two_digits == '19':
            one_digit = 'J'
        elif first_two_digits == '20':
            one_digit = 'K'
        elif first_two_digits == '21':
            one_digit = 'L'
        elif first_two_digits == '22':
            one_digit = 'M'
        elif first_two_digits == '23':
            one_digit = 'N'
        elif first_two_digits == '24':
            one_digit = 'O'
        elif first_two_digits == '25':
            one_digit = 'P'
        elif first_two_digits == '26':
            one_digit = 'Q'
        elif first_two_digits == '27':
            one_digit = 'R'
        elif first_two_digits == '28':
            one_digit = 'S'
        elif first_two_digits == '29':
            one_digit = 'T'
        elif first_two_digits == '30':
            one_digit = 'U'
        elif first_two_digits == '31':
            one_digit = 'V'
        elif first_two_digits == '32':
            one_digit = 'W'
        elif first_two_digits == '33':
            one_digit = 'X'
        elif first_two_digits == '34':
            one_digit = 'Y'
        elif first_two_digits == '35':
            one_digit = 'Z'
        elif first_two_digits == '36':
            one_digit = 'a'
        elif first_two_digits == '37':
            one_digit = 'b'
        elif first_two_digits == '38':
            one_digit = 'c'
        elif first_two_digits == '39':
            one_digit = 'd'
        elif first_two_digits == '40':
            one_digit = 'e'
        elif first_two_digits == '41':
            one_digit = 'f'
        elif first_two_digits == '42':
            one_digit = 'g'
        elif first_two_digits == '43':
            one_digit = 'h'
        elif first_two_digits == '44':
            one_digit = 'i'
        elif first_two_digits == '45':
            one_digit = 'j'
        elif first_two_digits == '46':
            one_digit = 'k'
        elif first_two_digits == '47':
            one_digit = 'l'
        elif first_two_digits == '48':
            one_digit = 'm'
        elif first_two_digits == '49':
            one_digit = 'n'
        elif first_two_digits == '50':
            one_digit = 'o'
        elif first_two_digits == '51':
            one_digit = 'p'
        elif first_two_digits == '52':
            one_digit = 'q'
        elif first_two_digits == '53':
            one_digit = 'r'
        elif first_two_digits == '54':
            one_digit = 's'
        elif first_two_digits == '55':
            one_digit = 't'
        elif first_two_digits == '56':
            one_digit = 'u'
        elif first_two_digits == '57':
            one_digit = 'v'
        elif first_two_digits == '58':
            one_digit = 'w'
        elif first_two_digits == '59':
            one_digit = 'x'
        elif first_two_digits == '60':
            one_digit = 'y'
        elif first_two_digits == '61':
            one_digit = 'z'
        else:
            raise Exception('unrecognized sbdb designation format',des1)
        des2 = one_digit + last_four_digits
    if len(des1) > 6:
        century = des1[0:2]
        year = des1[2:4]
        space = des1[4]
        letter1 = des1[5]
        letter2 = des1[6]
        if len(des1) == 7:
            number = 0
        elif len(des1) == 8:
            number = des1[7]
        elif len(des1) == 9:
            number = des1[7:9]
        elif len(des1) == 10:
            number = des1[7:10]
        else:
            raise Exception('unrecognized sbdb designation format',des1)
        number = int(number)
        if space != ' ':
            raise Exception('unrecognized sbdb designation format',des1)
        if century == '18':
            century = 'I'
        elif century == '19':
            century = 'J'
        elif century == '20':
            century = 'K'
        else:
            raise Exception('unrecognized sbdb designation format',des1)
        if number == 0:
            code = '00'
        elif number < 10:
            code = '0' + str(number)
        elif number < 100:
            code = str(number)
        elif number < 620: # 619 is the max number with the MPC packing scheme
            numberstr = str(number)
            first_two_digits = numberstr[0:2]
            last_digit = numberstr[2]
            if first_two_digits == '10':
                one_digit = 'A'
            elif first_two_digits == '11':
                one_digit = 'B'
            elif first_two_digits == '12':
                one_digit = 'C'
            elif first_two_digits == '13':
                one_digit = 'D'
            elif first_two_digits == '14':
                one_digit = 'E'
            elif first_two_digits == '15':
                one_digit = 'F'
            elif first_two_digits == '16':
                one_digit = 'G'
            elif first_two_digits == '17':
                one_digit = 'H'
            elif first_two_digits == '18':
                one_digit = 'I'
            elif first_two_digits == '19':
                one_digit = 'J'
            elif first_two_digits == '20':
                one_digit = 'K'
            elif first_two_digits == '21':
                one_digit = 'L'
            elif first_two_digits == '22':
                one_digit = 'M'
            elif first_two_digits == '23':
                one_digit = 'N'
            elif first_two_digits == '24':
                one_digit = 'O'
            elif first_two_digits == '25':
                one_digit = 'P'
            elif first_two_digits == '26':
                one_digit = 'Q'
            elif first_two_digits == '27':
                one_digit = 'R'
            elif first_two_digits == '28':
                one_digit = 'S'
            elif first_two_digits == '29':
                one_digit = 'T'
            elif first_two_digits == '30':
                one_digit = 'U'
            elif first_two_digits == '31':
                one_digit = 'V'
            elif first_two_digits == '32':
                one_digit = 'W'
            elif first_two_digits == '33':
                one_digit = 'X'
            elif first_two_digits == '34':
                one_digit = 'Y'
            elif first_two_digits == '35':
                one_digit = 'Z'
            elif first_two_digits == '36':
                one_digit = 'a'
            elif first_two_digits == '37':
                one_digit = 'b'
            elif first_two_digits == '38':
                one_digit = 'c'
            elif first_two_digits == '39':
                one_digit = 'd'
            elif first_two_digits == '40':
                one_digit = 'e'
            elif first_two_digits == '41':
                one_digit = 'f'
            elif first_two_digits == '42':
                one_digit = 'g'
            elif first_two_digits == '43':
                one_digit = 'h'
            elif first_two_digits == '44':
                one_digit = 'i'
            elif first_two_digits == '45':
                one_digit = 'j'
            elif first_two_digits == '46':
                one_digit = 'k'
            elif first_two_digits == '47':
                one_digit = 'l'
            elif first_two_digits == '48':
                one_digit = 'm'
            elif first_two_digits == '49':
                one_digit = 'n'
            elif first_two_digits == '50':
                one_digit = 'o'
            elif first_two_digits == '51':
                one_digit = 'p'
            elif first_two_digits == '52':
                one_digit = 'q'
            elif first_two_digits == '53':
                one_digit = 'r'
            elif first_two_digits == '54':
                one_digit = 's'
            elif first_two_digits == '55':
                one_digit = 't'
            elif first_two_digits == '56':
                one_digit = 'u'
            elif first_two_digits == '57':
                one_digit = 'v'
            elif first_two_digits == '58':
                one_digit = 'w'
            elif first_two_digits == '59':
                one_digit = 'x'
            elif first_two_digits == '60':
                one_digit = 'y'
            elif first_two_digits == '61':
                one_digit = 'z'
            else:
                raise Exception('unrecognized sbdb designation format',des1)
            code = one_digit + last_digit
        else:
            raise Exception('unrecognized sbdb designation format',des1)
        des2 = century + year + letter1 + code + letter2
    return des2
#%% unpack designation
def unpack(designation):
    packed_designation = designation.lstrip()
    packed_designation = packed_designation.rstrip()
    Nt = len(packed_designation)
    if Nt == 5:
        onedigit = packed_designation[0]
        fourdigits = packed_designation[1:5]
        if onedigit == '0':
            onedigit = '0'
        elif onedigit == '1':
            onedigit = '1'
        elif onedigit == '2':
            onedigit = '2'
        elif onedigit == '3':
            onedigit = '3'
        elif onedigit == '4':
            onedigit = '4'
        elif onedigit == '5':
            onedigit = '5'
        elif onedigit == '6':
            onedigit = '6'
        elif onedigit == '7':
            onedigit = '7'
        elif onedigit == '8':
            onedigit = '8'
        elif onedigit == '9':
            onedigit = '9'
        elif onedigit == 'A':
            onedigit = '10'
        elif onedigit == 'B':
            onedigit = '11'
        elif onedigit == 'C':
            onedigit = '12'
        elif onedigit == 'D':
            onedigit = '13'
        elif onedigit == 'E':
            onedigit = '14'
        elif onedigit == 'F':
            onedigit = '15'
        elif onedigit == 'G':
            onedigit = '16'
        elif onedigit == 'H':
            onedigit = '17'
        elif onedigit == 'I':
            onedigit = '18'
        elif onedigit == 'J':
            onedigit = '19'
        elif onedigit == 'K':
            onedigit = '20'
        elif onedigit == 'L':
            onedigit = '21'
        elif onedigit == 'M':
            onedigit = '22'
        elif onedigit == 'N':
            onedigit = '23'
        elif onedigit == 'O':
            onedigit = '24'
        elif onedigit == 'P':
            onedigit = '25'
        elif onedigit == 'Q':
            onedigit = '26'
        elif onedigit == 'R':
            onedigit = '27'
        elif onedigit == 'S':
            onedigit = '28'
        elif onedigit == 'T':
            onedigit = '29'
        elif onedigit == 'U':
            onedigit = '30'
        elif onedigit == 'V':
            onedigit = '31'
        elif onedigit == 'W':
            onedigit = '32'
        elif onedigit == 'X':
            onedigit = '33'
        elif onedigit == 'Y':
            onedigit = '34'
        elif onedigit == 'Z':
            onedigit = '35'
        elif onedigit == 'a':
            onedigit = '36'
        elif onedigit == 'b':
            onedigit = '37'
        elif onedigit == 'c':
            onedigit = '38'
        elif onedigit == 'd':
            onedigit = '39'
        elif onedigit == 'e':
            onedigit = '40'
        elif onedigit == 'f':
            onedigit = '41'
        elif onedigit == 'g':
            onedigit = '42'
        elif onedigit == 'h':
            onedigit = '43'
        elif onedigit == 'i':
            onedigit = '44'
        elif onedigit == 'j':
            onedigit = '45'
        elif onedigit == 'k':
            onedigit = '46'
        elif onedigit == 'l':
            onedigit = '47'
        elif onedigit == 'm':
            onedigit = '48'
        elif onedigit == 'n':
            onedigit = '49'
        elif onedigit == 'o':
            onedigit = '50'
        elif onedigit == 'p':
            onedigit = '51'
        elif onedigit == 'q':
            onedigit = '52'
        elif onedigit == 'r':
            onedigit = '53'
        elif onedigit == 's':
            onedigit = '54'
        elif onedigit == 't':
            onedigit = '55'
        elif onedigit == 'u':
            onedigit = '56'
        elif onedigit == 'v':
            onedigit = '57'
        elif onedigit == 'w':
            onedigit = '58'
        elif onedigit == 'x':
            onedigit = '59'
        elif onedigit == 'y':
            onedigit = '60'
        elif onedigit == 'z':
            onedigit = '61'
        else:
            raise Exception('unrecognized packed designation format',packed_designation)
        unpacked_designation = onedigit + fourdigits
    elif Nt == 7:
        century = packed_designation[0]
        year = packed_designation[1:3]
        letter1 = packed_designation[3]
        code1 = packed_designation[4]
        code2 = packed_designation[5]
        letter2 = packed_designation[6]
        space = ' '
        if century == 'I':
            century = '18'
        elif century == 'J':
            century = '19'
        elif century == 'K':
            century = '20'
        else:
            raise Exception('unrecognized packed designation format',packed_designation)
        if code1 == '0' and code2 == '0':
            code1 = ''
            code2 = ''
        else:
            if code1 == '0':
                code1 = ''
            elif code1 == '1':
                code1 = '1'
            elif code1 == '2':
                code1 = '2'
            elif code1 == '3':
                code1 = '3'
            elif code1 == '4':
                code1 = '4'
            elif code1 == '5':
                code1 = '5'
            elif code1 == '6':
                code1 = '6'
            elif code1 == '7':
                code1 = '7'
            elif code1 == '8':
                code1 = '8'
            elif code1 == '9':
                code1 = '9'
            elif code1 == 'A':
                code1 = '10'
            elif code1 == 'B':
                code1 = '11'
            elif code1 == 'C':
                code1 = '12'
            elif code1 == 'D':
                code1 = '13'
            elif code1 == 'E':
                code1 = '14'
            elif code1 == 'F':
                code1 = '15'
            elif code1 == 'G':
                code1 = '16'
            elif code1 == 'H':
                code1 = '17'
            elif code1 == 'I':
                code1 = '18'
            elif code1 == 'J':
                code1 = '19'
            elif code1 == 'K':
                code1 = '20'
            elif code1 == 'L':
                code1 = '21'
            elif code1 == 'M':
                code1 = '22'
            elif code1 == 'N':
                code1 = '23'
            elif code1 == 'O':
                code1 = '24'
            elif code1 == 'P':
                code1 = '25'
            elif code1 == 'Q':
                code1 = '26'
            elif code1 == 'R':
                code1 = '27'
            elif code1 == 'S':
                code1 = '28'
            elif code1 == 'T':
                code1 = '29'
            elif code1 == 'U':
                code1 = '30'
            elif code1 == 'V':
                code1 = '31'
            elif code1 == 'W':
                code1 = '32'
            elif code1 == 'X':
                code1 = '33'
            elif code1 == 'Y':
                code1 = '34'
            elif code1 == 'Z':
                code1 = '35'
            elif code1 == 'a':
                code1 = '36'
            elif code1 == 'b':
                code1 = '37'
            elif code1 == 'c':
                code1 = '38'
            elif code1 == 'd':
                code1 = '39'
            elif code1 == 'e':
                code1 = '40'
            elif code1 == 'f':
                code1 = '41'
            elif code1 == 'g':
                code1 = '42'
            elif code1 == 'h':
                code1 = '43'
            elif code1 == 'i':
                code1 = '44'
            elif code1 == 'j':
                code1 = '45'
            elif code1 == 'k':
                code1 = '46'
            elif code1 == 'l':
                code1 = '47'
            elif code1 == 'm':
                code1 = '48'
            elif code1 == 'n':
                code1 = '49'
            elif code1 == 'o':
                code1 = '50'
            elif code1 == 'p':
                code1 = '51'
            elif code1 == 'q':
                code1 = '52'
            elif code1 == 'r':
                code1 = '53'
            elif code1 == 's':
                code1 = '54'
            elif code1 == 't':
                code1 = '55'
            elif code1 == 'u':
                code1 = '56'
            elif code1 == 'v':
                code1 = '57'
            elif code1 == 'w':
                code1 = '58'
            elif code1 == 'x':
                code1 = '59'
            elif code1 == 'y':
                code1 = '60'
            elif code1 == 'z':
                code1 = '61'
        unpacked_designation = century + year + space + letter1 + letter2 + \
            code1 + code2
    unpacked_designation = unpacked_designation.lstrip()
    unpacked_designation = unpacked_designation.rstrip()
    return unpacked_designation
#%%
import pandas as pd
import numpy as np
import time, os, json, urllib
from astroquery.jplhorizons import Horizons
dfmpc = pd.read_csv('/Users/iggymatheson/Documents_off_iCloud/mm22 mpcorb databases/MPCORB_20211202.csv',low_memory=False)
dfmpc = dfmpc.loc[dfmpc['Semimajor axis AU']>28]
dfmpc = dfmpc.loc[dfmpc['Semimajor axis AU']<200]
dfmpc.to_csv('MPCORB_20211202_reduced.csv',index=False)
sbdb_file=  'sbdb_query_results.csv'
df = pd.read_csv(sbdb_file,low_memory=False)
heliocentric_a_list = df['a'].tolist()
heliocentric_sigma_a_list = df['sigma_a'].tolist()
heliocentric_fractional_sigma_a_list = []
Nobj = df.shape[0]
for iobj in range(Nobj):
    a_obj = heliocentric_a_list[iobj]
    sigma_a_obj = heliocentric_sigma_a_list[iobj]
    if np.isnan(sigma_a_obj) or np.isnan(a_obj):
        heliocentric_fractional_sigma_a_list.append(999)
    else:
        heliocentric_fractional_sigma_a_list.append(np.abs(sigma_a_obj/a_obj))
print('done reading fractional semimajor axis uncertainty')
df['sigma_a/a'] = heliocentric_fractional_sigma_a_list
df = df.loc[df['sigma_a/a']<0.05]
dfmpc = pd.read_csv('MPCORB_20211202_reduced.csv',low_memory=False)
sbdb_des = df['pdes'].tolist()
mpc_des = dfmpc['Number or provisional designation in packed form'].tolist()
mpc_opp = dfmpc['Number of oppositions'].tolist()
sbdb_name = df['full_name'].tolist()
sbdb_des2 = []
sbdb_opp = []
Nobj_sbdb = df.shape[0]
Nobj_mpc = dfmpc.shape[0]
for iobj in range(Nobj_sbdb):
    des1 = sbdb_des[iobj]
    des2 = sbdb_des_to_packed_mpc_des(des1)
    sbdb_des2.append(des2)
    # if not des2 in mpc_des:
        # print(iobj,Nobj_sbdb,des1,des2)
        # time.sleep(5)
    if des2 in mpc_des:
        index = mpc_des.index(des2)
        sbdb_opp.append(mpc_opp[index])
    else:
        sbdb_opp.append(-999)
print('done reading opposition counts')
df['pdes2'] = sbdb_des2
df['Nopp'] = sbdb_opp
df = df.loc[df['Nopp']>=3]
# now retrieve Horizons barycentric elements for 2022 January 21 00:00:00, ie JD 2459600.50000
Nobj_sbdb = df.shape[0]
sbdb_des = df['pdes'].tolist()
packed_list = []
unpacked_list = []
for iobj in range(Nobj_sbdb):
    des1 = sbdb_des[iobj]
    des2 = sbdb_des_to_packed_mpc_des(des1)
    packed_list.append(des2)
    des3 = unpack(des2)
    unpacked_list.append(des3)
print('done unpacking designations')
M_list = []
Argperi_list = []
Node_list = []
Inc_list = []
e_list = []
n_list = []
a_list = []
q_list = []
period_list = []
center = '500@0' # center of the solar system
JD = 2459600.50000
for iobj in range(Nobj_sbdb):
    packed = packed_list[iobj]
    unpacked = unpacked_list[iobj]
    id_type = 'smallbody'
    if packed == 'D4340':
        unpacked = '9' # pluto barycenter shouldn't be in the list, but just in case
        id_type = 'majorbody'
    try:
        obj = Horizons(id=unpacked,location=center,epochs=JD,id_type=id_type)
        el = obj.elements()
        M = float(el['M'])
        w = float(el['w'])
        Omega = float(el['Omega'])
        incl = float(el['incl'])
        e = float(el['e'])
        n = float(el['n'])
        a = float(el['a'])
        q = float(el['q'])
        M = np.mod(M,360)
        w = np.mod(w,360)
        Omega = np.mod(Omega,360)
        period = float(el['P'])
        M_list.append(M)
        Argperi_list.append(w)
        Node_list.append(Omega)
        Inc_list.append(incl)
        e_list.append(e)
        n_list.append(n)
        a_list.append(a)
        q_list.append(q)
        period_list.append(period)
    except:
        raise Exception('Horizons returned an error for object',unpacked,iobj,Nobj_sbdb)
        M_list.append(-999)
        Argperi_list.append(-999)
        Node_list.append(-999)
        Inc_list.append(-999)
        e_list.append(-999)
        n_list.append(-999)
        a_list.append(-999)
        q_list.append(-999)
        period_list.append(-999)
    print('done with object',iobj,Nobj_sbdb)
print('done retrieving Horizons elements')
df['Mean anomaly degrees barycentric 2022-01-21 00:00:00'] = M_list
df['Argument of perihelion J2000.0 degrees barycentric 2022-01-21 00:00:00'] = Argperi_list
df['Longitude of ascending node J2000.0 degrees barycentric 2022-01-21 00:00:00'] = Node_list
df['Inclination to ecliptic J2000.0 degrees barycentric 2022-01-21 00:00:00'] = Inc_list
df['Orbital eccentricity barycentric 2022-01-21 00:00:00'] = e_list
df['Mean daily motion degrees per day barycentric 2022-01-21 00:00:00'] = n_list
df['Semimajor axis AU barycentric 2022-01-21 00:00:00'] = a_list
df['Perihelion distance AU barycentric 2022-01-21 00:00:00'] = q_list
df['Orbital period days barycentric 2022-01-21 00:00:00'] = period_list
df['Unpacked MPC designation'] = unpacked_list
df['Packed MPC designation'] = packed_list
# now retrieve a JSON file for each object from the Small Body Database API
for iobj in range(Nobj_sbdb):
    des = packed_list[iobj]
    print(iobj,Nobj_sbdb,des)
    url = 'https://ssd-api.jpl.nasa.gov/sbdb.api?sstr='+des+'&full-prec=True&cov=mat'
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    outfile = 'ssd_json_' + des + '.json'
    with open(outfile,'w') as of:
        json.dump(data,of)
print('done retrieving JSON files')
# now parse JSON files to retrieve covariance and add it to the dataframe
cov11_list = []
cov12_list = []
cov13_list = []
cov14_list = []
cov15_list = []
cov16_list = []
cov21_list = []
cov22_list = []
cov23_list = []
cov24_list = []
cov25_list = []
cov26_list = []
cov31_list = []
cov32_list = []
cov33_list = []
cov34_list = []
cov35_list = []
cov36_list = []
cov41_list = []
cov42_list = []
cov43_list = []
cov44_list = []
cov45_list = []
cov46_list = []
cov51_list = []
cov52_list = []
cov53_list = []
cov54_list = []
cov55_list = []
cov56_list = []
cov61_list = []
cov62_list = []
cov63_list = []
cov64_list = []
cov65_list = []
cov66_list = []
for iobj in range(Nobj_sbdb):
    des = packed_list[iobj]
    print(iobj,Nobj_sbdb,des)
    if des == 'D4340': # Pluto, orbit is so well known, of course we include it
        2-2
        cov11_list.append(0)
        cov12_list.append(0)
        cov13_list.append(0)
        cov14_list.append(0)
        cov15_list.append(0)
        cov16_list.append(0)
        cov21_list.append(0)
        cov22_list.append(0)
        cov23_list.append(0)
        cov24_list.append(0)
        cov25_list.append(0)
        cov26_list.append(0)
        cov31_list.append(0)
        cov32_list.append(0)
        cov33_list.append(0)
        cov34_list.append(0)
        cov35_list.append(0)
        cov36_list.append(0)
        cov41_list.append(0)
        cov42_list.append(0)
        cov43_list.append(0)
        cov44_list.append(0)
        cov45_list.append(0)
        cov46_list.append(0)
        cov51_list.append(0)
        cov52_list.append(0)
        cov53_list.append(0)
        cov54_list.append(0)
        cov55_list.append(0)
        cov56_list.append(0)
        cov61_list.append(0)
        cov62_list.append(0)
        cov63_list.append(0)
        cov64_list.append(0)
        cov65_list.append(0)
        cov66_list.append(0)
    else: # SSD query for covariance matrix (no results from that if query is Pluto)
        infile = 'ssd_json_' + des + '.json'
        with open(infile) as json_file:
            data = json.load(json_file)
            cov = data['orbit']['covariance']['data']
            cov11_list.append(cov[0][0])
            cov12_list.append(cov[0][1])
            cov13_list.append(cov[0][2])
            cov14_list.append(cov[0][3])
            cov15_list.append(cov[0][4])
            cov16_list.append(cov[0][5])
            cov21_list.append(cov[1][0])
            cov22_list.append(cov[1][1])
            cov23_list.append(cov[1][2])
            cov24_list.append(cov[1][3])
            cov25_list.append(cov[1][4])
            cov26_list.append(cov[1][5])
            cov31_list.append(cov[2][0])
            cov32_list.append(cov[2][1])
            cov33_list.append(cov[2][2])
            cov34_list.append(cov[2][3])
            cov35_list.append(cov[2][4])
            cov36_list.append(cov[2][5])
            cov41_list.append(cov[3][0])
            cov42_list.append(cov[3][1])
            cov43_list.append(cov[3][2])
            cov44_list.append(cov[3][3])
            cov45_list.append(cov[3][4])
            cov46_list.append(cov[3][5])
            cov51_list.append(cov[4][0])
            cov52_list.append(cov[4][1])
            cov53_list.append(cov[4][2])
            cov54_list.append(cov[4][3])
            cov55_list.append(cov[4][4])
            cov56_list.append(cov[4][5])
            cov61_list.append(cov[5][0])
            cov62_list.append(cov[5][1])
            cov63_list.append(cov[5][2])
            cov64_list.append(cov[5][3])
            cov65_list.append(cov[5][4])
            cov66_list.append(cov[5][5])
df['e_e']    = cov11_list
df['e_q']    = cov12_list
df['e_tp']   = cov13_list
df['e_node'] = cov14_list
df['e_peri'] = cov15_list
df['e_i']    = cov16_list
df['q_e']    = cov21_list
df['q_q']    = cov22_list
df['q_tp']   = cov23_list
df['q_node'] = cov24_list
df['q_peri'] = cov25_list
df['q_i']    = cov26_list
df['tp_e']    = cov31_list
df['tp_q']    = cov32_list
df['tp_tp']   = cov33_list
df['tp_node'] = cov34_list
df['tp_peri'] = cov35_list
df['tp_i']    = cov36_list
df['node_e']    = cov41_list
df['node_q']    = cov42_list
df['node_tp']   = cov43_list
df['node_node'] = cov44_list
df['node_peri'] = cov45_list
df['node_i']    = cov46_list
df['peri_e']    = cov51_list
df['peri_q']    = cov52_list
df['peri_tp']   = cov53_list
df['peri_node'] = cov54_list
df['peri_peri'] = cov55_list
df['peri_i']    = cov56_list
df['i_e']    = cov61_list
df['i_q']    = cov62_list
df['i_tp']   = cov63_list
df['i_node'] = cov64_list
df['i_peri'] = cov65_list
df['i_i']    = cov66_list
df.to_csv('sbdb_query_results_addcov.csv',index=False)
print('done with covariance')
# now generate clones for each object
Nclones = 300
df = pd.read_csv('sbdb_query_results_addcov.csv',low_memory=False)
Nobj = df.shape[0]
packed_list = df['Packed MPC designation'].tolist()
a_list = df['Semimajor axis AU barycentric 2022-01-21 00:00:00'].tolist()
e_list = df['Orbital eccentricity barycentric 2022-01-21 00:00:00'].tolist()
i_list = df['Inclination to ecliptic J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
w_list = df['Argument of perihelion J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
W_list = df['Longitude of ascending node J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
M_list = df['Mean anomaly degrees barycentric 2022-01-21 00:00:00'].tolist()
q_list = df['Perihelion distance AU barycentric 2022-01-21 00:00:00'].tolist()
n_list = df['Mean daily motion degrees per day barycentric 2022-01-21 00:00:00'].tolist()
tp_list = [] # time of pericenter passage formatted as Julian day - we'll just use time SINCE pericenter passage, ie pericenter passage is time zero
e_e_list       = df['e_e'].tolist()
e_q_list       = df['e_q'].tolist()
e_tp_list      = df['e_tp'].tolist()
e_node_list    = df['e_node'].tolist()
e_peri_list    = df['e_peri'].tolist()
e_i_list       = df['e_i'].tolist()
q_e_list       = df['q_e'].tolist()
q_q_list       = df['q_q'].tolist()
q_tp_list      = df['q_tp'].tolist()
q_node_list    = df['q_node'].tolist()
q_peri_list    = df['q_peri'].tolist()
q_i_list       = df['q_i'].tolist()
tp_e_list      = df['tp_e'].tolist()
tp_q_list      = df['tp_q'].tolist()
tp_tp_list     = df['tp_tp'].tolist()
tp_node_list   = df['tp_node'].tolist()
tp_peri_list   = df['tp_peri'].tolist()
tp_i_list      = df['tp_i'].tolist()
node_e_list    = df['node_e'].tolist()
node_q_list    = df['node_q'].tolist()
node_tp_list   = df['node_tp'].tolist()
node_node_list = df['node_node'].tolist()
node_peri_list = df['node_peri'].tolist()
node_i_list    = df['node_i'].tolist()
peri_e_list    = df['peri_e'].tolist()
peri_q_list    = df['peri_q'].tolist()
peri_tp_list   = df['peri_tp'].tolist()
peri_node_list = df['peri_node'].tolist()
peri_peri_list = df['peri_peri'].tolist()
peri_i_list    = df['peri_i'].tolist()
i_e_list       = df['i_e'].tolist()
i_q_list       = df['i_q'].tolist()
i_tp_list      = df['i_tp'].tolist()
i_node_list    = df['i_node'].tolist()
i_peri_list    = df['i_peri'].tolist()
i_i_list       = df['i_i'].tolist()
for iobj in range(Nobj):
    des = packed_list[iobj]
    a_nom = a_list[iobj]
    e_nom = e_list[iobj]
    i_nom = i_list[iobj]
    peri_nom = w_list[iobj]
    node_nom = W_list[iobj]
    M_nom = M_list[iobj]
    q_nom = q_list[iobj]
    n_nom = n_list[iobj]
    tp_nom = M_nom / n_nom # time since pericenter passage in days is mean anomaly in degrees divided by mean motion in degrees per day
    e_e = e_e_list[iobj]
    e_q = e_q_list[iobj]
    e_tp = e_tp_list[iobj]
    e_node = e_node_list[iobj]
    e_peri = e_peri_list[iobj]
    e_i = e_i_list[iobj]
    q_e = q_e_list[iobj]
    q_q = q_q_list[iobj]
    q_tp = q_tp_list[iobj]
    q_node = q_node_list[iobj]
    q_peri = q_peri_list[iobj]
    q_i = q_i_list[iobj]
    tp_e = tp_e_list[iobj]
    tp_q = tp_q_list[iobj]
    tp_tp = tp_tp_list[iobj]
    tp_node = tp_node_list[iobj]
    tp_peri = tp_peri_list[iobj]
    tp_i = tp_i_list[iobj]
    node_e = node_e_list[iobj]
    node_q = node_q_list[iobj]
    node_tp = node_tp_list[iobj]
    node_node = node_node_list[iobj]
    node_peri = node_peri_list[iobj]
    node_i = node_i_list[iobj]
    peri_e = peri_e_list[iobj]
    peri_q = peri_q_list[iobj]
    peri_tp = peri_tp_list[iobj]
    peri_node = peri_node_list[iobj]
    peri_peri = peri_peri_list[iobj]
    peri_i = peri_i_list[iobj]
    i_e = i_e_list[iobj]
    i_q = i_q_list[iobj]
    i_tp = i_tp_list[iobj]
    i_node = i_node_list[iobj]
    i_peri = i_peri_list[iobj]
    i_i = i_i_list[iobj]
    state_nom = np.array([e_nom,q_nom,tp_nom,node_nom,peri_nom,i_nom])
    e_list_here = []
    q_list_here = []
    tp_list_here = []
    node_list_here = []
    peri_list_here = []
    i_list_here = []
    e_list_here.append(e_nom)
    q_list_here.append(q_nom)
    tp_list_here.append(tp_nom)
    node_list_here.append(node_nom)
    peri_list_here.append(peri_nom)
    i_list_here.append(i_nom)
    cov = np.array([[e_e,   e_q,   e_tp,   e_node,   e_peri,   e_i],\
                    [q_e,   q_q,   q_tp,   q_node,   q_peri,   q_i],\
                    [tp_e,  tp_q,  tp_tp,  tp_node,  tp_peri,  tp_i],\
                    [node_e,node_q,node_tp,node_node,node_peri,node_i],\
                    [peri_e,peri_q,peri_tp,peri_node,peri_peri,peri_i],\
                    [i_e,   i_q,   i_tp,   i_node,   i_peri,   i_i]])
    for iclone in range(Nclones):
        try:
            L = np.linalg.cholesky(cov)
            clone_1 = np.random.normal(0,1,6) # 6-element random Gaussian mean 0 stdev 1
            clone_2 = np.dot(L,clone_1) # Gaussian perturbations from zero correlated according to cov matrix
            clone_3 = clone_2 + state_nom # add perturbation to nominal orbit
            e_clone = clone_3[0]
            if e_clone < 0:
                e_clone = -e_clone
            q_clone = clone_3[1]
            tp_clone = clone_3[2]
            node_clone = clone_3[3]
            peri_clone = clone_3[4]
            i_clone = clone_3[5]
        except:
            e_clone = 9e9
            q_clone = 9e9
            tp_clone = 9e9
            node_clone = 9e9
            peri_clone = 9e9
            i_clone = 9e9
        e_list_here.append(e_clone)
        q_list_here.append(q_clone)
        tp_list_here.append(tp_clone)
        node_list_here.append(node_clone)
        peri_list_here.append(peri_clone)
        i_list_here.append(i_clone)
        print(iobj,Nobj,iclone,Nclones)
    outfile = 'sbdb_query_results_clones_' + des + '.csv'
    dictionary = {'e':e_list_here,
                  'q':q_list_here,
                  'tp':tp_list_here,
                  'node':node_list_here,
                  'peri':peri_list_here,
                  'i':i_list_here}
    df_out = pd.DataFrame(dictionary)
    df_out.to_csv(outfile,index=False)
print('done generating clones')
