import argparse

__whatami__ = 'Match sources to the DECam candidates database '
__author__ = 'Igor Andreoni <andreoni@caltech.edu>'

def get_info_result(result):
    list_names = []
    list_info = []
    list_photometry = ["name | ra | dec | filter | mag | mag_err | mjd "]
    for r in result:
        name_cand, ra_cand, dec_cand, filter_cand, mag_cand, mag_err_cand, mjd_cand,\
            sg_gaia_cand, cluname_cand, cludistmpc_cand, cluseparcsec_cand, parallax_cand = \
            r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9], r[10], r[11]
        if (name_cand in list_names) == False:
            list_names.append(name_cand)
            if sg_gaia_cand == None or sg_gaia_cand == '':
                gaia_info = '?'
            if sg_gaia_cand == 0:
                gaia_info = 'No match'
            if sg_gaia_cand == 1:
                gaia_info = f'Star with parallax {parallax_cand}'
            if sg_gaia_cand == 3:
                gaia_info = f'Match in Gaia without parallax'
            list_info.append({'name': name_cand, 'RA':format(ra_cand, '.6f'), 'Dec':format(dec_cand, '.6f'), \
              'Gaia ': gaia_info, 'Nearest CLU galaxy name': cluname_cand, \
              'Nearest CLU galaxy distance (Mpc)': cludistmpc_cand,  \
              'Nearest CLU galaxy separation (arcsec)': cluseparcsec_cand\
              })
        list_photometry.append( f"{name_cand} | {format(ra_cand, '.6f')} | {format(dec_cand, '.6f')} | {filter_cand} | {format(mag_cand, '.3f')} | {format(mag_err_cand, '.3f')} | {format(mjd_cand, '.5f')}"  )

    return list_names, list_info, list_photometry


def query_coords(ra,dec,cursor,radius=0.1):
    ''' Query the database to search for matches at the given coordinates '''
    query = f"SELECT name,ra,dec,filter,mag,mag_err, mjd, sg_gaia, cluname, cludistmpc, cluseparcsec,parallax FROM object WHERE q3c_radial_query(ra, dec,{ra},{dec},{radius})"
    cursor.execute(query)
    result = cursor.fetchall()
    print('Query result:')
    if result == []:
            print(f'No match within {radius * 3600.} arcsec.')
            return
    else:
        list_names, list_info, list_photometry = get_info_result(result)
        print(f"\n Found {len(list_names)} match(es): {', '.join(list_names)}") 
        print("\n Photometry:")
        for p in list_photometry:
            print(p)
        print("\n Source(s) summary:")
        for i in list_info:
            for k in i.keys():
                 print(k, i[k])
            print()
    print('\n Done.')

    return

def query_name(name):
    ''' Query the database to search for matches to the given candidate name '''
    query = f"SELECT name,ra,dec,filter,mag,mag_err, mjd, sg_gaia, cluname, cludistmpc, cluseparcsec,parallax FROM object WHERE name = '{name}' "
    cursor.execute(query)
    result = cursor.fetchall()
    print('Query result:')
    if result == []:
            print(f'The candidate {name} is not present in the database')
            return
    else:
        list_names, list_info, list_photometry = get_info_result(result)
        print("\n Photometry:")
        for p in list_photometry:
            print(p)
        print("\n Source summary:")
        for i in list_info:
            for k in i.keys():
                 print(k, i[k])
            print()
    print('\n Done.')

    return


if __name__ == "__main__":

    #Parser 
    parser = argparse.ArgumentParser(description='Match coordinates or retrieve info of candidates stored on the DECam database.')

    parser.add_argument('--c', dest='radec', type=float, nargs=2, required=False, default=None, \
    help = 'Coordinates of the object you want to find. RA Dec (deg).  Example: 13.23422 -67.32677') 
    parser.add_argument('--r', dest='match_radius', type=float, required=False, default=1.5, \
    help = 'Match radius (arcsec)') 
    parser.add_argument('--n', dest='name', type=str, required=False, default=' ', \
    help = 'Name of the candidate that you want to find, for example DGabcd')

    #Entry checks
    args = parser.parse_args()
    if args.radec == None and args.name == ' ':
        print('>>> Please use --h to see the input required.')
        exit()

    import numpy as np
    from astropy.io import fits
    import psycopg2
    from astropy.io import ascii
    import pdb
    from astropy.table import Table
    import matplotlib.pyplot as plt

    #Read the secrets file and make the connection to the database
    info = ascii.read('./db_access.csv', format='csv')
    info_decam = info[info['db'] == 'decam']
    db = f"host={info_decam['host'][0]} port={info_decam['port'][0]} dbname={info_decam['dbname'][0]} user={info_decam['user'][0]} password={info_decam['password'][0]}"
    connection = psycopg2.connect(db)
    cursor = connection.cursor()

    if args.radec != None:
        ra, dec = args.radec[0], args.radec[1]
        radius = args.match_radius/3600.
        query_coords(ra,dec,cursor,radius=radius)

    if args.name != ' ':
        query_name(args.name)


