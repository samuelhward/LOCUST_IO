#wall_converter

import context
from classes.input_classes.wall import Wall as w


radii=['1.01','1.05','1.10','1.15','1.20','1.30','1.40','1.50']
transp_filenames='OMF157418.LIM'
ascot_filenames='input.wall_2d'
locust_filenames='LCFS_DIII-D.dat'

for radius in radii:

	filename=transp_filenames+'_'+radius
	wall=w('ID','ufile',filename)
	wall.dump_data('ASCOT',ascot_filenames+'_'+radius)
	wall.dump_data('LOCUST_2D',locust_filenames+'_'+radius)