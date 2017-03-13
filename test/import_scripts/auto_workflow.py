from subprocess import call
import sys
import os
from datetime import date, timedelta

#REQUIRED ARGUMENTS:  data type, year startm, year end, start month, 

#Create dates that can be fed into batch file
data_type = sys.argv[1]
year_start = int(sys.argv[2])
year_end = int(sys.argv[3])
month = int(sys.argv[4])

start_date = date(year_start,month,1)

if (data_type == 'ERA'):
  end_date = date((year_end+1),1,1)
elif (data_type == 'MERRA'):
  end_date = date(year_start, (month+1), 1)
else:
  print "Invalid data choice."
  sys.exit(1)

add_delta = timedelta(days=30)
mod_delta = timedelta(days=1)

month_start_list = []
month_end_list = []
while start_date != end_date:
  month_start_list = month_start_list + [start_date]
  current_date = start_date + add_delta
  while current_date.month != start_date.month:
    current_date = current_date - mod_delta
  month_end_list = month_end_list + [current_date]
  start_date = current_date + mod_delta
print "List is " + str(len(month_start_list)) + " long."

#Make folder to place files
directory = '/scratch2/scratchdirs/marielp/'

#Create strings of date ranges for batch import
for i in range(0, len(month_start_list)):
#date string
  sdate = month_start_list[i]
  edate = month_end_list[i]
  subdir = '%s%s_%s/'%(directory, data_type, sdate.strftime('%Y'))
  dir_check = 0
  dir_check = int(os.path.exists(subdir))
  if dir_check ==0:
    print "Making directory " + subdir
    call(['mkdir', subdir])
#name string
  file_it = '%s_%s'%(data_type, sdate.strftime("%Y_%m"))

#STEP 1: DATA IMPORT
##ERA DATA IMPORT: INITIALIZING FILE CONDITIONS
  if (data_type == 'ERA'):
    str_date = '%04d%02d%02d/to/%04d%02d%02d'%(sdate.year, sdate.month, sdate.day, edate.year, edate.month, edate.day)
    filename_grib = subdir + file_it + ".grib"
    print 'Grib file path is %s.'%(filename_grib)
    filename_nc = filename_grib.replace(".grib", ".nc")
    mod_filename = filename_nc.replace(".nc", "_mod.nc")
    print 'Final filename is %s.'%(mod_filename)

##MERRA DATA IMPORT: INITIALIZING FILE CONDITIONS
####NOTE: THIS ONLY CALLS ONE MONTH'S WORTH OF DATA AS OPPOSED TO ERA (1 YR)
  elif (data_type == 'MERRA'):
    mon = sdate.strftime('%b').lower() 
    str_date = '0z%d'%(sdate.day) + mon + '%04d'%(sdate.year) + ' 23z%d'%(edate.day) + mon + '%04d'%(edate.year)
    filename_nc_noext = subdir + file_it
    filename_nc = filename_nc_noext + '.nc4'
    mod_filename = filename_nc4.replace(".nc4", "_mod.nc4")

  file_exists_check = int(os.path.exists(mod_filename))
  if (file_exists_check != 0):
    print "File %s already exists."%(mod_filename)
  else:
    if (data_type == 'ERA'):
####import file
      import_grib = 1
      import_grib = call([sys.executable, '/global/homes/m/marielp/batch_import.py', str_date, filename_grib])
      if (int(import_grib) != 0):
        print "Failed to import. Exiting."
        os._exit(1)
      print "Finished importing grib file. Now converting to netCDF.\n"
####convert grib to netCDF
      grib2nc = 1
      grib2nc = call(['ncl_convert2nc', filename_grib, '-o', subdir])
      if (int(grib2nc) != 0):
        print "Failed to convert. Exiting."
        os._exit(1)
      print "Successfully converted. Removing grib file."
####Remove grib file
      remove = 1
      remove = call(['rm', filename_grib])
      if (int(remove) != 0):
        print "Could not find specified file. Exiting."
        os._exit(1)
    elif (data_type == 'MERRA'):
      import_nc = 1
      import_nc = call(['/global/project/projectdirs/vacet/Wehner/marielle/grads-2.0.1.oga.1/Contents/lats4d.sh', \
 '-i', 'http://goldsmr3.sci.gsfc.nasa.gov:80/dods/MAI6NPANA', '-v', \
 '-vars', 't u v', '-time', str_date, '-lat', '-90 90', '-lon', '-180 180', \
 '-levs', '150 200 250 300 350 400 450 500', \
 '-format', 'netcdf4', '-o', filename_nc_noext])
      if (int(import_nc) != 0):
        print('Failed to import.Exiting.')
        os._exit(1)
      print 'Successfully imported MERRA data.\n'
    print "Now formatting netCDF file to VisIt readable format.\n"
    file_format = 1
    if (data_type == 'ERA'):
      file_format = call([sys.executable, '/global/homes/m/marielp/format_nc_ERA.py', filename_nc])
    elif(data_type == 'MERRA'):
      file_format = call([sys.executable, '/global/homes/m/marielp/format_nc_MERRA.py', filename_nc])
    if (int(file_format != 0)):
      print "Failed to format. Exiting."
      os._exit(1)
    remove = call(['rm', filename_nc])
    check_files = 1
    check_files = call(['ls', mod_filename])
    if (check_files == 0):
      print "Successfully completed import and format. File is located at " + mod_filename + ".\n\n"
#STEP 2: CREATE BATCH FILE
  call([sys.executable, '/global/homes/m/marielp/batch_create.py', mod_filename, data_type, str(sdate.year), str(sdate.month)])
#  if (data_type == 'ERA'):
#    call([sys.executable, '/global/homes/m/marielp/ERA_batch_create.py', mod_filename, data_type, str(sdate.year), str(sdate.month)])
sys.exit(0)
