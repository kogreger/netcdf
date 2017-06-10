import netCDF4
import zipfile
import shutil
import sys
import urllib2


def main():

	path_to_zip_file = 'C:\\Users\kgreger\Downloads\S3A_OL_1_EFR____20170609T093528_20170609T093828_20170609T113914_0179_018_307_2160_MAR_O_NR_002.zip'
	directory_to_extract_to = 'C:\\Users\kgreger\Downloads\output'

	zip_ref = zipfile.ZipFile(path_to_zip_file, 'r')
	zip_ref.extractall(directory_to_extract_to)
	zip_ref.close()

	# this is where the magic happens

	#urllib.urlretrieve ("https://coda.eumetsat.int/odata/v1/Products('7b9577fa-6c2c-4fbf-80bd-c01319d344d0')/$value", "data.zip")

	data = Dataset("Oa01_radiance.nc", "r", format="NETCDF4")
	print data.data_model


	# this is where the magic ends

	shutil.rmtree(directory_to_extract_to)

	return 0


if __name__ == "__main__":
    retval = main()
    sys.exit( retval )
