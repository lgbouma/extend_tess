#
# bash instructions to make Luke-analog tess sky plots.
# if you execute this script, it will install the "tessmaps" package to your
# local python environment. it will download ~200Mb of files as well, to a
# temporary directory.
#

mkdir temp
cd temp
git clone https://github.com/lgbouma/tessmaps
git clone https://github.com/lgbouma/extend_tess

# install "tessmaps" package to your local python environment. this package has
# a focal plane geometry model used by the plotting script.

cd tessmaps
python setup.py install 

# extend_tess.src.plot_extmission_field_positions.py is the plotting script.
# its usage is described in its docstring. first, download 2 cached files to
# make the plot that made it into the proposal. (otherwise,
# plot_extmission_field_positions.py takes ~5 minutes to do this -- shorter if
# you tune it to give fewer points and a rattier plot).

cd ../extend_tess/data
wget "https://www.dropbox.com/s/led02qjev990wy4/idea_14_final_shifted_coords_observed_forproposal.csv?dl=0" -O idea_14_final_shifted_coords_observed_forproposal.csv
wget "https://www.dropbox.com/s/qepc50zc0d2j0tx/idea_14_final_shifted_coords_observed_merged_forproposal.csv?dl=0" -O idea_14_final_shifted_coords_observed_merged_forproposal.csv

cd ../src
python plot_extmission_field_positions.py
