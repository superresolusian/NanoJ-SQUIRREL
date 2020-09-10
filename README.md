# NanoJ-SQUIRREL

Hello! This is the new home for the NanoJ-SQUIRREL source code that used to be hosted on Bitbucket. If you use SQUIRREL in your research, then please cite our paper: https://doi.org/10.1038/nmeth.4605

## Installation guide
If you're not interested in viewing the code, and just want to use the SQUIRREL plugin, then the best way to do this is via the Fiji updater:
* In Fiji, go to Help > Update...
* Click 'Manage Update Sites' in the bottom left window of the ImageJ Updater window that appears after Fiji is done checking for updates
* In the Manage Update Sites window that appears, tick the boxes next to 'NanoJ-Core' and 'NanoJ-SQUIRREL' and then press 'Close'
* Press 'Apply changes' in the ImageJ Updater window, and then restart Fiji
* You should now have the menu options for 'NanoJ-Core' and 'NanoJ-SQUIRREL' in your Plugins menu ready to go :)

For further instructions I strongly advise reading the manual, which is downloadable in the list of files above.

## Test data
If you want to give SQUIRREL a test run on a small data set, you can find the raw data from Figure 3 of the paper in the 'test_data' folder above. This is a small crop from a widefield reference image and three super-resolution reconstructions of that same crop reconstructed via different algorithms.

## Known issues
The underlying architecture of NanoJ relies on some OpenCL libraries for GPU-based processing that seem to disagree with some graphics cards. If your graphics card isn't compatible, then unfortunately you will find this out fairly quickly as Fiji will crash onces you click on any of the NanoJ-Core or NanoJ-SQUIRREL menu options :( We're working on a general fix for this, but in the interim I have made a version of SQUIRREL that does not do anything on the GPU. This is downloadable from the 'Releases' bar on the right hand side of the page, with a usage/installation note. Ironically, this is actually faster than the GPU-enabled version thanks to some clever maths from @uschmidt83, so don't worry about a loss of speed!
