# Template for running Gabor MP Analysis
%INPUT_OUTPUT
Numb_inputs=1
Numb_outputs=1
Mode=parallel
%INPUT
Path=/Users/crfetsch/Documents/MATLAB/MP/data/test/ImportData_SIG/
Header_file=sig.hdr
Calibrate=1
Numb_points=4096
Shift_points=4096
%OUTPUT
Path=/Users/crfetsch/Documents/MATLAB/MP/data/test/GaborMP/
All_chans=2
Numb_chans=2
Start_chan=1
Start_chan_no=0
Header_file=book.hdr
Type=book
File_format=double
Name_template=mp#.bok
Max_len=600
Chans_per_file=-1
%GABOR_DECOMPOSITION
Max_Iterations=100
