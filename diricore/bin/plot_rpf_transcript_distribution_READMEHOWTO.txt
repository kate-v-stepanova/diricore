Usage is as follows (this command was in fact used to make the plot you
attached in the first email):

minreads=100;
python plot_rpf_transcript_distribution.py \
     -o
KR20151104.2864_Noco_Torin.rpf_transcript_coverage.m${minreads}.pdf \
     -m ${minreads} \
'2864_1,../KR20141110/data/transcript_coordinate_rpf_counts/KR20141215.286
4_U2OS_Myc_Torin_Noco.txcoord_counts.hdf5,Control,#425D80'
\
'2864_5,../KR20141110/data/transcript_coordinate_rpf_counts/KR20141215.286
4_U2OS_Myc_Torin_Noco.txcoord_counts.hdf5,Torin,#806842'
\
'2864_3,../KR20141110/data/transcript_coordinate_rpf_counts/KR20141215.286
4_U2OS_Myc_Torin_Noco.txcoord_counts.hdf5,Noco,#428042'
\

The '-o' argument is required, the minimum-number of reads argument
defaults to 100 if not specified. The samples are specified as list of
comma-separated tuples of:
'<sample_id>,<HDF5 file with transcriptome-mapped RPFs>,<human readable
sample name>,<color for the line in the plot>'
The data of <sample_id> must obviously be present in the HDF5 file (in
fact the map_rpfs_to_transcriptome_positions.py script is what generates
the transcriptome-mapped HDF5 file given genome-mapped BAM files as
input).
