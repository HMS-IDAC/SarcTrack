NOTE: a new version is available; please use that one instead: https://github.com/HMS-IDAC/SarcTrack2

# SarcTrack

Sarcomere Detection and Contraction Analysis

Step 1:
    'dwProcessFrame' to test parameters

Step 2:
    fix parameters in 'dwProcessVideo'

Step 3:
    'dwProcessFolder' on a folder with one or two movies

Step 4:
    'dwCheckResults' on outputs from step 3; adjust parameters and run again if necessary

Step 5:
    'dwProcessFolder' on folder with arbitrary number of movies

Step 6:
    'dwAggegateTables' on folders whose tables are to be aggregated

Sample dataset: https://www.dropbox.com/s/k1p65tnfeixp9q0/SarcTrackSampleVideos.zip?dl=0

To generate and analyze synthetic movies, as well as to understand the code, use dwGenerateSynthMovie.m, dwProcessSynthFrame.m, and dwProcessSynthMovie.m
