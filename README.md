# B0_in_Dstr_Dstr [04.09.2022]
The Search of the B0 decay in the Dstr and anti Dstr. The repository contains the Xb_frame file, the MySelector file and the code for creating graphics.

## The data used:
The whole BParking dataset.

## The files

### Xb_frame
File on c++ with the code for the reconstruction of the particle which decays on D*(2010) and D*(2010). The file is launched on crab. The data is save on cernbox.

### MySelector_DstDst_2018_v2
The handler of events created after the Xb_frame. The data is saved on the cernbox. The file is launched on the cms server.

### Graphs
The program for handling the events after the MySelector. Returns the histograms with mass distributions.