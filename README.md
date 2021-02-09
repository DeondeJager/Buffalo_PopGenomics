**High diversity, inbreeding and a dynamic Pleistocene demographic history revealed by African buffalo genomes**

Deon de Jager, Brigitte Glanzmann, Marlo Möller, Eileen Hoal, Paul van Helden, Cindy Harper, Paulette Bloomer  

This manuscript has been accepted for publication in Scientific Reports. I will post a link here as soon as it is available online.

# Summary
This repository contains the scripts used in the processing and analysis of the resequencing data for this study, as well as R scripts for producing some of the figures in the manuscript.

The scripts are organised into self-explanatory folders. Folders and scripts are numbered to more-or-less follow the same order as the Methods and Results section of the manuscript. 

Some analyses required more than one script: Related scripts for such analyses are grouped by number and letter, e.g. 04A and 04B.

For several analyses I had separate scripts for each sample (40 samples = 40 scripts) to run the analyses more efficiently on the high-performance computing cluster. For these analyses, I provided the script for only one sample (A_243_14) as an example and to avoid unnecessary clutter and repetition in the repository.

# Data
The raw sequence reads of all 40 buffalo genomes are available on the Sequence Read Archive (SRA) under the BioProject accession number PRJNA645266.

Data to produce the figures are also provided as comma-separated values (CSV) files in the `/03-Figures` folder. The exception to this are the PSMC output files used to produce Fig 4 and Supplementary Fig 5, since these files are >500 MB and GitHub is not a data repository, but they are available from me upon request. They can, of course, also be produced using the scripts provided in this repository.

Please feel free to contact me (email on profile page) if you have any questions or data requests.
