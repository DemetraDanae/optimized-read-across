# optimized-read-across

In the present study, an optimized read-across method for the toxicity prediction of engineered nanomaterials is developed. The proposed method combines the analogue and the grouping approach using one or more similarity criteria. The main advantage of the proposed algorithm is that there is no need of prior grouping hypothesis to proceed with the predictions, as it is produced as an outcome of the method. The whole procedure is automated and only the balance between the level of predictive accuracy and the number of predicted samples is user-defined. The number of similarity criteria used in the analogues selection is also user-defined.   

The repository contains the main code for the read-across optimization with one similarity criterion as well as the necessary functions. The presented code can be easily extended for two ore more similarity criteria. The code was developed using Matlab<sup>TM</sup> however, it is compatible with GNU Octave. In the GNU Octave version, subtle changes (marked as comments) have been made to the original code for compatibility reasons. The intereseted users can also comment in/out the commands concerning the dataset partitioning and external validation.

Additional information about the employed methods and an detailed case study, can be found in the corresponding <a href=https://pubs.rsc.org/en/content/articlelanding/2019/na/c9na00242a#!divAbstract>publication</a>: Varsou <i>et al</i>. (2019), "Read-across predictions of nanoparticle hazard endpoints: a mathematical optimization approach". 

<a href="https://zenodo.org/badge/latestdoi/168539471"><img src="https://zenodo.org/badge/168539471.svg" alt="DOI"></a>

# Dataset
The proposed read-across method is demonstrated on data derived from the publication of Walkey <i>et al.</i> (2014) in the field of nanoinformatics. The dataset dataset consists of 84 gold anionic and cationic nanoparticles with known cell association with human A549 cells (toxicity index), 40 physicochemical descriptors and 129 protein corona fingerprints (PCF, biological descriptors). This dataset was used by Varsou <i>et al.</i> (2017) in the <b>toxFlow</b> web application, which performs Gene Set Variation Analysis (GSVA) of omics data coupled with read-across toxicity prediction. The GSVA on the above PCF data identified 63 proteins as statistically significant.  In the file PCF giltered by GSVA.csv contains the 40 physicochemical descriptors and the 63 proteins produced by the GSVA filtering for the 84 gold nanoparticles.  

# License
This application is released under <a href="https://www.gnu.org/licenses/gpl.html"> GNU General Public License v.3</a>. 
```html
Copyright (C) 2019  Dimitra-Danai Varsou

This program is free software: you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
You should have received a copy of the GNU General Public License along with this program.  
If not, see here: http://www.gnu.org/licenses/.

Copyright (C) 2019  Dimitra-Danai Varsou

```
