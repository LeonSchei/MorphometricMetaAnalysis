Below is information on the files in this folder, the inputs to these functions, 
and how to run the processing described in Lefebvre et al. (2021).
Lefebvre, A, Her­rling, Zorndt, A, G, Becker, M, Krämer, K, Winter, C, 2021.
Mor­pho­logy of es­tu­ar­ine bed­forms, Weser Es­tu­ary, Ger­many. Earth Sur­face 
Pro­cesses and Land­forms https://doi.org/10.1002/esp.5243

There are three main functions to do the processing
get_crest_troughlines.m = calculates the position and properties of crestlines
tidal_steep_faces.m = calculates the position and properties of tidal steep faces
tidal_bedform_properties.m = calculates the properties of BEP bedforms

An example of how to run the processing is given in example_processing.m
Many thanks to to Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee, Standort Bremerhaven 
for allowing us to share the data

The main input is with a 3D bathymetry (xr,yr,zr) from a (constrained) tidal environment. 
It is assumed that 
- xr is longitude/crosswise position
- yr is the latitude/streamwise position with ebb direction in increasing streamwise values
- zr is the depth

For tidal_steep_faces.m, user-input threshold values are to find the steep faces
tidal_bedform_properties.m uses results from the two other functions

The software was developped for data from the Weser Estuary with a grid size fo 2 m. It will
probably have to be adapated for data with higher resolution!

If you have any questions feel welcome to contact me at alefebvre@marum.de
I would be happy to work together to improve or fix any errors. Any feedback will be welcome!

How to cite this material
https://doi.org/10.5281/zenodo.10715613

Alice Lefebvre (she/her)
MARUM - University of Bremen
https://www.marum.de/en/Dr-Alice-Lefebvre.html


+---------------+
|    LICENCE    |
+---------------+

Copyright (c) 2021 Alice Lefebvre

This program is free software: you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

If you use the software, please cite the paper
Lefebvre, A., Herrling, G., Becker, M., Zorndt, A., Krämer, K. & Winter, C. (2021).
Morphology of estuarine bedforms, Weser Estuary, Germany. Earth Surface
Processes and Landforms,1–15. https://doi.org/10.1002/esp.5243

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


