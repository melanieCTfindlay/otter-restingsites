# otter-restingsites
Supporting data and code for Findlay, Briers, Ingledew, White (2022) An evidence-based approach to identifying resting sites of Eurasian otter (Lutra lutra) from camera-trap and field-sign data

Supporting code and data for two analyses in the above paper. 

The first analysis looks at whether field signs of Eurasian otter differed between resting sites (dens) and non resting sites (structures visited by otters but not rested in). We analyse field data collected for approximately a year per site, from 6 otter resting sites and 20 non resting sites across the Tweed catchment in Scotland and use models to see if any single variable (presence of bedding, presence of a latrine, presence of spraint piles, presence of a run, count of spraint to 1m of the resting site/den entrance) has a strong association with resting sites. The code is in the file called Findlay, Briers, Ingledew, White 2022 and is in R. The data is within FIELDMEANSDATA.csv.

The second analysis simulates the camera-trapping data to find what theoretical duration of camera trapping would one have to undertake to have a 95% chance of detecting a day when an otter rests at a structure. Two approaches to the simulations are taken, the first uses the combined winter and following spring and simulates a single period of camera trapping during that period. The second approach simulates two equal periods of camera trapping, one in winter and one in spring to see if any economies can be made on camera trapping effort. The code is in the file called Findlay, Briers, Ingledew, White 2022 and is in R. The data is within site.period5ab.csv which can be used to run the single period duration code for two site periods, and the two equal sampling periods, again for two site periods.

