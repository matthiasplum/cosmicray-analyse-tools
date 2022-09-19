# cosmicray-analysis-tools

This python library should help to perform cosmic ray mass composition analysis using template pdf. Similar to https://arxiv.org/abs/1906.04317

##### Version 0.0.2
* extended likelihood template fit method
* unbinnend and binned fit method available
* event number per template are constraint to max number of data
* add setup.py

##### ToDO:
* Remove probfit dependecy as it's currently not 100% combatible with iminuit 2.16.0
* Add Numba support and general speedup code for unbinned dataset fits
* Test weighted event data sets 