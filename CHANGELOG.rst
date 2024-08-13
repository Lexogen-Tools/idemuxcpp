=======
History
=======

0.3.0 (2024-09-13)
------------------
* Change parallelization to async tasks
* Choose correction table with smallest compatible barcode set
- Remove parameter --writer-buffer-gb
+ Add parameter --gzip-block-size

0.2.0 (2023-11-07)
------------------
* Auto-detect barcode orientation (--auto-detect)
* Added additional mapping tables (for Aviti)

0.1.9 (2021-10-07)
------------------
* Limit open file handles proportional to number of writing threads
* New reader buffer and writer buffer parameter

0.1.8 (2021-03-02)
------------------
* added comma in correction count csv

0.1.7 (2021-28-01)
------------------
* bug fix: renamed i5 rc files so they get loaded properly

0.1.6 (2020-11-30)
------------------
* fixed correction counter indices
* set limit for pe mode
* updated sample sheet reader for windows files

0.1.5 (2020-11-18)
------------------
* increase maximal open file limit if required
* fix i1 only mode

0.1.4 (2020-11-09)
------------------

* added single end mode
* extended sample sheet
* allow mixed barcode lengths

0.1.3 (2020-11-03)
------------------

* bugfix zip file reader
* optional correction count output file

0.1.2 (2020-09-10)
------------------

* Second release.

0.1.0 (2020-09-10)
------------------

* First release.
