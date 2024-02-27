# SBNAF Infrared Database Main Script

Main script for the SBNAF Infrared Database. For more info see: https://ird.konkoly.hu/

## Description

The main file wich creates the input csv files for the [SBNAF Infared Database](https://ird.konkoly.hu/). For the input file format see 
the [release note](https://ird.konkoly.hu/data/releaseNotes/SBNAF_IRDB_public_release_note_2021February02.pdf).

## How to install

Just git clone the repository to your working directory:
```shell
$ git clone https://github.com/rszakats/SBNAF_IRDB
$ cd SBNAF_IRDB
```
Then:

```shell
virtualenv sbnaf
source sbnaf/bin/activate
pip install -r requirements.txt
```

## Running

```shell
python irdatabase.1.2.py --help
```

There are two command line arguments that must be provided:
--workdir: The location of your input files and the directory where you want your output files.
--inst: The instrument name. This can be: 
    - IRAS
    - HSO
    - AKARI
    - MSX
    - WISE

So for example:
```shell
python irdatabase.1.2.py --workdir=/data/irdata/ --inst=HSO
```

## Author

@author: R. Szakáts, [Konkoly Observatory](https://konkoly.hu/en), 2019-2024

## Credits

- [NASA JPL](https://ssd-api.jpl.nasa.gov/about/)
- [SBNAF](http://www.sbnaf.eu/)
- Vladimir Ignatev (progress.py)
- Victor Ali-Lagoa (BB_cc_*.py)

## License

This software is licensed under the [GNU GPL-3.0](LICENSE) License.