### Installation

```
python -m venv eos-venv
source eos-venv/bin/activate
pip install -e git+https://github.com/cra/simple-eos.git#egg=simple_eos
```

### Usage

After installation, `get_eos` console command should be available:

```
get_eos --datafile eos.dat
```

Input file eos\.dat should contain two columns: Volume/cell and Energy/cell. Volume is expected to be provided in Ang^3 and Energy in eV.

Something like this:

```
17.58	-.89579572E+01
18.77	-.95953780E+01
20.03	-.10042777E+02
21.33	-.10334275E+02
22.69	-.10495856E+02
30.37	-.10036497E+02
32.09	-.97898900E+01
```


