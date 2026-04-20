# Pynucastro test suite

Test suite to auto generate pynucastro comparisons to GridFire by copying the GridFire topology to pynucastro.

Note that this may take a while to run as each pynucastro network must run through JIT compilation. The JIT time is not counted to 
pynucastro evaluation time


To run use

```bash
python GridFireValidationSuite.py --suite HotStar_No_QSE CoolStar_No_QSE --pynucastro-compare --pync-engine="GraphEngine"
```

to see all options


```bash
python GridFireValidationSuite.py --help
```

Results will be saved in a directory as json files which you may parse to analyze 
