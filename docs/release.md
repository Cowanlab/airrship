# Release notes

We recommend using the latest available AIRRSHIP release. If you have already installed an earlier version then you may have to upgrade your installed version using pip. 

```bash
pip install airrship -U
```

## v0.1.0

Initial release.

## v0.1.1

Bug fix: the --shm_random option no longer produces an error.

## v0.1.2

Adjustments to ensure that proportional VDJ usage is accurately recreated from reference files.

## v0.1.3

Added option for including non-productive sequences in output.
Added --shm_multiplier option to control mutation rates. 

## v0.1.4

Yanked. Issue with Python 3.7 compatibility.

## v0.1.5

Included --species option to make simulating non-human species using user specified data easier. 
!!! warning
    Default values for --het have changed to 0 0 0. This differs from previous releases where the default behaviour was to use the greatest proportion of heterozygous positions possible in the inbuilt data.


