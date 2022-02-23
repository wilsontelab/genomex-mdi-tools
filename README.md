# wilsontelab genomex-mdi-tools

The [Michigan Data Interface](https://midataint.github.io/) (MDI) 
is a framework for developing, installing and running a variety of 
HPC data analysis pipelines and interactive R Shiny data visualization 
applications within a standardized design and implementation interface.

Data analysis in the MDI is separated into 
[two stages of code execution](https://midataint.github.io/docs/analysis-flow/) 
called Stage 1 HPC **pipelines** and Stage 2 web applications (i.e., **apps**).
Collectively, pipelines and apps are referred to as **tools**.

## Repository contents

This **genomex** repository contains public 
tools, modules, environments, and options 
for MDI pipelines and apps that perform genomic data analysis,
from the 
[Thomas Wilson laboratory](https://wilsonte-umich.github.io)
at the University of Michigan.

### Stage 1 pipelines to establish genomics resources

This repo first carries the following pipelines to prepare 
your computer with required genomics resources:

- **download** = obtain genome and sample data files from the internet
- **genmap** = run the genmap utility to create genome mappability files
- **prepareBins** = create a file of fixed-width genome bins with associated metadata

The use these pipelines, install this tool suite as described below.

### Shared modules, environments, and option families

This repo further carries a series of shared elements to support
other Wilson laboratory pipelines and apps. You may also use
any of them in your own tool suites by following patterns like:

```yml
# _config.yml, in your tool suite root directory
suite_dependencies:
  - wilsontelab/genomex-mdi-tools
```

```yml
# pipeline.yml
actions:
  action:
    condaFamilies:
      - genomex-mdi-tools/XXX
```

---
## Quick Start 1: suite-centric installation

In suite-centric mode, you will:
- clone this tool suite repository
- call its 'install.sh' script to create a suite-specific MDI installation
- call its 'genomex' utility to use its tools

### Install this tool suite

```bash
git clone https://github.com/wilsontelab/genomex-mdi-tools.git
cd genomex-mdi-tools
./install.sh
```

### Execute a Stage 1 pipeline from the command line

For help, call the 'genomex' utility with no arguments - 
add the suite directory to PATH to run the tool suite from any directory.

```bash
./genomex
```

---
## Quick Start 2: mdi-centric installation

In mdi-centric mode, you will:
- clone and install the MDI
- add this tool suite (and potentially others) to your configuration file
- re-install the MDI to add this tool suite to your MDI installation
- call the 'mdi' utility to use its tools

### Install the MDI framework

Please read the 'install.sh' menu options and the 
[MDI installer instructions](https://github.com/MiDataInt/mdi.git) to decide
which installation option is best for you. Briefly, choose option 1
if you will only run Stage 1 HPC pipelines from your installation.

```bash
git clone https://github.com/MiDataInt/mdi.git
cd mdi
./install.sh
```

### Add an alias to .bashrc (optional)

These commands will help you create a permanent named alias to the 'mdi'
target script in your new installation.

```bash
./mdi alias --help
./mdi alias --alias genomex # change the alias name if you'd like 
`./mdi alias --alias genomex --get` # activate the alias in the current shell too
genomex
```

### Add this tool suite to your MDI installation

Edit file 'mdi/config/suites.yml' as follows:

```yml
# mdi/config/suites.yml
suites:
    - wilsontelab/genomex-mdi-tools
```

and re-run <code>install.sh</code>, which will go faster
the second time. 

### Execute a Stage 1 pipeline from the command line

For help, call the 'mdi' utility with no arguments.

```bash
genomex
```
