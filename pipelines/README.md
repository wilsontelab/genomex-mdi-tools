---
published: false
---

## Stage 1 Pipelines overview

### Pipelines vs. programs

A **pipeline** (also called a workflow) is a single coordinated set
of data analysis instructions that can be called by name using the
mdi command line utility, e.g.,

```bash
mdi <pipeline> ...
```

Despite being called by a convenient and familiar command line utility, 
a pipeline is distinct from many standalone programs in
that a pipeline might have many inputs and outputs and typically makes 
calls to many installed programs. Thus, rather than being a program itself, 
a pipeline is a means of coordinating a reproducible set of calls to 
other programs, in a manner that can be easily repeated on new data sets.

### Pipeline actions

Each MDI pipeline must have one or more **actions** defined
in 'pipeline.yml' and organized into pipeline subfolders.
Actions are entered by users at the command line or in
'\<data\>.yml' configuration files.

```yml
# pipeline.yml (the pipeline's configuration file)
actions:
    actionName:
        ... # see other documentation for what goes here
```

```bash
# command line
mdi <pipeline> <actionName> ...
mdi myPipeline do ... # e.g., a single-action pipeline called 'myPipeline'
```

```yml
# <data>.yml (the specifications for analysis of a specific data set)
pipeline: myPipeline
options: 
    ... # see other documentation for what goes here
execute: # the list of actions to execute
    - do # 'do' is the standardized name for a single action
```

### Action steps

Each pipeline action might finally be executed as a series of 
sequential **steps**, organized into action subfolders.
Steps are executed in series from a single pipeline action call.
Their success can be independently monitored to allow a 
pipeline action to restart without repeating previously satisfied 
steps, e.g., if a server crashes during step 3, steps 1 and 2 
would not need to be repeated.

### Conda environments

All pipelines use conda to construct an appropriate execution
environment with proper versions of all required program
dependencies, for explicit version control, reproducibility
and portability. These might include any program called from 
the command line or a shell script to do data analysis work.

### Singularity containers

When implemented by a tool suite developer, any pipeline can
be wrapped into an optional Singularity container to provide complete
control over the entirety of the operating system and codebase
available to support the pipeline. The operating system, e.g.,
Ubuntu, and system libraries are specified in 'singularity.def',
while program dependencies are provided by conda environments
pre-installed into the container.

## Pipeline construction

Create one folder in '\<suite\>/pipelines' for each distinct data 
analysis pipeline carried in your suite. Each pipeline 
might define its own specific code and/or use code elements 
in the '\<suite\>/shared' folder.

Only one file is essential and must be present, called 'pipeline.yml'.

Optional files include (see the _template pipeline for usage):
- 'README.md', to document the pipeline 
- 'pipeline.pl', which can be used to set custom environment variables 
- 'singularity.def', if your pipeline will offer or require a container to run

### Pipeline structure and definition

Begin by editing file **pipeline.yml**, which is the configuration file that establishes your pipeline's identity, options, actions, etc. It dictates how users will provide information to your pipeline and where the pipeline will look for supporting scripts and definitions.

Edit file **README.md** to provide a summary of the purpose and usage of your new pipeline. You will usually want to link this file into your GitHub Pages documentation by adding Jekyll front matter.

#### Pipeline actions

Next, create a subfolder in your pipeline directory for each discrete **action**
defined in pipeline.yml. Many pipelines only require one action, which by convention is called 'do'. Alternatively, you might need multiple actions executed independently, e.g., a first action 'analyze' applied to individual samples followed by a second action 'compare' that integrates information from multiple samples.

By convention, the target script in an action folder is called '**Workflow.sh**' - 
it is the script that performs the work of the pipeline action.
There are no restrictions on exactly how Workflow.sh does its work. You may incorporate 
other code, make calls to programs, including calls to nested workflow managers such as snakemake, etc. Thus, one common pattern might be:

```bash
# pipelines/<pipeline>/<action>/Worflow.sh
TARGET_FILE=$DATA_NAME.XYZ
snakemake $SN_DRY_RUN $SN_FORCEALL \
    --cores $N_CPU \
    --snakefile $ACTION_DIR/Snakefile \
    --directory $TASK_DIR \
    $TARGET_FILE
checkPipe
```

### Output conventions

Data files written by a pipeline are always placed into a folder you specify using 
options '--output-dir' and '--data-name'. File names should always be prefixed with the value of option
'--data-name', such that pipeline output files follow the pattern:

```
<output-dir>/<data-name>/<data-name>.XXX
```

'--output-dir' and '--data-name' are thus universally required options, and their 
values must be the same between sequential actions applied to the same data. 

#### Environment variables

The job launcher sets a number of environment variables that are available
to your action script. 

First, all pipeline options are available, where an option named "--abc-def" 
sets environment variable "ABC_DEF", e.g., '--n-cpu' becomes N_CPU, etc.

Additionally, the launcher sets the following derivative environment variables,
which are useful for locating files and for other purposes:

| variable name | value | description |
|---------------|---------------|-------------|
| **TASK_DIR**          | $OUTPUT_DIR/$DATA_NAME | where output files should be placed |
| **DATA_FILE_PREFIX**  | $TASK_DIR/$DATA_NAME   | prefix to use for output file names |
| **PLOTS_DIR**         | $TASK_DIR/plots        | where output plots should be placed |
| **PLOT_PREFIX**       | $PLOTS_DIR/$DATA_NAME  | prefix to use for output plot names |
| **SUITE_NAME**        | | the name of the tool suite that carries the running pipeline |
| **PIPELINE_NAME**     | | the name of the running pipeline |
| **PIPELINE_ACTION**   | | the name of the running pipeline action being applied to $DATA_NAME |
| **TASK_PIPELINE_DIR** | $TASK_DIR/$PIPELINE_NAME            | status, code and log files specific to the running pipeline and task |
| **TASK_ACTION_DIR**   | $TASK_PIPELINE_DIR/$PIPELINE_ACTION | code and log files specific to the running action and task |
| **SUITES_DIR**        | $TASK_ACTION_DIR/suites                | the directory where working versions of accessible MDI code suites are found |
| **SUITE_DIR**         | $SUITES_DIR/$SUITE_NAME                | the working root directory of tool suite $SUITE_NAME |
| **PIPELINE_DIR**      | $SUITE_DIR/pipelines/$PIPELINE_NAME    | the root directory of pipeline $PIPELINE_NAME, which contains 'pipeline.yml' |
| **ACTION_DIR**        | $PIPELINE_DIR/$PIPELINE_ACTION         | the directory that contains the scripts for action $PIPELINE_ACTION, including 'Workflow.sh' |
| **SCRIPT_DIR**        | $ACTION_DIR                            | legacy synonym for above |
| **ACTION_SCRIPT**     |                                        | the primary action script, usually '$ACTION_DIR/Workflow.sh' (unless overridden) |
| **SCRIPT_TARGET**     | $ACTION_SCRIPT                         | legacy synonym for above |
| **MODULES_DIR**       | $SUITE_DIR/shared/modules              | the directory that contains all shared code modules in suite $SUITE_NAME |
| **LOGS_DIR**          | $TASK_ACTION_DIR/logs | where log files should be placed |
| **LOG_FILE_PREFIX**   | $LOGS_DIR/$DATA_NAME  | prefix to use for log file names |
| **TASK_LOG_FILE**     | $LOG_FILE_PREFIX.$PIPELINE_NAME.$PIPELINE_ACTION.task.log | the main log file for a running task; starts with job config YAML |
| **RAM_PER_CPU_INT**   | int($RAM_PER_CPU)                     | e.g., 1M becomes 1000000 |
| **TOTAL_RAM_INT**     | $RAM_PER_CPU_INT * $N_CPU             | RAM available to entire job, e.g., 1000000 * 4 = 4000000 |
| **TOTAL_RAM**         | $TOTAL_RAM_INT (as string)            | e.g., 4000000 becomes 4M |

## Pipeline execution

Once the MDI is installed, pipelines can be called as executables from the command line, 
either by (i) a direct call to the pipeline or (ii) using the manager 'submit' command to queue
a pipeline job via a cluster server's job scheduler.

```bash
mdi <pipeline> ...
mdi submit <data.yml> ...
```

Use the --help options to learn about the different ways to configure pipeline calls.

```bash
mdi --help
mdi <command> --help
mdi <pipeline> --help
mdi <pipeline> <action> --help
```

### Job configuration files (specifying pipeline options)

You may provide all options through the command line 
as you would for any typical program.  Use '--help' to
see the available options. 

However, we recommend instead writing a '\<data\>.yml'
configuration file to set options, and then providing
the path to that file to the pipeline. This makes
it easy to assemble jobs and to keep a history of what
was done, especially if you use our job manager.

Config files are valid YAML files, although the interpreter
we use to read them only processes a subset of YAML features.
[Learn more about YAML on the internet](https://www.google.com/search?q=yaml+basics), 
or just proceed, it is intuitive and easy.

### Config file templates

To get a template to help you write your config file use:

```bash
mdi <pipelineName> template --help
mdi <pipelineName> template -a -c
```

In general, the syntax is:

```yml
# <data>.yml
---
pipeline: [suiteName/]pipelineName[:suiteVersion]
variables:
    VAR_NAME: value
pipelineAction:
    optionFamily:
        optionName1: $VAR_NAME # a single keyed option value
        optionName2:
            - valueA # an array of option values, executed in parallel
            - valueB
execute:
    - pipelineAction
```

The 'suiteName' and 'version' components of the pipeline declaration are optional, however, including suiteName can improve clarity and ensure that
you are always using the tool you intend. If provided, the version designation should be either:
- a suite release tag of the form 'v0.0.0'
- 'latest', to use the most recent release tag [the default]
- 'pre-release', to use the development code at the tip of the main branch
- for developers, the name of a code branch in the git repository

As a convenience for when you get tired of have many files
with the same option values (e.g., a shared data directory), you may
also create a file called 'pipeline.yml' or '\<pipelineName\>.yml'
in the same directory as '\<data\>.yml'. Options will be read
from 'pipeline.yml' first, then '\<data\>.yml', then finally
from any values you specify on the command line, with the last 
value that is read taking precedence, i.e., options specified on the 
command line have the highest precedence.

### Common workflow actions and options

Each pipeline supports parallel processing via options '--n-cpu' and '--ram-per-cpu'.

Option '--dry-run' allows a test of the action and options configuration prior to actual execution.

Many pipelines are designed to have staged execution with success monitoring. Commands 
'status' and 'rollback' provide support for monitoring and manipulating pipeline
execution status, repeating steps, etc. Pipelines that fail part-way through can be 
re-launched from the failure point by re-calling the initial action.

## Stage 1 versioning

### Pipeline versions

Individual pipeline versioning is optional but recommended as it will
help users to confidently access legacy versions of your code to analyze 
their data according to some previous standard, e.g., to ensure consistency 
between older and newer data sets.

Declaring pipeline versions is simple: just add a proper semantic version
declaration to pipeline.yml and update it prior to committing new code. 
It is not necessary to create Git tags for pipeline versions.

```yml
# pipelines/<pipeline>/pipeline.yml
pipeline:
    name: myPipeline
    description: "Description of myPipeline"
    version: v0.0.0
```

### Required versions of tools suite dependencies

If your pipeline uses code modules from external tool suites, you may
wish to specify the required versions of those external suites.
This is useful if you don't wish to adjust your pipeline to account for a
breaking change made in an external tool suite.  Declare such version
requirements as follows, replacing 'suiteName' with the name of the
external tool suite.

```yml
# pipelines/<pipeline>/pipeline.yml
suiteVersions:
    suiteName: v0.0.0
```

If you do not provide a version for an external tool suite,
the latest version of that suite will be used.

If you only use pipeline code from within your own tool suite, the 
suiteVersions dictionary can be omitted from pipeline.yml.
