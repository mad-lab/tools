





Overview
--------


+ This is a software that can be used to analyze Tn-Seq datasets. It includes various statistical calculations of essentiality of genes or genomic regions (including conditional essentiality between 2 conditions). These methods were developed and tested as a collaboration between the Sassetti lab (UMass) and the Ioerger lab (Texas A). 

.. thumbnail:: http://saclab.tamu.edu/essentiality/transit/images/tutorial_glyc_ctrl.png
   :width: 70%
   :align: center

|

+ TRANSIT assumes you have already done pre-processing of raw sequencing files (.fastq) and extracted read counts (at each TA dinucleotide). Actually, the `current protocol <http://www.springer.com/biomed/human+genetics/book/978-1-4939-2397-7>`_ utilizes internal barcodes that can be used to reduce raw read counts to unique template counts, and this this is the intended input to TRANSIT. The input for TRANSIT consists of `.wig <http://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ files, which simply list the coordinate of each TA sites in the genome and the number of templates observed.

|

+ There are various methods available for pre-processing (converting .fastq files to .wig files). You might have your own scripts (if so, massage the data into .wig format), or you might get the scripts used in the Sassetti lab. For convenience, we are including a separate tool called `TPP <http://saclab.tamu.edu/tom/TPP.html>`_ (Tn-Seq Pre-Processor) with this distribution that encodes the way we process .fastq files in the Ioerger lab. It's a complicated process with many steps (removing transposon prefixes of reads, mapping into genome, identifying barcodes and reducing read counts to template counts).

|

+ Most of the analysis methods in TRANSIT require an **annotation** to know the gene coordinates and names. This is the top file input in the GUI window. The annotation has to be in a somewhat non-standard format called a ".prot_table". If you know what you are doing, it is easy to convert annotations for other organsims into .prot_table format. But for convenience, we are distributing the prot_tables for 3 common versions of the H37Rv genome: H37Rv.prot_table (NC_000962.2, from Stewart Cole), H37RvMA2.prot_table (sequenced version from the Sassetti lab), and H37RvBD.prot_table (sequenced by the Broad Institute). All of these are slightly different, and it is **critical** that you use the same annotation file as the reference genome sequence used for mapping the reads (during pre-processing).

|

+ There are 2 main types of essentiality analyses: individaul, comaparative. In individual analysis, the goal is to distinguish essential vs. non-essential in a single growth condition, and to assess the statistical significance of these calls. Two methods for this are the Gumbel method and the HMM. They are computationally distinct. The Gumbel method is looking for significant stretches of TA sites lacking insertions, whereas the HMM looks for regions where the mean read count is locally suppresed or increased. The HMM can detect 'growth-advantaged' and 'growth-defect' regions. The HMM is also a bit more robust on low-density datasets (insertion density 20-30%). But both methods have their merits and are complementary. For compararative analysis, TRANSIT uses 're-sampling', which is analogous to a permutation test, to determine if the sum of read counts differs significantly between two conditions. Hence this can be used to identify conditionally essential regions and quantify the statistical significance.

|

+ TRANSIT has been designed to handle multiple replicates. If you have two or more replicate dataset of the same libarary selected in the same condition, you can provide them, and more of the computational methods will do something reasonable with them.

|

+ For those methods that generate p-values, we often also calculate adjusted p-value (or 'q-values') which are corrected for multiple tests typically the Benjamini-Hochberg procedure. A typical threshold for significance would be q<0.05 (not p<0.05).


+ It is important to understand the GUI model that TRANSIT uses It allows you to load up datasets (.wig files), select them, choose an analysis method, set parameters, and start the computation. It will generate **output files** in your local directory with the results. These files can then be loaded into the interface and browser with custom displays and graphs. The interface has 3 main windows or sections: 'Control Samples', 'Experimental Samples', 'Results Files.' The first two are for loading input files ('Control Samples' would be like replicate datasets from a reference condition, like in vitro, rich media, etc.; 'Experimental Samples' would be where you would load replicates for a comparative conditions, like in vivo, or minimmal media, or low-iron, etc.) The 'Results Files' section is initially empty, but after a computation finishes, it will automatically be populated with the corresponding output file. See the 'Tutorial' section below in this documentation for an illustraion of the overall process for a typical work-flow.

|

+ TRANSIT incorporates many interesting ways of looking at your data.

|

    + Track view shows you a visual representation of the read counts at each site at a locus of interest (for selected datasets) somewhat like IGV.
    
.. thumbnail:: http://saclab.tamu.edu/essentiality/transit/images/track_view.png
   :width: 70%
   :align: center
|
    + Scatter plots can show the correlation of counts between 2 datasets.

.. thumbnail:: http://saclab.tamu.edu/essentiality/transit/images/scatter.png
   :width: 70%
   :align: center


|    
    + Volcano plots can be used to visualize the results of resampling and assess the distribution between over- and under-represented genes in condition B vs. condition A. In addition you can look at histogram of the re-sample distributions for each gene.

.. thumbnail:: http://saclab.tamu.edu/essentiality/transit/images/result_volcano_graph.png
   :width: 70%
   :align: center


.. thumbnail:: http://saclab.tamu.edu/essentiality/transit/images/result_table_histogram.png
   :width: 70%
   :align: center


+ Most of the methods take a few minutes to run. (it depends on parameters, CPU clock speed, etc., but the point is, a) these calculations are complex and not instaneous, but b) we have tried to implement it so that they don't take hours)


+ Note: in this version of TRANSIT, most of the methods are oriented toward gene-level analysis. There are methods for analyzing essentiality of arbitrary genomic regions (e.g. sliding windows, HMMs...). We plan to incorporate some of these in future versions.





|

Installation
------------
TRANSIT can be downloaded from the public GitHub server,
`http://github.com/mad-lab/transit <http://github.com/mad-lab/transit>`_. It is released under a GPL
License. It can be downloaded with git as follows:

::

    
    
    git clone https://github.com/mad-lab/transit/
    

TRANSIT is python-based You must have python installed (installed by
default on most systems). In addition, TRANSIT relies on some python
packages/libraries/modules that you might need to install. Below are
the list of requirements:


|

Requirements
~~~~~~~~~~~~
The following libraries/modules are required to run TRANSIT:


+ `Python 2.7 <http://www.python.org>`_
+ `Numpy 1.6.1+ <http://www.numpy.org/>`_
+ `Scipy 0.14.0+ <http://www.scipy.org/>`_
+ `matplotlib 1.1.1+ <http://matplotlib.org/users/installing.html>`_
+ `wxpython 2.8.0+ <http://www.wxpython.org/>`_ (for Mac OSX, use the **cocoa** version of wxPython)
+ `PIL (Python Imaging Library) <http://www.pythonware.com/products/pil/>`_ or Pillow.


Generally, these requirements are install using the appropriate
methods for your operating system, i.e. apt-get or yum for unix
machines, pip or easy_install for OSX, or binary installers on
Windows. Below more detailed instructions are provided.

|

Detailed Instructions: Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the requirements are available in default package sources in
most Linux distributions. The following commands will install python,
numpy, scipy, matplotlib on the Ubuntu or Fedora Linux distributions:

::

    
    #Ubuntu:
    sudo apt-get install python python-numpy python-scipy python-matplotlib python-wxgtk3.0
    
    #Fedora:
    sudo yum install python numpy scipy python-matplotlib python-wxgtk3.0


The final requirement left to install is Pillow. First you need
install pip which simplifies the process of installing certain python
modules like Pillow:


::

    
    #Ubuntu:
    sudo apt-get install pip
    
    #Fedora:
    sudo yum install pip


Next, using pip you must have a clean installation of Pillow, and the
desired libraries. You can achieve this through the following
commands:

::

    
    #Ubuntu:
    pip uninstall pillow
    pip uninstall Pillow
    sudo apt-get install libjpeg-dev zlib1g-dev
    pip install -I Pillow
    
    #Fedora:
    pip uninstall pillow
    pip uninstall Pillow
    sudo yum install install libjpeg-dev zlib1g-dev
    pip install -I Pillow


Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_. Download the source files:

::

    
    `http://sourceforge.net/projects/bio-bwa/files/`_


Extract the files:

::

    
    tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::

    
    cd bwa-0.7.12
    make



|

Detailed Instructions: OSX
~~~~~~~~~~~~~~~~~~~~~~~~~~
First, download and install the latest Python 2.7.x installation file from the official python website:


    
    `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


Next make sure you have pip installed. Pip can be installed through easy_install, which should come with OSX:

::

    
    sudo easy_install pip


Next install numpy, scipy, and matplotlib and pillow using pip:

::

    
    sudo pip install numpy
    sudo pip install scipy
    sudo pip install matplotlib
    sudo pip install pillow


Download and install the OSX binary of wxpython (cocoa version) for python 2.7:

::

    
    `http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg`_

Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_ . Download the source files:

::

    
    `http://sourceforge.net/projects/bio-bwa/files/`_


Extract the files:

::

    
    tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::

    
    cd bwa-0.7.12
    make




|

Detailed Instructions: Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, download and install the latest Python 2.7.x installation file
from the official python website:


::

    
    `http://www.python.org/downloads/`_


Next, you will need to install pip. If you are using python 2.7.9+
then pip will come pre-installed and included in the default script
directory (i.e. C:\Python27\Scripts ). If you are using python 2.7.8
or older, you will need to manually install pip by downloading and
running the `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`_ script:


::

    
    python.exe get-pip.py


Make sure that "wheel" is installed. This is necessary to allow you to
install .whl (wheel) files:

::

    
    pip.exe install wheel


Download the .whl files for all the requirements (Note: Make sure to
choose the files that match your Windows version i.e. 32/64 bit)

  + `numpy-1.9.2+mkl-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/numpy-1.9.2+mkl-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/numpy-1.9.2+mkl-cp27-none-win32.whl>`_


  + `scipy-0.15.1-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/scipy-0.15.1-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/scipy-0.15.1-cp27-none-win32.whl>`_


  + `matplotlib-1.4.3-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/matplotlib-1.4.3-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/matplotlib-1.4.3-cp27-none-win32.whl>`_


  + `Pillow-2.8.2-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/Pillow-2.8.2-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/Pillow-2.8.2-cp27-none-win32.whl>`_


  + `wxPython-3.0.2.0-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win32.whl>`_


  + `wxPython_common-3.0.2.0-py2-none-any.whl <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_






Source: These files were obtained from the `Unofficial Windows Binaries for Python Extension Packages by Christoph Gohlke, Laboratory for Fluorescence Dynamics, University of California, Irvine. <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_


Finally, install the files using pip:

::

    
    pip.exe install numpy-1.9.2+mkl-cp27-none-win_amd64.whl
    pip.exe install scipy-0.15.1-cp27-none-win_amd64.whl
    pip.exe install matplotlib-1.4.3-cp27-none-win_amd64.whl
    pip.exe install Pillow-2.8.1-cp27-none-win_amd64.whl
    pip.exe install wxPython-3.0.2.0-cp27-none-win_amd64.whl
    pip.exe install wxPython_common-3.0.2.0-py2-none-any.whl


Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_. We provide a windows executable (.exe) for Windows 64 bit:

`bwa-0.7.12_windows.zip <http://saclab.tamu.edu/essentiality/transit/bwa-0.7.12_windows.zip>`_






|

Running TRANSIT
---------------


|

GUI Mode
~~~~~~~~
To run TRANSIT in GUI mode (should be the same on Linux, Windows and MacOS), from the command line run:

::

    
    python PATH/src/transit.py

where PATH is the path to the TRANSIT installation directory. You might be able to double-click on icon for transit.py, if your OS associates .py files with python and automatically runs them. Note, because TRANSIT has a graphical user interface, if you are trying to run TRANSIT across a network, for example, running on a unix server but displaying on a desktop machine, you will probably need to use 'ssh -Y' and a local X11 client (like Xming or Cygwin/X on PCs).


|

Command line Mode
~~~~~~~~~~~~~~~~~
TRANSIT can also be run from the command line, without the GUI interface. This is convenient if you want to run many analyses in batch, as you can write a script that automatically runs that automatically runs TRANSIT from the command line. TRANSIT expects the user to specify which analysis method they wish to run. The user can choose from "gumbel", "hmm", or "resampling". By choosing a method, and adding the "-h" flag, you will get a list of all the necessary parameters and optional flags for the chosen method:

::

    python PATH/src/transit.py gumbel -h




|

Gumbel
``````

To run the Gumbel analysis from the command line, type "python PATH/src/transit.py gumbel" followed by the following arguments:


+----------------+----------------+----------------+----------------+----------------+
| Argument       | Type           | Description    | Default        | Example        |
+================+================+================+================+================+
| annotation     | Required       | Path to        |                | genomes/H37Rv. |
|                |                | annotation     |                | prot\_table    |
|                |                | file in        |                |                |
|                |                | .prot\_table   |                |                |
|                |                | format         |                |                |
+----------------+----------------+----------------+----------------+----------------+
| control\_files | Required       | Comma-separate |                | data/glycerol\ |
|                |                | d              |                | _reads\_rep1.w |
|                |                | list of paths  |                | ig,data/glycer |
|                |                | to the \*.wig  |                | ol\_reads\_rep |
|                |                | replicate      |                | 2.wig          |
|                |                | datasets       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| output\_file   | Required       | Name of the    |                | results/gumbel |
|                |                | output file    |                | \_glycerol.dat |
|                |                | with the       |                |                |
|                |                | results.       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -s SAMPLES     | Optional       | Number of      | 10000          | -s 20000       |
|                |                | samples to     |                |                |
|                |                | take.          |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -m MINREAD     | Optional       | Smallest       | 1              | -m 2           |
|                |                | read-count     |                |                |
|                |                | considered to  |                |                |
|                |                | be an          |                |                |
|                |                | insertion.     |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -b BURNIN      | Optional       | Burn in        | 500            | -b 100         |
|                |                | period, Skips  |                |                |
|                |                | this number of |                |                |
|                |                | samples before |                |                |
|                |                | getting        |                |                |
|                |                | estimates. See |                |                |
|                |                | documentation. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -t TRIM        | Optional       | Number of      | 1              | -t 2           |
|                |                | samples to     |                |                |
|                |                | trim. See      |                |                |
|                |                | documentation. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -r REP         | Optional       | How to handle  | Sum            | -r Mean        |
|                |                | replicates     |                |                |
|                |                | read-counts:   |                |                |
|                |                | 'Sum' or       |                |                |
|                |                | 'Mean'.        |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iN IGNOREN    | Optional       | Ignore TAs     | 5              | -iN 0          |
|                |                | occuring at X% |                |                |
|                |                | of the N       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iC IGNOREC    | Optional       | Ignore TAs     | 5              | -iC 10         |
|                |                | occuring at X% |                |                |
|                |                | of the C       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+



::

    python PATH/src/transit.py gumbel genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig test_console_gumbel.dat -s 20000 -b 1000




|

Tn5 Gaps
````````

To run the Tn5 Gaps analysis from the command line, type "python
PATH/src/transit.py tn5gaps" followed by the following arguments:

Argument Type Description Default Example annotation Required Path to
annotation file in .prot_table format genomes/Salmonella-
Ty2.prot_table control_files Required Comma-separated list of paths to
the \*.wig replicate datasets
data/salmonella_2122_rep1.wig,data/salmonella_2122_rep2.wig
output_file Required Name of the output file with the results.
results/test_console_tn5gaps.dat -m MINREAD Optional Smallest read-
count considered to be an insertion. 1 -m 2 -r REP Optional How to
handle replicates read-counts: 'Sum' or 'Mean'. Sum -r Sum

Example Tn5 Gaps command:

::

    python PATH/src/transit.py tn5gaps genomes/Salmonella-Ty2.prot_table data/salmonella_2122_rep1.wig,data/salmonella_2122_rep2.wig results/test_console_tn5gaps.dat -m 2 -r Sum





Example HMM command:

::

    python PATH/src/transit.py hmm genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig test_console_hmm.dat -r Sum


| 

Resampling
``````````

To run the Resampling analysis from the command line, type "python
PATH/src/transit.py resampling" followed by the following arguments:

+----------------+----------------+----------------+----------------+----------------+
| Argument       | Type           | Description    | Default        | Example        |
+================+================+================+================+================+
| annotation     | Required       | Path to        |                | genomes/H37Rv. |
|                |                | annotation     |                | prot\_table    |
|                |                | file in        |                |                |
|                |                | .prot\_table   |                |                |
|                |                | format         |                |                |
+----------------+----------------+----------------+----------------+----------------+
| control\_files | Required       | Comma-separate |                | data/glycerol\ |
|                |                | d              |                | _reads\_rep1.w |
|                |                | list of paths  |                | ig,data/glycer |
|                |                | to the \*.wig  |                | ol\_reads\_rep |
|                |                | replicate      |                | 2.wig          |
|                |                | datasets for   |                |                |
|                |                | the control    |                |                |
|                |                | condition      |                |                |
+----------------+----------------+----------------+----------------+----------------+
| exp\_files     | Required       | Comma-separate |                | data/cholester |
|                |                | d              |                | ol\_reads\_rep |
|                |                | list of paths  |                | 1.wig,data/cho |
|                |                | to the \*.wig  |                | lesterol\_read |
|                |                | replicate      |                | s\_rep2.wig    |
|                |                | datasets for   |                |                |
|                |                | the            |                |                |
|                |                | experimental   |                |                |
|                |                | condition      |                |                |
+----------------+----------------+----------------+----------------+----------------+
| output\_file   | Required       | Name of the    |                | results/gumbel |
|                |                | output file    |                | \_glycerol.dat |
|                |                | with the       |                |                |
|                |                | results.       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -s SAMPLES     | Optional       | Number of      | 10000          | -s 5000        |
|                |                | permutations   |                |                |
|                |                | performed.     |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -H             | Optional       | Creates        | Not set        | -H             |
|                |                | histograms of  |                |                |
|                |                | the            |                |                |
|                |                | permutations   |                |                |
|                |                | for all genes. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -a             | Optional       | Performs       | Not set        | -a             |
|                |                | adaptive       |                |                |
|                |                | appoximation   |                |                |
|                |                | to resampling. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -N             | Optional       | Select which   | nzmean         | -N nzmean      |
|                |                | normalizing    |                |                |
|                |                | procedure to   |                |                |
|                |                | use. Can       |                |                |
|                |                | choose between |                |                |
|                |                | 'TTR',         |                |                |
|                |                | 'nzmean',      |                |                |
|                |                | 'totreads',    |                |                |
|                |                | 'zinfnb',      |                |                |
|                |                | 'betageom',    |                |                |
|                |                | and 'nonorm'.  |                |                |
|                |                | See the        |                |                |
|                |                | parameters     |                |                |
|                |                | section for    |                |                |
|                |                | the            |                |                |
|                |                | `Re-sampling   |                |                |
|                |                | method <http:/ |                |                |
|                |                | /saclab.tamu.e |                |                |
|                |                | du/essentialit |                |                |
|                |                | y/transit/tran |                |                |
|                |                | sit.html#resam |                |                |
|                |                | pling>`__      |                |                |
|                |                | for a          |                |                |
|                |                | description of |                |                |
|                |                | these          |                |                |
|                |                | normalization  |                |                |
|                |                | options.       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iN IGNOREN    | Optional       | Ignore TAs     | 5              | -iN 0          |
|                |                | occuring at X% |                |                |
|                |                | of the N       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iC IGNOREC    | Optional       | Ignore TAs     | 5              | -iC 10         |
|                |                | occuring at X% |                |                |
|                |                | of the C       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+

Example Resampling command:

::

    python PATH/src/transit.py resampling genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig data/cholesterol_reads_rep1.wig,data/cholesterol_reads_rep2.wig,data/cholesterol_reads_rep3.wig test_console_resampling.dat -H -s 10000 -N nzmean

| 

--------------

|

Analysis Methods
----------------

|

Gumbel
~~~~~~

The Gumbel can be used to determine which genes are essential in a
single condition. It does a gene-by-gene analysis of the insertions at
TA sites with each gene, makes a call based on the longest consecutive
sequence of TA sites without insertion in the genes, calculates the
probabily of this using a Baysian model.

|

How does it work?
`````````````````

| For a formal description of how this method works, see our paper:
|  DeJesus, M.A., Zhang, Y.J., Sassettti, C.M., Rubin, E.J.,
  Sacchettini, J.C., and Ioerger, T.R. (2013).
| `Bayesian analysis of gene essentiality based on sequencing of transposon insertion libraries. <http://www.ncbi.nlm.nih.gov/pubmed/23361328>`_ *Bioinformatics*, 29(6):695-703.

|

Parameters
``````````

-  **Samples:** Gumbel uses Metropolis-Hastings (MH) to generate samples
   of posterior distributions. The default setting is to run the
   simulation for 10,000 iterations. This is usually enough to assure
   convergence of the sampler and to provide accurate estimates of
   posterior probabilities. Less iterations may work, but at the risk of
   lower accuracy.

-  **Burn-In:** Because the MH sampler many not have stabalized in the
   first few iterations, a "burn-in" period is defined. Samples obtained
   in this "burn-in" period are discarded, and do not count towards
   estimates.

-  **Trim:** The MH sampler produces Markov samples that are correlated.
   This parameter dictates how many samples must be attempted for every
   sampled obtained. Increasing this parameter will decrease the
   auto-correlation, at the cost of dramatically increasing the
   run-time. For most situations, this parameter should be left at the
   default of "1".

-  **Minimum Read:** The minimum read count that is considered a true
   read. Because the Gumbel method depends on determining gaps of TA
   sites lacking insertions, it may be suceptible to spurious reads
   (e.g. errors). The default value of 1 will consider all reads as true
   reads. A value of 2, for example, will ignore read counts of 1.

-  **Replicates:** Determines how to deal with replicates by averaging
   the read-counts or suming read counts accross datasets. This should
   not have an affect for the Gumbel method, aside from potentially
   affecting spurious reads.

|

Outputs and diagnostics
```````````````````````

The Gumbel method generates a tab-seperated output file at the location
chosen by the user. This file will automatically be loaded into the
Results Files section of the GUI, allowing you to display it as a table.
Alternatively, the file can be opened in a spreadsheet software like
Excel as a tab-separated file. The columns of the output file are
defined as follows:

+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Column Header   | Column Definition                                                                                                             |
+=================+===============================================================================================================================+
| ORF             | Gene ID.                                                                                                                      |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Name            | Name of the gene.                                                                                                             |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Description     | Gene description.                                                                                                             |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| k               | Number of Transposon Insertions Observed within the ORF.                                                                      |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| n               | Total Number of TA dinucleotides within the ORF.                                                                              |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| r               | Length of the Maximum Run of Non-Insertions observed.                                                                         |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| s               | Span of nucleotidies for the Maximum Run of Non-Insertions.                                                                   |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| zbar            | Posterior Probability of Essentiality.                                                                                        |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Call            | Essentiality call for the gene. Depends on FDR corrected thresholds. E=Essential U=Uncertain, NE=Non-Essential, S=too short   |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+

| 
|  Note: Technically, Bayesian models are used to calculate posterior
  probabilities, not p-values (which is a concept associated with the
  frequentist framework). However, we have implemented a method for
  computing the approximate false-discovery rate (FDR) that serves a
  similar purpose. This determines a threshold for significance on the
  posterior probabilities that is corrected for multiple tests. The
  actual thresholds used are reported in the headers of the output file
  (and are near 1 for essentials and near 0 for non-essentials). There
  can be many genes that score between the two thresholds (t1 < zbar <
  t2). This reflects intrinsic uncertainty associated with either low
  read counts, sparse insertion density, or small genes. If the
  insertion\_density is too low (< ~30%), the method may not work as
  well, and might indicate an unusually large number of Uncertain or
  Essential genes.

|

Run-time
````````

The Gumbel method takes on the order of 10 minutes for 10,000 samples.
Run-time is linearly proportional to the 'samples' parameter, or length
of MH sampling trajectory. Other notes: Gumbel can be run on multiple
replicates; replicate datasets will be automatically merged.





|

Tn Gaps
~~~~~~~

The Tn5 Gaps method can be used to determine which genes are essential
in a single condition for a TN5 dataset. It does an analysis of the
insertions at each site within the genome, makes a call for a given
gene based on the length of the most heavily overlapping run of sites
without insertions (gaps), calculates the probabily of this using a
Bayesian model.


|

How does it work?
`````````````````
This method of analysis is based on the original gumbel analysis
method. For a formal description of how the original method works, see
our paper:

Griffin, J.E., Gawronski, J.D., DeJesus, M.A., Ioerger, T.R., Akerley, B.J., Sassetti, C.M. (2011). 
`High-resolution phenotypic profiling defines genes essential for mycobacterial survival and cholesterol catabolism. <http://www.ncbi.nlm.nih.gov/pubmed/21980284>`_  *PLoS Pathogens*, 7(9):e1002251.

The Tn5 Gaps method modifies the original method in order to work on
TN5 datasets by analyzing non-insertion runs throughout the whole
genome, including non-coding regions, instead of within single genes.
In doing so, the expected maximum run length is calculated and a
p-value can be derived for every run. A gene is then classified by
using the p-value of the run with the largest number of nucleotides
overlapping with the gene.

This method was tested on a salmonella TN5 dataset presented in this
paper:

Langridge GC1, Phan MD, Turner DJ, Perkins TT, Parts L, Haase J,
Charles I, Maskell DJ, Peters SE, Dougan G, Wain J, Parkhill J, Turner
AK. (2009). `Simultaneous assay of every Salmonella Typhi gene using one million
transposon mutants. <http://www.ncbi.nlm.nih.gov/pubmed/19826075>`_ *Genome Res.* , 19(12):2308-16.

This data was downloaded from SRA (located `herei <http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP000051>`_) , and used to make
wig files (`base <http://orca1.tamu.edu/essentiality/transit/data/salmonella_base.wig>`_ and `bile <http://orca1.tamu.edu/essentiality/transit/data/salmonella_bile.wig>`_) and the following 4 baseline datasets
were merged to make a wig file: (IL2_2122_1,3,6,8). Our analysis
produced 415 genes with adjusted p-values less than 0.05, indicating
essentiality, and the analysis from the above paper produced 356
essential genes. Of these 356 essential genes, 344 overlap with the
output of our analysis.

|

Parameters
``````````


+ **Minimum Read:** The minimum read count that is considered a true read. Because the Gumbel method depends on determining gaps of TA sites lacking insertions, it may be suceptible to spurious reads (e.g. errors). The default value of 1 will consider all reads as true reads. A value of 2, for example, will ignore read counts of 1.


+ **Replicates:** Determines how to deal with replicates by averaging the read-counts or suming read counts accross datasets. This should not have an affect for the Gumbel method, aside from potentially affecting spurious reads.



|

Outputs and diagnostics
```````````````````````
The Tn5 Gaps method generates a tab-seperated output file at the
location chosen by the user. This file will automatically be loaded
into the Results Files section of the GUI, allowing you to display it
as a table. Alternatively, the file can be opened in a spreadsheet
software like Excel as a tab-separated file. The columns of the output
file are defined as follows:
Column Header Column Definition ORF Gene ID. Name Name of the gene.
Desc Gene description. k Number of Transposon Insertions Observed
within the ORF. n Total Number of TA dinucleotides within the ORF. r
Length of the Maximum Run of Non-Insertions observed. ovr The number
of nucleotides in the overlap with the longest run partially covering
the gene. lenovr The length of the above run with the largest overlap
with the gene. pval P-value calculated by the permutation test. padj
Adjusted p-value controlling for the FDR (Benjamini-Hochberg). call
Essentiality call for the gene. Depends on FDR corrected thresholds.
Essential or Non-Essential

Note: Technically, Bayesian models are used to calculate posterior
probabilities, not p-values (which is a concept associated with the
frequentist framework). However, we have implemented a method for
computing the approximate false-discovery rate (FDR) that serves a
similar purpose. This determines a threshold for significance on the
posterior probabilities that is corrected for multiple tests. The
actual thresholds used are reported in the headers of the output file
(and are near 1 for essentials and near 0 for non-essentials). There
can be many genes that score between the two thresholds (t1 < zbar <
t2). This reflects intrinsic uncertainty associated with either low
read counts, sparse insertion density, or small genes. If the
insertion_density is too low, the method may not work as well, and
might indicate an unusually large number of Uncertain or Essential
genes.

|

Run-time
````````
The Gumbel method takes on the order of 10 minutes.
Other notes: Gumbel can be run on multiple replicates; replicate
datasets will be automatically merged.







|

HMM
~~~

The HMM method can be used to determine the essentiality of the entire genome, as opposed to gene-level analysis of the other methods. It is capable of identifying regions that have unusually high or unusually low read counts (i.e. growth advantage or growth defect regions), in addition to the more common categories of essential and non-essential.

|

How does it work?
`````````````````

| For a formal description of how this method works, see our paper:
|  DeJesus, M.A., Ioerger, T.R. `A Hidden Markov Model for identifying essential and growth-defect regions in bacterial genomes from transposon insertion sequencing data. <http://www.ncbi.nlm.nih.gov/pubmed/24103077>`_ *BMC Bioinformatics.* 2013. 14:303

|

Parameters
``````````

The HMM method automatically estimates the necessary statistical
parameters from the datasets. You can change how the method handles
replicate datasets:

-  **Replicates:** Determines how the HMM deals with replicate datasets
   by either averaging the read-counts or suming read counts accross
   datasets. For regular datasets (i.e. mean-read count > 100) the
   recommended setting is to average read-counts together. For sparse
   datasets, it suming read-counts may produce more accurate results.

|

Output and Diagnostics
``````````````````````

| The HMM method outputs two files. The first file provides the most
  likely assignment of states for all the TA sites in the genome. Sites
  can belong to one of the following states: "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage). In
  addition, the output includes the probability of the particular site
  belonging to the given state. The columns of this file are defined as
  follows:

+------------+-----------------------------------------------------------------------------------------------------+
| Column #   | Column Definition                                                                                   |
+============+=====================================================================================================+
| 1          | Coordinate of TA site                                                                               |
+------------+-----------------------------------------------------------------------------------------------------+
| 2          | Observed Read Counts                                                                                |
+------------+-----------------------------------------------------------------------------------------------------+
| 3          | Probability for ES state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 4          | Probability for GD state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 5          | Probability for NE state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 6          | Probability for GA state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 7          | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+------------+-----------------------------------------------------------------------------------------------------+
| 8          | Gene(s) that share(s) the TA site.                                                                  |
+------------+-----------------------------------------------------------------------------------------------------+

| 
|  The second file provides a gene-level classification for all the
  genes in the genome. Genes are classified as "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage)
  depending on the number of sites within the gene that belong to those
  states.

+-------------------+-----------------------------------------------------------------------------------------------------+
| Column Header     | Column Definition                                                                                   |
+===================+=====================================================================================================+
| Orf               | Gene ID                                                                                             |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Name              | Gene Name                                                                                           |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Desc              | Gene Description                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------+
| N                 | Number of TA sites                                                                                  |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n0                | Number of sites labeled ES (Essential)                                                              |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n1                | Number of sites labeled GD (Growth-Defect)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n2                | Number of sites labeled NE (Non-Essential)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n3                | Number of sites labeled GA (Growth-Advantage)                                                       |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Avg. Insertions   | Mean insertion rate within the gene                                                                 |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Avg. Reads        | Mean read count within the gene                                                                     |
+-------------------+-----------------------------------------------------------------------------------------------------+
| State Call        | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+-------------------+-----------------------------------------------------------------------------------------------------+

| 
|  Note: Libraries that are too sparse (e.g. < 30%) or which contain
  very low read-counts may be problematic for the HMM method, causing it
  to label too many Growth-Defect genes.

|

Run-time
````````

| The HMM method takes less than 10 minutes to complete. The parameters
  of the method should not affect the running-time.

--------------

|

Re-sampling
~~~~~~~~~~~

The re-sampling method is a comparative analysis the allows that can be
used to determine conditional essentiality of genes. It is based on a
permutation test, and is capable of determining read-counts that are
significantly different across conditions.

|

How does it work?
`````````````````

This technique has yet to be formally published in the context of
differential essentiality analysis. Briefly, the read-counts at each
genes are determined for each replicate of each condition. The total
read-counts in condition A is substracted from the total read counts at
condition B, to obtain an observed difference in read counts. The TA
sites are then permuted for a given number of "samples". For each one of
these permutations, the difference is read-counts is determined. This
forms a null disttribution, from which a p-value is calculated for the
original, observed difference in read-counts.

|

Parameters
``````````

The resampling method is non-parametric, and therefore does not require
any parameters governing the distributions or the model. The following
parameters are available for the method:

-  **Samples:** The number of samples (permutations) to perform. The
   larger the number of samples, the more resolution the p-values
   calculated will have, at the expense of longer computation time. The
   re-sampling method runs on 10,000 samples by default.

-  **Output Histograms:**\ Determines whether to output .png images of
   the histograms obtained from resampling the difference in
   read-counts.

-  **Adaptive Resampling:** An optional "adaptive" version of resampling
   which accelerates the calculation by terminating early for genes
   which are likely not significant. This dramatically speeds up the
   computation at the cost of less accurate estimates for those genes
   that terminate early (i.e. deemed not significant). This option is
   OFF by default.

-  **Normalization Method:** Determines which normalization method to
   use when comparing datasets. Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences.

   -  **TTR**: Trimmed Total Reads (TTR), normalized by the total
      read-counts (like totreads), but trims top and bottom 5% of
      read-counts. This is the recommended normalization method for most
      cases as it has the beneffit of normalizing for difference in
      saturation in the context of resampling.
   -  **nzmean**: Normalizes datasets to have the same mean over the
      non-zero sites.
   -  **totreads**: Normalizes datasets by total read-counts, and scales
      them to have the same mean over all counts.
   -  **zinfnb**: Fits a zero-inflated negative binomial model, and then
      divides read-counts by the mean. The zero-inflated negative
      binomial model will treat some empty sites as belonging to the
      "true" negative binomial distribution responsible for read-counts
      while treating the others as "essential" (and thus not influencing
      its parameters).
   -  **quantile**: Normalizes datasets using the quantile normalization
      method described by `Bolstad et al.
      (2003) <http://www.ncbi.nlm.nih.gov/pubmed/12538238>`_. In this
      normalization procedure, datasets are sorted, an empirical
      distribution is estimated as the mean across the sorted datasets
      at each site, and then the original (unsorted) datasets are
      assigned values from the empirical distribution based on their
      quantiles.
   -  **betageom**: Normalizes the datasets to fit an "ideal" Geometric
      distribution with a variable probability parameter *p*. Specially
      useful for datasets that contain a large skew.
   -  **nonorm**: No normalization is performed.

|

Output and Diagnostics
``````````````````````

The re-sampling method ouputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+-----------------+-----------------------------------------------------------------+
| Column Header   | Column Definition                                               |
+=================+=================================================================+
| Orf             | Gene ID.                                                        |
+-----------------+-----------------------------------------------------------------+
| Name            | Name of the gene.                                               |
+-----------------+-----------------------------------------------------------------+
| Description     | Gene description.                                               |
+-----------------+-----------------------------------------------------------------+
| N               | Number of TA sites in the gene.                                 |
+-----------------+-----------------------------------------------------------------+
| TAs Hit         | Number of TA sites with at least one insertion.                 |
+-----------------+-----------------------------------------------------------------+
| Sum Rd 1        | Sum of read counts in condition 1.                              |
+-----------------+-----------------------------------------------------------------+
| Sum Rd 2        | Sum of read counts in condition 2.                              |
+-----------------+-----------------------------------------------------------------+
| Delta Rd        | Difference in the sum of read counts.                           |
+-----------------+-----------------------------------------------------------------+
| p-value         | P-value calculated by the permutation test.                     |
+-----------------+-----------------------------------------------------------------+
| p-adj.          | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+-----------------+-----------------------------------------------------------------+

| 

Run-time
````````

A typical run of the re-sampling method with 10,000 samples will take
around 45 minutes (with the histogram option ON). Using the adaptive
resampling option, the run-time is reduced to around 10 minutes.





