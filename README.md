

# gene_cluster_networks_and_genetic_dereplication
Repository for "Uncovering secondary metabolite evolution and biosynthesis using gene cluster networks and genetic dereplication"


## Important

You need to load the data provided under the data availability statement into a mysql database and write your credentials/ paths accordingly into config.txt

## Creating a virtual environment

```bash
virtualenv /home/seth/virtualenvs/secMetPipeline -p /usr/bin/python3.5

pip install biopython pandas ete3 mysqlclient
```

# R packages

* ape
* RColorBrewer
* igraph
* parallel
* ggraph
* ggplot2
* ggtree
* ggstance
* optparse
* Gviz (bioconductor)

# Test run

```bash
python3 MAIN.py -o nigri_set.txt -bibase publication_nigri_biblast_corrected -biFinal publication_smurf_bidir_hits_nigri_test -t smBasicTree.nwk -l run.log -od nigri_test
```

# Data

**Link to data is included in the manuscript under Data availability**

## JGI links

* [Aspergillus nidulans](https://genome.jgi.doe.gov/Aspnid1/Aspnid1.home.html)
* [Aspergillus flavus](https://genome.jgi.doe.gov/Aspfl1/Aspfl1.home.html)
* [Aspergillus saccharolyticus](https://genome.jgi.doe.gov/Aspsac1/Aspsac1.home.html)
* [Aspergillus luchuensis CBS 106.47](https://genome.jgi.doe.gov/Aspfo1/Aspfo1.home.html)
* [Aspergillus neoniger](https://genome.jgi.doe.gov/Aspneo1/Aspneo1.home.html)
* [Aspergillus vadensis](https://genome.jgi.doe.gov/Aspvad1/Aspvad1.home.html)
* [Aspergillus costaricaensis](https://genome.jgi.doe.gov/Aspcos1/Aspcos1.home.html)
* [Aspergillus tubingensis](https://genome.jgi.doe.gov/Asptu1/Asptu1.home.html)
* [Aspergillus brunneoviolaceus](https://genome.jgi.doe.gov/Aspbru1/Aspbru1.home.html)
* [Aspergillus aculeatinus](https://genome.jgi.doe.gov/Aspacu1/Aspacu1.home.html)
* [Aspergillus fijiensis](https://genome.jgi.doe.gov/Aspfij1/Aspfij1.home.html)
* [Aspergillus fumigatus Af293](https://genome.jgi.doe.gov/Aspfu1/Aspfu1.home.html)
* [Aspergillus niger NRRL3](https://genome.jgi.doe.gov/Aspni\_NRRL3\_1/Aspni\_NRRL3\_1.home.html)
* [Aspergillus niger ATCC 1015](https://genome.jgi.doe.gov/Aspni7/Aspni7.home.html)
* [Aspergillus phoenicis](https://genome.jgi.doe.gov/Aspph1/Aspph1.home.html)
* [Aspergillus niger ATCC 13496](https://genome.jgi.doe.gov/Aspni\_bvT\_1/Aspni\_bvT\_1.home.html)
* [Aspergillus brasiliensis](https://genome.jgi.doe.gov/Aspbr1/Aspbr1.home.html)
* [Aspergillus lacticoffeatus](https://genome.jgi.doe.gov/Asplac1/Asplac1.home.html)
* [Aspergillus eucalypticola](https://genome.jgi.doe.gov/Aspeuc1/Aspeuc1.home.html)
* [Aspergillus homomorphus](https://genome.jgi.doe.gov/Asphom1/Asphom1.home.html)
* [Aspergillus indologenus](https://genome.jgi.doe.gov/Aspind1/Aspind1.home.html)
* [Aspergillus oryzae](https://genome.jgi.doe.gov/Aspor1/Aspor1.home.html)
* [Penicillium chrysogenum](https://genome.jgi.doe.gov/Pench1/Pench1.home.html)
* [Aspergillus violaceofuscus](https://genome.jgi.doe.gov/Aspvio1/Aspvio1.home.html)
* [Aspergillus japonicus](https://genome.jgi.doe.gov/Aspjap1/Aspjap1.home.html)
* [Aspergillus uvarum](https://genome.jgi.doe.gov/Aspuva1/Aspuva1.home.html)
* [Aspergillus sclerotiicarbonarius](https://genome.jgi.doe.gov/Aspscle1/Aspscle1.home.html)
* [Aspergillus sclerotioniger](https://genome.jgi.doe.gov/Aspscl1/Aspscl1.home.html)
* [Aspergillus carbonarius](https://genome.jgi.doe.gov/Aspca3/Aspca3.home.html)
* [Aspergillus ellipticus](https://genome.jgi.doe.gov/Aspell1/Aspell1.home.html)
* [Aspergillus heteromorphus](https://genome.jgi.doe.gov/Asphet1/Asphet1.home.html)
* [Aspergillus welwitschiae](https://genome.jgi.doe.gov/Aspwel1/Aspwel1.home.html)
* [Aspergillus niger CBS 513.88](https://genome.jgi.doe.gov/Aspni\_DSM\_1/Aspni\_DSM\_1.home.html)
* [Aspergillus piperis](https://genome.jgi.doe.gov/Asppip1/Asppip1.home.html)
* [Aspergillus ibericus](https://genome.jgi.doe.gov/Aspibe1/Aspibe1.home.html)
* [Aspergillus luchuensis IFO 4308](https://genome.jgi.doe.gov/Aspka1\_1/Aspka1\_1.home.html)
* [Aspergillus aculeatus](https://genome.jgi.doe.gov/Aspac1/Aspac1.home.html)

Repository for "Uncovering secondary metabolite evolution and biosynthesis using gene cluster networks and genetic dereplication". Code and data are only intended for review and use is only permitted for reviewers. After successfull submission, a license will be added.

