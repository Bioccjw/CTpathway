# CTpathway
CTpathway is a CrossTalk-based pathway enrichment analysis method to identify disease risk pathways. Disregarding pathway crosstalk reduces the ability to detect risk pathways, especially those whose member genes are not dysregulated in the situations, such as cancer early stages or conditions, but that have essential contributions in transducing the signals during disease pathogenesis. Our new tool overcomes limitations in earlier methods by incorporating pathway crosstalk into our algorithm and calculating novel risk scores for genes in the pathways.

##	What can CTpathway do
* Calculating gene risk score
* Calculating pathway risk score
* Identifying significantly risk pathways
* Reducing pathway redundancy
* Visualizing pathway enrichment analysis results

### **Step 1. Install dependent R packages**

```
wget https://cran.r-project.org/src/contrib/Archive/fExtremes/fExtremes_3042.82.tar.gz
tar -zxvf fExtremes_3042.82.tar.gz
To install this packages, start "R" and enter:
install.packages("/fExtremes",repos=NULL,type="source")
```

### **Step 2. Install dependent Python packages**
```
pip3 install numpy==1.19.5
pip3 install scipy==1.5.4
pip3 install pandas==0.22.0
pip3 install rpy2==3.3.5
pip3 install cairosvg==2.5.2
pip3 install IPython==7.16.1
pip3 install chardet==5.0.0
pip3 install py4cytoscape==0.0.9
pip3 install networkx==2.3
pip3 install matplotlib==3.3.4
pip3 install markov_clustering==0.0.6.dev0
```

### **step3. Install and run cytoscape(3.8.2)**
```
(https://github.com/cytoscape/cytoscape/releases/3.8.2/)
tar -xvf cytoscape-unix-3.8.2.tar.gz
sh cytoscape.sh
```

## Usage
### Example

```
python3 CTpathway.py example/example_diff.txt symbol Both
```

Once the program has run successfully, a series of result files will appear in the "output" folder. The result files are as follows:   
	pathway_result.txt --Pathway enrichment result includes source, pathway risk score (PS), p_value, FDR, node_number and website for all candidate pathways.  
	significant_pathway_result.txt --Pathway enrichment result includes source, pathway risk score (PS), p_value, FDR, node_number and website for significant pathways.  
	risk_score_result.txt --Gene risk score (RS) for all genes.  
	cluster_info_\*.txt --Reducing pathway redundancy result for significant pathways includes significant pathway clustering information. \* ranges 0.1-0.9, which represents a cutoff of Jaccard similarity coefficient for shared genes among all significant pathway pairs used to construct a pathway similarity network. Then, MCL clustering algorithm is employed to absorb most redundancies into representative clusters.  
	bar.pdf --Pathway enrichment result visualization by bar plot shows top 20 significant pathways.  
	bubble.pdf --Pathway enrichment result visualization by bubble plot shows top 20 significant pathways.  
	network_\*.svg --Reducing pathway redundancy result visualization corresponding to cluster_info_\*.txt. CTpathway chooses the most significant (lowest FDR) pathway within each cluster to represent the cluster. To obtain a better visualization, CTpathway shows the top 20 non-redundant pathways or clusters with low FDR, if there are more than 20 clusters or pathways. For each cluster, the top 10 pathways with lower FDR are shown in the enrichment map if there are more than 10 pathways are within one cluster.  
	top20pathwayname_\*.svg --The names of representative pathway within each cluster corresponding to network_\*.svg.  

### Help Information
```
python3 CTpathway.py -h
```

##### **Positional parameters:**
```
file
	Gene differential expression profile, which should include gene id, one or both of log2FC (Fold Change) value and P value. Notice: input file without headline is commanded.
gene_id
	Gene id type corresponds to gene differential expression profile. Available options are "symbol" or "entrez_id".
data_type
	Data type corresponds to gene differential expression profile. Available options are "Log2_FC", "P_value" or "Both".
```

##### **Optional parameters:**
```
--pathway_database
	Optional. Select pathway database used in this analysis. Availale options are "KEGG", "Reactome", "PANTHER", "HumanCyc", "INOH", "NetPath", "PID", "WikiPathways" or "All". [default= "KEGG"].
--permutation_number
	Optional. Permutation number is used to identify significant pathways. [default= 1000]. Notice: the larger permutation number is, the longer time consumes.
--r_value
	Optional. Restart coefficient corresponds to multi-RWR (multiple random walk with restart) algorithm. Available options are "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8" or "0.9". [default= 0.7]. Notice: the larger value represents more dependency on differential expression and less dependency on Global pathway crosstalk map (GPCM), which is based on prior biological knowledge.
--significant_method
	Optional. Significant method is to choose P value or FDR to screen significant pathways. Availale options are "FDR" or "P_value". [default= "FDR"].
--significant_threshold
	Optional. Significant threshold is threshold of significant pathways. [default= 0.01].
```	
	
### Contact	
All comments, questions and suggestions: weijiang@nuaa.edu.cn		
