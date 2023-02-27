## GSEA exercises as part of the Transcriptomics Module (2022-2023).
### Conda environment installation (for exercises in R)
1. Clone this repository:
```
git clone https://github.com/bioinfo-lessons/gsea
```
2. Go to gsea folder:
```
cd gsea/
```
3. Install a new conda environment using the provided YAML:
```
mamba env create -f gsea.yaml
```

4. Activate the enviroment:
```
conda activate ISCIII_GSEA
```

5. Deactivate when finished:
```
conda deactivate
```

### GSEA installation (Linux)
1. Go to GSEA official [webpage](http://www.gsea-msigdb.org/gsea/downloads.jsp). Click on download GSEA_Linux_4.3.2.zip.

2. Save GSEA folder whenever you want (e.g. your home directory: `~/`). 

3. Uncompress the zip.

### GSEA server launch
1. Go to GSEA folder:

```
cd ~/GSEA_Linux_4.3.2/
```

2. Run the following command (java >= 11 required):

```
java --module-path=modules -Xmx4g @gsea.args --patch-module=jide.common=lib/jide-components-3.7.4.jar:lib/jide-dock-3.7.4.jar:lib/jide-grids-3.7.4.jar --module=org.gsea_msigdb.gsea/xapps.gsea.GSEA
```
