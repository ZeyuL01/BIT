10xGenomics PBMCs cell-type-specific TRs
========================================

We utilize SnapATAC2 to process single-cell ATAC-seq (scATAC-seq) data obtained from `10xGENOMICS <https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k>`_. . The detailed usage guide for SnapATAC2 can be found in the `SnapATAC2-manual <https://kzhang.org/SnapATAC2/tutorials/index.html>`_. We used **SnapATAC2 v2.6.0** for our analysis.

For improved efficiency, some steps require a GPU to accelerate model training.

.. code-block:: python

  import warnings
  warnings.filterwarnings("ignore")

  import snapatac2 as snap
  import anndata as ad
  import pandas as pd
  import scanpy as sc
  import scvi

  #Set seed
  scvi.settings.seed = 2
  snap.__version__

  print("i")

  reference = snap.read(snap.datasets.pbmc10k_multiome(), backed=None)
  atac=snap.read("./10XGENOMICS/pbmc.h5ad",backed=None)
  query = snap.pp.make_gene_matrix(atac, gene_anno=snap.genome.hg38)
  query.obs['cell_type'] = pd.NA
  data = ad.concat(
      [reference, query],
      join='inner',
      label='batch',
      keys=["reference", "query"],
      index_unique='_',
  )

  sc.pp.filter_genes(data, min_cells=5)
  sc.pp.highly_variable_genes(
      data,
      n_top_genes = 5000,
      flavor="seurat_v3",
      batch_key="batch",
      subset=True
  )

  scvi.model.SCVI.setup_anndata(data, batch_key="batch")
  vae = scvi.model.SCVI(
      data,
      n_layers=2,
      n_latent=30,
      gene_likelihood="nb",
      dispersion="gene-batch",
  )

  vae.train(max_epochs=2000, early_stopping=True)

  data.obs["celltype_scanvi"] = 'Unknown'
  ref_idx = data.obs['batch'] == "reference"
  data.obs["celltype_scanvi"][ref_idx] = data.obs['cell_type'][ref_idx]

  lvae = scvi.model.SCANVI.from_scvi_model(
      vae,
      adata=data,
      labels_key="celltype_scanvi",
      unlabeled_category="Unknown",
  )

  lvae.train(max_epochs=2000, n_samples_per_label=100)

  data.obs["C_scANVI"] = lvae.predict(data)
  data.obsm["X_scANVI"] = lvae.get_latent_representation(data)

  sc.pp.neighbors(data, use_rep="X_scANVI")
  sc.tl.umap(data)

  data.write("./10XGENOMICSpbmc10k.h5ad",compression="gzip")

  atac.obs['cell_type'] = data.obs.loc[atac.obs_names + '_query']['C_scANVI'].to_numpy()

  #cluster the similar cell types
  mapping_dict = {
    'CD4 Naive': 'CD4+ cell',
    'CD4 TCM': 'CD4+ cell',
    'CD4 TEM': 'CD4+ cell',
    'Treg':'CD4+ cell',
    'CD8 Naive': 'CD8+ cell',
    'CD8 TEM_1': 'CD8+ cell',
    'CD8 TEM_2': 'CD8+ cell',
    'Naive B': 'B cell',
    'Intermediate B': 'B cell',
    'Memory B':'B cell',
    'Plasma':'B cell',
    'NK': 'NK cell',
    'cDC': 'Dendritic cell',
    'pDC': 'Dendritic cell',
    'CD14 Mono':'Monocyte',
    'CD16 Mono':'Monocyte',
    'HSPC':'HSPC',
    'MAIT':'MAIT cell',
    'gdT':'gdT cell'
  }

  atac.write("./10XGENOMICS/pbmc10k_annotated.h5ad", compression="gzip")

  data.obs["Cell_Types"]=data.obs["cell_type"].map(mapping_dict)

  snap.tl.macs3(data, groupby='Cell_Types')
  peaks = snap.tl.merge_peaks(data.uns['macs3'], snap.genome.hg38)
  peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks['Peaks'])
  os.mkdir("./10XGENOMICS/csv/scATAC_Peaks_"+scvi.settings.seed)

  plt.rcParams['figure.dpi'] = 1000  # for displaying
  plt.rcParams['savefig.dpi'] = 1000  # for saving
  plt.rcParams['figure.figsize'] = [8, 8]  # width, height in inches

  snap.tl.umap(data,random_state=15)
  sc.pl.umap(data, color="Cell_Types",
          size=15,  # Increase point size for visibility
          alpha=0.9,  # Slightly transparent points to visualize density
          legend_fontsize=20,
          legend_fontweight='bold',
          frameon=True,  # Hide the frame for a cleaner look
          ncols=2,  # Organize legend into two columns
          show=False,
          save='umap_plot_Test.pdf'
      )

  gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)
  gene_matrix

  sc.pp.filter_genes(gene_matrix, min_cells= 5)
  sc.pp.normalize_total(gene_matrix)
  sc.pp.log1p(gene_matrix)

  sc.external.pp.magic(gene_matrix, solver="approximate")
  gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]

  gene_matrix.write("pbmc10k_gene_mat.h5ad", compression='gzip')
  sc.set_figure_params(scanpy=True, dpi=1000,dpi_save=1000,fontsize=24,figsize=[10,10])
  marker_genes = [] #Type marker genes here to color the umap plot
  gene_matrix=snap.read("pbmc10k_gene_mat.h5ad",backed=None)

  for i in range(len(marker_genes)):
      sc.pl.umap(gene_matrix, use_raw=False, color=marker_genes[i],
          size=15,  # Increase point size for visibility
          alpha=0.9,  # Slightly transparent points to visualize density
          frameon=True,  # Hide the frame for a cleaner look
          ncols=5,  # Organize legend into two columns
          show=False,
          save='umap_plot_Test_UMAP_PBMCs_'+marker_genes[i]+'.pdf',
          color_map="plasma"
          )


The previous pipeline can generate the below UMAP plots:

Colored by cell types:

.. image:: ../images/Examples/singlecell/Pic1.png

Colored by PAX5 gene activity:

.. image:: ../images/Examples/singlecell/Pic2.png

We also get the marker peaks of each cell type:

.. code-block:: python

  marker_peaks=snap.tl.marker_regions(peak_mat, groupby='Cell_Types', pvalue=0.01)
  for keys in marker_peaks.keys():
          elements=marker_peaks[keys]
          chromosomes = []
          starts = []
          ends = []
          for element in elements:
              # Split each element into chromosome, start, and end
              chromosome, positions = element.split(':')
              start, end = positions.split('-')
              # Append the results to the corresponding lists
              chromosomes.append(chromosome)
              starts.append(start)
              ends.append(end)
          df = pd.DataFrame({'Chrom': chromosomes,'Start': starts,'End': ends})
          df.to_csv("./10XGENOMICS/csv/scATAC_Peaks_"+scvi.settings.seed+"/"+keys.replace(" ","_")+".csv",index=False)
      print("./10XGENOMICS/csv/scATAC_Peaks_"+scvi.settings.seed+" Done")

cell-type-specific marker regions:

.. image:: ../images/Examples/singlecell/Pic3.png

We run BIT on each of the region set:

.. code-block:: R
  work_dir<-"./10XGENOMICS/csv/"
  work_files<-list.files(work_dir)
  output_path<-"./10XGENOMICS/bit/"
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  for(i in seq_along(work_files)){
    BIT(paste0(work_dir,work_files[i]), output_path=output_dir, format="csv", bin_width=1000, genome="hg38")
  }


We ause the SnapATAC2 internal implemented motif enrichment analysis method to derive the SnapATAC2 motif enrichment results:

.. code-block:: python

  os.makedirs('./10XGENOMICS/motifs/', exist_ok=True)
  motifs = snap.tl.motif_enrichment(
      motifs=snap.datasets.cis_bp(unique=True),
      regions=marker_peaks,
      genome_fasta=snap.genome.hg38,
  )
  for keys in motifs.keys():
                  elements=motifs[keys]
                  df=elements.to_pandas()
                  df=df.sort_values(by="adjusted p-value",ascending=True)
                  df.to_csv('./10XGENOMICS/motifs/pbmc10k_'+str(i)+"/"+keys+'_motifs.csv',index=False)


We plot the top 10 TRs in each cell type along with the 95% credible intervals:

