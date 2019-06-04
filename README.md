# Spatially Unbiased Infratentorial Template (SUIT) 

 SUIT is a Matlab toolbox dedicated to the analysis of imaging data of the human cerebellum. The toolbox contains a high-resolution atlas template of the human cerebellum and brainstem, based on the anatomy of 20 young healthy individuals. By using automated isolation and nonlinear normalization methods, a more accurate intersubject-alignment than current whole-brain methods can be achieved. The toolbox you to...

- automatically isolate cerebellar structures from the cerebral cortex based on an anatomical image.
- achieve accurate anatomical normalisation of cerebellar structures into atlas space using the Dartel algorithm.
- display the functional data on a surface-based flatmap representation.
- normalize focal cerebellar lesions for lesion-symptom mapping.
- use Voxel-based morphometry (VBM) to determine patterns of cerebellar degeneration or growth.
- use a probabilistic atlas of cerebellar anatomy to assign locations to different cerebellar lobules and deep cerebellar nuclei.
- use a collection of functional atlases based on task-based data and task-free connectivity to define regions of interest.
- improve normalization of the deep cerebellar nuclei using an ROI-driven normalization.

For more information on SUIT, please visit http://diedrichsenlab.org/imaging/suit.htm.

## Licence and Acknowledgements

The SUIT toolbox has been developed by J. Diedrichsen (joern.diedrichsen@googlemail.com),  C. Hernandez-Castillo, M. King, G. Prichard, N. Lally, T.Wiestler, J. Schlerf, D. Zhi, and many others. It is distributed under the [Creative Commons Attribution-NonCommercial 3.0 Unported License](http://creativecommons.org/licenses/by-nc/3.0/deed.en_US), meaning that it can be freely used for non-commercial purposes, as long as proper attribution in form of acknowledgments and links (for online use) or citations (in publications) are given. The relevant references are:

##### SUIT normalisation and template

- Diedrichsen, J. (2006). A spatially unbiased atlas template of the human cerebellum. *Neuroimage. 33(1)*, 127-138. 

##### Probabilistic atlas for cerebellar lobules

- Diedrichsen, J., Balsters, J. H., Flavell, J., Cussans, E., & Ramnani, N. (2009). A probabilistic atlas of the human cerebellum. *Neuroimage.46(1)*, 39-46.

##### Probabilistic atlas and normalisation for deep cerebellar nuclei

- Diedrichsen, J., Maderwald, S., Kuper, M., Thurling, M., Rabe, K., Gizewski, E. R., et al. (2011). Imaging the deep cerebellar nuclei: A probabilistic atlas and normalization procedure. *Neuroimage. 54(3)*, 1786-1794.

##### Surface-based representation and flatmap

- Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar data. *PLOSOne*. 

##### SUIT-N: Neonatal template, isolation and normalization

- Hernandez-Castillo, CR., Limperopoulos, C., & Diedrichsen, J. (2019). A representative template of the neonatal cerebellum. *Neuroimage. 184*, 450-454. 

##### Connectivity atlases

- Buckner, R. L., Krienen, F. M., Castellanos, A., Diaz, J. C., & Yeo, B. T. (2011). The organization of the human cerebellum estimated by intrinsic functional connectivity. *J Neurophysiol*, 106(5), 2322-2345.
- Ji, J. L., Spronk, M., Kulkarni, K., Repovs, G., Anticevic, A., & Cole, M. W. (2019). Mapping the human brain’s cortical-subcortical functional network organization. NeuroImage, 185, 35–57.

##### Task-based data sets

- King, M., Hernandez-Castillo, C.R., Poldrack, R. Ivry, R., Diedrichsen, J. (in press). A multi-domain task battery reveals functional boundaries in the human cerebellum. *Nature Neuroscience.* 
- Barch DM, Burgess GC, Harms MP, Petersen SE, Schlaggar BL, Corbetta M, Glasser MF, Curtiss S, Dixit S, Feldt C, Nolan D, Bryant E, Hartley T, Footer O, Bjork JM, Poldrack R, Smith S, Johansen-Berg H, Snyder AZ, Van Essen DC (2013) Function in the human connectome: task-fMRI and individual differences in behavior. *Neuroimage* 80:169-189.